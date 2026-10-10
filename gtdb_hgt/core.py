"""Provisional, interpretable similarity scores; never HGT probabilities."""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import math
from pathlib import Path
from statistics import median

from .io import BLAST_COLUMNS, canonical_accession, open_text, tax_key
from .quality import alignment_quality, best_bit_hit, hit_quality


def blast_rows(path: Path):
    with open_text(path) as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "qseqid":
                continue
            if len(parts) != len(BLAST_COLUMNS):
                raise ValueError(f"{path}:{line_no}: expected 14 DIAMOND columns, found {len(parts)}")
            row = dict(zip(BLAST_COLUMNS, parts))
            for key in ("bitscore", "evalue", "pident"):
                row[key] = float(row[key])
                if not math.isfinite(row[key]) or row[key] < 0:
                    raise ValueError(f"Invalid {key} at {path}:{line_no}")
            for key in ("qlen", "slen", "qstart", "qend", "sstart", "send"):
                row[key] = int(row[key])
            if row["qlen"] <= 0 or row["slen"] <= 0:
                raise ValueError(f"Invalid sequence length at {path}:{line_no}")
            row["qcov"] = (abs(row["qend"] - row["qstart"]) + 1) / row["qlen"]
            row["scov"] = (abs(row["send"] - row["sstart"]) + 1) / row["slen"]
            yield row


def self_scores(path: Path) -> dict[str, float]:
    result = {}
    for r in blast_rows(path):
        if r["qseqid"] == r["sseqid"]:
            result[r["qseqid"]] = max(result.get(r["qseqid"], 0), r["bitscore"])
    if not result:
        raise ValueError("No exact query-ID self alignments in self-hits file")
    return result


def collect_hits(path, taxa, query_tax, close_rank, donor_rank, exclude_rank,
                 query_genome=None, subject_map=None, evalue=1e-10,
                 min_query_cover=0.5, min_subject_cover=0.5, disk_store=None):
    """All eligible hits, collapsed to the best protein per GTDB species.

    Reference taxonomies must cover actual genomes in the search DB. No hit cap
    is imposed here. Retaining one species vote prevents copy/strain inflation.
    Memory scales with distinct query/species pairs, not raw protein hit rows.
    """
    subject_map = subject_map or {}
    close_key = tax_key(query_tax, close_rank)
    if close_key is None:
        raise ValueError(f"Query taxonomy is incomplete at {close_rank}")
    exclude_key = tax_key(query_tax, exclude_rank) if exclude_rank != "none" else None
    if exclude_rank != "none" and exclude_key is None:
        raise ValueError(f"Query taxonomy is incomplete at self-exclude rank {exclude_rank}")
    best = disk_store if disk_store is not None else defaultdict(dict)
    counts = defaultdict(Counter)
    audit = Counter()
    query_genome = canonical_accession(query_genome) if query_genome else None
    for row in blast_rows(path):
        q, sid = row["qseqid"], row["sseqid"]
        counts[q]["raw_hits"] += 1
        audit["raw_hits"] += 1
        if row["evalue"] > evalue or row["qcov"] < min_query_cover or row["scov"] < min_subject_cover:
            counts[q]["filtered_alignment"] += 1
            continue
        genome = subject_map.get(sid)
        if genome is None:
            genome = canonical_accession(sid.split("::", 1)[0]) if "::" in sid else None
        if genome is None or genome not in taxa:
            counts[q]["unknown_taxonomy"] += 1
            audit["unknown_taxonomy"] += 1
            continue
        tax = taxa[genome]
        if genome == query_genome or (exclude_key is not None and tax_key(tax, exclude_rank) == exclude_key):
            counts[q]["excluded_self"] += 1
            continue
        species = tax_key(tax, "species")
        donor = tax_key(tax, donor_rank)
        ref_close = tax_key(tax, close_rank)
        if species is None or donor is None or ref_close is None:
            counts[q]["incomplete_taxonomy"] += 1
            continue
        if row["bitscore"] <= 0:
            continue
        # Validate the quality inputs, but do not filter away comparison hits.
        alignment_quality(row["pident"], row["qcov"], row["scov"])
        counts[q]["eligible_protein_hits"] += 1
        row = {**row, "genome_id": genome, "tax": tax,
               "group": "close" if ref_close == close_key else "distal", "donor_key": donor}
        if disk_store is not None:
            disk_store.put(q, species, row)
        else:
            old = best[q].get(species)
            if old is None or (row["bitscore"], row["qcov"], row["scov"], sid) > (old["bitscore"], old["qcov"], old["scov"], old["sseqid"]):
                best[q][species] = row
    if disk_store is not None:
        disk_store.flush()
    return best, counts, dict(audit)


def score_genes(genes, best, counts, selfs, taxa, query_tax, close_rank="genus",
                donor_rank="genus", exclude_rank="species", query_genome=None,
                min_close_genomes=3, min_donor_species=2, missing_close="zero", min_hit_quality=0):
    if not math.isfinite(min_hit_quality) or not 0 <= min_hit_quality <= 1:
        raise ValueError("min-hit-quality must be between 0 and 1")
    close_key = tax_key(query_tax, close_rank)
    excluded = tax_key(query_tax, exclude_rank) if exclude_rank != "none" else None
    query_genome = canonical_accession(query_genome) if query_genome else None
    close_genomes = [g for g, t in taxa.items() if tax_key(t, close_rank) == close_key
                     and g != query_genome and (excluded is None or tax_key(t, exclude_rank) != excluded)]
    # Reference availability is distinct species, not copies of one strain.
    close_species_available = len({tax_key(taxa[g], "species") for g in close_genomes if tax_key(taxa[g], "species")})
    scored, summaries = [], []
    for gene in genes:
        q = gene["orf_id"]
        hits = list(best.get(q, {}).values())
        c = [h for h in hits if h["group"] == "close"]
        d = defaultdict(list)
        for h in hits:
            if h["group"] == "distal":
                d[h["donor_key"]].append(h)
        for group, group_hits in d.items():
            values = [h["bitscore"] for h in group_hits]
            summaries.append({"orf_id": q, "group": "distal", "taxon": ";".join(group),
                              "species_count": len(values), "median_bitscore": median(values),
                              "min_bitscore": min(values), "max_bitscore": max(values)})
        if c:
            values = [h["bitscore"] for h in c]
            summaries.append({"orf_id": q, "group": "close", "taxon": ";".join(close_key),
                              "species_count": len(values), "median_bitscore": median(values),
                              "min_bitscore": min(values), "max_bitscore": max(values)})
        # Best distal taxonomic bin by median; all eligible species in each bin
        # contribute, rather than the top N hits. Sparse bins remain visible.
        supported = [key for key in d if len(d[key]) >= min_donor_species]
        donor = max(supported or list(d), key=lambda key: (median(h["bitscore"] for h in d[key]), len(d[key]), key)) if d else None
        dh = d[donor] if donor else []
        top = best_bit_hit(hits)
        top_quality = hit_quality(top)
        donor_quality = hit_quality(best_bit_hit(dh))
        db = median(h["bitscore"] for h in dh) if dh else None
        cb = median(h["bitscore"] for h in c) if c else None
        normalizer = selfs.get(q)
        flags = []
        eligible = True
        if not normalizer or normalizer <= 0:
            flags.append("missing_self_score")
            eligible = False
        if close_species_available < min_close_genomes:
            flags.append("insufficient_close_reference_species")
            eligible = False
        if not dh:
            flags.append("no_distal_hit")
        elif len(dh) < min_donor_species:
            flags.append("sparse_distal_support")
            eligible = False
        if cb is None:
            flags.append("close_not_detected_above_filters")
            if missing_close == "na":
                eligible = False
        if not hits:
            flags.append("no_usable_hit")
            eligible = False
        if counts[q]["unknown_taxonomy"]:
            flags.append("some_hits_unmapped_taxonomy")
        # Distal absence is a valid non-positive comparison when close hits
        # exist. Neither side observed is always NA. Missing close zero is an
        # explicit censoring convention, never evidence of biological absence.
        score = None
        if normalizer and hits and (cb is not None or missing_close == "zero"):
            score = ((db or 0.0) - (cb or 0.0)) / normalizer
            if abs(score) > 1:
                flags.append("normalized_score_outside_unit_range")
        raw_score = score
        quality_status = "no_usable_hit" if not hits else "not_applied"
        if hits and min_hit_quality > 0:
            quality_status = "pass"
            if top_quality is None:
                quality_status = "unavailable"
                flags.append("missing_top_hit_quality")
            elif top_quality < min_hit_quality:
                quality_status = "low_top_hit_quality"
                flags.append(quality_status)
            elif score is not None and score > 0:
                if donor_quality is None:
                    quality_status = "unavailable"
                    flags.append("missing_donor_hit_quality")
                elif donor_quality < min_hit_quality:
                    quality_status = "low_donor_hit_quality"
                    flags.append(quality_status)
            if quality_status != "pass":
                score = None
                eligible = False
        scored.append({**gene, "self_bitscore": normalizer, "close_rank": close_rank,
                       "donor_rank": donor_rank, "close_reference_species": close_species_available,
                       "close_hit_species": len(c), "distal_hit_species": sum(len(h) for h in d.values()),
                       "close_median_bitscore": cb, "donor_median_bitscore": db,
                       "donor_species": len(dh), "donor_taxon": ";".join(donor) if donor else "",
                       "donor_label": donor[-1] if donor else "", "raw_orf_score": raw_score, "orf_score": score,
                       "top_hit_subject": top.get("sseqid", "") if top else "",
                       "top_hit_pident": top.get("pident") if top else None,
                       "top_hit_qcov": top.get("qcov") if top else None,
                       "top_hit_scov": top.get("scov") if top else None,
                       "top_hit_evalue": top.get("evalue") if top else None,
                       "top_hit_quality": top_quality, "donor_top_hit_quality": donor_quality,
                       "min_hit_quality": min_hit_quality, "quality_status": quality_status,
                       "eligible_for_window": eligible, "flags": ";".join(flags),
                       **{key: counts[q][key] for key in ("raw_hits", "eligible_protein_hits", "filtered_alignment", "unknown_taxonomy", "excluded_self")}})
    return scored, summaries


def smooth_and_regions(scored, window=5, min_valid=3, threshold=0.15,
                       max_gap_bp=10000, min_region_orfs=3):
    if window < 1 or window % 2 == 0:
        raise ValueError("Window size must be a positive odd number")
    if not 1 <= min_valid <= window:
        raise ValueError("min-valid must be between 1 and window")
    if min_region_orfs < 1:
        raise ValueError("min-region-orfs must be positive")
    # Split at contig boundaries and at large intergenic gaps. Circular origin
    # wrapping is not assumed. Edge windows are clipped and explicitly marked.
    blocks, block = [], []
    for row in sorted(scored, key=lambda r: (r["contig"], r["start"], r["end"], r["orf_id"])):
        if block and (row["contig"] != block[-1]["contig"] or row["start"] - block[-1]["end"] - 1 > max_gap_bp):
            blocks.append(block)
            block = []
        block.append(row)
    if block:
        blocks.append(block)
    radius = window // 2
    for block in blocks:
        for i, row in enumerate(block):
            neighborhood = block[max(0, i - radius):i + radius + 1]
            valid = [r for r in neighborhood if r["eligible_for_window"] and r["orf_score"] is not None]
            value = sum(r["orf_score"] for r in valid) / len(valid) if len(valid) >= min_valid else None
            positive = [r for r in valid if r["orf_score"] > 0 and r["donor_taxon"]]
            consensus = Counter(r["donor_taxon"] for r in positive)
            taxon, support = sorted(consensus.items(), key=lambda x: (-x[1], x[0]))[0] if consensus else ("", 0)
            row.update(window_score=value, window_orfs=len(neighborhood), window_valid_orfs=len(valid),
                       window_edge=len(neighborhood) < window, window_start=neighborhood[0]["start"],
                       window_end=max(r["end"] for r in neighborhood),
                       window_donor_taxon=taxon, window_donor_agreement=support / len(positive) if positive else None,
                       region_candidate=value is not None and value >= threshold)
    regions = []
    for block in blocks:
        run = []
        def emit():
            if len(run) < min_region_orfs:
                return
            signals = [r for r in run if r["eligible_for_window"] and r["orf_score"] is not None]
            donor_counts = Counter(r["donor_taxon"] for r in signals if r["orf_score"] > 0 and r["donor_taxon"])
            donor, count = sorted(donor_counts.items(), key=lambda x: (-x[1], x[0]))[0] if donor_counts else ("", 0)
            regions.append({"region_id": f"candidate_{len(regions) + 1}", "contig": run[0]["contig"],
                            "start": run[0]["start"], "end": max(r["end"] for r in run),
                            "n_orfs": len(run), "n_scored_orfs": len(signals),
                            "mean_window_score": sum(r["window_score"] for r in run) / len(run),
                            "mean_raw_score": sum(r["orf_score"] for r in signals) / len(signals) if signals else None,
                            "donor_taxon": donor, "donor_agreement": count / sum(donor_counts.values()) if donor_counts else None,
                            "orf_ids": ";".join(r["orf_id"] for r in run),
                            "boundary_note": "boundaries are smoothed candidate centers, not inferred transfer breakpoints"})
        for row in block:
            if row["region_candidate"]:
                run.append(row)
            else:
                emit()
                run = []
        emit()
    return scored, regions
