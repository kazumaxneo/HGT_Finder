"""One ORF, one vote: exploratory contig taxonomy, independent of HGT score."""
from collections import Counter, defaultdict
import html
import json
import math
from pathlib import Path
from statistics import median

from .io import RANKS, tax_key, tsv_write

DISPLAY_RANKS = ("phylum", "class", "order", "family", "genus")
RANK_LABELS = dict(phylum="門", **{"class": "綱"}, order="目", family="科", genus="属")

STATUS_LABELS = {
    "ambiguous": "曖昧：分類群が競合",
    "insufficient_support": "支持種数不足",
    "no_usable_hit": "評価不能：有効ヒットなし",
    "incomplete_taxonomy": "評価不能：表示階層の分類なし",
    "low_hit_quality": "低い配列一致：由来判定保留",
    "quality_unavailable": "評価不能：ヒット品質情報なし",
}


def _assign_bin(gene, bins, host_key, min_species, margin, method, has_hits=True):
    if method not in ("best-hit", "median"):
        raise ValueError("Sankey method must be best-hit or median")
    metric = "max" if method == "best-hit" else "median"
    # A sparse but strongest bin must not silently be replaced by a weaker bin.
    candidates = list(bins.items()) if method == "best-hit" else [
        (key, value) for key, value in bins.items() if value["species"] >= min_species]
    candidates.sort(key=lambda item: (-item[1][metric], -item[1]["species"], item[0]))
    winner = candidates[0] if candidates else None
    runner = candidates[1] if len(candidates) > 1 else None
    gap = (winner[1][metric] - runner[1][metric]) / winner[1][metric] if runner else None
    exact_best = [key for key, value in candidates if value[metric] == winner[1][metric]] if winner else []
    near_best = [key for key, value in candidates
                 if (winner[1][metric] - value[metric]) / winner[1][metric] <= margin] if winner else []
    if not has_hits:
        status = "no_usable_hit"
    elif not bins:
        status = "incomplete_taxonomy"
    elif not winner:
        status = "insufficient_support"
    elif method == "best-hit" and len(exact_best) > 1:
        status = "ambiguous"
    elif winner[1]["species"] < min_species:
        status = "insufficient_support"
    elif method == "median" and runner and gap <= margin:
        status = "ambiguous"
    else:
        status = "assigned"
    competitors = exact_best if method == "best-hit" else near_best
    host_competes = host_key in competitors
    ambiguity = ("host_competes" if host_competes else "other_competes") if status == "ambiguous" else ""
    assigned = winner[0] if status == "assigned" else None
    label = assigned[-1] if assigned else STATUS_LABELS[status]
    if status == "ambiguous":
        label = ("同点：宿主分類群も含む" if host_competes else "同点：宿主以外の分類群") if method == "best-hit" else (
            "曖昧：宿主分類群も競合" if host_competes else "曖昧：宿主以外が競合")
    host = bins.get(host_key)
    return {**gene, "assignment_status": status, "assignment_method": method,
            "taxon": ";".join(assigned) if assigned else "", "taxon_label": label,
            "host_taxon_match": assigned == host_key if assigned else None,
            "support_species": winner[1]["species"] if winner else 0,
            "median_bitscore": winner[1]["median"] if winner else None,
            "best_bitscore": winner[1]["max"] if winner else None,
            "ranking_bitscore": winner[1][metric] if winner else None,
            "leading_taxon": ";".join(winner[0]) if winner else "",
            "runner_up_taxon": ";".join(runner[0]) if runner else "",
            "relative_margin": gap,
            "near_tie": bool(status == "assigned" and runner and gap <= margin),
            "ambiguity_kind": ambiguity,
            "host_in_competition": host_competes,
            "host_in_near_best_set": host_key in near_best,
            "host_support_species": host["species"] if host else 0,
            "host_median_bitscore": host["median"] if host else None,
            "host_best_bitscore": host["max"] if host else None,
            "nearby_taxa": ";".join(key[-1] for key in near_best),
            "competing_taxa": ";".join(key[-1] for key in competitors)}


def assign_orfs(genes, best, query_tax, rank="genus", min_species=2, margin=.01, method="best-hit"):
    return assign_orfs_by_rank(genes, best, query_tax, (rank,), min_species, margin, method)[rank]


def assign_orfs_by_rank(genes, best, query_tax, ranks=DISPLAY_RANKS, min_species=2, margin=.01, method="best-hit"):
    """Compute each display rank from species-best hits, loading each ORF once."""
    hosts = {rank: tax_key(query_tax, rank) for rank in ranks}
    if not all(hosts.values()):
        raise ValueError("Query taxonomy is incomplete at a Sankey display rank")
    result = {rank: [] for rank in ranks}
    for gene in genes:
        hits = list(best.get(gene["orf_id"], {}).values())
        bins = {rank: defaultdict(list) for rank in ranks}
        for hit in hits:
            for rank in ranks:
                key = tax_key(hit["tax"], rank)
                if key:
                    bins[rank][key].append(hit["bitscore"])
        for rank in ranks:
            summaries = {key: {"median": median(values), "max": max(values), "species": len(values)}
                         for key, values in bins[rank].items()}
            result[rank].append(_assign_bin(gene, summaries, hosts[rank], min_species, margin, method, bool(hits)))
    return result


def _genus_summaries(path):
    """Yield one ORF's disjoint species bins at a time."""
    import csv
    q, bins = None, {}
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["orf_id"] != q:
                if q is not None:
                    yield q, bins
                q, bins = row["orf_id"], {}
            key = tuple(row["taxon"].split(";"))
            if len(key) != 6 or not all(v.startswith(prefix) and len(v)>3 for v, prefix in
                                      zip(key, ("d__", "p__", "c__", "o__", "f__", "g__"))) or key in bins:
                raise ValueError("Expected distinct full genus lineages in cached summaries")
            bins[key] = {"median": float(row["median_bitscore"]), "max": float(row["max_bitscore"]),
                         "species": int(row["species_count"])}
        if q is not None:
            yield q, bins


def assign_from_genus_summaries(genes, path, query_tax, min_species=2, margin=.01, method="best-hit"):
    return assign_cached_by_rank(genes, path, query_tax, ("genus",), min_species, margin, method)["genus"]


def assign_cached_by_rank(genes, path, query_tax, ranks=DISPLAY_RANKS, min_species=2, margin=.01, method="best-hit"):
    """Exact max and species counts at every rank from disjoint genus bins.

    Higher-rank medians cannot be recovered from genus medians. They stay NA
    in best-hit mode; cached median mode therefore supports genus only.
    """
    if method == "median" and any(rank != "genus" for rank in ranks):
        raise ValueError("Higher-rank median comparisons require original species hits")
    by_id = {g["orf_id"]: g for g in genes}
    hosts = {rank: tax_key(query_tax, rank) for rank in ranks}
    if not all(hosts.values()) or any(rank not in DISPLAY_RANKS for rank in ranks):
        raise ValueError("Invalid or incomplete cached Sankey rank")
    result = {rank: {} for rank in ranks}
    seen = set()
    for q, genus_bins in _genus_summaries(path):
        if q not in by_id or q in seen:
            raise ValueError("Unknown or noncontiguous query in genus summaries: " + q)
        seen.add(q)
        for rank in ranks:
            if rank == "genus":
                bins = genus_bins
            else:
                bins = {}
                depth = RANKS.index(rank) + 1
                for key, value in genus_bins.items():
                    parent = key[:depth]
                    group = bins.setdefault(parent, {"max": 0, "species": 0, "median": None})
                    group["max"] = max(group["max"], value["max"])
                    group["species"] += value["species"]
            result[rank][q] = _assign_bin(by_id[q], bins, hosts[rank], min_species, margin, method)
    return {rank: [result[rank].get(g["orf_id"]) or _assign_bin(g, {}, hosts[rank], min_species, margin, method, False)
                   for g in genes] for rank in ranks}


def attach_regions(assignments, scores, regions):
    scored = {r["orf_id"]: r for r in scores}
    memberships = defaultdict(list)
    for region in regions:
        for q in region["orf_ids"].split(";"):
            memberships[q].append(region["region_id"])
    result = []
    for r in assignments:
        s = scored[r["orf_id"]]
        row = {**r, "orf_score": s["orf_score"], "window_score": s["window_score"],
               "score_flags": s["flags"], "candidate_regions": ";".join(memberships[r["orf_id"]]),
               "in_candidate_region": bool(memberships[r["orf_id"]])}
        for key in ("raw_orf_score", "top_hit_subject", "top_hit_pident", "top_hit_qcov", "top_hit_scov",
                    "top_hit_evalue", "top_hit_quality", "donor_top_hit_quality", "min_hit_quality", "quality_status"):
            if key in s:
                row[key] = s[key]
        if "quality_status" in s:
            row.update(unfiltered_taxon=r["taxon"], unfiltered_taxon_label=r["taxon_label"],
                       unfiltered_assignment_status=r["assignment_status"],
                       unfiltered_host_taxon_match=r["host_taxon_match"],
                       unfiltered_near_tie=r["near_tie"],
                       unfiltered_ambiguity_kind=r["ambiguity_kind"])
        # A weak winning hit must not be replaced by a weaker taxonomy or host.
        # Low donor quality affects HGT screening, not a strong host assignment.
        if s.get("quality_status") == "low_top_hit_quality" or "missing_top_hit_quality" in s["flags"].split(";"):
            status = "low_hit_quality" if s.get("quality_status") == "low_top_hit_quality" else "quality_unavailable"
            row.update(unfiltered_taxon=r["taxon"], unfiltered_taxon_label=r["taxon_label"],
                       unfiltered_assignment_status=r["assignment_status"],
                       assignment_status=status, taxon="", taxon_label=STATUS_LABELS[status],
                       host_taxon_match=None, near_tie=False, ambiguity_kind="")
        result.append(row)
    return result


def contig_summary(rows):
    totals = Counter(r["contig"] for r in rows)
    groups = defaultdict(list)
    for row in rows:
        groups[(row["contig"], row["assignment_status"], row["taxon"], row["taxon_label"])].append(row)
    return [{"contig": key[0], "assignment_status": key[1], "taxon": key[2], "taxon_label": key[3],
             "orf_count": len(group), "contig_orfs": totals[key[0]],
             "fraction_of_contig_orfs": len(group) / totals[key[0]],
             "candidate_region_orfs": sum(r["in_candidate_region"] for r in group)}
            for key, group in sorted(groups.items())]


def _read_rank_cache(path):
    """Restore typed TSV values, including the assignment before quality masking."""
    import csv
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for row in rows:
        for key in ("start", "end", "strand", "protein_length", "support_species", "host_support_species"):
            if key in row:
                row[key] = int(row[key]) if row[key] not in ("", "NA") else None
        for key in ("median_bitscore", "best_bitscore", "ranking_bitscore", "relative_margin",
                    "host_median_bitscore", "host_best_bitscore", "orf_score", "window_score", "raw_orf_score",
                    "top_hit_pident", "top_hit_qcov", "top_hit_scov", "top_hit_evalue", "top_hit_quality",
                    "donor_top_hit_quality", "min_hit_quality"):
            if key in row:
                row[key] = float(row[key]) if row[key] not in ("", "NA") else None
        for key in ("host_taxon_match", "near_tie", "host_in_competition", "host_in_near_best_set",
                    "in_candidate_region", "unfiltered_host_taxon_match", "unfiltered_near_tie"):
            if key in row:
                if row[key] not in ("True", "False", "", "NA"):
                    raise ValueError(f"Invalid boolean in {path.name}: {key}")
                row[key] = {"True": True, "False": False, "": None, "NA": None}[row[key]]
        if row["assignment_status"] in ("low_hit_quality", "quality_unavailable"):
            keys = ("taxon", "taxon_label", "assignment_status", "host_taxon_match", "near_tie", "ambiguity_kind")
            if any("unfiltered_" + key not in row for key in keys):
                raise ValueError("Cached masked assignments lack original details; replot from the original run or summaries")
            row.update({key: row["unfiltered_" + key] for key in keys})
    return rows


def _load_quality_audit(path, source, scores):
    """Accept only a complete audit linked to these original rank assignments."""
    import csv
    from .io import digest
    from .quality import alignment_quality
    path = Path(path)
    metadata_path = path.with_name("summary.json")
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("complete") is not True or metadata.get("recovered_hits") != metadata.get("expected_hits"):
        raise ValueError("Quality audit is incomplete")
    hashes = metadata.get("input_hashes", {})
    if any(name not in hashes for name in ("sankey_orfs_genus.tsv", "sankey_orfs_phylum.tsv")):
        raise ValueError("Quality audit lacks source assignment hashes")
    for name, expected in hashes.items():
        if Path(name).name != name or digest(source / name) != expected:
            raise ValueError(f"Quality audit source mismatch: {name}; use its original run directory")
    with path.open() as f:
        audit_rows = list(csv.DictReader(f, delimiter="\t"))
    by_id = {row["orf_id"]: row for row in audit_rows}
    if len(by_id) != len(audit_rows) or set(by_id) != {r["orf_id"] for r in scores}:
        raise ValueError("Quality audit ORF identifiers do not match the run")
    recovered = 0
    for score in scores:
        audit = by_id[score["orf_id"]]
        if (audit["contig"], int(audit["start"]), int(audit["end"])) != (score["contig"], score["start"], score["end"]):
            raise ValueError("Quality audit ORF coordinates do not match the run")
        metrics = ("pident", "qcov", "scov", "evalue", "quality")
        present = [audit[k] not in ("", "NA") for k in metrics]
        if any(present) != all(present) or bool(audit["subject"]) != all(present):
            raise ValueError("Quality audit has incomplete hit metrics")
        score["top_hit_subject"] = audit["subject"]
        for key in metrics:
            score["top_hit_" + key] = float(audit[key]) if all(present) else None
        if all(present):
            recovered += 1
            quality = alignment_quality(score["top_hit_pident"], score["top_hit_qcov"], score["top_hit_scov"])
            if (not math.isfinite(score["top_hit_quality"]) or
                    not math.isclose(quality, score["top_hit_quality"], rel_tol=0, abs_tol=1e-12) or
                    not math.isfinite(score["top_hit_evalue"]) or score["top_hit_evalue"] < 0):
                raise ValueError("Quality audit contains invalid metrics or an inconsistent quality value")
        elif audit["genus_status"] != "no_usable_hit" or audit["phylum_status"] != "no_usable_hit":
            raise ValueError("Quality audit is missing a usable top hit")
    if recovered != metadata["expected_hits"]:
        raise ValueError("Quality audit hit count does not match its metadata")
    return {"audit_tsv": str(path.resolve()), "audit_tsv_sha256": digest(path),
            "audit_metadata_sha256": digest(metadata_path), "validated_source_hashes": hashes}


def replot_from_run(source, folder, method="best-hit", min_species=2, margin=.01,
                    max_contigs=20, max_taxa=12, *, min_hit_quality=None, quality_audit=None):
    """Use cached genus max/median values; retain original score files byte for byte."""
    import csv
    import shutil
    import time
    from datetime import datetime, timezone
    from .io import digest, dump_json
    started = time.perf_counter()
    source, folder = Path(source), Path(folder)
    if folder.exists() and any(folder.iterdir()):
        raise ValueError("Replot output must be empty")
    manifest = json.loads((source / "run.json").read_text())
    if quality_audit and min_hit_quality is None:
        min_hit_quality = .5
    if min_hit_quality is not None and (not math.isfinite(min_hit_quality) or not 0 <= min_hit_quality <= 1):
        raise ValueError("min-hit-quality must be finite and between 0 and 1")
    if any(manifest["parameters"][r] != "genus" for r in ("close_rank", "donor_rank")):
        raise ValueError("Cached replot requires close-rank and donor-rank genus; other ranks need species hits")
    with (source / "orf_scores.tsv").open() as f:
        scores = list(csv.DictReader(f, delimiter="\t"))
    for row in scores:
        for key in ("start", "end"):
            row[key] = int(row[key])
        for key in ("orf_score", "window_score", "raw_orf_score", "top_hit_pident", "top_hit_qcov", "top_hit_scov",
                    "top_hit_evalue", "top_hit_quality", "donor_top_hit_quality", "min_hit_quality"):
            if key in row:
                row[key] = float(row[key]) if row[key] not in ("", "NA") else None
    genes = [{k: row[k] for k in ("orf_id", "contig", "start", "end", "strand", "gene", "product", "protein_length")
              if k in row} for row in scores]
    with (source / "regions.tsv").open() as f:
        regions = list(csv.DictReader(f, delimiter="\t"))
    ranks = DISPLAY_RANKS if method == "best-hit" else ("genus",)
    cached_config = manifest.get("sankey", {})
    use_rank_cache = (min_hit_quality is not None and
                      cached_config.get("method") == method and cached_config.get("min_species") == min_species and
                      cached_config.get("margin") == margin and
                      all((source / f"sankey_orfs_{rank}.tsv").exists() for rank in ranks))
    assignment_hashes = {}
    if use_rank_cache:
        assignments = {}
        for rank in ranks:
            path = source / f"sankey_orfs_{rank}.tsv"
            values = _read_rank_cache(path)
            if [r["orf_id"] for r in values] != [r["orf_id"] for r in scores]:
                raise ValueError("Cached rank assignments do not match the ORF score order")
            if any((r["contig"], r["start"], r["end"]) != (s["contig"], s["start"], s["end"])
                   for r, s in zip(values, scores)):
                raise ValueError("Cached rank assignment coordinates do not match the scores")
            assignments[rank] = values
            assignment_hashes[path.name] = digest(path)
        # Display-only overlays keep the original score TSV byte for byte, so
        # their recovered alignment metrics live in the saved rank TSV instead.
        for s, r in zip(scores, assignments["genus"]):
            for key in ("top_hit_subject", "top_hit_pident", "top_hit_qcov", "top_hit_scov",
                        "top_hit_evalue", "top_hit_quality"):
                if key not in s and key in r:
                    s[key] = r[key]
    else:
        assignments = assign_cached_by_rank(genes, source / "taxon_summaries.tsv", manifest["query_taxonomy"],
                                           ranks, min_species, margin, method)
    quality_provenance = _load_quality_audit(quality_audit, source, scores) if quality_audit else {}
    if min_hit_quality is not None:
        if min_hit_quality > 0 and not all("top_hit_quality" in s for s in scores):
            raise ValueError("Saved run lacks top-hit quality. Supply --quality-audit top_hit_quality.tsv, or analyze saved hits first")
        genus_by_id = {r["orf_id"]: r for r in assignments["genus"]}
        for s in scores:
            quality = s.get("top_hit_quality")
            if quality is not None and (not math.isfinite(quality) or not 0 <= quality <= 1):
                raise ValueError("Invalid saved top-hit quality")
            s["min_hit_quality"] = min_hit_quality
            s["quality_status"] = ("not_applied" if min_hit_quality == 0 else
                                   "no_usable_hit" if genus_by_id[s["orf_id"]]["assignment_status"] == "no_usable_hit" else
                                   "unavailable" if quality is None else
                                   "low_top_hit_quality" if quality < min_hit_quality else "pass")
            # These copies are for display; the original score TSV stays unchanged.
            flags = [f for f in s["flags"].split(";") if f != "missing_top_hit_quality"]
            if s["quality_status"] == "unavailable":
                flags.append("missing_top_hit_quality")
            s["flags"] = ";".join(flags)
    rows_by_rank = {rank: attach_regions(values, scores, regions) for rank, values in assignments.items()}
    rows = rows_by_rank["genus"]
    manifest["sankey"] = dict(rank="genus", method=method, min_species=min_species,
                              margin=margin, max_contigs=max_contigs, max_taxa=max_taxa,
                              min_hit_quality=min_hit_quality if min_hit_quality is not None else manifest["parameters"].get("min_hit_quality", 0),
                              quality_display_only=min_hit_quality is not None,
                              interpretation="hit taxonomy view; not HGT probability or contamination fraction")
    manifest["sankey_replot"] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_run_directory": str(source.resolve()),
        "source_run_sha256": digest(source / "run.json"),
        "source_orf_scores_sha256": digest(source / "orf_scores.tsv"),
        "source_regions_sha256": digest(source / "regions.tsv"),
        "hgt_scores_recomputed": False,
        "summary_source": str(source.resolve()) if use_rank_cache else str((source / "taxon_summaries.tsv").resolve()),
        "assignment_hashes": assignment_hashes,
        "search_rerun": False}
    if min_hit_quality is not None:
        manifest["sankey_replot"]["quality_overlay"] = {**quality_provenance, "min_hit_quality": min_hit_quality,
            "formula": "pident/100 * qcov * scov", "scope": "taxonomy display only; original ORF/window scores and candidate regions retained"}
    folder.mkdir(parents=True, exist_ok=True)
    for name in ("orf_scores.tsv", "regions.tsv", "report.html"):
        if (source / name).exists():
            shutil.copy2(source / name, folder / name)
    write_sankey(folder, rows, manifest, regions, rows_by_rank)
    dump_json(folder / "run.json", manifest)
    summary = dict(rank_counts={rank: {"status_counts": dict(Counter(r["assignment_status"] for r in values)),
                                      "host_assigned": sum(r["host_taxon_match"] is True for r in values),
                                      "other_assigned": sum(r["assignment_status"] == "assigned" and r["host_taxon_match"] is False for r in values),
                                      "near_tie_assigned": sum(r["near_tie"] for r in values)}
                               for rank, values in rows_by_rank.items()},
                   status_counts=dict(Counter(r["assignment_status"] for r in rows)),
                   host_assigned=sum(r["host_taxon_match"] is True for r in rows),
                   near_tie_assigned=sum(r["near_tie"] for r in rows),
                   ambiguity_counts=dict(Counter(r["ambiguity_kind"] for r in rows if r["ambiguity_kind"])),
                   total_orfs=len(rows), candidate_orfs=sum(r["in_candidate_region"] for r in rows),
                   elapsed_seconds=time.perf_counter() - started,
                   hgt_scores_recomputed=False)
    summary["min_hit_quality"] = manifest["sankey"]["min_hit_quality"]
    dump_json(folder / "sankey_replot_summary.json", summary)
    return summary


def write_sankey(folder: Path, rows, manifest, regions=None, rows_by_rank=None):
    tsv_write(folder / "sankey_orfs.tsv", rows)
    tsv_write(folder / "contig_taxonomy.tsv", contig_summary(rows))
    cfg = manifest["sankey"]
    rows_by_rank = rows_by_rank or {cfg["rank"]: rows}
    cfg["available_ranks"] = list(rows_by_rank)
    common_fields = {"orf_id", "contig", "start", "end", "strand", "gene", "product", "protein_length",
                     "orf_score", "window_score", "score_flags", "candidate_regions", "in_candidate_region",
                     "raw_orf_score", "top_hit_subject", "top_hit_pident", "top_hit_qcov", "top_hit_scov",
                     "top_hit_evalue", "top_hit_quality", "donor_top_hit_quality", "min_hit_quality", "quality_status"}
    rank_assignments = {}
    for rank, values in rows_by_rank.items():
        if [r["orf_id"] for r in values] != [r["orf_id"] for r in rows]:
            raise ValueError("Rank assignment rows must have identical ORF order")
        tsv_write(folder / f"sankey_orfs_{rank}.tsv", values)
        tsv_write(folder / f"contig_taxonomy_{rank}.tsv", contig_summary(values))
        if rank != cfg["rank"]:
            rank_assignments[rank] = [{k: v for k, v in r.items() if k not in common_fields} for r in values]
    payload = {"rows": rows, "config": cfg, "rank_assignments": rank_assignments}
    rank_buttons = ''.join('<button type="button" data-rank="' + rank + '" aria-pressed="false">' +
                           RANK_LABELS[rank] + ' (' + rank + ')</button>' for rank in cfg["available_ranks"])
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")
    logic = Path(__file__).with_name("sankey.js").read_text()
    options = ''.join('<option value="' + html.escape(c, quote=True) + '">' + html.escape(c) + '</option>'
                      for c in sorted({r["contig"] for r in rows}))
    region_options = ''.join('<option value="region:' + html.escape(r["region_id"], quote=True) + '">' +
                             html.escape(f'{r["region_id"]} · {r["contig"]}:{r["start"]}–{r["end"]} ({r["n_orfs"]} ORF)') +
                             '</option>' for r in (regions or []))
    label = html.escape(manifest["database_label"])
    page = """<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HGT_Finder · Sankey</title><style>
body{font:15px/1.6 system-ui;background:#f5f8fb;color:#173140;margin:0}main{max-width:1450px;margin:auto;padding:24px}h1{margin-bottom:4px}.notice{background:#fff3d7;padding:12px 16px;border-left:4px solid #d29b26}.rank-controls{display:flex;gap:8px;flex-wrap:wrap;margin:16px 0}.rank-controls button[aria-pressed=true]{font-weight:bold;border:2px solid #173140}.controls{display:flex;gap:14px;flex-wrap:wrap;margin:20px 0}input,select,button{font:inherit;padding:5px}input[type=number]{width:120px}#diagram{background:white;border:1px solid #d4e1e8;overflow:auto}#diagram svg{display:block;width:100%;min-width:950px}#detail{min-height:28px;padding:10px;background:#eaf2f7}.scroll{max-height:440px;overflow:auto}table{color:black;background:white;border-collapse:collapse;width:100%;font-size:12px}td,th{border:1px solid black;padding:6px;text-align:left}th{position:sticky;top:0;background:white}a{color:#146b98}.muted{color:#536c7b}
</style></head><body><main><h1>Contig → ヒットの所属分類群</h1><p>HGT_Finder · サンキープロット</p>
<p><strong>参照DB：</strong>""" + label + """</p>
<div class="notice">分類の混在を探索する表示です。汚染率やGUNCの判定を計算したものではなく、HGTと汚染をこの図だけで区別できません。</div>
<div id="rank-controls" class="rank-controls" aria-label="表示する分類階層">""" + rank_buttons + """</div>
<p id="method"></p>
<div class="controls"><label>contig <select id="contig"><option value="">すべて</option>""" + options + """</select></label>
<label>対象 <select id="scope"><option value="all">全ORF</option><option value="candidate">全候補領域内のORF</option><option value="near">最高ヒット・別分類群と僅差</option>""" + region_options + """</select></label>
<label>開始bp <input id="from" type="number" min="1" placeholder="制限なし"></label><label>終了bp <input id="to" type="number" min="1" placeholder="制限なし"></label>
<button id="apply">表示を更新</button><button id="reset">全体に戻す</button><button id="export">SVG保存</button></div>
<p class="muted">座標範囲を指定すると、その範囲に全体が収まるORFを表示します。</p>
<p id="stats"></p><div id="diagram"></div><div id="detail">線にカーソルを置くと内訳を表示します。線をクリックすると下のORF一覧を絞れます。</div>
<p class="muted">線の太さ＝ORF数。各ORFは1票です。青緑は宿主と同じ表示分類群、橙系はその他の分類群、灰色は曖昧・評価不能です。小さな分類群・contigは「その他」にまとめますが、ORFは除外しません。</p>
<h2>ORF一覧</h2><p id="table-count"></p><button id="clear-flow">線による絞り込みを解除</button><div class="scroll"><table><thead><tr><th>ORF</th><th>contig</th><th>開始</th><th>終了</th><th>遺伝子</th><th>注釈</th><th>類似先</th><th>状態</th><th>最高ヒットの一致率 (%)</th><th>query coverage</th><th>subject coverage</th><th>ヒット品質</th><th>宿主分類群のヒット種数</th><th>宿主分類群の最高bit score</th><th>競合する分類群</th><th>ORFスコア</th><th>候補領域</th></tr></thead><tbody id="orf-body"></tbody></table></div>
<p><a href="report.html">スコアのグラフ</a> · <a id="orf-tsv" href="sankey_orfs.tsv">全ORFの分類・支持情報TSV</a> · <a id="contig-tsv" href="contig_taxonomy.tsv">contig別の集計TSV</a> · <a href="run.json">実行設定</a></p>
</main><script>const payload=""" + encoded + ";\n" + logic + "\n" + """
const el=id=>document.getElementById(id);let activeRank=payload.config.rank,rankRows=rowsForRank(payload,activeRank),colorMap=makeTaxonColors(rankRows),current=[],graph=null,flow=null;
const rankCache=new Map([[activeRank,rankRows]]);const rankNames={phylum:'門',class:'綱',order:'目',family:'科',genus:'属'};
const bestHit=payload.config.method==='best-hit';
function updateMethod(){
el('method').textContent='表示階層：'+rankNames[activeRank]+' ('+activeRank+')'+' ／ 最低支持：'+payload.config.min_species+'種 ／ 僅差の基準：'+(100*payload.config.margin)+'%。'+(bestHit?'有効ヒット中の最高bit scoreの所属分類群を表示します。別分類群との同点は灰色、僅差は所属先を表示し一覧に記録します。最低支持はその群の全有効ヒット種数です。分類の確実性やHGT確率を示すものではありません。':'従来方式：分類群ごとの種別最良bit scoreの中央値を比較し、僅差も曖昧にします。')+'自己除外・有意性・coverage条件はHGT解析と共通です。';
el('method').textContent+=(payload.config.min_hit_quality>0?' 品質の下限：'+payload.config.min_hit_quality+'（一致率/100 × query coverage × subject coverage）。低品質の最高ヒットは由来判定保留です。':' この保存結果には品質による除外を適用していません。');
if(payload.config.quality_display_only)el('method').textContent+=' 品質条件はこの分類表示に適用しています。候補領域とORF/windowスコアは元のHGT解析結果を保持しています。';
}
function setRank(rank){
 if(!payload.config.available_ranks.includes(rank))return;
 activeRank=rank;if(!rankCache.has(rank))rankCache.set(rank,rowsForRank(payload,rank));
 rankRows=rankCache.get(rank);colorMap=makeTaxonColors(rankRows);
 for(const button of el('rank-controls').querySelectorAll('[data-rank]'))button.setAttribute('aria-pressed',String(button.dataset.rank===rank));
 el('orf-tsv').href='sankey_orfs_'+rank+'.tsv';el('contig-tsv').href='contig_taxonomy_'+rank+'.tsv';
 updateMethod();update();
}
function table(){
 const rows=flow?current.filter(r=>graph.contigMap[r.contig]===flow.source&&graph.taxonMap[taxonId(r)]===flow.target):current;
 el('table-count').textContent=rows.length.toLocaleString()+' ORF（一覧は先頭200件。全件はTSVに保存）';
 el('orf-body').innerHTML=rows.slice(0,200).map(r=>'<tr>'+[r.orf_id,r.contig,r.start,r.end,r.gene,r.product,r.taxon_label,r.assignment_status==='assigned'&&r.near_tie?'最高ヒット・別分類群と僅差':r.assignment_status,r.top_hit_pident??'NA',r.top_hit_qcov==null?'NA':r.top_hit_qcov.toFixed(3),r.top_hit_scov==null?'NA':r.top_hit_scov.toFixed(3),r.top_hit_quality==null?'NA':r.top_hit_quality.toFixed(3),r.host_support_species,r.host_best_bitscore??'NA',r.nearby_taxa||'—',r.orf_score===null?'NA':r.orf_score.toFixed(3),r.candidate_regions||'—'].map(v=>'<td>'+esc(v)+'</td>').join('')+'</tr>').join('');
}
function update(){
 const lo=el('from').value===''?1:Number(el('from').value),hi=el('to').value===''?Infinity:Number(el('to').value);
 if(!Number.isFinite(lo)||Number.isNaN(hi)||lo<1||hi<lo){el('detail').textContent='開始・終了座標を確認してください。';return;}
 const scope=el('scope').value;
 current=rankRows.filter(r=>(!el('contig').value||r.contig===el('contig').value)&&(scope==='all'||(scope==='near'&&r.near_tie)||(scope==='candidate'&&r.in_candidate_region)||(scope.startsWith('region:')&&r.candidate_regions.split(';').includes(scope.slice(7))))&&r.start>=lo&&r.end<=hi);
 flow=null;graph=aggregateRows(current,payload.config.max_contigs,payload.config.max_taxa,colorMap);el('diagram').innerHTML=svgMarkup(graph);
 const assigned=current.filter(r=>r.assignment_status==='assigned').length;
 el('stats').textContent=current.length.toLocaleString()+' ORF ／ 分類群を選べたORF '+assigned.toLocaleString()+' ／ 曖昧・評価不能 '+(current.length-assigned).toLocaleString()+' ／ 選べたORFのうち別分類群と僅差 '+current.filter(r=>r.near_tie).length.toLocaleString();
 el('detail').textContent=current.length?'線にカーソルを置くと内訳を表示します。':'条件に合うORFがありません。';
 for(const path of el('diagram').querySelectorAll('[data-link]')){
  const link=graph.links[Number(path.dataset.link)];
  path.onmouseenter=()=>{el('detail').textContent=graph.nodes[link.source].label+' → '+graph.nodes[link.target].label+'：'+link.value+' ORF（表示対象の'+(100*link.value/graph.total).toFixed(1)+'%）';};
  path.onclick=()=>{flow=link;table();};
 }
 table();
}
el('apply').onclick=update;el('contig').onchange=update;el('scope').onchange=update;
el('reset').onclick=()=>{el('contig').value='';el('scope').value='all';el('from').value='';el('to').value='';update();};
el('clear-flow').onclick=()=>{flow=null;table();};
el('export').onclick=()=>{const url=URL.createObjectURL(new Blob([svgMarkup(graph)],{type:'image/svg+xml'}));const a=document.createElement('a');a.href=url;a.download='contig_taxonomy_'+activeRank+'_sankey.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
for(const button of el('rank-controls').querySelectorAll('[data-rank]'))button.onclick=()=>setRank(button.dataset.rank);
setRank(activeRank);
</script></body></html>"""
    (folder / "sankey.html").write_text(page)
