from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .core import collect_hits, score_genes, self_scores, smooth_and_regions
from .io import (RANKS, canonical_accession, digest, dump_json, lineage, prepare_genbank,
                 read_genes, read_subject_map, read_taxonomy, tsv_write)
from .report import write_report
from .search import SENSITIVITY_MODES, build_db, search, search_self
from .storage import SpeciesHitStore
from .sankey import assign_orfs_by_rank, DISPLAY_RANKS, attach_regions, write_sankey, replot_from_run
from .distribution import plot_from_run as plot_distributions


def add_sankey_options(p):
    p.add_argument("--sankey-method", choices=("best-hit", "median"), default="best-hit",
                   help="Best-hit taxonomy view (default), or legacy genus-median comparison")
    p.add_argument("--sankey-min-species", type=int, default=2)
    p.add_argument("--sankey-margin", type=float, default=0.01,
                   help="Relative score gap (default 0.01 = 1%%): flag close alternatives in best-hit mode; ambiguous in median mode")
    p.add_argument("--sankey-max-contigs", type=int, default=20)
    p.add_argument("--sankey-max-taxa", type=int, default=12)


def add_analysis(p):
    p.add_argument("--taxonomy", type=Path, required=True, help="GTDB taxonomy TSV/gz; actual genomes in the search DB")
    q = p.add_mutually_exclusive_group(required=True)
    q.add_argument("--query-taxonomy", help="Seven-rank GTDB lineage for the host genome")
    q.add_argument("--query-genome", help="Host genome accession present in taxonomy table")
    p.add_argument("--subject-map", type=Path, help="For existing DB IDs: TSV with sseqid, genome_id")
    p.add_argument("--reference-genomes", type=Path, help="One genome ID per line; restrict taxonomy to actual DB genomes")
    p.add_argument("--close-rank", choices=RANKS[1:-1], default="genus", help="Close-group boundary (default genus)")
    p.add_argument("--donor-rank", choices=RANKS[1:-1], default="genus", help="Distal hit aggregation rank (default genus)")
    p.add_argument("--self-exclude-rank", choices=("none", "species", "genus"), default="species")
    p.add_argument("--min-close-species", type=int, default=3, help="Minimum distinct close reference species for window eligibility")
    p.add_argument("--min-donor-species", type=int, default=2)
    p.add_argument("--missing-close", choices=("zero", "na"), default="zero",
                   help="Explicit close-hit censoring policy: zero (flagged) or NA")
    p.add_argument("--min-query-cover", type=float, default=0.5)
    p.add_argument("--min-subject-cover", type=float, default=0.5)
    p.add_argument("--min-hit-quality", type=float, default=0.5,
                   help="WAAFLE-inspired protein quality floor: identity/100 * query coverage * subject coverage (default 0.5, provisional; 0 disables). Keeps comparison hits; weak evidence becomes NA.")
    p.add_argument("--evalue", type=float, default=1e-10)
    p.add_argument("--window", type=int, default=5, help="Odd count of adjacent ORFs (default 5)")
    p.add_argument("--min-valid", type=int, default=3)
    p.add_argument("--max-gap-bp", type=int, default=10000)
    p.add_argument("--threshold", type=float, default=0.15, help="Provisional score threshold; not validated probability")
    p.add_argument("--min-region-orfs", type=int, default=3)
    p.add_argument("--database-label", default="GTDB reference database (release not verified)")
    p.add_argument("--sankey", action="store_true", help="Write offline interactive contig taxonomy Sankey plot")
    p.add_argument("--sankey-rank", choices=RANKS[1:-1], default="genus")
    add_sankey_options(p)
    p.add_argument("--hit-distributions", action="store_true", help="Write species-best distribution pages for all ORFs (may create large outputs)")
    p.add_argument("--distribution-orf", action="append", help="Generate a distribution only for this ORF; repeat for multiple ORFs")
    p.add_argument("--distribution-rank", choices=DISPLAY_RANKS, default="genus")
    p.add_argument("--distribution-max-taxa", type=int, default=12)


def add_search(p):
    p.add_argument("--diamond", default="diamond")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--sensitivity", choices=SENSITIVITY_MODES, default="very-sensitive",
                   help="DIAMOND search mode; very-sensitive is the sensitivity-preserving default; faster/fast favor close matches; default omits mode flags")
    p.add_argument("--max-target-seqs", type=int, default=0, help="0 = no target count cap (default); significance filters still apply")


def parser():
    p = argparse.ArgumentParser(description="HGT_Finder: ORF similarity distribution + rolling region screening")
    p.add_argument("--version", action="version", version=__version__)
    commands = p.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare", help="Extract translated ORFs and coordinates from GenBank")
    prep.add_argument("--genbank", required=True, type=Path)
    prep.add_argument("--out", required=True, type=Path)
    build = commands.add_parser("build-db", help="Build encoded DIAMOND DB from per-genome protein FAA files")
    build.add_argument("--proteins", required=True, type=Path, help="Directory or tar archive of per-genome FAA files")
    build.add_argument("--taxonomy", required=True, type=Path)
    build.add_argument("--out", required=True, type=Path)
    build.add_argument("--database-label", default="GTDB R232 (user-supplied source)")
    build.add_argument("--diamond", default="diamond")
    build.add_argument("--threads", type=int, default=4)
    scan = commands.add_parser("search", help="Run DIAMOND with unlimited target count by default")
    scan.add_argument("--query-faa", required=True, type=Path)
    scan.add_argument("--db", required=True, type=Path)
    scan.add_argument("--out", required=True, type=Path)
    scan.add_argument("--evalue", type=float, default=1e-10)
    add_search(scan)
    ana = commands.add_parser("analyze", help="Analyze saved 14-column DIAMOND hits with query coordinates")
    ana.add_argument("--hits", required=True, type=Path)
    ana.add_argument("--self-hits", required=True, type=Path)
    ana.add_argument("--genes", required=True, type=Path)
    ana.add_argument("--out", required=True, type=Path)
    ana.add_argument("--assume-unlimited-hits", action="store_true", help="Explicit attestation for external hit files without search manifest")
    add_analysis(ana)
    run = commands.add_parser("run", help="Complete local workflow from GenBank or FAA + genes TSV")
    q = run.add_mutually_exclusive_group(required=True)
    q.add_argument("--genbank", type=Path)
    q.add_argument("--query-faa", type=Path)
    run.add_argument("--genes", type=Path, help="Required with --query-faa")
    run.add_argument("--db", type=Path, required=True)
    run.add_argument("--out", type=Path, required=True)
    add_search(run)
    add_analysis(run)
    plot = commands.add_parser("sankey", help="Replot saved genus summaries without rerunning searches or HGT scores")
    plot.add_argument("--run-dir", required=True, type=Path)
    plot.add_argument("--out", required=True, type=Path)
    plot.add_argument("--min-hit-quality", type=float, default=0.5,
                      help="Plot-only top-hit quality floor (default 0.5; 0 disables); HGT scores and regions are retained")
    plot.add_argument("--quality-audit", type=Path,
                      help="Audited top_hit_quality.tsv for an older run without saved quality metrics")
    add_sankey_options(plot)
    dist = commands.add_parser("distribution", help="Plot individual species-best hits from a saved SQLite store; no search or rescoring")
    dist.add_argument("--run-dir", required=True, type=Path)
    dist.add_argument("--out", required=True, type=Path)
    selection = dist.add_mutually_exclusive_group(required=True)
    selection.add_argument("--orf-id", action="append", help="ORF ID; repeat to select several")
    selection.add_argument("--all-orfs", action="store_true", help="Generate separate pages for every ORF")
    dist.add_argument("--taxonomy", type=Path, help="Original taxonomy TSV for a legacy SQLite store")
    dist.add_argument("--subject-map", type=Path, help="Subject ID mapping for legacy non-encoded protein IDs")
    dist.add_argument("--rank", dest="distribution_rank", choices=DISPLAY_RANKS, default="genus")
    dist.add_argument("--max-taxa", dest="distribution_max_taxa", type=int, default=12)
    dist.add_argument("--min-hit-quality", type=float, default=.5, help="Quality marking/filter threshold (default .5); all comparison points shown initially")
    return p


def validate(a):
    if hasattr(a, "distribution_max_taxa") and a.distribution_max_taxa < 1:
        raise ValueError("Distribution display limit must be positive")
    if hasattr(a, "min_hit_quality") and (not math.isfinite(a.min_hit_quality) or not 0 <= a.min_hit_quality <= 1):
        raise ValueError("min-hit-quality must be finite and between 0 and 1")
    if hasattr(a, "threads") and a.threads < 1:
        raise ValueError("threads must be positive")
    if hasattr(a, "max_target_seqs") and a.max_target_seqs < 0:
        raise ValueError("max-target-seqs cannot be negative")
    if hasattr(a, "evalue") and (not math.isfinite(a.evalue) or a.evalue <= 0):
        raise ValueError("evalue must be finite and positive")
    if hasattr(a, "close_rank"):
        if a.self_exclude_rank != "none" and RANKS.index(a.self_exclude_rank) <= RANKS.index(a.close_rank):
            raise ValueError("Self-exclude rank must be narrower than close rank")
        if a.min_close_species < 1 or a.min_donor_species < 1:
            raise ValueError("Species count requirements must be positive")
        if not 0 <= a.min_query_cover <= 1 or not 0 <= a.min_subject_cover <= 1:
            raise ValueError("Coverage fractions must be between 0 and 1")
        if a.max_gap_bp < 0 or a.min_region_orfs < 1:
            raise ValueError("max-gap-bp must be nonnegative and min-region-orfs positive")
        if not math.isfinite(a.threshold):
            raise ValueError("threshold must be finite")
        if a.window < 1 or a.window % 2 == 0 or not 1 <= a.min_valid <= a.window:
            raise ValueError("Window must be positive odd; min-valid must be between 1 and window")
    if hasattr(a, "sankey_min_species"):
        if a.sankey_min_species < 1 or a.sankey_max_contigs < 1 or a.sankey_max_taxa < 1:
            raise ValueError("Sankey species requirement and display limits must be positive")
        if not math.isfinite(a.sankey_margin) or not 0 <= a.sankey_margin < 1:
            raise ValueError("Sankey margin must be finite and between 0 (inclusive) and 1 (exclusive)")


def analyze(a, folder: Path):
    taxa = read_taxonomy(a.taxonomy)
    query_tax = lineage(a.query_taxonomy) if a.query_taxonomy else taxa.get(canonical_accession(a.query_genome))
    if query_tax is None:
        raise ValueError("Host genome not found in taxonomy table")
    if a.reference_genomes:
        available = {canonical_accession(s) for s in a.reference_genomes.read_text().splitlines() if s.strip()}
        unknown = available - set(taxa)
        if unknown:
            raise ValueError(f"Reference genome list has missing taxonomy: {sorted(unknown)[:5]}")
        taxa = {g: t for g, t in taxa.items() if g in available}
    if not taxa:
        raise ValueError("No reference genomes after filtering")
    genes = read_genes(a.genes)
    if a.distribution_orf:
        unknown_orfs = set(a.distribution_orf) - {g["orf_id"] for g in genes}
        if unknown_orfs:
            raise ValueError(f"Unknown distribution ORF IDs: {sorted(unknown_orfs)[:5]}")
    selfs = self_scores(a.self_hits)
    subject_map = read_subject_map(a.subject_map)
    search_manifest = a.hits.with_suffix(a.hits.suffix + ".json")
    provenance = json.loads(search_manifest.read_text()) if search_manifest.exists() else {}
    unlimited = provenance.get("max_target_seqs") == 0 or getattr(a, "assume_unlimited_hits", False)
    missing_policy = a.missing_close
    warnings = []
    if not unlimited and missing_policy == "zero":
        missing_policy = "na"
        warnings.append("Missing-close policy changed to NA: unlimited hit output was not verified.")
    disk = SpeciesHitStore(folder / "species_hits.sqlite", taxa, a.close_rank, a.donor_rank, query_tax)
    try:
        best, counts, audit = collect_hits(a.hits, taxa, query_tax, a.close_rank, a.donor_rank,
                                         a.self_exclude_rank, a.query_genome, subject_map, a.evalue,
                                         a.min_query_cover, a.min_subject_cover, disk)
        scored, summaries = score_genes(genes, best, counts, selfs, taxa, query_tax,
                                       a.close_rank, a.donor_rank, a.self_exclude_rank, a.query_genome,
                                       a.min_close_species, a.min_donor_species, missing_policy, a.min_hit_quality)
        assignments = assign_orfs_by_rank(genes, best, query_tax, DISPLAY_RANKS,
                                         a.sankey_min_species, a.sankey_margin, a.sankey_method) if a.sankey else None
    finally:
        disk.close()
    usable = sum(r["eligible_protein_hits"] for r in scored)
    if usable == 0:
        raise ValueError("No taxonomy-mapped eligible hits for any ORF. Check ID mapping, taxonomy, self exclusions and coverage filters.")
    scored, regions = smooth_and_regions(scored, a.window, a.min_valid, a.threshold,
                                         a.max_gap_bp, a.min_region_orfs)
    tsv_write(folder / "orf_scores.tsv", scored)
    region_fields = ["region_id", "contig", "start", "end", "n_orfs", "n_scored_orfs", "mean_window_score", "mean_raw_score", "donor_taxon", "donor_agreement", "orf_ids", "boundary_note"]
    tsv_write(folder / "regions.tsv", regions, region_fields)
    tsv_write(folder / "taxon_summaries.tsv", summaries,
              ["orf_id", "group", "taxon", "species_count", "median_bitscore", "min_bitscore", "max_bitscore"])
    parameters = {k: getattr(a, k) for k in ("close_rank", "donor_rank", "self_exclude_rank", "min_close_species", "min_donor_species", "min_query_cover", "min_subject_cover", "min_hit_quality", "evalue", "window", "min_valid", "max_gap_bp", "threshold", "min_region_orfs")}
    parameters["missing_close"] = missing_policy
    manifest = {"version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
                "database_label": a.database_label, "query_taxonomy": query_tax,
                "parameters": parameters, "warnings": warnings, "search": provenance,
                "score_definition": "(best-supported distal taxon median bit score - close species median bit score) / query self bit score",
                "quality_definition": "protein adaptation of WAAFLE: (pident / 100) * qcov * scov; no contig-end correction or multi-hit per-site merging; provisional floor on top bit-score hit, also donor-bin top hit for positive ORF scores; comparison hits retained",
                "candidate_status": "provisional similarity screen; not HGT probability or validated event direction",
                "counts": {"orfs": len(scored), "window_eligible_orfs": sum(r["eligible_for_window"] for r in scored),
                           "quality_excluded_orfs": sum(r["quality_status"] in ("low_top_hit_quality", "low_donor_hit_quality", "unavailable") for r in scored),
                           "raw_hits": audit.get("raw_hits", 0), "unknown_taxonomy_hits": audit.get("unknown_taxonomy", 0), "regions": len(regions)},
                "inputs": {"hits": str(a.hits.resolve()), "self_hits": str(a.self_hits.resolve()),
                           "genes": str(a.genes.resolve()), "taxonomy": str(a.taxonomy.resolve()),
                           "taxonomy_sha256": digest(a.taxonomy), "genes_sha256": digest(a.genes)}}
    if a.subject_map:
        manifest["inputs"]["subject_map"] = str(a.subject_map.resolve())
    if a.hit_distributions or a.distribution_orf:
        dump_json(folder / "run.json", manifest)
        distribution_result = plot_distributions(folder, folder / "hit_distributions", a.distribution_orf,
                                                a.distribution_rank, a.distribution_max_taxa, a.min_hit_quality)
        manifest["hit_distributions"] = {"index": "hit_distributions/index.html",
                                         "orf_ids": [r["orf_id"] for r in distribution_result["entries"]]}
    if a.sankey:
        manifest["sankey"] = {"rank": a.sankey_rank, "min_species": a.sankey_min_species,
                              "margin": a.sankey_margin, "max_contigs": a.sankey_max_contigs,
                              "max_taxa": a.sankey_max_taxa,
                              "method": a.sankey_method,
                              "min_hit_quality": a.min_hit_quality,
                              "interpretation": "taxonomy of strongest filtered hit; close alternatives retained; not HGT probability, GUNC score or contamination fraction"}
        if "hit_distributions" in manifest:
            manifest["sankey"].update(distribution_index=manifest["hit_distributions"]["index"],
                                      distribution_orfs=manifest["hit_distributions"]["orf_ids"])
        rows_by_rank = {rank: attach_regions(values, scored, regions) for rank, values in assignments.items()}
        write_sankey(folder, rows_by_rank[a.sankey_rank], manifest, regions, rows_by_rank)
    dump_json(folder / "run.json", manifest)
    write_report(folder / "report.html", scored, regions, manifest)
    for warning in warnings:
        print("WARNING: " + warning, file=sys.stderr)
    print(json.dumps(manifest["counts"], ensure_ascii=False))
    print("Report: " + str((folder / "report.html").resolve()))


def fresh_output(path: Path):
    if path.exists() and any(path.iterdir()):
        raise ValueError("Output directory must be empty; use a new path for each parameter comparison")
    path.mkdir(parents=True, exist_ok=True)


def main(argv=None):
    a = parser().parse_args(argv)
    try:
        validate(a)
        if a.command == "prepare":
            fresh_output(a.out)
            print(json.dumps(prepare_genbank(a.genbank, a.out)))
        elif a.command == "build-db":
            print(json.dumps(build_db(a.proteins, a.taxonomy, a.out, a.diamond, a.threads, a.database_label)))
        elif a.command == "search":
            print(json.dumps(search(a.query_faa, a.db, a.out, a.diamond, a.threads, a.evalue, a.max_target_seqs, a.sensitivity)))
        elif a.command == "analyze":
            fresh_output(a.out)
            analyze(a, a.out)
        elif a.command == "sankey":
            fresh_output(a.out)
            print(json.dumps(replot_from_run(a.run_dir, a.out, a.sankey_method,
                             a.sankey_min_species, a.sankey_margin,
                             a.sankey_max_contigs, a.sankey_max_taxa,
                             min_hit_quality=a.min_hit_quality,
                             quality_audit=a.quality_audit), ensure_ascii=False))
            print("Sankey: " + str((a.out / "sankey.html").resolve()))
        elif a.command == "distribution":
            fresh_output(a.out)
            result = plot_distributions(a.run_dir, a.out, a.orf_id, a.distribution_rank,
                                        a.distribution_max_taxa, a.min_hit_quality, a.taxonomy, a.subject_map)
            print(json.dumps({k: result[k] for k in ("orfs", "species_hits", "search_rerun", "hgt_scores_recomputed")}, ensure_ascii=False))
            print("Distributions: " + str((a.out / "index.html").resolve()))
        elif a.command == "run":
            run_started = time.perf_counter()
            stages = {}
            def timed(name, function, *args):
                print(f"Starting stage: {name}", file=sys.stderr, flush=True)
                started = time.perf_counter()
                try:
                    return function(*args)
                finally:
                    stages[name] = time.perf_counter() - started
                    print(f"Elapsed stage: {name} = {stages[name]:.3f} seconds",
                          file=sys.stderr, flush=True)
            if a.query_faa and not a.genes:
                raise ValueError("--genes is required with --query-faa")
            fresh_output(a.out)
            if a.genbank:
                timed("orf_preparation", prepare_genbank, a.genbank, a.out / "input")
                a.query_faa, a.genes = a.out / "input/query.faa", a.out / "input/genes.tsv"
            a.hits, a.self_hits = a.out / "hits.tsv.gz", a.out / "self_hits.tsv"
            timed("reference_search", search, a.query_faa, a.db, a.hits, a.diamond, a.threads,
                  a.evalue, a.max_target_seqs, a.sensitivity)
            timed("self_comparison", search_self, a.query_faa, a.self_hits, a.diamond, a.threads, a.sensitivity)
            timed("scoring_and_reports", analyze, a, a.out)
            dump_json(a.out / "stage_timings.json",
                      {"clock": "time.perf_counter", "stages_seconds": stages,
                       "run_body_wall_seconds": time.perf_counter() - run_started,
                       "includes": "ORF preparation, searches, scoring, TSV and HTML/Sankey output",
                       "excludes": "interpreter startup, imports and database construction"})
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as e:
        print("ERROR: " + str(e), file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
