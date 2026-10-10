import csv
import json
import tempfile
import unittest
from pathlib import Path

from gtdb_hgt.cli import parser, validate
from gtdb_hgt.io import digest, lineage, tsv_write
from gtdb_hgt.sankey import DISPLAY_RANKS, assign_orfs_by_rank, attach_regions, replot_from_run, write_sankey


class SankeyQualityOverlayTests(unittest.TestCase):
    def fixture(self, root):
        source = root / "source"
        source.mkdir()
        host = lineage("d__Bacteria;p__P;c__C;o__O;f__F;g__Host;s__Host q")
        genes = [dict(orf_id=q, contig="c", start=i * 100 + 1, end=i * 100 + 100,
                      gene="", product="") for i, q in enumerate(("weak", "boundary", "strong", "nohit"))]
        hits = {g["orf_id"]: {str(i): dict(tax=lineage(
            f"d__Bacteria;p__P;c__C;o__O;f__F;g__Far;s__Far {i}"), bitscore=100)
            for i in range(2)} for g in genes[:3]}
        scores = [{**g, "orf_score": .2, "window_score": .25, "flags": "original"} for g in genes]
        regions = [dict(region_id="r1", contig="c", start=1, end=400, n_orfs=4,
                        orf_ids="weak;boundary;strong;nohit")]
        manifest = dict(database_label="synthetic", query_taxonomy=host,
                        parameters=dict(close_rank="genus", donor_rank="genus"),
                        sankey=dict(rank="genus", method="best-hit", min_species=2, margin=.01,
                                    max_contigs=20, max_taxa=12))
        ranked = {rank: attach_regions(rows, scores, regions) for rank, rows in
                  assign_orfs_by_rank(genes, hits, host).items()}
        tsv_write(source / "orf_scores.tsv", scores)
        tsv_write(source / "regions.tsv", regions)
        (source / "report.html").write_text("original report")
        write_sankey(source, ranked["genus"], manifest, regions, ranked)
        (source / "run.json").write_text(json.dumps(manifest))
        audit = root / "audit"
        audit.mkdir()
        audit_rows = []
        for g, quality in zip(genes, (.3, .5, .8, None)):
            audit_rows.append({**g, "subject": "subject" if quality is not None else "",
                "pident": quality * 100 if quality is not None else "", "qcov": 1 if quality is not None else "",
                "scov": 1 if quality is not None else "", "evalue": 1e-30 if quality is not None else "",
                "quality": quality if quality is not None else "",
                "genus_status": "assigned" if quality is not None else "no_usable_hit",
                "phylum_status": "assigned" if quality is not None else "no_usable_hit"})
        tsv_write(audit / "top_hit_quality.tsv", audit_rows)
        (audit / "summary.json").write_text(json.dumps(dict(complete=True, expected_hits=3, recovered_hits=3,
            input_hashes={f"sankey_orfs_{rank}.tsv": digest(source / f"sankey_orfs_{rank}.tsv")
                          for rank in ("genus", "phylum")})))
        return source, audit / "top_hit_quality.tsv", audit_rows

    def test_audit_masks_all_ranks_and_keeps_original_scores(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source, audit, _ = self.fixture(root)
            out = root / "plot"
            result = replot_from_run(source, out, min_hit_quality=.5, quality_audit=audit)
            self.assertEqual(result["min_hit_quality"], .5)
            for rank in DISPLAY_RANKS:
                with (out / f"sankey_orfs_{rank}.tsv").open() as f:
                    rows = list(csv.DictReader(f, delimiter="\t"))
                self.assertEqual([r["assignment_status"] for r in rows],
                                 ["low_hit_quality", "assigned", "assigned", "no_usable_hit"])
                self.assertEqual(rows[0]["top_hit_quality"], "0.3")
                self.assertEqual(rows[0]["orf_score"], "0.2")
            for name in ("orf_scores.tsv", "regions.tsv", "report.html"):
                self.assertEqual((out / name).read_bytes(), (source / name).read_bytes())
            manifest = json.loads((out / "run.json").read_text())
            self.assertNotIn("min_hit_quality", manifest["parameters"])
            self.assertTrue(manifest["sankey"]["quality_display_only"])
            self.assertFalse(manifest["sankey_replot"]["hgt_scores_recomputed"])
            # Lowering the floor on a saved plot restores the original taxonomy
            # and booleans, even though its original score TSV lacks quality.
            restored = replot_from_run(out, root / "restored", min_hit_quality=0)
            self.assertEqual(restored["status_counts"], {"assigned": 3, "no_usable_hit": 1})
            self.assertEqual(restored["rank_counts"]["genus"]["other_assigned"], 3)
            raised = replot_from_run(out, root / "raised", min_hit_quality=.8)
            self.assertEqual(raised["status_counts"]["low_hit_quality"], 2)
            self.assertEqual(raised["status_counts"]["assigned"], 1)

    def test_incomplete_mismatched_or_invalid_audit_is_rejected(self):
        for corruption in ("incomplete", "hash", "duplicate", "coordinates", "quality"):
            with self.subTest(corruption=corruption), tempfile.TemporaryDirectory() as d:
                root = Path(d)
                source, audit, rows = self.fixture(root)
                metadata = json.loads(audit.with_name("summary.json").read_text())
                if corruption == "incomplete":
                    metadata["complete"] = False
                elif corruption == "hash":
                    metadata["input_hashes"]["sankey_orfs_genus.tsv"] = "wrong"
                elif corruption == "duplicate":
                    rows.append(rows[0])
                elif corruption == "coordinates":
                    rows[0]["start"] = 2
                else:
                    rows[0]["quality"] = .9
                tsv_write(audit, rows)
                audit.with_name("summary.json").write_text(json.dumps(metadata))
                with self.assertRaises(ValueError):
                    replot_from_run(source, root / "out", min_hit_quality=.5, quality_audit=audit)

    def test_cli_default_and_missing_quality(self):
        base = ["sankey", "--run-dir", "source", "--out", "out"]
        self.assertEqual(parser().parse_args(base).min_hit_quality, .5)
        commands = next(a for a in parser()._actions if getattr(a, "choices", None))
        self.assertIn("default 0.5", commands.choices["sankey"].format_help())
        for value in ("nan", "inf", "-0.1", "1.1"):
            with self.assertRaises(ValueError):
                validate(parser().parse_args(base + ["--min-hit-quality", value]))
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source, _, _ = self.fixture(root)
            with self.assertRaisesRegex(ValueError, "lacks top-hit quality"):
                replot_from_run(source, root / "out", min_hit_quality=.5)


if __name__ == "__main__":
    unittest.main()
