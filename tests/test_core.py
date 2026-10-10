import csv
import tempfile
import unittest
from collections import Counter, defaultdict
from pathlib import Path

from gtdb_hgt.cli import parser, validate
from gtdb_hgt.core import collect_hits, score_genes, smooth_and_regions
from gtdb_hgt.io import lineage, read_taxonomy, tax_key
from gtdb_hgt.storage import SpeciesHitStore


def tax(genus, species, family="F"):
    return lineage(f"d__Bacteria;p__P;c__C;o__O;f__{family};g__{genus};s__{genus} {species}")


def gene(i, contig="c", start=None):
    a = i * 100 if start is None else start
    return {"orf_id": f"q{i}", "contig": contig, "start": a, "end": a + 80, "gene": "", "product": ""}


class ScoringTests(unittest.TestCase):
    def setUp(self):
        self.qt = tax("A", "host")
        self.taxa = {"self": self.qt, "c1": tax("A", "a"), "c2": tax("A", "b"), "c3": tax("A", "c"),
                     "d1": tax("B", "a", "F2"), "d2": tax("B", "b", "F2"), "d3": tax("X", "a", "F3")}

    def write_hits(self, path):
        with path.open("w") as f:
            # Duplicate protein copies and multiple genomes of the same species
            # must not inflate the species vote.
            for g, sid, bits in [("self", "p", 500), ("c1", "p", 350), ("c1", "paralog", 100),
                                 ("c2", "p", 340), ("c3", "p", 330), ("d1", "p", 450),
                                 ("d2", "p", 430), ("d3", "p", 490)]:
                f.write(f"q1\t{g}::{sid}\t80\t100\t0\t0\t1\t100\t1\t100\t1e-50\t{bits}\t100\t100\n")

    def test_medians_use_species_and_supported_distal_bin(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "hits.tsv"
            self.write_hits(path)
            best, counts, _ = collect_hits(path, self.taxa, self.qt, "genus", "genus", "species")
            rows, _ = score_genes([gene(1)], best, counts, {"q1": 500}, self.taxa, self.qt)
            r = rows[0]
            self.assertEqual(r["close_hit_species"], 3)
            self.assertEqual(r["close_median_bitscore"], 340)
            self.assertEqual(r["donor_median_bitscore"], 440)
            self.assertEqual(r["donor_label"], "g__B")
            self.assertAlmostEqual(r["orf_score"], .2)
            self.assertTrue(r["eligible_for_window"])
            self.assertEqual(r["excluded_self"], 1)

    def test_disk_and_memory_aggregation_agree(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "hits.tsv"
            self.write_hits(path)
            best, counts, _ = collect_hits(path, self.taxa, self.qt, "genus", "genus", "species")
            expected, _ = score_genes([gene(1)], best, counts, {"q1": 500}, self.taxa, self.qt)
            store = SpeciesHitStore(Path(d) / "store.sqlite", self.taxa, "genus", "genus", self.qt)
            try:
                b, c, _ = collect_hits(path, self.taxa, self.qt, "genus", "genus", "species", disk_store=store)
                actual, _ = score_genes([gene(1)], b, c, {"q1": 500}, self.taxa, self.qt)
                self.assertEqual(actual, expected)
            finally:
                store.close()

    def test_close_absence_explicit_policy_and_no_hit_na(self):
        hits = defaultdict(dict)
        for g in ["d1", "d2"]:
            t = self.taxa[g]
            hits["q1"][tax_key(t, "species")] = {"bitscore": 450, "group": "distal", "donor_key": tax_key(t, "genus")}
        counts = defaultdict(Counter)
        rows, _ = score_genes([gene(1), gene(2)], hits, counts, {"q1": 500, "q2": 500}, self.taxa, self.qt)
        self.assertAlmostEqual(rows[0]["orf_score"], .9)
        self.assertIn("close_not_detected_above_filters", rows[0]["flags"])
        self.assertIsNone(rows[1]["orf_score"])
        self.assertFalse(rows[1]["eligible_for_window"])
        rows, _ = score_genes([gene(1)], hits, counts, {"q1": 500}, self.taxa, self.qt, missing_close="na")
        self.assertIsNone(rows[0]["orf_score"])

    def test_family_boundary_changes_group(self):
        # Same family, different genus: distal under genus; close under family.
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "hits.tsv"
            t = {"b": tax("B", "b")}
            path.write_text("q1\tb::p\t80\t100\t0\t0\t1\t100\t1\t100\t1e-30\t400\t100\t100\n")
            b, _, _ = collect_hits(path, t, self.qt, "genus", "genus", "species")
            self.assertEqual(next(iter(b["q1"].values()))["group"], "distal")
            b, _, _ = collect_hits(path, t, self.qt, "family", "genus", "species")
            self.assertEqual(next(iter(b["q1"].values()))["group"], "close")

    def test_accession_prefix_and_taxonomy_headers(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "tax.tsv"
            p.write_text("accession\tgtdb_taxonomy\nRS_GCF_000001.1\td__Bacteria;p__P;c__C;o__O;f__F;g__A;s__A x\n")
            self.assertIn("GCF_000001.1", read_taxonomy(p))


class WindowTests(unittest.TestCase):
    def row(self, i, score, contig="c", start=None):
        return {**gene(i, contig, start), "orf_score": score, "eligible_for_window": score is not None,
                "donor_taxon": "g__B", "flags": ""}

    def test_window_not_cross_contig_or_gap_and_na_not_zero(self):
        rows = [self.row(1, .5), self.row(2, None), self.row(3, .4), self.row(4, .3),
                self.row(5, .9, start=50000), self.row(6, .9, "different")]
        out, _ = smooth_and_regions(rows, window=5, min_valid=3)
        first = next(r for r in out if r["orf_id"] == "q2")
        self.assertAlmostEqual(first["window_score"], .4)
        self.assertEqual(first["window_valid_orfs"], 3)
        for orf in ["q5", "q6"]:
            self.assertIsNone(next(r for r in out if r["orf_id"] == orf)["window_score"])

    def test_isolated_positive_not_region_and_consecutive_signal_is(self):
        values = [-.4, -.3, .5, -.35, -.25]
        rows, regions = smooth_and_regions([self.row(i+1, v) for i,v in enumerate(values)])
        self.assertAlmostEqual(rows[2]["window_score"], -.16)
        self.assertEqual(regions, [])
        rows, regions = smooth_and_regions([self.row(i+1, v) for i,v in enumerate([.5,.4,.02,.45,.35])])
        self.assertAlmostEqual(rows[2]["window_score"], .344)
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["n_orfs"], 5)

    def test_invalid_window_rejected(self):
        with self.assertRaises(ValueError):
            smooth_and_regions([], window=4)


class DefaultTests(unittest.TestCase):
    def test_user_defaults(self):
        a = parser().parse_args(["run", "--genbank", "x.gbk", "--db", "db.dmnd", "--taxonomy", "tax.tsv",
                                 "--query-genome", "GCA_1.1", "--out", "run"])
        validate(a)
        self.assertEqual(a.close_rank, "genus")
        self.assertEqual(a.window, 5)
        self.assertEqual(a.max_target_seqs, 0)


if __name__ == "__main__":
    unittest.main()
