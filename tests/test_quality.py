import math
import tempfile
import unittest
from pathlib import Path

from gtdb_hgt.cli import parser, validate
from gtdb_hgt.core import collect_hits, score_genes, smooth_and_regions
from gtdb_hgt.io import lineage, tax_key
from gtdb_hgt.quality import alignment_quality
from gtdb_hgt.sankey import assign_orfs_by_rank, attach_regions, DISPLAY_RANKS
from gtdb_hgt.storage import SpeciesHitStore


def tax(genus, species):
    return lineage(f"d__Bacteria;p__P;c__C;o__O;f__F;g__{genus};s__{genus} {species}")


class QualityTests(unittest.TestCase):
    def setUp(self):
        self.host = tax("Host", "self")
        self.taxa = {"self": self.host, **{f"c{i}": tax("Host", str(i)) for i in range(3)},
                     **{f"d{i}": tax("Far", str(i)) for i in range(2)}}
        self.gene = dict(orf_id="q", contig="ctg", start=100, end=399, gene="", product="")

    def hit(self, genome, bits, identity, qcov=1, scov=1):
        t = self.taxa[genome]
        return dict(genome_id=genome, tax=t, bitscore=bits, pident=identity,
                    qcov=qcov, scov=scov, evalue=1e-50, sseqid=genome+"::p",
                    group="close" if genome.startswith("c") else "distal", donor_key=tax_key(t,"genus"))

    def rows(self, identity=80, high_host=False, min_quality=.5):
        from collections import Counter, defaultdict
        hits = {"q": {f"c{i}": self.hit(f"c{i}", 500 if high_host and i==0 else 350, 95 if high_host and i==0 else 30)
                      for i in range(3)}}
        hits["q"].update({f"d{i}": self.hit(f"d{i}", 450, identity) for i in range(2)})
        scored, _ = score_genes([self.gene], hits, defaultdict(Counter), {"q":600}, self.taxa,
                                self.host, min_hit_quality=min_quality)
        return hits, scored

    def test_identity_and_both_coverages_affect_quality(self):
        self.assertAlmostEqual(alignment_quality(90,1,1), .9)
        self.assertAlmostEqual(alignment_quality(60,.98,.98), .57624)
        self.assertAlmostEqual(alignment_quality(33.1,.98,.97), .3146486)
        self.assertAlmostEqual(alignment_quality(94.2,1,.544), .512448)
        self.assertAlmostEqual(alignment_quality(90,.5,.5), .225)
        for args in [(101,1,1),(90,1.01,1),(90,math.nan,1),(-1,1,1)]:
            with self.assertRaises(ValueError):alignment_quality(*args)

    def test_weak_top_is_na_not_zero_and_hits_remain_for_comparison(self):
        hits, scored = self.rows(33.1)
        row = scored[0]
        self.assertAlmostEqual(row["raw_orf_score"], 100/600)
        self.assertIsNone(row["orf_score"])
        self.assertEqual(row["quality_status"], "low_top_hit_quality")
        self.assertFalse(row["eligible_for_window"])
        self.assertEqual(row["close_hit_species"],3)
        self.assertEqual(row["distal_hit_species"],2)
        smooth, regions=smooth_and_regions(scored,window=1,min_valid=1,min_region_orfs=1)
        ranked=assign_orfs_by_rank([self.gene],hits,self.host)
        for rank in DISPLAY_RANKS:
            plotted=attach_regions(ranked[rank],scored,[])[0]
            self.assertEqual(plotted["assignment_status"],"low_hit_quality")
            self.assertEqual(plotted["taxon"],"")
        self.assertIsNone(smooth[0]["window_score"])
        self.assertEqual(regions,[])

    def test_weak_donor_not_rescued_by_strong_host_top(self):
        _, scored=self.rows(33.1,high_host=True)
        r=scored[0]
        self.assertAlmostEqual(r["top_hit_quality"],.95)
        self.assertEqual(r["quality_status"],"low_donor_hit_quality")
        self.assertIsNone(r["orf_score"])

    def test_strong_hit_passes_and_zero_restores_legacy_score(self):
        _, strong=self.rows(80)
        self.assertEqual(strong[0]["quality_status"],"pass")
        self.assertTrue(strong[0]["eligible_for_window"])
        self.assertAlmostEqual(strong[0]["orf_score"],100/600)
        _, legacy=self.rows(33.1,min_quality=0)
        self.assertAlmostEqual(legacy[0]["orf_score"],100/600)
        self.assertEqual(legacy[0]["quality_status"],"not_applied")
        # Lower-bound equality is retained, not excluded.
        _, boundary=self.rows(50)
        self.assertEqual(boundary[0]["quality_status"],"pass")

    def test_sqlite_retains_quality_of_selected_species_best_hit(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"hits.tsv"
            path.write_text("q\td0::weak\t33.1\t100\t0\t0\t1\t100\t1\t100\t1e-50\t450\t100\t100\n"
                            "q\td0::stronger_identity\t90\t100\t0\t0\t1\t100\t1\t100\t1e-40\t440\t100\t100\n")
            mem, counts, _=collect_hits(path,self.taxa,self.host,"genus","genus","species")
            store=SpeciesHitStore(Path(d)/"hits.sqlite",self.taxa,"genus","genus",self.host)
            try:
                disk, dc, _=collect_hits(path,self.taxa,self.host,"genus","genus","species",disk_store=store)
                a,_=score_genes([self.gene],mem,counts,{"q":600},self.taxa,self.host,min_hit_quality=.5)
                b,_=score_genes([self.gene],disk,dc,{"q":600},self.taxa,self.host,min_hit_quality=.5)
                self.assertEqual(a,b)
                self.assertAlmostEqual(b[0]["top_hit_quality"],.331)
                self.assertEqual(b[0]["top_hit_subject"],"d0::weak")
            finally:store.close()

    def test_cli_floor_is_provisional_and_invalid_values_rejected(self):
        base=["analyze","--hits","x","--self-hits","x","--genes","x","--taxonomy","x","--query-genome","x","--out","x"]
        self.assertEqual(parser().parse_args(base).min_hit_quality,.5)
        for value in ["nan","-1","1.1"]:
            with self.assertRaises(ValueError):validate(parser().parse_args(base+["--min-hit-quality",value]))

    def test_analysis_and_cached_replot_preserve_quality_across_all_ranks(self):
        import contextlib, csv, io, json
        from gtdb_hgt.cli import main
        from gtdb_hgt.io import write_taxonomy
        from gtdb_hgt.sankey import replot_from_run
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            write_taxonomy(root/"tax.tsv",self.taxa)
            (root/"genes.tsv").write_text("orf_id\tcontig\tstart\tend\nq1\tctg\t1\t300\nq2\tctg\t400\t699\n")
            (root/"self.tsv").write_text("".join(f"{q}\t{q}\t100\t100\t0\t0\t1\t100\t1\t100\t0\t600\t100\t100\n" for q in ("q1","q2")))
            with (root/"hits.tsv").open("w") as f:
                for q,identity in [("q1",80),("q2",33.1)]:
                    for g in ["c0","c1","c2","d0","d1"]:
                        bits=350 if g.startswith("c") else 450
                        f.write(f"{q}\t{g}::p\t{identity}\t100\t0\t0\t1\t100\t1\t100\t1e-50\t{bits}\t100\t100\n")
            args=["analyze","--hits",str(root/"hits.tsv"),"--self-hits",str(root/"self.tsv"),
                  "--genes",str(root/"genes.tsv"),"--taxonomy",str(root/"tax.tsv"),"--query-genome","self",
                  "--out",str(root/"out"),"--assume-unlimited-hits","--sankey",
                  "--window","1","--min-valid","1","--min-region-orfs","1"]
            with contextlib.redirect_stdout(io.StringIO()):main(args)
            run=json.loads((root/"out/run.json").read_text())
            self.assertEqual(run["counts"]["quality_excluded_orfs"],1)
            self.assertEqual(run["counts"]["window_eligible_orfs"],1)
            for rank in DISPLAY_RANKS:
                with (root/f"out/sankey_orfs_{rank}.tsv").open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
                self.assertEqual(rows[1]["assignment_status"],"low_hit_quality")
                self.assertIn("unfiltered_taxon",rows[1])
                self.assertNotEqual(rows[1]["unfiltered_taxon"],"")
            replot_from_run(root/"out",root/"replot")
            self.assertEqual((root/"out/orf_scores.tsv").read_bytes(),(root/"replot/orf_scores.tsv").read_bytes())
            for rank in DISPLAY_RANKS:
                with (root/f"replot/sankey_orfs_{rank}.tsv").open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
                self.assertEqual(rows[1]["assignment_status"],"low_hit_quality")
                self.assertAlmostEqual(float(rows[1]["top_hit_quality"]),.331)
            self.assertIn("品質の下限",(root/"out/sankey.html").read_text())


if __name__ == "__main__":unittest.main()
