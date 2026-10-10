import csv
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

from gtdb_hgt.cli import main, parser, validate
from gtdb_hgt.distribution import SavedHitReader, plot_from_run, summarize_points
from gtdb_hgt.io import digest, lineage, tax_key, tsv_write, write_taxonomy
from gtdb_hgt.storage import SpeciesHitStore


class DistributionTests(unittest.TestCase):
    def fixture(self, root, legacy=False):
        source = root / "source"
        source.mkdir()
        host = lineage("d__Bacteria;p__P;c__C;o__O;f__F;g__Host;s__Host query")
        taxa = {f"GCA_{i:09d}.1": lineage(f"d__Bacteria;p__P;c__C;o__O;f__F;g__{genus};s__{genus} {i}")
                for i, genus in enumerate(("Host", "Host", "Far", "Far", "Far"))}
        taxonomy = root / "taxonomy.tsv"
        write_taxonomy(taxonomy, taxa)
        points = []
        for i, (genome, tax) in enumerate(taxa.items()):
            points.append(dict(genome_id=genome, tax=tax, sseqid=genome + "::protein", bitscore=(400,420,600,100,110)[i],
                               qcov=1., scov=1., pident=(80,80,30,50,80)[i], evalue=1e-30))
        if legacy:
            db = sqlite3.connect(source / "species_hits.sqlite")
            db.execute("CREATE TABLE hits(q TEXT,species INTEGER,genome INTEGER,subject TEXT,bits REAL,qcov REAL,scov REAL,PRIMARY KEY(q,species))")
            db.executemany("INSERT INTO hits VALUES(?,?,?,?,?,?,?)", [("q",i,i,p["sseqid"],p["bitscore"],1,1) for i,p in enumerate(points)])
            db.commit(); db.close()
        else:
            store = SpeciesHitStore(source / "species_hits.sqlite", taxa, "genus", "genus", host)
            for p in points:
                store.put("q", tax_key(p["tax"], "species"), p)
            # A lower paralog of the same species must not create a second point.
            store.put("q", tax_key(points[0]["tax"], "species"), {**points[0], "bitscore": 50, "pident": 20})
            store.close()
        manifest = dict(query_taxonomy=host, inputs=dict(taxonomy=str(taxonomy), taxonomy_sha256=digest(taxonomy)))
        (source / "run.json").write_text(json.dumps(manifest))
        tsv_write(source / "orf_scores.tsv", [dict(orf_id=q,contig="c",start=i*100+1,end=i*100+100,
                   gene="nifH" if q=="q" else "",product="test </script><script>bad</script>",orf_score=.2)
                   for i,q in enumerate(("q","nohit"))])
        return source, taxonomy, manifest

    def test_species_points_exact_medians_and_no_score_mutations(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source, taxonomy, manifest=self.fixture(root)
            original=(source/"orf_scores.tsv").read_bytes()
            reader=SavedHitReader(source,manifest)
            points=reader.get("q");reader.close()
            self.assertEqual(len(points),5)
            groups=summarize_points(points,manifest["query_taxonomy"],"genus")
            self.assertEqual([(g["species_count"],g["median_bitscore"],g["max_bitscore"]) for g in groups],[(2,410,420),(3,110,600)])
            self.assertEqual(groups[1]["below_quality_floor"],1)
            self.assertEqual(sum(g["species_count"] for g in summarize_points(points,manifest["query_taxonomy"],"phylum")),5)
            # Embedded references make distribution generation independent of the input taxonomy path.
            taxonomy.unlink()
            result=plot_from_run(source,root/"out",["q","nohit"])
            self.assertEqual(result["species_hits"],5)
            self.assertEqual(result["entries"][1]["species_hits"],0)
            page=(root/"out/orf_000001.html").read_text()
            self.assertNotIn('</script><script>bad</script>',page)
            self.assertIn('data-rank="phylum"',page)
            self.assertIn('data-rank="genus"',page)
            self.assertEqual((source/"orf_scores.tsv").read_bytes(),original)
            self.assertFalse(result["hgt_scores_recomputed"])
            with (root/"out/orf_000001_hits.tsv").open() as f:
                self.assertEqual(len(list(csv.DictReader(f,delimiter="\t"))),5)

    def test_legacy_encoded_ids_and_missing_quality(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source,taxonomy,manifest=self.fixture(root,legacy=True)
            reader=SavedHitReader(source,manifest)
            self.assertTrue(all(p["quality"] is None for p in reader.get("q")))
            reader.close()
            relocated=root/"relocated.tsv";taxonomy.rename(relocated)
            plot_from_run(source,root/"out",["q"],taxonomy=relocated)
            relocated.write_text(relocated.read_text().replace("g__Far","g__Different"))
            with self.assertRaisesRegex(ValueError,"Taxonomy hash"):
                SavedHitReader(source,manifest,taxonomy=relocated)

    def test_missing_store_unknown_orf_and_invalid_display_options(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source,_,manifest=self.fixture(root)
            with self.assertRaisesRegex(ValueError,"Unknown ORF"):
                plot_from_run(source,root/"unknown",["absent"])
            self.assertFalse((root/"unknown").exists())
            (source/"species_hits.sqlite").unlink()
            with self.assertRaisesRegex(ValueError,"medians cannot reconstruct"):
                SavedHitReader(source,manifest)
        base=["distribution","--run-dir","source","--out","out","--orf-id","q"]
        self.assertEqual(parser().parse_args(base).min_hit_quality,.5)
        for tail in (["--min-hit-quality","nan"],["--min-hit-quality","1.1"],["--max-taxa","0"]):
            with self.assertRaises(ValueError):validate(parser().parse_args(base+tail))

    def test_analyze_generates_links_and_keeps_comparison_hits(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/"taxonomy.tsv").write_text('g1\td__Bacteria;p__P;c__C;o__O;f__F;g__Host;s__Host 1\ng2\td__Bacteria;p__P;c__C;o__O;f__F;g__Far;s__Far 1\n')
            (root/"genes.tsv").write_text('orf_id\tcontig\tstart\tend\tgene\tproduct\nq\tc\t1\t100\tnifH\ttest\n')
            (root/"hits.tsv").write_text('q\tg1::a\t80\t100\t0\t0\t1\t100\t1\t100\t1e-20\t400\t100\t100\nq\tg2::a\t30\t100\t0\t0\t1\t100\t1\t100\t1e-20\t600\t100\t100\n')
            (root/"self.tsv").write_text('q\tq\t100\t100\t0\t0\t1\t100\t1\t100\t0\t700\t100\t100\n')
            main(['analyze','--hits',str(root/'hits.tsv'),'--self-hits',str(root/'self.tsv'),
                  '--genes',str(root/'genes.tsv'),'--taxonomy',str(root/'taxonomy.tsv'),
                  '--query-taxonomy','d__Bacteria;p__P;c__C;o__O;f__F;g__Host;s__Host query',
                  '--out',str(root/'out'),'--sankey','--distribution-orf','q'])
            index=root/'out/hit_distributions/index.html'
            self.assertTrue(index.is_file())
            report=(root/'out/report.html').read_text()
            self.assertIn('hit_distributions/index.html?orf=q',report)
            data=json.loads((root/'out/hit_distributions/distribution.json').read_text())
            self.assertEqual(data['species_hits'],2)
            self.assertEqual(data['entries'][0]['quality_known'],2)
            manifest=json.loads((root/'out/run.json').read_text())
            self.assertEqual(manifest['sankey']['distribution_orfs'],['q'])
            from gtdb_hgt.sankey import replot_from_run
            replot_from_run(root/'out', root/'replot')
            self.assertTrue((root/'replot/hit_distributions/orf_000001.html').is_file())
            self.assertIn('"distribution_index": "hit_distributions/index.html"',(root/'replot/sankey.html').read_text())


if __name__=='__main__':unittest.main()
