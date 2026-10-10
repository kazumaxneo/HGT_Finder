import csv
import json
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from gtdb_hgt.io import lineage
from gtdb_hgt.sankey import assign_orfs, attach_regions, contig_summary, write_sankey, assign_from_genus_summaries, replot_from_run, assign_cached_by_rank, assign_orfs_by_rank, DISPLAY_RANKS
from gtdb_hgt.cli import parser, validate


def tax(genus, species):
    return lineage(f"d__Bacteria;p__P;c__C;o__O;f__F;g__{genus};s__{genus} {species}")


def gene(q, contig="c"):
    return dict(orf_id=q, contig=contig, start=1, end=100, gene="", product="")


class SankeyTests(unittest.TestCase):
    def test_uses_close_and_far_medians_not_far_donor_label(self):
        hits = {"q": {"a": {"tax": tax("Host","a"), "bitscore": 450},
                      "b": {"tax": tax("Host","b"), "bitscore": 430},
                      "c": {"tax": tax("Far","a"), "bitscore": 340},
                      "d": {"tax": tax("Far","b"), "bitscore": 360}}}
        row = assign_orfs([gene("q")], hits, tax("Host","query"))[0]
        self.assertEqual(row["taxon_label"], "g__Host")
        self.assertTrue(row["host_taxon_match"])
        self.assertEqual(row["support_species"], 2)
        self.assertEqual(row["median_bitscore"], 440)

    def test_ambiguous_sparse_and_unmapped_stay_separate(self):
        hits = defaultdict(dict)
        for key, bits in [("A",400),("B",390)]:
            for i in range(2):
                hits["q1"][f"{key}{i}"] = {"tax":tax(key,str(i)), "bitscore":bits}
        hits["q2"]["only"] = {"tax":tax("A","a"),"bitscore":600}
        rows = assign_orfs([gene("q1"),gene("q2"),gene("q3")],hits,tax("A","host"),method="median",margin=.05)
        self.assertEqual([r["assignment_status"] for r in rows],
                         ["ambiguous","insufficient_support","no_usable_hit"])
        self.assertTrue(all(not r["taxon"] for r in rows))
        self.assertEqual(rows[0]["relative_margin"], .025)
        coarse = assign_orfs([gene("q1")],hits,tax("A","host"),rank="family")
        self.assertEqual(coarse[0]["assignment_status"],"assigned")
        self.assertEqual(coarse[0]["support_species"],4)

    def hits(self, values):
        return {"q": {f"{genus}{i}": {"tax": tax(genus,str(i)), "bitscore": bits}
                      for genus, scores in values.items() for i, bits in enumerate(scores)}}

    def test_best_hit_recovers_host_hidden_by_median(self):
        hits = self.hits({"Host": [100,60,60], "Far": [90,90]})
        host = tax("Host","query")
        best = assign_orfs([gene("q")], hits, host)[0]
        legacy = assign_orfs([gene("q")], hits, host, method="median")[0]
        self.assertEqual(best["taxon_label"], "g__Host")
        self.assertEqual(legacy["taxon_label"], "g__Far")
        # Changing the host label must not change the selected strongest hit.
        other = assign_orfs([gene("q")], hits, tax("Far","query"))[0]
        self.assertEqual(other["taxon"], best["taxon"])
        self.assertFalse(other["host_taxon_match"])

    def test_near_alternative_retained_and_exact_tie_unresolved(self):
        host = tax("Host","query")
        row = assign_orfs([gene("q")], self.hits({"Host":[100,80], "Far":[99,70]}),host)[0]
        self.assertEqual(row["assignment_status"], "assigned")
        self.assertTrue(row["near_tie"])
        self.assertIn("g__Far", row["nearby_taxa"])
        tied = assign_orfs([gene("q")], self.hits({"Host":[100,80], "Far":[100]}),host)[0]
        self.assertEqual(tied["assignment_status"], "ambiguous")
        self.assertEqual(tied["ambiguity_kind"], "host_competes")
        self.assertEqual(tied["taxon"], "")

    def test_one_percent_boundary_and_change_from_five_percent(self):
        host = tax("Host","query")
        hits = self.hits({"Host":[100,80], "Far":[98,70]})
        new = assign_orfs([gene("q")], hits, host)[0]
        old = assign_orfs([gene("q")], hits, host, margin=.05)[0]
        self.assertFalse(new["near_tie"])
        self.assertTrue(old["near_tie"])
        self.assertEqual(new["taxon"], old["taxon"])
        at_boundary = assign_orfs([gene("q")], self.hits({"Host":[100,80], "Far":[99,70]}),host)[0]
        self.assertTrue(at_boundary["near_tie"])
        self.assertEqual(at_boundary["relative_margin"], .01)
        self.assertEqual(parser().parse_args(["sankey","--run-dir","x","--out","y"]).sankey_margin,.01)

    def test_sparse_strongest_never_replaced_by_weaker_supported_bin(self):
        row = assign_orfs([gene("q")], self.hits({"Host":[100], "Far":[90,80]}),tax("Host","q"))[0]
        self.assertEqual(row["assignment_status"], "insufficient_support")
        self.assertEqual(row["best_bitscore"],100)
        self.assertEqual(row["support_species"],1)

    def test_cached_replot_matches_species_hits_and_preserves_scores(self):
        host=tax("Host","query")
        hits=self.hits({"Host":[100,60,60], "Far":[99,99], "Other":[100,50]})
        genes=[gene("q"),gene("nohit")]
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/"source"; source.mkdir()
            summaries=source/"taxon_summaries.tsv"
            with summaries.open("w") as f:
                writer=csv.DictWriter(f,fieldnames=["orf_id","taxon","species_count","median_bitscore","max_bitscore"],delimiter="\t")
                writer.writeheader()
                for genus,med,best,n in [("Host",60,100,3),("Far",99,99,2),("Other",75,100,2)]:
                    writer.writerow(dict(orf_id="q",taxon=";".join(list(tax(genus,"x").values())[:6]),species_count=n,median_bitscore=med,max_bitscore=best))
            for method in ("best-hit", "median"):
                self.assertEqual(assign_from_genus_summaries(genes,summaries,host,method=method),
                                 assign_orfs(genes,hits,host,method=method))
            (source/"run.json").write_text(json.dumps(dict(database_label="test",query_taxonomy=host,parameters=dict(close_rank="genus",donor_rank="genus"))))
            with (source/"orf_scores.tsv").open("w") as f:
                fields=list(genes[0])+["orf_score","window_score","flags"]
                writer=csv.DictWriter(f,fieldnames=fields,delimiter="\t");writer.writeheader()
                for g in genes: writer.writerow({**g,"orf_score":"","window_score":"0.2","flags":"test"})
            (source/"regions.tsv").write_text("region_id\tcontig\tstart\tend\tn_orfs\torf_ids\nr1\tc\t1\t100\t1\tq\n")
            out=Path(d)/"plot"
            result=replot_from_run(source,out)
            self.assertEqual(result["total_orfs"],2)
            self.assertEqual(result["candidate_orfs"],1)
            self.assertFalse(result["hgt_scores_recomputed"])
            for name in ("orf_scores.tsv","regions.tsv"):
                self.assertEqual((source/name).read_bytes(),(out/name).read_bytes())
            with self.assertRaises(ValueError): replot_from_run(source,out)

    def test_rank_tie_resolution_and_species_support_are_recomputed(self):
        host=tax("Host","q")
        hits=self.hits({"Host":[100], "Far":[100]})
        ranked=assign_orfs_by_rank([gene("q")],hits,host)
        self.assertEqual(ranked["genus"][0]["assignment_status"],"ambiguous")
        family=ranked["family"][0]
        self.assertEqual(family["assignment_status"],"assigned")
        self.assertEqual(family["support_species"],2)
        self.assertTrue(family["host_taxon_match"])
        self.assertIsNone(family["relative_margin"])
        # A tie spanning families still resolves at their shared order.
        hits["q"]["Far0"]["tax"]["family"]="f__Other"
        ranked=assign_orfs_by_rank([gene("q")],hits,host)
        self.assertEqual(ranked["family"][0]["assignment_status"],"ambiguous")
        self.assertEqual(ranked["order"][0]["assignment_status"],"assigned")

    def test_cached_higher_rank_matches_species_hits_without_fake_medians(self):
        host=tax("Host","q")
        genes=[gene("q"),gene("missing")]
        hits=self.hits({"Host":[100,70], "Far":[100], "Other":[99,70]})
        for hit in hits["q"].values():
            if hit["tax"]["genus"]=="g__Other":hit["tax"]["family"]="f__Other"
        with tempfile.TemporaryDirectory() as d:
            cache=Path(d)/"cache.tsv"
            cache.write_text("orf_id\ttaxon\tspecies_count\tmedian_bitscore\tmax_bitscore\n"+
                "q\t"+";".join(list(host.values())[:6])+"\t2\t85\t100\n"+
                "q\t"+";".join(list(tax("Far","a").values())[:6])+"\t1\t100\t100\n"+
                "q\td__Bacteria;p__P;c__C;o__O;f__Other;g__Other\t2\t84.5\t99\n")
            cached=assign_cached_by_rank(genes,cache,host)
            raw=assign_orfs_by_rank(genes,hits,host)
            self.assertEqual(cached["genus"],raw["genus"])
            for rank in DISPLAY_RANKS[:-1]:
                for a,b in zip(cached[rank],raw[rank]):
                    self.assertIsNone(a["median_bitscore"])
                    self.assertIsNone(a["host_median_bitscore"])
                    for key in a:
                        if key not in ("median_bitscore","host_median_bitscore"):
                            self.assertEqual(a[key],b[key],(rank,key))
            with self.assertRaises(ValueError):assign_cached_by_rank(genes,cache,host,method="median")

    def test_regions_use_emitted_regions_not_isolated_positive_window(self):
        rows = assign_orfs([gene("q1"),gene("q2")],{},tax("A","host"))
        scores = [{"orf_id":q,"orf_score":None,"window_score":.2,"flags":"","region_candidate":True} for q in ["q1","q2"]]
        joined = attach_regions(rows,scores,[{"region_id":"r1","orf_ids":"q1"}])
        self.assertTrue(joined[0]["in_candidate_region"])
        self.assertFalse(joined[1]["in_candidate_region"])
        summary = contig_summary(joined)
        self.assertEqual(sum(r["orf_count"] for r in summary),2)
        self.assertEqual(sum(r["fraction_of_contig_orfs"] for r in summary),1)
        self.assertEqual(sum(r["candidate_region_orfs"] for r in summary),1)

    def test_html_annotation_cannot_close_script(self):
        g = gene("q")
        g["product"] = "</script><script>unsafe()</script>"
        rows = attach_regions(assign_orfs([g],{},tax("A","host")),
                              [{"orf_id":"q","orf_score":None,"window_score":None,"flags":""}],[])
        manifest = {"database_label":"demo <unsafe>",
                    "sankey":{"rank":"genus","min_species":2,"margin":.05,"max_contigs":20,"max_taxa":12}}
        with tempfile.TemporaryDirectory() as d:
            write_sankey(Path(d),rows,manifest)
            page = (Path(d)/"sankey.html").read_text()
            self.assertEqual(page.count("</script>"),1)
            self.assertNotIn("<script>unsafe()",page)
            self.assertIn("demo &lt;unsafe&gt;",page)

    def test_cli_opt_in_and_validation(self):
        base=["analyze","--hits","x","--self-hits","x","--genes","x","--taxonomy","x","--query-genome","x","--out","x"]
        self.assertFalse(parser().parse_args(base).sankey)
        self.assertEqual(parser().parse_args(base).sankey_method,"best-hit")
        with self.assertRaises(ValueError):
            validate(parser().parse_args(["sankey","--run-dir","x","--out","y","--sankey-margin","nan"]))
        a=parser().parse_args(base+["--sankey","--sankey-rank","family"])
        validate(a)
        self.assertEqual(a.sankey_rank,"family")
        with self.assertRaises(ValueError):
            validate(parser().parse_args(base+["--sankey-margin","nan"]))


if __name__ == "__main__":
    unittest.main()
