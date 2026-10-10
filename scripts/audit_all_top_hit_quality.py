"""Recover identity for stored species-best hits; do not rerun HGT or search.

Read-only inputs: original species SQLite, saved rank assignments and raw hits.
Outputs are a new audit directory. Bit-score winners remain unchanged.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
from statistics import median
import time


def read_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def quality(identity, qc, sc):
    if not (0 <= identity <= 100 and 0 <= qc <= 1 and 0 <= sc <= 1):
        raise ValueError("Invalid quality inputs")
    return identity / 100 * qc * sc


def summarize(rows):
    valid = [r for r in rows if r.get("quality") is not None]
    scores = [r["quality"] for r in valid]
    bins = Counter()
    for s in scores:
        key = "<=0.3" if s <= .3 else "(0.3,0.4)" if s < .4 else "[0.4,0.5)" if s < .5 else "[0.5,0.6)" if s < .6 else ">=0.6"
        bins[key] += 1
    fields = ("orf_id", "product", "pident", "qcov", "scov", "quality", "genus_label", "phylum_label")
    sort = sorted(valid, key=lambda r: (r["quality"], r["orf_id"]))
    return dict(orfs=len(rows), recovered=len(valid), missing=len(rows)-len(valid),
                min=min(scores) if scores else None, median=median(scores) if scores else None,
                max=max(scores) if scores else None,
                bins={key:bins[key] for key in ("<=0.3", "(0.3,0.4)", "[0.4,0.5)", "[0.5,0.6)", ">=0.6")},
                retained_at={str(t):sum(s >= t for s in scores) for t in (.3,.4,.5,.55,.6)},
                lowest=[{k:r[k] for k in fields} for r in sort[:8]],
                highest=[{k:r[k] for k in fields} for r in sort[-8:]])


def audit(run, assignments, out):
    if out.exists():
        raise ValueError("Use a new audit output directory")
    out.mkdir(parents=True)
    started=time.perf_counter()
    genus=read_rows(assignments/"sankey_orfs_genus.tsv")
    phylum={r["orf_id"]:r for r in read_rows(assignments/"sankey_orfs_phylum.tsv")}
    if len(phylum)!=len(genus) or set(phylum)!={r["orf_id"] for r in genus}:
        raise ValueError("Rank assignments disagree on ORF IDs")
    manifest=json.loads((run/"run.json").read_text())
    records=[]; targets={}
    db=sqlite3.connect((run/"species_hits.sqlite").resolve().as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    last=time.perf_counter()
    for r in genus:
        p=phylum[r["orf_id"]]
        row=dict(orf_id=r["orf_id"],contig=r["contig"],start=int(r["start"]),end=int(r["end"]),
                 gene=r["gene"],product=r["product"],
                 genus_status=r["assignment_status"],genus_label=r["taxon_label"],
                 other_genus=r["assignment_status"]=="assigned" and r["host_taxon_match"]=="False",
                 host_genus=r["assignment_status"]=="assigned" and r["host_taxon_match"]=="True",
                 phylum_status=p["assignment_status"],phylum_label=p["taxon_label"],
                 other_phylum=p["assignment_status"]=="assigned" and p["host_taxon_match"]=="False",
                 host_phylum=p["assignment_status"]=="assigned" and p["host_taxon_match"]=="True",
                 subject="",bitscore=None,qcov=None,scov=None,pident=None,evalue=None,quality=None)
        if r["best_bitscore"] not in ("", "NA"):
            expected=float(r["best_bitscore"])
            if float(p["best_bitscore"])!=expected:
                raise ValueError("Rank best scores disagree")
            hit=db.execute("SELECT subject,bits,qcov,scov FROM hits WHERE q=? AND bits=? ORDER BY qcov DESC,scov DESC,subject DESC LIMIT 1",(r["orf_id"],expected)).fetchone()
            if hit is None:
                raise ValueError("Stored best hit missing: "+r["orf_id"])
            row.update(subject=hit[0],bitscore=hit[1],qcov=hit[2],scov=hit[3])
            targets.setdefault(r["orf_id"].encode(),{})[hit[0].encode()]=row
        records.append(row)
        if time.perf_counter()-last>20:
            print(json.dumps(dict(stage="select_stored_winners",orfs=len(records),elapsed=time.perf_counter()-started)),flush=True)
            last=time.perf_counter()
    db.close()
    expected=sum(len(v) for v in targets.values())
    print(json.dumps(dict(stage="scan_saved_raw_hits",expected=expected,elapsed=time.perf_counter()-started)),flush=True)
    found=0; lines=0; last=time.perf_counter()
    with gzip.open(run/"hits.tsv.gz","rb") as f:
        for line in f:
            lines+=1
            if lines%1000000==0 and time.perf_counter()-last>=20:
                print(json.dumps(dict(stage="scan_saved_raw_hits",raw_lines=lines,recovered=found,expected=expected,elapsed=time.perf_counter()-started)),flush=True)
                last=time.perf_counter()
            # Avoid decoding 681 million rows; parse only a selected subject.
            q,rest=line.split(b"\t",1)
            subjects=targets.get(q)
            if subjects is None:
                continue
            sid,rest=rest.split(b"\t",1)
            row=subjects.get(sid)
            if row is None:
                continue
            parts=rest.split(b"\t")
            # parts now begin with pident (original column 3).
            bits=float(parts[9])
            qc=(abs(int(parts[5])-int(parts[4]))+1)/int(parts[10])
            sc=(abs(int(parts[7])-int(parts[6]))+1)/int(parts[11])
            if abs(bits-row["bitscore"])>1e-9 or abs(qc-row["qcov"])>1e-9 or abs(sc-row["scov"])>1e-9:
                continue
            identity=float(parts[0])
            row.update(pident=identity,evalue=float(parts[8]),quality=quality(identity,qc,sc))
            found+=1
            del subjects[sid]
            if not subjects:
                del targets[q]
            if found==expected:
                break
    with (out/"top_hit_quality.tsv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(records)
    groups=dict(all_orfs=records,other_genus=[r for r in records if r["other_genus"]],
                host_genus=[r for r in records if r["host_genus"]],
                other_phylum=[r for r in records if r["other_phylum"]],
                same_phylum_other_genus=[r for r in records if r["other_genus"] and r["host_phylum"]],
                nif_27=[r for r in records if r["contig"]=="contig_1" and r["start"]>=291513 and r["end"]<=313504])
    result=dict(complete=found==expected,expected_hits=expected,recovered_hits=found,raw_lines_scanned=lines,
                elapsed_seconds=time.perf_counter()-started,
                formula="pident/100 * qcov * scov",ranking="same species-best bit-score winners, same exclusions and coverage filters; deterministic qcov/scov/subject tie break",
                groups={k:summarize(v) for k,v in groups.items()},
                source_run=str(run.resolve()),source_assignments=str(assignments.resolve()),
                source_parameters=manifest["parameters"],
                input_hashes={name:hashlib.sha256((assignments/name).read_bytes()).hexdigest() for name in ("sankey_orfs_genus.tsv","sankey_orfs_phylum.tsv")},
                audit_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                search_rerun=False,hgt_scores_recomputed=False,
                interpretation="Quality audit of the original global top hits, not new HGT calls. A pass is alignment quality only; donor-bin quality and other eligibility rules are not recomputed.")
    (out/"summary.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(dict(stage="complete",complete=result["complete"],recovered=found,expected=expected,elapsed=result["elapsed_seconds"],groups={k:{f:v[f] for f in ("orfs","recovered","bins","retained_at")} for k,v in result["groups"].items()})),flush=True)
    if not result["complete"]:
        raise ValueError("Some selected hits were not recovered; audit is incomplete")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir",type=Path,required=True)
    p.add_argument("--assignments",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();audit(a.run_dir,a.assignments,a.out)
