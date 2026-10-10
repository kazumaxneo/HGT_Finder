"""Read selected strongest host/other-phylum hits from the saved raw search."""
import argparse,csv,gzip,json,sqlite3,time
from pathlib import Path


def audit(run, taxonomy, out, queries, max_seconds):
    started=time.perf_counter()
    manifest=json.loads((run/'run.json').read_text())
    host=manifest['query_taxonomy']['phylum']
    taxa={}
    with taxonomy.open() as f:
        for line in f:
            accession,tax=line.strip().split('\t');taxa[accession]=tax.split(';')
    db=sqlite3.connect('file:'+str(run/'species_hits.sqlite')+'?mode=ro',uri=True)
    targets={};examples=[]
    for q in queries:
        candidates={}
        for sid,bits,qc,sc in db.execute('SELECT subject,bits,qcov,scov FROM hits WHERE q=? ORDER BY bits DESC',(q,)):
            tax=taxa[sid.split('::')[0]]
            category='host_phylum' if tax[1]==host else 'other_phylum'
            if category not in candidates:
                candidates[category]={'subject':sid,'bitscore':bits,'qcov':qc,'scov':sc,'taxonomy':';'.join(tax)}
            if len(candidates)==2:break
        examples.append({'orf_id':q,'hits':candidates})
        for value in candidates.values():targets[q,value['subject']]=value
    db.close()
    relevant={q.encode() for q in queries};found=set();lines=0;stopped='all_selected_hits_found'
    with gzip.open(run/'hits.tsv.gz','rb') as f:
        for line in f:
            lines+=1
            if lines%100000==0 and time.perf_counter()-started>max_seconds:
                stopped='time_limit';break
            if line.split(b'\t',1)[0] not in relevant:continue
            parts=line.decode().rstrip('\n').split('\t')
            key=tuple(parts[:2]);target=targets.get(key)
            if target is None or key in found:continue
            bits=float(parts[11]);qc=(abs(int(parts[7])-int(parts[6]))+1)/int(parts[12]);sc=(abs(int(parts[9])-int(parts[8]))+1)/int(parts[13])
            if abs(bits-target['bitscore'])>1e-9 or abs(qc-target['qcov'])>1e-9 or abs(sc-target['scov'])>1e-9:continue
            target.update(percent_identity=float(parts[2]),alignment_length=int(parts[3]),evalue=float(parts[10]),query_length=int(parts[12]),subject_length=int(parts[13]))
            found.add(key)
            if len(found)==len(targets):break
        else:stopped='end_of_file'
    result={'selected_queries':queries,'examples':examples,'expected_hits':len(targets),'retrieved_hits':len(found),'raw_lines_scanned':lines,'elapsed_seconds':time.perf_counter()-started,'stopped':stopped,'note':'Purposive examples, not a representative estimate of all other-lineage hits.'}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True);p.add_argument('--taxonomy',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--queries',nargs='+',required=True);p.add_argument('--max-seconds',type=float,default=180)
    a=p.parse_args();audit(a.run_dir,a.taxonomy,a.out,a.queries,a.max_seconds)
