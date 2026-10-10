"""Audit existing genus Sankey decisions without rerunning DIAMOND or HGT scoring."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path


def audit(folder):
    run = json.loads((folder / 'run.json').read_text())
    config = run['sankey']
    if config['rank'] != 'genus' or run['parameters']['donor_rank'] != 'genus' or run['parameters']['close_rank'] != 'genus':
        raise ValueError('This audit requires genus summaries and genus Sankey assignments')
    assignments = {r['orf_id']: r for r in csv.DictReader((folder / 'sankey_orfs.tsv').open(), delimiter='\t')}
    host = ';'.join(run['query_taxonomy'][rank] for rank in ['domain','phylum','class','order','family','genus'])
    minimum, margin = config['min_species'], config['margin']
    counts, examples, per_orf = Counter(), [], []

    def finish(q, bins):
        if not bins:
            return
        row = assignments[q]
        supported = sorted((b for b in bins if b['species'] >= minimum), key=lambda b: (-b['median'], -b['species'], b['taxon']))
        overall_max = max(b['max'] for b in bins)
        max_genera = sorted(b['taxon'] for b in bins if b['max'] == overall_max)
        host_bin = next((b for b in bins if b['taxon'] == host), None)
        contenders = [b for b in supported if b['median'] >= supported[0]['median'] * (1-margin)] if supported else []
        host_compatible = any(b['taxon'] == host for b in contenders)
        ranks = [b['taxon'].split(';') for b in contenders]
        shared = []
        for column in zip(*ranks):
            if len(set(column)) != 1:
                break
            shared.append(column[0])
        item = {'orf_id': q, 'original_status': row['assignment_status'], 'median_contending_genera': len(contenders),
                'host_in_median_competitors': host_compatible, 'median_competitors_lca': ';'.join(shared),
                'host_median_bitscore': host_bin['median'] if host_bin else None,
                'host_max_bitscore': host_bin['max'] if host_bin else None,
                'overall_max_bitscore': overall_max, 'global_best_genera': max_genera,
                'global_best_unique_host': max_genera == [host],
                'host_max_within_5pct_of_global_best': bool(host_bin and host_bin['max'] >= overall_max*(1-margin))}
        per_orf.append(item)
        counts['orfs_with_summary'] += 1
        if row['assignment_status'] != 'ambiguous':
            return
        counts['ambiguous'] += 1
        counts['ambiguous_host_in_all_median_competitors'] += host_compatible
        counts['ambiguous_unique_global_best_host_genus'] += item['global_best_unique_host']
        counts['ambiguous_global_best_includes_host_genus'] += host in max_genera
        counts['ambiguous_host_max_within_5pct_of_global_best'] += item['host_max_within_5pct_of_global_best']
        counts['ambiguous_all_median_competitors_same_host_phylum'] += len(shared)>1 and shared[1] == run['query_taxonomy']['phylum']
        if item['global_best_unique_host'] and row['leading_taxon'] != host and len(examples)<8:
            examples.append({**item, 'median_winner':row['leading_taxon'], 'median_winner_bitscore':float(row['median_bitscore'])})

    q, bins = None, []
    with (folder / 'taxon_summaries.tsv').open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        for r in reader:
            if r['orf_id'] != q:
                finish(q,bins)
                q,bins=r['orf_id'],[]
            bins.append({'taxon':r['taxon'],'species':int(r['species_count']),
                         'median':float(r['median_bitscore']),'max':float(r['max_bitscore'])})
        finish(q,bins)
    return {'counts':dict(counts),'examples':examples,'per_orf':per_orf,
            'note':'Maximum-hit taxonomy is a diagnostic, not a validated replacement assignment. All median competitors, not only the top two, were checked.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    result=audit(a.folder)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'counts':result['counts'],'examples':result['examples'][:3]},ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
