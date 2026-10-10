"""Read-only, species-best hit distributions. No search or HGT rescoring."""
import csv
import html
import json
import math
import sqlite3
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

from .io import RANKS, canonical_accession, digest, dump_json, read_subject_map, read_taxonomy, tax_key, tsv_write
from .quality import alignment_quality
from .sankey import DISPLAY_RANKS, RANK_LABELS


def quantile(values, fraction):
    values = sorted(values)
    position = (len(values) - 1) * fraction
    lo, hi = math.floor(position), math.ceil(position)
    return values[lo] + (values[hi] - values[lo]) * (position - lo)


def summarize_points(points, host, rank, floor=.5):
    groups = defaultdict(list)
    for p in points:
        groups[tuple(p["taxonomy"][:RANKS.index(rank) + 1])].append(p)
    host_key = tax_key(host, rank)
    result = []
    for key, values in groups.items():
        bits = [v["bitscore"] for v in values]
        result.append(dict(rank=rank, taxon=";".join(key), taxon_label=key[-1] or "分類未確定",
                           host_match=key == host_key, species_count=len(values), min_bitscore=min(bits),
                           q1_bitscore=quantile(bits, .25), median_bitscore=median(bits),
                           q3_bitscore=quantile(bits, .75), max_bitscore=max(bits),
                           quality_known=sum(v["quality"] is not None for v in values),
                           below_quality_floor=sum(v["quality"] is not None and v["quality"] < floor for v in values)))
    return sorted(result, key=lambda r: (not r["host_match"], -r["median_bitscore"], r["taxon"]))


class SavedHitReader:
    def __init__(self, run_dir, manifest, taxonomy=None, subject_map=None):
        path = Path(run_dir) / "species_hits.sqlite"
        if not path.is_file():
            raise ValueError("species_hits.sqlite is required: medians cannot reconstruct individual hit distributions")
        self.db = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
        self.db.execute("PRAGMA query_only=ON")
        self.columns = {r[1] for r in self.db.execute("PRAGMA table_info(hits)")}
        if not {"q", "species", "genome", "subject", "bits", "qcov", "scov"}.issubset(self.columns):
            self.close()
            raise ValueError("Unrecognized species hit store schema")
        self.embedded = self.db.execute("SELECT 1 FROM sqlite_master WHERE name='reference_taxonomy'").fetchone() is not None
        self.taxa = None
        self.subject_map = {}
        self.taxonomy_hash = None
        try:
            if not self.embedded:
                inputs = manifest.get("inputs", {})
                taxonomy = Path(taxonomy or inputs.get("taxonomy", ""))
                if not taxonomy.is_file():
                    raise ValueError("Legacy hit store needs the original taxonomy TSV: use --taxonomy for a relocated copy")
                self.taxonomy_hash = digest(taxonomy)
                if inputs.get("taxonomy_sha256") and self.taxonomy_hash != inputs["taxonomy_sha256"]:
                    raise ValueError("Taxonomy hash differs from the original analysis")
                self.taxa = read_taxonomy(taxonomy)
                self.subject_map = read_subject_map(subject_map or inputs.get("subject_map"))
        except Exception:
            self.close()
            raise

    def get(self, orf_id):
        optional = ["h." + k if k in self.columns else "NULL" for k in ("pident", "evalue")]
        select = "SELECT h.species,h.genome,h.subject,h.bits,h.qcov,h.scov," + ",".join(optional)
        if self.embedded:
            select += ",r.accession,r.taxonomy FROM hits h LEFT JOIN reference_taxonomy r ON h.genome=r.genome"
        else:
            select += " FROM hits h"
        records = self.db.execute(select + " WHERE h.q=? ORDER BY h.species", (orf_id,))
        points, species_seen = [], set()
        for record in records:
            _, _, subject, bits, qc, sc, identity, evalue = record[:8]
            if self.embedded:
                genome, encoded = record[8:]
                if encoded is None:
                    raise ValueError("Hit store has an unmapped genome index")
                tax = json.loads(encoded)
            else:
                genome = self.subject_map.get(subject)
                if genome is None and "::" in subject:
                    genome = canonical_accession(subject.split("::", 1)[0])
                if genome not in self.taxa:
                    raise ValueError("Cannot map saved subject to taxonomy; supply --subject-map for non-encoded IDs")
                tax = self.taxa[genome]
            species = tax_key(tax, "species")
            if species is None or species in species_seen:
                raise ValueError("Saved hits must contain exactly one best hit per complete GTDB species")
            species_seen.add(species)
            if not math.isfinite(bits) or bits <= 0 or not all(math.isfinite(v) and 0 <= v <= 1 for v in (qc, sc)):
                raise ValueError("Invalid bit score or coverage in saved hit store")
            quality = alignment_quality(identity, qc, sc) if identity is not None else None
            if evalue is not None and (not math.isfinite(evalue) or evalue < 0):
                raise ValueError("Invalid E-value in saved hit store")
            points.append(dict(genome_id=genome, subject=subject, taxonomy=[tax[r] for r in RANKS],
                               bitscore=bits, qcov=qc, scov=sc, pident=identity, evalue=evalue, quality=quality))
        return points

    def close(self):
        self.db.close()


def write_distribution_page(folder, name, gene, points, host, config):
    payload = dict(gene=gene, points=points, host=host, config=config)
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")
    buttons = ''.join(f'<button type="button" data-rank="{r}" aria-pressed="false">{RANK_LABELS[r]}</button>' for r in DISPLAY_RANKS)
    label = html.escape(gene["orf_id"])
    page = '''<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HGT_Finder · Hit distribution</title><style>
body{font:15px/1.6 system-ui;margin:0;background:#f6f8fa;color:#173140}main{max-width:1450px;margin:auto;padding:24px}button,input,select{font:inherit;padding:5px}button[aria-pressed=true]{font-weight:bold;border:2px solid #173140}.controls{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}input[type=number]{width:85px}.scroll{overflow:auto;background:white;border:1px solid #ccd7dd}svg{display:block;width:100%;height:auto;min-width:700px}.taxa{max-height:240px;overflow:auto;border:1px solid #ccd7dd;background:white;padding:10px;display:grid;grid-template-columns:repeat(auto-fill,minmax(260px,1fr));gap:6px}.taxa label[hidden]{display:none}.taxa input{margin-right:4px}table{border-collapse:collapse;background:white;color:black;width:100%;font-size:13px}th,td{border:1px solid black;padding:6px;text-align:left}.muted{color:#536c7b}a{color:#146b98}
</style></head><body><main><h1>分類群ごとのヒット分布</h1><p>HGT_Finder · ''' + label + '''</p>
<p id="annotation"></p><p>各GTDB参照種の最良ヒットをbit scoreの区間ごとに数えたヒストグラムです。有意性・coverage・自己除外条件を通過したヒットを使います。HGT確率や転移方向を示す図ではありません。</p>
<div id="rank-controls" class="controls">''' + buttons + '''</div>
<div class="controls"><label>自動選択する分類群数 <input id="max-taxa" type="number" min="1" value="''' + str(config["max_taxa"]) + '''"></label>
<label>表示形式 <select id="plot-mode"><option value="overlay">選んだ分類群を重ねる</option><option value="panels">分類群ごとに並べる</option></select></label>
<label>ビン数 <input id="bins" type="number" min="1" max="200" value="30"></label>
<label><input id="quality-filter" type="checkbox">品質下限以上のみ（未取得も除外）</label>
<label>品質下限 <input id="quality-floor" type="number" min="0" max="1" step="0.05" value="''' + str(config["min_hit_quality"]) + '''"></label>
<button id="apply">表示を更新</button><button id="export">SVG保存</button><button id="export-bins">表示中の度数TSV保存</button></div>
<h2>比較する分類群を選択</h2><p>階級を選び、比較する属・科などにチェックしてください。全分類群から選べます。例えば「宿主のみ」を押してから、比較したい別属を検索して追加できます。</p>
<div class="controls"><label>分類群を検索 <input id="taxon-search" type="search"></label><button id="select-host">宿主のみ</button><button id="clear-taxa">選択を解除</button><button id="reset-taxa">上位Nを選択</button></div><div id="taxon-choices" class="taxa"></div>
<p id="stats"></p><p id="quality-note" class="muted"></p><div id="diagram" class="scroll"></div>
<p class="muted">横軸はbit score、縦軸はGTDB参照種数です。選択した分類群を色分けして同じ座標に重ねます。分類群ごとのパネル表示にも切り替えられます。全分類群・全階級で同じ区間幅と横軸範囲を使います。重ね描きの縦軸は共通、パネル表示の縦軸の範囲はパネルごとに異なります。区間は下端を含み上端を含まず、最終区間のみ上端も含みます。選択リストは宿主群を先頭に固定し、その他は中央値の高い順です。青緑は宿主群、他群は凡例の色で示します。選択内容は階級ごとに保持します。検索は選択リストの表示だけに作用します。棒にカーソルを置くと区間と種数が見えます。絞り込みはこの図だけに作用し、元のHGTスコア・候補領域を変更しません。</p>
<h2>表示中の分類群</h2><table><thead><tr><th>分類群</th><th>宿主と同じ</th><th>支持種数</th><th>最小</th><th>25%</th><th>中央値</th><th>75%</th><th>最大</th><th>低品質</th><th>品質未取得</th></tr></thead><tbody id="summary"></tbody></table>
<p><a href="index.html">ORFの選択に戻る</a> · <a href="''' + name + '''_hits.tsv">全参照種のヒットTSV</a> · <a href="''' + name + '''_summary.tsv">全5階層の集計TSV（品質で絞る前）</a></p>
</main><script>const payload=''' + encoded + ';\n' + Path(__file__).with_name("distribution.js").read_text() + '\n</script></body></html>'
    (folder / (name + ".html")).write_text(page)
    rows = [{"orf_id": gene["orf_id"], **{k: v for k, v in p.items() if k != "taxonomy"},
             "gtdb_taxonomy": ";".join(p["taxonomy"])} for p in points]
    tsv_write(folder / (name + "_hits.tsv"), rows,
              ["orf_id", "genome_id", "subject", "bitscore", "qcov", "scov", "pident", "evalue", "quality", "gtdb_taxonomy"])
    summaries = [r for rank in DISPLAY_RANKS for r in summarize_points(points, host, rank, config["min_hit_quality"])]
    tsv_write(folder / (name + "_summary.tsv"), summaries,
              ["rank", "taxon", "taxon_label", "host_match", "species_count", "min_bitscore", "q1_bitscore",
               "median_bitscore", "q3_bitscore", "max_bitscore", "quality_known", "below_quality_floor"])


def plot_from_run(source, folder, orf_ids=None, rank="genus", max_taxa=12, min_hit_quality=.5,
                  taxonomy=None, subject_map=None):
    source, folder = Path(source), Path(folder)
    if rank not in DISPLAY_RANKS or max_taxa < 1 or not math.isfinite(min_hit_quality) or not 0 <= min_hit_quality <= 1:
        raise ValueError("Invalid distribution rank, display limit or quality floor")
    if folder.exists() and any(folder.iterdir()):
        raise ValueError("Distribution output must be empty")
    manifest = json.loads((source / "run.json").read_text())
    with (source / "orf_scores.tsv").open() as f:
        genes = list(csv.DictReader(f, delimiter="\t"))
    by_id = {g["orf_id"]: g for g in genes}
    if len(by_id) != len(genes):
        raise ValueError("Duplicate ORF IDs in score table")
    if orf_ids:
        unknown = set(orf_ids) - set(by_id)
        if unknown:
            raise ValueError(f"Unknown ORF IDs: {sorted(unknown)[:5]}")
        genes = [by_id[q] for q in dict.fromkeys(orf_ids)]
    if any(tax_key(manifest["query_taxonomy"], r) is None for r in DISPLAY_RANKS):
        raise ValueError("Host taxonomy must be complete through genus for rank switching")
    reader = SavedHitReader(source, manifest, taxonomy, subject_map)
    entries = []
    config = dict(rank=rank, max_taxa=max_taxa, min_hit_quality=min_hit_quality, plot_type="histogram", bins=30, plot_mode="overlay")
    folder.mkdir(parents=True, exist_ok=True)
    try:
        for i, gene in enumerate(genes, 1):
            points = reader.get(gene["orf_id"])
            name = f"orf_{i:06d}"
            write_distribution_page(folder, name, gene, points, manifest["query_taxonomy"], config)
            entries.append(dict(orf_id=gene["orf_id"], gene=gene.get("gene", ""), product=gene.get("product", ""),
                                species_hits=len(points), quality_known=sum(p["quality"] is not None for p in points),
                                page=name + ".html"))
            if i % 100 == 0:
                print(f"Hit distributions: {i}/{len(genes)} ORFs", flush=True)
    finally:
        reader.close()
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), config=config, orfs=len(entries),
                   species_hits=sum(e["species_hits"] for e in entries), entries=entries,
                   source_run_sha256=digest(source / "run.json"), source_scores_sha256=digest(source / "orf_scores.tsv"),
                   source_taxonomy_sha256=reader.taxonomy_hash, search_rerun=False, hgt_scores_recomputed=False,
                   definition="One species-best filtered hit per GTDB species; all comparison hits retained; display groups limited explicitly; quality filtering is optional and display-only")
    dump_json(folder / "distribution.json", summary)
    table = ''.join('<tr>' + ''.join('<td>' + html.escape(str(e[k])) + '</td>' for k in ("orf_id", "gene", "product", "species_hits", "quality_known")) +
                    '<td><a href="' + e["page"] + '">分布を見る</a></td></tr>' for e in entries)
    routes = json.dumps({e["orf_id"]: e["page"] for e in entries}, ensure_ascii=False).replace("<", "\\u003c").replace("&", "\\u0026")
    (folder / "index.html").write_text('''<!doctype html><html lang="ja"><meta charset="utf-8"><title>HGT_Finder · Hit distributions</title>
<style>body{font:15px/1.6 system-ui;padding:24px}table{border-collapse:collapse;color:black;background:white}td,th{border:1px solid black;padding:8px}input{font:inherit}</style>
<h1>ORFごとのヒット分布</h1><p>各GTDB参照種の最良ヒットを数えたヒストグラム。各図で門・綱・目・科・属を切り替えられます。</p>
<label>ORF・遺伝子を検索 <input id="filter"></label><p id="message"></p><table><thead><tr><th>ORF</th><th>遺伝子</th><th>注釈</th><th>支持種数</th><th>品質取得済み</th><th>図</th></tr></thead><tbody>''' + table + '''</tbody></table>
<p><a href="distribution.json">生成条件・出典</a></p><script>const routes=''' + routes + ''';const q=new URLSearchParams(location.search).get('orf');
if(q&&Object.hasOwnProperty.call(routes,q))location.replace(routes[q]);else if(q)document.getElementById('message').textContent='このORFの分布は生成されていません：'+q;
document.getElementById('filter').oninput=e=>{for(const r of document.querySelectorAll('tbody tr'))r.hidden=!r.textContent.toLowerCase().includes(e.target.value.toLowerCase());};</script></html>''')
    return summary
