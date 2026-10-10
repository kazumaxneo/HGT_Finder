"""Standalone local HTML: no CDN, remote assets, telemetry or uploads."""
import html
import json
from pathlib import Path
from urllib.parse import quote


def write_report(path: Path, scores, regions, manifest):
    fields = ["orf_id", "contig", "start", "end", "gene", "product", "orf_score", "window_score",
              "window_valid_orfs", "close_hit_species", "donor_species", "donor_label",
              "window_donor_agreement", "top_hit_pident", "top_hit_qcov", "top_hit_scov", "top_hit_quality", "quality_status", "flags"]
    headings = ["ORF", "contig", "開始", "終了", "遺伝子", "注釈", "ORFスコア", "移動平均",
                "有効ORF数", "近縁ヒット種数", "候補系統の種数", "類似先の候補系統",
                "window系統一致率", "最高ヒットの一致率 (%)", "query coverage", "subject coverage", "ヒット品質", "品質評価", "注意点"]
    body = []
    for r in scores:
        values = []
        for f in fields:
            v = r.get(f)
            values.append("NA" if v is None else f"{v:.3f}" if isinstance(v, float) else str(v))
        cells = [html.escape(v) for v in values]
        distributions = manifest.get("hit_distributions", {})
        if r["orf_id"] in distributions.get("orf_ids", []):
            href = distributions["index"] + "?orf=" + quote(r["orf_id"], safe="")
            cells[0] = '<a href="' + html.escape(href, quote=True) + '">' + cells[0] + '</a>'
        body.append('<tr data-id="' + html.escape(r["orf_id"], quote=True) + '">' +
                    ''.join('<td>' + v + '</td>' for v in cells) + '</tr>')
    compact = [{k: r.get(k) for k in ("orf_id", "contig", "start", "end", "gene", "product", "orf_score", "window_score", "flags", "donor_label")} for r in scores]
    # Escape HTML/script terminators in user-provided annotations.
    data = json.dumps(compact, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")
    region_rows = ''.join('<tr>' + ''.join('<td>' + html.escape(str(r.get(k, ""))) + '</td>' for k in
                          ("region_id", "contig", "start", "end", "n_orfs", "mean_window_score", "donor_taxon")) + '</tr>' for r in regions)
    options = ''.join('<option>' + html.escape(c) + '</option>' for c in sorted({r["contig"] for r in scores}))
    counts = manifest["counts"]
    label = html.escape(manifest.get("database_label", "Unspecified reference database"))
    page = '''<!doctype html><html lang="ja"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HGT_Finder local report</title>
<style>body{font:14px/1.6 system-ui;margin:0;color:#173140;background:#f7fafc}main{max-width:1450px;margin:auto;padding:24px}h1{margin-bottom:4px}header p{margin-top:0;color:#496573}.notice{padding:14px;background:#fff4db;border-left:4px solid #d39c23;margin:16px 0}.cards{display:flex;gap:18px;flex-wrap:wrap}.cards div{background:white;border:1px solid #d5e1e8;padding:12px 18px;border-radius:8px}.cards b{display:block;font-size:24px}.controls{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}input,select,button{padding:6px;font:inherit}canvas{display:block;width:100%;height:280px;background:white;border:1px solid #ccd9e1}table{border-collapse:collapse;width:100%;font-size:12px;background:white;color:black}td,th{border:1px solid black;padding:6px;text-align:left}th{position:sticky;top:0;background:white}.scroll{max-height:530px;overflow:auto}a{color:#14628f}.legend{color:#4b6271}#hover{min-height:50px;padding:8px;background:#eef4f7;white-space:pre-wrap}</style>
<main><header><h1>HGT_Finder · 0.1.2</h1><p>ORFの類似性分布と連続領域の探索</p></header>
<div class="notice">スコアはHGT確率ではありません。候補系統は類似配列を持つ参照系統で、歴史的ドナーを確定するものではありません。近縁ヒット未検出は、生物学的な遺伝子の不在を証明しません。</div>
<p><strong>参照DB：</strong>''' + label + '''</p>
<div class="cards"><div><b>''' + str(counts["orfs"]) + '''</b>ORF</div><div><b>''' + str(counts["window_eligible_orfs"]) + '''</b>windowに使用できるORF</div><div><b>''' + str(len(regions)) + '''</b>仮の閾値で選んだ候補領域</div></div>
<p>近縁の境界：''' + html.escape(manifest["parameters"]["close_rank"]) + ''' ／ 候補系統の集計単位：''' + html.escape(manifest["parameters"]["donor_rank"]) + ''' ／ window：''' + str(manifest["parameters"]["window"]) + ''' ORF ／ 仮の領域閾値：''' + str(manifest["parameters"]["threshold"]) + '''</p>
<div class="controls"><label>contig <select id="contig">''' + options + '''</select></label><label>開始bp <input id="from" type="number" min="1"></label><label>終了bp <input id="to" type="number" min="1"></label><button id="draw">範囲を表示</button><button id="reset">contig全体</button><label>ORF検索 <input id="filter" placeholder="nif / ORF ID"></label></div>
<canvas id="plot"></canvas><p class="legend">灰色の点：ORFスコア ／ 青い線：移動平均 ／ 橙の破線：候補領域の仮の閾値。NAは描画しません。</p><div id="hover">点の近くにカーソルを置くとORF情報を表示します。</div>
<h2>候補領域</h2><p>境界は移動平均で選ばれたORFの範囲です。転移の切断点を推定したものではありません。</p><div class="scroll"><table><tr><th>領域</th><th>contig</th><th>開始</th><th>終了</th><th>ORF数</th><th>平均windowスコア</th><th>候補系統</th></tr>''' + region_rows + '''</table></div>
<h2>全ORF</h2><div class="scroll"><table id="genes"><thead><tr>''' + ''.join('<th>' + h + '</th>' for h in headings) + '''</tr></thead><tbody>''' + ''.join(body) + '''</tbody></table></div>
<p><a href="orf_scores.tsv">ORFスコアTSV</a> · <a href="regions.tsv">候補領域TSV</a> · <a href="taxon_summaries.tsv">分類群別ヒット分布TSV</a> · <a href="run.json">パラメーター・入力・実行記録</a></p></main>
<script>
const data=''' + data + ''';
const threshold=''' + json.dumps(manifest["parameters"]["threshold"]) + ''';
const maxGap=''' + json.dumps(manifest['parameters']['max_gap_bp']) + ''';
const $=id=>document.getElementById(id);let displayed=[];let points=[];
function extent(){const a=data.filter(r=>r.contig===$('contig').value);$('from').value=a.reduce((v,r)=>Math.min(v,r.start),Infinity);$('to').value=a.reduce((v,r)=>Math.max(v,r.end),-Infinity);}
function draw(){const a=Number($('from').value),b=Number($('to').value);if(!Number.isFinite(a)||!Number.isFinite(b)||b<=a){$('hover').textContent='開始・終了座標を確認してください。';return;}displayed=data.filter(r=>r.contig===$('contig').value&&r.end>=a&&r.start<=b);const c=$('plot'),dpr=window.devicePixelRatio||1;c.width=c.clientWidth*dpr;c.height=280*dpr;const ctx=c.getContext('2d');ctx.scale(dpr,dpr);const w=c.clientWidth,h=280,left=55,right=w-15,top=20,bottom=h-35;const vals=displayed.flatMap(r=>[r.orf_score,r.window_score]).filter(v=>v!==null);const lo=vals.reduce((v,s)=>Math.min(v,s),-.2),hi=vals.reduce((v,s)=>Math.max(v,s),Math.max(.4,threshold));const x=v=>left+(v-a)/(b-a)*(right-left),y=v=>bottom-(v-lo)/(hi-lo)*(bottom-top);ctx.font='12px system-ui';ctx.fillStyle='#516879';for(let j=0;j<=4;j++){const v=lo+(hi-lo)*j/4;ctx.strokeStyle='#e3eaf0';ctx.beginPath();ctx.moveTo(left,y(v));ctx.lineTo(right,y(v));ctx.stroke();ctx.fillText(v.toFixed(2),4,y(v)+4);}ctx.fillText(a.toLocaleString(),left,bottom+23);ctx.fillText(b.toLocaleString()+' bp',Math.max(left,right-110),bottom+23);ctx.setLineDash([5,4]);ctx.strokeStyle='#c9901e';ctx.beginPath();ctx.moveTo(left,y(threshold));ctx.lineTo(right,y(threshold));ctx.stroke();ctx.setLineDash([]);points=[];for(const r of displayed){const xx=x((r.start+r.end)/2);if(r.orf_score!==null){ctx.fillStyle='#8298a8';ctx.beginPath();ctx.arc(xx,y(r.orf_score),2,0,7);ctx.fill();}points.push({r,x:xx});}ctx.strokeStyle='#1274aa';ctx.lineWidth=2;ctx.beginPath();let prev=null;for(const r of displayed){if(r.window_score===null){prev=null;continue;}const xx=x((r.start+r.end)/2),yy=y(r.window_score);if(prev&&r.start-prev.end-1<=maxGap)ctx.lineTo(xx,yy);else ctx.moveTo(xx,yy);prev=r;}ctx.stroke();}
$('draw').onclick=draw;$('reset').onclick=()=>{extent();draw();};$('contig').onchange=()=>{extent();draw();};window.addEventListener('resize',draw);
$('plot').onmousemove=e=>{const x=e.clientX-$('plot').getBoundingClientRect().left;const p=points.reduce((a,b)=>!a||Math.abs(b.x-x)<Math.abs(a.x-x)?b:a,null);if(p)$('hover').textContent=p.r.orf_id+' '+(p.r.gene||'')+' '+p.r.start+'–'+p.r.end+' bp\\n'+(p.r.product||'')+'\\nORF='+p.r.orf_score+' / window='+p.r.window_score+' / '+(p.r.donor_label||'')+' / '+(p.r.flags||'');};
$('filter').oninput=()=>{const q=$('filter').value.toLowerCase();for(const tr of $('genes').tBodies[0].rows)tr.hidden=!tr.textContent.toLowerCase().includes(q);};extent();draw();
</script></html>'''
    if manifest.get("sankey"):
        page = page.replace("<h2>候補領域</h2>",
                            '<p><a href="sankey.html">Contig → 分類群のサンキープロットを開く</a></p><h2>候補領域</h2>')
    if manifest.get("hit_distributions"):
        page = page.replace("<h2>候補領域</h2>", '<p><a href="' +
                            html.escape(manifest["hit_distributions"]["index"], quote=True) +
                            '">ORFごとの分類群別ヒット分布を開く（ORF IDからも選択できます）</a></p><h2>候補領域</h2>')
    path.write_text(page)
