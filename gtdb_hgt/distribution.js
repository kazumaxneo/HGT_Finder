/* Species-best histograms and exact summaries. No sampling or network requests. */
const distributionRanks={phylum:1,class:2,order:3,family:4,genus:5};
const distributionRankNames={phylum:'門',class:'綱',order:'目',family:'科',genus:'属'};
function distributionEscape(value){return String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function distributionQuantile(sorted,p){if(!sorted.length)return null;const pos=(sorted.length-1)*p,lo=Math.floor(pos),hi=Math.ceil(pos);return sorted[lo]+(sorted[hi]-sorted[lo])*(pos-lo);}
function distributionGroups(points,host,rank,floor=.5,filter=false,maxTaxa=12){
 if(!(rank in distributionRanks)||!Number.isFinite(floor)||floor<0||floor>1||!Number.isInteger(maxTaxa)||maxTaxa<1)throw Error('Invalid display settings');
 const depth=distributionRanks[rank]+1; // Host keys use full lineage, never only the leaf name.
 const allRanks=['domain','phylum','class','order','family','genus','species'];
 const hostTax=allRanks.slice(0,depth).map(r=>host[r]),hostKey=JSON.stringify(hostTax),groups=new Map();
 const accepted=points.filter(p=>!filter||(p.quality!==null&&p.quality>=floor));
 for(const p of accepted){const tax=p.taxonomy.slice(0,depth),key=JSON.stringify(tax);if(!groups.has(key))groups.set(key,{key,tax,label:tax.at(-1)||'分類未確定',host:key===hostKey,points:[]});groups.get(key).points.push(p);}
 if(!groups.has(hostKey))groups.set(hostKey,{key:hostKey,tax:hostTax,label:hostTax.at(-1),host:true,points:[]});
 for(const g of groups.values()){
  const values=g.points.map(p=>p.bitscore).sort((a,b)=>a-b);
  Object.assign(g,{n:values.length,min:values.length?values[0]:null,max:values.length?values.at(-1):null,q1:distributionQuantile(values,.25),median:distributionQuantile(values,.5),q3:distributionQuantile(values,.75),low:g.points.filter(p=>p.quality!==null&&p.quality<floor).length,unknown:g.points.filter(p=>p.quality===null).length});
 }
 const all=Array.from(groups.values()).sort((a,b)=>Number(b.host)-Number(a.host)||(b.median??-1)-(a.median??-1)||a.key.localeCompare(b.key));
 const shown=all.slice(0,maxTaxa),labels=new Map();for(const g of all)labels.set(g.label,(labels.get(g.label)||0)+1);
 for(const g of shown)g.displayLabel=labels.get(g.label)>1?g.label+' ['+(g.tax.at(-2)||'分類未確定')+']':g.label;
 return {groups:shown,totalGroups:all.length,total:points.length,accepted:accepted.length,shown:shown.reduce((n,g)=>n+g.n,0),hiddenGroups:all.length-shown.length,maxBits:points.reduce((n,p)=>Math.max(n,p.bitscore),1),floor,filtered:filter};
}
function distributionHistogram(model,binCount=30){
 if(!Number.isInteger(binCount)||binCount<1||binCount>200)throw Error('Bin count must be an integer from 1 to 200');
 const magnitude=10**Math.floor(Math.log10(model.maxBits));
 const upper=Math.ceil(model.maxBits/magnitude)*magnitude,binWidth=upper/binCount;
 const groups=model.groups.map(g=>{
  const counts=Array(binCount).fill(0);
  for(const p of g.points)counts[Math.min(binCount-1,Math.floor(p.bitscore/binWidth))]++;
  return {...g,counts};
 });
 return {...model,groups,binCount,binWidth,upper};
}
function distributionSvg(model){
 const hist=distributionHistogram(model,model.binCount??30),esc=distributionEscape;
 const cols=Math.min(3,hist.groups.length),panelW=400,panelH=270,width=cols*panelW,height=60+Math.ceil(hist.groups.length/cols)*panelH;
 const fmt=v=>v==null?'—':Number(v.toFixed(2)).toString();
 let svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="分類群別ヒストグラム：横軸bit score、縦軸GTDB参照種数"><rect width="100%" height="100%" fill="white"/><style>text{font-family:system-ui,sans-serif}</style><text x="20" y="24" font-size="16">Taxon hit histograms · ${hist.binCount} common bins</text><text x="20" y="46" font-size="12">x: bit score / y: reference species count (scale varies by panel)</text>`;
 hist.groups.forEach((g,i)=>{
  const ox=(i%cols)*panelW,oy=60+Math.floor(i/cols)*panelH,left=ox+58,right=ox+380,top=oy+50,bottom=oy+212;
  const maximum=Math.max(1,...g.counts),y=v=>bottom-v/maximum*(bottom-top),x=v=>left+v/hist.upper*(right-left),color=g.host?'#148c82':'#dd8730';
  svg+=`<text x="${ox+20}" y="${oy+20}" font-size="14" fill="${g.host?'#0b7068':'black'}">${esc(g.displayLabel.length>35?g.displayLabel.slice(0,34)+'…':g.displayLabel)}${g.host?' (host)':''}<title>${esc(g.tax.join(';'))}</title></text><text x="${ox+20}" y="${oy+38}" font-size="12">n=${g.n.toLocaleString()} / median=${fmt(g.median)}${g.n?'':' (no hits)'}</text>`;
  const ticks=[...new Set([0,Math.round(maximum/2),maximum])];
  for(const value of ticks)svg+=`<path d="M ${left} ${y(value)} H ${right}" stroke="#ddd"/><text x="${left-7}" y="${y(value)+4}" text-anchor="end" font-size="11">${value}</text>`;
  g.counts.forEach((count,j)=>{
   const low=j*hist.binWidth,high=(j+1)*hist.binWidth;
   svg+=`<rect data-bin="true" data-count="${count}" x="${x(low)}" y="${y(count)}" width="${(right-left)/hist.binCount}" height="${bottom-y(count)}" fill="${color}" stroke="white" stroke-width=".4"><title>${esc(g.displayLabel)}; [${fmt(low)}, ${fmt(high)}${j===hist.binCount-1?']':')'} bit; ${count} reference species</title></rect>`;
  });
  svg+=`<path d="M ${left} ${top} V ${bottom} H ${right}" stroke="black" fill="none"/>`;
  for(let t=0;t<=4;t++){const value=hist.upper*t/4;svg+=`<text x="${x(value)}" y="${bottom+17}" text-anchor="middle" font-size="11">${fmt(value)}</text>`;}
  svg+=`<text x="${(left+right)/2}" y="${bottom+38}" text-anchor="middle" font-size="12">bit score</text>`;
 });
 return svg+'</svg>';
}
function distributionBinsTsv(model){
 const hist=distributionHistogram(model,model.binCount??30);
 const rows=['taxon\tbin_start_bitscore\tbin_end_bitscore\tend_inclusive\tspecies_count'];
 for(const g of hist.groups)g.counts.forEach((n,i)=>rows.push([g.tax.join(';'),i*hist.binWidth,(i+1)*hist.binWidth,i===hist.binCount-1,n].join('\t')));
 return rows.join('\n')+'\n';
}
if(typeof module!=='undefined')module.exports={distributionGroups,distributionSvg,distributionHistogram,distributionBinsTsv,distributionQuantile};
if(typeof document!=='undefined'){
 const el=id=>document.getElementById(id);let activeRank=payload.config.rank,model=null;
 el('annotation').textContent=[payload.gene.gene,payload.gene.product,payload.gene.contig+':'+payload.gene.start+'–'+payload.gene.end].filter(Boolean).join(' · ');
 function updateDistribution(){
  const maximum=Number(el('max-taxa').value),floor=Number(el('quality-floor').value),bins=Number(el('bins').value);
  if(!el('max-taxa').value||!el('quality-floor').value||!Number.isInteger(maximum)||maximum<1||!Number.isFinite(floor)||floor<0||floor>1||!Number.isInteger(bins)||bins<1||bins>200){el('stats').textContent='分類群数は1以上の整数、品質下限は0–1、ビン数は1–200の整数を指定してください。';return;}
  model=distributionGroups(payload.points,payload.host,activeRank,floor,el('quality-filter').checked,maximum);
  model.binCount=bins;
  el('diagram').innerHTML=distributionSvg(model);
  for(const button of el('rank-controls').querySelectorAll('[data-rank]'))button.setAttribute('aria-pressed',String(button.dataset.rank===activeRank));
  el('stats').textContent='表示階層：'+distributionRankNames[activeRank]+' ／ 保存済み '+model.total.toLocaleString()+'参照種 ／ フィルター後 '+model.accepted.toLocaleString()+'種 ／ 図中 '+model.shown.toLocaleString()+'種 ／ '+model.groups.length+'/'+model.totalGroups+'分類群（非表示 '+model.hiddenGroups+'群）。宿主分類群を表示に含めます。';
  const unknown=payload.points.filter(p=>p.quality===null).length;
  el('quality-note').textContent='品質下限 '+floor+'。'+(model.filtered?'品質取得済みで下限以上のヒットから度数・中央値・支持数を再計算しています。':'元の比較に使った全有効ヒットから中央値・支持数を計算しています。低品質のヒットも保持します。')+(unknown?' '+unknown.toLocaleString()+'種は品質未取得です。古いDBの一致率を推定して補っていません。':'');
  const fmt=v=>v==null?'—':Number(v.toFixed(3)).toString();
  el('summary').innerHTML=model.groups.map(g=>'<tr>'+[g.displayLabel,g.host?'はい':'いいえ',g.n,fmt(g.min),fmt(g.q1),fmt(g.median),fmt(g.q3),fmt(g.max),g.low,g.unknown].map(v=>'<td>'+distributionEscape(v)+'</td>').join('')+'</tr>').join('');
 }
 for(const button of el('rank-controls').querySelectorAll('[data-rank]'))button.onclick=()=>{activeRank=button.dataset.rank;updateDistribution();};
 el('apply').onclick=updateDistribution;el('quality-filter').onchange=updateDistribution;
 function download(content,type,suffix){if(!model)return;const url=URL.createObjectURL(new Blob([content],{type})),a=document.createElement('a');a.href=url;a.download=payload.gene.orf_id+'_'+activeRank+suffix;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
 el('export').onclick=()=>{if(model)download(distributionSvg(model),'image/svg+xml','_hit_distribution.svg');};
 el('export-bins').onclick=()=>{if(model)download(distributionBinsTsv(model),'text/tab-separated-values','_histogram_bins.tsv');};
 updateDistribution();
}
