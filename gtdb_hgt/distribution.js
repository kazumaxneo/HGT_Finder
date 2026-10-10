/* Species-best points and exact summaries. No sampling or network requests. */
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
function distributionJitter(p){let hash=2166136261;for(const ch of p.genome_id+'|'+p.taxonomy[6]){hash^=ch.charCodeAt(0);hash=Math.imul(hash,16777619);}return ((hash>>>0)/4294967295-.5)*82;}
function distributionSvg(model){
 const esc=distributionEscape,column=150,left=72,top=40,bottom=465,width=left+column*model.groups.length+40,height=690,max=model.maxBits*1.05;
 const y=v=>bottom-v/max*(bottom-top),fmt=v=>v==null?'—':Number(v.toFixed(2)).toString();
 let svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-label="分類群ごとの参照種別bit score分布"><rect width="100%" height="100%" fill="white"/><style>text{font-family:system-ui,sans-serif}</style>`;
 for(let i=0;i<=5;i++){const value=max*i/5;svg+=`<path d="M ${left-8} ${y(value)} H ${width-15}" stroke="#ddd"/><text x="${left-14}" y="${y(value)+4}" text-anchor="end" font-size="12">${fmt(value)}</text>`;}
 svg+=`<text x="16" y="24" font-size="14">bit score</text>`;
 model.groups.forEach((g,i)=>{
  const x=left+column*i+column/2,color=g.host?'#148c82':'#dd8730';
  if(g.n){svg+=`<path d="M ${x} ${y(g.min)} V ${y(g.max)} M ${x-12} ${y(g.min)} H ${x+12} M ${x-12} ${y(g.max)} H ${x+12}" stroke="black" fill="none"/><rect x="${x-26}" y="${y(g.q3)}" width="52" height="${Math.max(0,y(g.q1)-y(g.q3))}" fill="none" stroke="black"/>`;}
  for(const p of g.points){
   const weak=p.quality===null||p.quality<model.floor;
   const title=[p.genome_id,p.taxonomy[6],p.subject,'bit score: '+p.bitscore,'query coverage: '+p.qcov,'subject coverage: '+p.scov,'identity (%): '+(p.pident??'未取得'),'E-value: '+(p.evalue??'未取得'),'quality: '+(p.quality??'未取得')].join('\n');
   svg+=`<circle data-point="true" cx="${x+distributionJitter(p)}" cy="${y(p.bitscore)}" r="2.4" fill="${weak?'white':color}" fill-opacity="${weak?'.8':'.45'}" stroke="${weak?'#737d85':color}" stroke-width=".6"><title>${esc(title)}</title></circle>`;
  }
  if(g.n)svg+=`<path d="M ${x-32} ${y(g.median)} H ${x+32}" stroke="black" stroke-width="3"><title>中央値: ${g.median}; n=${g.n}</title></path>`;
  svg+=`<text x="${x}" y="${bottom+21}" text-anchor="middle" font-size="12">n=${g.n}</text><text transform="translate(${x-35} ${bottom+48}) rotate(40)" font-size="13" fill="${g.host?'#0b7068':'black'}">${esc(g.displayLabel.length>31?g.displayLabel.slice(0,30)+'…':g.displayLabel)}<title>${esc(g.tax.join(';'))}</title></text>`;
 });
 return svg+'</svg>';
}
if(typeof module!=='undefined')module.exports={distributionGroups,distributionSvg,distributionQuantile};
if(typeof document!=='undefined'){
 const el=id=>document.getElementById(id);let activeRank=payload.config.rank,model=null;
 el('annotation').textContent=[payload.gene.gene,payload.gene.product,payload.gene.contig+':'+payload.gene.start+'–'+payload.gene.end].filter(Boolean).join(' · ');
 function updateDistribution(){
  const maximum=Number(el('max-taxa').value),floor=Number(el('quality-floor').value);
  if(!el('max-taxa').value||!el('quality-floor').value||!Number.isInteger(maximum)||maximum<1||!Number.isFinite(floor)||floor<0||floor>1){el('stats').textContent='分類群数は1以上の整数、品質下限は0–1を指定してください。';return;}
  model=distributionGroups(payload.points,payload.host,activeRank,floor,el('quality-filter').checked,maximum);
  el('diagram').innerHTML=distributionSvg(model);
  for(const button of el('rank-controls').querySelectorAll('[data-rank]'))button.setAttribute('aria-pressed',String(button.dataset.rank===activeRank));
  el('stats').textContent='表示階層：'+distributionRankNames[activeRank]+' ／ 保存済み '+model.total.toLocaleString()+'参照種 ／ フィルター後 '+model.accepted.toLocaleString()+'種 ／ 図中 '+model.shown.toLocaleString()+'種 ／ '+model.groups.length+'/'+model.totalGroups+'分類群（非表示 '+model.hiddenGroups+'群）。宿主分類群を表示に含めます。';
  const unknown=payload.points.filter(p=>p.quality===null).length;
  el('quality-note').textContent='品質下限 '+floor+'。'+(model.filtered?'品質取得済みで下限以上の点から中央値・支持数を再計算しています。':'元の比較に使った全有効ヒットから中央値・支持数を計算しています。低品質の点も保持します。')+(unknown?' '+unknown.toLocaleString()+'種は品質未取得です。古いDBの一致率を推定して補っていません。':'');
  const fmt=v=>v==null?'—':Number(v.toFixed(3)).toString();
  el('summary').innerHTML=model.groups.map(g=>'<tr>'+[g.displayLabel,g.host?'はい':'いいえ',g.n,fmt(g.min),fmt(g.q1),fmt(g.median),fmt(g.q3),fmt(g.max),g.low,g.unknown].map(v=>'<td>'+distributionEscape(v)+'</td>').join('')+'</tr>').join('');
 }
 for(const button of el('rank-controls').querySelectorAll('[data-rank]'))button.onclick=()=>{activeRank=button.dataset.rank;updateDistribution();};
 el('apply').onclick=updateDistribution;el('quality-filter').onchange=updateDistribution;
 el('export').onclick=()=>{if(!model)return;const url=URL.createObjectURL(new Blob([distributionSvg(model)],{type:'image/svg+xml'})),a=document.createElement('a');a.href=url;a.download=payload.gene.orf_id+'_'+activeRank+'_hit_distribution.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
 updateDistribution();
}
