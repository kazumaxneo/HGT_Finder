// Pure aggregation and SVG rendering, also used by the verification harness.
function rowsForRank(payload,rank){
 if(rank===payload.config.rank)return payload.rows;
 const assignments=payload.rank_assignments?.[rank];
 if(!assignments||assignments.length!==payload.rows.length)throw new Error('Missing or invalid rank assignments: '+rank);
 return payload.rows.map((row,i)=>({...row,...assignments[i]}));
}
function esc(value){return String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function taxonId(r){return r.assignment_status==='assigned'?'tax:'+r.taxon:'status:'+r.assignment_status+(r.assignment_status==='ambiguous'&&r.ambiguity_kind?':'+r.ambiguity_kind:'');}
function shortLabel(value,max){const chars=Array.from(value);return chars.length>max?chars.slice(0,max-1).join('')+'…':value;}
function makeTaxonColors(rows){
 const colors=['#d98328','#b05c8a','#bc6848','#8073b7','#a99532','#5286ac','#996334','#629347','#a46562','#587eb5'];
 const keys=[...new Set(rows.filter(r=>r.assignment_status==='assigned'&&r.host_taxon_match!==true).map(taxonId))].sort();
 const map=Object.create(null);keys.forEach((key,i)=>map[key]=colors[i%colors.length]);return map;
}
function aggregateRows(rows,maxContigs=20,maxTaxa=12,colorMap=null){
 colorMap=colorMap||makeTaxonColors(rows);
 const counts=new Map(),taxa=new Map();
 for(const r of rows){counts.set(r.contig,(counts.get(r.contig)||0)+1);const key=taxonId(r);
  if(!taxa.has(key))taxa.set(key,{key,label:r.taxon_label,count:0,host:r.host_taxon_match===true,status:r.assignment_status});
  taxa.get(key).count++;
 }
 const order=(a,b)=>b[1]-a[1]||a[0].localeCompare(b[0]);
 const ordered=[...counts].sort(order),keep=new Set(ordered.slice(0,maxContigs).map(v=>v[0]));
 const assigned=[...taxa.values()].filter(t=>t.status==='assigned').sort((a,b)=>Number(b.host)-Number(a.host)||b.count-a.count||a.key.localeCompare(b.key));
 const keptTaxa=new Set(assigned.slice(0,maxTaxa).map(t=>t.key));
 const nodes=[],contigMap=Object.create(null),taxonMap=Object.create(null);
 for(const [c,n] of ordered.slice(0,maxContigs)){contigMap[c]=nodes.length;nodes.push({label:c,count:n,side:0,color:'#527c96'});}
 if(ordered.length>maxContigs){const id=nodes.length;nodes.push({label:'その他のcontig（'+(ordered.length-maxContigs)+'）',count:0,side:0,color:'#859caa'});for(const [c,n] of ordered.slice(maxContigs)){contigMap[c]=id;nodes[id].count+=n;}}
 for(const t of assigned.filter(t=>keptTaxa.has(t.key))){taxonMap[t.key]=nodes.length;nodes.push({label:t.label,full:t.key.slice(4),count:t.count,side:1,color:t.host?'#148c82':colorMap[t.key]});}
 const extra=assigned.filter(t=>!keptTaxa.has(t.key));
 if(extra.length){const id=nodes.length;nodes.push({label:'その他の分類群（'+extra.length+'）',count:extra.reduce((s,t)=>s+t.count,0),side:1,color:'#a58b72'});for(const t of extra)taxonMap[t.key]=id;}
 for(const t of [...taxa.values()].filter(t=>t.status!=='assigned').sort((a,b)=>b.count-a.count||a.key.localeCompare(b.key))){taxonMap[t.key]=nodes.length;nodes.push({label:t.label,count:t.count,side:1,color:'#97a4ac'});}
 const links=new Map();
 for(const r of rows){const source=contigMap[r.contig],target=taxonMap[taxonId(r)],key=source+','+target;if(!links.has(key))links.set(key,{source,target,value:0});links.get(key).value++;}
 return {nodes,links:[...links.values()].sort((a,b)=>a.source-b.source||a.target-b.target),contigMap,taxonMap,total:rows.length};
}
function svgMarkup(graph){
 const width=1400,gap=14,top=64,bottom=30,left=260,right=985,nodeWidth=18;
 const leftNodes=graph.nodes.filter(n=>n.side===0),rightNodes=graph.nodes.filter(n=>n.side===1),n=Math.max(leftNodes.length,rightNodes.length,1);
 const height=Math.max(450,top+bottom+n*32+(n-1)*gap),available=height-top-bottom;
 if(!graph.total)return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 '+width+' 180"><rect width="100%" height="100%" fill="white"/><text x="30" y="90" font-family="sans-serif" font-size="18" fill="#496573">条件に合うORFがありません。</text></svg>';
 const scale=(available-(n-1)*gap)/graph.total,pos=new Map();
 for(const side of [0,1]){const entries=graph.nodes.map((node,id)=>({node,id})).filter(e=>e.node.side===side);let y=top+(available-(graph.total*scale+(entries.length-1)*gap))/2;for(const e of entries){pos.set(e.id,{x:side?right:left,y,h:e.node.count*scale,used:0});y+=e.node.count*scale+gap;}}
 let paths='';
 for(let i=0;i<graph.links.length;i++){const link=graph.links[i],a=pos.get(link.source),b=pos.get(link.target),h=link.value*scale,x0=a.x+nodeWidth,x1=b.x,mid=(x0+x1)/2,y0=a.y+a.used,y1=b.y+b.used;a.used+=h;b.used+=h;
  const d='M '+x0+' '+y0+' C '+mid+' '+y0+' '+mid+' '+y1+' '+x1+' '+y1+' L '+x1+' '+(y1+h)+' C '+mid+' '+(y1+h)+' '+mid+' '+(y0+h)+' '+x0+' '+(y0+h)+' Z';
  paths+='<path data-link="'+i+'" d="'+d+'" fill="'+graph.nodes[link.target].color+'" fill-opacity=".36"><title>'+esc(graph.nodes[link.source].label+' → '+graph.nodes[link.target].label+': '+link.value+' ORF')+'</title></path>';
 }
 let labels='';
 for(let id=0;id<graph.nodes.length;id++){const node=graph.nodes[id],p=pos.get(id),x=node.side?p.x+nodeWidth+10:p.x-10,anchor=node.side?'start':'end';labels+='<g><rect x="'+p.x+'" y="'+p.y+'" width="'+nodeWidth+'" height="'+p.h+'" fill="'+node.color+'"/><text x="'+x+'" y="'+(p.y+p.h/2)+'" dominant-baseline="middle" text-anchor="'+anchor+'" font-size="14" fill="#173140">'+esc(shortLabel(node.label,node.side?40:24))+' ('+node.count+')</text><title>'+esc((node.full||node.label)+': '+node.count+' ORF')+'</title></g>';}
 return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 '+width+' '+height+'" role="img" aria-label="Contigと類似先分類群のORF数を示すサンキープロット"><style>text{font-family:Hiragino Sans,Noto Sans CJK JP,sans-serif}path[data-link]{cursor:pointer}path[data-link]:hover{fill-opacity:.72}</style><rect width="100%" height="100%" fill="white"/><text x="'+left+'" y="30" font-size="16" fill="#496573">Contig</text><text x="'+right+'" y="30" font-size="16" fill="#496573">類似先分類群</text>'+paths+labels+'</svg>';
}
