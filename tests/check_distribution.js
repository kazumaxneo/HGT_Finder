// Check exact summaries and the handlers embedded in a generated HTML page.
const fs=require('fs'),vm=require('vm'),assert=require('assert'),path=require('path');
const file=process.argv[2];if(!file)throw Error('Pass a generated ORF distribution HTML');
const html=fs.readFileSync(file,'utf8'),script=html.match(/<script>([\s\S]*?)<\/script>/)[1];
const payload=JSON.parse(script.match(/^const payload=(.*);\n/)[1]);
const pure=require('../gtdb_hgt/distribution.js');
const empty=pure.distributionGroups([],payload.host,'genus',.5,false,12);
assert.equal(empty.shown,0);assert.equal(empty.groups[0].n,0);assert(!/\b(?:NaN|Infinity)\b/.test(pure.distributionSvg(empty)));
const rows=fs.readFileSync(file.replace(/\.html$/,'_summary.tsv'),'utf8').trim().split('\n').map(r=>r.split('\t'));
const header=rows.shift();const expected=rows.map(values=>Object.fromEntries(header.map((h,i)=>[h,values[i]])));
for(const rank of ['phylum','class','order','family','genus']){
 const model=pure.distributionGroups(payload.points,payload.host,rank,.5,false,1000000);
 assert.equal(model.shown,payload.points.length);
 for(const g of model.groups.filter(g=>g.n)){
  const r=expected.find(r=>r.rank===rank&&r.taxon===g.tax.join(';'));assert(r);
  assert.equal(g.n,Number(r.species_count));
  for(const [a,b] of [['min','min_bitscore'],['max','max_bitscore'],['q1','q1_bitscore'],['median','median_bitscore'],['q3','q3_bitscore']])assert(Math.abs(g[a]-Number(r[b]))<1e-9);
 }
}
// Half-open bins, maximum-inclusive final bin, zero/empty support and count conservation.
const boundary={...empty,maxBits:4,groups:[{...empty.groups[0],n:6,points:[0,.99,1,2,3,4].map(bitscore=>({bitscore}))}]};
assert.deepEqual(pure.distributionHistogram(boundary,4).groups[0].counts,[2,1,1,2]);
assert.throws(()=>pure.distributionHistogram(boundary,0));
assert.throws(()=>pure.distributionHistogram(boundary,2.5));
assert.equal(pure.distributionHistogram(empty,30).groups[0].counts.reduce((a,b)=>a+b,0),0);
for(const rank of ['phylum','class','order','family','genus']){
 const model=pure.distributionGroups(payload.points,payload.host,rank,.5,false,1000000);
 const hist=pure.distributionHistogram(model,30);
 assert.equal(hist.groups.reduce((n,g)=>n+g.counts.reduce((a,b)=>a+b,0),0),payload.points.length);
 for(const g of hist.groups)assert.equal(g.counts.reduce((a,b)=>a+b,0),g.n);
 const bins=pure.distributionBinsTsv(model).trim().split('\n').slice(1);
 assert.equal(bins.reduce((n,line)=>n+Number(line.split('\t').at(-1)),0),payload.points.length);
 assert.equal(hist.upper,pure.distributionHistogram(pure.distributionGroups(payload.points,payload.host,'genus',.5,true,12),30).upper);
}
// Equal genus leaf names in different lineages remain independently selectable.
const lineagePoints=[{bitscore:100,quality:.7,taxonomy:['d__Bacteria','p__One','c__One','o__One','f__One','g__Same','s__One']},{bitscore:200,quality:.7,taxonomy:['d__Bacteria','p__Two','c__Two','o__Two','f__Two','g__Same','s__Two']}];
const duplicateLeaf=pure.distributionGroups(lineagePoints,payload.host,'genus',.5,false,12);
const leafGroups=duplicateLeaf.allGroups.filter(g=>g.label==='g__Same');assert.equal(leafGroups.length,2);assert.notEqual(leafGroups[0].displayLabel,leafGroups[1].displayLabel);
const oneLineage=pure.distributionSelect(duplicateLeaf,[leafGroups[0].key]);assert.equal(oneLineage.shown,1);assert.equal(oneLineage.groups[0].points[0].bitscore,200);
const sample=payload.points.slice(0,3).map((p,i)=>({...p,quality:[.3,.5,null][i]}));
assert.equal(pure.distributionGroups(sample,payload.host,'genus',.5,true,100).accepted,1);
assert.equal(pure.distributionGroups(sample,payload.host,'genus',.5,false,100).accepted,3);
const buttons=['phylum','class','order','family','genus'].map(rank=>({dataset:{rank},setAttribute(k,v){this[k]=v;}}));
const elements=new Map();function el(id){if(!elements.has(id))elements.set(id,{value:'',textContent:'',_html:'',checked:false,labels:[],set innerHTML(html){this._html=html;this.labels=[...html.matchAll(/<label data-search="([^"]*)"/g)].map(m=>({dataset:{search:m[1]},hidden:false}));},get innerHTML(){return this._html;},querySelectorAll(){return id==='rank-controls'?buttons:this.labels;}});return elements.get(id);}
el('plot-mode').value='overlay';el('bins').value='30';el('max-taxa').value=String(payload.config.max_taxa);el('quality-floor').value=String(payload.config.min_hit_quality);
const downloads=[];const ctx={document:{getElementById:el,createElement(){return {click(){downloads.push(this.download);}};}},Blob,URL:{createObjectURL(){return 'blob:test';},revokeObjectURL(){}},setTimeout(fn){fn();},console};vm.createContext(ctx);vm.runInContext(script,ctx);
for(const button of buttons){
 button.onclick();assert.equal(button['aria-pressed'],'true');
 const model=pure.distributionGroups(payload.points,payload.host,button.dataset.rank,payload.config.min_hit_quality,false,payload.config.max_taxa);
 const svg=el('diagram').innerHTML;
 assert.equal([...svg.matchAll(/data-count="(\d+)"/g)].reduce((n,m)=>n+Number(m[1]),0),model.shown);
 assert(!svg.includes('<circle'));assert(svg.includes('x: bit score / y: reference species count'));
 assert.equal((svg.match(/data-bin="true"/g)||[]).length,model.groups.length*30);
 assert(!/\b(?:NaN|Infinity)\b/.test(svg));
 el('export').onclick();assert(downloads.at(-1).includes(button.dataset.rank));
 el('export-bins').onclick();assert(downloads.at(-1).endsWith('_histogram_bins.tsv'));assert(downloads.at(-1).includes(button.dataset.rank));
}
// Choose a group beyond the automatic display limit, by its full lineage.
const catalog=pure.distributionGroups(payload.points,payload.host,'genus',.5,false,payload.config.max_taxa);
const host=catalog.allGroups.find(g=>g.host),other=catalog.allGroups[payload.config.max_taxa]||catalog.allGroups.find(g=>!g.host);
assert(other);
const countInSvg=()=>[...el('diagram').innerHTML.matchAll(/data-count="(\d+)"/g)].reduce((n,m)=>n+Number(m[1]),0);
el('select-host').onclick();assert.equal(countInSvg(),host.n);
const otherIndex=catalog.allGroups.findIndex(g=>g.key===other.key);
el('taxon-choices').onchange({target:{dataset:{index:String(otherIndex)},checked:true}});
assert.equal(countInSvg(),host.n+other.n);assert.equal((el('diagram').innerHTML.match(/data-bin="true"/g)||[]).length,60);
const chosen=pure.distributionSelect(catalog,[host.key,other.key]);
assert.equal(chosen.shown,host.n+other.n);assert.equal(chosen.groups.length,2);
assert(!/\b(?:NaN|Infinity)\b/.test(pure.distributionSvg({...chosen,mode:'panels'})));
const chosenTsv=pure.distributionBinsTsv(chosen).trim().split('\n').slice(1);
assert.equal(chosenTsv.reduce((n,r)=>n+Number(r.split('\t').at(-1)),0),host.n+other.n);
assert.deepEqual(new Set(chosenTsv.map(r=>r.split('\t')[0])),new Set([host.tax.join(';'),other.tax.join(';')]));
// Selection survives rank round trips, search, quality filtering and display-mode changes.
buttons[0].onclick();buttons[4].onclick();assert.equal(countInSvg(),host.n+other.n);
el('taxon-search').value=other.label.toLowerCase();el('taxon-search').oninput();assert(el('taxon-choices').labels.some(g=>g.hidden));assert.equal(countInSvg(),host.n+other.n);
el('plot-mode').value='panels';el('plot-mode').onchange();assert.equal(countInSvg(),host.n+other.n);
el('plot-mode').value='overlay';el('plot-mode').onchange();assert.equal(countInSvg(),host.n+other.n);
el('quality-filter').checked=true;el('quality-filter').onchange();
const filtered=pure.distributionGroups(payload.points,payload.host,'genus',.5,true,12);
assert.equal(countInSvg(),pure.distributionSelect(filtered,[host.key,other.key],catalog).shown);
assert.equal((el('diagram').innerHTML.match(/data-bin="true"/g)||[]).length,60);
el('quality-filter').checked=false;el('quality-filter').onchange();assert.equal(countInSvg(),host.n+other.n);
el('clear-taxa').onclick();assert.equal(countInSvg(),0);assert(el('diagram').innerHTML.includes('Select taxa'));assert(!/\b(?:NaN|Infinity)\b/.test(el('diagram').innerHTML));
el('reset-taxa').onclick();assert.equal(countInSvg(),catalog.shown);
el('taxon-search').value='';el('taxon-search').oninput();
el('quality-filter').checked=true;el('quality-filter').onchange();
assert.equal([...el('diagram').innerHTML.matchAll(/data-count="(\d+)"/g)].reduce((n,m)=>n+Number(m[1]),0),pure.distributionGroups(payload.points,payload.host,'genus',payload.config.min_hit_quality,true,payload.config.max_taxa).shown);
el('quality-filter').checked=false;el('quality-filter').onchange();
el('bins').value='10';el('apply').onclick();assert.equal((el('diagram').innerHTML.match(/data-bin="true"/g)||[]).length,pure.distributionGroups(payload.points,payload.host,'genus',.5,false,payload.config.max_taxa).groups.length*10);
buttons[0].onclick();assert.equal(el('bins').value,'10');assert.equal(el('quality-filter').checked,false);
el('bins').value='0';el('apply').onclick();assert(el('stats').textContent.includes('1–200'));el('bins').value='30';
el('max-taxa').value='0';el('apply').onclick();assert(el('stats').textContent.includes('1以上'));
const full=pure.distributionGroups(payload.points,payload.host,'genus',.5,false,payload.config.max_taxa);
fs.writeFileSync(file.replace(/\.html$/,'_preview.svg'),pure.distributionSvg(full));
const comparison=pure.distributionSelect(full,[host.key,other.key]);fs.writeFileSync(file.replace(/\.html$/,'_comparison.svg'),pure.distributionSvg(comparison));
fs.writeFileSync(file.replace(/\.html$/,'_script_check.js'),script);
console.log(JSON.stringify({species_points:payload.points.length,exact_five_rank_summaries:true,rank_handlers:true,exact_histogram_counts:true,quality_filter:true,svg_export:true,individual_taxon_selection:true,selection_survives_rank_filter:true,browser_layout_tested:false}));
