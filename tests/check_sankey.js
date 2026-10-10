// Run with: node tests/check_sankey.js
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const context={};vm.createContext(context);
vm.runInContext(fs.readFileSync('gtdb_hgt/sankey.js','utf8'),context);
const {aggregateRows,svgMarkup}=context;
function verify(rows,maxC=20,maxT=12,colors=null){
 const graph=aggregateRows(rows,maxC,maxT,colors);
 assert.equal(graph.links.reduce((s,l)=>s+l.value,0),rows.length);
 for(let id=0;id<graph.nodes.length;id++){
  const node=graph.nodes[id];
  const mass=graph.links.filter(l=>node.side?l.target===id:l.source===id).reduce((s,l)=>s+l.value,0);
  assert.equal(mass,node.count);
 }
 const svg=svgMarkup(graph);
 assert(!svg.includes('NaN')&&!svg.includes('Infinity'));
 return {graph,svg};
}
const rows=Array.from({length:45},(_,i)=>({contig:'c'+i,taxon:'d__Bacteria;p__P;g__T'+i,
 taxon_label:'g__T'+i,assignment_status:'assigned',host_taxon_match:i===0}));
const grouped=verify(rows,3,4);
assert(grouped.graph.nodes.some(n=>n.label.startsWith('その他のcontig')));
assert(grouped.graph.nodes.some(n=>n.label.startsWith('その他の分類群')));
assert.equal(grouped.graph.total,45);
const colors=context.makeTaxonColors(rows);
const first=aggregateRows(rows,20,50,colors),subset=aggregateRows(rows.slice(0,4),20,12,colors);
for(const node of subset.nodes.filter(n=>n.full))assert.equal(node.color,first.nodes.find(n=>n.full===node.full).color);
const unknown=verify([{contig:'<script>&"bad',assignment_status:'no_usable_hit',taxon_label:'<unsafe>'}],1,1);
assert(unknown.svg.includes('&lt;script&gt;'));
assert(!unknown.svg.includes('<script>'));
assert(svgMarkup(aggregateRows([])).includes('条件に合うORF'));
const ties=verify([
 {contig:'c',assignment_status:'ambiguous',ambiguity_kind:'host_competes',taxon_label:'host tie'},
 {contig:'c',assignment_status:'ambiguous',ambiguity_kind:'other_competes',taxon_label:'other tie'}]);
assert.equal(ties.graph.nodes.filter(n=>n.side===1).length,2);
assert(ties.graph.nodes.filter(n=>n.side===1).every(n=>n.color==='#97a4ac'));
console.log('Sankey conservation, aggregation limits, empty selection and escaping: OK');
const path=process.argv[2]||'runs/pilot_sankey_genus/sankey.html';
const folder=require('path').dirname(path);
if(fs.existsSync(path)){
 const html=fs.readFileSync(path,'utf8');
 const payload=JSON.parse(html.split('const payload=')[1].split(';\n// Pure aggregation')[0]);
 const rankValidation={};
 for(const rank of payload.config.available_ranks||[payload.config.rank]){
  const ranked=context.rowsForRank(payload,rank),colors=context.makeTaxonColors(ranked);
  const graph=verify(ranked,payload.config.max_contigs,payload.config.max_taxa,colors);
  const nif=verify(ranked.filter(r=>r.start>=291513&&r.end<=313504),payload.config.max_contigs,payload.config.max_taxa,colors);
  rankValidation[rank]={all:graph.graph.total,nif:nif.graph.total,host:ranked.filter(r=>r.host_taxon_match===true).length,near_ties:ranked.filter(r=>r.near_tie).length};
  fs.writeFileSync(folder+'/sankey_'+rank+'_all.svg',graph.svg);
  fs.writeFileSync(folder+'/sankey_'+rank+'_nif.svg',nif.svg);
 }
 fs.writeFileSync(folder+'/rank_graph_validation.json',JSON.stringify(rankValidation,null,2)+'\n');
 const colors=context.makeTaxonColors(payload.rows);
 const all=verify(payload.rows,payload.config.max_contigs,payload.config.max_taxa,colors);
 fs.writeFileSync(folder+'/sankey_all.svg',all.svg);
 const candidates=verify(payload.rows.filter(r=>r.in_candidate_region),payload.config.max_contigs,payload.config.max_taxa,colors);
 fs.writeFileSync(folder+'/sankey_candidates.svg',candidates.svg);
 const nif=verify(payload.rows.filter(r=>r.start>=291513&&r.end<=313504),payload.config.max_contigs,payload.config.max_taxa,colors);
 fs.writeFileSync(folder+'/sankey_nif.svg',nif.svg);
 const counts={};
 for(const r of payload.rows)counts[r.assignment_status]=(counts[r.assignment_status]||0)+1;
 const validation={all:all.graph.total,candidates:candidates.graph.total,nif:nif.graph.total,status:counts,host:payload.rows.filter(r=>r.host_taxon_match===true).length,near_ties:payload.rows.filter(r=>r.near_tie).length,flow_conservation:true,svg_finite:true,browser_interaction_tested:false};
 console.log(JSON.stringify(validation));
 fs.writeFileSync(folder+'/sankey_validation.json',JSON.stringify(validation,null,2)+'\n');
 fs.writeFileSync(folder+'/sankey_script_check.js',html.split('<script>')[1].split('</script>')[0]);
}
