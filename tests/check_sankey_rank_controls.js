// Execute the generated page's actual button handlers with a minimal DOM stub.
// This checks state and filtering; it does not claim browser layout verification.
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const path=process.argv[2];
if(!path)throw new Error('Pass the generated sankey.html path');
const html=fs.readFileSync(path,'utf8'),script=html.split('<script>')[1].split('</script>')[0];
const payload=JSON.parse(html.split('const payload=')[1].split(';\n// Pure aggregation')[0]);
const ranks=payload.config.available_ranks;
assert.deepEqual(ranks,['phylum','class','order','family','genus']);
const buttons=ranks.map(rank=>({dataset:{rank},attributes:{},setAttribute(k,v){this.attributes[k]=v;}}));
const ids=['rank-controls','method','orf-tsv','contig-tsv','table-count','orf-body','from','to','scope','contig','detail','diagram','stats','apply','reset','export','clear-flow'];
const elements=Object.fromEntries(ids.map(id=>[id,{value:'',innerHTML:'',textContent:'',querySelectorAll(){return id==='rank-controls'?buttons:[];}}]));
elements.scope.value='all';
const downloads=[];
const context={document:{getElementById:id=>elements[id],createElement(){return {click(){downloads.push(this.download);}};}},Blob,
 URL:{createObjectURL(){return 'blob:test';},revokeObjectURL(){}},setTimeout:fn=>fn()};
vm.createContext(context);vm.runInContext(script,context);
const evaluate=s=>vm.runInContext(s,context);
const validations={};
for(const button of buttons){
 button.onclick();
 assert.equal(evaluate('activeRank'),button.dataset.rank);
 assert.equal(button.attributes['aria-pressed'],'true');
 assert.equal(buttons.filter(b=>b.attributes['aria-pressed']==='true').length,1);
 assert.equal(evaluate('graph.total'),payload.rows.length);
 assert(!elements.diagram.innerHTML.includes('NaN'));
 assert.equal(elements['orf-tsv'].href,'sankey_orfs_'+button.dataset.rank+'.tsv');
 elements.export.onclick();assert.equal(downloads.at(-1),'contig_taxonomy_'+button.dataset.rank+'_sankey.svg');
 const ranked=evaluate('rankRows');
 for(let i=0;i<ranked.length;i++)for(const k of ['orf_id','orf_score','window_score','candidate_regions'])
  assert.equal(ranked[i][k],payload.rows[i][k]);
 elements.contig.value='contig_1';elements.from.value='291513';elements.to.value='313504';elements.apply.onclick();
 assert.equal(evaluate('graph.total'),27);
 // Rank button clicks preserve the contig, coordinates and scope.
 button.onclick();assert.equal(evaluate('graph.total'),27);
 assert.equal(elements.from.value,'291513');assert.equal(elements.to.value,'313504');
 elements.scope.value='candidate';button.onclick();assert.equal(evaluate('graph.total'),27);
 validations[button.dataset.rank]={all_orfs:ranked.length,nif_orfs:27,host_nif:evaluate('current.filter(r=>r.host_taxon_match===true).length'),nif_status:evaluate('current.reduce((s,r)=>(s[r.assignment_status]=(s[r.assignment_status]||0)+1,s),{})')};
 elements.reset.onclick();
 assert.equal(evaluate('activeRank'),button.dataset.rank);
 assert.equal(evaluate('graph.total'),payload.rows.length);
 elements.scope.value='near';elements.scope.onchange();
 assert.equal(evaluate('graph.total'),ranked.filter(r=>r.near_tie).length);
 elements.reset.onclick();
}
const result={rank_button_handlers:true,filters_preserved:true,rank_specific_downloads:true,score_fields_preserved:true,browser_layout_tested:false,rank_counts:validations};
fs.writeFileSync(require('path').dirname(path)+'/rank_controls_validation.json',JSON.stringify(result,null,2)+'\n');
console.log(JSON.stringify(result));
