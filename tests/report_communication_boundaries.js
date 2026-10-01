(() => {
 const root=document.getElementById('interest-dotmap');
 const data=JSON.parse(document.getElementById('interest-dotmap-data').textContent);
 let passed=0;
 const check=(ok,name)=>{if(!ok)throw Error(name);passed++};
 check(!!root && root.open,'static matrix visible by default');
 check(!document.getElementById('interest-paths'),'four-column graph removed');
 check(!document.getElementById('interest-dotplot'),'duplicate old dotplot removed');
 check(!root.querySelector('select,button,input'),'matrix needs no interactive controls');
 check(!root.querySelector('.dotmap-reading').open,'detailed reading notes start collapsed');
 check(data.layout.font_sizes.tick===12 && data.layout.font_sizes.title===14,'declared readable font hierarchy');
 const compact=text=>text.replace(/\s+/g,'');
 for(let i=0;i<data.x_keys.length;i++){
   const label=compact(data.labels.columns[i]);
   check(data.x_keys[i].every(value=>label.includes(compact(value))),'full receiver-sensor identity on every x tick');
 }
 for(let i=0;i<data.y_keys.length;i++){
   const label=compact(data.labels.rows[i]);
   check(data.y_keys[i].every(value=>label.includes(compact(value))),'full sender-HMDB identity on every y tick');
 }
 const figures=[...root.querySelectorAll('.dotmap-figure')], ids=new Set();
 check(figures.length===1,'one static matrix contains every displayed event');
 check(JSON.stringify(data.x_fields)==='["receiver","sensor_gene"]','horizontal axis is receiver-sensor');
 check(JSON.stringify(data.y_fields)==='["sender","hmdb_id"]','vertical axis is sender-metabolite');
 check(data.labels.rotation===45 && data.labels.metabolite_label==='id','fixed 45 degree x labels and default HMDB IDs');
 check([...data.labels.columns,...data.labels.rows].every(label=>!label.includes('\n')),'single-line composite ticks');
 const box=figures[0],svg=box.querySelector('svg'),bounds=svg.getBoundingClientRect();
 check(bounds.width>0 && bounds.height>0,'matrix visible without interaction');
 check(box.scrollWidth<=box.clientWidth+2 && box.scrollHeight<=box.clientHeight+2,'full image without inner scrollbars');
 const xs=new Map(),ys=new Map();
 for(const record of data.records){
   const id=record.id, dot=data.dots.find(d=>d.id===id);
   check(!ids.has(id),'one marker per full identity');ids.add(id);
   const marker=document.getElementById('interest-dotmap-event-'+id);
   check(!!marker,'actual rendered event marker exists');
   check(JSON.stringify(data.x_keys[dot.x])===JSON.stringify([record.receiver,record.sensor_gene]),'correct receiver-sensor column');
   check(JSON.stringify(data.y_keys[dot.y])===JSON.stringify([record.sender,record.hmdb_id]),'correct sender-metabolite row');
   const use=marker.querySelector('use');
   const x=Number(use.getAttribute('x')),y=Number(use.getAttribute('y'));
   if(xs.has(dot.x))check(xs.get(dot.x)===x,'same composite column has same rendered x');else xs.set(dot.x,x);
   if(ys.has(dot.y))check(ys.get(dot.y)===y,'same composite row has same rendered y');else ys.set(dot.y,y);
 }
 check(new Set(xs.values()).size===data.x_keys.length,'different receiver-sensor columns remain distinct');
 check(new Set(ys.values()).size===data.y_keys.length,'different sender-metabolite rows remain distinct');
 const labels=[...svg.querySelectorAll('text,g[id*="-text_"]')];
 check(labels.length>0,'exported text or frozen glyph labels exist');
 for(const text of labels){
   const rect=text.getBoundingClientRect();
   check(rect.left>=bounds.left-2 && rect.right<=bounds.right+2 && rect.top>=bounds.top-2 && rect.bottom<=bounds.bottom+2,'full labels inside exported SVG');
 }
 check(ids.size===data.n_plotted && ids.size===data.records.length,'complete display scope');
 const table=JSON.parse(document.getElementById('interest-data').textContent);
 const tableIds=new Set(table.rows.map(r=>r.event_id));
 check(data.records.every(r=>tableIds.has(r.id)),'every dot corresponds to an original table record');
 for(const value of document.querySelectorAll('.overview-card strong'))check(value.scrollWidth<=value.clientWidth+2,'card value wraps within width');
 for(const group of document.querySelectorAll('.selection-expression-group')){
   const picker=group.querySelector('select');
   for(const option of picker.options){
     picker.value=option.value;picker.dispatchEvent(new Event('change',{bubbles:true}));
     const shown=[...group.querySelectorAll('[data-expression-index]')].filter(p=>!p.hidden);
     check(shown.length===1 && shown[0].dataset.expressionIndex===option.value,'exactly chosen expression object visible');
     check(!!shown[0].querySelector('svg'),'expression object contains figure');
   }
 }
 check(!window.PWNED,'labels remain inert');
 root.scrollIntoView();
 return JSON.stringify({status:'PASS',passed,records:data.records.length,figures:figures.length,columns:data.x_keys.length,rows:data.y_keys.length});
})();
