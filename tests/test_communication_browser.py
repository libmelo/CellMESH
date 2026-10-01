"""Optional automated WebKit boundaries; no synthetic human-review artifact."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import pandas as pd
import pytest

from cellmesh import CellMeshResult, run_cell_mesh
from cellmesh.reporting import render_html_report
from test_communication_features import path_frame
from test_input_archive import inputs


@pytest.mark.skipif(os.environ.get('CELLMESH_BROWSER_TESTS')!='1',reason='optional system WebKit browser')
def test_single_static_matrix_and_expression_object_switching(inputs,tmp_path):
    frame=path_frame()
    # Exercise colorbar, exact zero, missingness and word wrapping in a browser.
    # These fixtures are only for automated boundaries, never real-data acceptance.
    frame['fdr_sensor_type']=([.1,.001,0.,float('nan')]*27)[:len(frame)]
    frame['metabolite']='A long metabolite name with several meaningful words'
    result=CellMeshResult(frame,pd.DataFrame(),pd.DataFrame(),{})
    tables={'coverage':{'summary':pd.DataFrame(),'details':pd.DataFrame()},
            'explanation':{'receiver_components':pd.DataFrame()}}
    text=render_html_report(result,selection=result.filter_events(source='all'),plot_top_n=None,tables=tables)
    large=tmp_path/'large.html';large.write_text(text)
    adata,enzyme,sensor=inputs
    enzyme=pd.concat([enzyme,enzyme.assign(hmdb_id='HMDB2',metabolite='Other')],ignore_index=True)
    sensor=pd.concat([sensor,sensor.assign(hmdb_id='HMDB2',metabolite='Other')],ignore_index=True)
    result=run_cell_mesh(adata,enzyme,sensor,cell_type_key='ct',sample_key='sample',min_cells=1,n_perms=0,retain_input=True)
    objects=result.to_html(tmp_path/'objects.html',selection=result.filter_events(source='all'))
    static_checks=tmp_path/'static_checks.js'
    static_checks.write_text('''(() => {
      const root=document.getElementById('interest-dotmap');
      const data=JSON.parse(document.getElementById('interest-dotmap-data').textContent);
      const figures=[...root.querySelectorAll('.dotmap-figure')];
      if(figures.length!==1 || !root.open)throw Error('one complete matrix must be visible without page scripts');
      if(!figures[0].querySelector('svg').getBoundingClientRect().height)throw Error('static figure missing');
      for(const record of data.records){
        if(!document.getElementById('interest-dotmap-event-'+record.id))throw Error('event requires JavaScript');
      }
      return JSON.stringify({status:'PASS',javascript:'disabled on load',figures:figures.length});
    })();''')
    server=subprocess.Popen(['/usr/bin/broadwayd',':24'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    env=dict(os.environ,GDK_BACKEND='broadway',BROADWAY_DISPLAY=':24',LIBGL_ALWAYS_SOFTWARE='1')
    tests=Path(__file__).parent
    try:
        time.sleep(.5)
        for path in [large,objects]:
            for width in [1400,420]:
                done=subprocess.run(['/usr/bin/python3',str(tests/'report_browser_webkit.py'),str(path),
                    str(tests/'report_communication_boundaries.js'),'--width',str(width)],
                    env=env,capture_output=True,text=True,timeout=60)
                assert done.returncode==0,done.stdout+'\n'+done.stderr
                assert json.loads(done.stdout.strip().splitlines()[-1])['status']=='PASS'
            done=subprocess.run(['/usr/bin/python3',str(tests/'report_browser_webkit.py'),str(path),
                str(static_checks),'--no-javascript'],env=env,capture_output=True,text=True,timeout=60)
            assert done.returncode==0,done.stdout+'\n'+done.stderr
            assert json.loads(done.stdout.strip().splitlines()[-1])['status']=='PASS'
    finally:
        server.terminate();server.wait(timeout=5)
