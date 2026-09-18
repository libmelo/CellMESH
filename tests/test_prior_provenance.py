"""A11: distinguish biological source annotations from input-file provenance."""
import hashlib
import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import load_cell_mesh_database,run_cell_mesh
from cellmesh.database import normalize_enzyme_database
from cellmesh.prior_coverage import summarize_prior_coverage


@pytest.fixture
def priors():
    enzyme=pd.DataFrame(dict(metabolite=['M','M','N'],hmdb_id=['HMDB1','HMDB1','HMDB2'],
                             gene=['G','E','G'],role=['production','export','production'],reaction=['p','e','q']))
    sensor=pd.DataFrame(dict(metabolite=['M'],hmdb_id=['HMDB1'],sensor_gene=['R'],sensor_type=['Transporter']))
    return enzyme,sensor


def direction_frame(enzyme):
    return enzyme.rename(columns={'role':'direction'}).assign(direction=lambda x:x.direction.map({'production':'product','export':'exporter'}))


@pytest.mark.parametrize('direction',[False,True])
@pytest.mark.parametrize('file_input',[False,True])
@pytest.mark.parametrize('sources',['absent','missing','supplied'])
def test_row_sources_never_invented_and_user_annotations_preserved(priors,tmp_path,direction,file_input,sources):
    enzyme,sensor=priors
    if direction: enzyme=direction_frame(enzyme)
    if sources=='missing': enzyme['source']=None
    if sources=='supplied': enzyme['source']=['User study',None,'Curated evidence']
    before=enzyme.copy(deep=True)
    if file_input:
        value=tmp_path/'Enzyme2.3.csv';enzyme.to_csv(value,index=False)
    else: value=enzyme
    loaded,_=load_cell_mesh_database(value,sensor)
    if sources=='absent': assert 'source' not in loaded
    elif sources=='missing': assert loaded.source.isna().all()
    else:
        assert loaded.source.iloc[0]=='User study' and loaded.source.iloc[2]=='Curated evidence'
        assert pd.isna(loaded.source.iloc[1])
    meta=loaded.attrs['input_provenance']
    assert meta['raw_rows']==len(enzyme) and meta['normalized_rows']==len(loaded)
    if file_input:
        assert meta['input_kind']=='user_file'
        assert meta['filename_version']=='2.3'
        assert meta['sha256']==hashlib.sha256(value.read_bytes()).hexdigest()
    else:
        assert meta['input_kind']=='dataframe'
        assert all(meta[key] is None for key in ['filename','path','filename_version','sha256'])
    pd.testing.assert_frame_equal(enzyme,before)
    assert enzyme.attrs==before.attrs


@pytest.mark.parametrize('defaults',[(True,True),(True,False),(False,True),(False,False)])
def test_default_and_custom_inputs_are_recorded_independently(priors,tmp_path,monkeypatch,defaults):
    import cellmesh.database as database
    enzyme,sensor=priors
    enzyme_path=tmp_path/'Enzyme9.2.csv';sensor_path=tmp_path/'Interaction8.1.csv'
    direction_frame(enzyme).to_csv(enzyme_path,index=False);sensor.to_csv(sensor_path,index=False)
    monkeypatch.setattr(database,'_default_database_paths',lambda:(enzyme_path,sensor_path))
    e,s=load_cell_mesh_database(None if defaults[0] else enzyme,None if defaults[1] else sensor)
    for table,is_default,path,version in [(e,defaults[0],enzyme_path,'9.2'),(s,defaults[1],sensor_path,'8.1')]:
        meta=table.attrs['input_provenance']
        assert meta['input_kind']==('default_file' if is_default else 'dataframe')
        assert meta['filename_version']==(version if is_default else None)
        assert meta['sha256']==(hashlib.sha256(path.read_bytes()).hexdigest() if is_default else None)
    assert 'source' not in e


def test_reloaded_dataframe_does_not_claim_stale_file_identity(priors,tmp_path):
    enzyme,sensor=priors
    path=tmp_path/'Enzyme1.0.csv';enzyme.to_csv(path,index=False)
    e,s=load_cell_mesh_database(path,sensor)
    original_meta=e.attrs['input_provenance'].copy()
    e.loc[0,'gene']='CHANGED'
    again,_=load_cell_mesh_database(e,s)
    assert again.attrs['input_provenance']['input_kind']=='dataframe'
    assert again.attrs['input_provenance']['sha256'] is None
    assert e.attrs['input_provenance']==original_meta


@pytest.mark.parametrize('name',['custom.csv','Enzyme4.2.csv','Enzyme4.2.csv.gz'])
def test_file_hashes_describe_bytes_and_only_filename_versions(priors,tmp_path,name):
    enzyme,sensor=priors;path=tmp_path/name
    enzyme.to_csv(path,index=False)
    table,_=load_cell_mesh_database(path,sensor)
    meta=table.attrs['input_provenance']
    assert meta['filename']==name
    assert meta['filename_version']==('4.2' if name.startswith('Enzyme') else None)
    assert meta['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()


def test_file_changed_during_load_is_rejected(priors,tmp_path,monkeypatch):
    import cellmesh.database as database
    enzyme,sensor=priors;path=tmp_path/'enzyme.csv';enzyme.to_csv(path,index=False)
    original=database._read_prior_table
    def changed(value,**kwargs):
        out=original(value,**kwargs)
        with open(value,'a') as stream: stream.write('\n')
        return out
    monkeypatch.setattr(database,'_read_prior_table',changed)
    with pytest.raises(ValueError,match='changed while loading'):
        load_cell_mesh_database(path,sensor)


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('file_input',[False,True])
def test_provenance_changes_do_not_change_scores_and_export_is_traceable(priors,tmp_path,mode,file_input):
    enzyme,sensor=priors;enzyme=direction_frame(enzyme)
    adata=AnnData(np.array([[4,2,1],[0,3,7],[3,1,2],[1,0,5]],dtype=float),
                  obs=pd.DataFrame({'cell_type':['A','B','A','B'],'sample':['S1','S1','S2','S2']},index=['c1','c2','c3','c4']),
                  var=pd.DataFrame(index=['G','E','R']))
    options=dict(sample_mode=mode,sample_key='sample',min_cells=1,n_perms=5,random_state=7)
    previous=run_cell_mesh(adata,enzyme.assign(source='packaged_enzyme_test'),sensor,**options)
    if file_input:
        ep=tmp_path/'Enzyme1.0.csv';sp=tmp_path/'Interaction1.0.csv'
        enzyme.to_csv(ep,index=False);sensor.to_csv(sp,index=False)
    else: ep,sp=enzyme,sensor
    result=run_cell_mesh(adata,ep,sp,**options)
    for name in ['events','sender_scores','receiver_scores']:
        pd.testing.assert_frame_equal(getattr(result,name),getattr(previous,name),check_exact=True)
    result.to_csv(tmp_path/'result')
    exported=json.loads((tmp_path/'result.parameters.json').read_text())
    assert exported['prior_inputs']==result.parameters['prior_inputs']
    assert exported['prior_inputs']['enzyme']['input_kind']==('user_file' if file_input else 'dataframe')
    assert result.parameters['export_weight']==previous.parameters['export_weight']==.2


def test_coverage_uses_unique_normalized_ids_and_explicit_denominators(priors):
    enzyme,sensor=priors
    enzyme=pd.concat([enzyme,enzyme.iloc[[1]].assign(hmdb_id=' hmdb1 ',metabolite='Alias')],ignore_index=True)
    summary=summarize_prior_coverage(enzyme,sensor)
    assert summary['enzyme_all']=={'exporter_hmdb_count':1,'total_hmdb_count':2,'percent':50.0}
    assert summary['shared_with_interaction']=={'exporter_hmdb_count':1,'total_hmdb_count':1,'percent':100.0}
    no_shared=summarize_prior_coverage(enzyme,sensor.assign(hmdb_id='HMDB999'))
    assert no_shared['shared_with_interaction']=={'exporter_hmdb_count':0,'total_hmdb_count':0,'percent':None}


def test_coverage_cli_is_reproducible_from_explicit_files(priors,tmp_path):
    enzyme,sensor=priors;ep=tmp_path/'enzyme.csv';sp=tmp_path/'sensor.csv'
    enzyme.to_csv(ep,index=False);sensor.to_csv(sp,index=False)
    result=subprocess.run([sys.executable,'-m','cellmesh.prior_coverage','--enzyme',str(ep),
                           '--interaction',str(sp)],check=True,capture_output=True,text=True)
    assert json.loads(result.stdout)==summarize_prior_coverage(ep,sp)


def test_direct_normalization_does_not_invent_source(priors):
    enzyme,_=priors
    assert 'source' not in normalize_enzyme_database(direction_frame(enzyme))


@pytest.mark.parametrize('file_input',[False,True])
def test_interaction_source_alias_is_preserved_separately_from_file_identity(priors,tmp_path,file_input):
    enzyme,sensor=priors
    sensor=sensor.rename(columns={'hmdb_id':'HMDB_ID','metabolite':'standard_metName',
                                   'sensor_gene':'Gene_name','sensor_type':'Annotation'})
    sensor['Database source']='Literature source'
    value=tmp_path/'Interaction3.0.csv' if file_input else sensor
    if file_input: sensor.to_csv(value,index=False)
    _,loaded=load_cell_mesh_database(enzyme,value)
    assert loaded.source.tolist()==['Literature source']
    assert loaded.attrs['input_provenance']['input_kind']==('user_file' if file_input else 'dataframe')


def test_empty_coverage_denominator_is_unavailable(priors):
    enzyme,sensor=priors
    summary=summarize_prior_coverage(enzyme.iloc[:0],sensor)
    for scope in ['enzyme_all','shared_with_interaction']:
        assert summary[scope]=={'exporter_hmdb_count':0,'total_hmdb_count':0,'percent':None}


def test_equal_filenames_with_different_contents_have_different_hashes(priors,tmp_path):
    enzyme,sensor=priors
    path=tmp_path/'Enzyme1.0.csv';enzyme.to_csv(path,index=False)
    first,_=load_cell_mesh_database(path,sensor)
    enzyme.assign(source='New evidence').to_csv(path,index=False)
    second,_=load_cell_mesh_database(path,sensor)
    assert first.attrs['input_provenance']['filename_version']==second.attrs['input_provenance']['filename_version']
    assert first.attrs['input_provenance']['sha256']!=second.attrs['input_provenance']['sha256']
