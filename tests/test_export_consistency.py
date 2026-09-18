"""A09: a successful export describes exactly one run; invalid exports preserve it."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from test_result_export import _result


def snapshot(directory):
    return {p.name: p.read_bytes() for p in directory.iterdir() if p.is_file()}


@pytest.mark.parametrize('first', ['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('second', ['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('first_empty', [False,True])
@pytest.mark.parametrize('second_empty', [False,True])
def test_reexport_matches_current_file_set(tmp_path,first,second,first_empty,second_empty):
    prefix=tmp_path/'analysis'
    _result(first,production_prior=not first_empty).to_csv(prefix)
    (tmp_path/'analysis.notes.txt').write_text('user notes')
    (tmp_path/'other.sample_events.csv').write_text('other prefix')
    (tmp_path/'Enzyme.csv').write_text('database')
    result=_result(second,production_prior=not second_empty)
    result.to_csv(prefix)
    manifest=json.loads((tmp_path/'analysis.manifest.json').read_text())
    assert set(snapshot(tmp_path)) == set(manifest['files']) | {'analysis.notes.txt','other.sample_events.csv','Enzyme.csv'}
    for name,entry in manifest['tables'].items():
        frame=pd.read_csv(tmp_path/entry['file'])
        assert len(frame)==entry['rows'] and frame.columns.tolist()==entry['columns']
    assert (tmp_path/'analysis.notes.txt').read_text()=='user notes'
    assert (tmp_path/'other.sample_events.csv').read_text()=='other prefix'
    assert (tmp_path/'Enzyme.csv').read_text()=='database'


def test_legacy_export_without_manifest_and_untrusted_manifest_paths(tmp_path):
    prefix=tmp_path/'run'
    _result('sample_aware').to_csv(prefix)
    (tmp_path/'run.manifest.json').unlink()
    _result().to_csv(prefix)
    assert not (tmp_path/'run.sample_events.csv').exists()
    outside=tmp_path/'keep.txt';outside.write_text('keep')
    (tmp_path/'run.manifest.json').write_text(json.dumps({'files':['keep.txt','../escape']}))
    _result().to_csv(prefix)
    assert outside.read_text()=='keep'


@pytest.mark.parametrize('value',[float('nan'),float('inf'),-float('inf'),np.float32('nan'),np.float64('inf')])
@pytest.mark.parametrize('nested',[False,True])
def test_strict_json_failure_preserves_existing_export(tmp_path,value,nested):
    result=_result('sample_aware');result.to_csv(tmp_path/'run');before=snapshot(tmp_path)
    changed=_result();changed.parameters['custom']={'items':[value]} if nested else value
    with pytest.raises(ValueError,match='JSON'):
        changed.to_csv(tmp_path/'run')
    assert snapshot(tmp_path)==before
    assert not any(p.is_dir() for p in tmp_path.iterdir())


def test_valid_numpy_parameters_are_strict_json(tmp_path):
    result=_result();result.parameters['custom']=[np.int64(7),np.float32(.5),np.bool_(True)]
    result.to_csv(tmp_path/'run')
    def reject(value):
        raise AssertionError(value)
    parameters=json.loads((tmp_path/'run.parameters.json').read_text(),parse_constant=reject)
    assert parameters['custom']==[7,.5,True]


@pytest.mark.parametrize('failure',['header','serialize','replace','cleanup'])
def test_failure_keeps_old_file_contents(tmp_path,monkeypatch,failure):
    import cellmesh.core as core
    _result('sample_aware').to_csv(tmp_path/'run');before=snapshot(tmp_path)
    result=_result()
    if failure=='header':
        result.events.insert(0,'hmdb_id',result.events.hmdb_id,allow_duplicates=True)
    elif failure=='serialize':
        original=pd.DataFrame.to_csv
        calls=0
        def fail(*args,**kwargs):
            nonlocal calls
            calls+=1
            if calls==2: raise OSError('injected')
            return original(*args,**kwargs)
        monkeypatch.setattr(pd.DataFrame,'to_csv',fail)
    elif failure=='replace':
        original=core.os.replace
        calls=0
        def fail(*args,**kwargs):
            nonlocal calls
            calls+=1
            if calls==2: raise OSError('injected')
            return original(*args,**kwargs)
        monkeypatch.setattr(core.os,'replace',fail)
    else:
        original=Path.unlink
        failed=False
        def fail(path,*args,**kwargs):
            nonlocal failed
            if path.name=='run.sample_events.csv' and not failed:
                failed=True;raise OSError('injected')
            return original(path,*args,**kwargs)
        monkeypatch.setattr(Path,'unlink',fail)
    with pytest.raises((ValueError,OSError)):
        result.to_csv(tmp_path/'run')
    assert snapshot(tmp_path)==before
    assert not any(p.is_dir() for p in tmp_path.iterdir())


@pytest.mark.parametrize('target',['run.events.csv','run.sample_events.csv','run.manifest.json'])
def test_symlink_target_rejected_without_touching_destination(tmp_path,target):
    actual=tmp_path/'user.txt';actual.write_text('original')
    (tmp_path/target).symlink_to(actual)
    with pytest.raises(ValueError,match='symlink'):
        _result().to_csv(tmp_path/'run')
    assert actual.read_text()=='original'


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
def test_empty_event_schema_and_permutation_metadata(mode):
    full=_result(mode);empty=_result(mode,production_prior=False)
    assert empty.events.empty
    assert empty.events.columns.tolist()==full.events.columns.tolist()
    assert empty.events.attrs['n_perms_completed']==0
    assert empty.events.attrs['null_scores_stored'] is False
    assert full.events.attrs['n_perms_completed']==3


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('production_prior',[False,True])
@pytest.mark.parametrize('store_null_scores',[False,True])
def test_no_permutations_has_explicit_metadata_and_fixed_schema(mode,production_prior,store_null_scores):
    result=_result(mode,production_prior=production_prior,n_perms=0,store_null_scores=store_null_scores)
    reference=_result(mode,production_prior=production_prior)
    assert result.events.columns.tolist()==reference.events.columns.tolist()
    assert result.events.attrs['n_perms_completed']==0
    assert result.events.attrs['null_scores_stored']==(store_null_scores and mode=='sample_aware')
    assert result.events[['perm_pvalue','fdr_global','fdr_sensor_type']].isna().all().all()


@pytest.mark.parametrize('has_previous',[False,True])
def test_manifest_publish_failure_rolls_back_all_changes(tmp_path,monkeypatch,has_previous):
    import cellmesh.core as core
    if has_previous:
        _result('sample_aware').to_csv(tmp_path/'run')
    before=snapshot(tmp_path)
    original=core.os.replace
    def fail(source,dest):
        if Path(source).name=='run.manifest.json' and Path(source).parent.name!='backup':
            raise OSError('manifest publish failure')
        return original(source,dest)
    monkeypatch.setattr(core.os,'replace',fail)
    with pytest.raises(OSError,match='manifest publish failure'):
        _result().to_csv(tmp_path/'run')
    assert snapshot(tmp_path)==before
    assert not any(p.is_dir() for p in tmp_path.iterdir())


def test_absent_optional_qc_is_cleaned_and_standard_namespace_is_reserved(tmp_path):
    result=_result();result.to_csv(tmp_path/'run')
    result.celltype_qc=None
    # Even a hand-created file with this exact reserved name cannot be
    # distinguished from a legacy export; this boundary is documented.
    (tmp_path/'run.sample_events.csv').write_text('manually edited standard filename')
    result.to_csv(tmp_path/'run')
    assert not (tmp_path/'run.celltype_qc.csv').exists()
    assert not (tmp_path/'run.sample_events.csv').exists()


def test_unsupported_parameter_and_directory_targets_leave_export_untouched(tmp_path):
    result=_result();result.to_csv(tmp_path/'run');before=snapshot(tmp_path)
    result.parameters['bad']=object()
    with pytest.raises(TypeError,match='not JSON serializable'):
        result.to_csv(tmp_path/'run')
    assert snapshot(tmp_path)==before
    del result.parameters['bad']
    (tmp_path/'run.sample_events.csv').mkdir()
    with pytest.raises(ValueError,match='regular file'):
        result.to_csv(tmp_path/'run')
    assert snapshot(tmp_path)==before
