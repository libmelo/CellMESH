"""A10: optional Loom dependency diagnostics and real reader behavior."""
import builtins
import importlib.util
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellmesh import read_anndata,run_cell_mesh


def test_missing_loom_dependency_has_install_hint(monkeypatch,tmp_path):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name=='loompy':
            raise ModuleNotFoundError("No module named 'loompy'",name='loompy')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    with pytest.raises(ImportError,match=r'cellmesh\[loom\].*pip install loompy'):
        read_anndata(tmp_path/'missing.loom',mode='loom')


def test_broken_transitive_dependency_is_not_mislabeled(monkeypatch,tmp_path):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name=='loompy':
            raise ModuleNotFoundError("broken internal dependency",name='other_dependency')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    with pytest.raises(ModuleNotFoundError,match='broken internal dependency') as error:
        read_anndata(tmp_path/'missing.loom',mode='loom')
    assert error.value.name=='other_dependency'


@pytest.fixture
def loom_file(tmp_path):
    loompy=pytest.importorskip('loompy')
    matrix=np.array([[4,0,3,1],[1,5,2,8]],dtype=np.float32)
    path=tmp_path/'input.loom'
    loompy.create(str(path),matrix,{'Gene':np.array(['G','R'])},
                  {'CellID':np.array(['c1','c2','c3','c4']),
                   'cell_type':np.array(['A','B','A','B']),
                   'sample':np.array(['S1','S1','S2','S2'])})
    return path,matrix


@pytest.mark.parametrize('as_sparse',[False,True])
@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
def test_real_loom_read_and_scoring_matches_memory(loom_file,as_sparse,mode):
    from anndata import AnnData
    path,matrix=loom_file
    with warnings.catch_warnings():
        warnings.simplefilter('error',FutureWarning)
        actual=read_anndata(path,mode='loom',sparse=as_sparse)
    np.testing.assert_array_equal(actual.X.toarray() if sparse.issparse(actual.X) else actual.X,matrix.T)
    assert actual.var_names.tolist()==['G','R']
    assert actual.obs_names.tolist()==['c1','c2','c3','c4']
    expected=AnnData(matrix.T.copy(),obs=actual.obs.copy(),var=actual.var.copy())
    enzyme=pd.DataFrame(dict(metabolite=['M'],hmdb_id=['HMDB1'],gene=['G'],role=['production'],reaction=['P']))
    sensor=pd.DataFrame(dict(metabolite=['M'],hmdb_id=['HMDB1'],sensor_gene=['R'],sensor_type=['Transporter']))
    kwargs=dict(sample_key='sample',sample_mode=mode,min_cells=1,n_perms=5,random_state=7)
    left=run_cell_mesh(actual,enzyme,sensor,**kwargs)
    right=run_cell_mesh(expected,enzyme,sensor,**kwargs)
    pd.testing.assert_frame_equal(left.events,right.events,check_exact=True)


def test_reader_file_errors_are_preserved(loom_file,tmp_path):
    path,_=loom_file
    broken=tmp_path/'broken.loom';broken.write_text('not a Loom file')
    with pytest.raises(OSError):
        read_anndata(broken,mode='loom')
    with pytest.raises(OSError):
        read_anndata(tmp_path/'missing.loom',mode='loom')


def test_old_anndata_import_fallback(loom_file,monkeypatch):
    import anndata
    path,_=loom_file
    original=importlib.util.find_spec
    monkeypatch.setattr(importlib.util,'find_spec',lambda name: None if name=='anndata.io' else original(name))
    sentinel=object();calls=[]
    def legacy(path,**kwargs):
        calls.append((path,kwargs));return sentinel
    monkeypatch.setitem(anndata.__dict__,'read_loom',legacy)
    assert read_anndata(path,mode='loom',sparse=False) is sentinel
    assert calls==[(path,{'sparse':False})]
