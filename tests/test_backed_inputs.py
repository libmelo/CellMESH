"""A10: disk representation and gene ordering must not alter scoring semantics."""
import hashlib

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import read_anndata, run_cell_mesh
from cellmesh.preprocess import _slice_expression, _build_celltype_pseudobulk, _compute_celltype_expr_frac
from cellmesh.score import compute_metabolite_availability, compute_sensor_scores


@pytest.fixture
def data():
    x=np.array([[4,9,1,0],[0,2,4,8],[3,7,2,1],[1,0,5,3],[0,0,2,5],[2,6,0,0]],dtype=float)
    adata=AnnData(x,obs=pd.DataFrame({'cell_type':['A','B','A','B','B','A'],
                                    'sample':['S1','S1','S2','S2','S2','S3']},
                                   index=[f'c{i}' for i in range(6)]),
                  var=pd.DataFrame(index=['G1','G2','R','UNUSED']))
    enzyme=pd.DataFrame(dict(metabolite=['M','M'],hmdb_id=['HMDB1']*2,gene=['G1','G2'],
                             role=['production']*2,reaction=['P']*2))
    sensor=pd.DataFrame(dict(metabolite=['M'],hmdb_id=['HMDB1'],sensor_gene=['R'],sensor_type=['Transporter']))
    return adata,enzyme,sensor


def as_storage(adata,kind):
    if kind!='dense':
        adata.X=getattr(sparse,kind+'_matrix')(adata.X)
    return adata


def assert_results(actual,expected):
    for name in ['events','sender_scores','receiver_scores','sample_events','sample_sender_scores',
                 'sample_receiver_scores','sample_validation','celltype_qc']:
        a,b=getattr(actual,name),getattr(expected,name)
        if b is None:
            assert a is None
        else:
            pd.testing.assert_frame_equal(a,b,check_exact=True)
    for key in ['n_perms_completed','null_scores_stored']:
        assert actual.events.attrs[key]==expected.events.attrs[key]
    if 'sample_aware_null_scores' in expected.events.attrs:
        pd.testing.assert_frame_equal(actual.events.attrs['sample_aware_null_scores'],
                                      expected.events.attrs['sample_aware_null_scores'],check_exact=True)


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('kind',['dense','csr','csc'])
@pytest.mark.parametrize('reverse',[False,True])
@pytest.mark.parametrize('layer',[None,'counts'])
@pytest.mark.parametrize('n_perms',[0,5])
def test_backed_matches_memory(data,tmp_path,mode,kind,reverse,layer,n_perms):
    adata,enzyme,sensor=data
    if reverse:
        adata=adata[:,::-1].copy()
    as_storage(adata,kind)
    if layer:
        adata.layers[layer]=adata.X.copy()
        adata.X=adata.X*0 # chosen layer must supply values
    path=tmp_path/'data.h5ad';adata.write_h5ad(path)
    original=hashlib.sha256(path.read_bytes()).digest()
    options=dict(sample_mode=mode,sample_key='sample',min_cells=1,layer=layer,
                 n_perms=n_perms,random_state=3,store_null_scores=True)
    expected=run_cell_mesh(adata,enzyme,sensor,**options)
    backed=read_anndata(path,backed='r')
    try:
        actual=run_cell_mesh(backed,enzyme,sensor,**options)
        assert_results(actual,expected)
        assert backed.isbacked and backed.file.is_open
        pd.testing.assert_frame_equal(backed.obs,adata.obs)
        pd.testing.assert_frame_equal(backed.var,adata.var)
    finally:
        backed.file.close()
    assert hashlib.sha256(path.read_bytes()).digest()==original


@pytest.mark.parametrize('kind',['dense','csr','csc'])
def test_standalone_preprocessing_and_scoring(data,tmp_path,kind):
    adata,enzyme,sensor=data
    adata=as_storage(adata[:,::-1].copy(),kind)
    path=tmp_path/'input.h5ad';adata.write_h5ad(path)
    backed=read_anndata(path,backed='r')
    try:
        for function in [_build_celltype_pseudobulk,_compute_celltype_expr_frac]:
            pd.testing.assert_frame_equal(function(backed),function(adata),check_exact=True)
        a=compute_metabolite_availability(backed,enzyme,min_cells=1)
        b=compute_metabolite_availability(adata,enzyme,min_cells=1)
        pd.testing.assert_frame_equal(a['availability'],b['availability'],check_exact=True)
        pd.testing.assert_frame_equal(compute_sensor_scores(backed,sensor,min_cells=1),
                                      compute_sensor_scores(adata,sensor,min_cells=1),check_exact=True)
    finally:
        backed.file.close()


@pytest.mark.parametrize('kind',['dense','csr','csc'])
@pytest.mark.parametrize('rows',[np.array([4,0,4,2]),np.array([True,False,True,False,True,False]),
                                 slice(None,None,-1),np.array([],dtype=int)])
def test_disk_slice_restores_rows_columns_and_duplicates(data,tmp_path,kind,rows,monkeypatch):
    import cellmesh.preprocess as preprocess
    monkeypatch.setattr(preprocess,'_EXPRESSION_CHECK_BLOCK_SIZE',4)
    adata,_,_=data;adata=as_storage(adata,kind)
    path=tmp_path/'input.h5ad';adata.write_h5ad(path)
    backed=read_anndata(path,backed='r')
    try:
        cols=np.array([3,0,2,0])
        actual=_slice_expression(backed.X,rows,cols)
        expected=adata.X[rows,:][:,cols]
        if sparse.issparse(expected):
            assert sparse.issparse(actual)
            np.testing.assert_array_equal(actual.toarray(),expected.toarray())
        else:
            np.testing.assert_array_equal(actual,expected)
    finally:
        backed.file.close()


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
def test_backed_threaded_permutations_equal_serial(data,tmp_path,mode):
    adata,enzyme,sensor=data;adata=adata[:,::-1].copy()
    path=tmp_path/'data.h5ad';adata.write_h5ad(path)
    backed=read_anndata(path,backed='r')
    options=dict(sample_mode=mode,sample_key='sample',min_cells=1,n_perms=7,random_state=7,store_null_scores=True)
    try:
        assert_results(run_cell_mesh(backed,enzyme,sensor,n_jobs=2,**options),
                       run_cell_mesh(adata,enzyme,sensor,n_jobs=1,**options))
    finally:
        backed.file.close()


@pytest.mark.parametrize('kind',['dense','csr','csc'])
@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
def test_backed_violins_match_memory(data,tmp_path,kind,mode):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cellmesh import plot_metabolite_secretion_violin,plot_receptor_expression_violin
    adata,enzyme,sensor=data;adata=as_storage(adata[:,::-1].copy(),kind)
    path=tmp_path/'data.h5ad';adata.write_h5ad(path)
    result=run_cell_mesh(adata,enzyme,sensor,sample_mode=mode,sample_key='sample',min_cells=1,n_perms=0)
    backed=read_anndata(path,backed='r')
    try:
        for plot,extra in [(plot_metabolite_secretion_violin,{}),(plot_receptor_expression_violin,{'receptor_gene':'R'})]:
            actual=plot(result,backed,hmdb_id='HMDB1',**extra)
            expected=plot(result,adata,hmdb_id='HMDB1',**extra)
            pd.testing.assert_frame_equal(actual['plot_data'],expected['plot_data'],check_exact=True)
            pd.testing.assert_frame_equal(actual['summary'],expected['summary'],check_exact=True)
    finally:
        backed.file.close();plt.close('all')


@pytest.mark.parametrize('kind',['dense','csr','csc'])
@pytest.mark.parametrize('value',[-1,np.nan,np.inf])
def test_backed_invalid_relevant_expression_still_rejected(data,tmp_path,kind,value):
    adata,enzyme,sensor=data
    adata.X[0,0]=value
    adata=as_storage(adata[:,::-1].copy(),kind)
    path=tmp_path/'input.h5ad';adata.write_h5ad(path)
    backed=read_anndata(path,backed='r')
    try:
        with pytest.raises(ValueError,match='finite and non-negative'):
            run_cell_mesh(backed,enzyme,sensor,min_cells=1,n_perms=0)
    finally:
        backed.file.close()


@pytest.mark.parametrize('kind',['dense','csr','csc'])
def test_backed_sample_materialization_uses_only_requested_layer(data,tmp_path,kind):
    from cellmesh.preprocess import _sample_expression_adata
    adata,_,_=data;as_storage(adata,kind)
    adata.layers['counts']=adata.X.copy()
    adata.layers['ignored']=adata.X.copy()
    adata.raw=adata.copy()
    path=tmp_path/'input.h5ad';adata.write_h5ad(path)
    backed=read_anndata(path,backed='r')
    try:
        mask=backed.obs['sample'].eq('S2').to_numpy()
        actual=_sample_expression_adata(backed,mask,'counts')
        assert not actual.isbacked
        assert actual.X is None and actual.raw is None
        assert list(actual.layers)==['counts']
        pd.testing.assert_frame_equal(actual.obs,adata.obs.loc[mask])
        pd.testing.assert_frame_equal(actual.var,adata.var)
        matrix=actual.layers['counts'];expected=adata.layers['counts'][mask,:]
        if sparse.issparse(expected):
            assert sparse.issparse(matrix)
            np.testing.assert_array_equal(matrix.toarray(),expected.toarray())
        else: np.testing.assert_array_equal(matrix,expected)
    finally:
        backed.file.close()


def test_dense_disk_row_blocks_only_read_requested_gene_columns(data,tmp_path,monkeypatch):
    import h5py
    import cellmesh.preprocess as preprocess
    adata,_,_=data
    path=tmp_path/'matrix.h5'
    with h5py.File(path,'w') as file:
        source=file.create_dataset('X',data=adata.X)
        accesses=[]
        class Recorded:
            shape=source.shape
            def __getitem__(self,key):
                accesses.append(key);return source[key]
        monkeypatch.setattr(preprocess,'_EXPRESSION_CHECK_BLOCK_SIZE',4)
        actual=_slice_expression(Recorded(),np.array([5,0,3]),np.array([2,0]))
        np.testing.assert_array_equal(actual,adata.X[[5,0,3],:][:,[2,0]])
        assert len(accesses)==3
        for rows,columns in accesses:
            assert isinstance(rows,slice) and rows.stop-rows.start<=2
            np.testing.assert_array_equal(columns,[0,2])


@pytest.mark.parametrize('kind',['dense','csr','csc'])
@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('layer',[None,'counts'])
def test_backed_two_axis_reordered_views_match_memory(data,tmp_path,kind,mode,layer):
    adata,enzyme,sensor=data;as_storage(adata,kind)
    if layer: adata.layers[layer]=adata.X.copy()
    path=tmp_path/'view.h5ad';adata.write_h5ad(path)
    backed=read_anndata(path,backed='r')
    try:
        rows=[5,1,3,0];cols=[2,1,0,3]
        view=backed[rows,cols]
        expected=adata[rows,cols].copy()
        options=dict(sample_key='sample',sample_mode=mode,layer=layer,min_cells=1,n_perms=3,random_state=11)
        assert_results(run_cell_mesh(view,enzyme,sensor,**options),run_cell_mesh(expected,enzyme,sensor,**options))
        assert view.is_view and view.isbacked and backed.file.is_open
    finally: backed.file.close()
