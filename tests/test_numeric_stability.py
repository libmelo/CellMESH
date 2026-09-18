"""A04: stable equivalent formulae and explicit, non-fatal underflow reports."""
import json
import logging
import warnings
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import run_cell_mesh
from cellmesh._numerics import (
    _abundance_weights, _event_score, _positive_product, _positive_reference_scores,
    _reaction_activity, _sender_base,
)
from cellmesh.score import compute_metabolite_availability
from cellmesh.plotting import (
    _cell_reaction_scores, _plot_grouped_violins,
    plot_metabolite_secretion_violin, plot_receptor_expression_violin,
)


@pytest.fixture
def inputs():
    adata = AnnData(np.array([[1., 2.], [2., 0.], [3., 4.], [4., 6.]]),
                    obs=pd.DataFrame({'cell_type': ['A','B','A','B'],
                                      'sample': ['s1','s1','s2','s2']},
                                     index=['c1','c2','c3','c4']),
                    var=pd.DataFrame(index=['G','S']))
    enzyme = pd.DataFrame({'metabolite':['M'], 'hmdb_id':['HMDB0000001'],
                           'gene':['G'], 'role':['production'], 'reaction':['r1']})
    sensor = pd.DataFrame({'metabolite':['M'], 'hmdb_id':['HMDB0000001'],
                          'sensor_gene':['S'], 'sensor_type':['Transporter']})
    return adata, enzyme, sensor


@pytest.mark.parametrize('values', [[1e-20], [1e-200], [0,1e-20,2e-20], [4.,9.], [0.,4.,9.], [0.,0.]])
def test_reaction_activity_matches_high_precision_mathematical_definition(values):
    with localcontext() as ctx:
        ctx.prec = 250
        logs = [(Decimal(1) + Decimal(str(x))).ln() for x in values]
        expected = float((sum(logs) / Decimal(len(logs))).exp() - Decimal(1))
    actual = _reaction_activity(np.array([values]))[0]
    assert actual == pytest.approx(expected, rel=3e-14, abs=0)
    if any(x>0 for x in values):
        assert actual > 0


def test_shared_cellwise_reaction_activity_agrees_at_tiny_values():
    values = np.array([[0,1e-20], [1e-200,0], [4.,9.]])
    definitions = pd.DataFrame({'genes':[['G','H']], 'direction':['product'], 'reaction':['r1']})
    components, _ = _cell_reaction_scores(values, ['G','H'], definitions)
    np.testing.assert_array_equal(components.production_score, _reaction_activity(values))


@pytest.mark.parametrize('p,c', [(2e-200,0.), (0.,0.), (0.,.5), (.25,.75), (1.,0.)])
def test_sender_base_avoids_unnecessary_square_underflow(p,c):
    expected = 0 if p+c==0 else p*(p/(p+c))
    assert _sender_base(np.array([p]),np.array([c]))[0] == expected


@pytest.mark.parametrize('sender,receiver,expected', [(2e-160,2e-200,2e-180), (0.,1e-200,0.), (.25,1.,.5)])
def test_event_score_avoids_underflowing_intermediate_product(sender,receiver,expected):
    assert _event_score(sender,receiver) == pytest.approx(expected,rel=2e-15,abs=0)


@pytest.mark.parametrize('operation', ['weights','product','reaction','base','reference'])
def test_actual_underflow_logs_and_continues_even_under_strict_numpy_and_warning_filters(operation, caplog):
    operations = {
        'weights': lambda: _abundance_weights([.5],2000),
        'product': lambda: _positive_product(1e-200,1e-200,'test product'),
        'reaction': lambda: _reaction_activity(np.array([[np.nextafter(0.,1.),0.]])),
        'base': lambda: _sender_base(np.array([1e-200]),np.array([1.])),
        'reference': lambda: _positive_reference_scores([1e-200,1e200],'mean','test reference')[0][0],
    }
    # Obtain the subnormal before enabling strict numpy handling.
    np.nextafter(0.,1.)
    if operation == 'reaction':
        minimum = np.nextafter(0.,1.)
        operations[operation] = lambda: _reaction_activity(np.array([[minimum,0.]]))
    with caplog.at_level(logging.WARNING,logger='cellmesh._numerics'), warnings.catch_warnings():
        warnings.simplefilter('error')
        with np.errstate(all='raise'):
            value = operations[operation]()
    assert np.asarray(value).item() == 0
    assert 'underflow' in caplog.text.lower()


def test_legitimate_zeros_and_representable_tiny_values_do_not_report(caplog):
    minimum = np.nextafter(0.,1.)
    with caplog.at_level(logging.WARNING,logger='cellmesh._numerics'):
        assert _reaction_activity(np.array([[minimum]]))[0] == minimum
        assert _positive_product(0.,1e-200,'test') == 0
        assert _positive_product(1e-200,1e-6,'test') > 0
        assert _sender_base(np.array([0.]),np.array([0.]))[0] == 0
        assert _event_score(0.,0.) == 0
    assert not caplog.records


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('storage',['dense','csr','csc'])
@pytest.mark.parametrize('scale',[1e-20,1e-200])
@pytest.mark.parametrize('reference',['mean','median'])
def test_tiny_single_gene_production_preserves_scale_invariant_scores_and_permutation_counts(inputs,mode,storage,scale,reference):
    adata,enzyme,sensor=inputs
    options=dict(sample_key='sample',sample_mode=mode,min_cells=1,n_perms=9,
                 random_state=7,pce_reference=reference,store_null_scores=True)
    expected=run_cell_mesh(adata,enzyme,sensor,**options)
    adata.X[:,0] *= scale
    if storage!='dense':
        adata.X=getattr(sparse,storage+'_matrix')(adata.X)
    result=run_cell_mesh(adata,enzyme,sensor,**options)
    keys=['sender','receiver','hmdb_id','sensor_gene']
    left=result.events.set_index(keys).sort_index()
    right=expected.events.set_index(keys).sort_index()
    np.testing.assert_allclose(left.cell_mesh_score,right.cell_mesh_score,rtol=3e-14,atol=0)
    np.testing.assert_array_equal(left[['perm_pvalue','fdr_global','fdr_sensor_type']],
                                  right[['perm_pvalue','fdr_global','fdr_sensor_type']])
    if mode=='sample_aware':
        np.testing.assert_allclose(result.events.attrs['sample_aware_null_scores'],
                                   expected.events.attrs['sample_aware_null_scores'],rtol=3e-14,atol=0)
    else:
        assert result.availability_results['metadata'].production_status.tolist()==['supported']
    assert 'numerical_diagnostics' not in result.parameters


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
@pytest.mark.parametrize('n_jobs',[1,2])
def test_true_underflow_completes_all_permutations_and_exports_one_aggregated_notice(inputs,mode,n_jobs,caplog,tmp_path):
    adata,enzyme,sensor=inputs
    with caplog.at_level(logging.WARNING,logger='cellmesh._numerics'),warnings.catch_warnings():
        warnings.simplefilter('error')
        result=run_cell_mesh(adata,enzyme,sensor,sample_key='sample',sample_mode=mode,
                             n_perms=5,n_jobs=n_jobs,store_null_scores=True,
                             min_cells=1,sender_abundance_exponent=2000)
    assert len([r for r in caplog.records if r.name=='cellmesh._numerics'])==1
    assert (result.events.cell_mesh_score==0).all()
    assert (result.events.perm_pvalue==1).all()
    assert result.events.attrs['n_perms_completed']==5
    diagnostic=result.parameters['numerical_diagnostics']
    assert diagnostic['underflow_detected'] is True
    assert diagnostic['policy']=='report_and_continue'
    assert any(item['stage']=='abundance weights' for item in diagnostic['stages'])
    assert result.events.attrs['numerical_diagnostics']==diagnostic
    prefix=tmp_path/'underflow'
    result.to_csv(prefix)
    assert json.loads((tmp_path/'underflow.parameters.json').read_text())['numerical_diagnostics']==diagnostic
    subsequent=run_cell_mesh(adata,enzyme,sensor,min_cells=1)
    assert 'numerical_diagnostics' not in subsequent.parameters


def test_threaded_permutation_underflow_diagnostics_match_serial(inputs):
    adata,enzyme,sensor=inputs
    adata.X[:,1]=[1e-200,1e200,1e-200,1e200]
    options=dict(sample_key='sample',sample_mode='sample_aware',n_perms=9,
                 min_cells=1,random_state=17,store_null_scores=True)
    serial=run_cell_mesh(adata,enzyme,sensor,n_jobs=1,**options)
    parallel=run_cell_mesh(adata,enzyme,sensor,n_jobs=2,**options)
    assert serial.parameters['numerical_diagnostics']==parallel.parameters['numerical_diagnostics']
    assert any(row['stage'].startswith('receiver') for row in serial.parameters['numerical_diagnostics']['stages'])
    pd.testing.assert_frame_equal(serial.events,parallel.events,check_exact=True)


def test_concurrent_public_runs_do_not_share_reports(inputs):
    adata,enzyme,sensor=inputs
    def run(exponent):
        return run_cell_mesh(adata,enzyme,sensor,min_cells=1,sender_abundance_exponent=exponent)
    with ThreadPoolExecutor(max_workers=2) as pool:
        bad,normal=list(pool.map(run,[2000,1]))
    assert 'numerical_diagnostics' in bad.parameters
    assert 'numerical_diagnostics' not in normal.parameters


def test_standalone_availability_reports_underflow(inputs,caplog):
    adata,enzyme,_=inputs
    result=compute_metabolite_availability(adata,enzyme,sender_abundance_exponent=2000)
    assert result['numerical_diagnostics']['underflow_detected']
    assert len(caplog.records)==1


def _violin(values, ax=None):
    return _plot_grouped_violins(pd.DataFrame({'group':['A']*len(values),'value':values}),
                                 group_col='group',value_col='value',group_order=['A'],
                                 color_values=pd.Series({'A':.5}),color_value_name='score',
                                 cmap='viridis',vmin=None,vmax=None,ylabel='value',title='test',
                                 show_median=True,ax=ax)


@pytest.mark.parametrize('values',[[0.,1e-9,2e-9],[0.,1e-200,2e-200],[1.,2.,3.],[1.,np.nextafter(1.,2.)]])
def test_nonconstant_violins_preserve_real_range(values):
    import matplotlib.pyplot as plt
    result=_violin(values)
    try:
        assert len(result['violin_bodies'])==1
        assert not result['constant_artists']
        assert not result['density_fallbacks']
        y=result['violin_bodies'][0].get_paths()[0].vertices[:,1]
        assert y.min()==min(values)
        assert y.max()==max(values)
    finally:
        plt.close(result['fig'])


@pytest.mark.parametrize('values',[[0.],[0.,0.],[1e-200,1e-200]])
def test_true_constants_still_use_actual_value_line(values):
    import matplotlib.pyplot as plt
    result=_violin(values)
    try:
        assert not result['violin_bodies']
        assert len(result['constant_artists'])==1
        np.testing.assert_array_equal(result['constant_artists'][0].get_segments()[0][:,1],[values[0],values[0]])
    finally:
        plt.close(result['fig'])


def test_kde_failure_uses_original_points_and_reports_fallback(monkeypatch):
    import matplotlib.pyplot as plt
    def fail(*args,**kwargs):
        raise np.linalg.LinAlgError('synthetic singular density')
    monkeypatch.setattr('scipy.stats.gaussian_kde',fail)
    result=_violin([0.,1e-9,2e-9])
    try:
        assert not result['constant_artists']
        assert not result['violin_bodies']
        assert result['density_fallbacks']==[{'cell_type':'A','reason':'synthetic singular density'}]
        np.testing.assert_array_equal(result['fallback_artists'][0].get_offsets()[:,1],[0.,1e-9,2e-9])
        assert len(result['median_artists'])==1
    finally:
        plt.close(result['fig'])


@pytest.mark.parametrize('kind',['metabolite','receptor'])
def test_public_violins_retain_tiny_single_cell_variation(inputs,kind):
    import matplotlib.pyplot as plt
    adata,enzyme,sensor=inputs
    adata.X*=1e-20
    result=run_cell_mesh(adata,enzyme,sensor,min_cells=1)
    if kind=='metabolite':
        plot=plot_metabolite_secretion_violin(result,adata,hmdb_id='HMDB0000001')
    else:
        plot=plot_receptor_expression_violin(result,adata,hmdb_id='HMDB0000001',receptor_gene='S')
    try:
        assert plot['violin_bodies']
        assert not plot['constant_artists']
    finally:
        plt.close(plot['fig'])


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
def test_pseudobulk_underflow_is_reported_before_zero_production_state(inputs,mode):
    adata,enzyme,sensor=inputs
    adata.obs['sample']='s1'
    adata.X[:,0]=[np.nextafter(0.,1.),0,0,0]
    result=run_cell_mesh(adata,enzyme,sensor,sample_key='sample',sample_mode=mode,n_perms=3,min_cells=1)
    stages={row['stage'] for row in result.parameters['numerical_diagnostics']['stages']}
    assert 'sender pseudobulk expression' in stages
    assert result.events.attrs['n_perms_completed']==3


def test_standalone_receiver_dataframe_retains_diagnostic(inputs):
    from cellmesh.score import compute_sensor_scores
    adata,_,sensor=inputs
    adata.X[:,1]=[1e-200,1e200,1e-200,1e200]
    result=compute_sensor_scores(adata,sensor)
    assert result.attrs['numerical_diagnostics']['underflow_detected']
    assert np.isfinite(result.sensor_score).all()


def test_public_metabolite_plot_retains_underflow_diagnostic(inputs):
    import matplotlib.pyplot as plt
    adata,enzyme,sensor=inputs
    adata.X[:,0]=np.nextafter(0.,1.)
    adata.X[:,1]=0
    enzyme['gene']='G;S'
    result=run_cell_mesh(adata,enzyme,sensor,min_cells=1)
    plot=plot_metabolite_secretion_violin(result,adata,hmdb_id='HMDB0000001')
    try:
        assert plot['numerical_diagnostics']['underflow_detected']
        assert plot['constant_artists']
    finally:
        plt.close(plot['fig'])


def test_scalar_and_vector_event_paths_use_identical_stable_roots():
    sender=np.array([0.,1e-300,2e-160,.1,1.])
    receiver=np.array([1.,1e-100,2e-200,.3,.5])
    np.testing.assert_array_equal(_event_score(sender,receiver),
                                  [_event_score(s,r) for s,r in zip(sender,receiver)])
