"""B01: duplicate coordinates must be combined before overflow-prone operations."""
import hashlib

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import (read_anndata, run_cell_mesh, plot_metabolite_secretion_violin,
                      plot_receptor_expression_violin)
from cellmesh.preprocess import _compute_celltype_expr_frac


def case(dtype, kind):
    value = {'int8': 80, 'uint8': 128, 'int16': 20000, 'float32': 80, 'float64': 80}[dtype]
    # Include duplicates for both the enzyme and receptor, and explicit zeros.
    matrix = sparse.csr_matrix((np.array([value, value, value, value, 0, 2, 1] * 2, dtype=dtype),
                               np.array([0, 0, 1, 1, 1, 0, 1] * 2),
                               np.array([0, 5, 7, 12, 14])), shape=(4, 2))
    matrix = getattr(sparse, kind)(matrix)
    data = AnnData(matrix, obs=pd.DataFrame({'cell_type': ['A', 'B'] * 2,
                                            'sample': ['s1', 's1', 's2', 's2']}, index=list('abcd')),
                   var=pd.DataFrame(index=['G', 'R']))
    enzyme = pd.DataFrame(dict(metabolite=['M'], hmdb_id=['HMDB1'], gene=['G'],
                              reaction=['p'], role=['production']))
    sensor = pd.DataFrame(dict(metabolite=['M'], hmdb_id=['HMDB1'], sensor_gene=['R'],
                              sensor_type=['Transporter']))
    return data, enzyme, sensor


def safe_dense(matrix):
    # Independent reference: accumulate stored entries in float64 without sparse
    # duplicate reduction or a low-precision toarray() operation.
    coo = matrix.tocoo(copy=True)
    result = np.zeros(matrix.shape, dtype=np.float64)
    np.add.at(result, (coo.row, coo.col), coo.data.astype(np.float64))
    return result


@pytest.mark.parametrize('kind', ['csr_matrix', 'csc_matrix', 'csr_array', 'csc_array'])
@pytest.mark.parametrize('dtype', ['int8', 'uint8', 'int16', 'float32', 'float64'])
@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_duplicate_coordinates_match_dense(dtype, kind, mode):
    data, enzyme, sensor = case(dtype, kind)
    raw = data.X
    snapshot = [v.copy() for v in (raw.data, raw.indices, raw.indptr)]
    reference = AnnData(safe_dense(raw), obs=data.obs.copy(), var=data.var.copy())
    options = dict(min_cells=1, min_expr_frac=.1, n_perms=5, random_state=7,
                   sample_mode=mode, sample_key='sample', store_null_scores=True)
    actual = run_cell_mesh(data, enzyme, sensor, **options)
    expected = run_cell_mesh(reference, enzyme, sensor, **options)
    for key in ['events', 'sender_scores', 'receiver_scores', 'sample_events']:
        a, b = getattr(actual, key), getattr(expected, key)
        if b is not None:
            pd.testing.assert_frame_equal(a, b, check_exact=True)
    pd.testing.assert_frame_equal(_compute_celltype_expr_frac(data),
                                  _compute_celltype_expr_frac(reference), check_exact=True)
    if mode == 'sample_aware':
        pd.testing.assert_frame_equal(actual.events.attrs['sample_aware_null_scores'],
                                      expected.events.attrs['sample_aware_null_scores'], check_exact=True)
    for value, before in zip((raw.data, raw.indices, raw.indptr), snapshot):
        np.testing.assert_array_equal(value, before)
    assert not raw.has_canonical_format


@pytest.mark.parametrize('kind', ['csr_matrix', 'csc_matrix'])
@pytest.mark.parametrize('backed', [False, True])
@pytest.mark.parametrize('layer', [None, 'counts'])
@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_duplicate_layers_disk_and_violins(tmp_path, kind, backed, layer, mode):
    import matplotlib.pyplot as plt
    data, enzyme, sensor = case('int8', kind)
    reference = AnnData(safe_dense(data.X), obs=data.obs.copy(), var=data.var.copy())
    if layer:
        data.layers[layer] = data.X.copy()
        data.X = sparse.csr_matrix(data.shape, dtype=np.int8)
        reference.layers[layer] = reference.X.copy()
        reference.X = np.zeros(reference.shape)
    path = tmp_path / 'data.h5ad'
    data.write_h5ad(path)
    original = hashlib.sha256(path.read_bytes()).digest()
    target = read_anndata(path, backed='r') if backed else data
    options = dict(min_cells=1, min_expr_frac=.1, n_perms=3, n_jobs=2,
                   layer=layer, sample_mode=mode, sample_key='sample')
    try:
        actual = run_cell_mesh(target, enzyme, sensor, **options)
        expected = run_cell_mesh(reference, enzyme, sensor, **options)
        pd.testing.assert_frame_equal(actual.events, expected.events, check_exact=True)
        for plot, extra in [(plot_metabolite_secretion_violin, {}),
                            (plot_receptor_expression_violin, {'receptor_gene': 'R'})]:
            a = plot(actual, target, hmdb_id='HMDB1', **extra)
            b = plot(expected, reference, hmdb_id='HMDB1', **extra)
            pd.testing.assert_frame_equal(a['plot_data'], b['plot_data'], check_exact=True)
    finally:
        plt.close('all')
        if backed:
            target.file.close()
    assert hashlib.sha256(path.read_bytes()).digest() == original


@pytest.mark.parametrize('value', [-80, np.nan, np.inf])
def test_invalid_raw_duplicates_are_not_hidden(value):
    data, enzyme, sensor = case('float64', 'csr_matrix')
    data.X.data[0] = value
    with pytest.raises(ValueError, match='finite and non-negative'):
        run_cell_mesh(data, enzyme, sensor, min_cells=1)


@pytest.mark.parametrize('kind', ['csr_matrix', 'csc_matrix'])
def test_float32_duplicates_are_promoted_before_violin_densification(kind):
    import matplotlib.pyplot as plt
    data, enzyme, sensor = case('float32', kind)
    data.X.data[data.X.data == 80] = np.float32(2e38)
    reference = AnnData(safe_dense(data.X), obs=data.obs.copy(), var=data.var.copy())
    result = run_cell_mesh(data, enzyme, sensor, min_cells=1)
    try:
        plotted = plot_receptor_expression_violin(result, data, hmdb_id='HMDB1', receptor_gene='R')
        np.testing.assert_array_equal(plotted['plot_data']['receptor_expression'], reference.X[:, 1])
        assert np.isfinite(plotted['plot_data']['receptor_expression']).all()
    finally:
        plt.close('all')


def test_duplicate_reduction_overflow_is_still_rejected():
    data, enzyme, sensor = case('float64', 'csr_matrix')
    data.X.data[data.X.data == 80] = 1e308
    with pytest.raises(ValueError, match='finite and non-negative'):
        run_cell_mesh(data, enzyme, sensor, min_cells=1, n_perms=3)
