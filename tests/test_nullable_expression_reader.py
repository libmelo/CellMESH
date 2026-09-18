"""B05: pandas nullable numeric input must yield a usable numeric AnnData X."""
import numpy as np
import pandas as pd
import pytest

from cellmesh import read_anndata, run_cell_mesh


def priors():
    enzyme = pd.DataFrame(dict(metabolite=['M'], hmdb_id=['HMDB1'], gene=['G'], role=['production'], reaction=['p']))
    sensor = pd.DataFrame(dict(metabolite=['M'], hmdb_id=['HMDB1'], sensor_gene=['R'], sensor_type=['Transporter']))
    return enzyme, sensor


@pytest.mark.parametrize('mode,sep', [('csv', ','), ('tsv', '\t')])
@pytest.mark.parametrize('transpose', [False, True])
@pytest.mark.parametrize('compressed', [False, True])
@pytest.mark.parametrize('option', ['Int64', 'Float64', 'backend', 'mixed'])
@pytest.mark.parametrize('sample_mode', ['pooled_stratified', 'sample_aware'])
def test_nullable_reader_matches_native_scoring(tmp_path, mode, sep, transpose, compressed, option, sample_mode):
    frame = pd.DataFrame([[1, 2], [3, 4], [0, 3], [2, 0]], index=['01', '1', 'NA', 'c4'], columns=['G', 'R'])
    source = frame.T if transpose else frame
    path = tmp_path / ('input.' + mode + ('.gz' if compressed else ''))
    source.to_csv(path, sep=sep)
    if option == 'backend':
        options = {'dtype_backend': 'numpy_nullable'}
    elif option == 'mixed':
        options = {'dtype': {column: ('Int64' if i % 2 == 0 else np.float64)
                             for i, column in enumerate(source.columns)}}
    else:
        options = {'dtype': option}
    actual = read_anndata(path, mode=mode, transpose=transpose, **options)
    expected = read_anndata(path, mode=mode, transpose=transpose)
    assert actual.X.dtype.kind in 'iuf'
    assert actual.obs_names.tolist() == frame.index.tolist()
    assert actual.var_names.tolist() == frame.columns.tolist()
    np.testing.assert_array_equal(actual.X, expected.X)
    for data in (actual, expected):
        data.obs['cell_type'] = ['A', 'B', 'A', 'B']
        data.obs['sample'] = ['s1', 's1', 's2', 's2']
    options = dict(sample_key='sample', sample_mode=sample_mode, min_cells=1, n_perms=3, random_state=7)
    got, ref = [run_cell_mesh(data, *priors(), **options) for data in (actual, expected)]
    pd.testing.assert_frame_equal(got.events, ref.events, check_exact=True)


@pytest.mark.parametrize('dtype', ['Int64', 'Float32', 'Float64'])
@pytest.mark.parametrize('transpose', [False, True])
def test_missing_values_remain_missing_and_scoring_rejects(tmp_path, dtype, transpose):
    frame = pd.DataFrame([[1., 2.], [np.nan, 3.]], index=['a', 'b'], columns=['G', 'R'])
    path = tmp_path / 'missing.csv'
    (frame.T if transpose else frame).to_csv(path)
    data = read_anndata(path, mode='csv', dtype=dtype, transpose=transpose)
    assert data.X.dtype.kind == 'f'
    assert np.isnan(data.X[1, 0])
    data.obs['cell_type'] = ['A', 'B']
    with pytest.raises(ValueError, match='finite and non-negative'):
        run_cell_mesh(data, *priors(), min_cells=1)


@pytest.mark.parametrize('dtype,value', [('UInt64', 2**63 + 1), ('Int64', 2**53 + 1), ('Float32', .25)])
def test_homogeneous_nullable_dtype_and_values_are_preserved(tmp_path, dtype, value):
    path = tmp_path / 'exact.csv'
    path.write_text(f'cell,G,R\na,{value},{value}\n')
    data = read_anndata(path, mode='csv', dtype=dtype)
    assert data.X.dtype == pd.api.types.pandas_dtype(dtype).numpy_dtype
    assert data.X[0, 0] == value


@pytest.mark.parametrize('bad', ['bad', 'True', '1+2j'])
def test_invalid_text_is_not_coerced_to_missing_or_real(tmp_path, bad):
    path = tmp_path / 'bad.csv'
    path.write_text(f'cell,G,R\na,{bad},1\n')
    with pytest.raises((TypeError, ValueError)):
        read_anndata(path, mode='csv', dtype='Float64')


def test_nullable_boolean_not_promoted_to_numeric_expression(tmp_path):
    path = tmp_path / 'bool.csv'
    path.write_text('cell,G,R\na,True,False\n')
    data = read_anndata(path, mode='csv', dtype='boolean')
    data.obs['cell_type'] = ['A']
    with pytest.raises(ValueError, match='real values'):
        run_cell_mesh(data, *priors(), min_cells=1)
