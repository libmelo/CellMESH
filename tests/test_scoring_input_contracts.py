"""A05: original observations, rather than supplied caches, define coverage."""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import run_cell_mesh
from cellmesh.preprocess import _all_celltype_counts, _build_celltype_pseudobulk, _compute_celltype_expr_frac
from cellmesh.score import compute_metabolite_availability, compute_sensor_scores


@pytest.fixture
def inputs():
    adata = AnnData(
        np.array([[1., 2., 0.], [2., 2., 4.], [3., 8., 0.], [4., 8., 0.]]),
        obs=pd.DataFrame({'cell_type': pd.Categorical(['A', 'A', 'B', 'B'], categories=['A', 'B', 'UNUSED']),
                          'sample': ['s1', 's2', 's1', 's2']}, index=['c1', 'c2', 'c3', 'c4']),
        var=pd.DataFrame(index=['G', 'S', 'OTHER']),
    )
    enzyme = pd.DataFrame({'metabolite': ['M'], 'hmdb_id': ['HMDB0000001'],
                           'gene': ['G'], 'role': ['production'], 'reaction': ['r1']})
    sensor = pd.DataFrame({'metabolite': ['M'], 'hmdb_id': ['HMDB0000001'],
                          'sensor_gene': ['S'], 'sensor_type': ['Transporter']})
    return adata, enzyme, sensor


def summaries(adata):
    return dict(pseudobulk=_build_celltype_pseudobulk(adata),
                expr_frac=_compute_celltype_expr_frac(adata),
                cell_counts=_all_celltype_counts(adata))


@pytest.mark.parametrize('reference', ['mean', 'median'])
@pytest.mark.parametrize('gate', [None, 0.5])
@pytest.mark.parametrize('supplied', [('pseudobulk',), ('expr_frac',), ('cell_counts',),
                                      ('pseudobulk', 'expr_frac', 'cell_counts')])
def test_valid_reordered_normalized_caches_match_uncached_without_mutation(inputs, reference, gate, supplied):
    adata, _, sensor = inputs
    expected = compute_sensor_scores(adata, sensor, receiver_reference=reference, min_expr_frac=gate)
    cached = {key: value for key, value in summaries(adata).items() if key in supplied}
    for key, value in cached.items():
        value = value.iloc[::-1].copy()
        value.index = [' ' + str(x) + ' ' for x in value.index]
        if isinstance(value, pd.DataFrame):
            value = value.loc[:, ['S']]
            value.columns = [' S ']
        cached[key] = value
    before = {key: value.copy(deep=True) for key, value in cached.items()}
    result = compute_sensor_scores(adata, sensor, receiver_reference=reference, min_expr_frac=gate, **cached)
    pd.testing.assert_frame_equal(result, expected, check_exact=True)
    for key, value in cached.items():
        (pd.testing.assert_frame_equal if isinstance(value, pd.DataFrame) else pd.testing.assert_series_equal)(value, before[key])
    assert result.set_index('receiver').loc['A', 'sensor_score'] == pytest.approx(2 / 7)


@pytest.mark.parametrize('name', ['pseudobulk', 'expr_frac', 'cell_counts'])
@pytest.mark.parametrize('defect', ['missing_group', 'extra_group', 'duplicate', 'collision', 'blank', 'missing_label'])
def test_cache_cell_axis_must_match_actual_observations(inputs, name, defect):
    adata, _, sensor = inputs
    cached = summaries(adata)[name]
    if defect == 'missing_group':
        cached = cached.loc[['A']]
    elif defect == 'extra_group':
        cached = pd.concat([cached, cached.loc[['B']].rename(index={'B': 'UNUSED'})])
    elif defect == 'duplicate':
        cached = pd.concat([cached, cached.loc[['B']]])
    elif defect == 'collision':
        cached = pd.concat([cached, cached.loc[['B']].rename(index={'B': ' B '})])
    else:
        cached.index = ['A', ' ' if defect == 'blank' else None]
    with pytest.raises(ValueError, match=name):
        compute_sensor_scores(adata, sensor, **{name: cached})


@pytest.mark.parametrize('name', ['pseudobulk', 'expr_frac'])
@pytest.mark.parametrize('defect', ['missing_sensor', 'unknown_gene', 'duplicate', 'collision', 'blank', 'missing_label'])
def test_cache_gene_axis_cannot_hide_measured_sensor_or_add_unmeasured_gene(inputs, name, defect):
    adata, _, sensor = inputs
    cached = summaries(adata)[name]
    if defect == 'missing_sensor':
        cached = cached.drop(columns='S')
    elif defect == 'unknown_gene':
        cached['UNMEASURED'] = 1.
    elif defect in ('duplicate', 'collision'):
        extra = cached[['S']].rename(columns={'S': 'S' if defect == 'duplicate' else ' S '})
        cached = pd.concat([cached, extra], axis=1)
    else:
        cached.columns = ['G', ' ' if defect == 'blank' else None, 'OTHER']
    with pytest.raises(ValueError, match=name):
        compute_sensor_scores(adata, sensor, **{name: cached})


@pytest.mark.parametrize('name', ['pseudobulk', 'expr_frac'])
@pytest.mark.parametrize('value', [-1., np.nan, np.inf, -np.inf, 'bad', True, 1+2j, pd.NA])
@pytest.mark.parametrize('gate', [None, 0.1])
def test_invalid_summary_values_are_rejected_even_without_gate(inputs, name, value, gate):
    adata, _, sensor = inputs
    cached = summaries(adata)[name].astype(object)
    cached.loc['A', 'S'] = value
    with pytest.raises((TypeError, ValueError), match=name):
        compute_sensor_scores(adata, sensor, min_expr_frac=gate, **{name: cached})


def test_fraction_upper_bound_and_legitimate_zero_mean(inputs):
    adata, _, sensor = inputs
    adata.X[2:, 1] = 0
    cached = summaries(adata)
    result = compute_sensor_scores(adata, sensor, min_cells=100, **cached).set_index('receiver')
    assert result.loc['B', 'sensor_score'] == 0
    assert result.loc['B', 'sensor_expr_frac'] == 0
    assert not result.loc['B', 'receiver_passes_min_cells']
    cached['expr_frac'].loc['A', 'S'] = 1.01
    with pytest.raises(ValueError, match='expr_frac'):
        compute_sensor_scores(adata, sensor, min_expr_frac=None, **cached)


@pytest.mark.parametrize('value', [-1.5, 1.9, 0, -1, np.nan, np.inf, 1, 3, True, 2+0j, 'bad', pd.NA])
def test_counts_are_validated_before_integer_conversion_and_match_obs(inputs, value):
    adata, _, sensor = inputs
    counts = pd.Series({'A': value, 'B': 2}, dtype=object)
    with pytest.raises((TypeError, ValueError), match='cell_counts'):
        compute_sensor_scores(adata, sensor, cell_counts=counts)


def test_integral_float_counts_and_numeric_text_means_are_supported(inputs):
    adata, _, sensor = inputs
    cached = summaries(adata)
    cached['cell_counts'] = cached['cell_counts'].astype(float)
    cached['pseudobulk'] = cached['pseudobulk'].astype(str)
    expected = compute_sensor_scores(adata, sensor)
    pd.testing.assert_frame_equal(compute_sensor_scores(adata, sensor, **cached), expected, check_exact=True)


@pytest.mark.parametrize('name', ['pseudobulk', 'expr_frac', 'cell_counts'])
def test_summary_container_types_are_explicit(inputs, name):
    adata, _, sensor = inputs
    with pytest.raises(TypeError, match=name):
        compute_sensor_scores(adata, sensor, **{name: {'A': 1}})


def test_current_subset_not_unused_categories_defines_coverage(inputs):
    adata, _, sensor = inputs
    subset = adata[:2].copy()
    result = compute_sensor_scores(subset, sensor, **summaries(subset))
    assert result.receiver.tolist() == ['A']
    assert result.sensor_score.tolist() == [0.5]


def test_unrelated_nonfinite_expression_columns_do_not_invalidate_receiver_cache(inputs):
    adata, _, sensor = inputs
    cached = summaries(adata)
    for name in ['pseudobulk', 'expr_frac']:
        cached[name]['OTHER'] = np.inf
    pd.testing.assert_frame_equal(compute_sensor_scores(adata, sensor, **cached),
                                  compute_sensor_scores(adata, sensor), check_exact=True)


def test_cache_does_not_recompute_expression_summaries(inputs, monkeypatch):
    adata, _, sensor = inputs
    cached = summaries(adata)
    def unexpected(*args, **kwargs):
        raise AssertionError('Provided summaries must not be recomputed')
    monkeypatch.setattr('cellmesh.score._build_celltype_pseudobulk', unexpected)
    monkeypatch.setattr('cellmesh.score._compute_celltype_expr_frac', unexpected)
    result = compute_sensor_scores(adata, sensor, **cached)
    assert len(result) == 2


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
@pytest.mark.parametrize('n_perms', [0, 3])
@pytest.mark.parametrize('name', ['allow_self', 'store_null_scores'])
@pytest.mark.parametrize('value', ['False', 'True', 0, 1, None, np.nan])
def test_main_boolean_parameters_reject_truthiness(inputs, mode, n_perms, name, value):
    adata, enzyme, sensor = inputs
    with pytest.raises(TypeError, match=name):
        run_cell_mesh(adata, enzyme, sensor, sample_key='sample', sample_mode=mode,
                      n_perms=n_perms, **{name: value})


@pytest.mark.parametrize('value', ['False', 0, 1, None, np.nan])
def test_availability_boolean_parameter_rejects_truthiness(inputs, value):
    adata, enzyme, _ = inputs
    with pytest.raises(TypeError, match='return_intermediates'):
        compute_metabolite_availability(adata, enzyme, return_intermediates=value)


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
@pytest.mark.parametrize('n_perms', [0, 3])
@pytest.mark.parametrize('value', [np.nan, np.inf, -1, 1.5, 1., True, np.bool_(False), '42', None, [42], np.random.default_rng(42)])
def test_seed_validation_does_not_depend_on_permutation_execution(inputs, mode, n_perms, value):
    adata, enzyme, sensor = inputs
    with pytest.raises((TypeError, ValueError), match='random_state'):
        run_cell_mesh(adata, enzyme, sensor, sample_key='sample', sample_mode=mode,
                      n_perms=n_perms, random_state=value)


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
@pytest.mark.parametrize('n_perms', [0, 3])
def test_invalid_seed_rejected_before_empty_event_shortcut(inputs, mode, n_perms):
    adata, enzyme, sensor = inputs
    adata.obs['cell_type'] = 'A'
    with pytest.raises(TypeError, match='random_state'):
        run_cell_mesh(adata, enzyme, sensor, sample_key='sample', sample_mode=mode,
                      n_perms=n_perms, allow_self=False, random_state=np.nan)


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_numpy_parameters_match_python_and_parallel_seeded_results(inputs, mode):
    adata, enzyme, sensor = inputs
    options = dict(sample_key='sample', sample_mode=mode, n_perms=7, min_cells=1)
    expected = run_cell_mesh(adata, enzyme, sensor, random_state=42, allow_self=False,
                             store_null_scores=True, n_jobs=1, **options)
    actual = run_cell_mesh(adata, enzyme, sensor, random_state=np.int64(42),
                           allow_self=np.bool_(False), store_null_scores=np.bool_(True),
                           n_jobs=2, **options)
    pd.testing.assert_frame_equal(actual.events, expected.events, check_exact=True)
    assert actual.parameters == expected.parameters | {'n_jobs': 2}
    assert not actual.events.sender.eq(actual.events.receiver).any()
    if mode == 'sample_aware':
        pd.testing.assert_frame_equal(actual.events.attrs['sample_aware_null_scores'], expected.events.attrs['sample_aware_null_scores'])
    plain = compute_metabolite_availability(adata, enzyme, return_intermediates=False)
    numpy_bool = compute_metabolite_availability(adata, enzyme, return_intermediates=np.bool_(False))
    assert plain.keys() == numpy_bool.keys()


@pytest.mark.parametrize('name', ['pseudobulk', 'expr_frac', 'cell_counts'])
@pytest.mark.parametrize('kind', ['datetime', 'timedelta', 'categorical_bool'])
def test_nonreal_summary_dtypes_are_not_silently_encoded_as_numbers(inputs, name, kind):
    adata, _, sensor = inputs
    cached = summaries(adata)[name]
    values = {'datetime': pd.date_range('2026-01-01', periods=2),
              'timedelta': pd.to_timedelta([1, 2], unit='ns'),
              'categorical_bool': pd.Categorical([True, False])}[kind]
    if isinstance(cached, pd.DataFrame):
        cached['S'] = values
    else:
        cached = pd.Series(values, index=cached.index)
    with pytest.raises(TypeError, match=name):
        compute_sensor_scores(adata, sensor, **{name: cached})


@pytest.mark.parametrize('seed', [0, np.uint64(42), 2**128])
def test_nonnegative_integer_seed_contract_preserves_large_seeds(inputs, seed):
    adata, enzyme, sensor = inputs
    result = run_cell_mesh(adata, enzyme, sensor, n_perms=2, random_state=seed)
    assert result.parameters['random_state'] == int(seed)
    assert type(result.parameters['random_state']) is int


def test_no_measured_sensor_does_not_require_unmeasured_cache_column(inputs):
    adata, _, sensor = inputs
    sensor['sensor_gene'] = 'UNMEASURED'
    result = compute_sensor_scores(adata, sensor, **summaries(adata))
    assert result.empty
    assert 'sensor_score' in result.columns


@pytest.mark.parametrize('layer', [None, 'expression'])
def test_sample_specific_summary_coverage_and_missing_endpoint(inputs, layer):
    adata, enzyme, sensor = inputs
    # s1 has A and B; s2 has A only. Neither an unused categorical B in s2
    # nor its presence elsewhere may force an invented zero-expression row.
    adata = adata[:3].copy()
    if layer:
        adata.layers[layer] = adata.X.copy()
        adata.X[:] = np.nan
    result = run_cell_mesh(adata, enzyme, sensor, sample_key='sample',
                           sample_mode='sample_aware', layer=layer, min_cells=1, n_perms=3)
    rows = result.sample_events
    missing = rows['sample'].eq('s2') & (rows.sender.eq('B') | rows.receiver.eq('B'))
    assert missing.any()
    assert rows.loc[missing, 'cell_mesh_score'].isna().all()
