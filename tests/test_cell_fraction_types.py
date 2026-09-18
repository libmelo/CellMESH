"""B03: validate external fraction types before any float coercion."""
from datetime import date, timedelta
import warnings

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import compute_metabolite_availability
from cellmesh.score import compute_sensor_scores


@pytest.fixture
def inputs():
    data = AnnData(np.array([[1., 2.], [3., 4.]]),
                   obs=pd.DataFrame({'cell_type': ['A', 'B']}, index=['a', 'b']),
                   var=pd.DataFrame(index=['G', 'R']))
    enzyme = pd.DataFrame(dict(metabolite=['M'], hmdb_id=['HMDB1'], gene=['G'],
                              role=['production'], reaction=['p']))
    sensor = pd.DataFrame(dict(metabolite=['M'], hmdb_id=['HMDB1'], sensor_gene=['R'],
                              sensor_type=['Transporter']))
    return data, enzyme, sensor


def score(inputs, api, fractions):
    data, enzyme, sensor = inputs
    fn, prior = (compute_metabolite_availability, enzyme) if api == 'sender' else (compute_sensor_scores, sensor)
    return fn(data, prior, min_cells=1, cell_fractions=fractions)


@pytest.mark.parametrize('api', ['sender', 'receiver'])
@pytest.mark.parametrize('value', [True, False, np.bool_(True), .5+2j, np.complex128(.5),
                                   pd.Timestamp('2026-09-17'), np.datetime64('2026-09-17'),
                                   date(2026, 9, 17), timedelta(days=1), pd.Timedelta('1 day'),
                                   'bad', 'True', '', 'NA', np.nan, pd.NA, None,
                                   np.inf, -np.inf, 0, -.1, 1.1])
def test_invalid_fraction_scalars_rejected_without_lossy_warnings(inputs, api, value):
    fractions = pd.Series([value, .5], index=['A', 'B'], dtype=object)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        with pytest.raises((TypeError, ValueError), match='cell_fractions'):
            score(inputs, api, fractions)


@pytest.mark.parametrize('api', ['sender', 'receiver'])
@pytest.mark.parametrize('dtype', [complex, object])
def test_native_and_object_complex_series_are_rejected(inputs, api, dtype):
    fractions = pd.Series([.5+20j, .5+30j], index=['A', 'B'], dtype=dtype)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        with pytest.raises(TypeError, match='cell_fractions'):
            score(inputs, api, fractions)


@pytest.mark.parametrize('api', ['sender', 'receiver'])
def test_single_group_boolean_is_not_a_fraction(inputs, api):
    data, enzyme, sensor = inputs
    inputs = (data[:1].copy(), enzyme, sensor)
    with pytest.raises(TypeError, match='cell_fractions'):
        score(inputs, api, pd.Series([True], index=['A']))


@pytest.mark.parametrize('api', ['sender', 'receiver'])
@pytest.mark.parametrize('form', ['float', 'text', 'nullable', 'dict', 'extra_population'])
def test_legitimate_custom_fractions_preserve_values_and_caller(inputs, api, form):
    fractions = pd.Series([.75, .25], index=[' B ', ' A '])
    if form == 'text':
        fractions = pd.Series(['7.5e-1', '0.25'], index=fractions.index)
    elif form == 'nullable':
        fractions = fractions.astype('Float64')
    elif form == 'dict':
        fractions = fractions.to_dict()
    elif form == 'extra_population':
        fractions = pd.Series([.3, .1, .6], index=[' B ', ' A ', 'Other'])
    before = fractions.copy()
    expected_fractions = pd.Series([.1, .3] if form == 'extra_population' else [.25, .75], index=['A', 'B'])
    actual = score(inputs, api, fractions)
    expected = score(inputs, api, expected_fractions)
    if api == 'sender':
        for key in ['availability', 'P', 'C', 'E', 'metadata']:
            pd.testing.assert_frame_equal(actual[key], expected[key], check_exact=True)
        pd.testing.assert_series_equal(actual['cell_fractions'], expected['cell_fractions'], check_exact=True)
    else:
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
    if isinstance(fractions, pd.Series):
        pd.testing.assert_series_equal(fractions, before)
    else:
        assert fractions == before


@pytest.mark.parametrize('api', ['sender', 'receiver'])
@pytest.mark.parametrize('defect', ['missing', 'duplicate', 'blank', 'sum'])
def test_fraction_axis_and_sum_contracts_remain_enforced(inputs, api, defect):
    fractions = pd.Series([.5, .5], index=['A', 'B'])
    if defect == 'missing': fractions = fractions.iloc[:1]
    elif defect == 'duplicate': fractions.index = ['A', ' A ']
    elif defect == 'blank': fractions.index = ['A', ' ']
    else: fractions[:] = .8
    with pytest.raises(ValueError, match='cell_fractions'):
        score(inputs, api, fractions)
