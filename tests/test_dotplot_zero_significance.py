"""A07: zero is a distinct display cap, never a fabricated positive probability."""
import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt

from cellmesh import plot_event_dotplot
from test_dotplot_significance import events, _plot, _legend, _assert_visible_points


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


CASES = [
    [0, 0.01, 0.1, 0.1],
    [0, 1, 1, 1],
    [0, 0, 0, 0],
    [0, 0, np.nan, np.nan],
    [0, np.nextafter(0., 1.), 1e-300, np.nan],
    [np.nextafter(0., 1.), 1e-300, 0.1, 1],
    [0, 0.2, 0.2, np.nan],
    [0.1, 0.1, 0.1, 0.1],
    [np.nan, np.nan, np.nan, np.nan],
]


@pytest.mark.parametrize('supplied_ax', [False, True])
@pytest.mark.parametrize('statistic', ['fdr_sensor_type', 'perm_pvalue'])
@pytest.mark.parametrize('values', CASES)
def test_zero_positive_missing_sizes_and_truthful_legend(events, supplied_ax, statistic, values):
    if statistic == 'perm_pvalue':
        events = events.drop(columns='fdr_sensor_type')
    events[statistic] = values
    result = _plot(events, supplied_ax)
    sizes = _assert_visible_points(result, statistic)
    title, labels = _legend(result)
    has_zero = 0 in values
    positive = sorted(set(v for v in values if pd.notna(v) and v > 0))
    if has_zero:
        assert '0 (display cap)' in labels
        assert result['size_encoding']['zero_display_cap'] == 260
    else:
        assert not any('display cap' in label for label in labels)
        assert result['size_encoding']['zero_display_cap'] is None
    numeric_labels = [float(label) for label in labels if label not in ['0 (display cap)', 'Unavailable']]
    assert all(value in positive for value in numeric_labels)
    assert all(value > 0 for value in numeric_labels)
    assert ('Unavailable' in labels) == any(pd.isna(v) for v in values)
    for record, probability in enumerate(values):
        if pd.isna(probability):
            assert sizes[record] == 140
        elif probability == 0:
            assert sizes[record] == 260
        elif has_zero:
            assert 20 <= sizes[record] <= 212 < 260
    for i, left in enumerate(values):
        for j, right in enumerate(values):
            if 0 < left < right:
                assert sizes[i] > sizes[j]
    actual = result['plot_events']
    assert result['dot_sizes'].index.equals(actual.index)
    assert result['dot_sizes'].tolist() == [sizes[record] for record in actual.record]
    assert result['size_encoding']['statistic'] == statistic
    if not positive and has_zero:
        assert labels == (['0 (display cap)', 'Unavailable'] if any(pd.isna(v) for v in values)
                          else ['0 (display cap)'])
        assert result['size_encoding']['positive_area_range'] is None
    if has_zero or positive:
        assert title == ('FDR' if statistic == 'fdr_sensor_type' else 'p-value')


@pytest.mark.parametrize('supplied_ax', [False, True])
@pytest.mark.parametrize('statistic', ['fdr_sensor_type', 'perm_pvalue'])
@pytest.mark.parametrize('values', [[0, 0.01, 0.1, np.nan], [0, 0, 0, 0], [0.01, 0.1, 1, np.nan]])
def test_fixed_area_explicitly_disables_size_encoding(events, supplied_ax, statistic, values):
    if statistic == 'perm_pvalue':
        events = events.drop(columns='fdr_sensor_type')
    events[statistic] = pd.Series(values, dtype='Float64')
    result = _plot(events, supplied_ax, min_dot_size=64, max_dot_size=64)
    sizes = _assert_visible_points(result, statistic, fixed_size=64)
    assert set(sizes.values()) == {64}
    title, labels = _legend(result)
    assert 'Fixed size (no size encoding)' in title
    assert not any('display cap' in label for label in labels)
    assert ('0' in labels) == (0 in values)
    assert result['size_encoding']['fixed_size'] is True


@pytest.mark.parametrize('statistic', ['fdr_sensor_type', 'perm_pvalue'])
def test_zero_outside_selected_context_does_not_compress_positive_areas(events, statistic):
    if statistic == 'perm_pvalue':
        events = events.drop(columns='fdr_sensor_type')
    events[statistic] = [0.01, 0.1, 0, 0]
    plotted = _plot(events, sender_labels=['A'])
    assert set(plotted['dot_sizes']) == {20, 260}
    assert plotted['size_encoding']['zero_display_cap'] is None
    assert '0 (display cap)' not in _legend(plotted)[1]


@pytest.mark.parametrize('statistic', ['fdr_sensor_type', 'perm_pvalue'])
def test_positive_only_mapping_is_unchanged(events, statistic):
    if statistic == 'perm_pvalue':
        events = events.drop(columns='fdr_sensor_type')
    events[statistic] = [1e-8, 1e-4, 0.01, 1]
    plotted = _plot(events, min_dot_size=40, max_dot_size=180)
    probabilities = plotted['plot_events'][statistic].to_numpy()
    transformed = -np.log10(probabilities)
    expected = 40 + (transformed - transformed.min()) / np.ptp(transformed) * 140
    np.testing.assert_array_equal(plotted['dot_sizes'], expected)


def test_missing_fdr_does_not_fall_back_to_zero_pvalue(events):
    events['fdr_sensor_type'] = np.nan
    events['perm_pvalue'] = 0.0
    result = _plot(events)
    assert _legend(result) == ('Significance', ['Unavailable'])
    assert result['size_encoding']['zero_display_cap'] is None
    assert result['size_encoding']['statistic'] == 'fdr_sensor_type'


def test_custom_probability_column_legend_source(events):
    events = events.drop(columns='fdr_sensor_type').rename(columns={'perm_pvalue': 'external_p'})
    events['external_p'] = [0, 0.01, 0.1, 1]
    result = _plot(events, pvalue_col='external_p')
    assert _legend(result)[0] == 'p-value (external_p)'
    assert result['size_encoding']['statistic'] == 'external_p'


def test_close_positive_values_do_not_get_identical_legend_labels(events):
    values = [0, 0.10000000001, 0.10000000002, 0.10000000003]
    events['fdr_sensor_type'] = values
    result = _plot(events)
    _, labels = _legend(result)
    assert len(set(labels)) == len(labels) == 4
    assert set(float(label) for label in labels[:-1]) == set(values[1:])
