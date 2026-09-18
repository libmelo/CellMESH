"""B06: valid float16 inputs must survive plotting's multi-column sorts."""
import numpy as np
import pandas as pd
import pytest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from cellmesh import plot_communication_network, plot_event_dotplot
from cellmesh.plotting import _plot_numeric_series


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


def events():
    # Include duplicate identities to check selection, not merely successful rendering.
    return pd.DataFrame(dict(sender=['A'] * 5, receiver=['B', 'B', 'C', 'B', 'C'],
                             hmdb_id=['HMDB1', 'HMDB1', 'HMDB2', 'HMDB3', 'HMDB4'],
                             sensor_gene=['R'] * 5, cell_mesh_score=[.2, .7, .4, .6, np.nan],
                             perm_pvalue=[.04, .01, 0., np.nan, .02],
                             fdr_sensor_type=[.04, .02, 0., np.nan, .03],
                             passes_min_cells=[True] * 5))


@pytest.mark.parametrize('plot', [plot_event_dotplot, plot_communication_network])
@pytest.mark.parametrize('columns', [('cell_mesh_score',), ('perm_pvalue',), ('fdr_sensor_type',),
                                    ('cell_mesh_score', 'perm_pvalue', 'fdr_sensor_type')])
@pytest.mark.parametrize('filtered', [False, True])
def test_float16_matches_same_values_promoted_to_float64(plot, columns, filtered):
    data = events()
    data[list(columns)] = data[list(columns)].astype(np.float16)
    before = data.copy(deep=True)
    reference = data.copy(deep=True)
    reference[list(columns)] = reference[list(columns)].astype(np.float64)
    options = dict(max_fdr=.05 if filtered else None)
    actual, expected = plot(data, **options), plot(reference, **options)
    actual['fig'].canvas.draw()
    expected['fig'].canvas.draw()
    key = 'plot_events' if plot is plot_event_dotplot else 'filtered_events'
    pd.testing.assert_frame_equal(actual[key], expected[key], check_exact=True)
    pd.testing.assert_frame_equal(actual['missing_score_events'], expected['missing_score_events'])
    if plot is plot_event_dotplot:
        pd.testing.assert_series_equal(actual['dot_sizes'], expected['dot_sizes'], check_exact=True)
        assert actual['selected_event_keys'] == expected['selected_event_keys']
    else:
        pd.testing.assert_frame_equal(actual['edge_table'], expected['edge_table'])
    pd.testing.assert_frame_equal(data, before)


@pytest.mark.parametrize('value', [np.inf, -np.inf, -.1, 1.5])
def test_float16_invalid_values_still_rejected(value):
    data = events()
    data['cell_mesh_score'] = data['cell_mesh_score'].astype(np.float16)
    data.loc[0, 'cell_mesh_score'] = np.float16(value)
    with pytest.raises(ValueError, match='cell_mesh_score'):
        plot_event_dotplot(data)


@pytest.mark.parametrize('dtype', ['float32', 'float64', 'Float32', 'Float64', 'Int64'])
def test_other_valid_dtypes_keep_existing_contract(dtype):
    data = pd.Series([0, 1, None], dtype=dtype, name='value')
    pd.testing.assert_series_equal(_plot_numeric_series(data, 'value', upper=1.), data)
