"""B02: QC text must not silently hide otherwise valid events."""
import io

import numpy as np
import pandas as pd
import pytest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from cellmesh import (plot_significant_event_counts, plot_communication_network,
                      plot_event_dotplot, plot_sample_event_scores)

PLOTS = [plot_significant_event_counts, plot_communication_network,
         plot_event_dotplot, plot_sample_event_scores]


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


def events(flags):
    return pd.DataFrame(dict(sender=['A'] * 3, receiver=['B'] * 3,
                             hmdb_id=['HMDB1'] * 3, sensor_gene=['R'] * 3,
                             sample=['s1', 's2', 's3'], cell_mesh_score=[.6, .5, .4],
                             perm_pvalue=[.01] * 3, fdr_sensor_type=[.02] * 3,
                             passes_min_cells=flags))


def plot_events(plot, data, qc_only=True):
    extra = dict(sender='A', receiver='B', hmdb_id='HMDB1', sensor_gene='R') if plot is plot_sample_event_scores else {}
    return plot(data, qc_only=qc_only, **extra)


@pytest.mark.parametrize('plot', PLOTS)
@pytest.mark.parametrize('flags', [
    [True, False, None], [np.bool_(True), np.bool_(False), np.nan],
    pd.array([True, False, pd.NA], dtype='boolean'), ['True', 'False', None],
    [' tRuE ', ' FALSE ', pd.NA], [1, 0, np.nan], [1., 0., None], ['1', '0', None],
])
@pytest.mark.parametrize('qc_only', [True, False])
def test_qc_representations_have_identical_filtering(plot, flags, qc_only):
    frame = events(flags)
    original = frame.copy(deep=True)
    actual = plot_events(plot, frame, qc_only)
    expected = plot_events(plot, events(pd.array([True, False, pd.NA], dtype='boolean')), qc_only)
    assert actual['qc'] == expected['qc']
    key = ('sample_scores' if plot is plot_sample_event_scores else
           'plot_events' if plot is plot_event_dotplot else 'filtered_events')
    pd.testing.assert_frame_equal(actual[key], expected[key])
    pd.testing.assert_frame_equal(frame, original)
    if not qc_only:
        assert actual['qc_failed_events']['passes_min_cells'].isna().sum() == 1


@pytest.mark.parametrize('plot', PLOTS)
@pytest.mark.parametrize('qc_only', [True, False])
@pytest.mark.parametrize('value', ['bad', '', 'NA', 2, -1, .5, np.inf, 1+0j,
                                  [True], pd.Timestamp('2026-09-17')])
def test_bad_qc_fails_explicitly_even_without_filter(plot, qc_only, value):
    frame = events([True, False, None])
    frame['passes_min_cells'] = frame['passes_min_cells'].astype(object)
    frame.at[1, 'passes_min_cells'] = value
    with pytest.raises(ValueError, match='passes_min_cells.*position 1'):
        plot_events(plot, frame, qc_only)


@pytest.mark.parametrize('plot', PLOTS)
def test_qc_survives_csv_read_as_text(plot):
    original = events(pd.array([True, False, pd.NA], dtype='boolean'))
    reloaded = pd.read_csv(io.StringIO(original.to_csv(index=False)), dtype=str)
    actual = plot_events(plot, reloaded)
    expected = plot_events(plot, original)
    assert actual['qc'] == expected['qc']
    if plot is plot_significant_event_counts:
        pd.testing.assert_frame_equal(actual['counts'], expected['counts'])


@pytest.mark.parametrize('plot', PLOTS)
def test_absent_qc_retains_all_events(plot):
    plotted = plot_events(plot, events([True, False, None]).drop(columns='passes_min_cells'))
    assert plotted['qc']['events_after_qc'] == 3
    assert not plotted['qc']['qc_applied']
