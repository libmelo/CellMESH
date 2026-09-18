"""A08: measured sidebar boxes fit and leave caller-owned subplots alone."""
from io import BytesIO
import warnings

import numpy as np
import pandas as pd
import pytest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

from cellmesh import plot_event_dotplot


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


def events(statistic='fdr_sensor_type', state='positive', long=False, senders=1):
    values = {
        'positive': [0.10000000001, 0.10000000002, 0.10000000003, 0.2, 1],
        'zero': [0, 0.01, 0.1, 0.5, 1],
        'partial': [0, 0.01, 0.1, 0.5, np.nan],
        'missing': [np.nan] * 5,
        'fixed': [0, 0.01, 0.1, 0.5, np.nan],
    }[state]
    out = pd.DataFrame(dict(sender=['SenderType']*5,
                            receiver=['Receiver type with a longer label' if long else 'B']*5,
                            hmdb_id=[f'HMDB{i}' for i in range(5)],sensor_gene=['RECEPTOR']*5,
                            metabolite=[('A metabolite with a considerably longer name ' if long else 'M')+str(i)
                                        for i in range(5)],
                            cell_mesh_score=np.linspace(.1,.9,5)))
    out[statistic] = values
    if senders > 1:
        out = pd.concat([out.assign(sender=f'SenderType_{i}') for i in range(senders)], ignore_index=True)
    return out


def assert_layout(result, owned=None):
    fig = result['fig']
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [axis.get_tightbbox(renderer) for axis in result['axes']]
    boxes += [result['legend'].get_window_extent(renderer),
              result['colorbar'].ax.get_tightbbox(renderer),
              result['title_artist'].get_window_extent(renderer)]
    region = fig.bbox if owned is None else owned.transformed(fig.transFigure)
    for box in boxes:
        assert box.width > 0 and box.height > 0
        assert box.x0 >= region.x0 - 0.5 and box.y0 >= region.y0 - 0.5
        assert box.x1 <= region.x1 + 0.5 and box.y1 <= region.y1 + 0.5
    for i, first in enumerate(boxes):
        for second in boxes[i+1:]:
            assert not first.overlaps(second), (first.bounds, second.bounds)
    # Individual legend handles and text must not escape the legend allocation.
    legend_box = boxes[-3]
    for artist in [result['legend'].get_title(), *result['legend'].get_texts(),
                   *result['legend'].legend_handles]:
        box = artist.get_window_extent(renderer)
        assert box.x0 >= legend_box.x0 - 1 and box.x1 <= legend_box.x1 + 1
        assert box.y0 >= legend_box.y0 - 1 and box.y1 <= legend_box.y1 + 1


@pytest.mark.parametrize('statistic', ['fdr_sensor_type', 'perm_pvalue'])
@pytest.mark.parametrize('state', ['positive', 'zero', 'partial', 'missing', 'fixed'])
@pytest.mark.parametrize('mode', ['auto_single', 'auto_multiple', 'external'])
@pytest.mark.parametrize('long', [False, True])
def test_sidebar_and_labels_fit_all_branches(statistic, state, mode, long):
    data = events(statistic, state, long, 3 if mode == 'auto_multiple' else 1)
    before = data.copy(deep=True)
    kwargs = dict(top_n=None)
    if state == 'fixed':
        kwargs.update(min_dot_size=64, max_dot_size=64)
    owned = None
    if mode == 'external':
        _, ax = plt.subplots(figsize=(16, 8))
        owned = ax.get_position().frozen()
        kwargs['ax'] = ax
    result = plot_event_dotplot(data, **kwargs)
    assert_layout(result, owned)
    pd.testing.assert_frame_equal(data, before)
    assert len(result['plot_events']) == len(data)
    assert result['layout']['external_ax'] == (mode == 'external')


@pytest.mark.parametrize('engine', [None, 'tight', 'constrained'])
def test_external_subplot_preserves_other_axes_and_figure_title(engine):
    fig, axes = plt.subplots(1,2,figsize=(18,8),layout=engine)
    axes[1].plot([0,1], [1,0]); axes[1].set_title('Other panel')
    original_title = fig.suptitle('User-owned figure title')
    fig.canvas.draw()
    if engine is not None:
        fig.set_layout_engine('none')
    owned = axes[0].get_position().frozen()
    other_position = axes[1].get_position().frozen()
    other_children = list(axes[1].get_children())
    original_size = fig.get_size_inches().copy()
    result = plot_event_dotplot(events(state='partial'), ax=axes[0])
    assert_layout(result, owned)
    np.testing.assert_array_equal(axes[1].get_position().bounds, other_position.bounds)
    np.testing.assert_array_equal(fig.get_size_inches(), original_size)
    assert fig._suptitle is original_title and original_title.get_text() == 'User-owned figure title'
    assert axes[1].get_children() == other_children


@pytest.mark.parametrize('engine', ['tight', 'constrained'])
def test_active_external_layout_engine_errors_before_mutation(engine):
    fig, ax = plt.subplots(figsize=(12,7),layout=engine)
    fig.canvas.draw()
    position = ax.get_position().frozen(); children = list(ax.get_children()); original_axes = list(fig.axes)
    with pytest.raises(ValueError, match='finalized figure layout.*set_layout_engine'):
        plot_event_dotplot(events(), ax=ax)
    np.testing.assert_array_equal(ax.get_position().bounds, position.bounds)
    assert list(ax.get_children()) == children
    assert fig.axes == original_axes


def test_small_external_ax_gives_required_and_available_space():
    _, ax = plt.subplots(figsize=(3,2))
    with pytest.raises(ValueError, match=r'needs at least .* inches.*available.*Enlarge'):
        plot_event_dotplot(events(state='fixed',long=True), ax=ax, min_dot_size=64,max_dot_size=64)


@pytest.mark.parametrize('file_format', ['png','svg','pdf'])
@pytest.mark.parametrize('dpi', [100,200])
def test_export_at_multiple_dpi_keeps_layout(file_format,dpi):
    result = plot_event_dotplot(events('perm_pvalue','partial',True,2),top_n=None)
    result['fig'].set_dpi(dpi)
    assert_layout(result)
    output = BytesIO()
    with warnings.catch_warnings():
        warnings.simplefilter('error',RuntimeWarning)
        result['fig'].savefig(output,format=file_format,dpi=dpi)
    assert output.tell() > 1000
    assert_layout(result)


def test_large_markers_and_long_custom_titles_are_measured():
    data = events('external_adjusted_probability_column_name','partial')
    result = plot_event_dotplot(data, fdr_col='external_adjusted_probability_column_name',
                               min_dot_size=100,max_dot_size=2500,
                               score_label='A deliberately much longer colorbar title')
    assert_layout(result)


def test_auto_layout_ignores_global_automatic_layout_defaults():
    with matplotlib.rc_context({'figure.autolayout':True}):
        result=plot_event_dotplot(events(state='zero'))
        assert_layout(result)
