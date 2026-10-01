"""Public-entry and Figure-preservation contracts, using automated fixtures only."""
import json
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import pytest

import cellmesh
from cellmesh import plotting
from cellmesh import _communication_dotmap as matrix
from cellmesh._report_dotmap import communication_dotmap_html
from cellmesh.reporting import render_html_report
from test_communication_features import path_frame
from test_reporting_and_persistence import make_result


def payload(text, prefix):
    return json.loads(re.search(r'id="' + prefix + r'-dotmap-data">(.*?)</script>', text, re.S)[1])


def test_public_names_export_one_function_and_identity_remains_stable():
    from cellmesh.selection import event_id
    from cellmesh._report_browser import event_id as report_id
    assert cellmesh.plot_communication_dotmap is plotting.plot_communication_dotmap is matrix.plot_communication_dotmap
    assert event_id is report_id
    assert event_id(dict(sender=' A ', receiver='B', hmdb_id='hmdb0000001', sensor_gene='R')) == event_id(
        dict(sender='A', receiver='B', hmdb_id='HMDB0000001', sensor_gene='R'))


@pytest.mark.parametrize('values', [[.02] * 5, [.01, .03, 0, None, .02], [None] * 5])
def test_embedding_keeps_public_figure_encoding_and_custom_options(values, monkeypatch):
    frame = path_frame().iloc[:5].copy()
    frame['fdr_sensor_type'] = values
    frame.loc[0, 'cell_mesh_score'] = float('nan')
    view = plotting.plot_communication_dotmap(frame, cmap='magma', metabolite_label='name')
    view['ax'].set_title('PUBLIC FIGURE preserved')
    original_fonts = [(text.get_text(), text.get_fontsize()) for ax in view['fig'].axes for text in ax.texts]
    def forbidden(*args, **kwargs):
        pytest.fail('HTML embedding must not prepare or draw a second plot')
    monkeypatch.setattr(matrix, 'communication_scope', forbidden)
    monkeypatch.setattr(matrix, 'dotmap_encoding', forbidden)
    monkeypatch.setattr(plotting, 'plot_communication_dotmap', forbidden)
    try:
        with mpl.rc_context({'svg.fonttype': 'none'}):
            text, record = communication_dotmap_html(view)
            assert mpl.rcParams['svg.fonttype'] == 'none'
        data = payload(text, 'selection')
        assert 'PUBLIC FIGURE preserved' in text  # SVG retains the supplied title's glyphs/comment.
        assert '<text' not in text and 'DejaVuSans' in text
        assert 'role="img" aria-label="通讯事件矩阵点图"' in text
        assert data['encoding'] == view['encoding']
        assert data['labels'] == view['labels']
        assert data['layout'] == view['layout']
        assert data['dots'] == view['dot_table'].to_dict('records')
        assert record['n_plotted'] == len(view['plot_events'])
        assert '代谢物（名称）' in text
        assert original_fonts == [(t.get_text(), t.get_fontsize()) for ax in view['fig'].axes for t in ax.texts]
        assert plt.fignum_exists(view['fig'].number)  # Formatter does not take ownership.
    finally:
        plt.close(view['fig'])


def test_tiers_use_public_figures_once_and_share_their_scope_with_evidence(tmp_path, monkeypatch):
    from cellmesh import _report_pages, _selection_evidence
    calls, scopes, evidence = [], [], []
    draw, prepare = plotting.plot_communication_dotmap, matrix.communication_scope
    def prepare_spy(*args, **kwargs):
        scope = prepare(*args, **kwargs)
        scopes.append(scope)
        return scope
    def draw_spy(*args, **kwargs):
        view = draw(*args, **kwargs)
        view['ax'].set_title('Public entry ' + str(len(calls)))
        calls.append(view)
        return view
    def evidence_spy(result, scope, **kwargs):
        evidence.append(scope)
        return ''
    monkeypatch.setattr(matrix, 'communication_scope', prepare_spy)
    monkeypatch.setattr(plotting, 'plot_communication_dotmap', draw_spy)
    monkeypatch.setattr(_report_pages, 'selection_evidence_html', evidence_spy)
    monkeypatch.setattr(_selection_evidence, 'selection_evidence_html', evidence_spy)
    result = make_result(n_perms=0)
    filtered = result.filter_events(source='all')
    selected = result.select_events(metabolite_sensor_pairs=[('HMDB0000001', 'R')], source='all')
    files = result.save(tmp_path / 'a', filtered_events=filtered, selected_events=selected, plot_top_n=1)
    assert len(calls) == len(scopes) == len(evidence) == 2
    for view, scope, used in zip(calls, scopes, evidence):
        assert view['scope'] is scope is used
        assert not plt.fignum_exists(view['fig'].number)
    pair_page = next(p for key, p in files.items() if key.startswith('pair_report:'))
    for index, (page, prefix) in enumerate([(files['filtered_report'], 'interest'), (pair_page, 'pair')]):
        text = page.read_text()
        assert 'Public entry ' + str(index) in text
        data = payload(text, prefix)
        assert data['dots'] == calls[index]['dot_table'].to_dict('records')
        assert data['records'] == calls[index]['scope']['records']
    assert calls[0]['scope']['top_n'] == 1 and calls[1]['scope']['top_n'] is None
    assert 'dotmap-data' not in files['report'].read_text()


def test_report_closes_owned_figure_if_embedding_fails(monkeypatch):
    from cellmesh import _report_dotmap
    frame = path_frame().iloc[:2].copy()
    result = cellmesh.CellMeshResult(frame, pd.DataFrame(), pd.DataFrame(), {})
    before = set(plt.get_fignums())
    def fail(*args, **kwargs):
        raise OSError('export interrupted')
    monkeypatch.setattr(_report_dotmap, 'communication_dotmap_html', fail)
    with pytest.raises(OSError, match='export interrupted'):
        render_html_report(result, selection=result.filter_events(source='all'))
    assert set(plt.get_fignums()) == before
