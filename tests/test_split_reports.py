"""Separate delivery scopes, expression reuse and complete sample display."""
import copy
import json
import re

import h5py
import matplotlib.pyplot as plt
import pandas as pd
import pytest

from cellmesh import CellMeshResult
from cellmesh._report_evidence import expression_evidence_html, sample_scores_svg
from cellmesh.result_io import _read, _write
from test_input_archive import inputs, run
from test_report_layout import Structure
from test_reporting_and_persistence import make_result


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_shared_expression_evidence_and_report_modes(inputs, tmp_path, mode):
    result = run(inputs, sample_mode=mode, n_perms=0, retain_input=True)
    before = copy.deepcopy(result.events)
    events = result.events[['sender','receiver','hmdb_id','sensor_gene']].to_dict('records')
    pairs = list(result.events[['hmdb_id','sensor_gene']].drop_duplicates().itertuples(index=False, name=None))
    selected_pairs = result.select_events(metabolite_sensor_pairs=pairs, source='all')
    paths = result.save(tmp_path/'a', filtered_events=result.filter_events(source='all'), selected_events=selected_pairs, include_input=True)
    general, selected = (paths[key].read_text() for key in ['report', 'filtered_report'])
    assert len(paths) == 8 + len(pairs)
    assert 'interest-data' not in general and 'id="event-data"' not in selected
    assert 'report-coverage' in general and 'report-coverage' not in selected
    assert selected.count('class="expression-evidence"') == 0
    for key, path in paths.items():
        if not key.startswith('pair_report:'):
            continue
        text = path.read_text()
        assert text.count('class="selection-expression-group"') == 2
        assert '该表达图不可用' not in text
        assert ('id="pair-samples"' in text) == (mode == 'sample_aware')
        ids = Structure(text).ids
        assert len(ids) == len(set(ids))
    assert not plt.get_fignums()
    pd.testing.assert_frame_equal(result.events, before)
    loaded = CellMeshResult.load(paths['result'])
    selected = loaded.select_events(metabolite_sensor_pairs=pairs, source='all')
    without = selected.to_html(tmp_path/'without.html')
    page = next(p for key,p in without.items() if key.startswith('pair_report:'))
    assert '未提供细胞级表达输入' in page.read_text()
    supplied = selected.to_html(tmp_path/'supplied.html', adata=loaded.load_input())
    page = next(p for key,p in supplied.items() if key.startswith('pair_report:'))
    assert page.read_text().count('class="selection-expression-group"') == 2
    assert '该表达图不可用' not in page.read_text()



def test_selection_does_not_change_general_and_report_false_closes_both(tmp_path):
    r = make_result(n_perms=0)
    p = r.save(tmp_path/'a')
    original = {k: p[k].read_bytes() for k in ['report', 'positive_csv']}
    p = r.add_selection(r.filter_events(source='all'), detail_n=1)
    assert all(p[k].read_bytes() == value for k, value in original.items())
    assert set(p) == {'result', 'report', 'filtered_report', 'filtered_csv', 'positive_csv', 'navigation'}
    r.save(tmp_path/'a', report=False, selection=r.get_selection())
    assert not p['report'].exists() and not p['filtered_report'].exists()
    out = r.add_selection(r.filter_events(source='all'), replace=True)
    assert 'report' not in out and 'filtered_report' not in out
    assert not list(tmp_path.glob('*.html'))


def test_migrates_old_combined_report_once(tmp_path):
    r = make_result(n_perms=0)
    p = r.save(tmp_path/'old', selection=r.filter_events(source='all'), detail_n=0)
    p['report'].write_text('<html>old combined report</html>')
    p['filtered_report'].unlink()
    r.delivery.pop('report_layout')
    with h5py.File(p['result'], 'r+') as h:
        del h['metadata/delivery']
        _write(h['metadata'], 'delivery', r.delivery, include_null_scores=False)
    r.add_selection(r.filter_events(source='all'), replace=True)
    assert 'three-level-v3' in p['report'].read_text() and 'interest-data' not in p['report'].read_text()
    assert 'interest-data' in p['filtered_report'].read_text()
    before = p['report'].read_bytes()
    r.add_selection(r.filter_events(source='all'), replace=True)
    assert p['report'].read_bytes() == before


def test_all_sample_scores_include_na_zero_failed_qc_without_truncation():
    frame = pd.DataFrame({'sample': [f's{i:02d}' for i in range(15)],
                          'cell_mesh_score': [0., None] + [.2]*13,
                          'passes_min_cells': [True, False, None] + [True]*12})
    text = sample_scores_svg(frame)
    assert all(f's{i:02d}' in text for i in range(15))
    assert 'NA · 质控未通过' in text and '0 · 质控通过' in text and '质控未记录' in text
    assert '中位数（0.2）' in text


def test_failed_selected_render_publishes_nothing(tmp_path, monkeypatch):
    r = make_result(n_perms=0)
    p = r.save(tmp_path/'a')
    original = {k: path.read_bytes() for k, path in p.items()}
    import cellmesh._report_pages as evidence
    def fail(*args, **kwargs):
        raise RuntimeError('evidence render failed')
    monkeypatch.setattr(evidence, 'render_overview', fail)
    with pytest.raises(RuntimeError, match='evidence render failed'):
        r.add_selection(r.filter_events(source='all'))
    assert all(p[k].read_bytes() == value for k, value in original.items())
    assert not (tmp_path/'a.filtered_events.csv').exists()


def test_selected_paths_protected_and_empty_selection(tmp_path):
    r = make_result(n_perms=0)
    p = r.save(tmp_path/'a', selection=r.filter_events(source='all', sender_labels=[]))
    assert pd.read_csv(p['filtered_csv']).empty
    assert '没有事件满足全部筛选条件' in p['filtered_report'].read_text()
    with pytest.raises(ValueError, match='standard result delivery'):
        r.get_selection().to_csv(p['filtered_report'])
    with pytest.raises(ValueError, match='overlap'):
        r.get_selection().to_html(p['filtered_report'])


def test_data_only_append_does_not_require_explanation_capture(tmp_path):
    original = make_result(n_perms=0)
    r = CellMeshResult(events=original.events, sender_scores=original.sender_scores,
                       receiver_scores=original.receiver_scores, parameters=original.parameters)
    r.save(tmp_path/'data_only', report=False)
    out = r.add_selection(r.filter_events(source='all'))
    assert len(pd.read_csv(out['filtered_csv'])) == len(r.events)
    assert not list(tmp_path.glob('*.html'))


def test_zero_expression_colorbar_uses_nonnegative_range(inputs):
    from cellmesh._report_evidence import _figure_svg
    from cellmesh import plot_receptor_expression_violin
    inputs[0].X[:, 2] = 0
    r = run(inputs, n_perms=0)
    captured = {}
    def draw(**options):
        output = plot_receptor_expression_violin(**options)
        captured.update(output)
        return output
    _figure_svg(draw, 'zero', result=r, adata=inputs[0], hmdb_id='HMDB1', receptor_gene='R')
    norm = captured['colorbar'].mappable.norm
    assert (norm.vmin, norm.vmax) == (0, 1)
    assert captured['plot_data'].receptor_expression.eq(0).all()
    assert captured['ax'].get_ylim()[0] == 0
    assert not plt.get_fignums()
