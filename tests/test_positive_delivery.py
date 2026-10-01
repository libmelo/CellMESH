"""Two-level reporting: positive calls and independent interest selections."""
import copy
import json
import re

import h5py
import numpy as np
import pandas as pd
import pytest

from cellmesh import CellMeshResult
from cellmesh.positivity import positive_calls, resolve_positive
from cellmesh.reporting import render_html_report
from test_reporting_and_persistence import make_result, assert_nested_equal


@pytest.fixture
def result():
    result = make_result(n_perms=0)
    result.events['perm_pvalue'] = [.001, .03, .5, np.nan]
    result.events['fdr_sensor_type'] = [.01, .05, .7, np.nan]
    result.events['fdr_global'] = [.04, .08, .9, np.nan]
    result.events['passes_min_cells'] = [False, True, True, True]
    return result


def data(text, key='event-data'):
    return json.loads(re.search('id="' + key + '">(.*?)</script>', text, re.S)[1])


def test_positive_scope_and_defaults(result, tmp_path):
    before = copy.deepcopy(result.events)
    paths = result.save(tmp_path/'analysis')
    assert set(paths) == {'result', 'positive_csv', 'report', 'navigation'}
    assert not (tmp_path/'analysis.all_events.csv').exists()
    assert len(pd.read_csv(paths['positive_csv'])) == 2
    html = paths['report'].read_text()
    assert len(data(html)['rows']) == 2 and 'id="report-explanations"' not in html
    assert '阳性事件 CSV' in html
    assert 'analysis.all_events.csv' not in html
    assert result.delivery['positive_summary'] == dict(n_total=4, n_positive=2, n_not_positive=1,
                                                       n_not_evaluated=1, status='partial')
    assert_nested_equal(before, result.events)


def test_positive_csv_without_report_and_config_inheritance(result, tmp_path):
    paths = result.save(tmp_path/'p', report=False, positive_col='perm_pvalue', positive_threshold=.01)
    assert set(paths) == {'result', 'positive_csv'}
    assert len(pd.read_csv(paths['positive_csv'])) == 1
    restored = CellMeshResult.load(paths['result'])
    assert len(restored.filter_events().events) == 1
    assert restored.filter_events().criteria['positive']['column'] == 'perm_pvalue'
    assert restored.get_selection() is None


def test_selection_independent_from_general_scope(result, tmp_path):
    focus = result.filter_events(source='all')
    paths = result.save(tmp_path/'a', report=True, selection=focus, detail_n=1)
    html = paths['report'].read_text()
    assert len(data(html)['rows']) == 2
    selected_html = paths['filtered_report'].read_text()
    assert 'interest-data' not in html
    assert len(data(selected_html, 'interest-data')['rows']) == 4
    selected = pd.read_csv(paths['filtered_csv'])
    assert selected.positive_status.tolist() == ['positive', 'positive', 'not_positive', 'not_evaluated']
    assert len(data(html, 'selection-plot-data')[0]['plotted_event_indices']) == 2
    assert len(data(selected_html, 'interest-plot-data')[0]['plotted_event_indices']) == 4
    assert selected_html.count('class="event-explanation"') == 0
    # The static interest preview must also link to real sections without JS.
    assert 'href="#interest-explanation-' not in html
    interest_preview = selected_html.split('id="interest-browser"', 1)[1].split('id="interest-data"', 1)[0]
    assert 'result.add_event_reports' in interest_preview
    assert CellMeshResult.load(paths['result']).get_selection().criteria == focus.criteria


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_one_step_and_append_are_identical(mode, tmp_path):
    r = make_result(mode, n_perms=0)
    r.events['perm_pvalue'] = .01
    r.events['fdr_sensor_type'] = [.01, .02, .4, .8]
    original = copy.deepcopy(r)
    (tmp_path/'one').mkdir(); (tmp_path/'two').mkdir()
    one = r.save(tmp_path/'one/analysis', report=True, selection=r.filter_events(), detail_n=1)
    two = original.save(tmp_path/'two/analysis', detail_n=1)
    unchanged = {k: two[k].read_bytes() for k in ['positive_csv', 'report']}
    loaded = CellMeshResult.load(two['result'])
    assert loaded.get_selection() is None
    out = loaded.add_selection(loaded.filter_events())
    for key in ['report', 'filtered_report', 'filtered_csv', 'positive_csv']:
        assert one[key].read_bytes() == out[key].read_bytes(), key
    for key, value in unchanged.items():
        assert two[key].read_bytes() == value
    restored = CellMeshResult.load(out['result'])
    assert restored.delivery == r.delivery
    pd.testing.assert_frame_equal(restored.events, original.events, check_exact=True)
    with pytest.raises(ValueError, match='replace=True'):
        restored.add_selection(restored.filter_events())
    restored.add_selection(restored.filter_events(sender_labels=[]), replace=True, detail_n=0)
    assert pd.read_csv(out['filtered_csv']).empty
    assert 'id="report-explanations"' not in out['filtered_report'].read_text()


def test_no_statistics_differs_from_no_positive(tmp_path):
    r = make_result(n_perms=0)
    p = r.save(tmp_path/'none', report=True)
    assert pd.read_csv(p['positive_csv']).empty
    assert '尚未完成达标判定' in p['report'].read_text()
    assert r.filter_events().events.empty and len(r.filter_events(source='all').events) == len(r.events)
    r.events['fdr_sensor_type'] = 1.
    p = r.save(tmp_path/'negative', report=True)
    assert '已完成判定的事件中，没有事件达到本次标准' in p['report'].read_text()


def test_qc_is_explicit_and_pvalue_threshold_supported(result):
    assert len(result.filter_events().events) == 2
    assert len(result.filter_events(positive_qc_only=True).events) == 1
    assert len(result.filter_events(source='all', max_perm_pvalue=.001).events) == 1
    assert len(result.filter_events(positive_col='fdr_global').events) == 1


def test_mismatched_selection_criterion_rejected(result, tmp_path):
    focus = result.filter_events()
    with pytest.raises(ValueError, match='different positive criterion'):
        result.save(tmp_path/'bad', selection=focus, positive_threshold=.01)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize('options', [{'positive_col':'bad'}, {'positive_threshold':True},
                                  {'positive_threshold':np.nan}, {'positive_threshold':1.1},
                                  {'positive_qc_only':'yes'}])
def test_invalid_positive_options_before_writes(result, tmp_path, options):
    with pytest.raises((ValueError, TypeError)):
        result.save(tmp_path/'bad', **options)
    assert not list(tmp_path.iterdir())


def test_append_rollback_and_stale_delivery(result, tmp_path, monkeypatch):
    import cellmesh.result_io as io
    paths = result.save(tmp_path/'a', report=True)
    before = {k: p.read_bytes() for k,p in paths.items()}
    previous = copy.deepcopy(result.delivery)
    original = io.os.replace
    failed = False
    def fail(source, target):
        nonlocal failed
        if str(target).endswith('.report.html') and not failed:
            failed = True
            raise OSError('publish failure')
        return original(source,target)
    with monkeypatch.context() as context:
        context.setattr(io.os, 'replace', fail)
        with pytest.raises(OSError, match='publish failure'):
            result.add_selection(result.filter_events(), detail_n=0)
    assert result.delivery == previous
    for key, content in before.items():
        assert paths[key].read_bytes() == content
    assert not (tmp_path/'a.filtered_events.csv').exists()
    stale = CellMeshResult.load(paths['result'])
    result.add_selection(result.filter_events(), detail_n=0)
    with pytest.raises(ValueError, match='reload'):
        stale.add_selection(stale.filter_events(), detail_n=0)


def test_standalone_report_exports_working_csvs(result, tmp_path):
    path = result.to_html(tmp_path/'review.html', selection=result.filter_events(), detail_n=0)
    assert 'review.filtered_events.csv' in path.read_text()
    assert not (tmp_path/'review.positive_events.csv').exists()
    assert not (tmp_path/'review.all_events.csv').exists()
    assert len(pd.read_csv(tmp_path/'review.filtered_events.csv')) == 2


def test_append_rejects_modified_general_csv(result, tmp_path):
    paths = result.save(tmp_path/'a')
    paths['positive_csv'].write_text('different data')
    before = paths['result'].read_bytes()
    with pytest.raises(ValueError, match='general CSV differs'):
        result.add_selection(result.filter_events(), detail_n=0)
    assert paths['result'].read_bytes() == before


def test_append_none_limits_and_clear_optional_score_gate(result, tmp_path):
    paths = result.save(tmp_path/'a', positive_min_score=.9, detail_n=1)
    assert result.filter_events().criteria['positive']['min_score'] == .9
    assert result.filter_events(positive_min_score=None).criteria['positive']['min_score'] is None
    result.add_selection(result.filter_events(source='all'), detail_n=None, plot_top_n=None)
    assert paths['result'].exists()
    assert 'detail_n' not in result.delivery and result.delivery['plot_top_n'] is None
