"""Selected CSV, report rows and figure scopes remain consistent."""
import copy
from html.parser import HTMLParser
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from cellmesh import CellMeshResult
from cellmesh._report_plots import selection_plots_html
from cellmesh.reporting import render_html_report
from test_reporting_and_persistence import make_result, assert_nested_equal


def payload(text, key="selection-plot-data"):
    return json.loads(re.search(r'id="' + key + r'">(.*?)</script>', text, re.S)[1])


@pytest.fixture
def result():
    rows = [dict(sender=s, receiver=r, hmdb_id=f"HMDB{i:07d}", sensor_gene=f"R{i}",
                 metabolite='Same name <script>bad()</script>', cell_mesh_score=i/30,
                 fdr_sensor_type=(30-i)/30, passes_min_cells=False)
            for i in range(1, 26) for s, r in [("A", "B"), ("B", "A")]]
    return CellMeshResult(pd.DataFrame(rows), pd.DataFrame(), pd.DataFrame(), {})


def fake_tables():
    return {"coverage": {"summary": pd.DataFrame(), "details": pd.DataFrame()},
            "explanation": {"receiver_components": pd.DataFrame()}}


@pytest.mark.parametrize("top_n,n_pairs", [(20, 20), (3, 3), (None, 25)])
def test_three_figures_scope_ranking_no_filter_or_mutation(result, top_n, n_pairs):
    # Every row fails the legacy default FDR/QC filters. Selection is authoritative.
    before = copy.deepcopy(result)
    focus = result.filter_events(source="all", )
    text = render_html_report(result, selection=focus, detail_n=0,
                              plot_top_n=top_n, tables=fake_tables())
    records = {r['kind']: r for r in payload(text, 'interest-plot-data')}
    assert 'id="selection-plot-data"' not in text
    assert 'id="selection-dotplot"' not in text
    assert records['network']['edge_event_total'] == 50
    assert 'counts' not in records  # Equal nonzero counts add no comparison.
    for kind in ['network']:
        assert records[kind]['plotted_event_indices'] == list(range(50))
        assert records[kind]['n_display_excluded'] == 0
    dot = records['dotmap']
    assert dot['n_selected'] == 50 and dot['n_plotted'] == n_pairs*2
    assert dot['n_display_excluded'] == (25-n_pairs)*2
    expected_ids = {f'HMDB{i:07d}' for i in range(26-n_pairs, 26)}
    assert {pair[0] for pair in dot['plotted_pairs']} == expected_ids
    assert set(focus.events.iloc[dot['plotted_event_indices']].hmdb_id) == expected_ids
    data = json.loads(re.search(r'id="interest-data">(.*?)</script>', text, re.S)[1])
    assert len(data['rows']) == 50
    assert len(re.findall(r'<div class="selection-plot"><svg ', text)) == 1
    assert 'id="interest-dotmap"' in text and 'id="interest-paths"' not in text
    assert '<script>bad()' not in text
    assert_nested_equal(result.events, before.events)
    assert not plt.get_fignums()


@pytest.mark.parametrize('missing', [0, 3, 50])
def test_missing_scores_retained_in_table_and_counted_separately(result, missing):
    result.events.loc[:missing-1, 'cell_mesh_score'] = np.nan if missing else result.events.loc[:missing-1, 'cell_mesh_score']
    focus = result.filter_events(source="all", )
    text = selection_plots_html(focus)
    for record in payload(text):
        assert record['n_selected'] == 50
        assert record['n_missing_score'] == missing
        assert record['n_plotted'] + record['n_display_excluded'] + (0 if record['kind']=='dotmap' else missing) == 50
        if record['kind']!='dotmap':
            assert focus.events.iloc[record['plotted_event_indices']].cell_mesh_score.notna().all()
    assert len(focus.events) == 50
    if missing == 50:
        assert '<svg' in text and '全部缺失评分' in text


def test_empty_and_nonfiltered_reports(result):
    text = selection_plots_html(result.filter_events(source="all", sender_labels=[]))
    assert payload(text) == []
    assert all(r['n_selected'] == r['n_plotted'] == 0 for r in payload(text))
    assert '<svg' not in text and '没有事件满足全部筛选条件' in text
    plain = render_html_report(result, tables=fake_tables())
    assert 'id="selection-plots"' in plain
    assert 'id="report-explanations"' not in plain
    assert payload(plain) == []  # the single positive cell pair has a compact summary


@pytest.mark.parametrize('value', [True, False, 0, -1, 2.5, '20', float('nan')])
def test_invalid_limits_fail_before_writing(result, tmp_path, value):
    for action in [lambda: result.filter_events(source="all", ).to_html(tmp_path/'report.html', plot_top_n=value),
                   lambda: result.save(tmp_path/'analysis', plot_top_n=value)]:
        with pytest.raises(ValueError, match='plot_top_n'):
            action()
        assert not list(tmp_path.iterdir())


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_all_public_report_routes_preserve_full_exports(tmp_path, mode):
    result = make_result(mode, n_perms=0)
    focus = result.filter_events(source="all", sender_labels=['A'])
    before = copy.deepcopy(result)
    paths = result.save(tmp_path/'analysis', report=True, selection=focus, detail_n=1,
                        plot_top_n=1, include_null_scores=True, include_all_events=True)
    focus.to_csv(tmp_path/'chosen.csv')
    direct = focus.to_html(tmp_path/'direct.html', detail_n=1, plot_top_n=1)
    alternate = result.to_html(tmp_path/'alternate.html', selection=focus, detail_n=1, plot_top_n=None)
    for path in [paths['filtered_report'], direct, alternate]:
        records = payload(path.read_text(), "interest-plot-data")
        assert all(r['n_selected'] == len(focus.events) for r in records)
        assert path.read_text().count('class="event-explanation"') == 0
    assert len(pd.read_csv(paths['events_csv'])) == len(result.events)
    assert len(pd.read_csv(tmp_path/'chosen.csv')) == len(focus.events)
    assert not list(tmp_path.glob('*.png')) and not list(tmp_path.glob('*.pdf'))
    assert_nested_equal(result.events, before.events)
    assert_nested_equal(CellMeshResult.load(paths['result']).events, result.events)


def test_svg_ids_are_unique_and_figures_closed_on_failure(result, monkeypatch):
    class Parser(HTMLParser):
        def __init__(self, text):
            super().__init__(); self.ids = []; self.feed(text)
        def handle_starttag(self, tag, attrs):
            self.ids += [v for k,v in attrs if k == 'id']
    ids = Parser(selection_plots_html(result.filter_events(source="all", ))).ids
    assert len(ids) == len(set(ids))
    from matplotlib.figure import Figure
    def fail(*args, **kwargs): raise OSError('disk/render failure')
    monkeypatch.setattr(Figure, 'savefig', fail)
    with pytest.raises(OSError, match='render failure'):
        selection_plots_html(result.filter_events(source="all", ))
    assert not plt.get_fignums()


def test_svg_namespacing_preserves_arbitrary_label_text(result):
    label = 'id="literal" href="#literal" url(#literal)'
    result.events['metabolite'] = label
    # Default metabolite labels now use HMDB IDs; exercise visible cell labels.
    result.events['receiver'] = label
    text = selection_plots_html(result.filter_events(source="all", ), plot_top_n=1)
    assert payload(text, 'selection-dotmap-data')['records'][0]['metabolite'].startswith(label)
    assert all(part in text for part in label.split())
    assert 'id="selection-dotplot-literal"' not in text


def test_constant_event_counts_use_integer_colorbar_ticks(monkeypatch):
    from matplotlib.figure import Figure
    rows = [dict(sender='A', receiver=receiver, hmdb_id='HMDB0000001',
                 sensor_gene='R', cell_mesh_score=0.5, fdr_sensor_type=0.01)
            for receiver in ['B', 'C']]
    result = CellMeshResult(pd.DataFrame(rows), pd.DataFrame(), pd.DataFrame(), {})
    ticks = []
    original_savefig = Figure.savefig
    def capture(fig, *args, **kwargs):
        if len(fig.axes) == 2:
            ticks.extend(fig.axes[-1].get_yticks())
        return original_savefig(fig, *args, **kwargs)
    monkeypatch.setattr(Figure, 'savefig', capture)
    selection_plots_html(result.filter_events(), include_dotplot=False)
    assert ticks and all(float(value).is_integer() for value in ticks)
