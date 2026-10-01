"""Presentation compaction preserves scientific values, links and evidence."""
from html.parser import HTMLParser
import json
import re

import h5py
import pytest

from cellmesh.reporting import render_html_report
from cellmesh.result_io import _read
from test_reporting_and_persistence import make_result, SELECTOR


class Structure(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.stack = []
        self.h2 = []
        self.charts = []
        self.visible_svgs = 0
        self.open_explanations = []
        self.evidence_tables = []
        self.ids = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        data = dict(attrs)
        hidden = any(t == 'details' and 'open' not in a for t, a in self.stack)
        if tag == 'svg' and not hidden:
            self.visible_svgs += 1
        if 'id' in data:
            self.ids.append(data['id'])
        if data.get('class') == 'selection-chart':
            self.charts.append((data['id'], 'open' in data))
        if data.get('class') == 'event-explanation':
            self.open_explanations.append('open' in data)
        if 'data-table' in data:
            self.evidence_tables.append((data['data-table'], 'open' in data))
        if tag not in ['meta', 'input', 'br', 'hr', 'img', 'link']:
            self.stack.append((tag, data))

    def handle_endtag(self, tag):
        for i in range(len(self.stack)-1, -1, -1):
            if self.stack[i][0] == tag:
                del self.stack[i:]
                break

    def handle_data(self, data):
        if self.stack and self.stack[-1][0] == 'h2':
            self.h2.append(data)


@pytest.mark.parametrize('selected', [True, False])
def test_compact_sections_keep_all_explanations_and_evidence(selected):
    result = make_result(n_perms=0)
    opts = {'selection': result.filter_events(source="all", ), 'detail_n': 1} if selected else {'events': [SELECTOR]}
    text = render_html_report(result, **opts)
    view = Structure(text)
    assert len(view.ids) == len(set(view.ids))
    assert view.open_explanations == []
    if selected:
        assert view.h2 == ['筛选概览', '通讯特征', '所选事件列表', '筛选与导出']
        assert view.charts == [('interest-dotmap', True), ('interest-network', False),
                               ('interest-samples', False), ('interest-expression', False)]
        assert view.evidence_tables == []
        assert 'class="sample-component"' not in text
    else:
        assert view.h2 == ['事件概览', '评分组成', '样本表现', '表达证据与先验覆盖', '数据明细与来源']
        assert {name for name, _ in view.evidence_tables} == {
            'aggregate_event', 'sample_events', 'sender_components', 'receiver_components',
            'reference_context', 'reaction_details', 'gene_details', 'enzyme_prior', 'interaction_prior'}
        assert all(not open_ for _, open_ in view.evidence_tables)
        assert text.count('class="explanation-panels"') == 3
    assert 'id="explanation-directory"' not in text
    assert 'id="report-coverage"' not in text
    assert 'Coverage describes genes' not in text


def test_chinese_headers_leave_raw_columns_and_scientific_values_intact():
    result = make_result(n_perms=0)
    text = render_html_report(result)
    data = json.loads(re.search(r'id="event-data">(.*?)</script>', text, re.S)[1])
    assert data['columns'] == list(result.events.columns)
    assert data['labels'][data['columns'].index('sender')] == '发送细胞'
    for row, (_, expected) in zip(data['rows'], result.events.iterrows()):
        score = row['values'][data['columns'].index('cell_mesh_score')]
        assert float(score) == expected.cell_mesh_score
    assert 'title="sender" disabled>发送细胞' in text


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_single_event_compact_details_preserve_all_records_and_missingness(mode):
    from cellmesh._report_browser import _cell_text
    result = make_result(mode, n_perms=0)
    detail = result.explain_event(**SELECTOR)
    text = render_html_report(result, events=[SELECTOR])
    assert text.count('id="event-table-help"') == 1
    assert text.count('<noscript>') == 1
    assert '<details id="event-score-method">' in text
    assert '完整 CSV' not in text
    for name, frame in detail.items():
        if not hasattr(frame, 'columns'):
            continue
        match = re.search(r'id="detail-' + name.replace('_', '-') + '-data">(.*?)</script>', text, re.S)
        assert match, name
        data = json.loads(match[1])
        assert data['row_unit'] == '记录' and data['compact']
        assert data['columns'] == list(frame.columns)
        assert [row['values'] for row in data['rows']] == [
            [_cell_text(value) for value in row] for row in frame.itertuples(index=False, name=None)]
    assert ('id="event-sample"' in text) == (mode == 'sample_aware')
    if mode == 'sample_aware':
        assert '不表示统计显著' in text
        assert 'receiver_unobserved' in text  # Missing samples remain in the complete data.


@pytest.mark.parametrize('mode', ['pooled_stratified', 'sample_aware'])
def test_stored_report_tables_render_identically_to_fresh_tables(tmp_path, mode):
    result = make_result(mode, n_perms=0)
    paths = result.save(tmp_path/'analysis')
    with h5py.File(paths['result'], 'r') as handle:
        tables = {group: {name: _read(node) for name, node in handle[group].items()}
                  for group in ['explanation', 'coverage']}
    focus = result.filter_events(source="all", )
    cached = render_html_report(result, selection=focus, detail_n=2, tables=tables)
    fresh = render_html_report(result, selection=focus, detail_n=2)
    assert cached == fresh
