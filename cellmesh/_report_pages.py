"""Three independent offline report scopes, using captured scientific values."""
from __future__ import annotations

import html
import json
from urllib.parse import quote

import numpy as np
import pandas as pd

from ._report_browser import (BROWSER_CSS, BROWSER_JS, _script_json, event_browser,
                              event_id, page_entries)
from ._report_layout import (LAYOUT_CSS, COMPONENT_LABELS, PANEL_TITLES, REASON_LABELS,
                            TABLE_LABELS, overview_html, methods_html, csv_links,
                            selection_overview_html, metric_card, missing_criteria_html,
                            percent_text)
from .positivity import (INHERIT, STATUS_LABELS, STATISTIC_LABELS, positive_calls, positive_description,
                         positive_selection, selection_table, resolve_positive, delivery_record)
from ._selection_evidence import EVIDENCE_CSS, EVIDENCE_JS, selection_evidence_html


PAGE_CSS = """
svg{max-width:100%;height:auto}svg[viewBox]:not([width]){width:100%}
.sample-component[hidden]{display:none}.sample-controls{margin:16px 0}
.sample-controls select{font:inherit;padding:6px}.event-path{font-size:19px;font-weight:600}
.event-request pre{min-width:260px;white-space:pre-wrap}.report-id{font-family:monospace}
.sample-score-plot [data-sample-index]{cursor:pointer}
.sample-score-plot [aria-current=true]{fill:#c16b3c;font-weight:bold}
.event-stat-cards{grid-template-columns:repeat(auto-fit,minmax(min(280px,100%),1fr))}
.event-stat-cards .overview-card-primary strong{font-size:clamp(26px,3vw,34px)}
.component-inputs dl{display:grid;grid-template-columns:minmax(150px,1fr) 2fr;gap:5px 18px;margin:8px 0}
.component-inputs dt{color:#51616b}.component-inputs dd{margin:0;overflow-wrap:anywhere}
.event-source{overflow-wrap:anywhere}
"""
SAMPLE_JS = r"""
(() => {
 const select = document.getElementById('event-sample');
 if (!select) return;
 function show(value) {
   select.value = String(value);
   document.querySelectorAll('[data-component-sample]').forEach(panel => {
     panel.hidden = panel.dataset.componentSample !== select.value;
   });
   document.querySelectorAll('[data-sample-index]').forEach(node => {
     node.setAttribute('aria-current', String(node.dataset.sampleIndex === select.value));
   });
 }
 select.addEventListener('change', () => show(select.value));
 document.querySelectorAll('[data-sample-index]').forEach(node => {
   const choose = () => {show(node.dataset.sampleIndex); select.focus();};
   node.addEventListener('click', choose);
   node.addEventListener('keydown', event => {
     if (event.key === 'Enter' || event.key === ' ') {event.preventDefault(); choose();}
   });
 });
 show(select.value);
})();
"""


def document(title, sections, *, scope, nav=(), general_report=None, navigation_report=None):
    links = ''.join('<a href="' + html.escape(quote(str(path)), quote=True) + '">' + label + '</a> '
                    for path, label in [(general_report, '一级统计报告'), (navigation_report, '报告导航')]
                    if path)
    menu = '<nav class="report-nav" aria-label="报告导航">' + ''.join(
        f'<a href="#{key}">{label}</a>' for key, label in nav) + '</nav>' if nav else ''
    return ("<!doctype html><html lang='zh'><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width, initial-scale=1'>"
            '<meta name="cellmesh-report-layout" content="' + ('three-level-v3' if scope == 'general' else 'independent-scopes-v4') + '">'
            f'<meta name="cellmesh-report-scope" content="{scope}">'
            '<title>' + html.escape(title) + '</title><style>'
            'body{font:15px/1.6 system-ui,sans-serif;margin:clamp(12px,3vw,32px);color:#22313b}'
            'main{max-width:1200px;margin:auto}h1,h2,h3{color:#165d69}.table{overflow:auto}'
            'table{border-collapse:collapse;font-size:13px}th,td{padding:7px;border:1px solid #dde4e8;vertical-align:top}'
            'th{background:#eef5f6}details{margin:14px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere}a{color:#165d69}'
            + BROWSER_CSS + LAYOUT_CSS + PAGE_CSS + (EVIDENCE_CSS if scope in ('selection', 'metabolite_sensor') else '') + '</style></head><body><main><header><h1>' +
            html.escape(title) + '</h1>' + ('<p>' + links + '</p>' if links else '') + menu +
            '</header>' + ''.join(sections) + '</main><script id="event-browser-script">' +
            BROWSER_JS + '</script><script id="event-sample-script">' + SAMPLE_JS + '</script>' +
            ('<script id="communication-features-script">' + EVIDENCE_JS + '</script>' if scope in ('selection', 'metabolite_sensor') else '') + '</body></html>')


def full_table(frame, name, title):
    return event_browser(frame, [], prefix='detail-' + name.replace('_', '-'), title=title,
                         show_details=False, event_table=False, compact=True)


def sample_heatmap(selection, plot_top_n):
    """All original sample rows; bound only the event columns, without rescoring."""
    result = selection._result
    if result.parameters.get('sample_mode') != 'sample_aware' or len(selection.events) <= 1:
        return ''
    from .selection import EVENT_KEYS
    from .reporting import _units
    from matplotlib import colormaps
    from matplotlib.colors import to_hex
    chosen = selection.events if plot_top_n is None else selection.events.head(plot_top_n)
    keys = list(chosen[EVENT_KEYS].itertuples(index=False, name=None))
    samples = [sample for sample, _ in _units(result)]
    source = result.sample_events
    if source is None:
        return '<p>未保存样本事件评分，无法绘制样本热图。</p>'
    subset = source.loc[source.hmdb_id.isin(chosen.hmdb_id) & source.sensor_gene.isin(chosen.sensor_gene)]
    lookup = {tuple(row[k] for k in EVENT_KEYS) + (row['sample'],): row
              for row in subset.to_dict('records')}
    width, height = 200 + 42 * len(keys), 100 + 28 * len(samples)
    parts = ['<details class="selection-chart" id="interest-samples"><summary>样本 × 事件评分热图</summary>',
             f'<p>所选 {len(selection.events):,} 个事件；展示 {len(keys):,} 个事件、全部 {len(samples):,} 个样本。'
             '按当前事件排序取展示列；评分色标固定为 0–1，灰色为 NA，× 为质控未通过，? 为质控未记录。</p>',
             '<div class="selection-plot">',
             f'<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-label="样本事件评分" '
             f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    rows = []
    for j, key in enumerate(keys):
        parts.append(f'<text x="{180+42*j}" y="25" font-size="12">E{j+1}</text>')
    for i, sample in enumerate(samples):
        y = 42 + i * 28
        parts.append(f'<text x="4" y="{y+16}" font-size="12">{html.escape(str(sample))}</text>')
        scores, qcs = [], []
        for j, key in enumerate(keys):
            row = lookup.get(key + (sample,), {})
            value = row.get('cell_mesh_score', np.nan)
            qc = row.get('passes_min_cells', None)
            score = None if pd.isna(value) else float(value)
            qc = None if pd.isna(qc) else bool(qc)
            scores.append(score); qcs.append(qc)
            color = '#d9dfe3' if score is None else to_hex(colormaps['viridis'](score))
            label = 'NA' if score is None else f'{score:.6g}'
            text = 'NA' if score is None else '0' if score == 0 else ''
            text += '?' if qc is None else '×' if not qc else ''
            parts.append(f'<g><title>{html.escape(str(sample))} / E{j+1}: {label}; QC={qc}</title>'
                         f'<rect x="{170+42*j}" y="{y}" width="38" height="24" fill="{color}"/>'
                         f'<text x="{173+42*j}" y="{y+17}" fill="{"#22313b" if score is None else "white"}" '
                         f'font-size="11">{text}</text></g>')
        rows.append(dict(sample=sample, scores=scores, qc=qcs))
    y = height - 19
    for j in range(11):
        color = to_hex(colormaps['viridis'](j / 10))
        parts.append(f'<rect x="{170+12*j}" y="{y-12}" width="12" height="12" fill="{color}"/>')
    parts += [f'<text x="170" y="{y+14}" font-size="11">0</text><text x="294" y="{y+14}" font-size="11">1</text>',
              '</svg></div><details><summary>事件列标识</summary><ol>']
    parts += ['<li>' + html.escape(f'{key[0]} → {key[2]} / {key[3]} → {key[1]}') + '</li>' for key in keys]
    parts += ['</ol></details><script type="application/json" id="sample-heatmap-data">',
              _script_json(dict(events=[dict(zip(EVENT_KEYS, key)) for key in keys], rows=rows,
                                n_total=len(selection.events), n_plotted=len(keys))), '</script></details>']
    return ''.join(parts)


def render_overview(result, *, selection=None, tables=None, delivery=None, csv_paths=None,
                    plot_top_n=20, input_storage=None, event_links=None, general_report=None,
                    navigation_report=None, adata=None, selection_evidence_path=None):
    from .reporting import build_report_tables, _table, _input_archive_html
    from ._report_plots import selection_plots_html
    spec = delivery['positive']
    selected = selection is not None
    if not selected:
        tables = build_report_tables(result) if tables is None else tables
        positive = positive_selection(result, spec)
        frame = positive.events
        sections = [csv_links(csv_paths, ['events_csv', 'positive_csv']),
                    '<section class="report-section" id="report-overview"><h2>分析概览</h2>',
                    overview_html(result, frame, spec, delivery['positive_summary']), '</section>',
                    '<section class="report-section" id="report-communication"><h2>通讯概览</h2>',
                    selection_plots_html(positive, plot_top_n, fdr_col=spec['column'], include_dotplot=False),
                    '</section><section class="report-section" id="report-events"><h2>达标事件列表</h2>',
                    event_browser(frame, [], show_details=False), '</section>']
        nav = [('report-overview', '分析概览'), ('report-communication', '通讯概览'),
               ('event-browser', '达标事件列表'), ('report-methods', '数据与方法')]
    else:
        from .plotting import plot_communication_dotmap
        import matplotlib.pyplot as plt
        frame = selection_table(selection, spec)
        dotmap = plot_communication_dotmap(frame, top_n=plot_top_n, fdr_col=spec['column'])
        scope = dotmap['scope']
        try:
            plots = selection_plots_html(selection, plot_top_n, prefix='interest', scope_label='条件筛选事件',
                                         fdr_col=spec['column'], dotmap=dotmap, secondary=True)
        finally:
            plt.close(dotmap['fig'])
        from ._selection_evidence import saved_selection_evidence
        evidence = (saved_selection_evidence(selection_evidence_path, result, scope)
                    if selection_evidence_path is not None else None)
        if evidence is None:
            evidence = selection_evidence_html(result, scope, adata=adata)
        sections = [csv_links(csv_paths, ['filtered_csv']),
                    '<section class="report-section" id="report-overview"><h2>筛选概览</h2>',
                    selection_overview_html(result, selection, frame, spec), '</section>',
                    '<section class="report-section" id="report-comparison"><h2>通讯特征</h2>',
                    plots,
                    sample_heatmap(selection, plot_top_n), evidence,
                    '</section><section class="report-section" id="report-events"><h2>所选事件列表</h2>']
        frame['positive_status'] = frame.positive_status.map(STATUS_LABELS)
        entries = event_links or []
        sections += [event_browser(frame, entries, prefix='interest', title='所选事件列表'), '</section>']
        nav = [('report-overview', '筛选概览'), ('report-comparison', '通讯特征'),
               ('interest-browser', '所选事件列表'), ('report-methods', '筛选与导出')]
    sections.append(methods_html(result, tables or {}, selection, _input_archive_html(result, input_storage),
                                 _table, delivery, selected=selected))
    return document('CELL MESH 通讯筛选报告' if selected else 'CELL MESH 通讯分析报告', sections,
                    scope='selection' if selected else 'general', nav=nav, general_report=general_report,
                    navigation_report=navigation_report)


def _value(value):
    return 'NA' if pd.api.types.is_scalar(value) and pd.isna(value) else html.escape(str(value))


def _stat_value(value):
    if pd.api.types.is_scalar(value) and pd.isna(value):
        return '—'
    if isinstance(value, (float, np.floating)) and np.isfinite(value):
        return f'{value:.6g}'
    return str(value)


def _event_decision_html(result, aggregate, spec, status):
    if status == 'not_evaluated':
        return missing_criteria_html(result, aggregate, spec)
    if status == 'positive':
        return ''
    # positive_calls has established that all enabled inputs are known.
    row = aggregate.iloc[0]
    reasons = []
    if float(row[spec['column']]) > spec['threshold']:
        reasons.append(STATISTIC_LABELS[spec['column']] + '超过阈值')
    if spec['qc_only'] and not row['passes_min_cells']:
        reasons.append('未通过细胞数 QC')
    if spec['min_score'] is not None and float(row['cell_mesh_score']) < spec['min_score']:
        reasons.append('通讯评分低于门槛')
    return '<p class="overview-status">未达标原因：' + html.escape('；'.join(reasons)) + '。</p>'


def _event_sample_summary(row):
    """Use captured aggregate fields; do not introduce sample-level inference."""
    if row is None:
        return ''
    parts = []
    positive, computable = row.get('n_samples_positive'), row.get('n_samples_coobserved')
    if pd.notna(positive):
        text = f'正评分样本：{int(positive):,}'
        if pd.notna(computable) and computable > 0:
            text += f' / {int(computable):,}'
            prevalence = row.get('event_prevalence')
            if pd.notna(prevalence):
                text += '（' + percent_text(prevalence) + '）'
        parts.append(text)
    qc = row.get('n_samples_passing_min_cells')
    if pd.notna(qc):
        parts.append(f'通过细胞数 QC：{int(qc):,} 个样本')
    if not parts:
        return ''
    if pd.notna(positive):
        parts.append('正评分指评分 > 0，比例以可计算样本为分母，不表示统计显著')
    return '<p class="event-sample-summary">' + html.escape('；'.join(parts)) + '。</p>'


def _component_inputs(current):
    """Keep raw inputs available without repeating field names in the main text."""
    groups = [('sender_components', '发送端', [
        ('cell_fraction', '细胞比例'), ('abundance_weight', '丰度权重'),
        ('P_unweighted', '加权前生产分数'), ('P', '加权后生产分数'),
        ('P_ref', '生产参考值'), ('C_ref', '消耗参考值'), ('E_ref', '外排参考值')]),
        ('receiver_components', '接收端', [
        ('mean_expression', '平均表达'), ('expr_frac', '表达细胞比例'),
        ('reference', '表达参考值'), ('reference_method', '参考值计算方法'),
        ('gate_applied', '是否触发表达比例门控')])]
    parts = ['<details class="component-inputs"><summary>计算输入与参考值</summary>']
    for name, title, fields in groups:
        if current[name].empty:
            continue
        component = current[name].iloc[0]
        parts.append('<h4>' + title + '</h4><dl>')
        for field, label in fields:
            if field not in component:
                continue
            value = component[field]
            display = _stat_value(value)
            if field == 'reference_method':
                display = {'mean': '平均数', 'median': '中位数'}.get(value, display)
            elif field == 'gate_applied' and pd.notna(value):
                display = '是' if value else '否'
            parts.append('<dt title="' + field + '">' + label + '</dt><dd title="' +
                         html.escape(str(value), quote=True) + '">' + html.escape(display) + '</dd>')
        parts.append('</dl>')
    parts.append('</details>')
    return ''.join(parts)


def render_event(result, selector, *, tables=None, spec=None, adata=None,
                 general_report=None, navigation_report=None, evidence_cache=None):
    from .reporting import build_report_tables, _explain, _panel_rows, _score_svg
    from ._report_evidence import expression_evidence_html, sample_scores_svg
    tables = build_report_tables(result) if tables is None else tables
    spec = resolve_positive(result) if spec is None else spec
    complete = {k: v for k, v in selector.items() if k != 'sample'}
    detail = _explain(result, tables, **complete)
    key = event_id(complete)
    title = f'{complete["sender"]} → {complete["hmdb_id"]} / {complete["sensor_gene"]} → {complete["receiver"]}'
    aggregate = detail['aggregate_event']
    mode = result.parameters.get('sample_mode')
    sample_aware = mode == 'sample_aware'
    sample_events = detail['sample_events']
    total_samples = sample_events.get('sample', pd.Series(dtype=object)).nunique()
    row = None
    mode_label = {'sample_aware': '按样本分析', 'pooled_stratified': '合并分析'}.get(mode, str(mode))
    parts = ['<section class="report-section" id="event-overview"><h2>事件概览</h2>',
             '<p class="event-path">' + html.escape(title) + '</p>',
             '<p>分析模式：' + html.escape(mode_label) + '；达标标准：' + html.escape(positive_description(spec)) + '</p>']
    if aggregate.empty:
        parts.append('<p>未生成该事件：本页仅提供诊断，不构造评分或统计量。</p>')
    else:
        from types import SimpleNamespace
        status, _ = positive_calls(SimpleNamespace(events=aggregate), spec)
        row = aggregate.iloc[0]
        score, statistic = row.get('cell_mesh_score'), row.get(spec['column'])
        parts += ['<div class="overview-cards event-stat-cards">',
                  metric_card(_stat_value(score), '通讯评分', primary=True,
                              note='样本评分中位数' if sample_aware else '合并分析评分',
                              title=f'原值：{score}' if pd.notna(score) else '原值缺失'),
                  metric_card(_stat_value(statistic), '本次' + STATISTIC_LABELS[spec['column']],
                              title=f'原值：{statistic}' if pd.notna(statistic) else '该统计量缺失',
                              status=status.iloc[0])]
        if sample_aware:
            computable = row.get('n_samples_coobserved')
            denominator = f'{total_samples:,}' if total_samples else '未记录'
            parts.append(metric_card(_stat_value(computable) + ' / ' + denominator, '可计算样本',
                                     note='占全部已有样本'))
        parts.append('</div>')
        parts.append(_event_decision_html(result, aggregate, spec, status.iloc[0]))
    if sample_aware:
        parts.append('<p class="report-note">统计量针对汇总事件；下方评分组成可按样本切换。</p>')
    parts += ['</section><section class="report-section" id="event-components"><h2>评分组成</h2>',
              '<p>每个分析单元的通讯评分 = √发送端评分 × √接收端评分。</p>',
              '<details id="event-score-method"><summary>评分计算说明</summary>',
              '<p>P、C、E 分别为生产、消耗和外排分数。先乘细胞丰度权重，再用加权后的值计算参考值，'
              '以正值计算参考值，按 x / (x + x_ref) 标准化；无正值时参考值记为 NA，标准化分数为 0。</p>',
              '<p>基础可用性由标准化后的生产和消耗分数组合，随后乘外排调节系数得到发送端评分；'
              '接收端评分由受体表达标准化并经表达比例门控得到。'
              '各项原值、参考值及异常状态保存在计算输入和数据明细中。</p></details>']
    samples = sample_events['sample'].tolist()
    if sample_aware:
        parts.append('<label class="sample-controls">查看样本 <select id="event-sample">')
        initial = selector.get('sample', samples[0] if samples else None)
        for i, sample in enumerate(samples):
            parts.append(f'<option value="{i}"' + (' selected' if sample == initial else '') + '>' + html.escape(str(sample)) + '</option>')
        parts.append('</select></label>')
    else:
        initial = None
    for i, sample in enumerate(samples):
        current = {}
        for name in ('sender_components', 'receiver_components', 'sample_events'):
            frame = detail[name]
            current[name] = frame if frame.empty else frame.loc[frame['sample'].isna() if sample is None else frame['sample'].eq(sample)]
        parts.append(f'<div class="sample-component" data-component-sample="{i}"' +
                     (' hidden' if sample_aware and sample != initial else '') + '>')
        parts.append('<h3>' + ('样本 ' + html.escape(str(sample)) if sample_aware else '合并分析组成') + '</h3>')
        parts.append('<div class="explanation-panels">')
        for title_, _, values in list(_panel_rows(current))[:3]:
            translated = [(COMPONENT_LABELS.get(label.rsplit(' / ', 1)[-1], label), value, color) for label, value, color in values]
            parts.append(_score_svg(translated, PANEL_TITLES[title_]))
        parts.append('</div>')
        parts.append(_component_inputs(current))
        reasons = [REASON_LABELS.get(reason, reason) for row in current['sample_events'].get('reasons', []) for reason in row]
        if reasons:
            parts.append('<p>状态说明：' + html.escape('；'.join(dict.fromkeys(reasons))) + '</p>')
        parts.append('</div>')
    parts.append('</section>')
    if sample_aware:
        parts += ['<section class="report-section" id="event-samples"><h2>样本表现</h2>',
                  _event_sample_summary(row), sample_scores_svg(sample_events, linked=True), '</section>']
    cache = {} if evidence_cache is None else evidence_cache
    cache_key = tuple(complete[k] for k in ('sender', 'receiver', 'hmdb_id', 'sensor_gene'))
    if cache_key not in cache:
        cache[cache_key] = expression_evidence_html(result, [complete], adata=adata)[0]
    parts += ['<section class="report-section" id="event-evidence"><h2>表达证据与先验覆盖</h2>', cache[cache_key],
              '<h3>先验覆盖</h3><p>覆盖度表示先验基因在输入中的可用情况，不是事件置信度。</p>',
              '<p id="event-table-help" class="report-note">表格已内嵌全部记录，可离线搜索、排序和翻页；'
              '操作仅调整表格视图。悬停字段标题可查看原字段名。</p>'
              '<noscript><p>当前未启用 JavaScript，各表仅预览前 50 条记录；启用后可浏览全部记录。</p></noscript>',
              full_table(detail['coverage'], 'coverage', TABLE_LABELS['coverage']), '</section>',
              '<section class="report-section" id="event-tables"><h2>数据明细与来源</h2>']
    for label, names in [('样本与评分', ['aggregate_event', 'sample_events', 'sender_components', 'receiver_components']),
                         ('反应与基因', ['reaction_details', 'gene_details']),
                         ('先验与标准化参考', ['enzyme_prior', 'interaction_prior', 'reference_context'])]:
        parts.append('<details class="event-data-group"><summary>' + label + '</summary>')
        for name in names:
            parts += [f'<details data-table="{name}"><summary>{TABLE_LABELS[name]}</summary>',
                      full_table(detail[name], name, TABLE_LABELS[name]), '</details>']
        parts.append('</details>')
    from . import __version__
    parts += ['<p class="event-source">来源：' +
              html.escape(result._result_path.name if result._result_path is not None else '当前内存结果') +
              f' · CELL MESH {html.escape(__version__)} · 报告布局 v3<br>事件标识：<span class="report-id">' +
              key + '</span></p></section>']
    nav = [('event-overview', '事件概览'), ('event-components', '评分组成')]
    if sample_aware:
        nav.append(('event-samples', '样本表现'))
    nav += [('event-evidence', '表达与覆盖'), ('event-tables', '数据与来源')]
    return document('CELL MESH 单事件报告', parts, scope='event', nav=nav,
                    general_report=general_report, navigation_report=navigation_report)


def render_navigation(*, general=None, selection=None, entries=(), pairs=(), notes=()):
    """One portable entry for independent report scopes."""
    parts = ['<p>通讯筛选报告与代谢物–sensor 通讯报告可独立生成。</p>']
    parts += ['<p class="report-note">' + html.escape(note) + '</p>' for note in notes]
    parts.append('<ul>')
    for path, label in [(general, '统计筛选报告'), (selection, '通讯筛选报告')]:
        if path:
            parts.append('<li><a href="' + html.escape(quote(str(path)), quote=True) + '">' + label + '</a></li>')
    parts += ['</ul><h2>代谢物–sensor 通讯报告</h2>']
    if pairs:
        parts.append('<ul>')
        parts += ['<li><a href="' + html.escape(quote(entry['href']), quote=True) + '">' +
                  html.escape(entry['label']) + '</a></li>' for entry in pairs]
        parts.append('</ul>')
    else:
        parts.append('<p>未指定代谢物–sensor 通讯报告。</p>')
    legacy = [entry for entry in entries if entry.get('scope') == '单事件报告']
    if legacy:
        parts.append('<details><summary>历史单事件报告（沿用）</summary><ul>')
        parts += ['<li><a href="' + html.escape(quote(entry['href']), quote=True) + '">' +
                  html.escape(entry['label']) + '</a></li>' for entry in legacy]
        parts.append('</ul></details>')
    return document('CELL MESH 报告导航', parts, scope='navigation')
