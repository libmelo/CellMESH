"""Offline reports for explicitly requested metabolite–sensor pairs."""
from __future__ import annotations

import html
import tempfile
from pathlib import Path

import pandas as pd

from ._pair_selection import pair_id
from .selection import EVENT_KEYS


def _pair_rows(frame, pair):
    return frame.loc[frame.hmdb_id.eq(pair[0]) & frame.sensor_gene.eq(pair[1])].copy()


def _pair_tables(result, tables, pair, frame):
    """Index existing explanation rows once; retain all samples, no rescoring."""
    senders, receivers = set(frame.sender), set(frame.receiver)
    explanation = {}
    for name, table in tables['explanation'].items():
        subset = table
        if 'hmdb_id' in subset:
            subset = subset.loc[subset.hmdb_id.eq(pair[0])]
        if 'sensor_gene' in subset:
            subset = subset.loc[subset.sensor_gene.eq(pair[1])]
        if 'sender' in subset:
            subset = subset.loc[subset.sender.isin(senders)]
        if 'receiver' in subset:
            subset = subset.loc[subset.receiver.isin(receivers)]
        if name == 'reaction_details':
            subset = subset.loc[subset.cell_type.isin(senders)]
        if name == 'reference_context':
            subset = subset.loc[((subset.feature_kind != 'receiver') & subset.feature_id.eq(pair[0])) |
                                ((subset.feature_kind == 'receiver') & subset.feature_id.eq(pair[1]))]
        explanation[name] = subset
    needed = {pair[1]}
    for genes in explanation['reaction_details'].get('genes', []):
        needed.update(genes)
    genes = explanation['gene_details']
    explanation['gene_details'] = genes.loc[genes.cell_type.isin(senders | receivers) & genes.gene.isin(needed)]
    cov = tables['coverage']['details']
    cov = cov.loc[cov.hmdb_id.eq(pair[0]) & cov.cell_type.isin(senders)]
    return dict(explanation=explanation, coverage=dict(details=cov))


def render_pair(result, selection, *, tables, spec, plot_top_n=20, adata=None,
                csv_path=None, general_report=None, navigation_report=None):
    from ._report_pages import (document, full_table, sample_heatmap, _component_inputs,
                                _stat_value, _event_decision_html)
    from ._report_browser import event_browser, event_id, _script_json
    from ._report_layout import (selection_overview_html, filter_details_html, PANEL_TITLES,
                                 COMPONENT_LABELS, TABLE_LABELS, REASON_LABELS)
    from ._report_plots import selection_plots_html
    from .plotting import plot_communication_dotmap
    import matplotlib.pyplot as plt
    from ._selection_evidence import selection_evidence_html
    from .positivity import selection_table, STATUS_LABELS, positive_calls
    from .reporting import _explain, _panel_rows, _score_svg, _table
    from urllib.parse import quote
    pair = selection.pairs[0]
    frame = selection_table(selection, spec)
    key = pair_id(pair)
    generated = _pair_rows(result.events, pair)
    dotmap = plot_communication_dotmap(frame, top_n=None, fdr_col=spec['column'])
    scope = dotmap['scope']
    try:
        plots = selection_plots_html(selection, None, prefix='pair', scope_label='当前明确配对的所选事件',
                                     fdr_col=spec['column'], dotmap=dotmap, secondary=True)
    finally:
        plt.close(dotmap['fig'])
    sections = ['<section class="report-section" id="pair-overview"><h2>分子对概览</h2>',
                '<p class="event-path">' + html.escape(' / '.join(pair)) + '</p>']
    names = generated.get('metabolite', pd.Series(dtype=str)).dropna().astype(str).unique()
    if len(names):
        sections.append('<p>代谢物名称：' + html.escape('、'.join(names)) + '</p>')
    sections += [selection_overview_html(result, selection, frame, spec),
                 '<p>本页固定一个代谢物–sensor 组合，评分和统计量仍属于各条细胞通讯事件。</p>']
    if frame.empty:
        message = ('原结果未生成该分子对的通讯事件。' if generated.empty else
                   f'原结果中该分子对有 {len(generated):,} 条事件，但均未满足本次来源或附加筛选条件。')
        sections.append('<p class="pair-diagnostic">' + message + '未构造评分或统计量。</p>')
    sections.append('</section>')
    entries = [dict(id=event_id(row), label=f'{row["sender"]} → {row["receiver"]}',
                    scope='细胞对详情', identity=tuple(row[k] for k in EVENT_KEYS))
               for row in frame.to_dict('records')]
    # EVENT_KEYS and browser identities have the same endpoint ordering.
    shown = frame.copy()
    shown['positive_status'] = shown.positive_status.map(STATUS_LABELS)
    sections += ['<section class="report-section" id="pair-communication"><h2>通讯分布</h2>',
                 plots,
                 event_browser(shown, entries, prefix='pair-events', title='完整通讯事件列表'), '</section>']
    sections += ['<section class="report-section" id="pair-evidence"><h2>评分与表达依据</h2>',
                 '<p>发送端与接收端证据按细胞类型组织，同一输入不按通讯事件重复绘制。</p>',
                 '<details><summary>评分计算说明</summary><p>每个分析单元的通讯评分 = √发送端评分 × √接收端评分。'
                 'P、C、E 分数先乘细胞丰度权重，再以加权后的正值计算参考值，按 x / (x + x_ref) 标准化。'
                 '发送端包含生产、消耗与外排调节；接收端包含 sensor 表达标准化与表达比例门控。'
                 'sample_aware 的汇总评分为事件样本评分的中位数，不由各组成项的中位数重新计算。</p></details>',
                 selection_evidence_html(result, scope, adata=adata)]
    scoped = _pair_tables(result, tables, pair, frame)
    explanation = scoped['explanation']
    for name in ['sender_components', 'receiver_components']:
        sections += ['<details><summary>' + TABLE_LABELS[name] + '</summary>',
                     full_table(explanation[name], 'pair-' + name, TABLE_LABELS[name]), '</details>']
    sections.append('</section>')
    sample_aware = result.parameters.get('sample_mode') == 'sample_aware'
    if sample_aware:
        sections += ['<section class="report-section" id="pair-samples"><h2>样本表现</h2>',
                     '<p>热图包含全部已有样本；事件列的展示上限不影响事件详情或 CSV。</p>',
                     sample_heatmap(selection, plot_top_n), '</section>']
    sections += ['<section class="report-section" id="pair-details"><h2>单事件详情</h2>',
                 '<p>展开细胞对查看判定、全部样本的计算输入和评分分解。共同的表达与先验见本页相应部分。</p>']
    observed = _pair_rows(result.sample_events, pair) if result.sample_events is not None else generated
    for row in frame.to_dict('records'):
        identity = {name: row[name] for name in EVENT_KEYS}
        event_key = event_id(identity)
        detail = _explain(result, scoped, **identity, observed=observed, aggregate_events=generated)
        status = row['positive_status']
        sections += [f'<details class="event-explanation" id="{event_key}"><summary>' +
                     html.escape(f'{row["sender"]} → {row["receiver"]}') +
                     f' · 评分 {_stat_value(row.get("cell_mesh_score"))} · {STATUS_LABELS[status]}</summary>',
                     _event_decision_html(result, detail['aggregate_event'], spec, status),
                     full_table(detail['sample_events'], event_key + '-samples', '样本评分与诊断')]
        for i, sample_row in enumerate(detail['sample_events'].to_dict('records')):
            sample = sample_row['sample']
            current = {}
            for name in ('sender_components', 'receiver_components', 'sample_events'):
                table = detail[name]
                current[name] = table.loc[table['sample'].isna() if pd.isna(sample) else table['sample'].eq(sample)]
            sections.append('<details class="pair-sample-detail"><summary>' +
                            html.escape('合并分析组成' if pd.isna(sample) else '样本 ' + str(sample)) + '</summary>')
            sections.append('<div class="explanation-panels">')
            for title, _, values in list(_panel_rows(current))[:3]:
                translated = [(COMPONENT_LABELS.get(label.rsplit(' / ', 1)[-1], label), value, color)
                              for label, value, color in values]
                sections.append(_score_svg(translated, PANEL_TITLES[title]))
            sections += ['</div>', _component_inputs(current)]
            reasons = sample_row.get('reasons', [])
            if reasons:
                sections.append('<p>状态说明：' + html.escape('；'.join(REASON_LABELS.get(r, r) for r in reasons)) + '</p>')
            sections.append('</details>')
        sections.append('</details>')
    sections += ['</section><section class="report-section" id="pair-sources"><h2>数据来源与方法</h2>',
                 '<p>共同明细保留全部已有样本；按发送细胞、接收细胞或 cell_type 字段定位对应事件的输入。'
                 '先验覆盖表示输入可用情况，不是事件置信度。</p>']
    enzyme = result.priors['enzyme']
    interaction = result.priors['interaction']
    sources = dict(coverage=scoped['coverage']['details'],
                   **{name: explanation[name] for name in ['reaction_details', 'gene_details', 'reference_context']},
                   enzyme_prior=enzyme.loc[enzyme.hmdb_id.eq(pair[0])],
                   interaction_prior=_pair_rows(interaction, pair))
    for name, table in sources.items():
        sections += ['<details><summary>' + TABLE_LABELS[name] + '</summary>',
                     full_table(table, 'pair-' + name, TABLE_LABELS[name]), '</details>']
    if csv_path:
        sections.append('<p><a href="' + html.escape(quote(str(csv_path)), quote=True) + '">下载指定代谢物–sensor 事件 CSV</a>'
                        '：该表包含本次交付所有请求分子对的所选事件，可按 pair_id、hmdb_id 和 sensor_gene 筛选。</p>')
    from . import __version__
    source_name = result._result_path.name if result._result_path is not None else '当前内存结果'
    sections += [filter_details_html(selection, _table),
                 '<p class="report-note">表格完整内嵌，可离线搜索和分页；未启用 JavaScript 时仅预览前 50 条。'
                 '图形为静态图，原始统计量及其校正范围保持不变。</p>',
                 '<p class="event-source">来源：' + html.escape(source_name) +
                 ' · CELL MESH ' + html.escape(__version__) + ' · 报告布局 v4</p>',
                 '<p class="report-id">' + key + '</p>',
                 '<script type="application/json" id="pair-report-data">' + _script_json(dict(
                     pair_id=key, pair=list(pair), criteria=selection.criteria,
                     event_keys=frame[EVENT_KEYS].to_dict('records'), n_generated=len(generated),
                     source='captured_scientific_results')) + '</script></section>']
    nav = [('pair-overview', '分子对概览'), ('pair-communication', '通讯分布'), ('pair-evidence', '评分与表达依据')]
    if sample_aware:
        nav.append(('pair-samples', '样本表现'))
    nav += [('pair-details', '单事件详情'), ('pair-sources', '数据来源与方法')]
    return document('CELL MESH 代谢物–sensor 通讯报告', sections, scope='metabolite_sensor', nav=nav,
                    general_report=general_report, navigation_report=navigation_report)


def write_pair_bundle(result, path, selected, *, plot_top_n=20, adata=None, **positive_options):
    from ._pair_selection import pair_records
    from .positivity import resolve_positive, delivery_record, INHERIT
    from ._report_delivery import render_bundle
    from ._report_pages import render_navigation
    from ._report_plots import validate_plot_top_n
    from .result_io import (_validate_targets, protected_bundle_paths, _staged, _publish,
                            _write_delivery_csvs)
    from .reporting import build_report_tables
    path = Path(path)
    if path.suffix.lower() != '.html':
        raise ValueError('Report destination must have an .html suffix')
    plot_top_n = validate_plot_top_n(plot_top_n)
    spec = resolve_positive(result, selection=selected, **positive_options)
    delivery = delivery_record(result, spec, None, [], INHERIT, plot_top_n)
    delivery.update(report_layout=4, pair_selections=pair_records(selected, spec))
    directory = path.parent / (path.stem + '.reports')
    paths = dict(navigation=path, selected_csv=path.with_name(path.stem + '.metabolite_sensor_events.csv'))
    for record in delivery['pair_selections']:
        paths['pair_report:' + record['id']] = directory / 'metabolite_sensor' / (record['id'] + '.html')
    _validate_targets(paths.values())
    if result._result_path is not None:
        protected, report_dir = protected_bundle_paths(result)
        if any(p.resolve() in protected or p.resolve().is_relative_to(report_dir) for p in paths.values()):
            raise ValueError('Standalone export cannot overwrite managed delivery files; use add_event_reports()')
    contents = render_bundle(result, paths, delivery, None, tables=build_report_tables(result), adata=adata)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.cellmesh-pairs-', dir=path.parent) as temporary:
        staging = Path(temporary)
        _write_delivery_csvs(result, delivery, None, paths, staging)
        for key, content in contents.items():
            dest = _staged(staging, paths[key])
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(content, encoding='utf-8')
        _publish(staging, paths, list(paths.values()))
    return dict(paths, report=path)
