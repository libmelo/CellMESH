"""Compact report presentation; scientific tables and values stay unchanged."""
from __future__ import annotations

import html
import json

import pandas as pd


LAYOUT_CSS = """
.report-nav{display:flex;flex-wrap:wrap;gap:20px;padding:12px 0;border-bottom:1px solid #dde4e8}
.report-nav a{text-decoration:none;font-weight:600}
.report-section{margin:28px 0;scroll-margin-top:16px}
.report-meta,.report-note{color:#526773}.report-meta{margin:8px 0}
.overview-cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(240px,100%),1fr));gap:16px;margin:18px 0;max-width:960px}
.selection-overview-cards{max-width:440px}
.overview-card{background:#f8fafb;border:1px solid #dde8e8;border-radius:6px;padding:18px;min-width:0}
.overview-card-primary{background:#eef7f7;border-left:4px solid #165d69}
.overview-card strong{display:block;font-size:26px;color:#165d69;font-variant-numeric:tabular-nums;overflow-wrap:anywhere}
.overview-card-primary strong{font-size:34px}
.overview-card span{font-size:13px;color:#526773}
.overview-card .metric-label{display:block;font-size:14px;font-weight:600;margin-bottom:6px}
.overview-card .metric-note{display:block;margin-top:6px}
.overview-card .status-badge{display:inline-block;margin-top:8px}
.status-badge{font-size:13px;border-radius:4px;padding:2px 8px;background:#eef1f3;color:#526773}
.status-badge.positive{background:#e2f1ed;color:#24614d}
.status-badge.not_evaluated,.overview-status .status-missing{background:#fff2db;color:#795313}
.overview-notice{padding:10px 14px;border-left:3px solid #bb9146;background:#fff8eb}
.overview-status,.event-sample-summary{margin:12px 0;color:#526773}
.overview-status .status-missing{padding:2px 6px;border-radius:3px}
.overview-counts{font-size:13px;color:#526773}
.selection-chart{border:1px solid #dde4e8;border-radius:6px;padding:12px;margin:14px 0}
.selection-chart>summary{font-weight:600;color:#165d69}.chart-summary{font-weight:400;color:#526773;margin-left:12px}
.selection-plot{margin-top:12px}
.selection-plot svg{display:block;width:auto;max-width:100%;height:auto;margin:auto}
.plot-scope{font-size:13px;color:#526773}
.explanation-panels{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;align-items:start}
.explanation-panels svg{width:100%;max-width:none;height:auto}
.event-data-group{border-top:1px solid #dde4e8;padding-top:10px}.event-data-group>summary{font-weight:600}
.event-data-group details{margin-left:12px}.report-help li{margin:6px 0}
@media(max-width:760px){.explanation-panels{grid-template-columns:1fr}}
"""


PANEL_TITLES = {
    "P/C/E normalized scores": "P/C/E 标准化评分",
    "Sender: base and exporter adjustment": "发送端：基础可用性与外排调节",
    "Receiver: before / after expression gate": "接收端：表达比例门控前后",
    "Event scores by sample": "各样本的通讯评分",
}
COMPONENT_LABELS = {
    "P_score": "生产评分", "C_score": "消耗评分", "E_score": "外排评分",
    "base_availability": "基础可用性", "E_effective": "外排有效值",
    "E_factor": "外排调节系数", "sender_score": "发送端评分",
    "score_before_gate": "门控前评分", "sensor_score": "门控后评分",
    "cell_mesh_score": "通讯评分",
}
TABLE_LABELS = {
    "aggregate_event": "事件汇总", "sample_events": "各样本事件记录",
    "sender_components": "发送端评分组成", "receiver_components": "接收端评分组成",
    "reference_context": "标准化参考值", "coverage": "输入基因覆盖度",
    "reaction_details": "反应明细", "gene_details": "基因明细",
    "enzyme_prior": "酶—代谢物先验", "interaction_prior": "代谢物—受体先验",
}
REASON_LABELS = {
    "interaction_prior_missing": "缺少代谢物—受体先验",
    "self_communication_disabled": "本次分析未启用自通讯",
    "sender_unobserved": "未观察到发送细胞类型", "receiver_unobserved": "未观察到接收细胞类型",
    "sensor_gene_unavailable": "输入中缺少受体基因", "production_prior_missing": "缺少生产反应先验",
    "production_genes_unavailable": "输入中缺少生产反应基因",
    "receiver_expression_fraction_gate": "接收端表达比例未达到门槛",
    "receiver_measured_zero": "接收端平均表达为零", "zero_production_capacity": "生产能力为零",
    "zero_score_see_components_and_numerical_diagnostics": "评分为零，需结合组成与数值诊断查看",
    "event_not_generated": "未生成该事件",
}


def _json(value):
    return '<pre>' + html.escape(json.dumps(value, ensure_ascii=False, indent=2, default=str)) + '</pre>'


def selection_description(criteria):
    labels = {"sender_labels": "发送细胞", "receiver_labels": "接收细胞", "cell_pairs": "有向细胞对",
              "hmdb_ids": "代谢物", "sensor_genes": "受体"}
    parts = ["来源：一级达标事件" if criteria.get("source") == "positive" else "来源：全量事件"]
    for key, label in labels.items():
        values = criteria[key]
        if values is None:
            continue
        text = [f"{value[0]} → {value[1]}" if key == "cell_pairs" else str(value) for value in values]
        shown = "、".join(text[:5]) if text else "空列表"
        if len(text) > 5:
            shown += f"等 {len(text)} 项"
        parts.append(f"{label}：{shown}")
    if criteria["max_fdr"] is not None:
        scope = "BH 校正 p 值（全局）" if criteria["fdr_col"] == "fdr_global" else "BH 校正 p 值（按受体类型）"
        parts.append(f'{scope} ≤ {criteria["max_fdr"]:g}')
    if criteria.get("max_perm_pvalue") is not None:
        parts.append(f'置换 p 值 ≤ {criteria["max_perm_pvalue"]:g}')
    if criteria["min_cell_mesh_score"] is not None:
        parts.append(f'通讯评分 ≥ {criteria["min_cell_mesh_score"]:g}')
    if criteria["qc_only"]:
        parts.append("通过细胞数质控")
    return "；".join(parts) if parts else "未设置筛选条件"


def percent_text(value):
    percent = float(value) * 100
    if 0 < percent < 0.01:
        return '<0.01%'
    return f'{percent:.2f}'.rstrip('0').rstrip('.') + '%'


def status_badge(status, text=None):
    labels = {'positive': '达标', 'not_positive': '未达标', 'not_evaluated': '无法判定'}
    return (f'<span class="status-badge {html.escape(status, quote=True)}">' +
            html.escape(labels[status] if text is None else text) + '</span>')


def metric_card(value, label, *, note=None, title=None, primary=False, status=None, status_text=None):
    """Escape plain display values; original scientific precision stays in tables."""
    css = 'overview-card overview-card-primary' if primary else 'overview-card'
    tooltip = '' if title is None else ' title="' + html.escape(str(title), quote=True) + '"'
    parts = [f'<div class="{css}"><span class="metric-label">', html.escape(label),
             '</span><strong' + tooltip + '>', html.escape(str(value)), '</strong>']
    if note:
        parts += ['<span class="metric-note">', html.escape(note), '</span>']
    if status is not None:
        parts.append(status_badge(status, status_text))
    return ''.join(parts) + '</div>'


def status_counts_html(counts, *, selected=False):
    """Secondary counts never compete with the scope's main metric."""
    labels = [('n_positive', '达标'), ('n_not_positive', '未达标'), ('n_not_evaluated', '无法判定')]
    all_selected_pass = selected and counts['n_total'] == counts['n_positive']
    visible = []
    if not all_selected_pass:
        for key, label in labels:
            if counts[key] and (selected or key != 'n_positive'):
                css = ' class="status-missing"' if key == 'n_not_evaluated' else ''
                visible.append(f'<span{css}>{label} {counts[key]:,}</span>')
    parts = ['<p class="overview-status">' + ' · '.join(visible) + '</p>'] if visible else []
    if not selected:
        complete = [('n_total', '全部已生成事件'), *labels]
        parts.append('<details class="overview-counts"><summary>完整判定计数</summary><p>' +
                     '；'.join(f'{label}：{counts[key]:,}' for key, label in complete) + '</p></details>')
    return ''.join(parts)


def missing_criteria_html(result, frame, spec):
    """Explain missing inputs only for criteria that actually participate."""
    from .positivity import STATISTIC_LABELS
    reasons = []
    statistic = frame.get(spec['column'], pd.Series(index=frame.index, dtype=float))
    if statistic.isna().any():
        if statistic.isna().all() and result.parameters.get('n_perms') == 0:
            reasons.append('未执行置换检验（n_perms=0），所用 p 值或 BH 校正值缺失')
        else:
            reasons.append(STATISTIC_LABELS[spec['column']] + '缺失')
    for enabled, column, label in [(spec['qc_only'], 'passes_min_cells', '细胞数 QC 状态'),
                                   (spec['min_score'] is not None, 'cell_mesh_score', '通讯评分')]:
        if enabled and frame.get(column, pd.Series(index=frame.index, dtype=float)).isna().any():
            reasons.append(label + '缺失')
    if not reasons:
        return ''
    return '<p class="overview-notice">无法判定的原因：' + html.escape('；'.join(reasons)) + '。</p>'


def event_scope_html(frame, *, label, criteria=None):
    if frame.empty:
        return ''
    parts = [f'{label}涉及 {frame.sender.nunique()} 类发送细胞、{frame.receiver.nunique()} 类接收细胞，'
             f'共 {len(frame[["sender", "receiver"]].drop_duplicates()):,} 个有向细胞对']
    # A named singleton is already visible in the selection condition line.
    for key, column, name in [('hmdb_ids', 'hmdb_id', '代谢物'), ('sensor_genes', 'sensor_gene', '受体')]:
        values = (criteria or {}).get(key)
        if values is None or len(values) != 1:
            parts.append(f'{frame[column].nunique():,} 个{name}')
    return '<p class="overview-scope">' + '；'.join(parts) + '。</p>'


def overview_html(result, frame, spec, counts):
    from .positivity import positive_description
    mode = {"sample_aware": "按样本分析（sample-aware）", "pooled_stratified": "合并分析（pooled）"}.get(
        result.parameters.get("sample_mode"), str(result.parameters.get("sample_mode", "未记录")))
    parts = [f'<p class="report-meta">{html.escape(mode)}</p>']
    provenance = result.input_provenance or {}
    dimensions = []
    if 'n_obs' in provenance and 'n_vars' in provenance:
        dimensions.append(f'{provenance["n_obs"]:,} 个细胞 × {provenance["n_vars"]:,} 个基因')
    validation = result.sample_validation
    if validation is not None and 'sample' in validation:
        dimensions.append(f'{validation["sample"].nunique()} 个样本')
    if dimensions:
        parts.append('<p>输入规模：' + '；'.join(dimensions) + '。</p>')
    parts.append('<p><strong>达标标准：</strong>' + html.escape(positive_description(spec)) + '</p>')
    total, positive = counts['n_total'], counts['n_positive']
    undetermined = total > 0 and counts['n_not_evaluated'] == total
    ratio = ('占全部已生成事件的 ' + percent_text(positive / total)) if total and not undetermined else None
    parts += ['<div class="overview-cards">',
              metric_card('尚未判定' if undetermined else f'{positive:,}', '达标事件', note=ratio, primary=True),
              metric_card(f'{total:,}', '全部已生成事件'), '</div>', status_counts_html(counts)]
    if counts['n_not_evaluated']:
        parts.append(missing_criteria_html(result, result.events, spec))
    if counts['status'] == 'not_evaluated':
        parts.append('<p>尚未完成达标判定：缺少必要判定信息；空达标表不表示没有通讯。</p>')
    elif counts['status'] == 'empty':
        parts.append('<p>本结果没有通讯事件。</p>')
    elif not counts['n_positive']:
        parts.append('<p>已完成判定的事件中，没有事件达到本次标准。</p>')
    else:
        parts.append(event_scope_html(frame, label='达标事件'))
    parts.append('<p class="report-note">全量为已生成的通讯事件。达标（阳性）表示满足本次判定标准，不代表实验验证。'
                 '通用图表和主列表展示达标事件，不随 selection 改变。</p>')
    return ''.join(parts)


def selection_overview_html(result, selection, frame, spec):
    from .positivity import positive_description
    values = frame.positive_status.value_counts()
    counts = dict(n_total=len(frame), n_positive=int(values.get('positive', 0)),
                  n_not_positive=int(values.get('not_positive', 0)),
                  n_not_evaluated=int(values.get('not_evaluated', 0)))
    all_pass = len(frame) > 0 and counts['n_positive'] == len(frame)
    parts = ['<p>' + html.escape(selection_description(selection.criteria)) + '</p>',
             '<p>达标标准：' + html.escape(positive_description(spec)) + '</p>',
             '<div class="overview-cards selection-overview-cards">',
             metric_card(f'{len(frame):,}', '所选事件', primary=True,
                         status='positive' if all_pass else None, status_text='全部达标'), '</div>']
    if frame.empty:
        parts.append('<p>未选中事件，请结合来源范围和筛选条件查看。</p>')
    else:
        parts += [status_counts_html(counts, selected=True),
                  event_scope_html(frame, label='所选事件', criteria=selection.criteria)]
        if counts['n_not_evaluated']:
            parts.append(missing_criteria_html(result, frame, spec))
    return ''.join(parts)


def csv_links(links, names):
    from urllib.parse import quote
    labels = {'events_csv': '全部事件 CSV', 'positive_csv': '阳性事件 CSV', 'filtered_csv': '筛选事件 CSV', 'selected_csv': '指定代谢物–sensor 事件 CSV'}
    return '<p class="csv-downloads">' + ' · '.join(
        '<a download href="' + html.escape(quote(str(links[key]))) + '">' + labels[key] + '</a>'
        for key in names if key in (links or {})) + '</p>' if links else ''


def filter_details_html(selection, table):
    if selection is None:
        return ''
    summary = selection.summary
    order = ('统计值升序、评分降序、事件标识升序' if summary['fdr_available']
             else '评分降序、事件标识升序')
    return ('<details id="selection-audit"><summary>筛选明细与原始条件</summary>'
            f'<p>原始事件：{summary["n_total"]:,}；入选事件：{summary["n_selected"]:,}；'
            f'排除事件：{summary["n_excluded"]:,}。</p>'
            '<p>不同条件取交集，同一条件内取并集。下表按筛选顺序统计排除数，不重复计数。'
            '筛选沿用原始评分与 FDR。</p>' + _json(selection.criteria) +
            table(pd.DataFrame(summary['steps'])) + f'<p>初始排序：{order}；缺失值排在最后。</p></details>')


def methods_html(result, tables, selection, input_html, table, delivery=None, *, selected=False):
    delivery_display = dict(delivery or {})
    # General HTML must not change when an interest selection is appended.
    delivery_display = {key: value for key, value in delivery_display.items()
                        if key in (('positive', 'plot_top_n') if selected else ('positive',))}
    title = '筛选与导出' if selected else '数据与方法'
    parts = [f'<section class="report-section" id="report-methods"><h2>{title}</h2>',
             '<details id="report-reading"><summary>如何阅读评分与缺失值</summary><ul class="report-help">',
             '<li>NA 表示该项无法计算或未提供；0 表示已计算且数值为零。零表达不能证明生物学功能不存在。</li>',
             '<li>输入基因覆盖度描述先验基因在输入矩阵中的可用情况，不参与评分，也不代表反应功能或置信度。</li>',
             '<li>P、C、E 分别表示生产、消耗和外排。缺少外排先验时使用的中性有效值，不代表实测外排证据。</li>']
    if result.parameters.get('sample_mode') == 'sample_aware':
        parts += ['<li>跨样本通讯评分为各样本事件评分的中位数，不由组成项的中位数重新计算。</li>',
                  '<li>正评分样本比例描述事件在样本中的出现情况，不是显著率或重复验证率。</li>']
    parts += ['</ul></details>', '<details id="positive-method"><summary>阳性判定与交付设置</summary>'
              '<p>BH 列保存调整后的 p 值，用于控制错误发现率；原始置换 p 值未进行多重检验校正。'
              '筛选和追加不重新计算 p 值或 BH 校正值；缺失统计量不视为阴性。</p>' + _json(delivery_display) + '</details>',
              filter_details_html(selection, table)]
    if not selected:
        parts += ['<details id="report-coverage"><summary>全量输入覆盖度（不随筛选改变）</summary>',
              '<p>统计单元为“样本 × 已观察细胞类型 × 先验代谢物”。基因覆盖率的分母为保留反应的先验基因数；'
              '未提供相应先验时记为 NA。</p>', table(tables['coverage']['summary']),
              '<details><summary>P/C/E 覆盖明细</summary>' + table(tables['coverage']['details']) + '</details>',
              '<details><summary>接收端表达与基因覆盖</summary>' + table(tables['explanation']['receiver_components']) + '</details></details>',
              '<details id="report-input"><summary>输入数据与归档</summary>' + input_html + '</details>',
              '<details id="report-parameters"><summary>运行参数与先验来源</summary>' + _json(result.parameters) + '</details>']
    parts += ['<details id="report-export"><summary>完整数据与导出</summary>',
              ('<p>列表包含全部所选事件。' if selected else '<p>列表包含全部阳性事件。') + '一级全局明细表最多预览 100 行；单事件明细包含全部记录并支持分页；'
              '完整明细可通过 <code>coverage_report()</code> 或 <code>explain_event()</code> 查询。</p>',
              '<p><code>selection.to_csv(path)</code> 导出全部感兴趣事件及阳性状态。'
              '<code>plot_top_n</code> 只限制比较图的展示范围，不截断事件列表或 CSV。'
              '单事件报告由完整事件标识独立指定，可通过 <code>add_event_reports()</code> 追加。</p>',
              '<p>本报告可离线查看。折叠仅缩短页面，不减少内嵌数据或文件大小。</p></details></section>']
    content = ''.join(parts)
    if selected:
        content = content.replace('selection.to_csv(path)', 'filtered.to_csv(path)').replace(
            '单事件报告由完整事件标识独立指定，可通过 <code>add_event_reports()</code> 追加。',
            '代谢物–sensor 通讯报告由 <code>select_events(metabolite_sensor_pairs=...)</code> 明确指定，'
            '可通过 <code>add_event_reports()</code> 追加。')
    return content
