"""Embed one complete static communication matrix; no plotting JavaScript."""
from __future__ import annotations

import html
import io
import re


def figure_svg(fig, prefix, *, preserve_canvas=False):
    buffer = io.StringIO()
    # The matrix renderer measures its own canvas, including browser-font slack.
    # A subsequent tight crop would discard that margin around the BH legend.
    fig.savefig(buffer, format='svg', bbox_inches=fig.bbox_inches.frozen() if preserve_canvas else 'tight',
                metadata={'Date': None})
    content = buffer.getvalue()
    content = content[content.index('<svg'):]
    def namespace(match):
        tag = re.sub(r'\bid="([^"]+)"', lambda m: f'id="{prefix}-{m[1]}"', match[0])
        tag = re.sub(r'\bhref="#([^"]+)"', lambda m: f'href="#{prefix}-{m[1]}"', tag)
        return re.sub(r'url\(#([^)]+)\)', lambda m: f'url(#{prefix}-{m[1]})', tag)
    return re.sub(r'<[^>]+>', namespace, content)


def communication_dotmap_html(output, *, prefix='selection'):
    """Embed a public plotting result without selecting, encoding or drawing again.

    The caller owns the Figure and closes it after embedding/exporting it.
    """
    from ._communication_dotmap import statistic_caption, KEYS
    from ._report_browser import _script_json
    scope = output['scope']
    n_scope = len(scope['events'])
    positions = {key: i for i, key in enumerate(scope['frame'][KEYS].itertuples(index=False, name=None))}
    indices = [positions[key] for key in scope['events'][KEYS].itertuples(index=False, name=None)]
    record = dict(kind='dotmap', n_selected=scope['n_selected'], n_plotted=n_scope,
                  n_missing_score=int(scope['frame'].cell_mesh_score.isna().sum()),
                  n_missing_score_plotted=int(scope['events'].cell_mesh_score.isna().sum()),
                  n_display_excluded=scope['n_selected']-n_scope,
                  plotted_event_indices=sorted(indices), plotted_pairs=scope['pairs'],
                  plot_top_n=scope['top_n'], top_n_unit='hmdb_id + sensor_gene', n_figures=1 if n_scope else 0)
    if not n_scope:
        return '', record
    import matplotlib as mpl
    encoding = output['encoding']
    statistic_name, statistic_context = statistic_caption(scope['statistic'])
    statistic_label = html.escape(statistic_name)
    statistic_note = html.escape('（' + statistic_context + '）') if statistic_context else ''
    scope_text = f'所选 {scope["n_selected"]:,} 条 · 图中 {n_scope:,} 条'
    if record['n_display_excluded']:
        scope_text += f' · 超出组合展示上限 {record["n_display_excluded"]:,} 条'
    name_mode = output['labels']['metabolite_label'] == 'name'
    metabolite_caption = '名称' if name_mode else 'HMDB ID'
    parts = [f'<details class="selection-chart" id="{prefix}-dotmap" open>',
             '<summary>通讯事件矩阵点图' + f'<span class="chart-summary">{n_scope:,} 条事件</span></summary>',
             '<p class="plot-scope">' + scope_text + '</p>',
             '<p class="report-note">横轴为接收细胞—sensor，纵轴为发送细胞—代谢物（' + metabolite_caption + '）；每个点代表一个完整事件。'
             '点面积表示通讯评分，颜色表示 −log10(' + statistic_label + ')' + statistic_note + '。</p>']
    # Freeze the Figure's actual glyphs, as in a default standalone SVG export.
    # Browser font substitution or inherited page CSS must not restyle the plot.
    with mpl.rc_context({'svg.fonttype': 'path', 'svg.hashsalt': f'{prefix}-dotmap'}):
        svg = figure_svg(output['fig'], f'{prefix}-dotmap', preserve_canvas=True)
    svg = svg.replace('<svg ', '<svg role="img" aria-label="通讯事件矩阵点图" ', 1)
    parts.append(f'<div class="selection-plot dotmap-figure">{svg}</div>')
    notes = ['刻度使用单行复合标签；横轴倾斜 45°，纵轴水平。代谢物显示' + metabolite_caption + '，完整标识见事件列表。',
             '圆点保留最小可见面积；面积用于比较原评分，不表示实测通量。',
             '空白表示该位置没有本次展示的事件，不能据此认定不存在生物学通讯。']
    if encoding['missing_score']:
        notes.append('× 表示评分缺失，不作为零评分。')
    if encoding['missing_statistic']:
        notes.append('灰色表示统计值缺失。')
    if encoding['has_zero']:
        notes.append('统计值为 0 使用独立色块，不替换为人为设定的小数。')
    if scope['statistic'] == 'fdr_sensor_type':
        notes.append('BH 校正按 sensor_type 分别进行，请结合各类型的校正范围解读。')
    if scope['top_n'] is not None:
        notes.append('设定组合上限时，按最低统计值、最高评分选取；无有效评分的组合排在其后。')
    notes.append('矩阵按细胞类型和分子标识排列。图形不另加显著性或 QC 筛选；完整 selection 见下方列表与 CSV。')
    parts.append('<details class="dotmap-reading"><summary>读图说明与展示范围</summary>' +
                 ''.join('<p class="report-note">' + note + '</p>' for note in notes) + '</details>')

    data = dict(records=scope['records'], pairs=scope['pair_records'], encoding=encoding,
                n_selected=scope['n_selected'], n_plotted=n_scope,
                x_keys=output['x_keys'], y_keys=output['y_keys'],
                x_fields=['receiver', 'sensor_gene'], y_fields=['sender', 'hmdb_id'],
                dots=output['dot_table'].to_dict('records'), labels=output['labels'], layout=output['layout'])
    parts += [f'<script type="application/json" id="{prefix}-dotmap-data">', _script_json(data),
              '</script></details>']
    return ''.join(parts), record
