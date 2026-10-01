"""Bounded, shared expression evidence and complete sample-score presentation."""
from __future__ import annotations

import html
import io
import re

import pandas as pd


def sample_scores_svg(frame, *, linked=False):
    """All captured samples, with missingness and QC independent of score."""
    if frame.empty:
        return '<p>无样本记录。</p>'
    scores = pd.to_numeric(frame.cell_mesh_score, errors='raise')
    median = scores.median() if scores.notna().any() else None
    height = 90 + 29 * len(frame)
    parts = [f'<svg class="sample-score-plot" xmlns="http://www.w3.org/2000/svg" '
             f'role="img" aria-label="各样本的通讯评分" viewBox="0 0 920 {height}">',
             '<text x="5" y="22">全部样本；NA 为不可计算，0 为已计算且评分为零</text>']
    if median is not None:
        x = 360 + 320 * float(median)
        parts.append(f'<line x1="{x}" x2="{x}" y1="38" y2="{height-25}" '
                     'stroke="#c16b3c" stroke-dasharray="4 3"/>')
    for i, (_, row) in enumerate(frame.iterrows()):
        y = 58 + 29 * i
        score = scores.iloc[i]
        qc = row.get('passes_min_cells', None)
        status = '质控未记录' if pd.isna(qc) else '质控通过' if qc else '质控未通过'
        link = f'data-sample-index="{i}" tabindex="0" role="button"' if linked else ''
        parts.append(f'<a {link}><text x="5" y="{y}" font-size="13">{html.escape(str(row["sample"]))}</text></a>')
        if pd.notna(score):
            parts.append(f'<rect x="360" y="{y-13}" width="{320*float(score)}" height="16" fill="#277b89"/>')
        label = 'NA' if pd.isna(score) else f'{score:.6g}'
        parts.append(f'<text x="690" y="{y}" font-size="13">{label} · {status}</text>')
    parts.append('</svg>')
    note = '所有样本均不可计算。' if median is None else f'虚线：非缺失评分的中位数（{median:.6g}）。'
    parts.append('<p class="report-note">' + note + ('点击样本名称可切换评分组成。' if linked else '') + '</p>')
    return ''.join(parts)


def _figure_svg(draw, namespace, **options):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    before = set(plt.get_fignums())
    try:
        with mpl.rc_context({'text.usetex': False, 'text.parse_math': False,
                             'svg.fonttype': 'none', 'svg.hashsalt': namespace}):
            output = draw(**options)
            # The report provides the accurate evidence title and scope.
            output['ax'].set_title('')
            output['ax'].set_ylim(bottom=0)
            if output['color_range'] == (0., 0.):
                # A constant zero mean is not a tiny negative/positive range.
                # All colored artists represent zero; median markers stay white.
                from matplotlib.colors import Normalize
                mappable = output['colorbar'].mappable
                mappable.set_norm(Normalize(0, 1))
                color = mappable.cmap(0.)
                for artist in output['violin_bodies']:
                    artist.set_facecolor(color)
                for artist in output['constant_artists'] + output['fallback_artists']:
                    artist.set_color(color)
            buffer = io.StringIO()
            output['fig'].savefig(buffer, format='svg', bbox_inches='tight', metadata={'Date': None})
            value = buffer.getvalue()
            value = value[value.index('<svg'):]
            def tag(match):
                text = re.sub(r'\bid="([^"]+)"', lambda m: f'id="{namespace}-{m[1]}"', match[0])
                text = re.sub(r'\bhref="#([^"]+)"', lambda m: f'href="#{namespace}-{m[1]}"', text)
                return re.sub(r'url\(#([^)]+)\)', lambda m: f'url(#{namespace}-{m[1]})', text)
            return re.sub(r'<[^>]+>', tag, value)
    finally:
        for number in set(plt.get_fignums()) - before:
            plt.close(number)


def expression_evidence_html(result, selectors, *, adata=None):
    """One plot per metabolite/context across the explained cell types.

    Never open an input archive implicitly. An explicit AnnData takes precedence
    over the in-memory input retained at computation time.
    """
    if not selectors:
        return '', {}
    adata = result._input_adata if adata is None else adata
    if adata is None:
        return ('<div id="expression-evidence"><p class="report-note">缺少细胞级表达输入，无法展示表达分布。</p>'
                '<details><summary>补充表达输入</summary><p>使用 adata 参数提供原输入，或传入 '
                'result.load_input() 的返回值后重新生成报告。</p></details></div>', {})
    from .plotting import plot_metabolite_secretion_violin, plot_receptor_expression_violin
    from ._report_browser import _identity
    senders, receivers = {}, {}
    for selector in selectors:
        sender, receiver, hmdb, receptor = _identity(selector)
        senders.setdefault(hmdb, set()).add(sender)
        receivers.setdefault((hmdb, receptor), set()).add(receiver)
    charts, links = [], {}
    for kind, groups, draw in [('sender', senders, plot_metabolite_secretion_violin),
                               ('receiver', receivers, plot_receptor_expression_violin)]:
        for key, cells in groups.items():
            index = len(charts) + 1
            chart_id = f'expression-evidence-{index}'
            hmdb = key if kind == 'sender' else key[0]
            title = ('发送端生产相关表达' if kind == 'sender' else '接收端受体表达')
            title += ' · ' + hmdb + ('' if kind == 'sender' else ' / ' + key[1])
            options = dict(result=result, adata=adata, hmdb_id=hmdb)
            options.update({'sender_labels': sorted(cells), 'vmin': 0, 'vmax': 1} if kind == 'sender' else
                           {'receiver_labels': sorted(cells), 'receptor_gene': key[1]})
            try:
                chart = '<div class="selection-plot">' + _figure_svg(draw, chart_id, **options) + '</div>'
            except (KeyError, ValueError) as exc:
                # Unscored candidates can legitimately lack a measured gene or score.
                chart = '<p class="report-note">该表达图不可用：' + html.escape(str(exc)) + '</p>'
            caption = ('形状表示生产反应相关表达分布，颜色表示代谢物可用性评分。' if kind == 'sender' else
                       '形状表示受体表达分布，颜色表示平均表达。')
            charts.append(f'<details class="expression-evidence" id="{chart_id}"><summary>' +
                          html.escape(title) + '</summary><p>' + html.escape('、'.join(sorted(cells))) +
                          '：' + caption + '</p>' + chart + '</details>')
            links[(kind, key)] = chart_id
    event_links = {}
    for selector in selectors:
        identity = _identity(selector)
        _, _, hmdb, receptor = identity
        event_links[identity] = '<p>当前事件表达证据：' + ' · '.join(
            f'<a data-evidence-link href="#{links[key]}">{title}</a>' for key, title in
            [(('sender', hmdb), '发送端'), (('receiver', (hmdb, receptor)), '接收端')]) + '</p>'
    sample_note = ('按细胞类型合并输入中的全部样本，不随上方样本选择切换。'
                   if result.parameters.get('sample_mode') == 'sample_aware' else '')
    content = ('<div id="expression-evidence"><p class="report-note">' + sample_note +
               '表达证据不等于实测代谢物浓度或外排通量。</p>' + ''.join(charts) +
               '<details class="expression-reading"><summary>图形读法</summary><p>'
               '全零或常数分布显示水平线。发送端色标固定为 0–1；接收端均值全为零时，'
               '色标也使用 0–1。</p></details></div>')
    return content, event_links
