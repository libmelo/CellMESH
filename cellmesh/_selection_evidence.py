"""Deduplicated expression comparisons for a selected collection of events."""
from __future__ import annotations

import hashlib
import html
import json

from ._report_evidence import _figure_svg

START = '<!--cellmesh-selection-expression:start-->'
END = '<!--cellmesh-selection-expression:end-->'


def evidence_scope_key(result, scope):
    value = dict(version=1, events=[record['id'] for record in scope['records']],
                 pairs=scope['pairs'], top_n=scope['top_n'],
                 cell_type_key=result.parameters.get('cell_type_key'), layer=result.parameters.get('layer'))
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=True).encode()).hexdigest()


def saved_selection_evidence(path, result, scope):
    """Reuse only an unchanged selection view from the same saved report bundle."""
    if not path.is_file():
        return None
    source = path.read_text(encoding='utf-8')
    _, separator, remainder = source.partition(START)
    if not separator:
        return None
    block, separator, _ = remainder.partition(END)
    if not separator or f'data-evidence-scope="{evidence_scope_key(result, scope)}"' not in block:
        return None
    # Missing-input notices are not figures: a later explicit adata can fill them.
    if 'data-expression-index=' not in block:
        return None
    return START + block + END


def selection_evidence_html(result, scope, *, adata=None):
    """Use existing violin functions once per metabolite or metabolite/sensor."""
    if not scope['records']:
        return ''
    key = evidence_scope_key(result, scope)
    parts = [START, '<details class="selection-chart" id="interest-expression" '
             f'data-evidence-scope="{key}"><summary>发送端与接收端表达</summary>']
    source = '显式提供的 adata' if adata is not None else '计算时保留的内存输入'
    adata = getattr(result, '_input_adata', None) if adata is None else adata
    if adata is None:
        return ''.join(parts) + ('<p>未提供细胞级表达输入，未生成小提琴图。'
                                '可通过 adata 提供原输入或显式 load_input() 的返回值。</p></details>' + END)
    from .plotting import (plot_metabolite_secretion_violin, plot_receptor_expression_violin,
                           _result_adata_context)
    try:
        cell_key, layer, _, matrix = _result_adata_context(result, adata, cell_type_key=None, layer=None)
    except (KeyError, ValueError, TypeError) as exc:
        return ''.join(parts) + '<p>表达输入不可用：' + html.escape(str(exc)) + '</p></details>' + END
    expression_layer = 'X' if layer is None else f'layers[{layer!r}]'
    parts.append('<p>来源：' + html.escape(source) + f'；输入 {matrix.shape[0]:,} 个细胞 × {matrix.shape[1]:,} 个基因；'
                 '表达矩阵：' + html.escape(expression_layer) + '；细胞类型列：' + html.escape(cell_key) + '。</p>')
    parts.append('<p class="report-note">图形按代谢物或代谢物—受体组合合并，只使用本图范围内事件涉及的相应细胞类型。'
                 '发送端形状为生产反应相关表达分布，颜色为细胞类型代谢物可用性评分（0–1）；'
                 '接收端形状为受体表达分布，颜色为平均表达。表达证据不等于实测代谢物浓度或外排通量。</p>')
    if result.parameters.get('sample_mode') == 'sample_aware':
        parts.append('<p>本图按细胞类型汇集所提供输入中的各样本细胞，不是逐样本评分分布；样本间表现见样本 × 事件评分热图。</p>')
    senders, receivers, names = {}, {}, {}
    for record in scope['records']:
        hmdb, sensor = record['hmdb_id'], record['sensor_gene']
        senders.setdefault(hmdb, set()).add(record['sender'])
        receivers.setdefault((hmdb, sensor), set()).add(record['receiver'])
        names[hmdb] = record['metabolite']
    for role, groups, draw in [('sender', senders, plot_metabolite_secretion_violin),
                               ('receiver', receivers, plot_receptor_expression_violin)]:
        title = '发送端生产相关表达' if role == 'sender' else '接收端受体表达'
        parts.append(f'<section class="selection-expression-group" data-expression-role="{role}"><h3>{title}</h3>')
        options, panels = [], []
        for index, (identity, cells) in enumerate(groups.items()):
            hmdb = identity if role == 'sender' else identity[0]
            label = names[hmdb] + ('' if role == 'sender' else ' → ' + identity[1])
            options.append(f'<option value="{index}">{html.escape(label)}</option>')
            args = dict(result=result, adata=adata, hmdb_id=hmdb)
            args.update(dict(sender_labels=sorted(cells), vmin=0, vmax=1) if role == 'sender' else
                        dict(receiver_labels=sorted(cells), receptor_gene=identity[1]))
            namespace = f'interest-expression-{role}-{index}'
            try:
                chart = '<div class="selection-plot">' + _figure_svg(draw, namespace, **args) + '</div>'
            except (KeyError, ValueError) as exc:
                chart = '<p>该表达图不可用：' + html.escape(str(exc)) + '</p>'
            panels.append(f'<div data-expression-index="{index}"><p><strong>{html.escape(label)}</strong></p>'
                          '<p>细胞范围：' + html.escape('、'.join(sorted(cells))) + '。</p>' + chart + '</div>')
        parts.append('<label class="expression-picker" hidden>查看对象 <select>' + ''.join(options) + '</select></label>')
        parts.extend(panels)
        parts.append('</section>')
    parts.append(f'<p>去重后共 {len(senders):,} 张发送端图、{len(receivers):,} 张接收端图；'
                 f'覆盖本图范围内 {len(scope["records"]):,} 条事件。组合展示上限由 plot_top_n 控制。</p></details>')
    parts.append(END)
    return ''.join(parts)


EVIDENCE_JS = r"""
(() => {
 for (const group of document.querySelectorAll('.selection-expression-group')) {
   const picker=group.querySelector('.expression-picker'), select=picker.querySelector('select');
   const show=()=>group.querySelectorAll('[data-expression-index]').forEach(panel=>{
     panel.hidden=panel.dataset.expressionIndex!==select.value;
   });
   select.addEventListener('change',show);
   picker.hidden=select.options.length<2;
   show();
 }
})();
"""

EVIDENCE_CSS = """
.selection-expression-group{border-top:1px solid #dde4e8;margin-top:12px;padding-top:4px}
.expression-picker select{font:inherit;padding:6px;max-width:100%}
[data-expression-index][hidden]{display:none}
"""
