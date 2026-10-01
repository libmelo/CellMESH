"""Selection overview figures for portable HTML reports.

Use the public plotting functions so reports and standalone plots share the
same event identities, ranking and missing-score policy.
"""
from __future__ import annotations

import html
import json
import re
from ._communication_dotmap import validate_plot_top_n


def selection_plots_html(selection, plot_top_n=20, *, prefix="selection", scope_label="阳性事件",
                         fdr_col=None, include_dotplot=True, dotmap=None, secondary=False):
    from .selection import EVENT_KEYS
    from .plotting import (plot_communication_network,
                           plot_significant_event_counts)

    plot_top_n = validate_plot_top_n(plot_top_n)
    frame = selection.events
    positions = {key: i for i, key in enumerate(frame[EVENT_KEYS].itertuples(index=False, name=None))}
    n_selected = len(frame)
    n_missing = int(frame.cell_mesh_score.isna().sum())
    definitions = [("counts", "细胞对事件数（热图）", plot_significant_event_counts),
                   ("network", "通讯网络图", plot_communication_network)]
    omitted_counts = False
    if secondary:
        valid = frame.loc[frame.cell_mesh_score.notna()]
        counts = valid.groupby(['sender', 'receiver']).size()
        n_relations = len(valid[['hmdb_id', 'sensor_gene']].drop_duplicates())
        if n_relations <= 1 or counts.nunique() <= 1:
            definitions = [entry for entry in definitions if entry[0] != 'counts']
            omitted_counts = True
    description = '网络图展示细胞连接关系；计数热图在具有事件数量差异时显示' if secondary else '热图比较各细胞对的事件数；网络图展示连接关系'
    sections = [
        '<section id="selection-plots" aria-label="通讯概览图形">',
        '<p class="report-note">' + description + '。图形范围：' + html.escape(scope_label)
        + '。计数与网络不纳入缺失评分的事件。</p>',
    ]
    if omitted_counts and n_selected:
        sections.append('<p class="report-note">计数热图用于比较多个代谢物—受体组合在细胞对间的数量差异；当前可绘结果不满足此条件，省去该图。网络图仍按适用条件保留。</p>')
    n_pairs = len(frame.loc[frame.cell_mesh_score.notna(), ["sender", "receiver"]].drop_duplicates())
    if n_pairs <= 1:
        definitions = []
        if n_pairs == 1:
            row = frame.loc[frame.cell_mesh_score.notna()].iloc[0]
            sections.append('<p class="compact-pair">' + html.escape(str(row.sender)) + ' → ' +
                            html.escape(str(row.receiver)) + f'：{n_selected - n_missing:,} 条可绘事件。</p>')
    if n_selected <= 1 or n_selected == n_missing:
        definitions = []
        message = ('没有事件满足全部筛选条件。' if not n_selected else
                   '所选事件全部缺失评分，省去计数和网络；完整记录保留在列表和 CSV 中。' if n_selected == n_missing else
                   ('仅选择了一个事件，省去计数和网络。' if secondary else
                    '仅选择了一个事件，省去计数和网络；单事件报告需另行指定。'))
        sections.append('<p class="plot-empty">' + message + '</p>')
    records = []
    if include_dotplot:
        from .plotting import plot_communication_dotmap
        from ._report_dotmap import communication_dotmap_html
        import matplotlib.pyplot as plt
        output = dotmap if dotmap is not None else plot_communication_dotmap(
            selection, top_n=plot_top_n, fdr_col=fdr_col)
        try:
            matrix_html, matrix_record = communication_dotmap_html(output)
        finally:
            if dotmap is None:
                plt.close(output['fig'])
        if matrix_html:
            sections.insert(1, matrix_html)
            records.append(matrix_record)
    # Empty/fully missing selections can be reported without plotting extras.
    if n_selected > n_missing:
        try:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("Selection report plots require matplotlib; install cellmesh[report]") from exc
    for name, title, draw in definitions:
        record = dict(kind=name, n_selected=n_selected, n_missing_score=n_missing,
                      n_plotted=0, n_display_excluded=0, plotted_event_indices=[],
                      plot_top_n=None, top_n_unit="hmdb_id + sensor_gene")
        svg = ""
        if n_selected > n_missing:
            # Local settings avoid TeX/math interpretation of arbitrary labels,
            # preserve accessible text and keep SVG IDs separate across figures.
            with mpl.rc_context({"text.usetex": False, "text.parse_math": False,
                                 "svg.fonttype": "none", "svg.hashsalt": f"cellmesh-{name}"}):
                options = {}
                source = selection
                if fdr_col is not None:
                    options["fdr_col"] = fdr_col
                output = draw(source, **options)
                try:
                    plotted = output.get("plot_events", output.get("filtered_events"))
                    indices = [positions[key] for key in plotted[EVENT_KEYS].itertuples(index=False, name=None)]
                    record.update(n_plotted=len(plotted), plotted_event_indices=sorted(indices),
                                  n_missing_score=n_missing,
                                  n_display_excluded=n_selected - n_missing - len(plotted))
                    if name == "counts":
                        record["count_total"] = int(output["counts"].to_numpy().sum())
                        from matplotlib.ticker import MaxNLocator
                        output["fig"].axes[-1].yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=1))
                    if name == "network":
                        record["edge_event_total"] = int(output["edge_table"]["event_count"].sum())
                    output["ax"].set_title("")
                    from ._report_dotmap import figure_svg
                    svg = figure_svg(output['fig'], f'selection-{name}')
                finally:
                    plt.close(output["fig"])
        records.append(record)
        sections += [f'<details class="selection-chart" id="selection-{name}"' + (' open' if not secondary and not include_dotplot and name == 'counts' else '') + '>',
                     f'<summary>{html.escape(title)}<span class="chart-summary">{record["n_plotted"]:,} 条事件</span></summary>',
                     f'<p class="plot-scope">入选 {n_selected:,} 条 · 绘制 {record["n_plotted"]:,} 条 · '
                     f'缺失评分 {n_missing:,} 条 · 因展示上限未绘制 {record["n_display_excluded"]:,} 条</p>']
        if svg:
            sections.append('<div class="selection-plot">' + svg + '</div>')
        else:
            message = '没有事件满足全部筛选条件。' if not n_selected else '所选事件全部缺失评分，无法绘图；事件仍保留在 CSV 和列表中。'
            sections.append(f'<p class="plot-empty">{message}</p>')
        sections.append('</details>')
    payload = json.dumps(records, ensure_ascii=True, allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    sections += ['<script type="application/json" id="selection-plot-data">' + payload + '</script>', '</section>']
    content = ''.join(sections)
    if prefix != "selection":
        def namespace(match):
            tag = re.sub(r'\bid="selection-', f'id="{prefix}-', match[0])
            tag = re.sub(r'\bhref="#selection-', f'href="#{prefix}-', tag)
            return tag.replace('url(#selection-', f'url(#{prefix}-')
        content = re.sub(r'<[^>]+>', namespace, content)
    return content
