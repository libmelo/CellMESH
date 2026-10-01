"""Static event matrices: one full communication identity per occupied position."""
from __future__ import annotations

import math
import re
import statistics
from numbers import Integral

import pandas as pd

KEYS = ['sender', 'receiver', 'hmdb_id', 'sensor_gene']

def validate_plot_top_n(value):
    """Validate a shared display limit without depending on HTML reporting."""
    if value is not None and (isinstance(value, bool) or not isinstance(value, Integral) or value <= 0):
        raise ValueError("plot_top_n must be a positive integer or None")
    return None if value is None else int(value)


def communication_scope(result_or_events, *, top_n=None, statistic=None):
    """Rank pairs once; retain complete events, including missing scores, within them.

    Pairs with plottable scores follow dotplot ranking. Pairs with no scores
    follow them, so their identities can still be explored without inventing
    a numerical score. This is a display limit, never a statistical filter.
    """
    from .plotting import (_events_frame, _canonical_event_table, _plot_numeric_series,
                           _metabolite_display_names, _metabolite_label)
    from .selection import EventSelection, event_id
    top_n = validate_plot_top_n(top_n)
    if isinstance(result_or_events, EventSelection):
        frame = result_or_events.events
        statistic = statistic or result_or_events.criteria['positive']['column']
    else:
        frame = _events_frame(result_or_events)
        statistic = statistic or 'fdr_sensor_type'
    required = set(KEYS + ['cell_mesh_score'])
    if missing := required.difference(frame.columns):
        raise KeyError('Missing event columns: ' + ', '.join(sorted(missing)))
    if 'sample' in frame:
        raise ValueError('Communication dotmap require aggregate events, not sample rows')
    frame = _canonical_event_table(frame).reset_index(drop=True)
    if frame.duplicated(KEYS).any():
        raise ValueError('Communication dotmap require unique complete event identities')
    frame['cell_mesh_score'] = _plot_numeric_series(frame.cell_mesh_score, 'cell_mesh_score', upper=1)
    for column in dict.fromkeys([statistic, 'perm_pvalue', 'fdr_sensor_type', 'fdr_global']):
        if column in frame:
            frame[column] = _plot_numeric_series(frame[column], column, upper=1, probability=True)
    names = _metabolite_display_names(frame)
    pair_columns = ['hmdb_id', 'sensor_gene']
    valid = frame.loc[frame.cell_mesh_score.notna()]
    ranking = dict(best_score=('cell_mesh_score', 'max'))
    if statistic in frame:
        ranking['best_statistic'] = (statistic, 'min')
    order = (['best_statistic'] if statistic in frame else []) + ['best_score']
    ascending = ([True] if statistic in frame else []) + [False]
    ranked = valid.groupby(pair_columns, as_index=False, sort=True).agg(**ranking)
    ranked = ranked.sort_values(order + pair_columns, ascending=ascending + [True, True], na_position='last')
    pairs = list(ranked[pair_columns].itertuples(index=False, name=None))
    seen = set(pairs)
    pairs.extend(pair for pair in frame[pair_columns].drop_duplicates().itertuples(index=False, name=None)
                 if pair not in seen)
    n_pairs = len(pairs)
    if top_n is not None:
        pairs = pairs[:top_n]
    positions = {pair: i for i, pair in enumerate(pairs)}
    pair_values = list(frame[pair_columns].itertuples(index=False, name=None))
    selected = frame.loc[[pair in positions for pair in pair_values]].copy()
    # Preserve event ordering within a pair; identity is never an array index.
    selected['_pair_order'] = [positions[pair] for pair in selected[pair_columns].itertuples(index=False, name=None)]
    selected = selected.sort_values('_pair_order', kind='stable').drop(columns='_pair_order').reset_index(drop=True)
    records = []
    for row in selected.to_dict('records'):
        record = {key: str(row[key]) for key in KEYS}
        record.update(id=event_id(row), pair=positions[(row['hmdb_id'], row['sensor_gene'])],
                      metabolite=_metabolite_label(names.get(row['hmdb_id']), row['hmdb_id']))
        for key in ['cell_mesh_score', statistic, 'sensor_type', 'positive_status']:
            value = row.get(key)
            record[key] = None if value is None or pd.isna(value) else str(value)
        records.append(record)
    pair_records = [dict(hmdb_id=hmdb, sensor_gene=sensor,
                         label=_metabolite_label(names.get(hmdb), hmdb) + ' → ' + sensor)
                    for hmdb, sensor in pairs]
    return dict(frame=frame, events=selected, pairs=pairs, pair_records=pair_records,
                records=records, n_selected=len(frame), n_pairs=n_pairs, top_n=top_n, statistic=statistic)



def legend_numbers(values):
    """Compact ticks without merging distinct values or rounding a small value to zero."""
    values = [0. if value == 0 else value for value in values]
    for precision in range(3, 18):
        labels = [format(value, f'.{precision}g') for value in values]
        if len(set(labels)) == len(set(values)):
            return labels
    return labels


def statistic_caption(statistic, *, english=False):
    captions = {
        'fdr_sensor_type': ('BH 校正 p', '按 sensor 类型进行 BH 校正',
                            'BH-adjusted p', 'BH correction: within each sensor type'),
        'fdr_global': ('BH 校正 p', '对全部受检事件进行 BH 校正',
                       'BH-adjusted p', 'BH correction: all tested events'),
        'perm_pvalue': ('原始 p', '', 'raw p', ''),
    }
    values = captions.get(statistic, (statistic, '', statistic, ''))
    return values[2:] if english else values[:2]


def node_label_lines(label, max_units=28):
    """Wrap at words where possible; never drop the end of a scientific label."""
    width = lambda text: sum(2 if ord(char) > 127 else 1 for char in text)
    lines, line = [], ''
    for word in re.findall(r'\S+', label):
        if line and width(line + ' ' + word) <= max_units:
            line += ' ' + word
            continue
        if line:
            lines.append(line)
        line = ''
        for char in word:
            if line and width(line + char) > max_units:
                lines.append(line)
                line = ''
            line += char
    if line:
        lines.append(line)
    return lines or ['']


def matrix_axes(scope):
    """Build composite axes for one matrix, without aggregating or splitting events."""
    records = scope['records']
    return dict(x_keys=sorted({(r['receiver'], r['sensor_gene']) for r in records}),
                y_keys=sorted({(r['sender'], r['hmdb_id']) for r in records}))


def dotmap_encoding(scope, cmap='viridis'):
    """Observed-value scales for the complete matrix; never replace zero probabilities."""
    from matplotlib import colormaps
    colormaps[cmap]  # Validate even for a constant or empty scope.
    records, statistic = scope['records'], scope['statistic']
    scores = [float(r['cell_mesh_score']) for r in records if r['cell_mesh_score'] is not None]
    values = [float(r[statistic]) for r in records if r.get(statistic) is not None]
    logs = [-math.log10(v) for v in values if v > 0]
    return dict(statistic=statistic, cmap=cmap,
                score_range=[min(scores), max(scores)] if scores else None,
                score_ticks=sorted({min(scores), statistics.median(scores), max(scores)}) if scores else [],
                log_range=[min(logs), max(logs)] if logs else None,
                has_zero=0 in values, missing_statistic=len(values) < len(records),
                missing_score=len(scores) < len(records))


def dot_area(value, encoding):
    if value is None:
        return 60.
    low, high = encoding['score_range']
    fraction = .5 if low == high else (float(value) - low) / (high - low)
    return 24. + 216. * fraction


def dot_color(value, encoding):
    from matplotlib import colormaps
    from matplotlib.colors import to_hex
    if value is None:
        return '#888888'
    if float(value) == 0:
        return '#b2186b'
    low, high = encoding['log_range']
    fraction = .5 if low == high else (-math.log10(float(value)) - low) / (high - low)
    return to_hex(colormaps[encoding['cmap']](float(fraction)))


def dotmap_labels(scope, axes, *, metabolite_label='id'):
    """Single-line composite ticks; display names never replace HMDB identity keys."""
    if not isinstance(metabolite_label, str) or metabolite_label not in ('id', 'name'):
        raise ValueError("metabolite_label must be 'id' or 'name'")
    names = {r['hmdb_id']: r['metabolite'] for r in scope['records']}
    single_line = lambda value: ' '.join(str(value).split())
    for hmdb, label in names.items():
        suffix = f" ({hmdb})"
        names[hmdb] = single_line(label[:-len(suffix)] if label.endswith(suffix) else label)
    counts = {}
    for name in names.values():
        counts[name] = counts.get(name, 0) + 1
    columns = [f'{single_line(receiver)} – {single_line(sensor)}' for receiver, sensor in axes['x_keys']]
    rows = []
    for sender, hmdb in axes['y_keys']:
        label = hmdb if metabolite_label == 'id' else names[hmdb]
        if metabolite_label == 'name' and counts[label] > 1:
            label += f' ({hmdb})'
        rows.append(f'{single_line(sender)} – {label}')
    return dict(columns=columns, rows=rows, context='', rotation=45,
                metabolite_label=metabolite_label)


def _draw_dotmap_legend(fig, encoding, *, constant_probability=None):
    """Separate score circles, a continuous colorbar, and special-value swatches."""
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle

    limits = encoding['log_range']
    continuous = limits is not None and limits[0] != limits[1]
    stat_name, context = statistic_caption(encoding['statistic'], english=True)
    score_rows = len(encoding['score_ticks']) + int(encoding['missing_score'])
    heading = '−log10(' + stat_name + ')' if continuous else stat_name
    heading_lines = node_label_lines(heading, 28)
    context_lines = node_label_lines(context, 40) if context else []
    heading_height = .23 * len(heading_lines) + .19 * len(context_lines) + .16
    score_height = .26 + .30 * score_rows
    special_rows = int(encoding['has_zero']) + int(encoding['missing_statistic'])
    stat_height = heading_height + (1.8 if continuous else .62 if limits else 0.)
    stat_height += .35 * special_rows + .15
    width, height = 2.5, score_height + .20 + stat_height
    legend_ax = fig.add_axes([.78, .18, .18, .65], label='communication-dotmap-legend')
    legend_ax.set(xlim=(0, width), ylim=(0, height))
    legend_ax.set_axis_off()
    main = dict(fontsize=11, color='#25313b', va='center', fontfamily='sans-serif')
    y = height - .05
    legend_ax.text(0, y, 'Communication score', fontsize=11.5, weight='semibold', va='top',
                   color='#25313b', fontfamily='sans-serif')
    y -= .28
    for score, label in zip(encoding['score_ticks'], legend_numbers(encoding['score_ticks'])):
        legend_ax.scatter([.15], [y], s=[dot_area(score, encoding)], color='#59656f', linewidths=.4)
        legend_ax.text(.43, y, label, **main)
        y -= .30
    if encoding['missing_score']:
        legend_ax.scatter([.15], [y], s=[60], color='#59656f', marker='x', linewidths=1.1)
        legend_ax.text(.43, y, 'Score unavailable', **main)
        y -= .30
    y -= .15
    legend_ax.text(0, y, '\n'.join(heading_lines), fontsize=11.5, weight='semibold', va='top',
                   linespacing=1.2, color='#25313b', fontfamily='sans-serif')
    y -= .23 * len(heading_lines)
    if context_lines:
        legend_ax.text(0, y, '\n'.join(context_lines), fontsize=10, va='top',
                       color='#526773', fontfamily='sans-serif')
        y -= .19 * len(context_lines)
    y -= .18
    colorbar, bar_box = None, None
    ticks = []
    if continuous:
        # Keep a physical 0.28 × 1.8 inch bar, including when the matrix is very tall.
        bar_box = (0., y - 1.8, .28, 1.8)
        color_ax = fig.add_axes([.80, .20, .02, .25], label='communication-dotmap-colorbar')
        colorbar = fig.colorbar(ScalarMappable(norm=Normalize(*limits), cmap=encoding['cmap']), cax=color_ax)
        count = 5 if limits[1]-limits[0] >= 1 else 3
        ticks = [limits[0] + (limits[1]-limits[0])*i/(count-1) for i in range(count)]
        colorbar.set_ticks(ticks, labels=legend_numbers(ticks))
        colorbar.ax.tick_params(labelsize=11, length=3, width=.65, color='#526773', labelcolor='#25313b', pad=6)
        for tick in colorbar.ax.get_yticklabels():
            tick.set_fontfamily('sans-serif')
        colorbar.outline.set_linewidth(.65)
        colorbar.outline.set_edgecolor('#aab5bd')
        y -= 2.00
    elif limits is not None:
        probability = constant_probability
        y -= .10
        legend_ax.add_patch(Rectangle((0., y-.08), .28, .18,
                                      facecolor=dot_color(probability, encoding), edgecolor='#8d9aa3', linewidth=.5))
        constant = legend_ax.text(.43, y, '= ' + legend_numbers([probability])[0], **main)
        note = 'All nonzero values equal' if encoding['has_zero'] else 'All available values equal'
        legend_ax.text(0, y-.23, note, fontsize=10, va='top', color='#526773', fontfamily='sans-serif')
        y -= .58
    for active, color, text in [(encoding['has_zero'], '#b2186b', stat_name + ' = 0'),
                                 (encoding['missing_statistic'], '#888888', 'Statistic unavailable')]:
        if active:
            legend_ax.add_patch(Rectangle((0., y-.08), .28, .18, facecolor=color,
                                          edgecolor='#8d9aa3', linewidth=.5))
            legend_ax.text(.43, y, text, **main)
            y -= .35
    return dict(ax=legend_ax, colorbar=colorbar, bar_box=bar_box, width=width, height=height,
                constant_artist=constant if limits is not None and not continuous else None,
                kind='continuous' if continuous else 'constant' if limits else 'unavailable', ticks=ticks)


def _layout_dotmap(fig, ax, legend, *, auto, n_columns, n_rows):
    """Fit adjacent full labels; size the matrix independently of the legend."""
    from matplotlib.transforms import Bbox
    supplied = ax.get_position().frozen()
    sidebar_width = legend['width']

    def tick_spacing(ticks, renderer, *, vertical=False):
        # Measure protrusion on each side of the anchor, rather than assigning
        # the widest label's full width to every column. Only neighbors collide.
        extents = []
        for tick in ticks:
            box = tick.get_window_extent(renderer)
            x, y = tick.get_transform().transform(tick.get_position())
            extents.append((box.y0-y, box.y1-y) if vertical else (box.x0-x, box.x1-x))
        gaps = [(second[1]-first[0] if vertical else first[1]-second[0])/fig.dpi + .10
                for first, second in zip(extents, extents[1:])]
        return max([.44 if vertical else .50] + gaps)

    def place_legend(x, top):
        fw, fh = fig.get_size_inches()
        bottom = top - legend['height']
        legend['ax'].set_position([x/fw, bottom/fh, legend['width']/fw, legend['height']/fh])
        if legend['colorbar'] is not None:
            bx, by, bw, bh = legend['bar_box']
            legend['colorbar'].ax.set_position([(x+bx)/fw, (bottom+by)/fh, bw/fw, bh/fh])
    # Single-line labels at a fixed angle occupy parallel strips. Allocate
    # their perpendicular thickness, not their overlapping axis-aligned boxes.
    # This lets long labels extend below the plot without widening every column.
    for tick in ax.get_xticklabels():
        tick.set_rotation(0)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    thickness = max([0.] + [tick.get_window_extent(renderer).height/fig.dpi for tick in ax.get_xticklabels()])
    rotation = 45
    column_width = max(.50, (thickness+.10)/math.sin(math.radians(rotation)))
    for tick in ax.get_xticklabels():
        tick.set_rotation(rotation)
    row_height = tick_spacing(ax.get_yticklabels(), renderer, vertical=True)
    grid_width = max(2.2, n_columns*column_width)
    grid_height = max(1.6, n_rows*row_height)

    # Measure decorations after positioning; keep the sidebar's physical size
    # without stretching a short matrix to match the sidebar's height.
    for _ in range(3):
        fw, fh = fig.get_size_inches()
        position = ax.get_position()
        place_legend(position.x1*fw + .32, position.y1*fh)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        dpi = fig.dpi
        body = ax.get_window_extent(renderer)
        outer = ax.get_tightbbox(renderer)
        extras = [body.x0-outer.x0, body.y0-outer.y0, outer.x1-body.x1, outer.y1-body.y1]
        left, bottom, right, top = [max(value/dpi, 0.) for value in extras]
        sidebar_boxes = [legend['ax'].get_tightbbox(renderer)]
        if legend['colorbar'] is not None:
            sidebar_boxes.append(legend['colorbar'].ax.get_tightbbox(renderer))
        sidebar_box = Bbox.union(sidebar_boxes)
        sidebar_width = max(sidebar_width, sidebar_box.width/dpi)
        if not auto:
            ax.set_position(supplied)
            place_legend(supplied.x1*fw + right + .32, supplied.y1*fh)
            break
        lower_space = max(bottom, legend['height']-grid_height)
        # SVG viewers may substitute a wider sans-serif font than Agg measured.
        # Keep 26 extra points at the right edge for the full BH scope label.
        width = .72 + left + grid_width + right + .32 + sidebar_width
        height = .36 + lower_space + grid_height + top
        fig.set_size_inches(width, height, forward=True)
        ax.set_position([(.18+left)/width, (.18+lower_space)/height, grid_width/width, grid_height/height])
        place_legend(.18+left+grid_width+right+.32, .18+lower_space+grid_height)
    return dict(font_sizes=dict(title=14, axis=12, tick=12, legend=11, legend_title=11.5, note=10),
                legend_kind=legend['kind'], colorbar_ticks=legend['ticks'],
                automatic_canvas=auto, x_rotation=rotation)


def plot_communication_dotmap(result_or_events, *, top_n=None, fdr_col=None, cmap="viridis", ax=None,
                              metabolite_label="id"):
    """One static matrix: receiver-sensor columns × sender-metabolite rows.

    Every complete event occupies exactly one position. All metabolites and
    sensors in the display scope share the same axes and scales; the figure
    grows with axis size and is never automatically split or paginated.
    Dot area maps communication score; color maps -log10 of the chosen p/BH
    statistic. No implicit statistical or QC filter is applied. Missing scores
    use crosses, missing statistics gray, exact zero statistics a separate color.
    ``top_n`` limits HMDB-sensor pairs, retaining every event within those pairs.
    ``fdr_col`` defaults to the selection's positive criterion or fdr_sensor_type.
    Every tick shows a single-line composite identity, including shared molecules.
    ``metabolite_label='id'`` displays HMDB IDs (default); ``'name'`` displays
    metabolite names, falling back to IDs when absent and adding IDs to ambiguous
    names. Identity keys always use HMDB IDs. Labels use a 12-point font, with
    x labels fixed at 45 degrees and y labels horizontal. Matrix and legend
    heights are allocated independently.
    Varying statistics use a vertical continuous colorbar; constant values use a
    solid swatch labelled with the original probability. Score circles have a
    separate legend. Returned labels/layout, legend_ax and colorbar expose these
    decorations alongside fig, ax, plot_events, dot_table, encoding and x/y_keys.
    The returned scope records the exact display range, shared with report evidence.
    Use fig.savefig for a complete PNG/PDF/SVG export; reports embed this Figure.
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    scope = communication_scope(result_or_events, top_n=top_n, statistic=fdr_col)
    with mpl.rc_context({'text.usetex': False, 'text.parse_math': False}):
        encoding = dotmap_encoding(scope, cmap)
        axes = matrix_axes(scope)
        labels = dotmap_labels(scope, axes, metabolite_label=metabolite_label)
        auto = ax is None
        if auto:
            fig = plt.figure(figsize=(10.5, 6.3), layout='none')
            ax = fig.add_axes([.15, .25, .52, .55])
        else:
            fig = ax.figure
        ax.set_xlim(-.5, max(.5, len(labels['columns'])-.5))
        ax.set_ylim(max(.5, len(labels['rows'])-.5), -.5)
        ax.set_xticks(range(len(labels['columns'])), labels['columns'], rotation=labels['rotation'],
                      ha='right', rotation_mode='anchor', fontsize=12, fontfamily='sans-serif')
        ax.set_yticks(range(len(labels['rows'])), labels['rows'], fontsize=12, fontfamily='sans-serif')
        ax.set_xlabel('Receiver – sensor', fontsize=12, labelpad=12, fontfamily='sans-serif')
        ax.set_ylabel('Sender – metabolite' + (' (HMDB ID)' if metabolite_label == 'id' else ''),
                      fontsize=12, labelpad=12, fontfamily='sans-serif')
        ax.grid(color='#e5e9ed', linewidth=.65)
        ax.set_axisbelow(True)
        ax.tick_params(length=0, pad=9, labelcolor='#25313b')
        for spine in ax.spines.values():
            spine.set_visible(False)
        context_artist = None
        title = ax.set_title(f'Selected cell–cell communication events (n = {len(scope["records"])})',
                            fontsize=14, weight='semibold', pad=28, loc='center', fontfamily='sans-serif')
        xs = {key: i for i, key in enumerate(axes['x_keys'])}
        ys = {key: i for i, key in enumerate(axes['y_keys'])}
        dots = []
        for record in scope['records']:
            score, value = record['cell_mesh_score'], record.get(encoding['statistic'])
            x, y = xs[(record['receiver'], record['sensor_gene'])], ys[(record['sender'], record['hmdb_id'])]
            color, area = dot_color(value, encoding), dot_area(score, encoding)
            marker = 'x' if score is None else 'o'
            artist = ax.scatter([x], [y], s=[area], c=[color], marker=marker, zorder=3,
                                linewidths=1.1 if marker == 'x' else .4)
            artist.set_gid('event-' + record['id'])
            dots.append(dict(id=record['id'], x=x, y=y, area=area, color=color, marker=marker,
                             **{key: record[key] for key in KEYS}))
        if not scope['records']:
            ax.text(.5, .5, 'No selected events', ha='center', va='center', transform=ax.transAxes)
        probability = next((float(r[encoding['statistic']]) for r in scope['records']
                            if r.get(encoding['statistic']) is not None and float(r[encoding['statistic']]) > 0), None)
        legend = _draw_dotmap_legend(fig, encoding, constant_probability=probability)
        layout = _layout_dotmap(fig, ax, legend, auto=auto,
                                n_columns=len(axes['x_keys']), n_rows=len(axes['y_keys']))
        labels['rotation'] = layout['x_rotation']
        return dict(fig=fig, ax=ax, title_artist=title, context_artist=context_artist,
                    legend_ax=legend['ax'], colorbar=legend['colorbar'], layout=layout, labels=labels,
                    plot_events=scope['events'].copy(),
                    dot_table=pd.DataFrame(dots, columns=['id', 'x', 'y', 'area', 'color', 'marker'] + KEYS),
                    encoding=encoding, scope=scope, **axes,
                    selection=dict(n_selected=scope['n_selected'], n_in_scope=len(scope['records']),
                                   n_display_excluded=scope['n_selected']-len(scope['records'])))
