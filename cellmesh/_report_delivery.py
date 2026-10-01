"""Report paths and rendering for independent general, selection and event pages."""
from __future__ import annotations

import os
from pathlib import Path

from ._report_browser import event_id, page_entries


def bundle_paths(path, delivery, *, include_input=False):
    stem = path.name[:-len('.cellmesh.h5')]
    paths = {'result': path, 'positive_csv': path.with_name(stem + '.positive_events.csv')}
    if delivery.get('include_all_events', False):
        paths['events_csv'] = path.with_name(stem + '.all_events.csv')
    if delivery.get('selection') is not None:
        paths['filtered_csv'] = path.with_name(stem + ('.filtered_events.csv' if delivery.get('report_layout', 3) >= 4 else '.selected_events.csv'))
    if delivery.get('pair_selections') is not None:
        paths['selected_csv'] = path.with_name(stem + '.metabolite_sensor_events.csv')
    if include_input:
        paths['input'] = path.with_name(stem + '.input.h5ad')
    if delivery.get('report_enabled', True):
        paths['report'] = path.with_name(stem + '.report.html')
        paths['navigation'] = path.parent / (stem + '.reports') / 'index.html'
        if delivery.get('selection') is not None:
            paths['filtered_report'] = path.with_name(stem + ('.filtered.report.html' if delivery.get('report_layout', 3) >= 4 else '.selected.report.html'))
        for selector in delivery.get('events') or []:
            key = event_id(selector)
            paths['event_report:' + key] = paths['navigation'].parent / 'events' / (key + '.html')
        for record in delivery.get('pair_selections') or []:
            paths['pair_report:' + record['id']] = paths['navigation'].parent / 'metabolite_sensor' / (record['id'] + '.html')
    return paths


def relative(path, parent):
    return os.path.relpath(path, parent)


def entries_for(paths, delivery, parent):
    selectors = [selector for selector in delivery.get('events') or []
                 if 'event_report:' + event_id(selector) in paths]
    hrefs = {event_id(selector): relative(paths['event_report:' + event_id(selector)], parent)
             for selector in selectors}
    entries = page_entries(selectors, hrefs)
    for record in delivery.get('pair_selections') or []:
        key = 'pair_report:' + record['id']
        if key in paths:
            for event in record['keys']:
                entries.append(dict(id=record['id'], label=' / '.join(record['pair']),
                                    scope='代谢物–sensor 通讯报告',
                                    identity=tuple(event[k] for k in ('sender', 'receiver', 'hmdb_id', 'sensor_gene')),
                                    href=relative(paths[key], parent)))
    return entries


def render_bundle(result, paths, delivery, selection, *, requested=None, tables=None,
                  adata=None, input_storage=None, selection_evidence_path=None):
    from .reporting import render_html_report, build_report_tables
    from ._report_pages import render_event, render_navigation
    requested = set(paths) if requested is None else set(requested)
    contents = {}
    cache = {}
    needs_tables = 'report' in requested or any(key.startswith(('event_report:', 'pair_report:')) for key in requested)
    if needs_tables and tables is None:
        tables = build_report_tables(result)
    for key in ['report', 'filtered_report']:
        if key not in requested or key not in paths:
            continue
        target = paths[key]
        contents[key] = render_html_report(result, tables=tables, delivery=delivery,
            adata=adata, selection_evidence_path=selection_evidence_path if key == 'filtered_report' else None,
            selection=selection if key == 'filtered_report' else None,
            plot_top_n=delivery.get('plot_top_n', 20), input_storage=input_storage,
            csv_paths={k: relative(v, target.parent) for k, v in paths.items() if k.endswith('csv')},
            event_links=entries_for(paths, delivery, target.parent) if key == 'filtered_report' else None,
            general_report=relative(paths['report'], target.parent) if key == 'filtered_report' else None,
            navigation_report=relative(paths['navigation'], target.parent))
    for selector in delivery.get('events') or []:
        key = 'event_report:' + event_id(selector)
        if key not in requested or key not in paths:
            continue
        target = paths[key]
        contents[key] = render_event(result, selector, tables=tables, spec=delivery['positive'], adata=adata,
            general_report=relative(paths['report'], target.parent) if 'report' in paths else None,
            navigation_report=relative(paths['navigation'], target.parent) if 'navigation' in paths else None,
            evidence_cache=cache)
    from ._pair_selection import restore_pairs
    from ._pair_reports import render_pair
    for record in delivery.get('pair_selections') or []:
        key = 'pair_report:' + record['id']
        if key in requested:
            target = paths[key]
            focus = restore_pairs(result, [record])[0]
            contents[key] = render_pair(result, focus, tables=tables, spec=delivery['positive'],
                plot_top_n=delivery.get('plot_top_n', 20), adata=adata,
                csv_path=relative(paths['selected_csv'], target.parent),
                general_report=relative(paths['report'], target.parent) if 'report' in paths else None,
                navigation_report=relative(paths['navigation'], target.parent) if 'navigation' in paths else None)
    if 'navigation' in requested and 'navigation' in paths:
        parent = paths['navigation'].parent
        notes = []
        for key, label in [('report', '一级统计报告'), ('filtered_report', '通讯筛选报告')]:
            if key in paths:
                notes.append(label + ('：本次调用生成或更新。' if key in requested else
                                      '：沿用现有交付，本次调用未重新生成。'))
        pair_keys = {key for key in paths if key.startswith('pair_report:')}
        if pair_keys:
            updated = len(pair_keys & requested)
            notes.append(f'代谢物–sensor 通讯报告：本次调用生成或更新 {updated} 份，'
                         f'沿用 {len(pair_keys) - updated} 份。')
        contents['navigation'] = render_navigation(
            general=relative(paths['report'], parent) if 'report' in paths else None,
            selection=relative(paths['filtered_report'], parent) if 'filtered_report' in paths else None,
            entries=entries_for(paths, delivery, parent),
            pairs=[dict(id=r['id'], label=' / '.join(r['pair']),
                        href=relative(paths['pair_report:' + r['id']], parent))
                   for r in delivery.get('pair_selections') or []], notes=notes)
    return contents
