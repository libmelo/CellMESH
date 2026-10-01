"""Versioned, non-executable HDF5 persistence for CELL MESH results.

No pickle payloads or dynamic class imports are used. Table axes, nullable
dtypes, nested metadata and missing-value identities are encoded explicitly.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

import h5py
import numpy as np
import pandas as pd


from .positivity import INHERIT

FORMAT_VERSION = 1
RESULT_TABLES = ("events", "sender_scores", "receiver_scores", "celltype_qc")
SAMPLE_TABLES = ("sample_validation", "sample_events", "sample_sender_scores", "sample_receiver_scores")
NULL_KEY = "sample_aware_null_scores"


def _pack(value):
    """Tagged JSON only for object cells/scalars; numeric matrices stay binary."""
    if value is pd.NA:
        return ["pd.NA"]
    if value is pd.NaT:
        return ["NaT"]
    if value is None:
        return ["none"]
    if isinstance(value, np.generic):
        if isinstance(value, (np.datetime64, np.timedelta64)):
            return ["numpy_time", str(value.dtype), int(value.astype("int64"))]
        value = value.item()
    if isinstance(value, bool):
        return ["bool", value]
    if isinstance(value, int):
        return ["int", str(value)]
    if isinstance(value, float):
        return ["float", value.hex()]
    if isinstance(value, str):
        return ["str", value]
    if isinstance(value, pd.Timestamp):
        return ["timestamp", value.isoformat()]
    if isinstance(value, pd.Timedelta):
        return ["timedelta", int(value.value)]
    if isinstance(value, (dt.datetime, dt.date)):
        return ["datetime" if isinstance(value, dt.datetime) else "date", value.isoformat()]
    if isinstance(value, (list, tuple, set, frozenset)):
        return [type(value).__name__, [_pack(x) for x in value]]
    if isinstance(value, dict):
        return ["dict", [[_pack(k), _pack(v)] for k, v in value.items()]]
    raise TypeError(f"Unsupported result metadata/cell type: {type(value).__name__}")


def _unpack(value):
    tag, *data = value
    if tag == "pd.NA": return pd.NA
    if tag == "NaT": return pd.NaT
    if tag == "none": return None
    if tag == "str": return data[0]
    if tag == "bool": return bool(data[0])
    if tag == "int": return int(data[0])
    if tag == "float": return float.fromhex(data[0])
    if tag == "timestamp": return pd.Timestamp(data[0])
    if tag == "timedelta": return pd.Timedelta(data[0], unit="ns")
    if tag == "datetime": return dt.datetime.fromisoformat(data[0])
    if tag == "date": return dt.date.fromisoformat(data[0])
    if tag == "numpy_time": return np.array(data[1], dtype="int64").view(data[0])[()]
    if tag in {"list", "tuple", "set", "frozenset"}:
        cls = {"list": list, "tuple": tuple, "set": set, "frozenset": frozenset}[tag]
        return cls(_unpack(x) for x in data[0])
    if tag == "dict": return {_unpack(k): _unpack(v) for k, v in data[0]}
    raise ValueError(f"Unsupported scalar encoding {tag!r}")


def _json_write(group, name, value):
    # JSON escapes embedded NULs; HDF5 variable-length strings cannot contain them.
    payload = json.dumps(_pack(value), ensure_ascii=True, allow_nan=False).encode("utf-8")
    dataset = group.create_dataset(name, data=np.frombuffer(payload, dtype="uint8"),
                                   **({"compression": "gzip"} if len(payload) >= 1024 else {}))
    dataset.attrs["encoding"] = "tagged-json-utf8"


def _json_read(group, name):
    return _unpack(json.loads(group[name][()].tobytes().decode("utf-8")))


def _write(parent, name, value, *, include_null_scores, _column=False):
    group = parent.create_group(name)
    if isinstance(value, pd.DataFrame):
        group.attrs["kind"] = "dataframe"
        _write(group, "index", value.index, include_null_scores=include_null_scores)
        _write(group, "columns", value.columns, include_null_scores=include_null_scores)
        for i in range(value.shape[1]):
            # Series.attrs can inherit DataFrame.attrs; persist it once at table level.
            column = value.iloc[:, i].copy(deep=False)
            column.attrs = {}
            _write(group, f"column_{i}", column, include_null_scores=include_null_scores, _column=True)
        _write_attrs(group, value.attrs, include_null_scores)
    elif isinstance(value, pd.Series):
        group.attrs["kind"] = "series"
        group.attrs["dtype"] = str(value.dtype)
        if not _column:
            _json_write(group, "name", value.name)
            _write(group, "index", value.index, include_null_scores=include_null_scores)
        if isinstance(value.dtype, pd.CategoricalDtype):
            group.attrs["categorical"] = True
            group.attrs["ordered"] = value.cat.ordered
            _write(group, "categories", value.cat.categories, include_null_scores=include_null_scores)
            _write(group, "values", value.cat.codes.to_numpy(), include_null_scores=include_null_scores)
        else:
            if isinstance(value.dtype, pd.StringDtype):
                group.attrs["string_storage"] = value.dtype.storage
            _write(group, "values", value.to_numpy(), include_null_scores=include_null_scores)
        _write_attrs(group, value.attrs, include_null_scores)
    elif isinstance(value, pd.MultiIndex):
        group.attrs["kind"] = "multiindex"
        _json_write(group, "names", list(value.names))
        _json_write(group, "sortorder", value.sortorder)
        for i, (level, codes) in enumerate(zip(value.levels, value.codes)):
            _write(group, f"level_{i}", level, include_null_scores=include_null_scores)
            _write(group, f"codes_{i}", np.asarray(codes), include_null_scores=include_null_scores)
    elif isinstance(value, pd.RangeIndex):
        group.attrs["kind"] = "rangeindex"
        _json_write(group, "range", [value.start, value.stop, value.step, value.name])
    elif isinstance(value, pd.Index):
        group.attrs["kind"] = "index"
        group.attrs["index_class"] = type(value).__name__
        _write(group, "series", pd.Series(value, name=value.name), include_null_scores=include_null_scores)
    elif isinstance(value, np.ndarray):
        group.attrs["kind"] = "array"
        group.attrs["dtype"] = str(value.dtype)
        if value.dtype.kind in "biufc":
            group.create_dataset("values", data=value,
                                 **({"compression": "gzip"} if value.ndim and value.nbytes >= 1024 else {}))
        else:
            group.attrs["json_array"] = True
            group.attrs["shape"] = value.shape
            _json_write(group, "values", value.ravel().tolist())
    elif isinstance(value, dict):
        group.attrs["kind"] = "dict"
        _json_write(group, "keys", list(value.keys()))
        for i, item in enumerate(value.values()):
            _write(group, str(i), item, include_null_scores=include_null_scores)
    else:
        group.attrs["kind"] = "scalar"
        _json_write(group, "value", value)


def _write_attrs(group, attrs, include_null_scores):
    attrs = dict(attrs)
    # The potentially large null matrix has a single canonical storage location.
    attrs.pop(NULL_KEY, None)
    if not include_null_scores:
        if "null_scores_stored" in attrs:
            attrs["null_scores_stored"] = False
    if attrs:
        _write(group, "attrs", attrs, include_null_scores=include_null_scores)


def _read(group):
    kind = group.attrs["kind"]
    if kind == "dataframe":
        index, columns = _read(group["index"]), _read(group["columns"])
        series = [_read(group[f"column_{i}"]) for i in range(len(columns))]
        frame = pd.concat(series, axis=1) if series else pd.DataFrame(index=range(len(index)))
        frame.index, frame.columns = index, columns
        frame.attrs = _read(group["attrs"]) if "attrs" in group else {}
        return frame
    if kind == "series":
        values = _read(group["values"])
        if group.attrs.get("categorical", False):
            values = pd.Categorical.from_codes(values, _read(group["categories"]),
                                                ordered=bool(group.attrs["ordered"]))
            dtype = values.dtype
        else:
            dtype = group.attrs["dtype"]
            if "string_storage" in group.attrs:
                dtype = pd.StringDtype(storage=group.attrs["string_storage"])
        series = pd.Series(values, index=_read(group["index"]) if "index" in group else None,
                           name=_json_read(group, "name") if "name" in group else None, dtype=dtype)
        series.attrs = _read(group["attrs"]) if "attrs" in group else {}
        return series
    if kind == "multiindex":
        names = _json_read(group, "names")
        return pd.MultiIndex(levels=[_read(group[f"level_{i}"]) for i in range(len(names))],
                             codes=[_read(group[f"codes_{i}"]) for i in range(len(names))],
                             names=names, sortorder=_json_read(group, "sortorder"))
    if kind == "rangeindex":
        start, stop, step, name = _json_read(group, "range")
        return pd.RangeIndex(start, stop, step, name=name)
    if kind == "index":
        series = _read(group["series"])
        return pd.Index(series.array, name=series.name, tupleize_cols=False)
    if kind == "array":
        if group.attrs.get("json_array", False):
            items = _json_read(group, "values")
            # np.array(list_of_lists, dtype=object) expands nested cells; assign explicitly.
            array = np.empty(len(items), dtype=group.attrs["dtype"])
            for i, value in enumerate(items):
                array[i] = value
            return array.reshape(tuple(group.attrs["shape"]))
        return group["values"][()]
    if kind == "dict":
        return {key: _read(group[str(i)]) for i, key in enumerate(_json_read(group, "keys"))}
    if kind == "scalar":
        return _json_read(group, "value")
    raise ValueError(f"Unsupported HDF5 node kind {kind!r}")


def _paths(prefix):
    value = os.fspath(prefix)
    suffix = ".cellmesh.h5"
    stem = value[:-len(suffix)] if value.endswith(suffix) else value
    if not Path(stem).name:
        raise ValueError("Result prefix must have a filename component")
    return Path(stem + suffix), Path(stem + ".report.html")


def save_result(result, prefix, *, filtered_events=None, selected_events=None, report=True, events=None, include_null_scores=False,
                include_input=False, include_all_events=False, selection=None, detail_n=INHERIT, plot_top_n=20,
                positive_col=None, positive_threshold=None, positive_qc_only=None, positive_min_score=INHERIT,
                adata=None):
    """Stage independent report scopes and scientific outputs before publication."""
    from . import __version__
    from .reporting import build_report_tables
    from .input_archive import stage_input
    from .selection import report_selection
    from .positivity import resolve_positive, delivery_record
    from ._report_plots import validate_plot_top_n
    from ._report_delivery import bundle_paths, render_bundle
    from ._pair_selection import resolve_scopes, pair_records
    selection, selected_events = resolve_scopes(result, filtered_events, selected_events, selection)
    if events is not None:
        raise TypeError("events= is retired. Use select_events(metabolite_sensor_pairs=...) and selected_events=; "
                        "restrict sender/receiver with cell_pairs to preserve the original scope.")
    plot_top_n = validate_plot_top_n(plot_top_n)
    for name, value in [("report", report), ("include_null_scores", include_null_scores),
                        ("include_input", include_input), ("include_all_events", include_all_events)]:
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be boolean")
    spec = resolve_positive(result, positive_col=positive_col, positive_threshold=positive_threshold,
                            positive_qc_only=positive_qc_only, positive_min_score=positive_min_score, selection=selection if selection is not None else selected_events)
    if report or selection is not None or events is not None or len(result.events.columns):
        _, events = report_selection(result, selection, events, detail_n, positive=spec)
    delivery = delivery_record(result, spec, selection, events, INHERIT, plot_top_n)
    if selected_events is not None:
        delivery['pair_selections'] = pair_records(selected_events, spec)
    delivery.update(report_layout=4, report_enabled=bool(report), include_all_events=bool(include_all_events))
    path, _ = _paths(prefix)
    paths = bundle_paths(path, delivery, include_input=include_input)
    _validate_targets(paths.values())
    saved, saved_input = {}, None
    if path.exists():
        with h5py.File(path, 'r') as handle:
            saved = _read(handle['metadata/delivery']) if 'metadata/delivery' in handle else {}
            saved_input = _read(handle['metadata/input_archive']) if 'metadata/input_archive' in handle else None
        if saved and ((saved.get('selection') is not None and saved.get('selection') != delivery['selection']) or
                      (saved.get('events') and saved.get('events') != delivery['events']) or
                      (saved.get('pair_selections') is not None and saved.get('pair_selections') != delivery.get('pair_selections'))):
            raise ValueError("Existing independent reports require add_selection(..., replace=True), "
                             "add_event_reports(), or a new save prefix")
    managed = list(paths.values())
    if saved and not report:
        old_paths = bundle_paths(path, saved)
        managed.extend(value for key, value in old_paths.items()
                       if key in ('report', 'filtered_report', 'navigation') or key.startswith(('event_report:', 'pair_report:')))
    if saved_input and not include_input:
        managed.append(path.with_name(path.name[:-len('.cellmesh.h5')] + '.input.h5ad'))
    managed = list(dict.fromkeys(managed))
    _validate_targets(managed)
    if include_input and result._input_adata is None and not result.input_archive:
        raise ValueError("Input was not retained; run_cell_mesh(retain_input=True) is required")
    has_explanation = result.priors is not None and result.availability_results is not None
    if report and not has_explanation:
        raise ValueError("HTML reporting requires a result with captured priors and intermediates")
    tables = build_report_tables(result) if has_explanation else {}
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".cellmesh-save-", dir=path.parent) as temporary:
        staging = Path(temporary)
        input_archive = stage_input(result, staging / paths['input'].name) if include_input else None
        contents = render_bundle(result, paths, delivery, selection, tables=tables, adata=adata,
                                 input_storage={"archive": input_archive}) if report else {}
        with h5py.File(staging / path.name, "w") as handle:
            handle.attrs["format"] = "CELL MESH result"
            handle.attrs["schema_version"] = FORMAT_VERSION
            metadata = handle.create_group("metadata")
            _write(metadata, "package_version", __version__, include_null_scores=include_null_scores)
            _write(metadata, "parameters", result.parameters, include_null_scores=include_null_scores)
            _write(metadata, "delivery", delivery, include_null_scores=False)
            _write(metadata, "input_provenance", result.parameters.get("prior_inputs", {}),
                   include_null_scores=include_null_scores)
            # Additive schema-v1 metadata: old files need neither entry.
            _write(metadata, "expression_input", result.input_provenance,
                   include_null_scores=include_null_scores)
            _write(metadata, "input_archive", input_archive,
                   include_null_scores=include_null_scores)
            _write(metadata, "storage", {"raw_single_cell_matrix_included": False,
                   "input_companion_included": bool(include_input),
                   "include_null_scores": bool(include_null_scores), "explanation_available": has_explanation},
                   include_null_scores=include_null_scores)
            for group_name, names in [("results", RESULT_TABLES), ("sample_results", SAMPLE_TABLES)]:
                group = handle.create_group(group_name)
                for name in names:
                    _write(group, name, getattr(result, name), include_null_scores=include_null_scores)
            prior_group = handle.create_group("priors")
            prior_group.attrs["available"] = result.priors is not None
            if result.priors is not None:
                for prior_name in ["enzyme", "interaction"]:
                    _write(prior_group, prior_name, result.priors[prior_name], include_null_scores=include_null_scores)
            # Preserve the original result object layout as well as queryable report tables.
            _write(handle, "availability_results", result.availability_results,
                   include_null_scores=include_null_scores)
            for name, content in tables.items():
                group = handle.create_group(name)
                for table_name, frame in content.items():
                    _write(group, table_name, frame, include_null_scores=include_null_scores)
            optional = handle.create_group("optional")
            if include_null_scores and NULL_KEY in result.events.attrs:
                _write(optional, "permutation_null", result.events.attrs[NULL_KEY], include_null_scores=True)
        _write_delivery_csvs(result, delivery, selection, paths, staging)
        for key, content in contents.items():
            staged = _staged(staging, paths[key])
            staged.parent.mkdir(parents=True, exist_ok=True)
            staged.write_text(content, encoding='utf-8')
        _publish_bundle(staging, path, delivery, paths, managed)
    # Publish associations only after all requested files have succeeded. This
    # lets to_html()/load_input(backed='r') use the most recent saved bundle.
    result.input_archive = input_archive
    result.delivery = delivery
    result._result_path = path.absolute()
    return paths


def load_result(path):
    """Load a schema-versioned result; no code from the file is executed."""
    from .core import CellMeshResult
    with h5py.File(path, "r") as handle:
        if handle.attrs.get("format") != "CELL MESH result":
            raise ValueError("Not a CELL MESH result file")
        if handle.attrs.get("schema_version") != FORMAT_VERSION:
            raise ValueError("Unsupported CELL MESH result schema version")
        kwargs = {}
        for group_name, names in [("results", RESULT_TABLES), ("sample_results", SAMPLE_TABLES)]:
            for name in names:
                kwargs[name] = _read(handle[group_name][name])
        kwargs["parameters"] = _read(handle["metadata/parameters"])
        for field, key in [("input_provenance", "expression_input"), ("input_archive", "input_archive"), ("delivery", "delivery")]:
            if key in handle["metadata"]:
                kwargs[field] = _read(handle["metadata"][key])
        kwargs["availability_results"] = _read(handle["availability_results"])
        kwargs["priors"] = ({name: _read(handle["priors"][name]) for name in ["enzyme", "interaction"]}
                            if handle["priors"].attrs.get("available", False) else None)
        if "permutation_null" in handle["optional"]:
            kwargs["events"].attrs[NULL_KEY] = _read(handle["optional/permutation_null"])
    # Bind the relative companion filename to the loaded result's directory.
    # Loading results never opens or requires the companion file.
    kwargs["_result_path"] = Path(path).absolute()
    return CellMeshResult(**kwargs)


def _validate_targets(targets):
    for target in targets:
        if target.is_symlink() or (target.exists() and not target.is_file()):
            raise ValueError("Destinations must be regular files, not symlinks or directories")
        for parent in target.parents:
            if parent.is_symlink() or (parent.exists() and not parent.is_dir()):
                raise ValueError("Destination directories must not be symlinks or files")


def _staged(staging, target):
    return staging / target.absolute().relative_to(staging.parent.absolute())


def _publish(staging, paths, managed):
    """Publish nested files together, preserving normal I/O failure rollback."""
    _validate_targets(managed)
    backup = staging / 'backup'
    backup.mkdir()
    existing = {p for p in managed if p.exists()}
    for target in existing:
        destination = backup / target.absolute().relative_to(staging.parent.absolute())
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target, destination)
    changed, created_dirs = [], []
    try:
        for target in paths.values():
            missing = [parent for parent in target.parents if not parent.exists()]
            target.parent.mkdir(parents=True, exist_ok=True)
            created_dirs.extend(reversed(missing))
            os.replace(_staged(staging, target), target)
            changed.append(target)
        for target in managed:
            if target not in paths.values() and target in existing:
                target.unlink()
                changed.append(target)
    except Exception:
        for target in reversed(changed):
            if target in existing:
                original = backup / target.absolute().relative_to(staging.parent.absolute())
                os.replace(original, target)
            else:
                target.unlink(missing_ok=True)
        for parent in reversed(created_dirs):
            if parent.exists() and not any(parent.iterdir()):
                parent.rmdir()
        raise


def _publish_bundle(staging, path, delivery, paths, managed):
    """Publish reports and archive obsolete names in the same transaction."""
    from ._report_history import stage_legacy_archive
    archives, obsolete = stage_legacy_archive(staging, path, delivery, paths)
    targets = list(dict.fromkeys([*managed, *archives.values(), *obsolete]))
    _publish(staging, {**paths, **archives}, targets)
    return archives


def archive_legacy_reports(prefix):
    """Repair an already migrated bundle without loading or rewriting its results.

    Save/append perform this housekeeping automatically. This maintenance helper
    handles existing layout-4 bundles containing obsolete ``selected`` filenames.
    Originals are preserved in a ZIP under ``<prefix>.reports/history/``; current
    reports, CSVs, H5 and input archives are untouched. Repeated calls are no-ops.
    """
    path, _ = _paths(prefix)
    _validate_targets([path])
    with h5py.File(path, 'r') as handle:
        if (handle.attrs.get('format') != 'CELL MESH result' or
                handle.attrs.get('schema_version') != FORMAT_VERSION):
            raise ValueError('Not a supported CELL MESH result file')
        delivery = _read(handle['metadata/delivery']) if 'metadata/delivery' in handle else {}
    if delivery.get('report_layout', 0) < 4:
        raise ValueError('Migrate reports with add_filtered_report() before archiving legacy files')
    with tempfile.TemporaryDirectory(prefix='.cellmesh-history-', dir=path.parent) as directory:
        return _publish_bundle(Path(directory), path, delivery, {}, [])


def _write_delivery_csvs(result, delivery, selection, paths, staging):
    from .core import _result_table_for_csv
    from .positivity import positive_selection, selection_table
    if result.events.empty and not len(result.events.columns):
        from .selection import EVENT_KEYS
        positives = result.events.reindex(columns=EVENT_KEYS + ['cell_mesh_score', delivery['positive']['column']])
    else:
        positives = positive_selection(result, delivery['positive']).events
    frames = {"positive_csv": positives}
    if 'events_csv' in paths:
        frames['events_csv'] = result.events
    if selection is not None:
        frames['filtered_csv'] = selection_table(selection, delivery['positive'])
    if 'selected_csv' in paths:
        from ._pair_selection import pair_table
        frames['selected_csv'] = pair_table(result, delivery.get('pair_selections', []), delivery['positive'])
    for key, frame in frames.items():
        if key in paths:
            _result_table_for_csv(frame, key).to_csv(staging / paths[key].name, index=False, encoding='utf-8')


def _report_csv_paths(path, selection, *, include_all_events=False):
    stem = path.name[:-len('.report.html')] if path.name.endswith('.report.html') else path.stem
    if selection is not None:
        if stem.endswith(('.selected', '.filtered')):
            stem = stem.rsplit('.', 1)[0]
        return {'report': path, 'filtered_csv': path.with_name(stem + '.filtered_events.csv')}
    paths = {'report': path, 'positive_csv': path.with_name(stem + '.positive_events.csv')}
    if include_all_events:
        paths['events_csv'] = path.with_name(stem + '.all_events.csv')
    return paths


def protected_bundle_paths(result):
    from ._report_delivery import bundle_paths
    path = result._result_path
    stem = path.name[:-len('.cellmesh.h5')]
    protected = {value.resolve() for value in bundle_paths(path, result.delivery or {}, include_input=True).values()}
    protected.update(path.with_name(stem + suffix).resolve() for suffix in
                     ['.all_events.csv', '.positive_events.csv', '.selected_events.csv',
                      '.report.html', '.selected.report.html', '.filtered.report.html', '.filtered_events.csv',
                      '.metabolite_sensor_events.csv', '.cellmesh.h5', '.input.h5ad'])
    return protected, path.with_name(stem + '.reports').resolve()


def write_report_bundle(result, path, *, filtered_events=None, selected_events=None, events=None, selection=None, include_all_events=False,
                        detail_n=INHERIT, plot_top_n=20, adata=None, **positive_options):
    from .reporting import render_html_report, build_report_tables
    from .selection import report_selection
    from .positivity import resolve_positive, delivery_record
    from ._report_plots import validate_plot_top_n
    from ._report_delivery import relative, bundle_paths
    from ._pair_selection import resolve_scopes
    selection, selected_events = resolve_scopes(result, filtered_events, selected_events, selection)
    if events is not None:
        raise TypeError("events= is retired; use selected_events=result.select_events(metabolite_sensor_pairs=...)")
    if selected_events is not None:
        if selection is not None or include_all_events:
            raise ValueError("Use save() to combine general, filtered and selected outputs")
        from ._pair_reports import write_pair_bundle
        return write_pair_bundle(result, path, selected_events, plot_top_n=plot_top_n,
                                 adata=adata, **positive_options)
    path = Path(path)
    plot_top_n = validate_plot_top_n(plot_top_n)
    if not isinstance(include_all_events, (bool, np.bool_)):
        raise TypeError("include_all_events must be boolean")
    if result.priors is None or result.availability_results is None:
        raise ValueError("HTML reporting requires captured priors and intermediates")
    if path.suffix.lower() != '.html':
        raise ValueError("Report destination must have an .html suffix")
    spec = resolve_positive(result, selection=selection, **positive_options)
    report_selection(result, selection, None, detail_n, positive=spec)
    if include_all_events and selection is not None:
        raise ValueError("include_all_events=True requires a general HTML export; "
                         "use save() to include it alongside selection or event reports")
    delivery = delivery_record(result, spec, selection, [], INHERIT, plot_top_n)
    delivery.update(report_layout=4, report_enabled=True, include_all_events=bool(include_all_events))
    paths = _report_csv_paths(path, selection, include_all_events=include_all_events)
    general_report = None
    if result._result_path is not None:
        protected_paths = bundle_paths(result._result_path, result.delivery or {}, include_input=True)
        protected, report_dir = protected_bundle_paths(result)
        if any(value.resolve() in protected or value.resolve().is_relative_to(report_dir) for value in paths.values()):
            raise ValueError("Report companions overlap the saved bundle; use append methods or choose a new prefix")
        general = protected_paths.get('report')
        if general is not None and general.is_file():
            general_report = general
    _validate_targets(paths.values())
    tables = build_report_tables(result) if selection is None else None
    contents = {'report': render_html_report(result, selection=selection, plot_top_n=plot_top_n,
        tables=tables, delivery=delivery, adata=adata,
        general_report=relative(general_report, path.parent) if general_report else None,
        csv_paths={k: v.name for k, v in paths.items() if k.endswith('csv')})}
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.cellmesh-report-', dir=path.parent) as directory:
        staging = Path(directory)
        _write_delivery_csvs(result, delivery, selection, paths, staging)
        for key, content in contents.items():
            staged = _staged(staging, paths[key])
            staged.parent.mkdir(parents=True, exist_ok=True)
            staged.write_text(content, encoding='utf-8')
        _publish(staging, paths, list(paths.values()))
    return paths


def add_selection(result, selection, *, replace=False, detail_n=INHERIT, plot_top_n=INHERIT, adata=None):
    from .selection import EventSelection
    from ._pair_selection import resolve_scopes
    resolve_scopes(result, filtered_events=selection)
    return _append_reports(result, selection=selection, replace=replace, detail_n=detail_n,
                           plot_top_n=plot_top_n, adata=adata)


def add_event_reports(result, selected_events, *, replace=False, adata=None):
    from ._pair_selection import resolve_scopes
    resolve_scopes(result, selected_events=selected_events)
    return _append_reports(result, selected_events=selected_events, replace=replace, adata=adata)


def _append_reports(result, *, selection=None, selected_events=None, replace=False, detail_n=INHERIT,
                    plot_top_n=INHERIT, adata=None):
    from .selection import report_selection, normalize_event_requests
    from .positivity import resolve_positive, delivery_record, restore_selection
    from ._report_plots import validate_plot_top_n
    from ._report_delivery import bundle_paths, render_bundle
    from ._report_browser import event_id
    if not isinstance(replace, (bool, np.bool_)):
        raise TypeError("replace must be boolean")
    if result._result_path is None or not result._result_path.is_file():
        raise ValueError("Save the result before adding reports")
    path = result._result_path
    _validate_targets([path])
    from ._pair_selection import pair_records, restore_pairs
    requested = pair_records(selected_events, resolve_positive(result))
    with h5py.File(path, 'r') as handle:
        saved = _read(handle['metadata/delivery']) if 'delivery' in handle['metadata'] else None
        if saved != result.delivery:
            raise ValueError("Saved delivery changed; reload the result before updating")
        original = _read(handle['results/events'])
        if not pd.DataFrame(original).equals(pd.DataFrame(result.events)) or (
                _pack(_read(handle['metadata/parameters'])) != _pack(result.parameters)):
            raise ValueError("Saved scientific result differs; save a new bundle before updating")
        if selection is not None and (saved or {}).get('selection') is not None and not replace:
            raise ValueError("A selection already exists; pass replace=True to replace it")
        active_selection = selection if selection is not None else restore_selection(result)
        existing_events = normalize_event_requests(result, (saved or {}).get('events'))
        existing = {record['id']: record for record in (saved or {}).get('pair_selections') or []}
        if not replace and any(record['id'] in existing for record in requested):
            raise ValueError("A metabolite–sensor report already exists; pass replace=True to replace it")
        for record in requested:
            existing[record['id']] = record
        spec = resolve_positive(result)
        plot_top_n = validate_plot_top_n((saved or {}).get('plot_top_n', 20) if plot_top_n is INHERIT else plot_top_n)
        report_selection(result, active_selection, existing_events, detail_n, positive=spec)
        report_enabled = (saved or {}).get('report_enabled', True)
        if selected_events is not None and not report_enabled:
            raise ValueError("HTML output is disabled; save with report=True before adding event reports")
        delivery = delivery_record(result, spec, active_selection, existing_events, INHERIT, plot_top_n)
        if selected_events is not None or 'pair_selections' in (saved or {}):
            delivery['pair_selections'] = list(existing.values())
        delivery.update(report_layout=4, report_enabled=report_enabled,
                        include_all_events=bool((saved or {}).get('include_all_events', False)))
        all_paths = bundle_paths(path, delivery, include_input=bool(result.input_archive))
        paths = {'result': path}
        migration = (saved or {}).get('report_layout', 0) < 4
        if active_selection is not None and (selection is not None or migration):
            paths['filtered_csv'] = all_paths['filtered_csv']
        if selected_events is not None:
            paths['selected_csv'] = all_paths['selected_csv']
        if report_enabled:
            paths['navigation'] = all_paths['navigation']
            if active_selection is not None and (selection is not None or requested or migration):
                paths['filtered_report'] = all_paths['filtered_report']
            # Legacy bundles predate the explicit CSV choice. Refresh their
            # general links once, leaving any existing full CSV untouched.
            if (saved or {}).get('report_layout', 0) < 3 or 'include_all_events' not in (saved or {}) or not all_paths['report'].exists():
                paths['report'] = all_paths['report']
            for record in requested:
                key = 'pair_report:' + record['id']
                paths[key] = all_paths[key]
        needs_tables = 'report' in paths or any(key.startswith(('event_report:', 'pair_report:')) for key in paths)
        if needs_tables and not {'explanation', 'coverage'} <= set(handle):
            raise ValueError("Saved result has no explanation tables")
        tables = {group: {name: _read(node) for name, node in handle[group].items()}
                  for group in ['explanation', 'coverage']} if needs_tables else None
    _validate_targets(all_paths.values())
    from .core import _result_table_for_csv
    from .positivity import positive_selection
    for key, frame in [('events_csv', result.events), ('positive_csv', positive_selection(result, spec).events)]:
        if key not in all_paths:
            continue
        if not all_paths[key].exists():
            paths[key] = all_paths[key]
            continue
        expected = _result_table_for_csv(frame, key).to_csv(index=False).encode('utf-8')
        actual = hashlib.sha256()
        with all_paths[key].open('rb') as stream:
            for chunk in iter(lambda: stream.read(4 * 1024**2), b''):
                actual.update(chunk)
        if actual.digest() != hashlib.sha256(expected).digest():
            raise ValueError("Existing general CSV differs from this result/criterion; save a new bundle")
    evidence_path = None
    if selection is None and active_selection is not None and 'filtered_report' in paths:
        evidence_path = all_paths['filtered_report']
        if not evidence_path.exists():
            evidence_path = bundle_paths(path, saved).get('filtered_report', evidence_path)
    contents = render_bundle(result, all_paths, delivery, active_selection, requested=paths,
                             tables=tables, adata=adata, selection_evidence_path=evidence_path) if report_enabled else {}
    with tempfile.TemporaryDirectory(prefix='.cellmesh-update-', dir=path.parent) as directory:
        staging = Path(directory)
        shutil.copy2(path, staging / path.name)
        with h5py.File(staging / path.name, 'r+') as handle:
            if 'delivery' in handle['metadata']:
                del handle['metadata/delivery']
            _write(handle['metadata'], 'delivery', delivery, include_null_scores=False)
        _write_delivery_csvs(result, delivery, active_selection, paths, staging)
        for key, content in contents.items():
            staged = _staged(staging, paths[key])
            staged.parent.mkdir(parents=True, exist_ok=True)
            staged.write_text(content, encoding='utf-8')
        _publish_bundle(staging, path, delivery, paths, list(paths.values()))
    result.delivery = delivery
    return all_paths
