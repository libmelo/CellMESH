"""Reusable, non-scoring selections of aggregate communication events."""
from __future__ import annotations

import copy
import hashlib
import json
from numbers import Integral
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from .positivity import INHERIT
from .database import _normalize_hmdb_id

from .plotting import (_canonical_event_table, _plot_labels, _plot_number,
                       _plot_numeric_series, _plot_qc_series, _require_hmdb_id)


EVENT_KEYS = ["sender", "receiver", "hmdb_id", "sensor_gene"]


def event_identity(values):
    """Canonical complete identity, shared by selection, plots and reports."""
    return tuple(_normalize_hmdb_id(values[key]) if key == "hmdb_id" else
                 str(values[key]).strip() for key in EVENT_KEYS)


def event_id(values):
    """Stable across event order, plot layout and report pagination."""
    identity = json.dumps(event_identity(values), ensure_ascii=True, separators=(",", ":"))
    return "event_" + hashlib.sha256(identity.encode("utf-8")).hexdigest()[:24]


def _values(frame):
    # Exclude attrs, which can contain a very large permutation null matrix.
    return pd.DataFrame(frame).copy(deep=True)


def _labels(values, name, *, hmdb=False):
    if values is None:
        return None
    labels = _plot_labels(values, name)
    if hmdb:
        labels = [_require_hmdb_id(value) for value in labels]
        if len(set(labels)) != len(labels):
            raise ValueError(f"{name} contains duplicate normalized HMDB IDs")
    return labels


def _known_labels(result, work):
    cells = set(work.sender) | set(work.receiver)
    availability = result.availability_results or {}
    units = [availability, *availability.get("availability_by_sample", {}).values()]
    for unit in units:
        pb = unit.get("pseudobulk")
        if isinstance(pb, pd.DataFrame):
            cells.update(str(value).strip() for value in pb.index)
    metabolites, sensors = set(work.hmdb_id), set(work.sensor_gene)
    for table in (result.priors or {}).values():
        if "hmdb_id" in table:
            metabolites.update(_require_hmdb_id(value) for value in table.hmdb_id.dropna())
        if "sensor_gene" in table:
            sensors.update(str(value).strip() for value in table.sensor_gene.dropna())
    return {"sender_labels": cells, "receiver_labels": cells,
            "hmdb_ids": metabolites, "sensor_genes": sensors}


class EventSelection:
    """An event snapshot linked to its original result, with defensive accessors.

    Create with ``result.filter_events()``. Public tables/metadata are copies;
    edits to them do not change the selection. Changes to the source event table
    or scientific parameters invalidate the selection and require reselection.
    """

    def __init__(self, result, source, positions, criteria, summary):
        from .result_io import _pack
        self._result = result
        self._source = source
        self._positions = np.asarray(positions, dtype=int)
        self._criteria = copy.deepcopy(criteria)
        self._summary = copy.deepcopy(summary)
        self._parameters = _pack(result.parameters)

    def _validate(self, result=None):
        from .result_io import _pack
        if result is not None and result is not self._result:
            raise ValueError("selection belongs to another result; create a scope from this result")
        if (not self._source.equals(self._result.events)
                or self._parameters != _pack(self._result.parameters)):
            raise ValueError("Source result changed; recreate the scope with filter_events() or select_events()")

    @property
    def events(self):
        self._validate()
        return self._source.iloc[self._positions].copy().reset_index(drop=True)

    @property
    def criteria(self):
        return copy.deepcopy(self._criteria)

    @property
    def summary(self):
        self._validate()
        return copy.deepcopy(self._summary)

    def to_html(self, path, *, detail_n=INHERIT, plot_top_n=20, **positive_options):
        """Write selected-event HTML and its selected CSV companion.

        This writes tier 2 only. ``plot_top_n`` limits comparison figures;
        explicit select_events(metabolite_sensor_pairs=...) requests define tier 3.
        Pass ``adata`` for grouped expression violins, or retain input at scoring.
        """
        self._validate()
        return self._result.to_html(path, filtered_events=self, detail_n=detail_n, plot_top_n=plot_top_n,
                                    **positive_options)

    def to_csv(self, path):
        """Atomically write all selected rows to exactly one UTF-8 CSV file."""
        from .core import _result_table_for_csv
        from .positivity import selection_table
        frame = _result_table_for_csv(selection_table(self, self.criteria["positive"]), "selected_events")
        path = Path(path)
        if path.is_symlink() or (path.exists() and not path.is_file()):
            raise ValueError("Selected CSV destination must be a regular file, not a symlink or directory")
        if self._result._result_path is not None:
            from .result_io import protected_bundle_paths
            protected, report_dir = protected_bundle_paths(self._result)
            if path.resolve() in protected or path.resolve().is_relative_to(report_dir):
                raise ValueError("Selected CSV cannot overwrite a standard result delivery file")
        with tempfile.TemporaryDirectory(prefix=".cellmesh-selected-", dir=path.parent) as temporary:
            staged = Path(temporary) / "selected.csv"
            frame.to_csv(staged, index=False, encoding="utf-8")
            os.replace(staged, path)
        return path


def filter_events(result, *, sender_labels=None, receiver_labels=None, cell_pairs=None,
                  hmdb_ids=None, sensor_genes=None, max_fdr=None,
                  fdr_col="fdr_sensor_type", min_cell_mesh_score=None, qc_only=False,
                  source="positive", max_perm_pvalue=None, positive_col=None,
                  positive_threshold=None, positive_qc_only=None, positive_min_score=INHERIT):
    """AND across filters, OR within a list; None is unrestricted, [] matches none.

    Output contains unchanged aggregate event values. Missing FDR/score is kept
    unless its threshold is requested. FDR values retain their original testing
    family; they are never recomputed on a selected subset.
    """
    if fdr_col not in ("fdr_sensor_type", "fdr_global"):
        raise ValueError("fdr_col must be 'fdr_sensor_type' or 'fdr_global'")
    if not isinstance(qc_only, (bool, np.bool_)):
        raise TypeError("qc_only must be boolean")
    from .positivity import resolve_positive, positive_calls
    if source not in ("positive", "all"):
        raise ValueError("source must be 'positive' or 'all'")
    scope = source
    spec = resolve_positive(result, positive_col=positive_col, positive_threshold=positive_threshold,
                            positive_qc_only=positive_qc_only, positive_min_score=positive_min_score)
    criteria = dict(sender_labels=_labels(sender_labels, "sender_labels"),
                    receiver_labels=_labels(receiver_labels, "receiver_labels"),
                    cell_pairs=None, hmdb_ids=_labels(hmdb_ids, "hmdb_ids", hmdb=True),
                    sensor_genes=_labels(sensor_genes, "sensor_genes"),
                    max_fdr=None if max_fdr is None else _plot_number(max_fdr, "max_fdr", upper=1),
                    fdr_col=fdr_col,
                    min_cell_mesh_score=None if min_cell_mesh_score is None else
                    _plot_number(min_cell_mesh_score, "min_cell_mesh_score", upper=1),
                    qc_only=bool(qc_only), source=scope, positive=spec,
                    max_perm_pvalue=None if max_perm_pvalue is None else
                    _plot_number(max_perm_pvalue, "max_perm_pvalue", upper=1))
    if cell_pairs is not None:
        if isinstance(cell_pairs, (str, bytes)):
            raise TypeError("cell_pairs must be a sequence of (sender, receiver) pairs")
        pairs = []
        for pair in cell_pairs:
            if isinstance(pair, (str, bytes)) or not hasattr(pair, "__len__") or len(pair) != 2:
                raise ValueError("Each cell_pairs entry must contain exactly a sender and receiver")
            # Self pairs are valid, so normalize their two endpoints separately.
            pairs.append((_labels([pair[0]], "sender_labels")[0],
                          _labels([pair[1]], "receiver_labels")[0]))
        if len(set(pairs)) != len(pairs):
            raise ValueError("cell_pairs contains duplicate normalized pairs")
        criteria["cell_pairs"] = pairs
    source = _values(result.events)
    required = set(EVENT_KEYS + ["cell_mesh_score"])
    if missing := required.difference(source.columns):
        raise KeyError("Missing event columns: " + ", ".join(sorted(missing)))
    if "sample" in source:
        raise ValueError("filter_events() requires aggregate events, not sample-level rows")
    work = _canonical_event_table(source).reset_index(drop=True)
    if work.duplicated(EVENT_KEYS).any():
        raise ValueError("Event identities must be unique: sender, receiver, hmdb_id, sensor_gene")
    known = _known_labels(result, work)
    for name, choices in known.items():
        requested = criteria[name]
        if requested is not None and (unknown := set(requested).difference(choices)):
            raise ValueError(f"Unknown {name}: {sorted(unknown)}")
    if criteria["cell_pairs"] is not None:
        endpoints = {label for pair in criteria["cell_pairs"] for label in pair}
        if unknown := endpoints.difference(known["sender_labels"]):
            raise ValueError(f"Unknown cell_pairs endpoints: {sorted(unknown)}")

    work["cell_mesh_score"] = _plot_numeric_series(work.cell_mesh_score, "cell_mesh_score", upper=1)
    if fdr_col in work:
        work[fdr_col] = _plot_numeric_series(work[fdr_col], fdr_col, upper=1, probability=True)
    fdr_available = fdr_col in work and work[fdr_col].notna().any()
    if criteria["max_fdr"] is not None and not fdr_available:
        raise ValueError(f"No available {fdr_col}; calculate FDR or omit max_fdr")
    if criteria["max_perm_pvalue"] is not None:
        if "perm_pvalue" not in work or not work.perm_pvalue.notna().any():
            raise ValueError("No available perm_pvalue; calculate p-values or omit max_perm_pvalue")
        work["perm_pvalue"] = _plot_numeric_series(work.perm_pvalue, "perm_pvalue", upper=1, probability=True)
    if criteria["qc_only"] and "passes_min_cells" not in work:
        raise ValueError("qc_only=True requires passes_min_cells")
    if "passes_min_cells" in work:
        work["passes_min_cells"] = _plot_qc_series(work.passes_min_cells)
    keep = pd.Series(True, index=work.index)
    steps = []

    def apply(name, active, mask=None, missing=None):
        before = int(keep.sum())
        absent = int((keep & missing).sum()) if active and missing is not None else 0
        if active:
            keep.loc[:] = keep & mask.fillna(False)
        after = int(keep.sum())
        steps.append(dict(filter=name, applied=bool(active), before=before, remaining=after,
                          excluded=before - after, excluded_missing=absent))

    status, positive_summary = positive_calls(result, spec)
    apply("positive_source", scope == "positive", status.eq("positive"), status.eq("not_evaluated"))

    for name, column in [("sender_labels", "sender"), ("receiver_labels", "receiver"),
                         ("hmdb_ids", "hmdb_id"), ("sensor_genes", "sensor_gene")]:
        values = criteria[name]
        apply(name, values is not None, work[column].isin(values or []))
    pair_keys = pd.Series(list(work[["sender", "receiver"]].itertuples(index=False, name=None)),
                          index=work.index, dtype=object)
    apply("cell_pairs", criteria["cell_pairs"] is not None, pair_keys.isin(criteria["cell_pairs"] or []))
    apply("max_fdr", criteria["max_fdr"] is not None,
          work[fdr_col].le(criteria["max_fdr"]) if criteria["max_fdr"] is not None else None,
          work[fdr_col].isna() if fdr_col in work else None)
    apply("max_perm_pvalue", criteria["max_perm_pvalue"] is not None,
          work.perm_pvalue.le(criteria["max_perm_pvalue"]) if criteria["max_perm_pvalue"] is not None else None,
          work.perm_pvalue.isna() if "perm_pvalue" in work else None)
    apply("min_cell_mesh_score", criteria["min_cell_mesh_score"] is not None,
          work.cell_mesh_score.ge(criteria["min_cell_mesh_score"]) if criteria["min_cell_mesh_score"] is not None else None,
          work.cell_mesh_score.isna())
    apply("qc_only", criteria["qc_only"], work.passes_min_cells.eq(True)
          if "passes_min_cells" in work else None,
          work.passes_min_cells.isna() if "passes_min_cells" in work else None)
    rank_col = spec["column"] if scope == "positive" else "perm_pvalue" if max_perm_pvalue is not None else fdr_col
    rank_available = rank_col in work and work[rank_col].notna().any()
    if rank_available:
        work[rank_col] = _plot_numeric_series(work[rank_col], rank_col, upper=1, probability=True)
    order = ([rank_col] if rank_available else []) + ["cell_mesh_score"] + EVENT_KEYS
    ascending = ([True] if rank_available else []) + [False] + [True] * len(EVENT_KEYS)
    ranked = work.loc[keep].sort_values(order, ascending=ascending, na_position="last", kind="stable")
    summary = dict(n_total=len(source), n_selected=len(ranked), n_excluded=len(source) - len(ranked),
                   source=scope, positive_summary=positive_summary,
                   steps=steps, order_by=order, ascending=ascending, fdr_available=bool(fdr_available),
                   missing_fdr=int(work[fdr_col].isna().sum()) if fdr_col in work else len(work),
                   missing_score=int(work.cell_mesh_score.isna().sum()),
                   empty_reason=("No source events" if source.empty else "No events satisfy all conditions")
                   if ranked.empty else None)
    return EventSelection(result, source, ranked.index.to_numpy(), criteria, summary)


def normalize_event_requests(result, events):
    """Validate explicit, complete event identities without expanding a selection."""
    from collections.abc import Mapping, Sequence
    if events is None:
        return []
    if isinstance(events, (str, bytes, Mapping)) or not isinstance(events, Sequence):
        raise TypeError("events must be a sequence of complete event dictionaries")
    normalized, seen = [], set()
    for event in events:
        if not isinstance(event, Mapping):
            raise TypeError("Each event must be a dictionary with sender, hmdb_id, sensor_gene, receiver")
        missing = set(EVENT_KEYS) - set(event)
        if missing:
            raise ValueError("Incomplete event; missing fields: " + ", ".join(sorted(missing)))
        extra = set(event) - set(EVENT_KEYS) - {"sample"}
        if extra:
            raise ValueError("Unknown event fields: " + ", ".join(sorted(extra)))
        item = {}
        for key in EVENT_KEYS:
            value = event[key]
            if not isinstance(value, str) or not value.strip() or "\x00" in value:
                raise ValueError(f"{key} must be a non-empty string without NUL")
            item[key] = value.strip()
        item["hmdb_id"] = _require_hmdb_id(item["hmdb_id"])
        if pd.isna(item["hmdb_id"]):
            raise ValueError("hmdb_id must be a valid HMDB identifier")
        identity = tuple(item[key] for key in EVENT_KEYS)
        if identity in seen:
            raise ValueError("events contains duplicate normalized event identities")
        seen.add(identity)
        if event.get("sample") is not None:
            sample = str(event["sample"]).strip()
            from .reporting import _units
            if sample not in [name for name, _ in _units(result)]:
                raise ValueError(f"Unknown sample {sample!r}; pooled results have no sample-level decomposition")
            item["sample"] = sample
        normalized.append(item)
    return normalized


def deprecated_detail_n(value):
    if value is INHERIT:
        return
    if value is not None and (isinstance(value, (bool, np.bool_))
                             or not isinstance(value, Integral) or value < 0):
        raise ValueError("detail_n must be a non-negative integer or None")
    import warnings
    warnings.warn("detail_n is deprecated and has no effect. Request separate single-event "
                  "reports with select_events(metabolite_sensor_pairs=...) and add_event_reports().",
                  DeprecationWarning, stacklevel=3)


def report_selection(result, selection, events, detail_n=INHERIT, positive=None):
    """Validate independent tier-2 and tier-3 scopes; never infer detail requests."""
    deprecated_detail_n(detail_n)
    from .positivity import resolve_positive, positive_selection
    if selection is not None:
        if not isinstance(selection, EventSelection):
            raise TypeError("selection must be returned by filter_events()")
        selection._validate(result)
    spec = resolve_positive(result, selection=selection) if positive is None else positive
    if selection is not None and selection.criteria["source"] == "positive" and selection.criteria["positive"] != spec:
        raise ValueError("selection uses a different positive criterion; reselect with matching positive options")
    return positive_selection(result, spec).events, normalize_event_requests(result, events)
