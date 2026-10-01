"""Explicit molecular-pair requests, separate from independent list filters."""
from __future__ import annotations

import copy
import hashlib
import json

import numpy as np
import pandas as pd

from .selection import EVENT_KEYS, EventSelection, _labels, filter_events


def pair_id(pair):
    value = json.dumps(list(pair), ensure_ascii=True, separators=(",", ":"))
    return "pair_" + hashlib.sha256(value.encode()).hexdigest()[:24]


class SelectedEvents(EventSelection):
    """A validated, explicit pair request; even empty requested pairs persist."""

    @property
    def pairs(self):
        return [tuple(pair) for pair in self.criteria["metabolite_sensor_pairs"]]

    def to_html(self, path, **options):
        return self._result.to_html(path, selected_events=self, **options)


def select_events(result, *, metabolite_sensor_pairs=None, **filters):
    if metabolite_sensor_pairs is None or "hmdb_ids" in filters or "sensor_genes" in filters:
        raise TypeError("select_events() requires metabolite_sensor_pairs=[(HMDB, sensor), ...]. "
                        "For independent hmdb_ids/sensor_genes lists use filter_events().")
    if isinstance(metabolite_sensor_pairs, (str, bytes, dict)):
        raise TypeError("metabolite_sensor_pairs must be a sequence of (HMDB, sensor) pairs")
    pairs = []
    for item in metabolite_sensor_pairs:
        if isinstance(item, (str, bytes, dict)) or not hasattr(item, "__len__") or len(item) != 2:
            raise ValueError("Each metabolite_sensor_pairs entry needs exactly an HMDB ID and sensor")
        pair = (_labels([item[0]], "hmdb_ids", hmdb=True)[0],
                _labels([item[1]], "sensor_genes")[0])
        if pair in pairs:
            raise ValueError("Duplicate metabolite_sensor_pairs entry: " + repr(pair))
        pairs.append(pair)
    base = filter_events(result, hmdb_ids=list(dict.fromkeys(p[0] for p in pairs)),
                         sensor_genes=list(dict.fromkeys(p[1] for p in pairs)), **filters)
    frame = base.events
    keep = pd.MultiIndex.from_frame(frame[["hmdb_id", "sensor_gene"]]).isin(pairs)
    positions = base._positions[keep]
    criteria, summary = base.criteria, base.summary
    criteria.update(hmdb_ids=None, sensor_genes=None, metabolite_sensor_pairs=pairs)
    summary.update(n_selected=len(positions), n_excluded=summary["n_total"] - len(positions),
                   n_requested_pairs=len(pairs))
    summary["steps"].append(dict(filter="metabolite_sensor_pairs", applied=True,
                                before=len(frame), remaining=len(positions),
                                excluded=len(frame) - len(positions), excluded_missing=0))
    summary["empty_reason"] = "No events satisfy the requested pairs and conditions" if not len(positions) else None
    return SelectedEvents(result, base._source, positions, criteria, summary)


def resolve_scopes(result, filtered_events=None, selected_events=None, selection=None):
    """Reject mixed scopes before staging any output. Legacy selection is tier 2."""
    if selection is not None:
        import warnings
        if filtered_events is not None:
            raise ValueError("Use filtered_events only; do not also pass legacy selection")
        warnings.warn("selection= is deprecated; use filtered_events= from filter_events().",
                      FutureWarning, stacklevel=3)
        filtered_events = selection
    if filtered_events is not None:
        if not isinstance(filtered_events, EventSelection) or isinstance(filtered_events, SelectedEvents):
            raise TypeError("filtered_events must come from filter_events()")
        filtered_events._validate(result)
    if selected_events is not None:
        if not isinstance(selected_events, SelectedEvents):
            raise TypeError("selected_events must come from select_events(metabolite_sensor_pairs=...)")
        selected_events._validate(result)
    return filtered_events, selected_events


def pair_records(selected, spec):
    if selected is None:
        return []
    selected._validate()
    if selected.criteria["source"] == "positive" and selected.criteria["positive"] != spec:
        raise ValueError("selected_events uses a different positive criterion; reselect with matching options")
    records = []
    for pair in selected.pairs:
        criteria = selected.criteria
        criteria["metabolite_sensor_pairs"] = [pair]
        frame = selected.events
        frame = frame.loc[frame.hmdb_id.eq(pair[0]) & frame.sensor_gene.eq(pair[1])]
        records.append(dict(id=pair_id(pair), pair=list(pair), criteria=criteria,
                            keys=frame[EVENT_KEYS].to_dict("records")))
    return records


def restore_pairs(result, records):
    selections = []
    for record in records:
        criteria = copy.deepcopy(record["criteria"])
        spec = criteria.pop("positive")
        criteria.pop("hmdb_ids", None)
        criteria.pop("sensor_genes", None)
        selected = select_events(result, **criteria, positive_col=spec["column"],
                                 positive_threshold=spec["threshold"],
                                 positive_qc_only=spec["qc_only"], positive_min_score=spec["min_score"])
        if len(selected.pairs) != 1 or pair_id(selected.pairs[0]) != record["id"]:
            raise ValueError("Saved pair identity does not match its criteria")
        if selected.events[EVENT_KEYS].to_dict("records") != record["keys"]:
            raise ValueError("Saved pair selection no longer matches the scientific result")
        selections.append(selected)
    return selections


def pair_table(result, records, spec):
    from .positivity import selection_table
    frames = []
    for record, selected in zip(records, restore_pairs(result, records)):
        frame = selection_table(selected, spec)
        frame.insert(0, "pair_id", record["id"])
        frames.append(frame)
    if frames:
        return pd.concat(frames, ignore_index=True)
    frame = result.events.iloc[:0].copy()
    frame["positive_status"] = pd.Series(dtype=str)
    frame.insert(0, "pair_id", pd.Series(dtype=str))
    return frame
