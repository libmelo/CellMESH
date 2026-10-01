"""Delivery-level positive calls; never recompute scores or multiple testing."""
from __future__ import annotations

import copy

import numpy as np
import pandas as pd


DEFAULT_POSITIVE = dict(column="fdr_sensor_type", threshold=0.05, qc_only=False, min_score=None)
INHERIT = object()
STATISTIC_LABELS = {
    "fdr_sensor_type": "BH 校正 p 值（按受体类型）",
    "fdr_global": "BH 校正 p 值（全局）",
    "perm_pvalue": "置换 p 值（未校正）",
}
STATUS_LABELS = {"positive": "阳性", "not_positive": "未达标准", "not_evaluated": "无法判定"}


def resolve_positive(result, *, positive_col=None, positive_threshold=None,
                     positive_qc_only=None, positive_min_score=INHERIT, selection=None):
    from .plotting import _plot_number
    if selection is not None:
        from .selection import EventSelection
        if not isinstance(selection, EventSelection):
            raise TypeError("selection must be returned by this result's filter_events()")
    spec = copy.deepcopy((result.delivery or {}).get("positive", DEFAULT_POSITIVE))
    if not result.delivery and selection is not None:
        spec = copy.deepcopy(selection.criteria.get("positive", spec))
    for key, value in [("column", positive_col), ("threshold", positive_threshold),
                       ("qc_only", positive_qc_only), ("min_score", positive_min_score)]:
        if (key == "min_score" and value is not INHERIT) or (key != "min_score" and value is not None):
            spec[key] = value
    if spec["column"] not in STATISTIC_LABELS:
        raise ValueError("positive_col must be fdr_sensor_type, fdr_global or perm_pvalue")
    spec["threshold"] = _plot_number(spec["threshold"], "positive_threshold", upper=1)
    if not isinstance(spec["qc_only"], (bool, np.bool_)):
        raise TypeError("positive_qc_only must be boolean")
    spec["qc_only"] = bool(spec["qc_only"])
    if spec["min_score"] is not None:
        spec["min_score"] = _plot_number(spec["min_score"], "positive_min_score", upper=1)
    return spec


def positive_calls(result, spec):
    from .plotting import _plot_numeric_series, _plot_qc_series
    frame = result.events.reset_index(drop=True)
    values = (_plot_numeric_series(frame[spec["column"]], spec["column"], upper=1, probability=True)
              if spec["column"] in frame else pd.Series(np.nan, index=frame.index))
    known = values.notna()
    passed = values.le(spec["threshold"])
    if spec["qc_only"]:
        if "passes_min_cells" not in frame:
            raise ValueError("positive_qc_only=True requires passes_min_cells")
        qc = _plot_qc_series(frame.passes_min_cells)
        known &= qc.notna()
        passed &= qc.eq(True).fillna(False)
    if spec["min_score"] is not None:
        scores = _plot_numeric_series(frame.cell_mesh_score, "cell_mesh_score", upper=1)
        known &= scores.notna()
        passed &= scores.ge(spec["min_score"])
    status = pd.Series(np.where(~known, "not_evaluated", np.where(passed, "positive", "not_positive")),
                       index=frame.index, name="positive_status")
    counts = dict(n_total=len(frame), n_positive=int(status.eq("positive").sum()),
                  n_not_positive=int(status.eq("not_positive").sum()),
                  n_not_evaluated=int(status.eq("not_evaluated").sum()),
                  status="empty" if frame.empty else "not_evaluated" if not known.any() else
                         "partial" if not known.all() else "evaluated")
    return status, counts


def positive_description(spec):
    text = f'{STATISTIC_LABELS[spec["column"]]} ≤ {spec["threshold"]:g}'
    if spec["qc_only"]:
        text += "；通过细胞数质控"
    if spec["min_score"] is not None:
        text += f'；通讯评分 ≥ {spec["min_score"]:g}'
    return text


def positive_selection(result, spec):
    from .selection import filter_events
    return filter_events(result, source="positive", positive_col=spec["column"],
                         positive_threshold=spec["threshold"], positive_qc_only=spec["qc_only"],
                         positive_min_score=spec["min_score"])


def selection_table(selection, spec):
    """Add delivery status to a copy; scientific result/selection rows stay intact."""
    selection._validate()
    status, _ = positive_calls(selection._result, spec)
    frame = selection.events
    if "positive_status" in frame:
        raise ValueError("positive_status is reserved for delivery annotations")
    frame["positive_status"] = status.iloc[selection._positions].to_numpy()
    return frame


def delivery_record(result, spec, selection, events, detail_n, plot_top_n):
    from .selection import EVENT_KEYS
    _, counts = positive_calls(result, spec)
    selected = None
    if selection is not None:
        selection._validate(result)
        if selection.criteria["source"] == "positive" and selection.criteria["positive"] != spec:
            raise ValueError("selection uses a different positive criterion; reselect with the delivery criterion or source='all'")
        selected = dict(criteria=selection.criteria, summary=selection.summary,
                        keys=selection.events[EVENT_KEYS].to_dict("records"))
    return dict(version=1, positive=spec, positive_summary=counts, selection=selected,
                events=copy.deepcopy(events or []), plot_top_n=plot_top_n)


def restore_selection(result):
    record = (result.delivery or {}).get("selection")
    if record is None:
        return None
    from .selection import filter_events, EVENT_KEYS
    criteria = copy.deepcopy(record["criteria"])
    spec = criteria.pop("positive")
    focus = filter_events(result, **criteria, positive_col=spec["column"], positive_threshold=spec["threshold"],
                          positive_qc_only=spec["qc_only"], positive_min_score=spec["min_score"])
    if focus.events[EVENT_KEYS].to_dict("records") != record["keys"]:
        raise ValueError("Saved selection no longer matches result events")
    return focus
