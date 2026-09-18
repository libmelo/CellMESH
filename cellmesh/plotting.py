"""Canonical plotting utilities for CELL MESH.

Plotting functions import optional visualization dependencies, such as
matplotlib, inside the function body so the core package import remains fast.
"""
from __future__ import annotations

from numbers import Integral
from typing import Any, Optional

import pandas as pd

from ._numerics import _reaction_activity, _with_numerical_diagnostics


_MIN_CELLS_QC_COLUMN = "passes_min_cells"


# 防回退：代谢物名称只能展示。代谢物—受体身份只有 HMDB ID + sensor_gene；
# sender/receiver（样本表再加 sample）用于区分观测记录，名称不得参与匹配、
# 分组、去重或排序主键。不能为兼容旧名称调用而回退到名称查找。
# 回归测试：tests/test_visual_identity.py。


def _plot_numeric_series(
    values: pd.Series, name: str, *, upper: Optional[float] = None,
    lower: Optional[float] = 0.0, allow_missing: bool = True,
    probability: bool = False,
) -> pd.Series:
    """Validate a numeric view without turning malformed input into biological NA."""
    import numpy as np
    from numbers import Number

    # 防回退：必须在阈值、排序、去重之前转换。errors="coerce" 只用于定位错误，
    # 不能将新产生的 NA 当作原本不可计算；复数/布尔值也不能被隐式转换为得分。
    # 真正下溢产生的有限零值合法，A04 的 report-and-continue 策略保持不变。
    # 回归测试：tests/test_plotting_input_contracts.py。
    invalid_type = values.map(
        lambda value: not pd.api.types.is_scalar(value)
        or isinstance(value, (bool, np.bool_, complex, np.complexfloating))
        or (not pd.isna(value) and not isinstance(value, (str, Number)))
    ).to_numpy(dtype=bool)
    original_missing = values.isna().to_numpy()
    # Mask unsupported dtypes through object storage, so a native complex column
    # cannot emit ComplexWarning (or lose its imaginary part) before our error.
    convertible = values.astype(object) if invalid_type.any() else values
    numeric = pd.to_numeric(convertible.mask(invalid_type), errors="coerce")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        array = numeric.to_numpy(dtype=float, na_value=np.nan)
    invalid = invalid_type | (~original_missing & ~np.isfinite(array))
    if not allow_missing:
        invalid |= original_missing
    if lower is not None:
        invalid |= array < lower
    if upper is not None:
        invalid |= array > upper
    if invalid.any():
        positions = np.flatnonzero(invalid)[:5]
        records = ", ".join(
            f"position {i} (index {values.index[i]!r}): {values.iloc[i]!r}" for i in positions
        )
        contract = "finite probabilities in [0, 1]" if probability else "finite real values"
        if not probability:
            if lower is not None:
                contract += f" >= {lower:g}"
            if upper is not None:
                contract += f" and <= {upper:g}"
        if allow_missing:
            contract += " or missing values"
        raise ValueError(f"{name} must contain {contract}; invalid records: {records}")
    # pandas multi-column sorting builds categorical indexes, which cannot use
    # float16. Promote the validated working copy BEFORE ranking/deduplication;
    # this preserves the stored values/NA, not precision lost by earlier casting.
    # Regression: tests/test_plotting_float16.py.
    if values.dtype == np.dtype("float16"):
        return values.astype(np.float64)
    # Preserve other valid native/nullable columns, including their NA dtype.
    if pd.api.types.is_numeric_dtype(values.dtype):
        return values.copy()
    return pd.Series(array, index=values.index, name=values.name)


def _plot_number(value: Any, name: str, *, lower: Optional[float] = 0.0,
                 upper: Optional[float] = None, positive: bool = False) -> float:
    result = float(_plot_numeric_series(pd.Series([value], dtype=object), name,
                                       lower=lower, upper=upper, allow_missing=False).iloc[0])
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _plot_range(values: Any, name: str, *, positive: bool = False) -> tuple[float, float]:
    try:
        if isinstance(values, (str, bytes)) or len(values) != 2:
            raise ValueError
        low, high = (_plot_number(v, name, lower=None, positive=positive) for v in values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite (minimum, maximum) pair"
                         + (" with positive values" if positive else "")) from exc
    if high < low:
        raise ValueError(f"{name} maximum must be >= minimum")
    return low, high


def _plot_labels(values: Any, name: str) -> list[str]:
    from .preprocess import _normalized_label_series

    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be a sequence of labels, not a string")
    labels = _normalized_label_series(pd.Series(list(values), dtype=object), name,
                                      reject_collisions=False).tolist()
    if len(set(labels)) != len(labels):
        raise ValueError(f"{name} must not contain duplicates after normalization")
    return labels


def _plot_thresholds(
    score: Any, pvalue: Any, fdr: Any, score_col: str,
) -> tuple[Optional[float], Optional[float], Optional[float]]:
    return (
        None if score is None else _plot_number(score, "min_cell_mesh_score",
                                                upper=1.0 if score_col == "cell_mesh_score" else None),
        None if pvalue is None else _plot_number(pvalue, "max_perm_pvalue", upper=1.0),
        None if fdr is None else _plot_number(fdr, "max_fdr", upper=1.0),
    )


def _plot_event_numbers(events: pd.DataFrame, score_col: str,
                        probability_cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = events.copy()
    out[score_col] = _plot_numeric_series(out[score_col], score_col,
                                         upper=1.0 if score_col == "cell_mesh_score" else None)
    for column in dict.fromkeys(probability_cols):
        if column in out:
            out[column] = _plot_numeric_series(out[column], column, upper=1.0, probability=True)
    excluded = out.loc[out[score_col].isna()].copy()
    excluded["plot_exclusion_reason"] = "missing_score"
    return out.loc[out[score_col].notna()].copy(), excluded


def _require_hmdb_id(value: Any) -> str:
    from .database import _normalize_hmdb_id

    identifier = _normalize_hmdb_id(value)
    if pd.isna(identifier):
        raise ValueError(
            "hmdb_id must not be empty or missing; provide hmdb_id, not a metabolite name"
        )
    return identifier


def _sensor_gene_id(value: Any) -> str:
    from .preprocess import _normalized_label_series

    return _normalized_label_series(pd.Series([value]), "sensor_gene").iloc[0]


def _canonical_event_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize matching fields on a copy; names are optional display metadata."""
    from .database import _normalize_hmdb_series
    from .preprocess import _normalized_label_series

    missing = {"hmdb_id", "sensor_gene"}.difference(frame.columns)
    if missing:
        raise KeyError("Missing required event columns: " + ", ".join(sorted(missing)))
    out = frame.copy()
    out["hmdb_id"] = _normalize_hmdb_series(out["hmdb_id"])
    if out["hmdb_id"].isna().any():
        raise ValueError("Event hmdb_id must not be empty or missing")
    out["sensor_gene"] = _normalized_label_series(
        out["sensor_gene"], "sensor_gene", reject_collisions=False,
    )
    for context in ("sender", "receiver", "sample"):
        if context in out:
            out[context] = _normalized_label_series(out[context], context)
    if "metabolite" not in out:
        out["metabolite"] = None
    return out


def _event_record_keys(unique_keys: Optional[list[str]]) -> list[str]:
    keys = ["sender", "receiver", "hmdb_id", "sensor_gene"]
    if unique_keys is None:
        return keys
    supplied = set(unique_keys)
    if not set(keys).issubset(supplied) or supplied.difference(keys + ["sample"]):
        raise ValueError(
            "unique_keys must include sender, receiver, hmdb_id and sensor_gene; "
            "only sample may be added. metabolite is display-only"
        )
    return keys + (["sample"] if "sample" in supplied else [])


def _metabolite_label(name: Any, hmdb_id: str) -> str:
    if name is None or pd.isna(name) or not str(name).strip() or str(name).strip() == hmdb_id:
        return hmdb_id
    return f"{str(name).strip()} ({hmdb_id})"


def _metabolite_display_names(frame: pd.DataFrame) -> dict[str, Any]:
    names = {}
    for identifier, name in zip(frame["hmdb_id"], frame["metabolite"]):
        if identifier not in names or pd.isna(names[identifier]) or not str(names[identifier]).strip():
            names[identifier] = name
    return names


def _event_selector(key: Any) -> tuple[str, str]:
    if isinstance(key, dict):
        hmdb_id, gene = key["hmdb_id"], key["sensor_gene"]
    elif isinstance(key, (tuple, list)) and len(key) in (2, 3):
        # Legacy (name, ID, gene) tuples retain only the two identifying fields.
        hmdb_id, gene = key[-2:]
    else:
        raise ValueError(
            "event_keys requires (hmdb_id, sensor_gene) tuples or dictionaries; "
            "display labels cannot be matched"
        )
    return _require_hmdb_id(hmdb_id), _sensor_gene_id(gene)


def _validate_optional_positive_int(value: Any, *, name: str) -> Optional[int]:
    """Validate an optional positive integer plotting limit."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be a positive integer or None")
    if value <= 0:
        raise ValueError(f"{name} must be greater than zero")
    return int(value)


def _plot_qc_series(values: pd.Series) -> pd.Series:
    """Normalize explicit QC values without confusing text True with failure."""
    import numpy as np
    from numbers import Real

    if pd.api.types.is_bool_dtype(values.dtype):
        return values.copy()
    normalized, invalid = [], []
    for position, value in enumerate(values):
        if not pd.api.types.is_scalar(value) or isinstance(value, (complex, np.complexfloating)):
            invalid.append(position)
            normalized.append(pd.NA)
        elif pd.isna(value):
            normalized.append(pd.NA)
        elif isinstance(value, (bool, np.bool_)):
            normalized.append(bool(value))
        elif isinstance(value, str) and value.strip().lower() in {'true', 'false', '1', '0'}:
            normalized.append(value.strip().lower() in {'true', '1'})
        elif isinstance(value, Real) and value in (0, 1):
            normalized.append(bool(value))
        else:
            invalid.append(position)
            normalized.append(pd.NA)
    if invalid:
        records = ', '.join(f'position {i} (index {values.index[i]!r}): {values.iloc[i]!r}'
                            for i in invalid[:5])
        raise ValueError(f'{_MIN_CELLS_QC_COLUMN} must contain booleans, True/False text, '
                         f'0/1 or genuine missing values; invalid records: {records}')
    return pd.Series(normalized, index=values.index, name=values.name, dtype='boolean')


def _apply_min_cells_event_qc(
    events: pd.DataFrame,
    *,
    qc_only: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return the requested event view and transparent min-cells QC metadata."""
    if not isinstance(qc_only, bool):
        raise TypeError("qc_only must be a boolean")

    has_qc_column = _MIN_CELLS_QC_COLUMN in events.columns
    if has_qc_column:
        # Validate even with qc_only=False: QC summaries must not silently label
        # malformed/text True values as failures. Keep NA and caller data intact.
        # Regression: tests/test_plotting_qc_types.py.
        events = events.copy()
        events[_MIN_CELLS_QC_COLUMN] = _plot_qc_series(events[_MIN_CELLS_QC_COLUMN])
        passes_qc = events[_MIN_CELLS_QC_COLUMN].eq(True).fillna(False).astype(bool)
        failed_events = events.loc[~passes_qc].copy()
        passing_count = int(passes_qc.sum())
    else:
        passes_qc = pd.Series(True, index=events.index, dtype=bool)
        failed_events = events.iloc[0:0].copy()
        passing_count = len(events)

    qc_applied = qc_only and has_qc_column
    if qc_applied:
        visible_events = events.loc[passes_qc].copy()
        excluded_events = failed_events.copy()
    else:
        visible_events = events.copy()
        excluded_events = events.iloc[0:0].copy()

    qc_summary = {
        "qc_only": qc_only,
        "qc_column": _MIN_CELLS_QC_COLUMN if has_qc_column else None,
        "qc_column_present": has_qc_column,
        "qc_applied": qc_applied,
        "events_before_qc": len(events),
        "events_passing_qc": passing_count,
        "events_failing_qc": len(failed_events),
        "events_after_qc": len(visible_events),
        "events_excluded_by_qc": len(excluded_events),
    }
    return visible_events, failed_events, excluded_events, qc_summary


def _events_frame(result_or_events: Any) -> pd.DataFrame:
    if isinstance(result_or_events, pd.DataFrame):
        return result_or_events
    events = getattr(result_or_events, "events", None)
    if isinstance(events, pd.DataFrame):
        return events
    raise TypeError("Expected a CellMeshResult-like object with .events or a pandas DataFrame")


def _sample_events_frame(result_or_sample_events: Any) -> pd.DataFrame:
    if isinstance(result_or_sample_events, pd.DataFrame):
        return result_or_sample_events
    sample_events = getattr(result_or_sample_events, "sample_events", None)
    if isinstance(sample_events, pd.DataFrame):
        return sample_events
    availability = getattr(result_or_sample_events, "availability_results", None)
    if isinstance(availability, dict):
        sample_events = availability.get("sample_events")
        if isinstance(sample_events, pd.DataFrame):
            return sample_events
    raise TypeError(
        "Expected sample-level events as a pandas DataFrame or a "
        "CellMeshResult-like object with .sample_events"
    )


def plot_significant_event_counts(
    result_or_events: Any,
    *,
    min_cell_mesh_score: Optional[float] = None,
    max_perm_pvalue: Optional[float] = None,
    max_fdr: Optional[float] = 0.05,
    fdr_col: str = "fdr_sensor_type",
    qc_only: bool = True,
    unique_keys: Optional[list[str]] = None,
    cmap: str = "Blues",
    ax: Any = None,
    show_values: bool = True,
) -> dict[str, Any]:
    """
    Plot thresholded sender-to-receiver metabolite-sensor event counts.

    Counts are based on unique ``sender + receiver + hmdb_id + sensor_gene``
    combinations after optional score, permutation p-value, and
    FDR filters. By default, the plot summarizes events with
    ``fdr_sensor_type <= 0.05``.
    Set ``fdr_col`` to ``"fdr_global"`` or ``"fdr_sensor_type"`` in either
    inference mode to select the desired correction scope.
    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    excludes rows that do not explicitly pass that QC. Set ``qc_only=False``
    to display all rows. Results without the QC column retain their historical
    behavior.
    ``metabolite`` is optional display metadata. ``unique_keys`` must include
    the four identity/context fields above; only ``sample`` may be added to
    count sample-specific records separately.

    Numeric strings are converted before filtering, ranking and deduplication.
    Native scores and probabilities must be finite values in [0, 1] or genuine
    missing values; malformed text, infinity, complex values and booleans raise
    an error identifying the column and records. Missing scores are excluded
    and returned in ``missing_score_events`` with reason ``missing_score``.
    """
    min_cell_mesh_score, max_perm_pvalue, max_fdr = _plot_thresholds(
        min_cell_mesh_score, max_perm_pvalue, max_fdr, "cell_mesh_score",
    )
    source_events = _canonical_event_table(_events_frame(result_or_events))
    if source_events.empty:
        raise ValueError("No events available to plot")

    required = {"sender", "receiver", "hmdb_id", "sensor_gene", "cell_mesh_score"}
    if max_perm_pvalue is not None:
        required.add("perm_pvalue")
    if max_fdr is not None:
        required.add(fdr_col)
    missing = sorted(required.difference(source_events.columns))
    if missing:
        raise KeyError(f"Missing required event columns: {', '.join(missing)}")

    events, qc_failed_events, qc_excluded_events, qc_summary = _apply_min_cells_event_qc(
        source_events,
        qc_only=qc_only,
    )
    probability_cols = (["perm_pvalue"] if max_perm_pvalue is not None else [])
    if max_fdr is not None:
        probability_cols.append(fdr_col)
    filtered, missing_score_events = _plot_event_numbers(events, "cell_mesh_score", probability_cols)
    if min_cell_mesh_score is not None:
        filtered = filtered[filtered["cell_mesh_score"] >= min_cell_mesh_score]
    if max_perm_pvalue is not None:
        filtered = filtered[filtered["perm_pvalue"].notna()]
        filtered = filtered[filtered["perm_pvalue"] <= max_perm_pvalue]
    if max_fdr is not None:
        filtered = filtered[filtered[fdr_col].notna()]
        filtered = filtered[filtered[fdr_col] <= max_fdr]

    all_senders = pd.Index(sorted(source_events["sender"].unique()), name="sender")
    all_receivers = pd.Index(sorted(source_events["receiver"].unique()), name="receiver")

    unique_keys = _event_record_keys(unique_keys)
    missing_unique = sorted(set(unique_keys).difference(source_events.columns))
    if missing_unique:
        raise KeyError(f"Missing unique key columns: {', '.join(missing_unique)}")

    unique_pair_events = filtered.drop_duplicates(unique_keys).copy()
    if unique_pair_events.empty:
        counts = pd.DataFrame(0, index=all_senders, columns=all_receivers)
    else:
        counts = (
            unique_pair_events
            .groupby(["sender", "receiver"])
            .size()
            .unstack(fill_value=0)
            .reindex(index=all_senders, columns=all_receivers, fill_value=0)
        )

    sender_order = counts.sum(axis=1).sort_values(ascending=False).index
    receiver_order = counts.sum(axis=0).sort_values(ascending=False).index
    counts = counts.loc[sender_order, receiver_order]

    title_parts = []
    if min_cell_mesh_score is not None:
        title_parts.append(f"score >= {min_cell_mesh_score:g}")
    if max_perm_pvalue is not None:
        title_parts.append(f"p <= {max_perm_pvalue:g}")
    if max_fdr is not None:
        title_parts.append(f"FDR <= {max_fdr:g}")
    threshold_label = ", ".join(title_parts) if title_parts else "no thresholds"

    import numpy as np
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(
            figsize=(1.35 * len(counts.columns) + 3.0, 1.0 * len(counts.index) + 2.8)
        )
    fig = ax.figure
    image = ax.imshow(counts.to_numpy(dtype=float), cmap=cmap, aspect="auto")

    ax.set_xticks(np.arange(len(counts.columns)))
    ax.set_xticklabels(counts.columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(counts.index)))
    ax.set_yticklabels(counts.index)
    ax.set_xlabel("Receiver")
    ax.set_ylabel("Sender")
    ax.set_title(f"Significant metabolite-sensor event counts ({threshold_label})")

    if show_values:
        max_count = counts.to_numpy().max()
        text_threshold = max_count / 2 if max_count > 0 else 0
        for i, sender in enumerate(counts.index):
            for j, receiver in enumerate(counts.columns):
                count = int(counts.loc[sender, receiver])
                text_color = "white" if count > text_threshold else "black"
                ax.text(j, i, str(count), ha="center", va="center", color=text_color, fontsize=10)

    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Unique metabolite-sensor event count")
    fig.tight_layout()

    summary = (
        counts
        .stack()
        .rename("metabolite_sensor_count")
        .reset_index()
        .query("metabolite_sensor_count > 0")
        .sort_values("metabolite_sensor_count", ascending=False)
    )

    return {
        "fig": fig,
        "ax": ax,
        "counts": counts,
        "summary": summary,
        "filtered_events": filtered.copy(),
        "qc_failed_events": qc_failed_events,
        "qc_excluded_events": qc_excluded_events,
        "missing_score_events": missing_score_events,
        "qc": qc_summary,
        "thresholds": {
            "min_cell_mesh_score": min_cell_mesh_score,
            "max_perm_pvalue": max_perm_pvalue,
            "max_fdr": max_fdr,
            "fdr_col": fdr_col,
            "events_before_filter": len(source_events),
            "events_after_filter": len(filtered),
            "unique_events_after_filter": len(unique_pair_events),
            **qc_summary,
        },
    }


def plot_communication_network(
    result_or_events: Any,
    *,
    min_cell_mesh_score: Optional[float] = None,
    max_perm_pvalue: Optional[float] = None,
    max_fdr: Optional[float] = 0.05,
    fdr_col: str = "fdr_sensor_type",
    score_col: str = "cell_mesh_score",
    qc_only: bool = True,
    sender_labels: Optional[list[str]] = None,
    receiver_labels: Optional[list[str]] = None,
    unique_keys: Optional[list[str]] = None,
    node_size_by: Optional[str] = "connected_cell_type_count",
    edge_color_by: Optional[str] = "total_score",
    node_cmap: str = "tab20",
    edge_cmap: str = "viridis",
    fixed_edge_color: str = "#6b7280",
    node_size_range: tuple[float, float] = (500.0, 2200.0),
    edge_width_range: tuple[float, float] = (0.8, 6.0),
    node_value_range: Optional[tuple[float, float]] = None,
    edge_width_value_range: Optional[tuple[float, float]] = None,
    edge_color_value_range: Optional[tuple[float, float]] = None,
    edge_alpha: float = 0.72,
    curve: float = 0.18,
    node_label_size: float = 9.0,
    title: Optional[str] = None,
    figsize: Optional[tuple[float, float]] = None,
    show_legends: bool = True,
    show_colorbar: bool = True,
    ax: Any = None,
) -> dict[str, Any]:
    """
    Plot a directed circular network of cell-cell communication events.

    Events are deduplicated by sender, receiver, normalized HMDB ID and sensor
    gene, and then aggregated for each
    directed sender-receiver pair. Edge width represents the number of unique
    communication events. By default, edge color represents the sum of
    ``score_col`` and node area represents the number of distinct cell types
    connected to that node in either direction.
    Names are optional display metadata and never enter deduplication.
    ``unique_keys`` may add ``sample`` to the four required fields, but cannot
    replace them or add name/annotation fields.

    ``node_size_by`` accepts ``"connected_cell_type_count"``,
    ``"connection_count"``, ``"event_count"``, ``"total_score"``, or
    ``None`` for fixed-size nodes. ``edge_color_by`` accepts
    ``"total_score"``, ``"mean_score"``, ``"event_count"``, or ``None``
    for ``fixed_edge_color``. The returned node and edge tables contain the
    exact aggregate and visual values used in the plot.

    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    uses only rows that explicitly pass min-cells QC. Set ``qc_only=False`` to
    include all rows. Older event tables without this column are unchanged.

    Numeric strings are converted before filtering, ranking and deduplication.
    Native scores and probabilities must be finite values in [0, 1] or genuine
    missing values; malformed text, infinity, complex values and booleans raise
    an error identifying the column and records. Missing scores are excluded
    and returned in ``missing_score_events`` with reason ``missing_score``.
    Custom ``score_col`` values must be finite and non-negative (no upper bound).
    Cell-type selectors are stripped like source labels; normalized duplicates
    are rejected. Only the QC-retained, explicitly selected context is checked.
    """
    min_cell_mesh_score, max_perm_pvalue, max_fdr = _plot_thresholds(
        min_cell_mesh_score, max_perm_pvalue, max_fdr, score_col,
    )
    source_events = _canonical_event_table(_events_frame(result_or_events))
    if source_events.empty:
        raise ValueError("No events available to plot")

    required = {
        "sender",
        "receiver",
        "hmdb_id",
        "sensor_gene",
        score_col,
    }
    if max_perm_pvalue is not None:
        required.add("perm_pvalue")
    if max_fdr is not None:
        required.add(fdr_col)
    missing = sorted(required.difference(source_events.columns))
    if missing:
        raise KeyError(f"Missing required event columns: {', '.join(missing)}")

    events, qc_failed_events, qc_excluded_events, qc_summary = _apply_min_cells_event_qc(
        source_events,
        qc_only=qc_only,
    )

    valid_node_metrics = {
        "connected_cell_type_count",
        "connection_count",
        "event_count",
        "total_score",
        None,
    }
    valid_edge_color_metrics = {"total_score", "mean_score", "event_count", None}
    if node_size_by not in valid_node_metrics:
        raise ValueError(
            "node_size_by must be one of connected_cell_type_count, "
            "connection_count, event_count, total_score, or None"
        )
    if edge_color_by not in valid_edge_color_metrics:
        raise ValueError(
            "edge_color_by must be one of total_score, mean_score, event_count, or None"
        )

    node_size_range = _plot_range(node_size_range, "node_size_range", positive=True)
    edge_width_range = _plot_range(edge_width_range, "edge_width_range", positive=True)
    if node_value_range is not None:
        node_value_range = _plot_range(node_value_range, "node_value_range")
    if edge_width_value_range is not None:
        edge_width_value_range = _plot_range(edge_width_value_range, "edge_width_value_range")
    if edge_color_value_range is not None:
        edge_color_value_range = _plot_range(edge_color_value_range, "edge_color_value_range")
    edge_alpha = _plot_number(edge_alpha, "edge_alpha", upper=1.0)
    curve = _plot_number(curve, "curve")
    node_label_size = _plot_number(node_label_size, "node_label_size", positive=True)
    if figsize is not None:
        if isinstance(figsize, (str, bytes)) or len(figsize) != 2:
            raise ValueError("figsize must contain two positive finite values")
        figsize = tuple(_plot_number(v, "figsize", positive=True) for v in figsize)

    filtered = events.copy()
    if sender_labels is not None:
        filtered = filtered[filtered["sender"].isin(_plot_labels(sender_labels, "sender_labels"))]
    if receiver_labels is not None:
        filtered = filtered[filtered["receiver"].isin(_plot_labels(receiver_labels, "receiver_labels"))]
    filtered, missing_score_events = _plot_event_numbers(filtered, score_col, ["perm_pvalue", fdr_col])
    if min_cell_mesh_score is not None:
        filtered = filtered[filtered[score_col] >= min_cell_mesh_score]
    if max_perm_pvalue is not None:
        filtered = filtered[filtered["perm_pvalue"].notna()]
        filtered = filtered[filtered["perm_pvalue"] <= max_perm_pvalue]
    if max_fdr is not None:
        filtered = filtered[filtered[fdr_col].notna()]
        filtered = filtered[filtered[fdr_col] <= max_fdr]
    if filtered.empty:
        raise ValueError("No events remain after filtering")

    unique_keys = _event_record_keys(unique_keys)
    missing_unique = sorted(set(unique_keys).difference(source_events.columns))
    if missing_unique:
        raise KeyError(f"Missing unique key columns: {', '.join(missing_unique)}")

    sort_columns: list[str] = []
    sort_ascending: list[bool] = []
    if fdr_col in filtered.columns:
        sort_columns.append(fdr_col)
        sort_ascending.append(True)
    if "perm_pvalue" in filtered.columns:
        sort_columns.append("perm_pvalue")
        sort_ascending.append(True)
    sort_columns.append(score_col)
    sort_ascending.append(False)
    unique_events = (
        filtered
        .sort_values(sort_columns, ascending=sort_ascending, na_position="last")
        .drop_duplicates(unique_keys)
        .copy()
    )

    edge_table = (
        unique_events
        .groupby(["sender", "receiver"], as_index=False, sort=True)
        .agg(
            event_count=(score_col, "size"),
            total_score=(score_col, "sum"),
            mean_score=(score_col, "mean"),
        )
    )
    node_names = sorted(
        set(edge_table["sender"].astype(str)).union(edge_table["receiver"].astype(str))
    )
    if not node_names:
        raise ValueError("No sender or receiver cell types remain after filtering")

    node_rows = []
    for node in node_names:
        incoming_edges = edge_table[edge_table["receiver"].astype(str) == node]
        outgoing_edges = edge_table[edge_table["sender"].astype(str) == node]
        incident_edges = edge_table[
            (edge_table["sender"].astype(str) == node)
            | (edge_table["receiver"].astype(str) == node)
        ]
        incident_events = unique_events[
            (unique_events["sender"].astype(str) == node)
            | (unique_events["receiver"].astype(str) == node)
        ]
        connected_nodes = set(outgoing_edges["receiver"].astype(str)).union(
            incoming_edges["sender"].astype(str)
        )
        node_rows.append(
            {
                "cell_type": node,
                "connected_cell_type_count": len(connected_nodes),
                "connection_count": len(incident_edges),
                "event_count": len(incident_events),
                "incoming_event_count": int(incoming_edges["event_count"].sum()),
                "outgoing_event_count": int(outgoing_edges["event_count"].sum()),
                "total_score": float(incident_events[score_col].sum()),
            }
        )
    node_table = pd.DataFrame(node_rows)

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as path_effects
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch

    def _scale_values(
        values: Any,
        output_range: tuple[float, float],
        value_range: Optional[tuple[float, float]],
    ) -> tuple[Any, tuple[float, float]]:
        array = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError("Visual encoding values must be finite")
        if value_range is None:
            value_min = float(array.min())
            value_max = float(array.max())
        else:
            if len(value_range) != 2:
                raise ValueError("Visual value ranges must contain exactly two values")
            value_min, value_max = map(float, value_range)
            if not np.isfinite(value_min) or not np.isfinite(value_max):
                raise ValueError("Visual value ranges must be finite")
            if value_max < value_min:
                raise ValueError("Visual value range maximum must be >= minimum")
        output_min, output_max = output_range
        if value_max == value_min:
            scaled = np.full(array.shape, (output_min / 2.0 + output_max / 2.0))
        else:
            clipped = np.clip(array, value_min, value_max)
            fraction = (clipped - value_min) / (value_max - value_min)
            scaled = output_min + fraction * (output_max - output_min)
        return scaled, (value_min, value_max)

    if node_size_by is None:
        node_table["node_size_value"] = np.nan
        node_table["node_size"] = (node_size_range[0] / 2.0 + node_size_range[1] / 2.0)
        node_limits = None
    else:
        node_table["node_size_value"] = node_table[node_size_by].astype(float)
        node_sizes, node_limits = _scale_values(
            node_table["node_size_value"], node_size_range, node_value_range
        )
        node_table["node_size"] = node_sizes

    edge_widths, edge_width_limits = _scale_values(
        edge_table["event_count"], edge_width_range, edge_width_value_range
    )
    edge_table["edge_width"] = edge_widths

    edge_cmap_object = plt.get_cmap(edge_cmap)
    colorbar = None
    color_norm = None
    edge_color_limits = None
    if edge_color_by is None:
        edge_table["edge_color_value"] = np.nan
        edge_table["edge_color"] = fixed_edge_color
    else:
        color_values = edge_table[edge_color_by].to_numpy(dtype=float)
        if edge_color_value_range is None:
            color_min = float(color_values.min())
            color_max = float(color_values.max())
        else:
            if len(edge_color_value_range) != 2:
                raise ValueError("edge_color_value_range must contain exactly two values")
            color_min, color_max = map(float, edge_color_value_range)
            if not np.isfinite(color_min) or not np.isfinite(color_max):
                raise ValueError("edge_color_value_range must be finite")
            if color_max < color_min:
                raise ValueError("edge_color_value_range maximum must be >= minimum")
        edge_color_limits = (color_min, color_max)
        if color_max == color_min:
            color_delta = max(abs(color_min) * 0.05, 0.5)
            color_norm = Normalize(color_min - color_delta, color_max + color_delta)
        else:
            color_norm = Normalize(color_min, color_max, clip=True)
        edge_table["edge_color_value"] = color_values
        edge_table["edge_color"] = [
            tuple(edge_cmap_object(color_norm(value))) for value in color_values
        ]

    node_cmap_object = plt.get_cmap(node_cmap)
    denominator = max(len(node_table) - 1, 1)
    node_table["node_color"] = [
        tuple(node_cmap_object(index / denominator)) for index in range(len(node_table))
    ]
    angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, len(node_table), endpoint=False)
    node_table["position_x"] = np.cos(angles)
    node_table["position_y"] = np.sin(angles)

    created_figure = ax is None
    if ax is None:
        if figsize is None:
            figsize = (
                max(7.5, 0.55 * len(node_table) + 5.5),
                max(6.5, 0.35 * len(node_table) + 4.5),
            )
        _, ax = plt.subplots(figsize=figsize)
    fig = ax.figure
    positions = {
        row.cell_type: np.array([row.position_x, row.position_y], dtype=float)
        for row in node_table.itertuples()
    }
    node_sizes_by_name = dict(zip(node_table["cell_type"], node_table["node_size"]))

    edge_artists = []
    for row in edge_table.itertuples():
        sender = str(row.sender)
        receiver = str(row.receiver)
        start = positions[sender].copy()
        end = positions[receiver].copy()
        connection_style = f"arc3,rad={curve:g}"
        sender_shrink = max(9.0, np.sqrt(node_sizes_by_name[sender]) / 2.2)
        receiver_shrink = max(9.0, np.sqrt(node_sizes_by_name[receiver]) / 2.2)
        if sender == receiver:
            radial = start / max(float(np.linalg.norm(start)), 1.0)
            tangent = np.array([-radial[1], radial[0]])
            loop_radius = 0.08 + np.sqrt(node_sizes_by_name[sender]) / 450.0
            loop_center = start + radial * loop_radius * 0.62
            start = loop_center - tangent * loop_radius * 0.92
            end = loop_center + tangent * loop_radius * 0.92
            connection_style = "arc3,rad=1.45"
            sender_shrink = 0.0
            receiver_shrink = 0.0
        edge_color = row.edge_color
        artist = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            connectionstyle=connection_style,
            mutation_scale=max(9.0, 7.0 + 1.8 * row.edge_width),
            linewidth=row.edge_width,
            color=edge_color,
            alpha=edge_alpha,
            shrinkA=sender_shrink,
            shrinkB=receiver_shrink,
            capstyle="round",
            joinstyle="round",
            zorder=1,
        )
        ax.add_patch(artist)
        edge_artists.append(artist)

    node_artist = ax.scatter(
        node_table["position_x"],
        node_table["position_y"],
        s=node_table["node_size"],
        c=node_table["node_color"].tolist(),
        edgecolors="#1f2937",
        linewidths=1.0,
        zorder=3,
    )
    for row in node_table.itertuples():
        label_x = row.position_x * 1.20
        label_y = row.position_y * 1.20
        horizontal_alignment = (
            "left" if row.position_x > 0.12 else "right" if row.position_x < -0.12 else "center"
        )
        vertical_alignment = (
            "bottom" if row.position_y > 0.12 else "top" if row.position_y < -0.12 else "center"
        )
        text_artist = ax.text(
            label_x,
            label_y,
            row.cell_type,
            ha=horizontal_alignment,
            va=vertical_alignment,
            fontsize=node_label_size,
            color="#111827",
            zorder=4,
        )
        text_artist.set_path_effects(
            [path_effects.withStroke(linewidth=3.0, foreground="white")]
        )

    title_parts = []
    if min_cell_mesh_score is not None:
        title_parts.append(f"score >= {min_cell_mesh_score:g}")
    if max_perm_pvalue is not None:
        title_parts.append(f"p <= {max_perm_pvalue:g}")
    if max_fdr is not None:
        title_parts.append(f"{fdr_col} <= {max_fdr:g}")
    threshold_label = ", ".join(title_parts) if title_parts else "no thresholds"
    if title is None:
        title = f"Cell-cell communication network ({threshold_label})"
    ax.set_title(title, pad=16)
    ax.set_aspect("equal")
    # Keep outward-facing labels inside the axes even when a compact figure is
    # rendered with tight bounding boxes (as in Jupyter's inline backend).
    ax.set_xlim(-1.80, 1.80)
    ax.set_ylim(-1.43, 1.43)
    ax.axis("off")

    score_label = "Cellmesh_score" if score_col == "cell_mesh_score" else score_col
    node_metric_labels = {
        "connected_cell_type_count": "Connected cell types",
        "connection_count": "Directed cell-pair links",
        "event_count": "Incident communication events",
        "total_score": f"Incident total {score_label}",
    }
    edge_metric_labels = {
        "event_count": "Communication events",
        "total_score": f"Cell-pair total {score_label}",
        "mean_score": f"Cell-pair mean {score_label}",
    }
    legends: dict[str, Any] = {}

    def _representative_values(values: Any, integer: bool = False) -> list[float]:
        array = np.asarray(values, dtype=float)
        representatives = np.quantile(array, [0.0, 0.5, 1.0])
        if integer:
            representatives = np.rint(representatives)
        unique_representatives = []
        for value in representatives:
            numeric_value = float(value)
            if numeric_value not in unique_representatives:
                unique_representatives.append(numeric_value)
        return unique_representatives

    if show_legends and node_size_by is not None:
        node_legend_values = _representative_values(
            node_table["node_size_value"],
            integer=node_size_by != "total_score",
        )
        node_legend_sizes, _ = _scale_values(
            node_legend_values, node_size_range, node_value_range or node_limits
        )
        node_legend_sizes = np.clip(node_legend_sizes * 0.22, 35.0, 420.0)
        node_handles = [
            ax.scatter(
                [],
                [],
                s=size,
                color="#9ca3af",
                edgecolor="#1f2937",
                linewidth=0.8,
            )
            for size in node_legend_sizes
        ]
        node_legend = ax.legend(
            node_handles,
            [f"{value:g}" for value in node_legend_values],
            title=node_metric_labels[node_size_by],
            loc="lower left",
            frameon=False,
            labelspacing=1.2,
            borderaxespad=0.0,
        )
        ax.add_artist(node_legend)
        legends["node_size"] = node_legend

    if show_legends:
        edge_legend_values = _representative_values(edge_table["event_count"], integer=True)
        edge_legend_widths, _ = _scale_values(
            edge_legend_values,
            edge_width_range,
            edge_width_value_range or edge_width_limits,
        )
        edge_handles = [
            Line2D([0], [0], color="#4b5563", linewidth=width, alpha=edge_alpha)
            for width in edge_legend_widths
        ]
        edge_legend = ax.legend(
            edge_handles,
            [f"{value:g}" for value in edge_legend_values],
            title="Communication events",
            loc="lower right",
            frameon=False,
            borderaxespad=0.0,
        )
        legends["edge_width"] = edge_legend

    if show_colorbar and edge_color_by is not None:
        scalar_mappable = plt.cm.ScalarMappable(norm=color_norm, cmap=edge_cmap_object)
        scalar_mappable.set_array([])
        colorbar = fig.colorbar(
            scalar_mappable,
            ax=ax,
            orientation="horizontal",
            fraction=0.045,
            pad=0.025,
            shrink=0.58,
            aspect=30,
        )
        colorbar.set_label(edge_metric_labels[edge_color_by])

    if created_figure:
        fig.tight_layout()

    return {
        "fig": fig,
        "ax": ax,
        "node_table": node_table,
        "edge_table": edge_table,
        "filtered_events": filtered.copy(),
        "unique_events": unique_events,
        "qc_failed_events": qc_failed_events,
        "qc_excluded_events": qc_excluded_events,
        "missing_score_events": missing_score_events,
        "qc": qc_summary,
        "node_artist": node_artist,
        "edge_artists": edge_artists,
        "colorbar": colorbar,
        "legends": legends,
        "thresholds": {
            "min_cell_mesh_score": min_cell_mesh_score,
            "max_perm_pvalue": max_perm_pvalue,
            "max_fdr": max_fdr,
            "fdr_col": fdr_col,
            "events_before_filter": len(source_events),
            "events_after_filter": len(filtered),
            "unique_events_after_filter": len(unique_events),
            **qc_summary,
        },
        "encodings": {
            "node_size_by": node_size_by,
            "node_size_range": node_size_range,
            "node_value_range": node_limits,
            "edge_width_by": "event_count",
            "edge_width_range": edge_width_range,
            "edge_width_value_range": edge_width_limits,
            "edge_color_by": edge_color_by,
            "edge_color_value_range": edge_color_limits,
        },
    }



def _layout_dotplot(fig, axes, legend_ax, colorbar_ax, legend, title, *, region,
                    n_events: int, n_receivers: int) -> dict[str, Any]:
    """Allocate measured decorations and two separate sidebar boxes in inches."""
    import numpy as np
    from matplotlib.transforms import Bbox

    # 防回退：布局不能依赖 has_missing。零值说明、固定大小说明、长概率标签
    # 都必须计入实际 renderer 尺寸；图例和色条必须拥有独立区域。
    # 不调用整张图的 tight_layout/subplots_adjust，不重排调用者的其他子图。
    # 回归测试：tests/test_dotplot_layout.py。
    gap, edge, bar_width = 0.18, 0.10, 0.18
    margins = np.zeros((len(axes), 4))  # left, bottom, right, top
    bar_margins = np.zeros(4)
    required = None
    for attempt in range(6):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        dpi = fig.dpi

        def extras(axis):
            body = axis.get_window_extent(renderer)
            outer = axis.get_tightbbox(renderer)
            return np.maximum(0, np.array([body.x0 - outer.x0, body.y0 - outer.y0,
                                           outer.x1 - body.x1, outer.y1 - body.y1]) / dpi)

        margins = np.maximum(margins, np.array([extras(axis) for axis in axes]))
        bar_margins = np.maximum(bar_margins, extras(colorbar_ax))
        legend_box = legend.get_window_extent(renderer)
        legend_width, legend_height = legend_box.width / dpi, legend_box.height / dpi
        title_box = title.get_window_extent(renderer)
        title_width, title_height = title_box.width / dpi, title_box.height / dpi
        sidebar_width = max(legend_width, bar_width + bar_margins[0] + bar_margins[2])
        sidebar_height = legend_height + gap + 1.0 + bar_margins[1] + bar_margins[3]
        plot_min_width = max(0.75, 0.20 * n_receivers)
        plot_height = max(1.25, 0.30 * n_events) if region is None else 0.75
        min_width = (sum(margins[:, 0] + margins[:, 2]) + len(axes) * plot_min_width
                     + gap * len(axes) + sidebar_width + 2 * edge)
        min_width = max(min_width, title_width + 2 * edge)
        min_height = max(plot_height + max(margins[:, 1] + margins[:, 3]), sidebar_height)
        min_height += title_height + gap + 2 * edge
        required = (float(min_width), float(min_height))
        if region is None:
            target_width = min_width + len(axes) * (max(1.7, 0.48 * n_receivers) - plot_min_width)
            width, height = max(9.0, target_width), max(4.6, min_height)
            fig.set_size_inches(width, height, forward=True)
            x0, y0 = 0.0, 0.0
        else:
            fig_width, fig_height = fig.get_size_inches()
            x0, y0 = region.x0 * fig_width, region.y0 * fig_height
            width, height = region.width * fig_width, region.height * fig_height
            if width + 1e-8 < min_width or height + 1e-8 < min_height:
                raise ValueError(
                    f"Dotplot layout needs at least {min_width:.2f} x {min_height:.2f} inches "
                    f"inside the supplied ax; available {width:.2f} x {height:.2f}. "
                    "Enlarge the figure or allocate a larger subplot before plotting."
                )
        fig_width, fig_height = fig.get_size_inches()
        content_bottom = y0 + edge
        content_top = y0 + height - edge - title_height - gap
        content_height = content_top - content_bottom
        panel_width = (width - 2 * edge - sidebar_width - gap * len(axes)
                       - sum(margins[:, 0] + margins[:, 2])) / len(axes)
        cursor = x0 + edge
        for axis, (left, bottom, right, top) in zip(axes, margins):
            axis.set_position([(cursor + left) / fig_width, (content_bottom + bottom) / fig_height,
                               panel_width / fig_width, (content_height - bottom - top) / fig_height])
            cursor += left + panel_width + right + gap
        legend_ax.set_position([cursor / fig_width, (content_top - legend_height) / fig_height,
                                sidebar_width / fig_width, legend_height / fig_height])
        colorbar_ax.set_position([(cursor + bar_margins[0]) / fig_width,
                                 (content_bottom + bar_margins[1]) / fig_height,
                                 bar_width / fig_width,
                                 (content_height - legend_height - gap - bar_margins[1]
                                  - bar_margins[3]) / fig_height])
        title.set_position(((x0 + width / 2) / fig_width, (y0 + height - edge) / fig_height))
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        owned_box = Bbox.from_bounds(x0 * dpi, y0 * dpi, width * dpi, height * dpi)
        boxes = [axis.get_tightbbox(renderer) for axis in axes]
        boxes += [legend.get_window_extent(renderer), colorbar_ax.get_tightbbox(renderer),
                  title.get_window_extent(renderer)]
        contained = all(box.x0 >= owned_box.x0 - 0.5 and box.y0 >= owned_box.y0 - 0.5
                        and box.x1 <= owned_box.x1 + 0.5 and box.y1 <= owned_box.y1 + 0.5 for box in boxes)
        separate = all(not a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1:])
        if contained and separate:
            return {"required_size_inches": required, "available_size_inches": (width, height),
                    "external_ax": region is not None}
    raise ValueError("Dotplot layout could not fit all labels; enlarge the figure or subplot.")

def plot_event_dotplot(
    result_or_events: Any,
    *,
    top_n: Optional[int] = 10,
    event_keys: Optional[list[Any]] = None,
    min_cell_mesh_score: Optional[float] = None,
    max_perm_pvalue: Optional[float] = None,
    max_fdr: Optional[float] = None,
    sender_labels: Optional[list[str]] = None,
    receiver_labels: Optional[list[str]] = None,
    score_col: str = "cell_mesh_score",
    score_label: Optional[str] = None,
    fdr_col: str = "fdr_sensor_type",
    pvalue_col: str = "perm_pvalue",
    qc_only: bool = True,
    cmap: str = "plasma",
    ax: Any = None,
    min_dot_size: float = 20.0,
    max_dot_size: float = 260.0,
) -> dict[str, Any]:
    """
    Bubble plot for metabolite-sensor events across sender-receiver pairs.

    Rows are metabolite-to-sensor events, each subplot is one sender cell type,
    columns are receiver cell types, bubble color is ``score_col`` and bubble
    size follows ``-log10(fdr_col)`` for positive probabilities. Zero uses a
    separately labeled maximum area; when zeros are present, positive values
    use the lower 80% of the area interval. This reservation is display-only.
    Without zeros, the existing positive size mapping is retained. If the FDR
    column is absent, the p-value column supplies the size scale. Missing significance uses a fixed-size diamond;
    its color still represents the score. When no significance values are
    available, all points have fixed size and the legend reads "Significance:
    Unavailable" without a numeric scale. Both statistic branches label actual
    probabilities, without inventing a positive floor for zero. Equal requested
    size bounds explicitly disable size encoding. Returned ``dot_sizes`` align
    with ``plot_events``; ``size_encoding`` records the source column, display
    ranges and positive legend probabilities. An existing but missing FDR never
    falls back to an unadjusted p-value. Non-missing probabilities must be
    finite and in [0, 1]; input values, including NA, are preserved.
    If ``event_keys`` is not provided, the top
    ``top_n`` events are selected by lowest FDR and highest score across all
    sender-receiver pairs. Each row is identified by the full
    ``(hmdb_id, sensor_gene)`` tuple, independently of its label.

    ``event_keys`` accepts ``(hmdb_id, sensor_gene)`` tuples or dictionaries
    containing those fields. Legacy three-field tuples are accepted with their
    first (name) field ignored. Display-label strings cannot select events.
    Selectors retain their requested order, with repeated matches included once.
    Labels use the first available name per HMDB ID and always include the ID
    and sensor gene; missing names fall back to the ID. ``selected_events``
    contains these labels, while ``selected_event_keys`` contains the two-field
    tuples to use in another selection. Numeric scores and p-values are unchanged.
    When supplying an existing ``ax``, the filtered data must contain exactly
    one sender; use ``sender_labels`` to select it. Leave ``ax=None`` for the
    default multi-sender faceted layout.
    Layout measures the rendered legend, colorbar and axis labels; automatic
    figures expand to fit. Supplied axes reserve the sidebar inside their
    existing rectangle, without resizing the figure or other axes or replacing
    a caller's suptitle. Finish any automatic layout first with
    ``fig.canvas.draw(); fig.set_layout_engine('none')``. Active layout engines
    are rejected before modifying the axes. An undersized supplied rectangle
    raises an error giving the minimum required dimensions. Finalize figure
    dimensions before plotting; later global layout changes or shrinking the
    figure require recreating the plot. ``legend``, ``legend_ax``, ``colorbar``,
    ``title_artist`` and ``layout`` are returned for inspection and export.
    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    excludes rows that do not explicitly pass min-cells QC. Set
    ``qc_only=False`` to include all rows. Older event tables without this
    column retain their historical behavior.

    Numeric strings are converted before filtering, ranking and deduplication.
    Native scores and probabilities must be finite values in [0, 1] or genuine
    missing values; malformed text, infinity, complex values and booleans raise
    an error identifying the column and records. Missing scores are excluded
    and returned in ``missing_score_events`` with reason ``missing_score``.
    Custom ``score_col`` values must be finite and non-negative (no upper bound).
    Cell-type selectors are stripped like source labels; normalized duplicates
    are rejected. Only the QC-retained, explicitly selected context is checked.
    """
    min_dot_size, max_dot_size = _plot_range((min_dot_size, max_dot_size),
                                             "dot size range", positive=True)
    top_n = _validate_optional_positive_int(top_n, name="top_n")
    min_cell_mesh_score, max_perm_pvalue, max_fdr = _plot_thresholds(
        min_cell_mesh_score, max_perm_pvalue, max_fdr, score_col,
    )
    source_events = _canonical_event_table(_events_frame(result_or_events))
    if source_events.empty:
        raise ValueError("No events available to plot")

    required = {
        "sender",
        "receiver",
        "hmdb_id",
        "sensor_gene",
        score_col,
    }
    if max_perm_pvalue is not None or pvalue_col in source_events.columns:
        required.add(pvalue_col)
    if max_fdr is not None or fdr_col in source_events.columns:
        required.add(fdr_col)
    missing = sorted(required.difference(source_events.columns))
    if missing:
        raise KeyError(f"Missing required event columns: {', '.join(missing)}")

    events, qc_failed_events, qc_excluded_events, qc_summary = _apply_min_cells_event_qc(
        source_events,
        qc_only=qc_only,
    )
    filtered = events.copy()
    if sender_labels is not None:
        filtered = filtered[filtered["sender"].isin(_plot_labels(sender_labels, "sender_labels"))]
    if receiver_labels is not None:
        filtered = filtered[filtered["receiver"].isin(_plot_labels(receiver_labels, "receiver_labels"))]
    if event_keys is not None:
        requested_keys = {_event_selector(key) for key in event_keys}
        keys = pd.Series(list(filtered[["hmdb_id", "sensor_gene"]].itertuples(index=False, name=None)),
                         index=filtered.index, dtype=object)
        filtered = filtered.loc[keys.isin(requested_keys)].copy()
        if filtered.empty:
            raise ValueError("None of the requested event_keys were found after filtering")
    filtered, missing_score_events = _plot_event_numbers(filtered, score_col, [pvalue_col, fdr_col])
    if min_cell_mesh_score is not None:
        filtered = filtered[filtered[score_col] >= min_cell_mesh_score]
    if max_perm_pvalue is not None:
        filtered = filtered[filtered[pvalue_col].notna()]
        filtered = filtered[filtered[pvalue_col] <= max_perm_pvalue]
    if max_fdr is not None:
        filtered = filtered[filtered[fdr_col].notna()]
        filtered = filtered[filtered[fdr_col] <= max_fdr]
    if filtered.empty:
        raise ValueError("No events remain after filtering")
    if score_label is None:
        score_label = "Cellmesh_score" if score_col == "cell_mesh_score" else score_col

    event_columns = ["hmdb_id", "sensor_gene"]
    catalog = source_events[event_columns].drop_duplicates().copy()
    catalog["_event_id"] = list(catalog.itertuples(index=False, name=None))
    display_names = _metabolite_display_names(source_events)
    labels_by_id = {
        key: f"{_metabolite_label(display_names.get(key[0]), key[0])} -> {key[1]}"
        for key in catalog["_event_id"]
    }
    filtered["_event_id"] = list(filtered[event_columns].itertuples(index=False, name=None))
    filtered["_event_label"] = filtered["_event_id"].map(labels_by_id)
    filtered["_event_key"] = filtered["_event_id"].map(lambda key: f"{key[0]} | {key[1]}")

    if event_keys is not None:
        available_ids = set(filtered["_event_id"])
        requested_ids = [_event_selector(key) for key in event_keys]
        event_order = list(dict.fromkeys(key for key in requested_ids if key in available_ids))
        plot_events = filtered[filtered["_event_id"].isin(event_order)].copy()
        if plot_events.empty:
            raise ValueError("None of the requested event_keys were found after filtering")
    else:
        rank_cols = ["_event_id"]
        if fdr_col in filtered.columns:
            event_rank = (
                filtered
                .groupby(rank_cols, as_index=False)
                .agg(_best_fdr=(fdr_col, "min"), _best_score=(score_col, "max"))
                .sort_values(["_best_fdr", "_best_score"], ascending=[True, False], na_position="last")
            )
        else:
            event_rank = (
                filtered
                .groupby(rank_cols, as_index=False)
                .agg(_best_score=(score_col, "max"))
                .sort_values("_best_score", ascending=False)
            )
        if top_n is not None:
            event_rank = event_rank.head(top_n)
        plot_events = filtered[filtered["_event_id"].isin(event_rank["_event_id"])].copy()

    # Keep one event per row/sender/receiver if duplicate rows exist, prioritizing lower FDR and higher score.
    sort_cols = ["_event_id", "sender", "receiver"]
    ascending = [True, True, True]
    if fdr_col in plot_events.columns:
        sort_cols.append(fdr_col)
        ascending.append(True)
    sort_cols.append(score_col)
    ascending.append(False)
    plot_events = (
        plot_events
        .sort_values(sort_cols, ascending=ascending, na_position="last")
        .drop_duplicates(["_event_id", "sender", "receiver"], keep="first")
        .reset_index(drop=True)
    )

    sender_order = (
        plot_events.groupby("sender")[score_col]
        .max()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    receiver_order = (
        plot_events.groupby("receiver")[score_col]
        .max()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    if event_keys is None and fdr_col in plot_events.columns:
        event_order = (
            plot_events.groupby("_event_id")
            .agg(_best_fdr=(fdr_col, "min"), _best_score=(score_col, "max"))
            .sort_values(["_best_fdr", "_best_score"], ascending=[True, False], na_position="last")
            .index
            .tolist()
        )
    elif event_keys is None:
        event_order = (
            plot_events.groupby("_event_id")[score_col]
            .max()
            .sort_values(ascending=False)
            .index
            .tolist()
        )

    label_order = [labels_by_id[key] for key in event_order]
    plot_events["sender"] = pd.Categorical(plot_events["sender"], categories=sender_order, ordered=True)
    plot_events["receiver"] = pd.Categorical(
        plot_events["receiver"], categories=receiver_order, ordered=True
    )
    plot_events["_x"] = plot_events["receiver"].cat.codes
    plot_events["_y"] = pd.Categorical(
        plot_events["_event_id"], categories=event_order, ordered=True,
    ).codes

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    if fdr_col in plot_events.columns:
        significance_col = fdr_col
        size_label = "FDR"
    elif pvalue_col in plot_events.columns:
        significance_col = pvalue_col
        size_label = "p-value" if pvalue_col == "perm_pvalue" else f"p-value ({pvalue_col})"
    else:
        significance_col = None
        size_label = "Significance"

    probabilities = (
        pd.to_numeric(plot_events[significance_col], errors="raise").astype(float)
        if significance_col is not None
        else pd.Series(float("nan"), index=plot_events.index)
    )
    has_significance = probabilities.notna()
    has_missing = not has_significance.all()
    valid_probabilities = probabilities.loc[has_significance]
    if not (np.isfinite(valid_probabilities) & valid_probabilities.between(0, 1)).all():
        raise ValueError(
            f"{significance_col} must contain finite probabilities in [0, 1] or missing values"
        )

    # NA 表示显著性不可用，不等于 p/FDR=0 或 1。只用有效值计算大小范围，
    # 缺失记录以固定大小菱形展示，避免 NaN 点大小让已评分事件悄悄消失。
    # 全 NA 时不生成数值显著性图例，也不能用未校正 p 值填补缺失 FDR。
    # 回归测试：tests/test_dotplot_significance.py。
    # 防回退：0 是输入的零概率，不是最小正概率，也不是人为的 1e-12。
    # 只对正概率取 log10；零单独使用显示上限，NA 保留原有菱形。
    # 有零时正值使用面积区间的前 80%，上方空间仅作显示区分，不代表统计距离。
    # 回归测试：tests/test_dotplot_zero_significance.py。
    positive = valid_probabilities[valid_probabilities > 0]
    zero_mask = probabilities.eq(0)
    has_zero = bool(zero_mask.any())
    has_available = not valid_probabilities.empty
    significance = -np.log10(positive)
    sig_min = float(significance.min()) if not positive.empty else None
    sig_max = float(significance.max()) if not positive.empty else None
    fixed_size = min_dot_size == max_dot_size
    fixed_dot_size = min_dot_size / 2 + max_dot_size / 2
    positive_max_size = min_dot_size + (0.8 if has_zero else 1.0) * (max_dot_size - min_dot_size)
    positive_fixed_size = min_dot_size / 2 + positive_max_size / 2

    def positive_sizes(log_values):
        if sig_max > sig_min:
            return min_dot_size + (log_values - sig_min) / (sig_max - sig_min) * (
                positive_max_size - min_dot_size
            )
        return np.full(np.shape(log_values), positive_fixed_size)

    dot_sizes = pd.Series(fixed_dot_size, index=plot_events.index, dtype=float)
    if not positive.empty:
        dot_sizes.loc[positive.index] = positive_sizes(significance.to_numpy())
    dot_sizes.loc[zero_mask] = max_dot_size

    if ax is None:
        fig = plt.figure(figsize=(9.0, 4.6), layout="none")
        axes = [fig.add_axes([0.30, 0.22, 0.35, 0.55]) for _ in sender_order]
        owned_region = None
        title_artist = fig.suptitle("Metabolite-sensor communication events", va="top")
    else:
        if len(sender_order) != 1:
            raise ValueError(
                "A supplied ax can display exactly one sender; use "
                "sender_labels to select one sender or leave ax=None"
            )
        fig = ax.figure
        from matplotlib.layout_engine import PlaceHolderLayoutEngine

        engine = fig.get_layout_engine()
        if engine is not None and not isinstance(engine, PlaceHolderLayoutEngine):
            raise ValueError(
                "An external dotplot ax requires a finalized figure layout. "
                "First call fig.canvas.draw(), then fig.set_layout_engine('none'), "
                "then plot_event_dotplot(..., ax=ax). Automatic tight/constrained "
                "layout would otherwise reposition this plot and other subplots."
            )
        fig.canvas.draw()
        owned_region = ax.get_position().frozen()
        axes = [ax]
        # A local title must not overwrite the caller's figure-wide suptitle.
        title_artist = fig.text(0.5, 0.95, "Metabolite-sensor communication events",
                                ha="center", va="top", fontsize=12, in_layout=False)
    lax = fig.add_axes([0.75, 0.65, 0.20, 0.25])
    cax = fig.add_axes([0.78, 0.22, 0.025, 0.25])
    lax.set_in_layout(False)
    cax.set_in_layout(False)
    lax.axis("off")

    score_values = plot_events[score_col].astype(float)
    score_min = float(score_values.min())
    score_max = float(score_values.max())
    scatter = None
    for idx, sender in enumerate(sender_order):
        axis = axes[idx]
        panel = plot_events[plot_events["sender"] == sender]
        panel_has_significance = has_significance.loc[panel.index]
        for visible, marker in (
            (panel.loc[panel_has_significance], "o"),
            (panel.loc[~panel_has_significance], "D"),
        ):
            if visible.empty:
                continue
            scatter = axis.scatter(
                visible["_x"],
                visible["_y"],
                c=visible[score_col].astype(float),
                s=dot_sizes.loc[visible.index],
                marker=marker,
                cmap=cmap,
                vmin=score_min,
                vmax=score_max,
                edgecolors="black",
                linewidths=0.35,
            )

        axis.set_xticks(range(len(receiver_order)))
        axis.set_xticklabels(receiver_order, rotation=45, ha="right")
        axis.set_yticks(range(len(label_order)))
        if idx == 0:
            axis.set_yticklabels(label_order)
            axis.set_ylabel("Metabolite -> receptor gene")
        else:
            axis.set_yticklabels([])
            axis.tick_params(axis="y", length=0)
        axis.invert_yaxis()
        axis.set_xlabel("Receiver")
        strip_height = 0.095
        axis.add_patch(
            Rectangle(
                (0.0, 1.0),
                1.0,
                strip_height,
                transform=axis.transAxes,
                facecolor="white",
                edgecolor="black",
                linewidth=1.0,
                clip_on=False,
                zorder=4,
            )
        )
        axis.text(
            0.5,
            1.0 + strip_height / 2.0,
            str(sender),
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=10,
            zorder=5,
        )
        axis.grid(axis="both", color="#e6e6e6", linewidth=0.6)
        axis.set_axisbelow(True)
        axis.set_xlim(-0.5, len(receiver_order) - 0.5)
        axis.set_ylim(len(label_order) - 0.5, -0.5)

    if scatter is None:
        raise ValueError("No events available to plot after sender panel construction")

    cbar = fig.colorbar(scatter, cax=cax)
    cbar.ax.set_title(score_label, fontsize=9, pad=6)
    cbar.ax.tick_params(direction="in")

    handles = []
    labels = []
    # Use observed probabilities for legend keys. Inverting interpolated logs can
    # invent probabilities or round tiny positive values to zero again.
    legend_probabilities = []
    if not positive.empty:
        unique_positive = np.sort(positive.unique())[::-1]
        positions = np.linspace(0, len(unique_positive) - 1, min(3, len(unique_positive))).astype(int)
        legend_probabilities = unique_positive[positions].tolist()
        legend_sizes = positive_sizes(-np.log10(np.asarray(legend_probabilities)))
        handles = [
            Line2D([], [], linestyle="none", marker="o", markersize=np.sqrt(size),
                   color="black", markeredgewidth=0.3)
            for size in legend_sizes
        ]
        labels = [str(value) for value in legend_probabilities]
    if has_zero:
        handles.append(Line2D([], [], linestyle="none", marker="o", markersize=np.sqrt(max_dot_size),
                              color="black", markeredgewidth=0.3))
        labels.append("0" if fixed_size else "0 (display cap)")
    if has_missing:
        handles.append(
            Line2D([], [], linestyle="none", marker="D", markersize=np.sqrt(fixed_dot_size),
                   color="black", markeredgewidth=0.3)
        )
        labels.append("Unavailable")
    legend_title = size_label if has_available else "Significance"
    if fixed_size and has_available:
        legend_title += "\nFixed size (no size encoding)"
    max_marker_size = max((handle.get_markersize() for handle in handles), default=0.0)
    from matplotlib.font_manager import FontProperties

    legend_font_size = FontProperties(size=plt.rcParams["legend.fontsize"]).get_size_in_points()
    legend = lax.legend(
        handles, labels, title=legend_title, loc="upper left", bbox_to_anchor=(0, 1),
        borderaxespad=0, frameon=False,
        # Include marker radii in the measured box and give rows enough height
        # for the largest dot, including caller-specified marker/font sizes.
        borderpad=max(0.4, max_marker_size / (2 * legend_font_size) + 0.1),
        handleheight=max(0.7, max_marker_size / (0.65 * legend_font_size)),
        handlelength=max(2.0, max_marker_size / legend_font_size),
    )
    layout = _layout_dotplot(fig, axes, lax, cax, legend, title_artist, region=owned_region,
                             n_events=len(label_order), n_receivers=len(receiver_order))

    return {
        "fig": fig,
        "ax": axes[0],
        "axes": axes,
        "legend": legend,
        "legend_ax": lax,
        "colorbar": cbar,
        "title_artist": title_artist,
        "layout": layout,
        "dot_sizes": dot_sizes.copy(),
        "size_encoding": {
            "statistic": significance_col,
            "fixed_size": fixed_size,
            "zero_display_cap": max_dot_size if has_zero else None,
            "positive_area_range": (min_dot_size, positive_max_size) if not positive.empty else None,
            "positive_legend_probabilities": legend_probabilities,
        },
        "plot_events": plot_events.drop(columns=["_event_id", "_x", "_y"]),
        "qc_failed_events": qc_failed_events,
        "qc_excluded_events": qc_excluded_events,
        "missing_score_events": missing_score_events,
        "qc": qc_summary,
        "selected_events": label_order,
        "selected_event_keys": event_order,
        "sender_order": sender_order,
        "receiver_order": receiver_order,
        "thresholds": {
            "top_n": top_n,
            "event_keys": event_keys,
            "min_cell_mesh_score": min_cell_mesh_score,
            "max_perm_pvalue": max_perm_pvalue,
            "max_fdr": max_fdr,
            "events_before_filter": len(source_events),
            "events_after_filter": len(filtered),
            "events_plotted": len(plot_events),
            **qc_summary,
        },
    }


def plot_sample_event_scores(
    result_or_sample_events: Any,
    *,
    sender: str,
    receiver: str,
    metabolite: Optional[str] = None,
    sensor_gene: str,
    hmdb_id: Optional[str] = None,
    sensor_type: Optional[str] = None,
    score_col: str = "cell_mesh_score",
    qc_only: bool = True,
    sample_order: Optional[list[str]] = None,
    color: str = "black",
    missing_color: str = "#777777",
    show_median: bool = True,
    ax: Any = None,
) -> dict[str, Any]:
    """Plot one sample-aware event score across samples.

    A bar is drawn for each retained non-missing sample-level score. With
    ``qc_only=False``, samples where the event is not computable remain visible
    and are explicitly marked ``NA`` rather than being represented as zero.
    This can occur when an endpoint is unobserved or the metabolite lacks
    production availability in that sample. An explicit ``hmdb_id`` and
    ``sensor_gene`` select the relation; sender/receiver identify its direction,
    and sample identifies each record. ``metabolite`` is an optional display
    override and never filters samples. Identical alias records count once;
    conflicting records for one sample raise an error.

    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    shows only sample rows that explicitly pass min-cells QC. Set
    ``qc_only=False`` to include failed rows, including any NA sample scores.
    Older sample-event tables without this column are unchanged.

    Scores are validated before comparing duplicate sample records. Native
    scores lie in [0, 1]; custom scores are finite and non-negative. Numeric
    strings are accepted, but malformed text/Inf/complex/bool values cannot
    become NA. Sender, receiver and sample-order labels share source stripping;
    duplicate normalized sample-order entries are rejected.
    """
    sample_events = _canonical_event_table(_sample_events_frame(result_or_sample_events))
    required = {
        "sample",
        "sender",
        "receiver",
        "hmdb_id",
        "sensor_gene",
        score_col,
    }
    if sensor_type is not None:
        required.add("sensor_type")
    missing = sorted(required.difference(sample_events.columns))
    if missing:
        raise KeyError(f"Missing required sample event columns: {', '.join(missing)}")

    selected_hmdb = _require_hmdb_id(hmdb_id)
    sensor_gene = _sensor_gene_id(sensor_gene)
    selectors = {
        "sender": _plot_labels([sender], "sender")[0],
        "receiver": _plot_labels([receiver], "receiver")[0],
        "hmdb_id": selected_hmdb,
        "sensor_gene": sensor_gene,
    }
    if sensor_type is not None:
        selectors["sensor_type"] = sensor_type

    selected_source = sample_events
    for column, value in selectors.items():
        selected_source = selected_source[
            selected_source[column].astype(str) == str(value)
        ]
    if selected_source.empty:
        raise ValueError("The requested sample-aware event was not found")

    selected, qc_failed_events, qc_excluded_events, qc_summary = _apply_min_cells_event_qc(
        selected_source,
        qc_only=qc_only,
    )
    if selected.empty:
        raise ValueError(
            "The requested sample-aware event has no rows after min-cells QC; "
            "set qc_only=False to include failed rows"
        )

    selected = selected.copy()
    selected[score_col] = _plot_numeric_series(selected[score_col], score_col,
                                              upper=1.0 if score_col == "cell_mesh_score" else None)
    event_columns = ["sender", "receiver", "hmdb_id", "sensor_gene"]
    matched_events = selected[event_columns].drop_duplicates()
    if len(matched_events) != 1:
        raise ValueError("The requested fields match multiple events")
    # Identical sample records with different display names count once; an
    # inconsistent value for one sample remains an input error, never a median.
    selected = selected.drop_duplicates(subset=[column for column in selected if column != "metabolite"])
    if selected["sample"].astype(str).duplicated().any():
        raise ValueError("The selected event contains conflicting duplicate rows for at least one sample")

    selected = selected.copy()
    selected["sample"] = selected["sample"].astype(str)
    available_samples = selected["sample"].tolist()
    if sample_order is None:
        order = sorted(available_samples)
    else:
        requested_order = _plot_labels(sample_order, "sample_order")
        order = [sample for sample in requested_order if sample in available_samples]
        order.extend(sorted(set(available_samples).difference(order)))
    selected["sample"] = pd.Categorical(selected["sample"], categories=order, ordered=True)
    selected = selected.sort_values("sample").reset_index(drop=True)

    valid = selected[score_col].notna()
    if not valid.any():
        raise ValueError("The requested event has no non-missing sample-level scores")
    median_score = float(selected.loc[valid, score_col].median())

    import numpy as np
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(max(6.0, 0.58 * len(selected) + 2.5), 3.8))
    fig = ax.figure
    positions = np.arange(len(selected))
    valid_positions = positions[valid.to_numpy()]
    missing_positions = positions[~valid.to_numpy()]

    ax.bar(
        valid_positions,
        selected.loc[valid, score_col].to_numpy(dtype=float),
        width=0.68,
        color=color,
        label="Sample score",
    )
    if len(missing_positions):
        ax.scatter(
            missing_positions,
            np.zeros(len(missing_positions)),
            marker="x",
            s=38,
            color=missing_color,
            linewidths=1.2,
            label="NA (event not computable)",
            zorder=3,
        )
        for position in missing_positions:
            ax.annotate(
                "NA",
                (position, 0.0),
                xytext=(0, 7),
                textcoords="offset points",
                ha="center",
                va="bottom",
                color=missing_color,
                fontsize=8,
            )
    if show_median:
        ax.axhline(
            median_score,
            color="#D55E00",
            linestyle="--",
            linewidth=1.2,
            label=f"Median = {median_score:.3g}",
        )

    event = selected.iloc[0]
    display_name = (
        metabolite if metabolite is not None
        else _metabolite_display_names(selected_source).get(selected_hmdb)
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(selected["sample"].astype(str), rotation=45, ha="right")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Cellmesh_score" if score_col == "cell_mesh_score" else score_col)
    ax.set_title(
        f"{_metabolite_label(display_name, selected_hmdb)} -> {sensor_gene} | "
        f"{event['sender']} -> {event['receiver']}"
    )
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="#e6e6e6", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()

    return {
        "fig": fig,
        "ax": ax,
        "sample_scores": selected,
        "qc_failed_events": qc_failed_events,
        "qc_excluded_events": qc_excluded_events,
        "qc": qc_summary,
        "event_key": matched_events.iloc[0].to_dict(),
        "median_score": median_score,
        "n_samples_coobserved": int(valid.sum()),
        "missing_samples": selected.loc[~valid, "sample"].astype(str).tolist(),
    }


def _result_adata_context(
    result: Any,
    adata: Any,
    *,
    cell_type_key: Optional[str],
    layer: Optional[str],
) -> tuple[str, Optional[str], pd.Series, Any]:
    """Resolve the expression matrix and cell labels used by a result."""
    from .preprocess import _validated_celltype_labels, _expression_source

    parameters = getattr(result, "parameters", {})
    if not isinstance(parameters, dict):
        parameters = {}

    resolved_cell_type_key = cell_type_key or parameters.get("cell_type_key", "cell_type")
    resolved_layer = layer if layer is not None else parameters.get("layer")

    obs = getattr(adata, "obs", None)
    if not isinstance(obs, pd.DataFrame):
        raise TypeError("adata must provide a pandas DataFrame in .obs")
    if resolved_cell_type_key not in obs:
        raise KeyError(f"{resolved_cell_type_key!r} not found in adata.obs")

    if resolved_layer is None:
        matrix = _expression_source(adata) if hasattr(type(adata), "X") or "X" in vars(adata) else None
        if matrix is None:
            raise TypeError("adata must provide an expression matrix in .X")
    else:
        layers = getattr(adata, "layers", None)
        if layers is None or resolved_layer not in layers:
            raise KeyError(f"{resolved_layer!r} not found in adata.layers")
        matrix = layers[resolved_layer]

    # 绘图重读 AnnData 时也必须沿用评分的标识规则；原对象仍保留原始名称。
    labels = _validated_celltype_labels(adata, resolved_cell_type_key)
    if matrix.shape[0] != len(labels):
        raise ValueError("The expression matrix row count does not match adata.obs")
    return resolved_cell_type_key, resolved_layer, labels, matrix


def _availability_intermediates(result: Any) -> list[dict[str, Any]]:
    availability = getattr(result, "availability_results", None)
    if not isinstance(availability, dict):
        raise TypeError(
            "Expected a CellMeshResult-like object with .availability_results"
        )

    intermediates = [availability]
    by_sample = availability.get("availability_by_sample")
    if isinstance(by_sample, dict):
        intermediates.extend(
            value for value in by_sample.values() if isinstance(value, dict)
        )
    return intermediates


def _first_availability_frame(result: Any, key: str) -> pd.DataFrame:
    for intermediate in _availability_intermediates(result):
        value = intermediate.get(key)
        if isinstance(value, pd.DataFrame) and not value.empty:
            return value
    raise ValueError(f"CellMeshResult does not contain a non-empty {key!r} intermediate")


class _UnknownPlotHMDB(ValueError):
    pass


def _metabolite_index_ids(table: pd.DataFrame) -> pd.Series:
    from .database import _normalize_hmdb_series

    if isinstance(table.index, pd.MultiIndex) and table.index.nlevels >= 2:
        level = "hmdb_id" if "hmdb_id" in table.index.names else 1
        values = table.index.get_level_values(level)
    elif table.index.name == "hmdb_id":
        values = table.index
    else:
        raise KeyError("Metabolite tables require an hmdb_id index; names cannot be matched")
    identifiers = _normalize_hmdb_series(pd.Series(values, dtype=object))
    if identifiers.isna().any():
        raise ValueError("Metabolite table hmdb_id must not be empty or missing")
    return identifiers


def _select_metabolite_row(
    table: pd.DataFrame,
    *,
    hmdb_id: Optional[str],
    table_name: str,
    cell_types: Optional[list[str]] = None,
) -> tuple[pd.Series, str]:
    identifier = _require_hmdb_id(hmdb_id)
    matches = table.iloc[_metabolite_index_ids(table).eq(identifier).to_numpy()]
    if matches.empty:
        raise _UnknownPlotHMDB(f"HMDB ID {identifier!r} was not found in {table_name}")
    # Alias rows may repeat identical values. Conflicting values for one ID
    # cannot be resolved by selecting a display name.
    from .preprocess import _normalized_label_series

    matches = matches.copy()
    matches.columns = pd.Index(_normalized_label_series(
        pd.Series(matches.columns, dtype=object), f"{table_name} cell types",
        reject_collisions=False,
    ), name=matches.columns.name)
    if matches.columns.duplicated().any():
        raise ValueError(f"{table_name} contains duplicate cell-type columns")
    if cell_types is not None:
        requested = _plot_labels(cell_types, f"{table_name} cell types")
        matches = matches.loc[:, matches.columns.isin(requested)]
    matches = matches.apply(lambda column: _plot_numeric_series(
        column, f"{table_name}[{column.name!r}]", upper=1.0 if table_name == "sender_scores" else None,
    ))
    if len(matches.drop_duplicates()) != 1:
        raise ValueError(f"Conflicting values for HMDB ID {identifier!r} in {table_name}")
    return matches.iloc[0].copy(), identifier


def _select_reaction_definitions(
    result: Any,
    *,
    hmdb_id: Optional[str],
) -> pd.DataFrame:
    reaction_genes = _first_availability_frame(result, "reaction_genes").copy()
    from .database import _normalize_hmdb_series

    required = {"hmdb_id", "reaction", "direction", "genes"}
    missing = sorted(required.difference(reaction_genes.columns))
    if missing:
        raise KeyError("reaction_genes is missing required columns: " + ", ".join(missing))
    identifier = _require_hmdb_id(hmdb_id)
    selected = reaction_genes.loc[_normalize_hmdb_series(reaction_genes["hmdb_id"]).eq(identifier)].copy()
    from .preprocess import _normalized_label_series

    # 防回退：不能只规范用于比较的临时键。后续取表达和计算必须使用相同的
    # 去空白、去重基因列表，否则空白基因会漏匹配，重复基因会被重复加权。
    # 保留未测到的基因用于反应身份判断；表达提取仍按既定规则跳过未测基因。
    selected["genes"] = selected["genes"].map(lambda genes: sorted(set(
        _normalized_label_series(pd.Series(_as_reaction_sequence(genes), dtype=object),
                                 "reaction_genes.genes", reject_collisions=False).tolist()
    )))
    selected["_genes_key"] = selected["genes"].map(tuple)
    selected = selected.drop_duplicates(["reaction", "direction", "_genes_key"])
    if selected.duplicated(["reaction", "direction"]).any():
        raise ValueError(f"Conflicting reaction definitions for HMDB ID {identifier!r}")
    selected = selected.drop(columns="_genes_key")
    if selected.empty:
        raise ValueError("No reaction definitions were found for the requested metabolite")
    if not (selected["direction"] == "product").any():
        raise ValueError("The requested metabolite has no production reaction definition")
    return selected.reset_index(drop=True)


def _ordered_cell_types(
    score_values: pd.Series,
    labels: pd.Series,
    requested: Optional[list[str]],
    *,
    role: str,
) -> tuple[list[str], pd.Series, list[str]]:
    from .preprocess import _normalized_label_series

    scores = score_values.copy()
    scores.index = pd.Index(
        _normalized_label_series(
            pd.Series(scores.index, dtype="object"),
            f"The {role} score table cell types",
            reject_collisions=False,
        ),
        name=scores.index.name,
    )
    if scores.index.duplicated().any():
        raise ValueError(f"The {role} score table contains duplicate cell-type columns")

    if requested is not None:
        requested = _plot_labels(requested, f"{role}_labels")
        # Numeric validation concerns the requested context, not unrelated types.
        numeric_mask = scores.index.isin(requested)
    else:
        numeric_mask = scores.index.isin(labels.astype(str))
    scores = scores.loc[numeric_mask]
    scores = _plot_numeric_series(scores, f"{role} score", upper=1.0)
    present = set(labels.astype(str))
    available = [
        cell_type
        for cell_type in scores.index
        if cell_type in present and pd.notna(scores.loc[cell_type])
    ]
    if requested is None:
        order = available
    else:
        order = requested
        missing = [cell_type for cell_type in order if cell_type not in available]
        if missing:
            raise ValueError(
                f"Unknown or unavailable {role} cell types: {', '.join(missing)}"
            )
    if not order:
        raise ValueError(f"No {role} cell types are available to plot")
    missing_scores = scores.index[scores.isna() & scores.index.isin(present)].tolist()
    return order, scores.reindex(order), missing_scores


def _expression_columns(
    adata: Any,
    matrix: Any,
    genes: list[str],
    cell_mask: Any,
    *,
    allow_missing: bool = False,
    layer: Optional[str] = None,
) -> tuple[Any, list[str]]:
    import numpy as np
    from scipy import sparse
    from .preprocess import (_normalized_gene_names, _validate_expression_values,
                             _slice_expression, _canonical_sparse_expression)

    var_names = _normalized_gene_names(getattr(adata, "var_names", []))
    if len(var_names) != matrix.shape[1]:
        raise ValueError("The expression matrix column count does not match adata.var_names")

    unique_genes = list(dict.fromkeys(str(gene) for gene in genes))
    locations = []
    measured_genes = []
    for gene in unique_genes:
        matched = np.flatnonzero(var_names.to_numpy() == gene)
        if len(matched) == 0:
            if allow_missing:
                continue
            raise KeyError(f"{gene!r} was not found in adata.var_names")
        locations.append(int(matched[0]))
        measured_genes.append(gene)

    selected = _slice_expression(matrix, cell_mask, locations)
    _validate_expression_values(selected, layer=layer)
    selected = _canonical_sparse_expression(selected)
    _validate_expression_values(selected, layer=layer)
    if sparse.issparse(selected):
        selected = selected.toarray()
    selected = np.asarray(selected, dtype=float)
    if selected.ndim == 1:
        selected = selected.reshape(-1, 1)
    return selected, measured_genes


def _as_reaction_sequence(value: Any) -> list[Any]:
    if isinstance(value, (list, tuple)):
        return list(value)
    if hasattr(value, "tolist") and not isinstance(value, str):
        converted = value.tolist()
        return converted if isinstance(converted, list) else [converted]
    return [value]


def _cell_reaction_scores(
    expression: Any,
    expression_genes: list[str],
    reaction_definitions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    import numpy as np
    from ._numerics import _check_numeric_result

    gene_locations = {gene: position for position, gene in enumerate(expression_genes)}
    reaction_columns: dict[str, Any] = {}
    reaction_metadata = []
    direction_totals = {
        "product": np.zeros(expression.shape[0], dtype=float),
        "substrate": np.zeros(expression.shape[0], dtype=float),
        "exporter": np.zeros(expression.shape[0], dtype=float),
    }

    for row_number, row in reaction_definitions.iterrows():
        genes = [str(gene) for gene in _as_reaction_sequence(row["genes"])]

        valid = [position for position, gene in enumerate(genes) if gene in gene_locations]
        if not valid:
            reaction_score = np.zeros(expression.shape[0], dtype=float)
            used_genes: list[str] = []
        else:
            used_genes = [genes[position] for position in valid]
            locations = [gene_locations[gene] for gene in used_genes]
            reaction_score = _reaction_activity(expression[:, locations], stage="plot reaction activity")
            _check_numeric_result(reaction_score, "plot reaction activity")

        direction = str(row["direction"])
        reaction_name = str(row["reaction"])
        column_name = f"{direction}:{reaction_name}:{row_number}"
        reaction_columns[column_name] = reaction_score
        if direction in direction_totals:
            with np.errstate(over="ignore", invalid="ignore"):
                direction_totals[direction] += reaction_score
        reaction_metadata.append(
            {
                "reaction_score_column": column_name,
                "reaction": reaction_name,
                "direction": direction,
                "genes": used_genes,
            }
        )

    for direction, capacity in direction_totals.items():
        _check_numeric_result(capacity, f"plot {direction} capacities")

    components = pd.DataFrame(
        {
            "production_score": direction_totals["product"],
            "consumption_score": direction_totals["substrate"],
            "export_score": direction_totals["exporter"],
        }
    )
    return components, pd.DataFrame(reaction_metadata)


def _plot_grouped_violins(
    plot_data: pd.DataFrame,
    *,
    group_col: str,
    value_col: str,
    group_order: list[str],
    color_values: pd.Series,
    color_value_name: str,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    ylabel: str,
    title: str,
    show_median: bool,
    ax: Any,
) -> dict[str, Any]:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from scipy.stats import gaussian_kde

    clean = plot_data.copy()
    clean[value_col] = _plot_numeric_series(clean[value_col], value_col)
    missing_value_data = clean.loc[clean[value_col].isna()].copy()
    clean = clean.loc[clean[value_col].notna()].copy()
    if clean.empty:
        raise ValueError("No finite single-cell values are available to plot")

    summary_rows = []
    for group in group_order:
        values = clean.loc[clean[group_col] == group, value_col]
        if values.empty:
            raise ValueError(f"No finite values are available for cell type {group!r}")
        summary_rows.append(
            {
                group_col: group,
                "n_cells": int(len(values)),
                f"mean_{value_col}": float(values.mean()),
                f"median_{value_col}": float(values.median()),
                f"q1_{value_col}": float(values.quantile(0.25)),
                f"q3_{value_col}": float(values.quantile(0.75)),
            }
        )
    summary = pd.DataFrame(summary_rows)

    colors = _plot_numeric_series(color_values.reindex(group_order), color_value_name,
                                  upper=1.0 if color_value_name == "metabolite_availability" else None)
    if colors.isna().any():
        missing = colors.index[colors.isna()].tolist()
        raise ValueError(
            f"Missing {color_value_name} values for cell types: {', '.join(missing)}"
        )
    if color_value_name not in summary:
        summary[color_value_name] = colors.to_numpy(dtype=float)

    color_min = float(colors.min()) if vmin is None else _plot_number(vmin, "vmin", lower=None)
    color_max = float(colors.max()) if vmax is None else _plot_number(vmax, "vmax", lower=None)
    if not np.isfinite(color_min) or not np.isfinite(color_max) or color_min > color_max:
        raise ValueError("vmin and vmax must define a finite increasing color range")
    if color_min == color_max:
        padding = max(abs(color_min) * 0.05, 1e-12)
        norm = Normalize(vmin=color_min - padding, vmax=color_max + padding)
    else:
        norm = Normalize(vmin=color_min, vmax=color_max)
    color_map = plt.get_cmap(cmap)

    if ax is None:
        _, ax = plt.subplots(figsize=(max(6.0, 0.82 * len(group_order) + 2.4), 4.5))
    fig = ax.figure
    violin_bodies = []
    constant_artists = []
    fallback_artists = []
    density_fallbacks = []
    medians = []
    for position, group in enumerate(group_order):
        values = clean.loc[clean[group_col] == group, value_col].to_numpy(dtype=float)
        color = color_map(norm(float(colors.loc[group])))
        if len(values) >= 2 and not np.all(values == values[0]):
            # Absolute allclose tolerances turn [0, 1e-9, 2e-9] into a false
            # zero constant. Fit KDE on a unit range, then restore actual units;
            # bandwidth is scale equivariant, and tiny variance cannot underflow.
            try:
                with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
                    low, high = float(values.min()), float(values.max())
                    span = high - low
                    normalized = (values - low) / span
                    grid = np.linspace(0.0, 1.0, 100)
                    density = gaussian_kde(normalized)(grid)
                    coords = low + grid * span
                if not (np.isfinite(density).all() and density.max() > 0
                        and np.isfinite(coords).all()):
                    raise ValueError("density or coordinates are not finite")
                coords[0], coords[-1] = low, high
                stats = dict(coords=coords, vals=density, mean=float(np.mean(values)),
                             median=float(np.median(values)), min=low, max=high)
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                fallback_artists.append(ax.scatter(np.full(len(values), position), values,
                                                   s=16, color=color, alpha=0.7, zorder=3))
                density_fallbacks.append({"cell_type": str(group), "reason": str(exc)})
            else:
                violin = ax.violin([stats], positions=[position], widths=0.82,
                                   showmeans=False, showmedians=False, showextrema=False)
                body = violin["bodies"][0]
                body.set_facecolor(color)
                body.set_edgecolor("#666666")
                body.set_linewidth(0.6)
                body.set_alpha(1.0)
                violin_bodies.append(body)
        else:
            artist = ax.hlines(
                float(values[0]),
                position - 0.3,
                position + 0.3,
                colors=[color],
                linewidth=3.0,
                zorder=3,
            )
            constant_artists.append(artist)

        if show_median:
            median_artist = ax.scatter(
                [position],
                [float(np.median(values))],
                s=24,
                facecolor="white",
                edgecolor="black",
                linewidth=0.7,
                zorder=4,
            )
            medians.append(median_artist)

    ax.set_xticks(np.arange(len(group_order)))
    ax.set_xticklabels(group_order, rotation=45, ha="right")
    ax.set_xlabel("Cell type")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color="#e6e6e6", linewidth=0.6)
    ax.set_axisbelow(True)

    scalar_mappable = ScalarMappable(norm=norm, cmap=color_map)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(scalar_mappable, ax=ax, pad=0.025)
    colorbar.set_label(color_value_name.replace("_", " ").title())
    fig.tight_layout()

    return {
        "fig": fig,
        "ax": ax,
        "plot_data": clean,
        "missing_value_data": missing_value_data,
        "summary": summary,
        "cell_order": group_order,
        "colorbar": colorbar,
        "violin_bodies": violin_bodies,
        "constant_artists": constant_artists,
        "fallback_artists": fallback_artists,
        "density_fallbacks": density_fallbacks,
        "median_artists": medians,
        "color_range": (color_min, color_max),
    }


def _mean_metabolite_intermediate(
    result: Any,
    key: str,
    *,
    hmdb_id: Optional[str],
    cell_types: Optional[list[str]] = None,
) -> Optional[pd.Series]:
    rows = []
    for intermediate in _availability_intermediates(result):
        table = intermediate.get(key)
        if not isinstance(table, pd.DataFrame) or table.empty:
            continue
        try:
            row, _ = _select_metabolite_row(
                table,
                hmdb_id=hmdb_id,
                table_name=key,
                cell_types=cell_types,
            )
        except _UnknownPlotHMDB:
            continue
        row = _plot_numeric_series(row, key)
        row.index = row.index.astype(str)
        rows.append(row)
    if not rows:
        return None
    return pd.concat(rows, axis=1).mean(axis=1, skipna=True)


@_with_numerical_diagnostics
def plot_metabolite_secretion_violin(
    result: Any,
    adata: Any,
    *,
    metabolite: Optional[str] = None,
    hmdb_id: Optional[str] = None,
    sender_labels: Optional[list[str]] = None,
    cell_type_key: Optional[str] = None,
    layer: Optional[str] = None,
    cmap: str = "Purples",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_median: bool = True,
    ax: Any = None,
) -> dict[str, Any]:
    """Plot one metabolite's secretion evidence across sender cell types.

    The violin geometry shows the single-cell distribution of CELL MESH
    production-reaction expression: genes within a reaction are combined by
    the same equal-weight geometric mean used by the scoring algorithm, then all
    production reactions for the metabolite are summed. As in scoring, only
    measured genes contribute, and reactions with no measured genes score zero.
    Violin fill color
    represents the formal abundance-adjusted cell-type
    ``metabolite_availability`` score. That score incorporates the configured
    sender abundance exponent, normalized production/consumption capacity, and
    the configured bounded exporter modulation.

    An explicit ``hmdb_id`` is required. Names never select a metabolite or
    reaction definition; ``metabolite`` is only an optional display override.
    The title includes the name and HMDB ID. This sender-only plot does not
    require a sensor gene. Identical alias rows/reactions contribute once;
    conflicting values for one identity raise an error.

    ``adata`` is required because ``CellMeshResult`` intentionally stores
    cell-type summaries rather than a copy of the single-cell expression
    matrix. By default, ``cell_type_key`` and ``layer`` are recovered from the
    result parameters. Passing ``sender_labels`` both filters and orders the
    displayed sender cell types.

    Native scores/expression fractions are checked in [0, 1]; single-cell
    expression, reaction activity and P/C/E capacities are finite non-negative
    values without an upper bound. Malformed numeric inputs raise explicit
    errors. With automatic cell-type selection, genuine missing scores are
    omitted and listed in ``missing_score_cell_types``; explicitly requesting
    an unavailable type raises an error. ``plot_data`` contains the actual
    plotted values; ``missing_value_data`` records any missing cell values.
    """
    sender_scores = getattr(result, "sender_scores", None)
    if not isinstance(sender_scores, pd.DataFrame) or sender_scores.empty:
        raise TypeError("Expected a CellMeshResult-like object with non-empty .sender_scores")
    availability, selected_hmdb = _select_metabolite_row(
        sender_scores,
        hmdb_id=hmdb_id,
        table_name="sender_scores",
        cell_types=sender_labels,
    )
    if metabolite is None:
        metabolite = availability.name[0] if isinstance(availability.name, tuple) else None
    reaction_definitions = _select_reaction_definitions(
        result,
        hmdb_id=selected_hmdb,
    )

    resolved_key, resolved_layer, labels, matrix = _result_adata_context(
        result,
        adata,
        cell_type_key=cell_type_key,
        layer=layer,
    )
    sender_order, availability, missing_score_cell_types = _ordered_cell_types(
        availability,
        labels,
        sender_labels,
        role="sender",
    )
    cell_mask = labels.isin(sender_order).to_numpy()

    reaction_gene_names = []
    for genes in reaction_definitions["genes"]:
        reaction_gene_names.extend(str(gene) for gene in _as_reaction_sequence(genes))
    # 完整反应定义包含未测到的基因。保持反应选择结果，仅在表达取列时跳过
    # 未测基因，与评分的顺序一致；不可在反应去重/子集判断前删基因。
    expression, expression_genes = _expression_columns(
        adata,
        matrix,
        reaction_gene_names,
        cell_mask,
        allow_missing=True,
        layer=resolved_layer,
    )
    components, reaction_metadata = _cell_reaction_scores(
        expression,
        expression_genes,
        reaction_definitions,
    )

    selected_obs = adata.obs.loc[cell_mask]
    plot_data = pd.DataFrame(
        {
            "cell_id": selected_obs.index.astype(str),
            "sender": labels.loc[cell_mask].to_numpy(),
            "metabolite": metabolite,
            "hmdb_id": selected_hmdb,
        }
    )
    plot_data = pd.concat([plot_data.reset_index(drop=True), components], axis=1)

    plotted = _plot_grouped_violins(
        plot_data,
        group_col="sender",
        value_col="production_score",
        group_order=sender_order,
        color_values=availability,
        color_value_name="metabolite_availability",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        ylabel="Single-cell production reaction score",
        title=f"{_metabolite_label(metabolite, selected_hmdb)} secretion evidence in sender cell types",
        show_median=show_median,
        ax=ax,
    )

    summary = plotted["summary"].copy()
    intermediate_columns = {
        "P": "abundance_adjusted_production",
        "C": "abundance_adjusted_consumption",
        "E": "abundance_adjusted_export",
    }
    for key, column_name in intermediate_columns.items():
        values = _mean_metabolite_intermediate(
            result,
            key,
            hmdb_id=selected_hmdb,
            cell_types=sender_order,
        )
        if values is not None:
            summary[column_name] = summary["sender"].map(values)

    plotted.update(
        {
            "missing_score_cell_types": missing_score_cell_types,
            "summary": summary,
            "metabolite": metabolite,
            "hmdb_id": selected_hmdb,
            "reaction_definitions": reaction_definitions,
            "reaction_metadata": reaction_metadata,
            "cell_type_key": resolved_key,
            "layer": resolved_layer,
        }
    )
    return plotted


@_with_numerical_diagnostics
def plot_receptor_expression_violin(
    result: Any,
    adata: Any,
    *,
    receptor_gene: str,
    metabolite: Optional[str] = None,
    hmdb_id: Optional[str] = None,
    receiver_labels: Optional[list[str]] = None,
    cell_type_key: Optional[str] = None,
    layer: Optional[str] = None,
    cmap: str = "Reds",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_median: bool = True,
    ax: Any = None,
) -> dict[str, Any]:
    """Plot single-cell receptor expression across receiver cell types.

    Violin fill color represents mean expression in each receiver cell type,
    matching the visual encoding of MEBOCOST's ``violin_plot``. The returned
    summary also includes CELL MESH sensor scores and expression fractions.
    An explicit ``hmdb_id`` and ``receptor_gene`` select the context using the
    same ID normalization as scoring. ``metabolite`` is an optional display
    override, never a matching field. The title includes the metabolite name,
    HMDB ID and receptor gene. Name-only and gene-only context lookup are not
    supported. Empty/missing or unmatched IDs raise an error. Identical alias
    rows count once; conflicting receiver scores raise an error.
    Passing ``receiver_labels`` both filters and orders the displayed types.

    Native scores/expression fractions are checked in [0, 1]; single-cell
    expression, reaction activity and P/C/E capacities are finite non-negative
    values without an upper bound. Malformed numeric inputs raise explicit
    errors. With automatic cell-type selection, genuine missing scores are
    omitted and listed in ``missing_score_cell_types``; explicitly requesting
    an unavailable type raises an error. ``plot_data`` contains the actual
    plotted values; ``missing_value_data`` records any missing cell values.
    """
    receiver_scores = getattr(result, "receiver_scores", None)
    if not isinstance(receiver_scores, pd.DataFrame) or receiver_scores.empty:
        raise TypeError(
            "Expected a CellMeshResult-like object with non-empty .receiver_scores"
        )
    required = {"receiver", "sensor_gene", "sensor_score", "sensor_expr_frac"}
    missing = sorted(required.difference(receiver_scores.columns))
    if missing:
        raise KeyError(
            "receiver_scores is missing required columns: " + ", ".join(missing)
        )

    from .database import _normalize_hmdb_series
    from .preprocess import _normalized_label_series

    selected_hmdb = _require_hmdb_id(hmdb_id)
    receptor_gene = _sensor_gene_id(receptor_gene)
    if "hmdb_id" not in receiver_scores:
        raise KeyError("receiver_scores does not contain an 'hmdb_id' column")
    genes = _normalized_label_series(receiver_scores["sensor_gene"], "sensor_gene", reject_collisions=False)
    selected_scores = receiver_scores.loc[
        genes.eq(receptor_gene)
        & _normalize_hmdb_series(receiver_scores["hmdb_id"]).eq(selected_hmdb)
    ].copy()
    if selected_scores.empty:
        raise ValueError(f"Receptor {receptor_gene!r} with HMDB ID {selected_hmdb!r} was not found in receiver_scores")
    if metabolite is None and "metabolite" in selected_scores:
        candidates = selected_scores["metabolite"].dropna()
        metabolite = next((value for value in candidates if str(value).strip()), None)
    selected_scores["receiver"] = _normalized_label_series(selected_scores["receiver"], "receiver")
    if receiver_labels is not None:
        receiver_labels = _plot_labels(receiver_labels, "receiver_labels")
        selected_scores = selected_scores.loc[selected_scores["receiver"].isin(receiver_labels)].copy()
    for column in ("sensor_score", "sensor_expr_frac"):
        selected_scores[column] = _plot_numeric_series(selected_scores[column], column, upper=1.0,
                                                       probability=column == "sensor_expr_frac")
    comparable = selected_scores.assign(hmdb_id=selected_hmdb, sensor_gene=receptor_gene).drop(
        columns="metabolite", errors="ignore",
    )
    selected_scores = selected_scores.loc[~comparable.duplicated()]
    if selected_scores["receiver"].duplicated().any():
        raise ValueError(f"Conflicting receiver scores for HMDB ID {selected_hmdb!r} and receptor {receptor_gene!r}")
    receiver_metadata = selected_scores.set_index("receiver")[["sensor_score", "sensor_expr_frac"]]

    resolved_key, resolved_layer, labels, matrix = _result_adata_context(
        result,
        adata,
        cell_type_key=cell_type_key,
        layer=layer,
    )
    receiver_order, _, missing_score_cell_types = _ordered_cell_types(
        receiver_metadata["sensor_score"],
        labels,
        receiver_labels,
        role="receiver",
    )
    cell_mask = labels.isin(receiver_order).to_numpy()
    expression, _ = _expression_columns(
        adata,
        matrix,
        [str(receptor_gene)],
        cell_mask,
        layer=resolved_layer,
    )

    selected_obs = adata.obs.loc[cell_mask]
    plot_data = pd.DataFrame(
        {
            "cell_id": selected_obs.index.astype(str),
            "receiver": labels.loc[cell_mask].to_numpy(),
            "receptor_gene": str(receptor_gene),
            "receptor_expression": expression[:, 0],
        }
    )
    mean_expression = (
        plot_data.groupby("receiver", sort=False)["receptor_expression"].mean()
    )
    plotted = _plot_grouped_violins(
        plot_data,
        group_col="receiver",
        value_col="receptor_expression",
        group_order=receiver_order,
        color_values=mean_expression,
        color_value_name="mean_receptor_expression",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        ylabel="Single-cell receptor expression",
        title=f"{_metabolite_label(metabolite, selected_hmdb)} -> {receptor_gene} expression in receiver cell types",
        show_median=show_median,
        ax=ax,
    )
    summary = plotted["summary"].copy()
    summary["sensor_score"] = summary["receiver"].map(receiver_metadata["sensor_score"])
    summary["sensor_expr_frac"] = summary["receiver"].map(
        receiver_metadata["sensor_expr_frac"]
    )
    plotted.update(
        {
            "missing_score_cell_types": missing_score_cell_types,
            "summary": summary,
            "receptor_gene": str(receptor_gene),
            "selected_receiver_scores": selected_scores,
            "hmdb_id": selected_hmdb,
            "event_key": (selected_hmdb, receptor_gene),
            "cell_type_key": resolved_key,
            "layer": resolved_layer,
        }
    )
    return plotted


__all__ = [
    "plot_significant_event_counts",
    "plot_communication_network",
    "plot_event_dotplot",
    "plot_sample_event_scores",
    "plot_metabolite_secretion_violin",
    "plot_receptor_expression_violin",
]
