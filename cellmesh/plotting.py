"""Canonical plotting utilities for CELL MESH.

Plotting functions import optional visualization dependencies, such as
matplotlib, inside the function body so the core package import remains fast.
"""
from __future__ import annotations

from numbers import Integral
from typing import Any, Optional

import pandas as pd


_MIN_CELLS_QC_COLUMN = "passes_min_cells"


def _validate_optional_positive_int(value: Any, *, name: str) -> Optional[int]:
    """Validate an optional positive integer plotting limit."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be a positive integer or None")
    if value <= 0:
        raise ValueError(f"{name} must be greater than zero")
    return int(value)


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

    Counts are based on unique ``sender + receiver + metabolite + hmdb_id +
    sensor_gene`` combinations after optional score, permutation p-value, and
    FDR filters. By default, the plot summarizes events with
    ``fdr_sensor_type <= 0.05``.
    Set ``fdr_col`` to ``"fdr_global"`` or ``"fdr_sensor_type"`` in either
    inference mode to select the desired correction scope.
    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    excludes rows that do not explicitly pass that QC. Set ``qc_only=False``
    to display all rows. Results without the QC column retain their historical
    behavior.
    """
    source_events = _events_frame(result_or_events).copy()
    if source_events.empty:
        raise ValueError("No events available to plot")

    required = {"sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "cell_mesh_score"}
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
    filtered = events
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

    if unique_keys is None:
        unique_keys = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene"]
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

    Events are deduplicated by ``unique_keys`` and then aggregated for each
    directed sender-receiver pair. Edge width represents the number of unique
    communication events. By default, edge color represents the sum of
    ``score_col`` and node area represents the number of distinct cell types
    connected to that node in either direction.

    ``node_size_by`` accepts ``"connected_cell_type_count"``,
    ``"connection_count"``, ``"event_count"``, ``"total_score"``, or
    ``None`` for fixed-size nodes. ``edge_color_by`` accepts
    ``"total_score"``, ``"mean_score"``, ``"event_count"``, or ``None``
    for ``fixed_edge_color``. The returned node and edge tables contain the
    exact aggregate and visual values used in the plot.

    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    uses only rows that explicitly pass min-cells QC. Set ``qc_only=False`` to
    include all rows. Older event tables without this column are unchanged.
    """
    source_events = _events_frame(result_or_events).copy()
    if source_events.empty:
        raise ValueError("No events available to plot")

    required = {
        "sender",
        "receiver",
        "metabolite",
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

    def _validate_output_range(name: str, values: tuple[float, float]) -> None:
        if len(values) != 2 or values[0] <= 0 or values[1] < values[0]:
            raise ValueError(f"{name} must be a positive (minimum, maximum) pair")

    _validate_output_range("node_size_range", node_size_range)
    _validate_output_range("edge_width_range", edge_width_range)
    if not 0 <= edge_alpha <= 1:
        raise ValueError("edge_alpha must be between 0 and 1")
    if curve < 0:
        raise ValueError("curve must be non-negative")

    filtered = events.copy()
    if sender_labels is not None:
        filtered = filtered[filtered["sender"].isin(sender_labels)]
    if receiver_labels is not None:
        filtered = filtered[filtered["receiver"].isin(receiver_labels)]
    filtered[score_col] = pd.to_numeric(filtered[score_col], errors="coerce")
    filtered = filtered[filtered[score_col].notna()]
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

    if unique_keys is None:
        unique_keys = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene"]
    missing_unique = sorted(set(unique_keys).difference(source_events.columns))
    if missing_unique:
        raise KeyError(f"Missing unique key columns: {', '.join(missing_unique)}")
    if not {"sender", "receiver"}.issubset(unique_keys):
        raise ValueError("unique_keys must include sender and receiver")

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
            scaled = np.full(array.shape, (output_min + output_max) / 2.0)
        else:
            clipped = np.clip(array, value_min, value_max)
            fraction = (clipped - value_min) / (value_max - value_min)
            scaled = output_min + fraction * (output_max - output_min)
        return scaled, (value_min, value_max)

    if node_size_by is None:
        node_table["node_size_value"] = np.nan
        node_table["node_size"] = (node_size_range[0] + node_size_range[1]) / 2.0
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


def _event_label(df: pd.DataFrame) -> pd.Series:
    return df["metabolite"].astype(str) + " -> " + df["sensor_gene"].astype(str)


def _event_key_label(metabolite: Any, hmdb_id: Any, sensor_gene: Any) -> str:
    return f"{metabolite} | {hmdb_id} | {sensor_gene}"


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
    size is ``-log10(fdr_col)``. If ``event_keys`` is not provided, the top
    ``top_n`` events are selected by lowest FDR and highest score across all
    sender-receiver pairs.

    ``event_keys`` accepts labels produced by this function
    (``"metabolite -> sensor_gene"`` or ``"metabolite | hmdb_id | sensor_gene"``),
    ``(metabolite, hmdb_id,
    sensor_gene)`` tuples, or dictionaries containing those three fields.
    When supplying an existing ``ax``, the filtered data must contain exactly
    one sender; use ``sender_labels`` to select it. Leave ``ax=None`` for the
    default multi-sender faceted layout.
    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    excludes rows that do not explicitly pass min-cells QC. Set
    ``qc_only=False`` to include all rows. Older event tables without this
    column retain their historical behavior.
    """
    top_n = _validate_optional_positive_int(top_n, name="top_n")
    source_events = _events_frame(result_or_events).copy()
    if source_events.empty:
        raise ValueError("No events available to plot")

    required = {
        "sender",
        "receiver",
        "metabolite",
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
        filtered = filtered[filtered["sender"].isin(sender_labels)]
    if receiver_labels is not None:
        filtered = filtered[filtered["receiver"].isin(receiver_labels)]
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

    filtered = filtered.assign(
        _event_label=_event_label(filtered),
        _event_key=(
            filtered["metabolite"].astype(str)
            + " | "
            + filtered["hmdb_id"].astype(str)
            + " | "
            + filtered["sensor_gene"].astype(str)
        ),
    )

    if event_keys is not None:
        selected_labels: list[str] = []
        selected_keys: list[str] = []
        for key in event_keys:
            if isinstance(key, str):
                if " | " in key:
                    selected_keys.append(key)
                else:
                    selected_labels.append(key)
            elif isinstance(key, dict):
                selected_keys.append(
                    _event_key_label(key["metabolite"], key["hmdb_id"], key["sensor_gene"])
                )
            else:
                metabolite, hmdb_id, sensor_gene = key
                selected_keys.append(_event_key_label(metabolite, hmdb_id, sensor_gene))
        plot_events = filtered[
            filtered["_event_label"].isin(selected_labels)
            | filtered["_event_key"].isin(selected_keys)
        ].copy()
        if plot_events.empty:
            raise ValueError("None of the requested event_keys were found after filtering")
    else:
        rank_cols = ["_event_label"]
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
        plot_events = filtered[filtered["_event_label"].isin(event_rank["_event_label"])].copy()

    # Keep one event per row/sender/receiver if duplicate rows exist, prioritizing lower FDR and higher score.
    sort_cols = ["_event_label", "sender", "receiver"]
    ascending = [True, True, True]
    if fdr_col in plot_events.columns:
        sort_cols.append(fdr_col)
        ascending.append(True)
    sort_cols.append(score_col)
    ascending.append(False)
    plot_events = (
        plot_events
        .sort_values(sort_cols, ascending=ascending, na_position="last")
        .drop_duplicates(["_event_label", "sender", "receiver"], keep="first")
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
    if event_keys is not None:
        requested = selected_labels + filtered.loc[
            filtered["_event_key"].isin(selected_keys), "_event_label"
        ].drop_duplicates().tolist()
        label_order = [label for label in requested if label in set(plot_events["_event_label"])]
    elif fdr_col in plot_events.columns:
        label_order = (
            plot_events.groupby("_event_label")
            .agg(_best_fdr=(fdr_col, "min"), _best_score=(score_col, "max"))
            .sort_values(["_best_fdr", "_best_score"], ascending=[True, False], na_position="last")
            .index
            .tolist()
        )
    else:
        label_order = (
            plot_events.groupby("_event_label")[score_col]
            .max()
            .sort_values(ascending=False)
            .index
            .tolist()
        )

    plot_events["sender"] = pd.Categorical(plot_events["sender"], categories=sender_order, ordered=True)
    plot_events["receiver"] = pd.Categorical(
        plot_events["receiver"], categories=receiver_order, ordered=True
    )
    plot_events["_event_label"] = pd.Categorical(
        plot_events["_event_label"], categories=label_order, ordered=True
    )
    plot_events["_x"] = plot_events["receiver"].cat.codes
    plot_events["_y"] = plot_events["_event_label"].cat.codes

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Rectangle

    if fdr_col in plot_events.columns:
        fdr_values = plot_events[fdr_col].astype(float)
        positive = fdr_values[fdr_values > 0]
        floor = positive.min() if not positive.empty else 1e-12
        significance = -np.log10(fdr_values.clip(lower=floor))
        size_label = "FDR"
        size_legend_transform = lambda values: np.power(10.0, -values)
    else:
        p_values = plot_events[pvalue_col].astype(float)
        positive = p_values[p_values > 0]
        floor = positive.min() if not positive.empty else 1e-12
        significance = -np.log10(p_values.clip(lower=floor))
        size_label = f"-log10({pvalue_col})"
        size_legend_transform = None

    sig_min = float(np.nanmin(significance)) if len(significance) else 0.0
    sig_max = float(np.nanmax(significance)) if len(significance) else 0.0
    if sig_max > sig_min:
        dot_sizes = min_dot_size + (significance - sig_min) / (sig_max - sig_min) * (
            max_dot_size - min_dot_size
        )
    else:
        dot_sizes = pd.Series((min_dot_size + max_dot_size) / 2, index=plot_events.index)

    if ax is None:
        max_label_len = max(
            (len(str(label)) for label in plot_events["_event_label"].astype(str).unique()),
            default=0,
        )
        fig = plt.figure(
            figsize=(
                max(9.0, 2.15 * len(sender_order) + 3.4 + 0.035 * max_label_len),
                max(4.0, 0.42 * len(label_order) + 2.3),
            )
        )
        grid = GridSpec(
            2,
            len(sender_order) + 2,
            figure=fig,
            width_ratios=[1.0] * len(sender_order) + [0.20, 0.32],
            height_ratios=[0.42, 0.58],
            wspace=0.18,
            hspace=0.28,
        )
        axes = [fig.add_subplot(grid[:, i]) for i in range(len(sender_order))]
        lax = fig.add_subplot(grid[0, -1])
        cax = fig.add_subplot(grid[1, -1])
    else:
        if len(sender_order) != 1:
            raise ValueError(
                "A supplied ax can display exactly one sender; use "
                "sender_labels to select one sender or leave ax=None"
            )
        axes = [ax]
        fig = ax.figure
        cax = None
        lax = None

    score_values = plot_events[score_col].astype(float)
    score_min = float(score_values.min())
    score_max = float(score_values.max())
    scatter = None
    for idx, sender in enumerate(sender_order):
        axis = axes[idx]
        panel = plot_events[plot_events["sender"] == sender]
        panel_sizes = dot_sizes.loc[panel.index]
        scatter = axis.scatter(
            panel["_x"],
            panel["_y"],
            c=panel[score_col].astype(float),
            s=panel_sizes,
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

    if cax is not None:
        cbar = fig.colorbar(scatter, cax=cax)
    else:
        cbar = fig.colorbar(scatter, ax=axes)
    cbar.ax.set_title(score_label, fontsize=9, pad=6)
    cbar.ax.tick_params(direction="in")

    legend_values = np.linspace(sig_min, sig_max, num=3) if sig_max > sig_min else np.array([sig_max])
    legend_sizes = (
        min_dot_size + (legend_values - sig_min) / (sig_max - sig_min) * (max_dot_size - min_dot_size)
        if sig_max > sig_min
        else np.array([(min_dot_size + max_dot_size) / 2])
    )
    handles = [
        axes[0].scatter([], [], s=size, color="black", edgecolors="black", linewidths=0.3)
        for size in legend_sizes
    ]
    if size_legend_transform is not None:
        labels = [f"{value:.2g}" for value in size_legend_transform(legend_values)]
    else:
        labels = [f"{value:.2g}" for value in legend_values]
    if lax is not None:
        lax.axis("off")
        lax.legend(handles, labels, title=size_label, loc="upper left", frameon=False)
    else:
        axes[0].legend(handles, labels, title=size_label, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    fig.suptitle("Metabolite-sensor communication events", y=0.995)
    if ax is None:
        max_label_len = max((len(str(label)) for label in label_order), default=0)
        left_margin = min(0.52, max(0.22, 0.0105 * max_label_len))
        fig.subplots_adjust(left=left_margin, right=0.93, top=0.86, bottom=0.24)

    return {
        "fig": fig,
        "ax": axes[0],
        "axes": axes,
        "plot_events": plot_events.drop(columns=["_x", "_y"]),
        "qc_failed_events": qc_failed_events,
        "qc_excluded_events": qc_excluded_events,
        "qc": qc_summary,
        "selected_events": label_order,
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
    metabolite: str,
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
    production availability in that sample. If the supplied event fields match
    more than one HMDB identifier or sensor type, provide the corresponding
    optional argument to disambiguate the event.

    When ``passes_min_cells`` is available, ``qc_only=True`` (the default)
    shows only sample rows that explicitly pass min-cells QC. Set
    ``qc_only=False`` to include failed rows, including any NA sample scores.
    Older sample-event tables without this column are unchanged.
    """
    sample_events = _sample_events_frame(result_or_sample_events).copy()
    required = {
        "sample",
        "sender",
        "receiver",
        "metabolite",
        "hmdb_id",
        "sensor_gene",
        "sensor_type",
        score_col,
    }
    missing = sorted(required.difference(sample_events.columns))
    if missing:
        raise KeyError(f"Missing required sample event columns: {', '.join(missing)}")

    selectors = {
        "sender": sender,
        "receiver": receiver,
        "metabolite": metabolite,
        "sensor_gene": sensor_gene,
    }
    if hmdb_id is not None:
        selectors["hmdb_id"] = hmdb_id
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

    event_columns = [
        "sender",
        "receiver",
        "metabolite",
        "hmdb_id",
        "sensor_gene",
        "sensor_type",
    ]
    matched_events = selected[event_columns].drop_duplicates()
    if len(matched_events) != 1:
        raise ValueError(
            "The requested fields match multiple events; provide hmdb_id and/or "
            "sensor_type to select exactly one event"
        )
    if selected["sample"].astype(str).duplicated().any():
        raise ValueError("The selected event contains duplicate rows for at least one sample")

    selected = selected.copy()
    selected["sample"] = selected["sample"].astype(str)
    selected[score_col] = pd.to_numeric(selected[score_col], errors="coerce")
    available_samples = selected["sample"].tolist()
    if sample_order is None:
        order = sorted(available_samples)
    else:
        requested_order = [str(sample) for sample in sample_order]
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

    event = matched_events.iloc[0]
    ax.set_xticks(positions)
    ax.set_xticklabels(selected["sample"].astype(str), rotation=45, ha="right")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Cellmesh_score" if score_col == "cell_mesh_score" else score_col)
    ax.set_title(
        f"{event['metabolite']} -> {event['sensor_gene']} | "
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
        "event_key": event.to_dict(),
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
        matrix = getattr(adata, "X", None)
        if matrix is None:
            raise TypeError("adata must provide an expression matrix in .X")
    else:
        layers = getattr(adata, "layers", None)
        if layers is None or resolved_layer not in layers:
            raise KeyError(f"{resolved_layer!r} not found in adata.layers")
        matrix = layers[resolved_layer]

    labels = obs[resolved_cell_type_key]
    if labels.isna().any():
        raise ValueError(f"adata.obs[{resolved_cell_type_key!r}] contains missing values")
    labels = labels.astype(str)
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


def _index_matches_metabolite(
    index_value: Any,
    *,
    metabolite: str,
    hmdb_id: Optional[str],
) -> bool:
    if isinstance(index_value, tuple):
        candidate_metabolite = index_value[0]
        candidate_hmdb = index_value[1] if len(index_value) > 1 else None
    else:
        candidate_metabolite = index_value
        candidate_hmdb = None
    if str(candidate_metabolite) != str(metabolite):
        return False
    return hmdb_id is None or str(candidate_hmdb) == str(hmdb_id)


def _select_metabolite_row(
    table: pd.DataFrame,
    *,
    metabolite: str,
    hmdb_id: Optional[str],
    table_name: str,
) -> tuple[pd.Series, Optional[str]]:
    matches = [
        position
        for position, index_value in enumerate(table.index.tolist())
        if _index_matches_metabolite(
            index_value,
            metabolite=metabolite,
            hmdb_id=hmdb_id,
        )
    ]
    if not matches:
        suffix = "" if hmdb_id is None else f" and HMDB ID {hmdb_id!r}"
        raise ValueError(f"{metabolite!r}{suffix} was not found in {table_name}")

    matched_keys = [table.index[position] for position in matches]
    canonical_keys = {
        (
            str(key[0]),
            str(key[1]) if isinstance(key, tuple) and len(key) > 1 else None,
        )
        if isinstance(key, tuple)
        else (str(key), None)
        for key in matched_keys
    }
    if len(canonical_keys) != 1 or len(matches) != 1:
        raise ValueError(
            f"{metabolite!r} matches multiple rows in {table_name}; provide hmdb_id"
        )

    index_value = matched_keys[0]
    selected_hmdb = (
        str(index_value[1])
        if isinstance(index_value, tuple) and len(index_value) > 1
        else hmdb_id
    )
    return table.iloc[matches[0]].copy(), selected_hmdb


def _select_reaction_definitions(
    result: Any,
    *,
    metabolite: str,
    hmdb_id: Optional[str],
) -> pd.DataFrame:
    reaction_genes = _first_availability_frame(result, "reaction_genes").copy()
    required = {"metabolite", "hmdb_id", "reaction", "direction", "genes"}
    missing = sorted(required.difference(reaction_genes.columns))
    if missing:
        raise KeyError(
            "reaction_genes is missing required columns: " + ", ".join(missing)
        )

    selected = reaction_genes[
        reaction_genes["metabolite"].astype(str) == str(metabolite)
    ]
    if hmdb_id is not None:
        selected = selected[selected["hmdb_id"].astype(str) == str(hmdb_id)]
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
) -> tuple[list[str], pd.Series]:
    scores = pd.to_numeric(score_values, errors="coerce").copy()
    scores.index = scores.index.astype(str)
    if scores.index.duplicated().any():
        raise ValueError(f"The {role} score table contains duplicate cell-type columns")

    present = set(labels.astype(str))
    available = [
        cell_type
        for cell_type in scores.index
        if cell_type in present and pd.notna(scores.loc[cell_type])
    ]
    if requested is None:
        order = available
    else:
        order = [str(cell_type) for cell_type in requested]
        if len(order) != len(set(order)):
            raise ValueError(f"{role}_labels must not contain duplicates")
        missing = [cell_type for cell_type in order if cell_type not in available]
        if missing:
            raise ValueError(
                f"Unknown or unavailable {role} cell types: {', '.join(missing)}"
            )
    if not order:
        raise ValueError(f"No {role} cell types are available to plot")
    return order, scores.reindex(order)


def _expression_columns(
    adata: Any,
    matrix: Any,
    genes: list[str],
    cell_mask: Any,
) -> tuple[Any, list[str]]:
    import numpy as np
    from scipy import sparse

    var_names = pd.Index(getattr(adata, "var_names", [])).astype(str)
    if len(var_names) != matrix.shape[1]:
        raise ValueError("The expression matrix column count does not match adata.var_names")

    unique_genes = list(dict.fromkeys(str(gene) for gene in genes))
    locations = []
    for gene in unique_genes:
        matched = np.flatnonzero(var_names.to_numpy() == gene)
        if len(matched) == 0:
            raise KeyError(f"{gene!r} was not found in adata.var_names")
        if len(matched) > 1:
            raise ValueError(f"{gene!r} occurs more than once in adata.var_names")
        locations.append(int(matched[0]))

    selected = matrix[cell_mask, :][:, locations]
    if sparse.issparse(selected):
        selected = selected.toarray()
    selected = np.asarray(selected, dtype=float)
    if selected.ndim == 1:
        selected = selected.reshape(-1, 1)
    return selected, unique_genes


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

    if np.any(np.isfinite(expression) & (expression < 0)):
        raise ValueError(
            "Metabolite secretion violins require a non-negative expression matrix"
        )

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
            reaction_score = np.expm1(
                np.mean(np.log1p(expression[:, locations]), axis=1)
            )

        direction = str(row["direction"])
        reaction_name = str(row["reaction"])
        column_name = f"{direction}:{reaction_name}:{row_number}"
        reaction_columns[column_name] = reaction_score
        if direction in direction_totals:
            direction_totals[direction] += reaction_score
        reaction_metadata.append(
            {
                "reaction_score_column": column_name,
                "reaction": reaction_name,
                "direction": direction,
                "genes": used_genes,
            }
        )

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

    clean = plot_data.copy()
    clean[value_col] = pd.to_numeric(clean[value_col], errors="coerce")
    clean = clean[np.isfinite(clean[value_col])]
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

    colors = pd.to_numeric(color_values.reindex(group_order), errors="coerce")
    if colors.isna().any():
        missing = colors.index[colors.isna()].tolist()
        raise ValueError(
            f"Missing {color_value_name} values for cell types: {', '.join(missing)}"
        )
    if color_value_name not in summary:
        summary[color_value_name] = colors.to_numpy(dtype=float)

    color_min = float(colors.min()) if vmin is None else float(vmin)
    color_max = float(colors.max()) if vmax is None else float(vmax)
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
    medians = []
    for position, group in enumerate(group_order):
        values = clean.loc[clean[group_col] == group, value_col].to_numpy(dtype=float)
        color = color_map(norm(float(colors.loc[group])))
        if len(values) >= 2 and not np.allclose(values, values[0]):
            violin = ax.violinplot(
                [values],
                positions=[position],
                widths=0.82,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
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
        "summary": summary,
        "cell_order": group_order,
        "colorbar": colorbar,
        "violin_bodies": violin_bodies,
        "constant_artists": constant_artists,
        "median_artists": medians,
        "color_range": (color_min, color_max),
    }


def _mean_metabolite_intermediate(
    result: Any,
    key: str,
    *,
    metabolite: str,
    hmdb_id: Optional[str],
) -> Optional[pd.Series]:
    rows = []
    for intermediate in _availability_intermediates(result):
        table = intermediate.get(key)
        if not isinstance(table, pd.DataFrame) or table.empty:
            continue
        try:
            row, _ = _select_metabolite_row(
                table,
                metabolite=metabolite,
                hmdb_id=hmdb_id,
                table_name=key,
            )
        except ValueError:
            continue
        row = pd.to_numeric(row, errors="coerce")
        row.index = row.index.astype(str)
        rows.append(row)
    if not rows:
        return None
    return pd.concat(rows, axis=1).mean(axis=1, skipna=True)


def plot_metabolite_secretion_violin(
    result: Any,
    adata: Any,
    *,
    metabolite: str,
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
    production reactions for the metabolite are summed. Violin fill color
    represents the formal abundance-adjusted cell-type
    ``metabolite_availability`` score. That score incorporates the configured
    sender abundance exponent and normalized production/consumption capacity.
    Export capacity is available separately as support evidence and does not
    enter the formal sender score.

    ``adata`` is required because ``CellMeshResult`` intentionally stores
    cell-type summaries rather than a copy of the single-cell expression
    matrix. By default, ``cell_type_key`` and ``layer`` are recovered from the
    result parameters. Passing ``sender_labels`` both filters and orders the
    displayed sender cell types.
    """
    sender_scores = getattr(result, "sender_scores", None)
    if not isinstance(sender_scores, pd.DataFrame) or sender_scores.empty:
        raise TypeError("Expected a CellMeshResult-like object with non-empty .sender_scores")
    availability, selected_hmdb = _select_metabolite_row(
        sender_scores,
        metabolite=metabolite,
        hmdb_id=hmdb_id,
        table_name="sender_scores",
    )
    reaction_definitions = _select_reaction_definitions(
        result,
        metabolite=metabolite,
        hmdb_id=selected_hmdb,
    )

    resolved_key, resolved_layer, labels, matrix = _result_adata_context(
        result,
        adata,
        cell_type_key=cell_type_key,
        layer=layer,
    )
    sender_order, availability = _ordered_cell_types(
        availability,
        labels,
        sender_labels,
        role="sender",
    )
    cell_mask = labels.isin(sender_order).to_numpy()

    reaction_gene_names = []
    for genes in reaction_definitions["genes"]:
        reaction_gene_names.extend(str(gene) for gene in _as_reaction_sequence(genes))
    expression, expression_genes = _expression_columns(
        adata,
        matrix,
        reaction_gene_names,
        cell_mask,
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
            "metabolite": str(metabolite),
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
        title=f"{metabolite} secretion evidence in sender cell types",
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
            metabolite=metabolite,
            hmdb_id=selected_hmdb,
        )
        if values is not None:
            summary[column_name] = summary["sender"].map(values)

    plotted.update(
        {
            "plot_data": plot_data,
            "summary": summary,
            "metabolite": str(metabolite),
            "hmdb_id": selected_hmdb,
            "reaction_definitions": reaction_definitions,
            "reaction_metadata": reaction_metadata,
            "cell_type_key": resolved_key,
            "layer": resolved_layer,
        }
    )
    return plotted


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
    ``metabolite`` and ``hmdb_id`` may be supplied to restrict receiver
    metadata when the receptor participates in more than one interaction. If
    the requested metabolite name maps to multiple HMDB identifiers,
    ``hmdb_id`` is required to select exactly one metabolite context.
    Passing ``receiver_labels`` both filters and orders the displayed types.
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

    selected_scores = receiver_scores[
        receiver_scores["sensor_gene"].astype(str) == str(receptor_gene)
    ].copy()
    if metabolite is not None:
        context_columns = {"metabolite", "hmdb_id"}
        missing_context = sorted(context_columns.difference(selected_scores.columns))
        if missing_context:
            raise KeyError(
                "receiver_scores is missing required metabolite context columns: "
                + ", ".join(missing_context)
            )
        metabolite_contexts = (
            selected_scores[["metabolite", "hmdb_id"]]
            .drop_duplicates()
            .set_index(["metabolite", "hmdb_id"])
        )
        _, selected_hmdb = _select_metabolite_row(
            metabolite_contexts,
            metabolite=metabolite,
            hmdb_id=hmdb_id,
            table_name="receiver_scores",
        )
        selected_scores = selected_scores[
            (selected_scores["metabolite"].astype(str) == str(metabolite))
            & (selected_scores["hmdb_id"].astype(str) == str(selected_hmdb))
        ]
    elif hmdb_id is not None:
        if "hmdb_id" not in selected_scores:
            raise KeyError("receiver_scores does not contain an 'hmdb_id' column")
        selected_scores = selected_scores[
            selected_scores["hmdb_id"].astype(str) == str(hmdb_id)
        ]
    if selected_scores.empty:
        raise ValueError("The requested receptor was not found in receiver_scores")

    receiver_metadata = (
        selected_scores.assign(receiver=selected_scores["receiver"].astype(str))
        .groupby("receiver", sort=False, as_index=True)
        .agg(
            sensor_score=("sensor_score", "mean"),
            sensor_expr_frac=("sensor_expr_frac", "mean"),
        )
    )

    resolved_key, resolved_layer, labels, matrix = _result_adata_context(
        result,
        adata,
        cell_type_key=cell_type_key,
        layer=layer,
    )
    receiver_order, _ = _ordered_cell_types(
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
        title=f"{receptor_gene} expression in receiver cell types",
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
            "plot_data": plot_data,
            "summary": summary,
            "receptor_gene": str(receptor_gene),
            "selected_receiver_scores": selected_scores,
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
