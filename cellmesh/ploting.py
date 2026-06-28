"""
Plotting utilities for CELL MESH.

This module is intentionally lightweight for now. Plotting functions should
import optional visualization dependencies, such as matplotlib, inside the
function body so the core package import remains fast.
"""
from __future__ import annotations

from typing import Any, Optional

import pandas as pd


def _events_frame(result_or_events: Any) -> pd.DataFrame:
    if isinstance(result_or_events, pd.DataFrame):
        return result_or_events
    events = getattr(result_or_events, "events", None)
    if isinstance(events, pd.DataFrame):
        return events
    raise TypeError("Expected a CellMeshResult-like object with .events or a pandas DataFrame")


def plot_significant_event_counts(
    result_or_events: Any,
    *,
    min_cell_mesh_score: Optional[float] = None,
    max_perm_pvalue: Optional[float] = None,
    max_fdr: Optional[float] = 0.05,
    unique_keys: Optional[list[str]] = None,
    cmap: str = "Blues",
    ax: Any = None,
    show_values: bool = True,
    top_n_summary: int = 20,
) -> dict[str, Any]:
    """
    Plot thresholded sender-to-receiver metabolite-sensor event counts.

    Counts are based on unique ``sender + receiver + metabolite + hmdb_id +
    sensor_gene`` combinations after optional score, permutation p-value, and
    FDR filters. By default, the plot summarizes events with ``fdr <= 0.05``.
    """
    events = _events_frame(result_or_events).copy()
    if events.empty:
        raise ValueError("No events available to plot")

    required = {"sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "cell_mesh_score"}
    if max_perm_pvalue is not None:
        required.add("perm_pvalue")
    if max_fdr is not None:
        required.add("fdr")
    missing = sorted(required.difference(events.columns))
    if missing:
        raise KeyError(f"Missing required event columns: {', '.join(missing)}")

    filtered = events
    if min_cell_mesh_score is not None:
        filtered = filtered[filtered["cell_mesh_score"] >= min_cell_mesh_score]
    if max_perm_pvalue is not None:
        filtered = filtered[filtered["perm_pvalue"].notna()]
        filtered = filtered[filtered["perm_pvalue"] <= max_perm_pvalue]
    if max_fdr is not None:
        filtered = filtered[filtered["fdr"].notna()]
        filtered = filtered[filtered["fdr"] <= max_fdr]

    all_senders = pd.Index(sorted(events["sender"].unique()), name="sender")
    all_receivers = pd.Index(sorted(events["receiver"].unique()), name="receiver")

    if unique_keys is None:
        unique_keys = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene"]
    missing_unique = sorted(set(unique_keys).difference(events.columns))
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
        "summary": summary.head(top_n_summary) if top_n_summary is not None else summary,
        "all_summary": summary,
        "filtered_events": filtered.copy(),
        "thresholds": {
            "min_cell_mesh_score": min_cell_mesh_score,
            "max_perm_pvalue": max_perm_pvalue,
            "max_fdr": max_fdr,
            "events_before_filter": len(events),
            "events_after_filter": len(filtered),
            "unique_events_after_filter": len(unique_pair_events),
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
    fdr_col: str = "fdr",
    pvalue_col: str = "perm_pvalue",
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
    """
    events = _events_frame(result_or_events).copy()
    if events.empty:
        raise ValueError("No events available to plot")

    required = {
        "sender",
        "receiver",
        "metabolite",
        "hmdb_id",
        "sensor_gene",
        score_col,
    }
    if max_perm_pvalue is not None or pvalue_col in events.columns:
        required.add(pvalue_col)
    if max_fdr is not None or fdr_col in events.columns:
        required.add(fdr_col)
    missing = sorted(required.difference(events.columns))
    if missing:
        raise KeyError(f"Missing required event columns: {', '.join(missing)}")

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
        size_label = "FDR" if fdr_col.lower() == "fdr" else pvalue_col
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
        axes = [ax]
        fig = ax.figure
        cax = None
        lax = None
        sender_order = sender_order[:1]

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
        "selected_events": label_order,
        "sender_order": sender_order,
        "receiver_order": receiver_order,
        "thresholds": {
            "top_n": top_n,
            "event_keys": event_keys,
            "min_cell_mesh_score": min_cell_mesh_score,
            "max_perm_pvalue": max_perm_pvalue,
            "max_fdr": max_fdr,
            "events_before_filter": len(events),
            "events_after_filter": len(filtered),
            "events_plotted": len(plot_events),
        },
    }


__all__ = ["plot_significant_event_counts", "plot_event_dotplot"]
