"""Unavailable significance must not invent probabilities or hide scored events."""
import warnings

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

from cellmesh import plot_event_dotplot, run_cell_mesh


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture
def events():
    return pd.DataFrame({
        "record": range(4),
        "sender": ["A", "A", "B", "B"],
        "receiver": ["B", "B", "A", "A"],
        "hmdb_id": ["HMDB0000001", "HMDB0000002", "HMDB0000001", "HMDB0000003"],
        "sensor_gene": ["R", "Q", "R", "R"],
        "metabolite": ["Shared", "Shared", "Alias", "Third"],
        "cell_mesh_score": [0.8, 0.6, 0.7, 0.0],
        "perm_pvalue": [0.001, 0.01, 0.003, 0.05],
        "fdr_sensor_type": [0.01, 0.1, 0.03, 0.5],
        "passes_min_cells": [True] * 4,
    })


def _plot(events, supplied_ax=False, **kwargs):
    if supplied_ax:
        events = events.assign(sender="A")
        _, axis = plt.subplots(figsize=(9, 5))
        kwargs["ax"] = axis
    before = events.copy(deep=True)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = plot_event_dotplot(events, top_n=None, **kwargs)
        result["fig"].canvas.draw()
    pd.testing.assert_frame_equal(events, before)
    # The displayed table must retain exact probabilities/NA and scores.
    columns = [column for column in events if column not in ("sender", "receiver")]
    actual = result["plot_events"].set_index("record")
    expected = before.set_index("record").loc[actual.index]
    columns.remove("record")
    pd.testing.assert_frame_equal(actual[columns], expected[columns])
    return result


def _legend(result):
    legends = [axis.get_legend() for axis in result["fig"].axes if axis.get_legend() is not None]
    assert len(legends) == 1
    legend = legends[0]
    return legend.get_title().get_text(), [label.get_text() for label in legend.get_texts()]


def _assert_visible_points(result, statistic, fixed_size=140.0):
    """Check rendered positions, marker shapes, sizes and score colors."""
    plotted = result["plot_events"]
    sizes = {}
    seen = set()
    for sender, axis in zip(result["sender_order"], result["axes"]):
        for dots in axis.collections:
            offsets = dots.get_offsets()
            if not len(offsets):
                continue
            assert not np.ma.getmaskarray(offsets).any(), "A scored event was hidden"
            assert not np.ma.getmaskarray(dots.get_array()).any()
            assert np.isfinite(dots.get_sizes()).all()
            point_sizes = np.broadcast_to(dots.get_sizes(), (len(offsets),))
            for (x, y), size, score in zip(offsets, point_sizes, dots.get_array()):
                identifier, gene = result["selected_event_keys"][int(y)]
                receiver = result["receiver_order"][int(x)]
                rows = plotted.loc[
                    plotted.sender.eq(sender) & plotted.receiver.eq(receiver)
                    & plotted.hmdb_id.eq(identifier) & plotted.sensor_gene.eq(gene)
                ]
                assert len(rows) == 1
                row = rows.iloc[0]
                assert row.record not in seen
                seen.add(row.record)
                sizes[row.record] = size
                assert score == row.cell_mesh_score
                missing = statistic is None or pd.isna(row[statistic])
                marker = MarkerStyle("D" if missing else "o")
                path = marker.get_path().transformed(marker.get_transform())
                np.testing.assert_allclose(dots.get_paths()[0].vertices, path.vertices)
                if missing:
                    assert size == pytest.approx(fixed_size)
    assert seen == set(plotted.record), "Every retained scored record must be visible"
    return sizes


@pytest.mark.parametrize("supplied_ax", [False, True])
@pytest.mark.parametrize("source", ["fdr", "nullable_fdr", "pvalue", "no_columns", "fdr_missing_p_valid"])
def test_unavailable_significance_has_fixed_diamonds_and_no_numeric_legend(events, supplied_ax, source):
    if source == "no_columns":
        events = events.drop(columns=["perm_pvalue", "fdr_sensor_type"])
        statistic = None
    elif source == "pvalue":
        events = events.drop(columns="fdr_sensor_type").assign(perm_pvalue=np.nan)
        statistic = "perm_pvalue"
    else:
        events["fdr_sensor_type"] = (
            pd.Series([pd.NA] * len(events), dtype="Float64")
            if source == "nullable_fdr" else np.nan
        )
        if source != "fdr_missing_p_valid":
            events["perm_pvalue"] = np.nan
        statistic = "fdr_sensor_type"
    result = _plot(events, supplied_ax, min_dot_size=40, max_dot_size=160)
    _assert_visible_points(result, statistic, fixed_size=100)
    assert _legend(result) == ("Significance", ["Unavailable"])
    assert set(result["selected_event_keys"]) == {
        ("HMDB0000001", "R"), ("HMDB0000002", "Q"), ("HMDB0000003", "R"),
    }


@pytest.mark.parametrize("supplied_ax", [False, True])
@pytest.mark.parametrize("statistic", ["fdr_sensor_type", "perm_pvalue"])
@pytest.mark.parametrize("values", [
    [0.01, 0.2, np.nan, np.nan],  # One sender has only unavailable significance.
    [0.01, np.nan, np.nan, np.nan],  # Only one value cannot define a size range.
    [1.0, 1.0, np.nan, np.nan],  # p/FDR=1 is calculated, distinct from missing.
])
def test_partial_significance_keeps_missing_points_and_uses_only_available_sizes(events, supplied_ax, statistic, values):
    if statistic == "perm_pvalue":
        events = events.drop(columns="fdr_sensor_type")
    events[statistic] = pd.Series(values, dtype="Float64")
    result = _plot(events, supplied_ax)
    sizes = _assert_visible_points(result, statistic)
    title, labels = _legend(result)
    assert title == ("FDR" if statistic == "fdr_sensor_type" else "-log10(perm_pvalue)")
    assert labels[-1] == "Unavailable"
    assert all(np.isfinite(float(label)) for label in labels[:-1])
    if values[1] == 0.2:
        assert sizes[0] == 260
        assert sizes[1] == 20
    else:
        assert set(sizes.values()) == {140}


@pytest.mark.parametrize("statistic", ["fdr_sensor_type", "perm_pvalue"])
def test_available_significance_retains_circle_size_order_and_numeric_legend(events, statistic):
    if statistic == "perm_pvalue":
        events = events.drop(columns="fdr_sensor_type")
    result = _plot(events)
    sizes = _assert_visible_points(result, statistic)
    assert sizes[0] == 260
    assert sizes[3] == 20
    assert sizes[0] > sizes[2] > sizes[1] > sizes[3]
    title, labels = _legend(result)
    assert title == ("FDR" if statistic == "fdr_sensor_type" else "-log10(perm_pvalue)")
    assert len(labels) == 3
    assert all(np.isfinite(float(label)) for label in labels)


@pytest.mark.parametrize("statistic", ["fdr_sensor_type", "perm_pvalue"])
def test_zero_and_one_are_available_and_are_never_replaced_with_na(events, statistic):
    if statistic == "perm_pvalue":
        events = events.drop(columns="fdr_sensor_type")
    events[statistic] = [0.0, 1.0, np.nan, 0.1]
    result = _plot(events)
    _assert_visible_points(result, statistic)
    assert result["plot_events"][statistic].isna().sum() == 1


@pytest.mark.parametrize("statistic", ["fdr_sensor_type", "perm_pvalue"])
@pytest.mark.parametrize("invalid", [np.inf, -np.inf, -0.1, 1.1])
def test_invalid_probabilities_are_reported_instead_of_being_treated_as_missing(events, statistic, invalid):
    if statistic == "perm_pvalue":
        events = events.drop(columns="fdr_sensor_type")
    events.loc[0, statistic] = invalid
    with pytest.raises(ValueError, match=r"finite probabilities in \[0, 1\]"):
        _plot(events)


def test_explicit_filters_keep_excluding_unavailable_significance(events):
    events.loc[2:, "fdr_sensor_type"] = np.nan
    result = _plot(events, max_fdr=0.05)
    _assert_visible_points(result, "fdr_sensor_type")
    assert result["plot_events"].record.tolist() == [0]
    assert "Unavailable" not in _legend(result)[1]
    with pytest.raises(ValueError, match="No events remain after filtering"):
        _plot(events.assign(fdr_sensor_type=np.nan), max_fdr=0.05)
    with pytest.raises(KeyError, match="fdr_sensor_type"):
        _plot(events.drop(columns="fdr_sensor_type"), max_fdr=0.05)


def test_missingness_is_determined_after_selection_with_custom_probability_columns(events):
    events = events.rename(columns={"fdr_sensor_type": "adjusted_p", "perm_pvalue": "raw_p"})
    events.loc[events.sender.eq("B"), "adjusted_p"] = np.nan
    result = _plot(events, sender_labels=["B"], fdr_col="adjusted_p", pvalue_col="raw_p")
    _assert_visible_points(result, "adjusted_p")
    assert _legend(result) == ("Significance", ["Unavailable"])
    assert len(result["plot_events"]) == 2


@pytest.mark.parametrize("supplied_ax", [False, True])
@pytest.mark.parametrize("all_missing", [False, True])
def test_unavailable_legend_is_visible_and_separate_from_colorbar(events, supplied_ax, all_missing):
    events["metabolite"] = ["First metabolite", "Second metabolite", "First alias", "Third metabolite"]
    events["fdr_sensor_type"] = [np.nan] * 4 if all_missing else [0.01, 0.2, np.nan, np.nan]
    result = _plot(events, supplied_ax)
    fig = result["fig"]
    renderer = fig.canvas.get_renderer()
    legend = next(axis.get_legend() for axis in fig.axes if axis.get_legend() is not None)
    legend_box = legend.get_window_extent(renderer)
    assert legend_box.x0 >= fig.bbox.x0
    assert legend_box.x1 <= fig.bbox.x1
    assert legend_box.y0 >= fig.bbox.y0
    assert legend_box.y1 <= fig.bbox.y1
    colorbar = next(axis for axis in fig.axes if axis.get_title() == "Cellmesh_score")
    assert not legend_box.overlaps(colorbar.get_tightbbox(renderer))
    assert all(not legend_box.overlaps(axis.get_window_extent(renderer)) for axis in result["axes"])


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("n_perms", [0, 3])
def test_actual_scoring_results_render_with_and_without_permutations(sample_mode, n_perms):
    from anndata import AnnData

    data = AnnData(
        np.array([[8., 1.], [0., 2.], [6., 1.], [0., 2.]]),
        obs=pd.DataFrame({"cell_type": ["A", "B", "A", "B"], "sample": ["S1", "S1", "S2", "S2"]},
                         index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["PROD", "R"]),
    )
    enzyme = pd.DataFrame([dict(metabolite="Enzyme name", hmdb_id="HMDB0000001", gene="PROD",
                               role="production", reaction="P")])
    sensor = pd.DataFrame([dict(metabolite="Interaction alias", hmdb_id="HMDB0000001",
                               sensor_gene="R", sensor_type="Transporter")])
    result = run_cell_mesh(data, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                           min_cells=1, n_perms=n_perms)
    before = result.events.copy(deep=True)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        plotted = plot_event_dotplot(result, top_n=None)
        plotted["fig"].canvas.draw()
    pd.testing.assert_frame_equal(result.events, before)
    assert len(plotted["plot_events"]) == len(before) == 4
    assert (plotted["plot_events"].cell_mesh_score == 0).any()
    plotted["plot_events"]["record"] = range(len(plotted["plot_events"]))
    _assert_visible_points(plotted, "fdr_sensor_type")
    if n_perms == 0:
        assert _legend(plotted) == ("Significance", ["Unavailable"])
        assert plotted["plot_events"][["perm_pvalue", "fdr_sensor_type", "fdr_global"]].isna().all().all()
    else:
        assert _legend(plotted)[0] == "FDR"
        assert all(np.isfinite(float(label)) for label in _legend(plotted)[1])
