"""One immutable selection scope across CSV, reports and aggregate figures."""
import copy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from cellmesh import (CellMeshResult, EventSelection, plot_communication_network,
                      plot_event_dotplot, plot_significant_event_counts)
from cellmesh.selection import EVENT_KEYS, report_selection
from test_reporting_and_persistence import make_result, assert_nested_equal


@pytest.fixture
def result():
    rows = [dict(sender=s, receiver=r, hmdb_id=f"HMDB{i:07d}", sensor_gene=f"R{i}",
                 metabolite="same display name", sensor_type="Receptor", cell_mesh_score=.5,
                 fdr_sensor_type=.8, fdr_global=.9, perm_pvalue=.7, passes_min_cells=False)
            for i in range(1, 13) for s, r in [("A", "B"), ("B", "A")]]
    return CellMeshResult(pd.DataFrame(rows), pd.DataFrame(), pd.DataFrame(), {})


@pytest.fixture(scope="module")
def captured():
    return make_result(n_perms=0)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def keys(frame):
    return set(frame[EVENT_KEYS].itertuples(index=False, name=None))


def test_and_or_direction_and_normalized_identity(result):
    focus = result.filter_events(source="all", sender_labels=[" A ", "B"], receiver_labels=["A", "B"],
                                 cell_pairs=[(" A ", "B")], hmdb_ids=[" hmdb0000001 ", "HMDB0000002"],
                                 sensor_genes=["R1", "R2"])
    assert isinstance(focus, EventSelection)
    assert keys(focus.events) == {("A", "B", "HMDB0000001", "R1"),
                                 ("A", "B", "HMDB0000002", "R2")}
    # Directed pair and endpoint constraints intersect; neither overrides the other.
    assert result.filter_events(source="all", sender_labels=["B"], cell_pairs=[("A", "B")]).events.empty
    assert len(result.filter_events(source="all", cell_pairs=[("B", "A")]).events) == 12
    assert result.filter_events(source="all", cell_pairs=[("A", "A")]).events.empty
    assert len(result.filter_events(source="all", hmdb_ids=["HMDB0000001"]).events) == 2  # names never merge IDs


@pytest.mark.parametrize("name", ["sender_labels", "receiver_labels", "cell_pairs", "hmdb_ids", "sensor_genes"])
def test_empty_lists_match_none_and_none_unrestricted(result, name):
    assert result.filter_events(source="all", **{name: []}).events.empty
    assert len(result.filter_events(source="all", **{name: None}).events) == 24


@pytest.mark.parametrize("kwargs", [dict(sender_labels=["unknown"]), dict(receiver_labels=["unknown"]),
    dict(hmdb_ids=["HMDB9999999"]), dict(sensor_genes=["unknown"]), dict(cell_pairs=[("A", "unknown")]),
    dict(cell_pairs=[("A",)]), dict(cell_pairs=["AB"]), dict(cell_pairs=[("A", "B"), (" A ", "B")]),
    dict(hmdb_ids=["hmdb0000001", "HMDB0000001"]), dict(sender_labels=["A", " A "]),
    dict(hmdb_ids=["same display name"]), dict(sender_labels="A"), dict(cell_pairs="AB")])
def test_bad_identifiers_rejected(result, kwargs):
    with pytest.raises((ValueError, TypeError)):
        result.filter_events(source="all", **kwargs)


@pytest.mark.parametrize("kwargs", [dict(max_fdr=-.1), dict(max_fdr=1.1), dict(max_fdr=True),
    dict(max_fdr=np.inf), dict(min_cell_mesh_score=np.nan), dict(min_cell_mesh_score=1.1),
    dict(min_cell_mesh_score=True), dict(qc_only="true"), dict(fdr_col="perm_pvalue")])
def test_invalid_options_rejected(result, kwargs):
    with pytest.raises((ValueError, TypeError)):
        result.filter_events(source="all", **kwargs)


def test_prior_known_unscored_identifiers_are_valid_empty(captured):
    assert captured.filter_events(source="all", hmdb_ids=["HMDB0000002"]).events.empty
    assert captured.filter_events(source="all", sensor_genes=["ABSENT_R"]).events.empty


def test_thresholds_are_inclusive_na_is_conditional_and_counted(result):
    result.events.loc[:3, "fdr_sensor_type"] = [.05, .05, np.nan, .2]
    result.events.loc[:3, "cell_mesh_score"] = [.5, .49, .7, np.nan]
    result.events.loc[0, "passes_min_cells"] = True
    assert len(result.filter_events(source="all", ).events) == 24
    focus = result.filter_events(source="all", max_fdr=.05, min_cell_mesh_score=.5, qc_only=True)
    assert len(focus.events) == 1
    steps = {row["filter"]: row for row in focus.summary["steps"]}
    assert steps["max_fdr"]["excluded_missing"] == 1
    assert steps["min_cell_mesh_score"]["excluded"] == 1
    assert sum(row["excluded"] for row in steps.values()) == 23
    score_only = result.filter_events(source="all", min_cell_mesh_score=.5)
    assert score_only.events.fdr_sensor_type.isna().sum() == 1
    assert score_only.summary["missing_score"] == 1
    assert len(result.filter_events(source="all", max_fdr=.9, fdr_col="fdr_global").events) == 24


@pytest.mark.parametrize("missing_column", [False, True])
def test_unavailable_fdr_is_not_interpreted_as_significant(result, missing_column):
    if missing_column:
        result.events = result.events.drop(columns="fdr_sensor_type")
    else:
        result.events["fdr_sensor_type"] = np.nan
    assert len(result.filter_events(source="all", ).events) == 24
    with pytest.raises(ValueError, match="calculate FDR or omit"):
        result.filter_events(source="all", max_fdr=.05)


def test_stable_sorting_and_no_scientific_mutation(result):
    result.events["fdr_sensor_type"] = np.nan
    result.events.loc[3:5, "fdr_sensor_type"] = .01
    result.events.loc[5, "cell_mesh_score"] = .8
    result.events = result.events.sample(frac=1, random_state=8)
    result.events.attrs["null_scores"] = np.arange(48).reshape(24, 2)
    before = copy.deepcopy(result)
    expected = result.events.sort_values(["fdr_sensor_type", "cell_mesh_score"] + EVENT_KEYS,
                       ascending=[True, False] + [True] * 4, na_position="last").reset_index(drop=True)
    focus = result.filter_events(source="all", )
    pd.testing.assert_frame_equal(focus.events, expected, check_exact=True)
    assert focus.events.attrs == {}
    assert_nested_equal(result.events, before.events)
    # Neither public result tables nor metadata can change selection membership.
    frame = focus.events
    frame.loc[:, "cell_mesh_score"] = 0
    meta = focus.criteria
    meta["sender_labels"] = ["B"]
    summary = focus.summary
    summary["steps"][0]["remaining"] = 0
    pd.testing.assert_frame_equal(focus.events, expected, check_exact=True)
    assert focus.criteria["sender_labels"] is None
    assert focus.summary["steps"][0]["remaining"] == 24
    result.events["fdr_sensor_type"] = np.nan
    score_sort = result.filter_events(source="all", )
    assert score_sort.summary["order_by"][0] == "cell_mesh_score"
    assert score_sort.events.iloc[0].cell_mesh_score == .8


@pytest.mark.parametrize("mutation", ["score", "fdr", "params", "drop"])
def test_source_changes_require_reselection(result, mutation, tmp_path):
    focus = result.filter_events(source="all", )
    if mutation == "params":
        result.parameters["min_cells"] = 99
    elif mutation == "drop":
        result.events = result.events.iloc[:-1]
    else:
        result.events.loc[0, "cell_mesh_score" if mutation == "score" else "fdr_global"] = .1
    with pytest.raises(ValueError, match="Source result changed"):
        focus.to_csv(tmp_path / "selected.csv")
    assert not list(tmp_path.iterdir())


def test_cannot_attach_selection_to_another_result(result, tmp_path):
    other = copy.deepcopy(result)
    with pytest.raises(ValueError, match="another result"):
        other.save(tmp_path / "wrong", report=True, selection=result.filter_events(source="all", ))
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("empty", [False, True])
def test_filtered_csv_exact_scope_single_file_with_headers(result, tmp_path, empty):
    focus = result.filter_events(source="all", sender_labels=[] if empty else ["A"])
    path = focus.to_csv(tmp_path / "selected.csv")
    assert list(tmp_path.iterdir()) == [path]
    table = pd.read_csv(path, float_precision="round_trip")
    assert list(table.columns) == list(result.events.columns) + ["positive_status"]
    if not empty:
        pd.testing.assert_frame_equal(table.drop(columns="positive_status"), focus.events, check_exact=True)
    else:
        assert table.empty


def test_csv_atomic_failure_and_destination_checks(result, tmp_path, monkeypatch):
    path = tmp_path / "selected.csv"
    path.write_text("previous")
    def fail(*args, **kwargs):
        raise OSError("simulated failure")
    with monkeypatch.context() as patch:
        patch.setattr("cellmesh.selection.os.replace", fail)
        with pytest.raises(OSError, match="simulated failure"):
            result.filter_events(source="all", ).to_csv(path)
    assert path.read_text() == "previous"
    assert list(tmp_path.iterdir()) == [path]
    for destination in [tmp_path, tmp_path / "link.csv"]:
        if destination != tmp_path:
            destination.symlink_to(path)
        with pytest.raises(ValueError, match="regular file"):
            result.filter_events(source="all", ).to_csv(destination)
    files = result.save(tmp_path / "analysis", report=False)
    for destination in files.values():
        with pytest.raises(ValueError, match="standard result delivery"):
            result.filter_events(source="all", ).to_csv(destination)


@pytest.mark.parametrize("detail_n,expected", [(0, 0), (1, 0), (10, 0), (None, 0)])
def test_report_explanation_limit_full_persistence_and_focused_table(captured, tmp_path, detail_n, expected):
    result = copy.deepcopy(captured)
    before = copy.deepcopy(result)
    focus = result.filter_events(source="all", sender_labels=["A"])
    files = result.save(tmp_path / "analysis", report=True, selection=focus, detail_n=detail_n, include_all_events=True)
    html = files["filtered_report"].read_text()
    assert html.count('class="event-explanation"') == expected
    assert "原始事件：4；入选事件：2；排除事件：2" in html
    assert "全量输入覆盖度（不随筛选改变）" not in html
    assert "全量输入覆盖度（不随筛选改变）" in files["report"].read_text()
    if expected:
        assert '<details data-table="gene_details"><summary>基因明细</summary>' in html
    assert "<script>bad()" not in html
    pd.testing.assert_frame_equal(pd.read_csv(files["events_csv"]), result.events, check_dtype=False)
    restored = CellMeshResult.load(files["result"])
    pd.testing.assert_frame_equal(restored.events, result.events, check_exact=True)
    pd.testing.assert_frame_equal(restored.filter_events(source="all", sender_labels=["A"]).events, focus.events)
    assert_nested_equal(before.events, result.events)
    focus.to_csv(tmp_path / "selected.csv")
    assert len(pd.read_csv(tmp_path / "selected.csv")) == 2
    focus.to_html(tmp_path / "focus.html", detail_n=detail_n)
    assert (tmp_path / "focus.html").read_text().count('class="event-explanation"') == expected


def test_empty_selection_report_and_legacy_manual_details(captured, tmp_path):
    focus = captured.filter_events(source="all", hmdb_ids=[])
    path = focus.to_html(tmp_path / "empty.html")
    html = path.read_text()
    assert "入选事件：0" in html
    assert "没有事件满足全部筛选条件" in html
    assert "<h2>事件解释：" not in html
    selectors = captured.events[EVENT_KEYS].to_dict("records")
    # Old four-field exports must not silently expand to pair-wide reports.
    with pytest.raises(TypeError, match="retired"):
        captured.to_html(tmp_path / "manual.html", events=selectors, detail_n=0)
    path = captured.to_html(tmp_path / "plain.html")
    assert "<h2>事件解释：" not in path.read_text()


def test_default_detail_limit_only_limits_explanations(result):
    focus = result.filter_events(source="all", )
    table, selectors = report_selection(result, focus, None, 10)
    assert len(table) == 0 and len(selectors) == 0  # main table uses positive calls
    assert len(report_selection(result, focus, None, None)[1]) == 0


@pytest.mark.parametrize("kwargs", [dict(events=[{}]), dict(detail_n=-1), dict(detail_n=True), dict(detail_n=1.5)])
def test_invalid_report_requests_fail_before_writing(captured, tmp_path, kwargs):
    with pytest.raises((ValueError, TypeError)):
        captured.save(tmp_path / "invalid", report=True, selection=captured.filter_events(source="all", ), **kwargs)
    assert not list(tmp_path.iterdir())


def test_selection_without_report_writes_csv(captured, tmp_path):
    paths = captured.save(tmp_path / "analysis", report=False, selection=captured.filter_events(source="all"))
    assert "filtered_csv" in paths and "report" not in paths


@pytest.mark.parametrize("plot", [plot_event_dotplot, plot_communication_network, plot_significant_event_counts])
def test_selected_plot_has_no_implicit_fdr_qc_or_top_limit(result, plot):
    focus = result.filter_events(source="all")
    output = plot(focus)
    plotted = output.get("plot_events", output.get("filtered_events"))
    assert len(plotted) == len(focus.events) == 24
    assert set(map(tuple, plotted[EVENT_KEYS].to_numpy())) == set(map(tuple, focus.events[EVENT_KEYS].to_numpy()))
    assert output['selection']['n_plotted'] == 24


def test_explicit_dotplot_top_n_limits_metabolite_sensor_pairs_only(result, tmp_path):
    focus = result.filter_events(source="all", )
    output = plot_event_dotplot(focus, top_n=3)
    assert len(output["selected_event_keys"]) == 3
    assert len(output["plot_events"]) == 6  # two directed pairs per metabolite/sensor
    assert output["selection"]["n_selected"] == 24
    assert "top 3 metabolite-sensor pairs" in output["title_artist"].get_text()
    assert len(pd.read_csv(focus.to_csv(tmp_path / "selected.csv"))) == 24


@pytest.mark.parametrize("plot", [plot_event_dotplot, plot_communication_network, plot_significant_event_counts])
def test_selection_plot_missing_scores_counted_without_modifying_selection(result, plot):
    result.events.loc[0, "cell_mesh_score"] = np.nan
    result.events["fdr_sensor_type"] = np.nan
    focus = result.filter_events(source="all", )
    output = plot(focus)
    assert output["selection"]["n_selected"] == 24
    assert output["selection"]["n_plotted"] == 23
    assert output["selection"]["n_missing_score"] == 1
    assert len(focus.events) == 24


@pytest.mark.parametrize("kwargs", [dict(max_fdr=.05), dict(min_cell_mesh_score=.2),
    dict(qc_only=True), dict(sender_labels=["A"]), dict(event_keys=[("HMDB1", "R1")])])
def test_plot_cannot_silently_rescope_selection(result, kwargs):
    with pytest.raises(ValueError, match="scope|filter_events|select_events"):
        plot_event_dotplot(result.filter_events(source="all", ), **kwargs)


@pytest.mark.parametrize("plot", [plot_event_dotplot, plot_communication_network, plot_significant_event_counts])
def test_empty_selection_plot_clear_error(result, plot):
    with pytest.raises(ValueError, match="selection is empty"):
        plot(result.filter_events(source="all", hmdb_ids=[]))


def test_legacy_plot_defaults_remain(result):
    # Old direct-result API still applies QC/FDR and a 10-pair dotplot default.
    result.events["passes_min_cells"] = True
    result.events["fdr_sensor_type"] = .01
    result.events.loc[0, "fdr_sensor_type"] = .9
    network = plot_communication_network(result)
    assert len(network["filtered_events"]) == 23
    assert "selection" not in network
    dot = plot_event_dotplot(result)
    assert len(dot["selected_event_keys"]) == 10


def test_empty_source_has_valid_csv_and_summary(result, tmp_path):
    result.events = result.events.iloc[:0]
    focus = result.filter_events(source="all", )
    assert focus.events.empty
    assert focus.summary["empty_reason"] == "No source events"
    assert list(pd.read_csv(focus.to_csv(tmp_path / "empty.csv")).columns) == list(result.events.columns) + ["positive_status"]


@pytest.mark.parametrize("column,value", [("cell_mesh_score", "bad"), ("cell_mesh_score", np.inf),
    ("cell_mesh_score", True), ("fdr_sensor_type", -.1), ("passes_min_cells", "bad")])
def test_malformed_measurements_are_not_silently_missing(result, column, value):
    result.events[column] = result.events[column].astype(object)
    result.events.loc[0, column] = value
    with pytest.raises(ValueError, match=column):
        result.filter_events(source="all", )


def test_qc_text_is_normalized_for_selection_but_original_values_exported(result):
    result.events["passes_min_cells"] = pd.Series(["True", "False", pd.NA] * 8, dtype="string")
    focus = result.filter_events(source="all", qc_only=True)
    assert len(focus.events) == 8
    assert focus.events.passes_min_cells.eq("True").all()
    assert focus.summary["steps"][-1]["excluded_missing"] == 8
    result.events = result.events.drop(columns="passes_min_cells")
    assert len(result.filter_events(source="all", ).events) == 24
    with pytest.raises(ValueError, match="requires passes_min_cells"):
        result.filter_events(source="all", qc_only=True)


def test_numeric_string_values_and_float16_sort_without_rewriting_source(result):
    result.events["cell_mesh_score"] = result.events.cell_mesh_score.astype("float16")
    result.events["fdr_sensor_type"] = result.events.fdr_sensor_type.astype(str)
    focus = result.filter_events(source="all", max_fdr=.8, min_cell_mesh_score=.5)
    assert len(focus.events) == 24
    assert focus.events.cell_mesh_score.dtype == np.dtype("float16")
    assert focus.events.fdr_sensor_type.eq("0.8").all()


def test_sample_rows_and_duplicate_aggregate_keys_rejected(result):
    result.events["sample"] = "S1"
    with pytest.raises(ValueError, match="aggregate events"):
        result.filter_events(source="all", )
    result.events = result.events.drop(columns="sample")
    result.events = pd.concat([result.events, result.events.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="identities must be unique"):
        result.filter_events(source="all", )


def test_filtered_report_preview_does_not_change_counts_or_csv(captured, tmp_path):
    # The event browser must retain all rows beyond the former 100-row preview.
    result = copy.deepcopy(captured)
    row = pd.DataFrame(result.events.iloc[[0]])
    result.events = pd.concat([row.assign(sensor_gene=f"R{i}") for i in range(120)], ignore_index=True)
    focus = result.filter_events(source="all", )
    html = focus.to_html(tmp_path / "preview.html", detail_n=0).read_text()
    assert "入选事件：120" in html
    assert "共 120 条事件；静态预览前 50 条" in html
    assert '<span class="metric-label">所选事件</span><strong>120</strong>' in html
    assert len(pd.read_csv(focus.to_csv(tmp_path / "selected.csv"))) == 120


def test_focused_plot_rejects_invalid_qc_option(result):
    with pytest.raises(TypeError, match="qc_only must be boolean"):
        plot_event_dotplot(result.filter_events(source="all", ), qc_only=None)
