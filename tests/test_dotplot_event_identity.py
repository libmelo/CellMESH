"""Dotplot row identity must retain the HMDB ID through every plotting step."""
import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cellmesh import plot_event_dotplot, run_cell_mesh


KEY1 = ("Shared", "HMDB0000001", "R")
KEY2 = ("Shared", "HMDB0000002", "R")
KEY3 = ("Other", "HMDB0000003", "Q")


@pytest.fixture
def same_name_events():
    # The best HMDB differs between sender panels. Merging by the display
    # label would make a single plotted row silently switch metabolites.
    rows = []
    for sender, key, score, fdr in (
        ("A", KEY1, 0.8, 0.01),
        ("A", KEY2, 0.6, 0.03),
        ("B", KEY1, 0.4, 0.04),
        ("B", KEY2, 0.9, 0.005),
        ("A", KEY3, 0.7, 0.02),
    ):
        rows.append(dict(zip(("metabolite", "hmdb_id", "sensor_gene"), key),
                         sender=sender, receiver="C", sensor_type="Transporter",
                         cell_mesh_score=score, perm_pvalue=fdr / 2,
                         fdr_sensor_type=fdr, passes_min_cells=True))
    return pd.DataFrame(rows)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def _identities(keys):
    return [tuple(key[-2:]) for key in keys]


def _keys(frame):
    return set(frame[["metabolite", "hmdb_id", "sensor_gene"]].itertuples(index=False, name=None))


def test_dotplot_preserves_distinct_hmdb_rows_and_panel_coordinates(same_name_events):
    original = same_name_events.copy(deep=True)
    result = plot_event_dotplot(same_name_events, top_n=None)
    plotted = result["plot_events"]
    assert len(plotted) == 5
    assert result["selected_event_keys"] == _identities([KEY2, KEY1, KEY3])
    assert result["selected_events"] == [
        "Shared (HMDB0000002) -> R", "Shared (HMDB0000001) -> R", "Other (HMDB0000003) -> Q",
    ]
    assert plotted.loc[plotted.hmdb_id.eq(KEY1[1]), "_event_label"].unique().tolist() == [
        "Shared (HMDB0000001) -> R",
    ]
    for sender, axis in zip(result["sender_order"], result["axes"]):
        panel = plotted.loc[plotted.sender.eq(sender)]
        dots = axis.collections[0]
        expected_y = [result["selected_event_keys"].index(key[1:])
                      for key in panel[["metabolite", "hmdb_id", "sensor_gene"]].itertuples(index=False, name=None)]
        np.testing.assert_array_equal(dots.get_offsets()[:, 1], expected_y)
        np.testing.assert_array_equal(dots.get_array(), panel.cell_mesh_score)
    assert [tick.get_text() for tick in result["ax"].get_yticklabels()] == result["selected_events"]
    pd.testing.assert_frame_equal(same_name_events, original)


@pytest.mark.parametrize("top_n, expected", [(1, [KEY2]), (2, [KEY2, KEY1]), (None, [KEY2, KEY1, KEY3])])
@pytest.mark.parametrize("with_fdr", [True, False])
def test_top_n_ranks_full_event_identities(same_name_events, top_n, expected, with_fdr):
    events = same_name_events if with_fdr else same_name_events.drop(columns="fdr_sensor_type")
    result = plot_event_dotplot(events, top_n=top_n)
    assert result["selected_event_keys"] == _identities(expected)
    assert _keys(result["plot_events"]) == set(expected)


@pytest.mark.parametrize("selection, expected", [
    ([KEY2, KEY1], [KEY2, KEY1]),
    ([dict(zip(("metabolite", "hmdb_id", "sensor_gene"), KEY2))], [KEY2]),
    ([KEY2[1:], KEY1[1:]], [KEY2, KEY1]),
    ([KEY1[1:], KEY2[1:]], [KEY1, KEY2]),
    ([KEY2, KEY3[1:], KEY1, KEY2], [KEY2, KEY3, KEY1]),
])
def test_explicit_selectors_keep_every_requested_identity_in_order(same_name_events, selection, expected):
    result = plot_event_dotplot(same_name_events, event_keys=selection, top_n=1)
    assert result["selected_event_keys"] == _identities(expected)
    assert _keys(result["plot_events"]) == set(expected)
    assert len(result["selected_events"]) == len(expected)


def test_generated_identity_keys_can_be_selected_again(same_name_events):
    first = plot_event_dotplot(same_name_events, top_n=None)
    second = plot_event_dotplot(same_name_events, event_keys=first["selected_event_keys"])
    assert second["selected_event_keys"] == first["selected_event_keys"]
    assert _keys(second["plot_events"]) == _keys(first["plot_events"])


def test_only_true_duplicate_pair_events_are_deduplicated(same_name_events):
    duplicate = same_name_events.iloc[[0]].assign(cell_mesh_score=0.99, fdr_sensor_type=0.2)
    higher_score_tie = same_name_events.iloc[[0]].assign(cell_mesh_score=0.85)
    events = pd.concat([same_name_events, duplicate, higher_score_tie], ignore_index=True)
    result = plot_event_dotplot(events, top_n=None)
    plotted = result["plot_events"]
    assert len(plotted) == 5
    selected = plotted.loc[plotted.sender.eq("A") & plotted.hmdb_id.eq(KEY1[1])].iloc[0]
    assert selected.cell_mesh_score == 0.85
    assert selected.fdr_sensor_type == 0.01
    assert _keys(plotted) == {KEY1, KEY2, KEY3}


def test_thresholds_and_qc_do_not_mix_same_name_events(same_name_events):
    events = same_name_events.copy()
    events.loc[events.hmdb_id.eq(KEY2[1]), "passes_min_cells"] = False
    default = plot_event_dotplot(events, event_keys=[KEY1[1:], KEY2[1:]], max_fdr=0.025)
    assert _keys(default["plot_events"]) == {KEY1}
    assert default["selected_events"] == ["Shared (HMDB0000001) -> R"]
    all_events = plot_event_dotplot(events, event_keys=[KEY1[1:], KEY2[1:]], max_fdr=0.025, qc_only=False)
    assert _keys(all_events["plot_events"]) == {KEY1, KEY2}
    assert len(all_events["plot_events"]) == 2


def test_existing_axis_preserves_both_hmdb_rows(same_name_events):
    fig, axis = plt.subplots()
    result = plot_event_dotplot(same_name_events, event_keys=[KEY1, KEY2], sender_labels=["A"], ax=axis)
    assert result["ax"] is axis
    assert result["selected_event_keys"] == _identities([KEY1, KEY2])
    assert len(axis.collections[0].get_offsets()) == 2


def test_tuple_selectors_do_not_collapse_delimiter_containing_fields(same_name_events):
    # Both tuples serialize to the same legacy " | " string. Internal
    # identity and an explicit tuple selection must still distinguish them.
    first = ("M | HMDB0000001", "HMDB0000002", "R")
    second = ("M", "HMDB0000001", "HMDB0000002 | R")
    events = same_name_events.iloc[:2].copy()
    events.loc[:, ["metabolite", "hmdb_id", "sensor_gene"]] = [first, second]
    both = plot_event_dotplot(events, top_n=None)
    assert _keys(both["plot_events"]) == {first, second}
    only_second = plot_event_dotplot(events, event_keys=[second])
    assert _keys(only_second["plot_events"]) == {second}


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_dotplot_keeps_actual_scoring_results_for_both_hmdb_ids(sample_mode):
    from anndata import AnnData

    adata = AnnData(
        np.array([[8., 1., 1.], [0., 4., 2.], [6., 2., 1.], [0., 5., 2.]]),
        obs=pd.DataFrame({"cell_type": ["A", "B", "A", "B"], "sample": ["S1", "S1", "S2", "S2"]},
                         index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["P1", "P2", "R"]),
    )
    enzyme = pd.DataFrame([
        dict(metabolite="Shared", hmdb_id=f"HMDB{i:07d}", gene=f"P{i}", role="production", reaction=f"RX{i}")
        for i in (1, 2)
    ])
    sensor = pd.DataFrame([
        dict(metabolite="Shared", hmdb_id=f"HMDB{i:07d}", sensor_gene="R", sensor_type="Transporter")
        for i in (1, 2)
    ])
    result = run_cell_mesh(adata, enzyme, sensor, min_cells=1, n_perms=9,
                           sample_mode=sample_mode, sample_key="sample")
    original = result.events.copy(deep=True)
    plotted = plot_event_dotplot(result, sender_labels=["A"], receiver_labels=["B"], event_keys=[KEY1, KEY2])
    assert plotted["selected_event_keys"] == _identities([KEY1, KEY2])
    actual = plotted["plot_events"].set_index("hmdb_id")
    expected = original.loc[original.sender.eq("A") & original.receiver.eq("B")].set_index("hmdb_id")
    columns = ["cell_mesh_score", "perm_pvalue", "fdr_sensor_type", "fdr_global"]
    pd.testing.assert_frame_equal(actual[columns].sort_index(), expected[columns].sort_index())
    pd.testing.assert_frame_equal(result.events, original)
