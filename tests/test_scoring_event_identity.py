"""Names are annotations throughout scoring, aggregation, export and plotting."""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import compute_metabolite_availability, run_cell_mesh
from cellmesh.core import _make_cell_mesh_events
from cellmesh.score import compute_sensor_scores


H1, H2, H0 = "HMDB0000001", "HMDB0000002", "HMDB0000003"
KEY = ["sender", "receiver", "hmdb_id", "sensor_gene"]


@pytest.fixture
def identity_case():
    # S3 lacks B; ZERO is measured and zero everywhere. Neither state is missing
    # display metadata, and changing names must preserve this distinction.
    data = AnnData(
        np.array([[8, 2, 0, 1, 5], [6, 3, 0, 2, 3],
                  [1, 5, 0, 7, 1], [2, 6, 0, 8, 2],
                  [5, 1, 0, 2, 4], [7, 2, 0, 1, 6],
                  [2, 8, 0, 6, 2], [1, 7, 0, 9, 1],
                  [4, 3, 0, 3, 2], [6, 4, 0, 2, 3]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"] * 2 + ["A", "A"],
                          "sample": ["S1"] * 4 + ["S2"] * 4 + ["S3"] * 2},
                         index=[f"c{i}" for i in range(10)]),
        var=pd.DataFrame(index=["P1", "P2", "ZERO", "R1", "R2"]),
    )
    enzyme = pd.DataFrame([
        dict(metabolite="Enzyme name", hmdb_id=h, gene=g, role="production", reaction="P")
        for h, g in [(H1, "P1"), (H2, "P2"), (H0, "ZERO")]
    ])
    sensor = pd.DataFrame([
        dict(metabolite="Interaction name", hmdb_id=h, sensor_gene=g, sensor_type=t)
        for h, g, t in [(H1, "R1", "Transporter"), (H1, "R2", "Cell surface receptor"),
                        (H2, "R1", "Transporter"), (H0, "R1", "Transporter")]
    ])
    return data, enzyme, sensor


def _run(data, enzyme, sensor, mode="sample_aware", **kwargs):
    return run_cell_mesh(data, enzyme, sensor, sample_mode=mode, sample_key="sample",
                         min_cells=1, n_perms=5, random_state=7,
                         store_null_scores=mode == "sample_aware", **kwargs)


def _values(frame, keys):
    out = frame.drop(columns="metabolite", errors="ignore").set_index(keys).sort_index()
    out.attrs = {}
    return out


def _sender_values(frame):
    out = frame.copy()
    out.index = out.index.droplevel("metabolite")
    return out.sort_index()


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("csv_input", [False, True])
@pytest.mark.parametrize("enzyme_name,sensor_name", [
    (None, "Interaction name"), ("Enzyme name", None), (pd.NA, pd.NA),
    ("  ", ""), ("Enzyme alias", "Interaction alias"),
])
def test_name_changes_preserve_scores_inference_and_sample_states(
    identity_case, tmp_path, mode, csv_input, enzyme_name, sensor_name,
):
    data, enzyme, sensor = identity_case
    expected = _run(data, enzyme, sensor, mode)
    renamed_enzyme = enzyme.assign(metabolite=enzyme_name)
    renamed_sensor = sensor.assign(metabolite=sensor_name)
    original_enzyme, original_sensor = renamed_enzyme.copy(), renamed_sensor.copy()
    if csv_input:
        enzyme_path, sensor_path = tmp_path / "enzyme.csv", tmp_path / "sensor.csv"
        renamed_enzyme.to_csv(enzyme_path, index=False)
        renamed_sensor.to_csv(sensor_path, index=False)
        inputs = enzyme_path, sensor_path
    else:
        inputs = renamed_enzyme, renamed_sensor
    actual = _run(data, *inputs, mode, n_jobs=2 if csv_input else 1)

    assert len(actual.events) == 16  # Four relations, four sender/receiver pairs.
    assert not actual.events.duplicated(KEY).any()
    pd.testing.assert_frame_equal(_values(actual.events, KEY), _values(expected.events, KEY), check_exact=True)
    pd.testing.assert_frame_equal(_sender_values(actual.sender_scores), _sender_values(expected.sender_scores), check_exact=True)
    receiver_key = ["hmdb_id", "sensor_gene", "receiver"]
    pd.testing.assert_frame_equal(_values(actual.receiver_scores, receiver_key),
                                  _values(expected.receiver_scores, receiver_key), check_exact=True)
    assert actual.events.loc[actual.events.hmdb_id.eq(H0), "perm_pvalue"].eq(1).all()
    if mode == "sample_aware":
        pd.testing.assert_frame_equal(_values(actual.sample_events, KEY + ["sample"]),
                                      _values(expected.sample_events, KEY + ["sample"]), check_exact=True)
        pd.testing.assert_frame_equal(_sender_values(actual.sample_sender_scores),
                                      _sender_values(expected.sample_sender_scores), check_exact=True)
        pd.testing.assert_frame_equal(_values(actual.sample_receiver_scores, receiver_key + ["sample"]),
                                      _values(expected.sample_receiver_scores, receiver_key + ["sample"]), check_exact=True)
        missing = actual.sample_events.loc[actual.sample_events.cell_mesh_score.isna()]
        assert len(missing) == 12
        assert missing["sample"].eq("S3").all()
        assert missing[["metabolite", "sensor_type"]].notna().all().all()
        assert missing["passes_min_cells"].isna().all()
        null = actual.events.attrs["sample_aware_null_scores"]
        assert null.index.names == KEY
        assert null.index.is_unique
        pd.testing.assert_frame_equal(null.sort_index(),
                                      expected.events.attrs["sample_aware_null_scores"].sort_index(), check_exact=True)
    pd.testing.assert_frame_equal(renamed_enzyme, original_enzyme)
    pd.testing.assert_frame_equal(renamed_sensor, original_sensor)


@pytest.mark.parametrize("source", ["enzyme", "interaction", "identifier"])
def test_display_name_priority_uses_first_nonempty_name(identity_case, source):
    data, enzyme, sensor = identity_case
    enzyme["metabolite"] = None
    sensor["metabolite"] = None
    if source == "enzyme":
        # The first row is missing, but a later retained prior row has a name.
        enzyme = pd.concat([enzyme, enzyme.assign(metabolite="Later enzyme name")], ignore_index=True)
        sensor["metabolite"] = "Interaction fallback"
        names = {h: "Later enzyme name" for h in [H1, H2, H0]}
    elif source == "interaction":
        sensor["metabolite"] = "Interaction fallback"
        names = {h: "Interaction fallback" for h in [H1, H2, H0]}
    else:
        names = {h: h for h in [H1, H2, H0]}
    result = _run(data, enzyme, sensor)
    for table in [result.events, result.receiver_scores, result.sample_events,
                  result.sample_receiver_scores, result.sender_scores.reset_index(),
                  result.sample_sender_scores.reset_index()]:
        assert table["metabolite"].tolist() == table.hmdb_id.map(names).tolist()
    for unit in result.availability_results["availability_by_sample"].values():
        assert unit["reaction_genes"].metabolite.tolist() == unit["reaction_genes"].hmdb_id.map(names).tolist()


def test_standalone_scores_keep_missing_names_as_display_only(identity_case):
    data, enzyme, sensor = identity_case
    enzyme["metabolite"] = None
    sensor["metabolite"] = None
    availability = compute_metabolite_availability(data, enzyme, min_cells=1)
    assert availability["availability"].index.get_level_values("metabolite").tolist() == [H1, H2, H0]
    receiver = compute_sensor_scores(data, sensor, min_cells=1)
    assert receiver.metabolite.tolist() == receiver.hmdb_id.tolist()


def test_sample_aggregation_ignores_names_from_individual_scoring_units(identity_case, monkeypatch):
    import cellmesh.core as core

    data, enzyme, sensor = identity_case
    expected = _run(data, enzyme, sensor)
    compute = core._compute_availability_scores

    def with_sample_aliases(adata, *args, **kwargs):
        sender, receiver, intermediates = compute(adata, *args, **kwargs)
        sample = adata.obs["sample"].iloc[0]
        name = None if sample == "S1" else f"Alias in {sample}"
        sender = sender.copy()
        sender.index = pd.MultiIndex.from_arrays(
            [[name] * len(sender), sender.index.get_level_values("hmdb_id")],
            names=["metabolite", "hmdb_id"],
        )
        return sender, receiver.assign(metabolite=name), intermediates

    # Protect aggregation itself: future caching/optimization must not require
    # every scoring unit to supply the same non-missing display name.
    monkeypatch.setattr(core, "_compute_availability_scores", with_sample_aliases)
    actual = _run(data, enzyme, sensor)
    pd.testing.assert_frame_equal(_values(actual.events, KEY), _values(expected.events, KEY), check_exact=True)
    pd.testing.assert_frame_equal(_sender_values(actual.sender_scores), _sender_values(expected.sender_scores), check_exact=True)
    receiver_key = ["hmdb_id", "sensor_gene", "receiver"]
    pd.testing.assert_frame_equal(_values(actual.receiver_scores, receiver_key),
                                  _values(expected.receiver_scores, receiver_key), check_exact=True)
    assert actual.events.metabolite.eq("Enzyme name").all()
    pd.testing.assert_frame_equal(actual.events.attrs["sample_aware_null_scores"].sort_index(),
                                  expected.events.attrs["sample_aware_null_scores"].sort_index(), check_exact=True)


def test_legacy_name_only_sender_matching_is_rejected(identity_case):
    data, enzyme, sensor = identity_case
    receiver = compute_sensor_scores(data, sensor, min_cells=1)
    sender = pd.DataFrame({"A": [0.5]}, index=["Interaction name"])
    with pytest.raises(ValueError, match="hmdb_id"):
        _make_cell_mesh_events(sender, receiver, allow_self=True)


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
def test_missing_names_export_without_losing_identifiers_or_na(identity_case, tmp_path, mode):
    data, enzyme, sensor = identity_case
    result = _run(data, enzyme.assign(metabolite=None), sensor.assign(metabolite=None), mode)
    prefix = tmp_path / "result"
    result.to_csv(prefix)
    tables = ["events", "sender_scores", "receiver_scores"]
    if mode == "sample_aware":
        tables += ["sample_events", "sample_sender_scores", "sample_receiver_scores"]
    for name in tables:
        path = tmp_path / f"result.{name}.csv"
        header = path.read_text().splitlines()[0].split(",")
        assert header.count("metabolite") == header.count("hmdb_id") == 1
        saved = pd.read_csv(path, dtype={"hmdb_id": str})
        assert len(saved) == len(getattr(result, name)) > 0
        assert saved.metabolite.tolist() == saved.hmdb_id.tolist()
        if name == "sample_events":
            assert saved.cell_mesh_score.isna().sum() == 12
            assert saved.loc[saved.hmdb_id.eq(H0), "cell_mesh_score"].dropna().eq(0).all()


@pytest.mark.parametrize("kind", ["counts", "network", "dotplot", "samples", "sender", "receiver"])
def test_all_plots_accept_results_scored_with_missing_names(identity_case, kind):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cellmesh import (plot_significant_event_counts, plot_communication_network,
                          plot_event_dotplot, plot_sample_event_scores,
                          plot_metabolite_secretion_violin, plot_receptor_expression_violin)
    data, enzyme, sensor = identity_case
    result = _run(data, enzyme.assign(metabolite=None), sensor.assign(metabolite=None))
    try:
        if kind == "counts":
            plotted = plot_significant_event_counts(result, max_fdr=None)
            assert plotted["counts"].at["A", "B"] == 4
        elif kind == "network":
            plotted = plot_communication_network(result, max_fdr=None)
            assert plotted["edge_table"].set_index(["sender", "receiver"]).at[("A", "B"), "event_count"] == 4
        elif kind == "dotplot":
            plotted = plot_event_dotplot(result, top_n=None)
            assert set(plotted["selected_event_keys"]) == {(H1, "R1"), (H1, "R2"), (H2, "R1"), (H0, "R1")}
            assert all(f"{H1} ({H1})" not in label for label in plotted["selected_events"])
        elif kind == "samples":
            plotted = plot_sample_event_scores(result, sender="A", receiver="B", hmdb_id=H1,
                                               sensor_gene="R1", qc_only=False)
            assert plotted["missing_samples"] == ["S3"]
            assert plotted["n_samples_coobserved"] == 2
        elif kind == "sender":
            plotted = plot_metabolite_secretion_violin(result, data, hmdb_id=H1)
            assert plotted["summary"].mean_production_score.sum() > 0
        else:
            plotted = plot_receptor_expression_violin(result, data, hmdb_id=H1, receptor_gene="R1")
            assert plotted["event_key"] == (H1, "R1")
    finally:
        plt.close("all")
