"""All visual identity uses HMDB/gene; metabolite names only label the output."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cellmesh import (
    plot_communication_network, plot_event_dotplot, plot_metabolite_secretion_violin,
    plot_receptor_expression_violin, plot_sample_event_scores,
    plot_significant_event_counts, run_cell_mesh,
)


H1, H2 = "HMDB0000001", "HMDB0000002"


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture
def aliases():
    rows = []
    for name, hmdb, gene, receiver, score in (
        ("First", H1, "R", "B", 0.5),
        ("Alias", " hmdb0000001 ", " R ", "B", 0.5),
        ("First", H2, "R", "B", 0.6),
        ("First", H1, "OTHER", "B", 0.7),
        ("First", H1, "R", "C", 0.8),
    ):
        rows.append(dict(metabolite=name, hmdb_id=hmdb, sensor_gene=gene, sender="A", receiver=receiver,
                         sensor_type="Transporter", cell_mesh_score=score,
                         perm_pvalue=0.01, fdr_sensor_type=0.02, passes_min_cells=True))
    return pd.DataFrame(rows)


@pytest.mark.parametrize("kind", ["counts", "network", "dotplot"])
@pytest.mark.parametrize("names", ["aliases", "renamed", "absent"])
def test_event_plots_identify_relations_without_names(aliases, kind, names):
    events = aliases.copy()
    if names == "renamed":
        events.metabolite = [f"Display {i}" for i in range(len(events))]
    elif names == "absent":
        events = events.drop(columns="metabolite")
    original = events.copy(deep=True)
    if kind == "counts":
        result = plot_significant_event_counts(events, max_fdr=None)
        assert result["counts"].at["A", "B"] == 3
        assert result["counts"].at["A", "C"] == 1
    elif kind == "network":
        result = plot_communication_network(events, max_fdr=None)
        pairs = result["edge_table"].set_index(["sender", "receiver"])
        assert pairs.at[("A", "B"), "event_count"] == 3
        assert pairs.at[("A", "B"), "total_score"] == pytest.approx(1.8)
        assert pairs.at[("A", "C"), "event_count"] == 1
    else:
        result = plot_event_dotplot(events, top_n=None)
        assert set(result["selected_event_keys"]) == {(H1, "R"), (H1, "OTHER"), (H2, "R")}
        assert len(result["plot_events"]) == 4
        for (identifier, gene), label in zip(result["selected_event_keys"], result["selected_events"]):
            assert identifier in label and gene in label
    pd.testing.assert_frame_equal(events, original)


@pytest.mark.parametrize("plot", [plot_significant_event_counts, plot_communication_network])
def test_custom_keys_cannot_reintroduce_name_identity_or_omit_receptor(aliases, plot):
    with pytest.raises(ValueError, match="metabolite is display-only"):
        plot(aliases, unique_keys=["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene"])
    with pytest.raises(ValueError, match="sensor_gene"):
        plot(aliases, unique_keys=["sender", "receiver", "hmdb_id"])
    sample_events = pd.concat([aliases.assign(sample="S1"), aliases.assign(sample="S2")], ignore_index=True)
    plotted = plot(sample_events, unique_keys=["sender", "receiver", "hmdb_id", "sensor_gene", "sample"])
    assert plotted["thresholds"]["unique_events_after_filter"] == 8


def test_dotplot_top_n_and_legacy_structured_selection_ignore_name(aliases):
    selected = plot_event_dotplot(aliases, event_keys=[("Not in data", " hmdb0000001 ", " R ")])
    assert selected["selected_event_keys"] == [(H1, "R")]
    assert len(selected["plot_events"]) == 2
    by_dict = plot_event_dotplot(
        aliases, event_keys=[dict(hmdb_id=H1, sensor_gene="R", metabolite="Unused label")],
    )
    assert by_dict["selected_event_keys"] == selected["selected_event_keys"]
    top = plot_event_dotplot(aliases, top_n=1)
    assert top["selected_event_keys"] == [(H1, "R")]
    # A returned display label is not a lookup key; only ID-bearing structures match.
    with pytest.raises(ValueError, match="display labels cannot be matched"):
        plot_event_dotplot(aliases, event_keys=selected["selected_events"])


@pytest.fixture
def samples(aliases):
    rows = pd.concat([aliases.iloc[[0]]] * 3, ignore_index=True)
    rows["sample"] = ["S1", "S2", "S3"]
    rows["metabolite"] = ["One", "Two", "Three"]
    rows["cell_mesh_score"] = [0.0, 0.8, np.nan]
    rows["passes_min_cells"] = pd.Series([True, True, pd.NA], dtype="boolean")
    return rows


def test_sample_scores_include_all_id_matches_and_keep_zero_and_missing(samples):
    original = samples.copy(deep=True)
    duplicate = samples.iloc[[0]].assign(metabolite="Another alias")
    result = plot_sample_event_scores(pd.concat([samples, duplicate], ignore_index=True), sender="A", receiver="B",
                                      hmdb_id=" hmdb0000001 ", sensor_gene=" R ", metabolite="Display only", qc_only=False)
    assert result["sample_scores"]["sample"].astype(str).tolist() == ["S1", "S2", "S3"]
    assert result["median_score"] == 0.4
    assert result["n_samples_coobserved"] == 2
    assert result["missing_samples"] == ["S3"]
    assert "metabolite" not in result["event_key"]
    assert f"Display only ({H1}) -> R" in result["ax"].get_title()
    pd.testing.assert_frame_equal(samples, original)


def test_conflicting_sample_values_are_not_resolved_by_a_name(samples):
    conflict = samples.iloc[[0]].assign(metabolite="Alias", cell_mesh_score=0.2)
    with pytest.raises(ValueError, match="conflicting duplicate"):
        plot_sample_event_scores(pd.concat([samples, conflict], ignore_index=True), sender="A", receiver="B",
                                 hmdb_id=H1, sensor_gene="R", qc_only=False)


def test_sample_names_cannot_be_used_without_an_id(samples):
    with pytest.raises(ValueError, match="provide hmdb_id"):
        plot_sample_event_scores(samples, sender="A", receiver="B", metabolite="One", sensor_gene="R")


def test_sample_plot_does_not_require_display_or_sensor_type_metadata(samples):
    plotted = plot_sample_event_scores(
        samples.drop(columns=["metabolite", "sensor_type"]),
        sender="A", receiver="B", hmdb_id=H1, sensor_gene="R", qc_only=False,
    )
    assert plotted["median_score"] == 0.4
    assert plotted["missing_samples"] == ["S3"]
    assert plotted["ax"].get_title().startswith(f"{H1} -> R")


@pytest.fixture(params=["pooled_stratified", "sample_aware"])
def scored(request):
    data = AnnData(np.array([[8., 1.], [0., 2.], [6., 1.], [0., 2.]]),
                   obs=pd.DataFrame({"cell_type": ["A", "B", "A", "B"], "sample": ["S1", "S1", "S2", "S2"]},
                                    index=["c1", "c2", "c3", "c4"]),
                   var=pd.DataFrame(index=["PROD", "R"]))
    enzyme = pd.DataFrame([dict(metabolite="EnzymeName", hmdb_id=H1, gene="PROD", role="production", reaction="RP")])
    sensor = pd.DataFrame([dict(metabolite="InteractionName", hmdb_id=H1, sensor_gene="R", sensor_type="Transporter")])
    result = run_cell_mesh(data, enzyme, sensor, sample_mode=request.param, sample_key="sample", min_cells=1, n_perms=3)
    return data, result


@pytest.mark.parametrize("name", [None, "Display alias"])
def test_both_violins_locate_by_ids_and_label_the_metabolite(scored, name):
    data, result = scored
    before = result.events.copy(deep=True)
    sender = plot_metabolite_secretion_violin(result, data, hmdb_id=" hmdb0000001 ", metabolite=name)
    receiver = plot_receptor_expression_violin(result, data, hmdb_id=" hmdb0000001 ", receptor_gene=" R ", metabolite=name)
    assert sender["hmdb_id"] == receiver["hmdb_id"] == H1
    assert H1 in sender["ax"].get_title()
    assert H1 in receiver["ax"].get_title()
    assert receiver["event_key"] == (H1, "R")
    if name:
        assert name in sender["ax"].get_title() and name in receiver["ax"].get_title()
    pd.testing.assert_frame_equal(result.events, before)


def test_sender_name_alone_cannot_select_even_one_available_metabolite(scored):
    data, result = scored
    with pytest.raises(ValueError, match="provide hmdb_id"):
        plot_metabolite_secretion_violin(result, data, metabolite="EnzymeName")


def test_both_violins_work_without_metabolite_names(scored):
    data, result = scored
    modified = deepcopy(result)
    modified.sender_scores.index = modified.sender_scores.index.droplevel("metabolite")
    modified.receiver_scores = modified.receiver_scores.drop(columns="metabolite")
    units = modified.availability_results.get("availability_by_sample", {"pooled": modified.availability_results})
    for unit in units.values():
        unit["reaction_genes"] = unit["reaction_genes"].drop(columns="metabolite")
        for key in ("P", "C", "E"):
            unit[key].index = unit[key].index.droplevel("metabolite")
    sender = plot_metabolite_secretion_violin(modified, data, hmdb_id=H1)
    receiver = plot_receptor_expression_violin(modified, data, hmdb_id=H1, receptor_gene="R")
    assert sender["metabolite"] is None
    assert sender["hmdb_id"] == receiver["hmdb_id"] == H1
    assert sender["ax"].get_title().startswith(f"{H1} secretion")
    assert receiver["ax"].get_title().startswith(f"{H1} -> R")
    assert sender["summary"].mean_production_score.sum() == pytest.approx(7.0)


def test_sender_reaction_aliases_do_not_double_capacity(scored):
    data, result = scored
    reference = plot_metabolite_secretion_violin(result, data, hmdb_id=H1)
    modified = deepcopy(result)
    modified.sender_scores = pd.concat([modified.sender_scores, modified.sender_scores.rename(index={"EnzymeName": "Other name"}, level=0)])
    units = modified.availability_results.get("availability_by_sample", {"pooled": modified.availability_results})
    for unit in units.values():
        reactions = unit["reaction_genes"]
        unit["reaction_genes"] = pd.concat([reactions, reactions.assign(metabolite="Another name")], ignore_index=True)
    actual = plot_metabolite_secretion_violin(modified, data, hmdb_id=H1)
    np.testing.assert_array_equal(actual["plot_data"].production_score, reference["plot_data"].production_score)
    assert len(actual["reaction_definitions"]) == len(reference["reaction_definitions"])


def test_sender_conflicting_capacity_rows_cannot_be_chosen_by_name(scored):
    data, result = scored
    modified = deepcopy(result)
    conflict = (modified.sender_scores + 0.1).rename(index={"EnzymeName": "Other name"}, level=0)
    modified.sender_scores = pd.concat([modified.sender_scores, conflict])
    with pytest.raises(ValueError, match="Conflicting values"):
        plot_metabolite_secretion_violin(modified, data, hmdb_id=H1, metabolite="EnzymeName")


def test_conflicting_reaction_definitions_are_not_resolved_by_a_name(scored):
    data, result = scored
    modified = deepcopy(result)
    units = modified.availability_results.get("availability_by_sample", {"pooled": modified.availability_results})
    for unit in units.values():
        reactions = unit["reaction_genes"]
        conflict = reactions.iloc[[0]].assign(metabolite="Another name", genes=[("R",)])
        unit["reaction_genes"] = pd.concat([reactions, conflict], ignore_index=True)
    with pytest.raises(ValueError, match="Conflicting reaction definitions"):
        plot_metabolite_secretion_violin(modified, data, hmdb_id=H1, metabolite="EnzymeName")


def test_receiver_alias_rows_are_unique_and_conflicts_are_rejected(scored):
    data, result = scored
    modified = deepcopy(result)
    alias = modified.receiver_scores.assign(metabolite="Alias")
    modified.receiver_scores = pd.concat([modified.receiver_scores, alias], ignore_index=True)
    actual = plot_receptor_expression_violin(modified, data, hmdb_id=H1, receptor_gene="R")
    assert len(actual["selected_receiver_scores"]) == len(result.receiver_scores)
    modified.receiver_scores.loc[0, "sensor_score"] += 0.1
    with pytest.raises(ValueError, match="Conflicting receiver scores"):
        plot_receptor_expression_violin(modified, data, hmdb_id=H1, receptor_gene="R")
