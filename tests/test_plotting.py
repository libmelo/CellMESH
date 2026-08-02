from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cellmesh
from cellmesh.plotting import (
    plot_communication_network,
    plot_event_dotplot,
    plot_metabolite_secretion_violin,
    plot_receptor_expression_violin,
    plot_sample_event_scores,
    plot_significant_event_counts,
)


@pytest.fixture
def events():
    return pd.DataFrame(
        [
            {
                "sender": "A",
                "receiver": "B",
                "metabolite": "M1",
                "hmdb_id": "HMDB1",
                "sensor_gene": "G1",
                "sensor_type": "Transporter",
                "cell_mesh_score": 0.8,
                "perm_pvalue": 0.01,
                "fdr_sensor_type": 0.02,
                "fdr_global": 0.10,
            },
            {
                "sender": "A",
                "receiver": "B",
                "metabolite": "M2",
                "hmdb_id": "HMDB2",
                "sensor_gene": "G2",
                "sensor_type": "Cell surface receptor",
                "cell_mesh_score": 0.6,
                "perm_pvalue": 0.02,
                "fdr_sensor_type": 0.04,
                "fdr_global": 0.20,
            },
            {
                "sender": "B",
                "receiver": "A",
                "metabolite": "M3",
                "hmdb_id": "HMDB3",
                "sensor_gene": "G3",
                "sensor_type": "Transporter",
                "cell_mesh_score": 0.5,
                "perm_pvalue": 0.03,
                "fdr_sensor_type": 0.20,
                "fdr_global": 0.03,
            },
        ]
    )


@pytest.fixture
def qc_events(events):
    qc_events = events.copy()
    qc_events["passes_min_cells"] = [True, False, True]
    return qc_events


@pytest.fixture
def single_cell_violin_inputs():
    genes = ["P1", "P2", "C1", "E1", "G1"]
    expression = np.array(
        [
            [4.0, 4.0, 0.0, 2.0, 0.0],
            [3.0, 3.0, 0.0, 1.0, 1.0],
            [2.0, 2.0, 0.0, 1.0, 0.0],
            [1.0, 1.0, 2.0, 0.0, 5.0],
            [0.5, 0.5, 2.0, 0.0, 4.0],
            [0.0, 0.0, 2.0, 0.0, 3.0],
        ]
    )
    obs = pd.DataFrame(
        {"cluster": ["A", "A", "A", "B", "B", "B"]},
        index=[f"cell_{index}" for index in range(6)],
    )
    adata = SimpleNamespace(
        X=expression,
        layers={},
        var_names=pd.Index(genes),
        obs=obs,
    )

    metabolite_index = pd.MultiIndex.from_tuples(
        [("M1", "HMDB1")],
        names=["metabolite", "hmdb_id"],
    )
    sender_scores = pd.DataFrame([[0.8, 0.1]], index=metabolite_index, columns=["A", "B"])
    receiver_scores = pd.DataFrame(
        {
            "metabolite": ["M1", "M1"],
            "hmdb_id": ["HMDB1", "HMDB1"],
            "sensor_gene": ["G1", "G1"],
            "sensor_type": ["Cell surface receptor", "Cell surface receptor"],
            "receiver": ["A", "B"],
            "sensor_score": [0.0, 0.75],
            "sensor_expr_frac": [1 / 3, 1.0],
            "receiver_n_cells": [3, 3],
        }
    )
    reaction_genes = pd.DataFrame(
        {
            "metabolite": ["M1", "M1", "M1"],
            "hmdb_id": ["HMDB1", "HMDB1", "HMDB1"],
            "reaction": ["R_product", "R_consume", "R_export"],
            "direction": ["product", "substrate", "exporter"],
            "genes": [["P1", "P2"], ["C1"], ["E1"]],
        }
    )
    availability_results = {
        "reaction_genes": reaction_genes,
        "P": pd.DataFrame([[3.0, 0.5]], index=metabolite_index, columns=["A", "B"]),
        "C": pd.DataFrame([[0.0, 2.0]], index=metabolite_index, columns=["A", "B"]),
        "E": pd.DataFrame([[1.3, 0.0]], index=metabolite_index, columns=["A", "B"]),
    }
    result = SimpleNamespace(
        sender_scores=sender_scores,
        receiver_scores=receiver_scores,
        availability_results=availability_results,
        parameters={"cell_type_key": "cluster", "layer": None},
    )
    return result, adata


def test_plotting_functions_are_available_from_public_imports():
    canonical_functions = (
        plot_significant_event_counts,
        plot_communication_network,
        plot_event_dotplot,
        plot_sample_event_scores,
        plot_metabolite_secretion_violin,
        plot_receptor_expression_violin,
    )

    assert all(function.__module__ == "cellmesh.plotting" for function in canonical_functions)
    assert cellmesh.plot_significant_event_counts is plot_significant_event_counts
    assert cellmesh.plot_communication_network is plot_communication_network
    assert cellmesh.plot_event_dotplot is plot_event_dotplot
    assert cellmesh.plot_sample_event_scores is plot_sample_event_scores
    assert cellmesh.plot_metabolite_secretion_violin is plot_metabolite_secretion_violin
    assert cellmesh.plot_receptor_expression_violin is plot_receptor_expression_violin


def test_significant_event_counts_supports_selectable_fdr_column(events):
    sensor_type_result = plot_significant_event_counts(events, max_fdr=0.05)
    global_result = plot_significant_event_counts(
        events,
        max_fdr=0.05,
        fdr_col="fdr_global",
    )

    assert sensor_type_result["counts"].loc["A", "B"] == 2
    assert sensor_type_result["counts"].loc["B", "A"] == 0
    assert global_result["counts"].loc["A", "B"] == 0
    assert global_result["counts"].loc["B", "A"] == 1
    assert global_result["thresholds"]["fdr_col"] == "fdr_global"

    plt.close(sensor_type_result["fig"])
    plt.close(global_result["fig"])


def test_significant_event_counts_defaults_to_min_cells_qc_and_can_show_all(qc_events):
    qc_result = plot_significant_event_counts(qc_events, max_fdr=None)
    all_result = plot_significant_event_counts(
        qc_events,
        max_fdr=None,
        qc_only=False,
    )

    assert qc_result["counts"].loc["A", "B"] == 1
    assert all_result["counts"].loc["A", "B"] == 2
    assert qc_result["qc"]["qc_applied"] is True
    assert qc_result["qc"]["events_excluded_by_qc"] == 1
    assert qc_result["qc_excluded_events"]["metabolite"].tolist() == ["M2"]
    assert qc_result["qc_failed_events"]["metabolite"].tolist() == ["M2"]
    assert qc_result["thresholds"]["events_before_qc"] == 3
    assert qc_result["thresholds"]["events_after_qc"] == 2
    assert all_result["qc"]["qc_applied"] is False
    assert all_result["qc"]["events_excluded_by_qc"] == 0
    assert all_result["qc_excluded_events"].empty
    assert all_result["qc"]["events_failing_qc"] == 1

    plt.close(qc_result["fig"])
    plt.close(all_result["fig"])


def test_significant_event_counts_returns_complete_summary(events):
    many_events = pd.concat(
        [
            events.iloc[[0]].assign(
                sender=f"S{index}",
                receiver=f"R{index}",
                metabolite=f"M{index}",
                hmdb_id=f"HMDB{index}",
                sensor_gene=f"G{index}",
            )
            for index in range(21)
        ],
        ignore_index=True,
    )

    result = plot_significant_event_counts(many_events, max_fdr=None)

    assert len(result["summary"]) == 21
    assert "all_summary" not in result
    with pytest.raises(TypeError, match="top_n_summary"):
        plot_significant_event_counts(events, top_n_summary=1)
    plt.close(result["fig"])


def test_communication_network_aggregates_and_encodes_events(events):
    duplicate = events.iloc[[0]].copy()
    duplicate["cell_mesh_score"] = 0.2
    duplicate["fdr_sensor_type"] = 0.03
    self_event = events.iloc[[0]].copy()
    self_event["sender"] = "C"
    self_event["receiver"] = "C"
    self_event["metabolite"] = "M4"
    self_event["hmdb_id"] = "HMDB4"
    self_event["sensor_gene"] = "G4"
    self_event["cell_mesh_score"] = 0.4
    network_events = pd.concat([events, duplicate, self_event], ignore_index=True)

    result = plot_communication_network(network_events, max_fdr=None)
    edges = result["edge_table"].set_index(["sender", "receiver"])
    nodes = result["node_table"].set_index("cell_type")

    assert edges.loc[("A", "B"), "event_count"] == 2
    assert np.isclose(edges.loc[("A", "B"), "total_score"], 1.4)
    assert edges.loc[("A", "B"), "edge_width"] > edges.loc[("B", "A"), "edge_width"]
    assert edges.loc[("A", "B"), "edge_color_value"] == pytest.approx(1.4)
    assert nodes.loc["A", "connected_cell_type_count"] == 1
    assert nodes.loc["A", "connection_count"] == 2
    assert nodes.loc["C", "connection_count"] == 1
    assert len(result["unique_events"]) == 4
    assert len(result["edge_artists"]) == 3
    assert result["colorbar"] is not None
    assert result["encodings"]["edge_width_by"] == "event_count"
    assert result["encodings"]["edge_color_by"] == "total_score"

    plt.close(result["fig"])


def test_communication_network_defaults_to_min_cells_qc_and_can_show_all(qc_events):
    qc_result = plot_communication_network(qc_events, max_fdr=None)
    all_result = plot_communication_network(
        qc_events,
        max_fdr=None,
        qc_only=False,
    )

    qc_edges = qc_result["edge_table"].set_index(["sender", "receiver"])
    all_edges = all_result["edge_table"].set_index(["sender", "receiver"])
    assert qc_edges.loc[("A", "B"), "event_count"] == 1
    assert qc_edges.loc[("A", "B"), "total_score"] == pytest.approx(0.8)
    assert all_edges.loc[("A", "B"), "event_count"] == 2
    assert all_edges.loc[("A", "B"), "total_score"] == pytest.approx(1.4)
    assert qc_result["qc"]["events_excluded_by_qc"] == 1
    assert qc_result["qc_excluded_events"]["metabolite"].tolist() == ["M2"]
    assert all_result["qc_excluded_events"].empty

    plt.close(qc_result["fig"])
    plt.close(all_result["fig"])


def test_communication_network_default_fdr_filter(events):
    result = plot_communication_network(events)

    edges = result["edge_table"].set_index(["sender", "receiver"])
    assert list(edges.index) == [("A", "B")]
    assert edges.loc[("A", "B"), "event_count"] == 2
    assert result["thresholds"]["events_after_filter"] == 2
    plt.close(result["fig"])


def test_event_dotplot_smoke_and_returns_selected_data(events):
    result = plot_event_dotplot(events, top_n=2, max_fdr=None)

    assert len(result["selected_events"]) == 2
    assert not result["plot_events"].empty
    assert result["thresholds"]["events_plotted"] == len(result["plot_events"])
    plt.close(result["fig"])


def test_event_dotplot_existing_ax_requires_exactly_one_sender(events):
    fig, ax = plt.subplots()
    with pytest.raises(ValueError, match="exactly one sender"):
        plot_event_dotplot(events, top_n=None, max_fdr=None, ax=ax)
    plt.close(fig)


@pytest.mark.parametrize("top_n", [True, 1.5, "2"])
def test_event_dotplot_rejects_non_integer_top_n(events, top_n):
    with pytest.raises(TypeError, match="top_n"):
        plot_event_dotplot(events, top_n=top_n, max_fdr=None)


@pytest.mark.parametrize("top_n", [0, -1])
def test_event_dotplot_rejects_non_positive_top_n(events, top_n):
    with pytest.raises(ValueError, match="top_n"):
        plot_event_dotplot(events, top_n=top_n, max_fdr=None)

    fig, ax = plt.subplots()
    result = plot_event_dotplot(
        events,
        top_n=None,
        max_fdr=None,
        sender_labels=["A"],
        ax=ax,
    )
    assert result["sender_order"] == ["A"]
    assert set(result["plot_events"]["sender"]) == {"A"}
    plt.close(fig)


def test_event_dotplot_defaults_to_min_cells_qc_and_can_show_all(qc_events):
    qc_result = plot_event_dotplot(qc_events, top_n=None, max_fdr=None)
    all_result = plot_event_dotplot(
        qc_events,
        top_n=None,
        max_fdr=None,
        qc_only=False,
    )

    assert set(qc_result["plot_events"]["metabolite"]) == {"M1", "M3"}
    assert set(all_result["plot_events"]["metabolite"]) == {"M1", "M2", "M3"}
    assert qc_result["qc"]["events_excluded_by_qc"] == 1
    assert qc_result["qc_excluded_events"]["metabolite"].tolist() == ["M2"]
    assert all_result["qc_excluded_events"].empty

    plt.close(qc_result["fig"])
    plt.close(all_result["fig"])


def test_sample_event_scores_preserves_missing_samples_as_na():
    sample_events = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3"],
            "sender": ["A", "A", "A"],
            "receiver": ["B", "B", "B"],
            "metabolite": ["M1", "M1", "M1"],
            "hmdb_id": ["HMDB1", "HMDB1", "HMDB1"],
            "sensor_gene": ["G1", "G1", "G1"],
            "sensor_type": ["Transporter", "Transporter", "Transporter"],
            "cell_mesh_score": [0.4, np.nan, 0.8],
        }
    )
    result_object = SimpleNamespace(sample_events=sample_events)

    result = plot_sample_event_scores(
        result_object,
        sender="A",
        receiver="B",
        metabolite="M1",
        hmdb_id="HMDB1",
        sensor_gene="G1",
        sensor_type="Transporter",
    )

    assert result["n_samples_coobserved"] == 2
    assert result["missing_samples"] == ["S2"]
    assert np.isclose(result["median_score"], 0.6)
    assert pd.isna(result["sample_scores"].loc[1, "cell_mesh_score"])
    assert result["qc"]["qc_column_present"] is False
    assert result["qc"]["qc_applied"] is False
    assert result["qc_excluded_events"].empty
    plt.close(result["fig"])


def test_sample_event_scores_defaults_to_min_cells_qc_and_can_show_all():
    sample_events = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3"],
            "sender": ["A", "A", "A"],
            "receiver": ["B", "B", "B"],
            "metabolite": ["M1", "M1", "M1"],
            "hmdb_id": ["HMDB1", "HMDB1", "HMDB1"],
            "sensor_gene": ["G1", "G1", "G1"],
            "sensor_type": ["Transporter", "Transporter", "Transporter"],
            "cell_mesh_score": [0.4, 1.0, 0.8],
            "passes_min_cells": [True, False, True],
        }
    )
    kwargs = {
        "sender": "A",
        "receiver": "B",
        "metabolite": "M1",
        "hmdb_id": "HMDB1",
        "sensor_gene": "G1",
        "sensor_type": "Transporter",
    }

    qc_result = plot_sample_event_scores(sample_events, **kwargs)
    all_result = plot_sample_event_scores(sample_events, qc_only=False, **kwargs)

    assert qc_result["sample_scores"]["sample"].astype(str).tolist() == ["S1", "S3"]
    assert all_result["sample_scores"]["sample"].astype(str).tolist() == [
        "S1",
        "S2",
        "S3",
    ]
    assert qc_result["median_score"] == pytest.approx(0.6)
    assert all_result["median_score"] == pytest.approx(0.8)
    assert qc_result["qc"]["events_excluded_by_qc"] == 1
    assert qc_result["qc_excluded_events"]["sample"].tolist() == ["S2"]
    assert all_result["qc_excluded_events"].empty
    assert all_result["qc_failed_events"]["sample"].tolist() == ["S2"]

    plt.close(qc_result["fig"])
    plt.close(all_result["fig"])


def test_metabolite_secretion_violin_uses_cellwise_reactions_and_availability_color(
    single_cell_violin_inputs,
):
    result_object, adata = single_cell_violin_inputs

    result = plot_metabolite_secretion_violin(
        result_object,
        adata,
        metabolite="M1",
        sender_labels=["B", "A"],
    )

    summary = result["summary"].set_index("sender")
    assert result["cell_order"] == ["B", "A"]
    assert result["hmdb_id"] == "HMDB1"
    assert len(result["plot_data"]) == 6
    assert summary.loc["A", "mean_production_score"] == pytest.approx(3.0)
    assert summary.loc["B", "mean_production_score"] == pytest.approx(0.5)
    assert summary.loc["A", "metabolite_availability"] == pytest.approx(0.8)
    assert summary.loc["B", "abundance_adjusted_consumption"] == pytest.approx(2.0)
    assert set(result["reaction_metadata"]["direction"]) == {
        "product",
        "substrate",
        "exporter",
    }
    plt.close(result["fig"])


def test_receptor_expression_violin_reports_expression_and_sensor_metadata(
    single_cell_violin_inputs,
):
    result_object, adata = single_cell_violin_inputs

    result = plot_receptor_expression_violin(
        result_object,
        adata,
        receptor_gene="G1",
        metabolite="M1",
        receiver_labels=["B", "A"],
    )

    summary = result["summary"].set_index("receiver")
    assert result["cell_order"] == ["B", "A"]
    assert len(result["plot_data"]) == 6
    assert summary.loc["B", "mean_receptor_expression"] == pytest.approx(4.0)
    assert summary.loc["A", "mean_receptor_expression"] == pytest.approx(1 / 3)
    assert summary.loc["B", "sensor_score"] == pytest.approx(0.75)
    assert summary.loc["A", "sensor_expr_frac"] == pytest.approx(1 / 3)
    plt.close(result["fig"])


def test_receptor_expression_violin_requires_hmdb_for_ambiguous_metabolite(
    single_cell_violin_inputs,
):
    result_object, adata = single_cell_violin_inputs
    second_context = result_object.receiver_scores.copy()
    second_context["hmdb_id"] = "HMDB2"
    ambiguous_result = SimpleNamespace(
        **{
            **vars(result_object),
            "receiver_scores": pd.concat(
                [result_object.receiver_scores, second_context],
                ignore_index=True,
            ),
        }
    )

    with pytest.raises(ValueError, match="provide hmdb_id"):
        plot_receptor_expression_violin(
            ambiguous_result,
            adata,
            receptor_gene="G1",
            metabolite="M1",
        )

    resolved = plot_receptor_expression_violin(
        ambiguous_result,
        adata,
        receptor_gene="G1",
        metabolite="M1",
        hmdb_id="HMDB1",
    )
    assert resolved["selected_receiver_scores"]["hmdb_id"].unique().tolist() == [
        "HMDB1"
    ]
    plt.close(resolved["fig"])


def test_single_cell_violins_accept_sparse_real_cellmesh_result():
    anndata = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    genes = ["P1", "P2", "C1", "E1", "G1"]
    expression = sparse.csr_matrix(
        np.array(
            [
                [4.0, 4.0, 0.0, 2.0, 0.0],
                [3.0, 3.0, 0.0, 1.0, 1.0],
                [2.0, 2.0, 0.0, 1.0, 0.0],
                [1.0, 1.0, 2.0, 0.0, 5.0],
                [0.5, 0.5, 2.0, 0.0, 4.0],
                [0.0, 0.0, 2.0, 0.0, 3.0],
            ]
        )
    )
    adata = anndata.AnnData(
        X=expression,
        obs=pd.DataFrame({"cluster": ["A", "A", "A", "B", "B", "B"]}),
        var=pd.DataFrame(index=genes),
    )
    enzyme_prior = pd.DataFrame(
        {
            "metabolite": ["M1", "M1", "M1", "M1"],
            "hmdb_id": ["HMDB1"] * 4,
            "gene": ["P1", "P2", "C1", "E1"],
            "reaction": ["R_product", "R_product", "R_consume", "R_export"],
            "role": ["production", "production", "degradation", "export"],
        }
    )
    sensor_prior = pd.DataFrame(
        {
            "metabolite": ["M1"],
            "hmdb_id": ["HMDB1"],
            "sensor_gene": ["G1"],
            "sensor_type": ["Cell surface receptor"],
        }
    )
    result_object = cellmesh.run_cell_mesh(
        adata,
        enzyme_metabolite=enzyme_prior,
        metabolite_sensor=sensor_prior,
        cell_type_key="cluster",
        min_cells=1,
    )

    metabolite_result = plot_metabolite_secretion_violin(
        result_object,
        adata,
        metabolite="M1",
    )
    receptor_result = plot_receptor_expression_violin(
        result_object,
        adata,
        receptor_gene="G1",
    )

    assert len(metabolite_result["plot_data"]) == adata.n_obs
    assert len(receptor_result["plot_data"]) == adata.n_obs
    assert set(metabolite_result["summary"]["sender"]) == {"A", "B"}
    assert set(receptor_result["summary"]["receiver"]) == {"A", "B"}
    plt.close(metabolite_result["fig"])
    plt.close(receptor_result["fig"])
