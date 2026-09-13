"""Receptor plots must locate event context by HMDB, across prior-name aliases."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cellmesh import plot_receptor_expression_violin, run_cell_mesh


HMDB1 = "HMDB0000001"
HMDB2 = "HMDB0000002"


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture
def receiver_context():
    adata = AnnData(
        np.array([[1., 4.], [2., 3.], [3., 2.], [4., 1.]]),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}, index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["R", "OTHER"]),
    )
    rows = []
    for hmdb, gene, scores in ((HMDB1, "R", [0.2, 0.4]), (HMDB2, "R", [0.7, 0.8]),
                               (HMDB1, "OTHER", [0.9, 0.3])):
        for receiver, score in zip(("A", "B"), scores):
            rows.append(dict(metabolite="Shared", hmdb_id=hmdb, sensor_gene=gene,
                             sensor_type="Transporter", receiver=receiver,
                             sensor_score=score, sensor_expr_frac=1.0))
    result = SimpleNamespace(receiver_scores=pd.DataFrame(rows), parameters={"cell_type_key": "cell_type"})
    return result, adata


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("legacy_receiver_alias", [False, True])
def test_receptor_plot_accepts_actual_event_name_when_interaction_uses_an_alias(
    sample_mode, legacy_receiver_alias,
):
    adata = AnnData(
        np.array([[8., 1.], [0., 2.], [6., 1.], [0., 2.]]),
        obs=pd.DataFrame({"cell_type": ["A", "B", "A", "B"], "sample": ["S1", "S1", "S2", "S2"]},
                         index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["PROD", "R"]),
    )
    enzyme = pd.DataFrame([
        dict(metabolite="EnzymeName", hmdb_id=HMDB1, gene="PROD", role="production", reaction="RP"),
    ])
    sensor = pd.DataFrame([
        dict(metabolite="InteractionAlias", hmdb_id=HMDB1, sensor_gene="R", sensor_type="Transporter"),
    ])
    result = run_cell_mesh(adata, enzyme, sensor, min_cells=1, n_perms=5,
                           sample_mode=sample_mode, sample_key="sample")
    # Scoring now shares the preferred Enzyme display name across result tables,
    # while caller-owned Interaction metadata and old result aliases stay valid.
    assert result.receiver_scores.metabolite.eq("EnzymeName").all()
    assert sensor.metabolite.eq("InteractionAlias").all()
    if legacy_receiver_alias:
        result.receiver_scores["metabolite"] = "InteractionAlias"
    events_before = result.events.copy(deep=True)
    receiver_before = result.receiver_scores.copy(deep=True)
    obs_before = adata.obs.copy(deep=True)
    matrix_before = adata.X.copy()
    event = result.events.iloc[0]
    assert event.metabolite == "EnzymeName"
    plotted = plot_receptor_expression_violin(
        result, adata, receptor_gene=event.sensor_gene,
        metabolite=event.metabolite, hmdb_id=event.hmdb_id,
    )
    control = plot_receptor_expression_violin(result, adata, receptor_gene="R", hmdb_id=HMDB1)
    pd.testing.assert_frame_equal(plotted["summary"], control["summary"])
    pd.testing.assert_frame_equal(plotted["plot_data"], control["plot_data"])
    pd.testing.assert_frame_equal(plotted["selected_receiver_scores"], receiver_before)
    np.testing.assert_array_equal(plotted["plot_data"].receptor_expression, matrix_before[:, 1])
    pd.testing.assert_frame_equal(result.events, events_before)
    pd.testing.assert_frame_equal(result.receiver_scores, receiver_before)
    pd.testing.assert_frame_equal(adata.obs, obs_before)
    np.testing.assert_array_equal(adata.X, matrix_before)


@pytest.mark.parametrize("requested, stored", [
    (HMDB1, HMDB1), (" hmdb0000001 ", HMDB1), (HMDB1, " hMdB0000001 "),
])
@pytest.mark.parametrize("metabolite", [None, "EnzymeAlias"])
def test_explicit_hmdb_uses_shared_normalization_and_keeps_receptor_context(
    receiver_context, requested, stored, metabolite,
):
    result, adata = receiver_context
    result.receiver_scores.loc[result.receiver_scores.hmdb_id.eq(HMDB1), "hmdb_id"] = stored
    original = result.receiver_scores.copy(deep=True)
    plotted = plot_receptor_expression_violin(
        result, adata, receptor_gene="R", metabolite=metabolite, hmdb_id=requested,
    )
    expected = original.loc[original.hmdb_id.eq(stored) & original.sensor_gene.eq("R")]
    pd.testing.assert_frame_equal(plotted["selected_receiver_scores"], expected)
    summary = plotted["summary"].set_index("receiver")
    assert summary.loc["A", "sensor_score"] == 0.2
    assert summary.loc["B", "sensor_score"] == 0.4
    pd.testing.assert_frame_equal(result.receiver_scores, original)


@pytest.mark.parametrize("metabolite", [None, "Shared"])
def test_name_or_gene_alone_cannot_select_metabolite_context(receiver_context, metabolite):
    result, adata = receiver_context
    with pytest.raises(ValueError, match="provide hmdb_id"):
        plot_receptor_expression_violin(result, adata, receptor_gene="R", metabolite=metabolite)


@pytest.mark.parametrize("hmdb", [" ", np.nan, "null"])
def test_missing_explicit_hmdb_does_not_fall_back_to_a_name(receiver_context, hmdb):
    result, adata = receiver_context
    with pytest.raises(ValueError, match="hmdb_id must not be empty or missing"):
        plot_receptor_expression_violin(result, adata, receptor_gene="R", metabolite="Shared", hmdb_id=hmdb)


def test_unknown_hmdb_does_not_fall_back_to_a_name(receiver_context):
    result, adata = receiver_context
    with pytest.raises(ValueError, match="HMDB0000009.*not found"):
        plot_receptor_expression_violin(result, adata, receptor_gene="R", metabolite="Shared", hmdb_id="HMDB0000009")


def test_hmdb_match_does_not_replace_an_unknown_receptor(receiver_context):
    result, adata = receiver_context
    with pytest.raises(ValueError, match="MISSING.*not found"):
        plot_receptor_expression_violin(result, adata, receptor_gene="MISSING", hmdb_id=HMDB1)


def test_explicit_hmdb_requires_an_hmdb_column(receiver_context):
    result, adata = receiver_context
    result.receiver_scores = result.receiver_scores.drop(columns="hmdb_id")
    with pytest.raises(KeyError, match="hmdb_id"):
        plot_receptor_expression_violin(result, adata, receptor_gene="R", metabolite="Shared", hmdb_id=HMDB1)


def test_omitted_context_requires_an_hmdb(receiver_context):
    result, adata = receiver_context
    with pytest.raises(ValueError, match="provide hmdb_id"):
        plot_receptor_expression_violin(result, adata, receptor_gene="R")
