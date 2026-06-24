import numpy as np
import pandas as pd
from anndata import AnnData

from cellmesh.core import _make_cell_mesh_events
from cellmesh.score import (
    _contrast_PCE,
    bounded_median_contrast,
    compute_sensor_scores,
)


def test_bounded_median_contrast_handles_nan_zero_and_clipping():
    np.testing.assert_array_equal(
        bounded_median_contrast(np.array([np.nan, np.nan])),
        np.array([0.0, 0.0]),
    )
    np.testing.assert_array_equal(
        bounded_median_contrast(np.array([0.0, 0.0, 0.0])),
        np.array([0.0, 0.0, 0.0]),
    )
    result = bounded_median_contrast(np.array([0.0, 2.0, 4.0, np.nan]))
    np.testing.assert_allclose(result, np.array([-1.0, 0.0, 1.0 / 3.0, 0.0]))
    assert np.all(result >= -1.0)
    assert np.all(result <= 1.0)


def test_sender_contrasts_require_production_and_use_neutral_missing_priors():
    index = pd.MultiIndex.from_tuples([("Met", "HMDB0000001")])
    columns = ["A", "B", "C"]
    P = pd.DataFrame([[1.0, 2.0, 4.0]], index=index, columns=columns)
    C = pd.DataFrame([[10.0, 10.0, 10.0]], index=index, columns=columns)
    E = pd.DataFrame([[0.0, 1.0, 9.0]], index=index, columns=columns)

    product_only = pd.DataFrame(
        [{"metabolite": "Met", "hmdb_id": "HMDB0000001", "direction": "product"}]
    )
    result = _contrast_PCE(P, C, E, product_only)
    p_plus = result["P_contrast"].clip(lower=0.0)

    # Missing E/C priors are neutral, and production at or below median is zero.
    np.testing.assert_allclose(p_plus.loc[index[0]], [0.0, 0.0, 1.0 / 3.0])

    all_priors = pd.DataFrame(
        [
            {"metabolite": "Met", "hmdb_id": "HMDB0000001", "direction": "product"},
            {"metabolite": "Met", "hmdb_id": "HMDB0000001", "direction": "substrate"},
            {"metabolite": "Met", "hmdb_id": "HMDB0000001", "direction": "exporter"},
        ]
    )
    result = _contrast_PCE(P, C, E, all_priors)
    c_plus = result["C_contrast"].clip(lower=0.0)
    e_plus = result["E_contrast"].clip(lower=0.0)
    sender = (
        p_plus.loc[index[0]]
        * (1.0 + e_plus.loc[index[0]])
        * (1.0 - c_plus.loc[index[0]])
    )

    assert sender["A"] == 0.0
    assert sender["B"] == 0.0
    assert c_plus.loc[index[0], "C"] == 0.0
    assert sender["C"] > p_plus.loc[index[0], "C"]


def test_receiver_contrast_optional_gate_and_cell_counts():
    adata = AnnData(
        X=np.array(
            [
                [1.0],
                [1.0],
                [2.0],
                [2.0],
                [8.0],
                [0.0],
                [99.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B", "C", "C", "ineligible"]}
        ),
        var=pd.DataFrame(index=["SENSOR"]),
    )
    prior = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
            }
        ]
    )

    scores = compute_sensor_scores(adata, prior, min_cells=2, min_expr_frac=None)
    assert set(scores["receiver"]) == {"A", "B", "C"}
    assert set(scores["receiver_n_cells"]) == {2}
    assert scores.set_index("receiver").loc["C", "sensor_score"] > 0

    gated = compute_sensor_scores(adata, prior, min_cells=2, min_expr_frac=0.75)
    c_row = gated.set_index("receiver").loc["C"]
    assert c_row["sensor_expr_frac"] == 0.5
    assert c_row["sensor_score"] == 0.0


def test_events_record_sender_and_receiver_cell_counts():
    index = pd.MultiIndex.from_tuples(
        [("Met", "HMDB0000001")], names=["metabolite", "hmdb_id"]
    )
    sender = pd.DataFrame([[0.25, 0.5]], index=index, columns=["A", "B"])
    receiver = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
                "receiver": "B",
                "sensor_score": 0.5,
                "sensor_expr_frac": 1.0,
                "receiver_n_cells": 12,
            }
        ]
    )
    events = _make_cell_mesh_events(
        sender,
        receiver,
        allow_self=False,
        cell_counts=pd.Series({"A": 10, "B": 12}),
    )
    assert events.loc[0, "sender_n_cells"] == 10
    assert events.loc[0, "receiver_n_cells"] == 12
    assert np.isclose(events.loc[0, "cell_mesh_score"], np.sqrt(0.25 * 0.5))
