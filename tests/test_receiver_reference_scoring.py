import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import run_cell_mesh
from cellmesh.score import compute_sensor_scores


def _sensor_prior():
    return pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
            }
        ]
    )


def _receiver_adata(values):
    values = np.asarray(values, dtype=float)
    return AnnData(
        X=values[:, None],
        obs=pd.DataFrame(
            {"cell_type": [f"T{index}" for index in range(len(values))]}
        ),
        var=pd.DataFrame(index=["SENSOR"]),
    )


def test_receiver_positive_reference_saturation_defaults_to_positive_median():
    adata = _receiver_adata([0.0, 1.0, 2.0, 9.0])

    scores = compute_sensor_scores(
        adata,
        _sensor_prior(),
        min_cells=1,
        min_expr_frac=None,
    ).set_index("receiver")

    # Only [1, 2, 9] enter the reference, whose median is 2.
    assert scores.loc["T0", "sensor_score"] == 0.0
    assert scores.loc["T1", "sensor_score"] == pytest.approx(1.0 / 3.0)
    assert scores.loc["T2", "sensor_score"] == pytest.approx(0.5)
    assert scores.loc["T3", "sensor_score"] == pytest.approx(9.0 / 11.0)


def test_receiver_positive_reference_saturation_supports_positive_mean():
    adata = _receiver_adata([0.0, 1.0, 2.0, 9.0])

    scores = compute_sensor_scores(
        adata,
        _sensor_prior(),
        min_cells=1,
        min_expr_frac=None,
        receiver_reference="mean",
    ).set_index("receiver")

    # The positive arithmetic mean is 4.
    assert scores.loc["T0", "sensor_score"] == 0.0
    assert scores.loc["T1", "sensor_score"] == pytest.approx(1.0 / 5.0)
    assert scores.loc["T2", "sensor_score"] == pytest.approx(1.0 / 3.0)
    assert scores.loc["T3", "sensor_score"] == pytest.approx(9.0 / 13.0)


def test_receiver_all_zero_expression_returns_zero_scores():
    scores = compute_sensor_scores(
        _receiver_adata([0.0, 0.0, 0.0]),
        _sensor_prior(),
        min_cells=1,
        min_expr_frac=None,
    )

    assert len(scores) == 3
    np.testing.assert_array_equal(
        scores["sensor_score"].to_numpy(dtype=float),
        np.zeros(3),
    )


def test_receiver_single_positive_expression_maps_to_half():
    scores = compute_sensor_scores(
        _receiver_adata([0.0, 0.0, 5.0]),
        _sensor_prior(),
        min_cells=1,
        min_expr_frac=None,
    ).set_index("receiver")

    assert scores.loc["T0", "sensor_score"] == 0.0
    assert scores.loc["T1", "sensor_score"] == 0.0
    assert scores.loc["T2", "sensor_score"] == pytest.approx(0.5)


@pytest.mark.parametrize("invalid_expression", [-1.0, np.nan, np.inf])
def test_receiver_rejects_invalid_pseudobulk_expression(invalid_expression):
    with pytest.raises(
        ValueError,
        match="receiver pseudobulk expression must be finite and non-negative",
    ):
        compute_sensor_scores(
            _receiver_adata([1.0, invalid_expression]),
            _sensor_prior(),
            min_cells=1,
            min_expr_frac=None,
        )


@pytest.mark.parametrize(
    ("invalid", "error_type"),
    [(None, TypeError), ("", ValueError), ("average", ValueError), (1, TypeError)],
)
def test_compute_sensor_scores_validates_receiver_reference(invalid, error_type):
    with pytest.raises(error_type, match="receiver_reference"):
        compute_sensor_scores(
            _receiver_adata([1.0, 2.0]),
            _sensor_prior(),
            min_cells=1,
            min_expr_frac=None,
            receiver_reference=invalid,
        )


def test_pooled_run_propagates_receiver_reference_and_uses_requested_mean():
    frame = pd.DataFrame(
        {
            "cell_type": ["T0", "T1", "T2", "T3"],
            "PROD": [1.0, 1.0, 1.0, 1.0],
            "SENSOR": [0.0, 1.0, 2.0, 9.0],
        }
    )
    adata = AnnData(
        X=frame[["PROD", "SENSOR"]].to_numpy(dtype=float),
        obs=frame[["cell_type"]].copy(),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "gene": "PROD",
                "role": "production",
                "reaction": "prod",
            }
        ]
    )

    result = run_cell_mesh(
        adata,
        enzyme,
        _sensor_prior(),
        min_cells=1,
        min_expr_frac=None,
        receiver_reference=" MEAN ",
        n_perms=0,
    )

    assert result.parameters["receiver_reference"] == "mean"
    scores = result.receiver_scores.set_index("receiver")["sensor_score"]
    assert scores.loc["T0"] == 0.0
    assert scores.loc["T1"] == pytest.approx(1.0 / 5.0)
    assert scores.loc["T2"] == pytest.approx(1.0 / 3.0)
    assert scores.loc["T3"] == pytest.approx(9.0 / 13.0)


def test_sample_aware_run_computes_requested_reference_within_each_sample():
    rows = []
    for sample, scale in [("S1", 1.0), ("S2", 2.0)]:
        for cell_type, sensor_expression in zip(
            ["T1", "T2", "T3"],
            [1.0, 2.0, 9.0],
        ):
            rows.append(
                {
                    "sample": sample,
                    "cell_type": cell_type,
                    "PROD": 1.0,
                    "SENSOR": scale * sensor_expression,
                }
            )
    frame = pd.DataFrame(rows)
    adata = AnnData(
        X=frame[["PROD", "SENSOR"]].to_numpy(dtype=float),
        obs=frame[["sample", "cell_type"]].copy(),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "gene": "PROD",
                "role": "production",
                "reaction": "prod",
            }
        ]
    )

    result = run_cell_mesh(
        adata,
        enzyme,
        _sensor_prior(),
        sample_key="sample",
        sample_mode="sample_aware",
        min_cells=1,
        min_expr_frac=None,
        receiver_reference="mean",
        n_perms=0,
    )

    assert result.parameters["receiver_reference"] == "mean"
    assert result.parameters["sample_aware_component_aggregation"] == "median"
    assert result.parameters["sample_aware_event_score"] == "median_of_sample_level_cell_mesh_scores"
    assert {
        "metabolite_availability_median",
        "sensor_score_median",
        "sensor_expr_frac_median",
    }.issubset(result.events.columns)
    assert {
        "metabolite_availability",
        "sensor_expr_frac",
    }.isdisjoint(result.events.columns)
    assert np.allclose(
        result.events["cell_mesh_score"],
        result.events["event_score_median"],
    )
    expected_sender = result.sample_sender_scores.groupby(
        level=["metabolite", "hmdb_id"]
    ).median()
    pd.testing.assert_frame_equal(result.sender_scores, expected_sender)
    expected_receiver = (
        result.sample_receiver_scores.groupby(
            ["metabolite", "hmdb_id", "sensor_gene", "sensor_type", "receiver"]
        )["sensor_score"]
        .median()
        .sort_index()
    )
    actual_receiver = (
        result.receiver_scores.set_index(
            ["metabolite", "hmdb_id", "sensor_gene", "sensor_type", "receiver"]
        )["sensor_score"]
        .sort_index()
    )
    pd.testing.assert_series_equal(actual_receiver, expected_receiver)
    for sample_result in result.availability_results[
        "availability_by_sample"
    ].values():
        assert sample_result["receiver_reference"] == "mean"

    sample_scores = result.sample_receiver_scores.set_index(
        ["sample", "receiver"]
    )["sensor_score"]
    expected = {
        "T1": 1.0 / 5.0,
        "T2": 1.0 / 3.0,
        "T3": 9.0 / 13.0,
    }
    for sample in ("S1", "S2"):
        for receiver, score in expected.items():
            assert sample_scores.loc[(sample, receiver)] == pytest.approx(score)
