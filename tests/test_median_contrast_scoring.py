import numpy as np
import pandas as pd
from anndata import AnnData

import cellmesh.core as core
from cellmesh.core import _make_cell_mesh_events
from cellmesh import run_cell_mesh
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


def test_permutation_uses_only_originally_eligible_cell_types(monkeypatch):
    adata = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B", "Small"]}),
        var=pd.DataFrame(index=["G"]),
    )
    obs_events = pd.DataFrame(
        [
            {
                "sender": "A",
                "receiver": "B",
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
                "cell_mesh_score": 0.5,
            }
        ]
    )
    seen = []

    def fake_compute(adata_perm, enzyme_prior, sensor_prior, celltype_col, **kwargs):
        seen.append(
            {
                "n_obs": adata_perm.n_obs,
                "original_labels": set(adata_perm.obs["cell_type"].astype(str)),
                "perm_labels": set(adata_perm.obs[celltype_col].astype(str)),
            }
        )
        return pd.DataFrame(), pd.DataFrame(), {"cell_counts": pd.Series({"A": 2, "B": 2})}

    def fake_make(sender_scores, receiver_scores, allow_self, cell_counts):
        return obs_events.copy()

    monkeypatch.setattr(core, "_compute_availability_scores", fake_compute)
    monkeypatch.setattr(core, "_make_cell_mesh_events", fake_make)

    result = core._empirical_pvalues_by_sensor_type(
        obs_events,
        adata=adata,
        cell_type_key="cell_type",
        sample_key=None,
        layer=None,
        enzyme_prior=pd.DataFrame(),
        sensor_prior=pd.DataFrame(),
        n_perms=3,
        random_state=0,
        min_expr_frac=None,
        allow_self=True,
        availability_kwargs={"eps_num": 1e-12, "min_cells": 2},
    )

    assert len(seen) == 3
    assert all(item["n_obs"] == 4 for item in seen)
    assert all(item["original_labels"] == {"A", "B"} for item in seen)
    assert all(item["perm_labels"] <= {"A", "B"} for item in seen)
    assert result["perm_pvalue"].between(0, 1).all()


def _sample_mode_case():
    adata = AnnData(
        X=np.array(
            [
                [1.0, 1.0],
                [1.0, 1.0],
                [1.0, 1.0],
                [1.0, 1.0],
                [1.0, 1.0],
                [1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "sample": ["S1", "S1", "S1", "S2", "S2", "S2"],
                "cell_type": ["A", "A", "B", "A", "B", "B"],
            }
        ),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "gene": "PROD",
                "role": "production",
            }
        ]
    )
    sensor = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
            }
        ]
    )
    return adata, enzyme, sensor


def test_pooled_mode_min_cells_filters_by_total_cell_type_count():
    adata, enzyme, sensor = _sample_mode_case()
    result = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        cell_type_key="cell_type",
        sample_key="sample",
        min_cells=3,
        min_expr_frac=None,
        n_perms=0,
    )

    assert result.parameters["sample_mode"] == "pooled_stratified"
    assert set(result.sender_scores.columns) == {"A", "B"}
    assert set(result.receiver_scores["receiver"]) == {"A", "B"}


def test_sample_aware_min_cells_filters_each_sample_cell_type_unit():
    adata, enzyme, sensor = _sample_mode_case()
    result = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        cell_type_key="cell_type",
        sample_key="sample",
        sample_mode="sample_aware",
        min_cells=2,
        min_expr_frac=None,
        n_perms=0,
    )

    validation = result.availability_results["sample_validation"]
    flags = validation.set_index(["sample", "cell_type"])["eligible_in_sample"].to_dict()
    assert flags[("S1", "A")] is True
    assert flags[("S1", "B")] is False
    assert flags[("S2", "A")] is False
    assert flags[("S2", "B")] is True

    sample_sender = result.availability_results["sample_sender_scores"]
    assert pd.isna(sample_sender.loc[("S1", "Met", "HMDB0000001"), "B"])
    assert pd.isna(sample_sender.loc[("S2", "Met", "HMDB0000001"), "A"])
    assert not pd.isna(sample_sender.loc[("S1", "Met", "HMDB0000001"), "A"])
    assert not pd.isna(sample_sender.loc[("S2", "Met", "HMDB0000001"), "B"])

    events = result.events
    for col in [
        "event_score_median",
        "event_score_iqr",
        "n_samples_coobserved",
        "n_samples_positive",
        "event_prevalence",
        "inference_mode",
    ]:
        assert col in events.columns
    assert set(events["inference_mode"]) == {"sample_aware"}
    assert result.sample_validation is result.availability_results["sample_validation"]
    assert result.sample_sender_scores is result.availability_results["sample_sender_scores"]
    assert result.sample_receiver_scores is result.availability_results["sample_receiver_scores"]
    assert result.sample_events is result.availability_results["sample_events"]


def test_sample_aware_missing_sender_or_receiver_event_scores_are_na():
    adata = AnnData(
        X=np.ones((8, 2)),
        obs=pd.DataFrame(
            {
                "sample": ["S1", "S1", "S1", "S1", "S2", "S2", "S2", "S2"],
                "cell_type": ["A", "A", "B", "B", "A", "A", "Small", "Small"],
            }
        ),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme = pd.DataFrame(
        [{"metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": "PROD", "role": "production"}]
    )
    sensor = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
            }
        ]
    )

    result = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        cell_type_key="cell_type",
        sample_key="sample",
        sample_mode="sample_aware",
        min_cells=2,
        min_expr_frac=None,
        n_perms=0,
    )

    sample_scores = result.availability_results["sample_event_scores"]
    assert pd.isna(sample_scores.loc[("A", "B", "Met", "HMDB0000001", "SENSOR", "Transporter"), "S2"])
    assert pd.isna(sample_scores.loc[("A", "Small", "Met", "HMDB0000001", "SENSOR", "Transporter"), "S1"])
    assert not pd.isna(sample_scores.loc[("A", "Small", "Met", "HMDB0000001", "SENSOR", "Transporter"), "S2"])

    row = result.events[result.events["n_samples_coobserved"] == 1].iloc[0]
    assert row["n_samples_coobserved"] == 1
    assert pd.isna(row["event_score_iqr"])


def test_sample_aware_permutation_uses_fixed_within_sample_analysis_structure(monkeypatch):
    adata = AnnData(
        X=np.ones((10, 1)),
        obs=pd.DataFrame(
            {
                "sample": ["S1", "S1", "S1", "S1", "S1", "S2", "S2", "S2", "S2", "S2"],
                "cell_type": ["A", "A", "B", "B", "Small", "A", "A", "C", "C", "Drop"],
            }
        ),
        var=pd.DataFrame(index=["G"]),
    )
    obs_events = pd.DataFrame(
        [
            {
                "sender": "A",
                "receiver": "B",
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
                "cell_mesh_score": 0.5,
                "event_score_median": 0.5,
            },
            {
                "sender": "A",
                "receiver": "C",
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "sensor_gene": "SENSOR",
                "sensor_type": "Transporter",
                "cell_mesh_score": 0.25,
                "event_score_median": 0.25,
            },
        ]
    )
    seen = []

    def fake_compute(adata_perm, enzyme_prior, sensor_prior, cell_type_key, sample_key, **kwargs):
        counts = (
            adata_perm.obs.groupby([sample_key, cell_type_key], observed=True)
            .size()
            .astype(int)
            .to_dict()
        )
        seen.append(
            {
                "n_obs": adata_perm.n_obs,
                "original_labels": set(adata_perm.obs["cell_type"].astype(str)),
                "sample_labels": {
                    sample: set(frame[cell_type_key].astype(str))
                    for sample, frame in adata_perm.obs.groupby(sample_key, observed=True)
                },
                "counts": counts,
            }
        )
        if len(seen) == 1:
            return pd.DataFrame(), pd.DataFrame(), obs_events.iloc[[0]].assign(event_score_median=0.8), {}
        if len(seen) == 2:
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(columns=obs_events.columns), {}
        return pd.DataFrame(), pd.DataFrame(), obs_events.iloc[[1]].assign(event_score_median=0.9), {}

    monkeypatch.setattr(core, "_compute_sample_aware_scores", fake_compute)

    result = core._empirical_pvalues_by_sensor_type(
        obs_events,
        adata=adata,
        cell_type_key="cell_type",
        sample_key="sample",
        layer=None,
        enzyme_prior=pd.DataFrame(),
        sensor_prior=pd.DataFrame(),
        n_perms=3,
        random_state=0,
        min_expr_frac=None,
        allow_self=True,
        availability_kwargs={"eps_num": 1e-12, "min_cells": 2},
        sample_mode="sample_aware",
    )

    assert len(seen) == 3
    assert all(item["n_obs"] == 8 for item in seen)
    assert all(item["original_labels"] == {"A", "B", "C"} for item in seen)
    assert all(item["sample_labels"]["S1"] <= {"A", "B"} for item in seen)
    assert all(item["sample_labels"]["S2"] <= {"A", "C"} for item in seen)
    for item in seen:
        assert item["counts"] == {("S1", "A"): 2, ("S1", "B"): 2, ("S2", "A"): 2, ("S2", "C"): 2}

    null_scores = result.attrs["sample_aware_null_scores"]
    assert null_scores.shape == (2, 3)
    assert not null_scores.isna().any().any()
    assert result["perm_pvalue"].between(0, 1).all()
    assert set(result["permutation_mode"]) == {"within_sample_label_shuffle"}
    pd.testing.assert_series_equal(result["fdr"], result["fdr_sensor_type"], check_names=False)


def test_sample_aware_permutation_reports_pvalues_for_single_sample_events():
    adata, enzyme, sensor = _sample_mode_case()
    result = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        cell_type_key="cell_type",
        sample_key="sample",
        sample_mode="sample_aware",
        min_cells=2,
        min_expr_frac=None,
        n_perms=2,
        random_state=0,
    )

    row = result.events[result.events["n_samples_coobserved"] == 1].iloc[0]
    assert row["n_samples_coobserved"] == 1
    assert not pd.isna(row["perm_pvalue"])
    assert row["permutation_mode"] == "within_sample_label_shuffle"
    assert row["inference_mode"] == "sample_aware"
    assert "fdr_global" in result.events.columns
    assert "fdr_sensor_type" in result.events.columns
    pd.testing.assert_series_equal(result.events["fdr"], result.events["fdr_sensor_type"], check_names=False)
    assert result.events.attrs["sample_aware_null_scores"].shape[1] == 2


def test_sample_aware_requires_valid_sample_key():
    adata, enzyme, sensor = _sample_mode_case()
    with np.testing.assert_raises(ValueError):
        run_cell_mesh(
            adata,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            sample_mode="sample_aware",
            min_cells=2,
        )

    bad = adata.copy()
    bad.obs["sample"] = ["S1", "S1", "", "S2", "S2", "S2"]
    with np.testing.assert_raises(ValueError):
        run_cell_mesh(
            bad,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            sample_key="sample",
            sample_mode="sample_aware",
            min_cells=2,
        )
