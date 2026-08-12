import inspect

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import cellmesh.core as core
import cellmesh
from cellmesh.core import _make_cell_mesh_events
from cellmesh import run_cell_mesh
from cellmesh.score import (
    _score_PCE_by_reference,
    compute_metabolite_availability,
    compute_sensor_scores,
)
from cellmesh.config import METABOLITE_AVAILABILITY_DEFAULTS


def test_eps_num_is_removed_from_scoring_interfaces_and_defaults():
    for function in [
        run_cell_mesh,
        compute_metabolite_availability,
        compute_sensor_scores,
        _score_PCE_by_reference,
    ]:
        assert "eps_num" not in inspect.signature(function).parameters
    assert "eps_num" not in METABOLITE_AVAILABILITY_DEFAULTS
    assert not hasattr(cellmesh, "bounded_median_contrast")


def test_sender_positive_reference_supports_mean_default_and_median_option():
    index = pd.MultiIndex.from_tuples([("Met", "HMDB0000001")])
    columns = ["A", "B", "C"]
    P = pd.DataFrame([[1.0, 2.0, 9.0]], index=index, columns=columns)
    C = pd.DataFrame([[10.0, 10.0, 10.0]], index=index, columns=columns)
    E = pd.DataFrame([[0.0, 1.0, 9.0]], index=index, columns=columns)

    result = _score_PCE_by_reference(P, C, E)
    np.testing.assert_allclose(
        result["P_score"].loc[index[0]],
        [1.0 / 5.0, 1.0 / 3.0, 9.0 / 13.0],
    )
    assert result["P_ref"].loc[index[0]] == 4.0
    np.testing.assert_allclose(result["C_score"].loc[index[0]], 0.5)

    median_result = _score_PCE_by_reference(
        P,
        C,
        E,
        pce_reference="median",
    )
    np.testing.assert_allclose(
        median_result["P_score"].loc[index[0]],
        [1.0 / 3.0, 0.5, 9.0 / 11.0],
    )
    assert median_result["P_ref"].loc[index[0]] == 2.0
    np.testing.assert_allclose(
        result["E_score"].loc[index[0]],
        [0.0, 1.0 / 6.0, 9.0 / 14.0],
    )
    assert result["C_ref"].loc[index[0]] == 10.0
    assert result["E_ref"].loc[index[0]] == 5.0


def test_positive_reference_saturation_distinguishes_values_when_all_cell_median_is_zero():
    index = pd.MultiIndex.from_tuples([("Met", "HMDB0000001")])
    columns = ["A", "B", "C", "D", "E"]
    P = pd.DataFrame([[0.0, 0.0, 0.0, 1.0, 10.0]], index=index, columns=columns)
    zeros = pd.DataFrame(0.0, index=index, columns=columns)
    result = _score_PCE_by_reference(P, zeros, zeros)

    assert result["P_ref"].loc[index[0]] == 5.5
    np.testing.assert_allclose(
        result["P_score"].loc[index[0]],
        [0.0, 0.0, 0.0, 1.0 / 6.5, 10.0 / 15.5],
    )


def test_positive_reference_saturation_single_positive_and_nonfinite_rejection():
    index = pd.MultiIndex.from_tuples([("Met", "HMDB0000001")])
    columns = ["A", "B", "C"]
    zeros = pd.DataFrame(0.0, index=index, columns=columns)
    single_positive = pd.DataFrame(
        [[0.0, 0.0, 4.0]],
        index=index,
        columns=columns,
    )
    result = _score_PCE_by_reference(single_positive, zeros, zeros)
    assert result["P_ref"].loc[index[0]] == 4.0
    np.testing.assert_allclose(
        result["P_score"].loc[index[0]],
        [0.0, 0.0, 0.5],
    )
    for bad in (np.nan, np.inf):
        invalid = single_positive.copy()
        invalid.loc[index[0], "B"] = bad
        with pytest.raises(ValueError, match="finite and non-negative"):
            _score_PCE_by_reference(invalid, zeros, zeros)


def test_receiver_positive_reference_saturation_optional_gate_and_cell_counts():
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
    assert set(scores["receiver"]) == {"A", "B", "C", "ineligible"}
    assert set(scores["receiver_n_cells"]) == {1, 2}
    scores = scores.set_index("receiver")
    # min_cells is QC-only: the one-cell type enters the positive receiver
    # reference. Means [1, 2, 4, 99] have median reference 3.
    assert scores.loc["A", "sensor_score"] == pytest.approx(1.0 / 4.0)
    assert scores.loc["B", "sensor_score"] == pytest.approx(2.0 / 5.0)
    assert scores.loc["C", "sensor_score"] == pytest.approx(4.0 / 7.0)
    assert scores.loc["ineligible", "sensor_score"] == pytest.approx(99.0 / 102.0)
    assert scores["receiver_passes_min_cells"].to_dict() == {
        "A": True,
        "B": True,
        "C": True,
        "ineligible": False,
    }

    gated = compute_sensor_scores(adata, prior, min_cells=2, min_expr_frac=0.75)
    c_row = gated.set_index("receiver").loc["C"]
    assert c_row["sensor_expr_frac"] == 0.5
    assert c_row["sensor_score"] == 0.0


def test_receiver_score_uses_per_celltype_mean_sensor_expression():
    adata = AnnData(
        X=np.array([[2.0]] * 6 + [[4.0]] * 2),
        obs=pd.DataFrame(
            {"cell_type": ["Abundant"] * 6 + ["Rare"] * 2}
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

    scores = compute_sensor_scores(
        adata,
        prior,
        min_cells=2,
        min_expr_frac=None,
    ).set_index("receiver")

    # Receiver abundance remains available as metadata, but the score follows
    # the within-cell-type mean expression rather than population receptor mass.
    assert scores.loc["Abundant", "receiver_cell_fraction"] == 0.75
    assert scores.loc["Rare", "receiver_cell_fraction"] == 0.25
    # Equal weighting of cell types gives ref=median([2, 4])=3. Repeating
    # "Abundant" cells does not make its six cells count six times in the ref.
    assert scores.loc["Abundant", "sensor_score"] == pytest.approx(2.0 / 5.0)
    assert scores.loc["Rare", "sensor_score"] == pytest.approx(4.0 / 7.0)


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


def test_pooled_permutation_is_min_cells_invariant():
    adata, enzyme, sensor = _sample_mode_case()
    common = dict(
        adata=adata,
        enzyme_metabolite=enzyme,
        metabolite_sensor=sensor,
        cell_type_key="cell_type",
        sample_key="sample",
        sample_mode="pooled_stratified",
        n_perms=3,
        random_state=0,
    )
    low = run_cell_mesh(**common, min_cells=1)
    high = run_cell_mesh(**common, min_cells=100)
    keys = core.EVENT_KEY_COLUMNS
    left = low.events.sort_values(keys).reset_index(drop=True)
    right = high.events.sort_values(keys).reset_index(drop=True)
    np.testing.assert_allclose(
        left[["perm_pvalue", "fdr_global", "fdr_sensor_type"]],
        right[["perm_pvalue", "fdr_global", "fdr_sensor_type"]],
    )


def test_pooled_zero_observed_score_has_unit_permutation_pvalue():
    adata, enzyme, sensor = _sample_mode_case()
    baseline = run_cell_mesh(
        adata, enzyme, sensor, cell_type_key="cell_type", n_perms=0, min_cells=1
    )
    obs_events = baseline.events.iloc[[0]].copy()
    obs_events["cell_mesh_score"] = 0.0
    enzyme_prior, sensor_prior = core.validate_priors(
        enzyme, sensor, adata.var_names
    )
    result = core._empirical_pvalues_by_sensor_type(
        obs_events, adata=adata, cell_type_key="cell_type", sample_key=None,
        layer=None, enzyme_prior=enzyme_prior, sensor_prior=sensor_prior,
        n_perms=3, random_state=0, min_expr_frac=None, allow_self=True,
        availability_kwargs={"min_cells": 1},
    )

    assert result.loc[0, "perm_pvalue"] == 1.0
    assert result.loc[0, "fdr_global"] == 1.0
    assert result.loc[0, "fdr_sensor_type"] == 1.0


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
                "reaction": "prod",
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


def test_pooled_mode_min_cells_records_qc_without_filtering_cell_types():
    adata, enzyme, sensor = _sample_mode_case()
    result = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        cell_type_key="cell_type",
        sample_key="sample",
        min_cells=3,
        min_expr_frac=None,
        sender_abundance_exponent=0.5,
        n_perms=0,
    )

    assert result.parameters["sample_mode"] == "pooled_stratified"
    assert result.parameters["sender_abundance_exponent"] == 0.5
    assert result.parameters["pce_reference"] == "mean"
    assert result.parameters["receiver_reference"] == "median"
    assert result.parameters["export_weight"] == pytest.approx(0.2)
    assert result.availability_results["pce_reference"] == "mean"
    pd.testing.assert_series_equal(
        result.availability_results["sender_abundance_weights"],
        result.availability_results["cell_fractions"]
        .pow(0.5)
        .rename("sender_abundance_weight"),
    )
    assert set(result.sender_scores.columns) == {"A", "B"}
    assert set(result.receiver_scores["receiver"]) == {"A", "B"}
    assert result.celltype_qc["passes_min_cells"].all()
    assert result.receiver_scores["receiver_passes_min_cells"].all()
    assert result.events["sender_passes_min_cells"].all()
    assert result.events["receiver_passes_min_cells"].all()
    assert result.events["passes_min_cells"].all()


def test_prior_gene_unavailable_status_is_consistent_across_public_apis():
    adata = AnnData(
        X=np.array(
            [
                [1.0, 1.0],
                [1.0, 1.0],
                [2.0, 2.0],
                [2.0, 2.0],
            ]
        ),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}),
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
            },
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "gene": "CONS_NOT_MEASURED",
                "role": "degradation",
                "reaction": "cons",
            },
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
    idx = ("Met", "HMDB0000001")

    low_level = core.compute_metabolite_availability(
        adata,
        enzyme,
        min_cells=2,
    )
    high_level = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        cell_type_key="cell_type",
        min_cells=2,
        n_perms=0,
    ).availability_results

    for result in (low_level, high_level):
        assert (
            result["metadata"].loc[idx, "consumption_status"]
            == "prior_gene_unavailable"
        )
        assert "has_substrate" not in result["metadata"]
        assert "has_usable_substrate" not in result["metadata"]
        pd.testing.assert_frame_equal(result["base_availability"], result["P_score"])
        pd.testing.assert_frame_equal(
            result["availability"], result["base_availability"] * 0.9
        )


def test_all_zero_product_keeps_receiver_schema_and_intermediates():
    adata = AnnData(
        X=np.array(
            [
                [0.0, 1.0],
                [0.0, 1.0],
                [0.0, 4.0],
                [0.0, 4.0],
            ]
        ),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}),
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
        min_cells=2,
        n_perms=0,
    )

    assert result.events.empty
    assert result.sender_scores.shape == (0, 2)
    assert result.sender_scores.columns.tolist() == ["A", "B"]
    assert set(result.receiver_scores["receiver"]) == {"A", "B"}
    assert result.availability_results["pseudobulk"].shape == (2, 2)
    assert result.availability_results["expr_frac"].shape == (2, 2)
    assert not result.availability_results["reaction_genes"].empty


@pytest.mark.parametrize(
    ("invalid", "error_type"),
    [
        (-0.1, ValueError),
        (np.nan, ValueError),
        (np.inf, ValueError),
        ("invalid", TypeError),
        (True, TypeError),
    ],
)
def test_run_cell_mesh_validates_sender_abundance_exponent_before_scoring(
    invalid, error_type
):
    adata, enzyme, sensor = _sample_mode_case()
    with pytest.raises(error_type, match="sender_abundance_exponent"):
        run_cell_mesh(
            adata,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            sample_key="sample",
            sample_mode="sample_aware",
            min_cells=100,
            sender_abundance_exponent=invalid,
            n_perms=0,
        )


@pytest.mark.parametrize(
    ("invalid", "error_type"),
    [(None, TypeError), ("", ValueError), ("mode", ValueError), (1, TypeError)],
)
def test_run_cell_mesh_validates_pce_reference_before_scoring(invalid, error_type):
    adata, enzyme, sensor = _sample_mode_case()
    with pytest.raises(error_type, match="pce_reference"):
        run_cell_mesh(
            adata,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            sample_key="sample",
            sample_mode="sample_aware",
            min_cells=100,
            pce_reference=invalid,
            n_perms=0,
        )


@pytest.mark.parametrize(
    ("invalid", "error_type"),
    [
        (-0.1, ValueError),
        (1.1, ValueError),
        (np.nan, ValueError),
        (np.inf, ValueError),
        ("invalid", TypeError),
        (True, TypeError),
    ],
)
def test_run_cell_mesh_validates_export_weight_before_scoring(invalid, error_type):
    adata, enzyme, sensor = _sample_mode_case()
    with pytest.raises(error_type, match="export_weight"):
        run_cell_mesh(
            adata,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            sample_key="sample",
            sample_mode="sample_aware",
            min_cells=100,
            export_weight=invalid,
            n_perms=0,
        )


@pytest.mark.parametrize(
    ("invalid", "error_type"),
    [(None, TypeError), ("", ValueError), ("mode", ValueError), (1, TypeError)],
)
def test_run_cell_mesh_validates_receiver_reference_before_scoring(
    invalid, error_type
):
    adata, enzyme, sensor = _sample_mode_case()
    with pytest.raises(error_type, match="receiver_reference"):
        run_cell_mesh(
            adata,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            sample_key="sample",
            sample_mode="sample_aware",
            min_cells=100,
            receiver_reference=invalid,
            n_perms=0,
        )


@pytest.mark.parametrize(
    ("invalid", "error_type"),
    [(-1, ValueError), (1.5, TypeError), (True, TypeError), ("10", TypeError)],
)
def test_run_cell_mesh_validates_n_perms_before_scoring(invalid, error_type):
    adata, enzyme, sensor = _sample_mode_case()
    with pytest.raises(error_type, match="n_perms"):
        run_cell_mesh(
            adata,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            sample_key="sample",
            sample_mode="sample_aware",
            min_cells=1,
            n_perms=invalid,
        )


@pytest.mark.parametrize(
    ("invalid", "error_type"),
    [
        (-0.1, ValueError),
        (1.1, ValueError),
        (np.nan, ValueError),
        (np.inf, ValueError),
        (True, TypeError),
        ("0.5", TypeError),
    ],
)
def test_receiver_expression_fraction_gate_is_strictly_validated(
    invalid, error_type
):
    adata, enzyme, sensor = _sample_mode_case()
    with pytest.raises(error_type, match="min_expr_frac"):
        run_cell_mesh(
            adata,
            enzyme,
            sensor,
            cell_type_key="cell_type",
            min_cells=1,
            min_expr_frac=invalid,
            n_perms=0,
        )
    with pytest.raises(error_type, match="min_expr_frac"):
        compute_sensor_scores(
            adata,
            sensor,
            min_cells=1,
            min_expr_frac=invalid,
        )


def test_sample_aware_min_cells_marks_but_retains_each_sample_cell_type_unit():
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
    flags = validation.set_index(["sample", "cell_type"])["passes_min_cells"].to_dict()
    assert flags[("S1", "A")] is True
    assert flags[("S1", "B")] is False
    assert flags[("S2", "A")] is False
    assert flags[("S2", "B")] is True
    assert "eligible_in_sample" not in validation.columns
    fractions = validation.set_index(["sample", "cell_type"])["cell_fraction"]
    assert fractions.loc[("S1", "A")] == 2.0 / 3.0
    assert fractions.loc[("S1", "B")] == 1.0 / 3.0
    assert fractions.loc[("S2", "A")] == 1.0 / 3.0
    assert fractions.loc[("S2", "B")] == 2.0 / 3.0
    assert (
        result.availability_results["availability_by_sample"]["S1"][
            "cell_fractions"
        ].loc["A"]
        == 2.0 / 3.0
    )
    assert (
        result.availability_results["availability_by_sample"]["S2"][
            "cell_fractions"
        ].loc["B"]
        == 2.0 / 3.0
    )

    sample_sender = result.availability_results["sample_sender_scores"]
    assert not sample_sender.isna().any().any()

    sample_receiver = result.availability_results["sample_receiver_scores"]
    receiver_flags = sample_receiver.set_index(["sample", "receiver"])[
        "receiver_passes_min_cells"
    ].to_dict()
    assert receiver_flags == {
        ("S1", "A"): True,
        ("S1", "B"): False,
        ("S2", "A"): False,
        ("S2", "B"): True,
    }

    sample_events = result.availability_results["sample_events"]
    assert {
        "sender_passes_min_cells",
        "receiver_passes_min_cells",
        "passes_min_cells",
    }.issubset(sample_events.columns)
    event_flags = sample_events.set_index(["sample", "sender", "receiver"])[
        "passes_min_cells"
    ].to_dict()
    assert event_flags[("S1", "A", "A")] is True
    assert event_flags[("S1", "A", "B")] is False
    assert event_flags[("S2", "B", "B")] is True
    assert event_flags[("S2", "B", "A")] is False

    events = result.events
    assert "n_observed_samples" not in events.columns
    assert "n_valid_samples" not in events.columns
    for col in [
        "event_score_median",
        "event_score_iqr",
        "n_samples_coobserved",
        "n_samples_positive",
        "event_prevalence",
        "n_samples_passing_min_cells",
        "min_cells_pass_prevalence",
        "passes_min_cells",
        "inference_mode",
    ]:
        assert col in events.columns
    event_qc = events.set_index(["sender", "receiver"])
    assert event_qc.loc[("A", "A"), "n_samples_passing_min_cells"] == 1
    assert event_qc.loc[("A", "A"), "min_cells_pass_prevalence"] == 0.5
    assert event_qc.loc[("A", "A"), "passes_min_cells"]
    assert event_qc.loc[("A", "B"), "n_samples_passing_min_cells"] == 0
    assert event_qc.loc[("A", "B"), "min_cells_pass_prevalence"] == 0.0
    assert not event_qc.loc[("A", "B"), "passes_min_cells"]
    assert set(events["inference_mode"]) == {"sample_aware"}
    assert result.sample_validation is result.availability_results["sample_validation"]
    assert result.sample_sender_scores is result.availability_results["sample_sender_scores"]
    assert result.sample_receiver_scores is result.availability_results["sample_receiver_scores"]
    assert result.sample_events is result.availability_results["sample_events"]


def test_sample_aware_applies_configured_reference_with_nondefault_exponent():
    rows = []
    for sample, cell_type, n_cells, prod, cons, sensor_expr in [
        ("S1", "A", 4, 4.0, 1.0, 4.0),
        ("S1", "B", 1, 1.0, 0.0, 1.0),
        ("S2", "A", 1, 2.0, 0.0, 1.0),
        ("S2", "B", 4, 4.0, 2.0, 4.0),
    ]:
        rows.extend(
            {
                "sample": sample,
                "cell_type": cell_type,
                "PROD": prod,
                "CONS": cons,
                "SENSOR": sensor_expr,
            }
            for _ in range(n_cells)
        )
    frame = pd.DataFrame(rows)
    adata = AnnData(
        X=frame[["PROD", "CONS", "SENSOR"]].to_numpy(dtype=float),
        obs=frame[["sample", "cell_type"]].copy(),
        var=pd.DataFrame(index=["PROD", "CONS", "SENSOR"]),
    )
    enzyme = pd.DataFrame(
        [
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "gene": "PROD",
                "role": "production",
                "reaction": "prod",
            },
            {
                "metabolite": "Met",
                "hmdb_id": "HMDB0000001",
                "gene": "CONS",
                "role": "degradation",
                "reaction": "cons",
            },
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

    result = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        cell_type_key="cell_type",
        sample_key="sample",
        sample_mode="sample_aware",
        min_cells=1,
        min_expr_frac=None,
        sender_abundance_exponent=0.5,
        pce_reference="mean",
        receiver_reference="mean",
        n_perms=0,
    )
    idx = ("Met", "HMDB0000001")
    assert result.parameters["receiver_reference"] == "mean"
    assert result.availability_results["export_weight"] == pytest.approx(0.2)

    for sample in ("S1", "S2"):
        sample_result = result.availability_results["availability_by_sample"][sample]
        assert sample_result["pce_reference"] == "mean"
        assert sample_result["receiver_reference"] == "mean"
        expected_weights = sample_result["cell_fractions"].pow(0.5)
        pd.testing.assert_series_equal(
            sample_result["sender_abundance_weights"],
            expected_weights.rename("sender_abundance_weight"),
        )
        for raw_key, score_key, ref_key in [
            ("P", "P_score", "P_ref"),
            ("C", "C_score", "C_ref"),
        ]:
            raw = sample_result[raw_key].loc[idx]
            positive = raw > 0.0
            reference = raw.loc[positive].mean()
            expected_score = pd.Series(0.0, index=raw.index)
            expected_score.loc[positive] = (
                raw.loc[positive] / (raw.loc[positive] + reference)
            )
            assert sample_result[ref_key].loc[idx] == pytest.approx(reference)
            pd.testing.assert_series_equal(
                sample_result[score_key].loc[idx],
                expected_score,
                check_names=False,
            )

        p_score = sample_result["P_score"].loc[idx]
        c_score = sample_result["C_score"].loc[idx]
        expected_base = (
            p_score.pow(2)
            .div((p_score + c_score).where((p_score + c_score) > 0.0))
            .fillna(0.0)
        )
        pd.testing.assert_series_equal(
            sample_result["base_availability"].loc[idx],
            expected_base,
            check_names=False,
        )
        expected_sender = expected_base * 0.9
        pd.testing.assert_series_equal(
            sample_result["availability"].loc[idx],
            expected_sender,
            check_names=False,
        )


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
        [{"metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": "PROD", "role": "production", "reaction": "prod"}]
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

    sample_scores = result.sample_events.pivot_table(
        index=["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "sensor_type"],
        columns="sample",
        values="cell_mesh_score",
        aggfunc="mean",
        dropna=False,
    )
    assert pd.isna(sample_scores.loc[("A", "B", "Met", "HMDB0000001", "SENSOR", "Transporter"), "S2"])
    assert pd.isna(sample_scores.loc[("A", "Small", "Met", "HMDB0000001", "SENSOR", "Transporter"), "S1"])
    assert not pd.isna(sample_scores.loc[("A", "Small", "Met", "HMDB0000001", "SENSOR", "Transporter"), "S2"])

    row = result.events[result.events["n_samples_coobserved"] == 1].iloc[0]
    assert row["n_samples_coobserved"] == 1
    assert pd.isna(row["event_score_iqr"])


def test_sample_aware_permutation_reports_pvalues_for_all_retained_events():
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
        store_null_scores=True,
    )

    assert not result.events.empty
    assert (result.events["n_samples_coobserved"] == 2).all()
    assert result.events["perm_pvalue"].notna().all()
    assert result.events["perm_pvalue"].between(0.0, 1.0).all()
    assert set(result.events["permutation_mode"]) == {
        "within_sample_label_shuffle"
    }
    assert set(result.events["inference_mode"]) == {"sample_aware"}
    assert "fdr_global" in result.events.columns
    assert "fdr_sensor_type" in result.events.columns
    assert result.events["fdr_global"].notna().all()
    assert result.events["fdr_sensor_type"].notna().all()
    assert "fdr" not in result.events.columns
    assert {
        "n_samples_passing_min_cells",
        "min_cells_pass_prevalence",
        "passes_min_cells",
    }.issubset(result.events.columns)
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
