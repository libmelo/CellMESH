import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import compute_metabolite_availability, run_cell_mesh
from cellmesh.score import compute_sensor_scores


EVENT_KEY = [
    "sender",
    "receiver",
    "metabolite",
    "hmdb_id",
    "sensor_gene",
    "sensor_type",
]
EVENT_QC_COLUMNS = [
    "sender_passes_min_cells",
    "receiver_passes_min_cells",
    "passes_min_cells",
]
SAMPLE_AGGREGATE_QC_COLUMNS = EVENT_QC_COLUMNS + [
    "n_samples_passing_min_cells",
    "min_cells_pass_prevalence",
]


def _priors():
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
    return enzyme, sensor


def _pooled_case():
    rows = []
    for cell_type, n_cells, prod, sensor in [
        ("A", 3, 2.0, 1.0),
        ("B", 2, 3.0, 4.0),
        ("Tiny", 1, 5.0, 2.0),
    ]:
        rows.extend(
            {
                "cell_type": cell_type,
                "PROD": prod,
                "SENSOR": sensor,
            }
            for _ in range(n_cells)
        )
    frame = pd.DataFrame(rows)
    adata = AnnData(
        X=frame[["PROD", "SENSOR"]].to_numpy(dtype=float),
        obs=frame[["cell_type"]].copy(),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme, sensor = _priors()
    return adata, enzyme, sensor


def _sample_aware_case():
    rows = []
    for sample, cell_type, n_cells, prod, sensor in [
        ("S1", "A", 2, 4.0, 1.0),
        ("S1", "B", 1, 2.0, 6.0),
        ("S2", "A", 1, 3.0, 2.0),
        ("S2", "B", 2, 5.0, 8.0),
    ]:
        rows.extend(
            {
                "sample": sample,
                "cell_type": cell_type,
                "PROD": prod,
                "SENSOR": sensor,
            }
            for _ in range(n_cells)
        )
    frame = pd.DataFrame(rows)
    adata = AnnData(
        X=frame[["PROD", "SENSOR"]].to_numpy(dtype=float),
        obs=frame[["sample", "cell_type"]].copy(),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme, sensor = _priors()
    return adata, enzyme, sensor


def _sorted(frame, keys):
    return frame.sort_values(keys).reset_index(drop=True)


def _assert_pooled_numerics_equal(left, right):
    pd.testing.assert_frame_equal(left.sender_scores, right.sender_scores)

    receiver_keys = [
        "metabolite",
        "hmdb_id",
        "sensor_gene",
        "sensor_type",
        "receiver",
    ]
    pd.testing.assert_frame_equal(
        _sorted(
            left.receiver_scores.drop(columns=["receiver_passes_min_cells"]),
            receiver_keys,
        ),
        _sorted(
            right.receiver_scores.drop(columns=["receiver_passes_min_cells"]),
            receiver_keys,
        ),
    )
    pd.testing.assert_frame_equal(
        _sorted(left.events.drop(columns=EVENT_QC_COLUMNS), EVENT_KEY),
        _sorted(right.events.drop(columns=EVENT_QC_COLUMNS), EVENT_KEY),
    )

    for key in [
        "pseudobulk",
        "expr_frac",
        "P",
        "C",
        "E",
        "P_score",
        "C_score",
        "E_score",
        "availability",
    ]:
        pd.testing.assert_frame_equal(
            left.availability_results[key],
            right.availability_results[key],
        )
    for key in [
        "cell_counts",
        "cell_fractions",
        "P_ref",
        "C_ref",
        "E_ref",
    ]:
        pd.testing.assert_series_equal(
            left.availability_results[key],
            right.availability_results[key],
        )


def test_pooled_min_cells_is_qc_only_for_scores_references_and_inference():
    adata, enzyme, sensor = _pooled_case()
    common = dict(
        cell_type_key="cell_type",
        min_expr_frac=None,
        n_perms=4,
        random_state=7,
    )
    low = run_cell_mesh(adata, enzyme, sensor, min_cells=1, **common)
    high = run_cell_mesh(adata, enzyme, sensor, min_cells=3, **common)

    _assert_pooled_numerics_equal(low, high)

    assert list(high.availability_results["pseudobulk"].index) == ["A", "B", "Tiny"]
    assert list(high.availability_results["cell_fractions"].index) == ["A", "B", "Tiny"]
    assert list(high.sender_scores.columns) == ["A", "B", "Tiny"]
    assert set(high.receiver_scores["receiver"]) == {"A", "B", "Tiny"}
    assert set(high.events["sender"]) == {"A", "B", "Tiny"}
    assert set(high.events["receiver"]) == {"A", "B", "Tiny"}
    assert high.availability_results["cell_fractions"].sum() == pytest.approx(1.0)
    assert high.availability_results["P_ref"].loc[("Met", "HMDB0000001")] == pytest.approx(
        17.0 / 18.0
    )

    assert low.celltype_qc["passes_min_cells"].all()
    expected_pass = {"A": True, "B": False, "Tiny": False}
    assert high.celltype_qc["passes_min_cells"].to_dict() == expected_pass
    assert high.availability_results["celltype_qc"][
        "passes_min_cells"
    ].to_dict() == expected_pass
    assert high.receiver_scores.set_index("receiver")[
        "receiver_passes_min_cells"
    ].to_dict() == expected_pass

    for _, row in high.events.iterrows():
        assert bool(row["sender_passes_min_cells"]) is expected_pass[row["sender"]]
        assert bool(row["receiver_passes_min_cells"]) is expected_pass[row["receiver"]]
        assert bool(row["passes_min_cells"]) is (
            expected_pass[row["sender"]] and expected_pass[row["receiver"]]
        )


def _assert_sample_aware_numerics_equal(left, right):
    pd.testing.assert_frame_equal(left.sender_scores, right.sender_scores)
    pd.testing.assert_frame_equal(left.sample_sender_scores, right.sample_sender_scores)
    assert "sample_event_scores" not in left.availability_results
    assert "sample_event_scores" not in right.availability_results

    receiver_keys = [
        "metabolite",
        "hmdb_id",
        "sensor_gene",
        "sensor_type",
        "receiver",
    ]
    receiver_qc = [
        "receiver_passes_min_cells",
        "n_samples_passing_min_cells",
        "min_cells_pass_prevalence",
    ]
    pd.testing.assert_frame_equal(
        _sorted(left.receiver_scores.drop(columns=receiver_qc), receiver_keys),
        _sorted(right.receiver_scores.drop(columns=receiver_qc), receiver_keys),
    )

    sample_receiver_keys = [
        "sample",
        "metabolite",
        "hmdb_id",
        "sensor_gene",
        "sensor_type",
        "receiver",
    ]
    pd.testing.assert_frame_equal(
        _sorted(
            left.sample_receiver_scores.drop(
                columns=["receiver_passes_min_cells"]
            ),
            sample_receiver_keys,
        ),
        _sorted(
            right.sample_receiver_scores.drop(
                columns=["receiver_passes_min_cells"]
            ),
            sample_receiver_keys,
        ),
    )
    pd.testing.assert_frame_equal(
        _sorted(
            left.sample_events.drop(columns=EVENT_QC_COLUMNS),
            ["sample"] + EVENT_KEY,
        ),
        _sorted(
            right.sample_events.drop(columns=EVENT_QC_COLUMNS),
            ["sample"] + EVENT_KEY,
        ),
    )
    pd.testing.assert_frame_equal(
        _sorted(left.events.drop(columns=SAMPLE_AGGREGATE_QC_COLUMNS), EVENT_KEY),
        _sorted(right.events.drop(columns=SAMPLE_AGGREGATE_QC_COLUMNS), EVENT_KEY),
    )
    pd.testing.assert_frame_equal(
        left.events.attrs["sample_aware_null_scores"],
        right.events.attrs["sample_aware_null_scores"],
    )

    base_validation = ["sample", "cell_type", "n_cells", "cell_fraction"]
    pd.testing.assert_frame_equal(
        _sorted(left.sample_validation[base_validation], ["sample", "cell_type"]),
        _sorted(right.sample_validation[base_validation], ["sample", "cell_type"]),
    )

    for sample in ("S1", "S2"):
        for key in [
            "pseudobulk",
            "expr_frac",
            "P",
            "C",
            "E",
            "P_score",
            "C_score",
            "E_score",
            "availability",
        ]:
            pd.testing.assert_frame_equal(
                left.availability_results["availability_by_sample"][sample][key],
                right.availability_results["availability_by_sample"][sample][key],
            )
        for key in ["cell_counts", "cell_fractions", "P_ref", "C_ref", "E_ref"]:
            pd.testing.assert_series_equal(
                left.availability_results["availability_by_sample"][sample][key],
                right.availability_results["availability_by_sample"][sample][key],
            )


def test_sample_aware_min_cells_is_qc_only_for_scores_and_inference():
    adata, enzyme, sensor = _sample_aware_case()
    common = dict(
        cell_type_key="cell_type",
        sample_key="sample",
        sample_mode="sample_aware",
        min_expr_frac=None,
        n_perms=4,
        random_state=11,
    )
    low = run_cell_mesh(adata, enzyme, sensor, min_cells=1, **common)
    high = run_cell_mesh(adata, enzyme, sensor, min_cells=2, **common)

    _assert_sample_aware_numerics_equal(low, high)

    assert not high.sample_sender_scores.isna().any().any()
    assert set(high.sample_receiver_scores[["sample", "receiver"]].itertuples(index=False, name=None)) == {
        ("S1", "A"),
        ("S1", "B"),
        ("S2", "A"),
        ("S2", "B"),
    }
    assert len(high.sample_events) == 8
    assert (high.events["n_samples_coobserved"] == 2).all()

    expected_unit_pass = {
        ("S1", "A"): True,
        ("S1", "B"): False,
        ("S2", "A"): False,
        ("S2", "B"): True,
    }
    validation_flags = high.sample_validation.set_index(["sample", "cell_type"])[
        "passes_min_cells"
    ].to_dict()
    assert validation_flags == expected_unit_pass
    assert "eligible_in_sample" not in high.sample_validation.columns
    assert high.sample_receiver_scores.set_index(["sample", "receiver"])[
        "receiver_passes_min_cells"
    ].to_dict() == expected_unit_pass
    assert "n_valid_samples" not in high.receiver_scores.columns
    assert (high.receiver_scores["n_observed_samples"] == 2).all()
    assert (high.receiver_scores["n_samples_passing_min_cells"] == 1).all()
    assert (high.receiver_scores["min_cells_pass_prevalence"] == 0.5).all()

    for _, row in high.sample_events.iterrows():
        sender_pass = expected_unit_pass[(row["sample"], row["sender"])]
        receiver_pass = expected_unit_pass[(row["sample"], row["receiver"])]
        assert bool(row["sender_passes_min_cells"]) is sender_pass
        assert bool(row["receiver_passes_min_cells"]) is receiver_pass
        assert bool(row["passes_min_cells"]) is (sender_pass and receiver_pass)

    aggregate = high.events.set_index(["sender", "receiver"])
    expected_aggregate = {
        ("A", "A"): (1, 0.5, True),
        ("A", "B"): (0, 0.0, False),
        ("B", "A"): (0, 0.0, False),
        ("B", "B"): (1, 0.5, True),
    }
    for key, (n_passing, prevalence, passes) in expected_aggregate.items():
        row = aggregate.loc[key]
        assert row["n_samples_passing_min_cells"] == n_passing
        assert row["min_cells_pass_prevalence"] == prevalence
        assert bool(row["passes_min_cells"]) is passes

    assert low.sample_validation["passes_min_cells"].all()
    assert low.sample_events["passes_min_cells"].all()
    assert (low.events["n_samples_passing_min_cells"] == 2).all()
    assert (low.events["min_cells_pass_prevalence"] == 1.0).all()
    assert low.events["passes_min_cells"].all()


def test_sample_aware_empty_events_preserve_aggregate_and_fdr_schema():
    frame = pd.DataFrame(
        {
            "sample": ["S1", "S1", "S2", "S2"],
            "cell_type": ["A", "A", "A", "A"],
            "PROD": [1.0, 1.0, 1.0, 1.0],
            "SENSOR": [1.0, 1.0, 1.0, 1.0],
        }
    )
    adata = AnnData(
        X=frame[["PROD", "SENSOR"]].to_numpy(dtype=float),
        obs=frame[["sample", "cell_type"]].copy(),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme, sensor = _priors()

    # A single observed cell type plus allow_self=False yields no event while
    # sender/receiver calculations themselves remain valid.
    result = run_cell_mesh(
        adata,
        enzyme,
        sensor,
        sample_key="sample",
        sample_mode="sample_aware",
        min_cells=2,
        min_expr_frac=None,
        allow_self=False,
        n_perms=2,
        random_state=0,
    )

    assert result.events.empty
    required = {
        "sender",
        "receiver",
        "metabolite",
        "hmdb_id",
        "sensor_gene",
        "sensor_type",
            "metabolite_availability_median",
            "sensor_score_median",
            "sensor_expr_frac_median",
        "sender_n_cells",
        "receiver_n_cells",
        "sender_passes_min_cells",
        "receiver_passes_min_cells",
        "passes_min_cells",
        "cell_mesh_score",
        "event_score_median",
        "event_score_iqr",
        "n_samples_coobserved",
        "n_samples_positive",
        "event_prevalence",
        "n_samples_passing_min_cells",
        "min_cells_pass_prevalence",
        "inference_mode",
        "permutation_mode",
        "perm_pvalue",
        "fdr_global",
        "fdr_sensor_type",
    }
    assert required.issubset(result.events.columns), required.difference(
        result.events.columns
    )
    assert "n_observed_samples" not in result.events.columns
    assert "n_valid_samples" not in result.events.columns
    assert result.sample_events is not None and result.sample_events.empty
    assert "sample_event_scores" not in result.availability_results
    null_scores = result.events.attrs.get("sample_aware_null_scores")
    assert null_scores is not None and null_scores.empty


def _empty_adata():
    return AnnData(
        X=np.empty((0, 2), dtype=float),
        obs=pd.DataFrame(
            {
                "sample": pd.Series(dtype=str),
                "cell_type": pd.Series(dtype=str),
            }
        ),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )


def test_empty_anndata_raises_the_same_error_in_pooled_and_sample_aware_modes():
    enzyme, sensor = _priors()
    errors = []
    for mode, extra in [
        ("pooled_stratified", {}),
        ("sample_aware", {"sample_key": "sample"}),
    ]:
        with pytest.raises(ValueError) as exc_info:
            run_cell_mesh(
                _empty_adata(),
                enzyme,
                sensor,
                sample_mode=mode,
                min_cells=1,
                n_perms=0,
                **extra,
            )
        errors.append(str(exc_info.value))

    assert errors[0] == errors[1]
    assert errors[0]


INVALID_MIN_CELLS = [
    pytest.param(True, id="true"),
    pytest.param(False, id="false"),
    pytest.param(0, id="zero"),
    pytest.param(-1, id="negative"),
    pytest.param(1.0, id="integer-valued-float"),
    pytest.param(2.5, id="fractional-float"),
    pytest.param(np.nan, id="nan"),
    pytest.param(np.inf, id="infinity"),
    pytest.param("2", id="string"),
    pytest.param(None, id="none"),
]


def _call_min_cells_entrypoint(entrypoint, min_cells):
    adata, enzyme, sensor = _pooled_case()
    if entrypoint == "run_cell_mesh":
        return run_cell_mesh(
            adata,
            enzyme,
            sensor,
            min_cells=min_cells,
            min_expr_frac=None,
            n_perms=0,
        )
    if entrypoint == "compute_metabolite_availability":
        return compute_metabolite_availability(
            adata,
            enzyme,
            min_cells=min_cells,
        )
    if entrypoint == "compute_sensor_scores":
        return compute_sensor_scores(
            adata,
            sensor,
            min_cells=min_cells,
            min_expr_frac=None,
        )
    raise AssertionError(entrypoint)


@pytest.mark.parametrize(
    "entrypoint",
    [
        "run_cell_mesh",
        "compute_metabolite_availability",
        "compute_sensor_scores",
    ],
)
@pytest.mark.parametrize("invalid", INVALID_MIN_CELLS)
def test_min_cells_rejects_non_positive_or_non_integer_values(entrypoint, invalid):
    error_type = ValueError if isinstance(invalid, (int, np.integer)) and not isinstance(invalid, (bool, np.bool_)) else TypeError
    with pytest.raises(error_type, match="min_cells"):
        _call_min_cells_entrypoint(entrypoint, invalid)


@pytest.mark.parametrize(
    "entrypoint",
    [
        "run_cell_mesh",
        "compute_metabolite_availability",
        "compute_sensor_scores",
    ],
)
def test_min_cells_accepts_numpy_integer(entrypoint):
    result = _call_min_cells_entrypoint(entrypoint, np.int64(2))
    if entrypoint == "run_cell_mesh":
        assert result.parameters["min_cells"] == 2
    elif entrypoint == "compute_metabolite_availability":
        assert not result["availability"].empty
    else:
        assert not result.empty
