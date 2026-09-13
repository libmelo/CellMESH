import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import compute_metabolite_availability, load_cell_mesh_database, run_cell_mesh
from cellmesh.database import normalize_enzyme_database


@pytest.fixture
def priors():
    enzyme = pd.DataFrame({
        "metabolite": ["代谢物, A"] * 3 + ["M2"],
        "hmdb_id": ["HMDB0000001"] * 3 + ["HMDB0000002"],
        "gene": ["PROD", "CONS", "EXPORT", "PROD2"],
        "role": ["production", "degradation", "export", "production"],
        "reaction": ["001", "002", "003", "004"],
        "evidence_level": ["curated"] * 4,
        "source": ["user enzyme"] * 4,
        "reference": ["00123"] * 4,
        "confidence": [0.75] * 4,
    })
    sensor = pd.DataFrame({
        "metabolite": ["代谢物, A", "M2"],
        "hmdb_id": ["HMDB0000001", "HMDB0000002"],
        "sensor_gene": ["S1", "S2"],
        "sensor_type": ["Transporter", "Cell surface receptor"],
        "evidence_level": ["measured", "curated"],
        "source": ["user sensor"] * 2,
        "protein_name": ["Protein 1", "Protein 2"],
        "reference": ["00456"] * 2,
        "confidence": [0.5] * 2,
    })
    return enzyme, sensor


def _as_input(frame, path, csv_input):
    if csv_input:
        frame.to_csv(path, index=False)
        return path
    return frame


def test_normalized_csv_roundtrip_preserves_records_metadata_and_inputs(tmp_path, priors):
    enzyme, sensor = priors
    snapshots = [frame.copy(deep=True) for frame in priors]
    paths = [_as_input(frame, tmp_path / f"{i}.csv", True) for i, frame in enumerate(priors)]
    direct = load_cell_mesh_database(enzyme, sensor)
    restored = load_cell_mesh_database(*paths)

    for original, snapshot, expected, actual in zip(priors, snapshots, direct, restored):
        pd.testing.assert_frame_equal(actual, expected)
        pd.testing.assert_frame_equal(actual[original.columns], original)
        pd.testing.assert_frame_equal(original, snapshot)
    assert set(restored[0]["role"]) == {"production", "degradation", "export"}


@pytest.mark.parametrize("uppercase", [False, True])
@pytest.mark.parametrize("csv_input", [False, True])
@pytest.mark.parametrize("export_direction", ["exporter", "export"])
def test_direction_schemas_match_role_schema(
    tmp_path, priors, uppercase, csv_input, export_direction,
):
    enzyme, sensor = priors
    raw = enzyme.copy()
    raw["direction"] = raw.pop("role").map({
        "production": "product", "degradation": "substrate", "export": export_direction,
    })
    if uppercase:
        raw = raw.rename(columns={
            "metabolite": "standard_metName", "hmdb_id": "HMDB_ID",
            "gene": "Gene_name", "reaction": "Reactions", "direction": "Direction",
        })
    loaded, _ = load_cell_mesh_database(
        _as_input(raw, tmp_path / "enzyme.csv", csv_input), sensor,
    )
    pd.testing.assert_frame_equal(loaded[enzyme.columns], enzyme)


def test_direction_schema_still_expands_gene_evidence():
    raw = pd.DataFrame({
        "standard_metName": ["M1"], "HMDB_ID": ["HMDB0000001"],
        "Reactions": ["r1"], "Direction": ["product"],
        "Gene_name": ["PROD[reviewed]; PROD2[curated]"],
    })
    loaded = normalize_enzyme_database(raw)
    assert loaded["gene"].tolist() == ["PROD", "PROD2"]
    assert loaded["evidence_level"].tolist() == ["reviewed", "curated"]
    assert loaded["role"].tolist() == ["production", "production"]


@pytest.mark.parametrize("direction_column", ["Direction", "direction"])
@pytest.mark.parametrize("csv_input", [False, True])
@pytest.mark.parametrize("export_direction", ["exporter", "export"])
def test_role_direction_agreement_and_conflict(
    tmp_path, priors, direction_column, csv_input, export_direction,
):
    enzyme, sensor = priors
    enzyme[direction_column] = ["product", "substrate", export_direction, "product"]
    enzyme_input = _as_input(enzyme, tmp_path / "enzyme.csv", csv_input)
    loaded, _ = load_cell_mesh_database(enzyme_input, sensor)
    assert loaded["role"].tolist() == enzyme["role"].tolist()

    enzyme.loc[0, direction_column] = "substrate"
    enzyme_input = _as_input(enzyme, tmp_path / "enzyme.csv", csv_input)
    with pytest.raises(ValueError, match="conflicting role.*direction"):
        load_cell_mesh_database(enzyme_input, sensor)


@pytest.mark.parametrize("missing", ["role", "reaction", "sensor_gene"])
@pytest.mark.parametrize("csv_input", [False, True])
def test_missing_required_fields_fail_explicitly(tmp_path, priors, missing, csv_input):
    enzyme, sensor = priors
    enzyme = enzyme.drop(columns=[missing], errors="ignore")
    sensor = sensor.drop(columns=[missing], errors="ignore")
    inputs = [_as_input(frame, tmp_path / f"{i}.csv", csv_input)
              for i, frame in enumerate([enzyme, sensor])]
    with pytest.raises(ValueError, match=missing):
        load_cell_mesh_database(*inputs)


@pytest.mark.parametrize("raw_schema", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("csv_input", [False, True])
def test_sensor_conflicts_are_rejected_before_deduplication(
    tmp_path, priors, raw_schema, reverse, csv_input,
):
    enzyme, sensor = priors
    conflict = sensor.iloc[[0, 0]].copy()
    conflict.iloc[1, conflict.columns.get_loc("hmdb_id")] = " hmdb0000001 "
    conflict.iloc[1, conflict.columns.get_loc("sensor_type")] = "Cell surface receptor"
    if reverse:
        conflict = conflict.iloc[::-1]
    if raw_schema:
        conflict = conflict.rename(columns={
            "metabolite": "standard_metName", "hmdb_id": "HMDB_ID",
            "sensor_gene": "Gene_name", "sensor_type": "Annotation",
        })
    sensor_input = _as_input(conflict, tmp_path / "sensor.csv", csv_input)
    with pytest.raises(ValueError, match="conflicting sensor_type.*HMDB0000001.*S1"):
        load_cell_mesh_database(enzyme, sensor_input)


@pytest.mark.parametrize("csv_input", [False, True])
def test_same_type_duplicates_keep_first_metadata(tmp_path, priors, csv_input):
    enzyme, sensor = priors
    duplicated = sensor.iloc[[0, 0, 1]].copy()
    duplicated.iloc[1, duplicated.columns.get_loc("hmdb_id")] = " hmdb0000001 "
    duplicated.iloc[1, duplicated.columns.get_loc("source")] = "later source"
    _, loaded = load_cell_mesh_database(
        enzyme, _as_input(duplicated, tmp_path / "sensor.csv", csv_input),
    )
    pd.testing.assert_frame_equal(loaded[sensor.columns], sensor)


@pytest.fixture
def expression_data():
    return AnnData(
        X=np.array([
            [8, 1, 5, 2, 1, 3], [6, 2, 3, 0, 2, 1],
            [0, 4, 1, 8, 9, 2], [1, 5, 0, 6, 7, 1],
            [7, 1, 4, 3, 0, 7], [5, 3, 5, 2, 1, 5],
            [2, 6, 0, 7, 8, 3], [0, 2, 2, 9, 6, 1],
        ], dtype=float),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B"] * 2, "sample": ["D1"] * 4 + ["D2"] * 4},
            index=[f"cell_{i}" for i in range(8)],
        ),
        var=pd.DataFrame(index=["PROD", "CONS", "EXPORT", "PROD2", "S1", "S2"]),
    )


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("csv_prior", ["enzyme", "sensor", "both"])
@pytest.mark.parametrize("enzyme_schema", ["role", "direction_export"])
def test_csv_and_dataframe_analysis_results_match(
    tmp_path, priors, expression_data, sample_mode, csv_prior, enzyme_schema,
):
    adata = expression_data
    kwargs = {"sample_key": "sample", "sample_mode": sample_mode,
              "min_cells": 1, "n_perms": 9, "random_state": 11}
    expected = run_cell_mesh(adata, *priors, **kwargs)
    if enzyme_schema == "direction_export":
        enzyme = priors[0].copy()
        enzyme["direction"] = enzyme.pop("role").map({
            "production": "product", "degradation": "substrate", "export": " Export ",
        })
        priors = (enzyme, priors[1])
    inputs = [_as_input(frame, tmp_path / f"{name}.csv", csv_prior in (name, "both"))
              for name, frame in zip(["enzyme", "sensor"], priors)]
    actual = run_cell_mesh(adata, *inputs, **kwargs)
    assert not actual.events.empty
    assert {"perm_pvalue", "fdr_global", "fdr_sensor_type"}.issubset(actual.events)
    for name in ["events", "sender_scores", "receiver_scores",
                 "sample_sender_scores", "sample_receiver_scores", "sample_events"]:
        expected_table = getattr(expected, name)
        if expected_table is not None:
            pd.testing.assert_frame_equal(getattr(actual, name), expected_table, check_exact=True)

    expected_availability = expected.availability_results.get(
        "availability_by_sample", {"pooled": expected.availability_results},
    )
    actual_availability = actual.availability_results.get(
        "availability_by_sample", {"pooled": actual.availability_results},
    )
    for sample, reference in expected_availability.items():
        # An export row must reach E scoring, not fall back to missing-export evidence.
        assert reference["E"].to_numpy().max() > 0
        for name in ("E", "E_score", "E_effective", "E_factor"):
            pd.testing.assert_frame_equal(
                actual_availability[sample][name], reference[name], check_exact=True,
            )


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("schema", ["raw", "mixed", "blank_direction"])
def test_raw_direction_cannot_override_normalized_role(
    tmp_path, priors, expression_data, sample_mode, schema,
):
    enzyme, sensor = priors
    kwargs = {"sample_key": "sample", "sample_mode": sample_mode,
              "min_cells": 1, "n_perms": 9, "random_state": 11}
    expected = run_cell_mesh(expression_data, enzyme, sensor, **kwargs)
    alternative = enzyme.copy()
    alternative["direction"] = [" Product ", " SUBSTRATE ", " Exporter ", " Product "]
    if schema == "raw":
        alternative = alternative.drop(columns="role")
    elif schema == "blank_direction":
        alternative["direction"] = ["", None, "  ", np.nan]

    for csv_input in (False, True):
        actual = run_cell_mesh(
            expression_data,
            _as_input(alternative, tmp_path / "enzyme.csv", csv_input), sensor,
            **kwargs,
        )
        for name in ("events", "sender_scores", "receiver_scores", "sample_events"):
            expected_table = getattr(expected, name)
            if expected_table is not None:
                pd.testing.assert_frame_equal(getattr(actual, name), expected_table, check_exact=True)


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("schema", [
    "role", "export", "exporter", "transporter", "legacy_fields", "mixed", "blank_direction",
])
def test_public_availability_entrypoints_share_prior_normalization(
    priors, expression_data, sample_mode, schema,
):
    enzyme, sensor = priors
    alternative = enzyme.copy()
    if schema == "role":
        alternative["role"] = alternative["role"].map(lambda value: f" {value.upper()} ")
    else:
        export_direction = schema if schema in ("export", "exporter", "transporter") else "export"
        alternative["direction"] = ["product", "substrate", export_direction, "product"]
        if schema in ("mixed", "blank_direction"):
            if schema == "mixed":
                alternative["direction"] = alternative["direction"].map(lambda value: f" {value.upper()} ")
            else:
                alternative["direction"] = ["", None, "  ", np.nan]
        else:
            alternative = alternative.drop(columns="role")
        if schema == "legacy_fields":
            alternative = alternative.rename(columns={
                "metabolite": "standard_metName", "hmdb_id": "HMDB_ID",
                "gene": "Gene_name", "reaction": "Reactions", "direction": "Direction",
            })
    snapshot = alternative.copy(deep=True)
    kwargs = {"sample_key": "sample", "sample_mode": sample_mode, "min_cells": 1, "n_perms": 0}
    expected = run_cell_mesh(expression_data, enzyme, sensor, **kwargs).availability_results
    via_main = run_cell_mesh(expression_data, alternative, sensor, **kwargs).availability_results
    if sample_mode == "sample_aware":
        expected = expected["availability_by_sample"]
        via_main = via_main["availability_by_sample"]
        data_by_sample = {
            sample: expression_data[expression_data.obs["sample"].eq(sample)].copy()
            for sample in expected
        }
    else:
        expected = {"pooled": expected}
        via_main = {"pooled": via_main}
        data_by_sample = {"pooled": expression_data}

    for sample, reference in expected.items():
        direct = compute_metabolite_availability(data_by_sample[sample], alternative, min_cells=1)
        assert reference["E"].to_numpy().max() > 0
        for name in ("P", "C", "E", "P_score", "C_score", "E_score",
                     "E_effective", "E_factor", "availability"):
            pd.testing.assert_frame_equal(direct[name], reference[name], check_exact=True)
            pd.testing.assert_frame_equal(via_main[sample][name], reference[name], check_exact=True)
    pd.testing.assert_frame_equal(alternative, snapshot)


@pytest.mark.parametrize("direction_column", ["Direction", "direction"])
def test_both_public_entrypoints_reject_role_direction_conflicts(
    priors, expression_data, direction_column,
):
    enzyme, sensor = priors
    enzyme[direction_column] = ["product", "substrate", "product", "product"]
    with pytest.raises(ValueError, match="conflicting role and direction"):
        run_cell_mesh(expression_data, enzyme, sensor, min_cells=1, n_perms=0)
    with pytest.raises(ValueError, match="conflicting role and direction"):
        compute_metabolite_availability(expression_data, enzyme, min_cells=1)


@pytest.mark.parametrize("export_evidence", ["zero", "unmeasured", "absent"])
def test_public_entrypoints_agree_on_export_evidence_states(priors, expression_data, export_evidence):
    enzyme, sensor = priors
    if export_evidence == "zero":
        expression_data.X[:, expression_data.var_names == "EXPORT"] = 0
    elif export_evidence == "unmeasured":
        expression_data = expression_data[:, expression_data.var_names != "EXPORT"].copy()
    else:
        enzyme = enzyme.loc[enzyme["role"].ne("export")].copy()
    raw = enzyme.copy()
    raw["direction"] = raw.pop("role").map({
        "production": "product", "degradation": "substrate", "export": "export",
    })
    reference = run_cell_mesh(expression_data, enzyme, sensor, min_cells=1, n_perms=0).availability_results
    for candidate in (
        run_cell_mesh(expression_data, raw, sensor, min_cells=1, n_perms=0).availability_results,
        compute_metabolite_availability(expression_data, raw, min_cells=1),
    ):
        for name in ("E", "E_score", "E_effective", "E_factor", "availability"):
            pd.testing.assert_frame_equal(candidate[name], reference[name], check_exact=True)
        pd.testing.assert_series_equal(candidate["metadata"]["export_status"], reference["metadata"]["export_status"])
