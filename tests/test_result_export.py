import csv
import json

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import run_cell_mesh


def _result(sample_mode="pooled_stratified", zero_production=False, production_prior=True,
            n_perms=3, store_null_scores=False):
    expression = np.array([[5, 1], [3, 1], [1, 8], [8, 0], [0, 5], [1, 9], [4, 0]], dtype=float)
    if zero_production:
        expression[:, 0] = 0.0
    adata = AnnData(
        expression,
        obs=pd.DataFrame(
            {
                "cell_type": ["A", "A", "B", "A", "B", "B", "A"],
                "sample": ["S1", "S1", "S1", "S2", "S2", "S2", "S3"],
            },
            index=[f"cell_{i}" for i in range(7)],
        ),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme = pd.DataFrame([{
        "metabolite": "代谢物, A", "hmdb_id": "HMDB0000001", "gene": "PROD",
        "role": "production" if production_prior else "degradation", "reaction": "prod",
    }])
    sensor = pd.DataFrame([{
        "metabolite": "代谢物, A", "hmdb_id": "HMDB0000001",
        "sensor_gene": "SENSOR", "sensor_type": "Transporter",
    }])
    return run_cell_mesh(
        adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
        min_cells=2, n_perms=n_perms, random_state=np.int64(7), store_null_scores=store_null_scores,
    )


def _read_table(path):
    # Inspect the raw header: pandas can silently rename duplicate CSV columns.
    with path.open(encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream))
    assert len(header) == len(set(header))
    assert "index" not in header
    return pd.read_csv(path, float_precision="round_trip"), header


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_to_csv_roundtrip_preserves_identifiers_values_and_parameters(tmp_path, sample_mode):
    result = _result(sample_mode)
    snapshots = {
        name: getattr(result, name).copy(deep=True)
        for name in (
            "events", "sender_scores", "receiver_scores", "celltype_qc",
            "sample_validation", "sample_sender_scores", "sample_receiver_scores", "sample_events",
        )
        if getattr(result, name) is not None
    }
    prefix = tmp_path / "result"
    result.to_csv(prefix)

    expected_files = {f"result.{name}.csv" for name in snapshots}
    expected_files.update(["result.parameters.json", "result.manifest.json"])
    assert {path.name for path in tmp_path.iterdir()} == expected_files

    indexed_tables = {"sender_scores", "sample_sender_scores", "celltype_qc"}
    for name, expected in snapshots.items():
        path = tmp_path / f"result.{name}.csv"
        actual, header = _read_table(path)
        if name in indexed_tables:
            index_names = list(expected.index.names)
            assert header == index_names + expected.columns.tolist()
            actual = actual.set_index(index_names)
        else:
            assert header == expected.columns.tolist()
            # CSV has no dtype metadata; restore nullable booleans for NA units.
            for column in expected.columns:
                if isinstance(expected[column].dtype, pd.BooleanDtype):
                    actual[column] = pd.array(actual[column], dtype="boolean")
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False, check_exact=True)
        pd.testing.assert_frame_equal(getattr(result, name), expected)

    parameters = json.loads((tmp_path / "result.parameters.json").read_text(encoding="utf-8"))
    assert parameters == result.parameters
    assert parameters["random_state"] == 7
    if sample_mode == "sample_aware":
        assert result.sample_events["cell_mesh_score"].isna().any()


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_to_csv_empty_scores_still_export_identifier_headers(tmp_path, sample_mode):
    # Missing production evidence yields an empty sender table; measured zero
    # production is now a valid score and is covered separately below.
    result = _result(sample_mode, production_prior=False)
    assert result.sender_scores.empty
    result.to_csv(tmp_path / "empty")

    scores, header = _read_table(tmp_path / "empty.sender_scores.csv")
    assert scores.empty
    assert header[:2] == ["metabolite", "hmdb_id"]
    events, _ = _read_table(tmp_path / "empty.events.csv")
    assert events.empty
    if sample_mode == "sample_aware":
        sample_scores, header = _read_table(tmp_path / "empty.sample_sender_scores.csv")
        assert sample_scores.empty
        assert header == ["sample", "metabolite", "hmdb_id", "A", "B"]


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_to_csv_keeps_measured_zero_scores_distinct_from_missing_samples(tmp_path, sample_mode):
    result = _result(sample_mode, zero_production=True)
    result.to_csv(tmp_path / "zero")
    events, _ = _read_table(tmp_path / "zero.events.csv")
    assert not events.empty
    assert events["cell_mesh_score"].eq(0.0).all()
    assert events["perm_pvalue"].eq(1.0).all()
    if sample_mode == "sample_aware":
        sample_events, _ = _read_table(tmp_path / "zero.sample_events.csv")
        pair = sample_events.loc[sample_events["sender"].eq("A") & sample_events["receiver"].eq("B")]
        assert pair.set_index("sample").loc[["S1", "S2"], "cell_mesh_score"].eq(0.0).all()
        assert pd.isna(pair.set_index("sample").loc["S3", "cell_mesh_score"])


def test_to_csv_does_not_add_identifiers_twice_to_flattened_scores(tmp_path):
    result = _result()
    result.sender_scores = result.sender_scores.reset_index()
    result.to_csv(tmp_path / "flat")
    actual, header = _read_table(tmp_path / "flat.sender_scores.csv")
    assert header == ["metabolite", "hmdb_id", "A", "B"]
    pd.testing.assert_frame_equal(actual, result.sender_scores, check_exact=True)


@pytest.mark.parametrize("table", ["sender_scores", "events", "receiver_scores"])
def test_to_csv_rejects_duplicate_headers_before_writing_files(tmp_path, table):
    result = _result()
    if table == "sender_scores":
        result.sender_scores = result.sender_scores.rename(columns={"A": "hmdb_id"})
    else:
        frame = getattr(result, table)
        frame.insert(0, "hmdb_id", frame["hmdb_id"], allow_duplicates=True)

    with pytest.raises(ValueError, match=table):
        result.to_csv(tmp_path / "conflict")
    assert list(tmp_path.iterdir()) == []
