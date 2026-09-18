import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import compute_metabolite_availability, load_cell_mesh_database, run_cell_mesh
from cellmesh.database import normalize_enzyme_database


@pytest.fixture
def alias_case():
    enzyme = pd.DataFrame({
        "metabolite": ["M", "M"], "hmdb_id": ["HMDB0000001"] * 2,
        "reaction": ["001", "002"], "gene": ["G1;G2", "EXP"],
        "role": ["production", "export"], "evidence_level": ["curated"] * 2,
        "source": ["user prior"] * 2, "reference": ["00123"] * 2,
    })
    sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["HMDB0000001"],
                           "sensor_gene": ["S"], "sensor_type": ["Transporter"]})
    adata = AnnData(
        X=np.array([[4, 9, 1, 2], [4, 9, 3, 1], [1, 4, 2, 7], [1, 4, 1, 8]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"], "sample": ["D1"] * 4},
                         index=["a1", "a2", "b1", "b2"]),
        var=pd.DataFrame(index=["G1", "G2", "EXP", "S"]),
    )
    return enzyme, sensor, adata


def _input(frame, tmp_path, csv_input):
    if not csv_input:
        return frame
    path = tmp_path / "enzyme.csv"
    frame.to_csv(path, index=False)
    return path


@pytest.mark.parametrize("csv_input", [False, True])
def test_equivalent_enzyme_aliases_merge_and_drop_stale_columns(alias_case, tmp_path, csv_input):
    enzyme, sensor, _ = alias_case
    mixed = enzyme.assign(
        standard_metName=[" M ", "M"], HMDB_ID=[" hmdb0000001 ", "HMDB0000001"],
        Reactions=["001", " 002 "], Gene_name=["G2[Unknown]|G1[Enzyme]", "EXP[Transporter]"],
        Direction=[" Product ", " Exporter "], direction=["product", "transporter"],
    )
    snapshot = mixed.copy(deep=True)
    normalized, _ = load_cell_mesh_database(_input(mixed, tmp_path, csv_input), sensor)
    pd.testing.assert_frame_equal(normalized, normalize_enzyme_database(enzyme))
    pd.testing.assert_frame_equal(normalize_enzyme_database(normalized), normalized)
    pd.testing.assert_frame_equal(mixed, snapshot)


@pytest.mark.parametrize("field,alias", [("metabolite", "standard_metName"), ("hmdb_id", "HMDB_ID"),
                                       ("reaction", "Reactions"), ("gene", "Gene_name"),
                                       ("role", "direction")])
@pytest.mark.parametrize("empty", [None, "", "   "])
@pytest.mark.parametrize("csv_input", [False, True])
def test_blank_enzyme_aliases_are_filled_row_by_row(alias_case, tmp_path, field, alias, empty, csv_input):
    enzyme, sensor, _ = alias_case
    mixed = enzyme.copy()
    mixed[alias] = ["product", "exporter"] if field == "role" else mixed[field]
    mixed.loc[0, field] = empty
    mixed.loc[1, alias] = empty
    normalized, _ = load_cell_mesh_database(_input(mixed, tmp_path, csv_input), sensor)
    pd.testing.assert_frame_equal(normalized, normalize_enzyme_database(enzyme))


@pytest.mark.parametrize("field,alias,value", [
    ("metabolite", "standard_metName", "different M"),
    ("hmdb_id", "HMDB_ID", "HMDB0000002"),
    ("reaction", "Reactions", "999"),
    ("gene", "Gene_name", "G1;DIFFERENT"),
    ("role", "direction", "substrate"),
])
@pytest.mark.parametrize("csv_input", [False, True])
def test_conflicting_enzyme_aliases_are_rejected_before_gene_filtering(
    alias_case, tmp_path, field, alias, value, csv_input,
):
    enzyme, sensor, adata = alias_case
    mixed = enzyme.copy()
    mixed[alias] = ["product", "exporter"] if field == "role" else mixed[field]
    mixed.loc[0, alias] = value
    prior_input = _input(mixed, tmp_path, csv_input)
    for call in (
        lambda: load_cell_mesh_database(prior_input, sensor),
        lambda: run_cell_mesh(adata, prior_input, sensor, min_cells=1),
        lambda: compute_metabolite_availability(adata, mixed, min_cells=1),
    ):
        with pytest.raises(ValueError, match=rf"conflicting {field}.*{alias}.*rows.*0"):
            call()


@pytest.mark.parametrize("first,second,expected", [
    ("", "product", "production"), (None, "substrate", "degradation"),
    ("exporter", "", "export"), ("export", "transporter", "export"),
])
@pytest.mark.parametrize("csv_input", [False, True])
def test_direction_aliases_merge_without_role(alias_case, tmp_path, first, second, expected, csv_input):
    enzyme, sensor, _ = alias_case
    mixed = enzyme.iloc[[0]].drop(columns="role").assign(Direction=first, direction=second)
    normalized, _ = load_cell_mesh_database(_input(mixed, tmp_path, csv_input), sensor)
    assert normalized["role"].tolist() == [expected, expected]
    assert "Direction" not in normalized and "direction" not in normalized


@pytest.mark.parametrize("csv_input", [False, True])
def test_conflicting_direction_columns_without_role_raise(alias_case, tmp_path, csv_input):
    enzyme, sensor, _ = alias_case
    mixed = enzyme.drop(columns="role").assign(Direction="substrate", direction="product")
    with pytest.raises(ValueError, match="conflicting role.*Direction.*direction.*rows"):
        load_cell_mesh_database(_input(mixed, tmp_path, csv_input), sensor)


def test_gene_alias_merge_preserves_inline_evidence(alias_case):
    enzyme, _, _ = alias_case
    mixed = enzyme.iloc[[0]].drop(columns="evidence_level").assign(
        gene="G1;G2[reviewed]", Gene_name="G2[curated];G1[Enzyme]",
    )
    normalized = normalize_enzyme_database(mixed)
    assert normalized[["gene", "evidence_level"]].to_records(index=False).tolist() == [
        ("G1", "Enzyme"), ("G2", "reviewed"), ("G2", "curated"),
    ]
    assert "Gene_name" not in normalized
    pd.testing.assert_frame_equal(normalize_enzyme_database(normalized), normalized)


@pytest.mark.parametrize("csv_input", [False, True])
def test_duplicate_enzyme_column_names_are_rejected(alias_case, tmp_path, csv_input):
    enzyme, sensor, _ = alias_case
    duplicate = pd.concat([enzyme, enzyme[["gene"]]], axis=1)
    with pytest.raises(ValueError, match="duplicate column names.*gene"):
        load_cell_mesh_database(_input(duplicate, tmp_path, csv_input), sensor)


def test_gene_alias_order_does_not_change_existing_display_name(alias_case):
    enzyme, _, _ = alias_case
    enzyme["metabolite"] = " M "
    mixed = enzyme.assign(standard_metName="M")
    pd.testing.assert_frame_equal(normalize_enzyme_database(mixed), normalize_enzyme_database(enzyme))


def test_raw_direction_normalization_preserves_existing_column_order():
    raw = pd.DataFrame({"metabolite": ["M"], "HMDB_ID": ["HMDB0000001"],
                        "reaction": ["001"], "gene": ["G1[Enzyme]"],
                        "direction": ["product"], "reference": ["00123"]})
    normalized = normalize_enzyme_database(raw)
    assert normalized.columns.tolist() == [
        "metabolite", "hmdb_id", "reaction", "gene", "reference", "role", "evidence_level",
    ]
    pd.testing.assert_frame_equal(normalize_enzyme_database(normalized), normalized)


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_merged_alias_analysis_matches_canonical_input(alias_case, tmp_path, sample_mode):
    enzyme, sensor, adata = alias_case
    mixed = enzyme.assign(HMDB_ID=enzyme["hmdb_id"], Gene_name=["G2|G1", "EXP"],
                          Reactions=enzyme["reaction"], Direction=["product", "export"],
                          direction=["product", "exporter"])
    mixed.loc[0, "gene"] = None
    mixed.loc[1, "role"] = ""
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, min_cells=1,
                  n_perms=5, random_state=37, store_null_scores=True)
    reference = run_cell_mesh(adata, enzyme, sensor, **kwargs)
    actual = run_cell_mesh(adata, _input(mixed, tmp_path, True), sensor, **kwargs)
    pd.testing.assert_frame_equal(actual.events, reference.events, check_exact=True)
    pd.testing.assert_frame_equal(actual.sender_scores, reference.sender_scores, check_exact=True)
    if sample_mode == "sample_aware":
        pd.testing.assert_frame_equal(actual.events.attrs["sample_aware_null_scores"],
                                      reference.events.attrs["sample_aware_null_scores"], check_exact=True)
