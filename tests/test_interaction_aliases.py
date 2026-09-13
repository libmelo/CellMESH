import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import load_cell_mesh_database, run_cell_mesh
from cellmesh.database import normalize_interaction_database


@pytest.fixture
def interaction_case():
    enzyme = pd.DataFrame({
        "metabolite": ["M"], "hmdb_id": ["HMDB0000001"],
        "gene": ["E1"], "reaction": ["001"], "role": ["production"],
    })
    sensor = pd.DataFrame({
        "metabolite": ["M", "M"], "hmdb_id": ["HMDB0000001"] * 2,
        "sensor_gene": ["S1", "S2"],
        "sensor_type": ["Transporter", "Cell surface receptor"],
        "source": ["user database", "NA"], "protein_name": ["Protein 1", "Protein 2"],
        "reference": ["00123", "00456"], "evidence_level": ["curated", "reviewed"],
        "confidence": [0.75, 0.5], "ID": [1, 2],
    })
    adata = AnnData(
        X=np.array([[4, 2, 1], [8, 1, 2], [1, 7, 4], [2, 8, 5],
                    [6, 3, 2], [5, 2, 1], [2, 5, 6], [1, 4, 7]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"] * 2,
                          "sample": ["D1"] * 4 + ["D2"] * 4},
                         index=[f"cell{i}" for i in range(8)]),
        var=pd.DataFrame(index=["E1", "S1", "S2"]),
    )
    return enzyme, sensor, adata


def _input(frame, tmp_path, csv_input):
    if not csv_input:
        return frame
    path = tmp_path / "interaction.csv"
    frame.to_csv(path, index=False, encoding="utf-8-sig")
    return path


@pytest.mark.parametrize("field,alias", [
    ("metabolite", "standard_metName"), ("metabolite", "standard_metname"),
    ("hmdb_id", "HMDB_ID"),
    ("sensor_gene", "Gene_name"), ("sensor_gene", "gene_name"), ("sensor_gene", "gene"),
    ("sensor_type", "Annotation"), ("sensor_type", "annotation"),
    ("source", "Database source"), ("source", "database_source"),
    ("protein_name", "Protein_name"), ("reference", "Reference"),
])
@pytest.mark.parametrize("csv_input", [False, True])
def test_each_documented_interaction_alias_is_used(interaction_case, tmp_path, field, alias, csv_input):
    enzyme, sensor, _ = interaction_case
    alternative = sensor.rename(columns={field: alias})
    _, actual = load_cell_mesh_database(enzyme, _input(alternative, tmp_path, csv_input))
    pd.testing.assert_frame_equal(actual, sensor)


def _all_aliases(sensor):
    mixed = sensor.assign(
        standard_metName=" M ", standard_metname="M", HMDB_ID=" hmdb0000001 ",
        Gene_name=[" S1 ", "S2"], gene_name=["S1", " S2 "], gene=["S1", "S2"],
        Annotation=["transport protein", "surface receptor"],
        annotation=["TRANSPORTER", " CELL SURFACE RECEPTOR "],
        database_source=sensor["source"], Protein_name=sensor["protein_name"],
        Reference=sensor["reference"],
    )
    mixed["Database source"] = sensor["source"]
    return mixed


@pytest.mark.parametrize("csv_input", [False, True])
def test_equivalent_interaction_aliases_merge_without_changing_metadata(interaction_case, tmp_path, csv_input):
    enzyme, sensor, _ = interaction_case
    mixed = _all_aliases(sensor)
    snapshot = mixed.copy(deep=True)
    _, actual = load_cell_mesh_database(enzyme, _input(mixed, tmp_path, csv_input))
    pd.testing.assert_frame_equal(actual, sensor)
    pd.testing.assert_frame_equal(normalize_interaction_database(actual), actual)
    pd.testing.assert_frame_equal(mixed, snapshot)


@pytest.mark.parametrize("field,alias", [
    ("metabolite", "standard_metName"), ("hmdb_id", "HMDB_ID"),
    ("sensor_gene", "Gene_name"), ("sensor_type", "Annotation"),
    ("source", "Database source"), ("protein_name", "Protein_name"), ("reference", "Reference"),
])
@pytest.mark.parametrize("empty", [None, "", "   "])
@pytest.mark.parametrize("csv_input", [False, True])
def test_interaction_aliases_fill_blanks_in_both_directions(
    interaction_case, tmp_path, field, alias, empty, csv_input,
):
    enzyme, sensor, _ = interaction_case
    mixed = sensor.copy()
    mixed[alias] = mixed[field]
    mixed.loc[0, field] = empty
    mixed.loc[1, alias] = empty
    _, actual = load_cell_mesh_database(enzyme, _input(mixed, tmp_path, csv_input))
    pd.testing.assert_frame_equal(actual, sensor)


@pytest.mark.parametrize("field,alias,value", [
    ("metabolite", "standard_metName", "different metabolite"),
    ("hmdb_id", "HMDB_ID", "HMDB0000002"),
    ("sensor_gene", "Gene_name", "different_gene"),
    ("sensor_type", "Annotation", "Cell surface receptor"),
    ("source", "Database source", "different source"),
    ("protein_name", "Protein_name", "different protein"),
    ("reference", "Reference", "99999"),
])
@pytest.mark.parametrize("csv_input", [False, True])
def test_interaction_alias_conflicts_fail_before_expression_filtering(
    interaction_case, tmp_path, field, alias, value, csv_input,
):
    enzyme, sensor, adata = interaction_case
    mixed = sensor.copy()
    mixed.loc[0, "sensor_gene"] = "UNMEASURED"
    mixed[alias] = mixed[field]
    mixed.loc[0, alias] = value
    prior_input = _input(mixed, tmp_path, csv_input)
    for call in (
        lambda: normalize_interaction_database(mixed),
        lambda: load_cell_mesh_database(enzyme, prior_input),
        lambda: run_cell_mesh(adata, enzyme, prior_input, min_cells=1),
    ):
        with pytest.raises(ValueError, match=rf"conflicting {field}.*{alias}.*rows.*1 .*{value}"):
            call()


@pytest.mark.parametrize("first,second,third", [
    (None, "S1", " S1 "), ("S1", "", "S1"), ("S1", " S1 ", None),
])
@pytest.mark.parametrize("csv_input", [False, True])
def test_gene_aliases_merge_without_canonical_column(
    interaction_case, tmp_path, first, second, third, csv_input,
):
    enzyme, sensor, _ = interaction_case
    raw = sensor.iloc[[0]].drop(columns="sensor_gene").assign(
        Gene_name=first, gene_name=second, gene=third,
    )
    _, actual = load_cell_mesh_database(enzyme, _input(raw, tmp_path, csv_input))
    pd.testing.assert_frame_equal(actual[sensor.columns], sensor.iloc[[0]])


@pytest.mark.parametrize("with_canonical", [False, True])
def test_later_gene_alias_cannot_hide_a_conflict(interaction_case, with_canonical):
    _, sensor, _ = interaction_case
    mixed = sensor.assign(Gene_name=sensor["sensor_gene"], gene_name=sensor["sensor_gene"], gene="WRONG")
    if not with_canonical:
        mixed = mixed.drop(columns="sensor_gene")
    with pytest.raises(ValueError, match="conflicting sensor_gene.*Gene_name.*gene_name.*WRONG"):
        normalize_interaction_database(mixed)


@pytest.mark.parametrize("column", ["sensor_gene", "Gene_name", "Annotation", "reference", "ID"])
@pytest.mark.parametrize("csv_input", [False, True])
def test_duplicate_interaction_headers_are_rejected_even_when_values_agree(
    interaction_case, tmp_path, column, csv_input,
):
    enzyme, sensor, _ = interaction_case
    frame = sensor.assign(Gene_name=sensor["sensor_gene"], Annotation=sensor["sensor_type"])
    duplicate = pd.concat([frame, frame[[column]]], axis=1)
    with pytest.raises(ValueError, match=rf"metabolite_sensor has duplicate column names.*{column}"):
        load_cell_mesh_database(enzyme, _input(duplicate, tmp_path, csv_input))
    with pytest.raises(ValueError, match=rf"duplicate column names.*{column}"):
        normalize_interaction_database(duplicate)


@pytest.mark.parametrize("gene,alias", [("S1", "s1"), ("S1;S2", "S2;S1"), ("S1", "S1[reviewed]")])
def test_interaction_symbols_do_not_use_enzyme_gene_set_semantics(interaction_case, gene, alias):
    _, sensor, _ = interaction_case
    mixed = sensor.iloc[[0]].assign(sensor_gene=gene, Gene_name=alias)
    with pytest.raises(ValueError, match="conflicting sensor_gene"):
        normalize_interaction_database(mixed)


@pytest.mark.parametrize("csv_input", [False, True])
def test_raw_annotation_evidence_survives_semantic_type_merging(interaction_case, tmp_path, csv_input):
    enzyme, sensor, _ = interaction_case
    raw = sensor.drop(columns=["sensor_type", "evidence_level"]).assign(
        Annotation=["intracellular receptor", "surface receptor (curated)"],
        annotation=["nuclear receptor", "surface receptor (curated)"],
    )
    _, actual = load_cell_mesh_database(enzyme, _input(raw, tmp_path, csv_input))
    assert actual["sensor_type"].tolist() == ["Other receptor", "Cell surface receptor"]
    assert actual["evidence_level"].tolist() == [
        "intracellular receptor; nuclear receptor", "surface receptor (curated)",
    ]
    assert "Annotation" not in actual and "annotation" not in actual
    pd.testing.assert_frame_equal(normalize_interaction_database(actual), actual)


@pytest.mark.parametrize("with_canonical", [False, True])
def test_annotation_fills_missing_evidence_but_preserves_explicit_evidence(interaction_case, with_canonical):
    _, sensor, _ = interaction_case
    mixed = sensor.assign(Annotation=["transport protein", "surface receptor (curated)"])
    mixed.loc[1, "evidence_level"] = pd.NA
    if not with_canonical:
        mixed = mixed.drop(columns="sensor_type")
    actual = normalize_interaction_database(mixed)
    assert actual["evidence_level"].tolist() == ["curated", "surface receptor (curated)"]


def test_standard_type_keeps_raw_annotation_as_evidence(interaction_case):
    _, sensor, _ = interaction_case
    mixed = sensor.drop(columns="evidence_level").assign(Annotation=["transport protein", "surface receptor"])
    actual = normalize_interaction_database(mixed)
    assert actual["evidence_level"].tolist() == ["transport protein", "surface receptor"]


def test_explicit_other_receptor_does_not_count_as_missing(interaction_case):
    _, sensor, _ = interaction_case
    mixed = sensor.iloc[[0]].assign(sensor_type="Other receptor", Annotation="Transporter")
    with pytest.raises(ValueError, match="conflicting sensor_type"):
        normalize_interaction_database(mixed)


def test_all_empty_type_aliases_keep_existing_other_receptor_fallback(interaction_case):
    _, sensor, _ = interaction_case
    mixed = sensor.assign(sensor_type=None, Annotation=" ", annotation=pd.NA)
    actual = normalize_interaction_database(mixed)
    assert actual["sensor_type"].tolist() == ["Other receptor", "Other receptor"]


def test_missing_gene_and_hmdb_sentinels_can_be_filled(interaction_case):
    _, sensor, _ = interaction_case
    mixed = sensor.assign(Gene_name=sensor["sensor_gene"], HMDB_ID=sensor["hmdb_id"])
    mixed["sensor_gene"] = [" nan ", pd.NA]
    mixed["hmdb_id"] = ["None", " null "]
    pd.testing.assert_frame_equal(normalize_interaction_database(mixed), sensor)


def test_pair_type_conflict_is_checked_after_alias_fallback(interaction_case):
    _, sensor, _ = interaction_case
    mixed = sensor.assign(sensor_gene=["S1", None], Gene_name=[None, "S1"])
    with pytest.raises(ValueError, match="conflicting sensor_type.*HMDB0000001.*S1"):
        normalize_interaction_database(mixed)


def test_alias_conflicts_are_checked_before_blank_gene_rows_are_dropped(interaction_case):
    _, sensor, _ = interaction_case
    mixed = sensor.assign(sensor_gene=None, HMDB_ID="HMDB0000002")
    with pytest.raises(ValueError, match="conflicting hmdb_id"):
        normalize_interaction_database(mixed)


def test_duplicate_row_indices_do_not_mix_up_alias_values(interaction_case):
    _, sensor, _ = interaction_case
    mixed = sensor.assign(Gene_name=sensor["sensor_gene"], sensor_gene=[None, "S2"])
    mixed.index = ["same", "same"]
    pd.testing.assert_frame_equal(normalize_interaction_database(mixed), sensor)
    mixed.iloc[1, mixed.columns.get_loc("Gene_name")] = "WRONG"
    with pytest.raises(ValueError, match="rows.*2 .*index='same'.*WRONG"):
        normalize_interaction_database(mixed)


@pytest.mark.parametrize("csv_input", [False, True])
def test_empty_interaction_with_mixed_schema_remains_loadable(interaction_case, tmp_path, csv_input):
    enzyme, sensor, _ = interaction_case
    _, actual = load_cell_mesh_database(enzyme, _input(_all_aliases(sensor).iloc[:0], tmp_path, csv_input))
    assert actual.empty
    assert actual.columns.tolist() == sensor.columns.tolist()


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("csv_input", [False, True])
def test_interaction_alias_merge_preserves_scores_and_permutation_results(
    interaction_case, tmp_path, sample_mode, csv_input,
):
    enzyme, sensor, adata = interaction_case
    mixed = _all_aliases(sensor)
    mixed.loc[0, "sensor_gene"] = None
    mixed.loc[1, "hmdb_id"] = ""
    mixed.loc[0, "sensor_type"] = " "
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, min_cells=1,
                  n_perms=5, random_state=37, store_null_scores=True)
    expected = run_cell_mesh(adata, enzyme, sensor, **kwargs)
    actual = run_cell_mesh(adata, enzyme, _input(mixed, tmp_path, csv_input), **kwargs)
    for name in ("events", "sender_scores", "receiver_scores", "sample_events",
                 "sample_sender_scores", "sample_receiver_scores"):
        reference = getattr(expected, name)
        if reference is not None:
            pd.testing.assert_frame_equal(getattr(actual, name), reference, check_exact=True)
    if sample_mode == "sample_aware":
        pd.testing.assert_frame_equal(actual.events.attrs["sample_aware_null_scores"],
                                      expected.events.attrs["sample_aware_null_scores"], check_exact=True)
