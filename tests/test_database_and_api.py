import numpy as np
import pandas as pd
import pytest

from cellmesh import load_cell_mesh_database, read_anndata, run_cell_mesh
from cellmesh.database import (
    _default_database_paths,
    _find_versioned_database_files,
    _select_default_database_paths,
    normalize_interaction_database,
    validate_priors,
)


class FakeAnnData:
    def __init__(self, X, var_names, obs):
        self.X = X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)


def test_load_packaged_database():
    enzyme, sensor = load_cell_mesh_database()
    assert not enzyme.empty
    assert not sensor.empty
    assert {"metabolite", "hmdb_id", "gene", "role"}.issubset(enzyme.columns)
    assert {"metabolite", "hmdb_id", "sensor_gene", "sensor_type"}.issubset(sensor.columns)
    assert enzyme["hmdb_id"].notna().all()
    assert sensor["hmdb_id"].notna().all()
    assert set(sensor["sensor_type"]).issubset({"Cell surface receptor", "Transporter", "Other receptor"})
    assert "weight" not in sensor.columns


def test_default_database_uses_highest_packaged_version():
    enzyme_path, interaction_path = _default_database_paths()

    assert enzyme_path.name == "Enzyme1.39.csv"
    assert interaction_path.name == "Interaction1.41.csv"


def test_versioned_database_discovery_ignores_noncanonical_names():
    class Entry:
        def __init__(self, name):
            self.name = name

    entries = [
        Entry("Enzyme1.0.csv"),
        Entry("Enzyme1.2.csv"),
        Entry("Enzyme2.2.csv"),
        Entry("Enzyme1.0_all3.csv"),
        Entry("enzyme_test.csv"),
    ]

    matches = _find_versioned_database_files(entries, "Enzyme")

    assert sorted(matches) == [(1, 0), (1, 2), (2, 2)]
    assert matches[(2, 2)].name == "Enzyme2.2.csv"


def test_default_database_selection_uses_independent_highest_versions():
    class Entry:
        def __init__(self, name):
            self.name = name

    entries = [
        Entry("Enzyme1.2.csv"),
        Entry("Enzyme2.2.csv"),
        Entry("Interaction1.0.csv"),
        Entry("Interaction1.9.csv"),
    ]

    enzyme_path, interaction_path = _select_default_database_paths(
        entries,
        Entry("enzyme_test.csv"),
        Entry("interaction_test.csv"),
    )

    assert enzyme_path.name == "Enzyme2.2.csv"
    assert interaction_path.name == "Interaction1.9.csv"


def test_load_database_accepts_lowercase_internal_schema(tmp_path):
    enzyme_path = tmp_path / "enzyme.csv"
    interaction_path = tmp_path / "interaction.csv"
    pd.DataFrame(
        {
            "metabolite": ["M1", "M1", "M2"],
            "hmdb_id": ["HMDB00001", "HMDB00001", "HMDB00002"],
            "reaction": ["R_prod", "R_sub", "R_exp"],
            "gene": ["G_prod", "G_sub", "G_exp"],
            "direction": ["product", "substrate", "exporter"],
        }
    ).to_csv(enzyme_path, index=False)
    pd.DataFrame(
        {
            "metabolite": ["M1", "M2"],
            "hmdb_id": ["HMDB00001", "HMDB00002"],
            "sensor_gene": ["S1", "S2"],
            "sensor_type": ["transporter", "intracellular receptor"],
            "weight": [2.0, 1.0],
        }
    ).to_csv(interaction_path, index=False)

    enzyme, sensor = load_cell_mesh_database(str(enzyme_path), str(interaction_path))

    assert set(enzyme["role"]) == {"production", "degradation", "export"}
    assert "weight" not in enzyme.columns
    assert "weight" not in sensor.columns
    assert set(sensor["sensor_type"]) == {"Transporter", "Other receptor"}
    assert {"metabolite", "hmdb_id", "gene", "role", "reaction"}.issubset(enzyme.columns)
    assert {"metabolite", "hmdb_id", "sensor_gene", "sensor_type"}.issubset(sensor.columns)


def test_normalize_interaction_database_ignores_weight_and_keeps_first_duplicate():
    raw = pd.DataFrame(
        {
            "metabolite": ["first", "second", "other"],
            "hmdb_id": ["HMDB00001", "HMDB00001", "HMDB00002"],
            "sensor_gene": ["S1", "S1", "S2"],
            "sensor_type": ["Transporter", "Cell surface receptor", "Transporter"],
            # The old implementation sorted this column and retained "second".
            "weight": [1.0, 100.0, 50.0],
        }
    )

    sensor = normalize_interaction_database(raw)

    assert "weight" not in sensor.columns
    assert list(sensor["metabolite"]) == ["first", "other"]
    assert list(sensor["sensor_type"]) == ["Transporter", "Transporter"]


def test_validate_priors_drops_legacy_sensor_weight():
    enzyme = pd.DataFrame(
        {
            "metabolite": ["M1"],
            "hmdb_id": ["HMDB00001"],
            "gene": ["E1"],
            "role": ["production"],
            "reaction": ["prod"],
        }
    )
    sensor = pd.DataFrame(
        {
            "metabolite": ["M1"],
            "hmdb_id": ["HMDB00001"],
            "sensor_gene": ["S1"],
            "sensor_type": ["Transporter"],
            "weight": [99.0],
        }
    )

    _, validated_sensor = validate_priors(enzyme, sensor, ["E1", "S1"])

    assert "weight" not in validated_sensor.columns


def test_validate_priors_deduplicates_canonical_hmdb_sensor_pairs():
    enzyme = pd.DataFrame(
        {
            "metabolite": ["M1"],
            "hmdb_id": ["HMDB00001"],
            "gene": ["E1"],
            "role": ["production"],
            "reaction": ["prod"],
        }
    )
    sensor = pd.DataFrame(
        {
            "metabolite": ["first name", "alias name"],
            "hmdb_id": [" hmdb00001 ", "HMDB00001"],
            "sensor_gene": ["S1", "S1"],
            "sensor_type": ["Transporter", "Transporter"],
            "source": ["first", "second"],
        }
    )

    _, validated_sensor = validate_priors(enzyme, sensor, ["E1", "S1"])

    assert validated_sensor[["hmdb_id", "sensor_gene"]].to_dict("records") == [
        {"hmdb_id": "HMDB00001", "sensor_gene": "S1"}
    ]
    assert validated_sensor.loc[0, "metabolite"] == "first name"
    assert validated_sensor.loc[0, "source"] == "first"


def test_validate_priors_rejects_sensor_type_conflicts_for_same_pair():
    enzyme = pd.DataFrame(
        {
            "metabolite": ["M1"],
            "hmdb_id": ["HMDB00001"],
            "gene": ["E1"],
            "role": ["production"],
            "reaction": ["prod"],
        }
    )
    sensor = pd.DataFrame(
        {
            "metabolite": ["M1", "M1"],
            "hmdb_id": ["hmdb00001", "HMDB00001"],
            "sensor_gene": ["S1", "S1"],
            "sensor_type": ["Transporter", "Cell surface receptor"],
        }
    )

    with pytest.raises(
        ValueError,
        match="conflicting sensor_type.*HMDB00001.*S1",
    ):
        validate_priors(enzyme, sensor, ["E1", "S1"])


@pytest.mark.parametrize(
    "enzyme, message",
    [
        (
            pd.DataFrame(
                {
                    "metabolite": ["M"],
                    "hmdb_id": ["H"],
                    "gene": ["E"],
                    "role": ["production"],
                }
            ),
            "reaction",
        ),
        (
            pd.DataFrame(
                {
                    "metabolite": ["M"],
                    "hmdb_id": ["H"],
                    "gene": ["E"],
                    "role": ["production"],
                    "reaction": ["   "],
                }
            ),
            "reaction values must be non-empty",
        ),
    ],
)
def test_validate_priors_requires_nonempty_reaction(enzyme, message):
    sensor = pd.DataFrame(
        {
            "metabolite": ["M"],
            "hmdb_id": ["H"],
            "sensor_gene": ["S"],
            "sensor_type": ["Transporter"],
        }
    )

    with pytest.raises(ValueError, match=message):
        validate_priors(enzyme, sensor, ["E", "S"])


@pytest.mark.parametrize("enzyme_kind", ["default", "str", "path", "dataframe"])
@pytest.mark.parametrize("sensor_kind", ["default", "str", "path", "dataframe"])
def test_run_cell_mesh_with_packaged_database(enzyme_kind, sensor_kind, monkeypatch):
    enzyme, sensor = load_cell_mesh_database()
    enzyme_path, sensor_path = _default_database_paths()
    enzyme_input = {"default": None, "str": str(enzyme_path), "path": enzyme_path, "dataframe": enzyme}[enzyme_kind]
    sensor_input = {"default": None, "str": str(sensor_path), "path": sensor_path, "dataframe": sensor}[sensor_kind]
    if enzyme_kind != "default" and sensor_kind != "default":
        def unexpected_defaults():
            pytest.fail("Explicit inputs must not depend on packaged databases")

        monkeypatch.setattr("cellmesh.database._default_database_paths", unexpected_defaults)
    # Use a matched metabolite with production evidence so the synthetic sender
    # expression creates a directional A -> B event.
    enzyme_production = enzyme[enzyme["role"] == "production"]
    pairs = set(zip(enzyme_production["metabolite"], enzyme_production["hmdb_id"])).intersection(zip(sensor["metabolite"], sensor["hmdb_id"]))
    for met, hmdb_id in sorted(pairs):
        e_genes = enzyme_production.loc[
            (enzyme_production["metabolite"] == met) & (enzyme_production["hmdb_id"] == hmdb_id),
            "gene",
        ]
        s_genes = sensor.loc[
            (sensor["metabolite"] == met) & (sensor["hmdb_id"] == hmdb_id),
            "sensor_gene",
        ]
        candidates = [(e_gene, s_gene) for e_gene in e_genes for s_gene in s_genes if e_gene != s_gene]
        if candidates:
            e_gene, s_gene = candidates[0]
            break
    else:
        pytest.fail("No packaged metabolite has distinct production and sensor genes")
    genes = [e_gene, s_gene, "BACKGROUND"]
    X = np.array([
        [5, 0, 0], [4, 0, 0], [5, 0, 0],
        [0, 4, 0], [0, 5, 0], [0, 4, 0],
    ], dtype=float)
    adata = FakeAnnData(X, genes, {"cell_type": ["A", "A", "A", "B", "B", "B"]})
    expected = run_cell_mesh(
        adata, enzyme, sensor, cell_type_key="cell_type", min_cells=2, allow_self=False
    )
    kwargs = {}
    if enzyme_kind != "default":
        kwargs["enzyme_metabolite"] = enzyme_input
    if sensor_kind != "default":
        kwargs["metabolite_sensor"] = sensor_input
    res = run_cell_mesh(adata, cell_type_key="cell_type", min_cells=2, allow_self=False, **kwargs)
    pd.testing.assert_frame_equal(res.events, expected.events)
    assert not res.events.empty
    assert "cell_mesh_score" in res.events.columns
    assert res.events.iloc[0]["sender"] == "A"
    assert res.events.iloc[0]["receiver"] == "B"


@pytest.mark.parametrize("parameter", ["enzyme_metabolite", "metabolite_sensor"])
def test_run_cell_mesh_rejects_invalid_prior_type(parameter):
    adata = FakeAnnData(np.ones((2, 1)), ["G"], {"cell_type": ["A", "B"]})
    with pytest.raises(TypeError, match=parameter + " must be"):
        run_cell_mesh(adata, **{parameter: 42})


def test_load_database_preserves_dataframes_and_reports_missing_csv(tmp_path):
    enzyme = pd.DataFrame({"gene": ["G"]})
    sensor = pd.DataFrame({"sensor_gene": ["S"]})
    loaded_enzyme, loaded_sensor = load_cell_mesh_database(enzyme, sensor)
    assert loaded_enzyme is enzyme
    assert loaded_sensor is sensor
    with pytest.raises(FileNotFoundError):
        load_cell_mesh_database(tmp_path / "missing.csv", sensor)


def test_read_csv_with_metadata(tmp_path):
    pytest.importorskip("anndata")

    expr_path = tmp_path / "expr.csv"
    cell_meta_path = tmp_path / "cells.csv"
    gene_meta_path = tmp_path / "genes.csv"
    pd.DataFrame({"G1": [1, 0], "G2": [0, 2]}, index=["C1", "C2"]).to_csv(expr_path)
    pd.DataFrame({"cell_id": ["C1", "C2"], "cell_type": ["A", "B"]}).to_csv(cell_meta_path, index=False)
    pd.DataFrame({"gene_id": ["G1", "G2"], "symbol": ["Gene1", "Gene2"]}).to_csv(gene_meta_path, index=False)

    adata = read_anndata(
        expr_path,
        mode="csv",
        cell_meta_path=cell_meta_path,
        gene_meta_path=gene_meta_path,
        cell_id_col="cell_id",
    )

    assert list(adata.obs_names) == ["C1", "C2"]
    assert list(adata.var_names) == ["G1", "G2"]
    assert adata.obs.loc["C2", "cell_type"] == "B"
    assert adata.var.loc["G1", "symbol"] == "Gene1"


def test_read_mtx_with_names(tmp_path):
    pytest.importorskip("anndata")
    scipy_io = pytest.importorskip("scipy.io")
    sparse = pytest.importorskip("scipy.sparse")

    mtx_path = tmp_path / "matrix.mtx"
    genes_path = tmp_path / "features.tsv"
    barcodes_path = tmp_path / "barcodes.tsv"
    scipy_io.mmwrite(mtx_path, sparse.csr_matrix(np.array([[1, 0], [0, 2]], dtype=float)))
    pd.DataFrame([["gene_id_1", "G1"], ["gene_id_2", "G2"]]).to_csv(genes_path, sep="\t", header=False, index=False)
    pd.Series(["C1", "C2"]).to_csv(barcodes_path, sep="\t", header=False, index=False)

    adata = read_anndata(mtx_path, mode="mtx", genes_path=genes_path, barcodes_path=barcodes_path)

    assert list(adata.obs_names) == ["C1", "C2"]
    assert list(adata.var_names) == ["G1", "G2"]
    assert adata.X.shape == (2, 2)

    with pytest.raises(TypeError, match="unexpected_keyword"):
        read_anndata(
            mtx_path,
            mode="mtx",
            genes_path=genes_path,
            barcodes_path=barcodes_path,
            unexpected_keyword=True,
        )
