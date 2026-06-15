import numpy as np
import pandas as pd
import pytest

from cellmesh import load_cell_mesh_database, read_anndata, run_cell_mesh


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


def test_run_cell_mesh_with_packaged_database():
    enzyme, sensor = load_cell_mesh_database()
    # Use a matched metabolite with available enzyme and sensor genes.
    pairs = set(zip(enzyme["metabolite"], enzyme["hmdb_id"])).intersection(zip(sensor["metabolite"], sensor["hmdb_id"]))
    met, hmdb_id = sorted(pairs)[0]
    e_gene = enzyme.loc[(enzyme["metabolite"] == met) & (enzyme["hmdb_id"] == hmdb_id), "gene"].iloc[0]
    s_gene = sensor.loc[(sensor["metabolite"] == met) & (sensor["hmdb_id"] == hmdb_id), "sensor_gene"].iloc[0]
    genes = [e_gene, s_gene, "BACKGROUND"]
    X = np.array([
        [5, 0, 0], [4, 0, 0], [5, 0, 0],
        [0, 4, 0], [0, 5, 0], [0, 4, 0],
    ], dtype=float)
    adata = FakeAnnData(X, genes, {"cell_type": ["A", "A", "A", "B", "B", "B"]})
    res = run_cell_mesh(adata, cell_type_key="cell_type", min_cells=2, allow_self=False)
    assert not res.events.empty
    assert "cell_mesh_score" in res.events.columns
    assert res.events.iloc[0]["sender"] == "A"
    assert res.events.iloc[0]["receiver"] == "B"


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
