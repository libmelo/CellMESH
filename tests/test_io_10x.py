"""10X imports preserve raw identifiers and matrices across file layouts."""
import csv
import gzip

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse
from scipy.io import mmwrite

from cellmesh import read_anndata, run_cell_mesh


@pytest.fixture(autouse=True)
def scanpy_available(tmp_path, monkeypatch):
    # Scanpy/Numba needs a writable cache when imported in the sandbox.
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba"))
    pytest.importorskip("scanpy")


def _write_table(path, rows):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8-sig", newline="") as stream:
        csv.writer(stream, delimiter="\t").writerows(rows)


def _case(tmp_path, *, modern=True, compressed=False, prefix="", features=None, barcodes=None, matrix=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    suffix = ".gz" if compressed else ""
    if matrix is None:
        matrix = np.array([[1, 2, 3], [10, 20, 30]])
    if features is None:
        features = [["0001", "G"], ["0002", "G-1"]]
        if modern:
            features = [row + ["Gene Expression"] for row in features]
    if barcodes is None:
        barcodes = ["01", "1", "NA"]
    feature_path = tmp_path / (prefix + ("features.tsv" if modern else "genes.tsv") + suffix)
    barcode_path = tmp_path / (prefix + "barcodes.tsv" + suffix)
    matrix_path = tmp_path / (prefix + "matrix.mtx" + suffix)
    _write_table(feature_path, features)
    _write_table(barcode_path, [[value] for value in barcodes])
    opener = gzip.open if compressed else open
    with opener(matrix_path, "wb") as stream:
        mmwrite(stream, sparse.coo_matrix(matrix), symmetry="general")
    return feature_path, barcode_path, matrix_path


@pytest.mark.parametrize("modern", [False, True])
@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("var_names", ["gene_symbols", "gene_ids"])
@pytest.mark.parametrize("prefix", ["", "patient_"])
def test_layouts_preserve_values_literal_ids_and_legal_suffixes(tmp_path, modern, compressed, var_names, prefix):
    _case(tmp_path, modern=modern, compressed=compressed, prefix=prefix)
    data = read_anndata(tmp_path, mode="10x", var_names=var_names, prefix=prefix)
    assert sparse.issparse(data.X)
    np.testing.assert_array_equal(data.X.toarray(), [[1, 10], [2, 20], [3, 30]])
    assert data.obs_names.tolist() == ["01", "1", "NA"]
    if var_names == "gene_symbols":
        assert data.var_names.tolist() == ["G", "G-1"]
        assert data.var.gene_ids.tolist() == ["0001", "0002"]
    else:
        assert data.var_names.tolist() == ["0001", "0002"]
        assert data.var.gene_symbols.tolist() == ["G", "G-1"]


@pytest.mark.parametrize("modern", [False, True])
@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("padded", [False, True])
def test_duplicate_symbols_fail_with_original_ids_and_positions(tmp_path, modern, compressed, padded):
    features = [["id1", "G"], ["id2", " G " if padded else "G"]]
    if modern:
        features = [row + ["Gene Expression"] for row in features]
    paths = _case(tmp_path, modern=modern, compressed=compressed, features=features)
    with pytest.raises(ValueError) as error:
        read_anndata(tmp_path, mode="10x")
    message = str(error.value)
    for text in [str(paths[0]), "duplicate", "id1", "id2", "G", "1", "2"]:
        assert text in message
    # The same source is unambiguous when explicitly matched by unique gene IDs.
    data = read_anndata(tmp_path, mode="10x", var_names="gene_ids", make_unique=False)
    assert data.var_names.tolist() == ["id1", "id2"]
    np.testing.assert_array_equal(data.X.toarray(), [[1, 10], [2, 20], [3, 30]])


@pytest.mark.parametrize("axis", ["symbol", "gene_id", "barcode"])
@pytest.mark.parametrize("problem", ["duplicate", "empty"])
def test_selected_axes_reject_empty_or_duplicate_identifiers(tmp_path, axis, problem):
    features = [["id1", "G1", "Gene Expression"], ["id2", "G2", "Gene Expression"]]
    barcodes = ["C1", "C2", "C3"]
    value = "" if problem == "empty" else {"symbol": " G1 ", "gene_id": " id1 ", "barcode": " C1 "}[axis]
    if axis == "barcode":
        barcodes[1] = value
    else:
        features[1][0 if axis == "gene_id" else 1] = value
    _case(tmp_path, features=features, barcodes=barcodes)
    with pytest.raises(ValueError, match="empty|duplicate"):
        read_anndata(tmp_path, mode="10x", var_names="gene_ids" if axis == "gene_id" else "gene_symbols")


def test_gex_filter_validates_selected_features_before_constructing_anndata(tmp_path):
    features = [["id1", "G", "Gene Expression"], ["ab1", "G", "Antibody Capture"]]
    _case(tmp_path, features=features)
    data = read_anndata(tmp_path, mode="10x")
    assert data.var_names.tolist() == ["G"]
    np.testing.assert_array_equal(data.X.toarray(), [[1], [2], [3]])
    with pytest.raises(ValueError, match="duplicate"):
        read_anndata(tmp_path, mode="10x", gex_only=False)


@pytest.mark.parametrize("part", ["features", "barcodes"])
@pytest.mark.parametrize("problem", ["short", "wide", "count"])
def test_metadata_width_and_matrix_dimensions_are_checked(tmp_path, part, problem):
    features = [["id1", "G1", "Gene Expression"], ["id2", "G2", "Gene Expression"]]
    feature_path, barcode_path, _ = _case(tmp_path, features=features)
    path = feature_path if part == "features" else barcode_path
    rows = features if part == "features" else [["C1"], ["C2"], ["C3"]]
    if problem == "short":
        rows[-1] = rows[-1][:-1]
    elif problem == "wide":
        rows[-1] += ["EXTRA"]
    else:
        rows.pop()
    _write_table(path, rows)
    with pytest.raises(ValueError, match="expected|dimensions"):
        read_anndata(tmp_path, mode="10x")


def test_mixed_compression_missing_and_ambiguous_files(tmp_path):
    features, barcodes, matrix = _case(tmp_path)
    with features.open("rb") as source, gzip.open(str(features) + ".gz", "wb") as target:
        target.write(source.read())
    with pytest.raises(ValueError, match="ambiguous"):
        read_anndata(tmp_path, mode="10x")
    features.unlink()
    assert read_anndata(tmp_path, mode="10x").shape == (3, 2)
    matrix.unlink()
    with pytest.raises(FileNotFoundError, match="matrix"):
        read_anndata(tmp_path, mode="10x")


@pytest.mark.parametrize("kwargs", [{"make_unique": True}, {"make_unique": "False"},
                                   {"gex_only": "False"}, {"var_names": "bad"}, {"cache": "False"},
                                   {"unknown_argument": True}])
def test_options_cannot_reenable_silent_input_changes(tmp_path, kwargs):
    _case(tmp_path)
    with pytest.raises((ValueError, TypeError)):
        read_anndata(tmp_path, mode="10x", **kwargs)


def test_cache_does_not_bypass_current_raw_identifier_validation(tmp_path):
    import scanpy as sc
    paths = _case(tmp_path / "source", compressed=True)
    previous_cache = sc.settings.cachedir
    try:
        sc.settings.cachedir = tmp_path / "cache"
        data = read_anndata(tmp_path / "source", mode="10x", cache=True, cache_compression="gzip")
        assert data.var_names.tolist() == ["G", "G-1"]
        cached = read_anndata(tmp_path / "source", mode="10x", cache=True)
        np.testing.assert_array_equal(cached.X.toarray(), data.X.toarray())
        pd.testing.assert_frame_equal(cached.var, data.var)
        pd.testing.assert_frame_equal(cached.obs, data.obs)
        _write_table(paths[0], [["id1", "G", "Gene Expression"], ["id2", "G", "Gene Expression"]])
        with pytest.raises(ValueError, match="duplicate"):
            read_anndata(tmp_path / "source", mode="10x", cache=True)
    finally:
        sc.settings.cachedir = previous_cache


def test_no_gene_expression_features_returns_an_empty_selected_axis(tmp_path):
    _case(tmp_path, features=[["id1", "A", "Antibody Capture"], ["id2", "B", "Antibody Capture"]])
    data = read_anndata(tmp_path, mode="10x")
    assert data.shape == (3, 0)
    assert data.obs_names.tolist() == ["01", "1", "NA"]
    assert data.var.columns.tolist() == ["gene_ids", "feature_types"]
    assert read_anndata(tmp_path, mode="10x", gex_only=False).shape == (3, 2)


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("modern,compressed", [(False, False), (False, True), (True, False), (True, True)])
@pytest.mark.parametrize("n_perms", [0, 5])
def test_10x_scoring_matches_native_input(tmp_path, mode, modern, compressed, n_perms):
    matrix = np.array([[8, 1, 6, 0], [1, 4, 2, 5]], dtype=np.float32)
    features = [["id1", "PROD"], ["id2", "R"]]
    if modern:
        features = [row + ["Gene Expression"] for row in features]
    barcodes = ["C1", "C2", "C3", "C4"]
    _case(tmp_path, modern=modern, compressed=compressed, features=features, barcodes=barcodes, matrix=matrix)
    loaded = read_anndata(tmp_path, mode="10x")
    obs = pd.DataFrame({"cell_type": ["A", "B", "A", "B"], "sample": ["01", "01", "1", "1"]}, index=barcodes)
    loaded.obs = obs.copy()
    native = AnnData(matrix.T, obs=obs.copy(), var=pd.DataFrame(index=["PROD", "R"]))
    enzyme = pd.DataFrame([dict(metabolite="M", hmdb_id="HMDB1", gene="PROD", role="production", reaction="P")])
    sensor = pd.DataFrame([dict(metabolite="M", hmdb_id="HMDB1", sensor_gene="R", sensor_type="Transporter")])
    options = dict(sample_key="sample", sample_mode=mode, min_cells=1, n_perms=n_perms, random_state=7, store_null_scores=True)
    actual = run_cell_mesh(loaded, enzyme, sensor, **options)
    expected = run_cell_mesh(native, enzyme, sensor, **options)
    pd.testing.assert_frame_equal(actual.events, expected.events, check_exact=True)
    if mode == "sample_aware":
        pd.testing.assert_frame_equal(actual.events.attrs["sample_aware_null_scores"],
                                      expected.events.attrs["sample_aware_null_scores"], check_exact=True)
