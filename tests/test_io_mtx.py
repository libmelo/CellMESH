"""Read both Matrix Market storage formats without changing matrix semantics."""
import gzip

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse
from scipy.io import mmwrite

from cellmesh import read_anndata, run_cell_mesh


def _write_matrix(tmp_path, matrix, storage, *, compressed=False):
    path = tmp_path / (storage + ".mtx" + (".gz" if compressed else ""))
    value = sparse.csr_matrix(matrix) if storage == "coordinate" else matrix
    opener = gzip.open if compressed else open
    with opener(path, "wb") as stream:
        mmwrite(stream, value, symmetry="general")
    with opener(path, "rb") as stream:
        assert stream.readline().startswith(f"%%MatrixMarket matrix {storage} ".encode())
    return path


def _write_names(path, names):
    path.write_text("".join(name + "\n" for name in names), encoding="utf-8")
    return path


def _array(matrix):
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


@pytest.mark.parametrize("dtype", [np.int64, np.float64], ids=["integer", "real"])
@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
def test_mtx_formats_preserve_values_transpose_and_literal_names(tmp_path, dtype, compressed):
    matrix = np.array([[0, 4, 9], [3, 0, 7]], dtype=dtype)
    if dtype is np.float64:
        matrix /= 4
    genes = _write_names(tmp_path / "features.tsv", ["id1\t0007", "id2\t7"])
    barcodes = _write_names(tmp_path / "barcodes.tsv", ["01", "1", "NA"])
    loaded = []
    for storage in ("coordinate", "array"):
        path = _write_matrix(tmp_path, matrix, storage, compressed=compressed)
        data = read_anndata(path, mode="mtx", genes_path=genes, barcodes_path=barcodes)
        assert data.shape == (3, 2)
        assert data.X.dtype == matrix.dtype
        assert sparse.issparse(data.X) == (storage == "coordinate")
        if storage == "array":
            assert isinstance(data.X, np.ndarray)
        np.testing.assert_array_equal(_array(data.X), matrix.T)
        assert data.var_names.tolist() == ["0007", "7"]
        assert data.obs_names.tolist() == ["01", "1", "NA"]
        loaded.append(data)
    pd.testing.assert_frame_equal(loaded[0].obs, loaded[1].obs)
    pd.testing.assert_frame_equal(loaded[0].var, loaded[1].var)


@pytest.mark.parametrize("storage", ["coordinate", "array"])
@pytest.mark.parametrize("shape", [(1, 3), (3, 1)])
def test_mtx_single_axis_inputs_stay_two_dimensional_without_name_files(tmp_path, storage, shape):
    matrix = np.arange(3, dtype=float).reshape(shape)
    path = _write_matrix(tmp_path, matrix, storage)
    data = read_anndata(path, mode="mtx")
    assert data.shape == shape[::-1]
    assert data.X.ndim == 2
    np.testing.assert_array_equal(_array(data.X), matrix.T)
    assert data.obs_names.tolist() == [str(i) for i in range(shape[1])]
    assert data.var_names.tolist() == [str(i) for i in range(shape[0])]


@pytest.mark.parametrize("axis", ["genes", "barcodes"])
@pytest.mark.parametrize("problem", ["length", "duplicate", "empty", "empty_first"])
def test_array_mtx_retains_identifier_validation(tmp_path, axis, problem):
    matrix = np.array([[0., 4., 9.], [3., 0., 7.]])
    path = _write_matrix(tmp_path, matrix, "array")
    names = ["G1", "G2"] if axis == "genes" else ["C1", "C2", "C3"]
    if problem == "length":
        names.pop()
        error = rf"{axis}_path has .* entries, expected"
    elif problem == "duplicate":
        names[1] = names[0] + " "
        error = "Name file.*duplicate identifiers"
    elif problem == "empty":
        names[-1] = ""
        error = "Name file.*empty identifiers"
    else:
        # With no first-row fields, pandas can reject the name file before
        # label validation. Both parser and label errors must still reject it.
        names[0] = ""
        error = None
    names_path = _write_names(tmp_path / (axis + ".tsv"), names)
    with pytest.raises(ValueError, match=error):
        read_anndata(path, mode="mtx", **{axis + "_path": names_path})


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_both_mtx_formats_match_native_anndata_scores_and_permutations(tmp_path, sample_mode):
    matrix = np.array([[8., 0., 6., 0.], [1., 4., 2., 5.], [1., 2., 1., 2.]])
    obs = pd.DataFrame({"cell_type": ["A", "B", "A", "B"], "sample": ["01", "01", "1", "1"]},
                       index=["C1", "C2", "C3", "C4"])
    var = pd.DataFrame(index=["PROD", "CONS", "SENSOR"])
    native = AnnData(matrix.T, obs=obs, var=var)
    genes = _write_names(tmp_path / "genes.tsv", var.index.tolist())
    barcodes = _write_names(tmp_path / "barcodes.tsv", obs.index.tolist())
    enzyme = pd.DataFrame([
        dict(metabolite="M", hmdb_id="HMDB0000001", gene="PROD", role="production", reaction="RP"),
        dict(metabolite="M", hmdb_id="HMDB0000001", gene="CONS", role="degradation", reaction="RC"),
    ])
    sensor = pd.DataFrame([
        dict(metabolite="M", hmdb_id="HMDB0000001", sensor_gene="SENSOR", sensor_type="Transporter"),
    ])
    options = dict(sample_mode=sample_mode, sample_key="sample", min_cells=1, n_perms=5, random_state=17)
    expected = run_cell_mesh(native, enzyme, sensor, **options)
    for storage in ("coordinate", "array"):
        path = _write_matrix(tmp_path, matrix, storage)
        data = read_anndata(path, mode="mtx", genes_path=genes, barcodes_path=barcodes)
        data.obs = obs.copy()
        actual = run_cell_mesh(data, enzyme, sensor, **options)
        pd.testing.assert_frame_equal(actual.events, expected.events)
        pd.testing.assert_frame_equal(actual.sender_scores, expected.sender_scores)
        pd.testing.assert_frame_equal(actual.receiver_scores, expected.receiver_scores)
        if sample_mode == "sample_aware":
            pd.testing.assert_frame_equal(actual.sample_events, expected.sample_events)
