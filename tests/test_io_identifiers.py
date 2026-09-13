"""Text readers must preserve identity before pandas can infer or rename it."""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite

from cellmesh import read_anndata, run_cell_mesh


def _write_rows(path, rows, *, separator=",", compress=False):
    opener = gzip.open if compress else open
    with opener(path, "wt", encoding="utf-8-sig", newline="") as stream:
        csv.writer(stream, delimiter=separator).writerows(rows)
    return path


@pytest.mark.parametrize("mode,separator", [("csv", ","), ("tsv", "\t")])
@pytest.mark.parametrize("transpose", [False, True])
@pytest.mark.parametrize("spaced", [False, True])
def test_duplicate_expression_genes_fail_before_scoring(tmp_path, mode, separator, transpose, spaced):
    other = " PROD " if spaced else "PROD"
    rows = [["cell", "PROD", other, "SENSOR"], ["C1", 8, 0, 1], ["C2", 0, 8, 1]]
    if transpose:
        rows = list(map(list, zip(*rows)))
    path = _write_rows(tmp_path / f"expression.{mode}", rows, separator=separator)
    with pytest.raises(ValueError, match="duplicate identifiers.*PROD.*positions"):
        read_anndata(path, mode=mode, transpose=transpose)


@pytest.mark.parametrize("options", [{}, {"usecols": [0, 1]}, {"names": ["cell", "P1", "P2"], "header": 0}])
def test_raw_duplicate_header_cannot_be_hidden_by_selection_or_renaming(tmp_path, options):
    path = _write_rows(tmp_path / "expression.csv", [["cell", "PROD", "PROD"], ["C1", 8, 0]])
    with pytest.raises(ValueError, match="header.*duplicate identifiers.*PROD"):
        read_anndata(path, mode="csv", **options)


def test_literal_dot_suffix_is_a_valid_gene_name(tmp_path):
    path = _write_rows(tmp_path / "expression.csv", [["cell", "PROD", "PROD.1", "SENSOR"], ["C1", 8, 3, 1]])
    data = read_anndata(path, mode="csv")
    assert data.var_names.tolist() == ["PROD", "PROD.1", "SENSOR"]
    np.testing.assert_array_equal(data.X, [[8, 3, 1]])


@pytest.mark.parametrize("mode,separator", [("csv", ","), ("tsv", "\t")])
@pytest.mark.parametrize("transpose", [False, True])
@pytest.mark.parametrize("compress", [False, True])
@pytest.mark.parametrize("engine", ["c", "python"])
def test_expression_axis_text_survives_reading_with_numeric_values(
    tmp_path, mode, separator, transpose, compress, engine,
):
    cells = ["01", "1", "NA", "nan"]
    genes = ["0007", "7", "N/A", "NULL"]
    values = np.arange(16).reshape(4, 4)
    rows = [[""] + genes] + [[cell] + values[i].tolist() for i, cell in enumerate(cells)]
    if transpose:
        rows = list(map(list, zip(*rows)))
    suffix = f".{mode}" + (".gz" if compress else "")
    path = _write_rows(tmp_path / ("expression" + suffix), rows, separator=separator, compress=compress)
    data = read_anndata(path, mode=mode, transpose=transpose, dtype=np.float32, engine=engine)
    assert data.obs_names.tolist() == cells
    assert data.var_names.tolist() == genes
    assert data.X.dtype == np.float32
    np.testing.assert_array_equal(data.X, values)


@pytest.mark.parametrize("options,contents", [
    ({"header": 1}, "preamble\ncell,G1,G2\n01,1,2\n"),
    ({"skiprows": 1}, "preamble\ncell,G1,G2\n01,1,2\n"),
    ({"skiprows": [0], "header": 1}, "ignored\npreamble\ncell,G1,G2\n01,1,2\n"),
    ({"comment": "#"}, "# comment\n\ncell,G1,G2\n01,1,2\n"),
    ({"header": None, "names": ["cell", "G1", "G2"]}, "01,1,2\n"),
    ({"header": 0, "names": ["cell", "G1", "G2"]}, "old,P1,P2\n01,1,2\n"),
    ({"delimiter": ";"}, 'cell;G1;G2\n"01";1;2\n'),
    ({"sep": None, "engine": "python"}, 'cell;G1;G2\n"01";1;2\n'),
    ({"skipfooter": 1, "engine": "python"}, "cell,G1,G2\n01,1,2\nfooter\n"),
])
def test_common_parser_options_preserve_headers_and_identifiers(tmp_path, options, contents):
    path = tmp_path / "expression.csv"
    path.write_text(contents)
    data = read_anndata(path, mode="csv", **options)
    assert data.obs_names.tolist() == ["01"]
    assert data.var_names.tolist() == ["G1", "G2"]
    np.testing.assert_array_equal(data.X, [[1, 2]])


@pytest.mark.parametrize("usecols", [[0, 2], ["cell", "G2"], lambda name: name != "G1"])
def test_expression_column_selection_and_dtype_remain_available(tmp_path, usecols):
    path = _write_rows(tmp_path / "expression.csv", [["cell", "G1", "G2"], ["01", 1, 2]])
    data = read_anndata(path, mode="csv", usecols=usecols, dtype={"G2": np.float32})
    assert data.obs_names.tolist() == ["01"]
    assert data.var_names.tolist() == ["G2"]
    assert data.X.dtype == np.float32
    np.testing.assert_array_equal(data.X, [[2]])


def test_expression_converter_still_applies_to_values(tmp_path):
    path = _write_rows(tmp_path / "expression.csv", [["cell", "G1"], ["NA", "2"]])
    data = read_anndata(path, mode="csv", converters={"G1": lambda value: float(value) * 2})
    assert data.obs_names.tolist() == ["NA"]
    np.testing.assert_array_equal(data.X, [[4.]])
    with pytest.raises(ValueError, match="converter.*identifier"):
        read_anndata(path, mode="csv", converters={"cell": int})


@pytest.mark.parametrize("rows,match", [
    ([["cell", "G1", ""], ["C1", 1, 2]], "header.*empty identifiers"),
    ([["cell", "G1"], ["", 1]], "cell identifiers.*empty identifiers"),
    ([["cell", "G1"], ["C1", 1], [" C1 ", 2]], "cell identifiers.*duplicate identifiers"),
    ([["cell", "G1"], ["C1", 1, 2]], "more fields than the header"),
])
def test_ambiguous_or_empty_matrix_identifiers_are_rejected(tmp_path, rows, match):
    path = _write_rows(tmp_path / "expression.csv", rows)
    with pytest.raises(ValueError, match=match):
        read_anndata(path, mode="csv")


def test_metadata_text_and_axis_alignment_preserve_literal_names(tmp_path):
    expression = _write_rows(tmp_path / "expression.csv", [["", "01", "1", "NA"], ["01", 1, 2, 3], ["1", 4, 5, 6], ["NA", 7, 8, 9]])
    cells = _write_rows(tmp_path / "cells.csv", [
        ["cell_id", "cell_type", "donor_code", "age", "note"],
        ["NA", "nan", "NA", "20", "NA"], ["1", "1", "1", "08", ""],
        ["01", "01", "01", "07", "NULL"],
    ])
    genes = _write_rows(tmp_path / "genes.csv", [["gene_id", "reference"], ["NA", "NA"], ["1", "1"], ["01", "01"]])
    data = read_anndata(expression, mode="csv", cell_meta_path=cells, gene_meta_path=genes, cell_id_col="cell_id")
    assert data.obs_names.tolist() == ["01", "1", "NA"]
    assert data.obs.cell_type.tolist() == ["01", "1", "nan"]
    assert data.obs.donor_code.tolist() == ["01", "1", "NA"]
    assert data.obs.age.tolist() == ["07", "08", "20"]
    assert data.obs.age.astype(int).tolist() == [7, 8, 20]
    assert pd.isna(data.obs.loc["1", "note"])
    assert data.obs.loc["NA", "note"] == "NA"
    assert data.var.reference.tolist() == ["01", "1", "NA"]
    assert data.X.dtype.kind in "iuf"


@pytest.mark.parametrize("rows,match", [
    ([["cell_id", "cell_type", "cell_type"], ["C1", "A", "B"]], "header.*duplicate identifiers"),
    ([["cell_id", "cell_type", " cell_type "], ["C1", "A", "B"]], "header.*duplicate identifiers"),
    ([["cell_id", "cell_type"], ["C1", "A"], [" C1 ", "B"]], "row identifiers.*duplicate identifiers"),
    ([["cell_id", "cell_type"], ["", "A"]], "row identifiers.*empty identifiers"),
])
def test_metadata_conflicts_fail_before_alignment(tmp_path, rows, match):
    expression = _write_rows(tmp_path / "expression.csv", [["cell", "G1"], ["C1", 1]])
    cells = _write_rows(tmp_path / "cells.csv", rows)
    with pytest.raises(ValueError, match=match):
        read_anndata(expression, mode="csv", cell_meta_path=cells)


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
def test_literal_sample_ids_match_native_anndata_scores_and_permutations(tmp_path, mode):
    expression = _write_rows(tmp_path / "expression.csv", [
        ["cell", "PROD", "SENSOR"], ["C1", 8, 1], ["C2", 0, 1], ["C3", 0, 1], ["C4", 8, 1],
    ])
    cells = _write_rows(tmp_path / "cells.csv", [
        ["cell_id", "cell_type", "donor_code"], ["C1", "A", "01"], ["C2", "B", "01"],
        ["C3", "A", "1"], ["C4", "B", "1"],
    ])
    loaded = read_anndata(expression, mode="csv", cell_meta_path=cells)
    native = ad.AnnData(
        np.array([[8, 1], [0, 1], [0, 1], [8, 1]]),
        obs=pd.DataFrame({"cell_type": ["A", "B", "A", "B"], "donor_code": ["01", "01", "1", "1"]}, index=["C1", "C2", "C3", "C4"]),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme = pd.DataFrame([{"metabolite": "M", "hmdb_id": "HMDB0000001", "gene": "PROD", "role": "production", "reaction": "R1"}])
    sensor = pd.DataFrame([{"metabolite": "M", "hmdb_id": "HMDB0000001", "sensor_gene": "SENSOR", "sensor_type": "Transporter"}])
    options = dict(sample_key="donor_code", sample_mode=mode, min_cells=1, n_perms=5, random_state=42)
    actual = run_cell_mesh(loaded, enzyme, sensor, **options)
    expected = run_cell_mesh(native, enzyme, sensor, **options)
    pd.testing.assert_frame_equal(actual.events, expected.events)
    if mode == "sample_aware":
        assert list(actual.availability_results["availability_by_sample"]) == ["01", "1"]
        assert actual.events.n_samples_coobserved.eq(2).all()
        score = actual.events.query("sender == 'A' and receiver == 'B'").cell_mesh_score.iloc[0]
        assert score == pytest.approx(0.23717082451262844)


@pytest.mark.parametrize("extension,separator", [("tsv", "\t"), ("csv", ","), ("txt", "\t")])
@pytest.mark.parametrize("compress", [False, True])
def test_mtx_name_files_preserve_text_and_compressed_formats(tmp_path, extension, separator, compress):
    matrix = tmp_path / "matrix.mtx"
    mmwrite(matrix, sparse.csr_matrix(np.arange(9).reshape(3, 3)))
    suffix = "." + extension + (".gz" if compress else "")
    genes = _write_rows(tmp_path / ("features" + suffix), [["id1", "0007"], ["id2", "7"], ["id3", "NA"]], separator=separator, compress=compress)
    cells = _write_rows(tmp_path / ("barcodes" + suffix), [["01"], ["1"], ["NA"]], separator=separator, compress=compress)
    data = read_anndata(matrix, mode="mtx", genes_path=genes, barcodes_path=cells)
    assert data.obs_names.tolist() == ["01", "1", "NA"]
    assert data.var_names.tolist() == ["0007", "7", "NA"]
    np.testing.assert_array_equal(data.X.toarray(), np.arange(9).reshape(3, 3).T)


@pytest.mark.parametrize("names", [["G1", "G1"], ["G1", " G1 "], ["G1", ""]])
def test_mtx_empty_or_duplicate_names_are_not_renamed(tmp_path, names):
    matrix = tmp_path / "matrix.mtx"
    mmwrite(matrix, sparse.eye(2, format="csr"))
    genes = _write_rows(tmp_path / "genes.tsv", [[name] for name in names], separator="\t")
    with pytest.raises(ValueError, match="Name file.*(duplicate|empty) identifiers"):
        read_anndata(matrix, mode="mtx", genes_path=genes)
