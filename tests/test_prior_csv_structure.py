"""Malformed prior records must fail before pandas can shift or pad fields."""
import csv
import io

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import load_cell_mesh_database, run_cell_mesh


@pytest.fixture
def priors():
    enzyme = pd.DataFrame([dict(metabolite="M", hmdb_id="HMDB1", gene="G",
                                role="production", reaction="001")])
    sensor = pd.DataFrame([dict(metabolite="M", hmdb_id="HMDB1", sensor_gene="R",
                                sensor_type="Transporter")])
    return enzyme, sensor


def _csv(rows):
    stream = io.StringIO(newline="")
    csv.writer(stream).writerows(rows)
    return stream.getvalue()


@pytest.mark.parametrize("which", [0, 1], ids=["enzyme", "interaction"])
@pytest.mark.parametrize("delta", [-1, 1, 2])
@pytest.mark.parametrize("after_valid", [False, True])
def test_bad_record_width_reports_file_line_and_counts(tmp_path, priors, which, delta, after_valid):
    frame = priors[which]
    row = frame.iloc[0].tolist()
    malformed = row[:delta] if delta < 0 else row + ["EXTRA"] * delta
    rows = [frame.columns.tolist()] + ([row] if after_valid else []) + [malformed]
    path = tmp_path / "bad.csv"
    path.write_text(_csv(rows), encoding="utf-8-sig")
    inputs = list(priors)
    inputs[which] = path
    with pytest.raises(ValueError) as error:
        load_cell_mesh_database(*inputs)
    message = str(error.value)
    assert str(path) in message
    assert f"line {3 if after_valid else 2}" in message
    assert f"expected {len(row)} fields" in message
    assert f"got {len(malformed)}" in message


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("which", [0, 1])
def test_run_rejects_shifted_prior_before_scoring(tmp_path, priors, mode, which):
    data = AnnData(np.array([[4., 1.], [1., 4.]]),
                   obs=pd.DataFrame({"cell_type": ["A", "B"], "sample": ["S", "S"]}, index=["c1", "c2"]),
                   var=pd.DataFrame(index=["G", "R"]))
    frame = priors[which]
    path = tmp_path / "shifted.csv"
    path.write_text(_csv([frame.columns.tolist(), frame.iloc[0].tolist() + ["EXTRA"]]))
    inputs = list(priors)
    inputs[which] = path
    with pytest.raises(ValueError, match="expected .* fields.*got"):
        run_cell_mesh(data, *inputs, min_cells=1, sample_key="sample", sample_mode=mode, n_perms=3)


@pytest.mark.parametrize("suffix", [".csv", ".csv.gz", ".csv.bz2", ".csv.xz", ".csv.zip"])
def test_quoted_multiline_empty_and_literal_fields_survive_validation(tmp_path, priors, suffix):
    enzyme, sensor = [frame.copy() for frame in priors]
    for frame in [enzyme, sensor]:
        frame["metabolite"] = 'Name, with "quotes"\nand another line'
        frame["reference"] = "00123"
        frame["source"] = "NA"
        frame["evidence_level"] = None
    expected = load_cell_mesh_database(enzyme, sensor)
    paths = []
    for i, frame in enumerate([enzyme, sensor]):
        path = tmp_path / f"prior{i}{suffix}"
        frame.to_csv(path, index=False, encoding="utf-8-sig")
        paths.append(path)
    actual = load_cell_mesh_database(*paths)
    for got, want in zip(actual, expected):
        assert got.evidence_level.isna().all() and want.evidence_level.isna().all()
        pd.testing.assert_frame_equal(got.drop(columns="evidence_level"), want.drop(columns="evidence_level"))


def test_csv_physical_line_range_accounts_for_multiline_field(tmp_path, priors):
    path = tmp_path / "multiline.csv"
    path.write_text('metabolite,hmdb_id,sensor_gene,sensor_type\n"name\ncontinued",HMDB1,R,Transporter,EXTRA\n')
    with pytest.raises(ValueError, match=r"lines 2-3.*expected 4 fields.*got 5"):
        load_cell_mesh_database(priors[0], path)


@pytest.mark.parametrize("record", ['"unterminated,HMDB1,R,Transporter\n', '"closed"junk,HMDB1,R,Transporter\n'])
def test_malformed_quoting_reports_source_location(tmp_path, priors, record):
    path = tmp_path / "quote.csv"
    path.write_text("metabolite,hmdb_id,sensor_gene,sensor_type\n" + record)
    with pytest.raises(ValueError, match=r"quote.csv.*line 2.*CSV"):
        load_cell_mesh_database(priors[0], path)


def test_blank_lines_are_distinct_from_short_quoted_records(tmp_path, priors):
    path = tmp_path / "blank.csv"
    sensor = priors[1]
    text = "\n  \n" + _csv([sensor.columns.tolist(), sensor.iloc[0].tolist()]) + "\n\t\n"
    path.write_text(text)
    _, result = load_cell_mesh_database(priors[0], path)
    assert len(result) == 1
    path.write_text(text + '""\n')
    with pytest.raises(ValueError, match="expected 4 fields.*got 1"):
        load_cell_mesh_database(priors[0], path)


def test_long_quoted_metadata_preserves_csv_parser_configuration(tmp_path, priors):
    before = csv.field_size_limit()
    sensor = priors[1].assign(reference="long,quoted\n" * 15000)
    path = tmp_path / "long.csv"
    sensor.to_csv(path, index=False)
    _, result = load_cell_mesh_database(priors[0], path)
    assert result.reference.iloc[0] == sensor.reference.iloc[0]
    assert csv.field_size_limit() == before


@pytest.mark.parametrize("contents", ["", "\n  \n"])
def test_empty_prior_has_explicit_header_error(tmp_path, priors, contents):
    path = tmp_path / "empty.csv"
    path.write_text(contents)
    with pytest.raises(ValueError, match=r"empty.csv.*header"):
        load_cell_mesh_database(priors[0], path)


def test_header_only_prior_and_custom_dataframe_index_remain_valid(tmp_path, priors):
    path = tmp_path / "header.csv"
    priors[1].iloc[:0].to_csv(path, index=False)
    _, result = load_cell_mesh_database(priors[0], path)
    assert result.empty
    enzyme, sensor = [frame.set_axis(["custom-row"]) for frame in priors]
    actual = load_cell_mesh_database(enzyme, sensor)
    expected = load_cell_mesh_database(*priors)
    for got, want in zip(actual, expected):
        pd.testing.assert_frame_equal(got, want)


def test_nul_cannot_be_silently_truncated_by_the_final_parser(tmp_path, priors):
    path = tmp_path / "nul.csv"
    text = priors[1].assign(reference="first~second").to_csv(index=False)
    path.write_text(text.replace("~", "\0"))
    with pytest.raises(ValueError, match=r"nul.csv.*line 2.*NUL"):
        load_cell_mesh_database(priors[0], path)
