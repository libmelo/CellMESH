"""User-visible reporting/persistence contracts, including scientific missingness."""
import copy

import h5py
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import CellMeshResult, run_cell_mesh
from cellmesh.result_io import _read, _write


def make_result(mode="sample_aware", *, all_zero=False, production=True, n_perms=7, min_cells=2):
    x = np.array([[4, 9, 0, 1], [2, 4, 0, 0], [0, 1, 0, 8],
                  [3, 1, 0, 1], [0, 2, 0, 6], [5, 2, 0, 1]], dtype=float)
    if all_zero:
        x[:, :2] = 0
    adata = AnnData(x, obs=pd.DataFrame({"ct": ["A", "A", "B", "A", "B", "A"],
                "sample": ["S1", "S1", "S1", "S2", "S2", "S3"]}, index=list("abcdef")),
                var=pd.DataFrame(index=["G1", "G2", "ZERO", "R"]))
    enzyme = pd.DataFrame([
        {"hmdb_id": "HMDB0000001", "metabolite": "M<script>bad()</script>", "gene": g,
         "reaction": rxn, "role": role, "source": "curated", "reference": "doi:example"}
        for rxn, role, genes in [("complete", "production", ["G1", "G2"]),
                                  ("duplicate", "production", ["G2", "G1"]),
                                  ("subset", "production", ["G1"]),
                                  ("partial", "production", ["G1", "MISSING"]),
                                  ("zero_cons", "degradation", ["ZERO"])]
        for g in genes
    ] + [{"hmdb_id": "HMDB0000002", "metabolite": "No measured producer", "gene": "MISSING",
          "reaction": "unknown", "role": "production"}])
    sensor = pd.DataFrame([
        {"hmdb_id": hmdb, "metabolite": "alias", "sensor_gene": g, "sensor_type": "Receptor"}
        for hmdb, g in [("HMDB0000001", "R"), ("HMDB0000001", "ABSENT_R"),
                        ("HMDB0000002", "R"), ("HMDB0000003", "R")]])
    if not production:
        enzyme["role"] = "degradation"
    return run_cell_mesh(adata, enzyme, sensor, cell_type_key="ct", sample_key="sample", sample_mode=mode,
                         min_cells=min_cells, min_expr_frac=0.75, n_perms=n_perms, random_state=18,
                         store_null_scores=True)


SELECTOR = dict(hmdb_id="HMDB0000001", sensor_gene="R", sender="A", receiver="B")


def assert_nested_equal(left, right):
    if isinstance(left, pd.DataFrame):
        pd.testing.assert_frame_equal(left, right, check_exact=True)
        assert_nested_equal(left.attrs, right.attrs)
    elif isinstance(left, pd.Series):
        pd.testing.assert_series_equal(left, right, check_exact=True)
    elif isinstance(left, np.ndarray):
        np.testing.assert_array_equal(left, right)
    elif isinstance(left, dict):
        assert left.keys() == right.keys()
        for key in left:
            assert_nested_equal(left[key], right[key])
    elif isinstance(left, (list, tuple)):
        assert len(left) == len(right)
        for a, b in zip(left, right):
            assert_nested_equal(a, b)
    elif left is pd.NA or left is None:
        assert left is right
    elif isinstance(left, float) and np.isnan(left):
        assert np.isnan(right)
    else:
        assert left == right


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("zero", [False, True])
def test_save_load_retains_scientific_results_and_explanations(tmp_path, mode, zero):
    result = make_result(mode, all_zero=zero)
    before = copy.deepcopy(result)
    detail = result.explain_event(**SELECTOR)
    files = result.save(tmp_path / "analysis", report=True, selected_events=result.select_events(metabolite_sensor_pairs=[('HMDB0000001','R')], source='all', cell_pairs=[('A','B')]), include_null_scores=True)
    assert set(p.name for p in tmp_path.iterdir()) == {
        "analysis.cellmesh.h5", "analysis.positive_events.csv", "analysis.report.html", "analysis.reports", "analysis.metabolite_sensor_events.csv"}
    restored = CellMeshResult.load(files["result"])
    for name in ["events", "sender_scores", "receiver_scores", "sample_events", "sample_validation",
                 "sample_sender_scores", "sample_receiver_scores", "celltype_qc", "availability_results",
                 "parameters", "priors"]:
        assert_nested_equal(getattr(result, name), getattr(restored, name))
        assert_nested_equal(getattr(result, name), getattr(before, name))
    assert_nested_equal(detail, restored.explain_event(**SELECTOR))
    with h5py.File(files["result"]) as handle:
        assert {"metadata", "results", "sample_results", "explanation", "coverage", "diagnostics", "priors"} <= set(handle)
        assert handle.attrs["schema_version"] == 1
    if mode == "sample_aware":
        samples = detail["sample_events"].set_index("sample")
        assert pd.isna(samples.at["S3", "cell_mesh_score"])
        assert "receiver_unobserved" in samples.at["S3", "reasons"]
        if zero:
            assert samples.loc[["S1", "S2"], "cell_mesh_score"].eq(0).all()
        assert detail["aggregate_event"].iloc[0].cell_mesh_score == samples.cell_mesh_score.median()
    # Old CSV output remains an independent, usable interface after restoring.
    before.to_csv(tmp_path / "old")
    restored.to_csv(tmp_path / "new")
    for name in ["events", "sender_scores", "receiver_scores", "celltype_qc"]:
        assert (tmp_path / f"old.{name}.csv").read_bytes() == (tmp_path / f"new.{name}.csv").read_bytes()


def test_coverage_distinguishes_missing_zero_partial_and_deduplication():
    result = make_result()
    cov = result.coverage_report()["details"]
    unit = cov.loc[(cov["sample"] == "S1") & (cov.cell_type == "A") & (cov.hmdb_id == "HMDB0000001")].set_index("direction")
    assert unit.at["P", "n_prior_genes"] == 3
    assert unit.at["P", "n_input_genes"] == 2
    assert unit.at["P", "input_gene_fraction"] == pytest.approx(2 / 3)
    assert unit.at["P", "n_annotated_reactions"] == 4
    assert unit.at["P", "n_reactions"] == 2
    assert unit.at["P", "n_fully_covered"] == 1
    assert unit.at["P", "n_partially_covered"] == 1
    assert unit.at["C", "expression_state"] == "measured_zero"
    assert unit.at["E", "expression_state"] == "prior_missing"
    assert pd.isna(unit.at["E", "input_gene_fraction"])
    detail = result.explain_event(**SELECTOR, sample="S1")
    reactions = detail["reaction_details"].set_index("reaction")
    assert reactions.at["duplicate", "selection_status"] == "duplicate_gene_set"
    assert reactions.at["subset", "selection_status"] == "strict_subset"
    assert pd.isna(reactions.at["duplicate", "activity_unweighted"])
    genes = detail["gene_details"].query("cell_type == 'A'").set_index("gene")
    assert not genes.at["MISSING", "input_present"]
    assert pd.isna(genes.at["MISSING", "mean_expression"])
    assert genes.at["ZERO", "input_present"]
    assert genes.at["ZERO", "mean_expression"] == 0


def test_explanation_recovers_actual_scoring_steps():
    result = make_result()
    detail = result.explain_event(**SELECTOR)
    sender = detail["sender_components"].set_index("sample")
    receiver = detail["receiver_components"].set_index("sample")
    for sample in ["S1", "S2"]:
        s, r = sender.loc[sample], receiver.loc[sample]
        assert s.P == pytest.approx(s.P_unweighted * s.abundance_weight)
        assert s.P_score == pytest.approx(s.P / (s.P + s.P_ref))
        assert s.base_availability == pytest.approx(s.P_score ** 2 / (s.P_score + s.C_score))
        assert s.E_effective == 0.5  # Missing exporter prior is explicitly neutral.
        assert s.sender_score == pytest.approx(s.base_availability * s.E_factor)
        assert r.sensor_score == pytest.approx(r.mean_expression / (r.mean_expression + r.reference))
        observed = detail["sample_events"].set_index("sample").at[sample, "cell_mesh_score"]
        assert observed == pytest.approx(np.sqrt(s.sender_score * r.sensor_score))
    gate_detail = result.explain_event(**dict(SELECTOR, sender="B", receiver="A"), sample="S1")
    assert gate_detail["receiver_components"].iloc[0].gate_applied
    assert "receiver_expression_fraction_gate" in gate_detail["sample_events"].iloc[0].reasons


@pytest.mark.parametrize("selector,reason", [
    (dict(SELECTOR, sensor_gene="ABSENT_R"), "sensor_gene_unavailable"),
    (dict(SELECTOR, hmdb_id="HMDB0000002"), "production_genes_unavailable"),
    (dict(SELECTOR, hmdb_id="HMDB0000003"), "production_prior_missing"),
    (dict(SELECTOR, sender="UNKNOWN"), "sender_unobserved"),
    (dict(SELECTOR, sensor_gene="UNKNOWN"), "interaction_prior_missing"),
])
def test_unscored_candidates_can_be_explained(selector, reason):
    detail = make_result().explain_event(**selector, sample="S1")
    assert detail["sample_events"].iloc[0].evaluation_status == "not_evaluable"
    assert reason in detail["sample_events"].iloc[0].reasons


def test_null_storage_is_explicit_and_does_not_mutate_result(tmp_path):
    result = make_result()
    null = result.events.attrs["sample_aware_null_scores"].copy()
    output = result.save(tmp_path / "small.cellmesh.h5")
    restored = CellMeshResult.load(output["result"])
    assert "sample_aware_null_scores" not in restored.events.attrs
    assert restored.events.attrs["null_scores_stored"] is False
    assert restored.parameters["store_null_scores"] is True  # Original computation option.
    pd.testing.assert_frame_equal(result.events.attrs["sample_aware_null_scores"], null)
    with h5py.File(output["result"]) as handle:
        assert "permutation_null" not in handle["optional"]


def test_codec_handles_axes_nullable_values_nested_cells_and_attributes(tmp_path):
    frame = pd.DataFrame({
        "nullable": pd.Series([1, pd.NA, 0], dtype="Int64"),
        "boolean": pd.Series([True, pd.NA, False], dtype="boolean"),
        "string": pd.Series(["中文", pd.NA, ""], dtype="string"),
        "category": pd.Categorical(["B", None, "A"], categories=["A", "B", "unused"], ordered=True),
        "object": [["G1", "G2"], {"source": "x", "n": None}, ("zero", 0)],
        "time": pd.to_datetime(["2026-01-01", None, "2026-01-03"], utc=True),
    })
    frame.index = pd.MultiIndex(levels=[["S1", "S2", "unused"], ["A", "B"]],
                                codes=[[0, 0, 1], [0, 1, 0]], names=["sample", "type"])
    frame.columns.name = "feature"
    frame.attrs = {"count": np.int64(3), "nested": {"missing": pd.NA, "pair": (1, 2)},
                   "array": np.array([1., np.nan])}
    path = tmp_path / "codec.h5"
    with h5py.File(path, "w") as handle:
        _write(handle, "frame", frame, include_null_scores=True)
    with h5py.File(path) as handle:
        restored = _read(handle["frame"])
    assert_nested_equal(frame, restored)


def test_html_is_self_contained_escaped_and_reports_missing_values(tmp_path):
    result = make_result()
    path = result.to_html(tmp_path / "report.html", selected_events=result.select_events(metabolite_sensor_pairs=[('HMDB0000001','R')], source='all', cell_pairs=[('A','B')]))
    text = next(value for key,value in path.items() if key.startswith("pair_report:")).read_text()
    assert "<svg" in text
    assert "NA" in text
    assert "&lt;script&gt;bad()&lt;/script&gt;" in text
    assert "<script>" not in text
    assert '<script src=' not in text and '<img src=' not in text
    assert "n_prior_genes" in text and "P_ref" in text


def test_save_failure_keeps_previous_deliverables(tmp_path, monkeypatch):
    import cellmesh.result_io as io
    result = make_result()
    paths = result.save(tmp_path / "analysis", report=True)
    before = {k: p.read_bytes() for k, p in paths.items()}
    original = io.os.replace
    failed = False

    def fail(source, target):
        nonlocal failed
        if str(target).endswith(".report.html") and not failed:
            failed = True
            raise OSError("simulated publication failure")
        return original(source, target)

    monkeypatch.setattr(io.os, "replace", fail)
    with pytest.raises(OSError, match="publication failure"):
        result.save(tmp_path / "analysis", report=True)
    assert {k: p.read_bytes() for k, p in paths.items()} == before
    assert {p for p in tmp_path.rglob('*') if p.is_file()} == set(paths.values())


def test_rejects_bad_format_and_preserves_saved_file_on_invalid_metadata(tmp_path):
    result = make_result()
    path = result.save(tmp_path / "result")["result"]
    old = path.read_bytes()
    result.parameters["unsupported"] = object()
    with pytest.raises(TypeError, match="Unsupported"):
        result.save(path)
    assert path.read_bytes() == old
    with h5py.File(path, "r+") as handle:
        handle.attrs["schema_version"] = 999
    with pytest.raises(ValueError, match="schema version"):
        CellMeshResult.load(path)


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
def test_empty_event_result_is_serializable_and_explains_missing_production(tmp_path, mode):
    result = make_result(mode, n_perms=0, production=False)
    assert result.events.empty and result.sender_scores.empty
    path = result.save(tmp_path / "empty")["result"]
    restored = CellMeshResult.load(path)
    assert restored.events.empty
    detail = restored.explain_event(**dict(SELECTOR, hmdb_id="HMDB0000002"))
    assert all("production_prior_missing" in reasons for reasons in detail["sample_events"].reasons)


def test_save_removes_stale_report_without_touching_csv_or_other_files(tmp_path):
    result = make_result(n_perms=0)
    result.save(tmp_path / "analysis", report=True)
    csv = tmp_path / "analysis.events.csv"
    csv.write_text("existing export")
    result.save(tmp_path / "analysis", report=False)
    assert not (tmp_path / "analysis.report.html").exists()
    assert csv.read_text() == "existing export"


def test_symlink_destination_is_rejected(tmp_path):
    result = make_result(n_perms=0)
    other = tmp_path / "other"
    other.write_text("do not overwrite")
    (tmp_path / "analysis.cellmesh.h5").symlink_to(other)
    with pytest.raises(ValueError, match="symlinks"):
        result.save(tmp_path / "analysis")
    assert other.read_text() == "do not overwrite"


def test_result_without_explanation_can_still_be_saved(tmp_path):
    result = CellMeshResult(events=pd.DataFrame(), sender_scores=pd.DataFrame(),
                            receiver_scores=pd.DataFrame(), parameters={"source": "manual"})
    path = result.save(tmp_path / "manual", report=False)["result"]
    restored = CellMeshResult.load(path)
    assert restored.priors is None and restored.availability_results is None
    with pytest.raises(ValueError, match="requires"):
        restored.to_html(tmp_path / "report.html")


def test_optional_figure_uses_captured_values_and_preserves_na(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cellmesh import plot_event_explanation
    result = make_result(n_perms=0)
    fig, axes = plot_event_explanation(result, **SELECTOR)
    try:
        last = axes.flat[3]
        detail = result.explain_event(**SELECTOR)["sample_events"]
        np.testing.assert_allclose([bar.get_width() for bar in last.patches],
                                   detail.cell_mesh_score.dropna().to_numpy())
        assert any("NA" in text.get_text() for text in last.texts)
        fig.savefig(tmp_path / "event.png", dpi=80)
        assert (tmp_path / "event.png").stat().st_size > 1000
    finally:
        plt.close(fig)


@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("production", [False, True])
def test_all_events_csv_matches_full_result_and_legacy_export(tmp_path, mode, production):
    result = make_result(mode, n_perms=0, production=production, all_zero=True, min_cells=99)
    if not result.events.empty:
        result.events.loc[0, "metabolite"] = '代谢物,含"引号"\n换行'
        assert result.events.cell_mesh_score.eq(0).any()
        assert not result.events.passes_min_cells.all()
    files = result.save(tmp_path / "analysis.cellmesh.h5", report=True, selected_events=result.select_events(metabolite_sensor_pairs=[('HMDB0000001','R')], source='all', cell_pairs=[('A','B')]), include_all_events=True)
    csv = pd.read_csv(files["events_csv"], float_precision="round_trip")
    assert files["events_csv"].name == "analysis.all_events.csv"
    assert csv.columns.tolist() == result.events.columns.tolist()
    pd.testing.assert_frame_equal(csv, result.events, check_dtype=False, check_exact=True)
    # Missing statistical results stay missing, not zero; zero/QC-failing
    # events stay in the CSV even though only one HTML explanation was selected.
    assert csv.perm_pvalue.isna().all()
    assert csv.fdr_global.isna().all() and csv.fdr_sensor_type.isna().all()
    result.to_csv(tmp_path / "legacy")
    assert files["events_csv"].read_bytes() == (tmp_path / "legacy.events.csv").read_bytes()


def test_all_events_csv_keeps_nonsignificant_events(tmp_path):
    result = make_result(n_perms=7)
    assert result.events.fdr_sensor_type.gt(0.05).any()
    files = result.save(tmp_path / "analysis", report=True, selected_events=result.select_events(metabolite_sensor_pairs=[('HMDB0000001','R')], source='all', cell_pairs=[('A','B')]), include_all_events=True)
    csv = pd.read_csv(files["events_csv"], float_precision="round_trip")
    pd.testing.assert_frame_equal(csv, result.events, check_dtype=False, check_exact=True)


@pytest.mark.parametrize("phase", ["staging", "publication"])
def test_csv_failure_rolls_back_entire_delivery(tmp_path, monkeypatch, phase):
    from pathlib import Path
    import cellmesh.result_io as io
    result = make_result(n_perms=0)
    files = result.save(tmp_path / "analysis", report=True, include_all_events=True)
    before = {name: path.read_bytes() for name, path in files.items()}
    failed = False
    if phase == "staging":
        original = pd.DataFrame.to_csv
        def fail_csv(self, path, *args, **kwargs):
            if str(path).endswith(".all_events.csv"):
                Path(path).write_text("incomplete CSV")
                raise OSError("simulated CSV failure")
            return original(self, path, *args, **kwargs)
        monkeypatch.setattr(pd.DataFrame, "to_csv", fail_csv)
    else:
        original = io.os.replace
        def fail_replace(source, target):
            nonlocal failed
            if target == files["events_csv"] and not failed:
                failed = True
                raise OSError("simulated CSV failure")
            return original(source, target)
        monkeypatch.setattr(io.os, "replace", fail_replace)
    with pytest.raises(OSError, match="CSV failure"):
        result.save(tmp_path / "analysis", report=True, include_all_events=True)
    assert {name: path.read_bytes() for name, path in files.items()} == before
    assert {p for p in tmp_path.rglob('*') if p.is_file()} == set(files.values())


def test_all_events_csv_symlink_is_rejected(tmp_path):
    other = tmp_path / "external.csv"
    other.write_text("do not touch")
    (tmp_path / "analysis.all_events.csv").symlink_to(other)
    with pytest.raises(ValueError, match="symlinks"):
        make_result(n_perms=0).save(tmp_path / "analysis", include_all_events=True)
    assert other.read_text() == "do not touch"
