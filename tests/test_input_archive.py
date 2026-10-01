"""Optional input bundles must preserve scientific results and original inputs."""
import shutil

import h5py
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from scipy import sparse

from cellmesh import CellMeshResult, run_cell_mesh


@pytest.fixture
def inputs():
    matrix = np.array([[4, 9, 1, 0], [0, 2, 4, 8], [3, 7, 2, 1],
                       [1, 0, 5, 3], [0, 0, 2, 5], [2, 6, 0, 0]], dtype=np.float32)
    adata = AnnData(matrix, obs=pd.DataFrame({
        "ct": pd.Categorical(["A", "B", "A", "B", "B", "A"], categories=["A", "B", "unused"]),
        "sample": ["S1", "S1", "S2", "S2", "S2", "S3"],
        "annotation": ["甲", "乙", "甲", "乙", "乙", "甲"],
    }, index=list("abcdef")), var=pd.DataFrame({"kind": ["enzyme", "enzyme", "sensor", "unused"]},
                                             index=["G1", "G2", "R", "UNUSED"]))
    enzyme = pd.DataFrame(dict(metabolite=["M", "M"], hmdb_id=["HMDB1"] * 2,
                              gene=["G1", "G2"], role=["production"] * 2, reaction=["P"] * 2))
    sensor = pd.DataFrame(dict(metabolite=["M"], hmdb_id=["HMDB1"], sensor_gene=["R"],
                              sensor_type=["Receptor"]))
    return adata, enzyme, sensor


def run(inputs, **kwargs):
    options = dict(cell_type_key="ct", sample_key="sample", min_cells=1,
                   n_perms=5, random_state=19, store_null_scores=True)
    options.update(kwargs)
    return run_cell_mesh(*inputs, **options)


def assert_results(left, right):
    for name in ["events", "sender_scores", "receiver_scores", "celltype_qc",
                 "sample_events", "sample_validation", "sample_sender_scores", "sample_receiver_scores"]:
        a, b = getattr(left, name), getattr(right, name)
        if a is None:
            assert b is None
        else:
            pd.testing.assert_frame_equal(a, b, check_exact=True)
    if "sample_aware_null_scores" in left.events.attrs:
        pd.testing.assert_frame_equal(left.events.attrs["sample_aware_null_scores"],
                                      right.events.attrs["sample_aware_null_scores"], check_exact=True)


def assert_matrix(left, right):
    assert sparse.issparse(left) == sparse.issparse(right)
    assert left.dtype == right.dtype
    if sparse.issparse(left):
        assert left.format == right.format
        np.testing.assert_array_equal(left.data, right.data)
        np.testing.assert_array_equal(left.indices, right.indices)
        np.testing.assert_array_equal(left.indptr, right.indptr)
    else:
        np.testing.assert_array_equal(left, right)


@pytest.mark.parametrize("mode", ["sample_aware", "pooled_stratified"])
@pytest.mark.parametrize("kind", ["dense", "csr", "csc"])
@pytest.mark.parametrize("layer", [None, "chosen"])
def test_selected_input_roundtrip_and_exact_replay(inputs, tmp_path, mode, kind, layer):
    adata, enzyme, sensor = inputs
    if kind != "dense":
        adata.X = getattr(sparse, kind + "_matrix")(adata.X)
    expected_matrix = adata.X.copy()
    expected_obs, expected_var = adata.obs.copy(), adata.var.copy()
    if layer:
        adata.layers[layer] = adata.X.copy()
        adata.X = adata.X * 0
    adata.layers["unselected"] = adata.X.copy()
    adata.raw = adata.copy()
    adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))
    adata.uns["unrelated"] = "do not archive"
    baseline = run(inputs, sample_mode=mode, layer=layer)
    result = run(inputs, sample_mode=mode, layer=layer, retain_input=True)
    assert_results(baseline, result)
    assert baseline.parameters == result.parameters
    # Subsequent caller mutations cannot change the captured matrix or grouping.
    if layer:
        adata.layers[layer][:] = 0
    else:
        adata.X[:] = 0
    adata.obs["ct"] = "changed"
    adata.var_names = ["different" + str(i) for i in range(adata.n_vars)]
    paths = result.save(tmp_path / "analysis", include_input=True, report=True, include_null_scores=True)
    assert set(paths) == {"result", "positive_csv", "report", "input", "navigation"}
    restored = CellMeshResult.load(paths["result"])
    assert restored._input_adata is None
    archived = restored.load_input()
    assert_matrix(archived.X, expected_matrix)
    pd.testing.assert_frame_equal(archived.obs, expected_obs)
    pd.testing.assert_frame_equal(archived.var, expected_var)
    assert not archived.layers and archived.raw is None and not archived.obsm
    assert "unrelated" not in archived.uns
    assert restored.input_provenance["source_layer"] == layer
    assert restored.input_provenance["raw_counts_status"] == "not inferred"
    replay = run((archived, restored.priors["enzyme"], restored.priors["interaction"]), sample_mode=mode)
    assert_results(result, replay)
    assert_results(result, restored)
    assert "analysis.input.h5ad" in paths["report"].read_text()
    # Public in-memory access also returns an independent copy.
    copy = result.load_input()
    copy.X = copy.X * 0
    assert_matrix(result.load_input().X, expected_matrix)
    backed = result.load_input(backed="r")
    try:
        assert backed.isbacked and backed.file.is_open
        assert backed.shape == expected_matrix.shape
    finally:
        backed.file.close()


@pytest.mark.parametrize("kind", ["dense", "csr", "csc"])
@pytest.mark.parametrize("layer", [None, "chosen"])
def test_capture_backed_view_keeps_order_and_original_file(inputs, tmp_path, kind, layer):
    from cellmesh.input_archive import _sha256
    adata, enzyme, sensor = inputs
    if kind != "dense":
        adata.X = getattr(sparse, kind + "_matrix")(adata.X)
    if layer:
        adata.layers[layer] = adata.X.copy()
        adata.X = adata.X * 0
    path = tmp_path / "original.h5ad"
    adata.write_h5ad(path)
    digest = _sha256(path)
    backed = read_h5ad(path, backed="r")
    rows, cols = [4, 0, 2, 3], [3, 2, 0, 1]
    try:
        view = backed[rows, cols]
        expected = run((adata[rows, cols].copy(), enzyme, sensor), layer=layer)
        result = run((view, enzyme, sensor), layer=layer, retain_input=True)
        assert_results(expected, result)
    finally:
        backed.file.close()
    assert _sha256(path) == digest
    paths = result.save(tmp_path / "captured", include_input=True)
    restored = CellMeshResult.load(paths["result"]).load_input()
    assert restored.obs_names.tolist() == adata.obs_names[rows].tolist()
    assert restored.var_names.tolist() == adata.var_names[cols].tolist()
    source = adata.layers[layer] if layer else adata.X
    selected = source[rows, :][:, cols]
    assert sparse.issparse(restored.X) == sparse.issparse(selected)
    assert restored.X.dtype == selected.dtype
    # Backed noncontiguous CSC slicing uses bounded CSR row blocks internally.
    # Logical values/order and sparse storage must survive that layout change.
    if sparse.issparse(selected):
        np.testing.assert_array_equal(restored.X.toarray(), selected.toarray())
    else:
        np.testing.assert_array_equal(restored.X, selected)


def test_default_avoids_capture_and_archive_loading(inputs, tmp_path, monkeypatch):
    import cellmesh.input_archive as module
    def forbidden(*args, **kwargs):
        raise AssertionError("optional matrix should not be accessed")
    monkeypatch.setattr(module, "capture_input", forbidden)
    result = run(inputs)
    assert result._input_adata is None and result.input_provenance is None
    with pytest.raises(ValueError, match="retain_input"):
        result.save(tmp_path / "missing", include_input=True)
    paths = result.save(tmp_path / "ordinary")
    monkeypatch.setattr(module, "read_h5ad", forbidden)
    restored = CellMeshResult.load(paths["result"])
    assert restored.input_archive is None
    with pytest.raises(ValueError, match="No input archive"):
        restored.load_input()
    assert {p.name for p in tmp_path.iterdir()} == {"ordinary.cellmesh.h5",
                                                   "ordinary.positive_events.csv", "ordinary.report.html", "ordinary.reports"}


def test_retained_input_is_still_optional_at_save(inputs, tmp_path):
    result = run(inputs, retain_input=True)
    paths = result.save(tmp_path / "small", report=True)
    restored = CellMeshResult.load(paths["result"])
    assert restored.input_provenance is not None and restored.input_archive is None
    assert "未附带细胞级输入文件" in paths["report"].read_text()
    assert not (tmp_path / "small.input.h5ad").exists()
    with pytest.raises(ValueError, match="No input archive"):
        restored.load_input()


def test_old_schema_v1_files_remain_loadable(inputs, tmp_path):
    path = run(inputs).save(tmp_path / "legacy")["result"]
    with h5py.File(path, "r+") as handle:
        del handle["metadata/expression_input"]
        del handle["metadata/input_archive"]
    restored = CellMeshResult.load(path)
    assert restored.input_provenance is None and restored.input_archive is None
    restored.save(tmp_path / "new")


def test_relocation_resave_and_checksum_verification(inputs, tmp_path, monkeypatch):
    result = run(inputs, retain_input=True)
    paths = result.save(tmp_path / "analysis", include_input=True)
    moved = tmp_path / "moved"
    moved.mkdir()
    for p in paths.values():
        shutil.move(p, moved / p.name)
    # Relative associations follow the result, independently of working directory.
    monkeypatch.chdir(tmp_path.parent)
    restored = CellMeshResult.load(moved / paths["result"].name)
    assert restored.load_input().shape == inputs[0].shape
    copied = restored.save(moved / "renamed", include_input=True, report=True)
    assert copied["input"].read_bytes() == (moved / paths["input"].name).read_bytes()
    assert restored.input_archive["filename"] == "renamed.input.h5ad"
    assert "renamed.input.h5ad" in restored.to_html(moved / "extra.html").read_text()
    # Alter a value while keeping a valid AnnData file: only explicit input load
    # fails; result inspection remains usable.
    with h5py.File(copied["input"], "r+") as handle:
        handle["X"][0, 0] += 1
    damaged = CellMeshResult.load(copied["result"])
    with pytest.raises(ValueError, match="checksum mismatch"):
        damaged.load_input(backed="r")
    with pytest.raises(ValueError, match="checksum mismatch"):
        damaged.save(moved / "bad_copy", include_input=True)
    assert not (moved / "bad_copy.cellmesh.h5").exists()
    copied["input"].unlink()
    assert not CellMeshResult.load(copied["result"]).events.empty
    with pytest.raises(FileNotFoundError, match="missing"):
        damaged.load_input()


@pytest.mark.parametrize("failure", ["input", "report", "cleanup"])
def test_publication_failure_restores_all_files_and_associations(inputs, tmp_path, monkeypatch, failure):
    import cellmesh.result_io as module
    result = run(inputs, retain_input=True)
    paths = result.save(tmp_path / "analysis", include_input=True, report=True)
    before = {k: p.read_bytes() for k, p in paths.items()}
    old_archive = result.input_archive.copy()
    replace, unlink = module.os.replace, module.Path.unlink
    failed = False
    def fail_replace(source, target):
        nonlocal failed
        if failure in ("input", "report") and target == paths[failure] and not failed:
            failed = True
            raise OSError("simulated publication failure")
        return replace(source, target)
    def fail_unlink(target, *args, **kwargs):
        nonlocal failed
        if failure == "cleanup" and target == paths["input"] and not failed:
            failed = True
            raise OSError("simulated cleanup failure")
        return unlink(target, *args, **kwargs)
    monkeypatch.setattr(module.os, "replace", fail_replace)
    monkeypatch.setattr(module.Path, "unlink", fail_unlink)
    with pytest.raises(OSError, match="simulated"):
        result.save(tmp_path / "analysis", include_input=failure != "cleanup", report=failure != "cleanup")
    assert {k: p.read_bytes() for k, p in paths.items()} == before
    assert result.input_archive == old_archive


def test_staging_failure_and_stale_companion_cleanup(inputs, tmp_path, monkeypatch):
    result = run(inputs, retain_input=True)
    paths = result.save(tmp_path / "analysis", include_input=True, report=True)
    before = {k: p.read_bytes() for k, p in paths.items()}
    def fail(*args, **kwargs):
        raise OSError("simulated H5AD failure")
    with monkeypatch.context() as patch:
        patch.setattr(AnnData, "write_h5ad", fail)
        with pytest.raises(OSError, match="H5AD failure"):
            result.save(tmp_path / "analysis", include_input=True, report=True)
    assert {k: p.read_bytes() for k, p in paths.items()} == before
    unrelated = tmp_path / "analysis.events.csv"
    unrelated.write_text("keep existing CSV")
    result.save(tmp_path / "analysis", report=False)
    assert not paths["input"].exists() and not paths["report"].exists()
    assert unrelated.read_text() == "keep existing CSV"
    assert result.input_archive is None


def test_companion_symlink_and_invalid_metadata_paths(inputs, tmp_path):
    result = run(inputs, retain_input=True)
    other = tmp_path / "unrelated.h5ad"
    other.write_text("do not touch")
    (tmp_path / "analysis.input.h5ad").symlink_to(other)
    with pytest.raises(ValueError, match="symlinks"):
        result.save(tmp_path / "analysis", include_input=True)
    assert other.read_text() == "do not touch"
    paths = result.save(tmp_path / "valid", include_input=True)
    restored = CellMeshResult.load(paths["result"])
    for filename in ["../unrelated.input.h5ad", "/other.input.h5ad", "..\\other.input.h5ad"]:
        restored.input_archive["filename"] = filename
        with pytest.raises(ValueError, match="Invalid input archive filename"):
            restored.load_input()


@pytest.mark.parametrize("bad", ["yes", 1, None])
def test_new_boolean_options_are_strict(inputs, tmp_path, bad):
    with pytest.raises((TypeError, ValueError), match="retain_input"):
        run(inputs, retain_input=bad)
    result = run(inputs)
    with pytest.raises(TypeError, match="include_input"):
        result.save(tmp_path / "result", include_input=bad)


def test_input_backed_mode_is_read_only(inputs, tmp_path):
    result = run(inputs, retain_input=True)
    for mode in ["r+", "w", True]:
        with pytest.raises(ValueError, match="read-only"):
            result.load_input(backed=mode)


@pytest.mark.parametrize("backed", [None, "r"])
def test_archived_selected_layer_works_with_existing_violin_plots(inputs, tmp_path, backed):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cellmesh import plot_metabolite_secretion_violin, plot_receptor_expression_violin
    adata, _, _ = inputs
    adata.layers["chosen"] = adata.X.copy()
    adata.X[:] = 0
    result = run(inputs, layer="chosen", retain_input=True)
    paths = result.save(tmp_path / "analysis", include_input=True)
    restored = CellMeshResult.load(paths["result"])
    archived = restored.load_input(backed=backed)
    try:
        for function, options in [(plot_metabolite_secretion_violin, dict(hmdb_id="HMDB1")),
                                  (plot_receptor_expression_violin, dict(hmdb_id="HMDB1", receptor_gene="R"))]:
            original = function(result, adata, **options)
            loaded = function(restored, archived, **options)
            try:
                pd.testing.assert_frame_equal(original["plot_data"], loaded["plot_data"], check_exact=True)
            finally:
                plt.close(original["fig"])
                plt.close(loaded["fig"])
            # An explicit request must never silently fall back to X.
            with pytest.raises(KeyError, match="not found in adata.layers"):
                function(restored, archived, layer="chosen", **options)
        archived.uns.clear()
        with pytest.raises(KeyError, match="not found in adata.layers"):
            plot_receptor_expression_violin(restored, archived, hmdb_id="HMDB1", receptor_gene="R")
    finally:
        if archived.isbacked:
            archived.file.close()
