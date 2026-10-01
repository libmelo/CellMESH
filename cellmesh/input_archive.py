"""Optional standard AnnData snapshots, separate from lightweight results."""
from __future__ import annotations

import copy
import hashlib
from pathlib import Path
import shutil

from anndata import AnnData, read_h5ad
import pandas as pd
from scipy import sparse

from .preprocess import _expression_source, _slice_expression


def _copy_annotations(frame):
    # Strip AnnData's DataFrameView before copying: its deepcopy can follow
    # parent references and return annotations for the entire unsliced object.
    result = pd.DataFrame(frame).copy(deep=True)
    for column in result.select_dtypes(include="object"):
        result[column] = pd.Series([copy.deepcopy(x) for x in result[column]],
                                   index=result.index, dtype=object)
    result.attrs = copy.deepcopy(frame.attrs)
    return result


def capture_input(adata, *, layer, cell_type_key, sample_key):
    """Own the selected matrix and annotations before scoring, without densifying."""
    from . import __version__

    source = _expression_source(adata, layer)
    matrix = _slice_expression(source).copy()
    # AnnData/H5AD support compressed sparse matrices. Preserve compressed
    # coordinates (including duplicates), values and dtype, without summation.
    if sparse.issparse(matrix):
        if matrix.format == "csc":
            matrix = sparse.csc_matrix(matrix, copy=False)
        else:
            matrix = sparse.csr_matrix(matrix, copy=False)
    provenance = {
        "package_version": __version__,
        "source": "adata.X" if layer is None else f"adata.layers[{layer!r}]",
        "source_layer": layer,
        "stored_matrix": "X",
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "dtype": str(matrix.dtype),
        "matrix_format": matrix.format if sparse.issparse(matrix) else "dense",
        "cell_type_key": cell_type_key,
        "sample_key": sample_key,
        "processing_state": "as supplied; CELL MESH applies no expression transformation",
        "raw_counts_status": "not inferred",
        "scope": "all input cells and genes; selected matrix and obs/var only",
    }
    snapshot = AnnData(X=matrix, obs=_copy_annotations(adata.obs), var=_copy_annotations(adata.var))
    # Avoid optional None values for compatibility with older AnnData writers.
    snapshot.uns["cellmesh_input"] = {k: v for k, v in provenance.items() if v is not None}
    return snapshot, provenance


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _archive_path(result):
    if not result.input_archive or result._result_path is None:
        raise ValueError("No input archive is associated with this result; run with retain_input=True")
    name = result.input_archive.get("filename")
    # Only sibling filenames emitted by save() are valid. Do not interpret
    # stored arbitrary paths relative to the current working directory.
    if (not isinstance(name, str) or not name.endswith(".input.h5ad")
            or Path(name).name != name or "/" in name or "\\" in name):
        raise ValueError("Invalid input archive filename in result metadata")
    return result._result_path.parent / name


def _verify_archive(path, info):
    if not path.is_file():
        raise FileNotFoundError(f"Associated input archive is missing: {path}")
    if path.stat().st_size != info.get("size_bytes") or _sha256(path) != info.get("sha256"):
        raise ValueError(f"Input archive checksum mismatch: {path}")


def stage_input(result, destination):
    """Write/copy one companion into staging; never publish incomplete output."""
    if result._input_adata is not None:
        # No mutation through automatic string-to-category conversion.
        result._input_adata.write_h5ad(
            destination, compression="gzip", convert_strings_to_categoricals=False,
        )
    else:
        source = _archive_path(result)
        _verify_archive(source, result.input_archive)
        shutil.copy2(source, destination)
        # Check the actual staged copy, including possible changes during copy.
        _verify_archive(destination, result.input_archive)
    return {"filename": destination.name, "sha256": _sha256(destination),
            "size_bytes": destination.stat().st_size, "format": "AnnData H5AD",
            "stored_matrix": "X"}


def load_input(result, *, backed=None):
    """Explicit loading only; verify the archive without materializing its matrix."""
    if backed not in (None, "r"):
        raise ValueError("backed must be None or 'r'; input archives are read-only")
    if backed is None and result._input_adata is not None:
        snapshot = result._input_adata
        return AnnData(X=snapshot.X.copy(), obs=_copy_annotations(snapshot.obs),
                       var=_copy_annotations(snapshot.var), uns=copy.deepcopy(snapshot.uns))
    path = _archive_path(result)
    _verify_archive(path, result.input_archive)
    return read_h5ad(path, backed=backed)
