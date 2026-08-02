"""
Expression preprocessing helpers.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy import sparse


def _validated_celltype_labels(adata, celltype_col: str) -> pd.Series:
    """Return non-empty cell-type labels without converting missing values to text."""
    if celltype_col not in adata.obs:
        raise KeyError(f"{celltype_col!r} not found in adata.obs")
    values = adata.obs[celltype_col]
    text = values.astype("string").str.strip()
    if values.isna().any() or text.isna().any() or text.eq("").any():
        raise ValueError(f"adata.obs[{celltype_col!r}] must not contain NA or empty labels")
    return text.astype(str)


def _validated_gene_names(adata) -> pd.Index:
    """Return unique, non-empty expression-matrix gene names."""
    raw = pd.Index(adata.var_names)
    values = pd.Series(raw, dtype="object")
    text = values.astype("string").str.strip()
    if values.isna().any() or text.isna().any() or text.eq("").any():
        raise ValueError("adata.var_names must not contain NA or empty gene names")
    genes = pd.Index(text.astype(str))
    if genes.has_duplicates:
        duplicates = genes[genes.duplicated()].unique().tolist()
        raise ValueError(
            "adata.var_names must be unique; duplicated genes: "
            + ", ".join(duplicates[:10])
        )
    return genes


def _as_1d_array(x) -> np.ndarray:
    """Convert a sliced aggregation result to a flat dense vector."""
    if hasattr(x, "A1"):
        return x.A1
    return np.asarray(x).ravel()


def _all_celltype_counts(
    adata,
    celltype_col: str = "cell_type",
) -> pd.Series:
    """Return counts for every observed cell type used in score calculation."""
    labels = _validated_celltype_labels(adata, celltype_col)
    group_counts = labels.value_counts()
    if group_counts.empty:
        raise ValueError("No observed cell types are available for analysis")
    return group_counts.astype(int)


def _compute_celltype_fractions(
    adata,
    celltype_col: str = "cell_type",
) -> pd.Series:
    """Return every observed cell-type count divided by all annotated cells."""
    total_cells = int(len(adata.obs))
    if total_cells <= 0:
        raise ValueError("Cannot compute cell-type fractions from an empty AnnData object")
    counts = _all_celltype_counts(adata, celltype_col)
    fractions = counts.astype(float) / float(total_cells)
    fractions.name = "cell_fraction"
    return fractions


def _build_celltype_pseudobulk(
    adata,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
) -> pd.DataFrame:
    """
    构建细胞类型的 pseudobulk 表达矩阵
    """
    X = adata.layers[layer] if layer is not None else adata.X
    genes = _validated_gene_names(adata)
    labels = _validated_celltype_labels(adata, celltype_col)

    # Every observed type contributes to pseudobulk construction.
    valid_groups = _all_celltype_counts(adata, celltype_col).index.tolist()

    pseudobulk = []
    group_names = []
    for group in valid_groups:
        idx = labels.values == group
        pseudobulk.append(_as_1d_array(X[idx, :].mean(axis=0)))
        group_names.append(group)

    return pd.DataFrame(np.vstack(pseudobulk), index=group_names, columns=genes)


def _compute_celltype_expr_frac(
    adata,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
) -> pd.DataFrame:
    """
    计算每个基因在每个细胞类型中的表达比例（表达>0的细胞比例）
    """
    X = adata.layers[layer] if layer is not None else adata.X
    genes = _validated_gene_names(adata)
    labels = _validated_celltype_labels(adata, celltype_col)

    # Expression prevalence is calculated for every observed cell type.
    valid_groups = _all_celltype_counts(adata, celltype_col).index.tolist()

    expr_frac = []
    group_names = []
    for group in valid_groups:
        idx = labels.values == group
        n_cells = idx.sum()
        group_x = X[idx, :]
        if sparse.issparse(group_x):
            frac = group_x.getnnz(axis=0) / n_cells
        else:
            frac = (np.asarray(group_x) > 0).sum(axis=0) / n_cells
        expr_frac.append(_as_1d_array(frac))
        group_names.append(group)

    return pd.DataFrame(np.vstack(expr_frac), index=group_names, columns=genes)
