"""
Expression preprocessing helpers.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy import sparse

from .config import MIN_CELL_COUNT


def _as_1d_array(x) -> np.ndarray:
    """Convert a sliced aggregation result to a flat dense vector."""
    if hasattr(x, "A1"):
        return x.A1
    return np.asarray(x).ravel()


def _eligible_celltype_counts(
    adata,
    celltype_col: str = "cell_type",
    min_cells: int = MIN_CELL_COUNT,
) -> pd.Series:
    """Return cell counts for cell types eligible for formal analysis."""
    if celltype_col not in adata.obs:
        raise KeyError(f"{celltype_col!r} not found in adata.obs")

    labels = adata.obs[celltype_col].astype(str)
    group_counts = labels.value_counts()
    eligible = group_counts[group_counts >= min_cells]
    if eligible.empty:
        raise ValueError(f"No cell types with at least {min_cells} cells")
    return eligible.astype(int)


def _build_celltype_pseudobulk(
    adata,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_cells: int = MIN_CELL_COUNT,
) -> pd.DataFrame:
    """
    构建细胞类型的 pseudobulk 表达矩阵
    """
    X = adata.layers[layer] if layer is not None else adata.X
    genes = pd.Index(adata.var_names).astype(str)
    labels = adata.obs[celltype_col].astype(str)

    valid_groups = _eligible_celltype_counts(adata, celltype_col, min_cells).index.tolist()

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
    min_cells: int = MIN_CELL_COUNT,
) -> pd.DataFrame:
    """
    计算每个基因在每个细胞类型中的表达比例（表达>0的细胞比例）
    """
    X = adata.layers[layer] if layer is not None else adata.X
    genes = pd.Index(adata.var_names).astype(str)
    labels = adata.obs[celltype_col].astype(str)

    valid_groups = _eligible_celltype_counts(adata, celltype_col, min_cells).index.tolist()

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
