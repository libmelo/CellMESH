"""
Expression preprocessing helpers.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy import sparse


def _normalized_label_series(
    values: pd.Series,
    name: str,
    *,
    reject_collisions: bool = True,
) -> pd.Series:
    """Normalize identifiers without mutating input or silently merging labels."""
    text = values.astype("string").str.strip()
    if values.isna().any() or text.isna().any() or text.eq("").any():
        raise ValueError(f"{name} must not contain NA or empty labels")
    if reject_collisions:
        # Repeated cells in one group are expected. Only distinct observed raw
        # labels that collapse to the same text are ambiguous; unused categories
        # must not trigger a collision. This also catches integer 1 versus "1".
        # Include normalized text in the deduplication key: pandas considers
        # numeric 1 and 1.0 equal, but their string identifiers differ.
        distinct = pd.DataFrame({"raw": values, "text": text}).drop_duplicates()
        collisions = distinct.loc[distinct["text"].duplicated(), "text"].unique().tolist()
        if collisions:
            raise ValueError(
                f"{name} has labels that collide after normalization: {collisions}"
            )
    return text.astype(str)


def _validated_obs_labels(adata, key: str) -> pd.Series:
    """Return canonical labels for either cell types or samples."""
    if key not in adata.obs:
        raise KeyError(f"{key!r} not found in adata.obs")
    return _normalized_label_series(adata.obs[key], f"adata.obs[{key!r}]")


def _validated_celltype_labels(adata, celltype_col: str) -> pd.Series:
    """Return canonical cell-type labels shared by observed and permuted scoring."""
    return _validated_obs_labels(adata, celltype_col)


def _normalized_gene_names(var_names) -> pd.Index:
    """Normalize expression-gene identifiers while preserving column positions."""
    raw = pd.Index(var_names)
    values = pd.Series(raw, dtype="object")
    genes = pd.Index(_normalized_label_series(values, "adata.var_names", reject_collisions=False))
    if genes.has_duplicates:
        duplicates = genes[genes.duplicated()].unique().tolist()
        raise ValueError(
            "adata.var_names must be unique; duplicated genes: "
            + ", ".join(duplicates[:10])
        )
    return genes


def _validated_gene_names(adata) -> pd.Index:
    """Return unique, non-empty expression-matrix gene names."""
    return _normalized_gene_names(adata.var_names)


_EXPRESSION_CHECK_BLOCK_SIZE = 1_000_000


def _validate_expression_values(X, *, layer: Optional[str] = None) -> None:
    """Check selected expression values without densifying sparse matrices."""
    source = "adata.X" if layer is None else f"adata.layers[{layer!r}]"
    message = (
        "scoring expression must be finite and non-negative real values; "
        f"check {source}"
    )
    values = X.data if sparse.issparse(X) else np.asarray(X)
    if values.dtype.kind not in "biuf":
        raise ValueError(message)
    # Sparse implicit zeros are valid. Check stored values without modifying
    # sparse structure; bound temporary arrays for dense and sparse input alike.
    for start in range(0, values.size, _EXPRESSION_CHECK_BLOCK_SIZE):
        block = values.flat[start:start + _EXPRESSION_CHECK_BLOCK_SIZE]
        if np.any(~np.isfinite(block)) or np.any(block < 0.0):
            raise ValueError(message)


def _validate_scoring_expression(adata, genes, layer: Optional[str] = None) -> None:
    """Validate measured scoring genes in the selected layer before aggregation."""
    measured_genes = _validated_gene_names(adata)
    requested = list(dict.fromkeys(str(gene) for gene in genes))
    columns = measured_genes.get_indexer(requested)
    columns = columns[columns >= 0]
    source = adata.layers[layer] if layer is not None else adata.X
    if not len(columns):
        return

    # 必须检查聚合前的表达值，不能为提速改成仅检查 pseudobulk 或 P/C/E：
    # [-1, 3] 的均值为 1，负反应也可能被正反应抵消。聚合结果合法不代表
    # 输入合法，否则 n_perms=0 会放行，开启置换才报错。
    # 只检查实际参与评分的基因和所选表达层；不改值，不改先验反应定义。
    # 回归测试：tests/test_expression_validation.py。
    rows_per_block = max(1, _EXPRESSION_CHECK_BLOCK_SIZE // len(columns))
    for start in range(0, source.shape[0], rows_per_block):
        selected = source[start:start + rows_per_block, columns]
        _validate_expression_values(selected, layer=layer)


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


def _grouped_expression_mean(X, indicator: sparse.csr_matrix, counts) -> np.ndarray:
    """Compute group means with the same float64 accumulation in both scorers.

    CSR membership multiplication fixes the accumulation order across dense,
    CSR, and CSC expression inputs, including single-gene matrices. NumPy and
    sparse mean/sum reductions can otherwise use different summation orders.
    """
    # Observed pseudobulks also retain unrelated genes. Consumers check the
    # scoring columns before using them, so unrelated overflow does not reject
    # an otherwise valid analysis. Keep the same accumulation/division order.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        sums = indicator.astype(np.float64, copy=False) @ X
        sums = sums.toarray() if sparse.issparse(sums) else np.asarray(sums)
        denominators = np.asarray(counts, dtype=np.float64)[:, None]
        return np.asarray(sums, dtype=np.float64) / denominators


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
        group_x = X[idx, :]
        n_cells = group_x.shape[0]
        # Process one cell type at a time to bound expression temporaries.
        indicator = sparse.csr_matrix(
            (
                np.ones(n_cells, dtype=np.float64),
                np.arange(n_cells),
                np.array([0, n_cells]),
            ),
            shape=(1, n_cells),
        )
        mean = _grouped_expression_mean(group_x, indicator, [n_cells])
        pseudobulk.append(_as_1d_array(mean))
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
        # 表达比例必须按数值 > 0 计数，不能为提速替换为 getnnz()/nnz。
        # 稀疏矩阵可能显式存储零；按存储条目计数会抬高表达比例，
        # 改变 min_expr_frac 门槛，并使观测与置换的统计口径不一致。
        # 布尔比较和求和保持稀疏计算，仅将按基因汇总的计数转为向量。
        # 回归测试：tests/test_expression_fraction.py。
        positive_counts = _as_1d_array((group_x > 0).sum(axis=0))
        expr_frac.append(positive_counts / n_cells)
        group_names.append(group)

    return pd.DataFrame(np.vstack(expr_frac), index=group_names, columns=genes)
