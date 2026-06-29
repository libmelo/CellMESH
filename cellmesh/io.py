"""
Single-cell data loading helpers.
"""
from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd


def _read_metadata_table(path: Union[str, Path], id_col: Optional[str] = None) -> pd.DataFrame:
    metadata = pd.read_csv(path)
    if id_col is not None:
        if id_col not in metadata.columns:
            raise ValueError(f"{id_col!r} not found in {path}")
        metadata = metadata.set_index(id_col)
    else:
        metadata = metadata.set_index(metadata.columns[0])
    metadata.index = metadata.index.astype(str)
    return metadata


def _read_name_list(path: Union[str, Path], prefer_second_column: bool = False) -> list[str]:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else None
    table = pd.read_csv(path, sep=sep, engine="python", header=None, comment="#")
    column = 1 if prefer_second_column and table.shape[1] > 1 else 0
    return table.iloc[:, column].astype(str).tolist()


def read_anndata(
    path: Union[str, Path],
    mode: Literal["h5ad", "10x", "csv", "tsv", "loom", "mtx"] = "h5ad",
    **kwargs,
):
    """
    从各种格式读取 AnnData 对象
    """
    path = Path(path)
    if mode == "h5ad":
        return _read_h5ad(path, **kwargs)
    if mode == "10x":
        return _read_10x(path, **kwargs)
    if mode == "csv":
        return _read_csv(path, **kwargs)
    if mode == "tsv":
        return _read_tsv(path, **kwargs)
    if mode == "loom":
        return _read_loom(path, **kwargs)
    if mode == "mtx":
        return _read_mtx(path, **kwargs)
    raise ValueError(f"不支持的读取模式: {mode}")


def _read_h5ad(path: Path, **kwargs):
    try:
        import anndata
    except ImportError:
        raise ImportError("读取 h5ad 文件需要 anndata 包")

    return anndata.read_h5ad(path, **kwargs)


def _read_10x(path: Path, gex_only: bool = True, **kwargs):
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError("读取 10X 数据需要 scanpy 包")

    return sc.read_10x_mtx(path, gex_only=gex_only, **kwargs)


def _read_csv(
    path: Path,
    cell_meta_path: Optional[Union[str, Path]] = None,
    gene_meta_path: Optional[Union[str, Path]] = None,
    cell_id_col: Optional[str] = None,
    transpose: bool = False,
    **kwargs,
):
    try:
        import anndata
    except ImportError:
        raise ImportError("读取 CSV 文件需要 anndata 包")

    df = pd.read_csv(path, index_col=0, **kwargs)
    if transpose:
        df = df.T
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    obs = None
    if cell_meta_path is not None:
        obs = _read_metadata_table(cell_meta_path, id_col=cell_id_col).reindex(df.index)

    var = None
    if gene_meta_path is not None:
        var = _read_metadata_table(gene_meta_path).reindex(df.columns)

    return anndata.AnnData(df, obs=obs, var=var)


def _read_tsv(path: Path, **kwargs):
    kwargs.setdefault("sep", "\t")
    return _read_csv(path, **kwargs)


def _read_loom(path: Path, **kwargs):
    try:
        import anndata
    except ImportError:
        raise ImportError("读取 Loom 文件需要 anndata 包")

    return anndata.read_loom(path, **kwargs)


def _read_mtx(
    path: Path,
    genes_path: Optional[Union[str, Path]] = None,
    barcodes_path: Optional[Union[str, Path]] = None,
    **kwargs,
):
    try:
        import anndata
        from scipy.io import mmread
    except ImportError:
        raise ImportError("读取 mtx 文件需要 anndata 和 scipy 包")

    mat = mmread(path).tocsr()
    adata = anndata.AnnData(mat.T)
    if barcodes_path is not None:
        barcodes = _read_name_list(barcodes_path)
        if len(barcodes) != adata.n_obs:
            raise ValueError(f"barcodes_path has {len(barcodes)} entries, expected {adata.n_obs}")
        adata.obs_names = barcodes
    if genes_path is not None:
        genes = _read_name_list(genes_path, prefer_second_column=True)
        if len(genes) != adata.n_vars:
            raise ValueError(f"genes_path has {len(genes)} entries, expected {adata.n_vars}")
        adata.var_names = genes
    return adata


def read_example_data(dataset: Literal["tiny", "small", "medium"] = "tiny"):
    """
    读取示例数据用于测试
    """
    try:
        import anndata
    except ImportError:
        raise ImportError("需要 anndata 包")

    rng = np.random.default_rng(42)

    if dataset == "tiny":
        n_cells, n_genes = 50, 50
        cell_types = ["A", "B", "C"]
    elif dataset == "small":
        n_cells, n_genes = 200, 100
        cell_types = ["Neutrophil", "Neuron", "Microglia", "T_cell"]
    elif dataset == "medium":
        n_cells, n_genes = 500, 200
        cell_types = ["Neutrophil", "Neuron", "Microglia", "T_cell", "B_cell", "Macrophage"]
    else:
        raise ValueError(f"不支持的数据集: {dataset}")

    X = rng.poisson(0.1, size=(n_cells, n_genes)).astype(float)
    gene_names = [f"Gene{i+1}" for i in range(n_genes)]
    cell_type_labels = rng.choice(cell_types, size=n_cells)

    for i, ct in enumerate(cell_types):
        if i < 5:
            ct_idx = cell_type_labels == ct
            gene_idx = i * 5 + np.arange(3)
            gene_idx = gene_idx[gene_idx < n_genes]
            X[np.ix_(ct_idx, gene_idx)] += rng.poisson(2, size=(ct_idx.sum(), len(gene_idx)))

    return anndata.AnnData(
        X,
        var=pd.DataFrame(index=gene_names),
        obs=pd.DataFrame({"cell_type": cell_type_labels, "sample": ["Sample1"] * n_cells}),
    )
