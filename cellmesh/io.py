"""
Single-cell data loading helpers.
"""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd

from ._table_io import _delimited_records, _record_location, _reject_nul_text


def _checked_text_labels(values, *, context: str, allow_first_blank: bool = False) -> pd.Index:
    """Keep literal identifiers while rejecting empty or ambiguous normalized names."""
    labels = pd.Series(list(values), dtype="object")
    normalized = labels.astype("string").str.strip()
    empty = normalized.isna() | normalized.eq("")
    if allow_first_blank and len(empty):
        empty.iloc[0] = False
    if empty.any():
        positions = (np.flatnonzero(empty.to_numpy(dtype=bool)) + 1).tolist()
        raise ValueError(f"{context}: empty identifiers at positions {positions[:10]}")
    repeated = normalized.duplicated(keep=False)
    if repeated.any():
        details = []
        for name in normalized[repeated].drop_duplicates().head(10):
            positions = (np.flatnonzero(normalized.eq(name).to_numpy(dtype=bool)) + 1).tolist()
            details.append(f"{name!r} at positions {positions}")
        raise ValueError(f"{context}: duplicate identifiers after stripping whitespace: " + "; ".join(details))
    return pd.Index(labels.astype(str), name=getattr(values, "name", None))


def _csv_columns(path: Path, options: dict, *, context: str) -> pd.Index:
    """Check original headers before pandas can rename duplicates to .1/.2."""
    _reject_nul_text(path, **{key: options[key] for key in
                             ("encoding", "compression", "encoding_errors", "storage_options")
                             if key in options})
    if options.get("chunksize") is not None or options.get("iterator", False):
        raise ValueError("AnnData CSV/TSV loading requires a complete table, not an iterator")
    if options.get("index_col") is not None:
        raise ValueError("AnnData CSV/TSV loading uses the first selected column as its identifier")
    if options.get("engine") == "pyarrow":
        raise ValueError("AnnData CSV/TSV header validation requires engine='c' or engine='python'")
    if options.get("parse_dates") is not None and options.get("parse_dates") is not False:
        raise ValueError("AnnData CSV/TSV loading preserves text identifiers; parse_dates is not supported")

    header = options.get("header", "infer")
    names = options.get("names")
    if header == "infer":
        header = 0 if names is None else None
    if header is not None and (
        isinstance(header, (bool, np.bool_))
        or not isinstance(header, (int, np.integer))
        or header < 0
    ):
        raise ValueError("AnnData CSV/TSV loading requires one header row or header=None")

    # Keep delimiter, quoting, compression, encoding, comments, and skipped rows
    # identical to the data read. Only data conversion/selection is disabled.
    probe = dict(options)
    for key in ("names", "usecols", "dtype", "converters", "na_values", "dtype_backend"):
        probe.pop(key, None)
    probe.update(index_col=None, dtype=str, keep_default_na=False, na_filter=False,
                 parse_dates=False, skipfooter=0)
    if header is not None:
        # A preamble before header=1, for example, may have fewer fields than
        # the real header. Supply its width so pandas pads that preamble instead
        # of treating the header itself as an over-wide data row.
        width = len(pd.read_csv(path, **{**probe, "header": int(header), "nrows": 0}).columns)
        raw = pd.read_csv(path, **{
            **probe, "header": None, "names": list(range(width)), "nrows": int(header) + 1,
        })
        if len(raw) <= header:
            raise ValueError(f"{context}: header row {header} is missing")
        if not isinstance(raw.index, pd.RangeIndex):
            raise ValueError(f"{context}: use skiprows for preamble rows wider than the header")
        _checked_text_labels(raw.iloc[int(header)], context=f"{context} header", allow_first_blank=True)
    if names is not None:
        _checked_text_labels(names, context=f"{context} supplied names", allow_first_blank=True)

    schema = pd.read_csv(path, **{
        **probe, "header": options.get("header", "infer"), "names": names,
        "usecols": options.get("usecols"), "nrows": 0,
    })
    if len(schema.columns) < 1:
        raise ValueError(f"{context}: no identifier column is available")
    return schema.columns


def _read_text_matrix(path: Path, options: dict) -> pd.DataFrame:
    """Read numeric expression while keeping the axis identifier as literal text."""
    columns = _csv_columns(path, options, context=f"Expression file {path}")
    identifier = columns[0]
    converters = dict(options.get("converters") or {})
    if identifier in converters or 0 in converters:
        raise ValueError("A converter cannot replace the expression matrix identifier column")
    # The Python parser still applies NA-token matching after a text converter.
    # A one-item tuple protects the literal during parsing; unwrap only after
    # the numeric columns have been read, before constructing the actual index.
    converters[identifier] = lambda value: (value,)
    read_options = {**options, "index_col": None, "converters": converters}
    dtype = read_options.get("dtype")
    if dtype is not None:
        if isinstance(dtype, Mapping):
            read_options["dtype"] = {key: value for key, value in dtype.items()
                                     if key not in (identifier, 0)}
        else:
            read_options["dtype"] = {column: dtype for column in columns[1:]}

    # 先将标识作为普通文本列读取，再设为索引。直接 index_col=0 即使有
    # converter，pandas 仍可能把空表头下的 01 变成 1，或把字面 NA 变成缺失。
    # 表达量继续使用数值解析；不能为方便将整个表达矩阵先读成字符串。
    # 回归测试：tests/test_io_identifiers.py。
    frame = pd.read_csv(path, **read_options)
    if not isinstance(frame.index, pd.RangeIndex):
        raise ValueError(f"Expression file {path}: data rows contain more fields than the header")
    frame[identifier] = frame[identifier].map(lambda value: value[0])
    frame = frame.set_index(identifier)
    # AnnData's conversion of nullable numeric columns can produce object X,
    # even without missing values. Convert columns BEFORE transpose (which can
    # itself discard extension dtype information). Never coerce arbitrary text,
    # booleans or complex data into real expression. Keep integer precision when
    # no NA exists; genuine missing integers require a floating NaN representation.
    # Regression: tests/test_nullable_expression_reader.py.
    for column in frame:
        values = frame[column]
        dtype = values.dtype
        if isinstance(dtype, pd.api.extensions.ExtensionDtype) and dtype.kind in "iuf":
            numpy_dtype = getattr(dtype, "numpy_dtype", None)
            if numpy_dtype is None:
                raise ValueError(f"Unsupported expression dtype {dtype!r} in column {column!r}")
            if values.isna().any():
                target = np.float64 if dtype.kind in "iu" else numpy_dtype
                frame[column] = values.to_numpy(dtype=target, na_value=np.nan)
            else:
                frame[column] = values.to_numpy(dtype=numpy_dtype)
    return frame


def _read_metadata_table(path: Union[str, Path], id_col: Optional[str] = None) -> pd.DataFrame:
    path = Path(path)
    options = {"encoding": "utf-8-sig"}
    _csv_columns(path, options, context=f"Metadata file {path}")
    # 分组列名可由 run_cell_mesh() 任意指定，因此所有元数据默认保留文本。
    # 读取后再 astype(str) 无法恢复已合并的 01/1，也无法恢复被当作 NA 的文本。
    # 只有空字段是真正的缺失；计数、年龄等附加字段由调用者显式转为数值。
    metadata = pd.read_csv(path, **options, dtype=str, keep_default_na=False, na_values=[""])
    if not isinstance(metadata.index, pd.RangeIndex):
        raise ValueError(f"Metadata file {path}: data rows contain more fields than the header")
    if id_col is not None:
        if id_col not in metadata.columns:
            raise ValueError(f"{id_col!r} not found in {path}")
        metadata = metadata.set_index(id_col)
    else:
        metadata = metadata.set_index(metadata.columns[0])
    metadata.index = _checked_text_labels(metadata.index, context=f"Metadata file {path} row identifiers")
    return metadata


def _read_name_list(path: Union[str, Path], prefer_second_column: bool = False) -> list[str]:
    path = Path(path)
    _reject_nul_text(path, encoding="utf-8-sig")
    # Do not sniff a delimiter from a single-column barcode such as 001 or NA.
    # Inner suffixes also recognize compressed .csv.gz/.tsv.gz name files.
    sep = "," if ".csv" in [suffix.lower() for suffix in path.suffixes] else "\t"
    table = pd.read_csv(path, sep=sep, header=None, comment="#", encoding="utf-8-sig",
                        dtype=str, keep_default_na=False, na_values=[""], skip_blank_lines=False)
    column = 1 if prefer_second_column and table.shape[1] > 1 else 0
    return _checked_text_labels(table.iloc[:, column], context=f"Name file {path}").tolist()


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


def _10x_file(path: Path, prefix: str, stems: tuple[str, ...]) -> Path:
    candidates = [path / f"{prefix}{stem}{suffix}" for stem in stems for suffix in ("", ".gz")]
    matches = [candidate for candidate in candidates if candidate.is_file()]
    if not matches:
        raise FileNotFoundError("10X file not found; expected one of: " + ", ".join(map(str, candidates)))
    if len(matches) != 1:
        raise ValueError("10X directory has ambiguous files: " + ", ".join(map(str, matches)))
    return matches[0]


def _10x_annotations(path: Path, width: int, *, keep_locations=True) -> tuple[pd.DataFrame, list[str]]:
    rows, locations = [], []
    # A blank barcode/feature record occupies a matrix position; never skip it.
    with _delimited_records(path, delimiter="\t", skip_blank_lines=False) as records:
        for start, end, row in records:
            if len(row) != width:
                raise ValueError(
                    f"10X file {path}, {_record_location(start, end)}: "
                    f"expected {width} fields, got {len(row)}"
                )
            rows.append(row)
            if keep_locations:
                locations.append(_record_location(start, end))
    if not rows:
        raise ValueError(f"10X file {path}: expected at least one annotation record")
    return pd.DataFrame(rows, columns=range(width), dtype=object), locations


def _read_10x(
    path: Path, gex_only: bool = True, *, var_names: str = "gene_symbols",
    make_unique: bool = False, prefix: Optional[str] = None, cache: bool = False,
    **kwargs,
):
    """Read raw 10X labels before matrix loading; never manufacture gene IDs.

    Accept legacy genes.tsv and modern features.tsv, each plain or gzip, with
    optional prefixes. Scanpy still handles numerical matrix loading/caching;
    raw annotations are revalidated even when a matrix cache is used.
    """
    for name, value in (("gex_only", gex_only), ("make_unique", make_unique), ("cache", cache)):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"10X {name} must be a boolean")
    if make_unique:
        raise ValueError(
            "10X make_unique=True is not supported: automatic gene renaming hides "
            "ambiguous prior matches. Resolve duplicate symbols or use unique "
            "var_names='gene_ids' with matching prior gene identifiers."
        )
    if var_names not in ("gene_symbols", "gene_ids"):
        raise ValueError("10X var_names must be 'gene_symbols' or 'gene_ids'")
    if prefix is not None and not isinstance(prefix, str):
        raise TypeError("10X prefix must be a string or None")
    unexpected = set(kwargs) - {"cache_compression"}
    if unexpected:
        raise TypeError(f"Unsupported 10X reader arguments: {sorted(unexpected)}")
    prefix = "" if prefix is None else prefix
    matrix_path = _10x_file(path, prefix, ("matrix.mtx",))
    features_path = _10x_file(path, prefix, ("genes.tsv", "features.tsv"))
    barcodes_path = _10x_file(path, prefix, ("barcodes.tsv",))
    modern = features_path.name in (f"{prefix}features.tsv", f"{prefix}features.tsv.gz")
    features, locations = _10x_annotations(features_path, 3 if modern else 2)
    barcodes, _ = _10x_annotations(barcodes_path, 1, keep_locations=False)
    if modern and features[2].str.strip().eq("").any():
        raise ValueError(f"10X file {features_path}: feature types must not be empty")
    selected = features.loc[features[2].str.strip().eq("Gene Expression")] if modern and gex_only else features
    selected_column = 1 if var_names == "gene_symbols" else 0
    try:
        genes = _checked_text_labels(selected[selected_column], context=f"10X file {features_path} {var_names}")
    except ValueError as error:
        normalized = selected[selected_column].str.strip()
        bad = normalized.eq("") | normalized.duplicated(keep=False)
        details = [
            f"record {i + 1} ({locations[i]}): gene_id={row[0]!r}, symbol={row[1]!r}"
            for i, row in selected.loc[bad].head(10).iterrows()
        ]
        raise ValueError(f"{error}; raw feature records: " + "; ".join(details)) from error
    cells = _checked_text_labels(barcodes[0], context=f"10X file {barcodes_path} barcodes")
    # sc.read_10x_mtx defaults to make_unique=True and also infers label dtypes.
    # G/G -> G/G-1 would evade later duplicate checks, while 01/1/NA could become
    # 1/1/NaN. Keep raw validated labels; neither renaming nor summing is allowed.
    # 回归测试：tests/test_io_10x.py，包含所有布局及 cache 路径。
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError("读取 10X 数据需要 scanpy 包")

    data = sc.read(matrix_path, cache=bool(cache), **kwargs).T
    expected_shape = (len(barcodes), len(features))
    if data.shape != expected_shape:
        raise ValueError(
            f"10X matrix {matrix_path}: dimensions {data.shape} after transpose, "
            f"expected {expected_shape} from {barcodes_path} and {features_path}"
        )
    if len(selected) != len(features):
        data = data[:, selected.index.to_numpy()].copy()
    data.obs = pd.DataFrame(index=pd.Index(cells.to_numpy()))
    annotations = pd.DataFrame(index=pd.Index(genes.to_numpy()))
    other_column = "gene_ids" if var_names == "gene_symbols" else "gene_symbols"
    annotations[other_column] = selected[1 - selected_column].to_numpy()
    if modern:
        annotations["feature_types"] = selected[2].to_numpy()
    data.var = annotations
    return data


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

    df = _read_text_matrix(path, kwargs)
    if transpose:
        df = df.T
    df.index = _checked_text_labels(df.index, context=f"Expression file {path} cell identifiers")
    df.columns = _checked_text_labels(df.columns, context=f"Expression file {path} gene identifiers")

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
    except ImportError as error:
        raise ImportError("读取 Loom 文件需要 anndata 包") from error
    try:
        import loompy  # noqa: F401 -- optional reader dependency, checked explicitly
    except ModuleNotFoundError as error:
        if error.name != "loompy":
            raise
        raise ImportError(
            "读取 Loom 文件需要可选依赖 loompy；请安装 pip install 'cellmesh[loom]' "
            "或 pip install loompy"
        ) from error
    # anndata.io is the supported entry in newer AnnData; older releases keep
    # their top-level reader. Do not catch corrupt-file/reader errors as imports.
    from importlib.util import find_spec

    if find_spec("anndata.io") is not None:
        from anndata.io import read_loom
    else:
        read_loom = anndata.read_loom
    return read_loom(path, **kwargs)


def _read_mtx(
    path: Path,
    genes_path: Optional[Union[str, Path]] = None,
    barcodes_path: Optional[Union[str, Path]] = None,
):
    try:
        import anndata
        from scipy import sparse
        from scipy.io import mmread
    except ImportError:
        raise ImportError("读取 mtx 文件需要 anndata 和 scipy 包")

    mat = mmread(path)
    # Matrix Market 的 coordinate 返回稀疏矩阵，array 返回 NumPy 数组。
    # 不能假定 mmread() 的结果都有 .tocsr()；按类型处理，保留稀疏存储以
    # 避免大矩阵被展开，同时保留稠密数组的值和维度。转置和名称校验共用。
    # 回归测试：tests/test_io_mtx.py。
    mat = mat.tocsr() if sparse.issparse(mat) else np.asarray(mat)
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
