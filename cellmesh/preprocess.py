"""
预处理模块
包含数据加载、先验验证和代谢物可用性计算等功能
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any, Union, Literal

import numpy as np
from scipy.stats import gmean
import pandas as pd
import warnings

# 导入集中配置
from .config import (
    MIN_CELL_COUNT,
    VALID_ROLES,
    VALID_SENSOR_TYPES,
    METABOLITE_AVAILABILITY_DEFAULTS
)


def _to_dense(x):
    """将稀疏矩阵转换为稠密矩阵"""
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


# ==================== 数据读取函数 ====================
def read_anndata(
    path: Union[str, Path],
    mode: Literal["h5ad", "10x", "csv", "tsv", "loom", "mtx"] = "h5ad",
    **kwargs
):
    """
    从各种格式读取 AnnData 对象

    参数:
        path: 数据路径
        mode: 数据格式
        **kwargs: 传递给具体读取函数的参数

    返回:
        AnnData 对象
    """
    path = Path(path)
    if mode == "h5ad":
        return _read_h5ad(path, **kwargs)
    elif mode == "10x":
        return _read_10x(path, **kwargs)
    elif mode == "csv":
        return _read_csv(path, **kwargs)
    elif mode == "tsv":
        return _read_tsv(path, **kwargs)
    elif mode == "loom":
        return _read_loom(path, **kwargs)
    elif mode == "mtx":
        return _read_mtx(path, **kwargs)
    else:
        raise ValueError(f"不支持的读取模式: {mode}")


def _read_h5ad(path: Path, **kwargs):
    """读取 h5ad 格式文件"""
    try:
        import anndata
    except ImportError:
        raise ImportError("读取 h5ad 文件需要 anndata 包")

    adata = anndata.read_h5ad(path, **kwargs)
    return adata


def _read_10x(path: Path, gex_only: bool = True, **kwargs):
    """读取 10X 格式数据"""
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError("读取 10X 数据需要 scanpy 包")

    adata = sc.read_10x_mtx(path, gex_only=gex_only, **kwargs)
    return adata


def _read_csv(
    path: Path,
    cell_meta_path: Optional[Union[str, Path]] = None,
    gene_meta_path: Optional[Union[str, Path]] = None,
    cell_id_col: Optional[str] = None,
    transpose: bool = False,
    **kwargs
):
    """读取 CSV 格式表达矩阵"""
    try:
        import anndata
    except ImportError:
        raise ImportError("读取 CSV 文件需要 anndata 包")

    df = pd.read_csv(path, index_col=0, **kwargs)
    if transpose:
        df = df.T
    adata = anndata.AnnData(df)
    return adata


def _read_tsv(path: Path, **kwargs):
    """读取 TSV 格式表达矩阵"""
    kwargs.setdefault('sep', '\t')
    return _read_csv(path, **kwargs)


def _read_loom(path: Path, **kwargs):
    """读取 Loom 格式文件"""
    try:
        import anndata
    except ImportError:
        raise ImportError("读取 Loom 文件需要 anndata 包")

    adata = anndata.read_loom(path, **kwargs)
    return adata


def _read_mtx(
    path: Path,
    genes_path: Optional[Union[str, Path]] = None,
    barcodes_path: Optional[Union[str, Path]] = None,
    **kwargs
):
    """读取 Mtx 格式表达矩阵"""
    try:
        import anndata
        from scipy.io import mmread
    except ImportError:
        raise ImportError("读取 mtx 文件需要 anndata 和 scipy 包")

    mat = mmread(path).tocsr()
    adata = anndata.AnnData(mat.T)
    return adata


def read_example_data(dataset: Literal["tiny", "small", "medium"] = "tiny"):
    """
    读取示例数据用于测试

    参数:
        dataset: 数据集大小,可选 "tiny", "small", "medium"

    返回:
        示例 AnnData 对象
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

    adata = anndata.AnnData(
        X,
        var=pd.DataFrame(index=gene_names),
        obs=pd.DataFrame({
            "cell_type": cell_type_labels,
            "sample": ["Sample1"] * n_cells
        })
    )
    return adata


# ==================== 先验验证 ====================
def validate_priors(
    enzyme_metabolite: pd.DataFrame,
    metabolite_sensor: pd.DataFrame,
    var_names
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    验证并标准化先验数据库

    参数:
        enzyme_metabolite: 酶-代谢物关系表
        metabolite_sensor: 代谢物-传感器关系表
        var_names: 数据集中的基因名列表

    返回:
        标准化后的 (enzyme_prior, sensor_prior) 元组
    """
    genes = set(pd.Index(var_names).astype(str))

    enz = enzyme_metabolite.copy()
    required_enz = {"metabolite", "gene", "role"}
    missing = required_enz - set(enz.columns)
    if missing:
        raise ValueError(f"enzyme_metabolite is missing columns: {sorted(missing)}")

    enz["gene"] = enz["gene"].astype(str)
    enz["role"] = enz["role"].astype(str).str.lower()
    enz = enz[enz["role"].isin(VALID_ROLES)]
    enz = enz[enz["gene"].isin(genes)]

    if "weight" not in enz:
        enz["weight"] = 1.0
    enz["weight"] = pd.to_numeric(enz["weight"], errors="coerce").fillna(1.0)

    sen = metabolite_sensor.copy()
    required_sen = {"metabolite", "sensor_gene", "sensor_type"}
    missing = required_sen - set(sen.columns)
    if missing:
        raise ValueError(f"metabolite_sensor is missing columns: {sorted(missing)}")

    sen["sensor_gene"] = sen["sensor_gene"].astype(str)
    # 不再转换或过滤 sensor_type，保持原样
    sen = sen[sen["sensor_gene"].isin(genes)]

    if "weight" not in sen:
        sen["weight"] = 1.0
    sen["weight"] = pd.to_numeric(sen["weight"], errors="coerce").fillna(1.0)

    return enz.reset_index(drop=True), sen.reset_index(drop=True)


# ==================== 代谢物可用性计算模块 ====================
def _build_celltype_pseudobulk(
    adata,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_cells: int = MIN_CELL_COUNT
) -> pd.DataFrame:
    """
    构建细胞类型的 pseudobulk 表达矩阵

    参数:
        adata: AnnData 对象
        celltype_col: 细胞类型列名
        layer: 使用的层,None 表示使用 adata.X
        min_cells: 每个细胞类型的最小细胞数

    返回:
        细胞类型 × 基因 的表达矩阵
    """
    if celltype_col not in adata.obs:
        raise KeyError(f"{celltype_col!r} not found in adata.obs")

    X = adata.layers[layer] if layer is not None else adata.X
    X = _to_dense(X)
    genes = pd.Index(adata.var_names).astype(str)
    labels = adata.obs[celltype_col].astype(str)

    # 过滤细胞数不足的细胞类型
    group_counts = labels.value_counts()
    valid_groups = group_counts[group_counts >= min_cells].index.tolist()
    if not valid_groups:
        raise ValueError(f"No cell types with at least {min_cells} cells")

    # 计算每个细胞类型的平均表达
    pseudobulk = []
    group_names = []
    for group in valid_groups:
        idx = labels.values == group
        pseudobulk.append(X[idx, :].mean(axis=0))
        group_names.append(group)

    return pd.DataFrame(
        np.vstack(pseudobulk),
        index=group_names,
        columns=genes
    )


def _compute_celltype_expr_frac(
    adata,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_cells: int = MIN_CELL_COUNT
) -> pd.DataFrame:
    """
    计算每个基因在每个细胞类型中的表达比例（表达>0的细胞比例）

    参数:
        adata: AnnData 对象
        celltype_col: 细胞类型列名
        layer: 使用的层,None 表示使用 adata.X
        min_cells: 每个细胞类型的最小细胞数

    返回:
        细胞类型 × 基因 的表达比例矩阵
    """
    if celltype_col not in adata.obs:
        raise KeyError(f"{celltype_col!r} not found in adata.obs")

    X = adata.layers[layer] if layer is not None else adata.X
    X = _to_dense(X)
    genes = pd.Index(adata.var_names).astype(str)
    labels = adata.obs[celltype_col].astype(str)

    # 过滤细胞数不足的细胞类型
    group_counts = labels.value_counts()
    valid_groups = group_counts[group_counts >= min_cells].index.tolist()
    if not valid_groups:
        raise ValueError(f"No cell types with at least {min_cells} cells")

    # 计算每个细胞类型中每个基因的表达比例
    expr_frac = []
    group_names = []
    for group in valid_groups:
        idx = labels.values == group
        n_cells = idx.sum()
        expr_frac.append((X[idx, :] > 0).sum(axis=0) / n_cells)
        group_names.append(group)

    return pd.DataFrame(
        np.vstack(expr_frac),
        index=group_names,
        columns=genes
    )


def compute_sensor_scores(
    adata,
    sensor_prior: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    lower: float = METABOLITE_AVAILABILITY_DEFAULTS["lower"],
    upper: float = METABOLITE_AVAILABILITY_DEFAULTS["upper"],
    min_expr_frac: float = 0.05,
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
) -> pd.DataFrame:
    """
    计算 sensor scores: robust min-max 标准化的 sensor 基因表达

    参数:
        adata: AnnData 对象
        sensor_prior: sensor 先验表（来自 normalize_interaction_database）
        celltype_col: 细胞类型列名
        layer: 使用的层
        lower: robust minmax 的下限百分位数
        upper: robust minmax 的上限百分位数
        min_expr_frac: 最小表达比例阈值
        min_cells: 每个细胞类型的最小细胞数

    返回:
        sensor scores DataFrame，包含列：
        metabolite, hmdb_id, sensor_gene, sensor_type, receiver, 
        sensor_score, sensor_expr_frac
    """
    # 获取 pseudobulk 表达和表达比例
    pseudobulk = _build_celltype_pseudobulk(adata, celltype_col, layer, min_cells)
    expr_frac = _compute_celltype_expr_frac(adata, celltype_col, layer, min_cells)
    
    # 筛选在数据中存在的 sensor genes
    valid_genes = [g for g in sensor_prior["sensor_gene"].unique() if g in pseudobulk.columns]
    if not valid_genes:
        return pd.DataFrame(columns=[
            "metabolite", "hmdb_id", "sensor_gene", "sensor_type", 
            "receiver", "sensor_score", "sensor_expr_frac"
        ])
    
    # 对每个 sensor gene 计算 robust min-max 标准化
    sensor_gene_scores = {}
    for gene in valid_genes:
        expr_values = pseudobulk[gene].values
        norm_scores = robust_minmax(expr_values, lower=lower, upper=upper)
        sensor_gene_scores[gene] = pd.Series(norm_scores, index=pseudobulk.index)
    
    # 构建结果
    rows = []
    for _, row in sensor_prior.iterrows():
        gene = row["sensor_gene"]
        if gene not in valid_genes:
            continue
        
        for receiver in pseudobulk.index:
            frac = expr_frac.loc[receiver, gene]
            score = sensor_gene_scores[gene].loc[receiver] if frac >= min_expr_frac else 0.0
            
            rows.append({
                "metabolite": row["metabolite"],
                "hmdb_id": row["hmdb_id"],
                "sensor_gene": gene,
                "sensor_type": row["sensor_type"],
                "receiver": receiver,
                "sensor_score": score,
                "sensor_expr_frac": frac
            })
    
    return pd.DataFrame(rows)


def _parse_reaction_table(reaction_table: pd.DataFrame) -> pd.DataFrame:
    """
    解析反应表,处理多基因和重复反应

    参数:
        reaction_table: 包含 metabolite, hmdb_id, reaction, gene, direction, weight 列的 DataFrame

    返回:
        解析后的 DataFrame
    """
    df = reaction_table.copy()

    # 确保必需列存在
    required_cols = ['metabolite', 'gene', 'direction']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Reaction table missing required column: {col}")

    # 处理可选字段默认值
    if 'hmdb_id' not in df.columns:
        df['hmdb_id'] = np.nan
    if 'reaction' not in df.columns:
        df['reaction'] = 'unknown'
    if 'weight' not in df.columns:
        df['weight'] = 1.0

    # 标准化字段类型
    df['weight'] = pd.to_numeric(df['weight'], errors='coerce').fillna(1.0)
    df['hmdb_id'] = df['hmdb_id'].astype(object).where(pd.notna(df['hmdb_id']), np.nan)

    return df


def _get_reaction_gene_sets(reaction_table: pd.DataFrame) -> pd.DataFrame:
    """
    为每个反应提取基因集合,处理多基因情况,去重,支持各种分隔符

    参数:
        reaction_table: 解析后的反应表

    返回:
        每个反应对应的基因集合
    """
    df = reaction_table.copy()

    # 创建反应唯一标识
    def get_reaction_id(row):
        hmdb_str = str(row['hmdb_id']) if pd.notna(row['hmdb_id']) else 'nan'
        return f"{row['metabolite']}|{hmdb_str}|{row['reaction']}|{row['direction']}"

    df['_reaction_id'] = df.apply(get_reaction_id, axis=1)

    # 解析基因列(可能包含多个基因,用分号、逗号、竖线分隔)
    def parse_genes(gene_str):
        if pd.isna(gene_str) or gene_str == '':
            return []
        # 尝试多种分隔符
        for sep in [';', ',', '|']:
            if sep in gene_str:
                genes = [g.strip() for g in gene_str.split(sep)]
                # 去除可能的 "[Enzyme]"、"[Transport]" 等注释后缀
                genes = [g.split('[')[0].strip() for g in genes]
                return [g for g in genes if g]
        # 单个基因的情况
        gene = gene_str.split('[')[0].strip()
        return [gene] if gene else []

    df['_genes'] = df['gene'].apply(parse_genes)

    # 合并同一反应的基因和权重,去重基因
    reaction_dict = {}
    for _, row in df.iterrows():
        rid = row['_reaction_id']
        genes = row['_genes']
        weight = row['weight']

        if rid not in reaction_dict:
            reaction_dict[rid] = {
                'metabolite': row['metabolite'],
                'hmdb_id': row['hmdb_id'],
                'reaction': row['reaction'],
                'direction': row['direction'],
                'gene_weights': {}  # 存储 {gene: max_weight},避免重复基因
            }

        # 合并基因,同一基因多次出现取最大权重,并报警告
        for gene in genes:
            if gene in reaction_dict[rid]['gene_weights']:
                if weight > reaction_dict[rid]['gene_weights'][gene]:
                    reaction_dict[rid]['gene_weights'][gene] = weight
                    warnings.warn(f"Gene {gene} appeared multiple times in reaction {rid}, using maximum weight {weight}")
            else:
                reaction_dict[rid]['gene_weights'][gene] = weight

    # 转换为 DataFrame 格式
    result = []
    for rid, info in reaction_dict.items():
        result.append({
            'metabolite': info['metabolite'],
            'hmdb_id': info['hmdb_id'],
            'reaction': info['reaction'],
            'direction': info['direction'],
            'genes': list(info['gene_weights'].keys()),
            'weights': list(info['gene_weights'].values())
        })

    return pd.DataFrame(result)


def _compute_reaction_scores(
    pseudobulk: pd.DataFrame,
    reaction_genes: pd.DataFrame
) -> pd.DataFrame:
    """
    计算每个反应在每个细胞类型中的得分,考虑基因权重

    参数:
        pseudobulk: 细胞类型 × 基因 的表达矩阵
        reaction_genes: 反应-基因映射表,包含 genes 和 weights 列

    返回:
        反应 × 细胞类型 的得分矩阵
    """
    scores = []
    reaction_ids = []

    for _, row in reaction_genes.iterrows():
        # 创建反应标识
        hmdb_str = str(row['hmdb_id']) if pd.notna(row['hmdb_id']) else 'nan'
        rid = f"{row['metabolite']}|{hmdb_str}|{row['reaction']}|{row['direction']}"
        reaction_ids.append(rid)

        # 找到在数据集中存在的基因和对应的权重
        valid_gene_idx = [i for i, g in enumerate(row['genes']) if g in pseudobulk.columns]

        if not valid_gene_idx:
            # 没有有效基因,分数为 0
            scores.append(pd.Series(0.0, index=pseudobulk.index))
        else:
            valid_genes = [row['genes'][i] for i in valid_gene_idx]
            valid_weights = np.array([row['weights'][i] for i in valid_gene_idx])

            # 提取基因表达,乘以权重
            expr = pseudobulk[valid_genes].values * valid_weights.reshape(1, -1)

            # 计算加权几何平均数作为反应得分
            # 几何平均数公式:gmean(x + 1) - 1,避免 0 值影响
            geo_mean = gmean(expr + 1, axis=1) - 1

            scores.append(pd.Series(geo_mean, index=pseudobulk.index))

    return pd.DataFrame(scores, index=reaction_ids, columns=pseudobulk.index)


def _compute_PCE_matrices(
    reaction_scores: pd.DataFrame,
    reaction_genes: pd.DataFrame
) -> Dict[str, pd.DataFrame]:
    """
    计算 P (production), C (consumption), E (efflux) 矩阵

    参数:
        reaction_scores: 反应 × 细胞类型 的得分矩阵
        reaction_genes: 反应-基因映射表

    返回:
        包含 P, C, E 矩阵的字典
    """
    # 为每个反应添加方向信息和代谢物信息
    reaction_info = {}
    metabolites_dict = {}  # 使用字典来跟踪唯一的代谢物

    for _, row in reaction_genes.iterrows():
        hmdb_str = str(row['hmdb_id']) if pd.notna(row['hmdb_id']) else 'nan'
        rid = f"{row['metabolite']}|{hmdb_str}|{row['reaction']}|{row['direction']}"
        reaction_info[rid] = {
            'direction': row['direction'],
            'metabolite': row['metabolite'],
            'hmdb_id': row['hmdb_id']
        }
        # 使用字符串来作为唯一键，避免 np.nan != np.nan 的问题
        met_key = f"{row['metabolite']}|{hmdb_str}"
        if met_key not in metabolites_dict:
            metabolites_dict[met_key] = (row['metabolite'], row['hmdb_id'])

    # 初始化矩阵
    cell_types = reaction_scores.columns
    metabolites = list(metabolites_dict.values())
    # 使用唯一的代谢物对构建索引
    met_index = pd.MultiIndex.from_tuples(metabolites, names=['metabolite', 'hmdb_id'])

    P = pd.DataFrame(0.0, index=met_index, columns=cell_types)
    C = pd.DataFrame(0.0, index=met_index, columns=cell_types)
    E = pd.DataFrame(0.0, index=met_index, columns=cell_types)

    # 填充矩阵
    for rid in reaction_scores.index:
        info = reaction_info[rid]
        direction = info['direction']
        met = info['metabolite']
        hmdb = info['hmdb_id']
        hmdb_str = str(hmdb) if pd.notna(hmdb) else 'nan'
        met_key = f"{met}|{hmdb_str}"
        met_idx = metabolites_dict[met_key]

        if direction == 'product':
            P.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == 'substrate':
            C.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == 'exporter':
            E.loc[met_idx] += reaction_scores.loc[rid]

    return {'P': P, 'C': C, 'E': E}


def robust_minmax(
    x: np.ndarray,
    lower: float = 5,
    upper: float = 95,
    eps: float = 1e-8
) -> np.ndarray:
    """
    Robust min-max 标准化,结果范围在 [0, 1] 之间

    参数:
        x: 输入数组
        lower: 下限百分位数
        upper: 上限百分位数
        eps: 防止除零的小值

    返回:
        标准化后的数组,范围 [0, 1]
    """
    x = np.asarray(x, dtype=float)

    # 保存 NaN 的位置
    nan_mask = np.isnan(x)

    # 处理全 NaN 的情况
    if np.all(nan_mask):
        return np.zeros_like(x)

    # 计算百分位数(忽略 NaN)
    lo = np.nanpercentile(x, lower)
    hi = np.nanpercentile(x, upper)

    # 如果上下限接近,返回全 0
    if np.isclose(hi, lo, rtol=1e-5):
        result = np.zeros_like(x)
    else:
        # 截断并标准化到 [0, 1]
        x_clip = np.clip(x, lo, hi)
        result = (x_clip - lo) / (hi - lo + eps)
        # 确保严格在 0-1 范围内
        result = np.clip(result, 0.0, 1.0)

    # 将原来的 NaN 位置填充为 0
    result[nan_mask] = 0.0

    return result


def _safe_hmdb_compare(row: pd.Series, met: str, hmdb: Optional[str]) -> bool:
    """
    安全比较代谢物和 hmdb_id,处理 hmdb_id 为 NaN 的情况

    参数:
        row: reaction_genes 中的行
        met: 代谢物名称
        hmdb: 代谢物 HMDB ID

    返回:
        是否匹配
    """
    if row["metabolite"] != met:
        return False
    # hmdb 都是 NaN,视为匹配
    if pd.isna(row["hmdb_id"]) and pd.isna(hmdb):
        return True
    # 字符串相等才匹配
    return str(row["hmdb_id"]) == str(hmdb)


def _normalize_PCE(
    P: pd.DataFrame,
    C: pd.DataFrame,
    E: pd.DataFrame,
    reaction_genes: pd.DataFrame,
    lower: float = 5,
    upper: float = 95,
    missing_C_norm: float = 0.2,
    missing_E_norm: float = 0.5
) -> Dict[str, pd.DataFrame]:
    """
    对 P, C, E 矩阵进行标准化，结果都在 [0, 1] 范围内

    参数:
        P, C, E: 原始矩阵
        reaction_genes: 反应-基因映射表（用于判断是否存在 substrate/exporter reaction）
        lower, upper: robust minmax 的百分位数
        missing_C_norm: 缺失 C 时的默认值，范围 [0, 1]
        missing_E_norm: 缺失 E 时的默认值，范围 [0, 1]

    返回:
        包含标准化后矩阵的字典
    """
    # 确保默认值在合法范围内
    missing_C_norm = np.clip(missing_C_norm, 0.0, 1.0)
    missing_E_norm = np.clip(missing_E_norm, 0.0, 1.0)

    # 初始化结果 DataFrame
    P_norm = pd.DataFrame(index=P.index, columns=P.columns, dtype=float)
    C_norm = pd.DataFrame(missing_C_norm, index=P.index, columns=P.columns, dtype=float)
    E_norm = pd.DataFrame(missing_E_norm, index=P.index, columns=P.columns, dtype=float)

    # 对每个代谢物进行标准化
    for met_idx in P.index:
        met, hmdb = met_idx

        # P 标准化：范围 [0, 1]
        p_vals = P.loc[met_idx].values.flatten()
        P_norm.loc[met_idx] = robust_minmax(p_vals, lower=lower, upper=upper)

        # 检查是否有 substrate reaction
        has_substrate = any(
            _safe_hmdb_compare(row, met, hmdb) and row['direction'] == 'substrate'
            for _, row in reaction_genes.iterrows()
        )

        # C 标准化 - 只有当有 substrate reaction 时才计算，否则保持默认值
        if has_substrate and met_idx in C.index:
            c_vals = C.loc[met_idx].values.flatten()
            C_norm.loc[met_idx] = robust_minmax(c_vals, lower=lower, upper=upper)

        # 检查是否有 exporter reaction
        has_exporter = any(
            _safe_hmdb_compare(row, met, hmdb) and row['direction'] == 'exporter'
            for _, row in reaction_genes.iterrows()
        )

        # E 标准化 - 只有当有 exporter reaction 时才计算，否则保持默认值
        if has_exporter and met_idx in E.index:
            e_vals = E.loc[met_idx].values.flatten()
            E_norm.loc[met_idx] = robust_minmax(e_vals, lower=lower, upper=upper)

    return {'P_norm': P_norm, 'C_norm': C_norm, 'E_norm': E_norm}


def compute_metabolite_availability(
    adata,
    reaction_table: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    lower: float = METABOLITE_AVAILABILITY_DEFAULTS["lower"],
    upper: float = METABOLITE_AVAILABILITY_DEFAULTS["upper"],
    eps: float = METABOLITE_AVAILABILITY_DEFAULTS["eps"],
    beta: float = METABOLITE_AVAILABILITY_DEFAULTS["beta"],
    missing_C_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_C_norm"],
    missing_E_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_E_norm"],
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    return_intermediates: bool = True
) -> Dict[str, Any]:
    """
    计算代谢物 availability,结果范围在 [0, 1] 之间

    参数:
        adata: AnnData 对象
        reaction_table: 反应表(内部使用,由 enzyme_metabolite 转换而来)
        celltype_col: 细胞类型列名
        layer: 使用的层
        lower, upper: robust minmax 的百分位数
        eps: 公式中的小值,避免除零
        beta: 消耗项的指数权重
        missing_C_norm: 缺失 C 时的默认值,范围 [0, 1]
        missing_E_norm: 缺失 E 时的默认值,范围 [0, 1]
        min_cells: 每个细胞类型的最小细胞数
        return_intermediates: 是否返回中间结果

    返回:
        包含 availability 和中间结果的字典
        availability 值范围:[0, 1],值越高代表该细胞类型释放该代谢物的能力越强

    计算公式:
        availability = (P_norm) * ((1 - C_norm) ** beta) * (0.8 + 0.2 * E_norm)
        其中:
        - P_norm: 标准化后的产生能力,范围 [0, 1]
        - C_norm: 标准化后的消耗能力,范围 [0, 1]
        - E_norm: 标准化后的外排能力,范围 [0, 1]
        - 最终结果严格限制在 [0, 1] 范围内
    """
    # 确保 eps 和 beta 是正数
    eps = max(eps, 1e-8)
    beta = max(beta, 0.0)

    # 1. 构建 pseudobulk
    pseudobulk = _build_celltype_pseudobulk(
        adata,
        celltype_col=celltype_col,
        layer=layer,
        min_cells=min_cells
    )

    # 2. 解析反应表
    parsed_reactions = _parse_reaction_table(reaction_table)

    # 3. 获取反应-基因集合
    reaction_genes = _get_reaction_gene_sets(parsed_reactions)

    # 4. 计算反应得分
    reaction_scores = _compute_reaction_scores(pseudobulk, reaction_genes)

    # 5. 计算 P, C, E 矩阵
    PCE = _compute_PCE_matrices(reaction_scores, reaction_genes)
    P, C, E = PCE['P'], PCE['C'], PCE['E']

    # 6. 过滤没有 product reaction 的代谢物(即 P 全为 0 的代谢物)
    valid_mets = P.index[P.sum(axis=1) > 0]

    # 如果没有找到任何有效代谢物,返回空结果而不抛异常
    if len(valid_mets) == 0:
        result = {
            'availability': pd.DataFrame(),
            'metadata': pd.DataFrame()
        }
        if return_intermediates:
            result.update({
                'P': pd.DataFrame(),
                'C': pd.DataFrame(),
                'E': pd.DataFrame(),
                'P_norm': pd.DataFrame(),
                'C_norm': pd.DataFrame(),
                'E_norm': pd.DataFrame(),
                'pseudobulk': pd.DataFrame(),
                'reaction_genes': pd.DataFrame()
            })
        return result

    # 过滤所有矩阵，确保它们有相同的索引
    P = P.loc[valid_mets]
    # 确保 C 和 E 有相同的索引
    C = C.reindex(valid_mets, fill_value=0.0)
    E = E.reindex(valid_mets, fill_value=0.0)

    # 7. 标准化
    normalized = _normalize_PCE(
        P, C, E, reaction_genes,
        lower=lower,
        upper=upper,
        missing_C_norm=missing_C_norm,
        missing_E_norm=missing_E_norm
    )
    P_norm, C_norm, E_norm = normalized['P_norm'], normalized['C_norm'], normalized['E_norm']

    # 8. 计算最终 availability,确保结果在 [0, 1] 之间
    availability = pd.DataFrame(index=P_norm.index, columns=P_norm.columns, dtype=float)

    for met_idx in P_norm.index:
        p = P_norm.loc[met_idx].values.flatten()
        c = C_norm.loc[met_idx].values.flatten()
        e = E_norm.loc[met_idx].values.flatten()

        # 按公式计算
        avail = p * np.power((1 - c), beta) * (0.8 + 0.2 * e)

        # 严格限制在 0-1 范围内
        avail = np.clip(avail, 0.0, 1.0)

        availability.loc[met_idx] = avail

    # 9. 准备元数据
    metadata = pd.DataFrame(index=P_norm.index)
    metadata['has_product'] = True
    metadata['has_substrate'] = [
        any(
            _safe_hmdb_compare(row, met, hmdb) and row['direction'] == 'substrate'
            for _, row in reaction_genes.iterrows()
        )
        for met, hmdb in P_norm.index
    ]
    metadata['has_exporter'] = [
        any(
            _safe_hmdb_compare(row, met, hmdb) and row['direction'] == 'exporter'
            for _, row in reaction_genes.iterrows()
        )
        for met, hmdb in P_norm.index
    ]
    metadata['n_product_reactions'] = [
        sum(
            _safe_hmdb_compare(row, met, hmdb) and row['direction'] == 'product'
            for _, row in reaction_genes.iterrows()
        )
        for met, hmdb in P_norm.index
    ]
    metadata['n_substrate_reactions'] = [
        sum(
            _safe_hmdb_compare(row, met, hmdb) and row['direction'] == 'substrate'
            for _, row in reaction_genes.iterrows()
        )
        for met, hmdb in P_norm.index
    ]
    metadata['n_exporter_reactions'] = [
        sum(
            _safe_hmdb_compare(row, met, hmdb) and row['direction'] == 'exporter'
            for _, row in reaction_genes.iterrows()
        )
        for met, hmdb in P_norm.index
    ]

    # 10. 整理结果
    result = {
        'availability': availability.astype(float),
        'metadata': metadata
    }

    if return_intermediates:
        result.update({
            'P': P.astype(float),
            'C': C.astype(float),
            'E': E.astype(float),
            'P_norm': P_norm.astype(float),
            'C_norm': C_norm.astype(float),
            'E_norm': E_norm.astype(float),
            'pseudobulk': pseudobulk.astype(float),
            'reaction_genes': reaction_genes
        })

    return result
