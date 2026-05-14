from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any

import numpy as np
import pandas as pd


VALID_ROLES = {"production", "degradation", "export", "import", "usage"}
VALID_SENSOR_TYPES = {"surface_receptor", "transporter", "nuclear_receptor", "intracellular_sensor"}


@dataclass
class AggregatedExpression:
    mean_expr: pd.DataFrame
    expr_frac: pd.DataFrame


def _to_dense(x):
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def aggregate_expression(adata, groupby: str, layer: str | None = None, min_cells_per_group: int = 5) -> AggregatedExpression:
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")
    X = adata.layers[layer] if layer is not None else adata.X
    X = _to_dense(X)
    genes = pd.Index(adata.var_names).astype(str)
    obs = adata.obs.copy()
    labels = obs[groupby].astype(str)
    groups = labels.value_counts()
    keep_groups = groups[groups >= min_cells_per_group].index.tolist()
    if not keep_groups:
        raise ValueError("No cell groups remain after min_cells_per_group filtering")

    mean_rows = []
    frac_rows = []
    names = []
    for group in keep_groups:
        idx = np.flatnonzero(labels.values == group)
        sub = X[idx, :]
        mean_rows.append(sub.mean(axis=0))
        frac_rows.append((sub > 0).mean(axis=0))
        names.append(group)
    mean_expr = pd.DataFrame(np.vstack(mean_rows), index=names, columns=genes)
    expr_frac = pd.DataFrame(np.vstack(frac_rows), index=names, columns=genes)
    return AggregatedExpression(mean_expr=mean_expr, expr_frac=expr_frac)


def validate_priors(enzyme_metabolite: pd.DataFrame, metabolite_sensor: pd.DataFrame, var_names) -> tuple[pd.DataFrame, pd.DataFrame]:
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
    sen["sensor_type"] = sen["sensor_type"].astype(str).str.lower()
    sen = sen[sen["sensor_type"].isin(VALID_SENSOR_TYPES)]
    sen = sen[sen["sensor_gene"].isin(genes)]
    if "weight" not in sen:
        sen["weight"] = 1.0
    sen["weight"] = pd.to_numeric(sen["weight"], errors="coerce").fillna(1.0)
    return enz.reset_index(drop=True), sen.reset_index(drop=True)


# ==================== NEW METABOLITE AVAILABILITY MODULE ====================


def _build_celltype_pseudobulk(
    adata,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_cells: int = 1
) -> pd.DataFrame:
    """
    构建细胞类型的 pseudobulk 表达矩阵
    
    参数:
        adata: AnnData 对象
        celltype_col: 细胞类型列名
        layer: 使用的层，None 表示使用 adata.X
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


def _parse_reaction_table(reaction_table: pd.DataFrame) -> pd.DataFrame:
    """
    解析反应表，处理多基因和重复反应
    
    参数:
        reaction_table: 包含 metabolite, HMDB_ID, reaction, gene, direction 列的 DataFrame
    
    返回:
        解析后的 DataFrame
    """
    df = reaction_table.copy()
    
    # 标准化列名
    col_mapping = {
        'standard_metName': 'metabolite',
        'HMDB_ID': 'hmdb_id',
        'Reactions': 'reaction',
        'Gene_name': 'gene',
        'Direction': 'direction'
    }
    df = df.rename(columns={k: v for k, v in col_mapping.items() if k in df.columns})
    
    # 确保必需列存在
    required_cols = ['metabolite', 'gene', 'direction']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Reaction table missing required column: {col}")
    
    # 处理 hmdb_id
    if 'hmdb_id' not in df.columns:
        df['hmdb_id'] = np.nan
    
    # 处理 reaction 列
    if 'reaction' not in df.columns:
        df['reaction'] = 'unknown'
    
    return df


def _get_reaction_gene_sets(reaction_table: pd.DataFrame) -> pd.DataFrame:
    """
    为每个反应提取基因集合，处理多基因情况
    
    参数:
        reaction_table: 解析后的反应表
    
    返回:
        每个反应对应的基因集合
    """
    df = reaction_table.copy()
    
    # 创建反应唯一标识
    df['_reaction_id'] = df.apply(
        lambda row: f"{row.get('metabolite', '')}|{row.get('hmdb_id', '')}|{row.get('reaction', '')}|{row['direction']}",
        axis=1
    )
    
    # 解析基因列（可能包含多个基因，用分号、逗号等分隔）
    def parse_genes(gene_str):
        if pd.isna(gene_str) or gene_str == '':
            return []
        # 尝试多种分隔符
        for sep in [';', ',', '|']:
            if sep in gene_str:
                genes = [g.strip() for g in gene_str.split(sep)]
                # 去除可能的 [Enzyme] 等后缀
                genes = [g.split('[')[0].strip() for g in genes]
                return [g for g in genes if g]
        # 单个基因的情况
        gene = gene_str.split('[')[0].strip()
        return [gene] if gene else []
    
    df['_genes'] = df['gene'].apply(parse_genes)
    
    # 合并同一反应的基因。考虑改写，酶文件不应出现完全一样的行
    reaction_genes = {}
    for _, row in df.iterrows():
        rid = row['_reaction_id']
        if rid not in reaction_genes:
            reaction_genes[rid] = {
                'metabolite': row['metabolite'],
                'hmdb_id': row.get('hmdb_id', np.nan),
                'reaction': row.get('reaction', 'unknown'),
                'direction': row['direction'],
                'genes': set()
            }
        reaction_genes[rid]['genes'].update(row['_genes'])
    
    # 转换为 DataFrame
    result = pd.DataFrame([
        {
            'metabolite': r['metabolite'],
            'hmdb_id': r['hmdb_id'],
            'reaction': r['reaction'],
            'direction': r['direction'],
            'genes': list(r['genes'])
        }
        for r in reaction_genes.values()
    ])
    
    return result


def _compute_reaction_scores(
    pseudobulk: pd.DataFrame,
    reaction_genes: pd.DataFrame
) -> pd.DataFrame:
    """
    计算每个反应在每个细胞类型中的得分
    
    参数:
        pseudobulk: 细胞类型 × 基因 的表达矩阵
        reaction_genes: 反应-基因映射表
    
    返回:
        反应 × 细胞类型 的得分矩阵
    """
    scores = []
    reaction_ids = []
    
    for _, row in reaction_genes.iterrows():
        # 创建反应标识
        rid = f"{row['metabolite']}|{row.get('hmdb_id', '')}|{row['reaction']}|{row['direction']}"
        reaction_ids.append(rid)
        
        # 找到在数据集中存在的基因
        genes = [g for g in row['genes'] if g in pseudobulk.columns]
        
        if not genes:
            # 没有基因，分数为 0
            scores.append(pd.Series(0.0, index=pseudobulk.index))
        else:
            # 取基因表达的最大值
            gene_expr = pseudobulk[genes]
            max_expr = gene_expr.max(axis=1)
            scores.append(max_expr)
    
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
    # 为每个反应添加方向信息
    reaction_dir = {}
    metabolite_map = {}
    for _, row in reaction_genes.iterrows():
        rid = f"{row['metabolite']}|{row.get('hmdb_id', '')}|{row['reaction']}|{row['direction']}"
        reaction_dir[rid] = row['direction']
        metabolite_map[rid] = (row['metabolite'], row.get('hmdb_id', np.nan))
    
    # 收集所有代谢物
    metabolites = set()
    for met, hmdb in metabolite_map.values():
        metabolites.add((met, hmdb))
    
    # 初始化矩阵
    cell_types = reaction_scores.columns
    P = pd.DataFrame(0.0, index=[f"{m}|{h}" for m, h in metabolites], columns=cell_types)
    C = pd.DataFrame(0.0, index=[f"{m}|{h}" for m, h in metabolites], columns=cell_types)
    E = pd.DataFrame(0.0, index=[f"{m}|{h}" for m, h in metabolites], columns=cell_types)
    
    # 填充矩阵
    for rid in reaction_scores.index:
        direction = reaction_dir[rid]
        met, hmdb = metabolite_map[rid]
        met_idx = f"{met}|{hmdb}"
        
        if direction == 'product':
            P.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == 'substrate':
            C.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == 'transporter':
            E.loc[met_idx] += reaction_scores.loc[rid]
    
    # 重新设置索引，分离 metabolite 和 hmdb_id
    def split_index(idx):
        parts = idx.split('|')
        met = parts[0]
        hmdb = parts[1] if len(parts) > 1 else np.nan
        return pd.Series({'metabolite': met, 'hmdb_id': hmdb})
    
    for mat in [P, C, E]:
        idx_df = mat.index.to_series().apply(split_index)
        mat.index = pd.MultiIndex.from_frame(idx_df, names=['metabolite', 'hmdb_id'])
    
    return {'P': P, 'C': C, 'E': E}


def robust_minmax(
    x: np.ndarray,
    lower: float = 5,
    upper: float = 95,
    eps: float = 1e-8
) -> np.ndarray:
    """
    Robust min-max 标准化
    
    参数:
        x: 输入数组
        lower: 下限百分位数
        upper: 上限百分位数
        eps: 防止除零的小值
    
    返回:
        标准化后的数组
    """
    x = np.asarray(x, dtype=float)
    
    # 保存 NaN 的位置
    nan_mask = np.isnan(x)
    
    # 处理全 NaN 的情况
    if np.all(nan_mask):
        return np.zeros_like(x)
    
    # 计算百分位数（忽略 NaN）
    lo = np.nanpercentile(x, lower)
    hi = np.nanpercentile(x, upper)
    
    # 如果上下限接近，返回全 0
    if np.isclose(hi, lo, rtol=1e-5):
        result = np.zeros_like(x)
    else:
        # 截断并标准化
        x_clip = np.clip(x, lo, hi)
        result = (x_clip - lo) / (hi - lo + eps)
    
    # 将原来的 NaN 位置填充为 0
    result[nan_mask] = 0.0
    
    return result


def _normalize_PCE(
    P: pd.DataFrame,
    C: pd.DataFrame,
    E: pd.DataFrame,
    reaction_genes: pd.DataFrame,
    lower: float = 5,
    upper: float = 95,
    missing_C_norm: float = 0.41,
    missing_E_norm: float = 0.75
) -> Dict[str, pd.DataFrame]:
    """
    对 P, C, E 矩阵进行标准化
    
    参数:
        P, C, E: 原始矩阵
        reaction_genes: 反应-基因映射表（用于判断是否存在 substrate/transporter reaction）
        lower, upper: robust minmax 的百分位数
        missing_C_norm: 缺失 C 时的默认值
        missing_E_norm: 缺失 E 时的默认值
    
    返回:
        包含标准化后矩阵的字典
    """
    P_norm = P.copy()
    
    # 初始化 C_norm 和 E_norm，形状与 P_norm 相同，先填充默认值
    C_norm = pd.DataFrame(missing_C_norm, index=P_norm.index, columns=P_norm.columns)
    E_norm = pd.DataFrame(missing_E_norm, index=P_norm.index, columns=P_norm.columns)
    
    # 对每个代谢物进行标准化
    for met_idx in P_norm.index:
        met, hmdb = met_idx
        
        # P 标准化
        p_vals = P.loc[met_idx].values
        P_norm.loc[met_idx] = robust_minmax(p_vals, lower=lower, upper=upper)
        
        # 检查是否有 substrate reaction
        has_substrate = len(reaction_genes[
            (reaction_genes['metabolite'] == met) & 
            (reaction_genes.get('hmdb_id') == hmdb) & 
            (reaction_genes['direction'] == 'substrate')
        ]) > 0
        
        # C 标准化 - 只有当有 substrate reaction 时才计算
        if has_substrate and met_idx in C.index:
            c_vals = C.loc[met_idx].values
            C_norm.loc[met_idx] = robust_minmax(c_vals, lower=lower, upper=upper)
        
        # 检查是否有 transporter reaction
        has_transporter = len(reaction_genes[
            (reaction_genes['metabolite'] == met) & 
            (reaction_genes.get('hmdb_id') == hmdb) & 
            (reaction_genes['direction'] == 'transporter')
        ]) > 0
        
        # E 标准化 - 只有当有 transporter reaction 时才计算
        if has_transporter and met_idx in E.index:
            e_vals = E.loc[met_idx].values
            E_norm.loc[met_idx] = robust_minmax(e_vals, lower=lower, upper=upper)
    
    return {'P_norm': P_norm, 'C_norm': C_norm, 'E_norm': E_norm}


def compute_metabolite_availability(
    adata,
    reaction_table: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    lower: float = 5,
    upper: float = 95,
    eps: float = 0.05,
    beta: float = 0.5,
    missing_C_norm: float = 0.41,
    missing_E_norm: float = 0.75,
    min_cells: int = 1,
    return_intermediates: bool = True
) -> Dict[str, Any]:
    """
    计算代谢物 availability
    
    参数:
        adata: AnnData 对象
        reaction_table: 反应表
        celltype_col: 细胞类型列名
        layer: 使用的层
        lower, upper: robust minmax 的百分位数
        eps: 公式中的小值
        beta: 公式中的指数
        missing_C_norm: 缺失 C 时的默认值
        missing_E_norm: 缺失 E 时的默认值
        min_cells: 每个细胞类型的最小细胞数
        return_intermediates: 是否返回中间结果
    
    返回:
        包含 availability 和中间结果的字典
    """
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
    
    # 6. 过滤掉没有 product reaction 的代谢物
    valid_mets = P.index[P.sum(axis=1) > 0]
    P = P.loc[valid_mets]
    C = C.loc[valid_mets.intersection(C.index)] if not C.empty else C.head(0)
    E = E.loc[valid_mets.intersection(E.index)] if not E.empty else E.head(0)
    
    # 7. 标准化
    normalized = _normalize_PCE(
        P, C, E, reaction_genes,
        lower=lower,
        upper=upper,
        missing_C_norm=missing_C_norm,
        missing_E_norm=missing_E_norm
    )
    P_norm, C_norm, E_norm = normalized['P_norm'], normalized['C_norm'], normalized['E_norm']
    
    # 8. 计算最终 availability
    availability = pd.DataFrame(index=P_norm.index, columns=P_norm.columns)
    for met_idx in P_norm.index:
        p = P_norm.loc[met_idx].values
        c = C_norm.loc[met_idx].values
        e = E_norm.loc[met_idx].values
        
        avail = (p + eps) * np.power((1 - c + eps), beta) * (e + eps)
        availability.loc[met_idx] = avail
    
    # 9. 准备元数据
    metadata = pd.DataFrame(index=P_norm.index)
    metadata['has_product'] = True
    metadata['has_substrate'] = [len(reaction_genes[
        (reaction_genes['metabolite'] == met[0]) & 
        (reaction_genes['hmdb_id'] == met[1]) & 
        (reaction_genes['direction'] == 'substrate')
    ]) > 0 for met in P_norm.index]
    metadata['has_transporter'] = [len(reaction_genes[
        (reaction_genes['metabolite'] == met[0]) & 
        (reaction_genes['hmdb_id'] == met[1]) & 
        (reaction_genes['direction'] == 'transporter')
    ]) > 0 for met in P_norm.index]
    metadata['n_product_reactions'] = [len(reaction_genes[
        (reaction_genes['metabolite'] == met[0]) & 
        (reaction_genes['hmdb_id'] == met[1]) & 
        (reaction_genes['direction'] == 'product')
    ]) for met in P_norm.index]
    metadata['n_substrate_reactions'] = [len(reaction_genes[
        (reaction_genes['metabolite'] == met[0]) & 
        (reaction_genes['hmdb_id'] == met[1]) & 
        (reaction_genes['direction'] == 'substrate')
    ]) for met in P_norm.index]
    metadata['n_transporter_reactions'] = [len(reaction_genes[
        (reaction_genes['metabolite'] == met[0]) & 
        (reaction_genes['hmdb_id'] == met[1]) & 
        (reaction_genes['direction'] == 'transporter')
    ]) for met in P_norm.index]
    
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
