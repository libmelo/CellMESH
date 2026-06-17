"""
CELL MESH 核心算法模块
完全基于 metabolite availability 算法计算代谢物通信事件
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any

import numpy as np
import pandas as pd

# 导入集中配置
from .config import (
    MIN_EXPR_FRAC,
    METABOLITE_AVAILABILITY_DEFAULTS
)

from .database import load_cell_mesh_database, validate_priors
from .score import (
    compute_metabolite_availability,
    compute_sensor_scores
)


EVENT_COLUMNS = [
    "sender",
    "receiver",
    "metabolite",
    "hmdb_id",
    "sensor_gene",
    "sensor_type",
    "metabolite_availability",
    "sensor_score",
    "sensor_expr_frac",
    "cell_mesh_score",
]


@dataclass
class CellMeshResult:
    """CELL MESH 运行结果容器"""
    events: pd.DataFrame
    sender_scores: pd.DataFrame
    receiver_scores: pd.DataFrame
    parameters: dict
    availability_results: Optional[Dict[str, Any]] = None

    def to_csv(self, prefix: str) -> None:
        """
        将结果保存为 CSV 文件

        参数:
            prefix: 保存路径前缀
        """
        self.events.to_csv(f"{prefix}.events.csv", index=False)
        self.sender_scores.to_csv(f"{prefix}.sender_scores.csv", index=False)
        self.receiver_scores.to_csv(f"{prefix}.receiver_scores.csv", index=False)


def _bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg 错误发现率校正

    参数:
        pvalues: p 值数组

    返回:
        校正后的 FDR 值数组
    """
    p = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    if n == 0:
        return pvalues
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def _same_hmdb(left: object, right: object) -> bool:
    if pd.isna(left) or pd.isna(right):
        return False
    return str(left) == str(right)


def _compute_availability_scores(
    adata,
    enzyme_prior: pd.DataFrame,
    sensor_prior: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_expr_frac: float = MIN_EXPR_FRAC,
    lower: float = METABOLITE_AVAILABILITY_DEFAULTS["lower"],
    upper: float = METABOLITE_AVAILABILITY_DEFAULTS["upper"],
    eps: float = METABOLITE_AVAILABILITY_DEFAULTS["eps"],
    beta: float = METABOLITE_AVAILABILITY_DEFAULTS["beta"],
    missing_C_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_C_norm"],
    missing_E_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_E_norm"],
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
) -> tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    """
    基于 metabolite availability 计算 sender 和 receiver 得分

    参数:
        adata: AnnData 对象
        enzyme_prior: 经过 validate_priors 验证后的酶-代谢物先验表
        sensor_prior: 代谢物-传感器先验表
        celltype_col: 细胞类型列名
        layer: 表达层
        min_expr_frac: 最小表达比例阈值
        **kwargs: availability 计算参数

    返回:
        (sender_scores, receiver_scores, availability_results) 元组
    """
    # 计算代谢物 availability（保持不变）
    avail_results = compute_metabolite_availability(
        adata,
        enzyme_prior,
        celltype_col=celltype_col,
        layer=layer,
        lower=lower,
        upper=upper,
        eps=eps,
        beta=beta,
        missing_C_norm=missing_C_norm,
        missing_E_norm=missing_E_norm,
        min_cells=min_cells,
        return_intermediates=True
    )
    
    availability = avail_results['availability']
    
    # 如果没有可用代谢物，返回空结果
    if availability.empty:
        return pd.DataFrame(), pd.DataFrame(), avail_results
    
    # sender_scores is the metabolite availability matrix kept as a separate result.
    sender_scores = availability.copy()
    
    # 计算新的 sensor scores
    receiver_scores = compute_sensor_scores(
        adata,
        sensor_prior,
        celltype_col=celltype_col,
        layer=layer,
        lower=lower,
        upper=upper,
        min_expr_frac=min_expr_frac,
        min_cells=min_cells,
        pseudobulk=avail_results.get("pseudobulk"),
        expr_frac=avail_results.get("expr_frac"),
    )
    
    return sender_scores, receiver_scores, avail_results


def _make_cell_mesh_events(
    sender_scores: pd.DataFrame,
    receiver_scores: pd.DataFrame,
    allow_self: bool
) -> pd.DataFrame:
    """
    构建 CELL MESH 通信事件

    参数:
        sender_scores: 发送方得分(来自 availability)
        receiver_scores: 接收方得分
        allow_self: 是否允许自分泌

    返回:
        通信事件 DataFrame
    """
    # 如果任意一方得分是空,返回空结果
    if sender_scores.empty or receiver_scores.empty:
        return pd.DataFrame(columns=EVENT_COLUMNS)

    rows = []
    for _, rr in receiver_scores.iterrows():
        metabolite = rr["metabolite"]
        hmdb_id = rr.get("hmdb_id", np.nan)
        if isinstance(sender_scores.index, pd.MultiIndex):
            sender_matches = [
                idx for idx in sender_scores.index
                if idx[0] == metabolite and _same_hmdb(idx[1], hmdb_id)
            ]
        else:
            sender_matches = [metabolite] if metabolite in sender_scores.index else []

        if not sender_matches:
            continue
        for sender_idx in sender_matches:
            for sender, availability_value in sender_scores.loc[sender_idx].items():
                receiver = rr["receiver"]
                if (not allow_self) and sender == receiver:
                    continue

                availability = float(availability_value)
                sensor_score = float(rr["sensor_score"])
                cell_mesh_score = float(np.sqrt(availability * sensor_score))

                rows.append(
                    {
                        "sender": sender,
                        "receiver": receiver,
                        "metabolite": metabolite,
                        "hmdb_id": hmdb_id,
                        "sensor_gene": rr["sensor_gene"],
                        "sensor_type": rr["sensor_type"],
                        "metabolite_availability": availability,
                        "sensor_score": sensor_score,
                        "sensor_expr_frac": float(rr["sensor_expr_frac"]),
                        "cell_mesh_score": cell_mesh_score,
                    }
                )

    if not rows:
        return pd.DataFrame(columns=EVENT_COLUMNS)
    return pd.DataFrame(rows, columns=EVENT_COLUMNS).sort_values("cell_mesh_score", ascending=False).reset_index(drop=True)


def _permute_labels(
    labels: pd.Series,
    sample_labels: Optional[pd.Series],
    rng: np.random.Generator
) -> pd.Series:
    """
    置换标签用于置换检验

    参数:
        labels: 原始标签
        sample_labels: 样本标签(如果提供,置换将在样本内进行)
        rng: 随机数生成器

    返回:
        置换后的标签
    """
    vals = labels.astype(str).copy()
    out = vals.copy()
    if sample_labels is None:
        out[:] = rng.permutation(vals.values)
    else:
        sample_text = sample_labels.astype(str)
        for sample in sample_text.unique():
            idx = np.flatnonzero(sample_text.values == sample)
            out.iloc[idx] = rng.permutation(vals.iloc[idx].values)
    return out


def _empirical_pvalues_by_sensor_type(
    obs_events: pd.DataFrame,
    adata,
    cell_type_key: str,
    sample_key: Optional[str],
    layer: Optional[str],
    enzyme_prior: pd.DataFrame,
    sensor_prior: pd.DataFrame,
    n_perms: int,
    random_state: int,
    min_expr_frac: float,
    allow_self: bool,
    availability_kwargs: dict,
) -> pd.DataFrame:
    """
    计算经验 p 值(置换检验)，按 sensor type 分别计算 null 分布

    参数:
        obs_events: 观察到的事件
        adata: AnnData 对象
        cell_type_key: 细胞类型列名
        sample_key: 样本列名
        layer: 表达层
        enzyme_prior: 酶先验
        sensor_prior: 传感器先验
        n_perms: 置换次数
        random_state: 随机种子
        **kwargs: 其他参数

    返回:
        带有 p 值和 FDR 的事件 DataFrame（FDR 按 sensor type 分别校正）
    """
    if n_perms <= 0 or obs_events.empty:
        obs_events["perm_pvalue"] = np.nan
        obs_events["fdr"] = np.nan
        return obs_events

    key_cols = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "sensor_type"]
    obs_keys = obs_events[key_cols].astype(str).agg("|".join, axis=1)
    ge_counts = pd.Series(0, index=obs_keys.values, dtype=int)
    obs_score = pd.Series(obs_events["cell_mesh_score"].values, index=obs_keys.values)
    obs_sensor_types = pd.Series(obs_events["sensor_type"].values, index=obs_keys.values)

    rng = np.random.default_rng(random_state)
    original = adata.obs[cell_type_key].copy()
    sample_labels = adata.obs[sample_key].copy() if sample_key is not None else None
    perm_key = "_cell_mesh_perm_label"

    try:
        for perm_idx in range(n_perms):
            # 1. 打乱细胞类型标签
            adata.obs[perm_key] = _permute_labels(original, sample_labels, rng).values

            # 2. 重新计算 availability 和得分
            sender_perm, receiver_perm, _ = _compute_availability_scores(
                adata,
                enzyme_prior,
                sensor_prior,
                celltype_col=perm_key,
                layer=layer,
                min_expr_frac=min_expr_frac,
                **availability_kwargs
            )

            # 3. 构建置换事件
            events_perm = _make_cell_mesh_events(sender_perm, receiver_perm, allow_self=allow_self)
            if events_perm.empty:
                continue

            # 4. 比较得分计数（按 sensor type 分别进行比较）
            perm_scores = events_perm.assign(_key=events_perm[key_cols].astype(str).agg("|".join, axis=1)).set_index("_key")["cell_mesh_score"]
            perm_sensor_types = events_perm.assign(_key=events_perm[key_cols].astype(str).agg("|".join, axis=1)).set_index("_key")["sensor_type"]
            
            common = obs_score.index.intersection(perm_scores.index)
            for key in common:
                # 只在同一 sensor type 内比较
                if obs_sensor_types.loc[key] == perm_sensor_types.loc[key]:
                    if perm_scores.loc[key] >= obs_score.loc[key]:
                        ge_counts.loc[key] += 1

    finally:
        if perm_key in adata.obs:
            del adata.obs[perm_key]

    # 计算 p 值
    p = (ge_counts.loc[obs_keys.values].values + 1) / (n_perms + 1)
    out = obs_events.copy()
    out["perm_pvalue"] = p

    # 按 sensor type 分别计算 FDR
    out["fdr"] = np.nan
    for sensor_type in out["sensor_type"].unique():
        mask = out["sensor_type"] == sensor_type
        type_pvalues = out.loc[mask, "perm_pvalue"].values
        out.loc[mask, "fdr"] = _bh_fdr(type_pvalues)

    return out


def _confidence_tier(row: pd.Series) -> str:
    """
    确定事件的置信等级

    参数:
        row: 事件行

    返回:
        置信等级字符串
    """
    if pd.isna(row.get("fdr", np.nan)):
        if row["cell_mesh_score"] >= 0.5 and row.get("sensor_expr_frac", 0) >= 0.1:
            return "Tier2_no_permutation"
        return "Tier3_exploratory"

    if row["fdr"] <= 0.05 and row["cell_mesh_score"] >= 0.5 and row.get("sensor_expr_frac", 0) >= 0.1:
        return "Tier1_high"
    if row["fdr"] <= 0.1 and row["cell_mesh_score"] >= 0.25:
        return "Tier2_medium"

    return "Tier3_exploratory"


def run_cell_mesh(
    adata,
    enzyme_metabolite: Optional[pd.DataFrame] = None,
    metabolite_sensor: Optional[pd.DataFrame] = None,
    cell_type_key: str = "cell_type",
    sample_key: Optional[str] = None,
    layer: Optional[str] = None,
    min_expr_frac: float = MIN_EXPR_FRAC,
    allow_self: bool = True,
    n_perms: int = 0,
    random_state: int = 0,
    lower: float = METABOLITE_AVAILABILITY_DEFAULTS["lower"],
    upper: float = METABOLITE_AVAILABILITY_DEFAULTS["upper"],
    eps: float = METABOLITE_AVAILABILITY_DEFAULTS["eps"],
    beta: float = METABOLITE_AVAILABILITY_DEFAULTS["beta"],
    missing_C_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_C_norm"],
    missing_E_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_E_norm"],
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
) -> CellMeshResult:
    """
    运行 CELL MESH: Metabolite-mediated Event Scoring with Sensor Hierarchies.
    完全基于 metabolite availability 算法计算代谢物通信事件。

    参数:
        adata: AnnData 对象,包含单细胞表达数据
        enzyme_metabolite: 酶-代谢物关系先验表,默认使用内置数据库
            必需列:metabolite, gene, role
            可选列:hmdb_id, reaction, weight, evidence_level, source
            role 取值:production (产生)、degradation (降解)、export (外排)
        metabolite_sensor: 代谢物-传感器关系先验表,默认使用内置数据库
        cell_type_key: 细胞类型列名,默认为 "cell_type"
        sample_key: 样本列名,用于置换检验时的样本内置换
        layer: 使用的表达层,None 表示使用 adata.X
        min_expr_frac: 最小表达比例阈值,低于该值的基因表达视为不表达
        allow_self: 是否允许自分泌通信
        n_perms: 置换检验次数,0 表示不进行置换检验
        random_state: 随机种子
        lower: availability 标准化的下限百分位数,默认 5
        upper: availability 标准化的上限百分位数,默认 95
        eps: availability 计算中的小常数,避免除以零,默认 0.05
        beta: 消耗项的指数权重,默认 0.5
        missing_C_norm: 当代谢物没有消耗证据时的默认 C_norm 值,默认 0.2
        missing_E_norm: 当代谢物没有外排证据时的默认 E_norm 值,默认 0.5
        min_cells: 每个细胞类型的最小细胞数,低于该值的细胞类型会被过滤

    返回:
        CellMeshResult 对象,包含所有计算结果

    内部逻辑说明:
        1. enzyme_metabolite 作为标准 enzyme prior 直接传入 availability 计算，
           availability 内部负责 role 到 direction 的映射:
           - production → product → 进入 P (产生) 矩阵
           - degradation → substrate → 进入 C (消耗) 矩阵
           - export → exporter → 进入 E (外排) 矩阵
        2. metabolite_availability 完全来自 metabolite availability 计算:
           availability = P_norm * ((1 - C_norm) ** beta) * (0.8 + 0.2 * E_norm)
           结果范围在 [0, 1] 之间,值越高代表该细胞类型释放该代谢物的能力越强
        3. sensor_score 基于 robust min-max 标准化的 sensor 基因表达
        4. cell_mesh_score = sqrt(metabolite_availability * sensor_score)
    """
    # 验证输入
    if sample_key is not None and sample_key not in adata.obs:
        raise KeyError(f"{sample_key!r} not found in adata.obs")

    # 加载默认数据库
    if enzyme_metabolite is None or metabolite_sensor is None:
        default_enzyme, default_sensor = load_cell_mesh_database()
        enzyme_metabolite = default_enzyme if enzyme_metabolite is None else enzyme_metabolite
        metabolite_sensor = default_sensor if metabolite_sensor is None else metabolite_sensor

    # 验证先验
    enzyme_prior, sensor_prior = validate_priors(enzyme_metabolite, metabolite_sensor, adata.var_names)

    if enzyme_prior.empty:
        raise ValueError("No enzyme prior genes found in adata.var_names")
    if sensor_prior.empty:
        raise ValueError("No sensor genes found in adata.var_names")

    # 计算 availability 和得分
    availability_kwargs = {
        'lower': lower,
        'upper': upper,
        'eps': eps,
        'beta': beta,
        'missing_C_norm': missing_C_norm,
        'missing_E_norm': missing_E_norm,
        'min_cells': min_cells,
    }

    sender_scores, receiver_scores, availability_results = _compute_availability_scores(
        adata,
        enzyme_prior,
        sensor_prior,
        celltype_col=cell_type_key,
        layer=layer,
        min_expr_frac=min_expr_frac,
        **availability_kwargs
    )

    # 构建事件
    events = _make_cell_mesh_events(sender_scores, receiver_scores, allow_self=allow_self)

    # 计算显著性（按 sensor type 分别计算）
    events = _empirical_pvalues_by_sensor_type(
        events,
        adata=adata,
        cell_type_key=cell_type_key,
        sample_key=sample_key,
        layer=layer,
        enzyme_prior=enzyme_prior,
        sensor_prior=sensor_prior,
        n_perms=n_perms,
        random_state=random_state,
        min_expr_frac=min_expr_frac,
        allow_self=allow_self,
        availability_kwargs=availability_kwargs
    )

    # 计算置信等级
    if not events.empty:
        events["confidence_tier"] = events.apply(_confidence_tier, axis=1)
        events = events.sort_values(["fdr", "cell_mesh_score"], ascending=[True, False], na_position="last").reset_index(drop=True)
    elif "confidence_tier" not in events.columns:
        events["confidence_tier"] = pd.Series(dtype=object)

    # 整理参数
    parameters = {
        "method": "CELL MESH",
        "acronym": "Metabolite-mediated Event Scoring with Sensor Hierarchies",
        "algorithm": "metabolite availability + robust min-max sensor scoring",
        "cell_type_key": cell_type_key,
        "sample_key": sample_key,
        "layer": layer,
        "min_expr_frac": min_expr_frac,
        "allow_self": allow_self,
        "n_perms": n_perms,
        "random_state": random_state,
        **availability_kwargs
    }

    return CellMeshResult(
        events=events, sender_scores=sender_scores, receiver_scores=receiver_scores,
        parameters=parameters, availability_results=availability_results
    )
