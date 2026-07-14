"""
CELL MESH 核心算法模块
完全基于 metabolite availability 算法计算代谢物通信事件
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any, Literal

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
from .preprocess import _eligible_celltype_counts


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
    "sender_n_cells",
    "receiver_n_cells",
    "cell_mesh_score",
]

SAMPLE_MODES = {"pooled_stratified", "sample_aware"}
EVENT_KEY_COLUMNS = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "sensor_type"]


@dataclass
class CellMeshResult:
    """CELL MESH 运行结果容器"""
    events: pd.DataFrame
    sender_scores: pd.DataFrame
    receiver_scores: pd.DataFrame
    parameters: dict
    availability_results: Optional[Dict[str, Any]] = None
    sample_validation: Optional[pd.DataFrame] = None
    sample_sender_scores: Optional[pd.DataFrame] = None
    sample_receiver_scores: Optional[pd.DataFrame] = None
    sample_events: Optional[pd.DataFrame] = None

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
    min_expr_frac: Optional[float] = MIN_EXPR_FRAC,
    eps_num: float = METABOLITE_AVAILABILITY_DEFAULTS["eps_num"],
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
        min_expr_frac: 可选 receiver 表达比例 gate

    返回:
        (sender_scores, receiver_scores, availability_results) 元组
    """
    # Compute P/C/E and the median-contrast sender score.
    avail_results = compute_metabolite_availability(
        adata,
        enzyme_prior,
        celltype_col=celltype_col,
        layer=layer,
        eps_num=eps_num,
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
        min_expr_frac=min_expr_frac,
        min_cells=min_cells,
        eps_num=eps_num,
        pseudobulk=avail_results.get("pseudobulk"),
        expr_frac=avail_results.get("expr_frac"),
        cell_counts=avail_results.get("cell_counts"),
    )
    
    return sender_scores, receiver_scores, avail_results


def _validate_sample_key(adata, sample_key: Optional[str], sample_mode: str) -> None:
    if sample_mode not in SAMPLE_MODES:
        raise ValueError("sample_mode must be one of {'pooled_stratified', 'sample_aware'}")
    if sample_mode == "sample_aware" and sample_key is None:
        raise ValueError("sample_aware mode requires a valid sample_key")
    if sample_key is None:
        return
    if sample_key not in adata.obs:
        raise KeyError(f"{sample_key!r} not found in adata.obs")
    values = adata.obs[sample_key]
    text = values.astype(str).str.strip()
    if values.isna().any() or (text == "").any():
        raise ValueError(f"{sample_key!r} must not contain NA or empty strings")


def _sample_validation_table(adata, sample_key: str, cell_type_key: str, min_cells: int) -> pd.DataFrame:
    obs = pd.DataFrame(
        {
            "sample": adata.obs[sample_key].astype(str),
            "cell_type": adata.obs[cell_type_key].astype(str),
        }
    )
    validation = (
        obs.groupby(["sample", "cell_type"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    validation["eligible_in_sample"] = validation["n_cells"] >= min_cells
    eligible_counts = validation[validation["eligible_in_sample"]].groupby("sample", observed=True).size()
    valid_samples = validation[validation["eligible_in_sample"]].groupby("cell_type", observed=True).size()
    validation["n_eligible_celltypes_in_sample"] = (
        validation["sample"].map(eligible_counts).fillna(0).astype(int)
    )
    validation["n_valid_samples_for_celltype"] = (
        validation["cell_type"].map(valid_samples).fillna(0).astype(int)
    )
    return validation


def _sample_stat_iqr(values: pd.Series) -> float:
    x = pd.to_numeric(values, errors="coerce").dropna()
    if len(x) < 2:
        return np.nan
    return float(x.quantile(0.75) - x.quantile(0.25))


def _sample_positive_prevalence(values: pd.Series) -> float:
    x = pd.to_numeric(values, errors="coerce").dropna()
    if len(x) == 0:
        return np.nan
    return float((x > 0).sum() / len(x))


def _assign_fdr_columns(
    events: pd.DataFrame,
    pvalue_col: str = "perm_pvalue",
) -> pd.DataFrame:
    out = events.copy()
    out["fdr_global"] = np.nan
    valid = out[pvalue_col].notna()
    if valid.any():
        out.loc[valid, "fdr_global"] = _bh_fdr(out.loc[valid, pvalue_col].to_numpy(dtype=float))

    out["fdr_sensor_type"] = np.nan
    for sensor_type in out.loc[valid, "sensor_type"].unique():
        mask = valid & (out["sensor_type"] == sensor_type)
        out.loc[mask, "fdr_sensor_type"] = _bh_fdr(out.loc[mask, pvalue_col].to_numpy(dtype=float))
    out["fdr"] = out["fdr_sensor_type"]
    return out


def _sample_aware_empirical_pvalues(
    obs_events: pd.DataFrame,
    adata,
    cell_type_key: str,
    sample_key: str,
    layer: Optional[str],
    enzyme_prior: pd.DataFrame,
    sensor_prior: pd.DataFrame,
    n_perms: int,
    random_state: int,
    min_expr_frac: Optional[float],
    allow_self: bool,
    availability_kwargs: dict,
) -> pd.DataFrame:
    out = obs_events.copy()
    out["permutation_mode"] = "within_sample_label_shuffle"

    null_index = pd.MultiIndex.from_frame(out[EVENT_KEY_COLUMNS].astype(str))
    null_df = pd.DataFrame(index=null_index)

    if out.empty or n_perms <= 0:
        out["perm_pvalue"] = np.nan
        out = _assign_fdr_columns(out)
        out.attrs["sample_aware_null_scores"] = null_df
        return out

    score_col = "event_score_median" if "event_score_median" in out.columns else "cell_mesh_score"
    obs_scores = pd.Series(out[score_col].to_numpy(dtype=float), index=null_index)

    min_cells = availability_kwargs.get("min_cells", METABOLITE_AVAILABILITY_DEFAULTS["min_cells"])
    validation = _sample_validation_table(adata, sample_key, cell_type_key, min_cells)
    eligible_units = validation.loc[
        validation["eligible_in_sample"], ["sample", "cell_type"]
    ]
    eligible_pairs = set(map(tuple, eligible_units.astype(str).to_numpy()))
    eligible_mask = [
        (str(sample), str(cell_type)) in eligible_pairs
        for sample, cell_type in zip(adata.obs[sample_key], adata.obs[cell_type_key])
    ]
    adata_perm_base = adata[eligible_mask, :].copy()

    rng = np.random.default_rng(random_state)
    original = adata_perm_base.obs[cell_type_key].copy()
    sample_labels = adata_perm_base.obs[sample_key].copy()
    perm_key = "_cell_mesh_perm_label"
    null_columns = []

    try:
        for perm_idx in range(n_perms):
            adata_perm_base.obs[perm_key] = _permute_labels(original, sample_labels, rng).values
            _, _, events_perm, _ = _compute_sample_aware_scores(
                adata_perm_base,
                enzyme_prior,
                sensor_prior,
                cell_type_key=perm_key,
                sample_key=sample_key,
                layer=layer,
                min_expr_frac=min_expr_frac,
                allow_self=allow_self,
                availability_kwargs=availability_kwargs,
            )

            if events_perm.empty:
                perm_scores = pd.Series(0.0, index=null_index)
            else:
                perm_score_col = (
                    "event_score_median"
                    if "event_score_median" in events_perm.columns
                    else "cell_mesh_score"
                )
                perm_index = pd.MultiIndex.from_frame(events_perm[EVENT_KEY_COLUMNS].astype(str))
                perm_scores = pd.Series(
                    events_perm[perm_score_col].to_numpy(dtype=float),
                    index=perm_index,
                ).reindex(null_index).fillna(0.0)
            null_columns.append(perm_scores.to_numpy(dtype=float))
    finally:
        if perm_key in adata_perm_base.obs:
            del adata_perm_base.obs[perm_key]

    if null_columns:
        null_df = pd.DataFrame(
            np.column_stack(null_columns),
            index=null_index,
            columns=[f"perm_{i}" for i in range(n_perms)],
        )
    else:
        null_df = pd.DataFrame(index=null_index)

    pvalues = pd.Series(np.nan, index=null_index, dtype=float)
    valid_obs = obs_scores.notna()
    if valid_obs.any():
        ge_counts = null_df.ge(obs_scores, axis=0).sum(axis=1)
        pvalues.loc[valid_obs] = (ge_counts.loc[valid_obs].to_numpy(dtype=float) + 1.0) / (n_perms + 1.0)

    out["perm_pvalue"] = pvalues.to_numpy()
    out = _assign_fdr_columns(out)
    out["permutation_mode"] = "within_sample_label_shuffle"
    out.attrs["sample_aware_null_scores"] = null_df
    return out


def _compute_sample_aware_scores(
    adata,
    enzyme_prior: pd.DataFrame,
    sensor_prior: pd.DataFrame,
    cell_type_key: str,
    sample_key: str,
    layer: Optional[str],
    min_expr_frac: Optional[float],
    allow_self: bool,
    availability_kwargs: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    min_cells = availability_kwargs.get("min_cells", METABOLITE_AVAILABILITY_DEFAULTS["min_cells"])
    validation = _sample_validation_table(adata, sample_key, cell_type_key, min_cells)
    all_celltypes = pd.Index(sorted(validation["cell_type"].unique()))

    sample_sender = []
    sample_receiver = []
    sample_events = []
    availability_by_sample: Dict[str, Dict[str, Any]] = {}
    valid_sample_names = []

    for sample in sorted(validation["sample"].unique()):
        sample_validation = validation[validation["sample"] == sample]
        eligible_celltypes = sample_validation.loc[
            sample_validation["eligible_in_sample"], "cell_type"
        ].tolist()
        if not eligible_celltypes:
            continue
        valid_sample_names.append(sample)

        obs = adata.obs
        mask = (
            (obs[sample_key].astype(str).values == sample)
            & obs[cell_type_key].astype(str).isin(eligible_celltypes).values
        )
        adata_sample = adata[mask, :].copy()
        sender_scores, receiver_scores, availability_results = _compute_availability_scores(
            adata_sample,
            enzyme_prior,
            sensor_prior,
            celltype_col=cell_type_key,
            layer=layer,
            min_expr_frac=min_expr_frac,
            **availability_kwargs,
        )
        availability_by_sample[sample] = availability_results

        if not sender_scores.empty:
            sender_reindexed = sender_scores.reindex(columns=all_celltypes, fill_value=np.nan)
            sender_reindexed.index = pd.MultiIndex.from_tuples(
                [(sample, met, hmdb) for met, hmdb in sender_reindexed.index],
                names=["sample", "metabolite", "hmdb_id"],
            )
            sample_sender.append(sender_reindexed)

        if not receiver_scores.empty:
            receiver_scores = receiver_scores.copy()
            receiver_scores.insert(0, "sample", sample)
            sample_receiver.append(receiver_scores)

        events = _make_cell_mesh_events(
            sender_scores,
            receiver_scores,
            allow_self=allow_self,
            cell_counts=availability_results.get("cell_counts"),
        )
        if not events.empty:
            events.insert(0, "sample", sample)
            sample_events.append(events)

    sample_sender_scores = (
        pd.concat(sample_sender, axis=0)
        if sample_sender
        else pd.DataFrame(columns=all_celltypes)
    )
    sender_scores = (
        sample_sender_scores.groupby(level=["metabolite", "hmdb_id"]).mean()
        if not sample_sender_scores.empty
        else pd.DataFrame()
    )

    sample_receiver_scores = (
        pd.concat(sample_receiver, ignore_index=True)
        if sample_receiver
        else pd.DataFrame()
    )
    receiver_key = ["metabolite", "hmdb_id", "sensor_gene", "sensor_type", "receiver"]
    if not sample_receiver_scores.empty:
        receiver_scores = (
            sample_receiver_scores.groupby(receiver_key, as_index=False)
            .agg(
                sensor_score=("sensor_score", "mean"),
                sensor_expr_frac=("sensor_expr_frac", "mean"),
                receiver_n_cells=("receiver_n_cells", "sum"),
                n_valid_samples=("sample", "nunique"),
            )
        )
    else:
        receiver_scores = pd.DataFrame()

    sample_events_df = (
        pd.concat(sample_events, ignore_index=True)
        if sample_events
        else pd.DataFrame(columns=["sample"] + EVENT_COLUMNS)
    )
    event_key = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "sensor_type"]
    if not sample_events_df.empty:
        event_index = sample_events_df[event_key].drop_duplicates()
        sample_event_index = pd.MultiIndex.from_frame(
            event_index.merge(pd.DataFrame({"sample": valid_sample_names}), how="cross")[
                event_key + ["sample"]
            ]
        )
        sample_events_complete = (
            sample_events_df.set_index(event_key + ["sample"])
            .reindex(sample_event_index)
            .reset_index()
        )
        sample_events_df = sample_events_complete

        grouped_scores = sample_events_df.groupby(event_key)["cell_mesh_score"]
        score_summary = grouped_scores.agg(
            event_score_median="median",
            n_samples_coobserved="count",
            n_samples_positive=lambda s: int((pd.to_numeric(s, errors="coerce").dropna() > 0).sum()),
        ).reset_index()
        score_summary["event_score_iqr"] = grouped_scores.apply(_sample_stat_iqr).values
        score_summary["event_prevalence"] = grouped_scores.apply(_sample_positive_prevalence).values

        meta_summary = (
            sample_events_df.groupby(event_key, as_index=False)
            .agg(
                metabolite_availability=("metabolite_availability", "median"),
                sensor_score=("sensor_score", "median"),
                sensor_expr_frac=("sensor_expr_frac", "median"),
                sender_n_cells=("sender_n_cells", "sum"),
                receiver_n_cells=("receiver_n_cells", "sum"),
            )
        )
        events = (
            meta_summary.merge(score_summary, on=event_key, how="inner")
            .query("n_samples_coobserved > 0")
            .assign(
                cell_mesh_score=lambda df: df["event_score_median"],
                n_valid_samples=lambda df: df["n_samples_coobserved"],
                inference_mode="sample_aware",
            )
            .sort_values("cell_mesh_score", ascending=False, na_position="last")
            .reset_index(drop=True)
        )
        sample_event_scores = sample_events_df.pivot_table(
            index=event_key,
            columns="sample",
            values="cell_mesh_score",
            aggfunc="mean",
            dropna=False,
        )
    else:
        events = pd.DataFrame(columns=EVENT_COLUMNS + ["n_valid_samples"])
        sample_event_scores = pd.DataFrame()

    availability_results = {
        "sample_validation": validation,
        "sample_sender_scores": sample_sender_scores,
        "sample_receiver_scores": sample_receiver_scores,
        "sample_events": sample_events_df,
        "sample_event_scores": sample_event_scores,
        "availability_by_sample": availability_by_sample,
    }
    return sender_scores, receiver_scores, events, availability_results


def _make_cell_mesh_events(
    sender_scores: pd.DataFrame,
    receiver_scores: pd.DataFrame,
    allow_self: bool,
    cell_counts: Optional[pd.Series] = None,
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
                        "sender_n_cells": (
                            int(cell_counts.loc[sender])
                            if cell_counts is not None and sender in cell_counts.index
                            else np.nan
                        ),
                        "receiver_n_cells": int(rr["receiver_n_cells"]),
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
    min_expr_frac: Optional[float],
    allow_self: bool,
    availability_kwargs: dict,
    sample_mode: str = "pooled_stratified",
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
    if sample_mode == "sample_aware":
        return _sample_aware_empirical_pvalues(
            obs_events,
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
            availability_kwargs=availability_kwargs,
        )

    if n_perms <= 0 or obs_events.empty:
        obs_events["perm_pvalue"] = np.nan
        obs_events["fdr"] = np.nan
        return obs_events

    key_cols = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "sensor_type"]
    obs_keys = obs_events[key_cols].astype(str).agg("|".join, axis=1)
    ge_counts = pd.Series(0, index=obs_keys.values, dtype=int)
    obs_score = pd.Series(obs_events["cell_mesh_score"].values, index=obs_keys.values)
    obs_sensor_types = pd.Series(obs_events["sensor_type"].values, index=obs_keys.values)

    min_cells = availability_kwargs.get("min_cells", METABOLITE_AVAILABILITY_DEFAULTS["min_cells"])
    if sample_mode == "sample_aware":
        validation = _sample_validation_table(adata, sample_key, cell_type_key, min_cells)
        eligible_units = validation.loc[
            validation["eligible_in_sample"], ["sample", "cell_type"]
        ]
        eligible_pairs = set(map(tuple, eligible_units.to_numpy()))
        eligible_mask = [
            (str(sample), str(cell_type)) in eligible_pairs
            for sample, cell_type in zip(adata.obs[sample_key], adata.obs[cell_type_key])
        ]
    else:
        eligible_celltypes = _eligible_celltype_counts(adata, cell_type_key, min_cells).index
        eligible_mask = adata.obs[cell_type_key].astype(str).isin(eligible_celltypes).values
    adata_perm_base = adata[eligible_mask, :].copy()

    rng = np.random.default_rng(random_state)
    original = adata_perm_base.obs[cell_type_key].copy()
    sample_labels = adata_perm_base.obs[sample_key].copy() if sample_key is not None else None
    perm_key = "_cell_mesh_perm_label"

    try:
        for perm_idx in range(n_perms):
            # 1. 打乱 eligible 细胞的细胞类型标签；低于 min_cells 的 cell type 不进入 null。
            adata_perm_base.obs[perm_key] = _permute_labels(original, sample_labels, rng).values

            # 2. 重新计算 availability 和得分
            if sample_mode == "sample_aware":
                _, _, events_perm, _ = _compute_sample_aware_scores(
                    adata_perm_base,
                    enzyme_prior,
                    sensor_prior,
                    cell_type_key=perm_key,
                    sample_key=sample_key,
                    layer=layer,
                    min_expr_frac=min_expr_frac,
                    allow_self=allow_self,
                    availability_kwargs=availability_kwargs,
                )
            else:
                sender_perm, receiver_perm, availability_perm = _compute_availability_scores(
                    adata_perm_base,
                    enzyme_prior,
                    sensor_prior,
                    celltype_col=perm_key,
                    layer=layer,
                    min_expr_frac=min_expr_frac,
                    **availability_kwargs
                )

                # 3. 构建置换事件
                events_perm = _make_cell_mesh_events(
                    sender_perm,
                    receiver_perm,
                    allow_self=allow_self,
                    cell_counts=availability_perm.get("cell_counts"),
                )
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
        if perm_key in adata_perm_base.obs:
            del adata_perm_base.obs[perm_key]

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
    min_expr_frac: Optional[float] = MIN_EXPR_FRAC,
    allow_self: bool = True,
    n_perms: int = 0,
    random_state: int = 0,
    eps_num: float = METABOLITE_AVAILABILITY_DEFAULTS["eps_num"],
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    sample_mode: Literal["pooled_stratified", "sample_aware"] = "pooled_stratified",
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
        min_expr_frac: 可选 receiver 表达比例 gate；None 表示不启用
        allow_self: 是否允许自分泌通信
        n_perms: 置换检验次数,0 表示不进行置换检验
        random_state: 随机种子
        eps_num: bounded median contrast 的数值保护常数,默认 1e-12
        min_cells : int
            Minimum number of cells required to construct a pseudobulk expression unit.
            In ``pooled_stratified`` mode, this filters cell types using their total
            cell count across the full AnnData object. In ``sample_aware`` mode, this
            filters each ``(sample, cell type)`` unit independently. Units below this
            threshold are excluded from pseudobulk calculation and are represented as
            missing values (NA), not zero, in sample-level outputs.
        sample_mode : {"pooled_stratified", "sample_aware"}, default="pooled_stratified"
            ``pooled_stratified`` computes pooled cell-type pseudobulks across all
            cells. When ``sample_key`` is provided, it is used only to stratify label
            permutations within samples.

            ``sample_aware`` computes pseudobulks, P/C/E scores, sender scores,
            receiver scores, and event scores separately for each ``(sample, cell
            type)`` unit, then aggregates event scores across samples.

    返回:
        CellMeshResult 对象,包含所有计算结果

    内部逻辑说明:
        1. enzyme_metabolite 作为标准 enzyme prior 直接传入 availability 计算，
           availability 内部负责 role 到 direction 的映射:
           - production → product → 进入 P (产生) 矩阵
           - degradation → substrate → 进入 C (消耗) 矩阵
           - export → exporter → 进入 E (外排) 矩阵
        2. sender score 以 production 的正向 median contrast 为必要锚点；
           exporter 高于背景时加分，消耗型复合酶水平高于背景时惩罚。
           缺少 exporter 或 consumption prior 时对应 factor 为 1。
        3. sensor_score 是 sensor 表达相对 eligible cell-type median 的正向 contrast
        4. cell_mesh_score = sqrt(metabolite_availability * sensor_score)
    """
    # 验证输入
    _validate_sample_key(adata, sample_key, sample_mode)

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
        'eps_num': eps_num,
        'min_cells': min_cells,
    }

    if sample_mode == "sample_aware":
        sender_scores, receiver_scores, events, availability_results = _compute_sample_aware_scores(
            adata,
            enzyme_prior,
            sensor_prior,
            cell_type_key=cell_type_key,
            sample_key=sample_key,
            layer=layer,
            min_expr_frac=min_expr_frac,
            allow_self=allow_self,
            availability_kwargs=availability_kwargs,
        )
    else:
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
        events = _make_cell_mesh_events(
            sender_scores,
            receiver_scores,
            allow_self=allow_self,
            cell_counts=availability_results.get("cell_counts"),
        )

    # 计算显著性（按 sensor type 分别计算）
    if sample_mode == "pooled_stratified" and not events.empty:
        events["inference_mode"] = "pooled_stratified"

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
        availability_kwargs=availability_kwargs,
        sample_mode=sample_mode,
    )
    event_attrs = events.attrs.copy()

    # 计算置信等级
    if not events.empty:
        events["confidence_tier"] = events.apply(_confidence_tier, axis=1)
        events = events.sort_values(["fdr", "cell_mesh_score"], ascending=[True, False], na_position="last").reset_index(drop=True)
        events.attrs.update(event_attrs)
    elif "confidence_tier" not in events.columns:
        events["confidence_tier"] = pd.Series(dtype=object)
        events.attrs.update(event_attrs)

    # 整理参数
    parameters = {
        "method": "CELL MESH",
        "acronym": "Metabolite-mediated Event Scoring with Sensor Hierarchies",
        "algorithm": "bounded cell-type median contrast scoring",
        "cell_type_key": cell_type_key,
        "sample_key": sample_key,
        "sample_mode": sample_mode,
        "layer": layer,
        "min_expr_frac": min_expr_frac,
        "allow_self": allow_self,
        "n_perms": n_perms,
        "random_state": random_state,
        **availability_kwargs
    }

    return CellMeshResult(
        events=events,
        sender_scores=sender_scores,
        receiver_scores=receiver_scores,
        parameters=parameters,
        availability_results=availability_results,
        sample_validation=availability_results.get("sample_validation"),
        sample_sender_scores=availability_results.get("sample_sender_scores"),
        sample_receiver_scores=availability_results.get("sample_receiver_scores"),
        sample_events=availability_results.get("sample_events"),
    )
