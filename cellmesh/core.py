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
    METABOLITE_AVAILABILITY_DEFAULTS,
    MISSING_EXPORT_SCORE,
)

from .database import _normalize_hmdb_id, load_cell_mesh_database, validate_priors
from .score import (
    _build_prior_role_coverage,
    _validate_export_weight,
    _validate_min_expr_frac,
    _validate_min_cells,
    _validate_pce_reference,
    _validate_receiver_reference,
    _validate_sender_abundance_exponent,
    compute_metabolite_availability,
    compute_sensor_scores
)
from .preprocess import (
    _compute_celltype_fractions,
    _validated_celltype_labels,
    _validated_gene_names,
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
    "sender_n_cells",
    "receiver_n_cells",
    "sender_passes_min_cells",
    "receiver_passes_min_cells",
    "passes_min_cells",
    "cell_mesh_score",
]

SAMPLE_MODES = {"pooled_stratified", "sample_aware"}
EVENT_KEY_COLUMNS = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "sensor_type"]
SAMPLE_AWARE_EVENT_COLUMNS = EVENT_KEY_COLUMNS + [
    "metabolite_availability_median",
    "sensor_score_median",
    "sensor_expr_frac_median",
    "sender_n_cells",
    "receiver_n_cells",
    "sender_passes_min_cells",
    "receiver_passes_min_cells",
    "passes_min_cells",
    "cell_mesh_score",
    "event_score_median",
    "event_score_iqr",
    "n_samples_coobserved",
    "n_samples_positive",
    "event_prevalence",
    "n_samples_passing_min_cells",
    "min_cells_pass_prevalence",
    "inference_mode",
]


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
    celltype_qc: Optional[pd.DataFrame] = None

    def to_csv(self, prefix: str) -> None:
        """
        将结果保存为 CSV 文件

        参数:
            prefix: 保存路径前缀
        """
        self.events.to_csv(f"{prefix}.events.csv", index=False)
        self.sender_scores.to_csv(f"{prefix}.sender_scores.csv", index=False)
        self.receiver_scores.to_csv(f"{prefix}.receiver_scores.csv", index=False)
        if self.celltype_qc is not None:
            self.celltype_qc.to_csv(f"{prefix}.celltype_qc.csv")


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
    left_id = _normalize_hmdb_id(left)
    right_id = _normalize_hmdb_id(right)
    if pd.isna(left_id) or pd.isna(right_id):
        return False
    return left_id == right_id


def _compute_availability_scores(
    adata,
    enzyme_prior: pd.DataFrame,
    sensor_prior: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_expr_frac: Optional[float] = MIN_EXPR_FRAC,
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    cell_fractions: Optional[pd.Series] = None,
    sender_abundance_exponent: float = METABOLITE_AVAILABILITY_DEFAULTS[
        "sender_abundance_exponent"
    ],
    pce_reference: str = METABOLITE_AVAILABILITY_DEFAULTS["pce_reference"],
    export_weight: float = METABOLITE_AVAILABILITY_DEFAULTS["export_weight"],
    receiver_reference: str = METABOLITE_AVAILABILITY_DEFAULTS[
        "receiver_reference"
    ],
    _prior_role_coverage: Optional[Dict[tuple[str, str], bool]] = None,
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
    # Compute abundance-adjusted P/C/E and the continuous sender score.
    avail_results = compute_metabolite_availability(
        adata,
        enzyme_prior,
        celltype_col=celltype_col,
        layer=layer,
        min_cells=min_cells,
        return_intermediates=True,
        cell_fractions=cell_fractions,
        sender_abundance_exponent=sender_abundance_exponent,
        pce_reference=pce_reference,
        export_weight=export_weight,
        _prior_role_coverage=_prior_role_coverage,
    )
    
    availability = avail_results['availability']

    # sender_scores is the metabolite availability matrix kept as a separate result.
    sender_scores = availability
    
    # 计算新的 sensor scores
    receiver_scores = compute_sensor_scores(
        adata,
        sensor_prior,
        celltype_col=celltype_col,
        layer=layer,
        min_expr_frac=min_expr_frac,
        min_cells=min_cells,
        pseudobulk=avail_results.get("pseudobulk"),
        expr_frac=avail_results.get("expr_frac"),
        cell_counts=avail_results.get("cell_counts"),
        cell_fractions=avail_results.get("cell_fractions"),
        receiver_reference=receiver_reference,
    )
    avail_results["receiver_reference"] = receiver_reference
    
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


def _validate_n_perms(value: Any) -> int:
    """Return a non-negative integer permutation count."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError("n_perms must be an integer")
    count = int(value)
    if count < 0:
        raise ValueError("n_perms must be greater than or equal to 0")
    return count


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
    sample_totals = validation.groupby("sample", observed=True)["n_cells"].transform("sum")
    validation["cell_fraction"] = validation["n_cells"] / sample_totals
    validation["passes_min_cells"] = validation["n_cells"] >= min_cells
    passing_celltypes = validation[validation["passes_min_cells"]].groupby("sample", observed=True).size()
    passing_samples = validation[validation["passes_min_cells"]].groupby("cell_type", observed=True).size()
    observed_counts = validation.groupby("sample", observed=True).size()
    observed_samples = validation.groupby("cell_type", observed=True).size()
    validation["n_celltypes_passing_min_cells_in_sample"] = (
        validation["sample"].map(passing_celltypes).fillna(0).astype(int)
    )
    validation["n_samples_passing_min_cells_for_celltype"] = (
        validation["cell_type"].map(passing_samples).fillna(0).astype(int)
    )
    validation["n_observed_celltypes_in_sample"] = (
        validation["sample"].map(observed_counts).fillna(0).astype(int)
    )
    validation["n_samples_observed_for_celltype"] = (
        validation["cell_type"].map(observed_samples).fillna(0).astype(int)
    )
    return validation


def _sample_cell_fractions(validation: pd.DataFrame) -> Dict[str, pd.Series]:
    """Calculate per-sample cell-type fractions from all observed units."""
    fractions: Dict[str, pd.Series] = {}
    for sample, group in validation.groupby("sample", observed=True):
        values = group.set_index("cell_type")["cell_fraction"].astype(float)
        values.index = values.index.astype(str)
        values.name = "cell_fraction"
        fractions[str(sample)] = values
    return fractions


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


def _qc_bool_sum(values: pd.Series) -> int:
    """Count true QC flags while treating missing structural units as false."""
    return int(pd.Series(values, dtype="boolean").fillna(False).sum())


def _qc_bool_any(values: pd.Series) -> bool:
    """Return whether any observed unit passes QC."""
    return bool(pd.Series(values, dtype="boolean").fillna(False).any())


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
    original_cell_fractions = _sample_cell_fractions(validation)
    # ``min_cells`` is QC-only. The null uses every observed analysis cell and
    # therefore the same score/reference universe as the observed data.
    adata_perm_base = adata.copy()

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
                cell_fractions_by_sample=original_cell_fractions,
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
    cell_fractions_by_sample: Optional[Dict[str, pd.Series]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    min_cells = availability_kwargs.get("min_cells", METABOLITE_AVAILABILITY_DEFAULTS["min_cells"])
    validation = _sample_validation_table(adata, sample_key, cell_type_key, min_cells)
    if cell_fractions_by_sample is None:
        cell_fractions_by_sample = _sample_cell_fractions(validation)
    all_celltypes = pd.Index(sorted(validation["cell_type"].unique()))

    sample_sender = []
    sample_receiver = []
    sample_events = []
    availability_by_sample: Dict[str, Dict[str, Any]] = {}
    valid_sample_names = []

    for sample in sorted(validation["sample"].unique()):
        # Every observed (sample, cell type) unit participates in scoring.
        # ``min_cells`` is carried only as a QC annotation.
        valid_sample_names.append(sample)

        obs = adata.obs
        mask = obs[sample_key].astype(str).values == sample
        adata_sample = adata[mask, :].copy()
        sender_scores, receiver_scores, availability_results = _compute_availability_scores(
            adata_sample,
            enzyme_prior,
            sensor_prior,
            celltype_col=cell_type_key,
            layer=layer,
            min_expr_frac=min_expr_frac,
            cell_fractions=cell_fractions_by_sample[str(sample)],
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
            min_cells=min_cells,
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
        sample_sender_scores.groupby(level=["metabolite", "hmdb_id"]).median()
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
                sensor_score=("sensor_score", "median"),
                sensor_expr_frac=("sensor_expr_frac", "median"),
                receiver_n_cells=("receiver_n_cells", "sum"),
                receiver_cell_fraction=("receiver_cell_fraction", "median"),
                n_observed_samples=("sample", "nunique"),
                n_samples_passing_min_cells=(
                    "receiver_passes_min_cells",
                    _qc_bool_sum,
                ),
                receiver_passes_min_cells=(
                    "receiver_passes_min_cells",
                    _qc_bool_any,
                ),
            )
        )
        receiver_scores["min_cells_pass_prevalence"] = (
            receiver_scores["n_samples_passing_min_cells"]
            / receiver_scores["n_observed_samples"]
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
        # Preserve the distinction between an event that is not computable in
        # this sample (NA) and a computed event that fails cell-count QC (False).
        for qc_column in [
            "sender_passes_min_cells",
            "receiver_passes_min_cells",
            "passes_min_cells",
        ]:
            sample_events_complete[qc_column] = sample_events_complete[
                qc_column
            ].astype("boolean")
        sample_events_df = sample_events_complete

        grouped_scores = sample_events_df.groupby(event_key)["cell_mesh_score"]
        score_summary = grouped_scores.agg(
            event_score_median="median",
            n_samples_coobserved="count",
            n_samples_positive=lambda s: int((pd.to_numeric(s, errors="coerce").dropna() > 0).sum()),
        ).reset_index()
        score_summary["event_score_iqr"] = grouped_scores.apply(_sample_stat_iqr).values
        score_summary["event_prevalence"] = grouped_scores.apply(_sample_positive_prevalence).values
        qc_summary = (
            sample_events_df.groupby(event_key)["passes_min_cells"]
            .apply(_qc_bool_sum)
            .rename("n_samples_passing_min_cells")
            .reset_index()
        )
        score_summary = score_summary.merge(qc_summary, on=event_key, how="left")
        score_summary["min_cells_pass_prevalence"] = (
            score_summary["n_samples_passing_min_cells"]
            / score_summary["n_samples_coobserved"]
        )

        meta_summary = (
            sample_events_df.groupby(event_key, as_index=False)
            .agg(
                metabolite_availability_median=("metabolite_availability", "median"),
                sensor_score_median=("sensor_score", "median"),
                sensor_expr_frac_median=("sensor_expr_frac", "median"),
                sender_n_cells=("sender_n_cells", "sum"),
                receiver_n_cells=("receiver_n_cells", "sum"),
                sender_passes_min_cells=(
                    "sender_passes_min_cells",
                    _qc_bool_any,
                ),
                receiver_passes_min_cells=(
                    "receiver_passes_min_cells",
                    _qc_bool_any,
                ),
                passes_min_cells=(
                    "passes_min_cells",
                    _qc_bool_any,
                ),
            )
        )
        events = (
            meta_summary.merge(score_summary, on=event_key, how="inner")
            .query("n_samples_coobserved > 0")
            .assign(
                cell_mesh_score=lambda df: df["event_score_median"],
                inference_mode="sample_aware",
            )
            .sort_values("cell_mesh_score", ascending=False, na_position="last")
            .reset_index(drop=True)
        )
    else:
        events = pd.DataFrame(columns=SAMPLE_AWARE_EVENT_COLUMNS)

    availability_results = {
        "sample_validation": validation,
        "celltype_qc": validation.set_index(["sample", "cell_type"])[
            ["n_cells", "cell_fraction", "passes_min_cells"]
        ].copy(),
        "sample_sender_scores": sample_sender_scores,
        "sample_receiver_scores": sample_receiver_scores,
        "sample_events": sample_events_df,
        "availability_by_sample": availability_by_sample,
        "export_weight": availability_kwargs.get(
            "export_weight",
            METABOLITE_AVAILABILITY_DEFAULTS["export_weight"],
        ),
        "receiver_reference": availability_kwargs.get(
            "receiver_reference",
            METABOLITE_AVAILABILITY_DEFAULTS["receiver_reference"],
        ),
    }
    return sender_scores, receiver_scores, events, availability_results


def _make_cell_mesh_events(
    sender_scores: pd.DataFrame,
    receiver_scores: pd.DataFrame,
    allow_self: bool,
    cell_counts: Optional[pd.Series] = None,
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
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
                idx for idx in sender_scores.index if _same_hmdb(idx[1], hmdb_id)
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
                sender_n_cells = (
                    int(cell_counts.loc[sender])
                    if cell_counts is not None and sender in cell_counts.index
                    else np.nan
                )
                receiver_n_cells = int(rr["receiver_n_cells"])
                sender_passes_min_cells = bool(
                    pd.notna(sender_n_cells) and sender_n_cells >= min_cells
                )
                receiver_passes_min_cells = bool(
                    rr.get(
                        "receiver_passes_min_cells",
                        receiver_n_cells >= min_cells,
                    )
                )

                rows.append(
                    {
                        "sender": sender,
                        "receiver": receiver,
                        "metabolite": sender_idx[0] if isinstance(sender_idx, tuple) else metabolite,
                        "hmdb_id": sender_idx[1] if isinstance(sender_idx, tuple) else hmdb_id,
                        "sensor_gene": rr["sensor_gene"],
                        "sensor_type": rr["sensor_type"],
                        "metabolite_availability": availability,
                        "sensor_score": sensor_score,
                        "sensor_expr_frac": float(rr["sensor_expr_frac"]),
                        "sender_n_cells": sender_n_cells,
                        "receiver_n_cells": receiver_n_cells,
                        "sender_passes_min_cells": sender_passes_min_cells,
                        "receiver_passes_min_cells": receiver_passes_min_cells,
                        "passes_min_cells": bool(
                            sender_passes_min_cells and receiver_passes_min_cells
                        ),
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
        带有经验 p 值、global FDR 和 sensor-type FDR 的事件 DataFrame
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
        out = obs_events.copy()
        out["perm_pvalue"] = np.nan
        return _assign_fdr_columns(out)

    key_cols = ["sender", "receiver", "metabolite", "hmdb_id", "sensor_gene", "sensor_type"]
    obs_keys = obs_events[key_cols].astype(str).agg("|".join, axis=1)
    ge_counts = np.zeros(len(obs_events), dtype=int)
    obs_score = obs_events["cell_mesh_score"].to_numpy(dtype=float)

    min_cells = availability_kwargs.get("min_cells", METABOLITE_AVAILABILITY_DEFAULTS["min_cells"])
    original_cell_fractions = _compute_celltype_fractions(
        adata,
        cell_type_key,
    )
    # ``min_cells`` is QC-only. Permute all observed analysis cells so the null
    # uses the same cell-type universe as the observed scores.
    adata_perm_base = adata.copy()

    rng = np.random.default_rng(random_state)
    original = adata_perm_base.obs[cell_type_key].copy()
    sample_labels = adata_perm_base.obs[sample_key].copy() if sample_key is not None else None
    perm_key = "_cell_mesh_perm_label"

    try:
        for perm_idx in range(n_perms):
            # 1. 打乱全部分析细胞的细胞类型标签；min_cells 只更新 QC 标记。
            adata_perm_base.obs[perm_key] = _permute_labels(original, sample_labels, rng).values

            # 2. 重新计算 availability 和得分
            sender_perm, receiver_perm, availability_perm = _compute_availability_scores(
                adata_perm_base,
                enzyme_prior,
                sensor_prior,
                celltype_col=perm_key,
                layer=layer,
                min_expr_frac=min_expr_frac,
                cell_fractions=original_cell_fractions,
                **availability_kwargs
            )

            # 3. 构建置换事件
            events_perm = _make_cell_mesh_events(
                sender_perm,
                receiver_perm,
                allow_self=allow_self,
                cell_counts=availability_perm.get("cell_counts"),
                min_cells=min_cells,
            )

            # 4. Compare the exact event key. A key absent from a permutation
            # has null score zero, matching the sample-aware implementation.
            if events_perm.empty:
                perm_values = np.zeros(len(obs_events), dtype=float)
            else:
                perm_scores = (
                    events_perm.assign(
                        _key=events_perm[key_cols].astype(str).agg("|".join, axis=1)
                    )
                    .groupby("_key")["cell_mesh_score"]
                    .max()
                )
                perm_values = (
                    perm_scores.reindex(obs_keys.to_numpy())
                    .fillna(0.0)
                    .to_numpy(dtype=float)
                )
            ge_counts += perm_values >= obs_score

    finally:
        if perm_key in adata_perm_base.obs:
            del adata_perm_base.obs[perm_key]

    # 计算 p 值
    p = (ge_counts + 1) / (n_perms + 1)
    out = obs_events.copy()
    out["perm_pvalue"] = p

    return _assign_fdr_columns(out)


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
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    sample_mode: Literal["pooled_stratified", "sample_aware"] = "pooled_stratified",
    sender_abundance_exponent: float = METABOLITE_AVAILABILITY_DEFAULTS[
        "sender_abundance_exponent"
    ],
    pce_reference: Literal["mean", "median"] = METABOLITE_AVAILABILITY_DEFAULTS[
        "pce_reference"
    ],
    export_weight: float = METABOLITE_AVAILABILITY_DEFAULTS["export_weight"],
    receiver_reference: Literal["mean", "median"] = METABOLITE_AVAILABILITY_DEFAULTS[
        "receiver_reference"
    ],
) -> CellMeshResult:
    """
    运行 CELL MESH: Metabolite-mediated Event Scoring with Sensor Hierarchies.
    完全基于 metabolite availability 算法计算代谢物通信事件。

    参数:
        adata: AnnData 对象,包含单细胞表达数据
        enzyme_metabolite: 酶-代谢物关系先验表,默认使用内置数据库
            必需列:metabolite, hmdb_id, gene, role, reaction
            可选列:evidence_level, source
            role 取值:production (产生)、degradation (降解)、export (外排)
        metabolite_sensor: 代谢物-传感器关系先验表,默认使用内置数据库
        cell_type_key: 细胞类型列名,默认为 "cell_type"
        sample_key: 样本列名,用于置换检验时的样本内置换
        layer: 使用的表达层,None 表示使用 adata.X
        min_expr_frac: 可选 receiver 表达比例 gate；必须在 [0, 1]，None 表示不启用
        allow_self: 是否允许自分泌通信
        n_perms: 非负整数置换检验次数，0 表示不进行置换检验
        random_state: 随机种子
        min_cells : int
            Positive-integer cell-count QC threshold. Every observed cell type contributes to
            pseudobulk construction, fractions, references, scores, events,
            and permutations regardless of this value. In
            ``pooled_stratified`` mode the QC flag uses the full-data cell-type
            count; in ``sample_aware`` mode it uses each ``(sample, cell type)``
            count. Changing ``min_cells`` changes QC annotations only.
        sample_mode : {"pooled_stratified", "sample_aware"}, default="pooled_stratified"
            ``pooled_stratified`` computes pooled cell-type pseudobulks across all
            cells. When ``sample_key`` is provided, it is used only to stratify label
            permutations within samples.

            ``sample_aware`` computes pseudobulks, P/C/E scores, sender scores,
            receiver scores, and event scores separately for each ``(sample, cell
            type)`` unit, then aggregates event scores across samples.
        sender_abundance_exponent : float, default=1.0
            Exponent applied to sender cell fraction before P/C/E construction:
            ``reaction_activity * cell_fraction ** sender_abundance_exponent``.
            ``1`` preserves linear abundance adjustment and ``0`` disables it.
        pce_reference : {"mean", "median"}, default="mean"
            Statistic calculated over strictly positive abundance-adjusted
            capacities to define each metabolite/direction reference.
        export_weight : float, default=0.2
            Bounded contribution of normalized exporter evidence to the sender
            score. Must be in [0, 1]; 0 exactly disables exporter modulation.
        receiver_reference : {"mean", "median"}, default="median"
            Unweighted statistic calculated over strictly positive observed
            cell-type mean sensor expression to define the receiver reference.

    返回:
        CellMeshResult 对象,包含所有计算结果

    内部逻辑说明:
        1. enzyme_metabolite 作为标准 enzyme prior 直接传入 availability 计算，
           availability 内部负责 role 到 direction 的映射:
           - production → product → 进入 P (产生) 矩阵
           - degradation → substrate → 进入 C (消耗) 矩阵
           - export → exporter → 进入 E (外排) 矩阵
        2. 每个代谢物的 P/C/E 默认使用所有 observed cell types 中严格正值的
           算术均值作为 reference；pce_reference="median" 时改用中位数。
           随后统一转换为 x / (x + positive_reference)。
           每条 reaction activity 在 P/C/E normalization 前乘对应细胞类型中的
           cell_fraction ** sender_abundance_exponent；sample-aware 模式按样本
           分别计算 fraction。
           base_sender_score = P_score ** 2 / (P_score + C_score)；P 为必要锚点，
           C 作为平滑竞争项。正式 sender score 为
           base_sender_score * ((1-export_weight) + export_weight*E_effective)。
           E 可计算时 E_effective=E_score；外排先验缺失或其基因未测到时取固定
           中性值 0.5；先验可计算但表达全零时取 0。
           缺少 consumption prior 时 C_score = 0，因此 base_sender_score = P_score；
           正式 sender score 仍乘以上述 E_factor。
           prior 存在但基因未测到、或基因已测到但表达全零时，C_score 同样为
           0，并由 consumption_status 区分其证据状态。
        3. sensor_score 默认使用正表达 observed cell types 的等权中位数
           R_ref（receiver_reference="mean" 时使用等权算术均值），并计算
           R_score = R / (R + R_ref)。receiver cell fraction 不进入得分或
           reference；min_expr_frac 仍可作为表达比例 gate。
        4. 单样本/pooled 事件满足
           cell_mesh_score = sqrt(metabolite_availability * sensor_score)。
           sample-aware 聚合事件的 cell_mesh_score 是样本级事件分数的中位数；
           组件列分别命名为 metabolite_availability_median 和
           sensor_score_median，二者仅是描述性跨样本汇总。
        5. min_cells 只生成 cell-count QC 标记，不过滤上述计算；所有事件保留，
           可视化默认仅展示 passes_min_cells=True 的事件。
    """
    # 验证输入
    _validate_sample_key(adata, sample_key, sample_mode)
    min_cells = _validate_min_cells(min_cells)
    min_expr_frac = _validate_min_expr_frac(min_expr_frac)
    n_perms = _validate_n_perms(n_perms)
    if cell_type_key not in adata.obs:
        raise KeyError(f"{cell_type_key!r} not found in adata.obs")
    if len(adata.obs) == 0:
        raise ValueError("No observed cell types are available for analysis")
    _validated_celltype_labels(adata, cell_type_key)
    _validated_gene_names(adata)
    sender_abundance_exponent = _validate_sender_abundance_exponent(
        sender_abundance_exponent
    )
    pce_reference = _validate_pce_reference(pce_reference)
    export_weight = _validate_export_weight(export_weight)
    receiver_reference = _validate_receiver_reference(receiver_reference)

    # 加载默认数据库
    if enzyme_metabolite is None or metabolite_sensor is None:
        default_enzyme, default_sensor = load_cell_mesh_database()
        enzyme_metabolite = default_enzyme if enzyme_metabolite is None else enzyme_metabolite
        metabolite_sensor = default_sensor if metabolite_sensor is None else metabolite_sensor

    # Keep the unfiltered enzyme prior only for coverage-status reporting.
    # Numerical scoring still uses the validated, expression-compatible prior.
    prior_role_coverage = _build_prior_role_coverage(
        enzyme_metabolite,
        adata.var_names,
    )

    # 验证先验
    enzyme_prior, sensor_prior = validate_priors(enzyme_metabolite, metabolite_sensor, adata.var_names)

    if enzyme_prior.empty:
        raise ValueError("No enzyme prior genes found in adata.var_names")
    if sensor_prior.empty:
        raise ValueError("No sensor genes found in adata.var_names")

    # 计算 availability 和得分
    availability_kwargs = {
        'min_cells': min_cells,
        'sender_abundance_exponent': sender_abundance_exponent,
        'pce_reference': pce_reference,
        'export_weight': export_weight,
        'receiver_reference': receiver_reference,
        '_prior_role_coverage': prior_role_coverage,
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
            min_cells=min_cells,
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

    if not events.empty:
        events = events.sort_values(
            ["fdr_sensor_type", "cell_mesh_score"],
            ascending=[True, False],
            na_position="last",
        ).reset_index(drop=True)
    events.attrs.update(event_attrs)

    # 整理参数
    parameters = {
        "method": "CELL MESH",
        "acronym": "Metabolite-mediated Event Scoring with Sensor Hierarchies",
        "algorithm": "sender-abundance-adjusted positive-reference saturation scoring",
        "population_adjustment": "sender_cell_fraction_power_before_pce_normalization",
        "pce_normalization": "positive_reference_saturation",
        "pce_reference": pce_reference,
        "sender_formula": "base_sender_score*((1-export_weight)+export_weight*E_effective)",
        "base_sender_formula": "P_score^2/(P_score+C_score)",
        "export_in_sender_score": True,
        "export_weight": export_weight,
        "missing_export_score": MISSING_EXPORT_SCORE,
        "receiver_normalization": "positive_reference_saturation",
        "receiver_reference": receiver_reference,
        "receiver_formula": "R/(R+R_ref)",
        "receiver_abundance_adjustment": "none",
        "min_cells_role": "qc_annotation_only",
        "calculation_celltypes": "all_observed",
        "cell_fraction_denominator": (
            "within_sample_all_cells"
            if sample_mode == "sample_aware"
            else "all_adata_cells"
        ),
        "sample_aware_component_aggregation": "median",
        "sample_aware_event_score": "median_of_sample_level_cell_mesh_scores",
        "cell_type_key": cell_type_key,
        "sample_key": sample_key,
        "sample_mode": sample_mode,
        "layer": layer,
        "min_expr_frac": min_expr_frac,
        "allow_self": allow_self,
        "n_perms": n_perms,
        "random_state": random_state,
        "min_cells": min_cells,
        "sender_abundance_exponent": sender_abundance_exponent,
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
        celltype_qc=availability_results.get("celltype_qc"),
    )
