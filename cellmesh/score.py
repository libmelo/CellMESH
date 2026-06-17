"""
Scoring algorithms for metabolite availability and sensor activity.
"""
from __future__ import annotations

from typing import Any, Dict, Optional
import warnings

import numpy as np
import pandas as pd
from scipy.stats import gmean

from .config import (
    METABOLITE_AVAILABILITY_DEFAULTS,
    MIN_EXPR_FRAC,
    ROLE_TO_DIRECTION,
    VALID_ROLES,
)
from .database import _valid_hmdb_mask
from .preprocess import _build_celltype_pseudobulk, _compute_celltype_expr_frac


def robust_minmax(
    x: np.ndarray,
    lower: float = 5,
    upper: float = 95,
    eps: float = 1e-8,
) -> np.ndarray:
    """
    Robust min-max 标准化,结果范围在 [0, 1] 之间
    """
    x = np.asarray(x, dtype=float)
    nan_mask = np.isnan(x)

    if np.all(nan_mask):
        return np.zeros_like(x)

    lo = np.nanpercentile(x, lower)
    hi = np.nanpercentile(x, upper)

    if np.isclose(hi, lo, rtol=1e-5):
        result = np.zeros_like(x)
    else:
        x_clip = np.clip(x, lo, hi)
        result = (x_clip - lo) / (hi - lo + eps)
        result = np.clip(result, 0.0, 1.0)

    result[nan_mask] = 0.0
    return result


def compute_sensor_scores(
    adata,
    sensor_prior: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    lower: float = METABOLITE_AVAILABILITY_DEFAULTS["lower"],
    upper: float = METABOLITE_AVAILABILITY_DEFAULTS["upper"],
    min_expr_frac: float = MIN_EXPR_FRAC,
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    pseudobulk: Optional[pd.DataFrame] = None,
    expr_frac: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    计算 sensor scores: robust min-max 标准化的 sensor 基因表达
    """
    if pseudobulk is None:
        pseudobulk = _build_celltype_pseudobulk(adata, celltype_col, layer, min_cells)
    if expr_frac is None:
        expr_frac = _compute_celltype_expr_frac(adata, celltype_col, layer, min_cells)

    valid_genes = [g for g in sensor_prior["sensor_gene"].unique() if g in pseudobulk.columns]
    if not valid_genes:
        return pd.DataFrame(
            columns=[
                "metabolite",
                "hmdb_id",
                "sensor_gene",
                "sensor_type",
                "receiver",
                "sensor_score",
                "sensor_expr_frac",
            ]
        )

    sensor_gene_scores = {}
    for gene in valid_genes:
        expr_values = pseudobulk[gene].values
        norm_scores = robust_minmax(expr_values, lower=lower, upper=upper)
        sensor_gene_scores[gene] = pd.Series(norm_scores, index=pseudobulk.index)

    rows = []
    for _, row in sensor_prior.iterrows():
        gene = row["sensor_gene"]
        if gene not in valid_genes:
            continue

        for receiver in pseudobulk.index:
            frac = expr_frac.loc[receiver, gene]
            score = sensor_gene_scores[gene].loc[receiver] if frac >= min_expr_frac else 0.0

            rows.append(
                {
                    "metabolite": row["metabolite"],
                    "hmdb_id": row["hmdb_id"],
                    "sensor_gene": gene,
                    "sensor_type": row["sensor_type"],
                    "receiver": receiver,
                    "sensor_score": score,
                    "sensor_expr_frac": frac,
                }
            )

    return pd.DataFrame(rows)


def _normalize_enzyme_metabolite(enzyme_metabolite: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize enzyme-metabolite prior rows for availability scoring.
    """
    df = enzyme_metabolite.copy()

    if "hmdb_id" not in df.columns and "HMDB_ID" in df.columns:
        df["hmdb_id"] = df["HMDB_ID"]
    if "direction" not in df.columns:
        if "role" not in df.columns:
            raise ValueError("enzyme_metabolite must contain either 'role' or 'direction'")
        df["role"] = df["role"].astype(str).str.lower()
        df = df[df["role"].isin(VALID_ROLES)]
        df["direction"] = df["role"].map(ROLE_TO_DIRECTION)
    else:
        df["direction"] = df["direction"].replace({"transporter": "exporter"})

    required_cols = ["metabolite", "hmdb_id", "gene", "direction"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"enzyme_metabolite is missing required column: {col}")

    if "reaction" not in df.columns:
        df["reaction"] = "unknown"
    if "weight" not in df.columns:
        df["weight"] = 1.0

    df["weight"] = pd.to_numeric(df["weight"], errors="coerce").fillna(1.0)
    df["hmdb_id"] = df["hmdb_id"].astype(object).where(pd.notna(df["hmdb_id"]), np.nan)
    df = df[_valid_hmdb_mask(df["hmdb_id"])]
    df = df[df["direction"].isin({"product", "substrate", "exporter"})]

    return df[["metabolite", "hmdb_id", "reaction", "gene", "direction", "weight"]].reset_index(drop=True)


def _get_reaction_gene_sets(reaction_table: pd.DataFrame) -> pd.DataFrame:
    """
    为每个反应提取基因集合,处理多基因情况,去重,支持各种分隔符
    """
    df = reaction_table.copy()
    if df.empty:
        return pd.DataFrame(columns=["metabolite", "hmdb_id", "reaction", "direction", "genes", "weights"])

    def get_reaction_id(row):
        hmdb_str = str(row["hmdb_id"]) if pd.notna(row["hmdb_id"]) else "nan"
        return f"{row['metabolite']}|{hmdb_str}|{row['reaction']}|{row['direction']}"

    def parse_genes(gene_str):
        if pd.isna(gene_str) or gene_str == "":
            return []
        for sep in [";", ",", "|"]:
            if sep in gene_str:
                genes = [g.strip() for g in gene_str.split(sep)]
                genes = [g.split("[")[0].strip() for g in genes]
                return [g for g in genes if g]
        gene = gene_str.split("[")[0].strip()
        return [gene] if gene else []

    df["_reaction_id"] = df.apply(get_reaction_id, axis=1)
    df["_genes"] = df["gene"].apply(parse_genes)

    reaction_dict = {}
    for _, row in df.iterrows():
        rid = row["_reaction_id"]
        genes = row["_genes"]
        weight = row["weight"]

        if rid not in reaction_dict:
            reaction_dict[rid] = {
                "metabolite": row["metabolite"],
                "hmdb_id": row["hmdb_id"],
                "reaction": row["reaction"],
                "direction": row["direction"],
                "gene_weights": {},
            }

        for gene in genes:
            if gene in reaction_dict[rid]["gene_weights"]:
                if weight > reaction_dict[rid]["gene_weights"][gene]:
                    reaction_dict[rid]["gene_weights"][gene] = weight
                    warnings.warn(f"Gene {gene} appeared multiple times in reaction {rid}, using maximum weight {weight}")
            else:
                reaction_dict[rid]["gene_weights"][gene] = weight

    result = []
    for _, info in reaction_dict.items():
        result.append(
            {
                "metabolite": info["metabolite"],
                "hmdb_id": info["hmdb_id"],
                "reaction": info["reaction"],
                "direction": info["direction"],
                "genes": list(info["gene_weights"].keys()),
                "weights": list(info["gene_weights"].values()),
            }
        )

    return pd.DataFrame(result)


def _compute_reaction_scores(
    pseudobulk: pd.DataFrame,
    reaction_genes: pd.DataFrame,
) -> pd.DataFrame:
    """
    计算每个反应在每个细胞类型中的得分,考虑基因权重
    """
    scores = []
    reaction_ids = []

    for _, row in reaction_genes.iterrows():
        hmdb_str = str(row["hmdb_id"]) if pd.notna(row["hmdb_id"]) else "nan"
        rid = f"{row['metabolite']}|{hmdb_str}|{row['reaction']}|{row['direction']}"
        reaction_ids.append(rid)

        valid_gene_idx = [i for i, g in enumerate(row["genes"]) if g in pseudobulk.columns]

        if not valid_gene_idx:
            scores.append(pd.Series(0.0, index=pseudobulk.index))
        else:
            valid_genes = [row["genes"][i] for i in valid_gene_idx]
            valid_weights = np.array([row["weights"][i] for i in valid_gene_idx])

            # Weight scales each gene expression value before the ordinary
            # geometric mean; it is not a normalized weighted geometric mean.
            expr = pseudobulk[valid_genes].values * valid_weights.reshape(1, -1)
            geo_mean = gmean(expr + 1, axis=1) - 1
            scores.append(pd.Series(geo_mean, index=pseudobulk.index))

    return pd.DataFrame(scores, index=reaction_ids, columns=pseudobulk.index)


def _compute_PCE_matrices(
    reaction_scores: pd.DataFrame,
    reaction_genes: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    """
    计算 P (production), C (consumption), E (efflux) 矩阵
    """
    reaction_info = {}
    metabolites_dict = {}

    for _, row in reaction_genes.iterrows():
        hmdb_str = str(row["hmdb_id"]) if pd.notna(row["hmdb_id"]) else "nan"
        rid = f"{row['metabolite']}|{hmdb_str}|{row['reaction']}|{row['direction']}"
        reaction_info[rid] = {"direction": row["direction"], "metabolite": row["metabolite"], "hmdb_id": row["hmdb_id"]}
        met_key = f"{row['metabolite']}|{hmdb_str}"
        if met_key not in metabolites_dict:
            metabolites_dict[met_key] = (row["metabolite"], row["hmdb_id"])

    cell_types = reaction_scores.columns
    met_index = pd.MultiIndex.from_tuples(list(metabolites_dict.values()), names=["metabolite", "hmdb_id"])

    P = pd.DataFrame(0.0, index=met_index, columns=cell_types)
    C = pd.DataFrame(0.0, index=met_index, columns=cell_types)
    E = pd.DataFrame(0.0, index=met_index, columns=cell_types)

    for rid in reaction_scores.index:
        info = reaction_info[rid]
        direction = info["direction"]
        met = info["metabolite"]
        hmdb = info["hmdb_id"]
        hmdb_str = str(hmdb) if pd.notna(hmdb) else "nan"
        met_idx = metabolites_dict[f"{met}|{hmdb_str}"]

        if direction == "product":
            P.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == "substrate":
            C.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == "exporter":
            E.loc[met_idx] += reaction_scores.loc[rid]

    return {"P": P, "C": C, "E": E}


def _safe_hmdb_compare(row: pd.Series, met: str, hmdb: Optional[str]) -> bool:
    if row["metabolite"] != met:
        return False
    if pd.isna(row["hmdb_id"]) or pd.isna(hmdb):
        return False
    return str(row["hmdb_id"]) == str(hmdb)


def _normalize_PCE(
    P: pd.DataFrame,
    C: pd.DataFrame,
    E: pd.DataFrame,
    reaction_genes: pd.DataFrame,
    lower: float = 5,
    upper: float = 95,
    missing_C_norm: float = 0.2,
    missing_E_norm: float = 0.5,
    eps: float = 1e-8,
) -> Dict[str, pd.DataFrame]:
    """
    对 P, C, E 矩阵进行标准化，结果都在 [0, 1] 范围内
    """
    missing_C_norm = np.clip(missing_C_norm, 0.0, 1.0)
    missing_E_norm = np.clip(missing_E_norm, 0.0, 1.0)

    P_norm = pd.DataFrame(index=P.index, columns=P.columns, dtype=float)
    C_norm = pd.DataFrame(missing_C_norm, index=P.index, columns=P.columns, dtype=float)
    E_norm = pd.DataFrame(missing_E_norm, index=P.index, columns=P.columns, dtype=float)

    for met_idx in P.index:
        met, hmdb = met_idx

        p_vals = P.loc[met_idx].values.flatten()
        P_norm.loc[met_idx] = robust_minmax(p_vals, lower=lower, upper=upper, eps=eps)

        has_substrate = any(
            _safe_hmdb_compare(row, met, hmdb) and row["direction"] == "substrate"
            for _, row in reaction_genes.iterrows()
        )
        if has_substrate and met_idx in C.index:
            c_vals = C.loc[met_idx].values.flatten()
            C_norm.loc[met_idx] = robust_minmax(c_vals, lower=lower, upper=upper, eps=eps)

        has_exporter = any(
            _safe_hmdb_compare(row, met, hmdb) and row["direction"] == "exporter"
            for _, row in reaction_genes.iterrows()
        )
        if has_exporter and met_idx in E.index:
            e_vals = E.loc[met_idx].values.flatten()
            E_norm.loc[met_idx] = robust_minmax(e_vals, lower=lower, upper=upper, eps=eps)

    return {"P_norm": P_norm, "C_norm": C_norm, "E_norm": E_norm}


def compute_metabolite_availability(
    adata,
    enzyme_metabolite: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    lower: float = METABOLITE_AVAILABILITY_DEFAULTS["lower"],
    upper: float = METABOLITE_AVAILABILITY_DEFAULTS["upper"],
    eps: float = METABOLITE_AVAILABILITY_DEFAULTS["eps"],
    beta: float = METABOLITE_AVAILABILITY_DEFAULTS["beta"],
    missing_C_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_C_norm"],
    missing_E_norm: float = METABOLITE_AVAILABILITY_DEFAULTS["missing_E_norm"],
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    return_intermediates: bool = True,
) -> Dict[str, Any]:
    """
    计算代谢物 availability,结果范围在 [0, 1] 之间
    """
    eps = max(eps, 1e-8)
    beta = max(beta, 0.0)

    pseudobulk = _build_celltype_pseudobulk(adata, celltype_col=celltype_col, layer=layer, min_cells=min_cells)
    expr_frac = _compute_celltype_expr_frac(adata, celltype_col=celltype_col, layer=layer, min_cells=min_cells)

    parsed_reactions = _normalize_enzyme_metabolite(enzyme_metabolite)
    reaction_genes = _get_reaction_gene_sets(parsed_reactions)
    reaction_scores = _compute_reaction_scores(pseudobulk, reaction_genes)
    PCE = _compute_PCE_matrices(reaction_scores, reaction_genes)
    P, C, E = PCE["P"], PCE["C"], PCE["E"]

    valid_mets = P.index[P.sum(axis=1) > 0]
    if len(valid_mets) == 0:
        result = {"availability": pd.DataFrame(), "metadata": pd.DataFrame()}
        if return_intermediates:
            result.update(
                {
                    "P": pd.DataFrame(),
                    "C": pd.DataFrame(),
                    "E": pd.DataFrame(),
                    "P_norm": pd.DataFrame(),
                    "C_norm": pd.DataFrame(),
                    "E_norm": pd.DataFrame(),
                    "pseudobulk": pd.DataFrame(),
                    "expr_frac": pd.DataFrame(),
                    "reaction_genes": pd.DataFrame(),
                }
            )
        return result

    P = P.loc[valid_mets]
    C = C.reindex(valid_mets, fill_value=0.0)
    E = E.reindex(valid_mets, fill_value=0.0)

    normalized = _normalize_PCE(
        P,
        C,
        E,
        reaction_genes,
        lower=lower,
        upper=upper,
        missing_C_norm=missing_C_norm,
        missing_E_norm=missing_E_norm,
        eps=eps,
    )
    P_norm, C_norm, E_norm = normalized["P_norm"], normalized["C_norm"], normalized["E_norm"]

    availability = pd.DataFrame(index=P_norm.index, columns=P_norm.columns, dtype=float)
    for met_idx in P_norm.index:
        p = P_norm.loc[met_idx].values.flatten()
        c = C_norm.loc[met_idx].values.flatten()
        e = E_norm.loc[met_idx].values.flatten()
        avail = p * np.power((1 - c), beta) * (0.8 + 0.2 * e)
        availability.loc[met_idx] = np.clip(avail, 0.0, 1.0)

    metadata = pd.DataFrame(index=P_norm.index)
    metadata["has_product"] = True
    metadata["has_substrate"] = [
        any(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "substrate" for _, row in reaction_genes.iterrows())
        for met, hmdb in P_norm.index
    ]
    metadata["has_exporter"] = [
        any(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "exporter" for _, row in reaction_genes.iterrows())
        for met, hmdb in P_norm.index
    ]
    metadata["n_product_reactions"] = [
        sum(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "product" for _, row in reaction_genes.iterrows())
        for met, hmdb in P_norm.index
    ]
    metadata["n_substrate_reactions"] = [
        sum(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "substrate" for _, row in reaction_genes.iterrows())
        for met, hmdb in P_norm.index
    ]
    metadata["n_exporter_reactions"] = [
        sum(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "exporter" for _, row in reaction_genes.iterrows())
        for met, hmdb in P_norm.index
    ]

    result = {"availability": availability.astype(float), "metadata": metadata}
    if return_intermediates:
        result.update(
            {
                "P": P.astype(float),
                "C": C.astype(float),
                "E": E.astype(float),
                "P_norm": P_norm.astype(float),
                "C_norm": C_norm.astype(float),
                "E_norm": E_norm.astype(float),
                "pseudobulk": pseudobulk.astype(float),
                "expr_frac": expr_frac.astype(float),
                "reaction_genes": reaction_genes,
            }
        )

    return result
