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
from .preprocess import (
    _build_celltype_pseudobulk,
    _compute_celltype_expr_frac,
    _eligible_celltype_counts,
)


def bounded_median_contrast(
    values: np.ndarray,
    eps: float = 1e-12,
) -> np.ndarray:
    """Compute bounded relative deviation from the finite-value median."""
    x = np.asarray(values, dtype=float)
    result = np.zeros_like(x, dtype=float)
    finite = np.isfinite(x)
    if not finite.any():
        return result

    finite_values = x[finite]
    if np.all(finite_values == 0):
        return result
    if np.any(finite_values < 0):
        raise ValueError("bounded_median_contrast only supports non-negative values")

    baseline = float(np.median(finite_values))
    denominator = finite_values + baseline
    valid = denominator > eps
    contrasted = np.zeros_like(finite_values, dtype=float)
    contrasted[valid] = (finite_values[valid] - baseline) / denominator[valid]
    result[finite] = np.clip(contrasted, -1.0, 1.0)
    return result


def compute_sensor_scores(
    adata,
    sensor_prior: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_expr_frac: Optional[float] = MIN_EXPR_FRAC,
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    eps_num: float = METABOLITE_AVAILABILITY_DEFAULTS["eps_num"],
    pseudobulk: Optional[pd.DataFrame] = None,
    expr_frac: Optional[pd.DataFrame] = None,
    cell_counts: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """
    Compute receiver scores as positive median-relative sensor expression.
    """
    if pseudobulk is None:
        pseudobulk = _build_celltype_pseudobulk(adata, celltype_col, layer, min_cells)
    if expr_frac is None:
        expr_frac = _compute_celltype_expr_frac(adata, celltype_col, layer, min_cells)
    if cell_counts is None:
        cell_counts = _eligible_celltype_counts(adata, celltype_col, min_cells)

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
                "receiver_n_cells",
            ]
        )

    sensor_gene_scores = {}
    for gene in valid_genes:
        expr_values = pseudobulk[gene].values
        contrast = bounded_median_contrast(expr_values, eps=eps_num)
        sensor_gene_scores[gene] = pd.Series(np.maximum(0.0, contrast), index=pseudobulk.index)

    rows = []
    for _, row in sensor_prior.iterrows():
        gene = row["sensor_gene"]
        if gene not in valid_genes:
            continue

        for receiver in pseudobulk.index:
            frac = expr_frac.loc[receiver, gene]
            score = sensor_gene_scores[gene].loc[receiver]
            if min_expr_frac is not None and frac < min_expr_frac:
                score = 0.0

            rows.append(
                {
                    "metabolite": row["metabolite"],
                    "hmdb_id": row["hmdb_id"],
                    "sensor_gene": gene,
                    "sensor_type": row["sensor_type"],
                    "receiver": receiver,
                    "sensor_score": score,
                    "sensor_expr_frac": frac,
                    "receiver_n_cells": int(cell_counts.loc[receiver]),
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


def _contrast_PCE(
    P: pd.DataFrame,
    C: pd.DataFrame,
    E: pd.DataFrame,
    reaction_genes: pd.DataFrame,
    eps_num: float = 1e-12,
) -> Dict[str, pd.DataFrame]:
    """Compute P/C/E median contrasts and prior-aware sender factors."""
    P_contrast = pd.DataFrame(index=P.index, columns=P.columns, dtype=float)
    C_contrast = pd.DataFrame(0.0, index=P.index, columns=P.columns, dtype=float)
    E_contrast = pd.DataFrame(0.0, index=P.index, columns=P.columns, dtype=float)
    has_substrate = pd.Series(False, index=P.index, dtype=bool)
    has_exporter = pd.Series(False, index=P.index, dtype=bool)

    for met_idx in P.index:
        met, hmdb = met_idx
        P_contrast.loc[met_idx] = bounded_median_contrast(P.loc[met_idx].to_numpy(), eps=eps_num)

        has_substrate.loc[met_idx] = any(
            _safe_hmdb_compare(row, met, hmdb) and row["direction"] == "substrate"
            for _, row in reaction_genes.iterrows()
        )
        if has_substrate.loc[met_idx]:
            C_contrast.loc[met_idx] = bounded_median_contrast(C.loc[met_idx].to_numpy(), eps=eps_num)

        has_exporter.loc[met_idx] = any(
            _safe_hmdb_compare(row, met, hmdb) and row["direction"] == "exporter"
            for _, row in reaction_genes.iterrows()
        )
        if has_exporter.loc[met_idx]:
            E_contrast.loc[met_idx] = bounded_median_contrast(E.loc[met_idx].to_numpy(), eps=eps_num)

    return {
        "P_contrast": P_contrast,
        "C_contrast": C_contrast,
        "E_contrast": E_contrast,
        "has_substrate": has_substrate,
        "has_exporter": has_exporter,
    }


def compute_metabolite_availability(
    adata,
    enzyme_metabolite: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    eps_num: float = METABOLITE_AVAILABILITY_DEFAULTS["eps_num"],
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    return_intermediates: bool = True,
) -> Dict[str, Any]:
    """
    Compute sender scores from median-relative P/C/E context.
    """
    eps_num = max(float(eps_num), np.finfo(float).eps)

    cell_counts = _eligible_celltype_counts(adata, celltype_col, min_cells)
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
                    "P_contrast": pd.DataFrame(),
                    "C_contrast": pd.DataFrame(),
                    "E_contrast": pd.DataFrame(),
                    "P_plus": pd.DataFrame(),
                    "relative_consumption_support": pd.DataFrame(),
                    "E_plus": pd.DataFrame(),
                    "pseudobulk": pd.DataFrame(),
                    "expr_frac": pd.DataFrame(),
                    "cell_counts": cell_counts,
                    "reaction_genes": pd.DataFrame(),
                }
            )
        return result

    P = P.loc[valid_mets]
    C = C.reindex(valid_mets, fill_value=0.0)
    E = E.reindex(valid_mets, fill_value=0.0)

    contrasted = _contrast_PCE(P, C, E, reaction_genes, eps_num=eps_num)
    P_contrast = contrasted["P_contrast"]
    C_contrast = contrasted["C_contrast"]
    E_contrast = contrasted["E_contrast"]
    P_plus = P_contrast.clip(lower=0.0)
    C_plus = C_contrast.clip(lower=0.0)
    E_plus = E_contrast.clip(lower=0.0)

    availability = pd.DataFrame(index=P.index, columns=P.columns, dtype=float)
    for met_idx in P.index:
        exporter_factor = 1.0 + E_plus.loc[met_idx] if contrasted["has_exporter"].loc[met_idx] else 1.0
        consumption_factor = 1.0 - C_plus.loc[met_idx] if contrasted["has_substrate"].loc[met_idx] else 1.0
        availability.loc[met_idx] = P_plus.loc[met_idx] * exporter_factor * consumption_factor
    availability = availability.clip(lower=0.0)

    metadata = pd.DataFrame(index=P.index)
    metadata["has_product"] = True
    metadata["has_substrate"] = [
        any(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "substrate" for _, row in reaction_genes.iterrows())
        for met, hmdb in P.index
    ]
    metadata["has_exporter"] = [
        any(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "exporter" for _, row in reaction_genes.iterrows())
        for met, hmdb in P.index
    ]
    metadata["n_product_reactions"] = [
        sum(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "product" for _, row in reaction_genes.iterrows())
        for met, hmdb in P.index
    ]
    metadata["n_substrate_reactions"] = [
        sum(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "substrate" for _, row in reaction_genes.iterrows())
        for met, hmdb in P.index
    ]
    metadata["n_exporter_reactions"] = [
        sum(_safe_hmdb_compare(row, met, hmdb) and row["direction"] == "exporter" for _, row in reaction_genes.iterrows())
        for met, hmdb in P.index
    ]

    result = {"availability": availability.astype(float), "metadata": metadata}
    if return_intermediates:
        result.update(
            {
                "P": P.astype(float),
                "C": C.astype(float),
                "E": E.astype(float),
                "P_contrast": P_contrast.astype(float),
                "C_contrast": C_contrast.astype(float),
                "E_contrast": E_contrast.astype(float),
                "P_plus": P_plus.astype(float),
                "relative_consumption_support": C_plus.astype(float),
                "E_plus": E_plus.astype(float),
                "pseudobulk": pseudobulk.astype(float),
                "expr_frac": expr_frac.astype(float),
                "cell_counts": cell_counts,
                "reaction_genes": reaction_genes,
            }
        )

    return result
