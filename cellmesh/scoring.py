from __future__ import annotations

import numpy as np
import pandas as pd


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -50, 50)))


def zscore_by_gene(mean_expr: pd.DataFrame, eps: float = 1e-8) -> pd.DataFrame:
    mu = mean_expr.mean(axis=0)
    sd = mean_expr.std(axis=0).replace(0, np.nan)
    return ((mean_expr - mu) / (sd + eps)).fillna(0.0)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    w = np.asarray(weights, dtype=float)
    v = np.asarray(values, dtype=float)
    denom = np.sum(np.abs(w))
    if denom == 0:
        return 0.0
    return float(np.sum(v * w) / denom)


def softmin(values: pd.Series, weights: pd.Series, tau: float = 1.0) -> float:
    v = np.asarray(values, dtype=float) * np.asarray(weights, dtype=float)
    if len(v) == 0:
        return 0.0
    return float(-tau * np.log(np.mean(np.exp(-v / tau))))


def aggregate_role_score(
    z_expr: pd.DataFrame,
    expr_frac: pd.DataFrame,
    enzyme_prior: pd.DataFrame,
    role: str,
    min_expr_frac: float,
    agg: str = "weighted_mean",
) -> pd.DataFrame:
    mets = sorted(enzyme_prior["metabolite"].dropna().astype(str).unique())
    out = pd.DataFrame(0.0, index=mets, columns=z_expr.index)
    role_prior = enzyme_prior[enzyme_prior["role"] == role]
    for metabolite, sub in role_prior.groupby("metabolite"):
        metabolite = str(metabolite)
        for cell_group in z_expr.index:
            genes = [g for g in sub["gene"].astype(str) if g in z_expr.columns and expr_frac.loc[cell_group, g] >= min_expr_frac]
            if not genes:
                continue
            ss = sub[sub["gene"].astype(str).isin(genes)]
            vals = z_expr.loc[cell_group, ss["gene"].astype(str)]
            weights = ss["weight"].astype(float)
            out.loc[metabolite, cell_group] = softmin(vals, weights) if agg == "softmin" else weighted_mean(vals, weights)
    return out


def specificity(score_matrix: pd.DataFrame, eps: float = 1e-8) -> pd.DataFrame:
    shifted = score_matrix - score_matrix.min(axis=1).values[:, None]
    bg = shifted.mean(axis=1).replace(0, np.nan)
    spec = np.log((shifted.add(eps)).div(bg.add(eps), axis=0))
    return spec.replace([np.inf, -np.inf], 0).fillna(0.0)
