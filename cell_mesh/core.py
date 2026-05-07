from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
import pandas as pd

from .database import load_cell_mesh_database
from .preprocess import aggregate_expression, validate_priors
from .scoring import aggregate_role_score, sigmoid, specificity, zscore_by_gene

SensorType = Literal["surface_receptor", "transporter", "nuclear_receptor", "intracellular_sensor"]


@dataclass
class CellMeshResult:
    """Container returned by :func:`run_cell_mesh`."""

    events: pd.DataFrame
    sender_scores: pd.DataFrame
    receiver_scores: pd.DataFrame
    role_scores: dict[str, pd.DataFrame]
    parameters: dict

    def to_csv(self, prefix: str) -> None:
        self.events.to_csv(f"{prefix}.events.csv", index=False)
        self.sender_scores.to_csv(f"{prefix}.sender_scores.csv")
        self.receiver_scores.to_csv(f"{prefix}.receiver_scores.csv", index=False)


def _bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def _compute_cell_mesh_scores(
    mean_expr: pd.DataFrame,
    expr_frac: pd.DataFrame,
    enzyme_prior: pd.DataFrame,
    sensor_prior: pd.DataFrame,
    min_expr_frac: float,
    role_agg: str,
    alpha_prod: float,
    alpha_deg: float,
    alpha_export: float,
    alpha_specificity: float,
    beta_sensor: float,
    beta_specificity: float,
    beta_import: float,
    beta_usage: float,
    beta_compartment: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    z = zscore_by_gene(mean_expr)
    role_scores = {
        role: aggregate_role_score(z, expr_frac, enzyme_prior, role, min_expr_frac, agg=role_agg)
        for role in ["production", "degradation", "export", "import", "usage"]
    }
    prod = role_scores["production"]
    deg = role_scores["degradation"]
    export = role_scores["export"]
    sender_specificity = specificity(prod - deg + export)
    sender_raw = alpha_prod * prod - alpha_deg * deg + alpha_export * export + alpha_specificity * sender_specificity
    sender_scores = pd.DataFrame(sigmoid(sender_raw.values), index=sender_raw.index, columns=sender_raw.columns)

    sensor_rows = []
    for _, row in sensor_prior.iterrows():
        metabolite = str(row["metabolite"])
        sensor_gene = str(row["sensor_gene"])
        sensor_type = row["sensor_type"]
        prior_weight = float(row.get("weight", 1.0))
        if sensor_gene not in z.columns or metabolite not in sender_scores.index:
            continue
        sensor_expr = z[sensor_gene]
        sensor_frac = expr_frac[sensor_gene]
        receiver_specificity = sensor_expr - sensor_expr.mean()
        import_score = role_scores["import"].loc[metabolite] if metabolite in role_scores["import"].index else pd.Series(0.0, index=z.index)
        usage_score = role_scores["usage"].loc[metabolite] if metabolite in role_scores["usage"].index else pd.Series(0.0, index=z.index)

        for receiver in z.index:
            if sensor_frac.loc[receiver] < min_expr_frac:
                receiver_score = 0.0
            else:
                if sensor_type == "surface_receptor":
                    raw = beta_sensor * sensor_expr.loc[receiver] + beta_specificity * receiver_specificity.loc[receiver]
                elif sensor_type == "transporter":
                    raw = (
                        beta_sensor * sensor_expr.loc[receiver]
                        + beta_specificity * receiver_specificity.loc[receiver]
                        + beta_import * import_score.loc[receiver]
                        + beta_usage * usage_score.loc[receiver]
                    )
                else:
                    raw = (
                        beta_sensor * sensor_expr.loc[receiver]
                        + beta_specificity * receiver_specificity.loc[receiver]
                        + beta_import * import_score.loc[receiver]
                        + beta_compartment * usage_score.loc[receiver]
                    )
                receiver_score = float(sigmoid(raw))
            sensor_rows.append(
                {
                    "metabolite": metabolite,
                    "hmdb_id": row.get("hmdb_id", np.nan),
                    "sensor_gene": sensor_gene,
                    "sensor_type": sensor_type,
                    "receiver": receiver,
                    "receiver_score": receiver_score,
                    "prior_weight": prior_weight,
                    "sensor_expr_z": float(sensor_expr.loc[receiver]),
                    "sensor_expr_frac": float(sensor_frac.loc[receiver]),
                    "import_score": float(import_score.loc[receiver]),
                    "usage_score": float(usage_score.loc[receiver]),
                    "sensor_evidence_level": row.get("evidence_level", "unknown"),
                    "sensor_source": row.get("source", "unknown"),
                }
            )
    receiver_scores = pd.DataFrame(sensor_rows)
    return sender_scores, receiver_scores, role_scores


def _make_cell_mesh_events(sender_scores: pd.DataFrame, receiver_scores: pd.DataFrame, allow_self: bool) -> pd.DataFrame:
    rows = []
    for _, rr in receiver_scores.iterrows():
        metabolite = rr["metabolite"]
        if metabolite not in sender_scores.index:
            continue
        for sender, sender_score in sender_scores.loc[metabolite].items():
            receiver = rr["receiver"]
            if (not allow_self) and sender == receiver:
                continue
            event_score = float(sender_score) * float(rr["receiver_score"]) * float(rr["prior_weight"])
            rows.append(
                {
                    "sender": sender,
                    "receiver": receiver,
                    "metabolite": metabolite,
                    "hmdb_id": rr.get("hmdb_id", np.nan),
                    "sensor_gene": rr["sensor_gene"],
                    "sensor_type": rr["sensor_type"],
                    "sender_score": float(sender_score),
                    "receiver_score": float(rr["receiver_score"]),
                    "prior_weight": float(rr["prior_weight"]),
                    "cell_mesh_score": event_score,
                    "communication_score": event_score,
                    "sensor_expr_z": float(rr["sensor_expr_z"]),
                    "sensor_expr_frac": float(rr["sensor_expr_frac"]),
                    "receiver_import_score": float(rr["import_score"]),
                    "receiver_usage_score": float(rr["usage_score"]),
                    "sensor_evidence_level": rr.get("sensor_evidence_level", "unknown"),
                    "sensor_source": rr.get("sensor_source", "unknown"),
                }
            )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("cell_mesh_score", ascending=False).reset_index(drop=True)


def _permute_labels(labels: pd.Series, sample_labels: Optional[pd.Series], rng: np.random.Generator) -> pd.Series:
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


def _empirical_pvalues(
    obs_events: pd.DataFrame,
    adata,
    cell_type_key: str,
    sample_key: Optional[str],
    layer: Optional[str],
    enzyme_prior: pd.DataFrame,
    sensor_prior: pd.DataFrame,
    n_perms: int,
    random_state: int,
    min_cells_per_group: int,
    min_expr_frac: float,
    role_agg: str,
    allow_self: bool,
    score_kwargs: dict,
) -> pd.DataFrame:
    if n_perms <= 0 or obs_events.empty:
        obs_events["perm_pvalue"] = np.nan
        obs_events["fdr"] = np.nan
        return obs_events

    key_cols = ["sender", "receiver", "metabolite", "sensor_gene", "sensor_type"]
    obs_keys = obs_events[key_cols].astype(str).agg("|".join, axis=1)
    ge_counts = pd.Series(0, index=obs_keys.values, dtype=int)
    obs_score = pd.Series(obs_events["cell_mesh_score"].values, index=obs_keys.values)

    rng = np.random.default_rng(random_state)
    original = adata.obs[cell_type_key].copy()
    sample_labels = adata.obs[sample_key].copy() if sample_key is not None else None
    perm_key = "_cell_mesh_perm_label"
    try:
        for _ in range(n_perms):
            adata.obs[perm_key] = _permute_labels(original, sample_labels, rng).values
            agg = aggregate_expression(adata, groupby=perm_key, layer=layer, min_cells_per_group=min_cells_per_group)
            sender_perm, receiver_perm, _ = _compute_cell_mesh_scores(
                agg.mean_expr,
                agg.expr_frac,
                enzyme_prior,
                sensor_prior,
                min_expr_frac=min_expr_frac,
                role_agg=role_agg,
                **score_kwargs,
            )
            events_perm = _make_cell_mesh_events(sender_perm, receiver_perm, allow_self=allow_self)
            if events_perm.empty:
                continue
            perm_scores = events_perm.assign(_key=events_perm[key_cols].astype(str).agg("|".join, axis=1)).set_index("_key")["cell_mesh_score"]
            common = obs_score.index.intersection(perm_scores.index)
            ge_counts.loc[common] += (perm_scores.loc[common] >= obs_score.loc[common]).astype(int)
    finally:
        if perm_key in adata.obs:
            del adata.obs[perm_key]

    p = (ge_counts.loc[obs_keys.values].values + 1) / (n_perms + 1)
    out = obs_events.copy()
    out["perm_pvalue"] = p
    out["fdr"] = _bh_fdr(p)
    return out


def _confidence_tier(row: pd.Series) -> str:
    if pd.isna(row.get("fdr", np.nan)):
        if row["cell_mesh_score"] >= 0.5 and row["sensor_expr_frac"] >= 0.1:
            return "Tier2_no_permutation"
        return "Tier3_exploratory"
    if row["fdr"] <= 0.05 and row["cell_mesh_score"] >= 0.5 and row["sensor_expr_frac"] >= 0.1:
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
    min_cells_per_group: int = 5,
    min_expr_frac: float = 0.05,
    role_agg: Literal["weighted_mean", "softmin"] = "weighted_mean",
    allow_self: bool = True,
    n_perms: int = 0,
    random_state: int = 0,
    alpha_prod: float = 1.0,
    alpha_deg: float = 1.0,
    alpha_export: float = 0.5,
    alpha_specificity: float = 0.25,
    beta_sensor: float = 1.0,
    beta_specificity: float = 0.25,
    beta_import: float = 0.5,
    beta_usage: float = 0.25,
    beta_compartment: float = 0.25,
) -> CellMeshResult:
    """Run CELL MESH: Metabolite-mediated Event Scoring with Sensor Hierarchies.

    If `enzyme_metabolite` and `metabolite_sensor` are omitted, the packaged
    enzyme and metabolite-sensor interaction databases are used.
    """
    if sample_key is not None and sample_key not in adata.obs:
        raise KeyError(f"{sample_key!r} not found in adata.obs")
    if enzyme_metabolite is None or metabolite_sensor is None:
        default_enzyme, default_sensor = load_cell_mesh_database()
        enzyme_metabolite = default_enzyme if enzyme_metabolite is None else enzyme_metabolite
        metabolite_sensor = default_sensor if metabolite_sensor is None else metabolite_sensor

    enzyme_prior, sensor_prior = validate_priors(enzyme_metabolite, metabolite_sensor, adata.var_names)
    if enzyme_prior.empty:
        raise ValueError("No enzyme prior genes were found in adata.var_names")
    if sensor_prior.empty:
        raise ValueError("No sensor genes were found in adata.var_names")

    agg = aggregate_expression(adata, groupby=cell_type_key, layer=layer, min_cells_per_group=min_cells_per_group)
    score_kwargs = dict(
        alpha_prod=alpha_prod,
        alpha_deg=alpha_deg,
        alpha_export=alpha_export,
        alpha_specificity=alpha_specificity,
        beta_sensor=beta_sensor,
        beta_specificity=beta_specificity,
        beta_import=beta_import,
        beta_usage=beta_usage,
        beta_compartment=beta_compartment,
    )
    sender_scores, receiver_scores, role_scores = _compute_cell_mesh_scores(
        agg.mean_expr,
        agg.expr_frac,
        enzyme_prior,
        sensor_prior,
        min_expr_frac=min_expr_frac,
        role_agg=role_agg,
        **score_kwargs,
    )
    events = _make_cell_mesh_events(sender_scores, receiver_scores, allow_self=allow_self)
    events = _empirical_pvalues(
        events,
        adata=adata,
        cell_type_key=cell_type_key,
        sample_key=sample_key,
        layer=layer,
        enzyme_prior=enzyme_prior,
        sensor_prior=sensor_prior,
        n_perms=n_perms,
        random_state=random_state,
        min_cells_per_group=min_cells_per_group,
        min_expr_frac=min_expr_frac,
        role_agg=role_agg,
        allow_self=allow_self,
        score_kwargs=score_kwargs,
    )
    if not events.empty:
        events["confidence_tier"] = events.apply(_confidence_tier, axis=1)
        events = events.sort_values(["fdr", "cell_mesh_score"], ascending=[True, False], na_position="last").reset_index(drop=True)

    params = dict(
        method="CELL MESH",
        acronym="Metabolite-mediated Event Scoring with Sensor Hierarchies",
        cell_type_key=cell_type_key,
        sample_key=sample_key,
        layer=layer,
        min_cells_per_group=min_cells_per_group,
        min_expr_frac=min_expr_frac,
        role_agg=role_agg,
        allow_self=allow_self,
        n_perms=n_perms,
        random_state=random_state,
        **score_kwargs,
    )
    return CellMeshResult(events=events, sender_scores=sender_scores, receiver_scores=receiver_scores, role_scores=role_scores, parameters=params)


def run_metcomm(*args, **kwargs) -> CellMeshResult:
    """Backward-compatible alias for older prototype code. Prefer run_cell_mesh."""
    return run_cell_mesh(*args, **kwargs)
