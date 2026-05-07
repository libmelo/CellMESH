from __future__ import annotations

from dataclasses import dataclass

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
