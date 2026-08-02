"""
Scoring algorithms for metabolite availability and sensor activity.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from scipy.stats import gmean

from .config import (
    METABOLITE_AVAILABILITY_DEFAULTS,
    MIN_EXPR_FRAC,
    PCE_REFERENCE_METHODS,
    RECEIVER_REFERENCE_METHODS,
    ROLE_TO_DIRECTION,
    VALID_ROLES,
)
from .database import _normalize_hmdb_series, _valid_hmdb_mask
from .preprocess import (
    _all_celltype_counts,
    _build_celltype_pseudobulk,
    _compute_celltype_fractions,
    _compute_celltype_expr_frac,
)


def _resolve_cell_fractions(
    adata,
    pseudobulk: pd.DataFrame,
    *,
    celltype_col: str,
    cell_fractions: Optional[pd.Series],
) -> pd.Series:
    """Validate or derive population fractions for all observed cell types."""
    if cell_fractions is None:
        fractions = _compute_celltype_fractions(adata, celltype_col)
    else:
        fractions = pd.Series(cell_fractions, dtype=float).copy()
        fractions.index = fractions.index.astype(str)
        if fractions.index.has_duplicates:
            raise ValueError("cell_fractions must contain one value per cell type")
        supplied_values = pd.to_numeric(fractions, errors="coerce").to_numpy(dtype=float)
        if np.any(~np.isfinite(supplied_values)) or np.any(supplied_values <= 0.0):
            raise ValueError("cell_fractions must be finite and strictly positive")
        if float(supplied_values.sum()) > 1.0 + 1e-9:
            raise ValueError("cell_fractions must sum to at most 1")

    required = pd.Index(pseudobulk.index).astype(str)
    missing = required.difference(fractions.index)
    if len(missing):
        raise ValueError(
            "cell_fractions is missing observed cell types: "
            + ", ".join(missing.tolist())
        )

    fractions = pd.to_numeric(fractions.reindex(required), errors="coerce")
    fractions.index = pseudobulk.index
    values = fractions.to_numpy(dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("cell_fractions must be finite and strictly positive")
    fractions.name = "cell_fraction"
    return fractions.astype(float)


def _validate_sender_abundance_exponent(value: Any) -> float:
    """Return a normalized finite, non-negative sender abundance exponent."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise TypeError("sender_abundance_exponent must be a real number")
    exponent = float(value)
    if not np.isfinite(exponent) or exponent < 0.0:
        raise ValueError(
            "sender_abundance_exponent must be a finite non-negative number"
        )
    return exponent


def _validate_pce_reference(value: Any) -> str:
    """Return a normalized P/C/E positive-value reference statistic."""
    if not isinstance(value, str):
        raise TypeError("pce_reference must be a string")
    method = value.strip().lower()
    if method not in PCE_REFERENCE_METHODS:
        raise ValueError("pce_reference must be one of {'mean', 'median'}")
    return method


def _validate_receiver_reference(value: Any) -> str:
    """Return a normalized receiver positive-expression reference statistic."""
    if not isinstance(value, str):
        raise TypeError("receiver_reference must be a string")
    method = value.strip().lower()
    if method not in RECEIVER_REFERENCE_METHODS:
        raise ValueError("receiver_reference must be one of {'mean', 'median'}")
    return method


def _validate_min_cells(value: Any) -> int:
    """Return a positive integer cell-count QC threshold."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError("min_cells must be an integer")
    threshold = int(value)
    if threshold < 1:
        raise ValueError("min_cells must be greater than or equal to 1")
    return threshold


def _validate_min_expr_frac(value: Any) -> Optional[float]:
    """Return ``None`` or a finite receiver-expression fraction in [0, 1]."""
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise TypeError("min_expr_frac must be None or a real number")
    fraction = float(value)
    if not np.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
        raise ValueError("min_expr_frac must be between 0 and 1 inclusive")
    return fraction


def compute_sensor_scores(
    adata,
    sensor_prior: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_expr_frac: Optional[float] = MIN_EXPR_FRAC,
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    pseudobulk: Optional[pd.DataFrame] = None,
    expr_frac: Optional[pd.DataFrame] = None,
    cell_counts: Optional[pd.Series] = None,
    cell_fractions: Optional[pd.Series] = None,
    receiver_reference: str = METABOLITE_AVAILABILITY_DEFAULTS[
        "receiver_reference"
    ],
) -> pd.DataFrame:
    """
    Compute receiver scores from per-cell-type mean sensor expression.

    ``cell_fractions`` is retained as receiver metadata but does not enter the
    receiver score or its reference. ``min_cells`` only annotates receiver QC;
    every observed cell type contributes to pseudobulk construction and the
    reference. For each sensor gene, the default
    reference is the unweighted median across strictly positive observed
    cell-type means; ``receiver_reference="mean"`` selects their unweighted
    arithmetic mean. Positive expression is normalized continuously as
    ``R / (R + R_ref)``. Zero expression remains zero.
    """
    min_cells = _validate_min_cells(min_cells)
    min_expr_frac = _validate_min_expr_frac(min_expr_frac)
    receiver_reference = _validate_receiver_reference(receiver_reference)
    if pseudobulk is None:
        pseudobulk = _build_celltype_pseudobulk(adata, celltype_col, layer)
    if expr_frac is None:
        expr_frac = _compute_celltype_expr_frac(adata, celltype_col, layer)
    if cell_counts is None:
        cell_counts = _all_celltype_counts(adata, celltype_col)
    cell_fractions = _resolve_cell_fractions(
        adata,
        pseudobulk,
        celltype_col=celltype_col,
        cell_fractions=cell_fractions,
    )

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
                "receiver_cell_fraction",
                "receiver_passes_min_cells",
            ]
        )

    sensor_gene_scores = {}
    for gene in valid_genes:
        expr_values = pseudobulk[gene].to_numpy(dtype=float)
        if np.any(~np.isfinite(expr_values)) or np.any(expr_values < 0.0):
            raise ValueError(
                "receiver pseudobulk expression must be finite and non-negative"
            )

        positive = expr_values > 0.0
        scores = np.zeros_like(expr_values, dtype=float)
        if positive.any():
            positive_values = expr_values[positive]
            if receiver_reference == "mean":
                reference = float(np.mean(positive_values))
            else:
                reference = float(np.median(positive_values))
            scores[positive] = positive_values / (positive_values + reference)
        sensor_gene_scores[gene] = pd.Series(scores, index=pseudobulk.index)

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
                    "receiver_cell_fraction": float(cell_fractions.loc[receiver]),
                    "receiver_passes_min_cells": bool(
                        int(cell_counts.loc[receiver]) >= min_cells
                    ),
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
        raise ValueError("enzyme_metabolite is missing required column: reaction")
    invalid_reaction = df["reaction"].isna() | df["reaction"].astype(str).str.strip().eq("")
    if invalid_reaction.any():
        raise ValueError("enzyme_metabolite reaction values must be non-empty")
    df["reaction"] = df["reaction"].astype(str).str.strip()
    df["hmdb_id"] = _normalize_hmdb_series(df["hmdb_id"])
    df = df[_valid_hmdb_mask(df["hmdb_id"])]
    df = df[df["direction"].isin({"product", "substrate", "exporter"})]

    return df[["metabolite", "hmdb_id", "reaction", "gene", "direction"]].reset_index(drop=True)


def _get_reaction_gene_sets(reaction_table: pd.DataFrame) -> pd.DataFrame:
    """
    Group genes by canonical HMDB ID, reaction, and direction.

    The first metabolite name observed for an HMDB ID is retained only as
    display metadata. Gene symbols are deduplicated within each reaction while
    preserving their first-seen order.
    """
    df = reaction_table.copy()
    if df.empty:
        return pd.DataFrame(columns=["metabolite", "hmdb_id", "reaction", "direction", "genes"])

    canonical_names = (
        df.drop_duplicates("hmdb_id", keep="first")
        .set_index("hmdb_id")["metabolite"]
        .to_dict()
    )

    def get_reaction_id(row):
        return (str(row["hmdb_id"]), str(row["reaction"]), str(row["direction"]))

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

        if rid not in reaction_dict:
            reaction_dict[rid] = {
                "metabolite": canonical_names[row["hmdb_id"]],
                "hmdb_id": row["hmdb_id"],
                "reaction": row["reaction"],
                "direction": row["direction"],
                "genes": [],
            }

        for gene in genes:
            if gene not in reaction_dict[rid]["genes"]:
                reaction_dict[rid]["genes"].append(gene)

    result = []
    for _, info in reaction_dict.items():
        result.append(
            {
                "metabolite": info["metabolite"],
                "hmdb_id": info["hmdb_id"],
                "reaction": info["reaction"],
                "direction": info["direction"],
                "genes": info["genes"],
            }
        )

    return pd.DataFrame(result)


def _build_prior_role_coverage(
    enzyme_metabolite: pd.DataFrame,
    measured_genes,
) -> Dict[tuple[str, str], bool]:
    """Map each canonical HMDB/direction prior to measured-gene availability."""
    parsed = _normalize_enzyme_metabolite(enzyme_metabolite)
    reactions = _get_reaction_gene_sets(parsed)
    genes = set(pd.Index(measured_genes).astype(str))
    coverage: Dict[tuple[str, str], bool] = {}
    for _, row in reactions.iterrows():
        key = (
            str(row["hmdb_id"]),
            str(row["direction"]),
        )
        usable = any(str(gene) in genes for gene in row["genes"])
        coverage[key] = coverage.get(key, False) or usable
    return coverage


def _compute_reaction_scores(
    pseudobulk: pd.DataFrame,
    reaction_genes: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute equal-weight multi-gene activity for each HMDB/reaction/direction.
    """
    scores = []
    reaction_ids = []

    for _, row in reaction_genes.iterrows():
        rid = (str(row["hmdb_id"]), str(row["reaction"]), str(row["direction"]))
        reaction_ids.append(rid)

        valid_gene_idx = [i for i, g in enumerate(row["genes"]) if g in pseudobulk.columns]

        if not valid_gene_idx:
            scores.append(pd.Series(0.0, index=pseudobulk.index))
        else:
            valid_genes = [row["genes"][i] for i in valid_gene_idx]
            expr = pseudobulk[valid_genes].values
            geo_mean = gmean(expr + 1, axis=1) - 1
            scores.append(pd.Series(geo_mean, index=pseudobulk.index))

    reaction_index = pd.MultiIndex.from_tuples(
        reaction_ids,
        names=["hmdb_id", "reaction", "direction"],
    )
    return pd.DataFrame(scores, index=reaction_index, columns=pseudobulk.index)


def _compute_PCE_matrices(
    reaction_scores: pd.DataFrame,
    reaction_genes: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    """
    Sum reactions into one P/C/E row per canonical HMDB identifier.
    """
    reaction_info = {}
    metabolites_dict = {}

    for _, row in reaction_genes.iterrows():
        hmdb_str = str(row["hmdb_id"])
        rid = (hmdb_str, str(row["reaction"]), str(row["direction"]))
        reaction_info[rid] = {"direction": row["direction"], "hmdb_id": row["hmdb_id"]}
        if hmdb_str not in metabolites_dict:
            metabolites_dict[hmdb_str] = (row["metabolite"], row["hmdb_id"])

    cell_types = reaction_scores.columns
    met_index = pd.MultiIndex.from_tuples(list(metabolites_dict.values()), names=["metabolite", "hmdb_id"])

    P = pd.DataFrame(0.0, index=met_index, columns=cell_types)
    C = pd.DataFrame(0.0, index=met_index, columns=cell_types)
    E = pd.DataFrame(0.0, index=met_index, columns=cell_types)

    for rid in reaction_scores.index:
        info = reaction_info[rid]
        direction = info["direction"]
        hmdb = info["hmdb_id"]
        met_idx = metabolites_dict[str(hmdb)]

        if direction == "product":
            P.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == "substrate":
            C.loc[met_idx] += reaction_scores.loc[rid]
        elif direction == "exporter":
            E.loc[met_idx] += reaction_scores.loc[rid]

    return {"P": P, "C": C, "E": E}


def _safe_hmdb_compare(row: pd.Series, hmdb: Optional[str]) -> bool:
    if pd.isna(row["hmdb_id"]) or pd.isna(hmdb):
        return False
    return str(row["hmdb_id"]) == str(hmdb)


def _score_PCE_by_reference(
    P: pd.DataFrame,
    C: pd.DataFrame,
    E: pd.DataFrame,
    pce_reference: str = METABOLITE_AVAILABILITY_DEFAULTS["pce_reference"],
) -> Dict[str, Any]:
    """Normalize P/C/E with a configurable positive-value reference.

    For each metabolite and direction, strictly positive capacities are scaled
    as ``x / (x + reference_positive)``. The reference is their arithmetic
    mean by default and can instead be their median. Zero capacity remains
    zero.
    """
    pce_reference = _validate_pce_reference(pce_reference)
    P_score = pd.DataFrame(0.0, index=P.index, columns=P.columns, dtype=float)
    C_score = pd.DataFrame(0.0, index=P.index, columns=P.columns, dtype=float)
    E_score = pd.DataFrame(0.0, index=P.index, columns=P.columns, dtype=float)
    P_ref = pd.Series(np.nan, index=P.index, dtype=float, name="P_ref")
    C_ref = pd.Series(np.nan, index=P.index, dtype=float, name="C_ref")
    E_ref = pd.Series(np.nan, index=P.index, dtype=float, name="E_ref")
    def normalize(values: np.ndarray) -> tuple[np.ndarray, float]:
        x = np.asarray(values, dtype=float)
        score = np.zeros_like(x, dtype=float)
        if np.any(~np.isfinite(x)):
            raise ValueError("P/C/E capacities must be finite and non-negative")
        if np.any(x < 0.0):
            raise ValueError("P/C/E capacities must be finite and non-negative")
        positive = x > 0.0
        if not positive.any():
            return score, np.nan

        if pce_reference == "mean":
            reference = float(np.mean(x[positive]))
        else:
            reference = float(np.median(x[positive]))
        denominator = x[positive] + reference
        score[positive] = x[positive] / denominator
        return score, reference

    for met_idx in P.index:
        p_score, p_ref = normalize(P.loc[met_idx].to_numpy())
        P_score.loc[met_idx] = p_score
        P_ref.loc[met_idx] = p_ref

        c_score, c_ref = normalize(C.loc[met_idx].to_numpy())
        C_score.loc[met_idx] = c_score
        C_ref.loc[met_idx] = c_ref

        e_score, e_ref = normalize(E.loc[met_idx].to_numpy())
        E_score.loc[met_idx] = e_score
        E_ref.loc[met_idx] = e_ref

    return {
        "P_score": P_score,
        "C_score": C_score,
        "E_score": E_score,
        "P_ref": P_ref,
        "C_ref": C_ref,
        "E_ref": E_ref,
    }


def compute_metabolite_availability(
    adata,
    enzyme_metabolite: pd.DataFrame,
    celltype_col: str = "cell_type",
    layer: Optional[str] = None,
    min_cells: int = METABOLITE_AVAILABILITY_DEFAULTS["min_cells"],
    return_intermediates: bool = True,
    cell_fractions: Optional[pd.Series] = None,
    sender_abundance_exponent: float = METABOLITE_AVAILABILITY_DEFAULTS[
        "sender_abundance_exponent"
    ],
    pce_reference: str = METABOLITE_AVAILABILITY_DEFAULTS["pce_reference"],
    _prior_role_coverage: Optional[Dict[tuple[str, str], bool]] = None,
) -> Dict[str, Any]:
    """
    Compute sender scores from continuously normalized P/C/E capacities.

    Reaction activity is multiplied by ``cell_fraction`` raised to
    ``sender_abundance_exponent`` before P/C/E construction. The default of
    1.0 preserves linear population-capacity scoring, while 0.0 removes the
    sender abundance effect. Each P/C/E direction is normalized against the
    arithmetic mean of its strictly positive cell-type capacities by default;
    ``pce_reference="median"`` selects the positive median instead. The formal
    sender score is
    ``P_score ** 2 / (P_score + C_score)``; E is retained as support metadata
    and does not enter that score. ``min_cells`` only annotates cell-count QC;
    all observed cell types contribute to pseudobulk, fractions, references,
    and scores. ``_prior_role_coverage`` is an internal
    sidecar used by :func:`run_cell_mesh` to distinguish a genuinely missing
    C/E relation from a relation whose genes are unavailable in the expression
    matrix; it never contributes to P/C/E calculation.
    """
    min_cells = _validate_min_cells(min_cells)
    sender_abundance_exponent = _validate_sender_abundance_exponent(
        sender_abundance_exponent
    )
    pce_reference = _validate_pce_reference(pce_reference)

    # ``min_cells`` is a QC/reporting threshold only. All observed cell types
    # participate in pseudobulk construction, abundance adjustment, P/C/E
    # references, and formal scores.
    cell_counts = _all_celltype_counts(adata, celltype_col)
    pseudobulk = _build_celltype_pseudobulk(adata, celltype_col=celltype_col, layer=layer)
    expr_frac = _compute_celltype_expr_frac(adata, celltype_col=celltype_col, layer=layer)
    cell_fractions = _resolve_cell_fractions(
        adata,
        pseudobulk,
        celltype_col=celltype_col,
        cell_fractions=cell_fractions,
    )
    celltype_qc = pd.DataFrame(
        {
            "n_cells": cell_counts.astype(int),
            "cell_fraction": cell_fractions.reindex(cell_counts.index).astype(float),
            "passes_min_cells": (cell_counts >= min_cells).astype(bool),
        },
        index=cell_counts.index,
    )
    celltype_qc.index.name = celltype_col

    parsed_reactions = _normalize_enzyme_metabolite(enzyme_metabolite)
    reaction_genes = _get_reaction_gene_sets(parsed_reactions)
    prior_role_coverage = (
        _build_prior_role_coverage(enzyme_metabolite, pseudobulk.columns)
        if _prior_role_coverage is None
        else dict(_prior_role_coverage)
    )
    reaction_scores = _compute_reaction_scores(pseudobulk, reaction_genes)
    sender_abundance_weights = cell_fractions.pow(sender_abundance_exponent)
    sender_abundance_weights.name = "sender_abundance_weight"
    reaction_scores = reaction_scores.mul(sender_abundance_weights, axis="columns")
    PCE = _compute_PCE_matrices(reaction_scores, reaction_genes)
    P, C, E = PCE["P"], PCE["C"], PCE["E"]
    for name, matrix in PCE.items():
        values = matrix.to_numpy(dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError(
                f"{name} capacities must be finite and non-negative; "
                "check the selected expression layer and enzyme prior"
            )

    valid_mets = P.index[P.sum(axis=1) > 0]
    if len(valid_mets) == 0:
        empty_index = P.index[:0]
        empty_matrix = pd.DataFrame(
            index=empty_index,
            columns=P.columns,
            dtype=float,
        )
        metadata = pd.DataFrame(index=empty_index)
        for column in [
            "n_product_reactions",
            "n_substrate_reactions",
            "n_exporter_reactions",
        ]:
            metadata[column] = pd.Series(index=empty_index, dtype=int)
        for column in ["consumption_status", "export_status"]:
            metadata[column] = pd.Series(index=empty_index, dtype=object)

        result = {
            "availability": empty_matrix.copy(),
            "metadata": metadata,
            "pce_reference": pce_reference,
            "celltype_qc": celltype_qc,
        }
        if return_intermediates:
            result.update(
                {
                    "P": empty_matrix.copy(),
                    "C": empty_matrix.copy(),
                    "E": empty_matrix.copy(),
                    "P_score": empty_matrix.copy(),
                    "C_score": empty_matrix.copy(),
                    "E_score": empty_matrix.copy(),
                    "P_ref": pd.Series(index=empty_index, dtype=float, name="P_ref"),
                    "C_ref": pd.Series(index=empty_index, dtype=float, name="C_ref"),
                    "E_ref": pd.Series(index=empty_index, dtype=float, name="E_ref"),
                    "pseudobulk": pseudobulk.astype(float),
                    "expr_frac": expr_frac.astype(float),
                    "cell_counts": cell_counts,
                    "cell_fractions": cell_fractions,
                    "sender_abundance_weights": sender_abundance_weights,
                    "sender_abundance_exponent": sender_abundance_exponent,
                    "reaction_genes": reaction_genes,
                }
            )
        return result

    P = P.loc[valid_mets]
    C = C.reindex(valid_mets, fill_value=0.0)
    E = E.reindex(valid_mets, fill_value=0.0)

    scored = _score_PCE_by_reference(
        P,
        C,
        E,
        pce_reference=pce_reference,
    )
    P_score = scored["P_score"]
    C_score = scored["C_score"]
    E_score = scored["E_score"]
    denominator = P_score + C_score
    availability = (
        P_score.pow(2)
        .div(denominator.where(denominator > 0.0))
        .fillna(0.0)
        .clip(lower=0.0, upper=1.0)
    )

    def prior_state(
        met_idx: tuple[Any, Any],
        direction: str,
        capacity: pd.DataFrame,
    ) -> str:
        _, hmdb = met_idx
        key = (str(hmdb), direction)
        has_prior = key in prior_role_coverage
        has_usable_gene = bool(prior_role_coverage.get(key, False))
        if not has_prior:
            return "prior_missing"
        if not has_usable_gene:
            return "prior_gene_unavailable"
        if bool((capacity.loc[met_idx] > 0.0).any()):
            return "supported"
        return "prior_no_expression"

    consumption_states = {
        met_idx: prior_state(met_idx, "substrate", C)
        for met_idx in P.index
    }
    export_states = {
        met_idx: prior_state(met_idx, "exporter", E)
        for met_idx in P.index
    }

    metadata = pd.DataFrame(index=P.index)
    metadata["n_product_reactions"] = [
        sum(_safe_hmdb_compare(row, hmdb) and row["direction"] == "product" for _, row in reaction_genes.iterrows())
        for _, hmdb in P.index
    ]
    metadata["n_substrate_reactions"] = [
        sum(_safe_hmdb_compare(row, hmdb) and row["direction"] == "substrate" for _, row in reaction_genes.iterrows())
        for _, hmdb in P.index
    ]
    metadata["n_exporter_reactions"] = [
        sum(_safe_hmdb_compare(row, hmdb) and row["direction"] == "exporter" for _, row in reaction_genes.iterrows())
        for _, hmdb in P.index
    ]
    metadata["consumption_status"] = [consumption_states[met_idx] for met_idx in P.index]
    metadata["export_status"] = [export_states[met_idx] for met_idx in P.index]

    result = {
        "availability": availability.astype(float),
        "metadata": metadata,
        "pce_reference": pce_reference,
        "celltype_qc": celltype_qc,
    }
    if return_intermediates:
        result.update(
            {
                "P": P.astype(float),
                "C": C.astype(float),
                "E": E.astype(float),
                "P_score": P_score.astype(float),
                "C_score": C_score.astype(float),
                "E_score": E_score.astype(float),
                "P_ref": scored["P_ref"].astype(float),
                "C_ref": scored["C_ref"].astype(float),
                "E_ref": scored["E_ref"].astype(float),
                "pseudobulk": pseudobulk.astype(float),
                "expr_frac": expr_frac.astype(float),
                "cell_counts": cell_counts,
                "cell_fractions": cell_fractions,
                "sender_abundance_weights": sender_abundance_weights,
                "sender_abundance_exponent": sender_abundance_exponent,
                "reaction_genes": reaction_genes,
            }
        )

    return result
