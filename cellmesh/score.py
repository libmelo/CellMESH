"""
Scoring algorithms for metabolite availability and sensor activity.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from ._numerics import _check_numeric_result, _positive_reference_scores, _reaction_activity

from .config import (
    METABOLITE_AVAILABILITY_DEFAULTS,
    MISSING_EXPORT_SCORE,
    MIN_EXPR_FRAC,
    PCE_REFERENCE_METHODS,
    RECEIVER_REFERENCE_METHODS,
    ROLE_TO_DIRECTION,
)
from .database import (
    _normalize_hmdb_id, _prior_metabolite_names, _valid_hmdb_mask,
    normalize_enzyme_database,
)
from .preprocess import (
    _all_celltype_counts,
    _build_celltype_pseudobulk,
    _compute_celltype_fractions,
    _compute_celltype_expr_frac,
    _normalized_label_series,
    _validate_scoring_expression,
    _validated_gene_names,
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
        fractions.index = pd.Index(_normalized_label_series(
            pd.Series(fractions.index, dtype=object), "cell_fractions index",
            reject_collisions=False,
        ))
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


def _validate_export_weight(value: Any) -> float:
    """Return a normalized finite exporter-modulation weight in [0, 1]."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise TypeError("export_weight must be a real number")
    weight = float(value)
    if not np.isfinite(weight) or not 0.0 <= weight <= 1.0:
        raise ValueError("export_weight must be between 0 and 1 inclusive")
    return weight


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


def _validate_boolean(value: Any, name: str) -> bool:
    """Do not interpret strings such as 'False' through Python truthiness."""
    if not isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a boolean")
    return bool(value)


def _summary_axis(values, name: str) -> pd.Index:
    labels = pd.Index(_normalized_label_series(
        pd.Series(values, dtype=object), name, reject_collisions=False,
    ))
    if labels.has_duplicates:
        duplicates = labels[labels.duplicated()].unique().tolist()
        raise ValueError(f"{name} must be unique after normalization: {duplicates[:10]}")
    return labels


def _summary_coverage(actual, required, name: str, *, exact: bool = True) -> None:
    missing = required.difference(actual)
    extra = actual.difference(required) if exact else pd.Index([])
    if len(missing) or len(extra):
        raise ValueError(
            f"{name} does not match current AnnData: "
            f"missing={missing.tolist()[:10]}, unexpected={extra.tolist()[:10]}"
        )


def _summary_numbers(frame: pd.DataFrame, name: str) -> np.ndarray:
    # Never discard imaginary parts, interpret booleans as counts, or coerce
    # invalid text to missing data. Numeric text may be converted explicitly.
    for column in frame.columns:
        values = frame[column]
        if (pd.api.types.is_bool_dtype(values.dtype)
                or pd.api.types.is_complex_dtype(values.dtype)
                or pd.api.types.is_datetime64_any_dtype(values.dtype)
                or pd.api.types.is_timedelta64_dtype(values.dtype)
                or (not pd.api.types.is_numeric_dtype(values.dtype) and any(
                    isinstance(v, (bool, np.bool_, complex, np.complexfloating))
                    for v in values
                ))):
            raise TypeError(f"{name} must contain real numeric values, not booleans, complex values or dates")
    try:
        return frame.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float, na_value=np.nan)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain real numeric values") from exc


def _checked_receiver_summary(value, *, name, counts, genes, required_genes):
    if not isinstance(value, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame")
    frame = value.copy(deep=False)
    frame.index = _summary_axis(value.index, f"{name} index")
    frame.columns = _summary_axis(value.columns, f"{name} columns")
    _summary_coverage(frame.index, counts.index, f"{name} cell types")
    _summary_coverage(genes, frame.columns, f"{name} genes", exact=False)
    _summary_coverage(frame.columns, required_genes, f"{name} required sensor genes", exact=False)
    # Use the actual observed groups and measured prior genes, never the cache
    # itself, to define coverage. Dropping B changes A's positive reference;
    # filling B with zero also changes the scientific meaning of the input.
    # Unrelated expression columns may contain overflow and are not scored.
    selected = frame.loc[counts.index, required_genes]
    values = _summary_numbers(selected, name)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError(f"{name} must be finite and non-negative")
    if name == "expr_frac" and np.any(values > 1):
        raise ValueError("expr_frac must be between 0 and 1 inclusive")
    return pd.DataFrame(values, index=counts.index, columns=required_genes)


def _checked_receiver_counts(value, observed: pd.Series) -> pd.Series:
    if not isinstance(value, pd.Series):
        raise TypeError("cell_counts must be a pandas Series")
    counts = value.copy(deep=False)
    counts.index = _summary_axis(value.index, "cell_counts index")
    _summary_coverage(counts.index, observed.index, "cell_counts cell types")
    values = _summary_numbers(counts.reindex(observed.index).to_frame(), "cell_counts")[:, 0]
    if (np.any(~np.isfinite(values)) or np.any(values <= 0)
            or np.any(values != np.floor(values))):
        raise ValueError("cell_counts must contain finite positive integer counts")
    if np.any(values != observed.to_numpy()):
        raise ValueError("cell_counts must match actual observed cell counts in current AnnData")
    # Only return integer counts after validation; int(1.9) must never hide an
    # invalid input. Counts come from obs, not a full expression recomputation.
    return observed.copy()


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

    Optional summaries must cover exactly the cell types observed in this
    AnnData (the current sample for sample-aware calls) and every measured
    prior sensor gene. Axes are normalized and aligned without mutating inputs.
    Unused categorical levels are not observed groups. Relevant means must be
    finite and non-negative; fractions must be in [0, 1]; cell counts must match
    obs. These checks do not prove that cached values came from this X/layer;
    callers remain responsible for cache provenance.
    """
    min_cells = _validate_min_cells(min_cells)
    min_expr_frac = _validate_min_expr_frac(min_expr_frac)
    receiver_reference = _validate_receiver_reference(receiver_reference)
    # Cached summaries cannot prove that the underlying cell values are valid.
    _validate_scoring_expression(adata, sensor_prior["sensor_gene"], layer)
    observed_counts = _all_celltype_counts(adata, celltype_col)
    genes = _validated_gene_names(adata)
    valid_genes = pd.Index([g for g in sensor_prior["sensor_gene"].unique() if g in genes])
    if pseudobulk is None:
        pseudobulk = _build_celltype_pseudobulk(adata, celltype_col, layer)
    if expr_frac is None:
        expr_frac = _compute_celltype_expr_frac(adata, celltype_col, layer)
    if cell_counts is None:
        cell_counts = observed_counts
    else:
        cell_counts = _checked_receiver_counts(cell_counts, observed_counts)
    pseudobulk = _checked_receiver_summary(
        pseudobulk, name="pseudobulk", counts=observed_counts,
        genes=genes, required_genes=valid_genes,
    )
    expr_frac = _checked_receiver_summary(
        expr_frac, name="expr_frac", counts=observed_counts,
        genes=genes, required_genes=valid_genes,
    )
    cell_fractions = _resolve_cell_fractions(
        adata,
        pseudobulk,
        celltype_col=celltype_col,
        cell_fractions=cell_fractions,
    )

    if len(valid_genes) == 0:
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
        scores, _ = _positive_reference_scores(
            expr_values, receiver_reference, "receiver pseudobulk expression",
        )
        sensor_gene_scores[gene] = pd.Series(scores, index=pseudobulk.index)

    display_names = _prior_metabolite_names(sensor_prior)
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
                    "metabolite": display_names.get(_normalize_hmdb_id(row["hmdb_id"]), row["hmdb_id"]),
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
    # 与数据库读取共用标准化规则，不能另设方向别名映射或优先使用原始
    # direction，否则直接调用 availability 与 run_cell_mesh 会得到不同的 E。
    # 回归测试：test_public_availability_entrypoints_share_prior_normalization。
    df = normalize_enzyme_database(enzyme_metabolite)
    df = df.loc[_valid_hmdb_mask(df["hmdb_id"])].copy()
    df["direction"] = df["role"].map(ROLE_TO_DIRECTION)

    return df[["metabolite", "hmdb_id", "reaction", "gene", "direction"]].reset_index(drop=True)


def _get_reaction_gene_sets(reaction_table: pd.DataFrame) -> pd.DataFrame:
    """
    Group normalized single-gene rows by canonical HMDB ID, reaction, and direction.

    The first non-empty metabolite name for an HMDB ID, or the ID itself, is
    retained only as display metadata. Gene symbols are deduplicated within each reaction. For
    the same HMDB ID and direction, only maximal gene sets contribute: exact
    duplicates are retained once, strict subsets are omitted, and partially
    overlapping sets that are not subsets of one another remain intact.
    These are complete prior gene sets, including unmeasured genes. Expression
    availability is considered only when computing the retained reactions.
    """
    df = reaction_table.copy()
    if df.empty:
        return pd.DataFrame(columns=["metabolite", "hmdb_id", "reaction", "direction", "genes"])

    canonical_names = _prior_metabolite_names(df)

    def get_reaction_id(row):
        return (str(row["hmdb_id"]), str(row["reaction"]), str(row["direction"]))

    df["_reaction_id"] = df.apply(get_reaction_id, axis=1)

    # 多基因字段已由 normalize_enzyme_database 在表达基因匹配前统一拆分。
    # 此处不要再独立解析分隔符或注释，否则主流程、直接评分与置换会分歧。
    reaction_dict = {}
    for _, row in df.iterrows():
        rid = row["_reaction_id"]
        gene = row["gene"]

        if rid not in reaction_dict:
            reaction_dict[rid] = {
                "metabolite": canonical_names[row["hmdb_id"]],
                "hmdb_id": row["hmdb_id"],
                "reaction": row["reaction"],
                "direction": row["direction"],
                "genes": [],
            }

        if gene not in reaction_dict[rid]["genes"]:
            reaction_dict[rid]["genes"].append(gene)

    candidates = []
    seen_gene_sets: Dict[tuple[str, str], set[frozenset[str]]] = {}
    for _, info in reaction_dict.items():
        aggregation_key = (str(info["hmdb_id"]), str(info["direction"]))
        gene_set = frozenset(info["genes"])
        seen = seen_gene_sets.setdefault(aggregation_key, set())
        if gene_set in seen:
            continue
        seen.add(gene_set)
        candidates.append((aggregation_key, gene_set, info))

    gene_sets_by_group = {
        key: sets
        for key, sets in seen_gene_sets.items()
    }
    result = []
    for aggregation_key, gene_set, info in candidates:
        if any(
            gene_set < other_gene_set
            for other_gene_set in gene_sets_by_group[aggregation_key]
        ):
            continue
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


def _production_evaluable_mask(hmdb_ids, prior_role_coverage) -> np.ndarray:
    """Select production priors with at least one measured reaction gene."""
    # 观测与置换必须共用“生成基因可用”条件，不能为提速改回 P > 0。
    # 基因已测得但 P 全零是可计算的零分，必须参与样本中位数和阳性比例；
    # 反应基因全部不可用时，内部容量的占位 0 不能当作测得的零表达。
    # 回归测试：tests/test_production_evaluability.py。
    return np.asarray(
        [bool(prior_role_coverage.get((str(hmdb), "product"), False)) for hmdb in hmdb_ids],
        dtype=bool,
    )


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
            # Numerical placeholder only. Production evaluability is checked
            # separately before this reaction can contribute a scored event.
            scores.append(pd.Series(0.0, index=pseudobulk.index))
        else:
            # 去重已使用完整基因集合；此处才选可测得基因。缺失基因不当作
            # 零表达基因，完全不可测得的反应容量为 0，并保留其先验证据状态。
            valid_genes = [row["genes"][i] for i in valid_gene_idx]
            expr = pseudobulk[valid_genes].values
            geo_mean = _reaction_activity(expr)
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

        with np.errstate(over="ignore", invalid="ignore"):
            if direction == "product":
                P.loc[met_idx] += reaction_scores.loc[rid]
            elif direction == "substrate":
                C.loc[met_idx] += reaction_scores.loc[rid]
            elif direction == "exporter":
                E.loc[met_idx] += reaction_scores.loc[rid]

    return {"P": P, "C": C, "E": E}


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
        return _positive_reference_scores(values, pce_reference, "P/C/E capacities")

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
    export_weight: float = METABOLITE_AVAILABILITY_DEFAULTS["export_weight"],
    _prior_role_coverage: Optional[Dict[tuple[str, str], bool]] = None,
) -> Dict[str, Any]:
    """
    Compute sender scores from continuously normalized P/C/E capacities.

    Enzyme DataFrames use the same schema normalization and role/direction
    conflict checks as :func:`run_cell_mesh`. Raw field aliases and standard
    role tables are supported; export, exporter, and the legacy transporter
    direction all denote export evidence.

    Reaction activity is multiplied by ``cell_fraction`` raised to
    ``sender_abundance_exponent`` before P/C/E construction. The default of
    1.0 preserves linear population-capacity scoring, while 0.0 removes the
    sender abundance effect. Each P/C/E direction is normalized against the
    arithmetic mean of its strictly positive cell-type capacities by default;
    ``pce_reference="median"`` selects the positive median instead. The formal
    base sender score is ``P_score ** 2 / (P_score + C_score)``. Normalized
    exporter evidence then has a bounded effect:
    ``availability = base * ((1 - export_weight) + export_weight * E_effective)``.
    ``E_effective`` is the normalized E score when exporter evidence is
    evaluable, 0 for an evaluable prior with no expression, and the fixed
    neutral value 0.5 when the prior is missing or its genes are unavailable.
    The default ``export_weight=0.2`` therefore constrains the E factor to
    [0.8, 1.0], while 0 exactly recovers the base formula.
    ``min_cells`` only annotates cell-count QC;
    all observed cell types contribute to pseudobulk, fractions, references,
    and scores. ``_prior_role_coverage`` is an internal
    sidecar used by :func:`run_cell_mesh` to distinguish a missing relation from
    one whose genes are unavailable. Production requires at least one measured
    gene; evaluable zero capacities retain zero scores in both inference modes.
    Metadata stays aligned to score matrices. A separate production_diagnostics
    table covers all enzyme-prior metabolites, including unscorable ones, with
    production_status and production_evaluable explaining their inclusion.
    """
    min_cells = _validate_min_cells(min_cells)
    return_intermediates = _validate_boolean(return_intermediates, "return_intermediates")
    sender_abundance_exponent = _validate_sender_abundance_exponent(
        sender_abundance_exponent
    )
    pce_reference = _validate_pce_reference(pce_reference)
    export_weight = _validate_export_weight(export_weight)

    # ``min_cells`` is a QC/reporting threshold only. All observed cell types
    # participate in pseudobulk construction, abundance adjustment, P/C/E
    # references, and formal scores.
    cell_counts = _all_celltype_counts(adata, celltype_col)
    parsed_reactions = _normalize_enzyme_metabolite(enzyme_metabolite)
    reaction_genes = _get_reaction_gene_sets(parsed_reactions)
    _validate_scoring_expression(
        adata, (gene for genes in reaction_genes["genes"] for gene in genes), layer,
    )
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

    prior_role_coverage = (
        _build_prior_role_coverage(enzyme_metabolite, pseudobulk.columns)
        if _prior_role_coverage is None
        else dict(_prior_role_coverage)
    )
    reaction_scores = _compute_reaction_scores(pseudobulk, reaction_genes)
    sender_abundance_weights = cell_fractions.pow(sender_abundance_exponent)
    sender_abundance_weights.name = "sender_abundance_weight"
    reaction_scores = reaction_scores.mul(sender_abundance_weights, axis="columns")
    _check_numeric_result(reaction_scores.to_numpy(), "abundance-adjusted reaction activity")
    PCE = _compute_PCE_matrices(reaction_scores, reaction_genes)
    P, C, E = PCE["P"], PCE["C"], PCE["E"]
    for name, matrix in PCE.items():
        values = matrix.to_numpy(dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError(
                f"{name} capacities must be finite and non-negative; "
                "check the selected expression layer and enzyme prior"
            )

    def prior_state(met_idx, direction, capacity):
        key = (str(met_idx[1]), direction)
        if key not in prior_role_coverage:
            return "prior_missing"
        if not prior_role_coverage[key]:
            return "prior_gene_unavailable"
        if bool((capacity.loc[met_idx] > 0.0).any()):
            return "supported"
        return "prior_no_expression"

    # Keep an audit record even when production is unscorable. Otherwise a
    # removed P row loses the distinction between absent prior and absent genes.
    metadata = pd.DataFrame(index=P.index)
    reaction_counts = reaction_genes.groupby(["hmdb_id", "direction"], sort=False).size().to_dict()
    for column, direction in [
        ("n_product_reactions", "product"),
        ("n_substrate_reactions", "substrate"),
        ("n_exporter_reactions", "exporter"),
    ]:
        metadata[column] = np.asarray(
            [reaction_counts.get((str(hmdb), direction), 0) for _, hmdb in P.index],
            dtype=int,
        )
    for column, direction, capacity in [
        ("consumption_status", "substrate", C),
        ("export_status", "exporter", E),
        ("production_status", "product", P),
    ]:
        metadata[column] = pd.Series(
            [prior_state(met_idx, direction, capacity) for met_idx in P.index],
            index=P.index, dtype=object,
        )
    metadata["production_evaluable"] = _production_evaluable_mask(
        P.index.get_level_values("hmdb_id"), prior_role_coverage,
    )
    production_diagnostics = metadata[
        ["n_product_reactions", "production_status", "production_evaluable"]
    ].copy()
    valid_mets = P.index[metadata["production_evaluable"]]
    # Existing callers use metadata indices to select P/C/E rows. Keep that
    # alignment while the separate diagnostics retain unscorable metabolites.
    metadata = metadata.loc[valid_mets].copy()
    if len(valid_mets) == 0:
        empty_index = P.index[:0]
        empty_matrix = pd.DataFrame(
            index=empty_index,
            columns=P.columns,
            dtype=float,
        )
        result = {
            "availability": empty_matrix.copy(),
            "metadata": metadata,
            "production_diagnostics": production_diagnostics,
            "pce_reference": pce_reference,
            "export_weight": export_weight,
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
                    "base_availability": empty_matrix.copy(),
                    "E_effective": empty_matrix.copy(),
                    "E_factor": empty_matrix.copy(),
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
    _check_numeric_result(denominator.to_numpy(), "sender normalization denominator")
    base_availability = (
        P_score.pow(2)
        .div(denominator.where(denominator > 0.0))
        .fillna(0.0)
        .clip(lower=0.0, upper=1.0)
    )

    E_effective = pd.DataFrame(0.0, index=P.index, columns=P.columns, dtype=float)
    for met_idx in P.index:
        state = metadata.at[met_idx, "export_status"]
        if state == "supported":
            E_effective.loc[met_idx] = E_score.loc[met_idx]
        elif state in {"prior_missing", "prior_gene_unavailable"}:
            E_effective.loc[met_idx] = MISSING_EXPORT_SCORE
        # ``prior_no_expression`` remains zero: the prior is evaluable and
        # supplies direct evidence that exporter expression is absent.
    E_factor = (1.0 - export_weight) + export_weight * E_effective
    availability = base_availability * E_factor
    _check_numeric_result(availability.to_numpy(), "sender scores")
    availability = availability.clip(lower=0.0, upper=1.0)

    result = {
        "availability": availability.astype(float),
        "metadata": metadata,
        "production_diagnostics": production_diagnostics,
        "pce_reference": pce_reference,
        "export_weight": export_weight,
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
                "base_availability": base_availability.astype(float),
                "E_effective": E_effective.astype(float),
                "E_factor": E_factor.astype(float),
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
