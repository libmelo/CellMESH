"""Compiled numerical kernels for CELL MESH label permutations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy import sparse

from ._numerics import (
    _abundance_weights, _check_numeric_result, _event_score, _positive_product,
    _positive_reference_scores, _reaction_activity, _report_underflow, _sender_base,
)

from .config import MISSING_EXPORT_SCORE
from .database import _normalize_hmdb_id
from .preprocess import (
    _canonical_sparse_expression,
    _grouped_expression_mean,
    _slice_expression,
    _expression_source,
    _validate_expression_values,
    _validated_celltype_labels,
    _validated_gene_names,
    _validated_obs_labels,
)
from .score import (
    _build_prior_role_coverage,
    _get_reaction_gene_sets,
    _normalize_enzyme_metabolite,
    _production_evaluable_mask,
)


@dataclass(frozen=True)
class _GroupPlan:
    """Fixed cell groups and abundance weights for one scoring unit."""

    cell_indices: np.ndarray
    group_names: tuple[str, ...]
    counts: np.ndarray
    abundance_weights: np.ndarray


class CompiledPermutationScorer:
    """Score only observed event keys without rebuilding public result tables."""

    def __init__(
        self,
        adata,
        enzyme_prior: pd.DataFrame,
        sensor_prior: pd.DataFrame,
        observed_events: pd.DataFrame,
        *,
        cell_type_key: str,
        sample_key: Optional[str],
        sample_mode: str,
        layer: Optional[str],
        min_expr_frac: Optional[float],
        sender_abundance_exponent: float,
        pce_reference: str,
        export_weight: float,
        receiver_reference: str,
        prior_role_coverage: Optional[Dict[tuple[str, str], bool]],
    ) -> None:
        self.sample_mode = sample_mode
        self.min_expr_frac = min_expr_frac
        self.pce_reference = pce_reference
        self.export_weight = float(export_weight)
        self.receiver_reference = receiver_reference

        # As in observed scoring, deduplicate complete reaction gene sets before
        # selecting measured columns. Filtering genes first changes subset rules.
        parsed = _normalize_enzyme_metabolite(enzyme_prior)
        reaction_genes = _get_reaction_gene_sets(parsed)
        measured_genes = _validated_gene_names(adata)
        measured_set = set(measured_genes)
        sensor_genes = [
            str(gene)
            for gene in sensor_prior["sensor_gene"].astype(str).unique()
            if str(gene) in measured_set
        ]
        enzyme_genes = [
            str(gene)
            for genes in reaction_genes["genes"]
            for gene in genes
            if str(gene) in measured_set
        ]
        selected_genes = list(dict.fromkeys(enzyme_genes + sensor_genes))
        if not selected_genes:
            raise ValueError("No permutation-scoring genes are available")

        gene_positions = measured_genes.get_indexer(selected_genes)
        source = _expression_source(adata, layer)
        self.X = _slice_expression(source, columns=gene_positions)
        # Same raw-value rule as observed scoring; validate once at compilation,
        # never inside the repeated label-shuffling/scoring loop.
        _validate_expression_values(self.X, layer=layer)
        self.X = _canonical_sparse_expression(self.X)
        _validate_expression_values(self.X, layer=layer)
        self.gene_to_col = {gene: i for i, gene in enumerate(selected_genes)}

        self.reactions: list[tuple[int, str, np.ndarray]] = []
        self.metabolites: list[tuple[str, str]] = []
        self.hmdb_to_met: dict[str, int] = {}
        for _, row in reaction_genes.iterrows():
            hmdb = str(row["hmdb_id"])
            if hmdb not in self.hmdb_to_met:
                self.hmdb_to_met[hmdb] = len(self.metabolites)
                self.metabolites.append((str(row["metabolite"]), hmdb))
            cols = np.asarray(
                [self.gene_to_col[str(g)] for g in row["genes"] if str(g) in self.gene_to_col],
                dtype=int,
            )
            self.reactions.append(
                (self.hmdb_to_met[hmdb], str(row["direction"]), cols)
            )

        self.sensor_gene_to_col = {
            gene: self.gene_to_col[gene]
            for gene in sensor_genes
        }
        self.prior_role_coverage = (
            _build_prior_role_coverage(enzyme_prior, measured_genes)
            if prior_role_coverage is None
            else dict(prior_role_coverage)
        )
        self.production_evaluable = _production_evaluable_mask(
            (hmdb for _, hmdb in self.metabolites), self.prior_role_coverage,
        )

        self.event_met = np.asarray(
            [
                self.hmdb_to_met.get(str(_normalize_hmdb_id(value)), -1)
                for value in observed_events["hmdb_id"]
            ],
            dtype=int,
        )
        self.event_sensor_col = np.asarray(
            [self.sensor_gene_to_col.get(str(value), -1) for value in observed_events["sensor_gene"]],
            dtype=int,
        )
        self.event_sender = observed_events["sender"].astype(str).to_numpy()
        self.event_receiver = observed_events["receiver"].astype(str).to_numpy()
        self.n_events = len(observed_events)

        # 这里必须使用与观测评分相同的标准化标识，不能直接 astype(str)。
        # 若观测为 A 而置换仍为 " A "，事件匹配会失败并被当作零分，
        # 使本应为 1 的 p 值变为 1/(n_perms+1)。基因名和样本名同样共用
        # 标准化规则；仅在编译时处理一次，不在每次置换中重复清理。
        # 回归测试：tests/test_identifier_normalization.py。
        labels = _validated_celltype_labels(adata, cell_type_key).to_numpy()
        self.original_labels = labels
        samples = _validated_obs_labels(adata, sample_key).to_numpy() if sample_key is not None else None
        exponent = float(sender_abundance_exponent)
        if sample_mode == "sample_aware":
            if sample_key is None:
                raise ValueError("sample_aware permutation scoring requires sample_key")
            self.sample_indices = [
                np.flatnonzero(samples == sample)
                for sample in pd.unique(samples)
            ]
            self.group_plans = [
                self._make_group_plan(indices, labels, exponent)
                for indices in self.sample_indices
            ]
        else:
            all_indices = np.arange(len(labels), dtype=int)
            self.sample_indices = (
                [np.flatnonzero(samples == sample) for sample in pd.unique(samples)]
                if sample_key is not None
                else [all_indices]
            )
            self.group_plans = [self._make_group_plan(all_indices, labels, exponent)]

    @staticmethod
    def _make_group_plan(
        indices: np.ndarray,
        labels: np.ndarray,
        exponent: float,
    ) -> _GroupPlan:
        local_labels = labels[indices]
        counts_series = pd.Series(local_labels).value_counts()
        group_names = tuple(counts_series.index.astype(str))
        group_to_code = {name: i for i, name in enumerate(group_names)}
        codes = np.asarray([group_to_code[str(label)] for label in local_labels], dtype=int)
        counts = np.bincount(codes, minlength=len(group_names)).astype(float)
        fractions = counts / float(len(indices))
        abundance_weights = _abundance_weights(fractions, exponent, "permutation abundance weights")
        _check_numeric_result(abundance_weights, "permutation abundance weights")
        return _GroupPlan(
            cell_indices=np.asarray(indices, dtype=int),
            group_names=group_names,
            counts=counts,
            abundance_weights=abundance_weights,
        )

    def permute_codes(self, rng: np.random.Generator) -> np.ndarray:
        """Generate one label permutation in the legacy RNG order."""
        out = np.empty(len(self.original_labels), dtype=object)
        for indices in self.sample_indices:
            out[indices] = rng.permutation(self.original_labels[indices])
        return out.astype(str)

    @staticmethod
    def _normalize_rows(
        values: np.ndarray,
        method: str,
        *,
        evaluable: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Normalize eligible rows while retaining compiled metabolite indices."""
        scores = np.zeros_like(values, dtype=float)
        rows = range(values.shape[0]) if evaluable is None else np.flatnonzero(evaluable)
        for i in rows:
            scores[i], _ = _positive_reference_scores(
                values[i], method, "permutation P/C/E capacities",
            )
        return scores

    def _pseudobulk(
        self,
        plan: _GroupPlan,
        permuted_labels: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        local_labels = permuted_labels[plan.cell_indices]
        group_to_code = {name: i for i, name in enumerate(plan.group_names)}
        codes = np.asarray([group_to_code[str(label)] for label in local_labels], dtype=int)
        indicator = sparse.csr_matrix(
            (
                np.ones(len(codes), dtype=np.float64),
                (codes, np.arange(len(codes), dtype=int)),
            ),
            shape=(len(plan.group_names), len(codes)),
        )
        local_x = self.X[plan.cell_indices, :]
        pseudobulk = _grouped_expression_mean(local_x, indicator, plan.counts)
        # 原始值只需检查一次，但重新分组可能让求和溢出。因此每次置换的
        # 中间结果都必须检查，不能用观测评分通过校验来替代。
        _check_numeric_result(pseudobulk, "permutation pseudobulk expression")
        # Keep X > 0 aligned with observed expression fractions. Do not use
        # nnz/getnnz as a shortcut: stored zeros must not pass expression gates.
        positive = indicator @ (local_x > 0.0)
        positive = positive.toarray() if sparse.issparse(positive) else np.asarray(positive)
        expr_frac = np.asarray(positive, dtype=float) / plan.counts[:, None]
        _report_underflow((positive > 0) & (pseudobulk == 0), "permutation pseudobulk expression")
        return pseudobulk, expr_frac

    def _score_group(
        self,
        plan: _GroupPlan,
        permuted_labels: np.ndarray,
    ) -> tuple[np.ndarray, dict[int, np.ndarray], np.ndarray]:
        pseudobulk, expr_frac = self._pseudobulk(plan, permuted_labels)
        n_met = len(self.metabolites)
        n_groups = len(plan.group_names)
        capacities = {
            "product": np.zeros((n_met, n_groups), dtype=float),
            "substrate": np.zeros((n_met, n_groups), dtype=float),
            "exporter": np.zeros((n_met, n_groups), dtype=float),
        }
        for met_idx, direction, cols in self.reactions:
            if len(cols) == 0:
                continue
            activity = _reaction_activity(pseudobulk[:, cols], stage="permutation reaction activity")
            with np.errstate(over="ignore", invalid="ignore"):
                capacities[direction][met_idx] += _positive_product(
                    activity, plan.abundance_weights, "permutation abundance-adjusted reaction activity",
                )

        # Inf/NaN must be rejected before positive masks, where(), or clipping
        # can turn the failure into a plausible zero score.
        for direction, capacity in capacities.items():
            _check_numeric_result(capacity, f"permutation {direction} capacities")

        P = capacities["product"]
        C = capacities["substrate"]
        E = capacities["exporter"]
        # Gene availability is fixed during label shuffling. A computed all-zero
        # sample must enter the same median as observed scoring, not become NA.
        valid_met = self.production_evaluable
        # 必须在归一化前应用与观测相同的可计算性筛选；仅在提取观测事件时
        # 过滤太晚，会对已排除代谢物的 C/E 求无用参考值，甚至造成溢出。
        # 不能改为 P > 0；可计算的零值仍需参与。原始表达和容量校验仍在上游。
        # 保留全长数组，未计算行用 0 占位，并由 valid_met 排除，避免索引错位。
        # 回归测试：tests/test_permutation_evaluability.py。
        P_score = self._normalize_rows(P, self.pce_reference, evaluable=valid_met)
        C_score = self._normalize_rows(C, self.pce_reference, evaluable=valid_met)
        E_score = self._normalize_rows(E, self.pce_reference, evaluable=valid_met)
        base = _sender_base(P_score, C_score, "permutation sender normalization")
        E_effective = np.zeros_like(E_score)
        for met_idx, (_, hmdb) in enumerate(self.metabolites):
            key = (str(hmdb), "exporter")
            if key not in self.prior_role_coverage or not self.prior_role_coverage[key]:
                E_effective[met_idx] = MISSING_EXPORT_SCORE
            elif np.any(E[met_idx] > 0.0):
                E_effective[met_idx] = E_score[met_idx]
        factor = (1.0 - self.export_weight) + self.export_weight * E_effective
        availability = _positive_product(base, factor, "permutation sender scores")
        _check_numeric_result(availability, "permutation sender scores")
        availability = np.clip(availability, 0.0, 1.0)

        receiver_scores: dict[int, np.ndarray] = {}
        for gene_col in set(self.sensor_gene_to_col.values()):
            values = pseudobulk[:, gene_col]
            score, _ = _positive_reference_scores(
                values, self.receiver_reference, "permutation receiver pseudobulk expression",
            )
            if self.min_expr_frac is not None:
                score[expr_frac[:, gene_col] < self.min_expr_frac] = 0.0
            receiver_scores[gene_col] = score
        return availability, receiver_scores, valid_met

    def score(self, permuted_labels: np.ndarray) -> np.ndarray:
        """Return null scores aligned to the observed event rows."""
        if self.sample_mode == "sample_aware":
            sample_scores = np.full(
                (len(self.group_plans), self.n_events), np.nan, dtype=float
            )
            for sample_idx, plan in enumerate(self.group_plans):
                sample_scores[sample_idx] = self._score_events_for_group(
                    plan, permuted_labels, structural_missing=np.nan
                )
            result = np.zeros(self.n_events, dtype=float)
            computable = np.isfinite(sample_scores).any(axis=0)
            if computable.any():
                result[computable] = np.nanmedian(
                    sample_scores[:, computable], axis=0
                )
            return result
        return self._score_events_for_group(
            self.group_plans[0], permuted_labels, structural_missing=0.0
        )

    def _score_events_for_group(
        self,
        plan: _GroupPlan,
        permuted_labels: np.ndarray,
        *,
        structural_missing: float,
    ) -> np.ndarray:
        availability, receiver_scores, valid_met = self._score_group(
            plan, permuted_labels
        )
        group_to_idx = {name: i for i, name in enumerate(plan.group_names)}
        out = np.full(self.n_events, structural_missing, dtype=float)
        computed = np.zeros(self.n_events, dtype=bool)
        for event_idx in range(self.n_events):
            met_idx = self.event_met[event_idx]
            gene_col = self.event_sensor_col[event_idx]
            sender_idx = group_to_idx.get(self.event_sender[event_idx])
            receiver_idx = group_to_idx.get(self.event_receiver[event_idx])
            if (
                met_idx < 0
                or gene_col < 0
                or sender_idx is None
                or receiver_idx is None
                or not valid_met[met_idx]
                or gene_col not in receiver_scores
            ):
                continue
            out[event_idx] = _event_score(
                availability[met_idx, sender_idx], receiver_scores[gene_col][receiver_idx],
                "permutation event scores",
            )
            computed[event_idx] = True
        # Only deliberately uncomputable sample-event slots may contain NaN.
        # Check computed slots before sample nanmedian can skip numerical NaN.
        _check_numeric_result(out[computed], "permutation event scores")
        return out
