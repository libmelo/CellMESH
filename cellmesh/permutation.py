"""Compiled numerical kernels for CELL MESH label permutations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import gmean

from .config import MISSING_EXPORT_SCORE
from .database import _normalize_hmdb_id
from .score import (
    _build_prior_role_coverage,
    _get_reaction_gene_sets,
    _normalize_enzyme_metabolite,
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

        parsed = _normalize_enzyme_metabolite(enzyme_prior)
        reaction_genes = _get_reaction_gene_sets(parsed)
        measured_genes = pd.Index(adata.var_names).astype(str)
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
        source = adata.layers[layer] if layer is not None else adata.X
        self.X = source[:, gene_positions]
        self._validate_expression(self.X)
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

        labels = adata.obs[cell_type_key].astype(str).to_numpy()
        self.original_labels = labels
        exponent = float(sender_abundance_exponent)
        if sample_mode == "sample_aware":
            if sample_key is None:
                raise ValueError("sample_aware permutation scoring requires sample_key")
            samples = adata.obs[sample_key].astype(str).to_numpy()
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
                [np.flatnonzero(adata.obs[sample_key].astype(str).to_numpy() == sample)
                 for sample in pd.unique(adata.obs[sample_key].astype(str).to_numpy())]
                if sample_key is not None
                else [all_indices]
            )
            self.group_plans = [self._make_group_plan(all_indices, labels, exponent)]

    @staticmethod
    def _validate_expression(X: Any) -> None:
        values = X.data if sparse.issparse(X) else np.asarray(X)
        values = np.asarray(values, dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError(
                "permutation-scoring expression must be finite and non-negative"
            )

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
        return _GroupPlan(
            cell_indices=np.asarray(indices, dtype=int),
            group_names=group_names,
            counts=counts,
            abundance_weights=np.power(fractions, exponent),
        )

    def permute_codes(self, rng: np.random.Generator) -> np.ndarray:
        """Generate one label permutation in the legacy RNG order."""
        out = np.empty(len(self.original_labels), dtype=object)
        for indices in self.sample_indices:
            out[indices] = rng.permutation(self.original_labels[indices])
        return out.astype(str)

    @staticmethod
    def _normalize_rows(values: np.ndarray, method: str) -> np.ndarray:
        scores = np.zeros_like(values, dtype=float)
        for i in range(values.shape[0]):
            row = values[i]
            positive = row > 0.0
            if not positive.any():
                continue
            positive_values = row[positive]
            reference = (
                float(np.mean(positive_values))
                if method == "mean"
                else float(np.median(positive_values))
            )
            scores[i, positive] = positive_values / (positive_values + reference)
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
                np.ones(len(codes), dtype=float),
                (codes, np.arange(len(codes), dtype=int)),
            ),
            shape=(len(plan.group_names), len(codes)),
        )
        local_x = self.X[plan.cell_indices, :]
        sums = indicator @ local_x
        positive = indicator @ (local_x > 0.0)
        sums = sums.toarray() if sparse.issparse(sums) else np.asarray(sums)
        positive = positive.toarray() if sparse.issparse(positive) else np.asarray(positive)
        pseudobulk = np.asarray(sums, dtype=float) / plan.counts[:, None]
        expr_frac = np.asarray(positive, dtype=float) / plan.counts[:, None]
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
            activity = gmean(pseudobulk[:, cols] + 1.0, axis=1) - 1.0
            capacities[direction][met_idx] += activity * plan.abundance_weights

        P = capacities["product"]
        C = capacities["substrate"]
        E = capacities["exporter"]
        valid_met = P.sum(axis=1) > 0.0
        P_score = self._normalize_rows(P, self.pce_reference)
        C_score = self._normalize_rows(C, self.pce_reference)
        E_score = self._normalize_rows(E, self.pce_reference)
        denominator = P_score + C_score
        base = np.divide(
            P_score ** 2,
            denominator,
            out=np.zeros_like(P_score),
            where=denominator > 0.0,
        )
        E_effective = np.zeros_like(E_score)
        for met_idx, (_, hmdb) in enumerate(self.metabolites):
            key = (str(hmdb), "exporter")
            if key not in self.prior_role_coverage or not self.prior_role_coverage[key]:
                E_effective[met_idx] = MISSING_EXPORT_SCORE
            elif np.any(E[met_idx] > 0.0):
                E_effective[met_idx] = E_score[met_idx]
        factor = (1.0 - self.export_weight) + self.export_weight * E_effective
        availability = np.clip(base * factor, 0.0, 1.0)

        receiver_scores: dict[int, np.ndarray] = {}
        for gene_col in set(self.sensor_gene_to_col.values()):
            values = pseudobulk[:, gene_col]
            score = np.zeros(n_groups, dtype=float)
            positive_mask = values > 0.0
            if positive_mask.any():
                positive_values = values[positive_mask]
                reference = (
                    float(np.mean(positive_values))
                    if self.receiver_reference == "mean"
                    else float(np.median(positive_values))
                )
                score[positive_mask] = positive_values / (
                    positive_values + reference
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
            out[event_idx] = np.sqrt(
                availability[met_idx, sender_idx]
                * receiver_scores[gene_col][receiver_idx]
            )
        return out
