import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh.core import (
    EVENT_KEY_COLUMNS,
    _compute_availability_scores,
    _compute_sample_aware_scores,
    _make_cell_mesh_events,
    _sample_cell_fractions,
    _sample_validation_table,
    run_cell_mesh,
)
from cellmesh.database import validate_priors
from cellmesh.permutation import CompiledPermutationScorer
from cellmesh.preprocess import _compute_celltype_fractions
from cellmesh.score import _build_prior_role_coverage


def _case():
    frame = pd.DataFrame(
        [
            ("S1", "A", 8.0, 1.0, 0.0, 1.0),
            ("S1", "A", 5.0, 0.0, 1.0, 2.0),
            ("S1", "B", 1.0, 5.0, 3.0, 8.0),
            ("S1", "B", 2.0, 4.0, 2.0, 6.0),
            ("S2", "A", 7.0, 1.0, 0.0, 2.0),
            ("S2", "A", 4.0, 2.0, 1.0, 1.0),
            ("S2", "B", 1.0, 7.0, 4.0, 9.0),
            ("S2", "B", 0.0, 5.0, 3.0, 7.0),
        ],
        columns=["sample", "cell_type", "PROD", "CONS", "EXPORT", "SENSOR"],
    )
    adata = AnnData(
        X=frame[["PROD", "CONS", "EXPORT", "SENSOR"]].to_numpy(),
        obs=frame[["sample", "cell_type"]].copy(),
        var=pd.DataFrame(index=["PROD", "CONS", "EXPORT", "SENSOR"]),
    )
    enzyme = pd.DataFrame(
        [
            {"metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": "PROD", "role": "production", "reaction": "prod"},
            {"metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": "CONS", "role": "degradation", "reaction": "cons"},
            {"metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": "EXPORT", "role": "export", "reaction": "export"},
        ]
    )
    sensor = pd.DataFrame(
        [{"metabolite": "Met", "hmdb_id": "HMDB0000001", "sensor_gene": "SENSOR", "sensor_type": "Transporter"}]
    )
    enzyme_prior, sensor_prior = validate_priors(enzyme, sensor, adata.var_names)
    coverage = _build_prior_role_coverage(enzyme, adata.var_names)
    kwargs = {
        "min_cells": 1,
        "sender_abundance_exponent": 1.0,
        "pce_reference": "mean",
        "export_weight": 0.2,
        "receiver_reference": "median",
        "_prior_role_coverage": coverage,
    }
    return adata, enzyme, sensor, enzyme_prior, sensor_prior, kwargs


def _aligned_scores(events, observed):
    if events.empty:
        return np.zeros(len(observed), dtype=float)
    keys = events[EVENT_KEY_COLUMNS].astype(str).agg("|".join, axis=1)
    observed_keys = observed[EVENT_KEY_COLUMNS].astype(str).agg("|".join, axis=1)
    values = pd.Series(events["cell_mesh_score"].to_numpy(), index=keys)
    return values.reindex(observed_keys).fillna(0.0).to_numpy()


@pytest.mark.parametrize(
    ("use_sparse", "pce_reference", "receiver_reference"),
    [(False, "mean", "median"), (True, "median", "mean")],
)
def test_compiled_pooled_scores_match_full_recomputation(
    use_sparse, pce_reference, receiver_reference
):
    adata, _, _, enzyme, sensor, kwargs = _case()
    if use_sparse:
        adata.X = sparse.csr_matrix(adata.X)
    kwargs["pce_reference"] = pce_reference
    kwargs["receiver_reference"] = receiver_reference
    sender, receiver, availability = _compute_availability_scores(
        adata, enzyme, sensor, celltype_col="cell_type", min_expr_frac=None, **kwargs
    )
    observed = _make_cell_mesh_events(
        sender, receiver, allow_self=True,
        cell_counts=availability["cell_counts"], min_cells=1,
    )
    scorer = CompiledPermutationScorer(
        adata, enzyme, sensor, observed,
        cell_type_key="cell_type", sample_key="sample",
        sample_mode="pooled_stratified", layer=None, min_expr_frac=None,
        sender_abundance_exponent=1.0, pce_reference=pce_reference,
        export_weight=0.2, receiver_reference=receiver_reference,
        prior_role_coverage=kwargs["_prior_role_coverage"],
    )
    rng = np.random.default_rng(7)
    fractions = _compute_celltype_fractions(adata, "cell_type")
    for _ in range(4):
        labels = scorer.permute_codes(rng)
        perm = adata.copy()
        perm.obs["perm"] = labels
        s, r, a = _compute_availability_scores(
            perm, enzyme, sensor, celltype_col="perm", min_expr_frac=None,
            cell_fractions=fractions, **kwargs,
        )
        events = _make_cell_mesh_events(
            s, r, allow_self=True, cell_counts=a["cell_counts"], min_cells=1
        )
        np.testing.assert_allclose(
            scorer.score(labels), _aligned_scores(events, observed), atol=1e-12
        )


@pytest.mark.parametrize(
    ("use_sparse", "pce_reference", "receiver_reference"),
    [(False, "mean", "median"), (True, "median", "mean")],
)
def test_compiled_sample_aware_scores_match_full_recomputation(
    use_sparse, pce_reference, receiver_reference
):
    adata, _, _, enzyme, sensor, kwargs = _case()
    if use_sparse:
        adata.X = sparse.csr_matrix(adata.X)
    kwargs["pce_reference"] = pce_reference
    kwargs["receiver_reference"] = receiver_reference
    _, _, observed, _ = _compute_sample_aware_scores(
        adata, enzyme, sensor, cell_type_key="cell_type", sample_key="sample",
        layer=None, min_expr_frac=None, allow_self=True,
        availability_kwargs=kwargs,
    )
    scorer = CompiledPermutationScorer(
        adata, enzyme, sensor, observed,
        cell_type_key="cell_type", sample_key="sample", sample_mode="sample_aware",
        layer=None, min_expr_frac=None, sender_abundance_exponent=1.0,
        pce_reference=pce_reference, export_weight=0.2,
        receiver_reference=receiver_reference,
        prior_role_coverage=kwargs["_prior_role_coverage"],
    )
    validation = _sample_validation_table(adata, "sample", "cell_type", 1)
    fractions = _sample_cell_fractions(validation)
    rng = np.random.default_rng(11)
    for _ in range(3):
        labels = scorer.permute_codes(rng)
        perm = adata.copy()
        perm.obs["perm"] = labels
        _, _, events, _ = _compute_sample_aware_scores(
            perm, enzyme, sensor, cell_type_key="perm", sample_key="sample",
            layer=None, min_expr_frac=None, allow_self=True,
            availability_kwargs=kwargs, cell_fractions_by_sample=fractions,
        )
        np.testing.assert_allclose(
            scorer.score(labels), _aligned_scores(events, observed), atol=1e-12
        )


def test_compiled_permutations_are_reproducible_across_workers_and_null_storage():
    adata, enzyme, sensor, _, _, _ = _case()
    common = dict(
        adata=adata, enzyme_metabolite=enzyme, metabolite_sensor=sensor,
        cell_type_key="cell_type", sample_key="sample",
        sample_mode="sample_aware", min_cells=1, n_perms=6, random_state=3,
    )
    serial = run_cell_mesh(**common, n_jobs=1, store_null_scores=True)
    parallel = run_cell_mesh(**common, n_jobs=2, store_null_scores=False)
    key = EVENT_KEY_COLUMNS
    serial_events = serial.events.sort_values(key).reset_index(drop=True)
    parallel_events = parallel.events.sort_values(key).reset_index(drop=True)
    np.testing.assert_allclose(
        serial_events[["perm_pvalue", "fdr_global", "fdr_sensor_type"]],
        parallel_events[["perm_pvalue", "fdr_global", "fdr_sensor_type"]],
    )
    assert serial.events.attrs["sample_aware_null_scores"].shape[1] == 6
    assert parallel.events.attrs["sample_aware_null_scores"].empty
    assert parallel.events.attrs["n_perms_completed"] == 6


@pytest.mark.parametrize(
    ("value", "error"),
    [(0, ValueError), (-2, ValueError), (1.5, TypeError), (True, TypeError)],
)
def test_run_cell_mesh_validates_n_jobs(value, error):
    adata, enzyme, sensor, _, _, _ = _case()
    with pytest.raises(error, match="n_jobs"):
        run_cell_mesh(
            adata, enzyme, sensor, cell_type_key="cell_type", n_jobs=value
        )


def test_run_cell_mesh_validates_store_null_scores():
    adata, enzyme, sensor, _, _, _ = _case()
    with pytest.raises(TypeError, match="store_null_scores"):
        run_cell_mesh(
            adata, enzyme, sensor, cell_type_key="cell_type",
            store_null_scores=1,
        )
