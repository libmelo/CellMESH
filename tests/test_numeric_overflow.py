import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import compute_metabolite_availability, run_cell_mesh
from cellmesh.permutation import CompiledPermutationScorer
from cellmesh.score import _score_PCE_by_reference, compute_sensor_scores


def _case():
    adata = AnnData(
        X=np.ones((4, 3)),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"], "sample": ["D1"] * 4},
                         index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["PROD", "SENSOR", "UNUSED"]),
    )
    enzyme = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["HMDB0000001"],
                           "gene": ["PROD"], "role": ["production"], "reaction": ["p"]})
    sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["HMDB0000001"],
                           "sensor_gene": ["SENSOR"], "sensor_type": ["Transporter"]})
    return adata, enzyme, sensor


def _scorer(adata, enzyme, sensor, events, *, sample_mode="pooled_stratified", layer=None,
            sender_abundance_exponent=1.0):
    return CompiledPermutationScorer(
        adata, enzyme, sensor, events,
        cell_type_key="cell_type", sample_key="sample", sample_mode=sample_mode,
        layer=layer, min_expr_frac=None, sender_abundance_exponent=sender_abundance_exponent,
        pce_reference="mean", export_weight=0.2, receiver_reference="median", prior_role_coverage=None,
    )


@pytest.mark.parametrize("gene", ["PROD", "SENSOR"])
@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("layer", [None, "lognorm"])
def test_finite_observed_data_can_overflow_after_permutation(gene, storage, sample_mode, layer):
    adata, enzyme, sensor = _case()
    matrix = adata.X.copy()
    matrix[:, adata.var_names.get_loc(gene)] = [1e308, 0, 1e308, 0]
    assert np.isfinite(matrix).all() and (matrix >= 0).all()
    if storage != "dense":
        matrix = getattr(sparse, f"{storage}_matrix")(matrix)
    if layer:
        adata.layers[layer] = matrix
    else:
        adata.X = matrix
    observed = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                             layer=layer, min_cells=1).events
    assert np.isfinite(observed["cell_mesh_score"]).all()
    scorer = _scorer(adata, enzyme, sensor, observed, sample_mode=sample_mode, layer=layer)
    # This shuffle puts both 1e308 values in A, overflowing only after shuffling.
    with pytest.raises(ValueError, match="pseudobulk.*finite and non-negative"):
        scorer.score(np.array(["A", "B", "A", "B"]))


@pytest.mark.parametrize("gene", ["PROD", "SENSOR"])
def test_observed_and_direct_scoring_reject_aggregation_overflow(gene):
    adata, enzyme, sensor = _case()
    adata.X[:, adata.var_names.get_loc(gene)] = [1e308, 1e308, 0, 0]
    with pytest.raises(ValueError, match="pseudobulk.*finite and non-negative"):
        run_cell_mesh(adata, enzyme, sensor, min_cells=1)
    with pytest.raises(ValueError, match="pseudobulk.*finite and non-negative"):
        if gene == "PROD":
            compute_metabolite_availability(adata, enzyme, min_cells=1)
        else:
            compute_sensor_scores(adata, sensor, min_cells=1)


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("sample_mode,store_null", [
    ("pooled_stratified", False), ("sample_aware", False), ("sample_aware", True),
])
def test_permutation_overflow_aborts_public_inference(n_jobs, sample_mode, store_null):
    adata, enzyme, sensor = _case()
    adata.X[:, 0] = [1e308, 0, 1e308, 0]
    with pytest.raises(ValueError, match="pseudobulk.*finite and non-negative"):
        run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                      n_perms=10, random_state=0, n_jobs=n_jobs, store_null_scores=store_null, min_cells=1)


@pytest.mark.parametrize("role", ["production", "degradation", "export"])
def test_reaction_capacity_overflow_is_not_masked_as_zero(role):
    adata, enzyme, sensor = _case()
    observed = run_cell_mesh(adata, enzyme, sensor, min_cells=1).events
    adata = adata[[0, 2]].copy()
    adata.var_names = ["PROD", "SENSOR", "SECOND"]
    adata.X[:, [0, 2]] = 1e308
    enzyme["role"] = role
    enzyme = pd.concat([enzyme, enzyme.assign(gene="SECOND", reaction="second")], ignore_index=True)
    scorer = _scorer(adata, enzyme, sensor, observed, sender_abundance_exponent=0.0)
    with pytest.raises(ValueError, match="capacities.*finite and non-negative"):
        scorer.score(scorer.original_labels)
    with pytest.raises(ValueError, match="capacities.*finite and non-negative"):
        compute_metabolite_availability(adata, enzyme, min_cells=1, sender_abundance_exponent=0.0)


@pytest.mark.parametrize("method", ["mean", "median"])
@pytest.mark.parametrize("stage,values", [
    ("reference", [1e308, 1e308]), ("denominator", [1e308, 0.0]),
])
def test_positive_reference_overflow_is_rejected_in_both_paths(method, stage, values):
    index = pd.MultiIndex.from_tuples([("M", "HMDB0000001")], names=["metabolite", "hmdb_id"])
    capacities = pd.DataFrame([values], index=index, columns=["A", "B"])
    zeros = capacities * 0.0
    with pytest.raises(ValueError, match=rf"{stage}.*finite and non-negative"):
        _score_PCE_by_reference(capacities, zeros, zeros, pce_reference=method)
    with pytest.raises(ValueError, match=rf"{stage}.*finite and non-negative"):
        CompiledPermutationScorer._normalize_rows(np.array([values]), method)
    adata, _, sensor = _case()
    adata = adata[[0, 2]].copy()
    adata.X[:, 1] = values
    with pytest.raises(ValueError, match=rf"{stage}.*finite and non-negative"):
        compute_sensor_scores(adata, sensor, min_cells=1, receiver_reference=method)


@pytest.mark.parametrize("gene", ["PROD", "SENSOR"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_permutation_denominator_can_overflow_with_finite_means(gene, sample_mode):
    adata, enzyme, sensor = _case()
    adata.obs["cell_type"] = ["A", "A", "A", "B"]
    adata.X[:, adata.var_names.get_loc(gene)] = [1.5e308, 0, 0, 5e307]
    observed = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                             min_cells=1, sender_abundance_exponent=0.0).events
    scorer = _scorer(adata, enzyme, sensor, observed, sample_mode=sample_mode,
                     sender_abundance_exponent=0.0)
    # Means remain finite after shuffling, but 1.5e308 + its reference does not.
    with pytest.raises(ValueError, match="denominator.*finite and non-negative"):
        scorer.score(np.array(["B", "A", "A", "A"]))


def test_sample_median_does_not_treat_numerical_nan_as_structural_missing(monkeypatch):
    adata, enzyme, sensor = _case()
    adata.obs["sample"] = ["D1", "D1", "D2", "D2"]
    observed = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode="sample_aware",
                             min_cells=1).events
    scorer = _scorer(adata, enzyme, sensor, observed, sample_mode="sample_aware")
    original = scorer._score_group

    def numerical_nan(plan, labels):
        availability, receiver_scores, valid_met = original(plan, labels)
        if plan.cell_indices[0] == 0:
            for values in receiver_scores.values():
                values[0] = np.nan
        return availability, receiver_scores, valid_met

    monkeypatch.setattr(scorer, "_score_group", numerical_nan)
    with pytest.raises(ValueError, match="permutation event scores.*finite and non-negative"):
        scorer.score(scorer.original_labels)


@pytest.mark.parametrize("bad_score", [np.nan, np.inf, -1.0])
@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("sample_mode,store_null", [
    ("pooled_stratified", False), ("sample_aware", False), ("sample_aware", True),
])
def test_invalid_null_scores_never_enter_tail_counts(monkeypatch, bad_score, n_jobs, sample_mode, store_null):
    adata, enzyme, sensor = _case()

    def invalid_scores(self, labels):
        return np.full(self.n_events, bad_score)

    monkeypatch.setattr(CompiledPermutationScorer, "score", invalid_scores)
    with pytest.raises(ValueError, match="permutation.*scores.*finite and non-negative"):
        run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                      n_perms=3, n_jobs=n_jobs, store_null_scores=store_null, min_cells=1)


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_large_representable_values_keep_unit_pvalues(storage, sample_mode):
    adata, enzyme, sensor = _case()
    adata.X *= 1e100
    if storage != "dense":
        adata.X = getattr(sparse, f"{storage}_matrix")(adata.X)
    result = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                           n_perms=9, n_jobs=2, store_null_scores=True, min_cells=1)
    np.testing.assert_allclose(result.events["cell_mesh_score"], np.sqrt(0.45 * 0.5), rtol=1e-12)
    np.testing.assert_array_equal(result.events[["perm_pvalue", "fdr_global", "fdr_sensor_type"]],
                                   np.ones((4, 3)))


def test_unrelated_gene_overflow_does_not_block_scoring():
    adata, enzyme, sensor = _case()
    reference = run_cell_mesh(adata, enzyme, sensor, min_cells=1, n_perms=3)
    adata.X[:, 2] = 1e308
    result = run_cell_mesh(adata, enzyme, sensor, min_cells=1, n_perms=3)
    pd.testing.assert_frame_equal(result.events, reference.events, check_exact=True)


def test_reaction_plot_rejects_capacity_overflow():
    pytest.importorskip("matplotlib")
    from cellmesh import plot_metabolite_secretion_violin

    adata, enzyme, sensor = _case()
    enzyme = pd.concat([enzyme, enzyme.assign(gene="UNUSED", reaction="second")], ignore_index=True)
    result = run_cell_mesh(adata, enzyme, sensor, min_cells=1)
    adata.X[:, [0, 2]] = 1e308
    with pytest.raises(ValueError, match="capacities.*finite and non-negative"):
        plot_metabolite_secretion_violin(result, adata, metabolite="M", hmdb_id="HMDB0000001")
