"""Exclude unscorable references without weakening permutation numeric checks."""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import run_cell_mesh
from cellmesh.core import EVENT_KEY_COLUMNS
from cellmesh.permutation import CompiledPermutationScorer


EXCLUDED = "HMDB0000001"
SIGNAL = "HMDB0000002"
ZERO = "HMDB0000003"


def _case(production_state="prior_missing", competing_role="degradation"):
    # One cell per type prevents pseudobulk overflow. Only the unnecessary
    # reference of the excluded metabolite would overflow on these finite data.
    adata = AnnData(
        np.array([[1., 2., 1e308, 0.], [8., 9., 1e308, 0.],
                  [0., 1., 1e308, 0.], [3., 6., 1e308, 0.]]),
        obs=pd.DataFrame({"cell_type": ["A", "B", "C", "D"],
                          "sample": ["D1", "D1", "D2", "D2"]},
                         index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["PROD", "SENSOR", "HUGE", "ZERO"]),
    )
    rows = [
        {"metabolite": "excluded", "hmdb_id": EXCLUDED, "gene": "HUGE",
         "role": competing_role, "reaction": "competing"},
        {"metabolite": "signal", "hmdb_id": SIGNAL, "gene": "PROD",
         "role": "production", "reaction": "production"},
        {"metabolite": "zero", "hmdb_id": ZERO, "gene": "ZERO",
         "role": "production", "reaction": "zero_production"},
    ]
    if production_state == "prior_gene_unavailable":
        rows.append({"metabolite": "excluded", "hmdb_id": EXCLUDED, "gene": "ABSENT",
                     "role": "production", "reaction": "unmeasured_production"})
    enzyme = pd.DataFrame(rows)
    sensor = pd.DataFrame([
        {"metabolite": name, "hmdb_id": hmdb, "sensor_gene": "SENSOR", "sensor_type": "Transporter"}
        for name, hmdb in (("excluded", EXCLUDED), ("signal", SIGNAL), ("zero", ZERO))
    ])
    return adata, enzyme, sensor


def _options(sample_mode="pooled_stratified", method="mean"):
    return dict(sample_key="sample", sample_mode=sample_mode, min_cells=1,
                sender_abundance_exponent=0.0, pce_reference=method)


def _scorer(adata, enzyme, sensor, observed, sample_mode="pooled_stratified", method="mean"):
    return CompiledPermutationScorer(
        adata, enzyme, sensor, observed, cell_type_key="cell_type", sample_key="sample",
        sample_mode=sample_mode, layer=None, min_expr_frac=None,
        sender_abundance_exponent=0.0, pce_reference=method, export_weight=0.2,
        receiver_reference="median", prior_role_coverage=None,
    )


@pytest.mark.parametrize("production_state", ["prior_missing", "prior_gene_unavailable"])
@pytest.mark.parametrize("competing_role", ["degradation", "export"])
@pytest.mark.parametrize("method", ["mean", "median"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_unscorable_references_match_full_recomputation(
    production_state, competing_role, method, sample_mode,
):
    adata, enzyme, sensor = _case(production_state, competing_role)
    options = _options(sample_mode, method)
    observed = run_cell_mesh(adata, enzyme, sensor, **options)
    assert set(observed.events.hmdb_id) == {SIGNAL, ZERO}
    units = (observed.availability_results["availability_by_sample"] if sample_mode == "sample_aware"
             else {"pooled": observed.availability_results})
    for unit in units.values():
        diagnostics = unit["production_diagnostics"]
        assert diagnostics.at[("excluded", EXCLUDED), "production_status"] == production_state
        assert not diagnostics.at[("excluded", EXCLUDED), "production_evaluable"]

    scorer = _scorer(adata, enzyme, sensor, observed.events, sample_mode, method)
    assert not scorer.production_evaluable[scorer.hmdb_to_met[EXCLUDED]]
    assert scorer.production_evaluable[scorer.hmdb_to_met[ZERO]]
    keys = pd.MultiIndex.from_frame(observed.events[EVENT_KEY_COLUMNS])
    rng = np.random.default_rng(17)
    for labels in [scorer.original_labels, *[scorer.permute_codes(rng) for _ in range(5)]]:
        shuffled = adata.copy()
        shuffled.obs["cell_type"] = labels
        reference = run_cell_mesh(shuffled, enzyme, sensor, **options).events
        expected = reference.set_index(EVENT_KEY_COLUMNS).cell_mesh_score.reindex(keys)
        assert expected.notna().all()
        actual = scorer.score(labels)
        assert np.isfinite(actual).all() and (actual >= 0).all()
        np.testing.assert_allclose(actual, expected.to_numpy(), rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("store_null_scores", [False, True])
def test_excluded_extreme_capacity_does_not_change_inference(sample_mode, n_jobs, store_null_scores):
    adata, enzyme, sensor = _case()
    options = {**_options(sample_mode), "n_perms": 5, "random_state": 17,
               "n_jobs": n_jobs, "store_null_scores": store_null_scores}
    expected = run_cell_mesh(adata, enzyme.loc[enzyme.hmdb_id.ne(EXCLUDED)], sensor, **options)
    actual = run_cell_mesh(adata, enzyme, sensor, **options)
    pd.testing.assert_frame_equal(actual.events, expected.events, check_exact=True)
    zero_events = actual.events.loc[actual.events.hmdb_id.eq(ZERO)]
    assert not zero_events.empty
    assert zero_events.cell_mesh_score.eq(0.0).all()
    assert zero_events.perm_pvalue.eq(1.0).all()
    if sample_mode == "sample_aware":
        pd.testing.assert_frame_equal(actual.sample_events, expected.sample_events, check_exact=True)
        pd.testing.assert_frame_equal(actual.events.attrs["sample_aware_null_scores"],
                                      expected.events.attrs["sample_aware_null_scores"], check_exact=True)
        assert actual.sample_events.cell_mesh_score.isna().any()  # Missing endpoint samples stay NA.


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -1.0])
def test_excluded_metabolite_raw_expression_is_still_validated(bad_value):
    adata, enzyme, sensor = _case()
    observed = run_cell_mesh(adata, enzyme, sensor, **_options()).events
    adata.X[0, adata.var_names.get_loc("HUGE")] = bad_value
    with pytest.raises(ValueError, match="scoring expression.*finite and non-negative"):
        run_cell_mesh(adata, enzyme, sensor, **_options())
    with pytest.raises(ValueError, match="scoring expression.*finite and non-negative"):
        _scorer(adata, enzyme, sensor, observed)


@pytest.mark.parametrize("competing_role", ["degradation", "export"])
def test_excluded_metabolite_capacity_overflow_is_still_rejected(competing_role):
    adata, enzyme, sensor = _case(competing_role=competing_role)
    observed = run_cell_mesh(adata, enzyme, sensor, **_options()).events
    adata = AnnData(np.column_stack([adata.X, np.full(adata.n_obs, 1e308)]),
                    obs=adata.obs.copy(), var=pd.DataFrame(index=[*adata.var_names, "HUGE2"]))
    enzyme = pd.concat([enzyme, pd.DataFrame([
        {"metabolite": "excluded", "hmdb_id": EXCLUDED, "gene": "HUGE2",
         "role": competing_role, "reaction": "second_competing"},
    ])], ignore_index=True)
    with pytest.raises(ValueError, match="capacities.*finite and non-negative"):
        run_cell_mesh(adata, enzyme, sensor, **_options())
    scorer = _scorer(adata, enzyme, sensor, observed)
    with pytest.raises(ValueError, match="capacities.*finite and non-negative"):
        scorer.score(scorer.original_labels)


@pytest.mark.parametrize("method", ["mean", "median"])
def test_evaluable_zero_production_still_requires_valid_competing_reference(method):
    adata, enzyme, sensor = _case()
    observed = run_cell_mesh(adata, enzyme, sensor, **_options(method=method)).events
    enzyme = pd.concat([enzyme, pd.DataFrame([
        {"metabolite": "excluded", "hmdb_id": EXCLUDED, "gene": "ZERO",
         "role": "production", "reaction": "now_evaluable"},
    ])], ignore_index=True)
    with pytest.raises(ValueError, match="reference.*finite and non-negative"):
        run_cell_mesh(adata, enzyme, sensor, **_options(method=method))
    scorer = _scorer(adata, enzyme, sensor, observed, method=method)
    assert scorer.production_evaluable[scorer.hmdb_to_met[EXCLUDED]]
    with pytest.raises(ValueError, match="reference.*finite and non-negative"):
        scorer.score(scorer.original_labels)
