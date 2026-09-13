import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import compute_metabolite_availability, run_cell_mesh
from cellmesh.core import EVENT_KEY_COLUMNS
from cellmesh.permutation import CompiledPermutationScorer


HMDB = {name: f"HMDB{i:07d}" for i, name in enumerate(
    ("signal", "zero", "unavailable", "no_prior", "partial", "partial_zero"), start=1,
)}
EVALUABLE = {HMDB[name] for name in ("signal", "zero", "partial", "partial_zero")}


@pytest.fixture
def production_case():
    adata = AnnData(
        X=np.array([[1., 0., 2., 1.]] * 4 + [[0., 0., 2., 1.]] * 4),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"] * 2,
                          "sample": ["D1"] * 4 + ["D2"] * 4},
                         index=[f"c{i}" for i in range(8)]),
        var=pd.DataFrame(index=["PROD", "ZERO", "CONS", "SENSOR"]),
    )
    enzyme = pd.DataFrame([
        {"metabolite": name, "hmdb_id": HMDB[name], "gene": gene,
         "role": role, "reaction": name}
        for name, gene, role in (
            ("signal", "PROD", "production"),
            ("zero", "ZERO", "production"),
            ("unavailable", "UNMEASURED", "production"),
            ("no_prior", "CONS", "degradation"),
            ("partial", "PROD;UNMEASURED", "production"),
            ("partial_zero", "ZERO;UNMEASURED", "production"),
        )
    ])
    sensor = pd.DataFrame([
        {"metabolite": name, "hmdb_id": hmdb, "sensor_gene": "SENSOR", "sensor_type": "Transporter"}
        for name, hmdb in HMDB.items()
    ])
    return adata, enzyme, sensor


def _pair(frame, metabolite="signal"):
    return frame.loc[frame["metabolite"].eq(metabolite)
                     & frame["sender"].eq("A") & frame["receiver"].eq("B")]


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("layer", [None, "analysis"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_measured_zeros_are_evaluable_but_unavailable_production_is_excluded(
    production_case, storage, layer, sample_mode,
):
    adata, enzyme, sensor = production_case
    matrix = adata.X.copy()
    if storage != "dense":
        matrix = getattr(sparse, f"{storage}_matrix")(matrix)
    if layer is None:
        adata.X = matrix
    else:
        adata.layers[layer] = matrix
        adata.X = np.full(adata.shape, -1.0)  # The unselected layer must not matter.
    result = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                           layer=layer, min_cells=1)
    assert set(result.events["hmdb_id"]) == EVALUABLE
    units = (result.availability_results["availability_by_sample"] if sample_mode == "sample_aware"
             else {"pooled": result.availability_results})
    for sample, unit in units.items():
        metadata = unit["metadata"]
        diagnostics = unit["production_diagnostics"]
        assert set(diagnostics.index.get_level_values("hmdb_id")) == set(HMDB.values())
        pd.testing.assert_index_equal(metadata.index, unit["P"].index)
        pd.testing.assert_index_equal(metadata.index, unit["availability"].index)
        assert set(unit["P"].index.get_level_values("hmdb_id")) == EVALUABLE
        for name, expected in (("unavailable", "prior_gene_unavailable"),
                               ("no_prior", "prior_missing")):
            idx = (name, HMDB[name])
            assert diagnostics.at[idx, "production_status"] == expected
            assert not diagnostics.at[idx, "production_evaluable"]
            assert idx not in unit["P"].index
            assert idx not in unit["availability"].index
            assert unit["P"].reindex(diagnostics.index).loc[idx].isna().all()
        for name in ("zero", "partial_zero"):
            idx = (name, HMDB[name])
            assert metadata.at[idx, "production_status"] == "prior_no_expression"
            assert metadata.at[idx, "production_evaluable"]
            assert unit["P"].loc[idx].eq(0.0).all()
            assert unit["P_score"].loc[idx].eq(0.0).all()
            assert unit["availability"].loc[idx].eq(0.0).all()
            assert pd.isna(unit["P_ref"].loc[idx])  # No positive reference is needed for zero scores.
        expected_status = "prior_no_expression" if sample == "D2" else "supported"
        assert metadata.at[("signal", HMDB["signal"]), "production_status"] == expected_status
        # Partial measurement uses available genes; missing genes are not zeros
        # in the geometric mean, and complete reaction definitions are retained.
        pd.testing.assert_series_equal(unit["P"].loc[("signal", HMDB["signal"])],
                                       unit["P"].loc[("partial", HMDB["partial"])], check_names=False)
        partial = unit["reaction_genes"].loc[lambda df: df["hmdb_id"].eq(HMDB["partial"])]
        assert partial["genes"].tolist() == [["PROD", "UNMEASURED"]]
        assert metadata.at[("partial", HMDB["partial"]), "n_product_reactions"] == 1
    if sample_mode == "sample_aware":
        sample_pair = _pair(result.sample_events).set_index("sample")
        assert sample_pair.at["D1", "cell_mesh_score"] == pytest.approx(np.sqrt(0.45 * 0.5))
        assert sample_pair.at["D2", "cell_mesh_score"] == 0.0
        assert sample_pair["sender_n_cells"].eq(2).all()
        assert sample_pair["receiver_n_cells"].eq(2).all()
        pair = _pair(result.events).iloc[0]
        assert pair["cell_mesh_score"] == pytest.approx(np.sqrt(0.45 * 0.5) / 2)
        assert pair["metabolite_availability_median"] == pytest.approx(0.225)
        assert pair["n_samples_coobserved"] == 2
        assert pair["n_samples_positive"] == 1
        assert pair["event_prevalence"] == 0.5
        assert result.events["n_samples_coobserved"].eq(2).all()
        assert _pair(result.events, "zero")["event_prevalence"].eq(0.0).all()


@pytest.mark.parametrize("min_cells", [1, 3])
def test_missing_endpoint_stays_na_while_zero_samples_keep_counts_and_qc(production_case, min_cells):
    adata, enzyme, sensor = production_case
    extra_obs = pd.DataFrame({"sample": ["D3"] * 4, "cell_type": ["A", "A", "C", "C"]},
                             index=[f"d3_{i}" for i in range(4)])
    adata = AnnData(np.vstack([adata.X, [[0., 0., 2., 1.]] * 4]),
                    obs=pd.concat([adata.obs, extra_obs]), var=adata.var.copy())
    result = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode="sample_aware",
                           min_cells=min_cells)
    pair = _pair(result.sample_events).set_index("sample")
    assert pair.at["D2", "cell_mesh_score"] == 0.0
    assert pair.at["D2", "sender_n_cells"] == 2
    assert pair.at["D2", "receiver_n_cells"] == 2
    assert bool(pair.at["D2", "passes_min_cells"]) == (min_cells <= 2)
    assert pd.isna(pair.at["D3", "cell_mesh_score"])
    assert pd.isna(pair.at["D3", "passes_min_cells"])
    summary = _pair(result.events).iloc[0]
    assert summary["n_samples_coobserved"] == 2
    assert summary["n_samples_positive"] == 1
    assert summary["event_prevalence"] == 0.5
    assert summary["sender_n_cells"] == 4
    assert summary["receiver_n_cells"] == 4
    assert summary["min_cells_pass_prevalence"] == (1.0 if min_cells <= 2 else 0.0)
    # Missing cell types in sender summaries stay NA too.
    assert pd.isna(result.sample_sender_scores.at[("D3", "signal", HMDB["signal"]), "B"])


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("store_null_scores", [False, True])
@pytest.mark.parametrize("n_jobs", [1, 2])
def test_all_measured_zero_events_have_unit_permutation_pvalues(
    production_case, sample_mode, store_null_scores, n_jobs,
):
    adata, enzyme, sensor = production_case
    adata.X[:, 0] = 0.0
    result = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                           min_cells=1, n_perms=5, random_state=19, n_jobs=n_jobs,
                           store_null_scores=store_null_scores)
    assert set(result.events["hmdb_id"]) == EVALUABLE
    assert result.events["cell_mesh_score"].eq(0.0).all()
    for name in ("perm_pvalue", "fdr_global", "fdr_sensor_type"):
        assert result.events[name].eq(1.0).all()
    if sample_mode == "sample_aware":
        assert result.events["n_samples_coobserved"].eq(2).all()
        assert result.events["n_samples_positive"].eq(0).all()
        assert result.events["event_prevalence"].eq(0.0).all()
        null = result.events.attrs["sample_aware_null_scores"]
        assert null.shape == (len(result.events), 5 if store_null_scores else 0)
        if store_null_scores:
            assert null.eq(0.0).all().all()


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_zero_events_without_permutations_have_no_pvalue(production_case, sample_mode):
    adata, enzyme, sensor = production_case
    result = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                           min_cells=1, n_perms=0)
    zeros = result.events.loc[result.events["metabolite"].eq("zero")]
    assert not zeros.empty
    assert zeros["cell_mesh_score"].eq(0.0).all()
    assert zeros["perm_pvalue"].isna().all()


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_unavailable_production_keeps_diagnostics_without_creating_zero_events(production_case, sample_mode):
    adata, enzyme, sensor = production_case
    enzyme = enzyme.loc[enzyme["metabolite"].isin(["unavailable", "no_prior"])].copy()
    result = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                           min_cells=1, n_perms=3)
    assert result.events.empty
    assert result.sender_scores.empty
    units = (result.availability_results["availability_by_sample"] if sample_mode == "sample_aware"
             else {"pooled": result.availability_results})
    for unit in units.values():
        assert unit["P"].empty
        assert unit["metadata"].empty
        diagnostics = unit["production_diagnostics"]
        assert not diagnostics["production_evaluable"].any()
        assert diagnostics.at[("unavailable", HMDB["unavailable"]), "production_status"] == "prior_gene_unavailable"
        assert diagnostics.at[("no_prior", HMDB["no_prior"]), "production_status"] == "prior_missing"


@pytest.mark.parametrize("return_intermediates", [False, True])
def test_standalone_availability_preserves_zero_and_unavailable_states(production_case, return_intermediates):
    adata, enzyme, sensor = production_case
    direct = compute_metabolite_availability(adata, enzyme, min_cells=1,
                                              return_intermediates=return_intermediates)
    main = run_cell_mesh(adata, enzyme, sensor, min_cells=1).availability_results
    for name in ("availability", "metadata", "production_diagnostics"):
        pd.testing.assert_frame_equal(direct[name], main[name], check_exact=True)
    assert direct["availability"].loc[("zero", HMDB["zero"])].eq(0.0).all()
    assert not direct["production_diagnostics"].at[("unavailable", HMDB["unavailable"]), "production_evaluable"]


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("storage", ["dense", "csr"])
def test_zero_sample_permutations_match_full_recomputation(production_case, sample_mode, storage):
    adata, enzyme, sensor = production_case
    adata.X[:4, 0] = [8., 2., 1., 0.]
    adata.X[:, 3] = [1., 3., 7., 2., 4., 2., 8., 1.]
    if storage == "csr":
        adata.X = sparse.csr_matrix(adata.X)
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, min_cells=1,
                  receiver_reference="mean")
    observed = run_cell_mesh(adata, enzyme, sensor, **kwargs).events
    scorer = CompiledPermutationScorer(
        adata, enzyme, sensor, observed, cell_type_key="cell_type", sample_key="sample",
        sample_mode=sample_mode, layer=None, min_expr_frac=None,
        sender_abundance_exponent=1.0, pce_reference="mean", export_weight=0.2,
        receiver_reference="mean", prior_role_coverage=None,
    )
    keys = pd.MultiIndex.from_frame(observed[EVENT_KEY_COLUMNS])
    rng = np.random.default_rng(7)
    null_columns = []
    for labels in [scorer.original_labels, *[scorer.permute_codes(rng) for _ in range(5)]]:
        permuted = adata.copy()
        permuted.obs["cell_type"] = labels
        reference = run_cell_mesh(permuted, enzyme, sensor, **kwargs).events
        expected = reference.set_index(EVENT_KEY_COLUMNS)["cell_mesh_score"].reindex(keys)
        assert expected.notna().all()
        np.testing.assert_allclose(scorer.score(labels), expected.to_numpy(), rtol=1e-12, atol=1e-14)
        null_columns.append(expected.to_numpy())
    # Compare public p-values with independent full-recomputation tail counts.
    expected_null = np.column_stack(null_columns[1:])
    obs = observed["cell_mesh_score"].to_numpy()
    threshold = obs - 100 * np.finfo(float).eps * np.abs(obs)
    expected_p = (1 + (expected_null >= threshold[:, None]).sum(axis=1)) / 6
    serial = run_cell_mesh(adata, enzyme, sensor, **kwargs, n_perms=5, random_state=7,
                           store_null_scores=True, n_jobs=1)
    parallel = run_cell_mesh(adata, enzyme, sensor, **kwargs, n_perms=5, random_state=7,
                             store_null_scores=False, n_jobs=2)
    pd.testing.assert_frame_equal(serial.events, parallel.events, check_exact=True)
    np.testing.assert_array_equal(serial.events["perm_pvalue"].to_numpy(), expected_p)
    if sample_mode == "sample_aware":
        np.testing.assert_allclose(serial.events.attrs["sample_aware_null_scores"].to_numpy(),
                                   expected_null, rtol=1e-12, atol=1e-14)


def test_permutation_zero_sample_is_computed_not_structurally_missing(production_case):
    adata, enzyme, sensor = production_case
    observed = run_cell_mesh(adata, enzyme, sensor, sample_key="sample",
                             sample_mode="sample_aware", min_cells=1).events
    scorer = CompiledPermutationScorer(
        adata, enzyme, sensor, observed, cell_type_key="cell_type", sample_key="sample",
        sample_mode="sample_aware", layer=None, min_expr_frac=None,
        sender_abundance_exponent=1.0, pce_reference="mean", export_weight=0.2,
        receiver_reference="median", prior_role_coverage=None,
    )
    zero_plan = scorer.group_plans[1]
    scores = scorer._score_events_for_group(zero_plan, scorer.original_labels, structural_missing=np.nan)
    assert np.isfinite(scores).all()
    assert (scores == 0.0).all()
    assert scorer.production_evaluable[scorer.hmdb_to_met[HMDB["zero"]]]
    assert not scorer.production_evaluable[scorer.hmdb_to_met[HMDB["unavailable"]]]
    assert not scorer.production_evaluable[scorer.hmdb_to_met[HMDB["no_prior"]]]
