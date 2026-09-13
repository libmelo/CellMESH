import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import compute_metabolite_availability, run_cell_mesh
from cellmesh.core import EVENT_KEY_COLUMNS
from cellmesh.database import validate_priors
from cellmesh.permutation import CompiledPermutationScorer


@pytest.fixture
def missing_gene_case():
    adata = AnnData(
        X=np.array([[4, 9, 1], [4, 9, 2], [1, 4, 5], [1, 4, 7],
                    [8, 3, 2], [8, 3, 1], [2, 5, 7], [2, 5, 5]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"] * 2,
                          "sample": ["D1"] * 4 + ["D2"] * 4},
                         index=[f"c{i}" for i in range(8)]),
        var=pd.DataFrame(index=["G1", "G2", "S"]),
    )
    rows = []
    for role in ("production", "degradation", "export"):
        for reaction, genes in (("partial", ["MISSING", "G1"]),
                                ("pair", ["G1", "G2"]),
                                ("duplicate", ["G2", "G1"]),
                                ("subset", ["G1"])):
            for gene in genes:
                rows.append({"metabolite": "Canonical M" if gene == "MISSING" else "Alias M",
                             "hmdb_id": "HMDB0000001", "reaction": f"{role}_{reaction}",
                             "gene": gene, "role": role})
    enzyme = pd.DataFrame(rows)
    sensor = pd.DataFrame({"metabolite": ["Sensor M"], "hmdb_id": ["HMDB0000001"],
                           "sensor_gene": ["S"], "sensor_type": ["Transporter"]})
    return adata, enzyme, sensor


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("csv_input", [False, True])
@pytest.mark.parametrize("use_sparse", [False, True])
def test_complete_reaction_sets_survive_missing_genes(
    tmp_path, missing_gene_case, sample_mode, csv_input, use_sparse,
):
    adata, enzyme, sensor = missing_gene_case
    if use_sparse:
        adata.X = sparse.csr_matrix(adata.X)
    prior_input = enzyme
    if csv_input:
        prior_input = tmp_path / "enzyme.csv"
        enzyme.to_csv(prior_input, index=False)
    result = run_cell_mesh(adata, prior_input, sensor, sample_key="sample",
                           sample_mode=sample_mode, min_cells=1, n_perms=0)
    if sample_mode == "sample_aware":
        actual = result.availability_results["availability_by_sample"]
        units = {sample: adata[adata.obs["sample"].eq(sample)].copy() for sample in actual}
    else:
        actual = {"pooled": result.availability_results}
        units = {"pooled": adata}
    for unit, main in actual.items():
        data = units[unit]
        direct = compute_metabolite_availability(data, enzyme, min_cells=1)
        for name in ("P", "C", "E", "availability", "metadata", "reaction_genes"):
            pd.testing.assert_frame_equal(main[name], direct[name], check_exact=True)
        assert len(main["reaction_genes"]) == 6
        for name in ("n_product_reactions", "n_substrate_reactions", "n_exporter_reactions"):
            assert main["metadata"][name].tolist() == [2]
        matrix = data.X.toarray() if sparse.issparse(data.X) else data.X
        means = pd.DataFrame(matrix, columns=data.var_names).groupby(
            data.obs["cell_type"].to_numpy(), sort=False,
        ).mean()
        expected = (means["G1"] + np.sqrt((means["G1"] + 1) * (means["G2"] + 1)) - 1) * 0.5
        for name in ("P", "C", "E"):
            np.testing.assert_allclose(main[name].iloc[0], expected.reindex(main[name].columns), rtol=1e-12)


@pytest.mark.parametrize("role,status_column", [("degradation", "consumption_status"),
                                                 ("export", "export_status")])
@pytest.mark.parametrize("state", ["unavailable", "zero", "absent"])
def test_complete_prior_preserves_unavailable_reaction_metadata(missing_gene_case, role, status_column, state):
    adata, enzyme, sensor = missing_gene_case
    enzyme = enzyme.loc[enzyme["role"].eq("production")].copy()
    if state != "absent":
        extra = enzyme.iloc[[0]].copy()
        extra["role"] = role
        extra["reaction"] = "evidence"
        extra["gene"] = "UNMEASURED" if state == "unavailable" else "ZERO"
        enzyme = pd.concat([enzyme, extra], ignore_index=True)
    if state == "zero":
        adata = AnnData(X=np.column_stack([adata.X, np.zeros(len(adata))]),
                        obs=adata.obs.copy(), var=pd.DataFrame(index=[*adata.var_names, "ZERO"]))
    main = run_cell_mesh(adata, enzyme, sensor, min_cells=1).availability_results
    direct = compute_metabolite_availability(adata, enzyme, min_cells=1)
    for name in ("P", "C", "E", "E_effective", "availability", "metadata", "reaction_genes"):
        pd.testing.assert_frame_equal(main[name], direct[name], check_exact=True)
    expected_status = {"unavailable": "prior_gene_unavailable", "zero": "prior_no_expression",
                       "absent": "prior_missing"}[state]
    assert main["metadata"][status_column].tolist() == [expected_status]


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_missing_gene_permutations_match_full_recomputation(missing_gene_case, sample_mode):
    adata, enzyme, sensor = missing_gene_case
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, min_cells=1, receiver_reference="mean")
    observed = run_cell_mesh(adata, enzyme, sensor, **kwargs).events
    scorer = CompiledPermutationScorer(
        adata, enzyme, sensor, observed, cell_type_key="cell_type", sample_key="sample",
        sample_mode=sample_mode, layer=None, min_expr_frac=None,
        sender_abundance_exponent=1.0, pce_reference="mean", export_weight=0.2,
        receiver_reference="mean", prior_role_coverage=None,
    )
    keys = pd.MultiIndex.from_frame(observed[EVENT_KEY_COLUMNS])
    rng = np.random.default_rng(19)
    for labels in (scorer.original_labels, scorer.permute_codes(rng), scorer.permute_codes(rng)):
        permuted = adata.copy()
        permuted.obs["cell_type"] = labels
        reference = run_cell_mesh(permuted, enzyme, sensor, **kwargs).events
        scores = reference.set_index(EVENT_KEY_COLUMNS)["cell_mesh_score"].reindex(keys, fill_value=0.0)
        np.testing.assert_allclose(scorer.score(labels), scores, rtol=1e-12, atol=1e-14)
    serial = run_cell_mesh(adata, enzyme, sensor, **kwargs, n_perms=5, random_state=19,
                           n_jobs=1, store_null_scores=True)
    parallel = run_cell_mesh(adata, enzyme, sensor, **kwargs, n_perms=5, random_state=19,
                             n_jobs=2, store_null_scores=False)
    pd.testing.assert_frame_equal(serial.events, parallel.events, check_exact=True)


def test_validation_can_preserve_complete_enzyme_prior(missing_gene_case):
    adata, enzyme, sensor = missing_gene_case
    full, _ = validate_priors(enzyme, sensor, adata.var_names, filter_enzyme_genes=False)
    matched, _ = validate_priors(enzyme, sensor, adata.var_names)
    assert "MISSING" in full["gene"].values
    assert "MISSING" not in matched["gene"].values


def test_main_still_rejects_entirely_unmeasurable_enzyme_prior(missing_gene_case):
    adata, enzyme, sensor = missing_gene_case
    enzyme["gene"] = "MISSING"
    with pytest.raises(ValueError, match="No enzyme prior genes"):
        run_cell_mesh(adata, enzyme, sensor, min_cells=1)
