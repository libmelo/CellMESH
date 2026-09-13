import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import compute_metabolite_availability, run_cell_mesh
from cellmesh.core import EVENT_KEY_COLUMNS, _permute_labels
from cellmesh.database import validate_priors
from cellmesh.permutation import CompiledPermutationScorer


@pytest.fixture
def identifier_case():
    adata = AnnData(
        X=np.array([[8, 1, 3, 1], [6, 2, 1, 0], [1, 5, 0, 7], [0, 3, 2, 9],
                    [7, 2, 4, 0], [5, 1, 2, 1], [2, 4, 1, 8], [1, 6, 0, 6]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"] * 2,
                          "sample": ["D1"] * 4 + ["D2"] * 4},
                         index=[f"c{i}" for i in range(8)]),
        var=pd.DataFrame(index=["PROD", "CONS", "EXP", "SENSOR"]),
    )
    enzyme = pd.DataFrame({"metabolite": ["M"] * 3, "hmdb_id": ["HMDB0000001"] * 3,
                           "gene": ["PROD", "CONS", "EXP"],
                           "role": ["production", "degradation", "export"],
                           "reaction": ["p", "c", "e"]})
    sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["HMDB0000001"],
                           "sensor_gene": ["SENSOR"], "sensor_type": ["Transporter"]})
    return adata, enzyme, sensor


def _padded(adata, field="all"):
    out = adata.copy()
    for key in ("cell_type", "sample"):
        if field in (key, "all"):
            out.obs[key] = pd.Categorical([f" \t{value}\u3000" for value in out.obs[key]])
    if field in ("genes", "all"):
        out.var_names = [f" {value}\t" for value in out.var_names]
    return out


@pytest.mark.parametrize("sample_mode,sample_key", [
    ("pooled_stratified", None), ("pooled_stratified", "sample"), ("sample_aware", "sample"),
])
@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
def test_padded_constant_labels_have_unit_pvalues(identifier_case, sample_mode, sample_key, storage):
    adata, enzyme, sensor = identifier_case
    adata.X = np.full(adata.shape, 0.1, dtype=np.float32)
    if storage != "dense":
        adata.X = getattr(sparse, f"{storage}_matrix")(adata.X)
    data = _padded(adata, "cell_type")
    result = run_cell_mesh(data, enzyme, sensor, sample_key=sample_key, sample_mode=sample_mode,
                           min_cells=1, n_perms=99, random_state=3, n_jobs=2,
                           store_null_scores=True)
    assert len(result.events) == 4
    np.testing.assert_array_equal(
        result.events[["perm_pvalue", "fdr_global", "fdr_sensor_type"]], np.ones((4, 3)),
    )


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("field", ["cell_type", "sample", "genes", "all"])
@pytest.mark.parametrize("layer", [None, "expression"])
def test_identifier_whitespace_does_not_change_results_or_inputs(identifier_case, sample_mode, field, layer):
    adata, enzyme, sensor = identifier_case
    if layer:
        adata.layers[layer] = sparse.csr_matrix(adata.X.astype(np.float32))
        adata.X = np.zeros(adata.shape)
    data = _padded(adata, field)
    obs_before, var_before = data.obs.copy(deep=True), data.var.copy(deep=True)
    x_before = data.X.copy()
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, layer=layer, min_cells=1,
                  min_expr_frac=0.25, n_perms=5, random_state=21, store_null_scores=True)
    reference = run_cell_mesh(adata, enzyme, sensor, **kwargs)
    result = run_cell_mesh(data, enzyme, sensor, **kwargs)
    for name in ("events", "sender_scores", "receiver_scores", "celltype_qc", "sample_validation",
                 "sample_sender_scores", "sample_receiver_scores", "sample_events"):
        expected = getattr(reference, name)
        if expected is not None:
            pd.testing.assert_frame_equal(getattr(result, name), expected, check_exact=True)
    if sample_mode == "sample_aware":
        pd.testing.assert_frame_equal(result.events.attrs["sample_aware_null_scores"],
                                      reference.events.attrs["sample_aware_null_scores"], check_exact=True)
        availability = result.availability_results["availability_by_sample"]
        units = {key: data[data.obs["sample"].astype(str).str.strip().eq(key)].copy() for key in availability}
    else:
        availability = {"pooled": result.availability_results}
        units = {"pooled": data}
    for key, main in availability.items():
        direct = compute_metabolite_availability(units[key], enzyme, layer=layer, min_cells=1)
        for name in ("P", "C", "E", "availability", "metadata", "celltype_qc"):
            pd.testing.assert_frame_equal(main[name], direct[name], check_exact=True)
    pd.testing.assert_frame_equal(data.obs, obs_before)
    pd.testing.assert_frame_equal(data.var, var_before)
    np.testing.assert_array_equal(data.X, x_before)
    if layer:
        np.testing.assert_array_equal(data.layers[layer].data, adata.layers[layer].data)


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_compiled_and_legacy_permutations_use_normalized_identifiers(identifier_case, sample_mode):
    adata, enzyme, sensor = identifier_case
    data = _padded(adata)
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, min_cells=1, min_expr_frac=None)
    observed = run_cell_mesh(adata, enzyme, sensor, **kwargs).events
    scorer = CompiledPermutationScorer(
        data, enzyme, sensor, observed, cell_type_key="cell_type", sample_key="sample",
        sample_mode=sample_mode, layer=None, min_expr_frac=None,
        sender_abundance_exponent=1.0, pce_reference="mean", export_weight=0.2,
        receiver_reference="median", prior_role_coverage=None,
    )
    keys = pd.MultiIndex.from_frame(observed[EVENT_KEY_COLUMNS])
    compiled_rng, legacy_rng = np.random.default_rng(17), np.random.default_rng(17)
    for _ in range(3):
        labels = scorer.permute_codes(compiled_rng)
        legacy = _permute_labels(data.obs["cell_type"], data.obs["sample"], legacy_rng)
        np.testing.assert_array_equal(labels, legacy.to_numpy())
        permuted = data.copy()
        permuted.obs["cell_type"] = labels
        actual = run_cell_mesh(permuted, enzyme, sensor, **kwargs).events
        expected = actual.set_index(EVENT_KEY_COLUMNS)["cell_mesh_score"].reindex(keys, fill_value=0.0)
        np.testing.assert_allclose(scorer.score(labels), expected, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("field", ["cell_type", "sample", "genes"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_identifier_collisions_are_rejected(identifier_case, field, sample_mode):
    adata, enzyme, sensor = identifier_case
    if field == "genes":
        adata.var_names = ["PROD", " PROD ", "EXP", "SENSOR"]
        message = "var_names must be unique"
    else:
        adata.obs.loc[adata.obs.index[1], field] = f" {adata.obs[field].iloc[0]} "
        message = rf"{field}.*collid.*normalization"
    with pytest.raises(ValueError, match=message):
        run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode, min_cells=1)
    if field != "sample":
        with pytest.raises(ValueError, match=message):
            compute_metabolite_availability(adata, enzyme, min_cells=1)


@pytest.mark.parametrize("field", ["cell_type", "sample"])
@pytest.mark.parametrize("values", [[1, "1", 2, 2], [1, 1.0, "1.0", 2]])
def test_string_conversion_collisions_are_rejected(identifier_case, field, values):
    adata, enzyme, sensor = identifier_case
    adata.obs[field] = values * 2
    with pytest.raises(ValueError, match=rf"{field}.*collid.*normalization"):
        run_cell_mesh(adata, enzyme, sensor, sample_key="sample", min_cells=1)


def test_unused_categories_do_not_create_identifier_collisions(identifier_case):
    adata, enzyme, sensor = identifier_case
    data = _padded(adata)
    data.obs["cell_type"] = data.obs["cell_type"].cat.add_categories(["A", "B", "unused"])
    data.obs["sample"] = data.obs["sample"].cat.add_categories(["D1", "D2", "unused"])
    result = run_cell_mesh(data, enzyme, sensor, sample_key="sample", min_cells=1)
    reference = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", min_cells=1)
    pd.testing.assert_frame_equal(result.events, reference.events, check_exact=True)


def test_explicit_cell_fractions_use_normalized_labels(identifier_case):
    adata, enzyme, _ = identifier_case
    data = _padded(adata)
    direct = compute_metabolite_availability(data, enzyme, min_cells=1,
                                             cell_fractions=pd.Series([0.5, 0.5], index=[" A ", " B "]))
    reference = compute_metabolite_availability(adata, enzyme, min_cells=1)
    pd.testing.assert_frame_equal(direct["availability"], reference["availability"], check_exact=True)
    with pytest.raises(ValueError, match="one value per cell type"):
        compute_metabolite_availability(data, enzyme, min_cells=1,
                                         cell_fractions=pd.Series([0.3, 0.2, 0.5], index=["A", " A ", "B"]))


def test_standalone_prior_validation_normalizes_expression_gene_names(identifier_case):
    adata, enzyme, sensor = identifier_case
    reference = validate_priors(enzyme, sensor, adata.var_names)
    result = validate_priors(enzyme, sensor, _padded(adata).var_names)
    for expected, actual in zip(reference, result):
        pd.testing.assert_frame_equal(actual, expected)


@pytest.mark.parametrize("plot_kind", ["metabolite", "receptor"])
def test_expression_plots_use_the_same_identifiers(identifier_case, plot_kind):
    plt = pytest.importorskip("matplotlib.pyplot")
    from cellmesh.plotting import plot_metabolite_secretion_violin, plot_receptor_expression_violin

    adata, enzyme, sensor = identifier_case
    data = _padded(adata)
    result = run_cell_mesh(data, enzyme, sensor, min_cells=1)
    if plot_kind == "metabolite":
        plot = plot_metabolite_secretion_violin
        kwargs = dict(metabolite="M", hmdb_id="HMDB0000001", sender_labels=[" A ", " B "])
        reference_kwargs = dict(metabolite="M", hmdb_id="HMDB0000001", sender_labels=["A", "B"])
    else:
        plot = plot_receptor_expression_violin
        kwargs = dict(receptor_gene="SENSOR", hmdb_id="HMDB0000001", receiver_labels=[" A ", " B "])
        reference_kwargs = dict(receptor_gene="SENSOR", hmdb_id="HMDB0000001", receiver_labels=["A", "B"])
    reference = plot(result, adata, **reference_kwargs)
    try:
        actual = plot(result, data, **kwargs)
        try:
            pd.testing.assert_frame_equal(actual["plot_data"], reference["plot_data"])
            pd.testing.assert_frame_equal(actual["summary"], reference["summary"])
        finally:
            plt.close(actual["fig"])
    finally:
        plt.close(reference["fig"])


@pytest.mark.parametrize("production_gene", ["PROD;UNMEASURED", "UNMEASURED"])
def test_reaction_plot_matches_scoring_with_unmeasured_prior_genes(identifier_case, production_gene):
    plt = pytest.importorskip("matplotlib.pyplot")
    from cellmesh.plotting import plot_metabolite_secretion_violin

    adata, enzyme, sensor = identifier_case
    enzyme.loc[0, "gene"] = production_gene
    has_measured_gene = "PROD" in production_gene
    if not has_measured_gene:
        # A scored metabolite requires at least one measurable production
        # reaction. Keep one alongside the unavailable reaction under review.
        other_reaction = enzyme.iloc[[0]].assign(gene="PROD", reaction="p_measured")
        enzyme = pd.concat([enzyme, other_reaction], ignore_index=True)
    data = _padded(adata)
    result = run_cell_mesh(data, enzyme, sensor, min_cells=1)
    actual = plot_metabolite_secretion_violin(result, data, metabolite="M", hmdb_id="HMDB0000001")
    try:
        np.testing.assert_allclose(actual["plot_data"]["production_score"], adata.X[:, 0], rtol=1e-12)
        production = actual["reaction_metadata"].query("direction == 'product' and reaction == 'p'")
        assert production.iloc[0]["genes"] == (["PROD"] if has_measured_gene else [])
        definition = actual["reaction_definitions"].query("direction == 'product' and reaction == 'p'")
        assert "UNMEASURED" in definition.iloc[0]["genes"]
    finally:
        plt.close(actual["fig"])
