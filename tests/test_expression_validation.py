import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import compute_metabolite_availability, run_cell_mesh
from cellmesh.permutation import CompiledPermutationScorer
from cellmesh.preprocess import _build_celltype_pseudobulk, _compute_celltype_expr_frac
from cellmesh.score import compute_sensor_scores


ERROR = "scoring expression must be finite and non-negative"


@pytest.fixture
def expression_case():
    adata = AnnData(
        X=np.array([[8, 1, 3, 1, 0], [6, 2, 1, 0, 0], [1, 5, 0, 7, 0], [0, 3, 2, 9, 0],
                    [7, 2, 4, 0, 0], [5, 1, 2, 1, 0], [2, 4, 1, 8, 0], [1, 6, 0, 6, 0]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"] * 2,
                          "sample": ["D1"] * 4 + ["D2"] * 4},
                         index=[f"c{i}" for i in range(8)]),
        var=pd.DataFrame(index=["PROD", "CONS", "EXP", "SENSOR", "UNUSED"]),
    )
    enzyme = pd.DataFrame({"metabolite": ["M"] * 3, "hmdb_id": ["HMDB0000001"] * 3,
                           "gene": ["PROD", "CONS", "EXP"],
                           "role": ["production", "degradation", "export"],
                           "reaction": ["p", "c", "e"]})
    sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["HMDB0000001"],
                           "sensor_gene": ["SENSOR"], "sensor_type": ["Transporter"]})
    return adata, enzyme, sensor


def _stored(matrix, storage):
    if storage == "dense":
        return matrix.copy()
    # Store zeros explicitly to cover valid sparse input as well as bad values.
    rows, columns = np.indices(matrix.shape)
    out = sparse.coo_matrix((matrix.ravel(), (rows.ravel(), columns.ravel())), shape=matrix.shape)
    return out.asformat(storage)


def _compiled(adata, enzyme, sensor, observed):
    return CompiledPermutationScorer(
        adata, enzyme, sensor, observed,
        cell_type_key="cell_type", sample_key="sample", sample_mode="pooled_stratified",
        layer=None, min_expr_frac=None, sender_abundance_exponent=1.0,
        pce_reference="mean", export_weight=0.2, receiver_reference="median",
        prior_role_coverage=None,
    )


@pytest.mark.parametrize("gene", ["PROD", "SENSOR"])
@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("layer", [None, "lognorm"])
def test_negative_cells_cannot_hide_in_positive_means(expression_case, gene, storage, sample_mode, layer):
    adata, enzyme, sensor = expression_case
    matrix = adata.X.copy()
    matrix[:, adata.var_names.get_loc(gene)] = [-1, 3, 1, 1] * 2
    if layer:
        adata.layers[layer] = _stored(matrix, storage)
    else:
        adata.X = _stored(matrix, storage)
    errors = []
    for n_perms in (0, 3):
        with pytest.raises(ValueError, match=ERROR) as exc:
            run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                          layer=layer, min_cells=1, n_perms=n_perms)
        errors.append(str(exc.value))
    assert errors[0] == errors[1]
    assert ("adata.X" if layer is None else "adata.layers['lognorm']") in errors[0]


@pytest.mark.parametrize("entrypoint", ["main", "availability", "sensor"])
@pytest.mark.parametrize("bad_value,storage", [(np.nan, "csr"), (np.inf, "csc"), (-np.inf, "dense")])
def test_invalid_values_fail_before_aggregation(expression_case, monkeypatch, entrypoint, bad_value, storage):
    import cellmesh.score as score_module

    adata, enzyme, sensor = expression_case
    gene = "PROD" if entrypoint == "availability" else "SENSOR"
    adata.X[0, adata.var_names.get_loc(gene)] = bad_value
    adata.X = _stored(adata.X, storage)

    def unexpected_aggregation(*args, **kwargs):
        pytest.fail("Invalid expression reached pseudobulk aggregation")

    monkeypatch.setattr(score_module, "_build_celltype_pseudobulk", unexpected_aggregation)
    with pytest.raises(ValueError, match=ERROR):
        if entrypoint == "main":
            run_cell_mesh(adata, enzyme, sensor, min_cells=1)
        elif entrypoint == "availability":
            compute_metabolite_availability(adata, enzyme, min_cells=1)
        else:
            compute_sensor_scores(adata, sensor, min_cells=1)


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("reaction_layout", ["separate", "multi_gene"])
def test_positive_reaction_totals_cannot_hide_negative_genes(expression_case, storage, reaction_layout):
    adata, enzyme, sensor = expression_case
    adata.X[:, 0] = -0.5
    adata.X[:, 1] = 3.0
    enzyme = enzyme.iloc[:2].copy()
    enzyme["role"] = "production"
    if reaction_layout == "multi_gene":
        enzyme = enzyme.iloc[[0]].assign(gene="PROD;CONS")
    adata.X = _stored(adata.X, storage)
    # Both summing -0.5 + 3 and gmean(0.5, 4) - 1 yield positive P capacity.
    with pytest.raises(ValueError, match=ERROR):
        compute_metabolite_availability(adata, enzyme, min_cells=1)
    with pytest.raises(ValueError, match=ERROR):
        run_cell_mesh(adata, enzyme, sensor, min_cells=1)


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_only_used_genes_and_selected_layer_are_checked(expression_case, storage, sample_mode):
    adata, enzyme, sensor = expression_case
    reference = run_cell_mesh(adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                              min_cells=1, n_perms=5, random_state=12)
    matrix = adata.X.copy()
    matrix[:, -1] = [-1, np.nan, np.inf, -np.inf] * 2
    adata.layers["lognorm"] = _stored(matrix, storage)
    adata.X = np.full(adata.shape, np.nan)
    # The same helper must preserve canonical gene matching on an AnnData view.
    adata.var_names = [f" {gene} " for gene in adata.var_names]
    data = adata[:, :]
    obs_before, var_before = data.obs.copy(deep=True), data.var.copy(deep=True)
    result = run_cell_mesh(data, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
                          layer="lognorm", min_cells=1, n_perms=5, random_state=12)
    for name in ("events", "sender_scores", "receiver_scores", "celltype_qc"):
        pd.testing.assert_frame_equal(getattr(result, name), getattr(reference, name), check_exact=True)
    pd.testing.assert_frame_equal(data.obs, obs_before)
    pd.testing.assert_frame_equal(data.var, var_before)
    assert data.is_view
    np.testing.assert_array_equal(data.X, np.full(data.shape, np.nan))
    actual = data.layers["lognorm"]
    np.testing.assert_array_equal(actual.toarray() if sparse.issparse(actual) else actual, matrix)


@pytest.mark.parametrize("bad_value", [-1.0, np.nan, np.inf])
def test_precomputed_sensor_summaries_do_not_bypass_validation(expression_case, bad_value):
    adata, _, sensor = expression_case
    pseudobulk = _build_celltype_pseudobulk(adata)
    fractions = _compute_celltype_expr_frac(adata)
    adata.X[0, 3] = bad_value
    with pytest.raises(ValueError, match=ERROR):
        compute_sensor_scores(adata, sensor, pseudobulk=pseudobulk, expr_frac=fractions, min_cells=1)


@pytest.mark.parametrize("storage", ["dense", "csr", "csc"])
@pytest.mark.parametrize("bad_value", [-1.0, np.nan, np.inf])
def test_compiled_entrypoint_uses_the_same_validation(expression_case, storage, bad_value):
    adata, enzyme, sensor = expression_case
    observed = run_cell_mesh(adata, enzyme, sensor, min_cells=1).events
    adata.X[0, 0] = bad_value
    adata.X = _stored(adata.X, storage)
    with pytest.raises(ValueError, match=ERROR):
        _compiled(adata, enzyme, sensor, observed)


@pytest.mark.parametrize("plot_kind", ["metabolite", "receptor"])
@pytest.mark.parametrize("bad_value", [-1.0, np.nan, np.inf])
def test_expression_plots_reject_invalid_values(expression_case, plot_kind, bad_value):
    pytest.importorskip("matplotlib")
    from cellmesh import plot_metabolite_secretion_violin, plot_receptor_expression_violin

    adata, enzyme, sensor = expression_case
    result = run_cell_mesh(adata, enzyme, sensor, min_cells=1)
    adata.X[0, 0 if plot_kind == "metabolite" else 3] = bad_value
    with pytest.raises(ValueError, match=ERROR):
        if plot_kind == "metabolite":
            plot_metabolite_secretion_violin(result, adata, metabolite="M", hmdb_id="HMDB0000001")
        else:
            plot_receptor_expression_violin(result, adata, receptor_gene="SENSOR", hmdb_id="HMDB0000001")


@pytest.mark.parametrize("storage", ["csr", "csc"])
def test_sparse_validation_does_not_densify_or_modify_input(expression_case, monkeypatch, storage):
    from cellmesh.preprocess import _validate_scoring_expression

    adata, _, _ = expression_case
    adata.X = _stored(adata.X, storage)
    data, indices, indptr = adata.X.data.copy(), adata.X.indices.copy(), adata.X.indptr.copy()

    def unexpected_dense(*args, **kwargs):
        pytest.fail("Expression validation must not densify a sparse matrix")

    monkeypatch.setattr(type(adata.X), "toarray", unexpected_dense)
    _validate_scoring_expression(adata, ["PROD", "SENSOR"])
    np.testing.assert_array_equal(adata.X.data, data)
    np.testing.assert_array_equal(adata.X.indices, indices)
    np.testing.assert_array_equal(adata.X.indptr, indptr)
