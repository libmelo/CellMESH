import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import cellmesh.core as core
from cellmesh.preprocess import _build_celltype_pseudobulk


def _constant_case(dtype=np.float64, storage="dense", counts=(3, 7)):
    n_cells = sum(counts)
    expression = np.tile(np.array([0.1, 0.3], dtype=dtype), (n_cells, 1))
    if storage == "csr":
        expression = sparse.csr_matrix(expression)
    elif storage == "csc":
        expression = sparse.csc_matrix(expression)
    elif storage == "fortran":
        expression = np.asfortranarray(expression)
    adata = AnnData(
        expression,
        obs=pd.DataFrame(
            {
                "cell_type": ["A"] * counts[0] + ["B"] * counts[1],
                "sample": ["S1"] * n_cells,
            },
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(index=["PROD", "SENSOR"]),
    )
    enzyme = pd.DataFrame(
        [{
            "metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": "PROD",
            "role": "production", "reaction": "prod",
        }]
    )
    sensor = pd.DataFrame(
        [{
            "metabolite": "Met", "hmdb_id": "HMDB0000001",
            "sensor_gene": "SENSOR", "sensor_type": "Transporter",
        }]
    )
    return adata, enzyme, sensor


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("storage", ["dense", "fortran", "csr", "csc"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_constant_expression_has_unit_permutation_pvalues(dtype, storage, sample_mode):
    adata, enzyme, sensor = _constant_case(dtype, storage)
    result = core.run_cell_mesh(
        adata, enzyme, sensor, sample_key="sample", sample_mode=sample_mode,
        min_cells=1, n_perms=99, random_state=0,
    )

    # Shuffling identical cells preserves both the means and group abundances.
    assert len(result.events) == 4
    assert (result.events["cell_mesh_score"] > 0.0).all()
    np.testing.assert_array_equal(
        result.events[["perm_pvalue", "fdr_global", "fdr_sensor_type"]].to_numpy(),
        np.ones((4, 3)),
    )

    pseudobulk = _build_celltype_pseudobulk(adata)
    assert pseudobulk.to_numpy().dtype == np.float64
    expected_means = np.tile(np.array([0.1, 0.3], dtype=dtype).astype(np.float64), (2, 1))
    np.testing.assert_allclose(pseudobulk, expected_means, rtol=1e-15, atol=0.0)
    assert adata.X.dtype == dtype


@pytest.mark.parametrize(("storage", "shared_gene"), [("csc", False), ("dense", True)])
def test_large_constant_groups_use_consistent_accumulation(storage, shared_gene):
    adata, enzyme, sensor = _constant_case(np.float64, storage, counts=(3001, 10007))
    if shared_gene:
        adata = adata[:, ["PROD"]].copy()
        sensor["sensor_gene"] = "PROD"
    # Exercise a selected layer independently of the contents of X.
    adata.layers["expression"] = adata.X.copy()
    adata.X = np.zeros(adata.shape, dtype=np.float64)
    result = core.run_cell_mesh(
        adata, enzyme, sensor, layer="expression", min_cells=1, n_perms=3,
    )
    assert len(result.events) == 4
    assert (result.events["cell_mesh_score"] > 0.0).all()
    np.testing.assert_array_equal(result.events["perm_pvalue"], np.ones(4))


@pytest.mark.parametrize(
    ("sample_mode", "store_null_scores"),
    [("pooled_stratified", False), ("sample_aware", False), ("sample_aware", True)],
)
def test_permutation_counts_distinguish_roundoff_ties_from_lower_scores(
    sample_mode, store_null_scores, monkeypatch
):
    adata, enzyme, sensor = _constant_case()
    common = dict(
        adata=adata, enzyme_metabolite=enzyme, metabolite_sensor=sensor,
        sample_key="sample", sample_mode=sample_mode, min_cells=1,
        store_null_scores=store_null_scores,
    )
    observed = core.run_cell_mesh(**common, n_perms=0).events
    scores_by_pair = observed.set_index(["sender", "receiver"])["cell_mesh_score"]

    def null_scores(scorer, *, n_perms, **kwargs):
        assert n_perms == 3
        values = np.array([
            scores_by_pair.loc[(sender, receiver)]
            for sender, receiver in zip(scorer.event_sender, scorer.event_receiver)
        ])
        yield np.nextafter(values, -np.inf)  # Roundoff-equivalent tie.
        yield values * (1.0 - 1e-8)  # A resolved decrease must not count.
        yield values * (1.0 + 1e-8)

    monkeypatch.setattr(core, "_compiled_permutation_scores", null_scores)
    result = core.run_cell_mesh(**common, n_perms=3)
    np.testing.assert_array_equal(result.events["perm_pvalue"], np.full(4, 0.75))
    if store_null_scores:
        assert result.events.attrs["sample_aware_null_scores"].shape == (4, 3)


def test_permutation_tolerance_preserves_small_positive_scores():
    observed = np.array([0.5, 0.5, 1e-30, 1e-30, 1e-30, 0.0])
    null_scores = np.array([
        np.nextafter(0.5, 0.0),
        0.5 * (1.0 - 1e-8),
        np.nextafter(1e-30, 0.0),
        1e-30 * (1.0 - 1e-8),
        0.0,
        0.0,
    ])
    np.testing.assert_array_equal(
        null_scores >= core._permutation_score_thresholds(observed),
        [True, False, True, False, False, True],
    )
