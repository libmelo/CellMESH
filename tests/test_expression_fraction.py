import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from cellmesh import run_cell_mesh
from cellmesh.preprocess import _compute_celltype_expr_frac


def _case(storage, layer=None):
    # A expresses SENSOR in half its cells; B passes a 0.75 gate in both
    # samples (3/4 in S1, 4/4 in S2). ZERO and ALL cover the endpoints.
    values = np.column_stack([
        [9, 10, 11, 12, 1, 2, 3, 4] * 2,
        [1, 0, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 3, 3, 3, 3],
        np.zeros(16),
        np.ones(16),
    ]).astype(np.float64)
    if storage == "dense":
        expression = values.copy()
    elif storage.endswith("_explicit_zeros"):
        n_cells, n_genes = values.shape
        # Store EVERY entry, including zeros. A normal csr_matrix(values)
        # drops zeros and therefore cannot reproduce the getnnz regression.
        expression = sparse.csr_matrix(
            (
                values.ravel().copy(),
                np.tile(np.arange(n_genes), n_cells),
                np.arange(0, values.size + 1, n_genes),
            ),
            shape=values.shape,
        ).asformat(storage.split("_")[0])
        assert expression.nnz == values.size
        assert (expression.data == 0.0).any()
    else:
        expression = sparse.csr_matrix(values).asformat(storage)

    adata = AnnData(
        expression if layer is None else np.full(values.shape, 99.0),
        obs=pd.DataFrame(
            {
                "cell_type": (["A"] * 4 + ["B"] * 4) * 2,
                "sample": ["S1"] * 8 + ["S2"] * 8,
            },
            index=[f"cell_{i}" for i in range(16)],
        ),
        var=pd.DataFrame(index=["PROD", "SENSOR", "ZERO", "ALL"]),
    )
    if layer is not None:
        adata.layers[layer] = expression
    enzyme = pd.DataFrame([{
        "metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": "PROD",
        "role": "production", "reaction": "prod",
    }])
    sensor = pd.DataFrame([{
        "metabolite": "Met", "hmdb_id": "HMDB0000001",
        "sensor_gene": "SENSOR", "sensor_type": "Transporter",
    }])
    return adata, enzyme, sensor


@pytest.mark.parametrize(
    "storage", ["dense", "csr", "csc", "csr_explicit_zeros", "csc_explicit_zeros"]
)
@pytest.mark.parametrize("layer", [None, "expression"])
def test_expression_fractions_count_positive_cells_independently_of_storage(storage, layer):
    adata, _, _ = _case(storage, layer)
    expression = adata.layers[layer] if layer is not None else adata.X
    before = expression.copy()

    actual = _compute_celltype_expr_frac(adata, layer=layer).loc[["A", "B"]]
    expected = pd.DataFrame(
        [[1.0, 0.5, 0.0, 1.0], [1.0, 0.875, 0.0, 1.0]],
        index=["A", "B"], columns=["PROD", "SENSOR", "ZERO", "ALL"],
    )
    pd.testing.assert_frame_equal(actual, expected, check_exact=True)

    # Counting must not clean up or otherwise mutate the supplied matrix.
    if sparse.issparse(expression):
        for attribute in ("data", "indices", "indptr"):
            np.testing.assert_array_equal(getattr(expression, attribute), getattr(before, attribute))
    else:
        np.testing.assert_array_equal(expression, before)


@pytest.mark.parametrize("storage", ["csr", "csc", "csr_explicit_zeros", "csc_explicit_zeros"])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_sparse_storage_preserves_expression_gate_scores_and_permutation_results(storage, sample_mode):
    dense, enzyme, sensor = _case("dense")
    other, _, _ = _case(storage)
    kwargs = dict(
        sample_key="sample", sample_mode=sample_mode, min_cells=1,
        min_expr_frac=0.75, n_perms=19, random_state=0,
        store_null_scores=sample_mode == "sample_aware",
    )
    expected = run_cell_mesh(dense, enzyme, sensor, **kwargs)
    actual = run_cell_mesh(other, enzyme, sensor, **kwargs)

    for attribute in ("receiver_scores", "sender_scores", "events"):
        pd.testing.assert_frame_equal(
            getattr(actual, attribute), getattr(expected, attribute), check_exact=True,
        )
    receivers = actual.receiver_scores.set_index("receiver")
    assert receivers.loc["A", "sensor_expr_frac"] == 0.5
    assert receivers.loc["A", "sensor_score"] == 0.0
    assert receivers.loc["B", "sensor_score"] > 0.0
    assert (actual.events.loc[actual.events["receiver"] == "A", "cell_mesh_score"] == 0.0).all()
    assert actual.events["perm_pvalue"].notna().all()

    if sample_mode == "sample_aware":
        for attribute in ("sample_receiver_scores", "sample_sender_scores", "sample_events"):
            pd.testing.assert_frame_equal(
                getattr(actual, attribute), getattr(expected, attribute), check_exact=True,
            )
        pd.testing.assert_frame_equal(
            actual.events.attrs["sample_aware_null_scores"],
            expected.events.attrs["sample_aware_null_scores"],
            check_exact=True,
        )
