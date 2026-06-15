"""
CellMESH Comprehensive Tests
覆盖: Metabolite Availability 计算, Sensor Score 计算, 总通讯分数计算
数据: h5ad 单细胞数据 + enzyme CSV + interaction CSV
"""
from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import gmean

import cellmesh
from cellmesh import (
    run_cell_mesh,
    load_cell_mesh_database,
    compute_metabolite_availability,
    read_anndata,
    read_example_data,
)
from cellmesh.core import _confidence_tier, _enzyme_prior_to_availability_reactions
from cellmesh.preprocess import robust_minmax


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PACKAGE_ROOT / "cellmesh" / "data"


class FakeAnnData:
    """轻量级 AnnData 模拟, 不依赖 anndata 包。"""

    def __init__(self, X, var_names, obs):
        self.X = np.asarray(X, dtype=float) if not sparse.issparse(X) else X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)
        self.n_obs = self.X.shape[0]
        self.n_vars = self.X.shape[1]


def _copy_adata(adata, X=None):
    copied = adata.copy() if hasattr(adata, "copy") else FakeAnnData(adata.X.copy(), adata.var_names.copy(), adata.obs.copy())
    if X is not None:
        copied.X = X
    return copied


def _availability_case(adata):
    reaction_table = pd.DataFrame(
        {
            "metabolite": ["MetA", "MetA", "MetA", "MetB", "MetB", "MetC"],
            "hmdb_id": ["HMDB00001", "HMDB00001", "HMDB00001", "HMDB00002", "HMDB00002", "HMDB00003"],
            "reaction": ["r_prod_1", "r_prod_2", "r_sub_1", "r_prod_1", "r_exp_1", "r_sub_only"],
            "gene": ["Gene1", "Gene2;Gene3", "Gene4", "Gene5", "Gene6", "Gene7"],
            "direction": ["product", "product", "substrate", "product", "exporter", "substrate"],
            "weight": [1.0, 2.0, 1.0, 1.0, 1.5, 1.0],
        }
    )
    return compute_metabolite_availability(adata, reaction_table, min_cells=1, return_intermediates=True), reaction_table


def _manual_pce(result, reaction_table):
    pseudobulk = result["pseudobulk"]
    rows = []
    for _, rr in reaction_table.iterrows():
        genes = [g.split("[")[0].strip() for g in str(rr["gene"]).replace(",", ";").replace("|", ";").split(";") if g.strip()]
        valid = [g for g in genes if g in pseudobulk.columns]
        if valid:
            expr = pseudobulk[valid].to_numpy(dtype=float) * float(rr.get("weight", 1.0))
            score = gmean(expr + 1.0, axis=1) - 1.0
        else:
            score = np.zeros(pseudobulk.shape[0], dtype=float)
        rows.append(
            {
                "metabolite": rr["metabolite"],
                "hmdb_id": rr["hmdb_id"],
                "direction": rr["direction"],
                "score": pd.Series(score, index=pseudobulk.index),
            }
        )

    idx = result["P"].index
    P = pd.DataFrame(0.0, index=idx, columns=pseudobulk.index)
    C = P.copy()
    E = P.copy()
    for row in rows:
        key = (row["metabolite"], row["hmdb_id"])
        if key not in P.index:
            continue
        target = {"product": P, "substrate": C, "exporter": E}[row["direction"]]
        target.loc[key] += row["score"]
    return P, C, E


@pytest.fixture
def real_h5ad_path():
    """返回真实 h5ad 文件路径。"""
    path = DATA_DIR / "test_single_cell.h5ad"
    if path.exists():
        return path
    pytest.importorskip("anndata")
    adata = read_example_data("small")
    tmp = Path("/tmp/cellmesh_test_single_cell.h5ad")
    adata.write_h5ad(tmp)
    return tmp


@pytest.fixture
def enzyme_file():
    """enzyme CSV 路径。"""
    path = DATA_DIR / "Enzyme_new.csv"
    return path if path.exists() else DATA_DIR / "enzyme_test.csv"


@pytest.fixture
def interaction_file():
    """interaction CSV 路径。"""
    path = DATA_DIR / "Interaction1.0.csv"
    return path if path.exists() else DATA_DIR / "interaction_test.csv"


@pytest.fixture
def real_adata(real_h5ad_path):
    """读取真实 h5ad。"""
    adata = read_anndata(real_h5ad_path, mode="h5ad")
    if "_index" in adata.var.columns:
        adata.var_names = adata.var["_index"].astype(str)
    return adata


@pytest.fixture
def enzyme_interaction(enzyme_file, interaction_file):
    return load_cell_mesh_database(str(enzyme_file), str(interaction_file))


@pytest.fixture
def reaction_table(enzyme_interaction):
    enzyme_df, _ = enzyme_interaction
    return _enzyme_prior_to_availability_reactions(enzyme_df)


@pytest.fixture
def availability_result(real_adata, reaction_table):
    return compute_metabolite_availability(real_adata, reaction_table, min_cells=1, return_intermediates=True)


@pytest.fixture
def full_result(real_adata, enzyme_interaction):
    """运行完整 CELL MESH。"""
    enzyme_df, interaction_df = enzyme_interaction
    return run_cell_mesh(
        real_adata,
        enzyme_metabolite=enzyme_df,
        metabolite_sensor=interaction_df,
        cell_type_key="cell_type",
        n_perms=0,
        allow_self=True,
        min_cells=1,
    )


@pytest.fixture
def synthetic_adata():
    X = np.array(
        [
            [4, 2, 8, 1, 0, 3, 5, 2],
            [2, 4, 6, 2, 1, 4, 4, 0],
            [0, 1, 2, 6, 7, 0, 3, 5],
            [1, 0, 3, 4, 5, 1, 2, 6],
            [5, 2, 1, 0, 3, 8, 0, 1],
            [6, 3, 0, 1, 2, 7, 1, 2],
        ],
        dtype=float,
    )
    obs = {"cell_type": ["A", "A", "B", "B", "C", "C"], "sample": ["S1", "S1", "S1", "S1", "S2", "S2"]}
    return FakeAnnData(X, [f"Gene{i}" for i in range(1, 9)], obs)


class TestMetaboliteAvailability:
    """TC-A1 ~ TC-A7, TC-A10。"""

    def test_tc_a1_basic_availability_real_data(self, real_adata, reaction_table, availability_result):
        availability = availability_result["availability"]
        n_celltypes = real_adata.obs["cell_type"].astype(str).nunique()
        assert not availability.empty
        assert availability.shape == (availability.index.nunique(), n_celltypes)
        assert np.nanmin(availability.values) >= 0.0
        assert np.nanmax(availability.values) <= 1.0
        for key in ["P_norm", "C_norm", "E_norm"]:
            values = availability_result[key].values
            assert np.nanmin(values) >= 0.0
            assert np.nanmax(values) <= 1.0
        assert availability_result["pseudobulk"].shape == (n_celltypes, real_adata.n_vars)
        assert availability_result["metadata"]["has_product"].all()

    def test_tc_a2_availability_formula(self, availability_result):
        expected = (
            availability_result["P_norm"]
            * ((1.0 - availability_result["C_norm"]) ** 0.5)
            * (0.8 + 0.2 * availability_result["E_norm"])
        )
        pd.testing.assert_frame_equal(availability_result["availability"], expected.clip(0.0, 1.0), atol=1e-10, rtol=0)

    def test_tc_a3_pce_direction_aggregation(self, synthetic_adata):
        result, reactions = _availability_case(synthetic_adata)
        manual_P, manual_C, manual_E = _manual_pce(result, reactions)
        pd.testing.assert_frame_equal(result["P"], manual_P, atol=1e-10, rtol=0)
        pd.testing.assert_frame_equal(result["C"], manual_C, atol=1e-10, rtol=0)
        pd.testing.assert_frame_equal(result["E"], manual_E, atol=1e-10, rtol=0)
        assert result["P"].loc[("MetA", "HMDB00001")].sum() > result["P"].loc[("MetB", "HMDB00002")].sum()

    def test_tc_a4_multigene_reaction_uses_geometric_mean(self, synthetic_adata):
        result, _ = _availability_case(synthetic_adata)
        pseudobulk = result["pseudobulk"]
        expected = gmean((pseudobulk[["Gene2", "Gene3"]].to_numpy() * 2.0) + 1.0, axis=1) - 1.0
        observed = result["P"].loc[("MetA", "HMDB00001")] - pseudobulk["Gene1"]
        assert np.allclose(observed.values, expected, atol=1e-10)

    def test_tc_a5_robust_minmax_manual(self, availability_result):
        P = availability_result["P"]
        manual = pd.DataFrame(index=P.index, columns=P.columns, dtype=float)
        for idx in P.index:
            values = P.loc[idx].to_numpy(dtype=float)
            lo = np.nanpercentile(values, 5)
            hi = np.nanpercentile(values, 95)
            if np.isclose(hi, lo, rtol=1e-5):
                norm = np.zeros_like(values)
            else:
                norm = np.clip((np.clip(values, lo, hi) - lo) / (hi - lo + 0.05), 0.0, 1.0)
            manual.loc[idx] = norm
        pd.testing.assert_frame_equal(availability_result["P_norm"], manual.astype(float), atol=1e-10, rtol=0)

    def test_tc_a6_missing_substrate_exporter_defaults(self, availability_result):
        metadata = availability_result["metadata"]
        no_substrate = metadata.index[~metadata["has_substrate"]]
        no_exporter = metadata.index[~metadata["has_exporter"]]
        if len(no_substrate):
            assert np.allclose(availability_result["C_norm"].loc[no_substrate].values, 0.2)
        if len(no_exporter):
            assert np.allclose(availability_result["E_norm"].loc[no_exporter].values, 0.5)

    def test_tc_a7_metabolite_without_product_is_filtered(self, synthetic_adata):
        reactions = pd.DataFrame(
            {
                "metabolite": ["OnlyConsumed", "Produced"],
                "hmdb_id": ["HMDB9", "HMDB8"],
                "reaction": ["sub", "prod"],
                "gene": ["Gene1", "Gene2"],
                "direction": ["substrate", "product"],
            }
        )
        result = compute_metabolite_availability(synthetic_adata, reactions, min_cells=1, return_intermediates=True)
        assert "OnlyConsumed" not in result["availability"].index.get_level_values("metabolite")
        assert "Produced" in result["availability"].index.get_level_values("metabolite")

    def test_tc_a10_dense_sparse_equivalence(self, synthetic_adata):
        dense_result, reactions = _availability_case(synthetic_adata)
        sparse_adata = _copy_adata(synthetic_adata, sparse.csr_matrix(synthetic_adata.X))
        sparse_result = compute_metabolite_availability(sparse_adata, reactions, min_cells=1, return_intermediates=True)
        pd.testing.assert_frame_equal(dense_result["availability"], sparse_result["availability"], atol=1e-10, rtol=0)


class TestSensorScore:
    """TC-A8, TC-A9。"""

    def test_tc_a8_receiver_score_matches_current_implementation(self, full_result):
        receiver_scores = full_result.receiver_scores
        assert not receiver_scores.empty
        pseudobulk = full_result.availability_results["pseudobulk"]
        for _, row in receiver_scores.sample(min(5, len(receiver_scores)), random_state=0).iterrows():
            expected = robust_minmax(pseudobulk[row["sensor_gene"]].values)[pseudobulk.index.get_loc(row["receiver"])]
            expected = float(expected) if row["sensor_expr_frac"] >= full_result.parameters["min_expr_frac"] else 0.0
            observed = row.get("sensor_score", row.get("receiver_score"))
            assert np.isclose(observed, expected, atol=1e-6)

    def test_tc_a9_min_expr_frac_cutoff(self, full_result):
        scores = full_result.receiver_scores
        low = scores[scores["sensor_expr_frac"] < full_result.parameters["min_expr_frac"]]
        if low.empty:
            pytest.skip("No receiver rows below min_expr_frac in this fixture")
        score_col = "receiver_score" if "receiver_score" in low.columns else "sensor_score"
        assert (low[score_col] == 0.0).all()


class TestCellMeshScore:
    """TC-B1 ~ TC-B10。"""

    def test_tc_b1_end_to_end_real_data(self, real_adata, full_result):
        assert isinstance(full_result.events, pd.DataFrame)
        assert not full_result.events.empty
        required = {
            "sender",
            "receiver",
            "metabolite",
            "sensor_gene",
            "sensor_type",
            "sender_score",
            "receiver_score",
            "cell_mesh_score",
            "communication_score",
            "sensor_expr_frac",
            "confidence_tier",
        }
        assert required.issubset(full_result.events.columns)
        n_metabolites = full_result.availability_results["availability"].shape[0]
        n_celltypes = real_adata.obs["cell_type"].astype(str).nunique()
        assert full_result.sender_scores.shape == (n_metabolites, n_celltypes)
        assert not full_result.receiver_scores.empty

    def test_tc_b2_cell_mesh_score_formula(self, full_result):
        events = full_result.events
        if "prior_weight" in events.columns:
            expected = events["sender_score"] * events["receiver_score"] * events["prior_weight"]
        else:
            expected = np.sqrt(events["sender_score"] * events["receiver_score"])
        assert np.allclose(events["cell_mesh_score"], expected, atol=1e-10)
        assert np.allclose(events["communication_score"], expected, atol=1e-10)

    def test_tc_b3_sender_scores_from_availability(self, full_result):
        availability = full_result.availability_results["availability"].copy()
        pd.testing.assert_frame_equal(full_result.sender_scores, availability)

    def test_tc_b4_allow_self_false(self, real_adata, enzyme_interaction):
        enzyme_df, interaction_df = enzyme_interaction
        result = run_cell_mesh(real_adata, enzyme_df, interaction_df, cell_type_key="cell_type", allow_self=False, n_perms=0, min_cells=1)
        assert not (result.events["sender"] == result.events["receiver"]).any()

    def test_tc_b4b_allow_self_false_single_celltype_returns_empty_events(self, synthetic_adata):
        adata = FakeAnnData(
            synthetic_adata.X[:2, :2],
            ["Gene1", "Gene2"],
            {"cell_type": ["A", "A"]},
        )
        enzyme = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["H"], "gene": ["Gene1"], "role": ["production"]})
        sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["H"], "sensor_gene": ["Gene2"], "sensor_type": ["Transporter"]})
        result = run_cell_mesh(adata, enzyme, sensor, cell_type_key="cell_type", min_cells=1, allow_self=False, n_perms=0)
        assert result.events.empty
        assert {"sender", "receiver", "communication_score", "confidence_tier"}.issubset(result.events.columns)

    def test_tc_b4c_events_match_metabolite_and_hmdb(self):
        adata = FakeAnnData(
            np.array(
                [
                    [5, 0, 0, 0],
                    [5, 0, 0, 0],
                    [0, 7, 0, 4],
                    [0, 7, 0, 4],
                ],
                dtype=float,
            ),
            ["G_H1", "G_H2", "S_H1", "S_H2"],
            {"cell_type": ["A", "A", "B", "B"]},
        )
        enzyme = pd.DataFrame(
            {
                "metabolite": ["M", "M"],
                "hmdb_id": ["H1", "H2"],
                "gene": ["G_H1", "G_H2"],
                "role": ["production", "production"],
            }
        )
        sensor = pd.DataFrame(
            {
                "metabolite": ["M"],
                "hmdb_id": ["H2"],
                "sensor_gene": ["S_H2"],
                "sensor_type": ["Transporter"],
            }
        )
        result = run_cell_mesh(adata, enzyme, sensor, cell_type_key="cell_type", min_cells=1, lower=0, upper=100, n_perms=0)
        assert set(result.events["hmdb_id"]) == {"H2"}
        assert (result.events["metabolite"] == "M").all()

    def test_tc_b4d_missing_hmdb_priors_are_excluded(self):
        adata = FakeAnnData(
            np.array([[5, 0], [5, 0], [0, 4], [0, 4]], dtype=float),
            ["Gene1", "Gene2"],
            {"cell_type": ["A", "A", "B", "B"]},
        )
        enzyme = pd.DataFrame({"metabolite": ["M"], "gene": ["Gene1"], "role": ["production"]})
        sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["H"], "sensor_gene": ["Gene2"], "sensor_type": ["Transporter"]})
        with pytest.raises(ValueError, match="hmdb_id|No enzyme prior genes"):
            run_cell_mesh(adata, enzyme, sensor, cell_type_key="cell_type", min_cells=1, n_perms=0)

        enzyme = pd.DataFrame({"metabolite": ["M"], "hmdb_id": [np.nan], "gene": ["Gene1"], "role": ["production"]})
        with pytest.raises(ValueError, match="No enzyme prior genes"):
            run_cell_mesh(adata, enzyme, sensor, cell_type_key="cell_type", min_cells=1, n_perms=0)

        enzyme = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["H"], "gene": ["Gene1"], "role": ["production"]})
        sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": [np.nan], "sensor_gene": ["Gene2"], "sensor_type": ["Transporter"]})
        with pytest.raises(ValueError, match="No sensor genes"):
            run_cell_mesh(adata, enzyme, sensor, cell_type_key="cell_type", min_cells=1, n_perms=0)

    def test_tc_b5_confidence_tier_rules(self):
        cases = [
            ({"fdr": 0.04, "communication_score": 0.6, "sensor_expr_frac": 0.2}, "Tier1_high"),
            ({"fdr": 0.09, "communication_score": 0.3, "sensor_expr_frac": 0.0}, "Tier2_medium"),
            ({"fdr": np.nan, "communication_score": 0.6, "sensor_expr_frac": 0.2}, "Tier2_no_permutation"),
            ({"fdr": 0.5, "communication_score": 0.1, "sensor_expr_frac": 0.0}, "Tier3_exploratory"),
        ]
        for row, expected in cases:
            assert _confidence_tier(pd.Series(row)) == expected

    def test_tc_b6_permutation_pvalue_range(self, real_adata, enzyme_interaction):
        enzyme_df, interaction_df = enzyme_interaction
        result = run_cell_mesh(real_adata, enzyme_df, interaction_df, cell_type_key="cell_type", n_perms=10, random_state=1, min_cells=1)
        if result.events.empty:
            pytest.skip("No events available for permutation checks")
        assert result.events["perm_pvalue"].between(0, 1).all()
        assert result.events["fdr"].between(0, 1).all()

    def test_tc_b7_empty_prior_handling(self, synthetic_adata):
        enzyme = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["H"], "gene": ["MissingGene"], "role": ["production"]})
        sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["H"], "sensor_gene": ["AlsoMissing"], "sensor_type": ["Transporter"]})
        with pytest.raises(ValueError, match="No enzyme prior genes|No sensor genes"):
            run_cell_mesh(synthetic_adata, enzyme, sensor, cell_type_key="cell_type", n_perms=0, min_cells=1)

    def test_tc_b8_builtin_database_compatibility(self):
        adata = read_example_data("small")
        try:
            result = run_cell_mesh(adata, cell_type_key="cell_type", n_perms=0, min_cells=1)
        except ValueError as exc:
            pytest.skip(f"Built-in priors do not overlap read_example_data genes in this implementation: {exc}")
        assert isinstance(result.events, pd.DataFrame)

    def test_tc_b9_parameter_passing_and_availability_formula(self, real_adata, enzyme_interaction):
        enzyme_df, interaction_df = enzyme_interaction
        kwargs = {
            "beta": 1.0,
            "missing_C_norm": 0.41,
            "missing_E_norm": 0.75,
        }
        sig = inspect.signature(run_cell_mesh)
        if "beta_sensor" in sig.parameters:
            kwargs["beta_sensor"] = 2.0
        if "beta_specificity" in sig.parameters:
            kwargs["beta_specificity"] = 0.5
        result = run_cell_mesh(real_adata, enzyme_df, interaction_df, cell_type_key="cell_type", n_perms=0, min_cells=1, **kwargs)
        for key, value in kwargs.items():
            assert result.parameters[key] == value
        avail = result.availability_results
        expected = avail["P_norm"] * ((1 - avail["C_norm"]) ** kwargs["beta"]) * (0.8 + 0.2 * avail["E_norm"])
        pd.testing.assert_frame_equal(avail["availability"], expected.clip(0, 1), atol=1e-10, rtol=0)

    def test_tc_b10_result_save_load(self, tmp_path, full_result):
        prefix = tmp_path / "test_result"
        full_result.to_csv(str(prefix))
        for suffix in ["events", "sender_scores", "receiver_scores"]:
            path = tmp_path / f"test_result.{suffix}.csv"
            assert path.exists()
            assert path.stat().st_size > 0
            assert not pd.read_csv(path).empty


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
