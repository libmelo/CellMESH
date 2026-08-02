import numpy as np
import pandas as pd
import pytest
from scipy.stats import gmean

from cellmesh import (
    load_cell_mesh_database, 
    run_cell_mesh, 
    compute_metabolite_availability
)
from cellmesh.config import MIN_CELL_COUNT


class FakeAnnData:
    def __init__(self, X, var_names, obs):
        self.X = X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)


def test_compute_metabolite_availability():
    """测试新的 metabolite availability 计算功能"""
    enzyme_metabolite = pd.DataFrame({
        'metabolite': ['MetA', 'MetA', 'MetA', 'MetB', 'MetB', 'MetC'],
        'hmdb_id': ['HMDB001', 'HMDB001', 'HMDB001', 'HMDB002', 'HMDB002', 'HMDB003'],
        'reaction': ['R1', 'R2', 'R3', 'R4', 'R5', 'R6'],
        'gene': ['Gene1', 'Gene2;Gene3', 'Gene4', 'Gene5', 'Gene6', 'Gene7'],
        'role': ['production', 'production', 'degradation', 'production', 'export', 'production']
    })
    
    # 创建测试数据
    n_cells = 12
    n_genes = 7
    genes = [f'Gene{i+1}' for i in range(n_genes)]
    
    X = np.zeros((n_cells, n_genes), dtype=float)
    cell_types = ['A'] * 6 + ['B'] * 6
    
    # 在细胞类型 A 中高表达 Gene1, Gene2, Gene3
    X[0:6, 0] += 3  # Gene1
    X[0:6, 1] += 2  # Gene2
    X[0:6, 2] += 2  # Gene3
    
    # 在细胞类型 B 中高表达 Gene5, Gene6
    X[6:12, 4] += 3  # Gene5
    X[6:12, 5] += 2  # Gene6
    X[:, 6] = 1       # Gene7
    
    adata = FakeAnnData(X, genes, {"cell_type": cell_types})
    
    # Toy data has only 6 cells per group, so this unit test lowers min_cells
    # explicitly. The public default remains MIN_CELL_COUNT=100 and is tested
    # separately below.
    result = compute_metabolite_availability(
        adata,
        enzyme_metabolite,
        celltype_col='cell_type',
        min_cells=1,
        return_intermediates=True
    )
    
    # 检查结果
    assert 'availability' in result
    assert 'P' in result
    assert 'C' in result
    assert 'E' in result
    assert 'P_score' in result
    assert 'C_score' in result
    assert 'E_score' in result
    assert 'metadata' in result
    assert result["pce_reference"] == "mean"
    
    availability = result['availability']
    assert not availability.empty
    metA_avail = availability.loc[('MetA', 'HMDB001')]
    assert metA_avail.loc['A'] > metA_avail.loc['B']

    metc_idx = ('MetC', 'HMDB003')
    assert metc_idx in availability.index
    assert np.allclose(
        result['C_score'].loc[metc_idx].values,
        0.0,
    )
    assert np.allclose(
        result['E_score'].loc[metc_idx].values,
        0.0,
    )


def test_compute_availability_default_min_cells_is_qc_only():
    adata = FakeAnnData(
        np.array([[5], [5], [0], [0]], dtype=float),
        ["Gene1"],
        {"cell_type": ["A", "A", "B", "B"]},
    )
    enzyme_metabolite = pd.DataFrame(
        {
            "metabolite": ["M"],
            "hmdb_id": ["HMDB00001"],
            "gene": ["Gene1"],
            "role": ["production"],
            "reaction": ["prod"],
        }
    )

    result = compute_metabolite_availability(adata, enzyme_metabolite)

    # The public default remains 100, but it is a QC threshold rather than an
    # analysis filter. Both two-cell types remain in every numerical result.
    assert set(result["pseudobulk"].index) == {"A", "B"}
    assert set(result["cell_fractions"].index) == {"A", "B"}
    assert set(result["availability"].columns) == {"A", "B"}
    assert result["cell_counts"].to_dict() == {"A": 2, "B": 2}
    qc = result["celltype_qc"]
    assert qc["n_cells"].to_dict() == {"A": 2, "B": 2}
    assert qc["cell_fraction"].to_dict() == {"A": 0.5, "B": 0.5}
    assert not qc["passes_min_cells"].any()


@pytest.mark.parametrize("reaction", [None, "   "])
def test_compute_availability_requires_nonempty_reaction(reaction):
    adata = FakeAnnData(
        np.array([[5], [2]], dtype=float),
        ["Gene1"],
        {"cell_type": ["A", "B"]},
    )
    enzyme = pd.DataFrame(
        {
            "metabolite": ["M"],
            "hmdb_id": ["HMDB00001"],
            "gene": ["Gene1"],
            "role": ["production"],
            "reaction": [reaction],
        }
    )

    with pytest.raises(ValueError, match="reaction values must be non-empty"):
        compute_metabolite_availability(adata, enzyme, min_cells=1)


def test_reactions_group_by_hmdb_and_deduplicate_genes_across_name_variants():
    adata = FakeAnnData(
        np.array([[3.0, 3.0, 1.0], [1.0, 1.0, 3.0]]),
        ["G1", "G2", "G3"],
        {"cell_type": ["A", "B"]},
    )
    enzyme = pd.DataFrame(
        {
            "metabolite": ["Canonical name", "Alias name", "Alias name", "Alias name"],
            "hmdb_id": [" hmdb00001 ", "HMDB00001", "HMDB00001", "HMDB00001"],
            "gene": ["G1; G1", "G1", "G2", "G3"],
            "role": ["production"] * 4,
            "reaction": ["shared", "shared", "shared", "second"],
        }
    )

    result = compute_metabolite_availability(adata, enzyme, min_cells=1)

    assert result["P"].index.tolist() == [("Canonical name", "HMDB00001")]
    reactions = result["reaction_genes"]
    assert reactions["metabolite"].unique().tolist() == ["Canonical name"]
    assert reactions[["hmdb_id", "reaction", "direction"]].to_dict("records") == [
        {"hmdb_id": "HMDB00001", "reaction": "shared", "direction": "product"},
        {"hmdb_id": "HMDB00001", "reaction": "second", "direction": "product"},
    ]
    assert reactions.loc[reactions["reaction"] == "shared", "genes"].iloc[0] == ["G1", "G2"]
    expected_shared = np.sqrt((adata.X[:, 0] + 1.0) * (adata.X[:, 1] + 1.0)) - 1.0
    expected_p = (expected_shared + adata.X[:, 2]) * 0.5
    assert np.allclose(result["P"].iloc[0].to_numpy(dtype=float), expected_p)
    assert result["metadata"].iloc[0]["n_product_reactions"] == 2


def test_sender_pce_uses_configurable_cell_fraction_exponent_before_normalization():
    adata = FakeAnnData(
        np.full((7, 1), 2.0),
        ["PROD"],
        {
            "cell_type": [
                "Abundant",
                "Abundant",
                "Abundant",
                "Abundant",
                "Rare",
                "Rare",
                "Ineligible",
            ]
        },
    )
    enzyme = pd.DataFrame(
        {
            "metabolite": ["Met"],
            "hmdb_id": ["HMDB00001"],
            "reaction": ["prod"],
            "gene": ["PROD"],
            "role": ["production"],
        }
    )

    result = compute_metabolite_availability(
        adata,
        enzyme,
        min_cells=2,
        return_intermediates=True,
    )

    # All cell types have the same per-cell mean, but their population P values
    # differ because abundance uses each type's fraction of all seven cells.
    assert result["pseudobulk"].loc["Abundant", "PROD"] == pytest.approx(2.0)
    assert result["pseudobulk"].loc["Rare", "PROD"] == pytest.approx(2.0)
    assert result["pseudobulk"].loc["Ineligible", "PROD"] == pytest.approx(2.0)
    assert result["cell_fractions"].loc["Abundant"] == pytest.approx(4.0 / 7.0)
    assert result["cell_fractions"].loc["Rare"] == pytest.approx(2.0 / 7.0)
    assert result["cell_fractions"].loc["Ineligible"] == pytest.approx(1.0 / 7.0)
    assert result["cell_fractions"].sum() == pytest.approx(1.0)
    assert result["celltype_qc"]["passes_min_cells"].to_dict() == {
        "Abundant": True,
        "Rare": True,
        "Ineligible": False,
    }
    assert result["sender_abundance_exponent"] == 1.0
    pd.testing.assert_series_equal(
        result["sender_abundance_weights"],
        result["cell_fractions"].rename("sender_abundance_weight"),
    )
    assert result["P"].loc[("Met", "HMDB00001"), "Abundant"] == pytest.approx(8.0 / 7.0)
    assert result["P"].loc[("Met", "HMDB00001"), "Rare"] == pytest.approx(4.0 / 7.0)
    assert result["P"].loc[("Met", "HMDB00001"), "Ineligible"] == pytest.approx(2.0 / 7.0)
    assert result["P_ref"].loc[("Met", "HMDB00001")] == pytest.approx(2.0 / 3.0)
    assert result["P_score"].loc[("Met", "HMDB00001"), "Abundant"] == pytest.approx(
        12.0 / 19.0
    )
    assert result["P_score"].loc[("Met", "HMDB00001"), "Rare"] == pytest.approx(6.0 / 13.0)
    assert result["P_score"].loc[("Met", "HMDB00001"), "Ineligible"] == pytest.approx(0.3)
    pd.testing.assert_frame_equal(result["availability"], result["P_score"])

    tempered = compute_metabolite_availability(
        adata,
        enzyme,
        min_cells=2,
        sender_abundance_exponent=0.5,
    )
    assert tempered["sender_abundance_weights"].loc["Abundant"] == pytest.approx(
        np.sqrt(4.0 / 7.0)
    )
    assert tempered["P"].loc[("Met", "HMDB00001"), "Rare"] == pytest.approx(
        2.0 * np.sqrt(2.0 / 7.0)
    )

    unadjusted = compute_metabolite_availability(
        adata,
        enzyme,
        min_cells=2,
        sender_abundance_exponent=0.0,
    )
    assert np.allclose(unadjusted["sender_abundance_weights"], 1.0)
    assert np.allclose(unadjusted["P"].loc[("Met", "HMDB00001")], 2.0)
    assert np.allclose(unadjusted["P_score"].loc[("Met", "HMDB00001")], 0.5)
    assert np.allclose(unadjusted["availability"].loc[("Met", "HMDB00001")], 0.5)

    for invalid, error_type in [
        (-0.1, ValueError),
        (np.nan, ValueError),
        (np.inf, ValueError),
        ("invalid", TypeError),
        (False, TypeError),
    ]:
        with pytest.raises(error_type, match="sender_abundance_exponent"):
            compute_metabolite_availability(
                adata,
                enzyme,
                min_cells=2,
                sender_abundance_exponent=invalid,
            )


def test_sender_formula_excludes_export_and_tracks_missing_prior_status():
    adata = FakeAnnData(
        np.array(
            [
                [1.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0, 0.0],
                [2.0, 2.0, 1.0, 0.0, 0.0],
                [2.0, 2.0, 1.0, 0.0, 0.0],
                [4.0, 4.0, 10.0, 0.0, 0.0],
                [4.0, 4.0, 10.0, 0.0, 0.0],
            ]
        ),
        ["PROD", "CONS", "EXPORT", "CONS_ZERO", "EXPORT_ZERO"],
        {"cell_type": ["A", "A", "B", "B", "C", "C"]},
    )
    product = {
        "metabolite": "Met",
        "hmdb_id": "HMDB00001",
        "reaction": "prod",
        "gene": "PROD",
        "role": "production",
    }
    consumption = {
        "metabolite": "Met",
        "hmdb_id": "HMDB00001",
        "reaction": "cons",
        "gene": "CONS",
        "role": "degradation",
    }
    export = {
        "metabolite": "Met",
        "hmdb_id": "HMDB00001",
        "reaction": "export",
        "gene": "EXPORT",
        "role": "export",
    }

    product_only = compute_metabolite_availability(
        adata, pd.DataFrame([product]), min_cells=2
    )
    product_only_median = compute_metabolite_availability(
        adata,
        pd.DataFrame([product]),
        min_cells=2,
        pce_reference="median",
    )
    with_consumption = compute_metabolite_availability(
        adata, pd.DataFrame([product, consumption]), min_cells=2
    )
    with_consumption_and_export = compute_metabolite_availability(
        adata, pd.DataFrame([product, consumption, export]), min_cells=2
    )
    idx = ("Met", "HMDB00001")

    pd.testing.assert_frame_equal(
        product_only["availability"], product_only["P_score"]
    )
    assert product_only["pce_reference"] == "mean"
    assert product_only_median["pce_reference"] == "median"
    assert product_only["P_ref"].loc[idx] != product_only_median["P_ref"].loc[idx]
    for removed in [
        "P_plus",
        "relative_consumption_support",
        "E_plus",
        "P_contrast",
        "C_contrast",
        "E_contrast",
    ]:
        assert removed not in product_only
    assert product_only["metadata"].loc[idx, "consumption_status"] == "prior_missing"
    assert product_only["metadata"].loc[idx, "export_status"] == "prior_missing"

    denominator = (
        with_consumption["P_score"] + with_consumption["C_score"]
    )
    expected = (
        with_consumption["P_score"].pow(2)
        .div(denominator.where(denominator > 0.0))
        .fillna(0.0)
    )
    pd.testing.assert_frame_equal(with_consumption["availability"], expected)
    assert with_consumption["metadata"].loc[idx, "consumption_status"] == "supported"

    # E remains available as support evidence but cannot increase, decrease,
    # or cancel C in the formal sender score.
    assert with_consumption_and_export["E_score"].loc[idx].max() > 0.0
    pd.testing.assert_frame_equal(
        with_consumption_and_export["availability"],
        with_consumption["availability"],
    )
    assert with_consumption_and_export["metadata"].loc[idx, "export_status"] == "supported"
    assert "export_used_in_sender_score" not in with_consumption_and_export["metadata"]

    no_consumption_expression = compute_metabolite_availability(
        adata,
        pd.DataFrame(
            [
                product,
                {
                    **consumption,
                    "gene": "CONS_ZERO",
                    "reaction": "cons_zero",
                },
            ]
        ),
        min_cells=2,
    )
    assert (
        no_consumption_expression["metadata"].loc[idx, "consumption_status"]
        == "prior_no_expression"
    )
    pd.testing.assert_frame_equal(
        no_consumption_expression["availability"],
        no_consumption_expression["P_score"],
    )

    no_export_expression = compute_metabolite_availability(
        adata,
        pd.DataFrame(
            [
                product,
                {
                    **export,
                    "gene": "EXPORT_ZERO",
                    "reaction": "export_zero",
                },
            ]
        ),
        min_cells=2,
    )
    assert (
        no_export_expression["metadata"].loc[idx, "export_status"]
        == "prior_no_expression"
    )
    assert np.isnan(no_export_expression["E_ref"].loc[idx])
    assert np.allclose(no_export_expression["E_score"].loc[idx], 0.0)
    pd.testing.assert_frame_equal(
        no_export_expression["availability"],
        no_export_expression["P_score"],
    )

    unavailable_consumption = compute_metabolite_availability(
        adata,
        pd.DataFrame(
            [
                product,
                {
                    **consumption,
                    "gene": "NOT_MEASURED",
                    "reaction": "cons_unavailable",
                },
            ]
        ),
        min_cells=2,
    )
    assert (
        unavailable_consumption["metadata"].loc[idx, "consumption_status"]
        == "prior_gene_unavailable"
    )
    assert "has_substrate" not in unavailable_consumption["metadata"]
    assert "has_usable_substrate" not in unavailable_consumption["metadata"]
    pd.testing.assert_frame_equal(
        unavailable_consumption["availability"],
        unavailable_consumption["P_score"],
    )

    for invalid_reference, error_type in [
        (None, TypeError),
        ("", ValueError),
        ("average", ValueError),
        (1, TypeError),
    ]:
        with pytest.raises(error_type, match="pce_reference"):
            compute_metabolite_availability(
                adata,
                pd.DataFrame([product]),
                min_cells=2,
                pce_reference=invalid_reference,
            )


def test_compute_availability_excludes_missing_hmdb_rows():
    adata = FakeAnnData(
        np.array([[5, 0], [5, 0], [0, 3], [0, 3]], dtype=float),
        ["Gene1", "Gene2"],
        {"cell_type": ["A", "A", "B", "B"]},
    )
    enzyme_metabolite = pd.DataFrame(
        {
            "metabolite": ["MissingHMDB", "ValidHMDB"],
            "hmdb_id": [np.nan, "HMDB00001"],
            "gene": ["Gene1", "Gene2"],
            "role": ["production", "production"],
            "reaction": ["missing_prod", "valid_prod"],
        }
    )

    result = compute_metabolite_availability(adata, enzyme_metabolite, min_cells=1)

    assert ("MissingHMDB", np.nan) not in result["availability"].index
    assert result["availability"].index.tolist() == [("ValidHMDB", "HMDB00001")]
    

def test_sparse_vs_dense():
    """测试 dense 和 sparse 输入是否得到一致结果"""
    sparse = pytest.importorskip("scipy.sparse")

    enzyme_metabolite = pd.DataFrame({
        'metabolite': ['MetX', 'MetX', 'MetY'],
        'hmdb_id': ['HMDB001', 'HMDB001', 'HMDB002'],
        'reaction': ['R1', 'R2', 'R3'],
        'gene': ['GeneA', 'GeneB', 'GeneC'],
        'role': ['production', 'degradation', 'production']
    })
    
    X_dense = np.array(
        [
            [1, 0, 2, 0, 1],
            [2, 1, 0, 1, 0],
            [1, 2, 1, 0, 0],
            [0, 1, 3, 1, 1],
            [2, 0, 1, 0, 2],
            [0, 3, 1, 2, 0],
            [1, 2, 0, 1, 1],
            [0, 2, 2, 0, 1],
            [1, 3, 1, 1, 0],
            [0, 1, 0, 2, 2],
        ],
        dtype=float,
    )
    genes = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE']
    cell_types = ['Type1'] * 5 + ['Type2'] * 5
    
    # 创建 dense AnnData
    adata_dense = FakeAnnData(X_dense, genes, {"cell_type": cell_types})
    
    # 创建 sparse AnnData
    X_sparse = sparse.csr_matrix(X_dense)
    adata_sparse = FakeAnnData(X_sparse, genes, {"cell_type": cell_types})
    
    # 计算 dense 结果
    result_dense = compute_metabolite_availability(
        adata_dense, enzyme_metabolite, min_cells=1
    )
    
    # 计算 sparse 结果
    result_sparse = compute_metabolite_availability(
        adata_sparse, enzyme_metabolite, min_cells=1
    )
    
    # 比较结果
    avail_dense = result_dense['availability']
    avail_sparse = result_sparse['availability']
    
    pd.testing.assert_frame_equal(avail_dense, avail_sparse)


def test_boundary_cases():
    """测试边界情况"""
    enzyme1 = pd.DataFrame({
        'metabolite': ['NoProduct'],
        'hmdb_id': ['HMDB000'],
        'reaction': ['R1'],
        'gene': ['Gene1'],
        'role': ['degradation']
    })
    
    adata = FakeAnnData(
        np.array([[1, 2], [3, 4]]),
        ['Gene1', 'Gene2'],
        {'cell_type': ['A', 'B']}
    )
    
    result1 = compute_metabolite_availability(adata, enzyme1, min_cells=1)
    assert result1['availability'].empty, "没有 product reaction 的代谢物应该被跳过"
    
    # 2. 所有 product 基因都缺失 - 应该与没有 product reaction 一样被跳过
    enzyme2 = pd.DataFrame({
        'metabolite': ['MissingGenes'],
        'hmdb_id': ['HMDB001'],
        'reaction': ['R1'],
        'gene': ['NotPresent'],
        'role': ['production']
    })
    
    result2 = compute_metabolite_availability(adata, enzyme2, min_cells=1)
    # availability 应该为空（因为 P 全 0，被过滤掉了）
    assert result2['availability'].empty, "所有 product 基因都缺失的代谢物应该被跳过"
    assert result2["availability"].columns.tolist() == ["A", "B"]
    assert result2["pseudobulk"].shape == (2, 2)
    assert result2["expr_frac"].shape == (2, 2)
    assert not result2["reaction_genes"].empty

    # 2b. Product gene 存在但表达全零时，同样跳过 metabolite，同时保持
    # observed cell-type schema 和真实中间结果。
    adata_zero = FakeAnnData(
        np.zeros((4, 1), dtype=float),
        ["Gene1"],
        {"cell_type": ["A", "A", "B", "B"]},
    )
    enzyme_zero = pd.DataFrame(
        {
            "metabolite": ["ZeroProduct"],
            "hmdb_id": ["HMDB009"],
            "reaction": ["R0"],
            "gene": ["Gene1"],
            "role": ["production"],
        }
    )
    result_zero = compute_metabolite_availability(
        adata_zero,
        enzyme_zero,
        min_cells=1,
    )
    assert result_zero["availability"].shape == (0, 2)
    assert result_zero["availability"].columns.tolist() == ["A", "B"]
    assert result_zero["pseudobulk"].shape == (2, 1)
    assert result_zero["expr_frac"].shape == (2, 1)
    assert not result_zero["reaction_genes"].empty
    
    # 3. 同一 reaction 多个基因和多行 - 应该合并为一个基因集合后按普通几何均值聚合
    enzyme3 = pd.DataFrame({
        'metabolite': ['MultiGene', 'MultiGene'],
        'hmdb_id': ['HMDB002', 'HMDB002'],
        'reaction': ['R1', 'R1'],  # 同一 reaction
        'gene': ['HighExpr;LowExpr', 'AnotherGene'],
        'role': ['production', 'production']
    })
    
    adata3 = FakeAnnData(
        np.array([[10, 1, 0.5], [10, 1, 0.5]]),  # HighExpr = 10, LowExpr = 1, AnotherGene = 0.5
        ['HighExpr', 'LowExpr', 'AnotherGene'],
        {'cell_type': ['A', 'A']}
    )
    
    result3 = compute_metabolite_availability(adata3, enzyme3, min_cells=1)
    expected3 = gmean(np.array([10.0, 1.0, 0.5]) + 1.0) - 1.0
    assert result3['P'].loc[('MultiGene', 'HMDB002'), 'A'] == pytest.approx(expected3)
    
    # 4. 多个 reaction 同一代谢物 - 应该 sum
    enzyme4 = pd.DataFrame({
        'metabolite': ['MultiRx', 'MultiRx'],
        'hmdb_id': ['HMDB003', 'HMDB003'],
        'reaction': ['R1', 'R2'],
        'gene': ['Gene1', 'Gene2'],
        'role': ['production', 'production']
    })
    
    adata4 = FakeAnnData(
        np.array([[5, 3], [5, 3]]),
        ['Gene1', 'Gene2'],
        {'cell_type': ['A', 'A']}
    )
    
    result4 = compute_metabolite_availability(adata4, enzyme4, min_cells=1)
    assert result4['P'].loc[('MultiRx', 'HMDB003'), 'A'] == pytest.approx(8.0)


def test_run_cell_mesh_with_availability():
    """测试使用 availability 方法的完整 CELL MESH 流程"""
    enzyme, sensor = load_cell_mesh_database()
    
    # 创建一个简单的测试数据集。显式选择一个有 production evidence 且
    # enzyme/sensor gene 不重复的 (metabolite, hmdb_id)，避免大表中的同名
    # gene 让 synthetic var_names 出现重复。
    enzyme_production = enzyme[enzyme["role"] == "production"]
    pairs = set(zip(enzyme_production["metabolite"], enzyme_production["hmdb_id"])).intersection(zip(sensor["metabolite"], sensor["hmdb_id"]))
    for met, hmdb_id in sorted(pairs):
        e_genes = enzyme_production.loc[
            (enzyme_production["metabolite"] == met) & (enzyme_production["hmdb_id"] == hmdb_id),
            "gene",
        ]
        s_genes = sensor.loc[
            (sensor["metabolite"] == met) & (sensor["hmdb_id"] == hmdb_id),
            "sensor_gene",
        ]
        candidates = [(e_gene, s_gene) for e_gene in e_genes for s_gene in s_genes if e_gene != s_gene]
        if candidates:
            e_gene, s_gene = candidates[0]
            break
    else:
        pytest.skip("No common metabolites with distinct production and sensor genes in packaged database")

    enzyme_subset = enzyme[
        (enzyme["metabolite"] == met)
        & (enzyme["hmdb_id"] == hmdb_id)
        & (enzyme["gene"] == e_gene)
    ]
    sensor_subset = sensor[
        (sensor["metabolite"] == met)
        & (sensor["hmdb_id"] == hmdb_id)
        & (sensor["sensor_gene"] == s_gene)
    ]

    genes = [e_gene, s_gene, "BACKGROUND"]
    X = np.array([
        [5, 0, 0], [4, 0, 0], [5, 0, 0],
        [0, 4, 0], [0, 5, 0], [0, 4, 0],
    ], dtype=float)
    
    adata = FakeAnnData(X, genes, {"cell_type": ["A", "A", "A", "B", "B", "B"]})
    
    res = run_cell_mesh(
        adata,
        enzyme_metabolite=enzyme_subset,
        metabolite_sensor=sensor_subset,
        cell_type_key="cell_type",
        min_cells=2,
        allow_self=False,
    )
    
    assert not res.events.empty
    assert "cell_mesh_score" in res.events.columns
    assert res.availability_results is not None
