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
    assert 'P_contrast' in result
    assert 'C_contrast' in result
    assert 'E_contrast' in result
    assert 'metadata' in result
    
    availability = result['availability']
    assert not availability.empty
    metA_avail = availability.loc[('MetA', 'HMDB001')]
    assert metA_avail.loc['A'] > metA_avail.loc['B']

    metc_idx = ('MetC', 'HMDB003')
    assert metc_idx in availability.index
    assert np.allclose(
        result['relative_consumption_support'].loc[metc_idx].values,
        0.0,
    )
    assert np.allclose(
        result['E_plus'].loc[metc_idx].values,
        0.0,
    )


def test_compute_availability_default_min_cells_is_100():
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
        }
    )

    with pytest.raises(ValueError, match=f"at least {MIN_CELL_COUNT} cells"):
        compute_metabolite_availability(adata, enzyme_metabolite)


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
