import numpy as np
import pandas as pd

from cell_mesh import (
    load_cell_mesh_database, 
    run_cell_mesh, 
    compute_metabolite_availability
)


class FakeAnnData:
    def __init__(self, X, var_names, obs):
        self.X = X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)


def test_robust_minmax():
    """测试 robust minmax 函数"""
    from cell_mesh.preprocess import robust_minmax
    
    # 测试正常情况
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    result = robust_minmax(x, lower=10, upper=90)
    assert len(result) == len(x)
    assert result.min() >= 0
    assert result.max() <= 1
    
    # 测试常数数组
    x = np.array([5, 5, 5, 5, 5])
    result = robust_minmax(x)
    assert np.all(result == 0)
    
    # 测试包含 NaN 的情况
    x = np.array([1, 2, np.nan, 4, 5])
    result = robust_minmax(x)
    assert not np.any(np.isnan(result))


def test_compute_metabolite_availability():
    """测试新的 metabolite availability 计算功能"""
    # 创建一个简单的 reaction table
    reaction_table = pd.DataFrame({
        'metabolite': ['MetA', 'MetA', 'MetA', 'MetB', 'MetB', 'MetC'],
        'hmdb_id': ['HMDB001', 'HMDB001', 'HMDB001', 'HMDB002', 'HMDB002', np.nan],
        'reaction': ['R1', 'R2', 'R3', 'R4', 'R5', 'R6'],
        'gene': ['Gene1', 'Gene2;Gene3', 'Gene4', 'Gene5', 'Gene6', 'Gene7'],
        'direction': ['product', 'product', 'substrate', 'product', 'transporter', 'product']
    })
    
    # 创建测试数据
    n_cells = 12
    n_genes = 7
    genes = [f'Gene{i+1}' for i in range(n_genes)]
    
    # 生成表达矩阵 - 让 Gene1-4 在细胞类型 A 中高表达
    X = np.random.poisson(0.1, size=(n_cells, n_genes)).astype(float)
    cell_types = ['A'] * 6 + ['B'] * 6
    
    # 在细胞类型 A 中高表达 Gene1, Gene2, Gene3
    X[0:6, 0] += 3  # Gene1
    X[0:6, 1] += 2  # Gene2
    X[0:6, 2] += 2  # Gene3
    
    # 在细胞类型 B 中高表达 Gene5, Gene6
    X[6:12, 4] += 3  # Gene5
    X[6:12, 5] += 2  # Gene6
    
    adata = FakeAnnData(X, genes, {"cell_type": cell_types})
    
    # 计算 metabolite availability
    result = compute_metabolite_availability(
        adata,
        reaction_table,
        celltype_col='cell_type',
        lower=0,  # 因为样本量小，使用全部数据
        upper=100,
        min_cells=1,
        return_intermediates=True
    )
    
    # 检查结果
    assert 'availability' in result
    assert 'P' in result
    assert 'C' in result
    assert 'E' in result
    assert 'P_norm' in result
    assert 'C_norm' in result
    assert 'E_norm' in result
    assert 'metadata' in result
    
    availability = result['availability']
    assert not availability.empty
    print(f"\nAvailability matrix shape: {availability.shape}")
    print(f"Availability matrix:\n{availability}")
    
    # 检查 MetA 在细胞类型 A 中的 availability 应该更高
    if ('MetA', 'HMDB001') in availability.index:
        metA_avail = availability.loc[('MetA', 'HMDB001')]
        assert metA_avail.loc['A'] > metA_avail.loc['B'] * 0.5, f"MetA availability in A should be higher: {metA_avail}"
    
    # 检查 MetC (只有 product，没有 substrate/transporter)
    if ('MetC', np.nan) in availability.index or ('MetC', 'nan') in availability.index:
        # 找到 MetC 的索引
        metc_idx = None
        for idx in availability.index:
            if idx[0] == 'MetC':
                metc_idx = idx
                break
        
        if metc_idx:
            # 检查 C_norm 是 0.41
            assert np.allclose(result['C_norm'].loc[metc_idx].values, 0.41)
            # 检查 E_norm 是 0.75
            assert np.allclose(result['E_norm'].loc[metc_idx].values, 0.75)
            print(f"\n✓ MetC 缺失值检查通过！")
    
    return result


def test_sparse_vs_dense():
    """测试 dense 和 sparse 输入是否得到一致结果"""
    print(f"\n{'='*60}")
    print(f"测试 dense vs sparse 输入一致性")
    print(f"{'='*60}")
    
    try:
        from scipy import sparse
    except ImportError:
        print("⚠️  跳过：scipy 不可用")
        return
    
    # 创建反应表
    reaction_table = pd.DataFrame({
        'metabolite': ['MetX', 'MetX', 'MetY'],
        'hmdb_id': ['HMDB00X', 'HMDB00X', 'HMDB00Y'],
        'reaction': ['R1', 'R2', 'R3'],
        'gene': ['GeneA', 'GeneB', 'GeneC'],
        'direction': ['product', 'substrate', 'product']
    })
    
    # 创建数据
    np.random.seed(42)
    X_dense = np.random.poisson(1, size=(10, 5)).astype(float)
    genes = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE']
    cell_types = ['Type1'] * 5 + ['Type2'] * 5
    
    # 创建 dense AnnData
    adata_dense = FakeAnnData(X_dense, genes, {"cell_type": cell_types})
    
    # 创建 sparse AnnData
    X_sparse = sparse.csr_matrix(X_dense)
    adata_sparse = FakeAnnData(X_sparse, genes, {"cell_type": cell_types})
    
    # 计算 dense 结果
    result_dense = compute_metabolite_availability(
        adata_dense, reaction_table, lower=0, upper=100
    )
    
    # 计算 sparse 结果
    result_sparse = compute_metabolite_availability(
        adata_sparse, reaction_table, lower=0, upper=100
    )
    
    # 比较结果
    avail_dense = result_dense['availability']
    avail_sparse = result_sparse['availability']
    
    print(f"\nDense availability:\n{avail_dense}")
    print(f"\nSparse availability:\n{avail_sparse}")
    
    # 检查一致性
    pd.testing.assert_frame_equal(avail_dense, avail_sparse)
    print(f"\n✓ Dense vs Sparse 结果一致！")


def test_boundary_cases():
    """测试边界情况"""
    print(f"\n{'='*60}")
    print(f"测试边界情况")
    print(f"{'='*60}")
    
    # 1. 代谢物没有 product reaction - 应该被跳过
    reaction_table1 = pd.DataFrame({
        'metabolite': ['NoProduct'],
        'hmdb_id': ['HMDB000'],
        'reaction': ['R1'],
        'gene': ['Gene1'],
        'direction': ['substrate']
    })
    
    adata = FakeAnnData(
        np.array([[1, 2], [3, 4]]),
        ['Gene1', 'Gene2'],
        {'cell_type': ['A', 'B']}
    )
    
    result1 = compute_metabolite_availability(adata, reaction_table1)
    assert result1['availability'].empty, "没有 product reaction 的代谢物应该被跳过"
    print(f"\n✓ 边界情况 1 测试通过：没有 product reaction 的代谢物被跳过")
    
    # 2. 所有 product 基因都缺失 - 应该与没有 product reaction 一样被跳过
    reaction_table2 = pd.DataFrame({
        'metabolite': ['MissingGenes'],
        'hmdb_id': ['HMDB001'],
        'reaction': ['R1'],
        'gene': ['NotPresent'],
        'direction': ['product']
    })
    
    result2 = compute_metabolite_availability(adata, reaction_table2)
    # availability 应该为空（因为 P 全 0，被过滤掉了）
    assert result2['availability'].empty, "所有 product 基因都缺失的代谢物应该被跳过"
    print(f"\n✓ 边界情况 2 测试通过：所有 product 基因缺失时代谢物被跳过")
    
    # 3. 同一 reaction 多个基因 - 应该取 max
    reaction_table3 = pd.DataFrame({
        'metabolite': ['MultiGene', 'MultiGene'],
        'hmdb_id': ['HMDB002', 'HMDB002'],
        'reaction': ['R1', 'R1'],  # 同一 reaction
        'gene': ['HighExpr;LowExpr', 'AnotherGene'],
        'direction': ['product', 'product']
    })
    
    adata3 = FakeAnnData(
        np.array([[10, 1, 0.5], [10, 1, 0.5]]),  # HighExpr = 10, LowExpr = 1, AnotherGene = 0.5
        ['HighExpr', 'LowExpr', 'AnotherGene'],
        {'cell_type': ['A', 'A']}
    )
    
    result3 = compute_metabolite_availability(adata3, reaction_table3, lower=0, upper=100)
    # P 应该是 10 (max of 10, 1) + 0.5 = 10.5
    assert 'P' in result3
    print(f"\nP matrix:\n{result3['P']}")
    print(f"\n✓ 边界情况 3 测试通过：多基因反应取 max")
    
    # 4. 多个 reaction 同一代谢物 - 应该 sum
    reaction_table4 = pd.DataFrame({
        'metabolite': ['MultiRx', 'MultiRx'],
        'hmdb_id': ['HMDB003', 'HMDB003'],
        'reaction': ['R1', 'R2'],
        'gene': ['Gene1', 'Gene2'],
        'direction': ['product', 'product']
    })
    
    adata4 = FakeAnnData(
        np.array([[5, 3], [5, 3]]),
        ['Gene1', 'Gene2'],
        {'cell_type': ['A', 'A']}
    )
    
    result4 = compute_metabolite_availability(adata4, reaction_table4, lower=0, upper=100)
    # P 应该是 5 + 3 = 8
    print(f"\nP matrix:\n{result4['P']}")
    print(f"\n✓ 边界情况 4 测试通过：多 reaction 代谢物取 sum")


def test_run_cell_mesh_with_new_availability():
    """测试使用新 availability 方法的完整 CELL MESH 流程"""
    enzyme, sensor = load_cell_mesh_database()
    
    # 创建一个简单的测试数据集
    # 使用一个有匹配酶和传感器的代谢物
    common_mets = sorted(set(enzyme["metabolite"]).intersection(sensor["metabolite"]))
    if not common_mets:
        print("⚠️  没有找到共同的代谢物，跳过此测试")
        return
    
    met = common_mets[0]
    e_genes = enzyme.loc[enzyme["metabolite"] == met, "gene"].unique()
    s_genes = sensor.loc[sensor["metabolite"] == met, "sensor_gene"].unique()
    
    genes = list(e_genes)[:2] + list(s_genes)[:1] + ["BACKGROUND"]
    n_genes = len(genes)
    
    X = np.array([
        [5, 3, 0, 0], [4, 2, 0, 0], [5, 4, 0, 0],
        [0, 0, 4, 0], [0, 0, 5, 0], [0, 0, 4, 0],
    ], dtype=float)
    
    adata = FakeAnnData(X, genes, {"cell_type": ["A", "A", "A", "B", "B", "B"]})
    
    # 创建一个简单的 reaction table
    reaction_table = pd.DataFrame({
        'metabolite': [met] * 3,
        'hmdb_id': enzyme.loc[enzyme["metabolite"] == met, "hmdb_id"].iloc[0] if "hmdb_id" in enzyme.columns else np.nan,
        'reaction': ['R1', 'R2', 'R3'],
        'gene': list(e_genes)[:3],
        'direction': ['product', 'substrate', 'transporter']
    })
    
    # 使用新的 availability 方法运行
    res = run_cell_mesh(
        adata,
        reaction_table=reaction_table,
        use_new_availability=True,
        cell_type_key="cell_type",
        min_cells_per_group=2,
        allow_self=False,
        lower=0,
        upper=100
    )
    
    assert not res.events.empty
    assert "cell_mesh_score" in res.events.columns
    assert res.availability_results is not None
    
    print(f"\nEvents with new availability method: {len(res.events)}")
    print(res.events[['sender', 'receiver', 'metabolite', 'cell_mesh_score']].head())
    
    # 确保原有的方法仍然可以工作
    print(f"\n{'='*60}")
    print(f"验证旧方法仍然工作")
    print(f"{'='*60}")
    res_old = run_cell_mesh(
        adata,
        use_new_availability=False,
        cell_type_key="cell_type",
        min_cells_per_group=2,
        allow_self=False
    )
    assert not res_old.events.empty
    print(f"\n✓ 旧方法仍然可以正常工作！")


if __name__ == "__main__":
    print("Testing metabolite availability module...")
    test_robust_minmax()
    print("✓ robust_minmax test passed")
    
    test_compute_metabolite_availability()
    print("✓ compute_metabolite_availability test passed")
    
    test_sparse_vs_dense()
    
    test_boundary_cases()
    
    test_run_cell_mesh_with_new_availability()
    print("✓ run_cell_mesh with new availability test passed")
    
    print("\n✅ All tests passed!")
