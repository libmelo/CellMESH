"""
示例：使用新的 metabolite availability 方法

这个示例演示了如何：
1. 单独使用 compute_metabolite_availability 函数
2. 在 run_cell_mesh 中使用新方法
"""

import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '../')

from cell_mesh import (
    load_cell_mesh_database, 
    run_cell_mesh, 
    compute_metabolite_availability
)


class FakeAnnData:
    """简单的 AnnData 模拟类"""
    def __init__(self, X, var_names, obs):
        self.X = X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)


def example_1_compute_availability():
    """示例1：单独计算 metabolite availability"""
    print("=" * 60)
    print("示例1：单独计算 metabolite availability")
    print("=" * 60)
    
    # 创建一个示例 reaction table
    reaction_table = pd.DataFrame({
        'metabolite': ['Glucose', 'Glucose', 'Lactate', 'Lactate', 'ATP'],
        'hmdb_id': ['HMDB00122', 'HMDB00122', 'HMDB00190', 'HMDB00190', 'HMDB00538'],
        'reaction': ['Glycolysis1', 'Glycolysis2', 'LactateDehydro', 'MCT1', 'ATP_synthase'],
        'gene': ['HK1;HK2', 'PFKM', 'LDHA', 'SLC16A1', 'ATP5A1'],
        'direction': ['product', 'product', 'product', 'transporter', 'product']
    })
    
    print("\nReaction table:")
    print(reaction_table)
    
    # 创建测试数据
    n_cells = 15
    n_genes = 5
    genes = ['HK1', 'HK2', 'PFKM', 'LDHA', 'SLC16A1', 'ATP5A1'][:n_genes]
    
    # 生成表达矩阵 - 模拟两种细胞类型
    np.random.seed(42)
    X = np.random.poisson(0.5, size=(n_cells, n_genes)).astype(float)
    cell_types = ['Hepatocyte'] * 7 + ['Myocyte'] * 8
    
    # Hepatocytes 高表达糖酵解基因
    X[0:7, 0] += 4  # HK1
    X[0:7, 1] += 3  # HK2
    X[0:7, 2] += 5  # PFKM
    
    # Myocytes 高表达乳酸相关基因
    X[7:15, 3] += 4  # LDHA
    X[7:15, 4] += 3  # SLC16A1
    
    adata = FakeAnnData(X, genes, {"cell_type": cell_types})
    
    print(f"\nAnnData shape: {adata.X.shape}")
    print(f"Cell types: {cell_types}")
    
    # 计算 metabolite availability
    result = compute_metabolite_availability(
        adata,
        reaction_table,
        celltype_col='cell_type',
        lower=0,  # 因为样本量小
        upper=100,
        min_cells=1,
        return_intermediates=True
    )
    
    print(f"\n✅ Metabolite availability 计算完成")
    
    # 显示结果
    availability = result['availability']
    print(f"\nAvailability matrix (metabolites × cell types):")
    print(availability)
    
    print(f"\nP (production) scores:")
    print(result['P'])
    
    print(f"\nMetadata:")
    print(result['metadata'])
    
    return result


def example_2_run_cellmesh_with_new_availability():
    """示例2：在 CELL MESH 中使用新方法"""
    print("\n" + "=" * 60)
    print("示例2：在 CELL MESH 中使用新方法")
    print("=" * 60)
    
    # 加载数据库
    enzyme, sensor = load_cell_mesh_database()
    print(f"\nLoaded enzyme table: {len(enzyme)} rows")
    print(f"Loaded sensor table: {len(sensor)} rows")
    
    # 找到一个有匹配酶和传感器的代谢物
    common_mets = sorted(set(enzyme["metabolite"]).intersection(sensor["metabolite"]))
    if not common_mets:
        print("❌ 没有找到共同的代谢物")
        return
    
    met = common_mets[0]
    print(f"\n使用代谢物: {met}")
    
    # 获取相关基因
    e_genes = enzyme.loc[enzyme["metabolite"] == met, "gene"].unique()[:3]
    s_genes = sensor.loc[sensor["metabolite"] == met, "sensor_gene"].unique()[:2]
    
    genes = list(e_genes) + list(s_genes) + ["BACKGROUND"]
    n_genes = len(genes)
    
    print(f"酶基因: {list(e_genes)}")
    print(f"传感器基因: {list(s_genes)}")
    
    # 创建测试数据
    X = np.array([
        [5, 3, 2, 0, 0, 0], 
        [4, 2, 3, 0, 0, 0], 
        [5, 4, 1, 0, 0, 0],
        [0, 0, 0, 4, 3, 0], 
        [0, 0, 0, 5, 2, 0], 
        [0, 0, 0, 4, 3, 0],
    ], dtype=float)
    
    adata = FakeAnnData(X, genes, {"cell_type": ["Sender", "Sender", "Sender", "Receiver", "Receiver", "Receiver"]})
    
    # 创建 reaction table
    reaction_table = pd.DataFrame({
        'metabolite': [met] * len(e_genes),
        'hmdb_id': enzyme.loc[enzyme["metabolite"] == met, "hmdb_id"].iloc[0] if "hmdb_id" in enzyme.columns else np.nan,
        'reaction': [f"R{i+1}" for i in range(len(e_genes))],
        'gene': list(e_genes),
        'direction': ['product'] * len(e_genes)
    })
    
    print(f"\nReaction table:")
    print(reaction_table)
    
    # 使用新方法运行 CELL MESH
    print(f"\n运行 CELL MESH (新方法)...")
    res = run_cell_mesh(
        adata,
        enzyme_metabolite=enzyme,
        metabolite_sensor=sensor,
        reaction_table=reaction_table,
        use_new_availability=True,
        cell_type_key="cell_type",
        min_cells_per_group=2,
        allow_self=False,
        lower=0,
        upper=100
    )
    
    print(f"\n✅ CELL MESH 完成")
    print(f"\n找到 {len(res.events)} 个事件")
    
    if not res.events.empty:
        print(f"\nTop 5 events:")
        print(res.events[['sender', 'receiver', 'metabolite', 'cell_mesh_score', 'confidence_tier']].head())
    
    # 也用旧方法运行，进行对比
    print(f"\n" + "-" * 60)
    print(f"对比：用旧方法运行 CELL MESH")
    print("-" * 60)
    
    res_old = run_cell_mesh(
        adata,
        enzyme_metabolite=enzyme,
        metabolite_sensor=sensor,
        use_new_availability=False,
        cell_type_key="cell_type",
        min_cells_per_group=2,
        allow_self=False
    )
    
    print(f"\n旧方法找到 {len(res_old.events)} 个事件")


if __name__ == "__main__":
    print("🧪 CELL MESH - Metabolite Availability 示例")
    print("=" * 60)
    
    # 运行示例
    example_1_compute_availability()
    example_2_run_cellmesh_with_new_availability()
    
    print("\n" + "=" * 60)
    print("🎉 所有示例完成！")
    print("=" * 60)
