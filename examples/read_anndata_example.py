#!/usr/bin/env python3
"""
CELL MESH: AnnData 读取功能示例脚本

展示如何使用多种格式读取单细胞数据
"""

import sys
from pathlib import Path

# 添加包路径
sys.path.insert(0, str(Path(__file__).parent.parent))

import cell_mesh


def example_1_read_example_data():
    """示例 1: 读取内置示例数据"""
    print("=" * 60)
    print("示例 1: 读取内置示例数据")
    print("=" * 60)
    
    # 读取小数据集
    print("\n--- 读取 small 示例数据 ---")
    adata = cell_mesh.read_example_data(dataset="small")
    
    print(f"\n数据统计:")
    print(f"  细胞数: {adata.n_obs}")
    print(f"  基因数: {adata.n_vars}")
    print(f"  细胞类型: {adata.obs['cell_type'].value_counts().to_dict()}")
    
    return adata


def example_2_read_from_csv():
    """示例 2: 从 CSV 文件读取"""
    print("\n" + "=" * 60)
    print("示例 2: 从 CSV 文件读取 (使用前先生成测试数据)")
    print("=" * 60)
    
    import tempfile
    import pandas as pd
    import numpy as np
    
    # 创建临时 CSV 文件
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # 生成示例数据
        rng = np.random.default_rng(42)
        n_cells, n_genes = 100, 50
        cell_types = rng.choice(["A", "B", "C"], size=n_cells)
        
        # 表达矩阵
        expr_matrix = pd.DataFrame(
            rng.poisson(0.1, size=(n_cells, n_genes)),
            index=[f"Cell{i+1}" for i in range(n_cells)],
            columns=[f"Gene{i+1}" for i in range(n_genes)]
        )
        
        # 细胞元数据
        cell_meta = pd.DataFrame({
            "cell_id": [f"Cell{i+1}" for i in range(n_cells)],
            "cell_type": cell_types,
            "sample": "Sample1"
        }, index=[f"Cell{i+1}" for i in range(n_cells)])
        
        # 保存文件
        expr_path = tmpdir / "expression.csv"
        meta_path = tmpdir / "cell_metadata.csv"
        
        expr_matrix.to_csv(expr_path)
        cell_meta.to_csv(meta_path)
        
        print(f"\n已生成临时文件:")
        print(f"  表达矩阵: {expr_path}")
        print(f"  细胞元数据: {meta_path}")
        
        # 读取 CSV
        print(f"\n正在使用 read_anndata 读取 CSV...")
        adata = cell_mesh.read_anndata(
            expr_path,
            mode="csv",
            cell_meta_path=meta_path,
            cell_id_col="cell_id"
        )
        
        print(f"\n成功读取!")
        print(f"  细胞数: {adata.n_obs}")
        print(f"  基因数: {adata.n_vars}")
        
        return adata


def example_3_comprehensive_workflow():
    """示例 3: 完整工作流程 - 读取数据 + 运行 CELL MESH"""
    print("\n" + "=" * 60)
    print("示例 3: 完整工作流程 - 读取 + 运行 CELL MESH")
    print("=" * 60)
    
    # 读取示例数据
    print("\n步骤 1: 读取示例数据")
    adata = cell_mesh.read_example_data(dataset="small")
    
    # 先为示例数据添加数据库中的基因，这样才能运行 CELL MESH
    print("\n步骤 2: 为示例数据添加数据库基因 (为了演示 CELL MESH)")
    
    # 读取数据库，获取基因名
    enzyme_db, sensor_db = cell_mesh.load_cell_mesh_database()
    db_genes = list(set(enzyme_db['gene']).union(set(sensor_db['sensor_gene'])))[:20]
    
    # 为示例数据添加这些基因
    import anndata
    import numpy as np
    import pandas as pd
    
    rng = np.random.default_rng(42)
    
    # 扩展基因
    all_genes = list(adata.var_names) + db_genes
    new_X = np.zeros((adata.n_obs, len(all_genes)))
    new_X[:, :adata.n_vars] = adata.X
    
    # 为新基因添加一些表达
    cell_types = adata.obs['cell_type'].unique()
    for i, ct in enumerate(cell_types):
        ct_idx = adata.obs['cell_type'] == ct
        gene_idx = adata.n_vars + i * 3 + np.arange(2)
        gene_idx = gene_idx[gene_idx < len(all_genes)]
        new_X[np.ix_(ct_idx, gene_idx)] += rng.poisson(1, size=(ct_idx.sum(), len(gene_idx)))
    
    # 创建新的 AnnData
    adata = anndata.AnnData(
        new_X,
        obs=adata.obs,
        var=pd.DataFrame(index=all_genes)
    )
    
    print(f"扩展后的数据: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
    
    # 运行 CELL MESH
    print("\n步骤 3: 运行 CELL MESH 分析")
    result = cell_mesh.run_cell_mesh(
        adata,
        cell_type_key="cell_type",
        n_perms=0,
        allow_self=False,
        min_expr_frac=0.01
    )
    
    print(f"\n分析完成!")
    print(f"检测到 {len(result.events)} 个通讯事件")
    
    if not result.events.empty:
        print(f"\nTop 5 事件:")
        print(result.events[
            ['sender', 'receiver', 'metabolite', 'sensor_gene', 'cell_mesh_score']
        ].head())
    
    return result


def main():
    print("CELL MESH: AnnData 读取功能示例")
    print("=" * 60)
    
    try:
        # 示例 1: 内置示例数据
        adata1 = example_1_read_example_data()
        
        # 示例 2: CSV 读取
        try:
            adata2 = example_2_read_from_csv()
        except Exception as e:
            print(f"\n示例 2 遇到问题 (可能缺少包): {e}")
            print("继续其他示例...")
        
        # 示例 3: 完整工作流
        result = example_3_comprehensive_workflow()
        
        print("\n" + "=" * 60)
        print("所有示例完成!")
        print("=" * 60)
        
    except ImportError as e:
        print(f"\n缺少必要的包: {e}")
        print("请安装: pip install anndata scanpy scipy")
    except Exception as e:
        print(f"\n示例运行出错: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
