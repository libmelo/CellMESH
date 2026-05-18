"""
冒烟测试脚本，验证重构后的 CELL MESH 是否正常工作
"""
import numpy as np
import pandas as pd
import anndata
from cell_mesh import run_cell_mesh, read_example_data

print("=" * 60)
print("CELL MESH 冒烟测试")
print("=" * 60)

# 1. 生成测试数据
print("\n1. 加载示例数据...")
adata = read_example_data("tiny")
print(f"  数据形状: {adata.shape}")
print(f"  细胞类型: {list(adata.obs['cell_type'].unique())}")

# 2. 构建测试 prior
print("\n2. 构建测试 prior 表...")

# 酶-代谢物 prior，包含所有 role 类型
enzyme_metabolite = pd.DataFrame([
    # ATP 相关
    {"metabolite": "ATP", "gene": "Gene1", "role": "production", "weight": 1.0},
    {"metabolite": "ATP", "gene": "Gene2", "role": "degradation", "weight": 0.8},
    {"metabolite": "ATP", "gene": "Gene3", "role": "usage", "weight": 0.7},
    {"metabolite": "ATP", "gene": "Gene4", "role": "export", "weight": 1.2},
    {"metabolite": "ATP", "gene": "Gene5", "role": "import", "weight": 0.9},
    # 谷氨酸相关
    {"metabolite": "Glutamate", "gene": "Gene6", "role": "production", "weight": 1.0},
    {"metabolite": "Glutamate", "gene": "Gene7", "role": "export", "weight": 1.0},
])

# 代谢物-传感器 prior
metabolite_sensor = pd.DataFrame([
    {"metabolite": "ATP", "sensor_gene": "Gene8", "sensor_type": "surface_receptor", "weight": 1.0},
    {"metabolite": "Glutamate", "sensor_gene": "Gene9", "sensor_type": "surface_receptor", "weight": 1.0},
])

print(f"  酶 prior 行数: {len(enzyme_metabolite)}")
print(f"  传感器 prior 行数: {len(metabolite_sensor)}")

# 3. 运行 CELL MESH，带置换检验
print("\n3. 运行 CELL MESH (n_perms=2)...")
res = run_cell_mesh(
    adata,
    enzyme_metabolite=enzyme_metabolite,
    metabolite_sensor=metabolite_sensor,
    cell_type_key="cell_type",
    n_perms=2,
    min_expr_frac=0.0,
    min_cells=1,
)

print(f"  事件总数: {len(res.events)}")
print(f"  sender_scores 形状: {res.sender_scores.shape}")
print(f"  receiver_scores 行数: {len(res.receiver_scores)}")
print(f"  availability_results 包含的键: {list(res.availability_results.keys())}")

# 4. 检查结果
print("\n4. 结果检查:")

# 检查 availability 结果是否在 0-1 之间
avail = res.availability_results['availability']
if not avail.empty:
    min_val = avail.min().min()
    max_val = avail.max().max()
    print(f"  availability 范围: [{min_val:.4f}, {max_val:.4f}] (预期在 0-1 之间)")
    assert 0 <= min_val <= 1 and 0 <= max_val <= 1, "availability 结果超出 0-1 范围"

# 检查事件结果
if not res.events.empty:
    print(f"  事件前 5 行:")
    print(res.events[['sender', 'receiver', 'metabolite', 'cell_mesh_score', 'perm_pvalue', 'fdr']].head())
    
    # 检查 score 范围
    min_score = res.events['cell_mesh_score'].min()
    max_score = res.events['cell_mesh_score'].max()
    print(f"  cell_mesh_score 范围: [{min_score:.4f}, {max_score:.4f}] (预期在 0-1 之间)")
    assert 0 <= min_score <= 1 and 0 <= max_score <= 1, "cell_mesh_score 超出 0-1 范围"

# 5. 检查参数是否正确
print("\n5. 参数检查:")
print(f"  包含的参数键: {list(res.parameters.keys())}")
# 确认没有旧参数
old_params = ['use_new_availability', 'role_agg', 'min_cells_per_group', 'alpha_prod', 'alpha_deg', 'alpha_export', 'alpha_specificity']
for p in old_params:
    assert p not in res.parameters, f"旧参数 {p} 不应出现在结果中"
print("  ✓ 所有旧参数已成功移除")

# 6. 检查元数据
if 'metadata' in res.availability_results and not res.availability_results['metadata'].empty:
    print("\n6. 元数据检查:")
    metadata = res.availability_results['metadata']
    print("  代谢物元数据:")
    print(metadata)
    # 检查 role 映射是否正确
    assert metadata.loc[('ATP', np.nan), 'has_substrate'] == True, "ATP 应该有 substrate 反应"
    assert metadata.loc[('ATP', np.nan), 'has_exporter'] == True, "ATP 应该有 exporter 反应"
    print("  ✓ role 映射正确")

print("\n" + "=" * 60)
print("✅ 所有冒烟测试通过！重构成功！")
print("=" * 60)
