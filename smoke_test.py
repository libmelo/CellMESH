"""
冒烟测试脚本，验证重构后的 CELL MESH 是否正常工作
"""
import numpy as np
import pandas as pd
import anndata
from cellmesh import run_cell_mesh, read_example_data

print("=" * 60)
print("CELL MESH 冒烟测试 - 新 Sensor 评分算法")
print("=" * 60)

# 1. 生成测试数据
print("\n1. 加载示例数据...")
adata = read_example_data("tiny")
print(f"  数据形状: {adata.shape}")
print(f"  细胞类型: {list(adata.obs['cell_type'].unique())}")
print(f"  基因名: {list(adata.var_names[:10])}...")

# 2. 构建测试 prior
print("\n2. 构建测试 prior 表...")

# 酶-代谢物 prior，只包含 3 个 roles
enzyme_metabolite = pd.DataFrame([
    # ATP 相关
    {"metabolite": "ATP", "hmdb_id": "HMDB00001", "gene": "Gene1", "role": "production", "weight": 1.0},
    {"metabolite": "ATP", "hmdb_id": "HMDB00001", "gene": "Gene2", "role": "degradation", "weight": 0.8},
    {"metabolite": "ATP", "hmdb_id": "HMDB00001", "gene": "Gene4", "role": "export", "weight": 1.2},
    # 谷氨酸相关
    {"metabolite": "Glutamate", "hmdb_id": "HMDB00002", "gene": "Gene6", "role": "production", "weight": 1.0},
    {"metabolite": "Glutamate", "hmdb_id": "HMDB00002", "gene": "Gene7", "role": "export", "weight": 1.0},
])

# 代谢物-传感器 prior，使用 3 种 sensor types
metabolite_sensor = pd.DataFrame([
    {"metabolite": "ATP", "hmdb_id": "HMDB00001", "sensor_gene": "Gene8", "sensor_type": "Cell surface receptor", "weight": 1.0},
    {"metabolite": "ATP", "hmdb_id": "HMDB00001", "sensor_gene": "Gene9", "sensor_type": "Transporter", "weight": 1.0},
    {"metabolite": "Glutamate", "hmdb_id": "HMDB00002", "sensor_gene": "Gene10", "sensor_type": "Other receptor", "weight": 1.0},
])

print(f"  酶 prior 行数: {len(enzyme_metabolite)}")
print(f"  传感器 prior 行数: {len(metabolite_sensor)}")
print(f"  传感器类型: {list(metabolite_sensor['sensor_type'].unique())}")

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
    display_cols = ['sender', 'receiver', 'metabolite', 'sensor_gene', 'sensor_type', 
                   'cell_mesh_score', 'sensor_expr_frac', 'perm_pvalue', 'fdr']
    available_cols = [c for c in display_cols if c in res.events.columns]
    print(res.events[available_cols].head())
    
    # 检查 score 范围
    min_score = res.events['cell_mesh_score'].min()
    max_score = res.events['cell_mesh_score'].max()
    print(f"  cell_mesh_score 范围: [{min_score:.4f}, {max_score:.4f}] (预期在 0-1 之间)")
    assert 0 <= min_score <= 1 and 0 <= max_score <= 1, "cell_mesh_score 超出 0-1 范围"
    
    # 检查 sensor_expr_frac 存在
    assert 'sensor_expr_frac' in res.events.columns, "缺少 sensor_expr_frac 列"
    assert 'sensor_type' in res.events.columns, "缺少 sensor_type 列"
    
    # 检查是否包含所需的 sensor types
    sensor_types = res.events['sensor_type'].unique()
    print(f"  检测到的 sensor types: {list(sensor_types)}")
    
    # 验证 cell_mesh_score 是几何均值
    # 检查一些样本
    sample = res.events.iloc[0]
    avail_score = sample['metabolite_availability']
    sensor_score = sample['sensor_score']
    expected = np.sqrt(avail_score * sensor_score)
    actual = sample['cell_mesh_score']
    print(f"  几何均值验证: sqrt({avail_score:.4f} * {sensor_score:.4f}) = {expected:.4f} (实际: {actual:.4f})")
    assert np.isclose(expected, actual), "cell_mesh_score 不是几何均值"

# 5. 检查参数是否正确
print("\n5. 参数检查:")
print(f"  包含的参数键: {list(res.parameters.keys())}")
# 确认没有旧参数
old_params = ['beta_sensor', 'beta_specificity']
for p in old_params:
    assert p not in res.parameters, f"旧参数 {p} 不应出现在结果中"
print("  ✓ 所有旧参数已成功移除")

# 6. 检查 receiver scores 格式
if not res.receiver_scores.empty:
    print("\n6. Receiver scores 检查:")
    print(f"  列名: {list(res.receiver_scores.columns)}")
    assert 'sensor_score' in res.receiver_scores.columns, "缺少 sensor_score 列"
    assert 'sensor_expr_frac' in res.receiver_scores.columns, "缺少 sensor_expr_frac 列"
    print("  ✓ receiver_scores 格式正确")

print("\n" + "=" * 60)
print("✅ 所有冒烟测试通过！新 Sensor 评分算法成功！")
print("=" * 60)
