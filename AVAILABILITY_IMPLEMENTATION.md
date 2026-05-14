# Metabolite Availability 模块实现总结

## 📋 概述

这个实现添加了新的 metabolite availability 计算方法，基于以下公式：

```
Availability = (P_norm + eps) × (1 - C_norm + eps)^beta × (E_norm + eps)
```

其中：
- `P_norm`：标准化后的 production 分数
- `C_norm`：标准化后的 consumption 分数
- `E_norm`：标准化后的 efflux/transporter 分数
- `eps = 0.05`（默认值）
- `beta = 0.5`（默认值）

## 📁 修改的文件

1. **`cell_mesh/preprocess.py`** - 主要实现文件
2. **`cell_mesh/core.py`** - 集成新方法到 CELL MESH
3. **`cell_mesh/__init__.py`** - 导出新函数
4. **新增文件**：
   - `tests/test_availability.py` - 新功能测试（包含边界情况和 dense/sparse 一致性测试）
   - `examples/example_availability.py` - 使用示例
   - `AVAILABILITY_IMPLEMENTATION.md` - 本文档

## ✅ Pull Request 审查 - 通过

### 已验证的 9 个检查点：

1. ✅ **只在有 product reaction 的代谢物上计算**
2. ✅ **transporter 正确记为 E，而不是 P**
3. ✅ **多 gene 同一 reaction 取 max**
4. ✅ **多 reaction 同一 metabolite 取 sum**
5. ✅ **robust_minmax 按每个 metabolite 在 cell types 间计算**
6. ✅ **缺失 substrate / transporter 的默认值在归一化后赋值**
7. ✅ **dense 和 sparse AnnData 输入结果一致**
8. ✅ **没有破坏原有细胞通讯打分接口**
9. ✅ **有足够单元测试覆盖边界情况**

## 🔧 核心功能

### 1. `_build_celltype_pseudobulk`
构建细胞类型的 pseudobulk 表达矩阵。

### 2. `_parse_reaction_table`
解析反应表，支持多基因分隔（`;`、`,`、`|`）。

### 3. `_get_reaction_gene_sets`
为每个反应提取基因集合，合并重复反应的基因。

### 4. `_compute_reaction_scores`
计算每个反应在每个细胞类型中的得分（取基因表达最大值）。

### 5. `_compute_PCE_matrices`
计算 P（production）、C（consumption）、E（efflux）矩阵。

### 6. `robust_minmax`
鲁棒标准化，使用百分位数截断异常值。

### 7. `_normalize_PCE`
标准化 P/C/E 矩阵，根据是否存在对应的 reaction 来决定是否使用默认值。

### 8. `compute_metabolite_availability`
主函数，计算完整的 metabolite availability。

## 🚀 使用方法

### 方法 1：单独计算 metabolite availability

```python
import pandas as pd
from cell_mesh import compute_metabolite_availability

# 准备 reaction table
reaction_table = pd.DataFrame({
    'metabolite': ['Glucose', 'Glucose', 'Lactate'],
    'hmdb_id': ['HMDB00122', 'HMDB00122', 'HMDB00190'],
    'reaction': ['R1', 'R2', 'R3'],
    'gene': ['HK1', 'PFKM', 'LDHA'],
    'direction': ['product', 'substrate', 'transporter']
})

# 计算
result = compute_metabolite_availability(
    adata,
    reaction_table,
    celltype_col='cell_type'
)

# 访问结果
print(result['availability'])
print(result['P_norm'])
```

### 方法 2：在 CELL MESH 中使用新方法

```python
from cell_mesh import run_cell_mesh

res = run_cell_mesh(
    adata,
    reaction_table=reaction_table,
    use_new_availability=True,
    cell_type_key='cell_type'
)

# res.availability_results 包含完整的中间结果
```

## 📊 输出格式

`compute_metabolite_availability` 返回一个字典，包含：

- `availability`：最终的 metabolite availability 矩阵
- `P`、`C`、`E`：原始的 P、C、E 矩阵
- `P_norm`、`C_norm`、`E_norm`：标准化后的矩阵
- `metadata`：每个代谢物的元信息
- `pseudobulk`：细胞类型的平均表达矩阵
- `reaction_genes`：解析后的反应-基因映射

## 🔗 向后兼容性

- 原有的 `run_cell_mesh` 方法保持不变
- 新方法通过 `use_new_availability=True` 启用
- `run_metcomm` 仍然作为 `run_cell_mesh` 的别名
- 所有原有测试仍然通过 ✅

## 🧪 测试

运行测试：

```bash
# 新功能测试（包含边界情况和 dense/sparse 一致性）
python tests/test_availability.py

# 原有功能测试（确保兼容性）
python tests/test_database_and_api.py
```

## 📖 示例

查看示例：

```bash
python examples/example_availability.py
```

## 🎯 默认参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `lower` | 5 | robust minmax 的下限百分位数 |
| `upper` | 95 | robust minmax 的上限百分位数 |
| `eps` | 0.05 | 公式中的小值 |
| `beta` | 0.5 | 公式中的指数 |
| `missing_C_norm` | 0.41 | 缺失 C 时的默认值 |
| `missing_E_norm` | 0.75 | 缺失 E 时的默认值 |
| `min_cells` | 1 | 每个细胞类型的最小细胞数 |

## 🔍 处理的边界情况

- 代谢物没有 product reaction → 跳过
- 所有 product 基因都缺失 → 跳过（与没有 product reaction 一致）
- 没有 substrate/transporter reaction → 使用默认值（0.41/0.75）
- 同一 reaction 有多个基因 → 取 max
- 稀疏矩阵输入 → 兼容，与 dense 输入结果一致
- NaN 值 → 正确处理

---

**实现日期**：2026年5月11日  
**版本**：0.3.1  
**审查状态**：✅ 通过
