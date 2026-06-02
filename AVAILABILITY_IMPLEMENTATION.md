
# Metabolite Availability 模块实现总结 (v0.4.0)

## 📋 概述

CELL MESH v0.4.0 **完全基于 metabolite availability 算法**，公式已更新为：

```
Availability = P_norm × (1 - C_norm)^beta × (0.8 + 0.2 × E_norm)
```

其中：
- `P_norm`：标准化后的 production 分数（范围 [0, 1]）
- `C_norm`：标准化后的 consumption 分数（范围 [0, 1]）
- `E_norm`：标准化后的 efflux/transporter 分数（范围 [0, 1]）
- `beta = 0.5`（默认值）
- 最终结果严格限制在 `[0, 1]` 范围内

## 📁 修改的文件

1. **`cellmesh/config.py`** - 新增配置文件
2. **`cellmesh/preprocess.py`** - 主要实现重构
3. **`cellmesh/core.py`** - 核心算法重构
4. **`cellmesh/scoring.py`** - 简化评分函数
5. **`cellmesh/__init__.py`** - 更新导出
6. **新增文件**：
   - `smoke_test.py` - 冒烟测试
   - `TEST_REPORT.md` - 测试报告

## ✅ 功能特性

### 1. 集中配置管理

所有参数统一在 `config.py` 中管理：

```python
from cellmesh.config import METABOLITE_AVAILABILITY_DEFAULTS, ROLE_TO_DIRECTION
```

### 2. Role 到 Direction 映射

| Role | Direction | Matrix |
|------|-----------|--------|
| `production` | `product` | P |
| `degradation` | `substrate` | C |
| `export` | `exporter` | E |

### 3. 加权几何平均数

Reaction 分数计算使用加权几何平均数：
- 支持 `weight` 列
- 重复基因取最大权重并报警告
- 公式：`gmean(expr * weight + 1) - 1`

### 4. 多基因支持

支持多种分隔符：
- `;`（分号）
- `,`（逗号）
- `|`（竖线）

自动去除 `[Enzyme]`、`[Transport]` 等注释后缀。

## 🚀 使用方法

### 方法 1：单独计算 metabolite availability

```python
import pandas as pd
from cellmesh import compute_metabolite_availability

# 准备 enzyme_metabolite prior（不是 reaction table！）
enzyme_metabolite = pd.DataFrame({
    'metabolite': ['ATP', 'ATP', 'ATP', 'Glutamate'],
    'gene': ['Gene1', 'Gene2', 'Gene4', 'Gene6'],
    'role': ['production', 'degradation', 'export', 'production'],
    'weight': [1.0, 0.8, 1.2, 1.0]
})

# 转换为内部 reaction 表
from cellmesh.core import _enzyme_prior_to_availability_reactions
reaction_table = _enzyme_prior_to_availability_reactions(enzyme_metabolite)

# 计算
result = compute_metabolite_availability(
    adata,
    reaction_table,
    celltype_col='cell_type'
)

# 访问结果
print(result['availability'])  # 范围 [0, 1]
print(result['P_norm'])
print(result['C_norm'])
print(result['E_norm'])
```

### 方法 2：完整 CELL MESH 流程

```python
from cellmesh import run_cell_mesh

res = run_cell_mesh(
    adata,
    enzyme_metabolite=enzyme_metabolite,  # 直接使用 enzyme prior
    metabolite_sensor=metabolite_sensor,
    cell_type_key='cell_type',
    # availability 参数（可选，有默认值）
    lower=5,
    upper=95,
    eps=0.05,
    beta=0.5,
    missing_C_norm=0.2,
    missing_E_norm=0.5,
    min_cells=1
)

# res.sender_scores = availability
# res.availability_results 包含所有中间结果
```

## 📊 输出格式

`compute_metabolite_availability` 返回字典：

| 键 | 说明 |
|----|------|
| `availability` | 最终的 metabolite availability 矩阵（范围 [0, 1]） |
| `P` / `C` / `E` | 原始的 P/C/E 矩阵 |
| `P_norm` / `C_norm` / `E_norm` | 标准化后的矩阵（范围 [0, 1]） |
| `metadata` | 每个代谢物的元信息（包含 has_product / has_substrate / has_exporter） |
| `pseudobulk` | 细胞类型的平均表达矩阵 |
| `reaction_genes` | 解析后的反应-基因映射 |

## 🎯 默认参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `lower` | 5 | robust minmax 的下限百分位数 |
| `upper` | 95 | robust minmax 的上限百分位数 |
| `eps` | 0.05 | （保留但公式已简化） |
| `beta` | 0.5 | 消耗项的指数权重 |
| `missing_C_norm` | 0.2 | 缺失 substrate 时的默认值 |
| `missing_E_norm` | 0.5 | 缺失 exporter 时的默认值 |
| `min_cells` | 1 | 每个细胞类型的最小细胞数 |

## 🔍 处理的边界情况

- ✅ 代谢物没有 product reaction → 跳过（返回空结果）
- ✅ 所有 product 基因都缺失 → 跳过
- ✅ 没有 substrate reaction → 使用 `missing_C_norm`
- ✅ 没有 exporter reaction → 使用 `missing_E_norm`
- ✅ 同一 reaction 有多个基因 → 加权几何平均
- ✅ 稀疏矩阵输入 → 兼容
- ✅ NaN 值 → 正确处理
- ✅ 重复基因 → 取最大权重并报警告

## 🧪 测试

```bash
# 冒烟测试
python smoke_test.py

# 运行所有测试
cd tests
pytest
```

---

**更新日期**：2026年5月18日  
**版本**：0.4.0  
**状态**：✅ 完成
