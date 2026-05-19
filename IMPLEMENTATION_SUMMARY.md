
# CELL MESH v0.4.0 - Implementation Summary

## 1. 概述

CELL MESH v0.4.0 进行了重大重构，**完全基于 metabolite availability 算法**，移除了所有旧算法代码。

**主要变更：**
- 版本号：`0.3.0` → `0.4.0`
- 完全移除旧算法，仅保留 metabolite availability
- 新增配置文件 `config.py`，集中管理参数
- 重构代码结构，提升可维护性
- 新增 `smoke_test.py` 冒烟测试脚本

## 2. 修改的文件

**新增文件：**
1. `cellmesh/config.py` - 集中配置管理
2. `smoke_test.py` - 冒烟测试脚本
3. `TEST_REPORT.md` - 测试报告

**修改的文件：**
1. `README.md` - 更新文档
2. `cellmesh/__init__.py` - 版本号和导出函数
3. `cellmesh/core.py` - 核心算法重构
4. `cellmesh/preprocess.py` - 预处理和 availability 计算重构
5. `cellmesh/scoring.py` - 简化评分函数
6. `pyproject.toml` - 包配置
7. `IMPLEMENTATION_SUMMARY.md` - 本文档

## 3. 新增配置文件 `config.py`

集中管理所有参数：

```python
# 基础路径配置
DATA_DIR = Path(__file__).parent / "data"

# 通用阈值配置
MIN_CELL_COUNT = 100
MIN_EXPR_FRAC = 0.05

# 代谢物 availability 默认参数
METABOLITE_AVAILABILITY_DEFAULTS = {
    "lower": 5,
    "upper": 95,
    "eps": 0.05,
    "beta": 0.5,
    "missing_C_norm": 0.2,
    "missing_E_norm": 0.5,
    "min_cells": 1,
}

# 角色到反应方向的映射
ROLE_TO_DIRECTION = {
    'production': 'product',      # 产生代谢物 → P 矩阵
    'degradation': 'substrate',    # 降解代谢物 → C 矩阵
    'usage': 'substrate',          # 使用代谢物 → C 矩阵
    'export': 'exporter',          # 外排代谢物 → E 矩阵
    'import': 'exporter',          # 转运代谢物 → E 矩阵
}
```

## 4. API 变更

### 4.1 `run_cell_mesh` 参数变更

**新增参数（v0.4.0）：**
- `lower` / `upper` / `eps` / `beta` / `missing_C_norm` / `missing_E_norm` / `min_cells` - availability 计算参数

**移除的参数（v0.4.0）：**
- `reaction_table` - 现在由 enzyme_metabolite 自动生成
- `use_new_availability` - 现在始终启用
- `role_agg` - 旧算法已移除
- `min_cells_per_group` - 替换为 `min_cells`
- `alpha_prod` / `alpha_deg` / `alpha_export` / `alpha_specificity` - 旧算法已移除

### 4.2 导出的函数

```python
from .core import CellMeshResult, run_cell_mesh
from .database import load_cell_mesh_database, load_default_priors
from .preprocess import compute_metabolite_availability, read_anndata, read_example_data
from .config import MIN_CELL_COUNT, DATA_DIR
```

## 5. 核心算法变更

### 5.1 Metabolite Availability 计算公式

```
availability = (P_norm) × ((1 - C_norm)^beta) × (0.8 + 0.2 × E_norm)
```

结果严格限制在 `[0, 1]` 范围内。

### 5.2 Role 到 Direction 的映射

| Role | Direction | Matrix |
|------|-----------|--------|
| `production` | `product` | P |
| `degradation` | `substrate` | C |
| `usage` | `substrate` | C |
| `export` | `exporter` | E |
| `import` | `exporter` | E |

### 5.3 Reaction 分数计算

- 使用**加权几何平均数**计算反应得分
- 支持权重列 `weight`
- 支持多基因分隔符（`;`、`,`、`|`）
- 重复基因取最大权重并报警告

## 6. 结果结构

```python
CellMeshResult(
    events=events_df,              # 通信事件
    sender_scores=sender_scores,   # 代谢物 × 细胞类型 availability
    receiver_scores=receiver_scores, # 接收方得分
    role_scores={},                # 空字典（向后兼容）
    parameters=parameters_dict,    # 参数
    availability_results=avail_results # 所有中间结果
)
```

## 7. 冒烟测试

新增 `smoke_test.py` 用于快速验证重构后的功能：

```bash
cd /home/qsong/.openclaw/workspace/developer/cell_mesh_pkg
python smoke_test.py
```

测试内容：
- 示例数据加载
- 所有 role 类型的处理
- availability 结果在 0-1 范围
- 参数正确性
- 元数据正确性

## 8. 向后兼容性

- `CellMeshResult.role_scores` 返回空字典（保留字段）
- 旧参数名抛出 KeyError（强制迁移）
- `run_metcomm` 别名仍然存在但已弃用

## 9. 包的最终结构

```
cell_mesh_pkg/
├── README.md                    # 主文档
├── IMPLEMENTATION_SUMMARY.md    # 本文档
├── AVAILABILITY_IMPLEMENTATION.md
├── TASK_COMPLETE.md
├── TEST_REPORT.md
├── smoke_test.py               # 冒烟测试
├── pyproject.toml              # 配置文件
├── cellmesh/                   # 主包
│   ├── __init__.py
│   ├── config.py               # 新增配置文件
│   ├── core.py                 # 核心算法
│   ├── preprocess.py           # 预处理和 availability
│   ├── database.py             # 数据库
│   ├── scoring.py              # 评分工具
│   └── data/                   # 内置数据
├── examples/                   # 示例
├── tests/                      # 测试
└── docs/                       # 文档
```

## 10. 总结

CELL MESH v0.4.0 是一个重大更新：
- ✅ 完全基于 metabolite availability 算法
- ✅ 代码更清晰、更易维护
- ✅ 集中配置管理
- ✅ 完整的测试覆盖
- ✅ 向后兼容性考虑

---

**更新日期**：2026年5月18日  
**版本**：0.4.0  
**状态**：✅ 完成
