# CELL MESH: AnnData 读取功能文档

## 概述

CELL MESH 现在支持多种格式的单细胞数据读取，使用统一的 `read_anndata()` 接口。

## 接口

### 主要函数

#### `read_anndata(path, mode, **kwargs)`

主函数，支持多种读取模式。`mode` 默认是 `"h5ad"`。

**参数:**
- `path`: 文件路径或目录路径
- `mode`: 读取模式，支持以下选项：
  - `"h5ad"`: 从 .h5ad 文件读取 (AnnData/Scanpy 格式)
  - `"10x"`: 从 10X Genomics 输出目录读取
  - `"csv"`: 从 CSV 文件读取 (细胞 × 基因矩阵)
  - `"tsv"`: 从 TSV 文件读取
  - `"loom"`: 从 Loom 文件读取
  - `"mtx"`: 从 Matrix Market 文件读取
- `**kwargs`: 传递给对应读取函数的参数。`mtx` 模式只接受
  `genes_path` 和 `barcodes_path`，未知参数会直接报错，不会被静默忽略。

**返回:** AnnData 对象

#### `read_example_data(dataset)`

快速生成内置示例数据，用于测试。

**参数:**
- `dataset`: 数据集大小，可选：
  - `"tiny"`: 50 细胞 × 50 基因（默认）
  - `"small"`: 200 细胞 × 100 基因
  - `"medium"`: 500 细胞 × 200 基因

## 使用示例

### 1. 读取 .h5ad 文件

```python
import cellmesh

adata = cellmesh.read_anndata("path/to/data.h5ad", mode="h5ad")
```

### 2. 读取 10X Genomics 数据

```python
adata = cellmesh.read_anndata("path/to/10x_directory", mode="10x")
```

### 3. 从 CSV 读取

```python
# 简单读取
adata = cellmesh.read_anndata("expression.csv", mode="csv")

# 读取时同时加载元数据
adata = cellmesh.read_anndata(
    "expression.csv",
    mode="csv",
    cell_meta_path="cell_metadata.csv",
    gene_meta_path="gene_metadata.csv",
    cell_id_col="cell_id"
)
```

### 4. 使用内置示例数据

```python
# 读取小示例数据
adata = cellmesh.read_example_data(dataset="small")

# 读取中等大小数据
adata = cellmesh.read_example_data(dataset="medium")
```

### 5. 完整工作流程

```python
import cellmesh

# 步骤 1: 读取 walkthrough 的固定表达矩阵
adata = cellmesh.read_anndata(
    cellmesh.DATA_DIR / "test_single_cell.h5ad",
    mode="h5ad",
)
# 该固定旧测试文件把原始基因名保存在 var["_index"] 中。
if "_index" in adata.var.columns:
    adata.var_names = adata.var["_index"].astype(str)

# 步骤 2: 加载与该表达矩阵配套的固定先验
enzyme_prior, sensor_prior = cellmesh.load_cell_mesh_database(
    enzyme_file=cellmesh.DATA_DIR / "Enzyme_new.csv",
    interaction_file=cellmesh.DATA_DIR / "Interaction1.0.csv",
)

# 步骤 3: 运行 CELL MESH
result = cellmesh.run_cell_mesh(
    adata,
    enzyme_metabolite=enzyme_prior,
    metabolite_sensor=sensor_prior,
    cell_type_key="cell_type",
    n_perms=0,
)

# 步骤 4: 查看结果
print(result.events.head())
```

`read_example_data()` 使用通用的 `Gene1`、`Gene2` 等名称，适合验证 AnnData
读取和对象结构，但不保证与默认生物学先验存在基因交集。

## 模式详细说明

### h5ad 模式

读取 AnnData 的 `.h5ad` 标准格式。

**需要:** `anndata` 包，已包含在 CELL MESH 核心依赖中。
**安装:** `pip install -e .`

### 10x 模式

读取 10X Genomics Cell Ranger 的输出目录。目录应包含：
- `matrix.mtx.gz` 或 `matrix.mtx`
- `genes.tsv.gz` / `genes.tsv` / `features.tsv.gz` / `features.tsv`
- `barcodes.tsv.gz` / `barcodes.tsv`

**需要:** `scanpy` 包
**安装:** `pip install -e ".[scanpy]"`

### CSV/TSV 模式

读取文本格式的表达矩阵。

**CSV 额外参数:**
- `cell_meta_path`: 细胞元数据 CSV 文件
- `gene_meta_path`: 基因元数据 CSV 文件
- `cell_id_col`: 细胞 ID 列名
- `transpose`: 是否转置矩阵 (默认: False)

### Loom 模式

读取 Loom 格式文件。

**需要:** `anndata` 包，已包含在 CELL MESH 核心依赖中。

### mtx 模式

读取 Matrix Market 格式文件。

**MTX 额外参数:**
- `genes_path`: 基因列表文件路径
- `barcodes_path`: 细胞 barcode 文件路径

MTX 矩阵按基因 × 细胞读取，并在构建 AnnData 时转为细胞 × 基因。提供的
基因名或 barcode 数量必须与矩阵维度一致。

**需要:** `anndata` 和 `scipy` 包，均已包含在 CELL MESH 核心依赖中。
**安装:** `pip install -e .`

## 示例 notebook

在 `examples/` 目录下提供完整 walkthrough：
- `examples/cellmesh_comprehensive_walkthrough.ipynb`

运行前安装 notebook 依赖：
```bash
pip install -e ".[notebook]"
```

该 notebook 为保证示例计算可复现，使用包内固定测试文件：
- `cellmesh/data/test_single_cell.h5ad`
- `cellmesh/data/Enzyme_new.csv`
- `cellmesh/data/Interaction1.0.csv`

这些是 walkthrough 固定输入，不代表默认数据库版本。当前默认数据库由
`load_cell_mesh_database()` 自动选择，为 `Enzyme2.0.csv` 和
`Interaction4.0.csv`。

## API 参考

### 已导出的函数

在 `__init__.py` 中已导出以下函数，可直接使用：

```python
from cellmesh import (
    CellMeshResult,
    run_cell_mesh,
    load_cell_mesh_database,
    load_default_priors,
    read_anndata,
    read_example_data,
)
```

## 注意事项

1. **依赖包**: 某些模式需要额外的包，安装提示会自动显示
2. **细胞类型**: 确保 `adata.obs` 中有 `cell_type` 列，或在运行 `run_cell_mesh()` 时指定 `cell_type_key`
3. **矩阵方向**: CSV/TSV 模式默认假设是细胞 × 基因，如果是基因 × 细胞，设置 `transpose=True`

## 版本说明

本文档对应 CELL MESH `0.4.0` 的公开读取接口。包的完整当前 API 和算法
说明以根目录 `README.md` 及 `docs/METHODS.md` 为准。
