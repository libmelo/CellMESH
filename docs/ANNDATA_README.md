# CELL MESH: AnnData 读取功能文档

## 概述

CELL MESH 现在支持多种格式的单细胞数据读取，使用统一的 `read_anndata()` 接口。

## 新功能

### 主要函数

#### `read_anndata(path, mode, **kwargs)`

主函数，支持多种读取模式。

**参数:**
- `path`: 文件路径或目录路径
- `mode`: 读取模式，支持以下选项：
  - `"h5ad"`: 从 .h5ad 文件读取 (Scanpy 格式)
  - `"10x"`: 从 10X Genomics 输出目录读取
  - `"csv"`: 从 CSV 文件读取 (细胞 × 基因矩阵)
  - `"tsv"`: 从 TSV 文件读取
  - `"loom"`: 从 Loom 文件读取
  - `"mtx"`: 从 Matrix Market 文件读取
- `**kwargs`: 传递给底层读取函数的额外参数

**返回:** AnnData 对象

#### `read_example_data(dataset)`

快速生成内置示例数据，用于测试。

**参数:**
- `dataset`: 数据集大小，可选:
  - `"tiny`: 50细胞 × 50基因
  - `"small`: 200细胞 × 100基因 (默认)
  - `"medium`: 500细胞 × 200基因

## 使用示例

### 1. 读取 .h5ad 文件

```python
import cell_mesh

adata = cell_mesh.read_anndata("path/to/data.h5ad", mode="h5ad")
```

### 2. 读取 10X Genomics 数据

```python
adata = cell_mesh.read_anndata("path/to/10x_directory", mode="10x")
```

### 3. 从 CSV 读取

```python
# 简单读取
adata = cell_mesh.read_anndata("expression.csv", mode="csv")

# 读取时同时加载元数据
adata = cell_mesh.read_anndata(
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
adata = cell_mesh.read_example_data(dataset="small")

# 读取中等大小数据
adata = cell_mesh.read_example_data(dataset="medium")
```

### 5. 完整工作流程

```python
import cell_mesh

# 步骤 1: 读取数据
adata = cell_mesh.read_example_data(dataset="small")

# 步骤 2: 运行 CELL MESH
result = cell_mesh.run_cell_mesh(
    adata,
    cell_type_key="cell_type",
    n_perms=100
)

# 步骤 3: 查看结果
print(result.events.head())
```

## 模式详细说明

### h5ad 模式

使用 Scanpy 的 AnnData 标准格式，最推荐的格式。

**需要:** `anndata` 包
**安装:** `pip install anndata`

### 10x 模式

读取 10X Genomics Cell Ranger 的输出目录。目录应包含：
- `matrix.mtx.gz` 或 `matrix.mtx`
- `genes.tsv.gz` / `genes.tsv` / `features.tsv.gz` / `features.tsv`
- `barcodes.tsv.gz` / `barcodes.tsv`

**需要:** `scanpy` 包
**安装:** `pip install scanpy`

### CSV/TSV 模式

读取文本格式的表达矩阵。

**CSV 额外参数:**
- `cell_meta_path`: 细胞元数据 CSV 文件
- `gene_meta_path`: 基因元数据 CSV 文件
- `cell_id_col`: 细胞 ID 列名
- `transpose`: 是否转置矩阵 (默认: False)

### Loom 模式

读取 Loom 格式文件。

**需要:** `anndata` 包

### mtx 模式

读取 Matrix Market 格式文件。

**MTX 额外参数:**
- `genes_path`: 基因列表文件路径
- `barcodes_path`: 细胞 barcode 文件路径

**需要:** `anndata` 和 `scipy` 包
**安装:** `pip install anndata scipy`

## 示例脚本

在 `examples/` 目录下提供完整的示例脚本：
- `examples/read_anndata_example.py` - 完整示例

运行示例：
```bash
cd /path/to/cell_mesh_pkg
python examples/read_anndata_example.py
```

## API 参考

### 已导出的函数

在 `__init__.py` 中已导出以下函数，可直接使用：

```python
from cell_mesh import (
    CellMeshResult,
    run_cell_mesh,
    run_metcomm,
    load_cell_mesh_database,
    load_default_priors,
    read_anndata,          # 新增
    read_example_data,     # 新增
)
```

## 注意事项

1. **依赖包**: 某些模式需要额外的包，安装提示会自动显示
2. **细胞类型**: 确保 `adata.obs` 中有 `cell_type` 列，或在运行 `run_cell_mesh()` 时指定 `cell_type_key`
3. **矩阵方向**: CSV/TSV 模式默认假设是细胞 × 基因，如果是基因 × 细胞，设置 `transpose=True`

## 更新日志

### v0.2.0 (新增)
- 添加 `read_anndata()` 函数，支持 6 种读取模式
- 添加 `read_example_data()` 函数，生成测试数据
- 更新 `__init__.py` 导出新函数
- 添加完整文档和示例
