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

### H5AD 的 backed 模式

```python
adata = cellmesh.read_anndata("data.h5ad", mode="h5ad", backed="r")
try:
    result = cellmesh.run_cell_mesh(
        adata, enzyme_prior, sensor_prior,
        sample_key="sample", sample_mode="sample_aware", min_cells=1,
    )
finally:
    adata.file.close()
```

`backed="r"` 以只读方式保持文件打开，主要表达矩阵按需读取；普通默认模式将
数据加载到内存。文件由调用者负责关闭，评分/绘图不会关闭或修改原文件。

评分前检查、独立评分、pseudobulk、表达比例、编译置换及两种小提琴图共用
兼容磁盘的切片规则：按递增、唯一位置读取基因列，再恢复请求顺序和重复项；
细胞选择按行块读取并恢复顺序。dense、CSR、CSC 以及重排细胞/基因后的 backed
视图均有回归测试。视图适配使用 AnnData 保存的视图到原对象的位置映射，避免按
可能重复或被改写的细胞名称重新匹配；升级 AnnData 后应运行这些兼容性测试。

sample-aware 对当前样本所选表达层构建内存对象，保留所有基因定义和 obs/var
对应关系，不额外复制 raw 或未选择的层。完整先验的反应去重规则保持不变。

**内存边界：** 支持 backed 输入不等于全程磁盘计算。分组汇总会加载当前细胞类型
的表达数据；sample-aware 会加载当前样本的所选矩阵；编译置换会保留所有细胞中
参与评分的基因子矩阵及计算中间量。稀疏矩阵保持稀疏，不会因兼容修订统一转为
稠密。AnnData 对 layers、obs/var 等内容的加载方式仍由其读取实现决定，不能假定
所有表达层都留在磁盘。内存需求取决于当前样本、细胞类型和评分基因的规模。

普通/backed 对照仅用于小型回归数据，日常分析不会自动运行两遍。测试比较事件
标识、得分、p/FDR、零与 NA，以及无置换/有置换、串行/并行和两种分析模式。

CSR/CSC 若含有重复的 `(cell, gene)` 存储坐标，计算会在稀疏副本中先提升为
float64 再合并，避免原始整数类型在比较或转为稠密时溢出。观测、编译置换及
小提琴取值共用该规则；backed CSC 分块拼接也保留原始条目供校验。原始相关
表达值先校验，负值不能通过重复坐标的正负抵消被隐藏。调用者的矩阵和文件不变，
没有重复坐标的规范矩阵无需此副本；程序不会统一将真实表达数据转换为 int8。

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

读取 10X Genomics Cell Ranger 的输出目录。三个文件均无表头，文本文件使用 UTF-8，
可带 BOM；各文件可独立选择未压缩或 gzip 形式：

| 文件 | 内容及规范 |
|---|---|
| `matrix.mtx` 或 `matrix.mtx.gz` | 基因/特征 × 细胞的 Matrix Market 矩阵，读入后转为细胞 × 特征 |
| legacy `genes.tsv` 或 `genes.tsv.gz` | 每条记录严格两列：基因 ID、基因符号 |
| modern `features.tsv` 或 `features.tsv.gz` | 每条记录严格三列：基因 ID、基因符号、特征类型 |
| `barcodes.tsv` 或 `barcodes.tsv.gz` | 每条记录严格一列：细胞 barcode |

目录中每类文件必须唯一。压缩和未压缩版本并存，或同时存在 `genes.tsv` 和
`features.tsv` 时会报告歧义；程序不猜测应采用哪一份。支持 `prefix="patient_"`
读取 `patient_matrix.mtx` 等带统一前缀的文件。注释记录数必须与原始矩阵维度一致，
空记录不跳过；格式错误报告文件及原始行位置。

**标识和选择参数：**

- `var_names="gene_symbols"`（默认）使用第二列基因符号作为 `adata.var_names`，
  第一列保存在 `adata.var["gene_ids"]`。选择 `var_names="gene_ids"` 时反过来，
  符号保存在 `adata.var["gene_symbols"]`。Enzyme 的 `gene` 和 Interaction 的
  `sensor_gene` 必须与所选标识体系一致；读取器不自动转换两种体系。
- 原始标识保持文本，`01`、`1`、`NA` 不合并。保留的基因轴和细胞轴必须非空，
  首尾空白去除后不能重复。重复基因错误同时列出原始 ID、符号及记录位置。
  不自动合并表达列或添加后缀；`make_unique=False` 可显式传入，`True` 会报错。
  原本合法且唯一的 `G-1` 等名称照常保留。
- `gex_only=True`（默认）在 modern 文件中只保留 `Gene Expression` 特征，
  特征类型保存在 `adata.var["feature_types"]`。基因重复检查作用于保留的特征；
  `gex_only=False` 时全部特征都必须满足所选标识的唯一性要求。legacy 文件不按类型筛选。
- `cache`、`cache_compression` 保留 Scanpy 的数值矩阵缓存行为；原始文件仍须存在，
  基因和 barcode 注释每次重新读取、校验。替换原始矩阵后应关闭或清理旧缓存。
  `gex_only`、`make_unique`、`cache` 接受布尔值，不接受字符串 `"False"`。

```python
adata = cellmesh.read_anndata("path/to/10x_directory", mode="10x")
# 仅在先验基因字段也使用这些 ID 时选择 gene_ids：
adata_by_id = cellmesh.read_anndata(
    "path/to/10x_directory", mode="10x", var_names="gene_ids", make_unique=False,
)
```

**需要:** `scanpy` 包
**安装:** `pip install -e ".[scanpy]"`

### CSV/TSV 模式

读取文本格式的表达矩阵。

**CSV 额外参数:**
- `cell_meta_path`: 细胞元数据 CSV 文件
- `gene_meta_path`: 基因元数据 CSV 文件
- `cell_id_col`: 细胞 ID 列名
- `transpose`: 是否转置矩阵 (默认: False)

**标识与重复检查：**

- 第一列为行标识，其余列为表达量；默认细胞 × 基因，`transpose=True` 时按基因 × 细胞解释。
  表头和最终轴标识均检查重复，首尾空格去除后重名也会报错，并指出名称及位置。
- 原始表头在 pandas 自动改名之前检查。`PROD,PROD` 会报错；原本就叫 `PROD.1` 的合法列保留。
  `usecols` 或显式 `names` 不会掩盖原始表头冲突。空的左上角表头允许存在，
  实际细胞/基因标识必须非空且唯一，不自动合并或添加编号。
- 细胞/基因标识从读取时保留文本，`01` 与 `1` 不同，字面 `NA`、`nan` 等不会被转换成缺失。
  表达矩阵正文仍使用数值解析；`dtype` 和表达量列的 `converters` 可以继续使用。
- 细胞和基因元数据的所有字段默认保留文本，因为后续分析可以选任意列作为样本或细胞类型。
  元数据空字段仍是缺失值。年龄、计数等附加字段需要时显式转为数值；例如：

```python
import pandas as pd

adata.obs["age"] = pd.to_numeric(adata.obs["age"], errors="raise")
```

元数据表头只允许首列留空，实际行标识必须非空；表头和行标识中的重复、
首尾空格后的重名同样拒绝。元数据按原始文本标识对齐；未匹配的记录仍保留缺失，
不会补成零。细胞/基因名中的首尾空格可由后续评分规范化，但导入时两份文件的对齐标识应一致。

**读取参数与兼容范围：**

CSV/TSV 使用 pandas 的 C 或 Python 解析器，支持单行表头（或 `header=None`）、
`names`、分隔符、引号、编码、压缩、注释、跳过行、`nrows`、`usecols` 等常用参数。
非表格前言可用 `skiprows` 跳过。完整读入后再建立 AnnData，不接受分块迭代结果、
多层表头、`parse_dates` 或标识列转换器；第一列固定用于行标识。
这些限制会明确报错，不会静默改变标识。

表达文本、细胞/基因元数据和 MTX 名称文件在 pandas 解析前检查实际 NUL 字符
（`0x00`），发现后报告文件和物理行号。检查按选定编码解码并支持压缩文件，
UTF-16/32 编码本身的零字节不会被误判；字面文本 `\0` 也不是 NUL。
整个文件都检查，包括被 `usecols`、`skiprows`、`nrows` 排除的内容，避免解析器
先截断字段，再将改写后的基因名、分组或表达值当作合法输入。

`dtype="Int64"`、`dtype="Float64"`、`dtype_backend="numpy_nullable"` 等 pandas
可空数值选项会在转置和建立 AnnData 前转换为 NumPy 数值列。无缺失整数列
保留整数类型，缺失整数使用浮点 NaN，缺失值不会补零。普通 NumPy float32/64
读取方式保持原行为；转换不用于掩盖非法文本、布尔或复数。实际参与评分的
缺失/非有限表达仍由评分校验拒绝。

**旧结果处理：**

旧版本可能把重复基因表头改成 `.1`，或把样本 `01`、`1` 都读成整数 `1`。
如果已发生这种信息丢失，应从原始文本文件重新读取并重跑分析；
对旧 AnnData 再执行 `astype(str)` 无法恢复原始名称。

### Loom 模式

需要核心依赖 `anndata` 以及可选依赖 `loompy`：

```bash
pip install 'cellmesh[loom]'
# 本地项目可使用：pip install -e '.[loom]'
```

```python
adata = cellmesh.read_anndata("data.loom", mode="loom", sparse=True)
```

Loom 使用 AnnData 对应版本的 reader，优先 `anndata.io.read_loom`，老版本回退
到原入口；额外参数原样传递。缺少 loompy 会给出安装提示，依赖内部错误、文件
损坏、文件不存在等保持原始异常，不会全部改写成安装依赖错误。真实 Loom 文件
的方向、基因/细胞标识、元数据和 dense/稀疏读取均有测试，并与同数据内存评分比较。
Loom 读取不会因本修订变为 backed 模式；backed 示例针对 H5AD。

### mtx 模式

读取 Matrix Market 格式文件，支持 `coordinate`（稀疏存储）和
`array`（稠密存储），以及 `.mtx.gz` 压缩文件。这里区分的是文件存储格式。

**MTX 额外参数:**
- `genes_path`: 基因列表文件路径
- `barcodes_path`: 细胞 barcode 文件路径

MTX 矩阵按基因 × 细胞读取，并在构建 AnnData 时转为细胞 × 基因。提供的
基因名或 barcode 数量必须与矩阵维度一致。两种格式沿用同一转置和名称
校验规则，保留原始表达值：`coordinate` 读取后的 `adata.X` 保持稀疏，
`array` 则保持二维 NumPy 数组，单个基因或单个细胞也不压缩维度。

名称文件无表头；基因文件有多列时使用第二列，barcode 使用第一列。
`.csv`（含 `.csv.gz` 等压缩形式）按逗号读取，其余名称文件按制表符读取，
单列文件不猜测分隔符。名称保留文本和前导零，字面 `NA` 等保持原值；
空名称、重复名称和去除首尾空格后重名都会明确报错。

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
`load_cell_mesh_database()` 自动选择，为 `Enzyme1.39.csv` 和
`Interaction1.41.csv`。

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
