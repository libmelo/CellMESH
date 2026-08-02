# CELL MESH 当前实现状态

本文件记录当前包的实现入口，避免继续保留早期开发任务中的旧包名和已删除接口。
完整算法说明见 `README.md`、`docs/METHODS.md` 和
`IMPLEMENTATION_SUMMARY.md`。

## 主要入口

- 包名：`cellmesh`
- 主分析函数：`cellmesh.run_cell_mesh()`
- 数据读取函数：`cellmesh.read_anndata()`
- 示例数据函数：`cellmesh.read_example_data()`
- 结果类型：`cellmesh.CellMeshResult`

历史文档中出现的 `cell_mesh` 包名和 `run_metcomm()` 别名均不是当前公开
API。

## AnnData 读取

读取实现位于 `cellmesh/io.py`，支持以下模式：

| 模式 | 输入 | 主要依赖 |
|---|---|---|
| `h5ad` | AnnData `.h5ad` 文件 | `anndata` |
| `10x` | 10X Genomics 矩阵目录 | `scanpy` 可选依赖 |
| `csv` | 细胞乘基因 CSV 矩阵 | `anndata`、`pandas` |
| `tsv` | 细胞乘基因 TSV 矩阵 | `anndata`、`pandas` |
| `loom` | Loom 文件 | `anndata` |
| `mtx` | Matrix Market 文件及可选名称文件 | `anndata`、`scipy` |

`read_example_data()` 支持 `tiny`、`small` 和 `medium`，默认值是
`tiny`。

```python
import cellmesh

adata = cellmesh.read_example_data("tiny")
print(adata.shape)
```

`read_example_data()` 生成的是读取接口测试数据，不保证其通用 `Gene1` 等基因
名称与默认生物学先验重叠。完整 CELL MESH 分析应使用包含先验基因的数据，
参见下列 workflow notebook。

更完整的参数说明和读取示例见 `docs/ANNDATA_README.md`；完整分析示例见
`examples/cellmesh_comprehensive_walkthrough.ipynb`、
`examples/demo_HNSC_cellmesh_workflow.ipynb` 和
`examples/demo_HNSC_sample_aware_workflow.ipynb`。

## 当前版本

当前包版本由 `cellmesh.__version__` 和 `pyproject.toml` 共同定义为
`0.4.0`。测试状态应以当前测试命令的实际输出为准，不在本状态文件中保存
容易过期的通过数量。
