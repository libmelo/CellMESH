# 常规输入与分析流程复审（2026-09-17）

**结论：本轮未确认新的待修订缺陷。** 审查对象为本地 `feature-modified` 工作树，包含 A01–A11、B01–B06 的既有修订。B07 按用户决定忽略，不进入待修订清单。本轮未修改生产代码或测试实现，只更新审查记录；复现脚本、合成数据与临时导出放在 `/tmp/cellmesh_review_regular_20260917/`。

## 范围

代码复核包括输入解析、先验规范化、pseudobulk/表达比例、生成可计算性、P/C/E 与受体评分、样本汇总、置换、FDR、CSV 导出和六种可视化。沿用已确认的模型规则：完整反应基因集合先去重、名称仅用于展示、实测零保留、结构性 NA 不补零、`min_cells` 仅作 QC、真实下溢报告后继续。

本轮新增运行验证重点为普通 scRNA-seq 计数矩阵和实际操作中的格式衔接，没有继续扩展 NumPy MaskedArray 等已决定忽略的输入范围。

## 新增对照结果

| 检查 | 本轮结果 |
|---|---|
| 表达载体与切片 | **48 组通过**：dense/CSR/CSC × int32/float32 × X/指定 layer × 两种推断模式 × 原对象/嵌套视图。视图含细胞反序、基因重排与细胞子集；与同一行列顺序的 float64 内存矩阵比较。 |
| 观测与置换 | 每组原始标签及同一随机序列的 7 次置换，各次均调用完整主流程重算，再与编译路径比较；合计 **384 个向量、18,432 个事件分数，最大差异 0**。样本内细胞类型计数保持不变。 |
| p 值 | 用完整重算的零分布独立累计尾部计数，沿用已记录的 100-epsilon 相对浮点平局容差与加一公式；**差异 0**。串/并行和保存/不保存零分布的运行对照通过。 |
| FDR | 以当前事件 p 值调用 SciPy 的 `false_discovery_control(method='bh')`，分别对照全局和按传感器类型分组的结果；**4,608 个值通过**。 |
| 缺失、零与 QC | 包含样本缺失端点、未使用的分类水平、零生成、生成先验缺失/基因未测得，以及未通过 QC 的细胞类型。可计算零事件保留且 p=1；不可计算生成不进入事件。样本中位数、可计算样本数和所存零分布均核对通过。 |
| 导出与绘图 | **17 项检查通过，覆盖全部六种函数**。事件计数图、网络图和 dotplot 对照内存结果与导出后重读表；两种小提琴核对单细胞数据；样本图确认缺失端点显示 NA。所有图均执行实际画布渲染。 |
| 输入不变 | 对照前后表达值和 obs 一致。选择 layer 时未将 X 中的无关负值带入计算。 |

这里的 48 组是有控制的组合对照，并非 48 套真实数据。编译路径与完整重算共用部分底层公式，因此这种一致性证据不能单独证明整个模型的生物学正确性；它验证的是当前实现路径、统计汇总和输出之间的对应关系。

## 回归测试

本轮实际重新执行以下 7 个测试文件，结果为 **220 passed、52 warnings，39.53 秒**。52 条均为既有 AnnData 索引转字符串提示，没有失败或跳过。

```bash
PYTHONPATH=/tmp/cellmesh-a10-deps MPLCONFIGDIR=/tmp/cellmesh-mpl NUMBA_CACHE_DIR=/tmp/cellmesh-numba \
  /home/qsong/miniconda3/envs/cellmesh/bin/python -m pytest -q \
  tests/test_compiled_permutation.py tests/test_permutation_precision.py \
  tests/test_permutation_evaluability.py tests/test_scoring_event_identity.py \
  tests/test_production_evaluability.py tests/test_export_consistency.py \
  tests/test_min_cells_qc.py
```

日志：`/tmp/cellmesh_review_regular_20260917_tests.log`。新增探针运行 63.18 秒，脚本为 `/tmp/cellmesh_review_regular_20260917/review_probe.py`，日志与结果 JSON 在同一目录。

上一轮同生产代码的全量基线为 **2493 passed、1 skipped、222 warnings**（`/tmp/b04-b06-full.log`）。本轮没有重新运行全量，不将它计为本次新增验证。

本轮证据、环境和源码 SHA-256 见 [机器可读记录](REVIEW_REGULAR_INPUTS_2026-09-17_EVIDENCE.json)。结论限于当前代码、本机环境和上述范围；未开展所有依赖版本兼容性、大数据性能或数据库生物学关系真实性的验证。数据库按既有安排保留本地。
