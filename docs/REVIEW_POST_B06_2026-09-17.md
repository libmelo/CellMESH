# B01–B06 修订后的审查（2026-09-17）

本轮对象为本地 `feature-modified` 当前工作树。B01–B06 已修订；此次仅检查、复现并记录，确认新增 **B07 一项**。既有 A01–A11 与 B01–B06 记录保持各自范围。

**后续决定（2026-09-17）：用户确认忽略 B07，不修订，也不列入待修订清单。** NumPy MaskedArray 属于特殊输入，mask 并非常规 scRNA-seq 缺失标记；没有证据表明实际数据受此问题影响。下文保留当时的复现及方案作为历史记录，其中 P1 表示触发后的数值影响，不代表当前项目的修订优先级。后续审查优先覆盖常规 scRNA-seq 输入与实际分析路径。

## 验证范围与结果

- 复核表达读取与校验、分组均值/表达比例、独立评分、观测/编译置换、样本汇总、FDR、导出和绘图之间的衔接，重点补充此前未覆盖的 NumPy 掩码数组。
- **40 组完整路径对照全部通过**：CSV、可空数值 CSV、gzip TSV、MTX、backed CSC 五种入口 × 两种推断模式 × 无置换/5 次置换 × 无 gate/0.75 gate。与同数据原生 AnnData 比较事件得分与 p 值；检查 CSV 导出、manifest、文本标识重新读取和计数图渲染。包含 `01`、`1`、字面 `NA`、逗号、引号、换行、缺失样本端点及 QC 未通过事件。
- 普通数组与掩码数组另做同序列置换对照；掩码数组复现 B07。普通数组在两种模式、X/指定 layer 下均无得分或 p 值差异。
- 引用上一轮同生产代码工作树的全量基线：**2493 passed、1 skipped、222 warnings，262.92 秒**，日志 `/tmp/b04-b06-full.log`；本轮未重新运行全量测试，没有把基线作为新增边界已通过的证据。
- 环境为 Python 3.10.20、NumPy 2.2.6、pandas 2.3.3、SciPy 1.15.3、AnnData 0.11.4、Matplotlib 3.10.9。结论限于已检查代码和上述输入；未声称穷尽所有依赖版本、后端及数据组合。

源码 SHA-256、环境、探针结果见 [机器可读证据](REVIEW_POST_B06_2026-09-17_EVIDENCE.json)。所有复现数据均为合成数据，脚本位于本机 `/tmp/cellmesh_review_post_b06_20260917/`。

## B07：掩码表达的缺失标记被丢弃，底层值进入评分与绘图（P1）

**触发条件。** 调用者将 `numpy.ma.MaskedArray` 放入 AnnData 的 X 或所选 layer，且参与评分/绘图的位置被标记为 masked，底层仍保存有限非负数值。`mask=True` 表示该位置不可用，不代表真实零表达。常规 dense、CSR、CSC 输入不会自动产生这种掩码；本问题不等同于它们的显式零或 NaN。

**位置与原因。**

- [preprocess.py](../cellmesh/preprocess.py) `_validate_expression_values()` 第 95 行先 `np.asarray(X)`，导致缺失标记消失，只检查底层数据。
- 同文件 `_grouped_expression_mean()` 使用稀疏成员矩阵乘表达矩阵，在本环境中底层值进入均值；观测表达比例的 masked reduction 又会跳过掩码位置。这使同一数据的均值与表达比例使用不同口径。
- [permutation.py](../cellmesh/permutation.py) `_pseudobulk()` 的稀疏矩阵乘比较结果也忽略掩码，置换的表达比例与观测不同。
- [plotting.py](../cellmesh/plotting.py) `_expression_columns()` 第 2041 行再次 `np.asarray(selected, dtype=float)`，使两种小提琴图画出被掩码的底层值。
- [preprocess.py](../cellmesh/preprocess.py) `_expression_source()` 第 125–130 行只对 backed X 视图特殊处理。在当前 AnnData 中，普通视图的 `.X` 和 `.layers[...]` 可先变成没有掩码的 `ArrayView`；只修改数值校验函数不足以覆盖视图。sample-aware 的普通对象复制路径对 masked layer 也可丢失掩码。

**合成复现。** 四个细胞、两种细胞类型 A/A/B/B、一个样本，表达基因为生成基因 G 和受体 R：

```python
expression = numpy.ma.array(
    [[1., 20.], [1., 2.], [4., 3.], [4., 8.]],
    mask=[[False, True], [False, False], [False, False], [False, False]],
)
```

A 的第一个 R 值为缺失，20 只是被掩码的底层值。当前流程仍计算 A 的 R 均值为 `(20+2)/2 = 11`；观测表达比例为 0.5，但编译路径为 1。设置 `min_expr_frac=0.75` 后，即使不改变标签，同一事件也可在观测中为 0、在编译路径中为约 0.562569。

使用 `n_perms=19, random_state=2`，在完全相同的标签置换序列下：

| 事件 | 每次调用完整观测路径重算后累计的 p | 当前编译置换得到的 p |
|---|---:|---:|
| B → B，HMDB1 / R | 0.15 | 0.35 |
| A → B，HMDB1 / R | 0.55 | 1.00 |

这组差异在 pooled 的 X/layer 和 sample-aware 的 X 均复现。sample-aware 的 masked layer 在复制时已丢标记，两条路径可能一致地读取底层值，因此“对照一致”不能证明缺失处理正确。

以上 p 值仅用于证明计算路径不一致，**两列均不能作为包含未处理掩码缺失的数据的有效推断结果**。修订应在计算前拒绝相关缺失，不应把其中一列当作正确结果来逼近。

独立 availability 也复现底层值泄漏：保持相同 mask 和可见表达，只将被掩码的 G 底层值从 1 改成 1000，A 的 sender availability 从约 **0.257143** 变为 **0.598406**。两种小提琴图会将被掩码的底层 20 作为真实 production score 或 receptor expression 绘制。

**影响范围。** 主入口的两种推断模式、独立 sender/receiver 评分、编译置换及两种表达小提琴图；可能改变 P/C/E、生成状态、参考值、表达 gate、事件分数和 p/FDR。只有无关基因被掩码的对照没有改变本次相关事件，但修订仍须保持“只检查所选层和实际相关表达”的既有规则。

**建议修订方案。**

1. 在任何 `np.asarray()`、分组、复制或点积前检查相关表达的 mask。发现 `mask=True` 即明确报错，指出来源层及缺失位置；不按底层值或补零继续，也不临时改变统计分母。
2. 普通 AnnData 视图按原对象的行列位置读取所选 X/layer，保留掩码直到校验完成；覆盖行/列重排和切片。不能只检查已经丢掉 mask 的 `ArrayView`。
3. 主入口、独立评分、编译构建和两种小提琴共用这一检查。全 False 的 mask 可以按普通数组计算；无关基因或未选择层的掩码不应阻断相关计算。
4. 回归覆盖 X/layer、原对象/视图、两种推断模式、无/有置换、底层有限值/负值/NaN、全 False mask、无关位置，以及输入对象未被修改。有效普通输入与既有结果保持一致。

本轮未确认第二项独立的新缺陷。常规完整读取、计算、导出和重绘对照已通过，不需要因本次掩码发现推定所有现有数据均受影响。B07 按上述用户决定忽略，代码保持现状。
