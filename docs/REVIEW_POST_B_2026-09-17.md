# B01–B03 修订后的补充审查（2026-09-17）

对象：本地 `feature-modified` 当前工作树，包含此前 A01–A11 及 B01–B03 修订。初次审查检查生产代码、运行独立复现并记录结论；随后按用户授权依次修订 B04–B06。以下问题描述与机器证据保留修订前状态，实施与验收结果见文末。

## 验证范围

- 阅读输入解析、先验规范化、独立评分、主入口、分样本汇总、编译置换、FDR、结果导出和六种绘图函数；复核 B01–B03 的共用路径。
- 新随机种子的 **64 组配置**：比较逐次完整重算与编译置换，合计 **384 个向量、18,360 个事件分数**。覆盖 dense/CSR/CSC、float32/64、两种推断模式、样本缺失端点、零表达/未测基因、完整反应去重、选择 layer、名称空白、mean/median 参考值、表达比例 gate、丰度指数、外排权重及串/并行。增加 `1e-200` 至 `1e100` 的跨基因量级差异。最大绝对差 **1.1102230246251565e-16**，**p 值差异为 0**。真正下溢仍给出诊断并继续。
- CSV/TSV、C/Python 引擎及 transpose 对照中，包含逗号、引号、字段内换行、中文、空白、`01`/`1`/`NA` 的合法标识均保持原值。重复 DataFrame 行索引、float32/nullable Float64 绘图，以及 categorical/string QC 的专项探针未发现新故障。
- 上一轮同一生产代码工作树的全量测试为 **2325 passed、1 skipped、213 warnings，249.65 秒**；本轮没有修改生产代码，因此引用该基线，没有把它记为重新运行。此次新发现的边界不在原有测试中。
- 当前环境：Python 3.10.20、NumPy 2.2.6、pandas 2.3.3、SciPy 1.15.3、AnnData 0.11.4、Matplotlib 3.10.9。没有验证所有依赖版本、全部可选解析器、大数据性能或数据库的生物学关系真实性。

本轮确认以下 **3 项**。它们均有明确触发条件；本轮未复现常规输入下新的主评分或置换 p 值不一致。源码 SHA-256、环境和合成探针结果见 [机器可读证据](REVIEW_POST_B_2026-09-17_EVIDENCE.json)。脚本及日志位于本机 `/tmp/cellmesh_review_post_b_20260917/`。

## B04：文本输入的 NUL 字符可能被静默截断（P2，优先处理）

**位置。** [io.py](../cellmesh/io.py) `_csv_columns()` 第 69–72 行、`_read_text_matrix()` 第 114 行、`_read_metadata_table()` 第 128 行及 `_read_name_list()` 第 146 行。

**问题与原因。** 这些路径在原始文本检查前就使用 pandas 默认 C 引擎。该引擎可将字段在真实 NUL 字符（字节 `0x00`）处截断；后续标识与数值校验只能看到已经被改写的字段。表达文件的表头“原始检查”也使用同一解析器，因而不能发现这种变化。这不是字面字符串 `"\\0"`，也不是正常缺失值 NA。

已复现：

| 文件原始内容（用 `\x00` 表示实际 NUL） | 默认读取后的值 |
|---|---|
| 基因名 `G\x00suffix` | `G` |
| 细胞 ID `a\x00suffix` | `a` |
| 表达字段 `1\x00junk` | 数值 `1` |
| 元数据分组 `A\x00suffix` | `A` |
| MTX 配套名称条目 `G\x00suffix` | `G` |

CSV/TSV 表达文件均复现；前三类输入还成功进入 `run_cell_mesh()` 并输出事件分数和 p 值。基因名截断可能导致错误匹配先验，细胞 ID/分组截断可能导致错误对齐或合组，数值截断可能把非法表达当作合法数值。只影响含 NUL 的异常文本；普通合法 CSV/TSV 不受此触发条件影响。Python 引擎在本环境下报错而不是静默截断。先验 CSV 与 10X 注释已有 `_delimited_records()` 的 NUL 检查。

**最小复现。**

```python
from pathlib import Path
from cellmesh import read_anndata

path = Path("nul.csv")
path.write_bytes(b"cell,G\x00suffix,R\na,1,2\nb,3,4\n")
adata = read_anndata(path, mode="csv")
assert adata.var_names.tolist() == ["G", "R"]  # 当前错误：未拒绝被截断的标识
```

补充范围确认：gzip CSV/TSV 同样截断；通过公开 `read_anndata(..., cell_meta_path=...)` 和 MTX 的 `genes_path` / `barcodes_path` 也复现。因此不是只存在于直接调用私有辅助函数的情形。

**方案。** 在首次 pandas 解析前，对解码后的原始文本执行 NUL 检查，报出文件和行位置。覆盖表达 CSV/TSV、元数据及 MTX 名称表；保持调用者选择的编码、压缩及读取选项。不能通过切换引擎、事后检查解析结果或删除 NUL 后继续来代替。验收应覆盖表头/正文、名称/数值、压缩文件，以及合法引号和多行字段，避免破坏既有读取能力。

## B05：可空数值读取选项返回 object 表达矩阵，正常数据随后无法评分（P2）

**位置。** [io.py](../cellmesh/io.py) `_read_text_matrix()` 第 114–118 行、`_read_csv()` 第 317 行；最终由 [preprocess.py](../cellmesh/preprocess.py) `_validate_expression_values()` 第 96–97 行拒绝。

**问题与原因。** `read_anndata(..., mode="csv", dtype_backend="numpy_nullable")`，或 `dtype="Float64"` / `dtype="Int64"`，会得到 pandas 的可空数值列。直接交给 `AnnData(df, ...)` 后，`adata.X.dtype` 在本环境下变成 `object`。即使所有值都是有限非负数、没有缺失值，主流程仍报“scoring expression must be finite and non-negative real values”。这里的大写 `Float64` 是 pandas 扩展类型，与正常工作的 `numpy.float64` 不同。

**影响范围。** 使用上述可选 dtype 参数的 CSV/TSV 读取链路，后续两种推断模式均无法评分；默认读取以及 `dtype=numpy.float32/float64` 不触发。没有证据表明它会产生错误 p 值，表现为读取成功后评分中断、报错原因不够准确。

最小输入为 `cell,G,R` 表头和 `a,1,2` / `b,3,4` 两行。默认读取与显式 NumPy float32 可评分；可空 dtype 读取后不能评分。仅在合成对照中把合法 object 数字显式转回 float64，两个可空类型均恢复与默认输入完全一致的事件得分和 p 值。

补充 gzip CSV/TSV × transpose 开/关 × 三类可空类型选项，共 **12 组**均得到 object 矩阵并在评分时被拒绝。

**方案。** 在读取边界将受支持的 pandas 可空数值列转换为 AnnData 支持的 NumPy 数值矩阵；缺失值保留为 NaN，非法文本、布尔、复数等不要通过无条件强制转换掩盖。保留细胞/基因轴和 transpose 行为。若某个后端无法安全支持，应在读取时明确拒绝该选项。验收覆盖无缺失/有缺失、混合数值列、CSV/TSV、transpose，以及读取后两种模式的实际评分。

## B06：float16 外部结果在网络图和 dotplot 中排序失败（P3，兼容性）

**位置。** [plotting.py](../cellmesh/plotting.py) `_plot_numeric_series()` 第 71–75 行、网络图第 655 行、dotplot 第 1281 行及后续排序。

**问题与原因。** 共用数值检查确认值合法后，对所有原生数值 dtype 直接保留原类型，包括 float16。pandas 在这些多列排序路径内部构建分类索引，而当前环境不支持 float16 索引，于是报 `NotImplementedError: float16 indexes are not supported`。

**影响范围。** 外部 DataFrame 或用户将结果列压缩为 float16 后绘图：网络图的 score、p 或 FDR 列可触发；dotplot 的 score 或 FDR 列可触发。已在取消显著性阈值的对照中排除“没有显著事件”的干扰。主流程原生输出为 float64，不受此触发条件影响；本次 float16 受体小提琴及计数图探针未复现同一错误。

**方案。** 对通过严格数值检查的 float16 绘图工作副本提升为 float32/64，再筛选、排序和聚合；不改写调用者的数据。提升仅解决运算兼容性，不能恢复用户此前压缩损失的精度。对照已验证：把现有 float16 值原样提升到 float64 后，两图均可完成绘制。验收比较排序、去重、事件选择及图形渲染，并保留 float32/64、nullable 类型和真正 NA 的既有行为。

## 建议次序与边界

修订按 B04、B05、B06 次序执行。B04 涉及异常文件的静默内容改变；B05、B06 属于特定选项或数值类型下的兼容性。此次结果不推翻 B01–B03 的已验证修订，也不把测试通过解释为所有可能输入均无问题。

## B04–B06 实施记录

三项均先以新增测试复现失败，再修改实现，并保留源码防回退注释。

- **B04：** 共用 `_reject_nul_text()` 在首次 pandas 解析前逐行检查解码文本，使用调用者的编码、错误处理和压缩配置；表达表、元数据及 MTX 名称表全部接入。新增 `test_text_nul_inputs.py` **66 项**，覆盖 CSV/TSV、C/Python 引擎、UTF-8/16/32、gzip、表头/正文/配套文件、字段内换行和读取筛选参数。结合既有读取/先验测试 **237 项通过**。
- **B05：** `_read_text_matrix()` 在转置前把合法 pandas 可空数值列转换为 NumPy 数值列；无缺失整数不经 float64，缺失值保留 NaN。新增 `test_nullable_expression_reader.py` **77 项**，覆盖两种推断模式、CSV/TSV、压缩/转置、混合 dtype、整数精度、真正缺失和非法输入；与既有标识读取测试 **138 项通过**。
- **B06：** `_plot_numeric_series()` 在严格校验后将 float16 工作副本提升为 float64，使排序/去重能正常执行；其余合法 dtype 保持原有契约。新增 `test_plotting_float16.py` **25 项**，对照相同 float16 值显式提升后的事件选择、去重、点大小和图形绘制，并检查调用者数据不变及非法值仍报错。结合既有绘图输入、QC 和显著性用例，**564 项通过**。

README 和读取说明已补充输入规范。本轮新增 **168 项**回归测试；全量结果为 **2493 passed、1 skipped、222 warnings，262.92 秒**，包括四本 notebook 的执行和图形测试。原有 213 条 AnnData 警告保留，新增 9 条来自新增混合整数/浮点或缺失整数测试中 AnnData 将矩阵统一为 float64 的提示；这些用例均验证数值和缺失语义正确。全量日志为本机 `/tmp/b04-b06-full.log`，`git diff --check` 通过。修订前证据 JSON 保持不变。

复验命令（本机 Loom 可选依赖路径）：

```bash
PYTHONPATH=/tmp/cellmesh-a10-deps MPLCONFIGDIR=/tmp/cellmesh-mpl NUMBA_CACHE_DIR=/tmp/cellmesh-numba /home/qsong/miniconda3/envs/cellmesh/bin/python -m pytest -q
```

B04–B06 在上述范围内完成修订，修改保留在本地工作树。
