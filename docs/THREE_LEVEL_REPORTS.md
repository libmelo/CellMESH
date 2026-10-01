# 报告与交付接口

一级默认生成；二级、三级使用不同入口，互不依赖。所有筛选沿用原始事件评分、p 值和 BH 校正值，不重新计算算法或统计量。

| 范围 | 定义入口 | `save()` 参数 | 报告 |
|---|---|---|---|
| 全局 | 无需额外指定 | 无 | 一级统计筛选报告，保持原有内容 |
| 条件筛选 | `filter_events(...)` | `filtered_events=...` | 一份通讯筛选报告 |
| 明确配对 | `select_events(metabolite_sensor_pairs=...)` | `selected_events=...` | 每个明确指定的代谢物–sensor 对一份通讯报告 |

## 二级：条件筛选

```python
filtered = result.filter_events(
    hmdb_ids=["HMDB0000961", "HMDB0006465"],
    sensor_genes=["CYSLTR1", "SLCO2A1"],
    # sender_labels=[...], receiver_labels=[...], cell_pairs=[(sender, receiver)],
)
```

两个分子列表独立匹配，允许多对多，只纳入结果中实际存在的事件，不按列表位置配对。不同条件取交集，同一列表内取并集。`None` 不限制，空列表匹配零条。即使两个列表各只有一项，或结果只剩一个事件，仍属于二级。

结构：筛选概览 → 通讯特征 → 完整事件列表 → 筛选与导出。列表链接已经生成且包含对应事件的三级报告；未请求的三级不自动生成。

## 三级：明确配对

```python
selected = result.select_events(
    metabolite_sensor_pairs=[
        ("HMDB0000961", "CYSLTR1"),
        ("HMDB0006465", "SLCO2A1"),
    ],
    # sender_labels=[...], receiver_labels=[...], cell_pairs=[(sender, receiver)],
)
```

配对列表为必填条件；不接受独立的 `hmdb_ids` 或 `sensor_genes` 列表。附加细胞、统计、评分与 QC 条件取交集。一个分子对的多个细胞通讯事件集中在同一页；指定多个分子对则分别生成，不受绘图上限影响。

已知标识组成但未生成事件的分子对，以及事件被条件全部排除的分子对，均保留诊断页面。未知标识、格式错误及重复配对在写文件前报错。空配对列表生成空的指定事件 CSV，不生成三级页面。

结构：

1. 分子对概览：身份、实际条件、事件数量、涉及细胞及统计状态。
2. 通讯分布：静态矩阵点图、通讯网络和完整事件表。
3. 评分与表达依据：按发送端/接收端组织并复用共同证据。
4. 样本表现：sample_aware 样本×事件评分热图；pooled 省略。
5. 单事件详情：按细胞对折叠查看全部样本的评分分解、计算输入和诊断。
6. 数据来源与方法：共同反应、基因、先验、覆盖和标准化参考表，每类输入只内嵌一次。

不存在额外的“分子对总体评分”或“分子对总体 p 值”。完整事件列表、详情和 CSV 不按绘图上限截断。共同明细可按 sender、receiver 或 cell_type 定位；覆盖表示输入可用情况，不是事件置信度。

## 默认统计范围

两个入口默认 `source="positive"`，采用保存时的阳性标准；未保存时默认 `fdr_sensor_type <= 0.05`，不隐式增加 QC 或评分门槛。`source="all"` 可纳入未达标或无法判定的已生成事件。`max_fdr`、`max_perm_pvalue`、`min_cell_mesh_score` 和 `qc_only` 为附加条件。

阳性标准通过保存接口的 `positive_col`、`positive_threshold`、`positive_qc_only`、`positive_min_score` 指定。使用原始 p 值的例子：`positive_col="perm_pvalue", positive_threshold=0.01`。阳性来源的对象必须与交付采用相同标准，否则需重新筛选；报告记录实际标准及原始条件。

## 保存、加载和追加

```python
files = result.save(
    "analysis",
    filtered_events=filtered,
    selected_events=selected,
    adata=adata,
)
# 只传 filtered_events 或 selected_events 即独立输出该范围。
# 都不传，仅输出默认全局交付。
```

也可以先输出通用结果，再追加：

```python
files = result.save("analysis")
files = result.add_filtered_report(filtered, adata=adata)
files = result.add_event_reports(selected, adata=adata)
# 两个追加入口可交换顺序，不隐式继承彼此的条件。
```

同一交付保留一个二级集合，三级可持续追加不同分子对。重复的二级或分子对需显式 `replace=True`。替换一个分子对不改变其他分子对页面；一级、阳性 CSV、科学结果保持。普通 I/O 失败回滚；同一前缀不支持并发写入。

`get_filtered_events()` 恢复二级对象；`get_selected_events()` 返回按分子对组织的对象列表，保留每对独立的条件。新前缀的 `save()` 不自动恢复旧报告请求。

```python
filtered.to_csv("filtered.csv")
selected.to_csv("selected.csv")
filtered.to_html("filtered.html", adata=adata)
selected.to_html("molecular_pairs.html", adata=adata)
```

明确配对的独立 HTML 导出返回导航、各配对页面及合并 CSV 的路径字典。完整验收采用 `save()`/追加接口的统一交付，不用临时独立页面代替。

## 统一交付结构

```text
analysis.cellmesh.h5
analysis.positive_events.csv
analysis.report.html
analysis.filtered_events.csv              # 请求二级
analysis.filtered.report.html             # 请求二级
analysis.metabolite_sensor_events.csv      # 请求三级，所有请求配对的所选事件合并表
analysis.all_events.csv                    # 仅 include_all_events=True
analysis.input.h5ad                        # 仅 include_input=True
analysis.reports/
    index.html
    metabolite_sensor/
        pair_<稳定编号>.html               # 每个明确配对一份
    history/
        legacy-selection-<摘要>.zip        # 仅迁移存在旧名文件的交付时生成
```

返回路径键为 `result`、`positive_csv`、`report`、`navigation`，以及适用的 `filtered_csv`、`filtered_report`、`selected_csv`、`pair_report:pair_<编号>`、`events_csv`、`input`。配对编号仅由规范化的 HMDB ID 与 sensor 产生；细胞筛选改变时需替换同一分子对页面。

三级合并 CSV 包含 `pair_id`、完整四元事件标识、原始科学字段和 `positive_status`。空请求仍记录在 H5 交付元信息和诊断页中，不伪造 CSV 事件行。独立对象的 `to_csv()` 保留完整事件标识，便于自行分组。

`report=True` 默认；`report=False` 关闭全部 HTML，仍输出阳性 CSV 和已请求的两类筛选 CSV。全事件 CSV 默认关闭，每次 `save()` 都需显式开启；追加继承原开关。H5 始终保存全量事件。输入 H5AD 默认关闭；表达图使用显式 `adata` 或计算时保留的输入，加载结果不会隐式打开输入归档。

## 图形和表格

通讯矩阵的独立绘图、二级报告和三级报告均调用 `plot_communication_dotmap()`。`cellmesh` 与 `cellmesh.plotting` 导出同一个函数实现，原来的 `draw_communication_dotmap()` 和 `draw_dotmap()` 私有入口已合并移除。范围整理、颜色编码、图例和布局仍由各自的辅助函数负责。

公开函数返回 Figure 及原有字段，并增加 `scope`，记录实际展示的事件和分子组合。报告直接嵌入该 Figure，使用返回的坐标、编码、标签和展示范围组织说明及小提琴图；HTML 层不重新绘图。SVG 保留 Figure 的字体轮廓，避免浏览器字体替换和页面 CSS 改变排版。图形文字和事件身份同时保留在内嵌数据及事件表中。

- 通讯矩阵为单张静态图：横轴 receiver–sensor，纵轴 sender–metabolite；代谢物默认 HMDB ID，可在独立绘图中用 `metabolite_label="name"` 切换。横轴单行标签旋转 45°，纵轴保持水平。
- 点面积表示原评分，颜色表示所用 p/BH 的 −log10。多值使用连续色条，常值使用矩形色块；零概率、统计缺失和评分缺失单独标示。校正范围的文字仅限定 BH 校正。
- 保留细胞通讯网络；计数热图仅在多分子组合且细胞对计数存在差异时展示，三级不重复绘制该退化图形。
- 表达小提琴图按发送端代谢物、接收端代谢物–sensor 组织，共同输入只绘一次；它们描述表达证据，不代表实测浓度或通量。
- 二级 `plot_top_n=20` 限制图形展示的分子组合及样本热图的事件列；三级始终绘制当前配对全部符合条件的事件，样本热图仍可限制事件列。所有样本行、事件列表、详细数据和 CSV 不截断。
- 图片完整展开，不使用内部图片滑块。表格完整内嵌，支持离线搜索、排序、分页；未启用 JavaScript 时仅静态预览前 50 条。

## 从旧接口迁移

| 旧用法 | 新用法 |
|---|---|
| `select_events(hmdb_ids=..., sensor_genes=...)` | `filter_events(hmdb_ids=..., sensor_genes=...)` |
| `save(selection=focus)` | `save(filtered_events=filtered)` |
| `add_selection(focus)` | `add_filtered_report(filtered)` |
| `save(events=[完整四元组])` | `select_events(metabolite_sensor_pairs=..., cell_pairs=...)` 后传入 `selected_events` |
| `add_event_reports([完整四元组])` | `add_event_reports(selected)` |

旧四元组请求若涉及多个分子对且各自细胞范围不同，应按分子对分别创建对象并追加，不能将所有 `cell_pairs` 混入同一个请求。

旧式 `select_events()` 和 `events=` 明确报出迁移提示，不静默扩大范围。旧 `selection=`、`add_selection()` 暂作二级兼容别名并发出弃用提示。旧 H5 可读取；原四元组仍可用 `explain_event()` 查询，历史单事件页面仅在导航的历史区沿用。迁移旧交付时二级改用 `filtered` 文件名，新三级放入 `metabolite_sensor/`；一级保持原位置和内容。

保存或追加报告时，已被替代的 `<prefix>.selected.report.html` 和 `<prefix>.selected_events.csv` 自动归入 `<prefix>.reports/history/legacy-selection-<摘要>.zip`，当前目录只保留新命名文件。归档与结果更新在同一事务中完成，写入失败恢复原文件；相同内容复用已核验的归档。ZIP 保留原文件字节和 SHA-256 清单，属于历史快照，旧 HTML 的相对链接仍基于原目录，不作为当前报告入口。

对于已迁移为当前布局、但残留旧名文件的交付，可单独整理，无需读取完整科学数据或重绘报告：

```python
from cellmesh.result_io import archive_legacy_reports
archives = archive_legacy_reports("output/analysis")
```

重复执行不会新增归档；新报告或 CSV 缺失时会报错并保留旧文件。旧名文件对应的 Jupyter checkpoint 一并归档，其它 checkpoint 保留。旧布局交付应先用保存或追加接口迁移。该维护调用不改写 H5、当前 HTML、当前 CSV 或输入 H5AD。

## 本轮验收（2026-10-01）

人工审查入口为包外 `../cell_mesh_validation/report_scopes_2026-10-01/report_scopes_review.ipynb`，已保存实际执行输出。各模式的完整交付仍从 `analysis.reports/index.html` 进入；源交付 `tier3_text_review_2026-10-01` 保持原样。

沿用 GSE176078 的完整原始计算结果和输入归档（100,064 个细胞、29,733 个基因、26 个样本），未重算算法。`n_perms=1000`、`random_state=0`、原有细胞/QC 参数和 `fdr_sensor_type <= 0.05` 阳性标准保持不变。两种模式各含 246,120 条汇总事件；一级报告、阳性 CSV、输入 H5AD、科学数据及原始源文件的完整性核对通过。

| 模式 | 二级事件数 | 三级明确配对及事件数 | 三级样本明细 |
|---|---:|---|---:|
| pooled | 28 | HMDB0000961 / CYSLTR1：28；HMDB0006465 / SLCO2A1：6 | 34 条合并分析记录 |
| sample_aware | 26 | HMDB0000961 / CYSLTR1：26；HMDB0001999 / SLC8A1：7 | 858 条样本记录 |

二级沿用原独立列表条件，三级用新接口明确指定分子对。每个配对的完整事件、原评分/p/BH、样本明细和 CSV 已对照结果对象核验，未按展示上限截断。默认未生成全事件 CSV。两种模式共 8 个当前页面通过浏览器检查，截图保存在验收目录的 `checks/`。

自动回归为 260 项通过、1 项可选浏览器检查跳过；该可选检查随后单独启用并通过，静态图浏览器检查也通过。最后的导航标注修订另有 32 项相关回归通过（与前述测试重叠，不累加计数）。桌面、窄屏和禁用 JavaScript 的静态图检查均通过；本轮同时修复了统计图例右缘的裁切问题。

## 绘图入口统一验收（2026-10-01）

后续统一绘图调用的验收位于 [dotmap_unified_review.ipynb](../../cell_mesh_validation/dotmap_unified_2026-10-01/dotmap_unified_review.ipynb)。沿用上一节完整交付的原数据、参数和筛选条件，一级报告沿用，二级与全部三级报告通过包的追加接口重新生成。源文件、H5 科学数据、计算参数和 CSV 完整性核对通过。

- 独立图修订前后逐像素一致；此前 52 条事件、4 个 BH 值的图形范围完整保留。
- [52 事件对照报告](../../cell_mesh_validation/dotmap_unified_2026-10-01/comparison_52.report.html)与同一 H5 事件表的独立图逐项比较：事件、坐标、点面积、颜色、标签和图例一致。SVG 使用 Figure 的字体轮廓，实际浏览器截图已检查。
- 两种模式的二级及每对三级矩阵均核实来自公开 `plot_communication_dotmap()` 的返回结果；26 + 7 条 sample_aware 事件的全部 858 条样本记录保留。
- 266 项回归通过；默认跳过的两项浏览器检查另行启用并通过。真实数据的 8 个交付页面及 52 事件对照报告共 9 页通过浏览器检查。
- 初次对照将 CSV 解码值与 H5 原值直接比较，发现约 1e-16 的浮点读取差异。失败记录保存在 `checks/history/`；随后改用同一 H5 输入严格比较，未修改科学数据或放宽相等断言。

统一交付入口：[pooled](../../cell_mesh_validation/dotmap_unified_2026-10-01/pooled/analysis.reports/index.html)、[sample_aware](../../cell_mesh_validation/dotmap_unified_2026-10-01/sample_aware/analysis.reports/index.html)。执行日志、截图、公开入口调用记录与代码指纹均保存在该包外验收目录中。

## 旧文件归档核验（2026-10-01）

[delivery_cleanup_review.ipynb](../../cell_mesh_validation/dotmap_unified_2026-10-01/delivery_cleanup_review.ipynb) 在当前统一交付目录运行包内维护函数。pooled 的旧 HTML、CSV 及其 HTML checkpoint 共 3 个文件，以及 sample_aware 的旧 HTML、CSV 共 2 个文件，已原样归入各自 `analysis.reports/history/` 下的 ZIP。归档内容和原文件 SHA-256 一致，当前目录不再保留这些旧名文件。

本轮只整理文件；计算结果、各级报告和图形均沿用。两个模式其余 26 个文件的完整 SHA-256 不变，包括 H5、输入 H5AD、当前 CSV、报告及导航。14 个现存报告页面的 221 个本地链接与锚点全部通过核验，重复运行无变更。原执行 Notebook、原验收记录和历史源目录保留。

保存/追加/持久化相关自动回归 80 项通过；加入 checkpoint 归档后，最终归档专项 12 项通过（与前述测试重叠，不累加），覆盖旧版本迁移、既有残留整理、重复执行、归档冲突、符号链接保护及失败回滚。核验摘要和本轮代码指纹保存在 `checks/delivery_cleanup/summary.json`。
