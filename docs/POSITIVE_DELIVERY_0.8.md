# 阳性结果与感兴趣事件：0.8 交付说明

> 历史说明：本文记录 0.8 的交付行为。当前版本采用独立的三级报告，selection 不再自动生成事件解释；请以[三级报告使用说明](THREE_LEVEL_REPORTS.md)为准。

> 全量 CSV 默认值已修正：当前仅在显式设置 `include_all_events=True` 时导出。下文“每次 save 均输出全量 CSV”为旧行为，不能作为当前调用依据；修正尚待统一测试，见[修订记录](ALL_EVENTS_CSV_REVISION_2026-09-29.md)。


`run_cell_mesh()` 计算科学结果；`save()` 按阳性标准导出通用结果；
`selection` 单独指定感兴趣事件。保存、筛选、追加均不重新计算评分、置换或 BH 校正。
下面的示例是独立流程；首次交付示例从尚未保存交付设置的 `result` 开始。
输出目录需要预先创建。

## 默认行为

| 项目 | 行为 |
|---|---|
| 阳性标准 | `fdr_sensor_type <= 0.05`，按受体类型进行 BH 校正 |
| 可选统计量 | `fdr_global` 或 `perm_pvalue`；一次使用一种 |
| QC、评分门槛 | 默认均不附加；可显式指定 `positive_qc_only`、`positive_min_score` |
| HTML 报告 | `report=True`，默认输出；可显式设置 `report=False` 关闭 |
| 输入归档 | `include_input=False`；开启时要求计算阶段 `retain_input=True` |
| 缺失统计量 | 无法判定，保留全量结果；空阳性表不等于不存在通讯 |
| `select_events()` | 默认 `source="positive"`；`source="all"` 可包括未达标准或无法判定事件 |
| 解释数量 | `detail_n=10`；0 不生成面板，`None` 解释全部感兴趣事件 |
| 点图数量 | 兴趣事件点图 `plot_top_n=20` 个代谢物—受体组合；`None` 不限 |

FDR 字段保留兼容名称，含义是 BH 校正后的 p 值。使用原始 p 值时，报告明确标注未校正。
“阳性”只表示满足本次标准，不表示已得到实验验证。全量数量以 `result.events` 为准，不推算理论候选组合。

## 一次导出

```python
focus = result.select_events(sender_labels=["Myeloid"])
files = result.save("analysis", selection=focus)
```

生成：

```text
analysis.cellmesh.h5
analysis.all_events.csv
analysis.positive_events.csv
analysis.selected_events.csv       # 指定 selection 时
analysis.report.html               # 默认通用报告
analysis.selected.report.html      # 指定 selection 时的独立报告；report=False 关闭两份 HTML
analysis.input.h5ad                # include_input=True
```

首次调用 `result.save("analysis")` 默认交付 H5、全量 CSV、阳性 CSV 和 HTML 四个文件。
返回字典默认包含 `result`、`events_csv`、`positive_csv`、`report`，按需增加 `selected_csv`、`selected_report`、`input`。
两份通用 CSV 不依赖 selection，也不依赖是否生成 HTML。感兴趣 CSV 自动导出，不必另外调用 `to_csv()`。
全量 CSV 在每次 `save()` 成功时写出，包括未达阳性标准、无法判定和未被 selection 选中的已生成事件。
仅运行 `run_cell_mesh()` 不会自动写文件；后续 `add_selection()` 保留已有通用 CSV，旧交付缺失时补齐。
该 CSV 增加 `positive_status`：`positive`、`not_positive`、`not_evaluated`；原始结果表和评分值不改写。

## 调整阳性标准

```python
# 通用结果按原始置换 p 值 <= 0.01 判定
result.save("pvalue_analysis", report=True,
            positive_col="perm_pvalue", positive_threshold=0.01)
focus = result.select_events(sender_labels=["Myeloid"])  # 继承已保存标准
result.add_selection(focus)

# 首次导出时也可以让 selection 显式携带该标准
focus = result.select_events(positive_col="fdr_global", positive_threshold=0.05,
                             sender_labels=["Myeloid"])
# 在尚无交付设置的新结果上，save() 省略阳性参数时采用 selection 的标准。
```

省略的阳性参数继承已保存设置；尚无设置时使用默认标准。
`positive_qc_only=False` 可关闭已启用的 QC 门槛，显式 `positive_min_score=None` 可清除评分门槛。
若阳性来源的 selection 与 `save()` 的标准不一致，会要求重新筛选，不隐式改变事件范围。
`source="all"` 不受阳性门槛限制；报告和交付 CSV 根据本次交付标准标注它的阳性状态。

## 查阅后追加

```python
result.save("analysis")  # 默认生成 HTML；无事件解释部分

from cellmesh import CellMeshResult
loaded = CellMeshResult.load("analysis.cellmesh.h5")
focus = loaded.select_events(sender_labels=["Myeloid"])
loaded.add_selection(focus)

# 同一个交付包目前维护一组感兴趣事件，替换需显式声明。
other = loaded.select_events(source="all", sensor_genes=["SLC8A1"])
loaded.add_selection(other, replace=True, detail_n=10)

restored = CellMeshResult.load("analysis.cellmesh.h5")
saved_focus = restored.get_selection()
```

追加继承阳性标准和省略的展示参数；显式 `detail_n=None` 或 `plot_top_n=None` 仍表示不限。
通用 HTML、全量 CSV、阳性 CSV、表达矩阵不重写；更新所选 CSV、所选 HTML 和 H5 交付元信息。
旧版混合 HTML 在首次追加时迁移为两份报告，此后保留通用报告。`report=False` 的交付追加时仍不生成 HTML。
更新采用临时文件与失败回滚；保存包已被其他对象更新、通用 CSV 被修改或源结果变化时拒绝更新。
已有解释对象要求 `replace=True`。同前缀 `save()` 省略 selection 会保留已保存的感兴趣事件；
若要替换，请使用追加接口。旧 H5 没有交付设置时采用默认阳性标准，首次追加补齐缺失的通用 CSV。

## HTML 结构与独立导出

通用报告固定四部分：分析概览、通讯概览、阳性事件列表、数据与方法。
概览记录所有已生成事件数量及阳性状态；图形和列表展示阳性事件。
通讯概览仅含计数热图与网络图，不含点图或逐事件解释。

所选报告固定四部分：所选事件概览、所选事件列表、事件解释与表达证据、筛选与方法。
计数热图默认展开，网络与代谢物—受体点图折叠；图形展开后完整显示，不用内部滑块。
列表包含全部所选事件，不再追加阳性过滤；解释入口整合进列表，无重复解释索引。
详细解释默认前 10 个事件，发送端的两个面板归为一组、接收端单列。
sample_aware 的跨样本评分图展示全部样本，标示 NA、零分、QC 与中位数；pooled 不画该面板。
组成图仍最多展示前 10 个样本单元，明细表预览最多 100 行；完整数据在 H5/API 中。
全局覆盖度与完整参数仅在通用报告中展示，所选报告保留必要统计口径和筛选过程。

小提琴图使用显式 `adata` 或计算时保留的输入；相同代谢物/受体上下文跨事件复用。
只为已有详细解释的事件生成，默认折叠，缺少输入时简短说明。
发送端图是生产相关表达分布，颜色为代谢物可用性（色标 0–1）；接收端图是受体表达分布，颜色为均值。
全零接收端使用 0–1 色标；全零或常数表达分布显示水平线。
这不是实测代谢物浓度或通量。sample_aware 小提琴图按细胞类型合并所提供输入的细胞，非逐样本图。
只有 H5 或输入归档路径时不会自动读取大型表达矩阵；需要时显式提供：

```python
loaded.add_selection(focus, adata=loaded.load_input())  # 首次追加
# 已有 selection 时须另加 replace=True
focus.to_html("my_review.report.html", adata=loaded.load_input())
```

```python
result.to_html("general_review.report.html")  # 通用 HTML + 全量/阳性 CSV
focus.to_html("my_review.report.html")        # 所选 HTML + 所选 CSV
focus.to_csv("my_interests.csv")              # 仅所选 CSV
```

独立导出不生成 H5、不修改原交付设置，需使用独立前缀。
标准交付文件受到覆盖保护，原包的追加请使用 `add_selection()`。
HTML 和 CSV 需要一起移动，链接使用相对路径。有已保存通用报告时，所选报告提供返回链接。
旧 `events=[...]` 仍可解释指定样本或未生成候选；`save()` 将其解释置于独立所选 HTML。
该兼容接口不凭空生成候选事件 CSV；新流程建议使用 `select_events()`。

本轮具体改动与验证范围见 [报告拆分记录](REPORT_SPLIT_2026-09-28.md)。
