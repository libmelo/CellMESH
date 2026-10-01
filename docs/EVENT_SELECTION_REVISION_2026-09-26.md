# 事件筛选与交付联动修订（2026-09-26）

本地版本：0.7.0。实现统一 `select_events()`、筛选 CSV、报告与聚合绘图联动。
本次不修改主线评分、标准化、置换或 FDR 计算。

## 前后对比

| 项目 | 原交付 | 本次修订 |
| --- | --- | --- |
| 选择事件 | 手工填写四个身份字段，或自行操作 DataFrame | `res.select_events()` 统一按细胞、定向细胞对、HMDB、受体、FDR、评分与 QC 筛选 |
| 条件组合 | 由调用者自行实现 | 条件之间取交集，列表内部取并集；`cell_pairs` 进一步约束端点列表 |
| 筛选说明 | 无统一记录 | `focus.criteria`、`focus.summary` 和 HTML 保存条件、逐步剩余数量与缺失排除数量 |
| CSV | `.all_events.csv` 全量表 | 全量表保留；`focus.to_csv(path)` 单独输出全部筛选事件 |
| HTML 范围 | 全量列表及手工指定事件解释 | 支持 `selection=focus`，提供所选概览、列表与解释 |
| 解释数量 | 每个手工指定事件均解释 | 筛选报告默认解释前 10 条；`detail_n=0` 不解释，`None` 全解释；旧手工指定行为保留 |
| 深层解释表 | 表格均默认展开 | 基因、反应、参考值、覆盖和先验表默认折叠 |
| 聚合绘图 | 各自存在默认 FDR/QC/top-N 规则 | 传入 `focus` 后使用统一范围，不再隐含筛选；直接传入原结果仍保持旧默认 |
| 点图显示限制 | 默认前 10 个代谢物—受体组合 | 对 `focus` 默认全显示；显式 `top_n` 仅作显示限制，标注选择/绘制数量 |
| 空结果 | 缺少统一操作路径 | 筛选 CSV 有表头，HTML 显示零事件，绘图明确报错，不放宽阈值 |
| 全量结果 H5 | 完整结果与解释中间量 | 继续保持全量；保存筛选报告不会缩减结果文件 |
| 表达矩阵 H5AD | 可选伴随文件 | 沿用原有 `retain_input=True` / `include_input=True` |

## 用法

```python
focus = res.select_events(
    sender_labels=["A", "C"], receiver_labels=["B"],
    cell_pairs=[("A", "B"), ("C", "B")],
    max_fdr=0.05, min_cell_mesh_score=0.1, qc_only=True,
)
focus.to_csv("analysis.selected_events.csv")
files = res.save("analysis", report=True, selection=focus, detail_n=10)

from cellmesh import plot_event_dotplot, plot_communication_network
plot_event_dotplot(focus)
plot_communication_network(focus)
```

也可单独 `focus.to_html("focus.html", detail_n=None)`。
原默认 `save(report=False, include_input=False)` 保持不变；筛选 CSV 按需输出。
默认 `select_events()` 不设评分/FDR/QC 门槛；示例阈值由用户根据分析目的指定。

筛选保留原值，不对子集重新计算 FDR。未设置对应阈值时保留缺失统计量；
设置阈值时排除并计数。FDR 全缺失却要求 FDR 阈值时报错。
筛选对象绑定原结果；改变事件表或科学参数后必须重新筛选。
完整规则见 [RESULT_FILES_AND_REPORTS.md](RESULT_FILES_AND_REPORTS.md#reusable-event-selection-07)。

## 示例与实现位置

- [示例交付说明](../examples/selection_demo_output/README.md)：全量 9 条、筛选 2 条、报告解释 1 条。
- [筛选 CSV](../examples/selection_demo_output/analysis.selected_events.csv)
- [筛选报告](../examples/selection_demo_output/analysis.report.html)
- [筛选点图](../examples/selection_demo_output/selected_events.png)
- [可运行脚本](../examples/export_explain_report.py)：加 `--selection-demo`。
- [筛选实现](../cellmesh/selection.py)：条件、身份验证、排序、元数据、CSV 导出。
- [结果对象](../cellmesh/core.py)：公开 `select_events()` 入口和保存/报告参数。
- [报告](../cellmesh/reporting.py)：所选范围、通讯概览和解释数量。
- [持久化](../cellmesh/result_io.py)：全量交付与可选筛选报告。
- [绘图](../cellmesh/plotting.py)：选择对象适配和显示范围标注。
- [新增回归测试](../tests/test_event_selection.py)：筛选与输出一致性、错误边界、旧接口兼容。

## 验证

- 全套测试：2627 项通过、7 项跳过、230 条警告，运行约 5 分钟；无失败。
- 全套测试启动后追加的预览边界与 QC 参数检查，以及报告最终调整：筛选/报告定向复验 109 项通过。
- 网络图标题最终调整：相关筛选图与旧默认兼容性复验 4 项通过。
- 筛选与绘图相关测试：614 项通过。
- 原报告、输入归档和导出测试：122 项通过。
- 示例已验证 HDF5/全量 CSV/筛选 CSV 数值一致，伴随输入重新运行后事件完全一致。
- 示例点图已人工检查，选择与绘制数量均为 2，标签与图例完整。
- `git diff --check` 通过；修订说明与示例说明中的本地文件链接已核对。

这些改动保留在本地工作区，尚未提交或推送到 GitHub。
