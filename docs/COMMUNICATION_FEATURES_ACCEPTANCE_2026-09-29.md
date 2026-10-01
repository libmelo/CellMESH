# 通讯特征、概览卡片及全量 CSV 统一验收

> 当前展示已改为静态矩阵点图，四层图停用。本文为历史记录；见 [最新修订与验收](DOTMAP_REVISION_2026-09-30.md)。

状态：本轮执行与核对完成，自动检查通过，待用户人工审查。真实数据 Notebook 于 2026-09-29 完成；2026-09-30 复核全部代码单元已执行、无错误输出、交付文件完整、当前包代码与验收快照一致，并检查保存的报告截图。

本记录保留该次执行的历史结果。随后用户要求参照 MEBOCOST In[13] 修订四层图，新实现与独立验证见 [2026-09-30 FlowPlot 修订记录](FLOWPLOT_REVISION_2026-09-30.md)；本页代码一致性及测试数字不代表后续修订的验证结果。

## 数据与运行入口

本次复用原完整计算结果，重新生成交付，没有重新计算评分或置换，也没有抽样或修改表达、统计量、QC、标签或 selection。

- 完整数据：GSE176078，100,064 个细胞、29,733 个基因、26 位患者、8 类细胞。
- 原计算结果：`cell_mesh_validation/positive_delivery_2026-09-28/`。
- 原参数：1,000 次置换、随机种子 0；其他参数逐项与原 H5 核对。
- 原阳性标准：`fdr_sensor_type <= 0.05`，不附加 QC 或最低评分门槛。
- 原 selection：`hmdb_ids=['HMDB0000961'], sensor_genes=['CYSLTR1']`；pooled 28 条，sample_aware 26 条。
- 三级报告沿用前次 selection 内、外各一条完整事件。
- 主交付沿用 `include_input=True` 附带完整输入，未指定 `include_all_events`，应默认关闭；pooled 另以同一完整结果验证显式开启。

[已完整执行并保存输出的审查 Notebook](../../cell_mesh_validation/communication_features_2026-09-29/communication_features_review.ipynb) 位于包外。报告由包的 `save()`、`to_html()`、`add_selection()`、`add_event_reports()` 生成；辅助 Python/JavaScript 只核对数据、操作浏览器和截图。

## 最终检查结果

| 项目 | 本轮结果 |
|---|---|
| 完整 pytest 回归 | 2,725 通过，7 跳过 |
| 补充概览卡片回归 | 6 通过 |
| 补充 WebKit 边界回归 | 1 通过，覆盖 105 条路径、跨页定位、多对象切换及 1,400/420 像素窗口 |
| 原输入、先验与原计算文件哈希 | 与原留存哈希一致 |
| pooled 默认交付与追加 | 通过；全量 CSV 默认关闭，已有表达图和所选 CSV 在追加后保留 |
| pooled 显式全量 CSV | 与原全量 CSV 逐字节一致 |
| pooled 离线浏览器 | 一级、二级、三级断言通过 |
| sample_aware 交付 | 通过；原 selection 26 条，两端表达图各一张，保留 26 个样本的评分记录 |
| 最终科学 H5 全节点比对 | 通过；pooled 180,316 个节点、sample_aware 4,664,968 个节点 |
| 真实数据 Notebook | 8 个代码单元完整执行，0 个错误输出，212 项数据/交付检查通过 |
| 两模式离线浏览器 | 271 项断言通过；三级报告导航、路径/列表联动、表达图、展开无内部滚动条 |
| 当前代码一致性 | 所有包 Python 文件与生成报告时的 SHA-256 快照一致 |
| 人工图形审查 | 待用户审查 |

自动化边界测试中的合成输入只用于测试程序边界，未用于真实数据报告，也不替代真实数据验收。真实案例少于 50 条且仅一个代谢物—受体组合，因此不通过改变真实 selection 来制造翻页或多对象效果。

## 本轮修正

- 将旧测试中“默认全量 CSV”的断言改为显式开启；默认行为单独验证为关闭。
- 同步二级结构、图形折叠状态、条件计数热图、概览卡片和内嵌脚本断言。
- 修复二级报告末节标题与导航名称不一致，统一为“筛选与导出”。
- 增加路径完整性、评分缺失、无隐藏显著性/QC 筛选、表达图复用与替换、CSV 开关继承及旧交付迁移的回归检查。

## 审查入口

- [pooled 二级报告](../../cell_mesh_validation/communication_features_2026-09-29/pooled/analysis.selected.report.html)
- [sample_aware 二级报告](../../cell_mesh_validation/communication_features_2026-09-29/sample_aware/analysis.selected.report.html)
- [pooled 报告导航](../../cell_mesh_validation/communication_features_2026-09-29/pooled/analysis.reports/index.html)；[sample_aware 报告导航](../../cell_mesh_validation/communication_features_2026-09-29/sample_aware/analysis.reports/index.html)。
- [机器可读结果汇总](../../cell_mesh_validation/communication_features_2026-09-29/summary.json)：`status=PASS`，`human_review=pending`。
- 执行过程：`execution.log`；自动回归日志：`checks/cellmesh-acceptance-*.log`。
- Notebook 最后统一生成并展示真实报告截图；历史验收目录未改写。

## 科学结果与文件范围

| 模式 | 全部事件 | 阳性事件 | selection 事件 |
|---|---:|---:|---:|
| pooled | 246,120 | 11,607 | 28 |
| sample_aware | 246,120 | 6,053 | 26 |

两种模式各输出 9 个文件：完整结果 H5、完整输入 H5AD、阳性 CSV、所选 CSV、一级 HTML、二级 HTML、导航 HTML、2 个独立单事件 HTML。H5AD 是本次沿用 `include_input=True` 的选择，不是修改包默认值。主交付目录均无 `.all_events.csv`；显式开启的独立检查文件位于 `checks/explicit_all/`。

H5 核对仅排除交付配置与包版本元数据；其余参数、科学表格、先验、输入归档元数据和解释中间量均逐节点一致。原输入、先验、原 H5、原输入归档及原全量/阳性 CSV 哈希在执行前后不变。阳性 CSV 和输入 H5AD 与原文件逐字节一致。

## 截图和实际限制

- [pooled 四层图](../../cell_mesh_validation/communication_features_2026-09-29/checks/pooled.paths.png)、[sample_aware 四层图](../../cell_mesh_validation/communication_features_2026-09-29/checks/sample_aware.paths.png)。
- [pooled 表达图](../../cell_mesh_validation/communication_features_2026-09-29/checks/pooled.expression.png)、[sample_aware 表达图](../../cell_mesh_validation/communication_features_2026-09-29/checks/sample_aware.expression.png)。
- 截图保存当前浏览器视口，页面下方部分需在 HTML 中继续向下阅读；HTML 图形没有内部滚动条。Notebook 已内嵌这些截图。
- 主真实 selection 为 28/26 条、一个代谢物—受体组合。分页、多对象和混合统计状态由纯自动化边界测试覆盖，并未更换真实数据案例。
- 本轮完整 Notebook 约耗时 70 分钟，主要用于完整 sample_aware H5 的读取、解释表构建、保存与全节点比对。未执行评分或置换计算；持久化耗时未作为本轮功能修改范围。
- 浏览器日志包含旧 WebKit API 的弃用提示及无 GPU 环境提示；本轮使用软件渲染，所有浏览器断言与截图均完成。7 项 pytest 跳过不计为通过。
