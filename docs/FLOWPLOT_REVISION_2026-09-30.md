# 四层通讯关系图：MEBOCOST FlowPlot 参考修订

> 后续修订：用户决定停用四层图，改为参考 MEBOCOST In[12] 的静态矩阵点图。本文保留为历史记录；当前实现和验收见 [矩阵点图修订](DOTMAP_REVISION_2026-09-30.md)。

状态：实现、相关自动回归、完整真实数据 Notebook 及浏览器核对已完成。最终截图已检查，待用户人工审查。

以上状态及下文首轮数字对应最初的 FlowPlot 修订。用户随后要求提高文字可读性、区分颜色和线宽图例；最新实现与验证见本文末尾的“文字与图例可读性修订”，首轮代码快照不代表后续版本。

## 参考与改动

用户指定的 [MEBOCOST Demo](https://github.com/kaifuchenlab/MEBOCOST/blob/main/Demo_one_sample_mCCC_analysis.ipynb) 中，绘图 In[13] 调用 `mebo_obj.FlowPlot(...)`。本次参考其圆点、直线和数值编码，修订包内独立绘图及二级 HTML。实现位置为 `cellmesh/_communication_paths.py`、`cellmesh/_report_paths.py`；公开入口仍为 `plot_communication_paths()`，无新增必填参数。

| 内容 | 修订后 |
|---|---|
| 四层布局 | sender → metabolite → sensor → receiver，圆点节点，相邻列直线连接；节点较少的列居中 |
| 节点颜色 | 区分发送细胞、代谢物、受体/传感器、接收细胞四种角色 |
| 节点面积 | 当前页涉及该节点的完整事件数，经面积尺度映射；共享事件逐条计数 |
| 连线宽度 | 原 `cell_mesh_score` 经线性尺度映射 |
| 连线颜色 | 当前选择采用的统计列取 −log10 后映射；独立绘图可用 `fdr_col` 指定 |
| 统计值相同 | 同一颜色，图例明确说明相同原值 |
| 零值与缺失 | 统计值为零单独着色，不替换原值；统计缺失为灰色虚线，评分缺失为固定宽度虚线 |
| 尺度 | 取完整图形范围的最小/最大值，翻页和搜索不重新归一化；节点计数仍对应当前页 |
| 原事件定位 | 每条路径保留完整事件编号，支持节点高亮、定位事件、列表及已有三级报告链接 |

共享节点和线段可以重叠，完整事件以事件列表和 CSV 为准。图中关系不能任意重组为新的 sender–receiver 事件；连线不代表实测通量。采用 `fdr_sensor_type` 时，统计强弱在同一 sensor_type 内解读。

每页 50 条、组合筛选、搜索、`plot_top_n` 范围、两端表达图和三级链接保持原有规则。图形展示不追加显著性或 QC 过滤，不更改阳性判定、BH 校正或科学计算。报告图内嵌在 HTML 中，未增加默认图片文件。

## 真实数据审查入口

- [执行 Notebook](../../cell_mesh_validation/flowplot_revision_2026-09-30/flowplot_review.ipynb)
- [pooled 二级报告](../../cell_mesh_validation/flowplot_revision_2026-09-30/pooled/analysis.selected.report.html)
- [sample_aware 二级报告](../../cell_mesh_validation/flowplot_revision_2026-09-30/sample_aware/analysis.selected.report.html)

沿用此前完整 GSE176078 真实数据：100,064 个细胞、29,733 个基因、26 位患者、8 类细胞。复用完整计算结果与 1,000 次置换、随机种子 0 等原参数，阳性标准仍为 `fdr_sensor_type <= 0.05`。原 selection 为 `HMDB0000961–CYSLTR1`，pooled 28 条、sample_aware 26 条；两模式完整事件表均有 246,120 条。

Notebook 将上轮完整交付复制到独立目录，核对复制前后哈希；显式加载原完整输入，再通过包的 `result.add_selection(..., replace=True, adata=adata)` 更新二级报告。一级、三级、输入归档和 CSV 复用原文件，H5 仅更新交付配置。未重新计算评分或置换，未抽样、更改参数或调整 selection。

独立 PNG 仅用于本次审查，在 Notebook 中显式调用 `plot_communication_paths()` 后保存，不是默认交付的新文件。原案例少于 50 条、仅一个代谢物—受体组合，分页、多组合、零值和缺失分支由自动化边界测试覆盖，不修改真实案例制造这些情形。

## 本次验证

相关回归 62 项通过；WebKit 浏览器边界回归 1 项通过，覆盖 105 条事件、翻页、跨页定位、对象切换及 1,400/420 像素视口。上述数字为此次实际执行结果；未将上一轮完整回归数字计入本次。

实际截图检查发现通用 SVG 样式会拉宽图例符号、挤压数值，已为图例指定固定宽高，并增加浏览器尺寸与文字溢出断言；修正后重新执行相关回归及真实数据 Notebook。

真实 Notebook 核对完整事件表、参数、原始评分/BH、源交付文件哈希，以及新图的圆点、直线、线宽/线色、事件联动、表达图及三级链接。检查生成报告期间包代码哈希不变。此次未重复上轮耗时的 H5 全节点比较；使用追加接口更新报告，并逐项核对完整事件表与参数。

最终 Notebook 的 4 个代码单元均已执行，0 个错误输出；61 项数据与交付检查、207 项真实报告浏览器断言全部通过。两个模式的报告截图均已检查，图例符号及数值完整显示，没有内部滚动条。包 Python 文件与生成报告时的 SHA-256 快照一致。

- [最终汇总](../../cell_mesh_validation/flowplot_revision_2026-09-30/summary.json)：`status=PASS`，`human_review=pending`。
- [pooled 报告截图](../../cell_mesh_validation/flowplot_revision_2026-09-30/checks/pooled.flow.png)、[sample_aware 报告截图](../../cell_mesh_validation/flowplot_revision_2026-09-30/checks/sample_aware.flow.png)。
- 执行日志：`execution.log`；回归日志：`checks/cellmesh-flow-tests-final.log`、`checks/cellmesh-flow-browser-tests.log`；浏览器断言：`checks/browser_results.json`。

上述截图已内嵌于执行 Notebook。人工审查结论留待用户确认。

## 文字与图例可读性修订

用户确认计划后实施。HTML 主标签增至约 16 px、HMDB 副标签约 13 px、列标题约 18 px；代谢物名称和编号独立成行，长标签按词换行，必要时拆分超长词，保留完整标签。行高和图形高度随标签调整，交互图根据容器实际宽度维持文字尺寸。静态图主标签为 11.5 pt，副标签为 9.5 pt，列标题为 13 pt。

颜色图例改为矩形连续色条和最小/中间/最大刻度；非零统计值相同则用单个实心色块，同时给出原统计值和 −log10。零值与缺失值分别使用独立色块。灰色短线仅用于评分线宽，空心圆用于整数事件数。图例通常使用三位有效数字，必要时增加精度区分相邻刻度，小数极小时采用科学计数法。只调整显示格式，事件原始数值和详情精度不变。详细读图说明折叠展示。

本轮相关回归 63 项通过，WebKit 自动边界测试 1 项通过，覆盖连续色条、常量/零值/缺失值、长标签及 1,400/420 像素窗口。新增断言检查实际渲染字号、标签边界、标签之间的遮挡、图例形状和尺寸；纯自动化边界数据不用于人工验收。

真实数据审查在同一 Notebook 分两阶段执行：先检查 pooled 的静态图和实际 HTML 截图，再运行 sample_aware，避免图形问题在耗时的完整 H5 加载之后才被发现。沿用首轮完整数据、参数、阳性标准和原 selection，仅通过包的 `add_selection()` 重新生成二级报告，独立图由 `plot_communication_paths()` 生成；历史目录保持不变。

- [本轮 Notebook](../../cell_mesh_validation/flowplot_typography_2026-09-30/typography_review.ipynb)
- [pooled 二级报告](../../cell_mesh_validation/flowplot_typography_2026-09-30/pooled/analysis.selected.report.html)、[sample_aware 二级报告](../../cell_mesh_validation/flowplot_typography_2026-09-30/sample_aware/analysis.selected.report.html)
- [pooled 截图](../../cell_mesh_validation/flowplot_typography_2026-09-30/checks/pooled.flow.png)、[sample_aware 截图](../../cell_mesh_validation/flowplot_typography_2026-09-30/checks/sample_aware.flow.png)

最终状态：两种模式的真实数据与截图核对均通过。Notebook 的 4 个代码单元分两阶段执行并保存输出，0 个错误；61 项数据与交付检查、279 项浏览器断言通过。原全部事件各 246,120 条，selection 仍为 pooled 28 条、sample_aware 26 条；参数、原评分/BH、CSV、输入归档和历史文件核对一致。两份报告的包代码快照与最终代码一致。

[本轮结果汇总](../../cell_mesh_validation/flowplot_typography_2026-09-30/summary.json) 为 `status=PASS`、`human_review=pending`。执行日志分别为 `execution_pooled.log`、`execution_sample_aware.log`；自动回归日志保存在 `checks/cellmesh-typography-tests.log` 和 `checks/cellmesh-typography-browser.log`。人工审查待用户完成。
