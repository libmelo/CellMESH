# 静态通讯事件矩阵点图修订

> 本页为上一轮矩阵分图的历史记录。用户指出坐标组合及拆图不符合要求，已按 [复合坐标轴修正](DOTMAP_AXES_REVISION_2026-09-30.md) 改为 sender–metabolite × receiver–sensor 的单张矩阵。本页测试数字不代表修正后的验收。

用户指定参考 [MEBOCOST Demo In[12] 的 `commu_dotmap()`](https://github.com/kaifuchenlab/MEBOCOST/blob/main/Demo_one_sample_mCCC_analysis.ipynb)，停用四层通讯关系图。本次不修改算法、原评分、p/BH 校正值、筛选标准或默认交付文件种类。

## 最终展示

- 二级报告开头为“通讯事件矩阵点图”，合并原四层图和报告内旧点图的位置；通讯网络图及其他适用图形保留。
- 横轴为发送 → 接收细胞对，纵轴为代谢物—sensor，保留代谢物名称和 HMDB 编号。每个位置唯一对应完整事件，没有共享线段、事件合并或评分平均。
- 与参考示例一致，点面积映射通讯评分，颜色映射所选统计量的 `−log10`。使用顺序色图 `viridis`，保留最小可见点面积。评分和统计范围在整个展示集合内固定，所有分图共用。
- 统计值相同时使用单色块并注明原值；真实零统计值单独着色，不补伪数值；缺失统计值用灰色，缺失评分用 ×。零评分仍有可见点。
- 每页最多 16 个细胞对列、20 个代谢物—sensor 行，空分块跳过；全部静态分图在 HTML 中连续内嵌，无选择器、悬停依赖、交互翻页或绘图 JavaScript，无图片内部滚动条。
- `plot_top_n` 仍仅限制代谢物—sensor 组合；所保留组合内的所有完整事件各画一次，不另加显著性或 QC 过滤。列表和 CSV 不受绘图上限限制。空白不代表不存在生物学通讯。
- 网络与计数图继续不纳入缺失评分的事件；矩阵保留这些事件并明确标记。原事件列表搜索和独立三级报告链接保留。

## 调用和迁移

```python
from cellmesh import plot_communication_dotmap

view = plot_communication_dotmap(focus, page=1)
view["fig"].savefig("communication_dotmap.svg", bbox_inches="tight")
# n_pages 给出静态页数；完整独立导出需依次请求所有页。

result.save("analysis", selection=focus, adata=adata)
# 或先输出通用结果，再追加/替换 selection：
result.add_selection(focus, replace=True, adata=adata)
```

原临时接口 `plot_communication_paths()` 及四层图实现已移除，改用 `plot_communication_dotmap()`。`plot_event_dotplot()` 仍保留为独立的 sender 分面图，沿用其原来的颜色和面积编码；二级报告不再重复调用它。

核心文件：`cellmesh/_communication_dotmap.py`（范围、分图、尺度、静态绘图）、`cellmesh/_report_dotmap.py`（离线 SVG 内嵌）、`cellmesh/_report_plots.py`（报告整合）、`cellmesh/plotting.py`（公开接口）。

## 本轮验证

自动化：102 项相关回归测试通过，另 1 项 WebKit 浏览器测试通过。覆盖分页不重不漏、完整身份、固定尺度、缺失/零/常量统计、原始 p 值、转义、完整图片显示、窄屏、禁用页面 JavaScript 后的静态展示、表达图切换和报告追加。合成数据只用于这些自动化边界测试。

真实验收入口：[已执行 Notebook](../../cell_mesh_validation/dotmap_revision_2026-09-30/dotmap_review.ipynb)。4 个代码单元均有执行输出，0 个错误；执行计数反映 pooled / sample_aware 分阶段运行。

沿用完整 GSE176078（100,064 细胞 × 29,733 基因）、原先验、1000 次置换、随机种子 0、原阳性标准 `fdr_sensor_type ≤ 0.05`，以及原 `HMDB0000961–CYSLTR1` selection。**复用完整原计算结果，重新生成交付**；没有重算评分或置换。Notebook 分两阶段实际执行，先 pooled 后 sample_aware，每个模式输出静态图、报告截图和断言结果。独立 PNG/PDF 是 Notebook 显式导出，不是默认交付新增文件。

两模式均已完成，合计 **67 项数据与交付检查、248 项真实报告浏览器断言通过**，两份报告均与当前包代码 SHA-256 匹配。报告截图和静态 PNG 已检查；最终人工审查由用户进行。

| 模式 | 全事件 | 所选事件 | 静态分图 | 浏览器断言 |
|---|---:|---:|---:|---:|
| pooled | 246,120 | 28 | 2 | 127 |
| sample_aware | 246,120 | 26 | 2 | 121 |

- [pooled 报告](../../cell_mesh_validation/dotmap_revision_2026-09-30/pooled/analysis.selected.report.html)
- [sample_aware 报告](../../cell_mesh_validation/dotmap_revision_2026-09-30/sample_aware/analysis.selected.report.html)
- [验收汇总](../../cell_mesh_validation/dotmap_revision_2026-09-30/summary.json)
- 独立静态图位于各模式目录：`communication_dotmap_1.png`、`communication_dotmap_2.png`、两页 `communication_dotmap.pdf`。

范围说明：核对源文件 SHA-256、原参数、完整 events 表、选择事件逐项数值、报告内事件集合和独立图尺度，以及 CSV、输入、一级与三级交付未变。本轮不重复执行之前的完整 H5 所有内部节点对比；也不把历史测试计入本轮通过数。
