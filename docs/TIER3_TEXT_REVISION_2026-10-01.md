# 三级报告文字与结构修订

按用户确认方案，保留事件判定与全部明细，精简默认正文和重复说明。公开调用方式及三级报告与 selection 的独立关系不变。

| 部分 | 修订 |
|---|---|
| 事件概览 | 使用中文分析模式；保留路径、评分、统计量、实际达标标准及异常原因；事件编号移至来源 |
| 评分组成 | 正文只保留两端评分的几何平均；P/C/E、丰度加权和标准化方法折叠；计算输入改为中文字段列表，悬停可查看原值和字段名 |
| 样本表现 | 合并正评分比例及 QC 说明；保留 NA 与 0 的区别、中位数虚线及全部样本；明确正评分不等于显著 |
| 表达与覆盖 | 删除重复小标题；表达范围只说明一次；逐图保留形状和颜色含义；特殊读图规则折叠；覆盖度仍只表示输入可用情况 |
| 数据与来源 | 表格操作提示只出现一次；基因、反应等表格按“记录”计数；保留完整字段、记录和精度；去掉未实际导出的 CSV 承诺 |

`_report_pages.py` 负责三级结构，`_report_evidence.py` 负责表达/样本文案，`_report_browser.py` 为完整明细提供紧凑表格。一级、二级表格仍保留其事件计数与操作提示，共用分页状态更新为简洁表述。

## 真实数据验收

原文案验收读取 `dotmap_axes_2026-09-30/` 的完整 H5 与原表达输入，沿用 pooled、sample_aware 各两个原单事件请求，以独立 `to_html(events=...)` 生成三级报告。该执行记录现归档于 `checks/history/tier3_text_review.ipynb`。

用户随后确认恢复并保持统一交付结构。当前入口为包外 `../cell_mesh_validation/tier3_text_review_2026-10-01/tier3_text_review.ipynb`，实际执行目录整理、字节校验和导航检查。每种模式集中保留完整归档、一级、二级和新版三级报告，后续在完整交付副本上使用保存或追加接口。

检查包括：原参数与源文件 SHA256、每页十张表的全部原字段和记录、有效报告链接、样本切换、表格搜索/排序/分页、无 JavaScript 静态预览和真实浏览器截图。人工验收内容在 Notebook 中执行并保留输出。

两种模式各 2 份报告已重新生成，四份报告的 40 张内嵌明细表与原报告逐项一致；每种模式仍有 246,120 个完整事件，原参数和源文件校验通过。沿用完整 GSE176078 输入（100,064 细胞 × 29,733 基因）、1000 次置换的既有结果、随机种子 0 及原达标标准，未重新执行计算。

真实浏览器交互断言共 367 条通过（pooled 152、sample_aware 215），包括完整记录检索、排序控制、分页和 26 个样本的切换；四页无 JavaScript 静态预览通过。已查看默认页面、展开计算输入、表达图、完整样本图及展开明细截图。最终人工审查仍待用户进行。

原文案验收 Notebook 共 7 个代码单元保留实际执行输出、0 个错误。主体与补充截图分两阶段执行，因此执行计数不连续；日志分别为 `execution.log` 和 `execution_supplement.log`。这些属于已归档的历史验证，不计为目录整理时重新执行。

- [执行 Notebook](../../cell_mesh_validation/tier3_text_review_2026-10-01/tier3_text_review.ipynb)
- [pooled 统一导航](../../cell_mesh_validation/tier3_text_review_2026-10-01/pooled/analysis.reports/index.html)
- [sample_aware 统一导航](../../cell_mesh_validation/tier3_text_review_2026-10-01/sample_aware/analysis.reports/index.html)
- [验收汇总](../../cell_mesh_validation/tier3_text_review_2026-10-01/summary.json)
- [样本表现截图](../../cell_mesh_validation/tier3_text_review_2026-10-01/checks/sample_aware.samples.png)

## 统一交付恢复

每种模式的一级 `analysis.report.html`、二级 `analysis.selected.report.html` 位于模式根目录；三级位于 `analysis.reports/events/`；统一导航为 `analysis.reports/index.html`。旧的 `events.html` 和 `events.reports/` 独立导出布局已移除，相关 Jupyter 自动检查点已归档保留。

H5、H5AD、CSV、一级和二级报告从原完整交付复制，逐文件 SHA256 与已验收基线一致；四份新版三级正文、SVG 和内嵌数据均保留，只修正跳转链接。导航由包内 `render_navigation()` 生成，注明一级/二级沿用及三级采用本轮文字修订版。全部报告、CSV 及导航链接位于同一模式目录内，原历史交付未被修改。

当前目录整理 Notebook 的 4 个代码单元均已执行、0 个错误；十个 HTML 页面的 58 条浏览器断言通过。另有 12 项三级报告回归测试通过。本轮只验证目录与链接，不重算、不新增全量 CSV，也未重新生成一级、二级图形。首次整理因发现自动检查点而停止，随后先归档检查点再完成迁移，首次执行记录保存在 `checks/history/unification_attempt1.*`。

## 补充回归发现

119 项相关 Python 回归通过（首轮 117 项，更新两处沿用旧文案/旧标签的测试断言后，定向重跑 2 项通过）。共用事件表的空表和 10,001 条记录浏览器检查通过。

额外执行的二级矩阵图浏览器边界测试有 1 项失败：超宽矩阵的 `BH correction: within each sensor type` 图例文字在桌面视口超出 SVG 右边界约 4.17 px，窄视口约 2.61 px。该图形模块本轮未改动；此问题单独保留，不能将本轮记录表述为全包或二级图形验收全部通过。实际日志及边界诊断位于包外验收目录的 `checks/`。
