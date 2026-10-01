# 报告验收交付约定

沿用上级目录关于包外 Notebook 和真实数据验收的约定。

用户于 2026-10-01 确认统一交付，并进一步批准独立的 `filter_events()` / `select_events(metabolite_sensor_pairs=...)` 入口及二级文件名改为 `filtered`。

- 完整验收按 `pooled/`、`sample_aware/` 分别保存，沿用包内 `bundle_paths()` 的命名：一级 `analysis.report.html`、二级 `analysis.filtered.report.html`、统一导航 `analysis.reports/index.html`、三级 `analysis.reports/metabolite_sensor/pair_*.html`。
- 只修订其中一个层级时，仍保留当前模式已生成的其它层级及实际启用的归档和 CSV；报告链接指向当前交付目录，不跨历史验收目录查找一级或二级页面。无需生成原来未请求的层级或可选输出。
- 导航应说明本轮更新与沿用的报告。整理路径不等于重新生成内容，不得把沿用的报告记为本轮重新运行。
- 后续修订在完整交付副本上使用包的保存或追加接口。独立 `selected.to_html(...)` 可用于临时查看，不能取代完整验收的统一交付入口。
- 人工审查 Notebook 展示当前统一入口，保留可追溯的历史执行记录；迁移文件时同步检查导航、CSV、二级至三级链接和内嵌数据。
- 当前交付根目录不得混放被替代的旧名 `*.selected.report.html` / `*.selected_events.csv`。保存和追加接口自动将它们完整归档至 `<prefix>.reports/history/legacy-selection-<摘要>.zip`；已迁移交付可调用 `cellmesh.result_io.archive_legacy_reports(prefix)` 整理。验收同时检查残留旧文件、归档内容及当前报告链接；保留历史不等于将旧报告继续放在当前结果旁。
