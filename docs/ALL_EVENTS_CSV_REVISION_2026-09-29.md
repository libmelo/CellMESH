# 全量事件 CSV 按需导出

状态：修订及本轮自动验收已完成，待用户审查。完整自动回归与开关边界测试通过；pooled/sample_aware 默认均不生成全量 CSV，显式 pooled 全量 CSV 与原完整结果一致，最终科学 H5 核对通过。详见 [统一验收记录](COMMUNICATION_FEATURES_ACCEPTANCE_2026-09-29.md)。

用户重申：只有明确指定参数，才生成全量事件 CSV。此前代码及部分文档仍是默认导出，现予以修正。

## 调用与文件

```python
files = result.save("analysis")
# analysis.cellmesh.h5：完整科学结果
# analysis.positive_events.csv：达标事件，默认 fdr_sensor_type <= 0.05
# analysis.report.html：一级报告
# analysis.reports/index.html：报告导航

files = result.save("with_all", include_all_events=True)
# 额外生成 with_all.all_events.csv；files["events_csv"] 返回其路径

result.to_html("general.html")  # 一级 HTML + 阳性 CSV
result.to_html("full.html", include_all_events=True)  # 再增加全量 CSV
```

- `save()`、一级 `to_html()` 的 `include_all_events` 默认均为 `False`，只接受布尔值；加载旧结果或重新保存时也不隐式开启。
- `report=False` 只控制 HTML，阳性 CSV 仍输出；全量 CSV 仍由显式开关决定。
- 指定 selection 增加所选 CSV 和二级报告；指定 events 增加独立三级页面。二者均不会自动开启全量 CSV。
- 独立 `to_html()` 若同时指定 selection/非空 events 和 `include_all_events=True`，会在写文件前报错，因为该调用没有一级报告。需要一次交付多个层级与全量 CSV 时使用 `save()`。
- H5 的完整科学结果、一级概览的全部事件数量、阳性判定与评分算法均不因该开关变化。
- 返回字典与一级 HTML 只包含本次实际交付的 CSV 路径和下载链接。
- 显式调用旧入口 `to_csv(prefix)` 的多表导出行为保持不变。

## 追加与兼容

- H5 交付元信息记录 `include_all_events`；`add_selection()`、`add_event_reports()` 继承该设置。关闭时不会补写全量 CSV；此前明确开启过的交付可在追加时恢复缺失的全量 CSV。
- 旧交付没有此字段时按关闭处理。首次追加刷新一级 HTML 的下载链接，保留已有全量 CSV；后续追加沿用记录的开关。
- 重新 `save()` 默认关闭，不沿用旧值。未开启时，不删除或覆盖同名旧全量 CSV，不将它返回或链接为本次交付；目录里的旧文件不能据此视为本次新生成结果。
- 已有独立报告的覆盖保护、CSV 内容核对与普通 I/O 失败回滚继续沿用。

## 统一验证范围

1. 使用原完整真实数据与科学参数，在包外 `.ipynb` 中验证两种模式的默认文件、报告链接与全部事件计数。
2. 验证显式开启时全量 CSV 与完整事件表一致，并检查 `report=False`、空阳性、selection 和独立 events 的组合。
3. 验证保存/加载后的开关记录，默认关闭与显式开启交付的追加行为，旧交付迁移，以及同前缀已有 CSV 的保留。
4. 核对 H5 科学数据、计算参数、阳性标准与源文件未变；明确记录复用原计算结果、重新生成交付。
5. 同步调整旧自动化断言中“默认必有 events_csv”的预期，补充相关边界检查并统一执行。本轮已更新断言并执行测试及 Notebook；最终完成状态见统一验收记录。

此前验收目录和历史示例输出保留原状，不代表本次代码行为。
