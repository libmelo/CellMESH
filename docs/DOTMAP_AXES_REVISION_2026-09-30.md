# 单张矩阵与复合坐标轴修正

> 后续字号、排版与色条修订已进入代码，但依用户要求暂不校验。本文的验收记录对应此前坐标轴修正版，不能作为当前可读性修订的验收。见 [可读性修订](DOTMAP_READABILITY_2026-09-30.md)。

用户指出上一轮坐标轴和自动拆图没有实现所需展示。本轮明确采用：**横轴 sender–metabolite，纵轴 receiver–sensor，整个展示范围绘在同一张静态图里**。

## 实现

- 横轴内部键为 `(sender, hmdb_id)`，标签显示发送细胞、代谢物名称、HMDB 编号。
- 纵轴内部键为 `(receiver, sensor_gene)`，标签显示接收细胞和 sensor。
- 每个完整事件对应唯一 `(x, y)`，多个代谢物和 sensor 共用同一个坐标系；不合并事件，不平均评分，不因代谢物同名而合并 HMDB。
- 移除自动分页及按组合分图逻辑；画布随复合轴数量扩大。报告只内嵌一张完整 SVG，无交互控件或绘图脚本依赖。
- 保留点面积映射评分、颜色映射所选 p/BH 的 `−log10`。常量统计值、零值、缺失值的区分沿用前版。
- `plot_top_n` 仍仅限制组合展示范围；完整 selection CSV、事件列表、阳性标准、网络图及其它模块不变。

参考来源：[MEBOCOST Demo In[12] 的 `commu_dotmap()`](https://github.com/kaifuchenlab/MEBOCOST/blob/main/Demo_one_sample_mCCC_analysis.ipynb)。该参考用于点面积/颜色编码；其原示例采用细胞对与代谢物–sensor 轴。本轮采用的 sender–metabolite / receiver–sensor 复合轴是用户明确要求的改编，不将其误述为参考 Notebook 的原始轴定义。

## 调用

```python
from cellmesh import plot_communication_dotmap

view = plot_communication_dotmap(focus)
view["fig"].savefig("communication_dotmap.pdf", bbox_inches="tight")
# x_keys: [(sender, hmdb_id), ...]
# y_keys: [(receiver, sensor_gene), ...]
# dot_table: 每个事件的完整身份、x/y 索引、面积、颜色及符号。

result.save("analysis", selection=focus, adata=adata)
# 或更新已有交付中的二级报告：
result.add_selection(focus, replace=True, adata=adata)
```

移除上一轮临时引入的 `page`、`max_cell_pairs`、`max_pairs` 绘图参数，以及 `n_pages`、`cell_pairs` 返回字段。改用 `x_keys`、`y_keys` 查看复合轴身份。函数名称仍为 `plot_communication_dotmap()`，一轮调用返回完整单图。

## 本轮检查

102 项相关回归测试和 1 项 WebKit 浏览器测试通过。新增行为检查覆盖多个代谢物与多个 sensor 共用单图、明确给定的预期坐标、同名代谢物的 HMDB 区分、每个事件恰好出现一次，以及浏览器中实际点的位置与复合轴的对应关系。关闭页面 JavaScript 后仍完整显示。合成数据仅用于自动化边界检查。

真实数据验收入口：[执行 Notebook](../../cell_mesh_validation/dotmap_axes_2026-09-30/dotmap_axes_review.ipynb)。沿用完整 GSE176078（100,064 细胞 × 29,733 基因）、原先验、1000 次置换、随机种子 0、`fdr_sensor_type ≤ 0.05` 以及原 `HMDB0000961–CYSLTR1` selection。**复用完整计算结果，重新生成交付**，不重新计算评分或置换，不为展示效果更换 selection。

两模式均已完成。Notebook 的 4 个代码单元均保存实际执行输出，0 个错误；分阶段执行计数为 1、2、2、3。合计 **75 项数据与交付检查、424 项真实报告浏览器断言通过**。两模式报告和静态 PNG 已检查，最终人工审查待用户进行。

| 模式 | 完整事件表 | 所选事件点 | 矩阵 | 浏览器断言 |
|---|---:|---:|---|---:|
| pooled | 246,120 | 28 | 1 张，8 列 × 4 行 | 218 |
| sample_aware | 246,120 | 26 | 1 张，8 列 × 4 行 | 206 |

- [pooled 报告](../../cell_mesh_validation/dotmap_axes_2026-09-30/pooled/analysis.selected.report.html)、[静态图](../../cell_mesh_validation/dotmap_axes_2026-09-30/pooled/communication_dotmap.png)
- [sample_aware 报告](../../cell_mesh_validation/dotmap_axes_2026-09-30/sample_aware/analysis.selected.report.html)、[静态图](../../cell_mesh_validation/dotmap_axes_2026-09-30/sample_aware/communication_dotmap.png)
- [验收汇总](../../cell_mesh_validation/dotmap_axes_2026-09-30/summary.json)

每模式的 PNG 和单页 PDF 是 Notebook 显式导出，默认交付文件种类不增加。源码哈希与两份实际生成报告一致，原历史交付文件均未改写。

本轮核对参数、完整 events 表、逐点身份及原始统计值、源文件 SHA-256、保存后事件表、保留的 CSV/输入/一级和三级文件，以及报告与独立图的复合轴和尺度。没有重复执行完整 H5 所有内部节点的对比，也不把历史通过记录算作本轮检查。
