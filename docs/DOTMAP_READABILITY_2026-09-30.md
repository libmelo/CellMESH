# 通讯矩阵点图：可读性修订

状态：按最新确认，已交换横纵轴、改用单行标签，横轴固定 45°、纵轴保持水平，代谢物默认显示 HMDB ID，可通过绘图参数改用名称。沿用上一轮批准的同一组 52 条真实事件，**本轮只重新绘制一张图片供人工核验**。未运行 pytest、浏览器测试或完整报告验收，未重新生成 HTML/H5/CSV。待其它修订完成后一起校验。

2026-10-01 追加：标题润色为 `Selected cell–cell communication events (n = …)`，相对矩阵居中，标题间距由 14 pt 增至 28 pt；删除评分图例中的 `Dot area`，同时收紧对应空白。BH 图例明确标注校正范围，避免把 `within sensor type` 误解为颜色或评分归一化的范围。

## 原因与修正

上一版把全图共有的代谢物和 sensor 移至图头，虽然内部坐标键不变，但刻度只剩细胞类型，未落实直接展示完整复合身份的约定。旧的新增测试也错误地要求省略共有代谢物，现已修正这些断言，但尚未执行。

- 每个横轴刻度使用单行 `receiver – sensor`；每个纵轴刻度默认使用单行 `sender – HMDB ID`。全图共有的身份也完整显示，不用图头替代。
- `plot_communication_dotmap(..., metabolite_label="id")` 为默认行为；传入 `"name"` 使用代谢物名称，缺名时回退 ID，同名不同 HMDB 时附加 ID 区分。名称选择只影响显示，报告默认使用 ID。
- 一张静态矩阵包含完整展示范围：横轴 `(receiver, sensor_gene)`，纵轴 `(sender, hmdb_id)`。每个完整事件对应一个点，不拆图、不分页、不合并或平均事件。
- 原评分、统计值、面积/颜色映射规则和过滤逻辑均不因排版变化而改变；本次预览的 selection 变更单独记录如下。

## 字体与排版

| 对象 | 当前设置 |
|---|---|
| 字体 | 统一无衬线字体 |
| 图标题 | 14 pt，半粗体 |
| 坐标轴标题 / 刻度 | 12 pt |
| 图例标题 / 数值 | 11.5 pt / 11 pt |
| 次要说明 | 10 pt |

- 复合标签均为单行，横轴固定 45°，纵轴水平；不自动换行或改成 90°。按倾斜标签的垂直厚度安排列距，长标签向下延伸，避免按外接矩形宽度逐列扩张。
- 矩阵高度由行数及标签尺寸决定，图例单独安排空间，不再强制把少行矩阵拉高到图例高度。
- 自动画布容纳标题、完整刻度和右侧图例，不靠不断缩小字体容纳更多对象。传入外部 `ax` 时保留其位置，由调用者预留周边空间。
- HTML 仍按页面宽度等比缩放；极大矩阵可使用独立 SVG/PDF 阅读。

## 图例与报告说明

- 点面积映射原通讯评分，评分图例使用中性灰色圆点，仅注明 `Communication score`。
- 多个不同的正 p/BH 值使用连续矩形色条，宽 0.28 英寸、长 1.8 英寸；标明 `−log10(p/BH)`，显示 3–5 个刻度，并说明原始 p 或 BH 校正范围。
- 图例使用 `BH correction: within each sensor type` / `BH correction: all tested events` 明确限定 p 值校正范围；原始 p 值的多值标题为 `−log10(raw p)`，不附加 BH 说明。
- 单一正值使用矩形色块，标注捕获的原始概率，不制造梯度。零值单独着色并注明原始 p/BH 为 0，统计缺失使用灰色，评分缺失保留 ×。
- 图例数字通常保留三位有效数字，必要时增加精度区分相近刻度；小概率不格式化成 0。
- 报告图前只保留范围、轴含义和面积/颜色编码；详细说明位于图下默认折叠的“读图说明与展示范围”。报告代码与说明已更新，报告文件本轮未重新生成。

连续矩形色条参考 [MEBOCOST In[12] 的 commu_dotmap](https://github.com/kaifuchenlab/MEBOCOST/blob/main/Demo_one_sample_mCCC_analysis.ipynb)。本包使用用户指定的 receiver—sensor / sender—metabolite 复合轴及 viridis 顺序色图。

## 本轮真实预览

复用原 GSE176078 pooled 完整计算的阳性事件表（11,607 条，`fdr_sensor_type <= 0.05`）；不重新计算。来源为 `cell_mesh_validation/dotmap_axes_2026-09-30/pooled/analysis.positive_events.csv`。

沿用上一轮 selection：`hmdb_id ∈ {HMDB0000696, HMDB0001414}` 且 `sensor_gene ∈ {SLC7A5, SLC3A2}`，保留该阳性表中满足条件的全部事件；发送/接收细胞不限，不按评分截取。

| 组合 | 事件数 | 不同 BH 值数 |
|---|---:|---:|
| L-Methionine（HMDB0000696）—SLC7A5 | 24 | 4 |
| Putrescine（HMDB0001414）—SLC3A2 | 28 | 3 |

共 52 条事件、换轴后为 12 列 × 13 行；sensor_type 均为 Transporter。BH 原值为 0.01085141427439716、0.021475298894653733、0.03210936972994851、0.042606993006993006。Notebook 记录来源哈希、代码哈希、选择代码、实际执行输出及图片；本次图片使用默认 ID 标签。

- [执行 Notebook](../../cell_mesh_validation/dotmap_title_legend_preview_2026-10-01/pooled_title_legend_preview.ipynb)
- [静态预览 PNG](../../cell_mesh_validation/dotmap_title_legend_preview_2026-10-01/pooled_title_legend_dotmap.png)

原 28 条事件、52 条事件复合标签及换轴预览保留在各自旧目录。本轮图片尚待用户核验，不作为完整验收通过记录；名称模式和边界情况尚未运行验证。

## 后续统一校验清单（本轮未执行）

1. 核对每个点的完整身份、评分、统计值、坐标与源结果一致。
2. 单/多代谢物、单/多 sensor 的单行刻度均显示完整复合身份；ID/名称模式、缺名回退及同名 HMDB 区分正确，名称切换不改变坐标、点面积或颜色。
3. 相邻标签、标题和图例不重叠；独立图及 HTML 缩放后的排版清晰且完整。
4. 多值色条与点颜色对应；常量、极小正值、零值、缺失值及其组合正确标注。
5. 少行矩阵不因图例拉高；完整报告、追加 selection、外部 ax、SVG/PDF 导出及离线浏览器行为符合约定。

相关实现：`cellmesh/_communication_dotmap.py`、`cellmesh/_report_dotmap.py`、`cellmesh/plotting.py`。待执行的自动化检查已同步至 `tests/test_communication_features.py` 和 `tests/report_communication_boundaries.js`。
