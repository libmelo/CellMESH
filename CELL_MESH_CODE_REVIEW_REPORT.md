# CELL MESH 包测试补充与代码审查报告

> 历史审查记录。本报告中的公式与 API 描述对应 median-contrast 改造前的实现，
> 不属于当前方法文档。现行方法以 `README.md`、`docs/METHODS.md` 和
> `AVAILABILITY_IMPLEMENTATION.md` 为准。

生成时间：2026-06-14  
审查对象：`developer/cell_mesh_pkg`  
Conda 环境：`cellmesh`  
包版本：`cellmesh 0.4.0`

---

## 1. 环境安装状态

已检查 `cellmesh` conda 环境中的安装情况：

- Python 路径：`/home/qsong/miniconda3/envs/cellmesh/bin/python`
- 可导入模块：`cellmesh`
- 模块路径：`/home/qsong/.openclaw/workspace/developer/cell_mesh_pkg/cellmesh/__init__.py`
- 包版本：`0.4.0`
- 安装方式：editable install，`direct_url.json` 指向当前源码目录 `file:///home/qsong/.openclaw/workspace/developer/cell_mesh_pkg`

结论：环境中安装的 `cellmesh` 已指向当前源码目录。为确保安装元数据与源码同步，已执行：

```bash
conda run -n cellmesh python -m pip install -e .
```

执行结果：成功重新安装 `cellmesh-0.4.0` editable 版本。

---

## 2. 本次测试补充

修改文件：

- `examples/test_cellmesh_comprehensive.ipynb`

未改动包源码文件，仅补充 notebook 中的展示和断言。

### 2.1 新增 pseudobulk 聚合一致性测试

新增 section：`新增测试 1：pseudobulk 聚合计算前后一致性`

测试逻辑：

1. 从 `adata.X` 取原始单细胞表达矩阵。
2. 按 `adata.obs['cell_type']` 手工分组。
3. 对每个 cell type 计算每个基因的均值，得到 `manual_pseudobulk`。
4. 与 `compute_metabolite_availability(..., return_intermediates=True)` 返回的 `avail['pseudobulk']` 对齐后比较。
5. 展示：
   - 聚合前单细胞表达矩阵片段；
   - 手动聚合后的 pseudobulk 片段；
   - 函数返回的 pseudobulk 片段；
   - 两者差值片段。

实际验证结果：

- `max abs diff: 0.0`
- `PASS pseudobulk aggregation consistency`

说明：手工聚合与函数内部 pseudobulk 数值完全一致。由于原始 `adata.X` 为 float32，而函数返回结果最终为 float64，断言使用 `check_dtype=False`，只验证数值一致性。

### 2.2 增强 TC-A2 展示：availability 公式前后数据

原 TC-A2 只做公式断言。本次增强为：

- 展示计算前输入片段：`P_norm` / `C_norm` / `E_norm`；
- 手工计算 `manual_availability`；
- 展示函数输出 `function_availability`；
- 展示 `abs_diff`；
- 使用 `pd.testing.assert_frame_equal` 对完整矩阵做严格数值断言。

公式：

```python
availability = P_norm * ((1 - C_norm) ** 0.5) * (0.8 + 0.2 * E_norm)
```

实际验证结果：

- `max abs diff: 0.0`
- `PASS TC-A2`

### 2.3 增强 TC-B9 展示：参数传递与自定义 availability 前后数据

原 TC-B9 验证参数传入和自定义公式。本次增强为：

- 展示默认参数片段：`beta`, `missing_C_norm`, `missing_E_norm`；
- 展示自定义参数：`beta=1.0`, `missing_C_norm=0.41`, `missing_E_norm=0.75`；
- 展示自定义参数实际进入 `result.parameters` 的值；
- 展示自定义参数下的 `P_norm` / `C_norm` / `E_norm` 输入片段；
- 手工计算 `manual_custom_availability`；
- 展示函数输出 `function_custom_availability`；
- 展示默认参数结果 vs 自定义参数结果，以及差值 `delta_custom_minus_default`。

自定义公式：

```python
availability = P_norm * ((1 - C_norm) ** beta) * (0.8 + 0.2 * E_norm)
```

实际验证结果：

- `max abs diff: 0.0`
- `PASS TC-B9`

---

## 3. 验证结果

### 3.1 pytest 验证

执行：

```bash
conda run -n cellmesh pytest -q tests/test_cellmesh_comprehensive.py
```

结果：

```text
18 passed, 2 skipped, 105 warnings in 16.79s
```

警告主要来自 `anndata` 读取旧格式 h5ad 文件的 `OldFormatWarning`，不是本次修改引入的测试失败。

### 3.2 notebook 逐 cell 轻量执行验证

环境中没有安装 `jupyter-nbconvert`，因此未使用 nbconvert。改用轻量执行脚本逐 cell 执行 notebook 代码，并跳过 Jupyter magic 行 `%matplotlib inline`。

结果：

```text
NOTEBOOK_OK
```

关键新增测试输出：

- pseudobulk 聚合：`max abs diff: 0.0`
- TC-A2：`max abs diff: 0.0`
- TC-B9：`max abs diff: 0.0`

---

## 4. 包结构与功能流分析

### 4.1 顶层结构

核心包目录：`cellmesh/`

主要模块：

- `cellmesh/__init__.py`：对外 API 汇总，暴露 `run_cell_mesh`, `compute_metabolite_availability`, `load_cell_mesh_database` 等。
- `cellmesh/config.py`：集中配置路径、默认参数、合法 role / sensor type 常量。
- `cellmesh/database.py`：加载和标准化内置或用户提供的 enzyme / interaction 数据库。
- `cellmesh/preprocess.py`：数据读取、先验验证、pseudobulk、sensor score、metabolite availability 计算。
- `cellmesh/core.py`：顶层 `run_cell_mesh`，事件构建、置换检验、置信等级、结果封装。
- `cellmesh/scoring.py`：少量通用评分函数，如 `sigmoid`, `zscore_by_gene`。

辅助目录：

- `cellmesh/data/`：内置 CSV 和测试 h5ad 数据。
- `examples/`：notebook 和脚本示例。
- `tests/`：pytest 测试。
- `docs/`：方法和 AnnData 说明文档。

### 4.2 主计算流程

`run_cell_mesh()` 是主入口，位于 `cellmesh/core.py`。

整体流程如下：

1. 若未提供先验，调用 `load_cell_mesh_database()` 加载默认 enzyme / interaction 数据。
2. 调用 `validate_priors()` 过滤并标准化先验，只保留表达矩阵中存在的 gene。
3. 调用 `_enzyme_prior_to_availability_reactions()` 将 enzyme prior 的 `role` 映射到内部 reaction `direction`。
4. 调用 `_compute_availability_scores()`：
   - `compute_metabolite_availability()` 计算 sender-side metabolite availability；
   - `compute_sensor_scores()` 计算 receiver-side sensor score。
5. 调用 `_make_cell_mesh_events()` 构造 sender × receiver × metabolite × sensor_gene 通信事件。
6. 可选调用 `_empirical_pvalues_by_sensor_type()` 做置换检验。
7. 调用 `_confidence_tier()` 给事件打置信等级。
8. 返回 `CellMeshResult`。

### 4.3 metabolite availability 流程

位于 `cellmesh/preprocess.py`。

关键步骤：

1. `_build_celltype_pseudobulk()`：按 cell type 聚合表达均值，得到 cell type × gene 矩阵。
2. `_parse_reaction_table()`：补齐 `hmdb_id`, `reaction`, `weight` 等列。
3. `_get_reaction_gene_sets()`：按 metabolite + hmdb_id + reaction + direction 合并反应基因集合。
4. `_compute_reaction_scores()`：对每个 reaction 使用几何均值计算反应得分。
5. `_compute_PCE_matrices()`：按 `direction` 聚合到：
   - `product` → P；
   - `substrate` → C；
   - `exporter` → E。
6. `_normalize_PCE()`：按代谢物维度 robust min-max 标准化 P/C/E。
7. `compute_metabolite_availability()`：按公式计算 availability：

```python
availability = P_norm * ((1 - C_norm) ** beta) * (0.8 + 0.2 * E_norm)
```

### 4.4 receiver / sensor score 流程

`compute_sensor_scores()` 中：

1. 再次构建 pseudobulk。
2. 计算每个 cell type 中每个基因表达比例 `expr_frac`。
3. 对每个 sensor gene 在 cell type 间做 robust min-max。
4. 若表达比例低于 `min_expr_frac`，对应 cell type 的 sensor score 置 0。
5. 输出 receiver-level 表格，包含 `metabolite`, `hmdb_id`, `sensor_gene`, `sensor_type`, `receiver`, `sensor_score`, `sensor_expr_frac`。

### 4.5 event score 流程

`_make_cell_mesh_events()` 中：

```python
communication_score = sqrt(metabolite_availability * sensor_score)
```

同时保留别名字段：

- `sender_score` = `metabolite_availability`
- `receiver_score` = `sensor_score`
- `cell_mesh_score` = `communication_score`

---

## 5. 主要代码审查发现

### 5.1 冗余或未使用设计/参数

#### 5.1.1 `eps` 在 availability 主公式中未实际使用

`compute_metabolite_availability()` 接受 `eps` 参数，并执行：

```python
eps = max(eps, 1e-8)
```

但后续公式中没有使用该局部 `eps`。P/C/E 标准化调用的是 `robust_minmax()` 的默认 `eps=1e-8`，没有将 `compute_metabolite_availability(eps=...)` 传入。

影响：

- 用户传入 `eps` 不会改变结果；
- `run_cell_mesh(..., eps=...)` 也只记录到 `parameters`，但对计算无实质影响。

建议：

- 如果设计上 `eps` 应控制 robust min-max 分母，应将其传入 `robust_minmax()` / `_normalize_PCE()`；
- 如果不需要暴露，建议从 public API 参数中移除或在文档中说明当前无效。

#### 5.1.2 `ROLE_TO_DIRECTION` 常量未被核心转换函数使用

`config.py` 定义了：

```python
ROLE_TO_DIRECTION = {
    'production': 'product',
    'degradation': 'substrate',
    'export': 'exporter',
}
```

但 `_enzyme_prior_to_availability_reactions()` 内部重新定义了同样映射：

```python
role_to_dir = {
    'production': 'product',
    'degradation': 'substrate',
    'export': 'exporter'
}
```

影响：

- 配置集中化没有完全落实；
- 若未来修改 role 映射，容易出现配置和实现不一致。

建议：使用 `config.ROLE_TO_DIRECTION`，删除函数内重复字典。

#### 5.1.3 `VALID_SENSOR_TYPES` 导入但未使用

`preprocess.py` 从 `config.py` 导入 `VALID_SENSOR_TYPES`，但 `validate_priors()` 中没有验证或过滤 sensor type，注释明确写着“保持原样”。

影响：

- 常量存在但逻辑不使用，读者会误解 sensor type 会被强校验；
- 与 `database.py` 的 `_normalize_sensor_type()` 有一定重叠。

建议：

- 若希望允许任意 sensor type，删除导入或把常量改为文档用途；
- 若希望限制三类 sensor type，则在 `validate_priors()` 中加入校验/归一化。

#### 5.1.4 `MIN_CELL_COUNT=100` 与实际默认 `min_cells=1` 不一致

`config.py` 顶层定义：

```python
MIN_CELL_COUNT: int = 100
```

但 public availability 默认参数中：

```python
METABOLITE_AVAILABILITY_DEFAULTS['min_cells'] = 1
```

`compute_metabolite_availability()` 和 `run_cell_mesh()` 默认使用 `METABOLITE_AVAILABILITY_DEFAULTS['min_cells']`，因此实际默认是 1。

影响：

- 配置含义不清：用户看到 `MIN_CELL_COUNT=100` 可能误以为默认会过滤少于 100 个细胞的 cell type；
- 内部 helper `_build_celltype_pseudobulk()` 默认 `MIN_CELL_COUNT=100`，但 public API 默认传入 1，导致直接调用 helper 与通过主 API 行为不同。

建议：统一默认值来源，至少在文档中明确 helper 默认值和 public API 默认值的差异。

#### 5.1.5 `MIN_EXPR_FRAC` 常量未被默认参数使用

`config.py` 定义：

```python
MIN_EXPR_FRAC: float = 0.05
```

但 `run_cell_mesh()` 和 `compute_sensor_scores()` 中直接写死 `min_expr_frac: float = 0.05`。

影响：配置集中化不完整。

建议：将默认值改为 `MIN_EXPR_FRAC`。

#### 5.1.6 `role_scores` 保留为空 dict

`CellMeshResult` 中仍有：

```python
role_scores: dict[str, pd.DataFrame]
```

但 `run_cell_mesh()` 返回时固定为：

```python
role_scores={}
```

注释说明“旧算法已删除”。

影响：

- API 向后兼容，但当前无实际内容；
- 新用户可能误以为仍有 role-level score 可用。

建议：

- 如果需要向后兼容，保留但在 README/API 文档明确标记为 deprecated；
- 未来主版本可移除。

#### 5.1.7 `scoring.py` 中函数目前只在测试兼容逻辑中使用

`sigmoid()` 和 `zscore_by_gene()` 位于 `cellmesh/scoring.py`，但当前主流程的 sensor score 已改为 robust min-max，不再使用 z-score/sigmoid 路径。

影响：

- 可能是旧算法遗留；
- 用户可能误解主算法使用 sigmoid/zscore。

建议：

- 若保留作为工具函数，应在文档中标为 utility；
- 若旧算法不会恢复，可考虑移除或移动到 legacy 模块。

#### 5.1.8 若干 import 冗余

静态检查发现潜在未使用 import：

- `cellmesh/config.py`：`Literal`
- `cellmesh/core.py`：`Literal`, `Union`, `warnings`, `MIN_CELL_COUNT`, `ROLE_TO_DIRECTION`, `load_default_priors`
- `cellmesh/preprocess.py`：`dataclass`, `VALID_SENSOR_TYPES`

建议：清理未使用 import，降低维护噪音。

### 5.2 可能影响正确性或用户预期的问题

#### 5.2.1 interaction 数据库规范化不支持 lower-case schema

`normalize_interaction_database()` 目前只读取：

```python
Gene_name, standard_metName, HMDB_ID, Annotation, Database source, Protein_name, Reference
```

不像 `normalize_enzyme_database()` 那样兼容 lower-case schema。

影响：如果用户提供已经接近内部格式的 interaction 表，例如 `metabolite`, `hmdb_id`, `sensor_gene`, `sensor_type`，该函数可能读不到 gene，返回空表。

建议：增加对内部 schema / lower-case schema 的兼容，或文档明确用户文件必须使用原始上传格式。

#### 5.2.2 `_to_dense()` 会把稀疏矩阵整体转稠密

`_build_celltype_pseudobulk()` 和 `_compute_celltype_expr_frac()` 都会调用 `_to_dense()`，将稀疏表达矩阵完整转换为 dense。

影响：

- 对大型 scRNA-seq 数据可能产生高内存占用；
- 与单细胞数据常见 sparse matrix 使用习惯不匹配。

建议：对 sparse 矩阵使用按 group 切片并保持 sparse 的均值/表达比例计算，避免全量 densify。

#### 5.2.3 pseudobulk 计算在 availability 与 sensor score 中重复执行

`_compute_availability_scores()` 调用 `compute_metabolite_availability()` 生成一次 pseudobulk，随后 `compute_sensor_scores()` 又重新构建 pseudobulk 和 expr_frac。

影响：

- 对大型数据有重复计算成本；
- 两处使用相同 celltype/layer/min_cells 参数，结果应相同。

建议：允许 `compute_sensor_scores()` 接收预计算的 `pseudobulk` / `expr_frac`，或在 `_compute_availability_scores()` 中缓存并复用。

#### 5.2.4 weighted geometric mean 目前不是标准“加权几何均值”

`_compute_reaction_scores()` 中：

```python
expr = pseudobulk[valid_genes].values * valid_weights.reshape(1, -1)
geo_mean = gmean(expr + 1, axis=1) - 1
```

这相当于先把表达量乘以权重，再做普通几何均值；不是数学意义上的 weighted geometric mean（通常是 `exp(sum(w * log(x)) / sum(w))`）。

影响：

- 注释“加权几何平均数”可能不准确；
- weight 的含义是表达缩放，而不是几何均值权重。

建议：明确 weight 语义。如果要做真正的 weighted geometric mean，需要调整公式；如果当前公式符合预期，建议改注释和方法文档。

#### 5.2.5 multi-gene reaction 的权重可能被统一应用

当 gene 字段中含 `Gene2;Gene3` 且该行 weight=2.0 时，当前逻辑对 Gene2 和 Gene3 都用同一 weight。若未来需要每个基因不同权重，当前 schema 不支持。

建议：若数据库中 weight 是 reaction-level 权重，当前逻辑可以接受；若是 gene-level 权重，应扩展 schema。

#### 5.2.6 permutation p-value 可能按完整事件 key 过于严格

置换检验中 key 包含：

```python
sender, receiver, metabolite, sensor_gene, sensor_type
```

置换后只有完全相同 key 的事件才比较。由于置换改变 cell label 分配，这种比较更像“同一事件的 label-shuffle null”，不是按 sensor type 或 metabolite 生成整体 null distribution。

影响：统计解释需要谨慎。

建议：文档说明当前 permutation 的 null 定义；如要更稳健，可考虑按 sensor_type 或 metabolite-sensor 层级聚合 null 分布。

#### 5.2.7 `load_cell_mesh_database()` 默认加载测试数据库

默认路径为：

```python
enzyme_test.csv
interaction_test.csv
```

但包里同时有 `Enzyme_new.csv` / `Interaction1.0.csv` 等更完整文件。

影响：用户直接调用 `run_cell_mesh(adata)` 时使用的不是 notebook 里用的完整数据库，而是测试数据库。

建议：确认默认数据库策略。若 `enzyme_test.csv` / `interaction_test.csv` 仅用于测试，不建议作为生产默认。

#### 5.2.8 `DATA_DIR.mkdir(exist_ok=True)` 在 import 时创建目录

`config.py` import 时会执行：

```python
DATA_DIR.mkdir(exist_ok=True)
```

影响：包 import 有文件系统副作用。对于 site-packages 或只读环境可能不理想。

建议：包数据目录通常不应在 import 时创建。若需要写入缓存，应使用用户 cache 目录。

### 5.3 文档/API 一致性问题

#### 5.3.1 文档中存在旧算法痕迹

当前代码注释多次说明“旧算法已删除”，但仍保留：

- `role_scores`
- `scoring.py` 的 zscore/sigmoid 工具
- 测试中对旧字段 `sensor_expr_z` / `beta_sensor` / `beta_specificity` 的兼容分支

建议：统一文档：当前算法是 availability + robust min-max sensor scoring；旧算法相关内容标为 legacy 或移除。

#### 5.3.2 `load_cell_mesh_database()` docstring 中引用 `cell_mesh.run_cell_mesh`

当前包名为 `cellmesh`，docstring 写的是 `cell_mesh.run_cell_mesh`。

建议：修正为 `cellmesh.run_cell_mesh`。

---

## 6. 建议优先级

### 高优先级

1. 明确并修复 `eps` 参数是否应生效。
2. 统一默认数据库：确认是否应从 `enzyme_test.csv` / `interaction_test.csv` 改为完整数据库。
3. 避免稀疏矩阵全量 densify，至少在文档中说明当前内存限制。
4. 统一 `min_cells` 默认值，避免 `MIN_CELL_COUNT=100` 与 public API 默认 `1` 混淆。

### 中优先级

1. 复用 pseudobulk，减少重复计算。
2. 清理旧算法遗留：`role_scores`, `scoring.py`, 测试兼容分支。
3. 使用 `ROLE_TO_DIRECTION` 常量，移除重复映射。
4. 支持 interaction lower-case schema。

### 低优先级

1. 清理未使用 import。
2. 修正文档中的包名、legacy 描述。
3. 避免 import 时创建数据目录。

---

## 7. 本次未改动项

按照要求，本次未改动包源码逻辑文件：

- `cellmesh/config.py`
- `cellmesh/core.py`
- `cellmesh/database.py`
- `cellmesh/preprocess.py`
- `cellmesh/scoring.py`

本次只修改/新增：

- 修改 `examples/test_cellmesh_comprehensive.ipynb`，补充展示与断言；
- 新增本报告 `CELL_MESH_CODE_REVIEW_REPORT.md`。

注意：审查开始前工作区已有若干未提交修改和未跟踪文件；本次未对这些源码修改做回滚或覆盖。
