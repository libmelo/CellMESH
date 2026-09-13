# 数据库列名与输入规范

本规范适用于 `load_cell_mesh_database()` 和 `run_cell_mesh()` 接收的 CSV 或
DataFrame。项目根目录 [README](../../README.md#enzyme-table-columns) 同时列出列名映射。
下表列出程序实际识别的标准列名和兼容列名；每个必需字段提供其中一列即可。

测试阶段数据库内容不上传 GitHub。本目录仅将此说明文件纳入版本管理，CSV、
AnnData 文件及归档数据继续忽略。默认加载器只扫描本目录，分别选择数字版本号最高的
`Enzyme<版本>.csv` 和 `Interaction<版本>.csv`，两者版本不必相同，也不按修改时间排序。
没有相应版本文件时分别回退到 `enzyme_test.csv` 或 `interaction_test.csv`。
`archive/`、笔记本检查点、`Enzyme_new.csv` 等名称不参与自动版本选择；需要时显式传入路径。

## 通用文件规则

1. CSV 使用 UTF-8，可带 BOM；第一行是列名。带逗号、引号或换行的单元格应按 CSV
   规则引用，建议使用 pandas 的 `to_csv(index=False)` 保存。
   每条逻辑记录的字段数必须与表头一致，空值也须保留相应空字段，不能省略列。
   超宽行、短行、非法引号或 NUL 字符会报告文件和原始行位置；不会猜测隐式索引、
   截断多余字段或跳过坏记录。引用字段内的换行属于同一记录，记录外的空白行忽略。
   CSV 压缩输入沿用 pandas 的压缩识别，包含 gzip、bz2、xz、zip 等支持的格式。
2. 列名必须按表中拼写提供，区分大小写和空格。例如 `Gene_name` 是已知列名，
   `GENE_NAME`、` Gene_name ` 不会自动转换。
3. 必需字段应当是表格列，不能只放在 DataFrame 索引中。缺少必需列会报错。
   “必需列”与“每一行都有效”不同；空值的具体处理见下文。
   代谢物名称列必须提供，但名称值允许缺失或空白；名称不参与身份匹配、分组或筛选。
4. 完全相同的列名重复时直接报错，即使内容相同。CSV 在 pandas 自动将重名列改成
   `.1` 之前检查；DataFrame 也检查重复列名。
5. 不同列名可能代表同一字段，例如 `sensor_gene` 与 `Gene_name`。同义列允许并存，
   但必须逐行一致：等价值合并，空值由其他列补齐，非空冲突报错。检查发生在删行、
   去重和表达基因匹配之前。标准列存在不代表可以跳过其他列。
6. 字段无冲突时，优先保留有效的标准列值，再按表中顺序取有效兼容列值。
   合并后移除已使用的兼容列；其他元数据保留。Enzyme 基因证据的合并规则见下文。
7. HMDB ID 去除首尾空白并转成大写后匹配；代谢物名称用于展示，不作为两库关联键。
   两库可以为同一 HMDB 使用不同展示名称。建议填写正式 HMDB 编号；当前代码不会
   核验其是否实际存在于 HMDB 数据库。
   评分展示名从规范化先验中选取：优先 Enzyme 的首个非空名称，其次 Interaction，
   均缺失时使用 HMDB ID；独立评分入口使用其可用先验。事件身份为
   `hmdb_id + sensor_gene + sender + receiver`，样本记录另加 `sample`。
   `metabolite` 和 `sensor_type` 不参与事件身份或样本分组；类型仍用于分类型 FDR。
8. 已识别的文本字段按字符串读取，保留 `00123` 的前导零和字面值 `NA`。
   CSV 空单元格作为缺失值；DataFrame 的标识符也建议预先存为字符串。
   HMDB 字段中的空白及 `nan`、`none`、`null`（不区分大小写）视为缺失。
9. 两库中的旧版或自定义 `weight` 均被移除，不参与评分。来源、蛋白名称、文献和
   证据信息也不作为数值权重。

若表达矩阵从 10X 目录读取，默认使用基因符号；Enzyme 的 `gene` 和 Interaction 的
`sensor_gene` 均须与所选表达矩阵标识一致。选择 `var_names="gene_ids"` 时，两类
先验也须使用对应基因 ID。读取器不自动编号、合并重复基因或转换标识体系。
详细文件布局和参数见 [AnnData 读取规范](../../docs/ANNDATA_README.md#10x-模式)。

## Interaction：代谢物与传感器

| 字段含义 | 标准列名 | 兼容列名（按优先顺序） | 必需 |
|---|---|---|---|
| 代谢物名称 | `metabolite` | `standard_metName`、`standard_metname` | 是 |
| HMDB 编号 | `hmdb_id` | `HMDB_ID` | 是 |
| 传感器基因 | `sensor_gene` | `Gene_name`、`gene_name`、`gene` | 是 |
| 传感器类型 | `sensor_type` | `Annotation`、`annotation` | 是 |
| 数据来源 | `source` | `Database source`、`database_source` | 否 |
| 蛋白名称 | `protein_name` | `Protein_name` | 否 |
| 参考文献 | `reference` | `Reference` | 否 |
| 证据信息 | `evidence_level` | 无 | 否 |

当前默认 `Interaction1.41.csv` 使用 `Gene_name`，加载后改名为 `sensor_gene` 并用于
表达基因匹配。`ID`、`Interaction_mode`、`STITCH_evidence` 等额外列保留为元数据，
它们不决定传感器基因或类型。

### 值的规范与冲突判定

- 每行填写一个传感器基因，并与 `adata.var_names` 中的基因符号匹配。基因值去除
  首尾空白，但区分大小写；`S1` 和 `s1` 不视为同一基因。多个传感器关系应拆成多行。
  Interaction 不会按 `;`、`,`、`|` 拆分基因，也不会解析 Enzyme 式 `[evidence]` 后缀。
- 同义基因列都为空（包括缺失值、空白、字面值 `nan`）时，该行被删除。
  只有标准列为空时，会先从兼容列补齐。运行时再排除 HMDB 缺失或基因未被测得的行。
- 代谢物名称、来源、蛋白名称、文献的同义列按去除首尾空白后的文本比较，区分大小写。
  比较一致后保留优先列的原始展示文本；不会查询名称同义词库，也不会自动合并不同文献。
- 类型推荐直接填写下列三个标准值。类型文本去除首尾空白后按不区分大小写的规则归类：

| 类型文本 | 标准类型 |
|---|---|
| 含 `cell surface` 或 `surface receptor` | `Cell surface receptor` |
| 不符合上一项，但含 `transport` | `Transporter` |
| 其余文本，或全部类型列为空 | `Other receptor` |

同时存在多个类型列时，先补空值，再比较归类结果。明确写出的 `Other receptor` 是有效
类型，不能用另一列的 `Transporter` 覆盖。未知文本归为 `Other receptor` 是现有兼容行为，
建议明确填写标准类型，避免拼写错误影响分类型 FDR。

`Annotation` / `annotation` 的原文用于补齐空缺的 `evidence_level`；已有非空证据优先保留。
没有显式证据时，同类但不同的注释原文按上述别名顺序以 `; ` 连接；去除首尾空白后相同的
原文只保留一次。之后移除原始类型别名，标准类型由 `sensor_type` 单独承载。

同一标准 HMDB ID 与同一传感器基因的重复行，类型相同时保留首条元数据；类型不同时报错。
这一跨行检查与同一行的同义列检查都必须通过。

### 最小示例

标准格式：

```csv
metabolite,hmdb_id,sensor_gene,sensor_type
M,HMDB0000001,S1,Transporter
```

等效原始格式：

```csv
standard_metName,HMDB_ID,Gene_name,Annotation
M,HMDB0000001,S1,Transporter
```

同义列并存时：

| `sensor_gene` | `Gene_name` | 处理 |
|---|---|---|
| `S1` | `S1` | 合并为 `sensor_gene=S1` |
| 空 | `S1` | 补齐为 `sensor_gene=S1` |
| `S1` | `S2` | 报错，列出字段、列名、值和数据行位置 |

错误中的数据行位置从 1 开始，不包含表头；同时提供 DataFrame 索引标签。
带多行单元格的 CSV 数据行位置不是文件的物理行号。

## Enzyme：酶、反应与代谢物

| 字段含义 | 标准列名 | 兼容列名 | 必需 |
|---|---|---|---|
| 代谢物名称 | `metabolite` | `standard_metName` | 是 |
| HMDB 编号 | `hmdb_id` | `HMDB_ID` | 是 |
| 酶基因 | `gene` | `Gene_name` | 是 |
| 反应编号 | `reaction` | `Reactions` | 是 |
| 反应角色 | `role` | `Direction`、`direction`（值映射见下表） | 是 |
| 证据信息 | `evidence_level` | 无 | 否 |
| 数据来源 | `source` | 无 | 否 |

Enzyme 只使用上表的列名映射。例如，Interaction 接受的 `gene_name` 不属于 Enzyme 的
兼容基因列名。文献及其他元数据原样保留，Enzyme 不另将 `Reference` 等元数据列重命名。

### 角色与基因规范

| `role` 标准值 | `Direction` / `direction` 原始值 | 容量矩阵 |
|---|---|---|
| `production` | `product` | P：生成 |
| `degradation` | `substrate` | C：消耗 |
| `export` | `exporter`、`export`、旧版 `transporter` | E：外排 |

- 角色和方向值去除首尾空白并转成小写后处理。`role` 应使用标准值，方向列应使用
  原始值；同时提供时按角色含义检查一致性。空值可补齐，冲突报错。
- 填写非空反应编号；缺失或仅含空白的反应编号会报错。反应编号按文本保留前导零。
  空基因条目被忽略，不支持的角色被过滤；应按上表提供有效值。
- Enzyme 的 `gene` / `Gene_name` 可以包含多个基因，以 `;`、`,`、`|` 分隔；
  允许 `G1[reviewed];G2[curated]`。方括号中的分隔符属于证据原文。
  加载时先拆成单基因行，再进行表达基因匹配，保留原反应编号和元数据。
- 多个同义基因列按解析后的基因集合比较，忽略排列和分隔符差异。基因符号仍区分大小写。
  同义列中的基因证据合并保留；已提供 `evidence_level` 列时沿用该列。
- 反应以标准 HMDB ID、反应编号和方向共同分组。不同反应编号不要误用为同一个反应。
  同一代谢物和方向下按完整基因集合去重、剔除严格子集，再用测得的基因计算活性。

### 生成能力 P：零值与缺失

P 表示生成能力；置换检验的统计学 p 值另称 `perm_pvalue`。
在完整反应集合去重后，以是否有生成基因存在于规范化的 `adata.var_names` 中判断可评估性，
不能用 `P > 0` 判断。只有部分基因可用时仍按现有规则用可用基因计算，不把缺失基因填零
加入几何均值。

| 情况 | `production_status` | `production_evaluable` | 处理 |
|---|---|---|---|
| 无生成先验 | `prior_missing` | False | 不生成发送方得分或事件 |
| 有生成先验，但生成基因全部不可用 | `prior_gene_unavailable` | False | 不生成发送方得分或事件 |
| 有可用生成基因，计算得到的 P 全零 | `prior_no_expression` | True | 保留零分，计入样本汇总 |
| 有可用生成基因，至少一种细胞类型的 P 为正 | `supported` | True | 按各细胞类型的实际结果评分，保留其中的零分 |

观测与置换均先校验原始表达及 P/C/E 容量，再按 `production_evaluable` 筛选参与
归一化的代谢物。无生成先验或生成基因全部不可用的行不计算参考值，避免在提取观测事件
之前进行无用的 C/E 归一化。可计算的零值继续保留，参考值仍使用该评分单元的全部
观测细胞类型；事件集合不会因置换扩展。

已评分代谢物的状态位于 `metadata` 中，该表与评分矩阵保持索引一致。独立的
`production_diagnostics` 表含生成反应数和上述两个状态列，覆盖所有 Enzyme 先验代谢物，
包括无法评分的代谢物。状态按整个评分单元汇总：pooled 模式是全部细胞，sample_aware 模式是单个样本。
“supported”不表示其中每个发送细胞类型都为正。

无法评估的反应内部可能用数值 0 占位，但不能把它当作已测得的零表达。展示完整 P 表时，
可用 `unit["P"].reindex(unit["production_diagnostics"].index)` 补出被排除的行；这些行显示为 NA。
如果整份 Enzyme 先验没有任何基因匹配表达矩阵，主流程仍直接报错。

样本中发送和接收细胞类型都存在且可以评分时，生成能力全零对应事件分数 0，计入中位数、
可评估样本数和阳性比例的分母。缺少某端细胞类型的样本保留 NA，不把 NA 填成零。
观测与置换使用同一规则。两种分析模式均保留可评估的全零代谢物；执行置换时零分事件
的 `perm_pvalue=1`，未执行置换时仍为 NA。保留这些事件可能改变 FDR 校正范围。

基因可用性依据输入矩阵的基因列。如果上游已把某些未测量数据填成 0，仅凭矩阵无法恢复
缺失原因，需要额外保留测量可用性信息；程序不会根据零表达推断基因未测量。

### 最小示例

标准格式：

```csv
metabolite,hmdb_id,gene,reaction,role
M,HMDB0000001,G1;G2,001,production
```

等效原始格式：

```csv
standard_metName,HMDB_ID,Gene_name,Reactions,Direction
M,HMDB0000001,G1;G2,001,product
```

## 维护要求

列名映射以 `cellmesh/database.py` 中的 `_ENZYME_ALIASES`、`_SENSOR_ALIASES` 及角色映射
为准。修改识别规则时，应同时更新本文件、项目根目录 README 和对应回归测试。
不能为提速跳过同义列检查，或先丢弃空标准列的行后才查兼容列；这会再次丢失有效记录。
Interaction 回归测试位于 `tests/test_interaction_aliases.py`，Enzyme 回归测试位于
`tests/test_enzyme_aliases.py` 和 `tests/test_enzyme_multigene.py`。
