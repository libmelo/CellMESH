# CellMesh

Metabolite-mediated Event Scoring with Sensor Hierarchies

**CELL MESH** = **Metabolite-mediated Event Scoring with Sensor Hierarchies**.

CELL MESH is a Python package for inferring metabolite-mediated cell-cell communication events from single-cell expression data, **fully based on the metabolite availability algorithm**.

## Algorithm Overview

### 1. Enzyme Prior Normalization
The enzyme-metabolite prior table is the canonical enzyme input. `run_cell_mesh()`
validates this prior once with `validate_priors()` and passes it directly into
the availability calculation. Internally, enzyme roles are normalized to the
availability directions used to build the P/C/E matrices:

| Role in prior | Direction | Matrix | Meaning |
|---|---|---|---|
| `production` | `product` | P | Metabolite production ability |
| `degradation` | `substrate` | C | Metabolite consumption-enzyme ability proxy |
| `export` | `exporter` | E | Metabolite efflux/transport ability |

### 2. Sender-abundance-aware P/C/E
CELL MESH first computes mean expression and reaction activity within each cell
type. It then adjusts each per-cell reaction activity by sender abundance before
constructing P/C/E:

```text
cell_fraction(t) = n_cells(t) / n_cells(all cell types)
sender_abundance_weight(t) = cell_fraction(t) ** sender_abundance_exponent
adjusted_reaction_activity(t) = reaction_activity(t) * sender_abundance_weight(t)
```

`sender_abundance_exponent` defaults to `1.0`, preserving the current linear
population-capacity model. Setting it to `0.0` removes sender abundance from
the score, values between 0 and 1 provide a sublinear adjustment, and values
above 1 strengthen the abundance effect.

In `sample_aware` mode, fractions are calculated independently within each
sample. Every observed cell type contributes to the fraction denominator and
to formal pseudobulk scoring. `min_cells` is a QC threshold only: units below
it remain in P/C/E construction, references, scores, events, and permutations,
but are marked as not passing min-cells QC. In pooled mode, fractions use all
cells in the full AnnData object and the same QC-only rule applies.

Production, consumption, and exporter reaction capacities are summed into the
P/C/E matrices after reaction-level maximal-gene-set filtering. Within one HMDB
ID and one direction, identical complete gene sets contribute once, and a
reaction whose gene set is a strict subset of another reaction is omitted.
Incomparable sets, including sets that overlap without either containing the
other, are both retained intact. For each metabolite and each direction
\(X \in \{P,C,E\}\),
CELL MESH calculates a reference from the strictly positive capacities across
observed cell types and applies a continuous saturation transform:

```text
X_ref(m) = mean({X(m,t) | X(m,t) > 0})             default
X_ref(m) = median({X(m,t) | X(m,t) > 0})           pce_reference="median"
X_score(m,t) = 0                                  if X(m,t) = 0
X_score(m,t) = X(m,t) / (X(m,t) + X_ref(m))       if X(m,t) > 0
```

The reference is not cell-count weighted because sender abundance has already
entered each capacity. The default positive arithmetic mean lets a cell type
with a genuinely large abundance-adjusted population contribution influence
the reference. The positive median remains available as a more outlier-robust,
equal-cell-type alternative. Restricting either statistic to positive values
prevents a zero-inflated cell-type panel from collapsing the reference to zero.
A capacity equal to the selected reference scores `0.5`; smaller positive
values remain between `0` and `0.5`, and larger values remain between `0.5`
and `1`.

### 3. Sender Score

The base sender score is:

```text
base_sender_score(m,t) = P_score(m,t)^2 / (P_score(m,t) + C_score(m,t))
```

The result is defined as zero when the denominator is zero. Production is the
required anchor; increasing `C_score` smoothly reduces the score without a hard
reference cutoff. If the consumption prior is missing, or its genes have no
observed expression, `C_score = 0` and the sender score reduces to `P_score`.
`consumption_status` distinguishes `prior_missing` (no database relation),
`prior_gene_unavailable` (relation present but none of its genes were measured),
`prior_no_expression` (measured genes present but all capacities are zero), and
`supported`. The first three all give `C_score = 0`, but retain distinct
evidence semantics.

Production is evaluable when at least one gene from its retained reaction
definitions exists in the expression matrix. A computed `P=0` remains a zero
sender/event score, including when every cell type in a sample has zero
production. Missing production genes do not supply measured zero evidence;
those metabolites are excluded from scoring and identified in metadata.

Exporter evidence supplies a bounded, non-increasing modulation:

```text
E_factor(m,t) = (1 - export_weight) + export_weight * E_effective(m,t)
sender_score(m,t) = base_sender_score(m,t) * E_factor(m,t)
```

For a supported exporter prior, `E_effective = E_score`. When the exporter
prior is missing or none of its genes were measured, `E_effective = 0.5` is a
fixed neutral-evidence convention. When measured exporter genes are present
but all capacities are zero, `E_effective = 0`. With the default
`export_weight=0.2`, `E_factor` is restricted to `[0.8, 1.0]`; setting the
weight to zero exactly recovers the base score. `export_weight=0.2` is the
current configurable model setting, not a value automatically derived from
database coverage. Database coverage is not inserted into individual scores.

For the local **Enzyme1.39.csv / Interaction1.41.csv** snapshot checked on
2026-09-16, exporter-prior coverage is **76 / 1095 = 6.9406%** of all unique
Enzyme HMDB IDs, or **56 / 346 = 16.1850%** when the denominator is restricted
to HMDB IDs shared by both tables. These are distinct denominators, calculated
from normalized complete priors without expression filtering. They do not
measure how many exporters are evaluable in a particular expression dataset.
File hashes and counts are recorded in the
[coverage snapshot](docs/PRIOR_COVERAGE_2026-09-16.json). This replaces the
unsupported 17.1% statement; neither the weight nor the formula changes.

Recompute coverage for the current defaults or explicit files:

```bash
python -m cellmesh.prior_coverage
python -m cellmesh.prior_coverage --enzyme cellmesh/data/Enzyme1.39.csv --interaction cellmesh/data/Interaction1.41.csv
```

The JSON report includes both numerators/denominators, file identities and
SHA-256 checksums. A zero denominator has `percent: null`, not a fabricated
zero coverage. Compare hashes before comparing a later run with the snapshot;
a matching filename alone does not establish identical contents.

Raw C is an expression-derived proxy for metabolite-consuming enzyme capacity;
it is not a direct measurement of extracellular clearance flux.

### 4. Sensor Score Calculation

For each sensor gene, CELL MESH first calculates its per-cell-type pseudobulk
mean expression. The reference is calculated across strictly positive observed
cell types, with each cell type receiving equal weight:

```text
R_ref(g) = median({R(g,t) | R(g,t) > 0})           default
R_ref(g) = mean({R(g,t) | R(g,t) > 0})             receiver_reference="mean"
R_score(g,t) = 0                                   if R(g,t) = 0
R_score(g,t) = R(g,t) / (R(g,t) + R_ref(g))        if R(g,t) > 0
```

A positive expression value equal to the selected reference scores `0.5`;
smaller positive values remain between `0` and `0.5`, and larger values remain
between `0.5` and `1`. In `sample_aware` mode, the reference is calculated
independently within each sample. When `min_expr_frac` is not `None`, receiver
units below that optional gate are assigned score zero.

Receiver cell abundance does not enter `sensor_score`. `receiver_cell_fraction`
is retained only as descriptive metadata and does not weight the reference.

### 5. Communication Score
Communication score is the geometric mean of metabolite availability and sensor score:
```
sample_level_cell_mesh_score = sqrt(metabolite_availability * sensor_score)
```
Events are joined by a stripped, case-normalized `hmdb_id`. Metabolite names are
display metadata and may differ between enzyme and sensor priors; the sender-side
name is used in the event table. Prior rows without `hmdb_id` are excluded before
event construction.

### 6. Sensor Types
Sensor types are read from `sensor_type` or its `Annotation`/`annotation` aliases
and normalized to `Cell surface receptor`, `Transporter`, or `Other receptor`.
Raw annotations provide provenance through `evidence_level`; explicitly supplied
evidence takes precedence. The current default database is `Interaction1.41.csv`.

Both inference modes report global FDR and sensor-type-specific FDR.

## Included Data
The package includes built-in prior databases and a small AnnData example file:

- `cellmesh/data/Enzyme1.39.csv`: current default enzyme/reaction/metabolite table
- `cellmesh/data/Interaction1.41.csv`: current default metabolite-sensor table
- `cellmesh/data/enzyme_test.csv`: legacy small enzyme prior used for compatibility tests
- `cellmesh/data/interaction_test.csv`: legacy small sensor prior used for compatibility tests
- `cellmesh/data/Enzyme_new.csv`: walkthrough/test enzyme prior
- `cellmesh/data/test_single_cell.h5ad`: walkthrough/test single-cell data

`load_cell_mesh_database()` automatically loads the highest packaged enzyme file
named `Enzyme<version>.csv` and the highest packaged interaction file named
`Interaction<version>.csv` if file paths are not explicitly provided. The two
version numbers do not need to match: if the highest enzyme prior is
`Enzyme1.2.csv` and the highest interaction prior is `Interaction2.2.csv`, those
two files are selected together. With the current packaged files, the default is
`Enzyme1.39.csv` and `Interaction1.41.csv`. Historical erroneous databases
`Enzyme2.0.csv` and `Interaction4.0.csv` are retained in `cellmesh/data/archive/`;
the loader scans only the data directory itself and does not select archived files.
If no versioned files are available,
the loader falls back to the legacy `enzyme_test.csv` / `interaction_test.csv`
pair. The comprehensive walkthrough notebook uses `Enzyme_new.csv`,
`Interaction1.0.csv`, and `test_single_cell.h5ad` so its calculations are fully
reproducible from packaged files.

Recognized database column names are listed in [the input tables below](#enzyme-table-columns).
The database directory also contains a self-contained Chinese
[数据库列名与输入规范](cellmesh/data/README.md), with both mappings, value rules,
and CSV examples. During testing, database contents remain excluded from Git;
the directory's README is versioned separately.

## Install

From a local checkout:

```bash
pip install -e .
```

For running the walkthrough notebook:

```bash
pip install -e ".[notebook]"
```

The four [example notebooks](examples/) should be run from a fresh kernel using
**Restart Kernel and Run All**. Saved outputs are cleared so they cannot be
mistaken for results from a different code/database version. The HNSC examples
use measured genes to select demo HMDB IDs, then take **all enzyme records** for
those IDs from the complete prior. Do not prefilter enzyme genes before passing
these records to the algorithm: this can change reaction subset relationships.
The comprehensive trace uses `validate_priors(..., filter_enzyme_genes=False)`
to match the main workflow. An empty significance selection is valid: the pooled
HNSC network reports it without relaxing thresholds. Later exploratory plots
explicitly use their own selection and do not imply statistical significance.

The toy example verifies production diagnostics and bounded exporter modulation.
Its inputs are generated in the notebook. The other examples need the local
HNSC/default prior or fixed walkthrough files described above; database contents
remain excluded from Git during testing. To check the example input contracts,
empty/nonempty network paths, and run all four notebooks in fresh Python
processes (including figure rendering), use:

```bash
python -m pytest -q tests/test_notebook_contracts.py
```

This checks sequential Python execution, not the Jupyter browser interface.
Notebook tests explicitly skip when IPython or their required local fixtures are absent.

For development and tests:

```bash
pip install -e ".[dev]"
pytest -q
```

10X directory loading uses Scanpy and can be enabled with:

```bash
pip install -e ".[scanpy]"
```

Loom loading also requires `loompy`; enable it with `pip install -e ".[loom]"`
(or `pip install 'cellmesh[loom]'` for an installed package). Missing dependencies
produce an installation hint; invalid-file errors remain distinguishable.

### Prior input provenance

Row-level `source` is evidence supplied by the data provider. Existing values
are retained; missing Enzyme sources are not filled with a guessed database
name based on the use of `direction` columns.

`load_cell_mesh_database()` separately records each input in
`table.attrs["input_provenance"]`; `run_cell_mesh()` stores both records under
`result.parameters["prior_inputs"]`, also exported in `.parameters.json`:

- `input_kind`: `default_file`, `user_file` or `dataframe`.
- File inputs: `filename`, resolved `path`, `filename_version` when the basename
  follows the Enzyme/Interaction version pattern, and `sha256` of file bytes
  (compressed bytes for compressed CSVs). A filename version is a naming claim,
  not independent verification of a database release. Content changes during
  loading are rejected.
- `raw_rows` and `normalized_rows`: counts before/after normalization and gene
  expansion, before expression-gene filtering.
- DataFrame inputs have null file/path/version/hash fields. Previously loaded
  DataFrames may have been edited, so a new load records `dataframe` instead of
  copying stale file claims from incoming attrs. User data and row-level source
  annotations are not mutated. Retain earlier provenance separately if needed.

These records describe the current input and do not certify its biological
accuracy. They do not enter scoring or change `export_weight`.

## Basic Usage

For expression-file import, see the [AnnData reading guide](docs/ANNDATA_README.md).
H5AD inputs opened with `backed="r"` are supported by both inference modes,
standalone scoring and violin plots, including dense/CSR/CSC storage and
reordered backed views. Shared slicing reads sorted positions and restores the
requested order; sample-aware materializes the current sample's chosen matrix.
The caller retains ownership of the open file and should close `adata.file`
after use. This is not a fully out-of-core algorithm: group/sample matrices,
permutation-scoring gene subsets and intermediate results still require memory;
AnnData may load layers into memory at read time. See the reading guide for the
memory boundary and file-lifetime example.

For CSR/CSC inputs with duplicate stored `(cell, gene)` coordinates, calculation
uses a sparse copy promoted to float64 before combining the entries. This avoids
overflow in the original dtype during expression-fraction comparisons or violin
extraction. Raw relevant values are checked before combination, so negative
entries cannot be hidden by cancellation. Caller matrices and files remain
unchanged; ordinary canonical matrices do not require this copy.

For CSV/TSV inputs, `read_anndata()` checks original headers before pandas
can rename duplicate genes, and preserves axis identifiers such as `01`, `1`,
and literal `NA` as distinct text. Expression values remain numeric. Cell and
gene metadata are read as text by default, including custom sample/cell-type
columns; empty fields remain missing. Convert additional numeric metadata
explicitly when needed. MTX gene/barcode lists preserve the same literal names.
Expression text, metadata and MTX name files are scanned for actual NUL
characters before parsing, using the selected encoding and compression. A NUL
raises an error with the file and physical line, even in a record excluded by
`usecols`, `skiprows` or `nrows`; it is never silently removed or truncated.
Pandas nullable numeric options such as `dtype="Float64"`, `dtype="Int64"` and
`dtype_backend="numpy_nullable"` are converted to NumPy numeric columns before
AnnData construction and transposition. Missing values remain NaN and undergo
the usual scoring validation. Homogeneous integer columns without missing
values retain their integer dtype.
The MTX reader supports both sparse `coordinate` and dense `array` storage,
including `.mtx.gz` files. Both preserve expression values and transpose from
gene-by-cell to cell-by-gene; sparse inputs stay sparse and dense inputs stay
NumPy arrays.
10X imports also preserve original gene/barcode text and reject empty or
duplicate selected identifiers before constructing the labeled matrix. Both
legacy `genes.tsv` and modern `features.tsv` layouts support plain or gzip
files and optional prefixes. `var_names="gene_symbols"` is the default;
`"gene_ids"` requires priors using the same gene IDs. Automatic gene renaming
is disabled, and explicit `make_unique=True` is rejected. See the reading guide
for feature filtering, file ambiguity, and cache behavior.
Duplicate or empty axis identifiers are rejected, including duplicates after
stripping whitespace. Results affected by older automatic renaming or numeric
label inference require rereading the original files and rerunning analysis.

Both `enzyme_metabolite` and `metabolite_sensor` default to `None`.
`run_cell_mesh(adata)` calls `load_cell_mesh_database()` to load and normalize
the highest version of each database in `cellmesh/data` automatically. Versions
are compared numerically, not by file modification time. If you supply only one
prior table, the other is filled from the packaged database.

Each prior also accepts a CSV path (`str` or `pathlib.Path`). CSV files are read
into DataFrames, then both input forms use the same schema normalization and
sensor conflict checks. Raw database fields and normalized prior fields are
accepted. Input DataFrames are copied, and existing provenance metadata is
preserved. Paths, DataFrames, and `None` can be mixed independently:

```python
res = run_cell_mesh(
    adata,
    enzyme_metabolite="/path/to/Enzyme.csv",
    metabolite_sensor="/path/to/Interaction.csv",
    cell_type_key="cell_type",
)
```

```python
from cellmesh import run_cell_mesh

res = run_cell_mesh(
    adata,
    cell_type_key="cell_type",
    sample_key="sample",      # optional, for within-sample permutation
    sample_mode="pooled_stratified",
    layer="lognorm",          # optional, use specific expression layer
    n_perms=1000,             # optional, number of permutations for p-value calculation
    n_jobs=1,                 # use -1 for all available CPUs
    store_null_scores=False,  # online counts avoid event × permutation storage
    min_cells=100,            # QC-only threshold; does not filter calculations
    min_expr_frac=0.10,       # optional, minimum expression fraction for sensor genes
    pce_reference="mean",     # "mean" (default) or "median"
    receiver_reference="median",  # "median" (default) or "mean"
    allow_self=True,          # optional, allow self-communication events
)

# View results
res.events.head()  # All communication events
res.sender_scores  # Metabolite × cell type availability matrix
res.receiver_scores  # Receiver/sensor scores
res.availability_results  # All intermediate calculation results
```

Expression values used by enzyme or sensor scoring must be real, finite, and
non-negative. CELL MESH checks these cell-level values in the selected `layer`
(or `adata.X`) before aggregation, including when `n_perms=0`. Standalone
availability/sensor scoring, compiled permutations, and expression plots use
the same rule. Genes outside the scoring priors and unselected layers do not
enter this check. Values are not clipped or replaced. Sparse validation retains
sparse storage and checks stored values in bounded blocks.

Finite input can still overflow during aggregation. Both observed and permuted
scoring check relevant pseudobulks, reaction activities/capacities, normalization
references and denominators before invalid intermediates can become zero scores.
Each permutation's event scores are checked again before p-value counting.
Numerical failures stop the analysis with an error naming the stage; they are
not replaced by zero or treated as structural sample missingness. Normal
structural missing values and defined zero scores retain their existing meaning.

Earlier observation-only runs could accept negative values when averaging or
reaction aggregation hid them. For example, `[-1, 3]` has a positive mean even
though it contains an invalid negative expression value. Analyses affected by
this case should be rerun with a valid non-negative expression layer.

Cell-type labels, sample labels (when `sample_key` is supplied), and
`adata.var_names` are converted to text and stripped of leading/trailing
whitespace for all internal matching, scoring, permutation, QC, and expression
plotting. The input AnnData object is not modified. Missing or empty identifiers
are rejected before scoring. Gene names must remain unique after normalization;
distinct observed cell-type or sample labels that become identical (for example,
`"A"` and `" A "`, or numeric `1` and text `"1"`) also raise an error instead of
being silently merged. Repeated cells with the same label are valid; unused
categorical levels do not create collisions. Explicit `cell_fractions` indices
and expression-plot cell-type selections use the same whitespace rule.

Explicit `cell_fractions` for standalone scoring accept real numbers and numeric
text, including scientific notation. Booleans, complex values, dates, missing
values and infinities are rejected before any lossy conversion. Fractions must
be strictly positive, sum to at most 1 (tolerance `1e-9`), have unique normalized
indices and cover all required observed cell types. Valid custom fractions are
preserved rather than replaced with proportions derived from cell counts.

Earlier versions could compare stripped observed labels with unstripped
permutation labels, producing false significance. Rerun analyses affected by
leading/trailing label whitespace to regenerate p-values and FDR.

## Visualization

Visualization functions are available from the top-level package and from the
`cellmesh.plotting` implementation module.

All six plotting functions validate numeric inputs before the operations that
use them. Event tables are checked after min-cells QC and explicit context
selection, but **before numeric thresholds, ranking and duplicate resolution**.
The count plot checks the probabilities used by its thresholds; network and
dot plots also check present probability columns used for ranking or display.
The supplied data are not modified.

| Plot input | Accepted values |
|---|---|
| Native `cell_mesh_score`, sender availability, receiver `sensor_score` | Finite real values in `[0, 1]`, or genuine missing values |
| p/FDR and `sensor_expr_frac` | Finite real values in `[0, 1]`, or genuine missing values |
| Custom event/sample `score_col` | Finite, non-negative real values; values above 1 are allowed |
| Single-cell expression, reaction activity, P/C/E capacities | Finite, non-negative real values; no `[0, 1]` upper bound |
| Thresholds | Finite values in the corresponding score/probability range |
| Dot/node size and edge-width output ranges | Positive finite minimum and maximum, in order; equal values give a fixed size |
| Color/value ranges | Finite bounds in order; equal bounds are allowed |

Numeric strings such as `"1e-2"` are converted to numbers before comparison, so
`0.01` ranks ahead of `0.5`. Malformed text (including literal `"NaN"` or empty
strings), infinity, complex values and booleans raise an error naming the
column and invalid record positions/indices. Genuine missing values are
`None`, `numpy.nan` or `pandas.NA`, rather than arbitrary text placeholders.
Valid native and nullable numeric dtypes are preserved, except that float16
working copies are promoted to float64 for pandas sorting compatibility.

Missingness remains distinct from zero:

- Count, network and dot plots omit missing-score records and return them in
  `missing_score_events`, with `plot_exclusion_reason="missing_score"`.
- Dot plots retain missing significance as an unavailable-significance marker;
  missing FDR is never filled with an unadjusted p-value.
- Sample plots retain the existing NA markers and exclude only missing scores
  from the median. Computed zeros remain in the median. Min-cells QC still
  applies; use `qc_only=False` to include the sample rows that fail it.
- Violin plots with automatic cell-type selection list omitted missing-score
  types in `missing_score_cell_types`. Explicitly requesting an unavailable
  type raises an error. `plot_data` contains plotted cell values, and
  `missing_value_data` records any missing cell values. If no usable scores
  remain, network/dot/sample/violin plots raise their no-data error; the count
  plot can return a zero-count matrix.

Dot-plot significance sizes distinguish zero, positive and unavailable values:

- Positive probabilities use `-log10(p/FDR)`. Without zeros, their existing area
  mapping is unchanged. If zeros are present in the selected plot, positives
  use the lower 80% of the requested area interval, and zeros use its maximum.
  This reserved area is a display convention, not a statistical distance.
- With default sizes, `[0, 0.01, 0.1]` produces areas `[260, 212, 20]`;
  `[0, 1, 1]` produces `[260, 116, 116]`. All-zero probabilities use area 260
  and a single `0 (display cap)` legend entry, without inventing `1e-12`.
- Both FDR and p-value legends label actual probabilities selected from the
  displayed data (up to three distinct positive values plus zero/NA keys).
  The p-value legend now reads `p-value` and labels raw probabilities, rather
  than labeling transformed `-log10` values. Sizes still use the log scale.
- NA retains its fixed-area diamond and `Unavailable` key. An existing NA FDR
  never falls back to a p-value. Equal requested size bounds retain fixed sizes
  and explicitly label that size encoding is disabled.
- Returned `dot_sizes` has the same index as `plot_events`. `size_encoding`
  records the statistic column, fixed-size setting, zero display cap, positive
  area range, and actual positive probabilities used in the legend. Input
  probabilities are not changed.

Dot plots measure the actual legend, marker, title, colorbar and axis-label
sizes before allocating space. The significance legend and score colorbar have
separate regions for every missingness state, including all available values,
zero-value keys and fixed-size explanations. Automatically created figures grow
to fit long labels or large legend markers. Returned `legend`, `legend_ax`,
`colorbar`, `title_artist` and `layout` expose the resulting layout;
`layout` includes required and available sizes in inches.

For a supplied `ax`, the plot, its labels, local title and sidebar all fit inside
that axes' existing rectangle. The function does not resize the figure,
reposition other subplots or replace a caller's figure-wide title. Allocate
sufficient subplot space; an undersized rectangle raises an error stating the
required and available dimensions. If the figure uses automatic `tight_layout`
or `constrained_layout`, finalize that layout **before** calling the dot plot:

```python
fig.canvas.draw()
fig.set_layout_engine("none")
# Now call plot_event_dotplot(..., ax=ax) in the finalized subplot.
```

An active automatic layout engine is rejected before modifying the supplied
axes, since it would otherwise rearrange the measured plot and other subplots
on redraw or export. The function never disables the caller's engine silently.
Finalize figure dimensions before plotting; do not subsequently shrink the
figure or run global layout commands on the finished plot. Recreate the plot
if those settings need to change. Export at a different DPI is supported;
PNG, SVG and PDF export are covered by regression checks.

Cell-type selectors and sample-order labels use the same whitespace stripping
as source labels, and duplicate normalized selector entries are rejected.
HMDB ID and receptor gene remain the matching identifiers; metabolite names
only affect display. Reaction-gene lists are stripped and deduplicated **once
for both identity comparison and actual expression calculation**. Complete
reaction definitions retain unmeasured genes for identity checks; expression
extraction still uses measured genes only.

Actual floating-point underflow retains the
[report-and-continue policy](#numerical-precision-and-underflow): its finite zero
is valid input to these plots and is not rejected by the new checks.

```python
from cellmesh import (
    plot_communication_network,
    plot_event_dotplot,
    plot_metabolite_secretion_violin,
    plot_receptor_expression_violin,
    plot_sample_event_scores,
    plot_significant_event_counts,
)

# Sender x receiver heatmap of unique events passing the selected thresholds.
overview = plot_significant_event_counts(
    res,
    max_fdr=0.05,
    fdr_col="fdr_global",  # or use "fdr_sensor_type"
    qc_only=True,          # default: display only events passing min-cells QC
)

# Directed circular network of significant cell-cell communication. Edge width
# is the number of unique communication events, edge color is their summed
# Cellmesh_score, and node area is the number of connected cell types.
network = plot_communication_network(
    res,
    max_fdr=0.05,
    fdr_col="fdr_global",
    qc_only=True,
)

# Alternative encodings can emphasize communication burden at each cell type
# or mean event strength for each directed cell pair.
network_by_burden = plot_communication_network(
    res,
    max_fdr=0.05,
    fdr_col="fdr_global",
    qc_only=True,
    node_size_by="event_count",
    edge_color_by="mean_score",
)

# Sender-faceted bubble plot. Color represents score and bubble area follows
# -log10(FDR), while the size legend displays the corresponding raw FDR values.
dotplot = plot_event_dotplot(
    res,
    top_n=10,
    fdr_col="fdr_global",
    qc_only=True,
)

# For sample-aware results, the default shows sample rows where both endpoints
# pass min-cells QC. Set qc_only=False to inspect all calculated rows. In that
# view, NA means the sample-level event was not computable (for example, an
# endpoint was unobserved). Measured zero production remains score 0; a low
# cell count alone never creates NA.
event = res.events.sort_values("fdr_global").iloc[0]
sample_plot = plot_sample_event_scores(
    res,
    sender=event["sender"],
    receiver=event["receiver"],
    metabolite=event["metabolite"],
    hmdb_id=event["hmdb_id"],
    sensor_gene=event["sensor_gene"],
    sensor_type=event["sensor_type"],
    qc_only=True,
)

# MEBOCOST-style single-cell violins. The metabolite violin shape shows the
# distribution of production-reaction enzyme activity in sender cells; its
# fill color is the formal CELL MESH metabolite availability score based on
# normalized production, consumption, and bounded exporter modulation.
secretion_plot = plot_metabolite_secretion_violin(
    res,
    adata,
    metabolite=event["metabolite"],
    hmdb_id=event["hmdb_id"],
)

# Receptor expression is taken directly from the same AnnData layer used by
# run_cell_mesh(), and violin fill color represents mean expression per
# receiver cell type.
receptor_plot = plot_receptor_expression_violin(
    res,
    adata,
    receptor_gene=event["sensor_gene"],
    metabolite=event["metabolite"],
    hmdb_id=event["hmdb_id"],
)
```

Each plotting function returns a dictionary containing the Matplotlib figure,
axes, and the exact filtered data used to construct the plot. The two violin
functions require the original `adata` because `CellMeshResult` stores
cell-type summaries rather than duplicating the single-cell expression matrix.

All six plotting functions treat metabolite names as **display metadata only**.
Names never participate in matching, grouping, deduplication, or ranking keys.
The identity rules are:

| Plot | Identity and record context |
|---|---|
| Event counts, communication network, event dotplot | `(hmdb_id, sensor_gene)` identifies a relation; `sender` and `receiver` distinguish directed records. |
| Sample event scores | The same relation and direction, with `sample` distinguishing records. |
| Metabolite secretion violin | Explicit `hmdb_id` only; production capacity does not depend on a receptor. |
| Receptor expression violin | Explicit `hmdb_id` plus `receptor_gene` (the result's `sensor_gene`). |

HMDB IDs use the same trimming and uppercase normalization as scoring; sensor
genes are trimmed with case preserved. A missing or unknown target ID raises
an error. Both violins and the sample plot require the explicit ID even if a
name is unique. Their optional `metabolite` argument only overrides the display
name. Thus an event's Enzyme name can differ from the Interaction name without
changing which receiver context is selected. Original result tables remain
unchanged.

Metabolite labels include `Name (HMDB ID)`, and relation labels add `-> GENE`.
If no name is available, the ID is displayed alone. Event labels use the first
non-empty name for each ID in the input table. Same-ID aliases represent one
identity; different IDs or sensor genes stay separate. Identical alias rows
count once. Conflicting values for one HMDB row, receiver context, or sample
record raise an error instead of using a name to choose a value. Counts,
network, and dotplot retain their existing threshold and duplicate-priority
rules. For counts and network, `unique_keys` must include `sender`, `receiver`,
`hmdb_id`, and `sensor_gene`; only `sample` may be added to count each sample
separately.

For the dotplot, `event_keys` accepts `(hmdb_id, sensor_gene)` tuples or
dictionaries containing those two fields, for example:

```python
dotplot = plot_event_dotplot(res, event_keys=[("HMDB0001509", "LTB4R2")])
```

Legacy `(name, hmdb_id, sensor_gene)` tuples remain accepted, but their name is
ignored. Display-label strings are no longer accepted as selectors.
`selected_events` returns display labels; `selected_event_keys` returns the
two-field tuples to reuse for selection. Requested order is retained and
repeated identities appear once. Update old name-only calls to pass IDs, then
replot affected results; event scores, permutation p-values, and FDR require
no recalculation.

Dotplot color always represents the selected score. Significance controls
marker size only when its value is available:

| Significance in the displayed rows | Markers and legend |
|---|---|
| All available | Circles sized by significance, with the numeric legend. |
| Partly missing | Available values size the circles; missing values use fixed-size diamonds labeled `Unavailable`. |
| All missing, including `n_perms=0`, or both probability columns absent | Fixed-size diamonds and `Significance: Unavailable`; no numeric significance legend. |

The fixed size is the midpoint of `min_dot_size` and `max_dot_size`. Missing
values stay NA in the returned data, and measured zero scores remain visible.
The FDR column, when present, supplies the size scale; missing FDR never falls
back to an unadjusted p-value. Only an absent FDR column enables the existing
p-value size scale. Non-missing displayed probabilities must be finite and in
`[0, 1]`; both 0 and 1 are available values. Explicit p-value/FDR thresholds
still exclude missing values, raising an error if no events remain or a
requested filter column is absent. Replot existing results to apply this
display fix; no scoring or permutation rerun is needed.

The four event-based plots expose their min-cells QC decision in the returned
`qc` dictionary. Their default `qc_only=True` changes only what is displayed;
it does not change the stored numerical scores or inference results. Set
`qc_only=False` to visualize all calculated events or sample rows.

For external tables, `passes_min_cells` accepts booleans, numeric 0/1 and text
`True`/`False`/`1`/`0` (case-insensitive, surrounding whitespace ignored).
Genuine missing values remain missing and do not pass the default QC filter.
Other values, including empty strings and literal `"NA"`, raise an error even
with `qc_only=False`. Normalization uses a copy and preserves the caller's table.
Tables without this optional column retain the existing behavior without QC
filtering.

Validated float16 plotting columns are promoted to float64 on a working copy
before sorting or ranking. This supports compressed external result tables
without changing their stored values or the caller's data; it cannot recover
precision already lost when the table was converted to float16.

## Numerical precision and underflow

Observation, compiled permutations and single-cell reaction plots share the
stable equivalent calculation `expm1(mean(log1p(expression)))`; a one-gene
reaction returns that gene's mean expression directly. Sender base scores use
`P_score * (P_score / (P_score + C_score))`, with the existing zero-denominator
rule. Events use `sqrt(sender_score) * sqrt(receiver_score)`. These forms avoid
losing a representable positive result through adding/subtracting 1, squaring
first, or multiplying before taking a square root. They preserve the existing
reaction model and positive-reference definition.

Actual positive-to-zero underflow is **reported and computation continues**.
There is no small-value cutoff and no arbitrary epsilon replacement. Checks
cover relevant pseudobulk means, reaction activity, abundance powers/products,
positive-reference normalization, sender base/exporter multiplication, and
event multiplication. Finite positive values remain usable however small;
mathematically zero inputs do not trigger an underflow notice. Negative values,
NaN, Inf and overflow still raise the existing errors.

Each affected `run_cell_mesh()` invocation emits one aggregated WARNING-level
logging message, including nested sample scoring and permutation workers. It
does not use Python warnings that could become exceptions under `-W error`.
The diagnostic is stored in `result.parameters["numerical_diagnostics"]` and
`result.events.attrs["numerical_diagnostics"]`, and is included in the existing
parameters JSON export. It contains `underflow_detected`,
`policy="report_and_continue"` and a `stages` list with `stage` and `n_values`.
Counts are occurrences across operations and permutations, not unique events.
Runs without detected underflow do not add this key.

Standalone availability and violin results expose the same diagnostic key in
their returned dictionaries; standalone receiver scores use DataFrame attrs.
Supplied receiver summaries that have positive expression fractions but zero
means are also reported as a possible mean-underflow inconsistency; the caller
remains responsible for cache provenance.

Rounded zeros continue through the existing scoring and permutation rules;
permutations are not skipped and the requested number of draws is retained.
This continuation does **not** recover unrepresentable positive values.
Affected references, `production_status` and other evidence states, scores and
p/FDR can differ from higher-precision arithmetic. In particular, a
`prior_no_expression` state associated with underflow is not proof of absent
expression. Consult the diagnostic stages before interpreting these results.
Structural sample NA remains distinct from computed zeros.

Both violin functions recognize a constant only when all values are exactly
equal. Nonconstant values are scaled to a unit range for density estimation and
then plotted in their original units. If density estimation fails, original
points replace the density and the existing `show_median` option still applies.
The returned `density_fallbacks` lists affected groups/reasons, and
`fallback_artists` contains the scatter artists. True constants retain their
actual-value horizontal line.

## Full API Reference

### `run_cell_mesh()` Parameters

| Parameter | Default | Description |
|---|---|---|
| `adata` | *required* | AnnData object containing single-cell expression data |
| `enzyme_metabolite` | `None` | Raw or normalized enzyme-metabolite DataFrame or CSV path (str/Path). None loads the latest built-in database |
| `metabolite_sensor` | `None` | Raw or normalized metabolite-sensor DataFrame or CSV path (str/Path). None loads the latest built-in database |
| `cell_type_key` | `"cell_type"` | Column name in adata.obs containing cell type annotations |
| `sample_key` | `None` | Column name in adata.obs containing sample annotations |
| `sample_mode` | `"pooled_stratified"` | `"pooled_stratified"` computes pooled cell-type pseudobulks; `"sample_aware"` computes and scores `(sample, cell type)` units separately before aggregating event scores across samples |
| `layer` | `None` | Name of expression layer to use. If None, uses adata.X |
| `min_expr_frac` | `None` | Optional receiver expression-fraction gate; must be in `[0, 1]` |
| `allow_self` | `True` | Python/NumPy boolean; whether to allow self-communication events (sender == receiver). Strings such as `"False"` are rejected |
| `n_perms` | `0` | Non-negative integer permutation count. 0 = no permutation |
| `random_state` | `0` | Non-negative Python/NumPy integer seed, validated even with zero permutations or empty events. Floats, booleans, strings, `None`, arrays and generator objects are rejected |
| `n_jobs` | `1` | Permutation worker threads; `-1` uses all available CPUs |
| `store_null_scores` | `False` | Python/NumPy boolean; store the full sample-aware event × permutation null matrix; p-values do not require it |
| `min_cells` | `100` | Positive-integer, QC-only cell-count threshold; observed units remain in every numerical calculation regardless of this flag |
| `sender_abundance_exponent` | `1.0` | Finite non-negative exponent applied to sender cell fraction before P/C/E construction; `0` disables sender abundance adjustment |
| `pce_reference` | `"mean"` | Strictly-positive P/C/E reference statistic; accepts `"mean"` or `"median"` |
| `export_weight` | `0.2` | Bounded exporter modulation weight in `[0, 1]`; `0` exactly disables the E contribution |
| `receiver_reference` | `"median"` | Strictly-positive, unweighted receiver-expression reference statistic; accepts `"median"` or `"mean"` |

`min_cells : int >= 1`
    QC-only cell-count threshold. In ``pooled_stratified`` mode it flags each
    observed cell type using its total count across the full AnnData object. In
    ``sample_aware`` mode it flags each observed ``(sample, cell type)`` unit
    independently. The threshold does not change pseudobulks, cell fractions,
    P/C/E values or references, sender/receiver/event scores, event coverage,
    or permutation inference. A unit below the threshold is still calculated
    and receives ``passes_min_cells=False``. NA in sample-level results means
    the event was not computable in that sample (for example, an endpoint was
    unobserved), not merely that its cell count was low. Evaluable zero
    production remains score 0 and contributes to the sample summaries.

`sample_mode : {"pooled_stratified", "sample_aware"}, default="pooled_stratified"`
    ``pooled_stratified`` computes pooled cell-type pseudobulks across all
    cells. When ``sample_key`` is provided, it is used only to stratify label
    permutations within samples.

    ``sample_aware`` computes pseudobulks, P/C/E scores, sender scores,
    receiver scores, and event scores separately for each ``(sample, cell
    type)`` unit, then aggregates event scores across samples.

`sender_abundance_exponent : float, default=1.0`
    Sender reaction activity is multiplied by
    ``cell_fraction ** sender_abundance_exponent`` before P/C/E construction.
    Use ``0`` for no abundance adjustment and a value between 0 and 1 for a
    sublinear abundance effect; values above 1 amplify abundance differences.

`pce_reference : {"mean", "median"}, default="mean"`
    Statistic applied to strictly positive abundance-adjusted P/C/E capacities
    for each metabolite and direction. ``"mean"`` represents the average
    active population contribution; ``"median"`` provides a more
    outlier-robust equal-cell-type reference.

`receiver_reference : {"median", "mean"}, default="median"`
    Statistic applied to strictly positive per-cell-type mean sensor
    expression. Observed cell types are equally weighted, so receiver cell
    count and fraction do not enter either the reference or score.

### Optional receiver summaries

`cellmesh.score.compute_sensor_scores()` can reuse `pseudobulk` (a DataFrame of
cell-type mean gene expression), `expr_frac` (a DataFrame of fractions expressing
each gene), and `cell_counts` (a Series of observed cell counts). Omit a summary
to compute it from the supplied AnnData and selected layer.

Supplied summaries follow these rules:

- Rows must contain exactly the cell types actually observed in the current
  AnnData, regardless of `min_cells` or all-zero expression. In sample-aware
  scoring this means the current sample, not the union across samples. Unused
  categorical levels are not required. Missing or extra groups raise errors;
  they are never silently dropped or filled with zeros.
- Both DataFrames must contain every sensor gene shared by the prior and
  `adata.var_names`. Additional genes from that expression matrix are allowed;
  genes absent from the matrix are rejected. Extra unrelated columns do not
  enter numerical validation or receiver scoring.
- Axis labels use the existing text/whitespace normalization. Empty/missing
  labels and duplicates after normalization raise errors. Different row/column
  orders are aligned by identifier without modifying the supplied objects.
- Relevant means must be real, finite and non-negative; expression fractions
  must be real, finite and in `[0,1]`, even when the expression gate is disabled.
  Numeric text may be converted; invalid strings, booleans, complex values,
  dates and missing numerical values are rejected.
- Counts must be finite positive integer values equal to the actual counts in
  `adata.obs`. Integral values such as `2.0` are accepted; values such as `1.9`
  are rejected before integer conversion. Numeric text follows the same value
  checks. Existing checks for optional `cell_fractions` continue to apply.

For example, if A and B have mean expression 2 and 8, supplying only A would
incorrectly change A's receiver score from `2/(2+5)` to `2/(2+2)`. This input
now reports missing B. Supply a complete summary or let the function compute it.
Validation checks observations and the relevant summary values without
recomputing supplied expression summaries. The caller must ensure the cached
values were derived from the same cells, expression data and layer; structural
and range checks cannot establish that provenance. Original cell-level
expression validation still runs when summaries are supplied.

`compute_metabolite_availability(return_intermediates=...)` also requires a
Python/NumPy boolean; strings and integer substitutes are rejected.

## Inspect the Packaged Database

```python
from cellmesh import load_cell_mesh_database

enzyme_metabolite, metabolite_sensor = load_cell_mesh_database()

print(enzyme_metabolite.head())
print(metabolite_sensor.head())
```

### Enzyme Table Columns

| Meaning | Standard column | Accepted aliases | Required field |
|---|---|---|---|
| Metabolite display name | `metabolite` | `standard_metName` | Yes |
| HMDB identifier | `hmdb_id` | `HMDB_ID` | Yes |
| Enzyme gene(s) | `gene` | `Gene_name` | Yes |
| Reaction identifier | `reaction` | `Reactions` | Yes |
| Reaction role | `role` | `Direction`, `direction` (value mapping below) | Yes |
| Evidence | `evidence_level` | None | No |
| Source | `source` | None | No |

Each required field needs one recognized column; supplying both a standard name
and its aliases is optional. Column names match the exact spellings above,
including case and spaces; arbitrary capitalization or padded headers are not
accepted aliases. Optional references and other custom columns are retained as
metadata. Enzyme metadata columns have no additional alias mapping.

The columns `metabolite`, `hmdb_id`, `gene`, `role`, and `reaction` are required.
The metabolite name may be missing or blank; it is only used for display.
Every reaction identifier must be non-empty; CELL MESH never infers a shared
enzyme complex from rows whose reaction identity is unknown.

Both Enzyme and Interaction CSVs validate every logical record before parsing:
its field count must equal the header width. Empty values need explicit empty
fields; short or extra-wide records, invalid quoting, and NUL characters raise
errors with the file and original line location. Quoted commas/newlines,
escaped quotes, literal identifiers, BOM, and supported CSV compression remain
valid. Blank lines outside quoted records are ignored. DataFrame inputs retain
their existing column/schema checks and may have a custom row index.

Standard `role` tables can be saved with `to_csv(index=False)` and loaded again
through `load_cell_mesh_database()` or `run_cell_mesh()`. Raw
`Direction`/`direction` tables use `product -> production`,
`substrate -> degradation`, and `exporter`/`export`/legacy `transporter` all map
to `role="export"`. If both role and direction are supplied, non-empty
direction values must agree with populated role values. Enzyme aliases are
checked before any rows are filtered: equivalent values are merged, blank
values are filled from non-empty aliases, and conflicts report the affected
rows and columns. This also applies to simultaneous `Direction`/`direction`
columns without `role`, and to `gene`/`Gene_name`, `hmdb_id`/`HMDB_ID`,
`reaction`/`Reactions`, and `metabolite`/`standard_metName`. Duplicate column
names and missing required columns raise an explicit error. After conversion,
the role is carried solely by the canonical `role` column so raw direction
values cannot override it during scoring. Reaction IDs and other textual
identifiers are read as text, preserving leading zeros. Existing
`source`, `evidence_level`, references, and custom metadata columns are retained.

In both `role` and `Direction`/`direction` tables, `gene` can contain one symbol
or multiple symbols separated by `;`, `,`, or `|`, such as
`G1[Enzyme]; G2[Enzyme]`. Database loading expands these fields to one gene per
row **before** matching expression genes, retaining the same reaction identity
and provenance. Separators inside `[evidence]` remain part of the annotation.
Existing `evidence_level` values take precedence over inline annotations;
otherwise inline evidence from equivalent gene aliases is retained, including
distinct annotations for the same gene. Consumed alias columns are removed so
they cannot contradict expanded rows on a later load. Empty gene entries are ignored. Repeated
normalization preserves the expanded records. Both public scoring entrypoints
and the permutation scorer consume this shared normalization.

Reaction genes are grouped by `canonical_hmdb_id + reaction + direction`, not
by metabolite name. Exact symbols are deduplicated within a reaction to define
its gene set and are combined by the equal-weight geometric mean. Reactions in
the same canonical HMDB ID and direction are compared by their complete gene
sets, ignoring gene order. Identical sets contribute once. Any strict subset is
omitted in favor of its superset, so only inclusion-maximal gene sets remain.
Sets that are not subsets of one another retain all their genes and all
contribute. The retained reaction activities are summed into one P, C, or E
value. Reaction-count metadata counts these maximal contributing gene sets.
Display names prefer the first non-empty Enzyme name for an HMDB ID, then the
first non-empty Interaction name, then the HMDB ID itself. Names are selected
from normalized scoring priors; standalone scoring uses its available prior.
Changing a name never changes event identity or sample aggregation.

Reaction gene-set comparisons use the **complete prior**, including genes not
measured in the expression matrix. Only after deduplication does each retained
reaction use its measured genes to calculate activity. Unmeasured genes are
excluded from the geometric mean; measured zero-expression genes remain in it.
A reaction with no measured genes has zero capacity and retains its prior
identity for counts and evidence-status reporting. Both public entrypoints and
compiled permutations use this order. Results from earlier versions may change
when filtering missing genes previously changed reaction subset relationships.

`compute_metabolite_availability()` also accepts this standard
`enzyme_metabolite` schema directly, as well as legacy direction-style inputs.
It shares enzyme normalization with `run_cell_mesh()`, including field aliases,
direction values, and conflict checks, instead of maintaining a separate
direction parser. Enzyme genes
within one reaction contribute equally; legacy/custom enzyme `weight` columns
are ignored and removed during prior validation.

### Production Evaluability and Zero Scores

Production availability is determined from the complete retained reaction
definitions and normalized `adata.var_names`, independently of the cell labels
or expression magnitude. Both observed and permuted scoring use this rule:

| Production evidence | `production_status` | `production_evaluable` | Scoring behavior |
|---|---|---|---|
| No production relation | `prior_missing` | False | No sender/event score |
| Production relations exist, but none of their genes are in the matrix | `prior_gene_unavailable` | False | No sender/event score |
| Some production genes are available, all calculated P capacities are zero | `prior_no_expression` | True | Keep zero sender/event scores |
| Some production genes are available and a cell type has positive P | `supported` | True | Score each cell type, including its zeros |

Both paths apply this eligibility mask before P/C/E reference normalization.
References still include all observed cell types in each scoring unit. Relevant
raw expression and computed capacities are validated before this selection;
production-ineligible rows skip reference calculation, rather than being
normalized and discarded only when extracting the observed events.

Partial gene availability retains the existing calculation from measured genes;
missing genes do not enter the reaction's geometric mean as zeros. Entirely
unmeasured reactions use internal numerical placeholders, but their zeros do
not establish production evaluability. If no enzyme prior gene matches the
matrix at all, `run_cell_mesh()` still raises an error before analysis.

These statuses summarize each metabolite across cell types in one scoring unit
(pooled data or one sample). For scored metabolites they appear in `metadata`,
whose index remains aligned to P/C/E and availability. The separate
`production_diagnostics` table includes `n_product_reactions` and both status
columns for **all Enzyme-prior metabolites**, including ones excluded from score
matrices. Each sample's availability result contains its own diagnostics. To
display unavailable production as NA alongside measured zeros:

```python
unit = res.availability_results  # pooled mode
# For sample_aware: unit = res.availability_results["availability_by_sample"]["D1"]
production_status = unit["production_diagnostics"][["production_status", "production_evaluable"]]
P_display = unit["P"].reindex(production_status.index)  # excluded rows become NA
```

In `sample_aware` mode, evaluable zero samples participate in medians, IQR,
`n_samples_coobserved`, cell-count/QC summaries, and the denominator of
`event_prevalence`. An absent sender or receiver cell type remains NA and is
excluded from that event's sample denominator. For example, scores
`[0.474342, 0, NA]` give a median of `0.237171`, two co-observed samples, and
positive prevalence `0.5`. Zero production never implies absent endpoint cells.

Both inference modes retain measured all-zero metabolites as zero-score events.
When permutations are run, a zero observed score has `perm_pvalue=1`; with
`n_perms=0`, p-values remain NA. Previous versions discarded all-zero P rows,
which could omit evaluable samples from summaries. Rerun affected analyses:
sample medians/prevalence can change, and adding evaluable zero events expands
the FDR correction family in either mode.

Availability currently follows the shared gene columns in the supplied AnnData
object. If upstream processing filled unavailable measurements with zeros,
the matrix alone cannot recover their missingness; that requires separate
measurement-availability information. Zeros are never used to infer missing genes.

### Sensor Table Columns

| Meaning | Standard column | Accepted aliases | Required field |
|---|---|---|---|
| Metabolite display name | `metabolite` | `standard_metName`, `standard_metname` | Yes |
| HMDB identifier | `hmdb_id` | `HMDB_ID` | Yes |
| Sensor gene | `sensor_gene` | `Gene_name`, `gene_name`, `gene` | Yes |
| Sensor type | `sensor_type` | `Annotation`, `annotation` | Yes |
| Source | `source` | `Database source`, `database_source` | No |
| Protein name | `protein_name` | `Protein_name` | No |
| Reference | `reference` | `Reference` | No |
| Evidence | `evidence_level` | None | No |

The first four fields are required; provide one recognized column for each.
Column names are matched exactly. The current `Interaction1.41.csv` supplies
`Gene_name`, which is renamed to `sensor_gene` and used for expression matching.
Interaction rows represent one sensor gene each. Multiple sensors must be
separate rows: Interaction gene fields are not split on `;`, `,`, or `|`, and do
not parse Enzyme-style `[evidence]` suffixes.

Both `load_cell_mesh_database()` and `run_cell_mesh()` apply these rules to CSV
and DataFrame inputs before gene filtering or pair deduplication:

1. Exact duplicate column names raise an error, even when their values agree.
   CSV headers are checked before pandas can rename duplicates to `.1`.
2. All coexisting standard/alias columns are compared row by row, including
   aliases without a standard column. Equivalent values merge; missing values
   are filled from populated aliases. Conflicting populated values raise an
   error identifying the field, columns, values, and 1-based data-row positions
   (plus DataFrame index labels).
3. HMDB values are stripped and compared in uppercase. Genes, metabolite names,
   and source/protein/reference metadata are compared as stripped, case-sensitive
   text. Sensor types are compared by their normalized category. Missing types
   are filled before applying the `Other receptor` fallback.
4. Equivalent values retain the populated standard value, or the first populated
   alias in the order listed above. This preserves existing display/provenance
   text; gene, HMDB, and type values then undergo their normal normalization.
   Consumed alias columns are removed. Other metadata, such as `ID`,
   `Interaction_mode`, and `STITCH_evidence`, is retained.
5. `Annotation`/`annotation` originals fill missing `evidence_level` values.
   Distinct annotation texts with the same normalized category are joined with
   `; ` in alias order; duplicate texts after trimming are kept once. Populated
   `evidence_level` values take precedence. Evidence never weights a score.

For example, `sensor_gene=S1` with `Gene_name=S1` merges; a blank `sensor_gene`
with `Gene_name=S1` is filled; `sensor_gene=S1` with `Gene_name=S2` raises.
These are different-name aliases, so checking duplicate headers alone cannot
replace the row-wise comparison. Both gene aliases empty still remove the row;
records without an HMDB identifier or a measured sensor gene are excluded during
runtime validation. Use the recommended three sensor categories explicitly.

CSV inputs use UTF-8 (a BOM is accepted), with column names in the first row.
Recognized text fields preserve leading zeros and literal `NA` values; actual
blank cells remain missing. For DataFrames, required fields must be columns,
not only index levels. Save normalized tables with `to_csv(index=False)`.

The packaged interaction database does not define a quantitative sensor
`weight`. Legacy/custom sensor `weight` columns are ignored and removed during
prior validation; receiver scoring depends only on expression and its selected
positive-value reference.

Database loading and runtime prior validation guarantee one row per
`canonical_hmdb_id + sensor_gene`. Repeated rows with the same sensor type keep
their first metadata record and cannot duplicate communication events. If the
same HMDB/sensor pair has multiple sensor types, loading raises an error before
deduplication for both CSV and DataFrame inputs, including pairs whose genes
are absent from a later expression dataset. Silently choosing a type by row
order would make sensor-type-specific FDR ambiguous.

## Supported Sensor Types

After normalization from `sensor_type`, `Annotation`, or `annotation`:

- `Cell surface receptor`
- `Transporter`
- `Other receptor`

Type values are stripped and matched without case sensitivity. Text containing
`cell surface` or `surface receptor` maps to `Cell surface receptor`; otherwise,
text containing `transport` maps to `Transporter`. Other or missing values map
to `Other receptor`. Simultaneous nonempty aliases must agree after this mapping;
an explicit `Other receptor` is a populated type, not a blank eligible for fallback.

## Main Outputs

### `res.events`
Contains one row per communication event. Important columns:
- `sender`, `receiver`: Cell type pair
- `metabolite`, `hmdb_id`: Metabolite information
- `sensor_gene`, `sensor_type`: Sensor information (sensor_type is one of "Cell surface receptor", "Transporter", "Other receptor")
- In pooled output, `metabolite_availability`, `sensor_score`, and
  `sensor_expr_frac` are the event components calculated from pooled cell-type
  summaries, and `cell_mesh_score` is their geometric-mean score.
- In sample-aware aggregate output, the corresponding descriptive component
  fields are `metabolite_availability_median`, `sensor_score_median`, and
  `sensor_expr_frac_median`.
- `sender_n_cells`, `receiver_n_cells`: Observed sender/receiver group sizes
- `sender_passes_min_cells`, `receiver_passes_min_cells`: Endpoint-specific
  min-cells QC flags
- `passes_min_cells`: Pair-level display QC. In pooled mode both endpoints must
  pass. In sample-aware aggregate output it is true when at least one
  co-observed sample has both endpoints passing.
- `cell_mesh_score`: pooled/sample-level geometric-mean score in pooled output;
  in sample-aware aggregate output it equals `event_score_median`, the median
  of sample-level geometric-mean scores. It is not recomputed from the two
  separately aggregated component medians.
- `perm_pvalue`: One-sided empirical p-value from label permutation
- `fdr_global`: Benjamini-Hochberg correction across all events
- `fdr_sensor_type`: correction separately within each sensor type
- In `sample_aware` mode, `permutation_mode` records
  `"within_sample_label_shuffle"`.
- In `sample_aware` mode, `n_samples_passing_min_cells` counts co-observed
  samples where both endpoints pass and `min_cells_pass_prevalence` divides
  that count by `n_samples_coobserved`. These are QC summaries only.

### Other Outputs
- `res.sender_scores`: `(metabolite, hmdb_id)` × cell type matrix of availability
  scores; sample-aware mode groups by HMDB ID for the across-sample median, then
  reattaches the name in the display index.
- `res.receiver_scores`: Table of sensor scores per metabolite-sensor-cell type
  combination; sample-aware component columns use across-sample medians grouped
  by `hmdb_id + sensor_gene + receiver`. Names and sensor types are annotations.
- `res.availability_results`: Dictionary containing abundance-adjusted raw
  `P`/`C`/`E` capacities, continuous `P_score`/`C_score`/`E_score` matrices,
  per-metabolite `P_ref`/`C_ref`/`E_ref` positive-value references,
  `pce_reference`, `receiver_reference`, per-cell pseudobulk means,
  `cell_counts`, `cell_fractions`, `celltype_qc`, `sender_abundance_weights`,
  `sender_abundance_exponent`, reaction counts, `consumption_status`,
  `export_status`, `production_status`, `production_evaluable`, and other
  intermediate results. `metadata` stays aligned to score matrices;
  `production_diagnostics` also includes unscorable Enzyme-prior metabolites.
  Only the formal continuous
  score matrices are exposed; obsolete signed contrasts and score aliases are
  no longer returned.
- In pooled mode, `res.celltype_qc` and
  `res.availability_results["celltype_qc"]` report `n_cells`, `cell_fraction`,
  and `passes_min_cells` for every observed cell type.
- In `sample_aware` mode, `res.sample_validation` reports `sample`,
  `cell_type`, `n_cells`, `cell_fraction`, `passes_min_cells`, observed-unit
  counts, and explicitly named min-cells pass-count summaries. Sample event
  rows carry the same three endpoint/pair flags as pooled events. NA indicates
  that the sample-level event was not computable (for example, an endpoint was
  unobserved), never merely that an observed unit failed min-cells QC or had
  evaluable zero production. Computed zeros keep their cell counts and QC flags.

### Export Results

Call `res.to_csv(prefix)` explicitly to export results; `run_cell_mesh()` does
not write these files automatically. The parent directory must already exist.

```python
import json
import pandas as pd

res.to_csv("analysis")
sender_scores = pd.read_csv(
    "analysis.sender_scores.csv",
    dtype={"metabolite": "string", "hmdb_id": "string"},
    float_precision="round_trip",
).set_index(["metabolite", "hmdb_id"])
with open("analysis.parameters.json", encoding="utf-8") as stream:
    parameters = json.load(stream)
```

Exports include `<prefix>.events.csv`, `.sender_scores.csv`,
`.receiver_scores.csv`, `.parameters.json`, and `.manifest.json`. When present, `.celltype_qc.csv`,
`.sample_validation.csv`, `.sample_sender_scores.csv`,
`.sample_receiver_scores.csv`, and `.sample_events.csv` are also written.
These exact standard filenames form the exporter's reserved namespace:
existing files are overwritten, and standard optional CSVs absent from the
current result are removed. This includes older exports without a manifest.
For example, sample-aware → pooled re-export with the same prefix removes the
four stale sample-level CSVs. An absent optional QC table is handled likewise.
Other prefixes, nonstandard suffixes (such as `analysis.notes.txt`) and database
files are not scanned or cleaned. Do not store hand-created files under the
same exact standard result names: their origin cannot be inferred from the
name, so they are subject to the same replacement/removal rules.

`<prefix>.manifest.json` lists this export's relative filenames (including
itself), each table's columns and row count, and a format version. Cleanup uses
only the exporter's fixed filename list, never paths from an existing manifest.
Symbolic links and directories at standard target paths are rejected.

All headers and parameters are prevalidated, and all CSV/JSON contents are
written to a temporary directory beside the destination before replacing any
result. Parameters use strict JSON: nested NaN/Infinity values and unsupported
objects are rejected; valid NumPy scalar parameters remain supported. Genuine
missing values in result tables still become empty CSV fields, while numeric
zero remains zero. The manifest is published last, after replacements and
stale-file cleanup. Ordinary write/cleanup failures trigger rollback of changed
files from temporary backups. Prevalidation or staging failure leaves the
previous export untouched. This is not a crash-safe multi-file transaction:
do not write concurrently to the same prefix or read while it is being updated.

Sender score CSVs contain `metabolite` and `hmdb_id` as ordinary columns;
sample sender scores also contain `sample`. These identifiers are written once,
even if the score table was already converted with `reset_index()`.
Event and receiver tables retain their existing columns, and QC index labels
are preserved as columns. Conflicting column names raise an error before any
files are written. Empty score tables retain identifier headers. Empty event
tables retain the same columns as nonempty tables within the same inference
mode, including `inference_mode`. Both modes set `events.attrs["n_perms_completed"]`
to zero when permutations are skipped or there are no events, and explicitly
set `null_scores_stored` (sample-aware follows the requested storage setting;
pooled does not store null matrices).

This exports the result tables and parameters, including sample-level missing
values as empty CSV fields. CSV does not store pandas dtype metadata;
`availability_results` intermediates, DataFrame attributes, and stored
permutation null matrices are not included.

## Notes
- **Transcriptomics-only**: CELL MESH estimates metabolite availability using expression proxies and prior knowledge. Direct metabolomics, spatial data, or perturbation experiments should be used to validate predictions.
- **Sender abundance adjustment**: P/C/E reaction activity is multiplied by
  `cell_fraction ** sender_abundance_exponent` before positive-reference
  saturation normalization
- **P/C/E reference**: Strictly-positive arithmetic mean by default; set
  `pce_reference="median"` for the robust median option
- **Min-cells QC**: `min_cells` never filters observed units from calculation or
  inference. It supplies QC flags, while event visualizations use
  `qc_only=True` by default to display only passing rows
- **Export evidence**: Normalized E applies the bounded factor
  `(1-export_weight) + export_weight*E_effective`; missing or unavailable
  exporter evidence uses fixed `E_effective=0.5`, while measured all-zero
  exporter capacity uses `0`
- **Sensor scoring**: Uses `R / (R + R_ref)` with the strictly-positive,
  unweighted cell-type median as the default reference; set
  `receiver_reference="mean"` for the arithmetic-mean option. Receiver
  abundance is not part of the score or reference
- **Communication score**: Geometric mean ensures both sender and receiver have meaningful scores
- **Multiple-testing correction**: Both inference modes report `fdr_global`
  across all events and `fdr_sensor_type` corrected separately within each
  sensor type.
- **Permutation null**: Empirical p-values compare each observed event key (`sender`, `receiver`, `hmdb_id`, `sensor_gene`) against the same key after cell-type label permutation. Names and sensor types are annotations, excluded from event matching and null indices; sensor type still determines the stratified FDR group. Static priors and relevant expression columns are compiled once; permutation scoring evaluates only retained observed keys and accumulates exceedance counts online. `n_jobs` controls deterministic batch parallelism. In sample-aware mode, the full null matrix is stored only with `store_null_scores=True`, with the four identity fields as its index levels.
