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
P/C/E matrices. For each metabolite and each direction \(X \in \{P,C,E\}\),
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

The sender score is:

```text
sender_score(m,t) = P_score(m,t)^2 / (P_score(m,t) + C_score(m,t))
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

`E_score` is still calculated and reported as export-support evidence, together
with the same four-state `export_status`, but it does not enter the formal
sender score. This avoids assigning an unsupported numerical boost when
exporter coverage in the prior is sparse, and keeps the normalized export
evidence bounded below one.

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
Sensor types are read from the `Annotation` column of the selected metabolite-sensor
database. CELL MESH preserves the database labels instead of imposing a fixed type
catalogue; the current default database is `Interaction4.0.csv`.

Both inference modes report global FDR and sensor-type-specific FDR.

## Included Data
The package includes built-in prior databases and a small AnnData example file:

- `cellmesh/data/Enzyme2.0.csv`: current default enzyme/reaction/metabolite table
- `cellmesh/data/Interaction4.0.csv`: current default metabolite-sensor table
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
`Enzyme2.0.csv` and `Interaction4.0.csv`. If no versioned files are available,
the loader falls back to the legacy `enzyme_test.csv` / `interaction_test.csv`
pair. The comprehensive walkthrough notebook uses `Enzyme_new.csv`,
`Interaction1.0.csv`, and `test_single_cell.h5ad` so its calculations are fully
reproducible from packaged files.

## Install

From a local checkout:

```bash
pip install -e .
```

For running the walkthrough notebook:

```bash
pip install -e ".[notebook]"
```

For development and tests:

```bash
pip install -e ".[dev]"
pytest -q
```

10X directory loading uses Scanpy and can be enabled with:

```bash
pip install -e ".[scanpy]"
```

## Basic Usage

```python
from cellmesh import run_cell_mesh

res = run_cell_mesh(
    adata,
    cell_type_key="cell_type",
    sample_key="sample",      # optional, for within-sample permutation
    sample_mode="pooled_stratified",
    layer="lognorm",          # optional, use specific expression layer
    n_perms=1000,             # optional, number of permutations for p-value calculation
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

`adata.obs[cell_type_key]` must contain a non-empty label for every cell, and
`adata.var_names` must be non-empty and unique. Invalid labels are rejected
before pseudobulk, abundance, or reference calculations begin.

## Visualization

Visualization functions are available from the top-level package and from the
`cellmesh.plotting` implementation module.

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
# endpoint was unobserved or no production availability was available); a low
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
# normalized production and consumption. Export is returned as support-only
# evidence and does not change this color.
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
The four event-based plots expose their min-cells QC decision in the returned
`qc` dictionary. Their default `qc_only=True` changes only what is displayed;
it does not change the stored numerical scores or inference results. Set
`qc_only=False` to visualize all calculated events or sample rows.

## Full API Reference

### `run_cell_mesh()` Parameters

| Parameter | Default | Description |
|---|---|---|
| `adata` | *required* | AnnData object containing single-cell expression data |
| `enzyme_metabolite` | `None` | Enzyme-metabolite prior table. If None, uses built-in database |
| `metabolite_sensor` | `None` | Metabolite-sensor prior table. If None, uses built-in database |
| `cell_type_key` | `"cell_type"` | Column name in adata.obs containing cell type annotations |
| `sample_key` | `None` | Column name in adata.obs containing sample annotations |
| `sample_mode` | `"pooled_stratified"` | `"pooled_stratified"` computes pooled cell-type pseudobulks; `"sample_aware"` computes and scores `(sample, cell type)` units separately before aggregating event scores across samples |
| `layer` | `None` | Name of expression layer to use. If None, uses adata.X |
| `min_expr_frac` | `None` | Optional receiver expression-fraction gate; must be in `[0, 1]` |
| `allow_self` | `True` | Whether to allow self-communication events (sender == receiver) |
| `n_perms` | `0` | Non-negative integer permutation count. 0 = no permutation |
| `random_state` | `0` | Random seed for reproducibility |
| `min_cells` | `100` | Positive-integer, QC-only cell-count threshold; observed units remain in every numerical calculation regardless of this flag |
| `sender_abundance_exponent` | `1.0` | Finite non-negative exponent applied to sender cell fraction before P/C/E construction; `0` disables sender abundance adjustment |
| `pce_reference` | `"mean"` | Strictly-positive P/C/E reference statistic; accepts `"mean"` or `"median"` |
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
    unobserved or no production availability was available), not merely that
    its cell count was low.

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

## Inspect the Packaged Database

```python
from cellmesh import load_cell_mesh_database

enzyme_metabolite, metabolite_sensor = load_cell_mesh_database()

print(enzyme_metabolite.head())
print(metabolite_sensor.head())
```

### Enzyme Table Columns
```
metabolite, hmdb_id, gene, role, evidence_level, source, reaction
```

`metabolite`, `hmdb_id`, `gene`, `role`, and `reaction` are required. Every
reaction identifier must be non-empty; CELL MESH never infers a shared enzyme
complex from rows whose reaction identity is unknown.

Reaction genes are grouped by `canonical_hmdb_id + reaction + direction`, not
by metabolite name. Within that group, exact gene symbols are deduplicated in
first-seen order and combined by the equal-weight geometric mean. Multiple
reactions in the same direction are then summed into a single P, C, or E value
for the HMDB ID. The first enzyme-prior metabolite name observed for an HMDB ID
is retained only as display metadata.

`compute_metabolite_availability()` also accepts this standard
`enzyme_metabolite` schema directly. Legacy direction-style inputs are only
kept as a compatibility path for low-level availability tests. Enzyme genes
within one reaction contribute equally; legacy/custom enzyme `weight` columns
are ignored and removed during prior validation.

### Sensor Table Columns
```
ID, HMDB_ID, standard_metName, Gene_name, Protein_name, Annotation, Database source, Reference
```

The packaged interaction database does not define a quantitative sensor
`weight`. Legacy/custom sensor `weight` columns are ignored and removed during
prior validation; receiver scoring depends only on expression and its selected
positive-value reference.

Runtime prior validation also guarantees one row per
`canonical_hmdb_id + sensor_gene`. Repeated rows with the same sensor type keep
their first metadata record and cannot duplicate communication events. If the
same HMDB/sensor pair has multiple sensor types, validation raises an error
because silently choosing one would make sensor-type-specific FDR ambiguous.

## Supported Sensor Types

From the `Annotation` column in the packaged interaction CSV:
- `Cell surface receptor`
- `Transporter`
- `Other receptor`

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
  scores; sample-aware mode uses the across-sample median.
- `res.receiver_scores`: Table of sensor scores per metabolite-sensor-cell type
  combination; sample-aware component columns use across-sample medians.
- `res.availability_results`: Dictionary containing abundance-adjusted raw
  `P`/`C`/`E` capacities, continuous `P_score`/`C_score`/`E_score` matrices,
  per-metabolite `P_ref`/`C_ref`/`E_ref` positive-value references,
  `pce_reference`, `receiver_reference`, per-cell pseudobulk means,
  `cell_counts`, `cell_fractions`, `celltype_qc`, `sender_abundance_weights`,
  `sender_abundance_exponent`, reaction counts, `consumption_status`,
  `export_status`, and other intermediate results. Only the formal continuous
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
  unobserved or no production availability was available), never merely that
  an observed unit failed min-cells QC.

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
- **Export evidence**: E is normalized and reported but does not enter the
  formal sender score
- **Sensor scoring**: Uses `R / (R + R_ref)` with the strictly-positive,
  unweighted cell-type median as the default reference; set
  `receiver_reference="mean"` for the arithmetic-mean option. Receiver
  abundance is not part of the score or reference
- **Communication score**: Geometric mean ensures both sender and receiver have meaningful scores
- **Multiple-testing correction**: Both inference modes report `fdr_global`
  across all events and `fdr_sensor_type` corrected separately within each
  sensor type.
- **Permutation null**: Empirical p-values compare each observed full event key (`sender`, `receiver`, `metabolite`, `hmdb_id`, `sensor_gene`, `sensor_type`) against the same key after cell-type label permutation, with FDR stratified by sensor type.
