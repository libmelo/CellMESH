# CELL MESH methods

**CELL MESH** stands for **Metabolite-mediated Event Scoring with Sensor Hierarchies**.

CELL MESH infers candidate metabolite-mediated cell-cell communication events from single-cell expression data. Each event is defined as:

\[
(c_s, c_r, m, s)
\]

where \(c_s\) is the sender cell group, \(c_r\) is the receiver cell group, \(m\) is a metabolite, and \(s\) is a sensor gene.

## Packaged database

This version packages versioned prior CSV files plus fixed walkthrough/test files:

1. `Enzyme1.0.csv`: versioned metabolite-enzyme-reaction table
2. `Interaction1.0.csv`: versioned metabolite-sensor interaction table
3. `enzyme_test.csv`: legacy small enzyme prior used for compatibility tests
4. `interaction_test.csv`: legacy small sensor prior used for compatibility tests
5. `Enzyme_new.csv`: walkthrough/test enzyme prior
6. `test_single_cell.h5ad`: walkthrough/test single-cell data

At runtime, `load_cell_mesh_database()` independently selects the highest
packaged enzyme file named `Enzyme<version>.csv` and the highest packaged
interaction file named `Interaction<version>.csv`, then normalizes the selected
prior files into:

- `enzyme_metabolite`
- `metabolite_sensor`

Both priors must provide `hmdb_id`. Records without `hmdb_id` are excluded before scoring, and sender availability is matched to receiver sensors by exact `(metabolite, hmdb_id)`.

The normalized `enzyme_metabolite` table is the canonical enzyme prior passed through `run_cell_mesh()`. It uses three roles:

- `production`
- `degradation`
- `export`

Packaged enzyme files that contain reaction directions are normalized into these roles during `load_cell_mesh_database()`:

- `product` -> `production`
- `substrate` -> `degradation`
- `exporter` -> `export`

Availability scoring then maps those roles back to internal P/C/E directions. This is an internal normalization step, not a separate public reaction-table API.

The interaction table maps annotations to the three supported sensor classes:

- `Cell surface receptor`
- `Transporter`
- `Other receptor`

## Expression aggregation

For each cell group \(c\) and gene \(g\), CELL MESH computes mean expression:

\[
\bar{x}_{g,c} = \frac{1}{|C_c|}\sum_{i \in C_c} x_{i,g}
\]

and expression fraction:

\[
\phi_{g,c} = \frac{1}{|C_c|}\sum_{i \in C_c} I(x_{i,g} > 0)
\]

For each reaction, CELL MESH scales each valid gene expression value by its gene weight and then applies an ordinary geometric mean:

\[
\operatorname{rxn}_{r,c} =
\operatorname{gmean}_{g \in G_r}(w_g \bar{x}_{g,c} + 1) - 1
\]

This weight semantics is a per-gene expression scaling step, not the normalized weighted geometric mean \(\exp(\sum w \log x / \sum w)\).

## Metabolite availability

Only cell types with at least `min_cells` cells are eligible. Reaction scores
are aggregated into unchanged production \(P\), consumption \(C\), and efflux
\(E\) matrices. For each non-negative vector across eligible cell types:

\[
D(x_c;b)=\frac{x_c-b}{x_c+b}, \qquad b=\operatorname{median}_c(x_c)
\]

with a numerical zero-denominator guard. Let \(p^+,c^+,e^+\) be the positive
parts of the P/C/E contrasts. The sender score is:

\[
A_{m,c}=p^+_{m,c}F^E_{m,c}F^C_{m,c}
\]

where \(F^E=1+e^+\) when an exporter prior exists and 1 otherwise, while
\(F^C=1-c^+\) when a consumption/substrate prior exists and 1 otherwise.
Production above the median is therefore required. The C term represents
relative consumption support or turnover/consumption context, not necessarily
true extracellular clearance flux.

## Receiver score

Receiver-side sensor score is the positive bounded median contrast of
pseudobulk sensor expression across eligible cell types. If `min_expr_frac` is
not `None`, it acts as an optional expression-fraction gate.

\[
R_{m,s,c_r} \in [0, 1]
\]

## CELL MESH event score

\[
\text{CELL\_MESH}_{c_s \to c_r}^{m,s} =
\sqrt{A_{m,c_s} R_{m,s,c_r}}
\]

In the output table this is stored as `cell_mesh_score`.

## Permutation testing

CELL MESH optionally estimates empirical p-values by shuffling cell group labels. If `sample_key` is provided, labels are shuffled within sample.

\[
p = \frac{1 + \sum_b I(C_b \ge C_{obs})}{B + 1}
\]

The current null is defined at the full event-key level: `sender`, `receiver`, `metabolite`, `hmdb_id`, `sensor_gene`, and `sensor_type` must match between observed and permuted events before scores are compared. P-values are adjusted using Benjamini-Hochberg FDR separately within each sensor type.
