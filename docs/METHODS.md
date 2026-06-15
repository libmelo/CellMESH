# CELL MESH methods

**CELL MESH** stands for **Metabolite-mediated Event Scoring with Sensor Hierarchies**.

CELL MESH infers candidate metabolite-mediated cell-cell communication events from single-cell expression data. Each event is defined as:

\[
(c_s, c_r, m, s)
\]

where \(c_s\) is the sender cell group, \(c_r\) is the receiver cell group, \(m\) is a metabolite, and \(s\) is a sensor gene.

## Packaged database

This version packages two uploaded CSV files:

1. `enzyme_test.csv`: metabolite-enzyme-reaction table
2. `interaction_test.csv`: metabolite-sensor interaction table

At runtime, `load_cell_mesh_database()` normalizes these files into:

- `enzyme_metabolite`
- `metabolite_sensor`

Both priors must provide `hmdb_id`. Records without `hmdb_id` are excluded before scoring, and sender availability is matched to receiver sensors by exact `(metabolite, hmdb_id)`.

The enzyme table maps reaction directions into the three enzyme roles used by the availability model:

- `product` -> `production`
- `substrate` -> `degradation`
- `exporter` -> `export`

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

Reaction scores are aggregated into three metabolite matrices: production \(P\), consumption \(C\), and efflux \(E\). Each matrix is robust min-max normalized across cell groups for each metabolite. The `eps` parameter is the denominator stabilizer in this normalization.

\[
P^{norm}_{m,c}, C^{norm}_{m,c}, E^{norm}_{m,c} \in [0, 1]
\]

The sender-side metabolite availability is:

\[
A_{m,c_s} =
P^{norm}_{m,c_s}
\times (1 - C^{norm}_{m,c_s})^\beta
\times (0.8 + 0.2 E^{norm}_{m,c_s})
\]

## Receiver score

Receiver-side sensor score uses robust min-max normalized pseudobulk sensor gene expression across cell groups. If the sensor expression fraction in the receiver group is below `min_expr_frac`, the score is set to 0.

\[
R_{m,s,c_r} \in [0, 1]
\]

## CELL MESH event score

\[
\text{CELL\_MESH}_{c_s \to c_r}^{m,s} =
\sqrt{A_{m,c_s} R_{m,s,c_r}}
\]

In the output table this is stored as `cell_mesh_score`. The column `communication_score` is also retained as a backward-compatible alias.

## Permutation testing

CELL MESH optionally estimates empirical p-values by shuffling cell group labels. If `sample_key` is provided, labels are shuffled within sample.

\[
p = \frac{1 + \sum_b I(C_b \ge C_{obs})}{B + 1}
\]

The current null is defined at the full event-key level: `sender`, `receiver`, `metabolite`, `hmdb_id`, `sensor_gene`, and `sensor_type` must match between observed and permuted events before scores are compared. P-values are adjusted using Benjamini-Hochberg FDR separately within each sensor type.
