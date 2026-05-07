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

The enzyme table maps `Direction = product` to `production` and `Direction = substrate` to `degradation`. The interaction table maps annotations to CELL MESH sensor classes:

- Cell surface receptor -> `surface_receptor`
- Other receptor -> `surface_receptor`
- Transporter -> `transporter`
- Nuclear receptor -> `nuclear_receptor`
- Intracellular sensor -> `intracellular_sensor`

## Expression aggregation

For each cell group \(c\) and gene \(g\), CELL MESH computes mean expression:

\[
\bar{x}_{g,c} = \frac{1}{|C_c|}\sum_{i \in C_c} x_{i,g}
\]

and expression fraction:

\[
\phi_{g,c} = \frac{1}{|C_c|}\sum_{i \in C_c} I(x_{i,g} > 0)
\]

Mean expression is z-scored across cell groups:

\[
z_{g,c} = \frac{\bar{x}_{g,c} - \mu_g}{\sigma_g + \epsilon}
\]

## Metabolite role scores

For metabolite \(m\), cell group \(c\), and role \(r\), CELL MESH computes:

\[
A_{m,c}^{(r)} = \operatorname{Agg}_{g \in G_{m,r}}(w_g z_{g,c})
\]

where \(r\) can be `production`, `degradation`, `export`, `import`, or `usage`.

The default aggregation is a weighted mean:

\[
A_{m,c}^{(r)} = \frac{\sum_g w_g z_{g,c}}{\sum_g |w_g|}
\]

A `softmin` option is available for bottleneck-like pathways.

## Sender score

CELL MESH models sender-side metabolite availability as:

\[
\tilde{S}_{m,c_s} =
\alpha_1 A_{m,c_s}^{prod}
- \alpha_2 A_{m,c_s}^{deg}
+ \alpha_3 A_{m,c_s}^{export}
+ \alpha_4 Q_{m,c_s}^{send}
\]

where \(Q\) is a cell-group specificity term. The bounded sender score is:

\[
S_{m,c_s} = \sigma(\tilde{S}_{m,c_s})
\]

## Receiver score with sensor hierarchy

CELL MESH differs from a unified metabolite-sensor scoring scheme by using different receiver models for different sensor classes.

### Surface receptor

\[
R^{surface}_{m,s,c_r} = \sigma(\beta_1 z_{s,c_r} + \beta_2 Q_{s,c_r}^{recv})
\]

### Transporter

\[
R^{transporter}_{m,s,c_r} = \sigma(
\beta_1 z_{s,c_r}
+ \beta_2 Q_{s,c_r}^{recv}
+ \beta_3 A_{m,c_r}^{import}
+ \beta_4 A_{m,c_r}^{usage}
)
\]

### Nuclear/intracellular sensor

\[
R^{intracellular}_{m,s,c_r} = \sigma(
\beta_1 z_{s,c_r}
+ \beta_2 Q_{s,c_r}^{recv}
+ \beta_3 A_{m,c_r}^{import}
+ \beta_5 A_{m,c_r}^{usage}
)
\]

## CELL MESH event score

The final event score is:

\[
\text{CELL\_MESH}_{c_s \to c_r}^{m,s} =
S_{m,c_s} \times R_{m,s,c_r}^{type(s)} \times W_{m,s}^{prior}
\]

In the output table this is stored as `cell_mesh_score`. The column `communication_score` is also retained as a backward-compatible alias.

## Permutation testing

CELL MESH optionally estimates empirical p-values by shuffling cell group labels. If `sample_key` is provided, labels are shuffled within sample.

\[
p = \frac{1 + \sum_b I(C_b \ge C_{obs})}{B + 1}
\]

P-values are adjusted using Benjamini-Hochberg FDR.
