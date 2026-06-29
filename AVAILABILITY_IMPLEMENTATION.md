# Sender Score Implementation

Only cell types with `n_cells >= min_cells` are eligible. P/C/E reaction scores
retain the existing reaction-gene aggregation and reaction summation logic.

For each metabolite and each P/C/E vector across eligible cell types:

\[
b=\operatorname{median}(x),\qquad
D(x_c;b)=
\begin{cases}
(x_c-b)/(x_c+b), & x_c+b>\varepsilon_{num}\\
0, & \text{otherwise}
\end{cases}
\]

`eps_num` defaults to `1e-12` and is numerical protection only. NaN values map
to zero contrast; all-NaN and all-zero vectors return all zeros. Contrasts are
clipped to `[-1, 1]`.

The sender score is:

\[
S^{sender}_{m,c}=p^+_{m,c}F^E_{m,c}F^C_{m,c}
\]

- \(p^+=\max(0,d^P)\), so production above the cell-type median is required.
- \(F^E=1+e^+\) when an exporter prior exists, otherwise 1.
- \(F^C=1-c^+\) when a consumption/substrate prior exists, otherwise 1.
- Missing exporter or consumption priors are neutral.

Raw C is the expression-derived proxy for metabolite consumption-enzyme
ability. The C-derived output is named `relative_consumption_support` because
it retains only the positive deviation of that proxy above the eligible
cell-type median. It is not measured extracellular clearance flux.

Intermediates include `P`, `C`, `E`, `P_contrast`, `C_contrast`, `E_contrast`,
`P_plus`, `E_plus`, `relative_consumption_support`, `pseudobulk`, `expr_frac`,
and `cell_counts`.
