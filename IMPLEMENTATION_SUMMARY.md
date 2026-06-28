# Implementation Summary

CELL MESH filters to cell types with at least `min_cells` cells before
pseudobulk, P/C/E, sender, and receiver calculations.

- Reaction gene aggregation and P/C/E reaction summation are unchanged.
- P/C/E sender components use bounded deviation from the equal-weight
  eligible-cell-type median.
- Positive production contrast is the required sender anchor.
- Exporter positive contrast adds support only when an exporter prior exists.
- The C matrix represents the expression-derived level of
  metabolite-consuming enzyme complexes.
- Its positive median deviation penalizes only when a consumption/substrate
  prior exists.
- Missing exporter or consumption priors use neutral factor 1.
- Receiver scores use positive bounded median contrast of sensor pseudobulk.
- `min_expr_frac=None` is the default; a numeric value enables an optional
  receiver-only gate.
- Event score is `sqrt(sender_score * receiver_score)`.
- Events record `sender_n_cells` and `receiver_n_cells`.

The derived C component is exposed as `relative_consumption_support`: it is the
positive deviation of the consumption ability proxy above the eligible
cell-type median, not a direct extracellular clearance-flux measurement.
