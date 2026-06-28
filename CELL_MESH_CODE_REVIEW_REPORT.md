# CELL MESH Historical Review Notice

This file previously documented the scoring implementation that existed before
the bounded median-contrast migration. Those formulae, parameters, and output
fields have been removed to avoid presenting obsolete behavior as current
documentation.

The active implementation is documented in:

- `README.md`
- `docs/METHODS.md`
- `AVAILABILITY_IMPLEMENTATION.md`
- `IMPLEMENTATION_SUMMARY.md`

The current method filters eligible cell types using `min_cells`, preserves the
reaction-gene aggregation and P/C/E summation logic, and uses bounded
cell-type-median contrasts for sender and receiver scores. Missing exporter or
consumption priors are neutral. Raw C represents the expression-derived level
of metabolite-consuming enzyme complexes; its median contrast is not measured
extracellular clearance flux.
