# CELL MESH Code Review Status

The historical bounded-median-contrast implementation has been removed. The
active implementation uses sender-abundance-adjusted positive-reference
saturation scoring and is documented in:

- `README.md`
- `docs/METHODS.md`
- `AVAILABILITY_IMPLEMENTATION.md`
- `IMPLEMENTATION_SUMMARY.md`

Removed historical behavior includes enzyme/sensor numerical weights,
`min_cells`-based calculation filtering, signed P/C/E contrasts, unbounded exporter boosts,
duplicate score aliases, heuristic confidence tiers, and the ambiguous `fdr`
alias. E normalization now enters only through a bounded, non-increasing factor
whose default range is `[0.8, 1.0]`.

The current input and aggregation safeguards also require a non-empty reaction
identifier, reject missing/empty cell-type labels and duplicate gene names,
join sender and receiver evidence by canonicalized HMDB ID, and distinguish
sample-aware component medians from the median sample-level event score. Within
each canonical HMDB ID and P/C/E direction, identical reaction gene sets
contribute once and strict subsets are omitted, leaving all inclusion-maximal
gene sets for P/C/E aggregation.
