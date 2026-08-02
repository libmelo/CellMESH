# CELL MESH Code Review Status

The historical bounded-median-contrast implementation has been removed. The
active implementation uses sender-abundance-adjusted positive-reference
saturation scoring and is documented in:

- `README.md`
- `docs/METHODS.md`
- `AVAILABILITY_IMPLEMENTATION.md`
- `IMPLEMENTATION_SUMMARY.md`

Removed historical behavior includes enzyme/sensor numerical weights,
`min_cells`-based calculation filtering, signed P/C/E contrasts, exporter boosts,
duplicate score aliases, heuristic confidence tiers, and the ambiguous `fdr`
alias. E normalization remains available as support evidence and does not enter
the formal sender score.

The current input and aggregation safeguards also require a non-empty reaction
identifier, reject missing/empty cell-type labels and duplicate gene names,
join sender and receiver evidence by canonicalized HMDB ID, and distinguish
sample-aware component medians from the median sample-level event score.
