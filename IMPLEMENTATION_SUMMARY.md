# Implementation Summary

- Every observed cell type participates in pseudobulk, reference, score, event,
  and permutation calculations; `min_cells` produces QC annotations only.
- Equal-weight multi-gene reaction activity uses
  `geometric_mean(expression + 1) - 1`.
- `reaction` is required and non-empty; unknown reactions are never collapsed
  into an inferred multi-gene complex.
- Reaction grouping uses canonical HMDB ID, reaction, and direction. For one
  HMDB ID/direction, identical gene sets contribute once and strict subsets are
  omitted; all inclusion-maximal sets contribute. Metabolite name is
  display-only.
- Sender reaction activity is multiplied by
  `cell_fraction ** sender_abundance_exponent` before P/C/E construction.
- P/C/E use positive-reference saturation `X / (X + X_ref)`, with the positive
  mean as default and positive median as an option.
- The base sender score is `P_score^2 / (P_score + C_score)`.
- E applies the bounded factor `(1-export_weight) + export_weight*E_effective`.
  The default weight is 0.2; missing/gene-unavailable evidence uses fixed
  `E_effective=0.5`, while measured all-zero exporter capacity uses 0.
- Receiver score is `R / (R + R_ref)` and does not use receiver abundance.
- Pooled and sample-level event score is `sqrt(sender_score * receiver_score)`.
- Sample-aware mode scores each sample independently and aggregates event scores
  using median, IQR, and prevalence summaries. Aggregate component fields have
  explicit `_median` names and top-level component tables use sample medians.
- Sender and receiver priors join on canonicalized HMDB ID, not metabolite name.
- Cell-type labels must be non-missing/non-empty and expression gene names must
  be unique.
- Both modes report `fdr_global` and `fdr_sensor_type`; no ambiguous `fdr` alias
  or heuristic confidence tier is produced.
- Availability metadata retains only reaction counts plus
  `consumption_status` and `export_status`.
