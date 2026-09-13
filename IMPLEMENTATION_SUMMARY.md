# Implementation Summary

See [README.md](README.md) for usage and output fields,
[methods](docs/METHODS.md) for calculation rules, and the
[review report](CELL_MESH_CODE_REVIEW_REPORT.md) for resolved findings and validation.

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
- Reaction selection uses complete prior gene sets before matching expression
  columns. Each retained reaction then uses its measured genes, including their
  zeros; missing genes are not inserted as zeros into the geometric mean.
- Sender reaction activity is multiplied by
  `cell_fraction ** sender_abundance_exponent` before P/C/E construction.
- P/C/E use positive-reference saturation `X / (X + X_ref)`, with the positive
  mean as default and positive median as an option.
- The base sender score is `P_score^2 / (P_score + C_score)`.
- E applies the bounded factor `(1-export_weight) + export_weight*E_effective`.
  The default weight is 0.2; missing/gene-unavailable evidence uses fixed
  `E_effective=0.5`, while measured all-zero exporter capacity uses 0.
- Receiver score is `R / (R + R_ref)` and does not use receiver abundance.
- Supplied receiver summaries must cover every actually observed cell type and
  every measured prior sensor gene. Normalized axes must be unique; missing and
  extra groups are rejected before reference calculation. Relevant means and
  fractions are checked for valid ranges, and counts must equal actual obs
  counts before integer conversion. Valid summaries are aligned without
  mutation or recomputation; callers remain responsible for cache provenance.
- Calculation booleans accept Python/NumPy booleans only. `random_state` accepts
  non-negative Python/NumPy integers only, including when no permutations run
  or the event table is empty. NumPy parameters are stored as Python scalars.
- Pooled and sample-level event score is `sqrt(sender_score * receiver_score)`.
- Sample-aware mode scores each sample independently and aggregates event scores
  using median, IQR, and prevalence summaries. Aggregate component fields have
  explicit `_median` names and top-level component tables use sample medians.
- Production is evaluable when at least one gene in the retained production
  reactions occurs in the expression matrix. This eligibility rule is shared by
  observed and permuted scoring and does not depend on positive P values.
  Both paths apply this mask before P/C/E reference normalization, retaining
  all observed cell types and the upstream expression/capacity checks.
  Evaluable zero scores remain in both inference modes. In sample-aware mode,
  they enter medians, IQR, co-observed sample counts, and prevalence denominators;
  an absent sender or receiver cell type leaves the sample event NA.
- Sender and receiver priors join on canonicalized HMDB ID, not metabolite name.
- Sample aggregation and stored null indices use
  `sender + receiver + hmdb_id + sensor_gene`, with `sample` added for sample
  records. Sender summaries group by HMDB ID; receiver summaries group by
  HMDB ID, sensor gene, and receiver. Names and validated sensor types are
  attached separately. Display names prefer a non-empty Enzyme name, then
  Interaction, then the ID; name-only sender matching is rejected.
- CSV and DataFrame priors share field-alias normalization and conflict checks;
  duplicate headers are rejected. Prior CSV logical-record widths and quoting
  are checked before pandas can infer an index or pad incomplete records.
  Enzyme multi-gene fields are expanded before
  expression matching. Interaction genes use single symbols. Accepted columns
  and values are listed in the [database guide](cellmesh/data/README.md).
- 10X loading validates raw feature/barcode text and matrix dimensions before
  attaching labels. Legacy/modern layouts support plain/gzip files, prefixes,
  feature filtering, and explicit symbol/ID selection. Duplicate selected
  identifiers raise errors; automatic suffixes are never added. Matrix caching
  retains raw annotation validation on every load.
- Gene, cell-type, and sample identifiers are normalized consistently. Missing
  or empty labels, duplicate normalized gene names, and distinct group labels
  that normalize to the same text are rejected.
- Text expression import checks raw headers before automatic duplicate renaming.
  Axis names, metadata, and MTX name lists retain literal text, including leading
  zeros and `NA`; expression values remain numeric. Metadata fields default to
  text so custom grouping columns retain their identity. See the
  [AnnData reading guide](docs/ANNDATA_README.md) for conversion and parser rules.
- Relevant cell-level expression must be real, finite, and non-negative in the
  selected layer. Shared numerical checks also reject invalid computed means,
  reaction activities, capacities, positive references, denominators, and scores.
- Observed and compiled group means use the same float64 accumulation. Expression
  fractions count values strictly greater than zero, including for sparse inputs.
- Both modes report `fdr_global` and `fdr_sensor_type`; no ambiguous `fdr` alias
  or heuristic confidence tier is produced.
- Permutation inference compiles static priors and relevant expression columns
  once, scores only observed event keys, and accumulates exceedance counts
  online. Deterministic worker batches are supported; sample-aware null matrices
  are opt-in via `store_null_scores=True`.
- Permutation tail counts include numerical ties using a float64 relative
  tolerance. When permutations run, zero-score events have `perm_pvalue=1`;
  with `n_perms=0`, p-values remain NA.
- Availability `metadata` stays aligned to the score matrices and contains
  reaction counts, `consumption_status`, `export_status`, `production_status`,
  and the Boolean `production_evaluable`. The production states are
  `prior_missing`, `prior_gene_unavailable`, `prior_no_expression`, and
  `supported`; the latter two are evaluable.
- A separate `production_diagnostics` table, indexed by `(metabolite, hmdb_id)`,
  contains `n_product_reactions`, `production_status`, and `production_evaluable`
  for all Enzyme-prior metabolites, including those excluded from scoring. Each
  sample has its own diagnostics in sample-aware mode. See
  [production evaluability and zero scores](README.md#production-evaluability-and-zero-scores)
  for state definitions, access examples, and the limits of gene availability.
- CSV export preserves identifier fields once, rejects conflicting headers
  before writing, and keeps measured zeros distinct from missing values.

Analyses affected by earlier positive-P filtering should be rerun. Retaining
evaluable zero events can change sample summaries and expands the FDR correction
family in either inference mode. The review report also records the effects of
earlier input, reaction-selection, and permutation fixes.
