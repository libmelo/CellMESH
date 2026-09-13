# Sender Score Implementation

This document describes sender calculations and their diagnostic outputs.
See [README.md](README.md#production-evaluability-and-zero-scores) for API usage
and [methods](docs/METHODS.md) for the complete scoring and inference rules.

Every observed cell type participates in pseudobulk construction, cell-fraction
adjustment, P/C/E references, scores, events, and permutations. `min_cells` is
used only to annotate cell-count QC.

Reaction activity is adjusted before P/C/E construction:

```text
measured_genes(r) = genes(r) intersect normalized adata.var_names
reaction_activity(r,t)
    = geometric_mean({expression(g,t) + 1 | g in measured_genes(r)}) - 1

adjusted_reaction_activity(r,t)
    = reaction_activity(r,t) * cell_fraction(t) ** sender_abundance_exponent
```

Here `expression(g,t)` is mean expression in cell type `t` from the selected
layer. The geometric mean retains measured zero-expression genes. A reaction
with no measured genes uses an internal zero capacity; that placeholder does
not establish measured zero production or make the metabolite evaluable.

Genes within a reaction contribute equally. The enzyme prior has no numerical
`weight`; legacy/custom enzyme `weight` columns are ignored and removed during
prior validation.
The `reaction` field is required and must be non-empty; rows without a known
reaction identity are rejected instead of being merged into an artificial
multi-gene complex.
Reaction grouping uses `canonical_hmdb_id + reaction + direction`; metabolite
name is display metadata. Exact gene symbols define each reaction's gene set.
Within the same HMDB ID and direction, complete gene sets are compared
independent of gene order. Identical sets contribute once, and strict subsets
are omitted in favor of their supersets. Incomparable sets, including overlaps
where neither contains the other, retain their full gene sets and all
contribute. These inclusion-maximal reaction gene sets are summed into one
HMDB-level capacity, and reaction-count metadata counts the contributing sets.

Reaction selection uses complete prior gene sets, including unmeasured genes.
Only after selecting the maximal sets does each reaction use its measured genes
for activity. Partial availability retains this existing calculation; missing
genes are not added as zeros or removed before comparing reaction sets.

## Production evaluability

P denotes production capacity; the statistical permutation p-value is
`perm_pvalue`. A metabolite's `production_evaluable` flag is true when at least
one gene in its retained production reactions exists in the normalized
expression gene columns. It does not require positive expression or complete
coverage of all production genes. Observed and compiled permutation scoring
share this rule; replacing it with a `P > 0` filter would discard valid zeros.
Both paths apply the mask before P/C/E reference normalization, after validating
raw expression and computed capacities. Excluded metabolites do not have their
references calculated; evaluable metabolites continue to use all observed cell
types in the scoring unit, including those with zero capacity.

`production_status` explains the evidence in each scoring unit:

| State | Production evidence | `production_evaluable` | Scoring behavior |
|---|---|---|---|
| `prior_missing` | No production relation | False | Exclude sender/event scores |
| `prior_gene_unavailable` | Production relations exist, but none of their genes are measured | False | Exclude sender/event scores |
| `prior_no_expression` | Some production genes are measured, but all calculated P capacities are zero | True | Keep zero sender/event scores |
| `supported` | Some production genes are measured and at least one cell type has positive P | True | Score every observed cell type, including its zeros |

One scoring unit is the pooled data or one sample in sample-aware mode. These
states summarize a metabolite across that unit's cell types; `supported` does
not imply positive P in every cell type or statistical significance. If no
enzyme-prior gene at all matches the matrix, `run_cell_mesh()` still raises an
error before analysis.

## Reference normalization and sender score

For each metabolite and each direction `X` in `P`, `C`, and `E`, the reference
is calculated from strictly positive adjusted capacities:

```text
X_ref(m) = mean({X(m,t) | X(m,t) > 0})              default
X_ref(m) = median({X(m,t) | X(m,t) > 0})            optional
X_score(m,t) = 0                                    if X(m,t) = 0
X_score(m,t) = X(m,t) / (X(m,t) + X_ref(m))         otherwise
```

Because both `X(m,t)` and `X_ref(m)` are strictly positive in the division
branch, their sum is positive and no epsilon parameter is required.
With no positive capacity, the reference remains NaN and the normalized scores
are zero. For an evaluable all-zero P row, `P_ref=NaN` therefore accompanies
valid zero `P_score` and sender scores. Relevant input expression must be real,
finite, and non-negative. Computed means, reactions, capacities, positive
references, and division denominators are checked too: finite inputs can still
overflow during arithmetic, and such failures raise an error rather than
becoming zero scores.

The base sender score and bounded exporter modulation are:

```text
base_sender_score(m,t) = P_score(m,t)^2 / (P_score(m,t) + C_score(m,t))
E_factor(m,t) = (1 - export_weight) + export_weight * E_effective(m,t)
sender_score(m,t) = base_sender_score(m,t) * E_factor(m,t)
```

The base result is zero when the denominator is zero. `E_effective` equals
`E_score` for supported exporter evidence, 0.5 for missing or gene-unavailable
exporter priors, and 0 for measured all-zero exporter capacity. The default
`export_weight=0.2` bounds `E_factor` to `[0.8, 1.0]`; a zero weight recovers
the base result exactly.

## Diagnostic outputs and sample summaries

Intermediates include raw `P`, `C`, and `E`, their positive-reference scores
and references, `base_availability`, `E_effective`, `E_factor`, pseudobulk
expression, expression fractions, cell counts, cell fractions, abundance
weights, and QC annotations. Obsolete signed contrasts and duplicate score
aliases are not returned.

Availability `metadata` has the same `(metabolite, hmdb_id)` index as P/C/E and
sender scores. It contains reaction counts, `consumption_status`, `export_status`,
`production_status`, and `production_evaluable` for the scored metabolites.
The separate `production_diagnostics` table uses the same index fields but
covers all Enzyme-prior metabolites, including ones excluded from scoring:

| Column | Meaning |
|---|---|
| `n_product_reactions` | Number of retained maximal production reaction gene sets, including unmeasured sets |
| `production_status` | One of the four production evidence states above |
| `production_evaluable` | Whether production has at least one measured gene |

The diagnostics do not include metabolites found only in Interaction and absent
from the Enzyme prior. Each sample's availability result contains its own table.
Unavailable production is excluded from the score matrices; to show it as NA
alongside measured zeros, reindex P for display:

```python
unit = res.availability_results  # pooled mode
# In sample_aware mode, select one sample instead:
# unit = res.availability_results["availability_by_sample"]["D1"]
diagnostics = unit["production_diagnostics"]
P_display = unit["P"].reindex(diagnostics.index)
```

In sample-aware mode, a computable zero event participates in medians, IQR,
`n_samples_coobserved`, cell-count/QC summaries, and the `event_prevalence`
denominator. A sample missing either endpoint cell type retains NA for that
event and is excluded from its denominator. For example, scores
`[0.474342, 0, NA]` yield median `0.237171`, two co-observed samples, and
positive prevalence `0.5`.

Both inference modes retain evaluable all-zero metabolites as zero-score events
when matching sensor evidence exists. With permutations, a zero observed event
score has `perm_pvalue=1`; with `n_perms=0`, p-values remain NA. Previous versions
dropped all-zero P rows. Rerun affected analyses because retaining these events
can change sample summaries and expands the FDR correction family in either mode.

Gene availability follows the shared columns in the supplied AnnData object.
If upstream processing filled unmeasured data with zeros, expression values
alone cannot recover that missingness; separate measurement-availability
information is required. Zero expression never determines gene availability.
