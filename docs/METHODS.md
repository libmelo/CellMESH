# CELL MESH methods

CELL MESH (Metabolite-mediated Event Scoring with Sensor Hierarchies) infers
candidate metabolite-mediated communication events from single-cell expression.
Each event is defined by a sender cell type, receiver cell type, canonical
HMDB identifier, and sensor gene. Sample-level records additionally use sample
identity. Metabolite names are display annotations; sensor type is a validated
annotation used for stratified FDR, not an additional event identity field.

## Prior databases

`load_cell_mesh_database()` independently selects the highest packaged
`Enzyme<version>.csv` and `Interaction<version>.csv`. The current defaults are
`Enzyme1.39.csv` and `Interaction1.41.csv`.

CSV and DataFrame inputs share schema normalization, also used when calling
`compute_metabolite_availability()` directly. Enzyme tables may use the
standard `role` schema or raw `Direction`/`direction` fields. Raw `export`,
`exporter`, and the legacy `transporter` value map to `role="export"` and feed
the E capacity; role/direction conflicts and missing required columns are
rejected. Existing provenance
metadata is preserved, and CSV identifier fields are read as text so standard
priors can be saved and loaded without losing reaction IDs or role records.

All simultaneous enzyme field aliases are compared before filtering any rows.
Equivalent values are merged, blank values are filled from non-empty aliases,
and conflicting values report the affected rows and columns. Direction aliases
are compared by role meaning even when `role` is absent. Gene aliases are
compared as sets of parsed symbols, preserving inline evidence from equivalent
lists. Consumed aliases are removed after merging. Exact duplicate column names
are rejected for both DataFrames and CSVs before CSV header renaming can hide them.

Prior CSVs also validate the width of every logical record before numerical or
schema parsing. Short/extra-wide records, malformed quoting, and NUL characters
raise file/line errors; missing cells must be represented by empty fields.
Quoted separators and newlines are preserved, and blank lines outside records
are ignored. An inferred row index cannot replace a missing field.

Multi-gene `gene` fields use the same parser for both enzyme schemas. Symbols
separated by `;`, `,`, or `|` are expanded before expression-gene matching,
including mixed separators and optional `[evidence]` annotations. Separators
within annotations are preserved. Each expanded row retains its reaction key
and existing provenance; inline evidence is used when no `evidence_level`
column was supplied. Empty entries are ignored. Expansion changes the table
representation only: repeated genes in the same reaction are still deduplicated
and contribute once. Observed and compiled permutation scoring both consume
normalized single-gene records instead of independently parsing the input.

The enzyme prior requires the `metabolite`, `hmdb_id`, `gene`, `role`, and
`reaction` columns. The name column may contain missing or blank values;
valid identifiers, genes, roles, and reaction IDs follow the checks above.
It uses `production`, `degradation`, and `export` roles,
mapped internally to P, C, and E. Genes within one reaction contribute equally; neither
the enzyme nor interaction database defines a numerical weight. Records without
a valid HMDB identifier are excluded, and sender and receiver evidence is joined
by a stripped, case-normalized HMDB identifier. Metabolite names may differ
between the two priors and do not participate in the join.

The reaction key is `canonical_hmdb_id + reaction + direction`. Exact symbols
are deduplicated within a reaction to define its gene set. Across reactions for
the same canonical HMDB ID and direction, complete gene sets are compared
independent of gene order. Identical sets contribute once. A reaction is omitted
when its gene set is a strict subset of another reaction gene set, leaving only
inclusion-maximal sets. Incomparable sets—including overlapping sets for which
neither contains the other—remain intact and both equal-weight geometric means
contribute. The retained reactions are summed into one HMDB-level P, C, or E
capacity. Reported reaction counts refer to these maximal contributing gene
sets. For display, the first non-empty Enzyme name for each HMDB ID is preferred,
followed by the first non-empty Interaction name, then the HMDB ID itself.
Names are chosen from normalized scoring priors. Standalone scoring uses the
available prior and the same ID fallback. Names never determine matching,
sample grouping, or which records survive aggregation.

Complete prior gene sets, including unmeasured genes, determine duplicate and
subset relationships in both public scoring entrypoints and compiled
permutations. Expression-gene matching happens only after reaction selection:
the geometric mean uses measured genes, retaining their zero expression values,
and a reaction with no measured genes receives zero capacity. Such a reaction
retains its prior identity in reaction counts and evidence-status reporting.
This zero is a numerical placeholder when no reaction genes were measured;
it does not count as evidence of measured zero production. A metabolite is
production-evaluable if any gene from its retained production reaction sets
exists in the expression matrix. Measured zero production remains evaluable.
The shared coverage rule, rather than a positive-P filter, selects metabolites
for observed scoring and every permutation in both inference modes.
Both paths apply that rule before calculating P/C/E references or normalized
scores. Production-ineligible rows do not enter reference calculation, while
evaluable zero rows and all observed cell types remain. Raw-expression and
computed-capacity checks still run before selection; only the unnecessary
normalization of excluded rows is skipped. Compiled event indices remain fixed.
Missing-gene filtering must not turn incomparable reactions into subsets.
The main pipeline requests complete enzyme records from runtime validation;
standalone `validate_priors()` callers retain the historical measured-gene
filter unless they pass `filter_enzyme_genes=False`.

Each scoring unit's `production_diagnostics` table includes all Enzyme-prior
metabolites. `metadata` keeps its existing alignment to the score matrices.
`production_evaluable` distinguishes eligible scores from missing production
evidence. `production_status` is `prior_missing` when no production relation
exists, `prior_gene_unavailable` when no production gene is measured,
`prior_no_expression` for evaluable all-zero capacities, or `supported` when
any cell type has positive capacity. The latter two states retain every
observed cell type's score, including zeros. Unavailable metabolites are absent
from the numerical score matrices and event tables; reindexing P to the
diagnostics table for display yields NA for those rows. Partial gene availability still uses
the measured genes without inserting missing genes as zeros.

Sensors are normalized to `Cell surface receptor`, `Transporter`, or
`Other receptor`. Evidence level, source, protein name, and references remain
provenance metadata and do not multiply the numerical score.
All coexisting Interaction field aliases are compared before filtering or
deduplication. HMDB aliases use stripped uppercase identifiers, sensor-type
aliases use the normalized category, and other fields use stripped,
case-sensitive text. Equivalent values merge, blanks are filled, and conflicts
report data-row positions, index labels, column names, and values. Duplicate
DataFrame column names and original CSV headers are rejected. Consumed aliases
are removed; raw `Annotation`/`annotation` text fills missing `evidence_level`
values, with explicit evidence taking precedence. Interaction genes are single
symbols and do not use the enzyme multi-gene parser. Full column mappings and
input conventions are listed in the [database directory guide](../cellmesh/data/README.md).
After HMDB canonicalization, database loading checks sensor-type conflicts
before deduplication and before expression-gene filtering, for CSV and DataFrame
inputs alike. Same-type duplicates keep the first metadata record; conflicting
types for one `canonical_hmdb_id + sensor_gene` pair raise an error. Runtime
prior validation reuses this check for its expression-compatible records.

## Expression aggregation and cell abundance

Optional receiver pseudobulk and expression-fraction summaries must cover every
cell type actually observed in the supplied AnnData and every measured sensor
gene in the prior. Coverage is checked against obs/var before computing a
positive reference; it is not inferred from the supplied summary. Missing or
extra groups and duplicate/empty normalized axes are rejected. Unused
categorical levels are excluded, while low-count and all-zero observed groups
remain. Each sample-aware scoring unit uses its own observed groups.
Relevant summary means must be finite and non-negative and fractions must be
finite and in [0, 1]. Supplied counts must be positive integer values equal to
actual obs counts; fractional values cannot be truncated. Valid summaries are
aligned by identifiers without modifying or recomputing them. These checks
cannot establish that a cache was calculated from the current X/layer; the
caller is responsible for that provenance. See the
[summary input contract](../README.md#optional-receiver-summaries).

Every measured gene used by enzyme or sensor scoring must have real, finite,
non-negative cell-level expression in the selected expression layer. The main
pipeline validates both sets before constructing any pseudobulk. Direct scoring
entrypoints also validate their relevant genes, including sensor calls with
precomputed summaries. Compiled permutation scoring checks its selected matrix
once during compilation; the check is not repeated for each permutation.
Expression plots validate the cells and genes they select. The shared check
does not alter values, validate unrelated genes/layers, densify sparse matrices,
or filter complete reaction definitions. Temporary arrays are bounded in size.

Checking only aggregated values is insufficient: `[-1, 3]` averages to `1`,
and a negative reaction can be hidden by a positive reaction in the summed
capacity. Such inputs must fail with the same rule whether or not permutations
are requested. Existing checks on receiver pseudobulks and final P/C/E
capacities remain in place to detect invalid supplied summaries or numerical
results after aggregation.

Finite input also does not guarantee finite floating-point intermediates. Both
scoring paths check relevant pseudobulk means, geometric reaction activities,
summed P/C/E capacities, positive references, and normalization denominators.
The same checked reaction and positive-reference implementations serve
observed and compiled scoring. Reaction activity uses the equivalent
`expm1(mean(log1p(x)))` expression, with one-gene reactions returning x directly;
this preserves tiny positive x that would disappear in `gmean(x+1)-1`.
A reference or denominator can overflow even when its input capacities
are finite; the denominator is checked before division can turn that failure
into a finite zero. Non-finite or negative results raise an error naming the affected stage.
Actual positive-to-zero underflow is reported without aborting, as requested
for the current analysis workflow. Relevant pseudobulk means, reactions,
abundance powers/products, positive-reference normalization, sender base and
exporter products, and event scores have explicit checks. Representable tiny
values are not thresholded. A run aggregates stage/occurrence counts across
observation and all permutation workers and stores `numerical_diagnostics` in
result parameters (and hence the parameters JSON export). One logging notice
describes the loss and continued use of rounded zeros. The same diagnostics
are available from standalone scoring and the violin functions; details are
in the [numerical input/output policy](../README.md#numerical-precision-and-underflow).
No permutation draws are discarded because of underflow. Continued inference
uses the rounded scores and can differ from higher-precision calculations:
zeros can change positive-reference membership, production/evidence states,
scores and p/FDR. A state associated with reported underflow must not be
interpreted as proof of absent biological expression. Negative/non-finite
values and overflow remain errors; structural sample missingness remains NA.

The existence of positive production is tested directly, without an unnecessary
potentially overflowing sum over cell types.

Expression-gene names, cell-type labels, and supplied sample labels share one
identifier normalization rule: convert to text and strip leading/trailing
whitespace, without changing case or internal whitespace. Actual missing values
and empty normalized text are rejected; literal text such as `"NA"` is retained.
Expression columns must have unique normalized gene names. Distinct observed
cell-type or sample labels that normalize to the same text are rejected rather
than merged; repeated labels and unused categorical levels are handled normally.
Normalization preserves row/column positions and does not modify the caller's
AnnData object.

Prior matching, pseudobulks, expression fractions, cell counts, sample masks,
coverage status, and compiled permutations all use these canonical identifiers.
Permutation compilation normalizes identifiers once before repeated shuffling;
raw labels must not be reintroduced for speed, as failed event-key matching can
turn null scores into zeros and change a correct p-value of 1 to
`1 / (n_perms + 1)`. Explicit cell-fraction indices and expression plotting use
the same identifiers. Reaction-expression plots retain the selected complete
reaction definitions, calculate activity from measured genes, and assign zero
when a reaction has no measured genes.

For each observed cell type `t` and gene `g`, CELL MESH calculates the per-cell
mean expression `mean_expression(g,t)` and expression fraction
`expression_fraction(g,t)`.

Expression fraction is the number of cells with expression strictly greater
than zero divided by the number of cells in the group. Both observed and
permuted scoring use this definition. Explicitly stored zeros in sparse
matrices do not count as expression, so storage format does not change the
expression fractions or the `min_expr_frac` gate.

For reaction `r`, equal-weight multi-gene activity is:

```text
reaction_activity(r,t)
    = geometric_mean({mean_expression(g,t) + 1 | g belongs to r}) - 1
```

Sender abundance is applied before P/C/E construction:

```text
abundance_weight(t) = cell_fraction(t) ** sender_abundance_exponent
adjusted_reaction_activity(r,t) = reaction_activity(r,t) * abundance_weight(t)
```

The default exponent is 1.0. Receiver abundance does not enter receiver scoring.
`min_cells` is QC-only: every observed cell type remains in pseudobulks,
fractions, references, scores, events, and permutations.

## Sender score

Adjusted reaction activities are summed into production P, consumption C, and
export E capacities. For each metabolite and direction X:

```text
X_ref(m) = mean({X(m,t) | X(m,t) > 0})              default
X_ref(m) = median({X(m,t) | X(m,t) > 0})            optional
X_score(m,t) = 0                                    if X(m,t) = 0
X_score(m,t) = X(m,t) / (X(m,t) + X_ref(m))         otherwise
```

The base sender score and bounded exporter modulation are:

```text
base_sender_score(m,t) = P_score(m,t)^2 / (P_score(m,t) + C_score(m,t))
E_factor(m,t) = (1 - export_weight) + export_weight * E_effective(m,t)
sender_score(m,t) = base_sender_score(m,t) * E_factor(m,t)
```

The base is zero when its denominator is zero. For supported exporter evidence,
`E_effective=E_score`; missing or gene-unavailable exporter priors use the fixed
neutral value `0.5`, and measured all-zero exporter capacity uses `0`. The
default `export_weight=0.2` restricts `E_factor` to `[0.8, 1.0]`, and zero
weight exactly recovers the base formula. Metadata keeps `consumption_status`
and `export_status` to distinguish these evidence states.

The implementation evaluates the base as `P_score * (P_score / denominator)`
instead of squaring P first, and event scores as
`sqrt(sender_score) * sqrt(receiver_score)`. These are mathematically equivalent
and avoid unnecessary underflow of intermediate products. Ordinary finite
results may differ only by floating-point rounding; the observation and null
paths share the stable helpers and retain the existing tail-tie rule.

## Receiver and event scores

For sensor gene `g`, receiver expression R is normalized across strictly
positive observed cell-type means:

```text
R_ref(g) = median({R(g,t) | R(g,t) > 0})            default
R_ref(g) = mean({R(g,t) | R(g,t) > 0})              optional
receiver_score(g,t) = 0                             if R(g,t) = 0
receiver_score(g,t) = R(g,t) / (R(g,t) + R_ref(g)) otherwise
```

`min_expr_frac` optionally gates receiver scores. The event score is:

```text
cell_mesh_score = sqrt(sender_score * receiver_score)
```

No heuristic confidence tier is assigned by the core algorithm.

## Sample-aware mode

Sender summaries group by HMDB ID; receiver summaries group by HMDB ID, sensor
gene, and receiver. Event statistics and sample completion use
`sender + receiver + hmdb_id + sensor_gene`, adding `sample` for sample-level
rows. Names and sensor types are attached separately, including on structurally
missing sample rows. The sender tables retain their `(metabolite, hmdb_id)`
display index for compatibility and export; its name level is not a grouping
key. Stored permutation-null rows use the four-field event identity index,
without names or sensor types.

`sample_aware` mode computes fractions, pseudobulks, references, and scores
independently within each sample. The aggregate `cell_mesh_score` and
`event_score_median` are the median of sample-level event scores. Descriptive
component summaries are reported separately as
`metabolite_availability_median`, `sensor_score_median`, and
`sensor_expr_frac_median`; their geometric mean is not the aggregate event
score. Top-level sender and receiver component summaries also use the
across-sample median. IQR, positive prevalence, and the number of co-observed
samples describe cross-sample support. A missing
sample-event value means the event was structurally uncomputable in that sample;
failing `min_cells` changes only its QC flag.

Evaluable zero-production samples retain zero event scores, endpoint cell
counts, and QC flags. Their zeros enter the event median and IQR, the
`n_samples_coobserved` count, and the denominator of `event_prevalence`.
Samples without an endpoint remain NA and do not enter that event's denominator.
For example, `[0.474342, 0, NA]` has median `0.237171`, two co-observed samples,
and positive prevalence `0.5`. Sender component medians follow the same rule.
The code never fills structural NA with zero to achieve this result.

Gene availability is assessed from normalized `adata.var_names`, shared by all
samples in the supplied matrix. Zeros introduced upstream to represent missing
measurements cannot be distinguished from measured zeros without additional
availability information; the scorer does not infer that information from zero
expression. An entirely unmeasurable enzyme prior still fails input validation.

## Permutation inference

The public entrypoint validates `random_state` as a non-negative Python/NumPy
integer before scoring, including zero-permutation and empty-event runs.
`allow_self` and `store_null_scores` require Python/NumPy booleans; strings such
as `"False"` cannot silently alter the event family or storage behavior.

Both modes use label permutation; when a sample key is supplied, labels are
shuffled within samples. Every observed analysis cell participates in the null.
Before permutation, CELL MESH compiles the validated reaction/sensor mappings,
restricts the numerical kernel to prior-relevant genes, and encodes the retained
observed event keys. Each permutation recomputes the same P/C/E and receiver
formulas from permuted cell labels, but does not rebuild public metadata or a
full sender-receiver event table. In sample-aware scoring, each null summary
includes evaluable zero samples and excludes samples without an endpoint,
using the same evaluability rule as observation. Missing event keys retain the
existing null-score fallback of zero.
For B permutations, the one-sided empirical p-value is:

```text
p = (1 + number_of_null_scores_at_least_observed) / (B + 1)
```

Observed and permuted pseudobulks share a float64 membership-matrix
aggregation kernel, summing expression before dividing by cell counts,
including for float32 input. Tail comparisons count
numerical ties using `null_score >= observed_score - tolerance`, where
`tolerance = 100 * float64_epsilon * abs(observed_score)`. This relative-only
tolerance follows [SciPy's permutation-test convention](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html) and avoids interpreting
roundoff differences as evidence of communication. Zero null scores still do
not count as ties with positive observed scores. Stored null scores retain
their original values; the tolerance affects only comparison counts.

Both modes return:

- `fdr_global`: Benjamini-Hochberg correction across all events.
- `fdr_sensor_type`: correction separately within each sensor type.

The ambiguous historical `fdr` alias is not returned.

Measured all-zero metabolites now yield zero-score events in both modes.
A zero observed event score has `perm_pvalue=1` when permutations are run;
`n_perms=0` still leaves p-values unavailable. Older results that dropped all-zero
P rows should be recomputed: sample summaries can change, and newly retained
zero events also change the FDR correction family. Unavailable production never
adds a hypothesis or a zero to the sample denominator.

Exceedance counts are accumulated online. In sample-aware mode the full
event-by-permutation null matrix is omitted by default and can be requested with
`store_null_scores=True`; this does not change p-values or FDR. Permutation
batches may be evaluated with `n_jobs > 1` (or `-1` for all CPUs). Label
assignments are generated in deterministic permutation-index order, so results
do not depend on worker completion order.

Each shuffle changes the grouping and therefore requires checks on newly
computed intermediates, even though raw expression values were validated once.
Computed sample-event values are checked before sample medians, separately from
deliberate structural missingness. Final permutation scores must be finite and
non-negative before they enter exceedance counts or stored null matrices. This
applies to serial and parallel execution and both null-storage settings. An
invalid score stops inference instead of being skipped or counted as a value
below the observation.
