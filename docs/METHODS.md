# CELL MESH methods

CELL MESH (Metabolite-mediated Event Scoring with Sensor Hierarchies) infers
candidate metabolite-mediated communication events from single-cell expression.
Each event is defined by a sender cell type, receiver cell type, metabolite,
HMDB identifier, and sensor gene.

## Prior databases

`load_cell_mesh_database()` independently selects the highest packaged
`Enzyme<version>.csv` and `Interaction<version>.csv`. The current defaults are
`Enzyme1.39.csv` and `Interaction1.41.csv`.

The enzyme prior requires non-empty `metabolite`, `hmdb_id`, `gene`, `role`, and
`reaction` fields. It uses `production`, `degradation`, and `export` roles,
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
sets. The first enzyme-prior metabolite name for each HMDB ID is retained only
for display.

Sensors are normalized to `Cell surface receptor`, `Transporter`, or
`Other receptor`. Evidence level, source, protein name, and references remain
provenance metadata and do not multiply the numerical score.
After HMDB canonicalization and expression-gene filtering, runtime validation
enforces one sensor record per `canonical_hmdb_id + sensor_gene`. Same-type
duplicates keep the first metadata record; conflicting sensor types for one
pair are rejected before communication events are constructed.

## Expression aggregation and cell abundance

For each observed cell type `t` and gene `g`, CELL MESH calculates the per-cell
mean expression `mean_expression(g,t)` and expression fraction
`expression_fraction(g,t)`.

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

## Permutation inference

Both modes use label permutation; when a sample key is supplied, labels are
shuffled within samples. Every observed analysis cell participates in the null.
Before permutation, CELL MESH compiles the validated reaction/sensor mappings,
restricts the numerical kernel to prior-relevant genes, and encodes the retained
observed event keys. Each permutation recomputes the same P/C/E and receiver
formulas from permuted cell labels, but does not rebuild public metadata or a
full sender-receiver event table. Missing event keys retain null score zero.
For B permutations, the one-sided empirical p-value is:

```text
p = (1 + number_of_null_scores_at_least_observed) / (B + 1)
```

Both modes return:

- `fdr_global`: Benjamini-Hochberg correction across all events.
- `fdr_sensor_type`: correction separately within each sensor type.

The ambiguous historical `fdr` alias is not returned.

Exceedance counts are accumulated online. In sample-aware mode the full
event-by-permutation null matrix is omitted by default and can be requested with
`store_null_scores=True`; this does not change p-values or FDR. Permutation
batches may be evaluated with `n_jobs > 1` (or `-1` for all CPUs). Label
assignments are generated in deterministic permutation-index order, so results
do not depend on worker completion order.
