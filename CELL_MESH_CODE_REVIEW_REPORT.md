# CELL MESH Code Review Status

**Latest review (2026-09-17):** A01–A11 and B01–B06 remain resolved. B07 was
explicitly set aside by the user and is not a pending revision. The subsequent
[regular-input review](docs/REVIEW_REGULAR_INPUTS_2026-09-17.md) confirmed no new
defect: 48 numerical configurations, 17 plotting checks, and 220 selected
regression tests passed. No production code changed during that review. The
dated entries below retain their historical findings and implementation records.

The fixes in the resolved table below remain implemented. The systematic audit
conducted on 2026-09-12–2026-09-13 confirmed **11 grouped open revision items** covering
sample aggregation, input parsing, examples, numerical stability, optional
inputs, plotting, export, format support, and provenance. The complete current
inventory, causes, reproductions, and proposed acceptance criteria are in
[the systematic audit](docs/SYSTEMATIC_AUDIT_2026-09-12.md), with
[machine-readable evidence](docs/audit_2026_09_12_evidence.json).

The audit changed documentation only. Its baseline was 1019 passed, 1 skipped,
and 213 warnings; 32 additional configurations compared 6120 observed/full-
recomputation and compiled-permutation event values with no differences. The
two previously recorded open dotplot display issues below are included as A07
and A08 in the consolidated inventory.

A01 was fixed on 2026-09-13.
Sample aggregation and stored null indices now use HMDB ID, sensor gene,
sender, and receiver, with sample added for sample-level records. Names and
sensor types are attached separately; the name-only matching fallback is
removed. Display names prefer a non-empty Enzyme name, then Interaction, then
the ID. The 34 new regression cases and updated receptor-context tests passed;
the final full suite was **1055 passed, 1 skipped, 213 warnings in 135.10 s**.
The original missing-name example now retains 8 sample events, 4 aggregate
events, 1 sender row, and 2 receiver rows under all four name conditions.

A02 was fixed on 2026-09-13, together with the 10X format-support subitem of
A10. Prior CSV files now undergo
logical-record width and quoting validation before pandas parsing, retaining
valid multiline fields, explicit empty fields, and compressed inputs. The 10X
reader validates original selected gene identifiers and barcodes before any
automatic renaming can occur, supports explicit gene-symbol/gene-ID selection,
and handles legacy/modern plain or gzip layouts, prefixes, and mixed compression.
Matrix caching still uses Scanpy; current raw annotations are checked on every
call. The 92 new regression cases cover malformed inputs and scoring/permutation
equivalence in both inference modes. The final full suite was **1147 passed,
1 skipped, 213 warnings in 138.50 s**, with the same skip and AnnData warnings
as the A01 baseline. A10 remains open for Loom dependencies/imports and backed
AnnData scoring; the audit evidence JSON remains the original pre-fix snapshot.

A05 was fixed on 2026-09-13.
Receiver summaries are checked against the current AnnData's actually observed
cell types and measured prior sensor genes. Normalized axes must be unique,
relevant means and fractions must have valid ranges, and supplied counts must
equal actual obs counts before integer conversion. Valid cached expression
summaries are aligned without mutation or recomputation; cache provenance
remains the caller's responsibility. Calculation booleans accept Python/NumPy
booleans only; the non-negative integer seed is validated even for zero
permutations and empty events. The 216 new regression cases cover cache
integrity, sample-local coverage, zero expression, invalid parameters and
seeded serial/parallel equivalence. The final full suite was **1363 passed,
1 skipped, 213 warnings in 147.98 s**, retaining the original skip and AnnData
warnings. A09's exporter-level strict JSON check remains open.

A04 was fixed on 2026-09-13, leaving 7 grouped items open at that stage.
Reaction activity uses shared log1p/expm1 arithmetic, with exact one-gene
evaluation; sender base divides before multiplying and event scores multiply
square roots. Per the user's explicit decision, actual underflow is reported
and computation continues with rounded zeros, without discarding permutation
draws. One aggregated logging notice and exportable `numerical_diagnostics`
record stage/occurrence counts, including threaded work. These zeros can affect
references, evidence states, scores and p/FDR and do not prove absent expression.
Other numeric validity and overflow errors remain in place. Both violin
functions use exact constant detection and scale-stable density estimation,
falling back to original-value points with returned diagnostics if KDE fails.
The 67 new regression cases and visual checks passed. Full validation was
**1430 passed, 1 skipped, 213 warnings in 149.90 s**; the skip and AnnData
warnings are unchanged from A05. The historical audit evidence JSON is unchanged.

A06 was fixed on 2026-09-14, leaving 6 grouped items open at that stage
(A03 and A07–A11, with A10 limited to Loom/backed AnnData). All six plotting
entrances now use numeric validation before the operations consuming those
values, consistent selector normalization, and explicit missing-score outputs.
Native scores/probabilities/fractions are bounded in [0, 1]; custom scores and
raw expression/capacities are finite and non-negative without an upper bound.
Invalid text, infinity, complex and boolean values cannot masquerade as NA;
native complex columns also fail without a preceding cast warning. Reaction
lists are normalized for actual calculation as well as duplicate comparison.
Display ranges are finite and ordered. The A04 underflow policy is unchanged.
The 288 new regression cases passed; final full validation was **1718 passed,
1 skipped, 213 warnings in 153.03 s**, retaining the existing AnnData
warnings and skip. README and public docstrings describe the input/output
contract; A07/A08's probability-zero display and layout work remains separate.

A07 was subsequently fixed on 2026-09-14, leaving 5 grouped items open at that stage (A03 and A08–A11). Zero probabilities have an explicit maximum-area display
cap; positive probabilities use the lower 80% of the area interval when zeros
are present. Without zeros, the previous size mapping is unchanged. Both FDR
and p-value legends label actual probabilities, including the smallest positive
float, without inventing a floor for zero. Fixed-size plots explicitly disable
size encoding, and NA keeps its unavailable marker. Returned `dot_sizes` and
`size_encoding` make the display mapping inspectable without changing the data.
The 55 new cases and rendered checks passed; full validation was **1773 passed,
1 skipped, 213 warnings in 155.34 s**. A08's layout work remains open.

A08 was fixed on 2026-09-14, leaving 4 grouped items open at that stage
(A03 and A09–A11). Dotplot legends and colorbars use separate regions sized from
rendered text, marker and axis extents, independent of missingness. Automatic
figures expand to fit. External axes fit all decorations inside their existing
rectangle, preserve other subplots and the caller's suptitle, and report required
space if undersized. Active automatic layout engines are rejected before axes
mutation; callers first draw their overall layout, then set the engine to
`none`. The function does not silently change the caller's layout engine.
The 74 new layout cases cover missingness/zero/fixed-size branches, both
statistics, long labels, large markers, multiple senders, external subplot
ownership and PNG/SVG/PDF at two DPIs. Rendered examples were also inspected.
Full validation was **1847 passed, 1 skipped, 213 warnings in 176.43 s**;
the original skip and AnnData warnings remain unchanged.

A09 was fixed on 2026-09-16, leaving 3 grouped items open at that stage
(A03, A10 and A11). Exports stage all CSV/strict JSON contents before updating
reserved result filenames, remove obsolete optional standard files, and publish
a manifest last. Ordinary update failures attempt rollback using backups;
prevalidation/staging failures preserve previous exports. Cleanup never scans
arbitrary prefix matches or trusts paths from an old manifest. Exact standard
names are reserved even if hand-created; symlink/directory targets are rejected.
This does not guarantee crash-safe multi-file atomicity or concurrent writes.
Nested non-finite JSON parameters are rejected while valid NumPy scalars remain
supported. Empty pooled events retain inference_mode and skipped permutation
paths carry explicit completion/storage metadata. CSV NA/zero and identifier
semantics are unchanged. The 49 new cases passed; final full validation was
**1896 passed, 1 skipped, 213 warnings in 184.24 s**, retaining the
original skip and AnnData warnings.

A10's remaining Loom/backed work was fixed on 2026-09-16, leaving 2 grouped
revision items open at that stage (A03 and A11). Loom declares its optional dependency and uses
the current AnnData reader when available, with accurate missing-dependency
hints and unchanged file/internal-dependency errors. Shared expression slicing
sorts disk indices and restores requested row/gene order, preserving sparse
storage. It covers validation, summaries, compiled permutations and violins,
including reordered backed views. Sample-aware materializes only the current
sample's chosen matrix plus obs/var. Files remain unchanged and caller-owned.
Documentation states the memory boundary: this is not a fully out-of-core
algorithm and AnnData may load layers eagerly. View mapping is isolated and
regression-tested against the installed AnnData version.
The 96 backed and 8 Loom cases passed; temporary Loom dependencies under `/tmp`
were enabled for full validation: **2000 passed, 1 skipped, 213 warnings
in 201.20 s**. No new Loom tests were skipped; the original skip and warnings
are unchanged. A10's previously completed 10X changes remain covered.

A11 was fixed on 2026-09-16, leaving only A03 open at that stage. Missing direction-format
Enzyme sources are no longer assigned an invented test-library name. Row-level
sources remain distinct from loader provenance in table attrs and exported
`parameters.prior_inputs`. File records carry input kind, filename/path,
filename-derived version, byte SHA-256 and raw/normalized row counts. DataFrames
have no claimed file identity; stale incoming file attrs are not inherited.
The coverage CLI and hash-identified snapshot reproduce 76/1095 (6.9406%) for all
Enzyme HMDB IDs and 56/346 (16.1850%) for the intersection with Interaction.
Neither statistic filters by expression. The configurable export weight 0.2 and
all scoring formulas remain unchanged; the unsupported 17.1% justification was
removed. The 32 new cases passed; full validation was **2032 passed,
1 skipped, 213 warnings in 210.22 s**, retaining the existing skip and
AnnData warnings. Database contents remain local and ignored by Git.

A03 was fixed on 2026-09-17; **all 11 groups in this audit are now closed**.
The HNSC examples retain complete enzyme reactions for selected HMDB IDs, and
the comprehensive trace explicitly preserves these priors and asserts agreement
with the main workflow. Empty thresholded networks report the empty state without
relaxing thresholds. The toy example checks current production metadata and
explains bounded exporter modulation. Old notebook outputs and execution metadata
were cleared. Eleven new regression cases cover actual example preparation,
both inference modes, empty/nonempty network selections, and all four notebooks
in fresh Python processes with figure rendering (not browser-level Jupyter testing).
Full validation: **2043 passed, 1 skipped, 213 warnings in 225.66s**, with the existing skip and AnnData warnings
unchanged. No scoring formula or database contents were changed for A03.

The supplementary review on 2026-09-17 identified B01–B03 beyond the original
A01–A11 scope. These have now been revised: duplicate sparse coordinates are
promoted before reduction, external plotting QC values are validated and
normalized, and external cell fractions are checked before float conversion.
The implementation adds 282 regression cases. See the
[supplementary review and validation record](docs/REVIEW_2026-09-17.md) for the
pre-fix evidence, affected inputs and post-fix results.

B04–B06 from the follow-up review have also been revised: decoded input text is
checked for NUL before parsing, nullable numeric expression columns are converted
before AnnData construction, and validated float16 plotting columns are promoted
before sorting. This adds 168 regression cases. See the
[follow-up review and implementation record](docs/REVIEW_POST_B_2026-09-17.md).

The historical bounded-median-contrast implementation has been removed. The
active implementation uses sender-abundance-adjusted positive-reference
saturation scoring and is documented in:

- [README.md](README.md): usage, input columns, and result fields.
- [docs/METHODS.md](docs/METHODS.md): calculation and inference rules.
- [AVAILABILITY_IMPLEMENTATION.md](AVAILABILITY_IMPLEMENTATION.md): sender
  calculations, production states, and diagnostic outputs.
- [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md): implementation overview.
- [cellmesh/data/README.md](cellmesh/data/README.md): accepted Enzyme and
  Interaction columns, aliases, values, and input conventions.

## Current scoring model

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
gene sets for P/C/E aggregation. These comparisons use the complete reaction
definitions before matching expression genes; each retained reaction then
calculates activity from its measured genes.

Permutation inference now uses a compiled null-only kernel: validated priors,
maximal reaction sets, relevant expression columns, group structure, and event
indices are prepared once. Permutations recompute the same numerical scores but
skip public intermediates and full event-table construction. Exceedance counts
are accumulated online, deterministic worker batches are supported, and the
sample-aware full null matrix is opt-in.

## Resolved findings

| Area | Current behavior | Regression coverage |
|---|---|---|
| A01: scoring and sample identity | Aggregate and align events by HMDB/gene and cell-type/sample context. Missing names and aliases cannot remove or split events; names and types remain annotations. | [Scoring identity](tests/test_scoring_event_identity.py), [receptor context](tests/test_receptor_context.py) |
| A02: original prior and 10X input integrity; A10: 10X formats | Reject malformed prior CSV logical records before implicit index shifts or padding. Validate original selected 10X gene and barcode labels without automatic suffixes; support legacy/modern plain and gzip layouts. | [Prior CSV structure](tests/test_prior_csv_structure.py), [10X import and scoring](tests/test_io_10x.py), [reading guide](docs/ANNDATA_README.md) |
| A05: receiver caches and calculation parameters | Validate cache coverage against current observations, relevant numerical ranges and actual integer cell counts. Reject ambiguous booleans and invalid seeds before zero-permutation or empty-event shortcuts. | [Scoring input contracts](tests/test_scoring_input_contracts.py), [summary contract](README.md#optional-receiver-summaries) |
| A04: tiny positive values and violin variability | Share stable reaction/base/event formulas. Report actual underflow and continue with exportable diagnostics, as requested. Preserve real violin variation; use original points if density estimation fails. | [Numerical stability](tests/test_numeric_stability.py), [overflow checks](tests/test_numeric_overflow.py), [precision policy](README.md#numerical-precision-and-underflow) |
| A06: plotting input contracts | Validate numeric views before thresholds/ranking/deduplication, preserve genuine NA, return missing-score exclusions, normalize selectors and actual reaction lists, and check finite ordered display ranges. | [Plotting input contracts](tests/test_plotting_input_contracts.py), [visualization contract](README.md#visualization) |
| A07: zero-probability dot sizes and legends | Separate zero display caps from positive log-scaled areas; show actual probability keys and preserve unavailable/fixed-size semantics. | [Zero significance](tests/test_dotplot_zero_significance.py), [visualization contract](README.md#visualization) |
| A08: dotplot layout | Allocate separate legend/colorbar regions from rendered dimensions. Expand automatic figures, preserve caller subplot ownership, and report insufficient external space or unfinished automatic layout. | [Dotplot layout](tests/test_dotplot_layout.py), [usage constraints](README.md#visualization) |
| A09: export consistency and empty results | Stage all files, strictly validate JSON, manage only reserved same-prefix filenames, publish a manifest last, and keep mode-specific empty schemas and permutation metadata consistent. | [Export consistency](tests/test_export_consistency.py), [export contract](README.md#export-results) |
| A10: Loom and backed AnnData | Declare Loom dependency and current reader; preserve file errors. Share disk-safe ordered slicing and materialize current sample matrices while preserving file ownership and sparse structure. | [Backed inputs](tests/test_backed_inputs.py), [Loom](tests/test_loom_reader.py), [reading guide](docs/ANNDATA_README.md) |
| A11: prior provenance and coverage | Preserve row sources without invention, separately export actual input identities, and reproduce explicitly defined coverage from hash-identified files without changing model weights. | [Provenance tests](tests/test_prior_provenance.py), [coverage snapshot](docs/PRIOR_COVERAGE_2026-09-16.json) |
| Observed/permutation precision | Both paths use the same float64 group-mean accumulation; tail counts include roundoff ties. Constant expression retains p-values of 1. | [Precision](tests/test_permutation_precision.py), [compiled scoring](tests/test_compiled_permutation.py) |
| Sparse expression fractions | Count expression values strictly greater than zero. Stored sparse zeros do not count as expressing cells. | [Expression fractions](tests/test_expression_fraction.py) |
| Result export | Preserve named index identifiers once, reject ambiguous headers before writing, and keep measured zeros distinct from missing values. | [CSV export](tests/test_result_export.py) |
| Prior schemas and aliases | CSV and DataFrame inputs share normalization. Equivalent aliases merge, blanks are filled, conflicts are reported, and duplicate original headers are rejected. Direction/role mappings and HMDB joins are consistent. | [Prior consistency](tests/test_prior_input_consistency.py), [Enzyme aliases](tests/test_enzyme_aliases.py), [Interaction aliases](tests/test_interaction_aliases.py) |
| Multi-gene Enzyme reactions | Expand gene lists before expression matching, preserve provenance, select maximal sets from complete reactions, then calculate from measured genes. Interaction genes remain single symbols. | [Multi-gene inputs](tests/test_enzyme_multigene.py), [reaction consistency](tests/test_enzyme_reaction_consistency.py); details in the [Enzyme review](docs/ENZYME_INPUT_REVIEW.md) |
| Expression and group identifiers | Normalize gene, cell-type, and sample labels consistently; reject missing labels, duplicate normalized genes, and group-label collisions. Reaction plots follow the measured-gene rule. | [Identifier normalization](tests/test_identifier_normalization.py) |
| Text expression import | Reject duplicate raw headers before pandas renaming. Preserve literal axis names, all metadata fields, and MTX names; keep expression values numeric. Leading-zero sample labels remain separate groups. | [Text import](tests/test_io_identifiers.py), [reading guide](docs/ANNDATA_README.md) |
| Numeric validity | Validate relevant cell-level expression in the selected layer and reject non-finite or negative numerical results before normalization or permutation tail counting. | [Expression validation](tests/test_expression_validation.py), [numerical overflow](tests/test_numeric_overflow.py) |
| Production zeros and missingness | Keep evaluable zero scores, exclude unavailable production, retain diagnostic reasons, and preserve sample NA when an endpoint is absent. Observed and compiled scoring share eligibility. | [Production evaluability](tests/test_production_evaluability.py), [CSV export](tests/test_result_export.py) |
| Permutation normalization range | Apply production eligibility before P/C/E reference normalization, matching observation. Excluded rows skip reference calculation; raw/capacity checks, evaluable zeros, and event indices are preserved. | [Permutation evaluability](tests/test_permutation_evaluability.py) |
| Dotplot event identity | Use `(hmdb_id, sensor_gene)` for ranking, selection, pair-level deduplication, and row coordinates. Names only label rows; all labels include HMDB IDs. | [Dotplot identity](tests/test_dotplot_event_identity.py), [visual identity](tests/test_visual_identity.py) |
| Matrix Market input | Accept sparse coordinate and dense array storage, including gzip files. Preserve values, matrix dimensionality, gene-by-cell transposition, and identifier checks; sparse inputs stay sparse. | [MTX formats](tests/test_io_mtx.py), [text identifiers](tests/test_io_identifiers.py) |
| Receptor plot context | An explicit normalized HMDB ID and receptor gene identify the receiver context across different Enzyme/Interaction names. Names only control display; original receiver metadata remains intact. | [Receptor context](tests/test_receptor_context.py), [plotting](tests/test_plotting.py) |
| Shared visualization identity | All six plotting functions exclude names from matching and identity. Same-ID aliases count once; distinct IDs/genes remain separate. Target plots require IDs and reject conflicting records. | [Visual identity](tests/test_visual_identity.py), [usage and migration](README.md#visualization) |
| Unavailable dotplot significance | Missing significance uses fixed-size diamonds without fabricating probabilities or hiding scored events. All-missing data has a categorical legend; numeric sizes use only available values. | [Dotplot significance](tests/test_dotplot_significance.py), [plotting](tests/test_plotting.py) |

## Production diagnostics and zero-score behavior

P is production capacity; `perm_pvalue` is the statistical permutation p-value.
`production_evaluable` depends on whether at least one gene in the retained
production reactions occurs in the expression matrix, independently of its
expression magnitude. Partial gene availability uses the measured genes under
the existing reaction formula; missing genes do not enter the geometric mean
as zeros.

`production_status` is `prior_missing` for no production relation,
`prior_gene_unavailable` for entirely unmeasured production genes,
`prior_no_expression` for evaluable all-zero P, or `supported` when at least one
cell type has positive P. The latter two states retain scores, including zeros.
States summarize one pooled scoring unit or one sample, rather than each
individual cell type.

`metadata` remains aligned to the score matrices and contains reaction counts,
`consumption_status`, `export_status`, `production_status`, and
`production_evaluable`. The separate `production_diagnostics` table is indexed
by `(metabolite, hmdb_id)` and contains `n_product_reactions`,
`production_status`, and `production_evaluable` for all Enzyme-prior metabolites,
including those excluded from scoring. Sample-aware results hold one such table
per sample. Unavailable production has no sender/event score; reindexing P to
the diagnostics for display yields NA for those rows. If the entire Enzyme prior
has no measured genes, the main analysis still raises an error.

Computable zero samples enter sample medians, IQR, co-observed counts, QC
summaries, and prevalence denominators. Samples missing a sender or receiver
cell type remain NA for the affected events. Both inference modes retain
evaluable all-zero events when matching sensor evidence exists; with permutations,
their `perm_pvalue` is 1, and without permutations it remains NA. See the
[production guide](README.md#production-evaluability-and-zero-scores) for examples.

Availability is inferred from the shared gene columns in the supplied AnnData.
Upstream zero-filled missing measurements require separate availability
information; the expression matrix alone cannot identify them as missing.

## Effects on existing analyses

- Rerun analyses affected by positive reactions or intermediate products being
  rounded to zero. Stable formulas preserve representable positive results.
  Review `numerical_diagnostics` when present: continued underflow runs retain
  floating-point zeros and may have affected references, states and inference.
  Replot tiny nonconstant single-cell distributions previously shown as a line.
- Rerun custom analyses whose supplied receiver summaries omitted observed
  groups or measured sensor genes, or contained invalid fractions/counts.
  Use complete summaries from the same X/layer or omit them for automatic
  computation. Replace string booleans such as `allow_self="False"` with actual
  booleans and rerun if they previously expanded the event/FDR family.
- Rerun analyses affected by mismatched observed/null accumulation, stored sparse
  zeros, identifier matching, or reaction selection after missing-gene filtering.
  These fixes can change scores, expression gates, permutation p-values, or FDR.
- Rerun analyses affected by the earlier positive-P row filter. Keeping evaluable
  zero events changes sample medians/prevalence when zero samples were omitted and
  expands the FDR correction family in either inference mode.
- Ambiguous prior aliases, duplicate headers, and invalid relevant expression now
  produce explicit errors. Correct those inputs rather than relying on silent
  column selection or numerical zero fallbacks.
- Re-export existing in-memory results when an older CSV omitted index-based
  metabolite/HMDB identifiers. The export fix preserves each identifier once.
- Reread original expression/metadata files and rerun analyses if older text
  import renamed duplicate genes or merged numeric-looking labels such as `01`
  and `1`. Metadata now defaults to text; convert additional numeric fields
  explicitly when needed.
- Replot affected visualizations if different names referred to the same HMDB
  ID, or multiple IDs shared a name. Previous plots could count aliases twice,
  omit matching samples, or put distinct relations on one dotplot row. Pass
  explicit IDs to target plots and reuse `selected_event_keys`, not display
  strings, for dotplot selection. Stored scores, permutation p-values, and FDR
  require no recalculation for these visualization fixes.

Database files remain local and excluded from GitHub during testing, as intended.
The database directory's schema guide is maintained separately from those files.

## Resolved follow-up: normalization of unscorable metabolites

The compiled scorer already returned only the observed event keys. Previously,
it normalized P/C/E for all compiled metabolites, while observed scoring first
excluded metabolites without evaluable production. With a finite C capacity near
`1e308` for an excluded metabolite, the observed run could finish but permutation
reference calculation could overflow and raise an error. Numeric checks rejected
the failure; no new event or accepted non-finite null score was produced.

Both paths now apply the same production-evaluability mask before P/C/E reference
normalization. Raw-expression and capacity checks remain upstream, and evaluable
zeros remain eligible. The compiled arrays retain their metabolite positions;
uncomputed rows use internal zero placeholders and stay excluded by the mask.
All observed cell types still participate in the references for eligible rows.

The regression cases cover missing production priors and entirely unmeasured
production genes, C/E capacities near `1e308`, both reference statistics and
inference modes, serial/parallel workers, and optional null storage. Compiled
scores match full recomputation; event scores, p-values, and FDR also match a
control with the excluded prior removed. Raw invalid values and capacity overflow
are still rejected, as is an invalid competing reference for an evaluable all-zero
production row. Analyses previously stopped by the unnecessary reference
calculation can be rerun.

## Resolved follow-up: dotplot event identity

The dotplot previously grouped and deduplicated by the readable
`metabolite -> sensor_gene` label. A full analysis with two same-name metabolites
and the same sensor retained both HMDB IDs, while automatic plotting and explicit
selection of both IDs displayed only one event per sender/receiver pair.

Ranking, selection, deduplication, and row coordinates now use the
`(hmdb_id, sensor_gene)` tuple. True duplicate pair events still prefer lower
FDR and then higher score. Display labels always include the HMDB ID. Tuples and
dictionaries match those two identifiers; legacy three-field tuples ignore their
name, and display-string selectors are rejected. Requested order is preserved
and repeated matches appear once. `selected_event_keys` returns the two-field
tuples alongside the `selected_events` display labels.

Regression cases cover same-name metabolites across sender panels, top-N
selection with and without FDR, supported selector forms, key reuse, duplicate
pair events, thresholds/QC, supplied axes, delimiter-containing fields, and actual
results from both inference modes. The plots retain the original scores,
p-values, and FDR without mutating the source event table.

## Resolved follow-up: Matrix Market array input

`read_anndata(..., mode="mtx")` previously called `.tocsr()` on every result of
`mmread()`. The common coordinate format returned a sparse matrix and worked,
but array-format files returned NumPy arrays and raised `AttributeError`.

The reader now checks the returned type before conversion. Coordinate inputs
retain the existing sparse path; array inputs retain their numeric arrays. Both
share the existing gene-by-cell transpose and identifier checks. No expression
values are filled, rounded, or otherwise transformed by the compatibility fix.
Previously rejected array files can now be read and analyzed.

Regression cases compare both formats using the same non-square integer/real
matrices, including gzip files, literal identifiers, single-gene/single-cell
shapes, and absent optional name files. Dense inputs also retain rejection of
incorrect name counts, duplicate names, and empty names. Analyses from both
formats match native AnnData scores, permutation p-values, and FDR in both
inference modes.

## Resolved follow-up: receptor context across prior-name aliases

Events use the Enzyme-side metabolite display name, while receiver scores retain
the Interaction-side name. These may differ for the same HMDB ID. Previously,
passing an event's name, ID, and receptor to `plot_receptor_expression_violin()`
could raise a not-found error even though the event and receptor were scored.
The failure was reproduced in both inference modes; selecting by ID alone worked.

An explicit ID is now required and uses the existing scoring normalization on
both the requested ID and receiver IDs. The receptor gene remains part of the
selection; the name only controls display. Missing or unmatched IDs raise an
error without falling back to a name or receptor-wide context. Original tables
and annotations are preserved. The regression cases cover actual aliased-prior
analyses, ID formats, same-name/different-ID contexts, missing or unknown IDs,
unknown receptors, and rejection of name-only/gene-only calls. Existing affected
results can be plotted again without rerunning scoring or permutations.

## Resolved follow-up: visualization identity independent of names

The earlier dotplot fix separated different IDs with the same name, but its
three-field key still treated different names for one ID as separate relations.
Counts and network deduplication also used names; sample selection could omit
samples whose name differed. Sender violin lookup and reaction selection could
likewise depend on a name. These inconsistencies could change what was shown
without changing the underlying scored identity.

All six public plotting functions now use normalized HMDB IDs and sensor genes
for relation identity. Sender/receiver fields distinguish directed records;
the sample plot additionally distinguishes samples. The production violin needs
only the HMDB ID. Lookup helpers no longer accept a name argument. Titles and
event labels use `Name (HMDB ID)` and include the gene for a relation. Names may
be absent or overridden without changing selection or numerical summaries.

Identical alias rows and repeated reaction definitions contribute once.
Conflicting sender values, reaction definitions, receiver contexts, or sample
records raise an error. Event counts, network, and dotplot keep their existing
filtering and duplicate-priority rules, using only identifiers as keys.
`unique_keys` cannot reintroduce name-based identity. Chinese code comments and
regression tests guard this rule against later optimization or compatibility
changes. README and notebook source examples document the explicit-ID API and
the two-field `selected_event_keys` return value.

## Resolved follow-up: unavailable dotplot significance

With `n_perms=0`, p-values and FDR correctly remain missing. The dotplot
previously checked only whether the FDR column existed before calculating its
size range and legend, producing two all-NA warnings and an `FDR: nan` legend.
Partially missing values also propagated NaN into dot sizes when the available
significance values differed: a reproducer returned three scored records but
rendered only two points.

Size ranges now use only available probabilities. Missing significance uses a
fixed-size diamond, whose color still represents the score. When every value is
missing, or both probability columns are absent, the plot uses these diamonds
and a `Significance: Unavailable` legend without a numeric significance scale.
With partially missing values, the numeric circle-size legend includes a
separate diamond key. Fixed size is the midpoint of the configured size bounds.
The missingness key is kept within the figure and separated from the colorbar
in both faceted layouts and a supplied axis.

An existing FDR column retains priority even when missing; unadjusted p-values
are not substituted for missing FDR. The historical p-value size scale is used
only when the FDR column is absent. Non-missing displayed probabilities must be
finite and in `[0, 1]`, with 0 and 1 distinct from missing. Explicit probability
filters retain their previous behavior, including excluding NA and reporting
an absent required column or an empty selection. Original scores, probabilities,
NA values, event identities, ranking and QC rules are preserved. Existing
results can be replotted without rerunning scoring or permutations.

## Validation record

The subsequent warning and edge-case review completed a full regression run
with **1019 passed, 1 skipped, and 213 warnings**. The warnings are the existing
AnnData string-index conversion and old H5AD encoding-metadata notices. No
scoring or plotting code was changed during this review. The two display
follow-ups below are outside the assertions covered by this passing suite.

The dotplot-significance revision added **44 regression cases**. The final
targeted run completed with **133 passed and 1 warning**, covering all new
significance tests, dotplot identity, shared visualization identity, plotting,
and receptor context. The warning is the existing AnnData string-index
conversion in a plotting fixture; the new cases introduced no warnings.
Coverage includes all/partial/no missingness, absent or custom probability
columns, nullable NA, valid 0/1 values, invalid probabilities, explicit filters,
both axis layouts, and actual results from both inference modes with and without
permutations. Rendered points retain their score colors, identity coordinates,
and original probabilities. Legend bounds and colorbar separation are checked
after drawing. PNG previews of all-missing, partially missing, and supplied-axis
plots were also visually inspected; `git diff --check` passed.

The visualization-identity revision added **32 regression cases**. Its focused
run completed with **268 passed and 1 warning**, covering visual identity,
dotplot selection, receptor context, all plotting tests, identifier and
expression validation, and numerical overflow.

The full regression run before the dotplot-significance revision, after the
visualization-identity, receptor-context, and Matrix Market changes, completed
with **975 passed, 1 skipped, and 213 warnings**. The warnings are the existing AnnData string-index
conversion and older H5AD encoding metadata warnings. The new visual identity
cases introduced no warnings.

All notebook code cells passed syntax checks. The affected code in three
example workflows also passed smoke checks with synthetic data whose Enzyme,
Interaction, and sample display names deliberately differ. These checks cover
ID-only prior selection, dotplot keys, both violins, sample scores, HMDB-indexed
intermediates, and heatmaps. Existing notebook outputs were preserved; the full
database-backed notebooks were not rerun. `git diff --check` also passed.

The receptor-context fix added 17 regression cases. Its targeted run completed
with **166 passed and 1 warning**, covering the new context tests, all plotting
tests, identifier and expression validation, and dotplot event identity. The
warning is the existing AnnData string-index conversion in a plotting fixture;
the new receptor-context cases introduced no warnings.

The Matrix Market compatibility fix added 18 regression cases. Its targeted run
completed with **80 passed** and no warnings, covering the new MTX cases, all text
identifier/import tests, and the existing MTX public-API test.

The earlier full regression run, after the initial dotplot identity fix and before
the Matrix Market compatibility and receptor-context changes, completed with
**908 passed, 1 skipped, and 213 warnings**. This fix added 19 regression cases.
The warnings came from existing AnnData index conversion and older H5AD metadata;
the new dotplot identity cases introduced no warnings.

The preceding permutation normalization-range milestone completed with
889 passed and 1 skipped and added 31 regression cases. Its focused run covering
normalization eligibility, production evaluability, numerical overflow, compiled
scoring, and permutation precision completed with **170 passed**.

The preceding text-input milestone completed with 858 passed and 1 skipped.
Those input fixes added 61 regression cases; their focused run with the two existing CSV/MTX loading checks
completed with **63 passed**. Cases cover duplicate/empty names, literal suffixes,
leading zeros, literal NA tokens, custom metadata labels, CSV/TSV and transposed
inputs, C/Python parsers, compression, parser options, MTX lists, and end-to-end
scores and permutations compared with native AnnData in both inference modes.

The earlier production-evaluability milestone completed with 797 passed and
1 skipped and added 35 cases. That coverage includes dense/CSR/CSC inputs,
selected expression layers, both inference modes, serial/parallel permutations,
optional null storage, and compiled scores versus full recomputation.

## Warning explanation: AnnData string-index conversion

The single warning in the focused plotting run is
`ImplicitModificationWarning: Transforming to str index.` It originates in
`test_single_cell_violins_accept_sparse_real_cellmesh_result`, whose `obs`
DataFrame uses its default integer index. The installed AnnData converts this
index to strings when constructing the object. The values of the expression
matrix and cell-type column are unchanged.

A separate reproducer compared implicit conversion with explicitly supplying
the same string indices. Explicit indices emitted no warning; the resulting
AnnData observations, expression values, CELL MESH scores, permutation p-values,
and FDR were exactly equal. Recommended cleanup: give this test explicit string
cell identifiers rather than suppressing warnings globally. The old-format
H5AD warnings in the full suite are a separate compatibility notice.

## Historical dotplot follow-ups: A07 and A08 resolved

1. **A07 historical finding, now resolved: zero FDR was mapped to the smallest positive FDR.** The previous
   log-transform floor is the minimum positive value in the displayed rows.
   For `[0, 0.01, 0.1]`, default dot sizes are `[260, 260, 20]`, and the numeric
   legend stops at `0.01`. For `[0, 1, 1]`, all sizes are `140`, and the legend
   shows only `1`. Returned probabilities remain correct. This affects external
   or rounded result tables containing zero; CELL MESH's add-one permutation
   p-values are strictly positive. Proposed fix: give zero an explicitly labeled
   display cap, preserve the positive-value ordering, and never label zero as a
   positive probability. Extend regression checks to size ordering and truthful
   zero legends; the existing 0/1 tests only verify availability and preservation.

2. **A08 historical finding, now resolved: complete-data layouts retained the older legend arrangement.** The previous
   adjustment applied only when significance was missing. With all values
   available, a supplied axis still places its numeric legend over the full-height
   colorbar. The default faceted p-value fallback can also place its longer
   `-log10(perm_pvalue)` title outside the figure. Both cases were reproduced;
   the missing-data layout remained separate and within bounds. Proposed fix:
   share the legend/colorbar layout rules across all missingness states and both
   statistic types, then check rendered bounds and overlap for every branch.

Both findings concern visualization, without changing stored event scores,
probabilities, identities, or inference. A07 and A08 are now resolved as recorded
above; the original reproductions are retained as historical context.

To reproduce the full check from the project root in the `cellmesh` environment:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/cellmesh-review-mpl \
PYTHONDONTWRITEBYTECODE=1 \
python -m pytest -q -p no:cacheprovider --tb=short
```

Earlier counts in the Enzyme review refer to their explicitly described
historical stages.
