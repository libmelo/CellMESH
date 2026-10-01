# Result files and explanatory reports

The current [three-level report guide](THREE_LEVEL_REPORTS.md) specifies independent
condition-filter and explicit metabolite–sensor reports, append methods and the complete output layout.

CELL MESH 0.5 adds a versioned `.cellmesh.h5` result file and an optional,
self-contained HTML report. These features describe the existing scoring model;
coverage never adds a score multiplier, changes the tested events, or changes FDR.
`run_cell_mesh()` still returns an in-memory `CellMeshResult`. It does not write
files automatically. Existing `to_csv(prefix)` filenames and table schemas are
unchanged.

CELL MESH 0.6 additionally offers an opt-in standard `.input.h5ad` companion
containing the selected cell-by-gene expression matrix and annotations. Default
memory use remains unchanged. The full, unfiltered `.all_events.csv` companion
requires explicit `include_all_events=True`; it is no longer a default output.

Since 0.8, every `save()` also writes `.positive_events.csv`, independently of
selection and `report=True`. Defaults are BH-adjusted p-value within sensor type
<= 0.05, without implicit QC/score gates. Missing statistics produce an empty
positive CSV and an explicit unavailable-assessment status. See the
[three-level delivery guide](THREE_LEVEL_REPORTS.md) for current defaults and appending.
HTML is now included by default (`report=True`); use `report=False` for data-only
saving. The positive CSV is always written; the full CSV has its own opt-in flag.

## Save and reload

```python
from cellmesh import CellMeshResult

# res was returned by run_cell_mesh(...).
files = res.save("analysis")
# analysis.cellmesh.h5, analysis.positive_events.csv, analysis.report.html,
# and analysis.reports/index.html.
restored = CellMeshResult.load(files["result"])

selector = dict(
    hmdb_id="HMDB0000001", sensor_gene="YOUR_SENSOR",
    sender="YOUR_SENDER", receiver="YOUR_RECEIVER",
)
# Replace these identifiers with an event from your own results.
selected = res.select_events(
    metabolite_sensor_pairs=[(selector["hmdb_id"], selector["sensor_gene"])],
    cell_pairs=[(selector["sender"], selector["receiver"])], source="all",
)
files = res.add_event_reports(selected)
# Also includes analysis.positive_events.csv and analysis.report.html.
```

Missing parent directories are created. A full `.cellmesh.h5` filename is also accepted.
`save()` returns a dictionary of output `Path` objects, including `result` and
`positive_csv`, plus `events_csv` only with `include_all_events=True`,
`filtered_csv` and `filtered_report` when a filter is supplied, `selected_csv` for explicit pairs, and
`report` and `navigation` by default, `pair_report:<stable_id>` for each explicit molecular pair, and `input` when requested. The standard output
names belong to the prefix: saving again replaces requested outputs, `report=False`
removes obsolete same-prefix `.report.html` and `.filtered.report.html`, and `include_input=False` removes
an obsolete same-prefix `.input.h5ad`. Other files, including the independent
legacy `to_csv()` outputs such as `.events.csv`, are untouched. Symbolic-link
and directory destinations are rejected.
An existing `.all_events.csv` is left untouched when `include_all_events=False`;
it is not returned or linked as part of the new delivery.
All requested outputs are staged before publication; ordinary publication
failures roll back prior files. This is not a crash-safe multi-file transaction and
does not support concurrent writers to one prefix.
After successful saving, the in-memory result's archive association refers to the
new bundle; failed saving leaves that association unchanged.

The result file stores expression summaries and relevant prior snapshots, not
the original cell-by-gene expression matrix. Receiver references and reaction
activities are captured from the scoring calculation. Reports do not independently
reimplement scoring. Detailed tables are materialized when requested or saved,
not during every permutation.

## Full communication-event CSV (0.6.1)

`save(..., include_all_events=True)` writes `prefix.all_events.csv` from **all rows and columns of
`result.events`**, in their original order, with UTF-8 encoding and no extra
DataFrame index column. This includes zero scores, nonsignificant events and
events failing cell-count QC. HTML event selectors and its preview limits do not
filter the CSV. No new score, p-value or FDR calculation is performed.

The table includes sender/receiver cell types, HMDB and sensor identities, scores,
permutation p-values, both FDR columns and QC annotations. In sample-aware mode,
one row is a cross-sample aggregate event; its additional sample-support columns
are included. Sample-level event rows remain in the HDF5 `sample_events` table
and can be exported with the existing `to_csv()` API if needed. This is a
cell-type communication table, not every individual cell-to-cell pair or a
Cartesian list of unscorable candidates.

Missing statistics (for example with `n_perms=0`) are written as empty fields;
measured zero scores remain zero. Empty computed results still produce column
headers. CSV preserves tabular values, while the HDF5 remains the source for
exact pandas types, attributes, intermediate quantities and explanations.

```python
files = res.save("analysis", include_all_events=True)
print(files["events_csv"])  # analysis.all_events.csv
```

The distinct `.all_events.csv` suffix leaves the existing multi-table
`to_csv(prefix)` export namespace and manifest unchanged. All requested files,
including this CSV, share staging, destination checks and ordinary-error rollback.
The default is False on every save or standalone general HTML export, including
loaded results. Full HDF5 content and report counts do not depend on this flag.

## Independent filter and explicit-pair scopes

`filter_events()` defines tier 2. Independent `hmdb_ids` and `sensor_genes`
lists match all existing combinations satisfying both, never positional pairs.
`select_events(metabolite_sensor_pairs=[(hmdb, sensor), ...])` defines tier 3:
one page per explicitly requested pair, containing its matched cell-pair events.
Both default to `source="positive"`; use `source="all"` to include nonpositive
or statistically unavailable generated events. Values and BH correction remain unchanged.

```python
filtered = res.filter_events(sender_labels=["A", "C"], receiver_labels=["B"])
selected = res.select_events(metabolite_sensor_pairs=[("HMDB0000961", "CYSLTR1")])
files = res.save("analysis", filtered_events=filtered, selected_events=selected)
filtered.to_csv("filtered.csv")
selected.to_csv("selected.csv")
```

The condition-filter report is `analysis.filtered.report.html`, with
`analysis.filtered_events.csv`. Explicit-pair pages are under
`analysis.reports/metabolite_sensor/`; their combined CSV is
`analysis.metabolite_sensor_events.csv`. The general report and
`analysis.reports/index.html` retain their locations. H5 always includes all
scientific events, while `.all_events.csv` remains opt-in.

Different conditions intersect; list members are alternatives. Directed cell
pairs intersect sender/receiver conditions. `None` is unrestricted; `[]` matches
nothing. Explicit empty pair requests still produce diagnostic pages, never
fabricated scientific rows. Unknown identifiers and duplicate requested pairs
are errors before publication.

`add_filtered_report(filtered)` and `add_event_reports(selected)` independently
append to an existing bundle, requiring `replace=True` for existing scopes.
General files and scientific data stay unchanged. Stored explanation tables are
reused; only delivery metadata is updated in the H5. Ordinary publication errors
roll back. Complete tables/CSV and requested pair pages are never truncated by
`plot_top_n`.

Both objects expose defensive `.events`, `.criteria`, `.summary` and `.to_csv()`.
Changes to the source result invalidate them. `get_filtered_events()` restores
the tier-2 filter; `get_selected_events()` returns independent stored pair scopes.
No prior scope is implicitly restored by a new `save()`.

`filtered.to_html(path)` writes a standalone filter report and its CSV;
`selected.to_html(path)` writes pair pages, their navigation and combined CSV.
Managed bundle paths are protected; use append methods for updates. Plotting
functions accept either object without silently adding filters. Scientific
missingness and zero values retain their existing interpretation.

See [the full report guide](THREE_LEVEL_REPORTS.md) for API examples, report
contents, graph encoding, file keys and migration from the old `select_events()`
and complete-event HTML interfaces. `explain_event()` remains available for
individual four-field identities.

## Optional cell-level input archive (0.6)

```python
from cellmesh import run_cell_mesh, CellMeshResult

res = run_cell_mesh(
    adata, enzyme_metabolite, metabolite_sensor,
    cell_type_key="cell_type", sample_key="sample",
    sample_mode="sample_aware", layer="normalized",
    retain_input=True,
)
files = res.save("analysis", include_input=True, report=True)
# analysis.cellmesh.h5, analysis.positive_events.csv, analysis.input.h5ad,
# analysis.report.html and analysis.reports/index.html

restored = CellMeshResult.load(files["result"])  # Does not read the input matrix.
input_adata = restored.load_input()             # Explicit in-memory AnnData.
# The selected expression matrix is now in input_adata.X, even if the original
# run used adata.layers["normalized"]. Use layer=None when rerunning this input.

disk_input = restored.load_input(backed="r")    # Read-only disk-backed X.
try:
    print(disk_input.shape)
finally:
    disk_input.file.close()
```

- `retain_input=False` is the run default: no snapshot is copied. Enabling it
  captures an independent matrix and obs/var annotations before scoring, and the
  run scores that snapshot. Later changes to the caller's AnnData cannot change
  the captured input. This costs additional memory proportional to the selected
  matrix and annotations, including for a disk-backed input.
- The snapshot includes **all supplied cells and genes**, in their original
  order, including genes not used by the current prior. Only the selected matrix
  is copied into `.X`. Values, dtype and sparse storage are preserved for dense,
  CSR and CSC inputs; backed slicing may change the compressed sparse layout.
  No other layers, `.raw`, embeddings, graphs or unrelated `.uns` are copied.
- `.obs` and `.var` preserve cell/sample/type and gene annotations; they must be
  serializable by AnnData. Unsupported annotations fail before publishing any
  requested output, leaving previous outputs intact. This is an analysis-input
  archive, not a backup of every attribute in the user's original AnnData.
- `input_provenance` and `.uns["cellmesh_input"]` describe the source expression
  layer, matrix dimensions, dtype, grouping keys and capture scope. Values are
  **as supplied**, with no new expression transformation. Their previous
  processing cannot be inferred: the archive does not claim to contain raw
  counts or reconstruct upstream normalization/filtering history.
- `include_input=False` is the save default, even for a run that retained input.
  Enabling it requires a captured snapshot or a loaded result associated with an
  existing archive. Missing input cannot be reconstructed from result summaries
  or old CSV exports. There is no automatic attachment of an unrelated AnnData.
- The result file's `metadata/expression_input` stores capture provenance;
  `metadata/input_archive` stores a sibling filename, SHA-256 and byte size.
  `metadata/storage.input_companion_included` distinguishes the companion from
  matrices embedded inside the result file (none). Move the result and its
  companion together; relative filenames are resolved beside the loaded result,
  independently of the working directory. `load_input()` checks the complete
  file checksum before opening it. Backed mode avoids loading X into RAM, but
  checksum verification still reads the file once. Changes made externally to
  the companion require a new analysis/archive, rather than silently reusing it.
- A missing or changed companion does not prevent loading scores or generating
  an explanation report. Explicit input loading or re-archiving fails with a
  clear error. HTML displays the **recorded** archive association and provenance;
  generating HTML alone does not verify companion availability or checksum.
- A loaded bundle can be saved to another prefix with `include_input=True`;
  its verified companion is copied without materializing X. Ordinary
  `save(prefix)` drops the companion association from the newly saved result
  and removes a stale same-prefix input file. A retained in-memory snapshot is
  still available through `load_input()`; its returned AnnData is a copy.
- `load_input(backed="r")` requires a successfully saved or loaded companion.
  Writable backed modes are rejected. Standard `anndata.read_h5ad(files["input"])`
  is also supported, but it bypasses CELL MESH's association/checksum checks.
- Existing `plot_metabolite_secretion_violin(result, input_adata, ...)` and
  `plot_receptor_expression_violin(result, input_adata, ...)` recognize the
  archived source-layer mapping and read `.X` automatically. Explicit `layer=`
  overrides still select that named layer and error if it is absent.

For a repeat run, use the saved priors, original scoring/permutation parameters
and `layer=None` (the archived matrix is in `.X`). Reproducibility also depends
on software versions and the numerical environment. The runnable example below
checks exact event equality on replay. Schema version 1 remains readable: the
new metadata fields are optional, and older result files load without them.

## File layout (schema version 1)

```text
analysis.cellmesh.h5
├── metadata/
│   ├── package_version
│   ├── parameters
│   ├── input_provenance       # prior inputs; expression layer is in parameters
│   ├── expression_input       # optional expression snapshot provenance (0.6)
│   ├── input_archive          # optional companion filename/checksum/size (0.6)
│   ├── delivery               # optional positive criterion and interests (0.8)
│   └── storage                # null-matrix policy, explanation availability
├── results/
│   ├── events
│   ├── sender_scores
│   ├── receiver_scores
│   └── celltype_qc
├── sample_results/            # entries are None when not applicable
│   ├── sample_validation
│   ├── sample_events
│   ├── sample_sender_scores
│   └── sample_receiver_scores
├── explanation/
│   ├── sender_components
│   ├── receiver_components
│   ├── reaction_details
│   ├── gene_details
│   └── reference_context
├── coverage/
│   ├── summary
│   └── details
├── diagnostics/
│   ├── production_status
│   └── event_status
├── priors/
│   ├── enzyme
│   └── interaction
├── availability_results      # original intermediate layout for existing consumers
└── optional/
    └── permutation_null      # only if explicitly included and already available
```

The root attributes identify the format and schema version. These are HDF5 groups
and datasets inside **one file**, not independently exported tables. Numeric arrays
use native HDF5 datasets; object cells and metadata use a restricted, tagged JSON
encoding. No pickle payloads are stored or executed. Table names, indexes,
MultiIndex levels, nullable pandas dtypes and supported table attributes are
preserved. Nested lists/dictionaries, timestamps and missing values are supported;
arbitrary Python class instances raise an error before publication.

The `availability_results` group preserves the original object layout, while the
`explanation` groups expose report tables. This intentionally adds some storage
overhead. Large null matrices have one canonical storage location. The result
file is specific to CELL MESH; it is not an AnnData `.h5ad` file.

## Query coverage and explain an event

```python
coverage = restored.coverage_report()
coverage["summary"]       # counts and explicit denominator
coverage["details"]       # per sample / observed type / metabolite / P,C,E

detail = restored.explain_event(**selector)
detail["sender_components"]
detail["receiver_components"]
detail["reaction_details"]
detail["gene_details"]
detail["reference_context"]
detail["sample_events"]
detail["enzyme_prior"]
detail["interaction_prior"]

# In sample-aware mode, optionally restrict the explanation to one sample.
one_sample = restored.explain_event(**selector, sample="S1")

# Original CSV export remains available. Detailed tables can be exported selectively.
restored.to_csv("legacy_analysis")
detail["reaction_details"].to_csv("selected_event_reactions.csv", index=False)
selected = restored.select_events(
    metabolite_sensor_pairs=[(selector["hmdb_id"], selector["sensor_gene"])],
    cell_pairs=[(selector["sender"], selector["receiver"])], source="all",
)
selected.to_html("selected_event.html")

# Optional separate PNG/PDF figure (requires matplotlib).
from cellmesh import plot_event_explanation
fig, axes = plot_event_explanation(restored, **selector)
fig.savefig("selected_event.png", dpi=200, bbox_inches="tight")
# fig.savefig("selected_event.pdf", bbox_inches="tight")
```

Identity uses HMDB + sensor gene + sender + receiver, with sample added for
sample-level records. Metabolite names are display annotations only. In pooled
mode, report tables use `sample=None`; a pooled result cannot be decomposed into
independent sample scores after saving. In sample-aware mode, `aggregate_event`
always contains the original all-sample summary, even when `sample=` selects one
sample for the detailed explanation. It is not recalculated from component medians.

Unscored candidates can also be queried: absent endpoints, unavailable production
genes, missing interaction priors, unmeasured sensors, and disabled self events
have explicit reasons. A zero score is evaluated, and is kept distinct from NA.
Receiver gating is reported separately from unavailable input. No Cartesian
product of all nonexistent events is stored: `diagnostics/event_status` covers
generated events and structural sample gaps; other candidates are diagnosed on demand.

## Interpretation of completeness

Coverage uses retained reaction gene sets after the same complete-gene-set
deduplication as scoring. The reaction table also lists omitted duplicate and
strict-subset reactions, including the retained reactions that represent them.

- Gene coverage is the number of unique prior genes present in the input matrix
  divided by the unique prior gene count, separately for P/C/E.
- A direction without a prior has an undefined coverage fraction (NA), not 0%.
- Full/partial/uncovered reaction counts use retained reactions as the denominator.
  The pre-deduplication annotated reaction count is reported separately.
- Zero expression is an observed numeric value, not an input-gene absence.
  Transcript dropout and numerical underflow can also produce zeros; a zero is
  not proof that a biological function is absent. Numerical diagnostics remain
  available in run parameters.
- Positive reaction counts describe pre-abundance activity in the particular
  sample and type. The weighted capacities and actual score remain separate.
- Summary counts use sample × observed type × union of valid prior HMDB IDs.
  They are not counts of independent donors or unique metabolites.
- Receptor panel coverage is recorded in `receiver_components.input_present`;
  prior entries with unmeasured receptors remain inspectable there.
- Coverage is an input-description metric, not reaction functionality, confidence,
  transport feasibility, a p-value, or a new quality threshold.

If absent genes were filled with zeros when datasets were combined, the original
assay coverage cannot be inferred from the combined matrix alone. The report
describes the input matrix actually supplied.

## HTML and optional data

General HTML includes global counts, positive communication overview, the complete
positive table and global methods. Tier 2 contains selection criteria, collection
comparisons, the complete selected table and filtering details. Tier 3 contains
one molecular pair's communication distribution, shared endpoint evidence,
all sample performance and collapsed cell-pair details. See the [current report guide](THREE_LEVEL_REPORTS.md).

The general overview contains count and network plots. Selection adds a
metabolite–sensor dotplot and, in sample-aware mode, a sample-by-event heatmap.
Explicit molecular-pair requests add separate pages, independently of the filter; none are
automatically chosen. All sample components can be accessed, and event detail
tables are fully embedded and paginated. Images expand without internal scrolling.

Expression figures use explicit `adata` or retained input and never implicitly
load an archive. They describe expression, not measured concentrations or flux.
No external scripts or server are needed. Keep the relative file layout when
moving the report bundle.

## Runnable example

```bash
python examples/export_explain_report.py --output examples/report_demo_output --include-input
# Optional selection CSV, focused HTML, and aggregate dotplot:
python examples/export_explain_report.py --output examples/selection_demo_output --include-input --selection-demo
# Larger offline browser example: 216 real computed synthetic events, 3 explicit molecular-pair reports.
python examples/browse_event_report.py --output examples/report_browser_demo_output
```

The example uses generated expression and explicitly constructed synthetic priors.
With `--include-input`, it adds the input companion, reloads both result and input,
checks exact event equality after repeating the analysis, and queries coverage
and one event. The report bundle also includes navigation and explicitly requested molecular-pair pages. It is a workflow
demonstration, not biological validation.
The selection variant keeps all nine events in HDF5/full CSV, selects A/C → B
(two events), explicitly requests one independent molecular-pair page and draws both in the PNG. Its two extra
deliveries are `analysis.filtered_events.csv` and `selected_events.png`.
