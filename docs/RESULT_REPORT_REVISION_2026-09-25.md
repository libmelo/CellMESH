# Result persistence and explanatory reporting revision

Date: 2026-09-25. Base commit: `3f70744ca2e933d40410779b13f4d550023cc213`.
Package version: 0.5.0. This revision is local work; it has not been pushed to GitHub.

## Delivered behavior

- `CellMeshResult.save(prefix, report=False, events=None, include_null_scores=False)`
  writes a versioned HDF5 result file; `report=True` additionally writes one
  self-contained HTML. Default save omits raw cell expression and null matrices.
- `CellMeshResult.load(path)` restores result tables, supported pandas dtypes,
  indexes, missingness, parameters, prior snapshots and availability intermediates.
- `coverage_report()` exposes descriptive counts and explicit denominators.
  `explain_event()` exposes actual sender/receiver components, reaction selection,
  gene expression summaries, references, source priors and candidate diagnostics.
- HTML embeds four component-chart panels and foldable tables. The optional
  `plot_event_explanation()` returns a Matplotlib figure for PNG/PDF export.
- `to_csv()` retains its existing schema, naming, identifiers and behavior.
  No changes were made to score formulas, event selection, reference rules,
  permutation computation or FDR. New coverage statistics do not filter events.

The implementation captures unweighted and weighted reaction scores and receiver
reference/gating values during observed scoring. Prior snapshots retain unmeasured
sensors, allowing missing-candidate explanations. Detailed report tables are
constructed on request, not during every permutation.

## Verification

1. Before editing, captured six fixtures from the existing result-export tests:
   pooled/sample-aware crossed with positive production, measured-zero production
   and absent production priors. After the final edits, event tables (including
   P/FDR), sender/receiver scores, sample events, validation and parameters match
   exactly. Baseline artifact: `/tmp/cellmesh_reporting_baseline.pkl`.
2. Initial affected-path suite: **140 passed**, 13 warnings, 29.82 seconds.
   Log: `/tmp/cellmesh_reporting_targeted.log`.
3. Full repository suite before the final storage/figure refinements:
   **2504 passed, 7 skipped**, 221 warnings, 274.33 seconds.
   Log: `/tmp/cellmesh_reporting_full.log`.
4. After final storage encoding, HTML cleanup, figure and public API changes,
   reran the new feature tests plus export, database/API, reaction consistency,
   receiver reference and production evaluability suites:
   **177 passed**, 13 warnings, 91.77 seconds.
   Log: `/tmp/cellmesh_reporting_final_targeted.log`.
5. Ran `examples/export_explain_report.py`, reloaded the generated file, checked
   events and queried coverage/one event. Produced HTML plus HDF5; rendered PNG
   and PDF from the reloaded result. PNG visually inspected.
6. `git diff --check` passed.

New tests cover exact result round trips, partial/missing/zero gene evidence,
reaction deduplication, score decomposition, sample endpoint gaps, unmeasured
receptors, excluded production, nullable/categorical/index metadata, explicit null
storage, HTML escaping, separate figures, save rollback, obsolete HTML cleanup,
symlink rejection and incompatible schema versions.

## Demonstration artifacts

- `/tmp/cellmesh-report-demo/analysis.cellmesh.h5` — synthetic example result.
- `/tmp/cellmesh-report-demo/analysis.report.html` — standalone report.
- `/tmp/cellmesh-report-demo/event_explanation.png` and `.pdf` — optional figure.
- `/tmp/cellmesh_reporting_demo.log` — example output.

The synthetic example result is 4,218,296 bytes; HTML is 59,405 bytes.

Also exercised the complete default priors (11,687 enzyme rows and 6,478
interaction rows) with 90 synthetic cells, three samples and three cell types.
This produced 63 events. Final-format save took 16.24 seconds, produced a
9,630,470-byte result and 251,154-byte HTML; load took 0.76 seconds, and event
tables matched exactly. Artifact prefix: `/tmp/cellmesh_full_prior_report`.
Timing is a single local smoke measurement, not a general performance benchmark.

The generated expression demonstrates software behavior, not biological accuracy.
The bundled `read_example_data()` uses synthetic `GeneN` names, so its example
uses explicitly matched synthetic priors rather than the default biological database.

## Limits and usage

- No original cell-by-gene matrix is embedded. Raw-cell plots still require AnnData.
- A legacy CSV export cannot recover intermediates that were never saved.
- HTML previews at most 100 rows per table and ten sample units per chart;
  complete tables are queryable in Python.
- Input gene coverage does not establish reaction function or transport feasibility.
- The container stores original availability intermediates alongside report tables
  to restore current consumers; these partly overlap and add storage overhead.
- Writes roll back ordinary I/O failures, but multi-file publication is not a
  crash-safe transaction and does not support concurrent writers to one prefix.
- The original three untracked 2026-09-22 audit artifacts were retained as-is.

See [format and usage documentation](RESULT_FILES_AND_REPORTS.md) and
[runnable example](../examples/export_explain_report.py).
