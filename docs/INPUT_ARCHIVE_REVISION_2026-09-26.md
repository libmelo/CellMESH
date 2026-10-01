# Optional input archive revision

Date: 2026-09-26. Package version: 0.6.0.
This extends the local 0.5 result/report revision; changes have not been committed
or pushed to GitHub.

## Delivered behavior

- `run_cell_mesh(..., retain_input=True)` captures an independent copy of the
  chosen expression matrix and obs/var before scoring, then scores that copy.
  The default does not capture or retain a cell-level matrix.
- `result.save(prefix, include_input=True)` adds `prefix.input.h5ad`, a standard
  AnnData archive containing all supplied cells and genes, with the selected
  expression layer stored in `.X`. Other layers, raw, embeddings, graphs and
  unrelated uns are omitted. Values remain as supplied; raw counts and upstream
  processing are not inferred. Sparse inputs remain sparse.
- Result HDF5 retains the original source layer in run parameters, plus optional
  capture metadata, a relative companion filename, SHA-256 and file byte size.
  The result HDF5 itself still contains no cell-level expression matrix.
- `CellMeshResult.load()` does not load or require the companion matrix.
  `load_input()` explicitly reads a verified archive or returns an independent
  copy of a retained snapshot. `load_input(backed="r")` supports read-only disk
  access to a saved archive; checksum verification reads the file in chunks.
- Bundle relocation follows relative companion paths. Re-saving a loaded
  result with input copies the verified companion without materializing X.
- HDF5, H5AD and HTML are staged before publication. Ordinary publication or
  stale-file cleanup errors roll back all affected output files and preserve the
  prior in-memory association. Successful save records the new association.
- Standard same-prefix optional files are removed when omitted on subsequent
  saves, preventing stale input/report files from appearing to describe new
  results. CSV exports and unrelated files are unaffected.
- HTML identifies the recorded input archive, expression source and matrix
  dimensions. HTML alone does not open or verify the companion.
- Existing sender/receptor expression violin plots recognize the archive's
  source-layer mapping and read `.X` automatically, while preserving explicit
  layer overrides. Their plotted cell-level values match those from the original
  chosen layer, for both in-memory and read-only disk-backed archive loading.
- New metadata entries are optional additions to schema version 1. Files from
  the 0.5 revision remain loadable without an input companion.

No scoring formula, tested-event selection, permutation calculation or FDR rule
has changed. Input capture is a storage option, not a scientific parameter.

## Verification

- Input archive, existing reporting/persistence and disk-backed input suites:
  **149 passed**, 8 warnings, 46.13 seconds.
  Log: `/tmp/cellmesh_input_archive_targeted.log`.
- Full repository suite before the final source-layer mapping adjustment for
  archived-input violin plots: **2540 passed, 7 skipped**, 230 warnings,
  291.47 seconds. Log: `/tmp/cellmesh_input_archive_full.log`.
- After that adjustment, reran input archive, plotting, visual identity, backed
  input and sparse duplicate-coordinate suites: **709 passed**, 9 warnings,
  69.72 seconds. Log: `/tmp/cellmesh_input_archive_final_targeted.log`.
- Exact event, sender/receiver, sample-result and null-score equality between
  ordinary runs, retained-input runs and archived-input replay, across dense,
  CSR and CSC matrices, both sample modes, and `.X`/selected-layer inputs.
- Original-object mutation isolation; correct annotation and matrix order for
  disk-backed views; unchanged source input files; sparse preservation;
  missing/corrupted companion behavior; relocation; legacy files; strict option
  validation; symlink/path checks; staging/publication/cleanup rollback.
- Eight additional exact capture/replay checks on noncanonical int8 sparse
  matrices and arrays (CSR/CSC, both sample modes) preserve duplicate coordinates,
  stored values and dtype without premature duplicate reduction or overflow.
- Runnable synthetic example regenerated with `--include-input`. Its archived
  **50 cells × 50 genes** input reproduces the original event table exactly.
  Log: `/tmp/cellmesh_input_archive_demo.log`.
- Final `git diff --check` passed.

## Inspectable example

Run from the repository root:

```bash
python examples/export_explain_report.py --output examples/report_demo_output --include-input
```

- [Result HDF5](../examples/report_demo_output/analysis.cellmesh.h5)
- [Cell-level input H5AD](../examples/report_demo_output/analysis.input.h5ad)
- [HTML report](../examples/report_demo_output/analysis.report.html)

The example contains synthetic expression and synthetic priors; it demonstrates
software behavior, not biological validation.

## Practical limits

- Enabling capture uses additional RAM for the selected matrix and annotations,
  including when starting from disk-backed input. Keeping it disabled preserves
  the previous memory behavior.
- An input snapshot is not a full backup of all original AnnData attributes.
  obs/var annotations must be supported by the installed AnnData H5AD writer.
- Replaying an archived selected layer requires `layer=None`, because archived
  expression is in `.X`. Reuse original scientific parameters and saved priors;
  numerical reproducibility also depends on software/environment versions.
- Ordinary error rollback is not crash-safe multi-file atomic publication.
  Concurrent writers to the same prefix are not supported.

See [usage and file schema](RESULT_FILES_AND_REPORTS.md).
