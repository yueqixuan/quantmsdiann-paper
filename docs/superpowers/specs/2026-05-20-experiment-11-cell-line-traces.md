# Experiment #11 — Cell-line traces for wallclock-vs-size scatter

**Status:** spec (not yet implemented)
**Date:** 2026-05-20
**Owner:** ypriverol@gmail.com
**Manuscript figures:** F2a, F2b (broadening), F2c (broadening)
**Brainstorm reference:** [§2.2 / §5.2](../../brainstorming.md#5-experiments)

## Goal

Collect `nextflow_trace.txt` artefacts for the three cell-line
reanalyses (PXD003539, PXD004701, PXD030304) so the scaling figures
can pool benchmarks + cell lines into a single "all runs, all
instruments" cloud — the actual story the manuscript wants to tell.

## Background

The published PRIDE artefacts under
`quantms-collections/absolute-expression-2.0/cell-lines/<PXD>/` do
**not** include `nextflow_trace.txt` for any cell-line run (verified
in the original benchmarks design doc). Consequently F2a today shows
benchmarks-only, weakening the "independent of instrument" framing
because the benchmarks span only 4 instruments while the cell lines
add a fifth (TripleTOF) and broader file-size range.

## Inputs

- The 3 cell-line SDRFs already cached under
  `data/<PXD>/<PXD>.sdrf.tsv`.
- A scratch SLURM / cluster allocation big enough to rerun the
  pipeline end-to-end on each cohort at one DIA-NN version (v2.5.0).
- The cached FASTA per cohort (the original reanalysis already
  downloaded them).

## Outputs

- Three new `nextflow_trace.txt` files (one per cohort), staged under
  `data/<PXD>/v2_5_0/nextflow_trace.txt` for consumption by F2a/F2b/F2c.
- Updated `data/parallelism_data.tsv` and `data/runtime_per_step.tsv`
  with cell-line rows; `analysis/figure_performance_runtime.py` and
  `figure_performance_trace.py` already loop over a `CELL_LINE_DATASETS`
  tuple — adding cell-line URLs to the fetch logic is the only code
  change.
- Re-rendered `parallelism_vs_wallclock.svg` and `runtime_per_step.svg`
  now showing all 23 analyses (20 benchmarks + 3 cell lines).

## Validation criteria

- Each new trace has > 0 `COMPLETED` rows.
- Total wallclock from the trace matches `pipeline_report.txt`
  duration to within 5 % (allowing for SLURM queue time differences).
- Instrument labelling in F2a still resolves each new point to a
  single instrument family (the existing `_majority_instrument()`
  helper handles mixed cohorts).

## Risks

- **Compute cost.** PXD030304 is 5798 files; a full rerun at v2.5.0
  is the largest single cohort in the paper outside experiment #12.
  Budget needs sign-off before scheduling.
- **Reproducibility of the original results.** We are reproducing
  reanalyses that were originally run at v1.8.1; the v2.5.0 IDs/quant
  will differ. Document the version used in the new traces; do not
  conflate with the original cell-line counts.
- **Storage.** The traces are small (~MB) but if we also re-publish
  `diann_report.parquet`, those are multi-GB. Decide before the run
  whether to keep the parquet or trace-only.

## Plan

1. Schedule the 3 reruns at v2.5.0 with `nextflow -with-trace -with-report`.
2. Stage the trace + pipeline_report files into `data/<PXD>/v2_5_0/`.
3. Add a `cell_line_url(dataset)` URL builder to
   `figure_performance_runtime.py` (the constant
   `CELL_LINE_BASE` already exists; just needs activating).
4. Re-run `python -m analysis.figure_performance_runtime` and
   `figure_performance_trace`; verify new dots appear in the SVGs.
5. Capture the new TSVs under `analysis/figures/performance/data/`.

## Out of scope

- Cross-version reruns on the cell lines. That's experiment #13.
- Resubmission of cell-line results to PRIDE. Separate workflow.
