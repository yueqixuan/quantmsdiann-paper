# Experiment #13 — Cell-line cross-version sweep

**Status:** spec (not yet implemented)
**Date:** 2026-05-20
**Owner:** ypriverol@gmail.com
**Manuscript figures:** F3.PXD003539, F3.PXD004701, F3.PXD030304, F3
**Brainstorm reference:** [§2.3 / §3 / §5.2](../../brainstorming.md#5-experiments)

## Goal

Carry the "evolution" claim (§2.3) from the benchmark side (F1a) onto
the use case (§3). Show that re-running the cell-line cohorts at
successive DIA-NN versions delivers monotonic ID/accuracy
improvements, with v2.5.0 substantially ahead of the original
`quantms`-1.8.1 baseline that PRIDE currently hosts.

## Background

Today each cell-line dataset is processed at only **one DIA-NN
version**. The "evolution" half of the manuscript is therefore
benchmark-only, which is the weakest framing — readers want to see
the gain on real cohorts.

Per the [§9.2 decision](../../brainstorming.md#9-decisions-taken-provisional),
we run endpoints first (v1.8.1 + v2.5.0) and add the three
intermediate versions only if compute permits.

## Inputs

- The 3 cell-line SDRFs and FASTAs already cached.
- A scratch cluster allocation sized for **6 reruns** (3 datasets × 2
  versions) at endpoint scope, or **15 reruns** at full-sweep scope.

## Outputs

- New `data/<PXD>/<version>/diann_report.pr_matrix.tsv` +
  `.pg_matrix.tsv` per (dataset, version) pair.
- Updated per-cohort figures (F3.PXDxxx) showing the v1.8.1 → v2.5.0
  delta on the same protein-count / gene-count / coverage axes the
  current figures use.
- Updated combined atlas (F3) gaining a "by-version" view.
- New TSV: `analysis/figures/combined/data/cell_line_version_progression.tsv`
  with rows `(dataset, version, n_precursors, n_proteins,
  n_unique_genes)`.

## Validation criteria

- v2.5.0 protein counts must be **≥** v1.8.1 counts on each cohort
  (this is the load-bearing finding the figure relies on; non-trivial
  failure would invalidate the framing).
- Total unique proteins across all 3 cohorts at v2.5.0 must be a
  strict superset of the original `quantms`-1.8.1 published numbers
  to within ±2 % (allowing for re-search variance).
- SDRF parsing must not drop runs vs the original reanalysis.

## Risks

- **Compute cost.** PXD030304 (5798 files) at v1.8.1 is expensive
  (1.8.1 is the slowest of the supported versions on this scale).
- **Reproducibility.** v1.8.1 reruns will diverge from the original
  `quantms` reanalysis (different pipeline versions, possibly
  different SDRF). Compare against the original *published* counts in
  the paper, not against an in-house v1.8.1 rerun.
- **Storage.** 6 (or 15) new pr_matrix.tsv files at ~10 MB each is
  trivial; pg_matrix is ~MB; full `diann_report.parquet` per analysis
  is ~GB and may or may not be retained — decide before runs.

## Plan

1. **Endpoint reruns (6 analyses).** v1.8.1 and v2.5.0 on each of the
   3 cohorts. Smallest spend that delivers the §2.3 claim on cell
   lines.
2. **Update per-cohort figures** to render a side-by-side comparison
   of v1.8.1 vs v2.5.0.
3. **Decision gate.** If endpoint deltas are large and monotonic,
   skip the intermediate versions; the §2.3 claim is already made.
4. **(Conditional) Full sweep.** Add v2.1.0 / v2.2.0 / v2.3.2 reruns
   only if a reviewer or co-author requests the version trajectory.

## Out of scope

- New cell-line datasets. The use case sticks to the three current
  cohorts.
- Re-publication of the new results to PRIDE — separate workflow.
- ProteoBench-style accuracy metrics for cell lines (no ground truth
  available).
