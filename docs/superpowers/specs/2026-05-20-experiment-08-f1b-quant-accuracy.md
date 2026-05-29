# Experiment #8 — F1b quant-accuracy (ProteoBench API integration)

**Status:** spec (not yet implemented)
**Date:** 2026-05-20
**Owner:** ypriverol@gmail.com
**Manuscript figures:** F1b (main), F1c (optional supp)
**Brainstorm reference:** [§4 / §5.2](../../brainstorming.md#5-experiments)

## Goal

Compute ProteoBench's full quant-accuracy metric set (epsilon, CV, ROC
AUC, per-species log2 fold-change) for our 20 benchmark analyses (4
datasets × 5 DIA-NN versions), so the manuscript can make an accuracy
claim alongside the ID-count claim that F1a already supports.

## Background

Each public ProteoBench submission JSON stores per-replicate-threshold
metrics:

- `mean_abs_epsilon_global`, `median_abs_epsilon_global`
- `mean_abs_epsilon_HUMAN/YEAST/ECOLI`, `_eq_species`, `_precision_*`
- `mean_log2_empirical_HUMAN/YEAST/ECOLI`
- `CV_median`, `CV_q75/q90/q95`
- `roc_auc`, `roc_auc_directional`
- `variance_epsilon_global`

We have the same precursor matrices (`diann_report.pr_matrix.tsv`) for
each of our 20 analyses but have never run them through ProteoBench's
metric computation. The metric formulas are non-trivial (species ratio
expectations + log2 transforms + multiple-testing corrections) and live
in the [`proteobench`](https://pypi.org/project/proteobench/) Python
package.

## Inputs

- `data/quantmsdiann_benchmarks/<dataset>/<version>/diann_report.pr_matrix.tsv`
  (×20, already on disk)
- The four ProteoBench module parse-settings INI files (bundled with
  the `proteobench` package; pinned to a known release)
- The per-module FASTA used by the benchmark (already cached locally)

## Outputs

- `analysis/figures/quantmsdiann_benchmarks/main_id_vs_epsilon.svg` —
  one panel per dataset; x = `nr_prec` at ≥3-rep, y =
  `median_abs_epsilon_global` at ≥3-rep; community PB points in grey
  cloud, quantmsdiann v1.8.1 → v2.5.0 in red joined by a thin line.
- `analysis/figures/quantmsdiann_benchmarks/supplementary/supp_per_species_log2.svg`
  (F1c) — per-species log2 fold-change strip plot.
- `analysis/figures/quantmsdiann_benchmarks/data/proteobench_metrics_per_analysis.tsv`
  — long-format table `(dataset, version, threshold, metric, value)`
  for the 20 analyses.
- Cache: `data/quantmsdiann_benchmarks/proteobench_metrics/<dataset>_<version>.json`
  — one JSON per analysis with the full ProteoBench metric dict.

## Validation criteria

- For each dataset, the published `nr_prec` at ≥3-rep for a *known*
  ProteoBench DIA-NN submission must be reproducible from its
  `pr_matrix.tsv` to within ±1 % (sanity check on the parse-settings
  pin).
- quantmsdiann ID counts already extracted in F1a must match the
  `nr_prec` value computed by the ProteoBench parser on the same
  matrix to within ±1 %.

## Risks

- **ProteoBench parser version drift.** Submissions in the public
  repos span `proteobench_version` 0.10 → 0.13. We pin one recent
  release (e.g. 0.13.2) for our 20 analyses and report it; older
  community points are not strictly comparable. Methods footnote.
- **Memory.** Per-module parse runs in-process — confirm <8 GB for the
  largest matrix (PXD062685 v2.5.0 pr_matrix is ~14.7 MB; should be
  fine).

## Plan (estimate)

- ~150 LOC: new `analysis/proteobench_metrics.py` (per-module parse
  invocation, cache-or-compute, long-format writer).
- ~80 LOC: new figure script `analysis/figure_id_vs_epsilon.py`
  reusing the existing dataset/version/colour conventions.
- 6 new tests covering: cache hit/miss, parse-settings selection per
  module, long-format reshape, dataset ordering.

## Out of scope

- Computing quant accuracy for the cell-line datasets (different ground
  truth; not a ProteoBench module).
- Re-running ProteoBench's full UI (we only call the metric layer).
