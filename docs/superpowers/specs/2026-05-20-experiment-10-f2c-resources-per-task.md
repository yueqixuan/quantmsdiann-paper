# Experiment #10 — F2c memory + CPU per task

**Status:** spec (not yet implemented)
**Date:** 2026-05-20
**Owner:** ypriverol@gmail.com
**Manuscript figures:** F2c
**Brainstorm reference:** [§2.2 / §5.2](../../brainstorming.md#5-experiments)

## Goal

Produce a per-step distribution of memory (peak RSS) and CPU
utilisation (%cpu) so a reader can size a cluster request before
running quantmsdiann.

## Background

`nextflow_trace.txt` already carries seven columns we read nothing of
today:

- `peak_rss`, `peak_vmem` — peak memory (bytes)
- `rss`, `vmem` — final memory at exit
- `%cpu`, `%mem` — average utilisation
- `rchar`, `wchar`, `read_bytes`, `write_bytes` — I/O totals

[`figure_performance_trace.py`](../../../analysis/figure_performance_trace.py)
parses the trace files already; we just need a second view over the
same data.

## Inputs

- `nextflow_trace.txt` from the 14 complete benchmark analyses (6 are
  truncated; F2a/F2b already document this and exclude them
  consistently).
- Optional: the same files from cell-line analyses once experiment #11
  lands. The current spec ships F2c benchmarks-only and broadens when
  #11 closes.

## Outputs

- `analysis/figures/performance/resources_per_step.svg` — two-panel
  figure: (left) box-plot of `peak_rss` per DIA-NN step, (right)
  box-plot of `%cpu` per step. Same step ordering as F2b for
  readability.
- `analysis/figures/performance/data/resources_per_step.tsv` —
  long-format `(dataset, version, step, peak_rss_bytes, pct_cpu,
  duration_s)` per task.

## Validation criteria

- Total task count after deduplication matches the count in F2b's
  `runtime_per_step.tsv` (every `COMPLETED` row should produce both a
  runtime point and a resource point).
- `peak_rss` for `FINAL_QUANTIFICATION` on PXD062685 (largest cohort)
  should land in the 20-40 GB range based on DIA-NN's documented
  envelope; values outside that range need investigating.

## Risks

- Truncated traces: the 6 known PXD062685 / PXD070049 partial traces
  carry valid `SAMPLESHEET_CHECK` + `SDRF_PARSING` resource rows; F2c
  should include those and surface the gap, not silently drop them
  (mirror F2b's behaviour).
- Unit handling: `peak_rss` in `nextflow_trace.txt` is bytes;
  human-readable axis labels need a fixed converter.

## Plan (estimate)

- ~80 LOC reusing
  [`analysis/figure_performance_trace.py`](../../../analysis/figure_performance_trace.py)
  parser, adding a new `collect_resource_rows()` + `render_resources_boxplot()`.
- 3-4 new tests on the resource-extraction layer (bytes parsing,
  empty-trace edge case, step-name unification).

## Out of scope

- I/O panels (`rchar`/`wchar`). Listed for future supp; skipped here
  because the benchmark inputs are too small to be representative.
- GPU utilisation. DIA-NN can use GPUs but our cluster runs are
  CPU-only; nothing to report.
