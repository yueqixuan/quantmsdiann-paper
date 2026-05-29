# Experiment #12 — Cluster queue-size scaling

**Status:** spec (not yet implemented)
**Date:** 2026-05-20
**Owner:** ypriverol@gmail.com
**Manuscript figures:** F2d
**Brainstorm reference:** [§2.2 / §5.2](../../brainstorming.md#5-experiments)

## Goal

Show that as the cluster gives quantmsdiann more concurrent slots, the
wallclock drops sub-linearly toward a floor set by the longest
non-parallel step (`FINAL_QUANTIFICATION`, which is single-job per
dataset). This is the headline scaling claim of §2.2.

## Hypothesis

For a fixed large cohort:

`wallclock(queueSize) = base + serial_floor + parallel_total / queueSize`

— monotonically decreasing with diminishing returns past a knee
determined by the longest serial path. Above ~500 slots the wallclock
should flatten.

## Design

- **Dataset.** A 2000-file synthetic cohort built by 10× SDRF
  replication of one of our cell-line cohorts. Cheap, reproducible,
  and lets a reader rerun. (See [§9.1](../../brainstorming.md#9-decisions-taken-provisional)
  decision.) Once the synthetic curve shows the expected shape, we
  validate against one real public DIA cohort of ≥1000 files (TBD;
  surveying PRIDE for candidates is part of this spec).
- **Fix every other variable.** DIA-NN v2.5.0, same fasta, same SLURM
  partition, same per-task thread count.
- **Sweep.** `executor.queueSize ∈ {50, 100, 200, 500, 1000}` = 5
  runs.
- **Capture per run.** `nextflow_trace.txt`,
  `pipeline_info/pipeline_report.txt`, total CPU-hours, peak
  concurrent tasks (already extractable by experiment #11's parser).

## Inputs

- A 200-file source cohort (candidate: a slice of PXD030304) with a
  clean SDRF.
- An SDRF-replication script that produces a 2000-row SDRF pointing at
  the same 200 physical files (acceptable provided file deduplication
  semantics are documented in methods).
- Cluster allocation sufficient for 5 multi-hour runs.

## Outputs

- `analysis/figures/performance/queue_size_sweep.svg` — log-log
  scatter of `queueSize` vs `wallclock`; expected Amdahl-bound
  overlaid as a dashed line.
- `analysis/figures/performance/data/queue_size_sweep.tsv` — rows
  `(queueSize, wallclock_s, cpu_hours, peak_concurrent_tasks,
  serial_floor_s)`.
- Five new trace files staged under `data/queue_size_sweep/q<N>/`.

## Validation criteria

- Wallclock at `queueSize = 50` must be ≥ wallclock at `queueSize =
  500` (monotonicity check; failure indicates cluster contention or
  misconfiguration).
- CPU-hours should be approximately constant across queue sizes (±
  10 %) — confirms we are sweeping the parallelism control knob, not
  changing the workload.
- `peak_concurrent_tasks` at each queueSize should approach but not
  exceed the configured limit.

## Risks

- **Cluster cost.** 5 × multi-hour runs on a 2000-file cohort is the
  single biggest compute spend in the paper. Budget needs sign-off.
- **The "right" control knob.** Confirm that Nextflow's
  `executor.queueSize` (not `process.maxForks` or the SLURM partition
  width) is what changes effective parallelism on our SLURM profile.
  A pilot at queueSize = 100 vs 1000 (without intermediate points)
  decides this before committing to the full sweep.
- **Synthetic feel.** Replicating the same 200 files 10× may not
  reflect real I/O contention. Validation against one real ≥1000-file
  cohort mitigates this; methods footnote acknowledges the synthetic
  step.

## Plan

1. **Pilot (cheap).** Two-point run at queueSize ∈ {100, 1000} to
   confirm the knob works and the curve shape is plausible.
2. **Pick a base cohort.** Decide synthetic-replication source.
3. **Full sweep.** 5 runs as designed.
4. **Real-cohort validation.** One additional run on a real ≥1000-file
   cohort at queueSize = 500.
5. **Render figure + TSV.** New `analysis/figure_queue_size_sweep.py`
   reusing the trace parser from experiment #11.

## Out of scope

- GPU back-end scaling.
- Multi-cohort sweeps (one cohort is enough for the headline claim).
- Cross-cluster comparisons (single SLURM profile).
