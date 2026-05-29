# quantmsdiann manuscript — brainstorming

**Status:** living brainstorming doc
**Last edited:** 2026-05-21
**Owner:** ypriverol@gmail.com

This is the source-of-truth outline for the quantmsdiann methods paper.
Each numbered experiment listed in §4 is refined into its own spec under
[`docs/superpowers/specs/`](superpowers/specs/) and tracked through to
implementation. The README and [MANIFEST.md](../MANIFEST.md) at the repo
root carry the technical layout; this doc carries the scientific story.

---

## Table of contents

- [§1 Paper at a glance](#1-paper-at-a-glance)
- [§2 The three claims](#2-the-three-claims)
- [§3 The use case](#3-the-use-case--cell-line-atlas)
- [§4 Figures index](#4-figures-index)
- [§5 Experiments — done / needed / optional](#5-experiments)
- [§6 Open experiment specs](#6-open-experiment-specs)
- [§7 Refinements applied](#7-refinements-applied-2026-05-20)
- [§8 Repository organization](#8-repository-organization)
- [§9 Decisions taken (provisional)](#9-decisions-taken-provisional)
- [§10 Next steps](#10-next-steps)

---

## 1. Paper at a glance

A methods paper introducing **quantmsdiann**, a Nextflow DIA pipeline
that wraps DIA-NN. Audience: proteomics methods readers (Nature Methods
/ MCP / JPR style) who need to be convinced of three things before they
adopt the tool:

1. **Equivalence** — quantmsdiann's identifications and quantification
   match what a careful user gets from the DIA-NN single-node GUI.
2. **Scaling** — the cost of running it grows slower than the data, and
   the cluster envelope (CPU, memory, wallclock) is predictable.
3. **Evolution** — quantmsdiann is the successor to the original
   `quantms` reanalyses. Every new DIA-NN release becomes a new variant
   of the pipeline, so the same input data can be reanalysed against
   any DIA-NN release without re-stitching a workflow.

The paper closes with a multi-dataset **cell-line atlas reanalysis**
showing concrete ID/accuracy gains over the original `quantms`
(DIA-NN 1.8.1) reanalyses.

---

## 2. The three claims

### 2.1 Equivalence — "you get the same results as DIA-NN GUI"

**Why it matters.** A Nextflow port introduces moving parts (SDRF
parsing, empirical-library assembly, fine-tuning hooks, multi-step
orchestration). Readers must not suspect that any of those drift from
what they would have obtained on a workstation.

**Vehicle.** The four community ProteoBench DIA modules. Each was
processed through **5 DIA-NN versions** (1.8.1 / 2.1.0 / 2.2.0 / 2.3.2
/ 2.5.0) = 20 analyses. ProteoBench is the trusted external arbiter
because community DIA-NN GUI submissions exist there, computed under
exactly the same ground truth.

**Figures:** F1a, F1b, F1c — all rendered. See §4.

### 2.2 Scaling — "the cost grows slower than the data, and is predictable"

This claim splits into four sub-arguments, each its own figure:

- **F2a — Wallclock vs total input size**, independent of instrument.
- **F2b — Per-step runtime distribution**, showing which DIA-NN steps
  dominate.
- **F2c — CPU and memory per task**, so a reader can size a cluster
  request before running.
- **F2d — Cluster queue-size scaling**, showing wallclock drops
  sub-linearly with cluster width on a single large cohort.

### 2.3 Evolution — "every new DIA-NN release is a new variant of the pipeline"

The original `quantms` DIA reanalyses on PRIDE were locked to DIA-NN
1.8.1. quantmsdiann tracks DIA-NN releases, so the same dataset can be
re-emitted at 1.8.1 / 2.1.0 / 2.2.0 / 2.3.2 / 2.5.0 with no workflow
surgery. Each release brings ID and accuracy improvements.

The version-progression panel of F1a already covers the benchmark side
of this claim. To carry it onto the cell lines, we need cross-version
reruns of PXD003539 / PXD004701 / PXD030304 — see experiment #13.

---

## 3. The use case — cell-line atlas

The paper's closing chapter. Already partly built:

- `figure_pxd003539_…` — NCI-60 PCT-SWATH (Guo 2019)
- `figure_pxd004701_…` — Sun 2023
- `figure_pxd030304_…` — ProCan 2022
- `figure_combined_cell_lines_atlas.py` — combined view (5 cohorts: NCI-60, ProCan, Sun, Tognetti, Wang). PXD017199 (Tognetti 2021) and PXD041421 (Wang 2023) are atlas-only; the Tognetti paper is a mass-cytometry study and the Wang dataset is a methodological batch-effect testbed, neither carries a DIA proteomics headline to compare against in a per-cohort figure.

**Filter policy.** Every count in this repo is target-only — rows with
`CONTAM_` / `Cont_` / `ENTRAP_` / `DECOY_` / `decoy_` tokens in their
Protein.Group are excluded under a conservative filter (any prefix → row
dropped, including mixed groups). The audit TSVs ship both unfiltered
and target-only counts. Design:
[docs/superpowers/specs/2026-05-21-contaminant-filter-and-pxd041421-design.md](superpowers/specs/2026-05-21-contaminant-filter-and-pxd041421-design.md).

Narrative:

1. Per-cohort ID lift vs the original published results.
2. ID/accuracy progression across DIA-NN versions (needs experiment
   #13).
3. Combined atlas: total proteome coverage across all **five** cohorts
   (PXD003539 / PXD030304 / PXD004701 + atlas-only PXD017199 /
   PXD041421), union of protein IDs, per-tissue distribution.

---

## 4. Figures index

Canonical map of every manuscript figure to its script, SVG output,
derived TSV(s), and status.

| ID | Claim | Script | SVG | Data (TSV) | Status |
|---|---|---|---|---|---|
| **F1a** | §2.1 ProteoBench parity | [`figure_quantmsdiann_benchmarks_vs_proteobench.py`](../analysis/figure_quantmsdiann_benchmarks_vs_proteobench.py) | `quantmsdiann_benchmarks/main_diann_quantmsdiann_parity.svg` + companion `main_benchmarks_*.svg`, `supp_vs_proteobench_*min3.svg` | `quantmsdiann_benchmarks/data/{counts,counts_min3,diann_quantmsdiann_parity_*,match_category_counts,median_nr_prec_by_version}.tsv` | ✅ done |
| **F1b** | §2.1 ID vs accuracy scatter | [`figure_id_vs_epsilon.py`](../analysis/figure_id_vs_epsilon.py) (consumes [`proteobench_metrics.py`](../analysis/proteobench_metrics.py)) | `quantmsdiann_benchmarks/main_id_vs_epsilon.svg` | `quantmsdiann_benchmarks/data/id_vs_epsilon_min3.tsv` | ✅ done |
| **F1c** | §2.1 per-species log2 strips (supp) | same as F1b | `quantmsdiann_benchmarks/supplementary/supp_per_species_log2_min3.svg` | `quantmsdiann_benchmarks/data/per_species_log2_min3.tsv` | ✅ done |
| **F2a** | §2.2 wallclock vs input size | [`figure_performance_runtime.py`](../analysis/figure_performance_runtime.py) | `performance/parallelism_vs_wallclock.svg` | `performance/data/parallelism_data.tsv` | ✅ done (uses cell-line `nextflow_report.html`, available on PRIDE). Caveat: 20 of 21 points have n_runs=6, weakens the "sub-linear" fit |
| **F2b** | §2.2 per-step runtime | [`figure_performance_trace.py`](../analysis/figure_performance_trace.py) | `performance/runtime_per_step.svg` | `performance/data/runtime_per_step.tsv` | ✅ benchmarks only; needs cell-line `nextflow_trace.txt` (#11) to broaden |
| **F2c** | §2.2 memory + CPU per task | [`figure_performance_trace.py`](../analysis/figure_performance_trace.py) | `performance/resources_per_step.svg` | `performance/data/resources_per_step.tsv` | ✅ done (benchmarks only; broadens with #11 like F2b) |
| **F2d** | §2.2 Cluster-node scaling | [`figure_queue_size_sweep.py`](../analysis/figure_queue_size_sweep.py) | `performance/queue_size_sweep.svg` (PXD071075 single-cell, 5 sweep points: 10/50/100/200/300 nodes) | `performance/data/queue_size_sweep.tsv` | ✅ done ([#12](superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md)): wallclock 38h → 2.4h, monotonically decreasing, knee at 200 nodes |
| **F3** | §3 cell-line atlas | [`figure_combined_cell_lines_atlas.py`](../analysis/figure_combined_cell_lines_atlas.py) | `combined/atlas_overlap.svg` (panels A/B/C — protein-accession UpSet + headline counts + detection histogram) + `combined/atlas_distribution.svg` (panels A/B/C — per-tissue cell lines + per-tissue unique proteins + breadth-vs-depth) | `combined/data/combined_counts.tsv` (audit rows prefixed by `atlas_overlap` / `atlas_distribution` for figure-name disambiguation) | ✅ done (split into two SVGs; each figure has independent A/B/C namespace; cell-line UpSet dropped as PXD030304-dominated; EA overlap dropped as redundant with PXD003539 supp) |
| **F3-version** | §2.3 + §3 cell-line cross-version progression | [`figure_cell_line_version_progression.py`](../analysis/figure_cell_line_version_progression.py) — **consumer only** | `combined/cell_line_version_progression.svg` — rendered only when cross-version matrices are staged; no placeholder SVG if data missing | `combined/data/cell_line_version_progression.tsv` | 🧪 consumer scaffolded; data-bound on [#13 spec](superpowers/specs/2026-05-20-experiment-13-cell-line-version-sweep.md) (cluster reruns) |
| **F3.PXD003539** | §3 NCI-60 reanalysis | [`figure_original_vs_quantmsdiann.py`](../analysis/figure_original_vs_quantmsdiann.py) | `PXD003539/*.svg` | `PXD003539/data/counts.tsv` | ✅ done (single version) |
| **F3.PXD004701** | §3 Sun-2023 reanalysis | [`figure_pxd004701_sun_vs_quantmsdiann.py`](../analysis/figure_pxd004701_sun_vs_quantmsdiann.py) | `PXD004701/*.svg` | `PXD004701/data/counts.tsv` | ✅ done (single version) |
| **F3.PXD030304** | §3 ProCan-2022 reanalysis | [`figure_pxd030304_procan_vs_quantmsdiann.py`](../analysis/figure_pxd030304_procan_vs_quantmsdiann.py) | `PXD030304/*.svg` | `PXD030304/data/counts.tsv` | ✅ done (single version) |

Legend: ✅ done, 📋 spec exists (not yet implemented), 🧪 consumer code scaffolded and tested but data-bound on cluster reruns, 🔧 in flight.

---

## 5. Experiments

### 5.1 Done

| # | Experiment | Section | Figure |
|---|---|---|---|
| 1 | 4 ProteoBench × 5 DIA-NN versions | §2.1 | F1a |
| 2 | Per-step trace analysis on benchmarks | §2.2 | F2b |
| 3 | Threads-vs-wallclock scatter | §2.2 | F2a |
| 4 | Per-dataset reanalysis (3 PXDs) | §3 | F3.PXD* |
| 5 | Combined cell-lines atlas | §3 | F3 |
| 6 | Matched-cohort parameter signature analysis | §2.1 | F1a |

### 5.2 Needed for the paper

| # | Experiment | Section | Figure | Status | Spec |
|---|---|---|---|---|---|
| 8 | Compute ProteoBench epsilon/CV/ROC for our 20 pr_matrix files | §2.1 | F1b, F1c | ✅ implemented + 20/20 cached | [link](superpowers/specs/2026-05-20-experiment-08-f1b-quant-accuracy.md) |
| 9 | Annotate parity figure with **signed** epsilon | §2.1 | F1a | ✅ done — see §7 | inline in §7 |
| 10 | Memory + CPU per task from existing traces | §2.2 | F2c | ✅ implemented | [link](superpowers/specs/2026-05-20-experiment-10-f2c-resources-per-task.md) |
| 11 | Cell-line `nextflow_trace.txt` collected for F2b/F2c broadening | §2.2 | F2b, F2c | 🧪 consumer ready; data-bound on cluster reruns (F2a already broad via `nextflow_report.html`) | [link](superpowers/specs/2026-05-20-experiment-11-cell-line-traces.md) |
| 12 | Cluster-node sweep (PXD071075 single-cell) | §2.2 | F2d | ✅ clean 5-point sweep (10/50/100/200/300 nodes); monotonic; 15× speedup at 300 nodes vs 10 | [link](superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md) |
| 13 | Cell-line cross-version sweep (3 datasets × ≥2 DIA-NN versions) | §2.3 + §3 | F3-version | 🧪 consumer ready; data-bound on cluster reruns | [link](superpowers/specs/2026-05-20-experiment-13-cell-line-version-sweep.md) |

### 5.3 Optional / supplementary

| # | Experiment | Section | Cost | Why |
|---|---|---|---|---|
| 14 | `--enable_fine_tuning true` on benchmarks at v2.5.0 | supp | 4 extra runs | DL fine-tuning knob (DIA-NN ≥ 2.3.2) |
| 15 | MBR via `FINAL_QUANTIFICATION { ext.args = '--reanalyse' }` on benchmarks at v2.5.0 | supp | 4 extra runs | closes residual gap to PB DIANN-predicted + QuantUMS + MBR cohort |
| 16 | `--light_models true` (DIA-NN ≥ 2.0) on benchmarks | supp | 4 extra runs | fast-screening knob |
| 17 | Pure-predicted-library mode (`skip_empirical_assembly`) | supp / discussion | requires pipeline PR | matches the ~95 % majority cohort head-to-head |

---

## 6. Open experiment specs

Each "needed" experiment in §5.2 has a one-page spec under
[`docs/superpowers/specs/`](superpowers/specs/). The spec is created
before the implementation plan; it captures inputs, outputs, validation
criteria, and risks.

| Experiment | Spec |
|---|---|
| #8 — F1b quant-accuracy (ProteoBench API integration) | [2026-05-20-experiment-08-f1b-quant-accuracy.md](superpowers/specs/2026-05-20-experiment-08-f1b-quant-accuracy.md) |
| #10 — F2c memory + CPU per task | [2026-05-20-experiment-10-f2c-resources-per-task.md](superpowers/specs/2026-05-20-experiment-10-f2c-resources-per-task.md) |
| #11 — Cell-line traces collection | [2026-05-20-experiment-11-cell-line-traces.md](superpowers/specs/2026-05-20-experiment-11-cell-line-traces.md) |
| #12 — Queue-size scaling sweep | [2026-05-20-experiment-12-queue-size-scaling.md](superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md) |
| #13 — Cell-line cross-version sweep | [2026-05-20-experiment-13-cell-line-version-sweep.md](superpowers/specs/2026-05-20-experiment-13-cell-line-version-sweep.md) |

---

## 7. Refinements applied (2026-05-20)

These three concrete improvements to existing artefacts are done:

1. ✅ **Sign-preserving epsilon** (experiment #9).
   [`diann_quantmsdiann_parity_epsilon.tsv`](../analysis/figures/quantmsdiann_benchmarks/data/diann_quantmsdiann_parity_epsilon.tsv)
   now stores `epsilon_frac = (qm − matched) / matched` (signed);
   the parity-figure annotation prints `ε = ±X.X %`. The flipped
   narrative is now visible on the SVG: at ≥3 replicates
   quantmsdiann **leads** the matched cohort on every comparable
   dataset by +7.6 % (Module 7) / +11.7 % (PXD062685) / +18.0 %
   (PXD049412). At ≥1 replicate it appears to trail (−12 % to
   −22 %), which is the well-documented DIA-NN 1.9.1 ≥1-rep
   anomaly in community submissions — the ≥3 threshold is the
   honest comparator per the Slack-driven correction in the
   benchmarks design doc.

2. ✅ **Per-version matched cohort.** The TSV gained a `qm_version`
   column. Each (dataset, threshold) pair now carries 6 rows: one
   aggregate (`qm_version = "all"`, median across the 5 versions)
   and five per-version rows where the matched cohort is recomputed
   against that specific quantmsdiann version's signature.

3. ✅ **Library-strategy colour key on the matched supp.** Bars in
   [`supp_vs_proteobench_matched_min3.svg`](../analysis/figures/quantmsdiann_benchmarks/supp_vs_proteobench_matched_min3.svg)
   (and `_min1`) are now coloured by `predictors_library` (empirical
   / DIANN-predicted / user-defined speclib / other tool). The
   `[exact|near|far]` text annotation is preserved on each bar for
   auditability but no longer drives colour. The legend is
   data-driven — it only lists library kinds that actually appear in
   the rendered bars.

---

## 8. Repository organization

Documented in full in [README.md](../README.md) and
[MANIFEST.md](../MANIFEST.md). Key conventions:

- Each figure family has one `analysis/figure_*.py` script.
- Outputs land under `analysis/figures/<group>/`:
  - `*.svg` — paper-ready vector figures (no PDF, no PNG)
  - `data/*.tsv` — auditable derived data
  - `supplementary/*.svg` — lower-threshold / reference variants
- Input caches live under the top-level `data/` (git-ignored).
- All 169 tests under `analysis/tests/` are inline-fixture-only — no
  network.

**Cleanup log (2026-05-20):**

| Item | Status |
|---|---|
| README rewritten around manuscript figure IDs | ✅ |
| `MANIFEST.md` for datasets | ✅ |
| `venn_protein_accessions.py` kept as shared utility | ✅ |
| `_min1` supps moved to `…/supplementary/` | ✅ |
| `__pycache__` / `.pytest_cache` git-ignored | ✅ |
| Test inventory verified (116 pass) | ✅ |
| SVG-only outputs (44 PDF + PNG files deleted; 19 render fns refactored) | ✅ |
| Per-analysis `data/` subdir for TSV outputs | ✅ |
| Figure legend: dropped empty "exact param match" swatch | ✅ |
| `CITATION.cff` stub | ⏳ deferred until preprint DOI |

---

## 9. Decisions taken (provisional)

Each call below is the default direction we work to until overturned.
Reasoning is recorded so the call survives a future re-read.

1. **Cluster scaling dataset → synthetic replication first, with one
   real-cohort validation.** A 2000-file synthetic cohort built by
   replicating one of our cell-line cohorts 10× in the SDRF is cheap,
   reproducible, and lets a reader rerun. Once it shows the expected
   sub-linear curve, we validate on one real public DIA cohort of
   ≥1000 files. Methods footnote acknowledges the synthetic step.
2. **Cell-line cross-version reruns → endpoints first, full sweep
   only if compute permits.** Rerun the 3 cell-line datasets at
   **v1.8.1** (original `quantms` baseline) and **v2.5.0** (current) =
   6 runs. The intermediate three versions (2.1.0 / 2.2.0 / 2.3.2 = 9
   more runs) are an extension that lands as a supp if budget allows.
3. **Quant-accuracy scope → main paper carries ID-vs-ε scatter and a
   ROC AUC table. Per-species fold-change strips and CV panel are
   supp.**
4. **MBR / fine-tuning → in scope as a single "exposed knobs"
   supplementary section.** Four extra runs each on v2.5.0 across the
   4 benchmark datasets. Pre-empts the inevitable reviewer comment.
5. **Empirical-vs-predicted-library framing → defend empirical
   assembly on its merits.** Empirical-library assembly is a design
   choice, not an oversight: it calibrates RT/IM/MS2 from the actual
   data. Discussion lists `skip_empirical_assembly` as a roadmap
   option, not a paper experiment. We do **not** commit to a pipeline
   PR pre-submission.

---

## 10. Next steps

§7 refinements + experiments #8 / #9 / #10 / #12 are all landed.
Remaining work is **cluster-bound** — the consumer code (figure
scripts, TSV writers, tests) is in the repo and verified on empty
inputs; the experiments themselves run on the SLURM cluster.
Recommended sequence:

1. **Experiment #11** — cell-line `nextflow_trace.txt` collection.
   The smallest cluster-bound experiment. Unblocks F2b/F2c broadening
   (F2a already broadened via `nextflow_report.html`). Consumers
   auto-detect traces under `data/<PXD>/pipeline_info/nextflow_trace.txt`
   via `has_cell_line_traces()`; no further code changes needed.
2. **Experiment #13** — cell-line cross-version sweep (endpoints).
   Stage `data/<PXD>/v1_8_1/diann_report.{pr,pg}_matrix.tsv` and
   `…/v2_5_0/…` and rerun
   [`figure_cell_line_version_progression.py`](../analysis/figure_cell_line_version_progression.py).
   This is the largest open gap in §2.3 — without it the evolution
   claim is benchmarks-only.
3. **F2a re-scaling** — switch x-axis from `n_runs` to `total input
   bytes` so the "sub-linear" framing reads correctly across
   instruments (currently 20/21 points sit at n=6, so the fit is fragile).
4. **Optional supps (#14 / #15 / #16)** — DIA-NN parameter
   experiments. Standalone runs; figure scaffolding TBD until results
   exist.
