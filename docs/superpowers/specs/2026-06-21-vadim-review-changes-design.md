# Design — Vadim review changes + unified rebuild script

_Date: 2026-06-21. Status: approved in principle (3 decisions confirmed); pending user review of this spec._

## Context

Vadim Demichev (author V.D., DIA-NN) reviewed the pre-submission manuscript and sent a
batch of changes. This spec consolidates them, grounded in the actual code, plus a
unified figure-reproduction script Vadim requested ("rebuild --all from original
sources"). The paper is about to be submitted; scientific correctness of the counting
rules is the highest-priority item.

Confirmed decisions (user, 2026-06-21):
1. **Implement the new filter rule fully, repo-wide.**
2. **Add PXD049412 to the single-cell figure** (already analyzed; data in repo).
3. **Recompute spatial DVP (PXD064049) from matrices** (make it script-reproducible).

Already applied before this spec (safe, additive):
- Conflict of Interest: "V.D. holds shares of Aptila Biotech. The other authors declare no competing interests."
- Acknowledgments: BMBF / MSCoreSys grant 161L0221 funding for V.D.

## Goal 1 — Counting rule (scientific correctness, IMPERATIVE)

Vadim's rule, verbatim intent:

1. **Per-run protein-group numbers** (distribution boxplots, per-run averages):
   filter on `PG.Q.Value` **only**. No run-specific precursor-level filters, no global
   filters of any kind, no contaminant/target-only filter, and **no positive-quantity
   filter — zero quantities must be counted**.
2. **Per-run precursor numbers**: filter on `Q.Value` **only**.
3. **Global numbers**: the only admissible filters are `Lib.PG.Q.Value` (protein
   groups) and `Lib.Q.Value` (precursors). No other global or run-specific filters.

### Current state (what violates the rule)

`analysis/count_report_ids.py::count_report` (lines 95–150):

| Metric | Vadim requires | Code currently does |
|---|---|---|
| per-run protein groups (`prot_avg`, `prot_complete`) | `PG.Q.Value` only, count zeros | `PG.Q.Value` + target-only contaminant filter |
| per-run precursors | `Q.Value` only | `Q.Value` + `Global.Q.Value` |
| global proteins | `Lib.PG.Q.Value` only | `Global.PG.Q.Value` + target-only |
| global precursors | `Lib.Q.Value` only | `Q.Value` + `Global.Q.Value` |

The contaminant/target-only filter is the previously-approved "Option A" work; Vadim's
rule explicitly forbids it for these comparisons. This reverses that decision.

### Changes

- Rewrite `count_report` to emit, per report:
  - `prec_perrun` — precursors per run, `Q.Value <= q` only (q = DIA-NN per-version
    `--qvalue`: 0.01 for 1.8.1, 0.05 for >= 2.5.0, unchanged).
  - `prot_perrun_avg`, `prot_complete` — protein groups per run / in every run,
    `PG.Q.Value <= 0.01` only, **no target filter, zeros counted** (do not drop rows
    with zero/empty `PG.Quantity`).
  - `prec_global` — global precursors, `Lib.Q.Value <= 0.01` only.
  - `prot_global` — global protein groups, `Lib.PG.Q.Value <= 0.01` only.
- Drop the `_tgt`/`_unf` dual columns and the `is_target` contaminant logic from the
  comparison path. (Keep the contaminant helper in `contaminant_filter.py` for any
  non-comparison audit use, but it is no longer applied to reported numbers.)
- `Decoy == 0` is retained (decoys are not real identifications; this is not one of the
  forbidden "filters" — it removes the FDR null model, standard for every count).
- **plexDIA** (per-channel counting, MSV000093870): the rule **does** extend here
  (user decision 2026-06-21). Strip per-channel precursor counts to `Q.Value` only — drop
  the `Channel.Q.Value` gate that the plexDIA path currently applies. This will move the
  plexDIA numbers; recompute and update Fig 4 panel + prose. The
  `plexdia-channel-fdr` memory is now superseded and must be retired.

### Feasibility / where it runs — NO cluster needed

- Reports carry the needed columns (`PG.Q.Value`, `Q.Value`, `Lib.PG.Q.Value`,
  `Lib.Q.Value` all present — verified in the local single-cell report header).
- **All full reports are on the public PRIDE FTP** (verified 2026-06-21). The benchmark
  cohorts' `diann_report.parquet` is at
  `quantmsdiann-benchmarks/proteobench/<dataset>/<version>/quant_tables/diann_report.parquet`
  (note: real path has NO `quantmsdiann_results` segment). Single-cell, cell-line,
  spatial, and phospho trees likewise hold `quant_tables/diann_report.parquet`.
- Therefore the recount runs **locally**: download each report from FTP (cache under
  `data/`) and run the rewritten `count_report_ids` to regenerate
  `data/quantmsdiann_benchmarks/report_counts.tsv` under the new rule. This replaces the
  one-time cluster stage that produced the current staged file, and folds directly into
  the `rebuild` data-prep stage (Goal 6). Some reports are multi-GB (e.g. ProCan), so the
  counter streams/reads only needed columns, as it already does.
- **Original deposited side**: only `*_pg_matrix.tsv` / `*_pr_matrix.tsv` exist (no
  q-value columns — they are already-filtered count matrices). These are counted as
  matrix rows as-is; the rule applies to OUR reanalysis side only. Captions must say so.

### Propagation

Every per-run / global number that changes flows into: Fig 2 (benchmarks), Fig 3
(reanalysis-recovery), Fig 4 (single-cell), supplementary figures, and all numeric
claims in abstract/results/captions/Table-class numbers. Each regenerated figure's
headline numbers get re-read into the prose. A before/after numbers table is produced as
an audit artifact (`docs/.../numbers_old_vs_new.md`) even though the user approved the
full (non-staged) path — it documents the change for the submission record.

## Goal 2 — Remove One-Tip (PXD044991), full purge

One-Tip is woven through:
- `analysis/figure_single_cell_combined.py` — panel A box/jitter group, panel B
  completeness (dashed series), panel C CV; the `ACC`/`ACC_SHORT`/`DS_STYLE` dicts and
  the solid/dashed legend.
- `analysis/figure_reanalysis_improvement.py` — one PG bar + one precursor bar.
- `analysis/make_single_cell_tables.py` — `ACC`/`FTP_DIR` entries; generated
  `mv_*.tsv`, `sc_totals.tsv` rows.
- `analysis/figure_performance_trace.py:600` — `("PXD044991", "Orbitrap Astral", 12)`.
- Prose: abstract (main.tex:75), results (149), Supplementary Note 1 (117), caption
  text (218–233, 247–254); Table S1 (supplementary.tex:92).

Changes: delete every One-Tip code path, dict entry, generated row, prose mention, and
table row. Fig single-cell panels become Astral + plexDIA + PXD049412 (see Goal 4).
Recompute the abstract headline gains (avg ~30%, +17–54%, ~32% precursors) since
dropping One-Tip and re-filtering will move them.

## Goal 3 — Figure changes

- **Swap Fig 3 <-> Fig 4 order**: reanalysis-recovery (`fig_reanalysis_improvement`)
  becomes **Fig 3** (main message); single-cell (`fig3_single_cell_combined`) becomes
  **Fig 4** (niche). Reorder the two `\begin{figure}` floats in main.tex (lines
  213–262), swap `\label`s consistently (`fig:reanalysis` vs `fig:single-cell`), update
  every `\ref` and the "Fig. 3/4" prose, rename the SVG/PDF targets if we want filenames
  to match numbering (optional — labels carry the reference, filenames can stay).
- **Single-cell panel A**: per-run protein groups recomputed `PG.Q.Value`-only; add
  **percent-change labels on top** of each version pair.
- **Single-cell panel B** (completeness): y-axis starts at **0.0**.
- **Reanalysis-recovery caption (now Fig 3)**: state explicitly that the 1.8.1 side
  **is** a quantmsdiann reanalysis (not the deposited original), and how each side is
  counted (reanalysis = report under the new rule; original = deposited matrix rows).
  Clarify panel (b) "original" provenance.

## Goal 4 — Single-cell representation (PXD049412) — BLOCKED on data

After One-Tip removal single-cell = HeLa Astral + plexDIA only. Intent: add **PXD049412**
to restore breadth.

**Blocker found (2026-06-21):** the PXD049412 data in this repo (and on the FTP under
`proteobench/PXD049412/`) is **only the ProteoBench Module 9 HYE benchmark subset** — 6
runs, the `200pg_50pg_H_Y` / `240pg_10pg_H_Y` dilution series. The repo SDRF
(`PXD049412.sdrf.tsv`, 6 rows) annotates only these. The **A549/H460 single cells**,
**blanks**, and **10x/20x carriers** Vadim refers to are in the *original* PXD049412
deposit but were **never reanalyzed** by quantmsdiann. So there is currently no per-cell
A549/H460 report to plot.

Sampling rule (to apply *if/when* the single-cell runs are reanalyzed): count only
A549/H460 **single cells**; **exclude blanks** and **exclude the 10x and 20x carrier**
runs (carriers are not cells). Numbers reported per cell, not per carrier.

**Decision (2026-06-21): reanalyze the PXD049412 single-cell runs.** This is a
prerequisite pipeline run on the raw files (A549/H460 single cells), the one piece NOT
available from the FTP — it goes through the quantmsdiann workflow (cluster; codon-cluster
skill), producing a `diann_report.parquet` that is then counted under the new rule and
added to Fig 4 (single-cell) + the reanalysis-recovery bars. Run selection for this
reanalysis and for counting: **A549 and H460 single cells only; exclude blank runs;
exclude 10x and 20x carrier runs** (carriers are not cells; numbers are per cell). The
SDRF must be (re)built to annotate the single-cell runs (the current 6-row SDRF only
covers the Module 9 HYE benchmark subset). Sequencing: this reanalysis blocks the
PXD049412 single-cell addition but NOT the rest of the changes — the filter rule,
One-Tip removal, figure swap, spatial DVP, and rebuild script proceed independently.

## Goal 5 — Spatial DVP (PXD064049) reproducibility

Numbers (2,573 / 2,838 / 13,147 / 19,780) are hand-staged in
`analysis/figures/PXD064049/data/cache/*_matrix.tsv` with no computing script. Write a
real counter (extend `figure_pxd064049_spatial_vs_quantmsdiann.py` or a small module)
that computes them from the cached pg/pr matrices under the new rule, so they are
script-reproducible. If the recomputed numbers differ from the staged ones, the
recomputed (reproducible) numbers win and prose updates to match.

## Goal 6 — Unified `rebuild --all` script

`scripts/rebuild.py` (+ thin `scripts/rebuild` shim) — one entry point reproducing every
figure from PRIDE-FTP deposited reports / published matrices (user's confirmed scope),
then building the PDFs.

Stages, dependency-ordered:
1. **Env check** — `numpy<2`, `rsvg-convert`, required Python deps; fail with the conda
   hint if missing.
2. **Data prep** — `make_single_cell_tables`, `make_phospho_tables`,
   `compute_cohort_consistency` (pull deposited reports from FTP into `data/`/caches).
3. **Figures** — every `analysis/figure_*.py` + `analysis/plexDIA/figure_*` in order;
   the atlas (Fig S13) runs under the numpy<2 env automatically.
4. **Assemble** — `cd paper && make figures && make pdf && make supplementary`.

Interface: `rebuild --all | --figures-only | --only <name> | --list | --no-pdf`. Each
stage logs figure->source provenance (consolidating the scattered docstring provenance
into one manifest). Idempotent: re-downloads from FTP only when caches are absent. The
upstream raw-file DIA-NN reanalysis is explicitly out of scope (documented as a separate
HPC step, not part of `rebuild`).

## Goal 7 — `methods.md` (living, single source of truth)

Vadim wants one human-readable markdown that captures what is otherwise only knowable by
prompting the AI. Create and **maintain** `methods.md` (repo root or `docs/`) with:

- **(i) Data inventory** — every dataset (PXD/MSV accession, class, instrument,
  acquisition, year, role: benchmark / reanalysis / scalability), what files exist where
  (FTP path + local cache path), and which DIA-NN versions were run.
- **Script inventory** — every `analysis/*.py` (and plexDIA/), one line each: what it
  computes, its inputs, its outputs (figure SVG / table TSV).
- **Data-flow specification** — raw data -> DIA-NN report (FTP) -> counter -> figure
  data TSV -> figure SVG -> PDF, traced per figure. A table mapping each manuscript
  figure/number to the exact script + data file that produces it.
- **(ii) Filter specification** — for every reported quantity, the exact filters applied
  (q-value columns + thresholds, quantity handling incl. whether zeros are counted,
  peptide-count rules, target/contaminant handling, decoy handling). This is the
  executable-truth companion to the prose; it must match the rewritten `count_report_ids`
  and Vadim's rule exactly.

`methods.md` is updated as part of every change in this spec (it is not a one-off). The
`rebuild` script references it as the provenance manifest source.

## Goal 8 — Multi-audit protocol (Vadim: "multiple audits are essential")

After the changes land and `methods.md` exists, run **three independent thorough audits**
(separate agent passes so each is uncoloured by the others):

- **(a) General audit** — data flow, scripts, reproducibility: does `rebuild --all`
  actually reproduce every figure; any orphaned/stale data; any number with no producing
  script (the PXD064049 problem class).
- **(b) Filter audit** — *only* the filter specification (ii): every reported quantity
  audited against Vadim's three rules + the plexDIA decision; flag any residual
  contaminant/target/Global.* filter or uncounted-zero.
- **(c) Paper-text audit** — `methods.md` + code vs the manuscript prose: every numeric
  claim, caption, and method statement checked for correctness and completeness against
  what the code actually does.

Each audit produces a findings list; findings are triaged and fixed; re-audit until
clean. (Implemented with subagents; not a one-shot.)

## Out of scope
- Re-running DIA-NN from raw vendor files (HPC reanalysis) — documented, not scripted in
  `rebuild`.
- Unrelated refactors of figure scripts beyond what the filter rule and One-Tip removal
  require.

## Verification
- `analysis/tests/` updated to assert the new filter rule (no target filter on per-run
  PG; `Lib.*` on global; zeros counted) and the One-Tip absence. All tests pass under the
  conda env (numpy<2).
- `rebuild --all` runs clean end-to-end; both PDFs build with no undefined refs/cites.
- Grep both PDFs to confirm no "One-Tip"/"PXD044991" remains.
- Audit table `numbers_old_vs_new.md` committed.

## Risks / open points
- The recount is fully FTP-driven and local (no cluster). It downloads multi-GB reports;
  sequencing: regenerate `report_counts.tsv` under the new rule before updating the
  prose numbers / figures that read it.
- Dropping the contaminant filter may surface contaminant protein groups in headline
  counts; this is Vadim's explicit intent (rawest comparison) but worth a sanity check
  that no cohort's counts become dominated by contaminants.
- Memory `contaminant-filter-policy` and `plexdia-channel-fdr` must be updated/retired to
  reflect the reversal so future sessions don't re-apply the old filter.
