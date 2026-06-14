# quantmsdiann manuscript — figure deck redesign

**Date:** 2026-06-13
**Status:** Draft for review

## Motivation

Two drivers:

1. **Lead with real heterogeneous single cell, not benchmarks.** The Astral single-cell
   datasets (PXD046357, PXD044991 HeLa One-Tip, PXD044991 mouse zygote/oocyte) have now
   been processed on both DIA-NN 1.8.1 and 2.5.1, and 2.5 shows a clear gain on real
   low-input data — especially in data completeness (proteins quantified in *every* cell).
   This is the story ProteoBench triplicates wash out.
2. **Reframe ProteoBench as an equivalence argument.** Instead of an internal DIA-NN
   version sweep, Figure 2 should show that quantmsdiann (distributed on the cluster)
   produces results inside the cloud of existing single-machine ProteoBench community
   submissions — i.e. "use quantmsdiann and get what one machine would." MBR is preserved
   by design: quantmsdiann runs a preliminary-analysis pass across all runs, assembles a
   combined empirical library, then does final quantification of each run against that
   shared library (`preliminary_analysis → assemble_empirical_library → final_quantification`).

## Methodology correction (applies to all ID counts)

All identification counts come from the **DIA-NN precursor report**
(`diann_report.parquet` for 2.x, `diann_report.tsv` for 1.8.1) via the `count_report`
logic in `analysis/count_report_ids.py` — **NOT** the `*_matrix.tsv` files. The matrices
are filtered at different `--qvalue` thresholds per version (0.01 for 1.8.1, 0.05 for
2.5.x), so matrix row counts are apples-to-oranges across versions. The report-based
filter is uniform: `Global.Q.Value ≤ 0.01` and `Global.PG.Q.Value ≤ 0.01`, run-specific
precursor `Q.Value ≤ 0.01` (1.8.1) / `≤ 0.05` (2.5.1), target-only (drop
Cont_/CONTAM_/ENTRAP_/DECOY_).

Headline single-cell metric: **average proteins per cell/run** (run-specific
`PG.Q.Value ≤ 0.01`, target-only — the DIA-NN-author method) and **complete-profile
proteins** (quantified in every run).

## Single-cell numbers (report-based, confirmed 2026-06-13)

| Dataset | precursors | proteins (global) | avg PG/cell | complete-profile |
|---|---|---|---|---|
| PXD046357 HeLa Astral — 1.8.1 → 2.5.1 | 19,365 → 23,286 (+20%) | 3,900 → 4,230 (+8%) | 3,582 → 4,198 (+17%) | 2,886 → 3,665 (+27%) |
| PXD044991 HeLa One-Tip — 1.8.1 → 2.5.1 | 10,627 → 15,210 (+43%) | 1,586 → 1,979 (+25%) | 691 → 906 (+31%) | 218 → 351 (+61%) |
| PXD044991 mouse zygote/oocyte | *pending re-runs into correct trees* | | | |

## Figure deck (4 main figures)

### Figure 1 — Architecture & scaling *(unchanged)*
`figure_f1_architecture_scaling.py`: workflow panel + queue-size sweep + parallelism vs
wall-clock. No change.

### Figure 2 — ProteoBench equivalence *(reframe; no new image asset)*
Reorient the outputs of `figure_quantmsdiann_benchmarks_vs_proteobench.py` so the lead
message is *distributed quantmsdiann ≈ single-machine*: quantmsdiann's marker sitting
inside the predicted-library DIA-NN community-submission cloud (precursors and, where
available, quantification accuracy). The script already fetches public ProteoBench
submissions and separates predicted- vs empirical-library DIA-NN entries — this is mostly
narrative, marker highlighting, and caption work.

Caption states the MBR mechanism explicitly (preliminary analysis → assembled empirical
library → final quantification) so the equivalence reads as by-design, with ProteoBench as
empirical confirmation. Acknowledge in one sentence that the community submissions serve as
the single-machine reference (no separate local/sequential control run).

### Figure 3 — Single-cell reanalysis *(NEW; 1 composite image asset)*
New script `analysis/figure_single_cell.py`, house style from `figure_style.py`, all
counts report-based.
- Three Astral panels, DIA-NN 1.8.1 (orange) vs 2.5.1 (blue): PXD046357 HeLa Astral,
  PXD044991 HeLa One-Tip, PXD044991 mouse zygote/oocyte. Lead metric avg proteins/cell +
  complete-profile; precursors secondary.
- One plexDIA panel: MSV000093870 deposited (Galatidou 2024) vs quantms.io reanalysis,
  reusing `analysis/plexDIA/figure_msv000093870_galatidou_vs_quantmsdiann.py` data/logic.
  This panel is intentionally on a different axis (deposited-vs-reanalysis); the caption
  sub-labels it as such.

### Figure 4 — Cell-line reanalyses + pan-cohort *(merge two existing figures)*
Merge `fig:cell-line-reanalysis` and `fig:pancohort`:
- **Remove** the plexDIA panel (now in Fig 3).
- Keep the best representative cell-line examples (from PXD003539, PXD030304, PXD064049,
  PXD049692) at good resolution.
- Add a condensed pan-cohort view (best examples) + a **deposited-vs-new evidence panel**.
- Move the full atlas and the complete per-dataset description to **supplementary**.

**Assumption (confirm on review):** the deposited-vs-new evidence panel shows, aggregated
across the cell-line datasets, the split of protein groups into **shared / deposited-only /
quantms.io-new** — i.e. what the reanalysis adds versus the original deposit.

## Supplementary changes
- Full pan-cohort atlas (`combined/atlas_main`) + complete per-dataset description moved
  from main text to supplementary.
- Existing supplementary panels retained.

## Files

**New**
- `analysis/figure_single_cell.py` — Figure 3 composite.

**Modify**
- `analysis/figure_quantmsdiann_benchmarks_vs_proteobench.py` — reframe Fig 2 lead/markers.
- The cell-line/atlas figure scripts (`figure_combined_cell_lines_atlas.py` and the
  per-dataset `figure_*_vs_quantmsdiann.py`) — Fig 4 merge + deposited-vs-new panel.
- `paper/sections/figures.tex` — new Fig 3 include; Fig 4 merge; renumber.
- `paper/sections/results.tex` — narrative for the new Fig 2 equivalence framing and the
  single-cell Fig 3; rewire figure references.
- `paper/supplementary.tex` — receive full atlas + descriptions.

## Out of scope
- plexDIA both-version pipeline run (MSV000093870 has no SDRF and needs the dedicated
  plexDIA launcher — separate task if ever wanted; the deposited-vs-reanalysis panel does
  not need it).
- Any pg_matrix-based counting.

## Data dependency
Figure 3 mouse-zygote panel is blocked on the in-flight re-runs landing in the correct
trees (1.8.1 → `_v1_8_1/` tree, 2.5.1 → default tree); recompute report-based once they
complete. HeLa Astral and One-Tip numbers are final.
