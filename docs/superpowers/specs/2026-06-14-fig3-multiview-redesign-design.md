# Fig 3 multi-view redesign + Fig 2 CV panel

**Date:** 2026-06-14
**Status:** Draft for review

## Motivation

Figs 2 and 3 lean almost entirely on count barplots (means of precursors/proteins).
Bars show *that* 2.5.1 Enterprise beats 1.8.1, not *how*. We have rich
within-experiment data (per-cell protein counts, per-protein cross-cell detection,
per-protein intensities, CVs) that can show the mechanism — match-between-runs /
data completeness, per-cell consistency, depth into low abundance, and quantitative
precision — far more compellingly than a bar of means.

Scope: **single-cell HeLa datasets, 1.8.1 vs 2.5.1 Enterprise** (per Vadim's
two-build showcase). Mouse-zygote stays out (KB-on-mouse caveat); the standard-2.5.1
three-way lives in its own supplement already.

## Fig 3 — single-cell multi-view (2x2), 1.8.1 vs 2.5.1 Enterprise

Flagship dataset for the detailed curves is **PXD046357 HeLa Astral** (deepest,
cleanest). Panel B also includes **PXD044991 HeLa One-Tip**.

- **A — Data-completeness curve** *(headline / MBR)*: x = number of cells (1..N),
  y = protein groups quantified in ≥ x cells; two curves (1.8.1, Enterprise).
  The whole curve lifts with 2.5; the right endpoint equals the complete-profile
  count. Data: per-(Run, Protein.Group) at PG.Q ≤ 1% (target-only) → histogram of
  "#cells each PG is seen in" → reverse-cumulative.
- **B — Per-cell distribution**: violin + jittered points of protein groups per
  cell across all cells, for **both HeLa datasets × version** (4 groups). Shows 2.5
  lifts every cell. Data: per-run unique PG at PG.Q ≤ 1% (target-only).
- **C — Rank-abundance / dynamic range**: protein groups ranked by mean intensity
  (log10), 1.8.1 vs Enterprise, HeLa Astral. Shows 2.5 extends into low abundance.
  Data: pg_matrix mean intensity per PG, sorted descending.
- **D — CV density**: distribution of per-protein CV across cells, 1.8.1 vs
  Enterprise, HeLa Astral. Shows precision (more IDs at comparable/tighter CV).
  Data: pg_matrix per-PG CV = sd/mean across cell columns (PGs in ≥3 cells).

**plexDIA** (deposited-vs-reanalysis) → moves to **supplementary** (different
comparison axis; doesn't belong in a version-impact multi-view).

Colours: 1.8.1 = Okabe-Ito sky_blue, 2.5.1 Enterprise = vermillion (house style).

## Fig 2 — add a CV panel

Keep the existing two panels (depth-vs-accuracy scatter; per-concentration
fold-change). **Add a CV panel**: community CV cloud + quantms.io marker per
version, using `CV_median` (and optionally `CV_q90`) already present in the
ProteoBench submission JSONs and our metrics cache. Same two comparable datasets
(Module 7, PXD062685), 1.8.1 vs Enterprise.

## Data / staging

All computable on the cluster from existing outputs:
- Per-run PG counts + completeness histogram → small staged TSVs.
- pg_matrix per-PG mean intensity (rank-abundance) and per-PG CV → staged as
  compact arrays/TSVs (subsample or bin the rank curve to keep files small).
- ProteoBench CV from submission JSONs (already local).

Staged under `data/single_cell/` and `data/quantmsdiann_benchmarks/`. Computed
once on the cluster (anaconda python + pyarrow), terse outputs only.

## Files
**New**
- `analysis/figure_single_cell_multiview.py` — Fig 3 (A–D). Replaces
  `figure_single_cell.py` as the main Fig 3 producer (old script retained for the
  supplementary bar/plexDIA panel, or repurposed).
- staged data: `data/single_cell/per_cell_counts.tsv`,
  `completeness.tsv`, `rank_abundance.tsv`, `cv.tsv`.

**Modify**
- `analysis/figure_proteobench_accuracy.py` — add the CV panel (3rd row/panel).
- `paper/sections/figures.tex`, `results.tex` — captions + the plexDIA→supp move.

## Out of scope
- mouse-zygote (dropped) and standard-2.5.1 in Fig 3 (the std three-way is its own
  supplement).
- Raw precursor counts as a cross-version metric (q-value confound).

## Open notes
- A/C/D on HeLa Astral flagship; B on both HeLa datasets (approved).
- Rank-abundance/CV need per-protein intensities — confirm the pg_matrix intensity
  columns are populated for both versions before finalizing those two panels.
