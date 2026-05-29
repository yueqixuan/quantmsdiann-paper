# Fig 2 accuracy-panel redesign — per-species fold-change + community comparison

_Date: 2026-05-29_

## Problem

Fig 2b is currently a 4-module scatter of precursors-quantified (x) vs median
`|ε|` (y) with the quantmsdiann 5-version trajectory overlaid. We established
that the **version-to-version accuracy differences are not relevant**:

| module | between-version median `|ε|` span | linear | within-version per-precursor SD |
|---|---|---|---|
| Astral | 0.016 log2 | ~1.1% | 0.54 |
| diaPASEF | 0.008 log2 | ~0.6% | 0.34 |
| single-cell | 0.013 log2 | ~0.9% | 0.63 |
| ZenoTOF | 0.007 log2 | ~0.5% | 0.36 |

The version-to-version median shift (~0.01 log2) is ~20–80× smaller than the
per-precursor scatter and ~10% of the across-community spread (Astral `|ε|`
0.15–0.33). The scatter's y-axis therefore resolves noise, and on the two
modules without community comparators (single-cell, ZenoTOF) it degenerates
into a lonely zig-zag exaggerated by axis auto-scaling. The real, robust
version effect is the +14–27% identification gain, already shown in Fig 2a.

## Design

Replace the scatter with two message-bearing panels; Fig 2 becomes three panels.

- **(a)** unchanged — `main_benchmarks_precursors` (precursors at 1% FDR per
  module × version; the identification-growth story).
- **(b) Per-species fold-change accuracy** — 2×2 grid (Astral, diaPASEF,
  single-cell, ZenoTOF). x = HYE species (Human / Yeast / E. coli); y = measured
  log2 fold-change. A dashed reference line at the **expected** log2 ratio per
  species (accuracy = distance to the line) and the **5 DIA-NN versions overlaid**
  as points per species (their tight clustering shows version-invariant
  accuracy). Carries: "accuracy is stable across versions" + "quantmsdiann
  recovers the known species ratios" (with the expected DIA ratio-compression at
  the extreme E. coli ratio shown honestly).
- **(c) quantmsdiann within the predicted-library community** — for the two
  modules with predicted-library DIA-NN comparators (Astral, diaPASEF): a box +
  strip of the community median `|ε|` distribution with the 5 quantmsdiann
  versions overlaid as a tight cluster. Carries: "quantmsdiann sits within / at
  the better end of the community." Single-cell and ZenoTOF are omitted from (c)
  (no predicted-library comparators).

## Data sources (all cached; offline)

- (b) quantmsdiann per-species measured log2: `mean_log2_empirical_<SPECIES>`
  from `data/quantmsdiann_benchmarks/proteobench_metrics/<ds>_<ver>.json`
  (existing `extract_qm_per_species_log2`).
- (b) expected ratios: validated against **`ParseSettingsBuilder(...).build_parser("DIA-NN").species_expected_ratio()`**.
  Authoritative A_vs_B ratios: HYE modules (Astral/diaPASEF/ZenoTOF) Human 1.0
  (log2 0), Yeast 2.0 (+1), E. coli 0.25 (−2); single-cell Human 1.2 (+0.263),
  Yeast 0.2 (−2.32). These **match the existing hardcoded
  `SPECIES_EXPECTED_LOG2_A_vs_B` values** — only the code *comment* was wrong
  (it claimed a "0.5 / E. coli −1" design). Keep the validated dict for offline
  rebuilds, fix the comment, and add a test asserting the dict matches
  `species_expected_ratio()`. The measured E. coli ~−1.6 vs expected −2.0 is the
  expected DIA ratio-compression toward 0.
- (c) community median `|ε|`: `median_abs_epsilon_global` per predicted-library
  DIA-NN submission from `data/quantmsdiann_benchmarks/proteobench/<module>.json`
  (reuse the predicted-library filter from the id_vs_epsilon work).

## Output & integration

- New combined SVG `main_accuracy.svg` (internal sub-labels **(b)** top 2×2,
  **(c)** bottom 1×2), replacing `main_id_vs_epsilon.svg` as Fig 2's second
  image. The old `render_id_vs_epsilon` scatter is retained for the supplement /
  removed from the main figure (decide during implementation; default: keep the
  per-species + community in the main, drop the scatter from the main figure).
- `paper/sections/figures.tex`: Fig 2 = (a) precursors + the combined accuracy
  image; rewrite the caption to describe (b) accuracy-stability/per-species and
  (c) community comparison, and to state the |ε| relevance caveat.
- Verify the figure + caption fit the 548 pt text height (Fig 2 previously
  overflowed; size panels accordingly).

## Tests

- Extend `analysis/tests/` with a smoke test that the new render writes a
  non-empty SVG from synthetic per-species + community inputs, and a unit test
  that the expected-ratio source matches ProteoBench's `species_expected_ratio()`.

## Out of scope

- Recomputing ProteoBench metrics (cache is authoritative).
- Changing the identification panel (a) or the supp vs-proteobench figures.
