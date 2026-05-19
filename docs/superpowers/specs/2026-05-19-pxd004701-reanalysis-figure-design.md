# PXD004701 two-way reanalysis figure (Sun et al. 2023 / quantmsdiann)

**Status:** Approved (design, 2026-05-19)
**Date:** 2026-05-19
**Owner:** ypriverol@gmail.com

## Goal

Produce a paper-ready figure set positioning the quantmsdiann reanalysis of
PXD004701 ("76 breast cancer cell-line proteomes by PCT-SWATH") against the
published 2023 analysis by Sun et al.
([Mol Cell Proteomics, doi:10.1016/j.mcpro.2023.100602](https://doi.org/10.1016/j.mcpro.2023.100602),
[PMC10392136](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10392136/)). Mirrors
the scope already shipped for PXD030304 (main 2-condition figure + per-tissue
supp + per-run completeness supp + counts.tsv).

## Inputs

### Reanalysis (quantmsdiann, DIA-NN)

From `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD004701/`:

- `quant_tables/diannsummary.log` — confirmed headline numbers:
  `Target precursors at 1% global q-value: 100,499`,
  `Protein groups with global q-value <= 0.01: 7,746`.
  Spectral library: 14,619 protein groups / 133,806 precursors.
- `quant_tables/diann_report.pr_matrix.tsv` — precursor matrix; peptide counts.
- `quant_tables/diann_report.pg_matrix.tsv` — protein-group matrix (per-run
  completeness).
- `quant_tables/diann_report.parquet` — 33 GB long-format report; streamed for
  the Sun-style consistency filter (see below).
- `quant_tables/diann_report.unique_genes_matrix.tsv` — gene-level cross-check.
- `sdrf/PXD004701.sdrf.tsv` — 300 rows, 76 cell lines.

Counting rules (locked, mirroring PXD030304):

- Protein groups (strict 1% FDR): from `diannsummary.log`
  `Protein groups with global q-value <= 0.01: 7,746`.
- Precursors: from `diannsummary.log`
  `Target precursors at 1% global q-value: 100,499`.
- Peptides: unique `Stripped.Sequence` in `pr_matrix.tsv` among rows with
  ≥1 non-NA per-run column (same logic as PXD003539/PXD030304).
- Protein groups (Sun-style consistency filter): see next section.

### Original (Sun et al. 2023, MCP)

Hardcoded headline constants from the paper:

- `SUN_PROTEINS = 6091` — proteins consistently identified at 1% Global.Q.Value
  with proteotypic peptides AND ≤90 % missing-value rate across all samples,
  after removing 2 PCA-outlier samples.
- `SUN_PEPTIDES = 90762` — proteotypic peptides under the same filter.
- `SUN_LIBRARY_PRECURSORS = 194899` — pan-human CAL library precursors
  ([PXD009597](https://www.ebi.ac.uk/pride/archive/projects/PXD009597); NOT
  identifications; included for TSV context only).
- `SUN_LIBRARY_PROTEINS = 10323` — pan-human CAL library SwissProt proteins.
- `SUN_TNBC = 39`, `SUN_NON_TNBC = 37` — breast cancer subtype split reported
  in the paper. The paper's 76 lines = 39 TNBC + 37 non-TNBC (no "normal-like"
  category at the paper level; non-tumorigenic mammary lines are folded into
  non-TNBC in their analysis).

### PMC supplement: not machine-readable

The PMC supplement is behind a proof-of-work CAPTCHA interstitial ("Preparing
to download…"). A direct HTTP fetch with browser UA + Referer still returns
the interstitial HTML, not the supplementary files. The exact per-cell-line
TNBC / non-TNBC assignment table is therefore not retrievable
programmatically. The script encodes a hardcoded `BC_SUBTYPES` dictionary
derived from the standard breast-cancer cell-line literature
([Heiser et al. 2012, PNAS](https://doi.org/10.1073/pnas.1018854108);
[Neve et al. 2006, Cancer Cell](https://doi.org/10.1016/j.ccr.2006.10.008);
[Lehmann et al. 2011, JCI](https://doi.org/10.1172/JCI45014);
[Cellosaurus](https://www.cellosaurus.org/) characteristics) and is the
source of truth for the per-subtype supp figure.

### TNBC mapping methodology and limits

The script's `BC_SUBTYPES` dict assigns each of the 76 SDRF cell lines to one
of: `TNBC`, `non-TNBC`, `normal-like`, `unknown`. Definitions:

- **TNBC**: ER- / PR- / HER2- as reported in the cited literature. 34 lines.
- **non-TNBC**: receptor-positive (ER+, PR+, and/or HER2+) tumorigenic. 37
  lines.
- **normal-like**: non-tumorigenic mammary-epithelial lines (184A1, 184B5,
  HBL100, MCF10A, MCF12A). 5 lines.
- **unknown**: lines where canonical receptor status cannot be resolved
  from the cited literature. 0 lines after we hardcode the most defensible
  call.

This sums to 34 + 37 + 5 = 76, exactly matching the SDRF. The paper's split
is 39 TNBC vs 37 non-TNBC; the residual disagreement (5-line difference vs
their 39 TNBC count) likely reflects:

1. The paper grouping non-tumorigenic mammary lines (184A1, 184B5, HBL100,
   MCF10A, MCF12A) into one of the two subtypes; we keep them as
   `normal-like` because they are not breast cancers.
2. Borderline calls for MDA-MB-453 (often listed TNBC despite AR+ behaviour;
   we keep TNBC), SKBR7 (ambiguous receptor status; classified non-TNBC
   because HER2 expression is weak-positive in Cellosaurus), HCC2185 (rare
   metaplastic carcinoma; classified TNBC per Lehmann 2011).
3. The paper not separately publishing the per-cell mapping list (PMC supp
   inaccessible programmatically; see "PMC supplement" above).

Documented in the supp figure caption: "subtype assignments are this paper's
classification of the 76 SDRF cell lines using standard breast-cancer
literature and may differ from Sun et al. on 1–3 borderline lines".

## Scope

In scope:

- **Main figure** — 2 conditions × 3 metrics. Sun et al. 2023 vs
  quantmsdiann (DIA-NN, this work), on (a) **Protein groups (consistency
  filter)** — 6,091 (paper) vs quantmsdiann after the Sun-style ≤90 %-missing
  filter; (b) **Proteotypic peptides** — 90,762 (paper) vs quantmsdiann
  unique `Stripped.Sequence` in `pr_matrix.tsv` after the same filter;
  (c) **Protein groups (strict 1 % FDR, no consistency)** — quantmsdiann only
  (7,746 from diannsummary.log); the Sun-side bar is the paper's pre-filter
  count of **8,952** proteins identified before the missing-value rule.
  Paper-ready (no title, no footer); SVG + PDF + PNG.
- **Supp figure A: per-subtype protein counts.** quantmsdiann-only grouped
  bar chart over the 3 BC subtypes (TNBC / non-TNBC / normal-like), each
  bar = union of `Protein.Group` across all MS runs of that subtype after
  the Sun-style consistency filter is applied within-subtype. A horizontal
  reference line at 6,091 marks the Sun et al. global consistency-filtered
  total (the paper does not report per-subtype protein counts, so a
  side-by-side per-subtype comparison is not possible).
- **Supp figure B: per-run completeness.** Per-run fraction of pg_matrix
  protein groups with non-NA value, identical to PXD030304's supp B.
  quantmsdiann only (Sun et al. publish no per-run completeness data).
- **counts.tsv** with auditable numbers (headline rows + per-subtype rows +
  inline methodology notes).

Out of scope:

- Streaming or processing PXD004701's raw .wiff files.
- TNBC-specific subtype-pattern reproduction (BL1/BL2/M/MSL/IM/LAR per
  Lehmann 2011 et seq.) — needs a per-cell expression-clustering analysis
  beyond this manuscript's scope.
- Programmatically fetching the PMC supplement (blocked behind a
  proof-of-work CAPTCHA; manual download is the workaround if ever
  needed).

## Joining the two sides

- **Reanalysis run → cell line:** SDRF `comment[data file]` (`*.wiff`)
  rewritten `.wiff` → `.mzML` for matrix-column lookup, then bare-basename
  stripping for `Run` lookup in the long-format parquet.
- **Cell line → subtype:** hardcoded `BC_SUBTYPES` dict in the script
  (covers all 76 SDRF cell lines).

## Outputs

All under `analysis/figures/PXD004701/`:

- `main_comparison.{pdf,png,svg}` — 2 conditions × 3 metrics.
- `supp_proteins_per_subtype.{pdf,png,svg}` — quantmsdiann per-subtype
  protein counts + 6,091 reference line.
- `supp_missing_values_per_run.{pdf,png,svg}` — quantmsdiann per-run
  completeness.
- `counts.tsv` — auditable totals + per-subtype numbers + methodology notes.

## Script architecture

Single new file
[`analysis/figure_pxd004701_sun_vs_quantmsdiann.py`](../../../analysis/figure_pxd004701_sun_vs_quantmsdiann.py),
mirroring the PXD030304 script's structure:

1. URL constants for the quantms-collections directory.
2. Hardcoded `BC_SUBTYPES` covering all 76 SDRF cell lines.
3. Reuses `download_if_missing` from `figure_original_vs_quantmsdiann` and
   `parse_diann_summary_log` from `figure_pxd030304_procan_vs_quantmsdiann`.
4. New parsers:
   - `parse_sdrf_data_file_to_cell_line(sdrf_path)` — same logic as the
     PXD030304 helper, kept local to PXD004701 for module isolation.
   - `proteins_per_subtype_quantmsdiann_consistency_filter(parquet_url,
     sdrf_path, subtype_dict, *, qvalue_cutoff=0.01,
     min_detection_fraction=0.10)` — pyarrow + fsspec streaming over the
     parquet's 4 needed columns, applies Sun's two-stage filter
     (`Proteotypic == 1 AND Global.Q.Value <= 0.01`, then drop any
     Protein.Group detected in <10 % of mapped samples), aggregates to one
     set per subtype, JSON-caches result.
5. `main()` orchestrates downloads, headline parsing, peptide counting,
   the streaming subtype-filter compute (cache to
   `data/PXD004701/diann_per_subtype_consistency_filter.json`), figure
   rendering, counts.tsv write.

Tests live under
[`analysis/tests/test_pxd004701_parsers.py`](../../../analysis/tests/test_pxd004701_parsers.py)
with small inline fixtures and never hit the network.

## Cross-checks (logged, non-gating)

- quantmsdiann `diannsummary.log` protein groups must equal **7,746**;
  precursors must equal **100,499**. Mismatch ⇒ warning.
- SDRF must contain **76** distinct cell lines and **300** rows.
- `BC_SUBTYPES` coverage assertion: every SDRF cell line classifies to one
  of TNBC / non-TNBC / normal-like / unknown — no silent drops.
- TNBC + non-TNBC + normal-like + unknown == 76.

## Open questions / follow-ups

- If a manual download of the PMC supplement ever happens, swap the
  hardcoded `BC_SUBTYPES` for the paper's canonical assignment and rerun
  the per-subtype figure.
- Whether to add a TNBC-subtype (BL1/BL2/M/MSL/IM/LAR) breakdown using the
  Lehmann 2011 framework. Default: skip; needs expression-clustering work
  that is out of scope here.
