# PXD030304 two-way reanalysis figure (ProCan-DepMapSanger / quantmsdiann)

**Status:** Approved (design, 2026-05-19)
**Date:** 2026-05-19
**Owner:** ypriverol@gmail.com

## Goal

Produce a paper-ready figure set positioning the quantmsdiann reanalysis of
PXD030304 ("Pan-cancer proteomic map of 949 human cell lines /
ProCan-DepMapSanger") against the published 2022 analysis. Mirrors the scope
already shipped for PXD003539 (main 2-condition figure + per-condition
gene/protein supp + protein-set Venn).

## Why we don't touch the 237 GB original output

The deposited `ProCan-DepMapSanger_DIANN_output.tsv` at PRIDE is 237.7 GB
(long-format DIA-NN report, 6,864 MS runs × ~50 columns). Streaming it once
to recompute peptide / protein counts is possible (HTTP supports byte-range)
but expensive and unnecessary. The Cancer Cell paper (Gonçalves et al. 2022,
PMC9387775, doi:10.1016/j.ccell.2022.06.010) deposited small post-processed
matrices on figshare ([10.6084/m9.figshare.19345397](https://doi.org/10.6084/m9.figshare.19345397))
that give us exactly what we need:

| File | Size | Purpose |
|---|---:|---|
| `protein_matrix_8498_averaged.txt` | 94 MB | 949 cell lines × 8,498 proteins (averaged across 3 replicates) |
| `protein_matrix_6692_averaged.txt` | 64 MB | same, restricted to the 6,692 ≥2-peptide "stringent" set |
| `peptide_counts_per_protein_per_sample.txt` | 143 MB | per-sample peptide counts per protein |
| `mapping_file_averaged.txt` | 82 KB | Cell_line → SIDM → Tissue_type (28 tissues) |
| `mapping_file_replicates.txt` | 989 KB | per-MS-run mapping, includes HEK293T QC runs |

These files are the canonical published representation of the ProCan
analysis; processing them for headline counts is the apples-to-apples
comparison the authors themselves use.

## Scope

In scope:

- **Main figure** — 2 conditions × 2 metrics: ProCan-DepMapSanger 2022 vs
  quantmsdiann (DIA-NN, this work), on (a) total protein groups and (b)
  protein groups supported by ≥2 unique peptides. ProCan's paper reports
  no identified-peptide count (just library size 144,578 precursors), so we
  substitute the ≥2-peptide metric for "peptides" — it's the closest
  apples-to-apples number ProCan publishes (`6,692` stringent proteins).
  Paper-ready (no title, no footer); SVG + PDF + PNG.
- **Supp figure A: per-tissue protein counts.** Grouped bar chart over the
  28 ProCan tissue categories, ProCan vs quantmsdiann. Both sides apply the
  same per-MS-run union semantics: a protein is detected in tissue T if it
  is non-NA in at least one MS run mapped to T. ProCan side reads the
  per-replicate matrix (`protein_matrix_8498_replicates.txt`, ~519 MB,
  chunked-read) joined to `mapping_file_replicates.txt` to map each MS run
  to its tissue; the 1,064 HEK293T QC runs (Tissue_type == "Control_HEK293T")
  are excluded so only the 5,800 cancer-tissue runs feed the per-tissue
  union. quantmsdiann side reads pr_matrix.tsv joined to the PXD030304
  SDRF, identical union logic.
- **Supp figure B: per-run completeness.** Same idea as PXD003539's Supp B
  (strict 1% per-run Q.Value FDR). Reanalysis read from `diann_report.parquet`
  (33 GB; chunked-read). For the original side, peptide counts per
  (protein × sample) are in `peptide_counts_per_protein_per_sample.txt`, so
  per-run completeness on the ProCan side is the fraction of stringent-set
  proteins with ≥1 peptide in that run.
- **Supp figure C: protein-group Venn.** Sets defined at the protein /
  protein-group level: ProCan 8,498 vs quantmsdiann 9,370. UniProt-accession
  normalisation matches the existing `venn_protein_accessions` helper.
- **counts.tsv** with auditable numbers (figure values + auxiliary context
  such as the 6,692-stringent count, paper-reported library size).

Out of scope:

- Streaming or processing the 237 GB raw DIA-NN long-format TSV.
- Quantitative agreement (CV, abundance correlation) comparisons.
- Drug-response / CRISPR association reproducibility (out of scope of the
  reanalysis manuscript; those live in the figshare deposit too).

## Inputs

### Reanalysis (quantmsdiann, DIA-NN 2.5)

From `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD030304/`:

- `quant_tables/diannsummary.log` — headline numbers (parsed for the global
  q-value protein group total). 253 KB.
- `quant_tables/diann_report.pr_matrix.tsv` — 2.16 GB; peptide counts.
- `quant_tables/diann_report.pg_matrix.tsv` — 199 MB; per-sample protein
  group quants (used for per-tissue counts).
- `quant_tables/diann_report.unique_genes_matrix.tsv` — 191 MB (for
  optional gene-level cross-check).
- `quant_tables/diann_report.parquet` — 33.5 GB; used only for the
  strict-per-run-FDR Supp B if rendered (chunked-read via pyarrow).
- `sdrf/PXD030304.sdrf.tsv` — 3.4 MB; maps `comment[data file]` →
  `characteristics[cell line]` for run-to-cell-line resolution.

Counting rules (locked, mirroring PXD003539):

- Protein groups: from `diannsummary.log`,
  `Protein groups with global q-value <= 0.01: 9370`.
- Precursors: from `diannsummary.log`,
  `Target precursors at 1% global q-value: 153644`.
- Peptides: unique `Stripped.Sequence` in `pr_matrix.tsv` among rows with
  ≥1 non-NA per-run column (same logic as PXD003539).

### Original (ProCan-DepMapSanger 2022)

Hardcoded headline constants from the published Cancer Cell paper:

- `PROCAN_PROTEINS = 8498` — total proteins quantified at proteotypic
  peptides Global.Q.Value ≤ 0.01 (paper, Results section + STAR Methods).
- `PROCAN_PROTEINS_STRINGENT = 6692` — proteins with ≥2 supporting peptides.
- `PROCAN_LIBRARY_PRECURSORS = 144578` — spectral library size (NOT
  identifications; included for TSV context only, no figure bar).
- `PROCAN_LIBRARY_PROTEINS = 12487` — spectral library protein count.
- `PROCAN_MS_RUNS = 6864` — runs used in the final dataset (paper).

Per-tissue / per-sample inputs from figshare 19345397:

- `mapping_file_averaged.txt` — 949 cell lines × 6 columns (Cell_line,
  SIDM, Project_Identifier, Tissue_type, Cancer_type, Cancer_subtype). The
  canonical tissue axis (28 categories: Lung, Haematopoietic and Lymphoid,
  Skin, Central Nervous System, Breast, Large Intestine, Ovary, Bone, etc.).
- `protein_matrix_8498_averaged.txt` — first column `Project_Identifier`,
  remaining 8,498 columns are protein groups like `P37108;SRP14_HUMAN`.
  Cell present (non-NA) means "protein detected in that cell line".

### Joining the two sides

- **Reanalysis run → tissue:** SDRF maps `comment[data file]` (`*.wiff`) to
  `characteristics[cell line]`; rewrite `.wiff` → `.mzML` for DIA-NN matrix
  column lookup, then look up tissue from figshare `mapping_file_averaged.txt`
  via the cell-line name. Cell-line spellings already match between the
  PXD030304 SDRF and figshare mapping (`BC-1`, `L-363`, ...); we still
  normalise (case-insensitive, strip whitespace) for safety.
- **Original cell line → tissue:** direct lookup in `mapping_file_averaged.txt`.

### 28-tissue canonical axis

The figshare `Tissue_type` column gives the published 28 tissues:

```
Lung (187), Haematopoietic and Lymphoid (163), Skin (54),
Central Nervous System (53), Breast (50), Large Intestine (47),
Ovary (41), Bone (40), Head and Neck (39), Esophagus (35),
Kidney (32), Peripheral Nervous System (30), Pancreas (30),
Stomach (27), Soft Tissue (22), Bladder (18), Thyroid (16),
Liver (16), Cervix (13), Endometrium (10), Prostate (8),
Other tissue (4), Biliary Tract (4), Vulva (3), Testis (3),
Placenta (2), Small Intestine (1), Adrenal Gland (1).
```

Counts in parentheses are cell lines per tissue in the figshare averaged
mapping. The Supp A figure orders tissues by descending sample count.

## Outputs

All written under `analysis/figures/PXD030304/`:

- `main_comparison.{pdf,png,svg}` — 2 conditions × 2 metrics.
- `supp_proteins_per_tissue.{pdf,png,svg}` — per-tissue protein counts.
- `supp_missing_values_per_run.{pdf,png,svg}` — per-run completeness
  (optional, only rendered when the parquet is downloaded).
- `supp_venn_protein_accessions.{pdf,png,svg}` — protein-group Venn.
- `counts.tsv` — auditable totals + per-tissue numbers.

## Script architecture

Single new file `analysis/figure_pxd030304_procan_vs_quantmsdiann.py`,
mirroring the PXD003539 script's structure:

1. URL constants for the quantms-collections directory + figshare file IDs
   (figshare URLs are stable `https://ndownloader.figshare.com/files/<id>`).
2. Reuses `download_if_missing` and the `Counts`-style dataclass via import
   from the existing `figure_original_vs_quantmsdiann` module (or local
   duplication if the import surface gets ugly).
3. New parsers:
   - `parse_procan_mapping(path)` → dict[Cell_line → Tissue_type].
   - `proteins_per_tissue_procan(matrix_path, mapping_path)` →
     dict[tissue → set[protein_id]].
   - `proteins_per_tissue_quantmsdiann(pg_matrix_path, sdrf_path,
     procan_mapping_path)` → dict[tissue → set[protein_id]].
4. `main()` orchestrates downloads, counters, figure rendering, TSV writing.

Tests live under `analysis/tests/test_pxd030304_*.py` with small inline
fixtures, matching the PXD003539 testing pattern.

## Cross-checks (logged, non-gating)

- quantmsdiann `diannsummary.log` protein groups must equal **9,370**;
  precursors must equal **153,644**. Mismatch ⇒ warning.
- ProCan protein matrix row count must equal **949**; column count must be
  ≥8,498 (allows extra metadata columns).
- ProCan tissue set must have **28** distinct values matching the figshare
  mapping.

## Finding: per-replicate ≡ averaged for "any-non-NA" detection

Switching ProCan's per-tissue input from the averaged matrix
(`protein_matrix_8498_averaged.txt`) to the per-MS-run matrix
(`protein_matrix_8498_replicates.txt`) produced **identical** per-tissue
counts (e.g. Lung 8,419; Adrenal Gland 4,947). The averaged matrix is
computed as `mean(non-NA)` per cell line, so a cell value is non-NA IFF at
least one of that cell line's replicates was non-NA — semantically the
same as a per-run union. The averaging does not introduce a lift.

The residual gap to quantmsdiann therefore comes from a different
filtering strictness, not from aggregation:

- **ProCan replicates matrix**: a per-run cell is non-NA whenever the
  global 8,498-protein post-curation filter (Global.Q.Value ≤ 0.01,
  proteotypic peptides) accepted it; the per-cell value can be present
  even if the local per-run identification was marginal.
- **quantmsdiann pr_matrix**: filtered at 1% precursor + 1% PG FDR per
  cell (DIA-NN's strictest per-cell filter; `--qvalue 0.01
  --matrix-qvalue 0.01`). Many proteins that pass global 1% FDR are NA
  in specific runs because their per-run identification didn't clear
  1% per-cell.

For the smaller tissues (1–4 cell lines) the per-cell strictness gap
shows up in the bars; in the largest tissues, quantmsdiann's global gain
(9,256 union vs 8,497 union) dominates and the per-tissue ranking flips.

## Open questions / follow-ups

- Should Supp A use the 8,498 set or the 6,692 stringent set for the
  ProCan side? Default: 8,498 (paper headline). Add the 6,692 set as a
  second pair of bars only if both fit cleanly on one figure.
- Whether to render the per-run completeness supp (requires 33 GB parquet
  download). Default: render it; deletable later if it duplicates the
  PXD003539 message.
