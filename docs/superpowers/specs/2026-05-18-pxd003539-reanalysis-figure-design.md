# PXD003539 reanalysis figure — original vs quantmsdiann

**Status:** Approved (design)
**Date:** 2026-05-18
**Owner:** ypriverol@gmail.com

## Goal

Produce a single paper-quality figure that demonstrates the value of reanalysing
the public DIA dataset PXD003539 ("NCI60 proteome by PCT-SWATH", Guo et al.,
iScience 2019) with the quantmsdiann pipeline: more peptides and more proteins
quantified than the original analysis reported.

This is the first reanalysis demonstration in the quantmsdiann manuscript repo.

> Unit note: the user framed this as "more precursors and proteins". The Guo
> 2019 paper reports numbers at the **peptide** level (22,554 proteotypic
> peptides), so for an apples-to-apples comparison we count unique peptide
> sequences from the DIA-NN `pr_matrix.tsv` (`Stripped.Sequence` column) rather
> than rows (precursors = peptide × charge × modification). The total precursor
> count is also reported in the auxiliary TSV but is not the headline figure.

## Scope

In scope:

- One figure: aggregate study-wide totals (peptides and protein groups),
  original vs quantmsdiann, as a grouped bar chart.
- One reusable Python script that downloads the DIA-NN report from the PRIDE
  FTP, computes counts under a defined rule, and renders the figure.
- A small TSV that records the underlying counts (so the figure is auditable).

Out of scope (explicitly deferred):

- Per-run / per-cell-line distributions of identifications.
- FDR-matched re-analysis of the original SWATH-expert output.
- Quantitative agreement (CV, abundance) comparisons.
- Any other PXD dataset.

## Inputs

### quantmsdiann (reanalysis)

Downloaded from PRIDE FTP:

```
https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/quant_tables/
```

Files used:

- `diann_report.pr_matrix.tsv` — precursor × run quantification matrix
  (filtered at 1% precursor and protein-group FDR by DIA-NN).
- `diann_report.pg_matrix.tsv` — protein-group × run quantification matrix
  (filtered at 5% protein q-value by DIA-NN, per `diann_config.cfg`).
- `diannsummary.log` — used to cross-check totals reported by DIA-NN
  (117,720 target precursors at 1% global q-value; 6,927 protein groups at 1%
  global q-value). Note: by config the `pr_matrix.tsv` is filtered at 1%
  precursor + 1% protein-group FDR, but `pg_matrix.tsv` is filtered at 5%
  protein q-value, so matrix protein counts can legitimately exceed 6,927.

Counting rule (locked):

- **Peptides (headline):** number of *unique* values of `Stripped.Sequence` in
  `pr_matrix.tsv` for which at least one sample column has a non-NA
  quantification across the 120 runs.
- **Protein groups (headline):** number of rows in `pg_matrix.tsv` with at
  least one non-NA quantification across the 120 runs.
- **Precursors (auxiliary, TSV only):** number of rows in `pr_matrix.tsv` with
  at least one non-NA quantification.

Metadata column anchors (locked, validated against the public files on
2026-05-18):

- `pr_matrix.tsv` metadata columns (everything before per-run sample columns):
  `Protein.Group`, `Protein.Ids`, `Protein.Names`, `Genes`,
  `First.Protein.Description`, `Proteotypic`, `Stripped.Sequence`,
  `Modified.Sequence`, `Precursor.Charge`, `Precursor.Id`.
- `pg_matrix.tsv` metadata columns: `Protein.Group`, `Protein.Names`, `Genes`,
  `First.Protein.Description`, `N.Sequences`, `N.Proteotypic.Sequences`.
- Sample columns are detected as "anything not in the metadata set". The
  script validates that the matrix header contains all expected metadata
  columns and raises a clear `ValueError` if any are missing (this is what
  catches a DIA-NN version change).

Cross-check (sanity, not gating): the matrix-derived precursor row count is
expected within ~1 % of the log total (117,720 target precursors); the
matrix-derived protein-group count is expected to be ≥ 6,927 and within ~25 %
of it (because `pg_matrix.tsv` is at 5 % protein q-value vs the log's 1 %).
The script logs the auxiliary precursor count and the two log totals into the
output TSV so the figure remains auditable.

### Original analysis (Guo et al., iScience 2019)

Source: Guo T, Luna A, Rajapakse VN, et al. *Quantitative Proteome Landscape of
the NCI-60 Cancer Cell Lines.* iScience. 2019;21:664–680.
doi: 10.1016/j.isci.2019.10.059. PubMed: 31733513.

Numbers are stored as hardcoded constants at the top of the script with an
inline comment naming the exact source.

**Baseline values (confirmed 2026-05-18 from PMC6889472 full text):**

- `ORIGINAL_PEPTIDES = 22554` — Guo 2019, Results, "DIA-Expert-Curated Results
  of the NCI-60 Proteome": *"Excluding proteins/peptides that were not
  technically reproducible resulted in 22,554 proteotypic peptides from 3,171
  proteins…"* This is the count of proteotypic peptide sequences in the final
  data matrix used by the paper for all downstream analyses.
- `ORIGINAL_PROTEINS = 3171` — same sentence as above, confirmed by the
  Introduction (*"a data matrix (120 proteomes vs. 3,171 proteins)"*) and the
  Discussion (*"3,171 proteins were included for further analyses"*).

For context, the paper also reports an initial OpenSWATH identification of
6,556 protein groups before DIA-expert curation, and a spectral library of
86,209 proteotypic peptides / 8,056 SwissProt proteins. These are *not* used
as headline baselines — the paper's reproducible quantified totals are
22,554 / 3,171.

## Outputs

- `analysis/figures/PXD003539_original_vs_quantmsdiann.pdf` — vector figure for
  inclusion in the manuscript.
- `analysis/figures/PXD003539_original_vs_quantmsdiann.png` — raster preview.
- `analysis/figures/PXD003539_counts.tsv` — tab-separated table with the four
  counts plotted, their source, and the date the script was run.

## Figure design

- Single-panel grouped bar chart, matplotlib.
- X axis: two groups, "Peptides" and "Protein groups".
- Within each group: two bars side-by-side, "Original (Guo 2019)" and
  "quantmsdiann".
- Y axis: count (linear). If the quantmsdiann/original ratio exceeds ~5× for
  either metric, switch to a log y-axis and note it in the caption.
- Value labels on top of each bar.
- Grey for the original bar, single accent colour for quantmsdiann.
- No grid, no chart-junk. Title omitted (figure is referenced from the paper
  caption).

## Script architecture

Single file: `analysis/figure_original_vs_quantmsdiann.py`.

Layout:

1. Constants block — original-paper numbers and their citations, PRIDE FTP
   URLs, output paths, cache directory (`data/PXD003539/`), and the two
   metadata-column lists for `pr_matrix.tsv` and `pg_matrix.tsv`.
2. `download_if_missing(url, dest)` — fetches a file once, idempotent.
3. `count_quantified_rows(tsv_path, metadata_cols, unique_by=None)` — reads
   the matrix, treats every column *not* in `metadata_cols` as a per-run
   sample column, identifies rows with at least one non-NA quantification,
   and returns either the count of those rows (default) or the count of
   distinct values of `unique_by` among them. Raises `ValueError` if the
   header is missing any expected metadata column.
4. `main()` — downloads, counts, writes TSV, renders figure.

Dependencies (`analysis/requirements.txt`): `pandas`, `matplotlib`, `requests`.

The script is runnable with `python analysis/figure_original_vs_quantmsdiann.py`
from the repo root and prints the four counts to stdout.

## Error handling

- Network: `download_if_missing` retries twice on failure and aborts with a
  clear message naming the missing URL.
- Schema: column-name detection raises an explicit `ValueError` listing the
  unexpected columns if the DIA-NN matrix layout has changed.
- Original numbers: the constants are confirmed and hardcoded (see "Baseline
  values" above). If they are ever changed to `None` during edits, the script
  refuses to render and prints the list of missing constants.

## Testing

Light by design (this is a one-figure analysis script, not a library):

- Run end-to-end on a fresh checkout: downloads the matrices, produces PDF +
  PNG + TSV, exits 0.
- Spot-check: matrix-derived precursor count is within ~1 % of 117,720;
  matrix-derived protein count is ≥ 6,927 and within ~25 % of it. Mismatches
  outside these bands surface in the output TSV (not a hard failure).
- Eyeball the figure: bars labelled correctly, quantmsdiann > original for
  both metrics.

## Open questions / follow-ups

- If the gap is dramatic (>10×), revisit whether log axis or a second panel
  (per-run distribution) is needed; tracked but not in scope here. With the
  confirmed numbers (22,554 peptides / 3,171 proteins original vs ~117,720
  precursors and ≥ 6,927 protein groups from quantmsdiann), the protein-group
  ratio is ~2× but the peptide ratio could exceed 5× depending on how many
  unique peptide sequences DIA-NN's `Stripped.Sequence` resolves to.
