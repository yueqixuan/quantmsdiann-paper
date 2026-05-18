# PXD003539 reanalysis figure — original vs quantmsdiann

**Status:** Approved (design)
**Date:** 2026-05-18
**Owner:** ypriverol@gmail.com

## Goal

Produce a single paper-quality figure that demonstrates the value of reanalysing
the public DIA dataset PXD003539 ("NCI60 proteome by PCT-SWATH", Guo et al.,
iScience 2019) with the quantmsdiann pipeline: more precursors and more proteins
quantified than the original analysis, at the **same FDR level (1%)**.

This is the first reanalysis demonstration in the quantmsdiann manuscript repo.

> FDR-matched comparison: we compare the OpenSWATH 1%-FDR output deposited in
> PRIDE (`feature_alignment_requant_matrix.tsv` under
> `https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD003539/`) against
> the DIA-NN 1%-FDR output from quantmsdiann. We deliberately do *not* compare
> against the paper's curated 22,554 / 3,171 numbers, which were further
> filtered by the manual DIA-expert system on top of the 1% FDR cutoff and
> would make the comparison unfair.

## Scope

In scope:

- One figure: aggregate study-wide totals (precursors and protein groups),
  original vs quantmsdiann, as a grouped bar chart, both at 1% FDR.
- One reusable Python script that downloads the DIA-NN report from the PRIDE
  FTP, computes counts under a defined rule, and renders the figure.
- A small TSV that records the underlying counts (so the figure is auditable).

Out of scope (explicitly deferred):

- Per-run / per-cell-line distributions of identifications.
- FDR-matched re-analysis of the original SWATH-expert output.
- Quantitative agreement (CV, abundance) comparisons.
- Any other PXD dataset.

## Inputs

### quantmsdiann (reanalysis, DIA-NN at 1% FDR)

Downloaded from PRIDE FTP:

```
https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/quant_tables/
```

Files used:

- `diann_report.pr_matrix.tsv` — precursor × run quantification matrix
  (filtered at 1% precursor and protein-group FDR by DIA-NN). Used for the
  precursor count.
- `diannsummary.log` — parsed for the 1% global q-value protein-group total
  (`Protein groups with global q-value <= 0.01: 6927`). This is the headline
  protein count, *not* the `pg_matrix.tsv` (which is at 5% protein q-value
  and would inflate the number unfairly).

Counting rule for quantmsdiann (locked):

- **Precursors:** number of rows in `pr_matrix.tsv` with at least one non-NA
  quantification across the 120 runs (expected ≈ 117,720, matching the log).
- **Protein groups:** parsed directly from `diannsummary.log` as the integer
  after `Protein groups with global q-value <= 0.01:` (6,927).

`pr_matrix.tsv` metadata columns (everything before per-run sample columns,
validated against the public file on 2026-05-18): `Protein.Group`,
`Protein.Ids`, `Protein.Names`, `Genes`, `First.Protein.Description`,
`Proteotypic`, `Stripped.Sequence`, `Modified.Sequence`, `Precursor.Charge`,
`Precursor.Id`. Sample columns are "anything not in the metadata set". The
matrix header is validated; missing metadata columns raise `ValueError`.

### Original analysis (OpenSWATH at 1% FDR, deposited PRIDE submission)

Downloaded from PRIDE FTP:

```
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD003539/feature_alignment_requant_matrix.tsv
```

This is the OpenSWATH + pyprophet + feature_alignment requantification matrix
deposited by Guo et al. with PXD003539 on 2020-06-01. Layout: one row per
precursor; columns are `Peptide`, `Protein`, then triplets of
`Intensity_<run>`, `RT_<run>`, `score_<run>` for each of the 120 runs (≈365
columns total, ≈50,406 rows). The file is ~184 MB.

Decoys: rows are decoys when the `Peptide` column starts with `DECOY_` *or*
the `Protein` column contains `DECOY` (case-insensitive). Decoys are excluded
from all reported counts.

Counting rule for OpenSWATH (locked):

- **Precursors:** number of target rows (non-decoy) with at least one non-NA
  `Intensity_<run>` column. Expected ≈ 48,374.
- **Protein groups:** number of unique values of the `Protein` column among
  target rows with ≥ 1 non-NA `Intensity_<run>`. The `Protein` field is a
  protein-group string of the form `<N>/sp|...|XXX_HUMAN/sp|...|YYY_HUMAN`
  where `<N>` is the group size; we count unique full strings, not split
  members. Expected ≈ 6,556 (matches Guo 2019 Results: "we identified 6,556
  protein groups").

Cross-check (logged, not gating): if either OpenSWATH count differs from the
expected value by more than ~1%, log a warning. These values are deterministic
on a fixed input file, so a mismatch likely means the file changed upstream.

### Reference citation

The original dataset is described in: Guo T, Luna A, Rajapakse VN, et al.
*Quantitative Proteome Landscape of the NCI-60 Cancer Cell Lines.* iScience.
2019;21:664–680. doi: 10.1016/j.isci.2019.10.059. PubMed: 31733513.

For context, the paper itself reports two different sets of numbers:

- Pre-curation (OpenSWATH 1% FDR): 6,556 protein groups across 48,374 target
  precursors. This is what is **deposited in PRIDE** and what this analysis
  uses as the fair-comparison baseline.
- Post-curation (DIA-expert curated subset published in figures): 22,554
  proteotypic peptides / 3,171 proteins. This is the paper's headline number
  but is much stricter than 1% FDR. It is **not** used as the baseline here;
  it is reported as auxiliary context only.

The spectral library used by the original study contained 86,209 proteotypic
peptides from 8,056 SwissProt proteins.

## Outputs

- `analysis/figures/PXD003539_original_vs_quantmsdiann.pdf` — vector figure for
  inclusion in the manuscript.
- `analysis/figures/PXD003539_original_vs_quantmsdiann.png` — raster preview.
- `analysis/figures/PXD003539_counts.tsv` — tab-separated table with the four
  counts plotted, their source, and the date the script was run.

## Figure design

- Single-panel grouped bar chart, matplotlib.
- X axis: two groups, "Precursors" and "Protein groups".
- Within each group: two bars side-by-side, "Original (OpenSWATH 1% FDR)" and
  "quantmsdiann (DIA-NN 1% FDR)".
- Y axis: count (linear). If the quantmsdiann/original ratio exceeds ~5× for
  either metric, switch to a log y-axis and note it in the caption.
- Value labels on top of each bar.
- Grey for the original bar, single accent colour for quantmsdiann.
- No grid, no chart-junk. Title omitted (figure is referenced from the paper
  caption).

## Script architecture

Single file: `analysis/figure_original_vs_quantmsdiann.py`.

Layout:

1. Constants block — PRIDE FTP URLs, output paths, cache directory
   (`data/PXD003539/`), and the `pr_matrix.tsv` metadata column list.
2. `download_if_missing(url, dest)` — fetches a file once, idempotent.
3. `count_quantified_rows(tsv_path, metadata_cols, unique_by=None)` — for
   DIA-NN matrices. Reads the matrix, treats every column *not* in
   `metadata_cols` as a per-run sample column, identifies rows with at least
   one non-NA quantification, and returns either the count of those rows
   (default) or the count of distinct `unique_by` values among them. Raises
   `ValueError` if the header is missing any expected metadata column.
4. `count_openswath_quantified(tsv_path)` — for the OpenSWATH
   `feature_alignment_requant_matrix.tsv`. Returns a `(precursors, proteins)`
   tuple of integer counts among target (non-decoy) rows with at least one
   non-NA `Intensity_<run>` column. Uses chunked reading (`pandas.read_csv`
   with `chunksize`) because the file is ~184 MB; we never need the
   intensities themselves, only whether each row has any.
5. `parse_summary_log(log_path)` — parses `diannsummary.log` for the line
   `Protein groups with global q-value <= 0.01: N` and returns N; raises
   `ValueError` if the line is not present.
6. `main()` — downloads, counts, writes TSV, renders figure.

Dependencies (`analysis/requirements.txt`): `pandas`, `matplotlib`, `requests`.

The script is runnable with `python analysis/figure_original_vs_quantmsdiann.py`
from the repo root and prints the four counts to stdout.

## Error handling

- Network: `download_if_missing` retries twice on failure and aborts with a
  clear message naming the missing URL.
- DIA-NN schema: column-name detection raises an explicit `ValueError` listing
  the missing columns if the matrix layout has changed.
- OpenSWATH schema: the matrix loader checks for `Peptide`, `Protein`, and at
  least one `Intensity_*` column, and raises `ValueError` otherwise.
- Summary-log parse: `parse_summary_log` raises `ValueError` if the expected
  protein-group total line is missing.

## Testing

Light by design (this is a one-figure analysis script, not a library):

- Unit tests cover the three counters (`count_quantified_rows` with and
  without `unique_by`, schema validation, `count_openswath_quantified` with a
  small fixture, and `parse_summary_log`).
- Run end-to-end on a fresh checkout: downloads the matrices, produces PDF +
  PNG + TSV, exits 0.
- Spot-check: DIA-NN matrix precursor count within ~1% of 117,720; DIA-NN log
  protein count = 6,927 exactly; OpenSWATH precursors within ~1% of 48,374;
  OpenSWATH protein groups within ~1% of 6,556.
- Eyeball the figure: bars labelled correctly, quantmsdiann > original for
  both metrics.

## Open questions / follow-ups

- Expected magnitudes (all at 1% FDR): precursors 48,374 → 117,720 (~2.4×);
  protein groups 6,556 → 6,927 (~1.06×). The protein ratio is small but the
  precursor ratio is meaningful. If we later want a more striking protein
  comparison, we could add per-cell-line counts in a follow-up figure (out of
  scope here).
- The 184 MB OpenSWATH matrix download is one-time; cached under
  `data/PXD003539/`. Subsequent runs are fast.
