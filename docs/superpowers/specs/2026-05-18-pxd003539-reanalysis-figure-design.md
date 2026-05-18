# PXD003539 three-way reanalysis figure (Guo / Walzer / quantmsdiann)

**Status:** Approved (design, revised 2026-05-18)
**Date:** 2026-05-18
**Owner:** ypriverol@gmail.com

## Goal

Produce a single paper-quality figure that positions the quantmsdiann
reanalysis of PXD003539 ("NCI60 proteome by PCT-SWATH") against the two prior
analyses published on this dataset, on peptide and protein-group counts at 1%
FDR:

1. **Guo 2019 (deposited matrix):** OpenSWATH + pyprophet, 1% FDR. Source:
   `feature_alignment_requant_matrix.tsv` under
   `https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD003539/`.
2. **Walzer 2022 (Scientific Data reanalysis):** Pan-human CAL library +
   OpenSWATH + pyprophet + MSstats top3, 1% FDR. Source: Walzer et al. 2022,
   doi:10.1038/s41597-022-01380-9, Supplementary Table S2.
3. **quantmsdiann (this work):** DIA-NN 2.5 with empirical library, 1% global
   q-value. Source: PRIDE quantms-collections directory.

This is the first reanalysis demonstration in the quantmsdiann manuscript
repo.

## Scope

In scope:

- **Main figure** — aggregate study-wide totals at 1% FDR, **two** conditions
  × two metrics (peptides and protein groups), original Guo 2019 OpenSWATH
  vs quantmsdiann. We compare against the *original submission* only here
  because that's the only set of raw intermediate files we still have
  programmatic access to (Walzer 2022's S3 bucket is decommissioned).
- **Supplementary figure A** — Walzer 2022 (E-PROT-73 in Expression Atlas) vs
  quantmsdiann at the **gene level, broken down by tissue**. The two
  pipelines are processed differently downstream (Walzer applies a gene-ID
  mapping + consistency filter; we don't), so this gene-by-tissue view is
  the only fair lens on what Walzer published.
- **Supplementary figure B** — missing values per MS run, Guo 2019 OpenSWATH
  vs quantmsdiann. Same matrices as the main figure; per-run completeness is
  expected to be substantially higher in quantmsdiann.
- One Python script that downloads inputs, computes counts, and renders all
  three figures.
- A TSV recording the underlying numbers (figure data, fully auditable).

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

Counting rules for quantmsdiann (locked):

- **Peptides:** number of unique `Stripped.Sequence` values in `pr_matrix.tsv`
  among rows with at least one non-NA quantification across the 120 runs.
- **Protein groups:** parsed directly from `diannsummary.log` as the integer
  after `Protein groups with global q-value <= 0.01:` (6,927).

Auxiliary (TSV only, not on figure): the total precursor row count (≈ 117,720,
matching the log).

`pr_matrix.tsv` metadata columns (everything before per-run sample columns,
validated against the public file on 2026-05-18): `Protein.Group`,
`Protein.Ids`, `Protein.Names`, `Genes`, `First.Protein.Description`,
`Proteotypic`, `Stripped.Sequence`, `Modified.Sequence`, `Precursor.Charge`,
`Precursor.Id`. Sample columns are "anything not in the metadata set". The
matrix header is validated; missing metadata columns raise `ValueError`.

### Guo 2019 deposited matrix (OpenSWATH at 1% FDR)

Downloaded from PRIDE FTP:

```
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD003539/feature_alignment_requant_matrix.tsv
```

This is the OpenSWATH + pyprophet + feature_alignment requantification matrix
deposited by Guo et al. with PXD003539 on 2020-06-01. Layout: one row per
precursor; columns are `Peptide`, `Protein`, then triplets of
`Intensity_<run>`, `RT_<run>`, `score_<run>` for each of the 120 runs (≈365
columns total, ≈50,406 rows). The file is ~184 MB.

The `Peptide` column has the form `<id>_<modified_sequence>_<charge>_run0`,
so a precursor row maps to a peptide sequence by stripping the leading
`<id>_` and trailing `_<charge>_run0`, then removing any `(UniMod:<n>)`
modification tags.

Decoys: rows are decoys when the `Peptide` column starts with `DECOY_` *or*
the `Protein` column contains `DECOY` (case-insensitive). Decoys are excluded
from all reported counts.

Counting rules for Guo deposited (locked):

- **Peptides:** number of unique *stripped* peptide sequences (after removing
  the `<id>_` prefix, the `_<charge>_run0` suffix, and any `(UniMod:<n>)`
  modification tags) among target rows with at least one non-NA
  `Intensity_<run>` column. Expected ≈ 40,592.
- **Protein groups:** number of unique values of the `Protein` column among
  target rows with ≥ 1 non-NA `Intensity_<run>`. The `Protein` field is a
  protein-group string of the form `<N>/sp|...|XXX_HUMAN/sp|...|YYY_HUMAN`
  where `<N>` is the group size; we count unique full strings, not split
  members. Expected ≈ 6,556 (matches Guo 2019 Results: "we identified 6,556
  protein groups").

Auxiliary: 48,374 target precursors (TSV only, not on figure).

Cross-check (logged, not gating): if any count differs from the expected
value by more than ~1%, log a warning. These values are deterministic on a
fixed input file, so a mismatch likely means the file changed upstream.

### Walzer 2022 reanalysis (CAL + OpenSWATH at 1% FDR, top3 inference)

Hardcoded constants (no download). The Walzer et al. 2022 *Scientific Data*
paper (doi:10.1038/s41597-022-01380-9) reanalysed PXD003539 using the
pan-human CAL spectral library + OpenSWATH + pyprophet + MSstats with the
'top3' protein-inference setting at 1% global peptide+protein FDR.
Supplementary Table S2 of that paper reports (confirmed 2026-05-18 from the
PMC9197839 supplementary PDF):

- `WALZER_PEPTIDES = 77014` — Supplementary Table S2, row `PXD003539`,
  column `Peptides`, 1% FDR, 'top3' inference, unfiltered.
- `WALZER_PROTEINS = 7097` — same row, column `Reanalysis proteins`, 1% FDR,
  'top3' inference, unfiltered.

For context (TSV only): with the paper's '50% per group' consistency filter
applied, the protein number drops to 6,867; with 'all' inference instead of
'top3' at 1% FDR it drops to 5,412.

Note: the Walzer paper's GitHub repo (PRIDE-reanalysis/DIA-reanalysis)
pointed to programmatically downloadable TRIC and MSstats intermediate
result files, but the underlying S3 bucket (`uk1s3.embassy.ebi.ac.uk/DIA-reanalysis`)
returns `NoSuchBucket` as of 2026-05-18. We rely on the supplementary table
numbers instead, which are stable.

### Reference citations

- **Guo 2019 (original publication):** Guo T, Luna A, Rajapakse VN, et al.
  *Quantitative Proteome Landscape of the NCI-60 Cancer Cell Lines.*
  iScience. 2019;21:664–680. doi: 10.1016/j.isci.2019.10.059. PubMed: 31733513.
- **Walzer 2022 (reanalysis):** Walzer M, García-Seisdedos D, Prakash A, et al.
  *Implementing the reuse of public DIA proteomics datasets: from the PRIDE
  database to Expression Atlas.* Sci Data. 2022;9(1):335.
  doi: 10.1038/s41597-022-01380-9. PubMed: 35701420.

For context, Guo 2019 reports two distinct numbers for the same data:

- Pre-curation (OpenSWATH 1% FDR): ≈40,592 unique peptides / 6,556 protein
  groups (≈ 48,374 target precursors). This is what is **deposited in PRIDE**
  and what this analysis uses as the Guo baseline.
- Post-curation (DIA-expert curated): 22,554 proteotypic peptides / 3,171
  proteins. This is the paper's headline number but reflects manual curation
  on top of 1% FDR. **Not** used as the baseline; reported in the TSV.

The spectral library used by Guo 2019 contained 86,209 proteotypic peptides
from 8,056 SwissProt proteins. Walzer 2022 used the pan-human CAL library
(139,449 proteotypic peptides / 10,316 proteins).

## Outputs

- `analysis/figures/PXD003539_main_comparison.pdf` / `.png` — main figure
  (2 conditions × 2 metrics: Guo 2019 OpenSWATH vs quantmsdiann; peptides
  and protein groups at 1% FDR).
- `analysis/figures/PXD003539_supp_walzer_genes_by_tissue.pdf` / `.png` —
  supplementary figure A (gene counts per tissue, Walzer E-PROT-73 vs
  quantmsdiann).
- `analysis/figures/PXD003539_supp_missing_values_per_run.pdf` / `.png` —
  supplementary figure B (per-run missing-value rate, Guo OpenSWATH vs
  quantmsdiann).
- `analysis/figures/PXD003539_counts.tsv` — tab-separated table with all
  figure values plus auxiliary context (Guo curated 22,554/3,171; Walzer
  raw Supplementary Table S2 totals; precursor counts).

## Figure designs

### Main figure: 2 conditions × 2 metrics

- Single-panel grouped bar chart.
- X axis: two groups, "Peptides" and "Protein groups".
- Within each group, two bars in this fixed order:
  1. "Guo 2019 (OpenSWATH)" — grey `#9e9e9e`
  2. "quantmsdiann (DIA-NN)" — accent blue `#1f77b4`
- Value labels on top of each bar.
- Y axis: count (linear); switch to log if max/min ratio within either
  metric > 5× (expected ratio ≈ 2.4×, comfortably linear).
- Footer: italic, "All counts at 1% FDR (Guo 2019: OpenSWATH+pyprophet; this
  work: DIA-NN)."

### Supplementary figure A: genes per tissue (Walzer vs quantmsdiann)

- Single-panel grouped bar chart.
- X axis: tissue groups (NCI-60 has 9: blood/leukemia, brain/CNS, breast,
  colon, lung, melanoma, ovary, prostate, renal). Order alphabetical.
- Within each tissue, two bars: Walzer 2022 (E-PROT-73) light-blue vs
  quantmsdiann (DIA-NN unique_genes_matrix) accent blue.
- A row is "detected" in a tissue if any sample in that tissue has a non-zero
  / non-NA value.
- Sample-to-tissue mapping:
  - Walzer (E-PROT-73): parse `E-PROT-73-configuration.xml` to map
    `g<N>` → assay (cell line) → disease label, then bucket disease labels
    into the 9 tissue categories (we'll use a small literal mapping inside
    the script).
  - quantmsdiann: parse `PXD003539.sdrf.tsv` to map `assay name` → MS run
    file → `characteristics[organism part]` (already at tissue granularity).
- Y axis: unique gene count.

### Supplementary figure B: missing values per run

- Single-panel, line / strip plot.
- X axis: MS run index, sorted by run name (the 120 files are alphabetical
  by acquisition date, so the order doubles as an acquisition timeline).
- Y axis: fraction of precursors with a non-NA quantification in that run.
- Two series:
  - Guo 2019 OpenSWATH: from `feature_alignment_requant_matrix.tsv`,
    fraction of target precursor rows with non-NA `Intensity_<run>`.
  - quantmsdiann (DIA-NN): from `diann_report.pr_matrix.tsv`, fraction of
    rows with non-NA value in each per-run column.
- Both denominators are the TOTAL number of precursor rows in the respective
  matrix (so the two series are not directly stack-comparable; what is
  meaningful is the *shape* of each curve and the per-pipeline floor).
- Optional secondary panel: missing-value *count* on the same x axis. We
  start with a single panel and add the count if the percentage view alone
  doesn't tell the story.

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
   `feature_alignment_requant_matrix.tsv`. Returns a
   `(precursors, peptides, proteins)` tuple of integer counts among target
   (non-decoy) rows with at least one non-NA `Intensity_<run>` column:
   - precursors: count of qualifying rows.
   - peptides: count of unique *stripped* peptide sequences (after removing
     the `<id>_` prefix, the `_<charge>_run0` suffix, and any `(UniMod:<n>)`
     modification tags) among qualifying rows.
   - proteins: count of unique `Protein` column values among qualifying rows.
   Uses chunked reading (`pandas.read_csv` with `chunksize`) because the file
   is ~184 MB; we never need the intensities themselves, only whether each
   row has any.
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
