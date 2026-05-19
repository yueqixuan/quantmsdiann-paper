# quantmsdiann × ProteoBench benchmarks figure design

**Status:** Approved (design, 2026-05-19)
**Date:** 2026-05-19
**Owner:** ypriverol@gmail.com

## Goal

Show what quantmsdiann buys you in the ProteoBench setting: one pipeline,
one click, four benchmark datasets each processed through **five** DIA-NN
versions (1.8.1, 2.1.0, 2.2.0, 2.3.2, 2.5.0). The point of the figure is
two-fold:

1. **Version progression within quantmsdiann** — how identifications change
   across DIA-NN releases on the same dataset.
2. **Where we land vs the community ProteoBench submissions** — the public
   ProteoBench page for each module already collects DIA-NN, AlphaDIA,
   FragPipe, Spectronaut, etc. submissions; we overlay our quantmsdiann
   numbers on top of that distribution.

## Inputs

### quantmsdiann (this work)

PRIDE FTP root:
`https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/proteobench/quantmsdiann_results/`

Layout:
`<dataset>/<version>/quant_tables/diann_report.{pr,pg}_matrix.tsv`

- Datasets: `PXD049412`, `PXD062685`, `PXD070049`, `ProteoBench_Module_7`
- DIA-NN versions per dataset: `v1_8_1`, `v2_1_0`, `v2_2_0`, `v2_3_2`,
  `v2_5_0` (5 × 4 = 20 analyses)
- Per-analysis counters:
  - **Precursors** = row count of `diann_report.pr_matrix.tsv` (one row
    per quantified precursor at 1% precursor + protein-group FDR).
  - **Protein groups** = row count of `diann_report.pg_matrix.tsv` (one
    row per protein group at 1% precursor + protein-group FDR).
- These two matrices replace `diannsummary.log` as the headline source
  because the log format is not uniform across DIA-NN versions: 1.8.1
  prints neither headline line, 2.1.0-2.3.2 print the protein-group line
  but not the precursor line, and 2.5.0 prints both. Matrix row counts
  are version-stable and match ProteoBench's `nr_prec` definition (which
  is computed from the same pr_matrix.tsv that we submit).
- The matrices are small for these benchmarks: pr_matrix.tsv is
  4.7-14.7 MB and pg_matrix.tsv is 520-800 KB per analysis. We do not
  download `diann_report.parquet` (multi-GB) — counting matrix rows is
  enough.

### Dataset ↔ ProteoBench module mapping

Verified against `docs/available-modules/active-modules/*.md` in the
`Proteobench/ProteoBench` repo:

| quantmsdiann folder | PRIDE PXD | ProteoBench module | module_id (code) | Results repo |
|---|---|---|---|---|
| `PXD049412` | PXD049412 | Module 9 — DIA single-cell | `quant_lfq_DIA_ion_singlecell` | `Proteobench/Results_quant_ion_DIA_singlecell` |
| `PXD062685` | PXD062685 | Module 5 — DIA diaPASEF | `quant_lfq_DIA_ion_diaPASEF` | `Proteobench/Results_quant_ion_DIA_diaPASEF` |
| `PXD070049` | PXD070049 | Module 10 — DIA ZenoTOF | `quant_lfq_DIA_ion_ZenoTOF` | `Proteobench/Results_quant_ion_DIA_ZenoTOF` |
| `ProteoBench_Module_7` | (not yet on PX) | Module 7 — DIA Astral 2Th | `quant_lfq_DIA_ion_Astral` | `Proteobench/Results_quant_ion_DIA_Astral` |

Module 7's raw files are still hosted only at proteobench.cubimed.rub.de;
this is the only one without a public PXD accession — hence the literal
`ProteoBench_Module_7` folder in our FTP layout.

### ProteoBench-published submissions

Each results repo is a flat directory of `<intermediate_hash>.json`
datapoints (one per public submission). Fields we use:

- `software_name`, `software_version` — tool/version label.
- `nr_prec` — total quantified precursor ions (top-level field; equals
  the `"1"` key in the `results` dict, i.e. precursors quantified in at
  least one run; this is ProteoBench's headline precursor count, which is
  what their UI plots on the y-axis).

**Important constraint:** the ProteoBench "ion" modules report only
precursor counts. There is no `nr_prot` / protein-group field in the
ion-level datapoints. Protein groups exist only for quantmsdiann here;
the head-to-head comparison vs ProteoBench is precursors-only.

We fetch the directory listing via the GitHub contents API and then each
`.json` from `raw.githubusercontent.com`. Cached under
`data/quantmsdiann_benchmarks/proteobench/<module>.json` as a single
consolidated list-of-dicts (idempotent).

## Scope

In scope:

- **Main overview figure** — per dataset, DIA-NN version progression of
  (a) precursors and (b) protein groups. 4 datasets × 2 metrics = 8
  panels, laid out as a 4×2 grid. Grouped bar chart with version on the
  x-axis. Paper-ready (no title, no footer); SVG + PDF + PNG.
- **vs-ProteoBench supp** — per dataset, a horizontal bar chart of
  precursor counts from all ProteoBench public submissions (grey) with
  the five quantmsdiann DIA-NN versions overlaid (coloured). One panel
  per dataset = 4 panels (vertical stack). Headline message: "we
  reproduce the community range, and we offer all five DIA-NN versions
  in one pipeline". Paper-ready; SVG + PDF + PNG.
- **counts.tsv** with the long-format `(dataset, tool, version,
  precursors, proteins)` table — auditable.

Out of scope:

- Quantitative-accuracy (median |epsilon|, ROC AUC) reproduction.
  ProteoBench publishes those; we do not recompute them.
- Protein-group head-to-head with ProteoBench: not possible at this time
  because the ion-level modules don't report protein counts.
- Submitting quantmsdiann's runs back to ProteoBench. The figure-rendering
  pipeline is offline; submission is a separate step.

## Joining the two sides

- Mapping is by **(PXD or module folder) → module_id**, computed at the
  top of the script as a literal dict. The four folders we encounter on
  the PRIDE FTP map deterministically to the four modules above.
- ProteoBench `software_name` strings vary in capitalisation (`DIA-NN`,
  `DIANN`, `Diann`); we normalise to `dia-nn` for highlighting our own
  version dots against community DIA-NN submissions.

## Outputs

All under `analysis/figures/quantmsdiann_benchmarks/`:

- `main_benchmarks_overview.{pdf,png,svg}` — 4×2 grid (dataset × metric),
  DIA-NN version progression.
- `supp_vs_proteobench.{pdf,png,svg}` — 4 panels, ProteoBench-community
  vs quantmsdiann precursor counts.
- `counts.tsv` — long format
  `dataset, source, tool, version, precursors, proteins, note`.

## Script architecture

Single new file
`analysis/figure_quantmsdiann_benchmarks_vs_proteobench.py`:

1. URL constants for the PRIDE FTP root, the four datasets/versions, and
   the four GitHub results repos.
2. Reuses `download_if_missing` and `SUMMARY_LOG_PROTEIN_LINE_RE` via
   import from `analysis.figure_original_vs_quantmsdiann`. Defines a
   local `SUMMARY_LOG_PRECURSOR_LINE_RE` mirroring the one in
   `analysis.figure_pxd030304_procan_vs_quantmsdiann`.
3. New helpers:
   - `count_matrix_data_rows(path)` — TSV line count minus the header;
     used uniformly across DIA-NN versions for the precursor and
     protein-group counts.
   - `parse_diann_summary_log(path)` — `(proteins, precursors)`; kept
     for downstream callers that have a v2.5.0-style log.
   - `fetch_proteobench_module(repo, dest)` — directory listing +
     per-file fetch, cached as a single JSON list under
     `data/quantmsdiann_benchmarks/proteobench/<module>.json`.
   - `parse_proteobench_datapoints(json_path)` — yields
     `(software, version, nr_prec)` tuples.
   - `normalise_software_name(name)` — case-insensitive collapse.
4. `main()` orchestrates downloads (pr_matrix + pg_matrix per analysis +
   per-module ProteoBench JSON list), aggregates the long-format
   DataFrame, renders the two figures, writes `counts.tsv`.

Tests under `analysis/tests/test_quantmsdiann_benchmarks.py` with small
inline fixtures: log parser, ProteoBench-JSON parser, normalisation,
dataset↔module mapper. No network in tests.

## Cross-checks (logged, non-gating)

- Every dataset must have exactly 5 version directories; missing ⇒ warn,
  proceed with what's available.
- pr_matrix.tsv and pg_matrix.tsv must both fetch + have > 0 data rows
  per analysis. Zero rows ⇒ warn.
- ProteoBench fetch: ≥1 datapoint per module; zero ⇒ warn (still render
  the quantmsdiann-only side).

## Open questions

- ProteoBench publishes precursor counts only at the ion level for these
  modules. Once a protein-level module exists for any of the four
  datasets we should regenerate the supp with both metrics.
- DIA-NN versions 1.8.1 and 2.x have meaningfully different defaults
  (libfree behaviour, MBR, normalisation). The figure shows raw headline
  numbers; we deliberately do not normalise.
