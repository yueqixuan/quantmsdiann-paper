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

## Slack-driven correction: ≥3-replicate threshold (2026-05-19)

After a Slack exchange with the ProteoBench maintainers (Marie Locard-Paulet
and Robbe Devreese) two corrections to the supp figure were agreed:

### 1. `predictors_library` is the canonical library-strategy field

Robbe confirmed that all of his DIA-NN submissions ran "right fasta + DIA-NN
predicted library" (i.e. within DIA-NN, FASTA-derived in-silico) and that
this is the policy for ~95% of current ProteoBench DIA points. The existing
`predictors_library` field — a per-component dict `{'RT': 'DIANN', 'IM':
'DIANN', 'MS2_int': 'DIANN'}` for predicted, `None` for empirical, the
literal string `"User defined speclib"` for externally uploaded — is
therefore the canonical declaration. No new extraction is required; we
document and keep `classify_predictors_library` as the single source of
truth for the colour coding in the supp figure.

### 2. `nr_prec` at the ≥1-replicate threshold is misleading for cross-version ranking

Robbe specifically called out DIA-NN 1.9.1's headline `nr_prec`:

> "you will notice that 2.3.0 will overtake 2.2 and 1.9.1 in n precursors at
> ≥3 replicate observations. Especially the effect of 1.9.1 being so high at
> 1 replicate observation is a very weird quirk that has been observed by a
> few people."

#### Field structure (verified)

ProteoBench submissions embed per-replicate-threshold counts under
`entry['results']` keyed by replicate count as a string. The exact shape,
inspected across all four modules' submissions:

```python
entry['results'] = {
    '1': { 'nr_prec': N, 'nr_prec_HUMAN': ..., 'mean_abs_epsilon_global': ... },
    '2': { 'nr_prec': N, ... },
    '3': { 'nr_prec': N, ... },
    '4': { 'nr_prec': N, ... },
    '5': { 'nr_prec': N, ... },
    '6': { 'nr_prec': N, ... },
}
```

`entry['results']['1']['nr_prec']` equals the top-level
`entry['nr_prec']`. The `nr_prec` at key `K` counts precursors quantified in
at least `K` of the six ProteoBench sample replicates. The figure now
renders both thresholds:

| File | Threshold | Default |
|---|---|---|
| `supp_vs_proteobench_min1.{pdf,png,svg}` | ≥1 replicate (legacy) | retained for comparison |
| `supp_vs_proteobench_min3.{pdf,png,svg}` | ≥3 replicates (Slack-corrected) | **new default**; matches Robbe's recommendation |

The corresponding quantmsdiann ≥3 count is computed from
`diann_report.pr_matrix.tsv` by counting rows with non-NA intensity in ≥3 of
the 6 sample columns (`count_pr_matrix_min_replicates(path, 3)`).

#### Verified ranking flips

`analysis/figures/quantmsdiann_benchmarks/median_nr_prec_by_version.tsv`
records the per-DIA-NN-version median `nr_prec` per module under both
thresholds. The ProteoBench-side ranking flips Robbe predicted are visible:

- **PXD062685 (timsTOF SCP, ProteoBench community submissions):**
  - ≥1: 1.9.1 median = **178,954** (highest, anomalously high)
  - ≥3: 1.9.1 median = **98,211** (drops to mid-pack; 2.0/2.2 now exceed it)
- **ProteoBench Module 7 (Astral, ProteoBench community submissions):**
  - ≥1: 1.9.1 median = **175,476** (highest)
  - ≥3: 1.9.1 median = **107,622** (drops below 1.7.16, 1.9.2, 2.0, 2.1.0,
    2.2.0, 2.3.0/2.3.2)
- The same pattern is consistent with the "1.9.1 quirk" Robbe described:
  unusually permissive at ≥1, in-line with neighbours at ≥3.

#### quantmsdiann ranking under the corrected threshold

quantmsdiann's per-version progression remains monotonic (or nearly so)
under ≥3, which is the desired behaviour:

| Dataset | quantmsdiann ≥3 medians (1.8.1 → 2.5.0) |
|---|---|
| PXD049412 (Astral) | 38,442 → 37,800 → 40,183 → 40,696 → 41,325 |
| PXD062685 (timsTOF SCP) | 91,688 → 100,861 → 110,281 → 111,805 → 111,825 |
| PXD070049 (ZenoTOF 7600) | 85,869 → 91,942 → 93,138 → 95,676 → 94,602 |
| ProteoBench Module 7 (Astral) | 104,301 → 107,227 → 111,497 → 113,958 → 113,270 |

The corrected supp figure is the version we ship in the paper; the ≥1 panel
is retained as a "what people usually see on ProteoBench's UI" reference
only.
