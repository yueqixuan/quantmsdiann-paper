# Reproducibility & data provenance

Every number and figure in the manuscript is derived from **public** data by a
script in this repository. The large DIA-NN reports themselves are *not* checked
in (they are 10s–100s of MB each); instead each figure script either downloads
them on first run and caches them, or reads a small derived table produced by a
generator that does. This file maps every figure to its source data and the
command that regenerates it.

## Data sources

| Source | Base URL / location | What |
|---|---|---|
| **Our reanalyses** (this pipeline) | `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/` | DIA-NN reports + matrices for every dataset we reprocessed, one tree per `DIA-NN` version (`v1_8_1`, `v2_5_1`, `v2_5_1_enterprise`), under `quant_tables/`. Subtrees: `single-cell/`, `proteobench/`, `PXD071075_cluster_sizes/`, plus per-PXD dirs (spatial, phospho). |
| **Original deposits** | `https://ftp.pride.ebi.ac.uk/pride/data/archive/<YYYY>/<MM>/<PXD>/` | The authors' originally deposited results we compare against (e.g. PXD003539 OpenSWATH matrix, 2020/06; PXD064049 `DIANN_results.zip`, 2025/07). |
| **quantms collections** | `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/` | Earlier curated reanalysis collections reused for some cohort figures (e.g. ProCan `proteo-search-data/ProCan/`). |
| **Pipeline** | `https://github.com/bigbio/quantmsdiann` | The Nextflow/DIA-NN workflow. |
| **SDRF tooling** | `https://github.com/bigbio/sdrf-pipelines` | `convert-diann` (SDRF → per-run DIA-NN config). |
| **Containers** | `https://github.com/bigbio/quantms-containers` | Per-`DIA-NN`-version Singularity/Docker recipes. |

## Counting conventions (the one place the rules live)

`analysis/count_report_ids.py` is the canonical counter; `analysis/contaminant_filter.py`
is the canonical target filter. All "our reanalysis" counts follow:

- **target-only**: drop any protein group with a `CONTAM_`/`Cont_`/`ENTRAP_`/`DECOY_`/`decoy_`
  token (`is_target_protein_group`).
- **protein groups**: unique `Protein.Group` at `Global.PG.Q.Value ≤ 0.01` (cross-run union).
- **precursors**: distinct `Precursor.Id` at run-specific `Q.Value` (0.01 for 1.8.1, 0.05 for
  ≥ 2.5.x — DIA-NN's per-version `--qvalue`) **and** `Global.Q.Value ≤ 0.01`.
- counts come from the **report** (`diann_report.parquet`/`.tsv`), never the `*_matrix.tsv`
  files (those bake in a version-dependent `--matrix-spec-q` and are not comparable across versions).

## Public availability (verified 2026-06-17)

Reproducibility from **public** data depends on the reanalysis reports being
deposited under `quantmsdiann-benchmarks/`. Verified by HTTP HEAD:

| Dataset / figure | Public on FTP? |
|---|---|
| single-cell PXD046357 / PXD044991 (Fig 3) | ✅ live |
| ProteoBench modules — Module_7, PXD062685, PXD049412, PXD070049 (Fig 2) | ✅ live |
| PXD071075 cluster-size sweep (Fig 1) | ✅ live |
| NCI-60 PXD003539 (Fig 4, via `quantms-collections`) | ✅ live |
| ProCan PXD030304 (Fig 4) — `quantms-collections/absolute-expression-2.0/cell-lines/PXD030304/` (figure source; also mirrored at `quantmsdiann-benchmarks/cell-lines/PXD030304/v2_5_1/`) | ✅ live |
| plexDIA MSV000093870 (Fig 3F) — `single-cell/MSV000093870/v2_5_1/` | ✅ live |
| spatial PXD064049 (Supp) — `spatial/PXD064049/v2_5_0/` | ✅ live |
| phospho PXD049692 / PXD034128 / PXD034623 (Supp) — `phospho/<ds>/v{2_5_1,2_5_1_enterprise}/` | ✅ live |

All reanalysis reports are now deposited under `quantmsdiann-benchmarks/`
(verified by HTTP HEAD). Each figure's generator pulls its report from the FTP:
single-cell → `make_single_cell_tables.py`; phospho → `make_phospho_tables.py`;
spatial/plexDIA → the figure scripts download directly. ProCan's per-tissue
generator is the one remaining item.

## Deposition checklist (to close the ⏳ rows)

To make the pending datasets reproducible from public data, upload their DIA-NN
reports to `quantmsdiann-benchmarks/` using the **same layout as single-cell**
(which already works), i.e. one tree per DIA-NN version, report under
`quant_tables/`:

```
quantmsdiann-benchmarks/<class>/<dataset>/v<version>/quant_tables/diann_report.{parquet,tsv}
```

| Dataset | Suggested path | Versions | Files needed |
|---|---|---|---|
| ProCan PXD030304 | `cell-lines/PXD030304/v<ver>/quant_tables/` | 1_8_1, 2_5_1, 2_5_1_enterprise | `diann_report.parquet` (or `.tsv` for 1.8.1) |
| plexDIA MSV000093870 | `single-cell/MSV000093870/v<ver>/quant_tables/` | 2_5_1, 2_5_1_enterprise | `diann_report.parquet` |
| spatial PXD064049 | `spatial/PXD064049/v<ver>/quant_tables/` | 1_8_1, 2_5_1_enterprise | `diann_report.parquet` |
| phospho PXD049692, PXD034128, PXD034623 | `phospho/<PXD>/v<ver>/quant_tables/` | 2_5_1, 2_5_1_enterprise | `diann_report.parquet` **and** `diann_report.site_report.parquet` (for class-I sites) |

Once these are live, the per-dataset generators (modeled on
`analysis/make_single_cell_tables.py`) read them directly from the FTP — no
hardcoded numbers. Until then, those figures read cluster-staged derived tables.

### Current cluster source paths (to upload from)

Under `/hps/nobackup/juan/pride/reanalysis/absolute-expression/`, version trees:
`_v1_8_1/` = 1.8.1, no-prefix = 2.5.1 standard, `_v2_5_1_enterprise/` = enterprise.
Each `quant_tables/` holds `diann_report.parquet` (+ `diann_report.site_report.parquet` for phospho).

| Dataset | 2.5.1 standard (no-prefix) | enterprise (`_v2_5_1_enterprise/`) |
|---|---|---|
| ProCan PXD030304 | `cell-lines-proteomes/PXD030304/results/quant_tables/` | (same subpath under `_v2_5_1_enterprise/`) |
| plexDIA MSV000093870 | `single-cell-proteomics/MSV000093870/results-plexDIA/quant_tables/` | ″ |
| spatial PXD064049 | `spatial-proteomics/PXD064049/results-CHP212-MYCN-DVP-diaPASEF-carbamidomethyl/quant_tables/` | ″ |
| phospho PXD049692 | `phospho-proteomics/PXD049692/results-NK-fibrin-IL15-phospho-diaPASEF/quant_tables/` | ″ |
| phospho PXD034128 | `phospho-proteomics/PXD034128/results-phospho-{biological-study,highspeed-DIA}/quant_tables/` | ″ |
| phospho PXD034623 | `phospho-proteomics/PXD034623/results-M2-Galectin1-phospho-DIA/quant_tables/` | ″ |

(1.8.1 reanalyses, where they exist, are under `_v1_8_1/<same subpath>`.)

## Per-figure provenance

| Figure | Numbers | Source (public) | Regenerate | Status |
|---|---|---|---|---|
| **Fig 1** pipeline/scaling | wall-clock per node-count & per dataset | `quantmsdiann-benchmarks/PXD071075_cluster_sizes/` (queue sweep) + per-dataset `nextflow_trace.txt` / `run_metadata.json` | `python -m analysis.figure_performance_trace`, `figure_queue_size_sweep` | trace/metadata staged from FTP — generator pending |
| **Fig 2** ProteoBench | depth vs accuracy, per-version IDs | ProteoBench community submissions + our `quantmsdiann-benchmarks/proteobench/` reports | `python -m analysis.figure_proteobench_accuracy` | community JSONs staged — generator pending |
| **Fig 3** single-cell | per-cell / completeness / dynamic range / CV / totals / plexDIA | `quantmsdiann-benchmarks/single-cell/{PXD046357,PXD044991_one-tip}/v{1_8_1,2_5_1_enterprise}/quant_tables/diann_report.{tsv,parquet}`; plexDIA `MSV000093870_plexDIA/` vs Galatidou 2024 matrix (`github.com/SlavovLab/single_cell_oocyte`) | **`python -m analysis.make_single_cell_tables`** then `figure_single_cell_combined` | **fully reproducible** ✅ |
| **Fig 4** cell-line atlas | per-cohort target-only protein groups, overlaps | `data/PXD00{3539,4701},PXD030304,...` reports (from `quantms-collections`/`quantmsdiann-benchmarks`) + per-cohort report JSONs | `python -m analysis.figure_combined_cell_lines_atlas` | matrices reproducible; per-tissue/per-subtype JSONs — generator pending |
| **Supp** phospho (PXD049692) | phosphopeptides / class-I sites, deposited vs reanalysis | our `quantmsdiann-benchmarks/.../PXD049692/.../diann_report.parquet` + deposited `*_PH_Report.tsv` | `python -m analysis.figure_phospho` | `phospho_counts.tsv` staged — generator pending |
| **Supp** spatial (PXD064049) | target-only protein groups, deposited vs reanalysis | `quantmsdiann-benchmarks/PXD064049-MYCN-DVP-diaPASEF/quant_tables/` + original `2025/07/PXD064049/DIANN_results.zip` | `python -m analysis.figure_pxd064049_spatial_vs_quantmsdiann` | script reads from FTP, but report ⏳ deposition pending |

## External published baselines (cited, not recomputed)

These are numbers reported by the original studies; we do not recompute them, we
cite them and compare against our reanalysis:

- **Guo 2019** (NCI-60, PXD003539) — deposited OpenSWATH matrix + paper text.
- **Walzer 2022** (PXD003539 reanalysis) — Supplementary Table S2.
- **Gonçalves 2022 / ProCan** (PXD030304) — paper protein/peptide counts.
- **Sun 2023** (PXD004701) — paper protein/peptide/library counts.
- **Galatidou 2024** (MSV000093870 plexDIA) — published proteins × cells matrix.

Each appears as a clearly-labelled constant in the relevant figure script with a
source comment.

---
*Generators that download from the FTP cache under `data/**/cache/` (git-ignored).
Run any `analysis/figure_*.py` or `analysis/make_*` module from the repo root.*
