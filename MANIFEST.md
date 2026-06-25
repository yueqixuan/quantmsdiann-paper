# Dataset manifest

Every dataset under `data/` with its source, instrument, file count, and
the manuscript figure(s) it feeds. Inputs are git-ignored and downloaded
on demand by each figure script; this file is the single source of truth
for what lives where.

## ProteoBench benchmark datasets (§2.1 — equivalence)

`data/quantmsdiann_benchmarks/<dataset>/v<version>/` — quantmsdiann
output per (dataset, DIA-NN version) pair. 4 datasets × 5 versions = 20
analyses. Cached artefacts per analysis: `diann_report.pr_matrix.tsv`,
`diann_report.pg_matrix.tsv`, `nextflow_trace.txt`, `pipeline_info/pipeline_report.txt`,
`quant_tables/diannsummary.log`.

| PXD / module | ProteoBench module | Instrument | Sample replicates | Source |
|---|---|---|---|---|
| `PXD049412` | Module 9 — DIA single-cell | Orbitrap Astral | 6 | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/proteobench/quantmsdiann_results/PXD049412/) |
| `PXD062685` | Module 5 — DIA diaPASEF | timsTOF SCP | 6 | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/proteobench/quantmsdiann_results/PXD062685/) |
| `PXD070049` | Module 10 — DIA ZenoTOF | ZenoTOF 7600 | 6 | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/proteobench/quantmsdiann_results/PXD070049/) |
| `ProteoBench_Module_7` | Module 7 — DIA Astral 2Th | Orbitrap Astral | 6 | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/proteobench/quantmsdiann_results/ProteoBench_Module_7/) (raw files at proteobench.cubimed.rub.de — no PXD accession) |

DIA-NN versions per dataset: `v1_8_1`, `v2_1_0`, `v2_2_0`, `v2_3_2`,
`v2_5_0`.

Consumed by: **F1a** (parity), **F1b** (planned ID vs ε scatter),
**F2a–F2b** (performance), **F2c** (planned resources panel).

## ProteoBench community submissions (§2.1 — equivalence)

`data/quantmsdiann_benchmarks/proteobench/<module>.json` — consolidated
list-of-dicts cached from the GitHub
[`Proteobench/Results_quant_ion_DIA_*`](https://github.com/Proteobench)
repositories. One JSON per module; ~50–600 KB. Idempotent fetch via
GitHub contents API.

| File | Submissions (total / DIA-NN) | GitHub repo |
|---|---|---|
| `PXD049412.json` | 12 / 10 | `Proteobench/Results_quant_ion_DIA_singlecell` |
| `PXD062685.json` | 32 / 24 | `Proteobench/Results_quant_ion_DIA_diaPASEF` |
| `PXD070049.json` | 4 / 3 | `Proteobench/Results_quant_ion_DIA_ZenoTOF` |
| `ProteoBench_Module_7.json` | 47 / 37 | `Proteobench/Results_quant_ion_DIA_Astral` |

Consumed by: **F1a**, **F1b** (planned).

## Cell-line reanalysis datasets (§3 — use case)

`data/<PXD>/` — quantmsdiann reanalysis outputs for the use-case
chapter. One DIA-NN version per dataset today; cross-version sweep is
experiment #13 in [docs/brainstorming.md](docs/brainstorming.md).

| PXD | Cohort | Instrument | Files (SDRF rows) | Original publication | quantmsdiann FTP |
|---|---|---|---|---|---|
| `PXD003539` | NCI-60 PCT-SWATH | TripleTOF 5600+ | 120 | Guo et al. 2019 (Nat. Med.) | [PRIDE](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/) |
| `PXD004701` | Sun 2023 — pan-cancer | Orbitrap | 300 | Sun et al. 2023 | [PRIDE](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD004701/) |
| `PXD030304` | ProCan-DepMapSanger 949 lines | Orbitrap | 5798 | Gonçalves et al. 2022 (Cancer Cell) | [PRIDE](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD030304/) |
| `PXD017199` | Tognetti 2021 — 67 breast lines (5-6 normal-like) | Q Exactive Plus | 206 | Tognetti et al. 2021 (Cell Systems, doi:10.1016/j.cels.2021.04.002) | [PRIDE](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD017199/) — atlas-only (no per-cohort figure: paper is mass-cytometry) |
| `PXD041421` | Wang 2023 — 2 lines (A549 lung + K562 leukemia) × 24 reps each, 2-batch design | timsTOF Pro (diaPASEF) | 48 | Wang et al. 2023 (batch-effect benchmark dataset) | [PRIDE FTP](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD041421/) — atlas-only (no per-cohort figure: methodological dataset) |

Cached artefacts vary per dataset (pr_matrix, pg_matrix, SDRF,
per-tissue protein JSONs). See the `atlas` stage and its helpers in
[`scripts/rebuild.py`](scripts/rebuild.py) (`python -m scripts.rebuild --only atlas`)
for the canonical layout.

Consumed by: **F3** (combined atlas), **F3.PXDxxx** (per-cohort
reanalysis comparisons).

## UniProt SwissProt FASTAs (§2.1 — equivalence)

`data/quantmsdiann_benchmarks/uniprot/` — SwissProt FASTA streams used by
[`proteobench_metrics.py`](analysis/proteobench_metrics.py) to rebuild a
UniProt-accession → species map. quantmsdiann's `Protein.Ids` cells are
bare accessions (`Q96P70`); ProteoBench's species-detection works by
substring match on `_HUMAN` / `_YEAST` / `_ECOLI` suffixes, so we
re-annotate accessions before passing them to ProteoBench's parser.

| File | Source | Used for |
|---|---|---|
| `human_reviewed.fasta` | `uniprot.org/uniprotkb/stream?reviewed:true+AND+organism_id:9606` | `_HUMAN` suffix annotation |
| `yeast_reviewed.fasta` | `… organism_id:559292` (S. cerevisiae S288C) | `_YEAST` suffix annotation |
| `ecoli_reviewed.fasta` | `… organism_id:83333` (E. coli K12) | `_ECOLI` suffix annotation |

Consumed by: **F1b**, **F1c**.

## ProteoBench metrics cache (§2.1 — equivalence)

`data/quantmsdiann_benchmarks/proteobench_metrics/<dataset>_<version>.json`
— per-analysis ProteoBench metric payload computed by
[`compute_proteobench_metrics`](analysis/proteobench_metrics.py). One
JSON per (dataset, DIA-NN version) pair, 20 total. Each payload carries
the same `results` schema as a published ProteoBench submission
(per-replicate-threshold dict keyed by `"1"…"6"`, with `nr_prec`,
`median_abs_epsilon_global`, `roc_auc`, `CV_median`, per-species log2
fold-change, etc.).

Consumed by: **F1b**, **F1c**.

## Auxiliary data

| File | Purpose | Consumed by |
|---|---|---|
| `data/E-PROT-73-query-results.tsv` | Expression Atlas (ArrayExpress) query output, used to cross-reference cell-line gene expression | `figure_pxd030304_procan_vs_quantmsdiann.py` |

## Caching contract

- All paths under the top-level `data/` are git-ignored (see `.gitignore`).
- Every figure script implements its own `download_if_missing(url, dest)`
  idempotent fetch: present-and-non-empty file ⇒ no network call.
- The single most important rule: **do not commit anything under
  the top-level `data/`** — those are cached inputs. The SVGs in
  `analysis/figures/<group>/` and derived TSVs in
  `analysis/figures/<group>/data/` are the only outputs that ship in
  git. (The two `data/` directories are deliberately distinct: top-level
  for input caches, per-figure for derived auditable tables.)

## Output layout reminder

```
analysis/figures/<group>/
  *.svg              # paper-ready vector figures
  data/              # auditable derived data
    *.tsv            # per-analysis tables consumed by the script and by reviewers
  supplementary/     # lower-threshold / reference variants of the same figures
```

## Nextflow queueSize sweep (§2.2 — F2d scaling)

`data/queue_size_sweep/q<N>/{nextflow_trace.txt, diannsummary.log,
run_metadata.json}` — nextflow trace + sweep metadata from the
**PXD071075 single-cell sweep** (Wang 2025, developing human brain
proteomics, Orbitrap Eclipse). Five sweep points at DIA-NN 2.5.0 with
`executor.queueSize` varied: 10 / 50 / 100 / 200 / 300. Each sweep
point's `run_metadata.json` confirms `queue_size == sweep_cores` —
each Nextflow task occupies one SLURM job in this configuration.

| Sweep point | queueSize | Source on PRIDE |
|---|---|---|
| q010 | 10 | `quantms-collections/quantmsdiann-benchmarks/PXD071075_cluster_sizes/v2_5_0_sweep_010cores/` |
| q050 | 50 | `… /v2_5_0_sweep_050cores/` |
| q100 | 100 | `… /v2_5_0_sweep_100cores/` |
| q200 | 200 | `… /v2_5_0_sweep_200cores/` |
| q300 | 300 | `… /v2_5_0_sweep_300cores/` |

Consumed by: **F2d** (queueSize scaling).

## Counting rule

Every identification count in this repo follows a single rule
([methods.md §1](methods.md)): exactly one admissible q-value filter per
quantity and nothing else. Per run, protein groups at `PG.Q.Value <= 0.01`
and precursors at `Q.Value <= 0.01`; global (dataset) totals at
`Lib.PG.Q.Value <= 0.01` and `Lib.Q.Value <= 0.01`. There is **no**
contaminant/target filter and **no** positive-quantity filter (zeros are
counted); decoys are dropped. Counts come from the DIA-NN report, never the
`*_matrix.tsv` files. The canonical primitive is `count_report` in
[`scripts/rebuild.py`](scripts/rebuild.py); a §4 guard flags datasets where
>1.2% of protein groups are multi-accession.

