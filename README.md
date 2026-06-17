# quantmsdiann reanalysis manuscript

Source repository for the **quantmsdiann methods paper**. quantmsdiann is the
SDRF-driven Nextflow/DIA-NN pipeline (<https://github.com/bigbio/quantmsdiann>);
this repo holds the LaTeX manuscript and the analysis scripts that generate its
figures.

- **Manuscript** (MCP / `elsarticle`): [`paper/`](paper/) — build with `cd paper && make pdf`.
- **Datasets**: see the table below; [MANIFEST.md](MANIFEST.md) adds instrument / file-count detail.
- **Reproducibility**: [REPRODUCIBILITY.md](REPRODUCIBILITY.md) maps every paper number → source URL → generator command.
- **Design docs**: [docs/superpowers/specs/](docs/superpowers/specs/).

## Datasets and data sources

Every figure number is computed from public data. **Our reanalyses** are deposited on the
PRIDE FTP under the benchmarks base (`$QB` below); the **original deposits** we compare against
live on PRIDE Archive / MassIVE. DIA-NN version trees are `v1_8_1`, `v2_5_1`, `v2_5_1_enterprise`,
each with `quant_tables/diann_report.{parquet,tsv}` (+ `…site_report.parquet` for phospho).

```
$QB = https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks
```

| Dataset (paper) | Accession | Our reanalysis (under `$QB/`) | Original deposit |
|---|---|---|---|
| HeLa Astral single-cell | PXD046357 | `single-cell/PXD046357/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD046357) |
| HeLa One-Tip single-cell | PXD044991 | `single-cell/PXD044991_one-tip/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD044991) |
| Oocyte plexDIA | MSV000093870 | `single-cell/MSV000093870/v*/` | [MassIVE](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000093870); matrix [SlavovLab/single_cell_oocyte](https://github.com/SlavovLab/single_cell_oocyte) |
| NCI-60 (Guo 2019) | PXD003539 | `quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/`¹ | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD003539) |
| ProCan-DepMapSanger | PXD030304 | `cell-lines/PXD030304/v2_5_1/`² | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD030304) |
| Sun breast PCT-SWATH | PXD004701 | (cell-line panel) | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD004701) |
| Tognetti breast | PXD017199 | (cell-line panel) | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD017199) |
| MultiPro batch testbed | PXD041421 | (cell-line panel) | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD041421) |
| Spatial DVP (MYCN) | PXD064049 | `spatial/PXD064049/v2_5_0/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD064049) |
| Phospho (NK Fe-NTA) | PXD049692 | `phospho/PXD049692/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD049692) |
| Phospho-enriched | PXD034128 | `phospho/PXD034128-{biological-study,highspeed-DIA}/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD034128) |
| Galectin-1 phospho | PXD034623 | `phospho/PXD034623/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD034623) |
| Scaling sweep | PXD071075 | `PXD071075_cluster_sizes/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD071075) |
| ProteoBench Module 7 (Astral) | — | `proteobench/ProteoBench_Module_7/v*/` | [ProteoBench](https://proteobench.cubimed.rub.de/) |
| ProteoBench Module 9 (SC Astral) | PXD049412 | `proteobench/PXD049412/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD049412) |
| ProteoBench Module 5 (diaPASEF) | PXD062685 | `proteobench/PXD062685/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD062685) |
| ProteoBench Module 10 (ZenoTOF) | PXD070049 | `proteobench/PXD070049/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD070049) |

¹ NCI-60 is served from the earlier `quantms-collections` deposition. ² ProCan is also mirrored
under `$QB/cell-lines/PXD030304/v2_5_1/`. Tooling: pipeline [bigbio/quantmsdiann](https://github.com/bigbio/quantmsdiann),
SDRF adapter [bigbio/sdrf-pipelines](https://github.com/bigbio/sdrf-pipelines) (`convert-diann`),
containers [bigbio/quantms-containers](https://github.com/bigbio/quantms-containers).

## Setup

```bash
conda env create -f environment.yml   # pins numpy<2, bundles rsvg-convert
conda activate quantmsdiann
```

(pip/venv alternative: `pip install -r analysis/requirements.txt`, plus a system
`rsvg-convert`/`librsvg` for the manuscript PDF build.)

## Usage

```bash
python -m analysis.figure_single_cell_combined   # one figure family per script
pytest analysis/tests/                           # test suite
```

Each script writes SVGs to `analysis/figures/<group>/` and derived tables to
`analysis/figures/<group>/data/`. Inputs are cached under `data/` (git-ignored,
re-downloaded on demand; see [MANIFEST.md](MANIFEST.md)).

All reported counts are target-only under a conservative contaminant / entrapment
/ decoy filter (drop a protein group if any accession token carries a `CONTAM_` /
`Cont_` / `ENTRAP_` / `DECOY_` prefix); canonical helper:
[`analysis/contaminant_filter.py`](analysis/contaminant_filter.py).
