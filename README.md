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

| Dataset (paper) | Accession | Role in paper | Our reanalysis (under `$QB/`) | Original deposit |
|---|---|---|---|---|
| HeLa Astral single-cell | PXD046357 | Fig 3 (single-cell) | `single-cell/PXD046357/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD046357) |
| HeLa One-Tip single-cell | PXD044991 | Support only³ | `single-cell/PXD044991_one-tip/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD044991) |
| Oocyte plexDIA | MSV000093870 | Support only³ | `single-cell/MSV000093870/v*/` | [MassIVE](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000093870); matrix [SlavovLab/single_cell_oocyte](https://github.com/SlavovLab/single_cell_oocyte) |
| NCI-60 (Guo 2019) | PXD003539 | Fig 4 (reanalysis + atlas) | `quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/`¹ | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD003539) |
| ProCan-DepMapSanger | PXD030304 | Fig 4 (reanalysis + atlas) | `cell-lines/PXD030304/v2_5_1/`² | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD030304) |
| Sun breast PCT-SWATH | PXD004701 | Fig 4 (reanalysis + atlas) | (cell-line panel) | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD004701) |
| Tognetti breast | PXD017199 | Fig 4 (atlas only) | (cell-line panel) | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD017199) |
| MultiPro batch testbed | PXD041421 | Fig 4 (atlas only) | (cell-line panel) | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD041421) |
| Spatial DVP (MYCN) | PXD064049 | Fig 4 (reanalysis) | `spatial/PXD064049/v2_5_1_enterprise/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD064049) |
| Phospho (NK Fe-NTA) | PXD049692 | Support only³ | `phospho/PXD049692/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD049692) |
| Phospho-enriched | PXD034128 | Support only³ | `phospho/PXD034128-{biological-study,highspeed-DIA}/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD034128) |
| Galectin-1 phospho | PXD034623 | Support only³ | `phospho/PXD034623/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD034623) |
| Scaling sweep | PXD071075 | Fig 1 (scaling) | `PXD071075_cluster_sizes/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD071075) |
| ProteoBench Module 7 (Astral) | — | Fig 2 (equivalence) | `proteobench/ProteoBench_Module_7/v*/` | [ProteoBench](https://proteobench.cubimed.rub.de/) |
| ProteoBench Module 9 (SC Astral) | PXD049412 | Fig 2 (equivalence); Fig 3 (A549/H460) | `proteobench/PXD049412/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD049412) |
| ProteoBench Module 5 (diaPASEF) | PXD062685 | Fig 2 (equivalence) | `proteobench/PXD062685/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD062685) |
| ProteoBench Module 10 (ZenoTOF) | PXD070049 | Fig 2 (equivalence) | `proteobench/PXD070049/v*/` | [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD070049) |

³ **Support only**: processed end-to-end to demonstrate workflow support, but **not**
benchmarked for identifications (DIA-NN counts are not comparable across versions for
these acquisition modes). They appear only in the runtime figures (Fig 2b / Supp Fig S3),
which show every modality finishes in minutes to hours.

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

All analysis logic lives in one self-contained script, `scripts/rebuild.py`,
exposed as named stages:

```bash
python -m scripts.rebuild --list                      # list every stage + what it produces
python -m scripts.rebuild --all                       # data prep -> all figures -> PDFs
python -m scripts.rebuild --only single_cell_combined # rebuild one figure
python -m scripts.rebuild --only paper_numbers        # re-aggregate every manuscript number
pytest tests/                                         # schema/regression tests
```

Each stage writes SVGs to `analysis/figures/<group>/` and derived tables to
`analysis/figures/<group>/data/` (or `data/`). Inputs are cached under `data/`
(git-ignored, re-downloaded on demand; see [MANIFEST.md](MANIFEST.md)). The PDF
build runs only when every prior stage succeeds.

All reported identification counts follow a single rule ([methods.md §1](methods.md)):
exactly one admissible q-value filter per quantity and nothing else. **No**
contaminant/target filter, no positive-quantity filter (zeros are counted), decoys
dropped. Per-run: protein groups at `PG.Q.Value <= 0.01`, precursors at
`Q.Value <= 0.01`. Global (dataset totals): protein groups at `Lib.PG.Q.Value <= 0.01`,
precursors at `Lib.Q.Value <= 0.01`. The counting primitive is `count_report` in
[`scripts/rebuild.py`](scripts/rebuild.py). plexDIA and phosphoproteomics are
processed end-to-end to demonstrate workflow support, but are **not** benchmarked
across DIA-NN versions (their identification counts are not comparable, so they
appear only in the runtime figures).
