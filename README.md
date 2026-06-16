# quantmsdiann reanalysis manuscript

Source repository for the **quantmsdiann methods paper**. quantmsdiann is the
SDRF-driven Nextflow/DIA-NN pipeline (<https://github.com/bigbio/quantmsdiann>);
this repo holds the LaTeX manuscript and the analysis scripts that generate its
figures.

- **Manuscript** (MCP / `elsarticle`): [`paper/`](paper/) — build with `cd paper && make pdf`.
- **Datasets**: [MANIFEST.md](MANIFEST.md) lists every PXD / module with instrument, file count, and source.
- **Design docs**: [docs/superpowers/specs/](docs/superpowers/specs/).

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
