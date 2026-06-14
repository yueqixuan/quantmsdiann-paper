# quantmsdiann reanalysis manuscript

Source repository for the **quantmsdiann methods paper**. quantmsdiann is
the Nextflow DIA pipeline that succeeds [`nf-core/quantms`](https://nf-co.re/quantms)
for DIA-only analyses and tracks DIA-NN releases (1.8.1 → 2.5.0).

For the full manuscript scope, claims, and experiment plan, see
[docs/brainstorming.md](docs/brainstorming.md). For the per-figure design
docs see [docs/superpowers/specs/](docs/superpowers/specs/). The LaTeX
first draft lives in [`paper/`](paper/) (MCP / `elsarticle` format).

## Manuscript figure → script mapping

| Figure | Claim | Script | SVG output(s) |
|---|---|---|---|
| **F2a — ProteoBench precursors by version** | quantmsdiann matches DIA-NN GUI on community benchmarks (1.8.1 / 2.5.1 / 2.5.1-enterprise) | [`figure_quantmsdiann_benchmarks_vs_proteobench.py`](analysis/figure_quantmsdiann_benchmarks_vs_proteobench.py) | `quantmsdiann_benchmarks/main_benchmarks_precursors.svg` (+ audit `data/counts.tsv`, `median_nr_prec_by_version.tsv`) |
| **F2b — ID vs accuracy scatter** | quantmsdiann matches accuracy too, not just IDs | [`figure_id_vs_epsilon.py`](analysis/figure_id_vs_epsilon.py) (consumes [`proteobench_metrics.py`](analysis/proteobench_metrics.py)) | `quantmsdiann_benchmarks/main_accuracy.svg` |
| **F1c — Per-species log2 fold-change (supp)** | per-species accuracy detail | same as F1b | `quantmsdiann_benchmarks/supplementary/supp_per_species_log2_min3.svg` |
| **F2a — Wallclock vs input size** | sub-linear scaling across instruments | [`figure_performance_runtime.py`](analysis/figure_performance_runtime.py) | `performance/parallelism_vs_wallclock.svg` |
| **F2b — Per-step runtime distribution** | DIA-NN steps dominate; glue is cheap | [`figure_performance_trace.py`](analysis/figure_performance_trace.py) | `performance/runtime_per_step.svg` |
| **F2c — Memory + CPU per task** | the cluster envelope is predictable | [`figure_performance_trace.py`](analysis/figure_performance_trace.py) | `performance/resources_per_step.svg` |
| **F2d — Cluster-node scaling** | wallclock vs cluster nodes (each Nextflow task occupies one node = `executor.queueSize`) on a real single-cell cohort (PXD071075, Wang 2025) | [`figure_queue_size_sweep.py`](analysis/figure_queue_size_sweep.py) | `performance/queue_size_sweep.svg` — 5 sweep points at 10/50/100/200/300 nodes, rendered from the PXD071075 sweep on PRIDE. Wallclock drops monotonically 37.7 h → 2.4 h with a knee at 200 nodes |
| **F3 — Cell-line atlas** | proteome coverage across NCI-60 / Sun / ProCan / Tognetti / Wang | [`figure_combined_cell_lines_atlas.py`](analysis/figure_combined_cell_lines_atlas.py) | `combined/atlas_overlap.svg` (panels A/D/G — cohort scope + protein UpSet + detection histogram), `combined/atlas_distribution.svg` (panels C/F/E — per-tissue distribution + breadth-vs-depth scatter). Panel B (cell-line UpSet) was removed because PXD030304's 947 lines dominated the inter-cohort intersections; Panel H (Expression Atlas overlap) was removed because it duplicates `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg` |
| **F3.PXDxxx — Per-cohort reanalysis** | per-dataset ID lift vs original publication | [`figure_original_vs_quantmsdiann.py`](analysis/figure_original_vs_quantmsdiann.py) (PXD003539), [`figure_pxd004701_sun_vs_quantmsdiann.py`](analysis/figure_pxd004701_sun_vs_quantmsdiann.py), [`figure_pxd030304_procan_vs_quantmsdiann.py`](analysis/figure_pxd030304_procan_vs_quantmsdiann.py) | `PXD003539/*.svg`, `PXD004701/*.svg`, `PXD030304/*.svg` |
| **F3.PXD017199 — Tognetti 2021** | atlas-only — no per-cohort figure (paper is mass-cytometry; SWATH-MS is supplementary) | covered by F3 atlas via [`figure_combined_cell_lines_atlas.py`](analysis/figure_combined_cell_lines_atlas.py) | — (atlas only) |
| **F3.PXD041421 — Wang 2023** | atlas-only — no per-cohort figure (batch-effect testbed; 2 cell lines, 48 runs, timsTOF Pro / diaPASEF) | covered by F3 atlas via [`figure_combined_cell_lines_atlas.py`](analysis/figure_combined_cell_lines_atlas.py) | — (atlas only) |

Shared utility module (not a figure): [`analysis/venn_protein_accessions.py`](analysis/venn_protein_accessions.py)
provides accession extraction + Venn rendering used by the per-cohort
reanalysis figures and the combined atlas.

**Contaminant / entrapment / decoy filter.** Every count in this repo
is target-only under a conservative filter (drop the row if any
semicolon-separated accession token has a `CONTAM_` / `Cont_` /
`ENTRAP_` / `DECOY_` / `decoy_` prefix). Canonical helper:
[`analysis/contaminant_filter.py`](analysis/contaminant_filter.py).
Per-figure audit TSVs carry paired unfiltered + target-only counts.
Design: [docs/superpowers/specs/2026-05-21-contaminant-filter-and-pxd041421-design.md](docs/superpowers/specs/2026-05-21-contaminant-filter-and-pxd041421-design.md).

## Datasets

See [MANIFEST.md](MANIFEST.md) for the full table of every PXD / module
under `data/` with instrument, file count, source FTP URL, and the
figures that consume it.

## Environment

The recommended route is the conda environment — it pins `numpy<2` (load-bearing:
numpy 2 crashes the atlas UpSet panel at render) and bundles `rsvg-convert` for the
SVG→PDF step, so a single environment reproduces **every** figure, the pan-cohort
atlas (Fig. 4) included:

```bash
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate quantmsdiann
```

A pip/venv path is also available (you must additionally install `librsvg`/`rsvg-convert`
yourself for the manuscript build):

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r analysis/requirements.txt   # numpy<2 is pinned here too
```

## Running a figure

```bash
python -m analysis.figure_quantmsdiann_benchmarks_vs_proteobench
```

Outputs land in `analysis/figures/<group>/`. Downloaded inputs are cached
under `data/<PXD-id>/` and are git-ignored (re-downloaded on demand; see
[MANIFEST.md](MANIFEST.md) for sources).

## Tests

```bash
pytest analysis/tests/
```

Tests are inline-fixture-only — no network — and cover parser, log-line,
matrix-row, and parameter-signature logic. Current count: 170 tests.

## Layout

```
analysis/
  figure_*.py                 # one figure family per script (see table)
  venn_protein_accessions.py  # shared accession + Venn utility
  tests/                      # pytest test suite
  figures/<group>/            # rendered SVGs (paper-ready, committed)
  figures/<group>/data/       # per-analysis derived data (TSVs; git-ignored, regenerated by the scripts)
  figures/<group>/supplementary/  # legacy / lower-threshold variants
data/                         # cached inputs (git-ignored, downloaded on demand)
environment.yml               # conda environment (numpy<2 + rsvg-convert) reproducing all figures
docs/
  brainstorming.md            # manuscript scope, claims, and experiment plan
  superpowers/
    specs/                    # per-figure design documents
    plans/                    # implementation plans for in-flight work
```

**Convention:** every figure script writes its rendered figures as SVG
to `analysis/figures/<group>/` and its derived TSV / JSON data tables to
`analysis/figures/<group>/data/`. PDF and PNG are not produced — open
the SVG directly in any browser or vector editor.
