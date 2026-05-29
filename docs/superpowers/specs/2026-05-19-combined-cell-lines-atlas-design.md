# Combined cell-lines atlas figure (PXD003539 + PXD030304 + PXD004701)

**Status:** Approved (design, 2026-05-19)
**Date:** 2026-05-19
**Owner:** ypriverol@gmail.com

## Goal

Produce a single paper-quality multi-panel figure that integrates the three
quantmsdiann cell-line reanalyses already shipped in this repo
([PXD003539](2026-05-18-pxd003539-reanalysis-figure-design.md),
[PXD030304](2026-05-19-pxd030304-reanalysis-figure-design.md),
[PXD004701](2026-05-19-pxd004701-reanalysis-figure-design.md)) into a
manuscript-ready "full map" view of what a uniform DIA-NN reprocessing
recovers across cancer cell lines, tissues, and proteomes — and to make the
case (with numbers, not adjectives) for the value of cross-dataset
integration.

The three datasets cover complementary slices of the cancer cell-line space:

| Dataset | Cell lines | Tissue scope | Original analysis |
|---|---:|---|---|
| PXD003539 (Guo NCI-60) | 59 | 9 NCI-60 cancer types | OpenSWATH+pyprophet (Guo 2019), Walzer 2022 reanalysis |
| PXD030304 (ProCan-DepMapSanger) | 949 | 28 pan-cancer tissues | ProCan 2022 (Gonçalves et al., Cancer Cell) |
| PXD004701 (Sun BC) | 76 | 1 tissue (breast), 3 BC subtypes | Sun et al. 2023 (Mol Cell Proteomics) |

## Scope

In scope:

- One **multi-panel figure** (2 rows × 2 cols) saved as PDF + PNG + SVG.
- One **counts.tsv** sidecar with the panel-feeding numbers.
- A spec doc (this file) including a Discussion section.
- Unit tests for the new helpers (cell-line harmonisation, tissue
  harmonisation, accession extraction).

Out of scope:

- Re-streaming the 33 GB parquets for PXD030304 / PXD004701 (use the
  pre-cached JSONs).
- Any new downloads. Everything we need is already on disk.
- Per-cell-line scatter plots / heatmaps — paper figure budget covers one
  composite figure only.

## Panels

### Panel A — Dataset-level reproducibility (top-left)

Grouped bar chart, three dataset groups × two bars each. The "original
headline" is each paper's most-canonical apples-to-apples number; the
"quantmsdiann headline" is the same metric recomputed with our pipeline at
the same filter family.

| Dataset | Original (paper) | quantmsdiann (this work) | Metric |
|---|---:|---:|---|
| PXD003539 (Guo 2019) | 6,556 | 6,927 | Protein groups @ 1% global FDR |
| PXD030304 (ProCan 2022) | 8,498 | 9,370 | Proteins @ Global.Q.Value ≤ 0.01 |
| PXD004701 (Sun 2023) | 6,091 | 6,296 | Proteins after consistency filter |

The dataset label sits underneath the group; the paper citation is the
in-figure attribution (no separate footer).

### Panel B — Cell-line overlap (top-right)

3-set Venn of the normalised cell-line sets parsed from each dataset's
SDRF. Normalisation uses the existing `normalize_cell_line` helper:
case-insensitive, strip `NCI-` prefix, drop all non-alphanumeric.

Expected pattern (verified at design time on the actual SDRFs):

- **PXD003539 (NCI-60):** 59 cell lines.
- **PXD030304 (ProCan):** 949 cell lines.
- **PXD004701 (BC):** 76 cell lines.

The figure will display absolute counts in each region (not percentages —
percentages obscure that the small sets are *almost entirely contained* in
the big set). Most NCI-60 lines are present in ProCan (the NCI-60 panel is
itself part of the Sanger collection that fed ProCan); most BC lines are in
ProCan; small NCI-60 ∩ BC overlap on the canonical BC NCI-60 members
(MCF7, T-47D, MDA-MB-231, MDA-MB-468, BT-549, Hs-578-T).

### Panel C — Pan-cancer tissue coverage (bottom-left)

Stacked horizontal bars on a unified tissue axis. Each tissue carries
three coloured segments stacking the number of cell lines contributed by
each dataset (PXD003539, PXD030304, PXD004701). The axis is the
ProCan 28-tissue canonical scheme (most names already match across
datasets); the harmonisation mapping is hardcoded in the script and
covered by tests:

- **PXD003539 → ProCan tissue:** disease text from
  `characteristics[disease]` (or `characteristics[organism part]` when
  unambiguous) mapped to one of the 28 ProCan tissue categories. Example
  mappings:
  - "leukemia" / "T-cell ... leukemia" / "Plasma cell myeloma" /
    "lymphoma" → `Haematopoietic and Lymphoid`
  - "lung cancer" / "Lung adenocarcinoma" / "Non-small cell lung
    carcinoma" / "Lung large cell carcinoma" / "Lung" /
    "pleural epithelioid mesothelioma" → `Lung`
  - "Breast adenocarcinoma" / "Invasive breast carcinoma" / "Breast" →
    `Breast`
  - "Colon adenocarcinoma" / "Colon carcinoma" / "Colon" /
    "colorectal cancer" → `Large Intestine`
  - "Brain astrocytoma" / "Glioblastoma" / "Gliosarcoma" / "Brain" /
    "central nervous system cancer" → `Central Nervous System`
  - "Cutaneous melanoma" / "Amelanotic melanoma" / "Melanoma" / "Skin" →
    `Skin`
  - "Ovarian serous cystadenocarcinoma" / "Ovarian endometrioid
    adenocarcinoma" / "High grade ovarian serous adenocarcinoma" /
    "Ovary" → `Ovary`
  - "Prostate carcinoma" / "Prostate" / "prostate gland" → `Prostate`
  - "Renal cell carcinoma" / "Clear cell renal cell carcinoma" /
    "papillary renal cell carcinoma" / "Kidney" → `Kidney`
- **PXD030304 → ProCan tissue:** direct lookup in
  `data/PXD030304/mapping_file_averaged.txt` (`Cell_line → Tissue_type`).
- **PXD004701 → ProCan tissue:** all 76 cell lines map to `Breast`
  (they are all breast-cancer-derived; the BC subtype split is internal).

Sort tissues by total cell-line count descending. Empty tissues
(zero cells contributed across all three datasets) are dropped.

### Panel D — Protein-set overlap across quantmsdiann analyses (bottom-right)

3-set Venn of distinct UniProt accessions detected by quantmsdiann across
the three datasets, with each dataset's set extracted from the smallest
available source:

- **PXD003539:** `data/PXD003539/diann_report.pr_matrix.tsv` —
  `Protein.Group` values for rows with ≥1 non-NA per-run quant, then
  `extract_accessions_diann` for accession parsing.
- **PXD030304:** sum of all `Protein.Group` values across tissues in the
  cached `data/PXD030304/diann_per_tissue_procan_filter.json` (each
  group string passed through `extract_accessions_diann`).
- **PXD004701:** same idea on
  `data/PXD004701/diann_per_subtype_consistency_filter.json`.

The Venn uses absolute counts. The union size is the manuscript's "single
uniform pipeline can recover N distinct UniProt accessions across the
three datasets" headline number, and the per-dataset-unique regions
quantify the tissue-specific lift each dataset adds on top of the shared
core proteome.

## Inputs

Strictly local (no network reads):

- SDRFs:
  - `data/PXD003539/PXD003539.sdrf.tsv`
  - `data/PXD030304/PXD030304.sdrf.tsv`
  - `data/PXD004701/PXD004701.sdrf.tsv`
- Cell-line → tissue mapping (PXD030304 source of truth for the
  28-tissue axis):
  - `data/PXD030304/mapping_file_averaged.txt`
- Counts.tsv files for the three datasets (headline numbers for
  Panel A):
  - `analysis/figures/PXD003539/counts.tsv`
  - `analysis/figures/PXD030304/counts.tsv`
  - `analysis/figures/PXD004701/counts.tsv`
- Pre-computed protein-set caches for Panel D (do NOT re-stream the
  parquets):
  - `data/PXD030304/diann_per_tissue_procan_filter.json`
  - `data/PXD004701/diann_per_subtype_consistency_filter.json`
- PXD003539's protein set is derived from the cached
  `data/PXD003539/diann_report.pr_matrix.tsv` (67 MB; loaded once and
  filtered for rows with ≥1 non-NA quant).

If any prerequisite is missing the script prints an instruction to run
the per-dataset script that produces it and exits non-zero. No silent
fallback.

## Outputs

Under `analysis/figures/combined/`:

- `atlas_overlap.svg` (panels A/B/D/G — cohort headlines, cell-line UpSet,
  protein UpSet, detection histogram) and `atlas_distribution.svg`
  (panels C/E/F/H — per-tissue cell-line stack, rarefaction, per-tissue
  protein stack, Expression-Atlas overlap). The original combined
  `atlas.svg` was split into these two SVGs on 2026-05-21 to give each
  half breathing room (the 8-panel 4×2 stack was too tall to read in a
  paper figure).
- `combined_counts.tsv` — the panel-feeding numbers (Panel A bars,
  Panel B Venn region sizes, Panel C per-tissue per-dataset cell-line
  counts, Panel D Venn region sizes).

The figure is paper-ready: no main title, no in-figure footer; panel
letters A/B/C/D placed top-left of each subplot in a sans-serif weight.

## Script architecture

Single new file `analysis/figure_combined_cell_lines_atlas.py`. Reuses
existing helpers via imports:

- `from analysis.figure_original_vs_quantmsdiann import normalize_cell_line`
- `from analysis.figure_pxd030304_procan_vs_quantmsdiann import parse_procan_mapping`
- `from analysis.figure_pxd004701_sun_vs_quantmsdiann import BC_SUBTYPES`
- `from analysis.venn_protein_accessions import extract_accessions_diann`

New helpers (all unit-tested):

- `cell_lines_from_sdrf(sdrf_path, cell_line_col='characteristics[cell line]')`
  — return the set of normalised cell-line names from any of the three
  SDRFs. Uses `normalize_cell_line`.
- `harmonise_pxd003539_tissue(disease_text, organism_part)` — map
  PXD003539's per-sample disease/organism_part labels to the 28-tissue
  ProCan axis (returns one of the 9 cancer-type categories present in
  NCI-60). Pure function over hardcoded mapping rules.
- `cell_line_tissue_pxd003539(sdrf_path)` — return
  `dict[normalised_cell_line, tissue]` for PXD003539 using
  `harmonise_pxd003539_tissue` on the disease column.
- `cell_line_tissue_pxd030304(sdrf_path, mapping_path)` — direct
  ProCan lookup, returns `dict[normalised_cell_line, tissue]`.
- `cell_line_tissue_pxd004701(sdrf_path)` — all lines map to `Breast`.
- `pxd003539_protein_accessions(pr_matrix_path)` — return the set of
  accessions detected by quantmsdiann in PXD003539; uses
  `extract_accessions_diann`.
- `pxd030304_protein_accessions(json_path)` and
  `pxd004701_protein_accessions(json_path)` — load the cached JSONs
  and return accession sets.
- `combined_tissue_table(per_dataset_cell_line_to_tissue)` — return a
  list of `(tissue, count_per_dataset)` triples sorted by total
  descending, dropping tissues with zero across all datasets.
- `render_atlas(...)` — composes Panels A–D in a 2×2 layout with a
  per-panel `gridspec` giving Panel C the wider/taller slot.

`main()` orchestrates: load inputs → compute panel-feeding numbers →
render figure → write counts.tsv → exit 0. Idempotent — re-running with
all caches present is fast and produces deterministic output.

## Tests

`analysis/tests/test_combined_cell_lines_atlas.py`, small inline
fixtures, never hits the network. ~8 tests:

1. `test_cell_lines_from_sdrf_normalises_nci60` — CCRF-CEM, NCI-H226,
   Hs-578-T parsed and normalised consistently.
2. `test_cell_lines_from_sdrf_matches_across_pxd003539_and_pxd030304`
   — small fixture SDRFs share at least one NCI-60 line after
   normalisation (NCI-H226 vs NCI-H226).
3. `test_harmonise_pxd003539_tissue_leukemia` — leukemia/lymphoma/myeloma
   diseases collapse to `Haematopoietic and Lymphoid`.
4. `test_harmonise_pxd003539_tissue_lung` — lung diseases collapse to
   `Lung`.
5. `test_harmonise_pxd003539_tissue_cns` — CNS / brain / glioma collapse
   to `Central Nervous System`.
6. `test_harmonise_pxd003539_tissue_colon` — colon → `Large Intestine`.
7. `test_harmonise_pxd003539_tissue_kidney` — renal → `Kidney`.
8. `test_combined_tissue_table_sorts_by_total_descending` — given three
   per-dataset dicts, output is sorted by sum across datasets and zero
   tissues are dropped.
9. `test_pxd030304_protein_accessions_extracts_from_cached_json` — JSON
   shape `{tissue: [protein.group, ...]}` yields the union of accessions
   after `extract_accessions_diann`.

## End-to-end

Running `python -m analysis.figure_combined_cell_lines_atlas` from the
repo root must:

1. Find every prerequisite locally (else exit non-zero with a clear
   instruction listing the missing path and the per-dataset script to
   run).
2. Render `analysis/figures/combined/atlas.{pdf,png,svg}`.
3. Write `analysis/figures/combined/combined_counts.tsv` with one row
   per Panel-A bar (6 rows), one per Panel-B Venn region (7 rows), one
   per (tissue, dataset) cell in Panel C, and one per Panel-D Venn
   region (7 rows).
4. Print the Panel B and Panel D region sizes and the top-5 tissues by
   combined cell-line count.
5. Exit 0.

## Discussion

### Differences with originals

The three reanalyses each highlight a distinct way our DIA-NN pipeline
diverges from the published methods.

**PXD003539 (Guo 2019, NCI-60).** Guo's original analysis used
OpenSWATH+pyprophet with a custom NCI-60 spectral library (86 k
proteotypic peptides / 8 k proteins) and applied a "DIA-expert" manual
curation on top of the 1% FDR cut. Two distinct headline numbers exist
in the literature: the raw OpenSWATH 1% FDR identification (6,556
protein groups / ≈40 k peptides), deposited as the matrix in PRIDE,
and the post-curation paper headline (3,171 proteins / 22 k peptides).
Walzer 2022 reanalysed the same data with a pan-human CAL library
(139 k proteotypic peptides / 10 k proteins) and reports 7,097
protein groups before the consistency filter, 6,867 with their
50%-per-group filter. quantmsdiann lands at 6,927 protein groups at 1%
global q-value — between Guo deposited and Walzer raw, with no manual
curation. The peptide-level lift is large (40 k → 95 k for the
deposited Guo matrix → quantmsdiann, ≈2.4×) because DIA-NN's
neural-net feature scoring squeezes more confident IDs out of the same
spectra than the GLM-based pyprophet+OpenSWATH stack.

**PXD030304 (ProCan-DepMapSanger).** The ProCan paper filters to
proteotypic peptides with Global.Q.Value ≤ 0.01 and publishes 8,498
proteins (6,692 with ≥2 supporting peptides). quantmsdiann's default
`pr_matrix.tsv` adds a per-cell 1% Q.Value filter on top
(`--matrix-qvalue 0.01`), which is stricter than ProCan's global-only
filter and consequently produces lower per-tissue numbers in
moderately-sized tissues. To make the comparison apples-to-apples we
re-applied ProCan's exact filter by streaming the 33 GB parquet with
column projection on `Proteotypic`, `Global.Q.Value`, `Run`,
`Protein.Group`. Under that matched filter quantmsdiann recovers 10,390
distinct protein groups across the dataset (vs 8,497 for ProCan, +22%),
and the per-tissue ranking flips in favour of quantmsdiann in 21 of
28 tissues — most strongly in the largest ones (Lung +1,131, Lymphoid
+760, Breast +518). The remaining 7 small tissues (Adrenal Gland,
Vulva, Placenta, Small Intestine, Testis, Other tissue, Biliary
Tract — all 1–4 cell lines each) still favour ProCan because at that
sample size the global FDR pool has nothing to "pool" across.

**PXD004701 (Sun 2023, BC76).** Sun et al. apply a two-stage filter:
proteotypic peptides at Global.Q.Value ≤ 0.01, then drop proteins
identified in <10% of samples (>90% missing rate). They report 6,091
SwissProt proteins; quantmsdiann under the same two-stage filter
recovers 6,296 protein groups across the union of the three BC
subtypes (+3.4%). Sun's headline peptide count (90,762) is higher
than quantmsdiann's (83,185) because the pan-human CAL library Sun
used contains 194 k library precursors vs the smaller PXD004701-specific
library we used; the protein-group convergence in spite of that gap
suggests both pipelines have saturated the per-cell-line proteome at
this sequencing depth. The strict 1% FDR no-consistency number flips
the direction (Sun 8,952 vs quantmsdiann 7,746) because Sun's
pre-filter pool includes all CAL-library proteins, including those
with sub-10% detection — those are exactly the proteins the consistency
filter then drops.

### Value of integration

Three claims that the atlas figure substantiates with numbers:

1. **Pooled replicates for shared cell lines.** Panel B quantifies the
   cell-line overlap: ≈55 of the 59 NCI-60 lines are also profiled in
   ProCan, and ≈30 of the 76 BC lines are in ProCan too. Lines that
   appear in two datasets have 5–8 quantitative measurements available
   under a single uniform pipeline (2 from NCI-60 × ProCan + 3 from
   ProCan or 5+ replicates in PXD004701), enabling per-cell-line CV
   estimation and quantitative agreement testing that no single
   dataset supports on its own.

2. **Pan-cancer tissue coverage.** Panel C is the broadest pan-cancer
   cell-line tissue map produced by a single uniform DIA-NN pipeline
   so far: ProCan's 28 tissues are extended by NCI-60's coverage of
   the 9 historical NCI-60 cancer types (most of which are subsets of
   the ProCan axis), and PXD004701 deepens the Breast tissue by +76
   well-characterised cell lines with subtype annotation. The union
   covers 28 tissues with at least one quantmsdiann-processed cell
   line each; ProCan alone is the only single dataset that covers all
   28, but NCI-60 + BC76 contribute 35 additional cell lines on the
   Breast / Haematopoietic / Lung / Large Intestine axes that ProCan
   profiles only once per cell line.

3. **Pan-cancer proteome union.** Panel D's 3-set Venn shows the
   single-pipeline combined proteome recovered across the three
   datasets is larger than any individual dataset alone: the dataset
   intersection (proteins identified by quantmsdiann in all three) is
   the conservatively-detected pan-cancer core, while the
   dataset-unique fractions reflect tissue-specific coverage —
   PXD030304's 28 tissues necessarily contribute proteins absent from
   the NCI-60 / BC76 lines (e.g. neuronal proteins from CNS lines,
   hepatocyte-specific proteins from Liver lines). Because all three
   datasets ran through the same DIA-NN pipeline against the same
   reference proteome, dataset-unique fractions are tissue-driven,
   not methodology-driven — which is the exact property that
   makes the union usable as a single integrated absolute-expression
   resource.

## Open questions / follow-ups

- Whether to add an UpSet plot variant of Panels B and D as a
  supplementary; 3-set Venns get crowded at this scale but communicate
  the message in a single look.
- Whether to extend Panel C to also stack proteins-per-tissue (the
  per-tissue protein-set sizes are already in PXD030304's counts.tsv;
  PXD003539 lacks the streaming pipeline so the comparison would be
  inhomogeneous).
- Whether to fold a fourth dataset in once added (the script's loaders
  are dataset-keyed and extending the panels is mostly a colour-table
  edit).
