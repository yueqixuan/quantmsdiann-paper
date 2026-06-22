# PXD049412 single-cell SDRF

SDRF-Proteomics annotation for the **single-cell** portion of PRIDE dataset
**PXD049412** (Brunner et al., *Challenging the Astral mass analyzer – up to 5300
proteins per single-cell …*, single-cell DIA on Orbitrap Astral). This is the
single-cell subset requested for reanalysis with the quantmsdiann / DIA-NN
workflow, distinct from the ProteoBench HYE dilution subset already annotated
elsewhere in the repo.

## File

- `PXD049412-singlecell.sdrf.tsv` — one full SDRF, 152 data rows.

A **single** SDRF was chosen (not one per cell line). It carries an explicit
`characteristics[run class]` column so downstream counting can keep only the
true single cells and drop carriers. Blanks/washes are excluded from the file
entirely (reviewer request: "run without blanks").

## Scope and run-file inventory

Selected from the deposit (357 raw files total). Classes included:

| run class    | cell line | cells/well | rows |
|--------------|-----------|-----------:|-----:|
| single-cell  | A549      | 1          | 126  |
| single-cell  | NCI-H460  | 1          | 20   |
| carrier-20x  | A549      | 20         | 3    |
| carrier-40x  | A549      | 40         | 3    |
| **total**    |           |            | **152** |

Filter to true single cells downstream with `characteristics[run class] ==
single-cell` (146 rows). Carriers are tagged `carrier-20x` / `carrier-40x`.

### A549 single cells (126) come from four prep batches
- `run69p2` (2024-06-14, Astral): 20 — same plate/run as the H460 cells.
- `2023-11-21` cellenONE plate (`SC_A549`): 66.
- `2023-12-11` Evo96 (`Evo96_A549`): 20.
- `2023-12-11` Armadillo (`Armadillo_A549`): 20.

Batch is recorded in `comment[sample preparation batch]`. The H460 cells exist
**only** in the `run69p2` (2024-06-14) plate, alongside 20 A549 cells from the
same plate — that 20+20 subset is the cleanest matched A549/H460 comparison if a
single-batch contrast is wanted.

### Carriers
The deposit's multi-cell "library/carrier" wells for A549 are **20x and 40x**
(20 and 40 pooled cells per well), per the sample-processing protocol ("for 20x
or 40x libraries, the respective cell-number was sorted into a single well").
There is **no literal 10x carrier** in the A549/H460 portion, and **no H460
carrier** at all. The reviewer's "10x and 20x carriers" is interpreted as
"include the carrier/library runs"; the actual available A549 carriers (20x,
40x) are included and clearly labelled. **Open question for the user:** confirm
whether 20x+40x is the intended carrier set, or whether only 20x should be kept.

### Excluded (present in deposit, not in this SDRF)
- All blanks / washes (50 files): `*_blank_*`, `*_wash_*`.
- HeLa, K562 single-cell / dilution runs.
- TE / hPSC (blastocyst-lineage) single cells and their 100x carriers.
- The HYE two-proteome ProteoBench dilution runs (annotated separately).

## Annotation decisions / assumptions

- **Column structure** modeled on the HeLa single-cell SDRF for PXD046357
  (`…/quantmsdiann-benchmarks/single-cell/PXD046357/v2_5_1/sdrf/`), extended with
  `characteristics[cell line]` and `characteristics[run class]`.
- **Cell lines** (CLO): A549 = `CLO:0001601`; NCI-H460 = `CLO:0008089`.
- **Disease**: A549 = lung adenocarcinoma (`EFO:0000571`); H460 = lung large cell
  carcinoma (`MONDO:0003050`).
- **organism part** = lung; **cell type** = epithelial cell (both lines are
  lung epithelial carcinoma lines).
- **Isolation protocol** = cellenONE (cellenONE robot, per sample-processing
  protocol).
- **Modifications**: per the data-processing protocol, **carbamidomethyl (C) is
  intentionally omitted for single cells** (no alkylation step was performed).
  Listed mods: Oxidation (M, variable) and Acetyl (protein N-term, variable).
  Enzyme = Trypsin. **Re-add Carbamidomethyl as a fixed mod only if you decide to
  alkylate-equivalent in the DIA-NN search; for faithful reanalysis leave it
  off.**
- **Instrument** = Orbitrap Astral (`MS:1003378`); acquisition =
  data-independent acquisition; dissociation = HCD.
- **Mass tolerances / collision energy** = `not available`. The original analysis
  was Spectronaut with dynamic tolerances; per-run NCE is not published in the
  PRIDE metadata. DIA-NN defaults (or your workflow's configured values) should
  be used — these fields are deliberately not invented.
- **`cell identifier`** is synthesized as
  `<line>_SC_<batch>_<plate-well>` from the filename's trailing plate-well token,
  giving a unique, traceable per-cell ID (required by the single-cell template).
- `biological replicate` = 1 and `technical replicate` = 1 for every run (each
  raw file is an independent single cell / well; no technical re-injections are
  declared in the deposit).

## File paths

- `comment[file uri]` points to the PRIDE archive:
  `https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD049412/<file>.raw`

## Validation

`parse_sdrf validate-sdrf --sdrf_file PXD049412-singlecell.sdrf.tsv`
(default `ms-proteomics` template) → **"Everything seems to be fine."** No
errors, no warnings.
