# Contaminant/entrapment/decoy filter + PXD041421 atlas integration — design

**Status:** spec (not yet implemented)
**Date:** 2026-05-21
**Owner:** ypriverol@gmail.com
**Manuscript figures touched:** F1a, F3 atlas (all 8 panels via cohort additions), F3.PXDxxx per-cohort
**Brainstorm reference:** [§3 use case](../../brainstorming.md#3-the-use-case--cell-line-atlas) (atlas) + [§2.1 equivalence](../../brainstorming.md#21-equivalence--you-get-the-same-results-as-dia-nn-gui) (benchmark headline counts)

This spec bundles two related changes that share a single implementation
pass:

1. **Apply a conservative contaminant/entrapment/decoy filter** across
   every count site in the codebase. Today our `extract_accessions_diann`
   strips prefixes instead of dropping prefixed rows, so contaminant /
   entrapment entries leak into our protein catalogs as fake target
   accessions. Quantified impact: 1.9 % (PXD003539) to 8.3 % (PXD030304)
   of cell-line cohort rows; ~1 % of benchmark precursor rows.
2. **Add PXD041421 (He Wang 2023)** as a 5th atlas cohort. Cell-line
   coverage is unchanged (A549, K562 already in PXD003539 / PXD030304),
   but the dataset adds timsTOF Pro / diaPASEF as a 3rd instrument
   family and a batch-effect testbed.

Both changes touch the same set of files (atlas, per-cohort scripts,
benchmark counts). Implementing together avoids a second sweep over the
same code.

---

## Part 1 — Conservative contaminant filter

### 1.1 Goal

Ensure every reported protein-group / accession / gene count across the
manuscript reflects **target-only** identifications: no contaminants,
no entrapment proteins, no decoys.

### 1.2 The bug today

`analysis/venn_protein_accessions.py:21`:

```python
_PREFIX_RE = re.compile(r"^(?:CONTAM_|ENTRAP_|DECOY_)+")
```

Two defects:

- **Strips, doesn't drop.** `extract_accessions_diann("CONTAM_HORSE_ALB")`
  returns `{"HORSE_ALB"}`, adding a phantom human protein accession to
  our union catalog. Pure-contaminant rows should return the empty set.
- **Misses the ProteoBench prefix.** Cell-line cohorts use
  `CONTAM_` (upper-case) per the entrap-contaminants FASTA; ProteoBench
  benchmarks use `Cont_` (capitalised, single underscore). The current
  regex matches only the cell-line convention, so ~1,200 `Cont_` rows
  per benchmark dataset leak into our F1a precursor counts.

Concrete impact at the row level. Unique-PG counts come from different
sources depending on what is cached on disk for each cohort: the
PXD003539 / PXD017199 / PXD041421 numbers count distinct Protein.Group
strings on `diann_report.pr_matrix.tsv`; PXD030304 and PXD004701 count
distinct strings across the cached per-tissue / per-subtype JSONs. The
per-tissue cache aggregates across the per-tissue consistency filter
and so has a slightly different denominator than pg_matrix.tsv — the
audit TSV will spell out both per cohort during implementation.

| Dataset | Source | Unique PG | CONTAM_/Cont_ | ENTRAP_ | DECOY_ | Bad % |
|---|---|---|---|---|---|---|
| PXD003539 (Guo) | pr_matrix | 6,875 | 90 | 39 | 0 | 1.9 % |
| PXD030304 (ProCan) | per-tissue JSON | 10,390 | 161 | 706 | 0 | **8.3 %** |
| PXD030304 (ProCan) | pg_matrix | 9,251 | (incl.) | (incl.) | 0 | 3.6 % |
| PXD004701 (Sun) | per-subtype JSON | 6,296 | 82 | 68 | 0 | 2.4 % |
| PXD017199 (Tognetti) | pr_matrix | 10,526 | 164 | 97 | 0 | 2.5 % |
| PXD041421 (Wang) | pr_matrix | 9,090 | 104 | 128 | 0 | 2.6 % |
| PXD049412 / PXD062685 / PXD070049 / Module 7 | pr_matrix | 45-122k each | ~1,200 each | 0 | 0 | 1-3 % |

### 1.3 Filter policy — **conservative**

A Protein.Group string passes the filter iff **none** of its
semicolon-separated tokens carry a prefix from the set
`{CONTAM_, Cont_, ENTRAP_, DECOY_, decoy_}`.

This is stricter than necessary: a mixed group like
`CONTAM_P02768;P02768` (contaminant entry sharing peptides with real
human albumin) will be **dropped**, even though one accession is a
real target. The justification: when DIA-NN's protein-grouping has
placed a target inside a contaminant-named group, the inference is
ambiguous; conservatively excluding such rows guards against
contamination of the target catalog. This is the user-chosen policy
(2026-05-21 brainstorm).

### 1.4 Affected prefixes — case-sensitive

Recognised exactly, no case-insensitive matching:

- `CONTAM_` — cell-line FASTA `Homo-sapiens-uniprot-reviewed-entrap-contaminants-202605.fasta`
- `Cont_` — ProteoBench `ProteoBenchFASTA_*` series
- `ENTRAP_` — cell-line FASTA entrapment proteins
- `DECOY_` — DIA-NN decoy convention (rarely seen post-FDR, defensive)
- `decoy_` — lowercase variant (defensive; not currently observed in our
  data, low cost to include)

### 1.5 The canonical helper

New module `analysis/contaminant_filter.py`:

```python
_FILTER_PREFIXES = ("CONTAM_", "Cont_", "ENTRAP_", "DECOY_", "decoy_")

def is_target_protein_group(pg_string: str | None) -> bool:
    """Return True iff every semicolon-separated accession in `pg_string`
    is a target (carries no contaminant/entrapment/decoy prefix). Empty
    / whitespace-only inputs return False (defensive)."""

def target_accessions(pg_string: str | None) -> set[str]:
    """Return the set of clean accessions if `is_target_protein_group`
    is True; otherwise the empty set. Reuses the existing _clean_token
    normalisation (sp|...|... → middle field, strip isoform suffix)."""
```

`venn_protein_accessions.extract_accessions_diann` becomes a thin
wrapper over `target_accessions` for backwards compatibility. The
existing `_PREFIX_RE` is replaced by `_FILTER_PREFIXES`.

### 1.6 Call sites that move to the new policy

**Cell-line side** (all use the entrap-contaminants FASTA):

- `figure_combined_cell_lines_atlas.py`:
  `pxd003539_protein_accessions`, `pxd017199_protein_accessions`,
  `pxd017199_accessions_per_cell_line`,
  `_accessions_from_json_cache` (filter Protein.Group at read time so
  the existing per-tissue / per-subtype JSON caches don't need to be
  regenerated).
- Same file's `DATASET_HEADLINES.diann_count`: recompute per-cohort
  protein-group counts from `diann_report.pg_matrix.tsv` after applying
  the filter, instead of trusting the unfiltered `diannsummary.log`
  headline. The diannsummary line stays in the audit TSV for context.
- `figure_original_vs_quantmsdiann.py` (PXD003539): every protein-group
  count in `counts.tsv` (currently four: post-filter, no-consistency,
  quantmsdiann strict, etc.).
- `figure_pxd004701_sun_vs_quantmsdiann.py`: ditto.
- `figure_pxd030304_procan_vs_quantmsdiann.py`: ditto.

**Benchmark side** (ProteoBench FASTAs, `Cont_` prefix):

- `figure_quantmsdiann_benchmarks_vs_proteobench.py`:
  `count_pr_matrix_rows` and `count_pr_matrix_min_replicates` switch
  to "rows whose Protein.Group passes the filter". F1a parity nr_prec
  values drop by ~1 %; community PB submissions are already filtered
  by their own `contaminant_flag = Cont_` in ProteoBench's parser, so
  the gap to the community baseline closes.
- F1a parity TSV (`diann_quantmsdiann_parity_long.tsv`,
  `diann_quantmsdiann_parity_epsilon.tsv`,
  `median_nr_prec_by_version.tsv`): recomputed with the filter.

**Out of scope (no change needed):**

- F1b / F1c (proteobench API path): the `convert_to_standard_format`
  step in ProteoBench's parser sets a `contaminant` column flag via
  `contaminant_flag = Cont_` and excludes those rows from the metric
  computation. Our melt → parser invocation gets the filter for free.
- F2a / F2b / F2c / F2d: performance figures that don't count proteins.
- F3 atlas Panel H: gene-symbol overlap with E-PROT-73; the
  `pr_matrix.tsv` Genes column will now be filtered (PG passes filter →
  the row's genes contribute). Coverage % is recomputed.

### 1.7 Audit trail

Every per-figure `data/<...>.tsv` gains paired rows reporting both the
pre-filter and post-filter count so a reviewer can audit the
contamination delta. Schema unchanged.

Example new rows in
`analysis/figures/combined/data/combined_counts.tsv` (PXD030304 values
verified against `diann_report.pg_matrix.tsv` at spec-write time —
final implementation may differ slightly if the on-disk matrix is
re-fetched):

```
Panel A | PXD030304 | quantmsdiann diann_count | unfiltered (pg_matrix rows pre-filter) | 9251  | DIA-NN protein groups at 1% global q-value, incl. CONTAM_/ENTRAP_
Panel A | PXD030304 | quantmsdiann diann_count | target-only (conservative filter)       | 8921  | drops 330 rows (3.6%) whose Protein.Group has CONTAM_/ENTRAP_ token
```

The headline `count` column on the figure carries the target-only
number.

### 1.8 Tests

New `analysis/tests/test_contaminant_filter.py`:

```
test_pure_contaminant_returns_empty
test_mixed_group_dropped_under_conservative
test_pure_target_passes
test_multi_target_preserved
test_lowercase_cont_recognised        # ProteoBench prefix
test_lowercase_decoy_recognised
test_is_target_protein_group_predicate
test_empty_or_whitespace_returns_false
```

Existing
`analysis/tests/test_venn_protein_accessions.py::test_extract_accessions_diann`
needs an update: `extract_accessions_diann("CONTAM_P02768-1;P02768") == {"P02768"}`
becomes `== set()`. The old behavior is the bug; the change is
intentional.

Existing per-cohort tests that pin specific count values may need
their numbers updated. Audit pass during implementation.

---

## Part 2 — PXD041421 atlas integration

### 2.1 Dataset facts

- **PXD041421** (He Wang, 2023 — TIQUEST batch-effect testbed)
- **48 DIA runs**, 2 cell lines: A549 (lung adenocarcinoma, 24 runs) +
  K562 (CML blast-phase, 24 runs)
- **2 batches** × **3 biological replicates** × **4 fractions** per
  cell line, explicit batch design (see run-column prefixes
  `CAD…` and `N…`)
- **Instrument:** timsTOF Pro / diaPASEF (1st diaPASEF in our atlas)
- **DIA-NN 2.5.0** headline (unfiltered): 152,779 precursors, 9,124
  protein groups at 1 % global q-value
- **FASTA:** same `Homo-sapiens-uniprot-reviewed-entrap-contaminants-202605.fasta`
  as the other cell-line cohorts → filter from Part 1 applies
- **PRIDE FTP:** [`/quantms-collections/absolute-expression-2.0/cell-lines/PXD041421/`](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD041421/)

### 2.2 Scope decision: atlas-only, no per-cohort figure

Same precedent as PXD017199 in [[pxd017199-atlas-only]]. The
PXD041421 deposit is a *methodological* dataset — its purpose is
batch-effect benchmarking, not deep-proteome discovery. The associated
paper carries no DIA proteomics headline number to compare against.
PXD041421 is therefore added as a 5th atlas cohort only; no
`figure_pxd041421_*.py` script.

### 2.3 What PXD041421 adds (and doesn't)

Honest assessment for the atlas reader:

- **Cell-line coverage**: A549 and K562 are already in PXD003539
  (NCI-60) and PXD030304 (ProCan). PXD041421 contributes **zero
  unique cell lines**. Panels B, C, E will show PXD041421 as a strict
  subset of existing cohorts.
- **Protein coverage**: deep diaPASEF on 24 reps may reveal proteins
  the other cohorts miss. Panel D (protein UpSet) and Panel G
  (detection histogram) will surface this.
- **Instrument family**: only timsTOF Pro / diaPASEF in the atlas;
  diversifies the manuscript's instrument-coverage claim.
- **Batch-effect insight**: out of scope for the atlas figure (atlas
  is about union proteome). Worth a footnote in the methods, not a
  panel.

### 2.4 Atlas plumbing changes

Mirror the PXD017199 integration (see [[pxd017199-atlas-only]]). The
changes are mechanical:

- `DATASET_COLORS["PXD041421"] = "#8da0cb"` (lavender — distinct from
  the existing 4 palette entries)
- `DATASET_LABELS["PXD041421"] = "PXD041421\n(Wang 2023)"`
- `DATASET_HEADLINES["PXD041421"]`: `paper_label=""`, `paper_count=0`
  (no paper bar — same convention as PXD017199), `diann_count` from
  the **filtered** pg_matrix.tsv count.
- `PXD041421_SDRF`, `PXD041421_PR_MATRIX` path constants.
- `pxd041421_protein_accessions(pr_matrix_path)` — filtered union.
- `pxd041421_accessions_per_cell_line(pr_matrix_path, sdrf_path)` —
  for Panel E rarefaction. Since only 2 cell lines × 24 reps, the
  rarefaction "curve" is effectively 2 points; render anyway for
  completeness, with a footnote.
- `cell_line_tissue_pxd041421(sdrf_path)`: A549 → "Lung Cancer", K562 →
  "Leukemia" (re-use the existing tissue labels already present in the
  unified-axis table).
- `ds_order` extended to 5: `[PXD003539, PXD030304, PXD004701,
  PXD017199, PXD041421]`.
- `_set_region_sizes` already generic over N sets (from the
  PXD017199 work) — 2^5 - 1 = 31 region keys.
- Panel B + D UpSet plots auto-scale to 5 sets.
- Panels C, F (stacked bars) gain a 5th colour.
- Panel E rarefaction gains a 5th curve (2 points; annotate "2 cell
  lines only" in the legend).
- Panel G detection histogram becomes 5 buckets (proteins detected in
  1, 2, 3, 4, 5 datasets).
- Panel H stays NCI-60-specific (PXD003539 ↔ E-PROT-73).

### 2.5 Risks and notes

- **Visual crowding**: 5 colours in the UpSet plot and stacked bars is
  near the limit of legibility. We may need to bump `figsize` from
  `(14, 22)` to `(14, 24)` and tighten panel padding.
- **Misleading "PXD041421-only" region**: because A549 and K562 are
  already in other cohorts, the "PXD041421-only" cell-line region will
  be exactly **zero**. Panel B UpSet will reflect this honestly. Panel
  D may show some unique proteins (deep diaPASEF coverage) — this is
  the dataset's true contribution and should be highlighted in the
  manuscript caption.
- **Rarefaction with 2 cell lines**: looks anaemic next to PXD003539's
  ~60-point curve. Acceptable — annotate.

---

## Implementation order

1. Land the contaminant filter (Part 1) first. Re-render every
   affected figure; record before/after numbers in the audit TSVs.
   Tests pass.
2. Add PXD041421 (Part 2) on top of the filtered code, so
   PXD041421's accession sets and headline counts are already
   computed under the new policy.
3. Update README + MANIFEST + brainstorming to reflect both changes.

Single PR. Two commits.

## Out of scope

- Pipeline-side filtering: no change to quantmsdiann (which is
  upstream of this repo). The pipeline continues to report
  contaminants in `diann_report.pg_matrix.tsv` and `pr_matrix.tsv`;
  this repo filters at the consumer level.
- Reconciling cell-line cohort headlines against the original
  publications' numbers. The "vs paper" bars in per-cohort figures
  keep the paper's original count; only the quantmsdiann side gets
  the filter. Methodological footnote in the figure caption explains
  the asymmetry.
- A batch-effect-focused per-cohort figure for PXD041421. The
  dataset's batch-effect properties are real but the atlas figure is
  about union proteome coverage; surfacing batch effects would need
  a separate figure outside §3 use-case scope.

## Validation criteria

After implementation:

- Every per-figure `data/<...>.tsv` carries both pre-filter and
  post-filter counts for every protein-group / precursor / accession /
  gene metric.
- Atlas Panel A `diann_count` values decrease by 2-8 % per cell-line
  cohort; ProteoBench benchmark nr_prec values decrease by ~1 %.
- F1a parity epsilons recompute correctly with the new counts.
- F3 atlas renders with 5 cohorts; PXD041421's "unique cell-line"
  region in Panel B UpSet is empty as expected.
- All existing tests pass; the test count grows by ~10 (filter tests +
  PXD041421-specific tests).
- No SVG carries a count that contradicts the audit TSV.
