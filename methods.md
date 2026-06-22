# methods.md — quantmsdiann manuscript: data, scripts, data flow, and filters

Single human-readable source of truth for how every number and figure in the
manuscript is produced. Maintained alongside the code: if a filter, dataset, or
script changes, this file changes in the same commit.

Reproduce everything with one command:

```bash
python -m scripts.rebuild --all      # data prep -> all figures -> PDFs
python -m scripts.rebuild --list     # show every stage + what it produces
```

---

## 1. Filter specification (IMPERATIVE — scientific correctness)

Per the 2026-06-21 review (V. Demichev). Every reported identification count
falls into one of two classes, each with EXACTLY one admissible filter and
nothing else — no contaminant/target filter, no positive-quantity filter (zero
quantities are counted), no extra q-value gates. Decoys (`Decoy == 1`) are
always dropped (that removes the FDR null model; it is not one of the forbidden
"filters").

| Quantity | Unit | Admissible filter | Nothing else |
|---|---|---|---|
| Protein groups, per run/cell | within one run | `PG.Q.Value <= 0.01` | no target, no global, zeros counted |
| Precursors, per run/cell | within one run | `Q.Value <= q` (run-specific `--qvalue`) | no `Global.Q.Value` |
| Protein groups, global (dataset total) | union | `Lib.PG.Q.Value <= 0.01` | no target, no `Global.PG.Q.Value` |
| Precursors, global (dataset total) | union | `Lib.Q.Value <= 0.01` | no `Q.Value`/`Global.Q.Value` |

Run-specific precursor `q` = DIA-NN's per-version `--qvalue`: **0.01** for 1.8.1,
**0.05** for >= 2.5.0 (the recommended operating point, applied by the pipeline).

A dataset TOTAL (union across runs) is a **global** number, so headline "total
precursors / protein groups" use `Lib.Q.Value` / `Lib.PG.Q.Value`. Per-cell and
per-run distributions, and >=3-run replicate counts, use `Q.Value` / `PG.Q.Value`.

Canonical implementation: the `count_report` function in
[`scripts/rebuild.py`](scripts/rebuild.py) (formerly `analysis/count_report_ids.py`).
Emitted keys: `prec_min1`, `prec_min3` (Q.Value replicate), `prec_global`
(Lib.Q.Value), `prot_global` (Lib.PG.Q.Value), `prot_perrun_avg`, `prot_complete`
(PG.Q.Value), `peptides`, `prot_2pep` (global rule).

### plexDIA (MSV000093870) — channel is the per-run unit

A single cell = one channel of one run. The per-cell precursor q-value is
therefore `Channel.Q.Value` (the channel-level analog of `Q.Value`; run-level
`Q.Value` passes a precursor in all 3 channels and is not a per-cell number).
Per-cell protein groups use `PG.Q.Value`. No quantity/contaminant filter.
Implemented in [`scripts/rebuild.py`](scripts/rebuild.py) (the `plexdia_per_cell`
stage; functions `load_channel_confident` / `per_cell_counts`, formerly
`analysis/plexDIA/figure_msv000093870_oocyte_plexdia.py`).
**Open flag for Vadim:** with `Channel.Q.Value` only, median precursors/cell is
~15.9k (proteins/cell ~1.47k) for these deep oocyte cells; this may warrant an
additional run-level `Q.Value <= 0.01` gate — pending confirmation.

### Phosphoproteomics
Phosphopeptides counted at the global precursor rule (`Lib.Q.Value <= 0.01`).
Phosphosite localisation uses the DIA-NN site report `Probability >= 0.99`
(class I) — a localisation-confidence axis, not one of the forbidden q-filters.

### Original / deposited side of reanalysis comparisons
Deposited datasets provide only `*_pg_matrix.tsv` / `*_pr_matrix.tsv` (already
q-filtered count matrices, no q-value columns). The "original" bar is the matrix
row count as-is; the rule governs the quantmsdiann (reanalysis) side.

## Questions for Vadim (need confirmation)

> These are the open decisions where the rule, as applied, has a judgement call.
> Each is implemented with a documented default; flag here if a different choice
> is wanted before submission.

- **[VADIM] "Total" = global rule.** A dataset union/total is treated as a GLOBAL
  number, so headline "total precursors / protein groups" use `Lib.Q.Value` /
  `Lib.PG.Q.Value` (not the per-run `Q.Value`/`PG.Q.Value` union). Per-cell/per-run
  distributions and >=3-run replicate counts use the per-run columns. Confirm this
  split.
- **[VADIM] plexDIA per-cell precursors.** Counted with `Channel.Q.Value <= 0.01`
  only (the channel is the per-run unit). Run-level `Q.Value`-only is wrong here
  (it gives ~16-19k precursors/cell). Current result: median ~15.9k precursors/cell
  and ~1.47k proteins/cell for these deep oocyte cells. **May want an added run-level
  `Q.Value <= 0.01` gate** — confirm whether `Channel.Q.Value` alone is acceptable.
- **[VADIM] plexDIA per-cell protein median dropped below the original.** Under the
  stricter `PG.Q.Value` per-cell gate, quantmsdiann median is 1,473 vs Galatidou
  1,784 per cell, while TOTAL protein groups rose to 3,328 vs 2,122 (+57%) and
  QC-matched r=0.97. Confirm this "more total, slightly lower per-cell median"
  framing is acceptable.
- **[VADIM] Reanalysis "original" bar for DIA-NN-originated deposits (4a).** For
  HeLa Astral (PXD046357) and Spatial DVP (PXD064049) the deposited analysis was
  itself DIA-NN 1.8.1, so the grey "original" bar is ambiguous: deposited/published
  count vs OUR 1.8.1 rerun. **Default = deposited/published** (HeLa Astral original
  3,903 PG / 19,365 prec; our 1.8.1 rerun would instead be 4,134 / 22,531). Confirm.
- **[VADIM] Benchmark "precursors" metric.** The ProteoBench benchmark headline
  "total precursors" uses the global `Lib.Q.Value` count (`prec_global`, +9-22%
  v1.8.1->ent), and the supplement (Note 5 / Fig S3) was made consistent with it.
  An alternative is the per-run-union "precursors identified in >=1 run"
  (`prec_min1`, Q.Value, +29-37%) — closer to ProteoBench's own precursor metric.
  Default = global `Lib.Q.Value`. Confirm which the benchmark should report.
- **[VADIM] Phosphopeptides + peptides under the global rule.** Phosphopeptides and
  the distinct-stripped-peptide / >=2-peptide-protein counts are gated on the global
  `Lib.Q.Value` / `Lib.PG.Q.Value`. Phosphosite localisation keeps `Probability >=
  0.99` (class I) as a localisation-confidence axis, not a forbidden q-filter.
  Confirm the peptide/site rules.

---

## 2. Data inventory

All reanalysis reports are on the public PRIDE FTP under
`databases/pride/resources/proteomes/quantmsdiann-benchmarks/`:
`proteobench/`, `single-cell/`, `cell-lines/`, `spatial/`, `phospho/`,
`PXD071075_cluster_sizes/`. Each cohort/version holds
`v<ver>/quant_tables/diann_report.parquet` (>= 2.x) or `diann_report.tsv`
(1.8.1), which carry the `PG.Q.Value`, `Q.Value`, `Lib.PG.Q.Value`,
`Lib.Q.Value` columns the rule needs. Local caches live under `data/` and
`analysis/figures/<DATASET>/data/` (both gitignored).

| Accession | Class | Role | Versions |
|---|---|---|---|
| ProteoBench_Module_7 | Astral DIA | ProteoBench M7 benchmark | 1.8.1 / 2.5.1 / 2.5.1-ent |
| PXD049412 | single-cell DIA | ProteoBench M9 benchmark (HYE subset) | 1.8.1 / 2.5.1 / 2.5.1-ent |
| PXD062685 | diaPASEF | ProteoBench M5 benchmark | 1.8.1 / 2.5.1 / 2.5.1-ent |
| PXD070049 | ZenoTOF | ProteoBench M10 benchmark | 1.8.1 / 2.5.1 / 2.5.1-ent |
| PXD046357 | single-cell HeLa Astral | reanalysis | 1.8.1 / 2.5.1-ent |
| PXD049412 | single-cell A549/H460 (Astral) | reanalysis (146 single cells; 20x/40x A549 library runs excluded; blanks/HeLa-dilution/HYE not used) | 1.8.1 / 2.5.1-ent |
| MSV000093870 | oocyte plexDIA | reanalysis | 2.5.0 |
| PXD003539 | NCI-60 bulk | reanalysis | 2.5.x |
| PXD030304 | ProCan bulk | reanalysis | 2.5.x |
| PXD004701 | Sun bulk | reanalysis | 2.5.x |
| PXD064049 | spatial DVP diaPASEF | reanalysis | 2.5.0 |
| PXD049692 / PXD034128 / PXD034623 | phospho | reanalysis (supp) | 2.5.1 / -ent |

_(PXD044991 "One-Tip" has been removed from the manuscript — see the design
spec. PXD049412 is a mixed deposit: its ProteoBench M9 HYE subset feeds the
benchmark (Fig 2), and its A549/H460 single cells (146 cells + 20x/40x A549
library runs) were reanalysed for the single-cell figure (Fig 4); the 50 blanks
and the HeLa low-input dilution series in the same deposit are not used.)_

---

## 3. Script inventory + data flow

**All analysis logic lives in ONE self-contained script, `scripts/rebuild.py`.**
The former `analysis/*.py` modules were inlined into it; each is now a named
**stage** (run a subset with `--only <stage>`, see them all with
`python -m scripts.rebuild --list`). There are no separate `analysis/*.py`
files anymore -- only `analysis/figures/` outputs + `analysis/requirements.txt`
remain. The counting primitive `count_report` (the §1 rule) is a function
inside the script, shared by every counting stage.

Flow: **FTP report -> `count_report` -> figure-data TSV -> figure SVG -> (rsvg) PDF -> manuscript**.

The stages below are listed in dependency order (data prep first, then figures,
then the numbers aggregator); `--list` prints the same registry live.

| Stage (`--only`) | Group | Produces | Filter path |
|---|---|---|---|
| `report_counts` | data | `data/quantmsdiann_benchmarks/report_counts.tsv` | downloads FTP reports, runs `count_report` |
| `reanalysis_pg_counts` | data | per-cohort `diann_report_protein_counts.json` caches | Lib rule; inputs for the bulk-cohort figures |
| `single_cell_tables` | data | `data/single_cell/mv_*.tsv`, `sc_totals.tsv` | per-cell `PG.Q.Value`; totals global |
| `phospho_tables` | data | phosphopeptide / phosphosite tables | `Lib.Q.Value`; site Prob>=0.99 |
| `benchmarks` | fig | Fig 2 benchmark panels, `counts.tsv` | reads `report_counts.tsv` |
| `queue_sweep` | fig | `queue_size_sweep.tsv` (feeds `fig2_validation`) | runtime trace |
| `fig2_validation` | fig | Fig 2 validation composite | composes sweep + accuracy |
| `id_vs_epsilon` | fig | ProteoBench id-vs-epsilon panel | reads `report_counts.tsv` |
| `proteobench_accuracy` | fig | ProteoBench accuracy panels | HYE fold-changes |
| `reanalysis_improvement` | fig | reanalysis-recovery figure | original matrix vs reanalysis report |
| `single_cell_combined` | fig | single-cell figure | `single_cell_tables` outputs |
| `plexdia_per_cell` | fig | plexDIA per-cell depth (MSV000093870) | `Channel.Q.Value` / `PG.Q.Value` |
| `plexdia_vs_galatidou` | fig | plexDIA deposited vs quantmsdiann | same per-cell counter |
| `pxd003539` | fig | PXD003539 (NCI-60) panels | matrix vs report |
| `pxd004701` | fig | PXD004701 (Sun) panels | matrix vs report |
| `pxd030304` | fig | PXD030304 (ProCan) panels | streams a 2GB matrix |
| `pxd064049_spatial` | fig | spatial DVP panels | matrix-based (see spec Goal 5) |
| `atlas` | fig | pan-cohort (Fig S13) | needs numpy<2 env |
| `phospho` | fig | phosphoproteomics supp figure | `Lib.Q.Value`; site Prob>=0.99 |
| `venn` | fig | protein-accession overlap (supp) | per-cohort |
| `performance_trace`, `mdc_cluster_runtime` | fig | Fig 1 / performance traces | runtime + resource traces |
| `paper_numbers` | num | `data/paper_numbers.tsv` + `paper/generated_numbers.tex` | aggregates ALL manuscript numbers |
| `paper/Makefile` | pdf | SVG->PDF + `main`/`supplementary` PDFs | rsvg-convert |

The orchestrator [`scripts/rebuild.py`](scripts/rebuild.py) runs all of the
above in dependency order; `--list` prints this table from the live stage
registry.

---

## 4. Audit status
See the design spec `docs/superpowers/specs/2026-06-21-vadim-review-changes-design.md`
(Goal 8) for the three-pass audit protocol (general / filters / paper-text) and
`docs/superpowers/specs/numbers_old_vs_new.md` for the benchmark old-vs-new diff.
