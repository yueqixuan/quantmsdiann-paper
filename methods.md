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

# CRITICAL DIRECTIVE

The following is essential to achieve uncompromising scientific correctness and transparency. Failure to do so is unacceptable.  

The manuscript involves comparisons of metrics obtained with different softwares, software versions and workflows. The comparisons must be audited for correctness. Comparisons are only allowed if the way particular metric is calculated is considered comparable between the softwares. While strict comparativeness cannot ever be guaranteed, given possibly vastly differnet FDR control accuracy or lack thereof in different softwares, the manuscript provides such comparison numbers to illustrate to the reader what numbers can be expected with DIA-NN 2.5.1 Enterprise used as part of quantmsdiann and how they may compare to previously published works. Only metrics explicitly designated as comparable are considered comparable. Any comparisons of numbers not explicitly permitted are not allowed in the manuscript. Any such comparison occurrence must be removed.  

The manuscript must explicitly indicate the following considerations, where appropriate and if applicable to the data being described:
* DIA-NN 2.5/2.5.1 is assumed to have largely correct FDR control; prior versions or other softwares may or may not have vastly poor FDR control, with possibly FDR at 10% or higher at reported 1% q-values. The proof of correct control by recent DIA-NN is beyond the scope of this manuscript, the verification of other softwares or prior versions is technically not possible on specific data considered and is likewise out of scope.  
* Numbers of localised sites in phosphoproteomics cannot be compared between DIA-NN 2.5/2.5.1 and any other softwares or versions. Reason: false localisation rates are known for a fact to be drastically different, at same reported confidence levels.
* An immunopeptidomics workflow in DIA-NN 2.5/2.5.1 would normally entail peptidoform-level confidence (achieved with Peptidoforms or Proteoforms mode in DIA-NN), making the data only partially comparable to outputs of other tools that cannot ensure peptidoform confidence. Thus, significant advantages of DIA-NN are therefore interpetable as 'likely' (a degree of uncertainty still present, see writing imperative below) lower bound for the true advantage of DIA-NN, while lack of such significant advantage, if occurs on particular data, is uninformative (not interpretable in any meaninful way) and must be flagged for refactor by the user.  

Compliance: third-party commercial tools Spectronaut, FragPipe, MSFragger, PEAKS or Chimerys must NOT be mentioned by name. If a comparison with a public dataset analysed by one those tools is performed, the dataset must be referred to as FirstAuthor et al., without indicating the software used.  

# Writing IMPERATIVE

This repostitory is a scientific manuscript draft. The manuscript is intended as a top quality highly rigorous scientific publication. The following is imperative:
* Any claim or statement made is strictly correct and substantiated with an appropriate reference to data being presented or a prior publication, unless obvious to a typical scientist within the field.
* Notwithstanding the above, speculation or opinion by the authors are allowed, if appropriately and transparently expressed as such and if properly grounded in presented data or prior knowledge. Examples of wording required include explicit indications of uncertainty (works like 'speculate', 'likely', 'suggests'). 
* The manuscript text is highly optimised. The overall structure fits the narrative and comprises sections and paragraphs that fit cohesively together. Every sentence serves specific purpose of communicating or explaining something to the reader, driven by the narrative: no 'water' in the text. Ambigous statements or something that can be misunderstood are strictly forbidden. Text that appears stylistically fitting, fluent and plausible, but could make the reader confused about what the authors actually mean must be actively purged from the manuscript, as it is incompatible with readability, transparency and scientific rigour.
* The above mandates the text to be concise. However, 'concise' must NOT be understood as being text short. Making text compact is a no-goal. Rather: ANY length is allowed, so long as the narrative is presented in appropriate level of detail and is optimised to be as concise as possible while (important!) maintaining that detail. 
* For clarity, concise text means no repetitions: enough to state a particular point once.
* The manuscript is written in reserved highly-precise and on-point language: NO colourful statements, overinterpretation, 'overselling' the paper etc - all these bad practices are strictly forbidden. Presented data must speak for itself, this manuscript presents solid research that does not require such bad presentation practices.

## 1. Filter specification (IMPERATIVE — scientific correctness)

Every reported identification count falls into one of two classes, each with EXACTLY one admissible filter and
nothing else, unless explicitly instructed in this document — no contaminant/target filter, no positive-quantity filter (zero
quantities are counted), no extra q-value gates. Decoys (`Decoy == 1`) are
always dropped (that removes the FDR null model; it is not one of the forbidden
"filters").

| Quantity | Unit | Admissible filter | 
|---|---|---|
| Protein groups, per run | within one run | `PG.Q.Value <= 0.01` |
| Precursors, per run | within one run | `Q.Value <= 0.01` |
| Protein groups, global (unique per dataset) | union | `Lib.PG.Q.Value <= 0.01` |
| Precursors, global (unique per dataset) | union | `Lib.Q.Value <= 0.01` |

A dataset TOTAL (union across runs) is a **global** number, so headline "total
precursors / protein groups" use `Lib.Q.Value` / `Lib.PG.Q.Value`. Per-cell and
per-run distributions, and >=3-run replicate counts, use `Q.Value` / `PG.Q.Value`.

Numbers comparisons between softwares based on the above are considered 'comparable' (see the above directive), except for plexDIA or phosphoproteomics, where they are not comparable (see below). It is further clarified that DIA-NN 1.8.1 pr_matrix.tsv report can be considered filtered at Q.Value <= 0.01 for per run comparisons and can be considered filtered at Lib.Q.Value <= 0.01 for global comparisons. Similar for pg_matrix.tsv generated by DIA-NN 1.8.1, just protein group level.

Canonical implementation: the `count_report` function in
[`scripts/rebuild.py`](scripts/rebuild.py) (formerly `analysis/count_report_ids.py`).
Emitted keys: `prec_min1`, `prec_min3` (Q.Value replicate), `prec_global`
(Lib.Q.Value), `prot_global` (Lib.PG.Q.Value), `prot_perrun_avg`, `prot_complete`
(PG.Q.Value), `peptides`, `prot_2pep` (global rule).

### plexDIA (MSV000093870) — channel is the per-run unit

A single cell = one channel of one run. The per-cell precursor q-value filter therefore involves `Channel.Q.Value` (the channel-level analog of `Q.Value`; run-level
`Q.Value` passes a precursor in all 3 channels and is not a per-cell number) in addition to `Q.Value` at 0.01.
Per-cell protein groups use `PG.Q.Value` which only applied to DIA-NN 2.5 and later output.
Implemented in [`scripts/rebuild.py`](scripts/rebuild.py) (the `plexdia_per_cell`
stage; functions `load_channel_confident` / `per_cell_counts`, formerly
`analysis/plexDIA/figure_msv000093870_oocyte_plexdia.py`).

CRITICAL: PG.Q.Value in DIA-NN versions before 2.0 does not reflect per-cell/channel confidence. Any comparisons of numbers obtained using this filter between versions is invalid and must not be included in this manuscript. Numbers of precursors at Q.Value <= 0.01 and Channnel.Q.Value <= 0.01, however, are considered comparable between DIA-NN 1.8.1 and 2.5/2.5.1.

### Phosphoproteomics
Phosphopeptides counted at the global precursor rule (`Lib.Q.Value <= 0.01`).
Phosphosite localisation uses the DIA-NN site report `Probability >= 0.75` — a localisation-confidence axis, not one of the forbidden q-filters. Since a single precursor mapping to multiple proteins gives rise to multiple sites, each identified as the protein accession combined with the site residue position within the protein, site counts must be reported only using first protein in the protein group (index within protein group is part of DIA-NN's site_report.parquet and can be used directly). This must be indicated to the reader appropriately.

### Original / deposited side of reanalysis comparisons
Deposited datasets provide only `*_pg_matrix.tsv` / `*_pr_matrix.tsv` (already
q-filtered count matrices, no q-value columns). The "original" bar is the matrix
row count as-is, respective 'virtual' equivalent filters are listed above; the rule governs the quantmsdiann (reanalysis) side.

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
| MSV000093870 | oocyte plexDIA | reanalysis | 2.5.1 |
| PXD003539 | NCI-60 bulk | reanalysis | 2.5.x |
| PXD030304 | ProCan bulk | reanalysis | 2.5.x |
| PXD004701 | Sun bulk | reanalysis | 2.5.x |
| PXD064049 | spatial DVP diaPASEF | reanalysis | 2.5.1-ent |
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
| `runtime_by_dataset` | fig | per-run wall-clock by dataset type (Supp Fig S3) | `parallelism_data.tsv` |
| `benchmarks` | fig | Fig 2 benchmark panels, `counts.tsv` | reads `report_counts.tsv` |
| `queue_sweep` | fig | `queue_size_sweep.tsv` (feeds `fig2_validation`) | runtime trace |
| `fig2_validation` | fig | Fig 2 validation composite | composes sweep + accuracy |
| `id_vs_epsilon` | fig | ProteoBench id-vs-epsilon panel | reads `report_counts.tsv` |
| `proteobench_accuracy` | fig | ProteoBench accuracy panels | HYE fold-changes |
| `reanalysis_improvement` | fig | reanalysis-recovery figure | original matrix vs reanalysis report |
| `single_cell_combined` | fig | single-cell figure | `single_cell_tables` outputs |
| `pxd003539` | fig | PXD003539 (NCI-60) panels | matrix vs report |
| `pxd004701` | fig | PXD004701 (Sun) panels | matrix vs report |
| `pxd030304` | fig | PXD030304 (ProCan) panels | streams a 2GB matrix |
| `pxd064049_spatial` | fig | spatial DVP panels | matrix-based (see spec Goal 5) |
| `atlas` | fig | pan-cohort (Fig S10) | needs numpy<2 env |
| `venn` | fig | protein-accession overlap (supp) | per-cohort |
| `performance_trace`, `mdc_cluster_runtime` | fig | Fig 1 / performance traces | runtime + resource traces |
| `paper_numbers` | num | `data/paper_numbers.tsv` + `paper/generated_numbers.tex` | aggregates ALL manuscript numbers |
| `paper/Makefile` | pdf | SVG->PDF + `main`/`supplementary` PDFs | rsvg-convert |

The orchestrator [`scripts/rebuild.py`](scripts/rebuild.py) runs all of the
above in dependency order; `--list` prints this table from the live stage
registry.

---

## 4. Audit targets

This section lists necessary checks for a structured audit. The audit must produce an audit artifact consisting of sections, each corresponding to an audit target. Each section lists an inventory of objects pertaining to the particular audit target, such as input data, scripts, individual algorithms and approaches, figures, figure legends, manuscript sections or paragraphs. The audit is first performed sequentially, section after section. The audit artifact is updated for each section and saved to disk. IMPERATIVE: the previous section audit must be saved to disk before commencing the next section, this rule must be followed under all circumstances and simultaneous work on multiple sections is forbidden. Second, an overall audit of the manuscript repository as a whole is performed, including previous findings. All findings are analysed together to produce a mitigations plan, if necessary.

If an audit finding cannot easily be fixed, e.g. conforming to a directive requrires completely removing a particular figure or refactoring it in a way that will change the narrative, the audit finding must be addressed by adding a pending item to the refactor targets below. This way the user has an opportunity to edit the task for this item, possibly after obtaining updated input data and repurposing the figure.

* Compliance with the critical directive on software comparisons.
* Compliance with filtering mandates.
* Compliance with manuscript writing mandates.
* Script correctness, maintainability and quality. Scripts must be absolutely minimal, with clean and easily auditable and maintainable code.
* Full data and scripts inventory is complete, correct and up-to-date in this methods.md.
* Audit for anything in the data/scripts where it is not obvious if it complies with the intent and not just the wording of mandatates in this methods.md. Could it be that the mandates were not formulated precise enough and an unintended scientifically dubious and possibly misleading benchmark/comparison gets included in the manuscript. Such instances - anything suspicious - must be flagged and explained in detail for the user to review.
* Protein inference sanity check. Any output from DIA-NN 1.8.1 or third-party software may contains inappropriately inferred protein groups: in this case the numbers of protein groups reported get significantly inflated, making them NOT comparable - such data must be directed to the user for reanalysis and refactor. Each script reading such data (e.g. DIA-NN 1.8.1 .tsv report of pg_matrix.tsv) must verify that among unique Protein.Group entries at 1% global q-value filter no more than 1.2% (known empirical cap) include multiple protein ids/accessions (multiple = Protein.Group string contains ';' symbol). 

## 5. Refactor targets

This section lists pending changes. Once a change is implemented, the list below must be updated accordingly.

* Reformat this methods.md document to use consistent text style, without artificial line breaks. Correct any typos.
* Rearrange figures order: the figure on quantmsdiann 2.5/2.5.1 comparison with public data must precede the single cell performance comparison with DIA-NN 1.8.1.
* Refactor phpsphoproteomics benchmarks to show only DIA-NN 2.5 or 2.5.1 numbers and runtimes but NOT comparing to other softwares.
