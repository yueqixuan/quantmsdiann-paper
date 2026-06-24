# Audit artifact — methods.md rule compliance (code + outputs)

Date: 2026-06-24. Scope: every rule in `methods.md` verified against `scripts/rebuild.py`
and the committed outputs (figure-data TSVs, `data/paper_numbers.tsv`,
`paper/main.tex`, `paper/supplementary.tex`). Numbers were **recomputed from the
on-disk caches** (plexDIA report, phospho reports, bulk cell-line matrices). Method:
8 audit clusters mapped to the §4 audit targets; every finding re-derived by an
independent skeptic (60 agents total). Refuted findings are excluded.

Legend: **CRITICAL** = wrong number, impermissible comparison, or non-reproducible
output reaching the manuscript. **MAJOR** = filter/rule deviation or mandated check
missing. **MINOR** = documentation/label staleness, no number affected.

---

## Target 1 — Compliance with the critical directive on software comparisons

**1.1 (CRITICAL) plexDIA protein-group comparison is impermissible.** methods.md §1
line 75 (frozen): "Any comparisons of [PG.Q.Value] numbers between versions is invalid
and must not be included … Numbers of precursors at Q.Value ≤ 0.01 and Channel.Q.Value
≤ 0.01, however, are considered comparable." Only **precursors** are comparable for
plexDIA; protein-group counts are not. The manuscript nonetheless makes an explicit
cross-version plexDIA protein-group comparison: main.tex:153 "3,328 protein groups …
versus 2,122 … (+57%)", "recovered 94% (2,001 of 2,122)", "+1,326 groups"; Fig 4a
(`main_galatidou_total_pg.svg`, `main_galatidou_comparison.svg`); the abstract's +57%
upper bound. Per the directive ("must be removed") this comparison must be removed or
restricted to precursors. Code: rebuild.py:7224-7263, 7199-7222.

**1.2 (MAJOR) phospho cross-software comparison sits against §1 line 58 + selective
reporting.** main.tex:163 / supp panel C: "23% more phosphopeptide backbones
(5,240 vs 4,254)" vs the deposited Spectronaut directDIA report. Line 58 is categorical:
phospho §1 numbers "are not comparable" between softwares. Line 24 (sites) is strictly
satisfied (backbones are not localized sites), but the comparability claim conflicts with
line 58, and the refactor target (line 187) says "show only DIA-NN 2.5/2.5.1 numbers …
NOT comparing to other softwares." Additionally **selective reporting**: the same
`PXD049692/counts.tsv` shows `phosphopeptidoforms_modified` deposited = 7,993 vs
quantmsdiann = 5,737 (quantmsdiann lower) — the manuscript reports only the backbone
metric where quantmsdiann wins. Naming the tool is user-permitted in supplementary;
the comparability claim and selective reporting are the issue.

**1.3 (MAJOR) FDR-control caveat absent (directive line 23).** The directive MANDATES
the manuscript explicitly state: DIA-NN 2.5/2.5.1 assumed to have largely correct FDR
control; prior versions/other softwares may have far worse true FDR (10%+) at the same
nominal 1%; verification out of scope. `grep` confirms no such statement anywhere in
main.tex or supplementary.tex. Every "+X%" recovery comparison (main.tex:145,151,153,
160,163; supp:341-347) is presented as if FDR-equivalent. Must add the caveat.

**1.4 (permitted, with caveat) bulk/spatial reanalysis count comparisons** (NCI-60,
ProCan, Sun, HeLa, spatial). §1 line 58 designates these §1-rule count comparisons as
comparable (non-plexDIA, non-phospho); original side = deposited matrix-row count
(lines 82-84). Permissible. Should be paired with the 1.3 FDR caveat and a one-line note
that protein-group grouping conventions differ between OpenSWATH/PCT-SWATH and DIA-NN.

---

## Target 2 — Compliance with filtering mandates (§1)

**Core counting (`count_report`, rebuild.py:530-575): fully compliant.** All nine emitted
keys recompute correctly from a cached report: prec_min1/min3 = Q.Value≤0.01 (per-run,
≥1/≥3 runs); prec_global = Lib.Q.Value≤0.01; prot_global = Lib.PG.Q.Value≤0.01;
prot_perrun_avg/prot_complete = PG.Q.Value≤0.01; peptides = distinct stripped among
Lib.Q.Value-passing; prot_2pep = global PG with ≥2 stripped peptides. Decoys dropped;
no contaminant/target/quantity filter (verified: 23 contaminant-prefixed PGs ARE counted,
proving no target filter). `is_target_protein_group`/`strip_known_prefix` are used only in
entrapment-FDR diagnostics and accession-overlap/atlas set logic, never in an
identification count. Constants flat 0.01 across versions. **No action.**

**2.1 (CRITICAL) plexDIA `qm_total` — global PG counted with NO q-value filter.**
rebuild.py:7227 `qm_total = confident['Protein.Group'].nunique()` on the decoy-dropped
frame. `REPORT_COLUMNS` (7309) does not even load `Lib.PG.Q.Value`. Recompute on cache:
no-filter = **4,217**; Channel.Q.Value≤0.01 = 4,137; PG.Q.Value≤0.01 = 2,621;
Lib.PG.Q.Value≤0.01 (§1 rule) = **2,500**. Committed `protein_groups_total` = 3,328
matches **none** of these. Decoy-drop-only counting is explicitly a §1 violation
(lines 43-45). (Mooted if 1.1 removes the PG comparison; otherwise must use Lib.PG.Q.Value.)

**2.2 (MAJOR) plexDIA per-cell precursors — missing the `Q.Value` gate.** methods.md
line 69: per-cell precursor filter is Channel.Q.Value≤0.01 "in addition to Q.Value at
0.01" (BOTH). rebuild.py:7350 applies Channel.Q.Value only; `Q.Value` is not loaded.
Recompute: median 17,695 (Channel only) vs **16,062** (both) — ~10% inflation.
Fix: add `Q.Value` to REPORT_COLUMNS; gate on both. Also supplementary.tex:125 documents
only Channel.Q.Value.

**2.3 (MAJOR) plexDIA `write_counts` unique_protein_groups/unique_precursors — no
filter + mislabeled.** rebuild.py:7389-7390 count nunique with no q-filter, labelled
`source='target, channel-confident'` although `load_channel_confident` explicitly DROPS
the target filter. Recompute: §1 globals are Lib.PG.Q.Value≤0.01 = 2,500 PG /
Lib.Q.Value≤0.01 = 29,433 precursors; committed counts.tsv = 3,328 / 33,133 (wrong filter
AND stale). Fix filter + relabel.

**2.4 (MAJOR) phospho sites not restricted to first protein in group.** methods.md line
79: "site counts must be reported only using first protein in the protein group (index …
in site_report.parquet)." rebuild.py:6854-6857 dedups `(Protein, Site)` over ALL proteins;
it never reads `Protein.Index.In.Group`. Committed `sites_classI`/`sites_all` equal the
all-protein recompute (confirmed). Inflation is small but systematic (e.g.
biological-study classI 16,745 all-prot vs 16,700 first-prot; highspeed classI 11,615 vs
11,536). Fix: read `Protein.Index.In.Group`, filter `== 0` before counting; add the
reader-facing note to the supp caption.

**2.5 (suspicious→flag) Sun "consistency filter" adds a detection-fraction gate.**
rebuild.py:5404-5472: Sun reanalysis (7,561) gates Proteotypic==1 AND Lib.Q.Value≤0.01
AND ≥10%-detection — three deviations from the §1 admissible global rule
(Lib.PG.Q.Value≤0.01 alone). It is deliberately like-for-like with Sun's own published
consistency filter, which §1 line 84 / Supp Note 3 arguably permit, but §1 does not
explicitly bless a "consistency filter" as admissible. Number recomputes correctly (7,561).
Flag for user: confirm §1 authorizes original-applied thresholds like-for-like.

**2.6 (suspicious→flag) spatial reanalysis side = matrix-row count, not §1 report rule.**
rebuild.py:6354-6357 computes the PXD064049 2.5.1-ent reanalysis PG (3,099) and precursors
(20,705) via `count_matrix_rows`, not the §1 Lib.PG.Q.Value/Lib.Q.Value report rule. No
report parquet is cached for spatial; the 2.5.x pg_matrix is baked at `--matrix-spec-q
0.05`, so its row count is NOT the Lib.PG.Q.Value≤0.01 number. Fig 4 caption claims "the
quantmsdiann side is counted from the precursor report (parquet)" — contradicted for
spatial. §1 line 58 grants matrix-row equivalence only to DIA-NN 1.8.1, not 2.5.x. Fix:
compute from report parquet, or amend caption + confirm §1 permits the 2.5.x equivalence.

---

## Target 3 — Compliance with manuscript writing mandates

**3.1 (CRITICAL) atlas accession numbers in main.tex:161 are stale.** Prose: "12,857
unique accessions, of which 5,506 (42.8%) … in all five … 1,528 (11.9%) in exactly four."
Generated `combined_counts.tsv` (Panel C): union 12,914; all-five 5,716 (44.3%); exactly-
four 1,472 (11.4%). Hardcoded literals drifted from the regenerated figure. Fix prose to
the generated values (or regenerate atlas if 12,857 was intended).

**3.2 (MINOR) "28 pan-cancer tissue categories" vs 29 generated.** main.tex:162 says 28;
atlas Panel B/A both render 29 distinct tissue labels (one is "Healthy (Non-cancer)").
Mapping 29→28 is undocumented. State which bin is excluded, or use 29.

**3.3 (covered in 1.1/1.2/1.3)** plexDIA +57%, phospho "23% more", and the missing FDR
caveat are also writing-mandate issues (overstated/non-permitted comparisons).

---

## Target 4 — Script correctness, maintainability, quality

`count_report` is clean and exactly implements §1. The plexDIA path
(`load_channel_confident` / `per_cell_counts` / `qm_total` / `write_counts`) is the
locus of the filter bugs (2.1-2.3): it loads a column set that makes the §1 global rule
impossible and counts several totals with no q-value gate. The phospho `_count`
(2.4) omits the first-protein restriction. Recommended: route every plexDIA/phospho count
through an explicit, single admissible filter, mirroring `count_report`'s discipline.

---

## Target 5 — Inventory completeness/correctness in methods.md (§4 line 177)

**5.1 (MINOR) §3 stage table lines 140 & 155 say "site Prob>=0.99"** — contradicts §1
line 79 and the code (0.75). Doc only; no number affected. Fix → 0.75.

**5.2 (MINOR) §2 line 105 plexDIA version "2.5.0"** — actual is 2.5.1 (FTP path
`…/MSV000093870/v2_5_1`, report header "DIA-NN 2.5.1 Academia", commit fb5db78). Fix → 2.5.1.

**5.3 (MINOR) §2 line 109 spatial PXD064049 version "2.5.0"** — actual is 2.5.1-enterprise
(FTP path `…/spatial/PXD064049/v2_5_1_enterprise`, figure legend). Fix → 2.5.1-ent.

**5.4 (flag) bulk cohort version provenance.** §2 lines 106-108 list NCI-60/ProCan/Sun as
"2.5.x" (wildcard). The reanalysis figure labels them 2.5.1-enterprise/2.5.1/2.5.1. The
skeptic established the count-producing fetches use `_CELL_LINE_QT_BASE` =
v2_5_1_enterprise (PXD003539) / v2_5_1 (PXD030304/PXD004701); the `CELL_LINE_ANALYSES`
`version='v2_5_0'` field (rebuild.py:4602) is dead metadata (only dict keys are read).
Recommend pinning §2 to the actual point releases and removing the dead `version` field.

---

## Target 6 — Intent (not just wording): suspicious/misleading comparisons

**6.1** The plexDIA PG comparison (1.1/2.1) is the clearest intent violation: it is the
abstract's headline +57% upper bound, is impermissible per §1 line 75, and is
non-reproducible (committed 3,328 reproduces from neither the current code nor any
admissible filter).

**6.2** The phospho backbone comparison (1.2) is selective: backbone metric reported
(quantmsdiann wins), peptidoform metric omitted (quantmsdiann loses 5,737 vs 7,993).

**6.3 (MAJOR) reanalysis_improvement.tsv is a hand-maintained static file.** No code
writes it (grep shows reads only). The entire Fig 4 + the abstract's +24%–+57% range and
~38% mean derive from hand-typed constants, defeating the §3 "FTP report → count_report →
figure-data TSV" mandate and `rebuild --all` reproducibility. The plexDIA row (3,328) is
non-reproducible; the HeLa original (3,903) is inconsistent with `sc_totals` (4,137) used
elsewhere; NCI-60 6,433 was hand-edited (git) and recomputes to 6,553 (>=2-pep) / 7,553
(Lib rule) — the source parquet is absent so 6,433 is currently unverifiable. Generate
this TSV from a stage.

---

## Target 7 — Protein-inference sanity (§4 line 179: ≤1.2% multi-accession)

The guard exists only inside `count_report` (rebuild.py:549-558). It is **missing** from
every matrix/accession-reading path that ingests a pg/pr_matrix:
`count_matrix_rows` (431), `_report_global_protein_groups` (747), `pxd003539/017199/041421
_protein_accessions` + per-cell-line (1007-1237), `pxd003539_gene_symbols` (1239),
`pxd030304/004701` JSON-accession paths (1030-1057, 1260-1280). **No cohort actually
breaches** the cap (measured fractions: PXD030304 0.323%, PXD004701 0.478%, PXD017199
0.342%, PXD041421 0.385%, PXD003539 0.441-0.452%), so no number is currently wrong — but
the mandated check is absent and `extract_accessions_diann` silently splits ';', exactly
the inflation the guard targets. Fix: add a shared `assert_protein_inference_ok(...)`
guard at every matrix/pr/pg ingestion.

---

## Whole-repository synthesis + mitigation plan

**Confirmed wrong/impermissible numbers reaching the manuscript:**
1. plexDIA PG comparison (+57%, 3,328 vs 2,122, 94%, 1,326) — impermissible per §1 L75
   AND non-reproducible. → remove or restrict to precursors (1.1/2.1/2.3/6.1).
2. atlas accessions 12,857/5,506/1,528 → 12,914/5,716/1,472 (3.1).
3. plexDIA per-cell precursor median inflated ~10% by missing Q.Value gate (2.2).
4. phospho site counts inflated <1% by all-protein counting (2.4).

**Mandated additions (directive/§4):**
5. FDR-control caveat (1.3).
6. ≤1.2% multi-accession guard on all matrix readers (Target 7) — no number changes.

**Provenance/reproducibility refactors:**
7. Generate reanalysis_improvement.tsv from a stage; reconcile HeLa baseline; re-stage
   PXD003539 parquet to verify NCI-60 6,433 (6.3).
8. spatial reanalysis from report parquet or caption fix (2.6); Sun consistency-filter
   admissibility (2.5).

**Doc fixes (no number):** §3 0.99→0.75; §2 versions 2.5.0→2.5.1/2.5.1-ent; "28"→document
or 29; dead `version` field.

**Decisions requiring the user** (narrative-affecting): scope of the plexDIA reframe; how
to handle the phospho-vs-Spectronaut comparison (disclose peptidoform metric / add
non-comparability caveat / remove); whether to regenerate reanalysis_improvement.tsv now.

**Update (resolved):** NCI-60 original 4,284 is now derived from the deposited OpenSWATH
matrix on PRIDE FTP (`feature_alignment_requant_matrix.tsv`, proteins with >=2 distinct
peptides among quantified rows; `count_openswath_proteins_min_peptides`), not a hand-typed
constant. The reanalysis-side count (6,812) likewise comes from the DIA-NN report parquet
on FTP. Spatial-cohort inference inflation (>1.2% multi-accession) remains an open
follow-up for re-analysis/caveat.
