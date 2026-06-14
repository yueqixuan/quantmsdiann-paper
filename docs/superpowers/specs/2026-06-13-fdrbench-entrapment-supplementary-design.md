# Supplementary: FDRBench entrapment validation of quantms.io identifications

**Date:** 2026-06-13
**Status:** Draft for review (and to run by Vadim before committing cluster runs)

## Motivation

The main figures show 2.5.1 Enterprise (with `--kb`) recovering substantially more
precursors/proteins than 1.8.1, especially on single-cell data. A reviewer will ask the
obvious question — *are the extra identifications real, or false transfers from MBR?*
Vadim's guidance is to show the gains transparently **with a clear note on caveats**; this
analysis turns that caveat into evidence by measuring the true false-discovery proportion
(FDP) behind the nominal 1% q-value, using the entrapment method (Wen & Noble, Nat Methods
2025). It also directly addresses the false-transfer concern raised by the Krijgsveld lab
data (Nat Commun 2024, 10.1038/s41467-024-52605-x).

This is **supplementary**, framed as a selected, transparent spot-check — NOT an exhaustive
FDR benchmark, and NOT a reproduction of the Krijgsveld DIA-ME experiment (that is a
wet-lab Matching-Enhancer design we cannot regenerate).

## Tool

FDRBench (Wen & Noble), already vendored at
`databases/fdrbench-0.0.4/fdrbench-0.0.4.jar`. Two phases:
1. **Build entrapment DB:** `java -jar fdrbench-0.0.4.jar -db <human.fasta> -fold 1
   -method shuffle -diann -o <entrapment.fasta>` (combined-entrapment; shuffled decoy-style
   entrapment sequences appended to the target DB). Optionally `-ms <foreign.fasta>` for
   multi-species (foreign-proteome) entrapment as a cross-check.
2. **Compute FDP:** after searching against the entrapment DB, `java -jar fdrbench-0.0.4.jar
   -i <diann_report> -diann -r <#entrapment/#target> -o <fdp.tsv>` at
   precursor/peptide/protein level.

## Design

**Datasets (selected — single-cell, where FDR control is hardest and gains largest):**
- PXD046357 (HeLa Astral single-cell)
- PXD044991 (HeLa One-Tip)
- (optional) one ProteoBench module as a standard-load reference.

**Entrapment DB:** human reviewed + contaminants + **1-fold shuffled entrapment**
(`-fold 1 -method shuffle`). Same DB reused across all conditions for comparability.

**Conditions:** DIA-NN **1.8.1** vs **2.5.1 Enterprise** (`--kb` on), each ideally
**MBR on vs off**. MBR-off isolates the false-transfer component (the Krijgsveld FTR idea).

**Readout:** FDRBench FDP at precursor / peptide / protein level, evaluated at the nominal
1% q-value. Supplementary panel: grouped bars of FDP (dataset × version × MBR × level),
with a dashed line at the nominal 1% to show over/under-control.

## Workflow / files

1. Build entrapment FASTA on the cluster (compute node; FDRBench needs Java).
2. Re-run quantms.io against the entrapment FASTA for each dataset × version (× MBR) into a
   dedicated `_entrapment` results tree (do NOT pollute the canonical trees — see the
   `_v1_8_1` / `_v2_5_1_enterprise` tree-separation lesson).
3. Run FDRBench FDP on each `diann_report` → small `fdp.tsv` staged into
   `data/entrapment/`.
4. New script `analysis/figure_entrapment_fdp.py` → `supplementary/supp_entrapment_fdp.svg`.
5. Caveat paragraph in the discussion (draft already prepared from the lit review:
   filter MBR on `Lib.Q.Value` not `Global.Q.Value`; report protein-level as the
   conservative readout; note gains partly reflect software/version differences).

## Expected result (honest)

Consistent with the literature: **protein-level FDP near-nominal**; **peptide/precursor-level
FDP likely inflated**, more so with MBR and at single-cell input. That is itself the
transparent caveat — reported, not hidden.

## Open items / risks (resolve before launching runs)

1. **MBR-off feasibility is UNconfirmed.** The pipeline exposes `skip_preliminary_analysis`
   (calibration only), not a clean MBR/empirical-library-reanalysis toggle. The Indiv-vs-MBR
   contrast may need a small pipeline change (skip `assemble_empirical_library` +
   requantification) or approximating Indiv via per-file independent runs. **If MBR-off is
   not feasible quickly, run MBR-on FDP only** (still shows whether the final gains are
   FDR-controlled, just not the Indiv-vs-MBR decomposition).
2. **Cluster cost:** ~6–8 DIA-NN runs against the entrapment DB. Queue *after* the Enterprise
   single-cell runs finish; mind the recently-hit disk pressure on `/hps/nobackup`.
3. **Confirm with Vadim** before investing the runs — he may prefer citing Krijgsveld +
   Wen/Noble over an in-house entrapment analysis. He also noted FDP interpretation is
   challenging (entrapment assumptions, protein vs peptide level differ).

## Citations
- Wen & Noble, entrapment FDR assessment / FDRBench (Nat Methods 2025; biorxiv
  2024.06.01.596967).
- Krijgsveld lab DIA-ME, Nat Commun 2024, 10.1038/s41467-024-52605-x (false-transfer
  reference).

## Out of scope
- Reproducing the Krijgsveld DIA-ME FPR/FTR/ROC figure (their wet-lab design).
- Any production / canonical-tree runs (entrapment runs live in a throwaway `_entrapment`
  tree).
