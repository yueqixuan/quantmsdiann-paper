#!/usr/bin/env python
"""PXD049692 (NK-cell Fe-NTA phosphoproteomics, diaPASEF on timsTOF HT)
reanalysis: quantmsdiann (DIA-NN 2.5.0, library-free) versus the originally
deposited Spectronaut directDIA analysis on the identical 10 runs.

!!! STALE — DO NOT USE IN THE MANUSCRIPT YET (2026-06-12) !!!
The staged quantmsdiann run was produced with sdrf-pipelines 0.1.4, whose
convert-diann emitted `--var-mod Phospho,79.966331,S,T,Y` (comma bug) so DIA-NN
searched **serine only** — the matrix is 100% pS, with all pT/pY missing.
Comparing that S-only count to the Spectronaut S/T/Y report is apples-to-oranges.
A re-run on the current dev (sdrf-pipelines 0.1.5 -> `Phospho,...,STY`) is in
flight; once it lands, REPLACE this with a phosphosite-level figure (class-I
localised pS/pT/pY) per the manuscript focus. See memory:
reference_sdrf_diann_multiresidue_mod_bug.


Both analyses are library-free / directDIA, so the comparison is like-for-like.
The deposited Spectronaut report (`*_PH_Report.tsv`) is a fragment-level report
without a PTM-localisation-probability column, so a class-I *site* head-to-head
is not possible; the fair, engine-independent metric is the number of distinct
phosphopeptides (stripped sequences carrying a phospho modification) at 1% FDR.
On the shared runs quantmsdiann recovers more phosphopeptide backbones
(4,784 vs 4,254, +12.5%); the original Spectronaut report carries more distinct
phosphopeptidoforms (7,993 vs 5,196), reflecting its site-variant reporting.

main_comparison.svg -- phosphopeptides at 1% FDR, original (grey) vs
quantmsdiann (blue); the bar style and (7,5) canvas match the other Fig.~3
panels.

Run:  PYTHONPATH=. python -m analysis.figure_pxd049692_phospho_vs_quantmsdiann
"""
from __future__ import annotations

import csv
import sys
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analysis import figure_style as fs
fs.apply_house_style()
import pandas as pd

from analysis.contaminant_filter import is_target_protein_group

REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "PXD049692"
CACHE_DIR = FIGURES_DIR / "data" / "cache"

ORIG_COLOUR = fs.COMPARISON["original"]
QM_COLOUR = fs.COMPARISON["quantmsdiann"]

_QB = ("https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
       "quantmsdiann-benchmarks/PXD049692")
# Deposited Spectronaut directDIA phospho report (fragment-level, ~270 MB).
_ORIG_PH = ("https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/06/PXD049692/"
            "20231023_090251_2023-10-21_KA1_PH_Report.tsv")


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or dest.stat().st_size == 0:
        with urllib.request.urlopen(url, timeout=900) as r:
            dest.write_bytes(r.read())
    return dest


def _qm_run_names() -> set[str]:
    sdrf = _download(
        f"{_QB}/sdrf/PXD049692-NK-fibrin-IL15-phospho-diaPASEF.sdrf.tsv",
        CACHE_DIR / "qm.sdrf.tsv")
    rows = list(csv.reader(open(sdrf), delimiter="\t"))
    fcol = next(i for i, h in enumerate(rows[0]) if "data file" in h.lower())
    return {r[fcol].replace(".d", "") for r in rows[1:] if len(r) > fcol}


def quantmsdiann_phosphopeptides() -> tuple[int, int]:
    """(stripped phosphopeptide backbones, distinct phosphopeptidoforms)."""
    pr = _download(f"{_QB}/quant_tables/diann_report.pr_matrix.tsv",
                   CACHE_DIR / "qm_pr_matrix.tsv")
    df = pd.read_csv(pr, sep="\t", dtype=str)
    df = df[df["Protein.Group"].map(is_target_protein_group)]
    ph = df[df["Modified.Sequence"].fillna("").str.contains("UniMod:21")]
    return int(ph["Stripped.Sequence"].nunique()), int(ph["Modified.Sequence"].nunique())


def original_phosphopeptides(qm_runs: set[str]) -> tuple[int, int]:
    """(stripped phosphopeptide backbones, distinct phosphopeptidoforms) in
    the deposited Spectronaut directDIA report, on the shared runs, at 1%
    precursor q-value."""
    report = _download(_ORIG_PH, CACHE_DIR / "orig_PH_Report.tsv")
    f = csv.reader(open(report, encoding="utf-8", errors="replace"), delimiter="\t")
    ix = {c: i for i, c in enumerate(next(f))}
    iFN, iMod, iStrip, iQ = (ix["R.FileName"], ix["EG.ModifiedSequence"],
                             ix["PEP.StrippedSequence"], ix["EG.Qvalue"])
    strp: set[str] = set()
    forms: set[str] = set()
    for r in f:
        if len(r) <= iQ or r[iFN] not in qm_runs:
            continue
        try:
            if float(r[iQ]) > 0.01:
                continue
        except ValueError:
            continue
        if "Phospho" in r[iMod]:
            strp.add(r[iStrip])
            forms.add(r[iMod])
    return len(strp), len(forms)


def render_main_comparison(orig: int, qm: int, svg_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar([0, 1], [orig, qm], width=0.55,
                  color=[ORIG_COLOUR, QM_COLOUR], edgecolor="#37474f")
    ax.set_xlim(-0.7, 1.7)
    for b, v in zip(bars, (orig, qm)):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center",
                va="bottom", fontsize=11)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Original\n(Spectronaut directDIA)",
                        "quantmsdiann\n(DIA-NN)"])
    ax.set_ylabel("Phosphopeptides (stripped sequence, 1% FDR)")
    ax.set_ylim(0, max(orig, qm) * 1.16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")  # SVG-only (repo convention)
    plt.close(fig)


def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    qm_runs = _qm_run_names()
    qm_strip, qm_forms = quantmsdiann_phosphopeptides()
    orig_strip, orig_forms = original_phosphopeptides(qm_runs)
    # The panel compares stripped-sequence phosphopeptide backbones (parity).
    # Phosphopeptidoform (modified-sequence) counts are recorded too: at that
    # site-resolved level the two analyses diverge (Spectronaut reports more
    # positional isomers), so the manuscript frames this as backbone-level,
    # not site-level, parity.
    render_main_comparison(orig_strip, qm_strip, FIGURES_DIR / "main_comparison.svg")
    (FIGURES_DIR / "counts.tsv").write_text(
        "metric\toriginal_spectronaut\tquantmsdiann_diann250\n"
        f"phosphopeptides_stripped\t{orig_strip}\t{qm_strip}\n"
        f"phosphopeptidoforms_modified\t{orig_forms}\t{qm_forms}\n")
    print(f"stripped phosphopeptides: original={orig_strip}  quantmsdiann={qm_strip}")
    print(f"phosphopeptidoforms:      original={orig_forms}  quantmsdiann={qm_forms}")
    print(f"wrote {FIGURES_DIR}/main_comparison.svg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
