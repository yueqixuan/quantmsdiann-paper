#!/usr/bin/env python
"""Graphical abstract for the quantmsdiann manuscript.

A simple three-stage banner conveying the paper's three messages:
  1. Metadata standards  -- diverse public DIA datasets, each described by one
     standardised SDRF, are the single input.
  2. Scalability         -- one SDRF-driven DIA-NN workflow parallelised across
     HPC/cloud (archive-scale cohorts in hours).
  3. Integration         -- every modality is analysed by the *same* workflow
     and exported to one harmonised, queryable (QPX/MSstats) output.

Run:  PYTHONPATH=. python -m analysis.figure_graphical_abstract
Output: analysis/figures/manuscript/graphical_abstract.{svg,pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "analysis" / "figures" / "manuscript"

QMBLUE = "#1565C0"
QMGREEN = "#2E7D32"
DARK = "#1A1A1A"
GREY = "#555555"
LIGHTBLUE = "#E8F1FB"
LIGHTGREEN = "#E8F3E9"
LIGHTGREY = "#EFEFEF"
FONT = "DejaVu Sans"


def _panel(ax, x, y, w, h, face, edge):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=1.2",
                       linewidth=1.2, edgecolor=edge, facecolor=face, zorder=1)
    ax.add_patch(p)


def _chip(ax, x, y, w, h, label, face, edge, fontsize=10.5, color=DARK, bold=False):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.8",
                       linewidth=1.0, edgecolor=edge, facecolor=face, zorder=2)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
            fontsize=fontsize, color=color, fontweight="bold" if bold else "normal",
            family=FONT, zorder=3)


def _arrow(ax, x0, x1, y):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                                 mutation_scale=22, linewidth=2.4,
                                 color=QMGREEN, zorder=4))


def build() -> Path:
    fig, ax = plt.subplots(figsize=(11.2, 4.7))
    ax.set_xlim(0, 112)
    ax.set_ylim(0, 47)
    ax.axis("off")

    # Title strip.
    ax.text(2, 44.2, "quantmsdiann", fontsize=20, fontweight="bold",
            color=QMBLUE, family=FONT)
    ax.text(2, 40.2, "one SDRF-driven DIA-NN workflow for archive-scale "
            "reanalysis of public DIA proteomics", fontsize=11, color=GREY,
            family=FONT)

    # ---- Stage 1: standardised input ----
    _panel(ax, 1.5, 5, 31, 31, LIGHTGREEN, QMGREEN)
    ax.text(17, 32.5, "Standardised input", ha="center", fontsize=12.5,
            fontweight="bold", color=QMGREEN, family=FONT)
    # SDRF sheet glyph
    _chip(ax, 9.5, 26.5, 15, 4.0, "SDRF metadata", "#ffffff", QMGREEN,
          fontsize=11, color=QMGREEN, bold=True)
    # dataset-type chips
    datasets = ["Bulk cell lines", "Single cell (plexDIA)", "Spatial (DVP)",
                "Phosphoproteomics", "timsTOF / Astral / Orbitrap / ZenoTOF"]
    cy = 22.0
    for d in datasets:
        _chip(ax, 3.5, cy, 27, 3.0, d, "#ffffff", "#bbbbbb", fontsize=9.5)
        cy -= 3.6
    ax.text(17, 6.3, "diverse public DIA datasets,\none SDRF each",
            ha="center", fontsize=9, color=GREY, family=FONT, style="italic")

    _arrow(ax, 33.5, 39.5, 20.5)

    # ---- Stage 2: one scalable workflow ----
    _panel(ax, 40, 5, 31, 31, LIGHTBLUE, QMBLUE)
    ax.text(55.5, 32.5, "One scalable workflow", ha="center", fontsize=12.5,
            fontweight="bold", color=QMBLUE, family=FONT)
    _chip(ax, 46, 26.5, 19, 4.2, "nf-core $\\cdot$ DIA-NN", "#ffffff", QMBLUE,
          fontsize=11, color=QMBLUE, bold=True)
    # parallel cluster nodes
    nx = 44.5
    for _ in range(6):
        _chip(ax, nx, 19.5, 3.0, 4.0, "", QMBLUE, QMBLUE)
        nx += 3.9
    ax.text(55.5, 16.2, "parallel \\texttt{DIA-NN} across HPC / cloud nodes"
            .replace("\\texttt{", "").replace("}", ""),
            ha="center", fontsize=9.5, color=DARK, family=FONT)
    _chip(ax, 43.5, 9.5, 24, 4.6, "2,300 files in 2.4 h\non 300 nodes",
          "#ffffff", QMBLUE, fontsize=10.5, color=QMBLUE, bold=True)

    _arrow(ax, 72, 78, 20.5)

    # ---- Stage 3: harmonised, integrated output ----
    _panel(ax, 78.5, 5, 31, 31, LIGHTGREEN, QMGREEN)
    ax.text(94, 32.5, "Harmonised output", ha="center", fontsize=12.5,
            fontweight="bold", color=QMGREEN, family=FONT)
    for i, lab in enumerate(["QPX (Parquet) archive", "MSstats input",
                             "pmultiqc QC report"]):
        _chip(ax, 81, 26.0 - i * 4.4, 26, 3.4, lab, "#ffffff", "#bbbbbb",
              fontsize=9.8)
    ax.text(94, 9.0, "one integrated, queryable,\nFAIR matrix across all "
            "modalities", ha="center", fontsize=9, color=GREY, family=FONT,
            style="italic")

    # bottom concept strip
    ax.text(56, 1.6, "Metadata standards   $\\bullet$   Scalability   "
            "$\\bullet$   Integration", ha="center", fontsize=11.5,
            fontweight="bold", color=DARK, family=FONT)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUT_DIR / "graphical_abstract"
    fig.tight_layout(pad=0.3)
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return stem.with_suffix(".svg")


if __name__ == "__main__":
    print("wrote", build())
