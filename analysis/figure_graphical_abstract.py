#!/usr/bin/env python
"""Graphical abstract for the quantmsdiann manuscript.

A visual three-act story conveying the paper's three messages:
  1. Metadata standards -- diverse public DIA datasets (bulk, single-cell,
     spatial, phospho), each described by one standardised SDRF (the green
     "spine" binding the dataset cards).
  2. Scalability        -- one nf-core/DIA-NN engine parallelised across HPC/
     cloud, with the *real* queue-size scaling curve embedded (37.7 h -> 2.4 h
     from 10 to 300 nodes on the 2,300-file single-cell cohort).
  3. Integration        -- every modality flows into one harmonised, queryable
     QPX/MSstats archive (the database cylinder).

Run:  PYTHONPATH=. python -m analysis.figure_graphical_abstract
Output: analysis/figures/manuscript/graphical_abstract.{svg,pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse, Rectangle

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "analysis" / "figures" / "manuscript"

QMBLUE = "#1565C0"
QMGREEN = "#2E7D32"
DARK = "#1A1A1A"
GREY = "#5A5A5A"
LIGHTBLUE = "#EAF2FB"
LIGHTGREEN = "#E9F4EA"
FONT = "DejaVu Sans"

# queue-size sweep (PXD071075): nodes -> wall-clock hours (real data).
SWEEP_NODES = [10, 50, 100, 200, 300]
SWEEP_HOURS = [37.7, 8.1, 4.8, 2.6, 2.4]

# distinct soft colours for the dataset "cards" (diversity of modalities).
CARD_COLOURS = ["#1e88e5", "#00897b", "#8e24aa", "#e8a000", "#6d4c41"]
CARD_LABELS = [
    "Bulk cell lines",
    "Single cell · plexDIA",
    "Spatial · DVP",
    "Phosphoproteomics",
    "timsTOF · Astral · Orbitrap · ZenoTOF",
]


def _round(ax, x, y, w, h, face, edge, lw=1.2, rounding=1.4, z=1, alpha=1.0):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rounding}",
        linewidth=lw, edgecolor=edge, facecolor=face, zorder=z, alpha=alpha))


def _cylinder(ax, cx, top, bot, w, face, edge):
    rx, ry = w / 2.0, 1.7
    # body
    ax.add_patch(Rectangle((cx - rx, bot), w, top - bot, facecolor=face,
                           edgecolor="none", zorder=2))
    ax.plot([cx - rx, cx - rx], [bot, top], color=edge, lw=1.4, zorder=3)
    ax.plot([cx + rx, cx + rx], [bot, top], color=edge, lw=1.4, zorder=3)
    # bottom + disk separators + top
    ax.add_patch(Ellipse((cx, bot), w, 2 * ry, facecolor=face, edgecolor=edge,
                         lw=1.4, zorder=2))
    for yy in (bot + (top - bot) * 0.34, bot + (top - bot) * 0.67):
        ax.add_patch(Ellipse((cx, yy), w, 2 * ry, facecolor="none",
                             edgecolor=edge, lw=0.7, alpha=0.5, zorder=4))
    ax.add_patch(Ellipse((cx, top), w, 2 * ry, facecolor="#ffffff",
                         edgecolor=edge, lw=1.4, zorder=5))


def _arrow(ax, x0, x1, y):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                                 mutation_scale=26, linewidth=3.0,
                                 color=QMGREEN, zorder=6))


def build() -> Path:
    fig, ax = plt.subplots(figsize=(11.6, 5.0))
    ax.set_xlim(0, 116)
    ax.set_ylim(0, 50)
    ax.axis("off")

    # ---- title ----
    ax.text(2.5, 46.6, "quantmsdiann", fontsize=20, fontweight="bold",
            color=QMBLUE, family=FONT)
    ax.text(2.5, 42.7, "one SDRF-driven DIA-NN workflow for archive-scale "
            "reanalysis of public DIA proteomics", fontsize=10.5, color=GREY,
            family=FONT)

    # ---- ACT 1: diverse datasets bound by one SDRF schema ----
    ax.text(20, 37.6, "Diverse DIA datasets, one SDRF", ha="center",
            fontsize=11.5, fontweight="bold", color=QMGREEN, family=FONT)
    # SDRF spine (the standardised-metadata "binding")
    _round(ax, 3.0, 9.5, 4.6, 24.5, QMGREEN, QMGREEN, rounding=1.2, z=2)
    ax.text(5.3, 22, "SDRF", rotation=90, ha="center", va="center",
            fontsize=12, fontweight="bold", color="#ffffff", family=FONT)
    # dataset cards fanning off the spine
    n = len(CARD_LABELS)
    cw, ch, top_y, step = 26.0, 3.9, 30.0, 5.0
    for i, (lab, col) in enumerate(zip(CARD_LABELS, CARD_COLOURS)):
        cy = top_y - i * step
        ax.plot([7.6, 9.5], [cy + ch / 2, cy + ch / 2], color=col, lw=1.3,
                zorder=1, alpha=0.8)
        _round(ax, 9.5, cy, cw, ch, "#ffffff", "#c7c7c7", lw=1.0, rounding=0.7, z=3)
        ax.add_patch(Rectangle((9.5, cy + 0.25), 1.5, ch - 0.5, facecolor=col,
                               edgecolor="none", zorder=4))
        ax.text(12.3, cy + ch / 2, lab, ha="left", va="center", fontsize=9.4,
                color=DARK, family=FONT)

    _arrow(ax, 37.0, 42.5, 21.5)

    # ---- ACT 2: one scalable engine (with the real scaling curve) ----
    _round(ax, 43.0, 8.5, 32.0, 30.0, LIGHTBLUE, QMBLUE, rounding=1.6, z=1)
    ax.text(59, 35.2, "One scalable workflow", ha="center", fontsize=11.5,
            fontweight="bold", color=QMBLUE, family=FONT)
    ax.text(59, 31.9, "nf-core $\\cdot$ DIA-NN, parallel on HPC / cloud",
            ha="center", fontsize=9.3, color=GREY, family=FONT)
    # embedded scaling curve (hand-plotted in data coords)
    bx0, bx1, by0, by1 = 48.0, 71.0, 13.0, 29.0
    hmin, hmax = min(SWEEP_HOURS), max(SWEEP_HOURS)
    def _px(node): return bx0 + (node - SWEEP_NODES[0]) / (SWEEP_NODES[-1] - SWEEP_NODES[0]) * (bx1 - bx0)
    def _py(h): return by0 + (h - hmin) / (hmax - hmin) * (by1 - by0)
    ax.plot([bx0, bx0], [by0, by1 + 1], color="#b0bec5", lw=1.0, zorder=2)
    ax.plot([bx0, bx1], [by0, by0], color="#b0bec5", lw=1.0, zorder=2)
    xs = [_px(nd) for nd in SWEEP_NODES]
    ys = [_py(h) for h in SWEEP_HOURS]
    ax.plot(xs, ys, "-o", color=QMBLUE, lw=2.4, ms=5, mfc=QMBLUE,
            mec="#ffffff", mew=0.8, zorder=4)
    ax.text(xs[0] + 0.6, ys[0] + 0.3, "37.7 h", fontsize=9, color=QMBLUE,
            fontweight="bold", ha="left", va="bottom", family=FONT)
    ax.text(xs[-1], ys[-1] + 1.3, "2.4 h", fontsize=9, color=QMBLUE,
            fontweight="bold", ha="right", va="bottom", family=FONT)
    ax.text(bx0 - 0.5, by0 - 1.3, "10", fontsize=8, color=GREY, ha="center", va="top")
    ax.text(bx1, by0 - 1.3, "300 nodes", fontsize=8, color=GREY, ha="right", va="top")
    ax.text((bx0 + bx1) / 2, by1 + 1.6, "wall-clock vs cluster width "
            "(2,300-file cohort)", ha="center", fontsize=8.2, color=GREY,
            family=FONT, style="italic")

    _arrow(ax, 75.5, 81.0, 21.5)

    # ---- ACT 3: one harmonised, integrated output ----
    ax.text(95.5, 37.6, "Harmonised output", ha="center", fontsize=11.5,
            fontweight="bold", color=QMGREEN, family=FONT)
    _cylinder(ax, 95.5, 31.5, 17.5, 19.0, LIGHTGREEN, QMGREEN)
    ax.text(95.5, 26.5, "QPX", ha="center", va="center", fontsize=15,
            fontweight="bold", color=QMGREEN, family=FONT)
    ax.text(95.5, 22.7, "Parquet archive", ha="center", va="center",
            fontsize=8.6, color=GREY, family=FONT)
    ax.text(95.5, 19.2, "+ MSstats $\\cdot$ pmultiqc QC", ha="center",
            va="center", fontsize=8.6, color=DARK, family=FONT)
    ax.text(95.5, 11.4, "one integrated, queryable matrix\nacross all "
            "modalities", ha="center", fontsize=8.8, color=GREY,
            family=FONT, style="italic")

    # ---- bottom pillars ----
    ax.text(58, 4.0, "Metadata standards      $\\bullet$      Scalability"
            "      $\\bullet$      Integration", ha="center", fontsize=12,
            fontweight="bold", color=DARK, family=FONT)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUT_DIR / "graphical_abstract"
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return stem.with_suffix(".svg")


if __name__ == "__main__":
    print("wrote", build())
