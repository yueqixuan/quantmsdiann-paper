#!/usr/bin/env python
"""Graphical abstract for the quantmsdiann manuscript (boxless, flowing design).

A circular convergence story conveying the three messages:
  1. Metadata standards -- diverse public DIA datasets converge through one
     standardised SDRF "gate".
  2. Scalability        -- into a single nf-core/DIA-NN engine hub ringed by
     cluster nodes, with the real queue-size scaling curve (37.7 h -> 2.4 h,
     10 -> 300 nodes on the 2,300-file cohort) beneath it.
  3. Integration        -- and out to one harmonised, queryable QPX archive.

Run:  PYTHONPATH=. python -m analysis.figure_graphical_abstract
Output: analysis/figures/manuscript/graphical_abstract.{svg,pdf,png}
"""
from __future__ import annotations

import math
from pathlib import Path as FsPath

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle, Ellipse, FancyArrowPatch

# Soft drop shadow for the main nodes (depth without boxes).
SHADOW = [pe.withSimplePatchShadow(offset=(1.8, -1.8), shadow_rgbFace="#7a7a7a",
                                   alpha=0.30)]

REPO_ROOT = FsPath(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "analysis" / "figures" / "manuscript"

QMBLUE = "#1565C0"
QMBLUE_L = "#dceaf8"  # light blue fill
QMGREEN = "#2E7D32"
QMGREEN_L = "#e3f1e4"
DARK = "#1A1A1A"
GREY = "#5A5A5A"
FONT = "DejaVu Sans"

SWEEP_NODES = [10, 50, 100, 200, 300]
SWEEP_HOURS = [37.7, 8.1, 4.8, 2.6, 2.4]

DATASETS = [
    ("Bulk cell lines", "#1e88e5"),
    ("Single-cell plexDIA", "#00897b"),
    ("Spatial DVP", "#8e24aa"),
    ("Phosphoproteomics", "#e8a000"),
    ("Astral / timsTOF / Orbitrap", "#6d4c41"),
]


def _stream(ax, p0, p1, color, lw=2.2, alpha=0.55):
    """A smooth cubic-Bezier flow line from p0 to p1."""
    dx = p1[0] - p0[0]
    c1 = (p0[0] + dx * 0.45, p0[1])
    c2 = (p0[0] + dx * 0.55, p1[1])
    path = Path([p0, c1, c2, p1],
                [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(PathPatch(path, fill=False, edgecolor=color, lw=lw,
                           alpha=alpha, zorder=2, capstyle="round"))


def build(out_name: str = "graphical_abstract", with_title: bool = True) -> FsPath:
    fig, ax = plt.subplots(figsize=(12.0, 5.4))
    ax.set_xlim(0, 122)
    ax.set_ylim(0, 54)
    ax.axis("off")

    # ---- title (standalone graphical abstract only; omitted for the in-text
    #      Figure 1a panel, where the paper title would be redundant) ----
    if with_title:
        ax.text(3, 50.4, "quantmsdiann", fontsize=21, fontweight="bold",
                color=QMBLUE, family=FONT)
        ax.text(3, 46.4, "one SDRF-driven DIA-NN workflow for archive-scale "
                "reanalysis of public DIA proteomics", fontsize=10.5,
                color=GREY, family=FONT)

    hub = (58.0, 30.0)
    hub_r = 8.5
    gate_x = 37.0

    # ---- diverse datasets converging through the SDRF gate ----
    ys = [40.5, 34.5, 28.5, 22.5, 16.5]
    for (lab, col), y in zip(DATASETS, ys):
        ax.add_patch(Circle((6.0, y), 1.15, facecolor=col, edgecolor="white",
                            lw=1.0, zorder=4)).set_path_effects(SHADOW)
        ax.text(8.6, y, lab, ha="left", va="center", fontsize=10.2,
                color=DARK, family=FONT)
        _stream(ax, (gate_x, y), (hub[0] - hub_r + 1, hub[1]), col)  # gate -> hub
        _stream(ax, (28.0, y), (gate_x, y), col, lw=2.0, alpha=0.5)  # label -> gate
    ax.text(20, 44.6, "Diverse public DIA datasets", ha="center",
            fontsize=11.5, fontweight="bold", color=DARK, family=FONT)

    # SDRF "gate": a soft vertical lens all streams pass through
    ax.add_patch(Ellipse((gate_x, 30), 6.5, 30, facecolor=QMGREEN_L,
                        edgecolor=QMGREEN, lw=1.6, zorder=3,
                        alpha=0.92)).set_path_effects(SHADOW)
    ax.text(gate_x, 30, "SDRF", ha="center", va="center",
            rotation=90, fontsize=14, fontweight="bold", color=QMGREEN,
            family=FONT, zorder=5)
    ax.text(gate_x, 11.5, "standardised\nmetadata", ha="center", va="center",
            fontsize=8.8, color=QMGREEN, family=FONT, style="italic")

    # ---- the engine hub, ringed by cluster nodes (scalability) ----
    for k in range(16):
        a = 2 * math.pi * k / 16
        ax.add_patch(Circle((hub[0] + (hub_r + 2.4) * math.cos(a),
                             hub[1] + (hub_r + 2.4) * math.sin(a)), 0.85,
                            facecolor=QMBLUE, edgecolor="white", lw=0.5,
                            alpha=0.55, zorder=3))
    ax.add_patch(Circle(hub, hub_r, facecolor=QMBLUE_L, edgecolor=QMBLUE,
                        lw=2.2, zorder=5)).set_path_effects(SHADOW)
    ax.text(hub[0], hub[1] + 2.3, "quantmsdiann", ha="center", va="center",
            fontsize=12.5, fontweight="bold", color=QMBLUE, family=FONT, zorder=6)
    ax.text(hub[0], hub[1] - 1.0, "nf-core", ha="center", va="center",
            fontsize=9.5, color=DARK, family=FONT, zorder=6)
    ax.text(hub[0], hub[1] - 3.4, "DIA-NN", ha="center", va="center",
            fontsize=9.5, color=DARK, family=FONT, zorder=6)
    ax.text(hub[0], 44.6, "One scalable workflow",
            ha="center", fontsize=11.5, fontweight="bold", color=QMBLUE,
            family=FONT)

    # ---- real scaling curve beneath the hub (no frame) ----
    bx0, bx1, by0, by1 = 47.0, 69.0, 6.5, 16.5
    hmin, hmax = min(SWEEP_HOURS), max(SWEEP_HOURS)
    px = lambda nd: bx0 + (nd - SWEEP_NODES[0]) / (SWEEP_NODES[-1] - SWEEP_NODES[0]) * (bx1 - bx0)
    py = lambda h: by0 + (h - hmin) / (hmax - hmin) * (by1 - by0)
    ax.plot([bx0, bx0], [by0, by1 + 0.5], color="#c0ccd6", lw=0.9, zorder=2)
    ax.plot([bx0, bx1], [by0, by0], color="#c0ccd6", lw=0.9, zorder=2)
    ax.plot([px(n) for n in SWEEP_NODES], [py(h) for h in SWEEP_HOURS],
            "-o", color=QMBLUE, lw=2.4, ms=5, mfc=QMBLUE, mec="white",
            mew=0.8, zorder=4)
    ax.text(px(10) + 0.6, py(37.7), "37.7 h", fontsize=8.8, color=QMBLUE,
            fontweight="bold", ha="left", va="bottom", family=FONT)
    ax.text(px(300), py(2.4) + 1.0, "2.4 h", fontsize=8.8, color=QMBLUE,
            fontweight="bold", ha="right", va="bottom", family=FONT)
    ax.text(bx0, by0 - 1.2, "10", fontsize=8, color=GREY, ha="center", va="top")
    ax.text(bx1, by0 - 1.2, "300 nodes", fontsize=8, color=GREY, ha="right", va="top")
    ax.text((bx0 + bx1) / 2, by1 + 1.4, "wall-clock vs cluster width "
            "(2,300 files)", ha="center", fontsize=8, color=GREY,
            family=FONT, style="italic")

    # ---- one harmonised output node ----
    _stream(ax, (hub[0] + hub_r - 1, hub[1]), (101.5, hub[1]), QMGREEN,
            lw=3.4, alpha=0.85)
    out = (107.0, hub[1])
    ax.add_patch(Circle(out, 7.0, facecolor=QMGREEN_L, edgecolor=QMGREEN,
                        lw=2.2, zorder=5)).set_path_effects(SHADOW)
    ax.text(out[0], out[1] + 1.6, "QPX", ha="center", va="center",
            fontsize=15, fontweight="bold", color=QMGREEN, family=FONT, zorder=6)
    ax.text(out[0], out[1] - 2.0, "harmonised\nmatrix", ha="center", va="center",
            fontsize=8.6, color=DARK, family=FONT, zorder=6)
    ax.text(out[0], 44.6, "Integrated output", ha="center",
            fontsize=11.5, fontweight="bold", color=QMGREEN, family=FONT)
    ax.text(out[0], out[1] - 9.5, "Parquet $\\cdot$ MSstats $\\cdot$ pmultiqc QC\n"
            "queryable across all modalities", ha="center", fontsize=8.6,
            color=GREY, family=FONT, style="italic")

    # ---- pillars ----
    ax.text(61, 1.8, "Metadata standards      $\\bullet$      Scalability"
            "      $\\bullet$      Integration", ha="center", fontsize=12,
            fontweight="bold", color=DARK, family=FONT)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUT_DIR / out_name
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return stem.with_suffix(".svg")


if __name__ == "__main__":
    # standalone graphical abstract (with paper title) ...
    print("wrote", build("graphical_abstract", with_title=True))
    # ... and the title-less variant used as main-text Figure 1a.
    print("wrote", build("overview_flow", with_title=False))
