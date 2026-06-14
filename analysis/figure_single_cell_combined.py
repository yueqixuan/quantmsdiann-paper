#!/usr/bin/env python
"""Fig 3 (combined) - single-cell reanalysis, DIA-NN 1.8.1 vs 2.5.1 Enterprise.

Merges the count overview and the mechanistic multi-view into one figure, with
dataset accessions on every panel. 2x3:
  A  Per-cell protein groups  (box + jitter, both HeLa datasets x version)
  B  Data-completeness curve  (HeLa Astral) -- the MBR headline
  C  plexDIA deposited vs quantms.io (MSV000093870) -- orthogonal reanalysis axis
  D  Relative rank-abundance  (HeLa Astral) -- depth, normalized to each run's top
  E  CV across cells          (HeLa Astral) -- quantitative precision
  F  legend

Data: data/single_cell/mv_{completeness,per_cell,rank_abundance,cv}.tsv.

Run:  python -m analysis.figure_single_cell_combined
Out:  analysis/figures/manuscript/fig3_single_cell_combined.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()

REPO = Path(__file__).resolve().parents[1]
D = REPO / "data" / "single_cell"
OUT = REPO / "analysis" / "figures" / "manuscript" / "fig3_single_cell_combined.svg"

VERS = ["1_8_1", "2_5_1_enterprise"]
VLAB = {"1_8_1": "1.8.1", "2_5_1_enterprise": "2.5.1 Enterprise"}
VCOL = {v: fs.VERSION_COLORS[v] for v in VERS}
ACC = {"HeLa Astral SC": "PXD046357", "HeLa One-Tip": "PXD044991"}
FLAG = "HeLa Astral SC"
FLAG_T = f"HeLa Astral SC\n({ACC[FLAG]})"
# plexDIA (Galatidou 2024) protein groups, deposited vs quantms.io reanalysis.
# Canonical report-based target-only counts (channel-confident), from
# analysis/figures/plexDIA/MSV000093870/comparison_counts.tsv.
PLEX_DEP, PLEX_QM = 2122, 2904


def _completeness(ax):
    df = pd.read_csv(D / "mv_completeness.tsv", sep="\t")
    df = df[df["dataset"] == FLAG]
    for v in VERS:
        s = df[df["version"] == v].sort_values("min_cells")
        ax.plot(s["min_cells"], s["n_proteins"], "-o", ms=4, color=VCOL[v], label=VLAB[v])
    ax.set_xlabel("quantified in ≥ N cells"); ax.set_ylabel("protein groups")
    ax.set_title(f"Data completeness — {FLAG_T}")
    fs.kfmt_axis(ax.yaxis); fs.despine(ax)


def _percell(ax):
    df = pd.read_csv(D / "mv_per_cell.tsv", sep="\t")
    datasets = list(dict.fromkeys(df["dataset"]))
    bw = 0.36; rng = np.random.default_rng(0)
    for gi, ds in enumerate(datasets):
        for k, v in enumerate(VERS):
            vals = df[(df["dataset"] == ds) & (df["version"] == v)]["pg_count"].values
            x = gi + (k - 0.5) * bw
            bp = ax.boxplot([vals], positions=[x], widths=bw * 0.85, patch_artist=True, showfliers=False)
            fs.style_boxplot(bp, color=VCOL[v])
            ax.scatter(x + rng.uniform(-0.07, 0.07, len(vals)), vals, s=12, color=VCOL[v], alpha=0.7, edgecolors="none", zorder=3)
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels([f"{d}\n({ACC[d]})" for d in datasets])
    ax.set_ylabel("protein groups / cell"); ax.set_title("Per-cell protein groups")
    fs.kfmt_axis(ax.yaxis); fs.despine(ax)


def _plexdia(ax):
    vals = [PLEX_DEP, PLEX_QM]; cols = [fs.COMPARISON["original"], fs.COMPARISON["quantmsdiann"]]
    bw = 0.38; pos = [-bw/2, bw/2]
    b = ax.bar(pos, vals, bw, color=cols, edgecolor="white", linewidth=0.6)
    pct = round(100*(vals[1]-vals[0])/vals[0])
    ax.annotate(f"+{pct}%", (b[1].get_x()+b[1].get_width()/2, vals[1]), textcoords="offset points",
                xytext=(0, 2), ha="center", va="bottom", fontsize=8, fontweight="bold", color=fs.COMPARISON["quantmsdiann"])
    ax.set_xticks(pos); ax.set_xticklabels(["Galatidou 2024\n(deposited)\n1.8.1b16", "quantms.io\n2.5.0"])
    ax.set_xlim(-0.6, 0.6); ax.set_ylabel("protein groups"); ax.set_title("plexDIA — MSV000093870")
    fs.kfmt_axis(ax.yaxis); fs.despine(ax)


def _rank(ax):
    df = pd.read_csv(D / "mv_rank_abundance.tsv", sep="\t")
    for v in VERS:
        s = df[df["version"] == v].sort_values("rank")
        y = s["log10_intensity"].values - s["log10_intensity"].max()  # relative to top protein
        ax.plot(s["rank"], y, "-", color=VCOL[v], label=VLAB[v])
    ax.set_xlabel("protein group rank"); ax.set_ylabel("log10 intensity (rel. to top)")
    ax.set_title(f"Dynamic range — {FLAG_T}")
    fs.despine(ax)


def _cv(ax):
    df = pd.read_csv(D / "mv_cv.tsv", sep="\t")
    bins = np.linspace(0, 1.0, 41)
    for v in VERS:
        cv = df[df["version"] == v]["cv"].clip(0, 1.0)
        ax.hist(cv, bins=bins, density=True, histtype="step", linewidth=1.6, color=VCOL[v], label=VLAB[v])
        ax.axvline(np.median(df[df["version"] == v]["cv"]), color=VCOL[v], linestyle=":", linewidth=1)
    ax.set_xlabel("CV across cells"); ax.set_ylabel("density")
    ax.set_title(f"Quantitative precision — {FLAG_T}")
    fs.despine(ax)


# Total counts (report-based, target-only): precursors and protein groups,
# 1.8.1 -> 2.5.1 Enterprise, for the two HeLa single-cell datasets.
TOTALS = {
    "Astral":  {"precursors": (19674, 24111), "proteins": (3903, 4574)},
    "One-Tip": {"precursors": (11618, 16534), "proteins": (1597, 2306)},
}


def _totals(ax):
    """Dual-axis totals: precursors (left axis) and protein groups (right axis),
    side by side, so the two very different scales are both readable."""
    bw = 0.36
    ax2 = ax.twinx()
    dsx = list(TOTALS)                       # [Astral, One-Tip]
    prec_x = {d: i for i, d in enumerate(dsx)}            # 0, 1
    prot_x = {d: i + len(dsx) + 0.6 for i, d in enumerate(dsx)}  # 2.6, 3.6
    def draw(target, xmap, metric):
        for k, v in enumerate(VERS):
            idx = 0 if v == "1_8_1" else 1
            for d in dsx:
                x = xmap[d] + (k - 0.5) * bw
                hi = TOTALS[d][metric][idx]
                target.bar(x, hi, bw, color=VCOL[v], edgecolor="white", linewidth=0.6)
                if v != "1_8_1":
                    lo = TOTALS[d][metric][0]
                    target.annotate(f"+{round(100*(hi-lo)/lo)}%", (x, hi), textcoords="offset points",
                                    xytext=(0, 2), ha="center", va="bottom", fontsize=7, fontweight="bold", color=VCOL[v])
    draw(ax, prec_x, "precursors")
    draw(ax2, prot_x, "proteins")
    ax.axvline(len(dsx) - 0.2, color="#cccccc", linewidth=0.8)
    ticks = [prec_x[d] for d in dsx] + [prot_x[d] for d in dsx]
    ax.set_xticks(ticks)
    ax.set_xticklabels([d.replace("HeLa ", "").replace(" SC", "") for d in dsx] * 2, fontsize=8)
    ax.set_xlim(-0.7, prot_x[dsx[-1]] + 0.7)
    ax.set_ylabel("precursors"); ax2.set_ylabel("protein groups")
    ax.set_ylim(0, 27000); ax2.set_ylim(0, 5400)
    fs.kfmt_axis(ax.yaxis); fs.kfmt_axis(ax2.yaxis)
    ax.set_title("Total identifications")
    ax.text(np.mean([prec_x[d] for d in dsx]), -0.16, "precursors", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=8.5, fontweight="bold")
    ax.text(np.mean([prot_x[d] for d in dsx]), -0.16, "protein groups", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=8.5, fontweight="bold")
    for sp in ("top",): ax.spines[sp].set_visible(False); ax2.spines[sp].set_visible(False)


def render(out: Path) -> Path:
    fig, ax = plt.subplots(2, 3, figsize=(13.5, 7.8))
    _totals(ax[0, 0]); _percell(ax[0, 1]); _completeness(ax[0, 2])
    _rank(ax[1, 0]); _cv(ax[1, 1]); _plexdia(ax[1, 2])
    handles = [Line2D([0], [0], color=VCOL[v], marker="o", linewidth=2, markersize=7, label=f"DIA-NN {VLAB[v]}") for v in VERS]
    handles += [Line2D([0], [0], marker="s", linestyle="none", markersize=9, markerfacecolor=fs.COMPARISON["original"],
                markeredgecolor="white", label="plexDIA deposited (Galatidou 1.8.1b16)"),
                Line2D([0], [0], marker="s", linestyle="none", markersize=9, markerfacecolor=fs.COMPARISON["quantmsdiann"],
                markeredgecolor="white", label="plexDIA quantms.io (2.5.0)")]
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.015), fontsize=9.5)
    for a, lab in zip([ax[0,0],ax[0,1],ax[0,2],ax[1,0],ax[1,1],ax[1,2]], "ABCDEF"):
        a.text(-0.14, 1.06, lab, transform=a.transAxes, fontsize=14, fontweight="bold", va="bottom", ha="right")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out); plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
