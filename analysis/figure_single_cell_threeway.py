#!/usr/bin/env python
"""Supplementary: single-cell three-way DIA-NN version comparison.

1.8.1 vs 2.5.1 (standard, free) vs 2.5.1 Enterprise, for the human HeLa
single-cell datasets. Transparency companion to the main Fig 3 (which shows
1.8.1 vs Enterprise only): shows that most of the gain comes from 1.8 -> 2.5,
with Enterprise adding a marginal further improvement.

Only the two fair, uniform-threshold metrics are shown (both at PG.Q <= 1%):
  A  avg protein groups per cell
  B  complete-profile protein groups (quantified in every cell)
Raw precursor counts are intentionally omitted because the pipeline's per-version
precursor q-value (0.01 for 1.8.1, 0.05 for 2.5.x) makes them non-comparable.

Run:  python -m analysis.figure_single_cell_threeway
Out:  analysis/figures/manuscript/supplementary/supp_single_cell_threeway.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()

REPO_ROOT = Path(__file__).resolve().parents[1]
COUNTS = REPO_ROOT / "data" / "single_cell" / "report_counts_threeway.tsv"
OUT = REPO_ROOT / "analysis" / "figures" / "manuscript" / "supplementary" / "supp_single_cell_threeway.svg"

VERSIONS = ["1_8_1", "2_5_1", "2_5_1_enterprise"]
VLABEL = {"1_8_1": "1.8.1", "2_5_1": "2.5.1 (free)", "2_5_1_enterprise": "2.5.1 Enterprise"}
VCOL = {v: fs.VERSION_COLORS[v] for v in VERSIONS}


def _val(df, ds, ver, col):
    s = df[(df["dataset"] == ds) & (df["version"] == ver)]
    return int(s[col].iloc[0]) if len(s) else 0


def _panel(ax, df, datasets, col, title, ylabel):
    bw = 0.26
    xs = range(len(datasets))
    for k, ver in enumerate(VERSIONS):
        vals = [_val(df, d, ver, col) for d in datasets]
        off = (k - 1) * bw
        ax.bar([x + off for x in xs], vals, bw, color=VCOL[ver],
               edgecolor="white", linewidth=0.6, label=VLABEL[ver])
    ax.set_xticks(list(xs))
    ax.set_xticklabels(datasets)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fs.kfmt_axis(ax.yaxis)
    fs.despine(ax)


def render(out: Path) -> Path:
    df = pd.read_csv(COUNTS, sep="\t")
    datasets = list(dict.fromkeys(df["dataset"]))
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.3))
    _panel(axes[0], df, datasets, "avg_pg_per_run", "Protein groups per cell", "PG / cell (mean)")
    _panel(axes[1], df, datasets, "complete_profile", "Complete-profile protein groups", "PG in every cell")
    handles = [Line2D([0], [0], marker="s", linestyle="none", markersize=9,
               markerfacecolor=VCOL[v], markeredgecolor="white",
               label=f"DIA-NN {VLABEL[v]}") for v in VERSIONS]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
