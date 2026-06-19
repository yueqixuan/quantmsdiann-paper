#!/usr/bin/env python
"""Reanalysis-improvement figure: original analysis vs quantmsdiann reanalysis.

Two panels of paired horizontal bars, one row per public DIA deposit:
  (a) Protein groups  — all seven deposits.
  (b) Precursors      — only the deposits whose ORIGINAL analysis was DIA-NN
      (HeLa Astral, One-Tip, spatial DVP); "precursor" is a DIA-NN concept, so
      OpenSWATH/PCT-SWATH/protein-only originals have no comparable count.

Original (grey) = deposit matrix / published headline; reanalysis (coloured by
the DIA-NN version used) = counted from the precursor report (parquet), per the
counting rule Vadim specified. Single-cell PG keeps the pg_matrix union (so the
match-between-runs gain stays visible); NCI-60/ProCan use >=2-peptide; Sun its
consistency filter; plexDIA channel-confident mTRAQ. Provenance per row is in
the data TSV.

Source: analysis/figures/reanalysis/data/reanalysis_improvement.tsv
Out:    analysis/figures/manuscript/fig_reanalysis_improvement.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "analysis" / "figures" / "reanalysis" / "data" / "reanalysis_improvement.tsv"
OUT = REPO / "analysis" / "figures" / "manuscript" / "fig_reanalysis_improvement.svg"

ORIG_COLOUR = "#9e9e9e"
# DIA-NN version palette: 2.5.0 -> 2.5.1 light -> dark blue; enterprise amber.
VERSION_COLOUR = {
    "2.5.0": "#90caf9",
    "2.5.1": "#1976d2",
    "2.5.1-enterprise": "#ff8f00",
}


def _panel(ax, df, orig_col, new_col, xlabel):
    """Paired horizontal bars (grey original / version-coloured reanalysis),
    sorted by gain (largest at top). `df` already filtered to rows with data."""
    df = df.copy()
    df["gain"] = df[new_col] / df[orig_col] - 1.0
    df = df.sort_values("gain", ascending=True).reset_index(drop=True)
    xmax = df[new_col].max()
    bar_h = 0.36
    for i, row in df.iterrows():
        vcol = VERSION_COLOUR.get(row["diann_version"], "#1976d2")
        ax.barh(i + bar_h / 2 + 0.02, row[orig_col], height=bar_h,
                color=ORIG_COLOUR, edgecolor="#37474f", linewidth=0.6, zorder=2)
        ax.barh(i - bar_h / 2 - 0.02, row[new_col], height=bar_h,
                color=vcol, edgecolor="#37474f", linewidth=0.6, zorder=2)
        ax.text(row[orig_col] + xmax * 0.008, i + bar_h / 2 + 0.02,
                f"{int(row[orig_col]):,}", va="center", ha="left", fontsize=8, color="#555555")
        ax.text(row[new_col] + xmax * 0.008, i - bar_h / 2 - 0.02,
                f"{int(row[new_col]):,}", va="center", ha="left", fontsize=8,
                fontweight="bold", color="#222222")
        ax.text(xmax * 1.17, i, f"+{row['gain'] * 100:.0f}%",
                va="center", ha="right", fontsize=10, fontweight="bold", color=vcol)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{r['label']}\n{r['original_engine']} → DIA-NN {r['diann_version']}"
                        for _, r in df.iterrows()], fontsize=8.5)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_xlim(0, xmax * 1.20)
    ax.set_ylim(-0.7, len(df) - 0.3)
    fs.despine(ax)
    ax.tick_params(axis="y", length=0)
    return df


def render(out: Path) -> Path:
    df = pd.read_csv(DATA, sep="\t")
    pg = df.dropna(subset=["reanalysis"])
    prec = df.dropna(subset=["orig_precursors", "new_precursors"]).copy()
    prec["orig_precursors"] = prec["orig_precursors"].astype(int)
    prec["new_precursors"] = prec["new_precursors"].astype(int)

    n_pg, n_pr = len(pg), len(prec)
    fig = plt.figure(figsize=(8.4, 0.66 * (n_pg + n_pr) + 2.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[n_pg, n_pr], hspace=0.32)
    ax_pg = fig.add_subplot(gs[0])
    ax_pr = fig.add_subplot(gs[1])

    _panel(ax_pg, pg, "original", "reanalysis", "Protein groups")
    _panel(ax_pr, prec, "orig_precursors", "new_precursors", "Precursors")
    ax_pg.set_title("(a) Protein groups — all reanalysed deposits",
                    loc="left", fontsize=11, fontweight="bold")
    ax_pr.set_title("(b) Precursors — deposits with a DIA-NN original",
                    loc="left", fontsize=11, fontweight="bold")

    handles = [Patch(facecolor=ORIG_COLOUR, edgecolor="#37474f", label="Original analysis")]
    for v, c in VERSION_COLOUR.items():
        if (df["diann_version"] == v).any():
            handles.append(Patch(facecolor=c, edgecolor="#37474f",
                                 label=f"quantmsdiann (DIA-NN {v})"))
    ax_pg.legend(handles=handles, loc="lower right", frameon=False, fontsize=8,
                 title="reanalysis", title_fontsize=8.5)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
