#!/usr/bin/env python
"""Fig 2 - quantmsdiann validation row (scaling + ProteoBench accuracy).

One row of three panels (the former Fig 1 b/c/d, now a standalone figure so the
workflow gets Fig 1 to itself):
  (a) wall-clock versus cluster nodes (PXD071075 single-cell sweep)
  (b) wall-clock to finish each reanalysis, one bar per dataset
  (c) ProteoBench quantification-accuracy concordance vs standalone DIA-NN

Reuses the existing per-panel renderers (composite/ax mode) so the numbers stay
identical to the standalone figures.

Out: analysis/figures/manuscript/fig2_validation.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()
from analysis.figure_queue_size_sweep import render_queue_size_sweep
from analysis.figure_performance_trace import render_parallelism_scatter
from analysis import figure_proteobench_accuracy as acc

REPO = Path(__file__).resolve().parents[1]
PERF = REPO / "analysis" / "figures" / "performance" / "data"
OUT = REPO / "analysis" / "figures" / "manuscript" / "fig2_validation.svg"


def render(out: Path) -> Path:
    from matplotlib.patches import Patch
    from analysis.figure_performance_trace import INSTRUMENT_COLOURS

    dq = pd.read_csv(PERF / "queue_size_sweep.tsv", sep="\t")
    dp = pd.read_csv(PERF / "parallelism_data.tsv", sep="\t")
    # Three EQUAL-width panels that each fill their cell (no forced square),
    # so a/b/c are the same size with no whitespace gaps.
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 4.6))
    render_queue_size_sweep(dq, ax=ax[0], composite=True)
    render_parallelism_scatter(dp, ax=ax[1], composite=True, show_legend=False,
                               short_labels=True)
    # panel c shares panel b's instrument colours (read from the bottom legend);
    # shape encodes the quantms-diann version, so c needs no colour legend.
    from analysis.figure_id_vs_epsilon import _COMMUNITY_COMPARATOR_DATASETS
    ds_colors = {}
    for ds in _COMMUNITY_COMPARATOR_DATASETS:
        inst = dp.loc[dp["dataset"] == ds, "instrument"]
        if len(inst):
            ds_colors[ds] = INSTRUMENT_COLOURS.get(inst.iloc[0], "#9e9e9e")
    acc.draw_strip(ax[2], compact=True, dataset_colors=ds_colors)
    for a, lab in zip(ax, "abc"):
        a.text(-0.06, 1.05, f"({lab})", transform=a.transAxes, fontsize=14,
               fontweight="bold", va="bottom", ha="left")
    # instrument legend (for panel b) as a figure-level strip at the bottom,
    # so it does not make panel b taller than a/c.
    insts = [i for i in dict.fromkeys(dp["instrument"]) if isinstance(i, str)]
    handles = [Patch(facecolor=INSTRUMENT_COLOURS.get(i, "#9e9e9e"), edgecolor="#222222", label=i)
               for i in insts]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7,
               frameon=False, title="Instrument / dataset (panels b, c)", title_fontsize=7.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.13, 1, 1), w_pad=0.4)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
