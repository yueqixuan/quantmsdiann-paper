#!/usr/bin/env python
"""Per-step wall-clock of one quantmsdiann run on an independent (non-EBI) cluster.

Demonstrates that quantmsdiann runs unchanged on a third-party SLURM cluster:
a 279-raw-file phosphoproteomics cohort processed by a different group at the
Max Delbruck Center (MDC Berlin) under the SLURM+Singularity profile (DIA-NN
2.5.0, quantmsdiann v2.1.0). Total wall-clock 7h49m; the per-file stages run in
parallel across the cluster, so the summed task time (hundreds of CPU-hours)
collapses to a few hours of wall-clock.

Privacy: this figure is built ONLY from an anonymised per-step duration table
(step name + seconds + memory + %CPU); raw file names, sample identifiers,
internal paths and contact details from the source run are never imported or
stored.

Source: analysis/figures/performance/data/mdc_step_durations.tsv
Out:    analysis/figures/performance/mdc_cluster_runtime.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "analysis" / "figures" / "performance" / "data" / "mdc_step_durations.tsv"
OUT = REPO / "analysis" / "figures" / "performance" / "mdc_cluster_runtime.svg"

# Non-sensitive run facts (public software versions + co-author institution).
N_FILES = 279
TOTAL_WALLCLOCK = "7 h 49 m"
CLUSTER = "Max Delbrück Center (MDC Berlin)"

# Pipeline order, and which stages run once-per-file (parallel) vs once-per-cohort.
STEP_ORDER = [
    "SAMPLESHEET_CHECK", "SDRF_PARSING", "INSILICO_LIBRARY_GENERATION",
    "PRELIMINARY_ANALYSIS", "ASSEMBLE_EMPIRICAL_LIBRARY", "INDIVIDUAL_ANALYSIS",
    "FINAL_QUANTIFICATION", "DIANN_MSSTATS", "SUMMARY_PIPELINE",
]
PER_FILE = {"PRELIMINARY_ANALYSIS", "INDIVIDUAL_ANALYSIS"}
LABEL = {
    "SAMPLESHEET_CHECK": "SDRF check",
    "SDRF_PARSING": "SDRF parsing",
    "INSILICO_LIBRARY_GENERATION": "In-silico library",
    "PRELIMINARY_ANALYSIS": "Preliminary analysis",
    "ASSEMBLE_EMPIRICAL_LIBRARY": "Empirical library",
    "INDIVIDUAL_ANALYSIS": "Individual analysis",
    "FINAL_QUANTIFICATION": "Final quantification",
    "DIANN_MSSTATS": "MSstats conversion",
    "SUMMARY_PIPELINE": "pmultiqc report",
}
PER_FILE_COL = "#d62728"   # red = per-file parallel stage (matches Fig 1 subway)
COLLECTIVE_COL = "#2e7d32"  # green = collective stage


def render(out: Path) -> Path:
    df = pd.read_csv(DATA, sep="\t")
    df = df[df["status"] == "COMPLETED"]
    df["minutes"] = df["duration_s"] / 60.0
    steps = [s for s in STEP_ORDER if s in set(df["step"])]
    rng = np.random.default_rng(0)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for y, step in enumerate(steps):
        vals = df.loc[df["step"] == step, "minutes"].values
        col = PER_FILE_COL if step in PER_FILE else COLLECTIVE_COL
        if step in PER_FILE and len(vals) > 1:
            bp = ax.boxplot([vals], positions=[y], vert=False, widths=0.55,
                            patch_artist=True, showfliers=False)
            for b in bp["boxes"]:
                b.set(facecolor=col, alpha=0.35, edgecolor=col, linewidth=1.0)
            for w in bp["whiskers"] + bp["caps"]:
                w.set(color=col, linewidth=1.0)
            for m in bp["medians"]:
                m.set(color=col, linewidth=1.4)
            ax.scatter(vals, y + rng.uniform(-0.16, 0.16, len(vals)), s=7,
                       color=col, alpha=0.5, edgecolors="none", zorder=3)
            ax.text(max(vals) * 1.02, y, f"n={len(vals)}", va="center", ha="left",
                    fontsize=7.5, color="#555555")
        else:
            v = float(vals[0]) if len(vals) else 0.0
            ax.barh(y, v, height=0.55, color=col, alpha=0.85, edgecolor=col)
            lab = f"{v * 60:.0f} s" if v < 1 else f"{v:.0f} min"
            ax.text(v * 1.08, y, lab, va="center", ha="left",
                    fontsize=7.5, color="#333333")

    ax.set_yticks(range(len(steps)))
    ax.set_yticklabels([LABEL.get(s, s) for s in steps], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Per-task wall-clock (minutes)", fontsize=10)
    ax.set_xscale("log")
    ax.set_title(f"quantmsdiann on an independent (non-EBI) SLURM cluster "
                 f"({N_FILES} raw files, total {TOTAL_WALLCLOCK})", fontsize=10.5)
    fs.despine(ax)
    handles = [Patch(facecolor=PER_FILE_COL, alpha=0.5, label="Per-file (parallel) stage"),
               Patch(facecolor=COLLECTIVE_COL, alpha=0.85, label="Collective (once per cohort) stage")]
    ax.legend(handles=handles, loc="center right", frameon=False, fontsize=8,
              bbox_to_anchor=(1.0, 0.42))
    ax.text(0.0, 1.13, f"{CLUSTER}; SLURM + Singularity; DIA-NN 2.5.0; quantmsdiann v2.1.0",
            transform=ax.transAxes, fontsize=8, color="#666666", va="bottom")
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
