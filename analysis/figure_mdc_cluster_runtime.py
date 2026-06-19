#!/usr/bin/env python
"""Per-step wall-clock of one quantmsdiann run on an independent (non-EBI) cluster.

Companion to the EBI per-step figure (`runtime_per_step.svg`): same horizontal
per-step box-plot style, rendered by the *same* function
(`figure_performance_trace.render_per_step_boxplot`) so panels (a) EBI and
(b) non-EBI in Supplementary Fig. S1 are directly comparable. Shows that
quantmsdiann runs unchanged on a third-party SLURM cluster — a 279-raw-file
cohort processed by a different group at the Max Delbruck Center (MDC Berlin)
under the SLURM+Singularity profile (DIA-NN 2.5.0, quantmsdiann v2.1.0; total
wall-clock 7h49m). The per-file stages run in parallel across the cluster, so
the summed task time collapses to a few hours of wall-clock; the run-level
facts live in the figure caption, matching the title-free house style of the
EBI panel.

Privacy: built ONLY from an anonymised per-step duration table (step name +
seconds + memory + %CPU); raw file names, sample identifiers, internal paths
and contact details from the source run are never imported or stored.

Source: analysis/figures/performance/data/mdc_step_durations.tsv
Out:    analysis/figures/performance/mdc_cluster_runtime.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd

from analysis.figure_performance_trace import render_per_step_boxplot

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "analysis" / "figures" / "performance" / "data" / "mdc_step_durations.tsv"
EBI_TSV = REPO / "analysis" / "figures" / "performance" / "data" / "runtime_per_step.tsv"
OUT = REPO / "analysis" / "figures" / "performance" / "mdc_cluster_runtime.svg"


def _ebi_panel_height() -> float:
    """Figure height of the EBI per-step panel (S1a) so this companion panel
    (S1b) renders at the same height when laid out side by side. Mirrors the
    auto-height rule in render_per_step_boxplot, keyed on the EBI step count."""
    n_ebi = len(pd.read_csv(EBI_TSV, sep="\t"))
    return max(3.5, 0.45 * n_ebi + 1.5)


def render(out: Path) -> Path:
    df = pd.read_csv(DATA, sep="\t")
    df = df[df["status"] == "COMPLETED"]
    # Per-step duration distributions, in seconds (per-file stages carry the
    # full n=279 distribution; collective stages contribute a single point).
    durations = {step: sub["duration_s"].tolist()
                 for step, sub in df.groupby("step")}
    # Same ordering convention as the EBI panel: descending median (slowest at top).
    order = sorted(durations, key=lambda s: pd.Series(durations[s]).median(),
                   reverse=True)
    summary = pd.DataFrame({"step": order})
    render_per_step_boxplot(durations, summary, out, fig_h=_ebi_panel_height())
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
