"""F2d — Nextflow queueSize scaling.

Per-run sweep on a single cohort (PXD071075 single-cell early
results, Wang 2025). The x-axis is **Nextflow `executor.queueSize`**
— the maximum number of concurrent Nextflow tasks. Each task is one
SLURM job (single-core process in this configuration), so queueSize
also bounds the total CPU cores in use; `run_metadata.json` in each
sweep dir reports `queue_size == sweep_cores`.

5 sweep points: 10 / 50 / 100 / 200 / 300. The earlier 4-point sweep
(10 / 20 / 100 / 200) had a q20 outlier from cluster contention; the
re-run replaces q20 with q50 and adds q300 for the high end.

When run with no data present, `main()` exits cleanly with a
"no input yet" message; the suite never fails just because the
sweep hasn't been collected yet.

Spec: docs/superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis.figure_performance_trace import (
    load_trace,
    peak_concurrent_tasks,
    trace_wallclock_seconds,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
SWEEP_DIR = REPO_ROOT / "data" / "queue_size_sweep"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "performance"

# Sweep points staged on PRIDE for the PXD071075 single-cell sweep
# (Wang 2025): 10, 50, 100, 200, 300. Each sweep point's `run_metadata.json`
# confirms `queue_size == sweep_cores` — each Nextflow task is allocated
# one SLURM job, so the sweep varies max concurrent tasks (and total
# cores in lockstep). The earlier 4-point sweep (10/20/100/200) had a
# q20 outlier from cluster contention; this re-run replaces q20 with q50
# and adds q300 for the high end.
#
# Spec: docs/superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md
DEFAULT_QUEUE_SIZES = (10, 50, 100, 200, 300)


def iter_sweep_traces(
    queue_sizes: Iterable[int] = DEFAULT_QUEUE_SIZES,
) -> Iterable[tuple[int, Path]]:
    """Yield (queue_size, trace_path) for every sweep point whose
    `nextflow_trace.txt` exists locally. Looks for both `q<N>/` and
    zero-padded `q<NNN>/` directory layouts (the PRIDE deposit uses
    zero-padded `v2_5_0_sweep_010cores` etc.; the staged local layout
    follows that convention). Points without a trace are silently
    skipped — `main()` decides whether the remaining set is sufficient
    to render the figure."""
    for q in queue_sizes:
        candidates = (
            SWEEP_DIR / f"q{q:03d}" / "nextflow_trace.txt",
            SWEEP_DIR / f"q{q}" / "nextflow_trace.txt",
        )
        for path in candidates:
            if path.exists():
                yield q, path
                break


def collect_sweep_rows(
    queue_sizes: Iterable[int] = DEFAULT_QUEUE_SIZES,
) -> pd.DataFrame:
    """Per-queueSize row: `(queue_size, wallclock_s, peak_concurrent,
    n_tasks)`. Returns an empty DataFrame with that schema when no
    sweep traces are present.

    The wallclock is computed from
    `max(submit+duration) − min(submit)` across all rows in the trace
    (FAILED retries included — they did occupy slots), matching the
    existing benchmark trace-wallclock semantics."""
    rows: list[dict] = []
    for q, path in iter_sweep_traces(queue_sizes):
        df = load_trace(path)
        wallclock_s = trace_wallclock_seconds(df)
        peak, _med = peak_concurrent_tasks(df)
        rows.append({
            "queue_size": q,
            "wallclock_s": float(wallclock_s),
            "peak_concurrent": int(peak),
            "n_tasks": int(len(df)),
        })
    return pd.DataFrame(
        rows,
        columns=["queue_size", "wallclock_s", "peak_concurrent", "n_tasks"],
    ).sort_values("queue_size").reset_index(drop=True)


def render_queue_size_sweep(
    df: pd.DataFrame,
    svg_path: Path | None = None,
    *,
    ax: plt.Axes | None = None,
    composite: bool = False,
) -> None:
    """Wallclock vs cluster-node count for the PXD071075 single-cell
    queueSize sweep. Log-log axes with explicit decimal tick labels
    (10/50/100/200/300 — no 10^1 / 10^2 power notation). Dot fill is
    proportional to queueSize so the eye reads the sweep direction at
    a glance. The wallclock value (hours) is printed next to each dot
    for direct readability.

    When `df` has fewer than 2 points the function still emits an SVG
    but with a "needs more sweep points" placeholder.

    Pass `ax` to draw into an existing axes (composite figures); omit
    `svg_path` in that mode."""
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(5.6, 4.4))
    if len(df) < 2:
        ax.text(
            0.5, 0.5,
            f"sweep is data-bound on experiment #12 — "
            f"need at least 2 sweep points, have {len(df)}",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color="#888888",
        )
        ax.set_axis_off()
        if own_fig:
            fig.tight_layout()
            assert svg_path is not None
            svg_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(svg_path, bbox_inches="tight")
            plt.close(fig)
        return

    df = df.sort_values("queue_size").reset_index(drop=True)
    xs = df["queue_size"].astype(float).values
    ys = (df["wallclock_s"] / 3600.0).values

    # Uniform dot colour — single brand-blue used for the line and
    # all markers. The sweep direction is already legible from the x
    # ordering plus the per-dot wallclock annotation; a colour
    # gradient added visual noise without information.
    dot_color = "#1976d2"
    dot_size = 90 if composite else 140
    line_width = 1.4 if composite else 1.8
    ann_size = 7.0 if composite else 8.5
    label_size = 9 if composite else 11
    tick_size = 8 if composite else 10
    ax.plot(
        xs, ys,
        color=dot_color, linewidth=line_width, alpha=0.55, zorder=2,
    )
    ax.scatter(
        xs, ys,
        s=dot_size, c=dot_color,
        edgecolors="#0d47a1", linewidths=0.8,
        zorder=3,
    )

    # Print the wallclock value next to each dot. Right-side anchor
    # so the labels sit just outside the markers without colliding.
    for xi, yi in zip(xs, ys):
        ax.annotate(
            f"{yi:.1f} h",
            xy=(xi, yi),
            xytext=(6, 6), textcoords="offset points",
            fontsize=ann_size, color="#1a237e", fontweight="bold",
            ha="left", va="bottom",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    # Explicit integer tick labels — no 10^1 / 10^2 power notation.
    from matplotlib.ticker import FixedLocator, FixedFormatter
    ax.xaxis.set_major_locator(FixedLocator(xs))
    ax.xaxis.set_major_formatter(FixedFormatter([str(int(x)) for x in xs]))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    # Y-axis: pick a handful of round-number ticks in the data range.
    y_ticks = []
    for cand in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0):
        if ys.min() / 1.4 <= cand <= ys.max() * 1.4:
            y_ticks.append(cand)
    if y_ticks:
        ax.yaxis.set_major_locator(FixedLocator(y_ticks))
        ax.yaxis.set_major_formatter(FixedFormatter(
            [f"{t:g}" for t in y_ticks]
        ))
    ax.yaxis.set_minor_locator(FixedLocator([]))

    ax.set_xlabel("Cluster nodes", fontsize=label_size)
    ax.set_ylabel("Workflow wallclock (hours)", fontsize=label_size)
    # No gridlines — keep the panel clean; the per-dot wallclock
    # annotations make every value readable without grid reference.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=tick_size)

    # Pad axis limits so labels at the edges don't get clipped.
    ax.set_xlim(xs.min() / 1.4, xs.max() * 1.8)
    ax.set_ylim(ys.min() / 1.6, ys.max() * 1.6)

    if own_fig:
        fig.tight_layout()
        assert svg_path is not None
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)


def write_sweep_tsv(df: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep="\t", index=False)


def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = FIGURES_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    df = collect_sweep_rows()
    if df.empty:
        print(
            "F2d (queue-size sweep) — no trace data under "
            f"{SWEEP_DIR}/q<N>/nextflow_trace.txt yet; skipping render. "
            "See docs/superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md."
        )
        return 0
    write_sweep_tsv(df, data_dir / "queue_size_sweep.tsv")
    render_queue_size_sweep(df, FIGURES_DIR / "queue_size_sweep.svg")
    print(
        f"F2d (queue-size sweep) rendered from {len(df)} sweep point(s) "
        f"({list(df['queue_size'])})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
