"""Performance plots derived from per-task nextflow_trace.txt files.

Two paper-ready panels are produced from the same set of 20
quantmsdiann x ProteoBench benchmark runs (4 datasets x 5 DIA-NN versions);
the 3 cell-line datasets are excluded because their nextflow_trace.txt files
are not published on PRIDE.

Plot 1 (`parallelism_vs_wallclock`)
-----------------------------------

One dot per analysis. X-axis = peak concurrent task count observed during
the workflow (max overlap of [submit, submit+duration] intervals across all
trace rows). Y-axis = total workflow wallclock in hours, defined as
`max(submit + duration) - min(submit)` across the same trace. Colour =
instrument family (re-uses `INSTRUMENT_COLOURS` from
`analysis.figure_performance_runtime`). MS run count is constant (6) across
all benchmarks so the dot size is held at a single readable value and the
size dimension is intentionally not used.

Plot 2 (`runtime_per_step`)
---------------------------

Horizontal box plot of per-task `duration` in seconds, one row per workflow
step, ordered by descending median. The step name is the last colon-
separated segment of `name` before the parenthesised input argument
(e.g. `SAMPLESHEET_CHECK`, `DIANN_RUN`-equivalents like `PRELIMINARY_ANALYSIS`
and `INDIVIDUAL_ANALYSIS`). All 20 traces are aggregated. FAILED tasks
(Nextflow retry attempts) are excluded from the per-step duration
distributions; they remain in the parallelism analysis because they did
occupy cluster slots. The x-axis is log scale when the observed range spans
more than two orders of magnitude.

Truncated traces
----------------

PRIDE publishes nextflow_trace.txt files for all 20 analyses but several
copies are truncated (PXD062685 x5 versions ship the first two rows only,
PXD070049/v2_3_2 ships header-only). For Plot 1 we still emit a TSV row for
every analysis but mark `complete = False` for those with <= 2 data rows and
omit them from the scatter (legend annotation explains the count). For
Plot 2 we ingest whatever rows exist from every trace (the partial PXD062685
rows contribute valid SAMPLESHEET_CHECK / SDRF_PARSING durations).
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis.figure_original_vs_quantmsdiann import download_if_missing
from analysis.figure_performance_runtime import (
    BENCHMARK_DATASETS,
    BENCHMARK_INSTRUMENT,
    BENCHMARK_BASE,
    DIANN_VERSIONS,
    INSTRUMENT_COLOURS,
    parse_duration_to_seconds,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "performance"

# Minimum number of data rows we require before we treat a trace as
# "complete enough" to compute peak-concurrent + wallclock. PRIDE publishes
# several truncated traces (PXD062685 *5 = header + 2 rows; PXD070049/v2_3_2
# = header only). All complete traces have >= 18 rows; the threshold is
# conservatively set to 5 to keep the partition obvious in the data.
MIN_ROWS_COMPLETE = 5


# ---------------------------------------------------------------------------
# Per-trace helpers
# ---------------------------------------------------------------------------

def extract_step_name(full_name: str) -> str:
    """Return the step identifier from a Nextflow `name` cell.

    The trace stores fully-qualified process paths like
    `BIGBIO_QUANTMSDIANN:QUANTMSDIANN:INPUT_CHECK:SAMPLESHEET_CHECK (PXD049412.sdrf.tsv)`.
    The step is the last `:`-separated segment, with any parenthesised input
    argument stripped. Rows without an argument (e.g. `SUMMARY_PIPELINE`)
    return the bare step name.
    """
    s = (full_name or "").strip()
    if not s:
        return ""
    # Strip parenthesised tag if present (last open-paren onward).
    paren = s.find(" (")
    if paren != -1:
        s = s[:paren]
    # Last `:`-separated segment.
    if ":" in s:
        s = s.rsplit(":", 1)[1]
    return s.strip()


def load_report_window_data(report_path: Path) -> pd.DataFrame:
    """Read a Nextflow-generated `nextflow_report.html` and return its
    embedded `window.data.trace` table as a normalised DataFrame.

    Nextflow embeds the full per-task trace inside the HTML report as a
    JavaScript object literal. We brace-balance the `window.data = {...}`
    assignment, repair its invalid `\\'` escapes (Nextflow over-escapes
    single quotes inside JSON strings), and lift the `trace` array.

    The report is strictly better than `nextflow_trace.txt` for our purposes:
    it's published for every dataset (under `pipeline_info/` for cell-line
    runs, alongside the dataset for benchmark runs), and it contains the
    full task list even when the sibling `trace.txt` is upstream-truncated.

    Returns columns: `step`, `status`, `submit` (epoch seconds, float),
    `duration_s` (float), `finish` (epoch seconds, float). Returns an empty
    DataFrame with those columns when the HTML cannot be parsed or the
    trace array is empty."""
    import json
    cols = ["step", "status", "submit", "duration_s", "finish"]
    try:
        html = report_path.read_text(encoding="utf-8")
    except (FileNotFoundError, OSError):
        return pd.DataFrame(columns=cols)
    start = html.find("window.data = ")
    if start < 0:
        return pd.DataFrame(columns=cols)
    i = html.find("{", start)
    depth, in_str, esc, end = 0, False, False, -1
    for j in range(i, len(html)):
        ch = html[j]
        if esc:
            esc = False
            continue
        if ch == "\\":
            esc = True
            continue
        if ch == '"':
            in_str = not in_str
            continue
        if in_str:
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    if end < 0:
        return pd.DataFrame(columns=cols)
    blob = html[i:end].replace("\\'", "'")
    try:
        data = json.loads(blob)
    except json.JSONDecodeError:
        return pd.DataFrame(columns=cols)
    trace = data.get("trace", [])
    if not trace:
        return pd.DataFrame(columns=cols)
    rows: list[dict] = []
    for rec in trace:
        # submit/duration are epoch-ms strings in the report's JSON; convert
        # to seconds. start/complete are also available but submit+duration
        # matches the trace.txt semantics we already use elsewhere.
        submit_ms = rec.get("submit")
        duration_ms = rec.get("duration")
        try:
            submit_s = float(submit_ms) / 1000.0 if submit_ms is not None else float("nan")
        except (TypeError, ValueError):
            submit_s = float("nan")
        try:
            duration_s = float(duration_ms) / 1000.0 if duration_ms is not None else 0.0
        except (TypeError, ValueError):
            duration_s = 0.0
        # Step name: prefer the `process` field (already the colon-qualified
        # path with no input-arg suffix), fall back to extracting from name.
        process = rec.get("process") or rec.get("name") or ""
        step = extract_step_name(process)
        rows.append({
            "step": step,
            "status": rec.get("status") or "",
            "submit": submit_s,
            "duration_s": duration_s,
            "finish": submit_s + duration_s,
        })
    return pd.DataFrame(rows, columns=cols)


def load_trace(trace_path: Path) -> pd.DataFrame:
    """Read a nextflow_trace.txt as a normalised DataFrame.

    Returns columns: `step`, `status`, `submit`, `duration_s`, `finish`.
    Drops rows with unparseable submit timestamps. Empty traces (header
    only) return an empty DataFrame with the expected columns. FAILED tasks
    (Nextflow retry attempts) are kept so that callers can decide whether
    to include them: the parallelism analysis counts them (they did occupy
    cluster slots while running), the per-step duration analysis filters to
    COMPLETED.
    """
    df = pd.read_csv(trace_path, sep="\t", dtype=str)
    if df.empty:
        return pd.DataFrame(
            columns=["step", "status", "submit", "duration_s", "finish"],
        )
    submit = pd.to_datetime(df["submit"], errors="coerce")
    dur_s = df["duration"].fillna("").map(parse_duration_to_seconds)
    step = df["name"].fillna("").map(extract_step_name)
    status = df.get("status", pd.Series([""] * len(df))).fillna("")
    out = pd.DataFrame({
        "step": step,
        "status": status,
        "submit": submit,
        "duration_s": dur_s,
    })
    out = out.dropna(subset=["submit"]).reset_index(drop=True)
    out["finish"] = out["submit"] + pd.to_timedelta(out["duration_s"], unit="s")
    return out


def peak_concurrent_tasks(df: pd.DataFrame) -> tuple[int, float]:
    """Compute (peak, median) concurrent-task counts from a normalised
    trace DataFrame produced by `load_trace`.

    Peak = max number of intervals [submit, finish] that overlap at any
    instant. Median = median concurrency *over the union of event points*
    (so it counts the concurrency seen at each distinct submit/finish
    instant). This biases the median toward periods with more activity,
    which is what we want for a "how parallel did this workflow get" metric.

    Returns (0, 0.0) for an empty frame.
    """
    if df.empty:
        return (0, 0.0)
    events: list[tuple[pd.Timestamp, int]] = []
    for _, row in df.iterrows():
        events.append((row["submit"], +1))
        events.append((row["finish"], -1))
    # Sort by time; tie-break on delta so finishes process before new starts
    # at the same instant (consistent with closed-on-left intervals).
    events.sort(key=lambda x: (x[0], x[1]))
    concurrency = 0
    peak = 0
    samples: list[int] = []
    for _, delta in events:
        concurrency += delta
        peak = max(peak, concurrency)
        if concurrency > 0:
            samples.append(concurrency)
    if not samples:
        return (peak, 0.0)
    med = float(pd.Series(samples).median())
    return (peak, med)


def trace_wallclock_seconds(df: pd.DataFrame) -> float:
    """`max(submit + duration) - min(submit)` in seconds. Returns 0.0 for an
    empty frame. Handles both pandas Timestamp/Timedelta (from the legacy
    nextflow_trace.txt loader) and plain epoch-second floats (from the new
    nextflow_report.html loader)."""
    if df.empty:
        return 0.0
    span = df["finish"].max() - df["submit"].min()
    if hasattr(span, "total_seconds"):
        return float(span.total_seconds())
    return float(span)


def aggregate_step_durations(
    traces: Iterable[pd.DataFrame],
    *,
    completed_only: bool = True,
) -> dict[str, list[float]]:
    """Aggregate per-task `duration_s` values across many traces, keyed by
    step name. Empty / missing-step rows are dropped. When `completed_only`
    is True (default) we drop FAILED rows so the duration distribution
    reflects task success times — important for the per-step box plot
    because PXD070049/v1_8_1 has ~14 failed SDRF_PARSING retry attempts
    that would otherwise dominate the distribution."""
    out: dict[str, list[float]] = {}
    for df in traces:
        if df.empty:
            continue
        sub = df
        if completed_only and "status" in df.columns:
            sub = df[df["status"] == "COMPLETED"]
        for step, dur in zip(sub["step"], sub["duration_s"]):
            if not step:
                continue
            out.setdefault(step, []).append(float(dur))
    return out


# ---------------------------------------------------------------------------
# Analysis enumeration
# ---------------------------------------------------------------------------

PRIDE_CELL_LINE_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines"
)

# Cell-line datasets: each has a single Nextflow run with a report under
# `pipeline_info/nextflow_report.html`. Versions are inferred from the
# diannsummary.log command line (all currently 2.5.0).
CELL_LINE_ANALYSES: dict[str, dict[str, str]] = {
    "PXD003539": {"instrument": "TripleTOF 5600", "version": "v2_5_0"},
    "PXD030304": {"instrument": "TripleTOF 6600", "version": "v2_5_0"},
    "PXD004701": {"instrument": "TripleTOF 5600", "version": "v2_5_0"},
}


def iter_analyses() -> Iterable[tuple[str, str, str, str, Path]]:
    """Yield (dataset, version, instrument, report_url, local_report_path)
    for every analysis we know about — 3 cell-line + 20 benchmark = 23 rows.
    Both groups serve `nextflow_report.html` but at different URL shapes:
    cell-lines under `<dataset>/pipeline_info/`; benchmarks under
    `<dataset>/<version>/`."""
    for dataset, info in CELL_LINE_ANALYSES.items():
        base = f"{PRIDE_CELL_LINE_BASE}/{dataset}"
        url = f"{base}/pipeline_info/nextflow_report.html"
        local = (DATA_DIR / dataset / "pipeline_info"
                 / "nextflow_report.html")
        yield dataset, info["version"], info["instrument"], url, local
    for dataset in BENCHMARK_DATASETS:
        for version in DIANN_VERSIONS:
            base = f"{BENCHMARK_BASE}/{dataset}/{version}"
            url = f"{base}/nextflow_report.html"
            local = (DATA_DIR / "quantmsdiann_benchmarks" / dataset
                     / version / "nextflow_report.html")
            yield dataset, version, BENCHMARK_INSTRUMENT.get(
                dataset, "unknown"
            ), url, local


# ---------------------------------------------------------------------------
# Row builders
# ---------------------------------------------------------------------------

def collect_parallelism_rows(*, fetch: bool = True) -> pd.DataFrame:
    """One row per analysis (3 cell-line + 20 benchmark = 23) reading the
    per-task trace from each `nextflow_report.html`'s embedded
    `window.data.trace` JSON. Reports below MIN_ROWS_COMPLETE rows (none
    expected in practice, but defensive against future truncations) are
    dropped."""
    rows: list[dict] = []
    for dataset, version, instrument, url, local in iter_analyses():
        if fetch:
            download_if_missing(url, local)
        df = load_report_window_data(local)
        n = int(len(df))
        if n < MIN_ROWS_COMPLETE:
            continue
        peak, med = peak_concurrent_tasks(df)
        wallclock = trace_wallclock_seconds(df)
        rows.append({
            "dataset": dataset,
            "version": version,
            "instrument": instrument,
            "n_tasks": n,
            "peak_concurrent": peak,
            "median_concurrent": med,
            "wallclock_seconds": wallclock,
            "source_file": str(local.relative_to(REPO_ROOT)),
        })
    return pd.DataFrame(rows)


def collect_step_runtime_rows(*, fetch: bool = True) -> tuple[
    dict[str, list[float]], pd.DataFrame,
]:
    """Aggregate per-step task durations across every analysis (3 cell-line
    + 20 benchmark = 23) by reading each `nextflow_report.html` trace.

    Returns (durations_by_step, summary_df) where summary_df has columns
    `step, n, p05_seconds, p25_seconds, p50_seconds, p75_seconds,
    p95_seconds, min_seconds, max_seconds`.
    """
    traces: list[pd.DataFrame] = []
    for _, _, _, url, local in iter_analyses():
        if fetch:
            download_if_missing(url, local)
        traces.append(load_report_window_data(local))
    durations = aggregate_step_durations(traces)

    summary: list[dict] = []
    for step, vals in durations.items():
        s = pd.Series(vals)
        summary.append({
            "step": step,
            "n": int(len(s)),
            "p05_seconds": float(s.quantile(0.05)),
            "p25_seconds": float(s.quantile(0.25)),
            "p50_seconds": float(s.quantile(0.50)),
            "p75_seconds": float(s.quantile(0.75)),
            "p95_seconds": float(s.quantile(0.95)),
            "min_seconds": float(s.min()),
            "max_seconds": float(s.max()),
        })
    df = pd.DataFrame(summary).sort_values("p50_seconds", ascending=False)
    return durations, df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_parallelism_scatter(
    df: pd.DataFrame,
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """One-panel scatter: x = peak concurrent tasks, y = wallclock hours,
    colour = instrument family. Truncated traces are already filtered out
    upstream by `collect_parallelism_rows`."""
    plot_df = df.copy()
    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    hours = plot_df["wallclock_seconds"] / 3600.0

    # Dot size held constant — n_runs is 6 for every benchmark so the size
    # dimension would be uninformative.
    DOT_SIZE = 140.0

    for instrument in sorted(plot_df["instrument"].unique()):
        mask = plot_df["instrument"] == instrument
        ax.scatter(
            plot_df.loc[mask, "peak_concurrent"],
            hours[mask],
            s=DOT_SIZE,
            c=INSTRUMENT_COLOURS.get(instrument, "#9e9e9e"),
            alpha=0.78,
            edgecolors="#222222",
            linewidths=0.6,
            label=instrument,
        )

    ax.set_xlabel("Peak concurrent tasks (Nextflow trace, log scale)",
                  fontsize=10)
    ax.set_ylabel("Workflow wallclock (hours)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)

    # Log x-axis: the benchmark analyses peak at 6-7 concurrent tasks while
    # the cell-line PXD004701 run hits ~300; linear would crush the
    # benchmarks against the y-axis.
    ax.set_xscale("log")
    xmin = max(1, int(plot_df["peak_concurrent"].min()) - 1)
    xmax = int(plot_df["peak_concurrent"].max()) * 1.5
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, max(hours) * 1.20 if len(hours) else 1.0)

    inst_legend = ax.legend(
        title="Instrument", loc="upper left", fontsize=8, title_fontsize=9,
        frameon=False,
    )
    ax.add_artist(inst_legend)

    # Per-instrument peak-concurrent median annotation (replaces the
    # informationless size legend from the threads scatter).
    by_inst = (
        plot_df.groupby("instrument")["peak_concurrent"]
        .median()
        .sort_values(ascending=False)
    )
    lines = ["Peak concurrent (median per instrument):"]
    for inst, val in by_inst.items():
        lines.append(f"  {inst}: {val:.0f}")
    # Truncated traces are no longer present in `df`; no in-figure note
    # needed. The dataset/version exclusions are documented in the spec.
    ax.text(
        0.98, 0.02, "\n".join(lines),
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=8, family="monospace",
        bbox=dict(facecolor="white", edgecolor="#cccccc", boxstyle="round,pad=0.4"),
    )

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def render_per_step_boxplot(
    durations: dict[str, list[float]],
    summary: pd.DataFrame,
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """Horizontal box plot of per-task durations, one row per step, ordered
    by descending median."""
    # Step order: descending median (same as summary).
    steps = summary["step"].tolist()
    data = [durations[s] for s in steps]

    # Figure height scales with the number of steps so labels don't clip.
    fig_h = max(3.5, 0.45 * len(steps) + 1.5)
    fig, ax = plt.subplots(figsize=(8.0, fig_h))

    # Box-plot with whiskers at the 5th/95th percentile, individual outlier
    # points kept visible.
    bp = ax.boxplot(
        data,
        vert=False,
        whis=(5, 95),
        showfliers=True,
        flierprops=dict(
            marker="o", markerfacecolor="#888888", markeredgecolor="none",
            markersize=3.5, alpha=0.6,
        ),
        medianprops=dict(color="#e53935", linewidth=1.6),
        boxprops=dict(color="#1e88e5", linewidth=1.2),
        whiskerprops=dict(color="#555555", linewidth=1.0),
        capprops=dict(color="#555555", linewidth=1.0),
        patch_artist=False,
    )
    del bp  # silence unused-var lint

    # Reverse y-tick labels so the slowest-median step ends up at the top.
    ax.set_yticks(range(1, len(steps) + 1))
    ax.set_yticklabels(steps, fontsize=9)
    ax.invert_yaxis()

    ax.set_xlabel("Task duration (s)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)

    # Log scale when the range exceeds two orders of magnitude.
    flat = [v for vals in data for v in vals if v > 0]
    if flat:
        ratio = max(flat) / max(min(flat), 1e-6)
        if ratio > 100:
            ax.set_xscale("log")
            # Keep the left edge above 0 so log scale is well-defined.
            ax.set_xlim(left=max(0.5, min(flat) * 0.7))

    ax.grid(True, axis="x", which="both", alpha=0.25, linestyle=":")

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# TSV writers
# ---------------------------------------------------------------------------

def write_parallelism_tsv(df: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "dataset", "version", "instrument", "n_tasks", "peak_concurrent",
        "median_concurrent", "wallclock_seconds", "source_file",
    ]
    out = df[cols].copy().sort_values(["dataset", "version"])
    out.to_csv(tsv_path, sep="\t", index=False)


def write_per_step_tsv(summary: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "step", "n", "p05_seconds", "p25_seconds", "p50_seconds",
        "p75_seconds", "p95_seconds", "min_seconds", "max_seconds",
    ]
    summary[cols].to_csv(tsv_path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    par_df = collect_parallelism_rows(fetch=True)
    write_parallelism_tsv(par_df, FIGURES_DIR / "parallelism_data.tsv")
    render_parallelism_scatter(
        par_df,
        FIGURES_DIR / "parallelism_vs_wallclock.pdf",
        FIGURES_DIR / "parallelism_vs_wallclock.png",
        FIGURES_DIR / "parallelism_vs_wallclock.svg",
    )

    durations, summary = collect_step_runtime_rows(fetch=True)
    write_per_step_tsv(summary, FIGURES_DIR / "runtime_per_step.tsv")
    render_per_step_boxplot(
        durations, summary,
        FIGURES_DIR / "runtime_per_step.pdf",
        FIGURES_DIR / "runtime_per_step.png",
        FIGURES_DIR / "runtime_per_step.svg",
    )

    print(
        f"Wrote {len(par_df)} parallelism rows (complete-trace analyses only) "
        f"and {len(summary)} step-summary rows."
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
