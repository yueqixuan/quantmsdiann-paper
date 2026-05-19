"""Performance scatter: workflow runtime vs DIA-NN threads.

One dot per analysis we have run end-to-end through the quantmsdiann nf-core
workflow:

  - 3 cell-line datasets (PXD003539, PXD030304, PXD004701) processed once
    with DIA-NN 2.5.0.
  - 4 ProteoBench benchmark datasets (PXD049412, PXD062685, PXD070049,
    ProteoBench_Module_7) processed once per DIA-NN release
    (1.8.1, 2.1.0, 2.2.0, 2.3.2, 2.5.0).

For each analysis we plot:

  x = DIA-NN `--threads N` from the quant_tables/diannsummary.log command line
  y = total Nextflow wallclock duration in hours, parsed from
      pipeline_info/pipeline_report.txt (the "duration: 3h 26m 38s" line)
  size = number of MS runs (count of `--f` arguments in the DIA-NN command
      line; matches SDRF row count for cell-line datasets)
  color = MS instrument family (TripleTOF 5600 / 6600, timsTOF SCP,
      Orbitrap Astral, ZenoTOF 7600)

Paper-ready output (no title, no footer) is written to
  analysis/figures/performance/runtime_vs_threads.{pdf,png,svg}

An auditable TSV (one row per dot) is written to
  analysis/figures/performance/runtime_data.tsv

We use pipeline_report.txt rather than aggregating nextflow_trace.txt because
some uploaded trace files are truncated (PXD062685 traces only contain the
SAMPLESHEET_CHECK + SDRF_PARSING rows; PXD070049/v2_3_2 is header-only),
whereas pipeline_report.txt always carries the run-level total. Both come
from the same Nextflow run, so the wallclock is identical when both are
intact.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis.figure_original_vs_quantmsdiann import download_if_missing

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "performance"

# ---------------------------------------------------------------------------
# URL maps (mirror the patterns already used by the per-dataset figure scripts)
# ---------------------------------------------------------------------------

CELL_LINE_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines"
)

BENCHMARK_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/proteobench/quantmsdiann_results"
)

CELL_LINE_DATASETS = ("PXD003539", "PXD004701", "PXD030304")
BENCHMARK_DATASETS = (
    "PXD049412", "PXD062685", "PXD070049", "ProteoBench_Module_7",
)
DIANN_VERSIONS = ("v1_8_1", "v2_1_0", "v2_2_0", "v2_3_2", "v2_5_0")

# Per-dataset instrument family. For cell-line datasets the SDRF
# comment[instrument] entry is read at runtime; this dict is only the fallback
# / authoritative value for benchmark datasets (where the ProteoBench module
# name carries the instrument identity that is not always preserved in the
# SDRF the workflow consumed). Verified against the diannsummary.log run file
# names per dataset (see `infer_instrument_from_runs`).
BENCHMARK_INSTRUMENT = {
    "PXD049412": "Orbitrap Astral",
    "PXD062685": "timsTOF SCP",
    "PXD070049": "ZenoTOF 7600",
    "ProteoBench_Module_7": "Orbitrap Astral",
}


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

# Nextflow's `duration:` field can be any subset of "<h>h <m>m <s>s" with
# optional milliseconds. Examples we observed in the wild:
#   "3h 26m 38s", "1h 8m 8s", "45m 13s", "13m 7s", "1m 9s", "45s", "2.5s".
_DURATION_PART_RE = re.compile(
    r"(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>h|m|s|ms)\b"
)

# `duration: 3h 26m 38s` line of pipeline_report.txt.
_PIPELINE_DURATION_RE = re.compile(
    r"duration:\s*([0-9hms\.\sM]+?)\s*\)"
)

# DIA-NN `--threads N` on the single-line command in diannsummary.log.
_THREADS_RE = re.compile(r"--threads\s+(\d+)\b")

# Each `--f <path>` on the same command line.
_F_ARG_RE = re.compile(r"--f\s+(\S+)")


def parse_duration_to_seconds(text: str) -> float:
    """Parse a Nextflow-style duration token like "3h 26m 38s" or "1m 9s"
    into seconds (float). Returns 0.0 for an empty string. Tokens may carry
    decimal seconds ("2.5s") and milliseconds ("250ms")."""
    s = text.strip()
    if not s:
        return 0.0
    total = 0.0
    matched = False
    for m in _DURATION_PART_RE.finditer(s):
        matched = True
        v = float(m.group("value"))
        unit = m.group("unit")
        if unit == "h":
            total += v * 3600.0
        elif unit == "m":
            total += v * 60.0
        elif unit == "s":
            total += v
        elif unit == "ms":
            total += v / 1000.0
    if not matched:
        raise ValueError(f"Unrecognised duration string: {text!r}")
    return total


def parse_pipeline_report_duration(report_path: Path) -> float:
    """Extract the Nextflow wallclock seconds from a pipeline_report.txt
    file. Raises ValueError if the `duration:` line cannot be found.

    Format (one line near the bottom):
      The workflow was completed at <iso8601> (duration: 3h 26m 38s)
    """
    with open(report_path, encoding="utf-8") as fh:
        for line in fh:
            m = _PIPELINE_DURATION_RE.search(line)
            if m:
                return parse_duration_to_seconds(m.group(1))
    raise ValueError(
        f"`duration: ... )` not found in {report_path}"
    )


def aggregate_nextflow_trace_duration(trace_path: Path) -> float:
    """Wallclock seconds inferred from a Nextflow trace.txt file by computing
    `max(submit + duration) - min(submit)` across rows.

    This is provided as a fallback / cross-check. We don't use it for the
    figure because several published traces in the benchmark set are
    truncated (header only or first two rows only); pipeline_report.txt is
    the authoritative source. Tested via fixtures.

    Returns 0.0 if the trace has zero data rows. Raises ValueError if a
    required column is missing on a non-empty trace.
    """
    df = pd.read_csv(trace_path, sep="\t", dtype=str)
    if df.empty:
        return 0.0
    for col in ("submit", "duration"):
        if col not in df.columns:
            raise ValueError(
                f"Trace {trace_path} missing required column {col!r}"
            )
    submit = pd.to_datetime(df["submit"], errors="coerce")
    dur_s = df["duration"].fillna("").map(parse_duration_to_seconds)
    finish = submit + pd.to_timedelta(dur_s, unit="s")
    if submit.isna().all():
        return 0.0
    span = finish.max() - submit.min()
    return float(span.total_seconds())


def parse_diann_command(log_path: Path) -> tuple[int, int, list[str]]:
    """Return (threads, n_files, run_basenames) parsed from a diannsummary.log
    file. The DIA-NN command is always on a single early line starting with
    `diann --lib ...`.

    Raises ValueError if the command line cannot be located or `--threads`
    is missing.
    """
    cmd: str | None = None
    with open(log_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("diann ") and "--threads" in line:
                cmd = line
                break
    if cmd is None:
        raise ValueError(
            f"DIA-NN command line not found in {log_path}"
        )
    m = _THREADS_RE.search(cmd)
    if not m:
        raise ValueError(f"--threads N flag not found in {log_path}")
    threads = int(m.group(1))
    runs = [m.group(1) for m in _F_ARG_RE.finditer(cmd)]
    return threads, len(runs), runs


def infer_instrument_from_runs(run_names: Iterable[str]) -> str | None:
    """Heuristic instrument-family classifier from run/raw file basenames.
    Used as a sanity check when SDRF is not available. Returns None if the
    name pattern is not recognised."""
    if not run_names:
        return None
    sample = " ".join(list(run_names)[:8]).lower()
    if "timstof" in sample or "ttscp" in sample or sample.endswith(".d"):
        return "timsTOF SCP"
    if "zenotof" in sample or "zeno" in sample:
        return "ZenoTOF 7600"
    if "astral" in sample:
        return "Orbitrap Astral"
    return None


def parse_sdrf_instrument(sdrf_path: Path) -> str | None:
    """Return the comment[instrument] value (NT= field stripped) from an SDRF.
    Returns None if the column or NT= value is absent. If multiple distinct
    instruments appear, returns the most common one (PXD030304 SDRF mixes
    6 instrument identifiers but all map to the same NT=TripleTOF 6600)."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str, nrows=10000)
    col = None
    for c in df.columns:
        if c.strip().lower() == "comment[instrument]":
            col = c
            break
    if col is None:
        return None
    vals = df[col].dropna().astype(str)
    if vals.empty:
        return None
    # NT=<name>;AC=<accession> -> <name>
    def _strip(v: str) -> str:
        if "NT=" in v:
            inner = v.split("NT=", 1)[1]
            return inner.split(";", 1)[0].strip()
        return v.strip()
    cleaned = vals.map(_strip)
    return cleaned.value_counts().idxmax()


# ---------------------------------------------------------------------------
# Row builder
# ---------------------------------------------------------------------------

def collect_runtime_rows(*, fetch: bool = True) -> pd.DataFrame:
    """Build the per-analysis runtime table by reading (and optionally
    fetching) the small set of text files we need.

    Columns:
      dataset, version, instrument, threads, n_runs, wallclock_seconds,
      source_file_for_duration, source_file_for_threads
    """
    rows: list[dict] = []

    for dataset in CELL_LINE_DATASETS:
        dset_dir = DATA_DIR / dataset
        report = dset_dir / "pipeline_info" / "pipeline_report.txt"
        log = dset_dir / "diannsummary.log"
        sdrf = dset_dir / f"{dataset}.sdrf.tsv"
        if fetch:
            download_if_missing(
                f"{CELL_LINE_BASE}/{dataset}/pipeline_info/pipeline_report.txt",
                report,
            )
            download_if_missing(
                f"{CELL_LINE_BASE}/{dataset}/quant_tables/diannsummary.log",
                log,
            )
            download_if_missing(
                f"{CELL_LINE_BASE}/{dataset}/sdrf/{dataset}.sdrf.tsv",
                sdrf,
            )
        secs = parse_pipeline_report_duration(report)
        threads, n_files, run_names = parse_diann_command(log)
        instrument = parse_sdrf_instrument(sdrf) if sdrf.exists() else None
        if instrument is None:
            instrument = infer_instrument_from_runs(run_names) or "unknown"
        rows.append({
            "dataset": dataset,
            "version": "v2_5_0",
            "instrument": instrument,
            "threads": threads,
            "n_runs": n_files,
            "wallclock_seconds": secs,
            "source_file_for_duration": str(report.relative_to(REPO_ROOT)),
            "source_file_for_threads": str(log.relative_to(REPO_ROOT)),
        })

    for dataset in BENCHMARK_DATASETS:
        for version in DIANN_VERSIONS:
            base = f"{BENCHMARK_BASE}/{dataset}/{version}"
            dset_dir = DATA_DIR / "quantmsdiann_benchmarks" / dataset / version
            report = dset_dir / "pipeline_info" / "pipeline_report.txt"
            log = dset_dir / "quant_tables" / "diannsummary.log"
            if fetch:
                download_if_missing(
                    f"{base}/pipeline_info/pipeline_report.txt", report,
                )
                download_if_missing(
                    f"{base}/quant_tables/diannsummary.log", log,
                )
            secs = parse_pipeline_report_duration(report)
            threads, n_files, run_names = parse_diann_command(log)
            instrument = (
                infer_instrument_from_runs(run_names)
                or BENCHMARK_INSTRUMENT.get(dataset, "unknown")
            )
            rows.append({
                "dataset": dataset,
                "version": version,
                "instrument": instrument,
                "threads": threads,
                "n_runs": n_files,
                "wallclock_seconds": secs,
                "source_file_for_duration": str(report.relative_to(REPO_ROOT)),
                "source_file_for_threads": str(log.relative_to(REPO_ROOT)),
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

# Stable instrument colour palette. Keys must match `BENCHMARK_INSTRUMENT`
# values and the SDRF NT= strings exactly so the same dot color is reused
# across datasets.
INSTRUMENT_COLOURS = {
    "TripleTOF 5600": "#fb8c00",   # orange
    "TripleTOF 6600": "#e53935",   # red
    "timsTOF SCP":   "#8e24aa",    # purple
    "Orbitrap Astral": "#1e88e5",  # blue
    "ZenoTOF 7600":  "#00897b",    # teal
    "unknown":       "#9e9e9e",    # grey
}


def render_runtime_scatter(
    df: pd.DataFrame,
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """One-panel scatter: x = threads, y = wallclock (hours), size = n_runs,
    color = instrument family. Paper-ready: no title, no footer."""
    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    # Size mapping: log-scaled so 6-run benchmarks and 1700-run PXD030304
    # are simultaneously legible. Marker area in points^2 = base + slope *
    # log10(n_runs). Tuned so n=6 -> ~60 pt^2 and n=1734 -> ~600 pt^2.
    def _size(n: int) -> float:
        import math
        return 60.0 + 220.0 * math.log10(max(n, 1))

    hours = df["wallclock_seconds"] / 3600.0
    sizes = df["n_runs"].map(_size)

    # Plot each instrument family separately so we get a clean colour legend.
    for instrument in sorted(df["instrument"].unique()):
        mask = df["instrument"] == instrument
        ax.scatter(
            df.loc[mask, "threads"],
            hours[mask],
            s=sizes[mask],
            c=INSTRUMENT_COLOURS.get(instrument, "#9e9e9e"),
            alpha=0.78,
            edgecolors="#222222",
            linewidths=0.6,
            label=instrument,
        )

    ax.set_xlabel("DIA-NN threads (--threads)", fontsize=10)
    ax.set_ylabel("Workflow wallclock (hours)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)

    # x-axis: pad a bit on both sides so dots aren't clipped.
    tmin = int(df["threads"].min()) - 2
    tmax = int(df["threads"].max()) + 4
    ax.set_xlim(max(0, tmin), tmax)
    ax.set_ylim(0, max(hours) * 1.18)

    # Two-part legend: instrument color + n_runs size key.
    inst_legend = ax.legend(
        title="Instrument", loc="upper left", fontsize=8, title_fontsize=9,
        frameon=False,
    )
    ax.add_artist(inst_legend)

    from matplotlib.lines import Line2D
    size_handles = []
    # Span the observed n_runs range (6 ProteoBench replicates -> 5800 cell
    # lines x replicates in PXD030304). Keys spaced ~one decade apart for a
    # legible dot-size key.
    for n_demo in (6, 120, 300, 5800):
        size_handles.append(Line2D(
            [0], [0], marker="o", color="w", markerfacecolor="#bdbdbd",
            markeredgecolor="#222222", markeredgewidth=0.6,
            markersize=(_size(n_demo) ** 0.5),
            label=f"{n_demo:,} MS runs",
        ))
    ax.legend(
        handles=size_handles, title="Dataset size", loc="lower right",
        fontsize=8, title_fontsize=9, frameon=False, labelspacing=1.2,
        borderpad=1.0, handletextpad=1.0,
    )

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def write_runtime_tsv(df: pd.DataFrame, tsv_path: Path) -> None:
    """Write an auditable wide-format TSV (one row per dot)."""
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "dataset", "version", "instrument", "threads", "n_runs",
        "wallclock_seconds", "source_file_for_duration",
        "source_file_for_threads",
    ]
    out = df[cols].copy().sort_values(["dataset", "version"])
    out.to_csv(tsv_path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df = collect_runtime_rows(fetch=True)
    write_runtime_tsv(df, FIGURES_DIR / "runtime_data.tsv")
    render_runtime_scatter(
        df,
        FIGURES_DIR / "runtime_vs_threads.pdf",
        FIGURES_DIR / "runtime_vs_threads.png",
        FIGURES_DIR / "runtime_vs_threads.svg",
    )
    print(f"Wrote {len(df)} rows to runtime_data.tsv")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
