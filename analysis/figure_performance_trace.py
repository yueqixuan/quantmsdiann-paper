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
    parse_pipeline_report_duration,
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


# ---------------------------------------------------------------------------
# Resource-column parsers (memory + CPU). Nextflow's nextflow_trace.txt
# stores these in human-readable units; the embedded report JSON does not
# carry them, so the resources figure can only read trace.txt.
# ---------------------------------------------------------------------------

_SIZE_UNITS = {
    "": 1,
    "B": 1,
    "KB": 1024,
    "MB": 1024 ** 2,
    "GB": 1024 ** 3,
    "TB": 1024 ** 4,
}


def parse_size_to_bytes(text: str) -> float:
    """Parse a Nextflow size string like `"165 MB"` / `"4.7 GB"` / `"3.2 KB"`
    to bytes. Empty / `"-"` / unparseable strings return NaN."""
    s = (text or "").strip()
    if not s or s in {"-", "0"}:
        return float("nan") if not s or s == "-" else 0.0
    parts = s.split()
    if len(parts) == 1:
        # Bare number (bytes).
        try:
            return float(parts[0])
        except ValueError:
            return float("nan")
    try:
        value = float(parts[0])
    except ValueError:
        return float("nan")
    unit = parts[1].strip().upper()
    factor = _SIZE_UNITS.get(unit)
    if factor is None:
        return float("nan")
    return value * factor


def parse_pct_cpu(text: str) -> float:
    """Parse a Nextflow `%cpu` cell like `"94.8%"` to a float (94.8). Empty
    / `"-"` strings return NaN; values are *not* capped at 100 because
    multi-thread tasks can legitimately exceed 100 %."""
    s = (text or "").strip().rstrip("%").strip()
    if not s or s == "-":
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def load_trace_resources(trace_path: Path) -> pd.DataFrame:
    """Read a `nextflow_trace.txt` and return per-task resource rows.

    Columns: `step`, `status`, `peak_rss_bytes`, `pct_cpu`,
    `duration_s`. Drops rows with empty step or status. FAILED rows
    are kept; callers filter as needed (the F2c box plot uses COMPLETED
    rows only — same convention as `aggregate_step_durations`).

    Empty / header-only traces return an empty DataFrame with the
    expected columns."""
    cols = ["step", "status", "peak_rss_bytes", "pct_cpu", "duration_s"]
    try:
        df = pd.read_csv(trace_path, sep="\t", dtype=str)
    except (FileNotFoundError, OSError, pd.errors.EmptyDataError):
        return pd.DataFrame(columns=cols)
    if df.empty:
        return pd.DataFrame(columns=cols)
    needed = {"name", "status", "%cpu", "peak_rss", "duration"}
    if not needed.issubset(df.columns):
        return pd.DataFrame(columns=cols)
    out = pd.DataFrame({
        "step": df["name"].fillna("").map(extract_step_name),
        "status": df["status"].fillna(""),
        "peak_rss_bytes": df["peak_rss"].fillna("").map(parse_size_to_bytes),
        "pct_cpu": df["%cpu"].fillna("").map(parse_pct_cpu),
        "duration_s": df["duration"].fillna("").map(parse_duration_to_seconds),
    })
    out = out[out["step"].astype(bool)].reset_index(drop=True)
    return out


def aggregate_step_resources(
    traces: Iterable[pd.DataFrame],
    *,
    completed_only: bool = True,
) -> dict[str, dict[str, list[float]]]:
    """Aggregate per-task `peak_rss_bytes` and `pct_cpu` across traces,
    keyed first by step, then by metric. Mirrors
    `aggregate_step_durations` but for resource columns. Returns a dict
    `{step: {"peak_rss_bytes": [...], "pct_cpu": [...]}}`. NaN values
    are dropped per-metric so a row missing one column still contributes
    the other."""
    out: dict[str, dict[str, list[float]]] = {}
    for df in traces:
        if df.empty:
            continue
        sub = df
        if completed_only and "status" in df.columns:
            sub = df[df["status"] == "COMPLETED"]
        for step, rss, cpu in zip(
            sub["step"], sub["peak_rss_bytes"], sub["pct_cpu"]
        ):
            if not step:
                continue
            bucket = out.setdefault(
                step, {"peak_rss_bytes": [], "pct_cpu": []}
            )
            if pd.notna(rss):
                bucket["peak_rss_bytes"].append(float(rss))
            if pd.notna(cpu):
                bucket["pct_cpu"].append(float(cpu))
    return out


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


_PARAMS_FILENAME_RE = __import__("re").compile(
    r"params_(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})\.json"
)
_PIPELINE_REPORT_COMPLETION_RE = __import__("re").compile(
    r"completed at\s+(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?)"
)


def list_params_timestamps(pipeline_info_url: str) -> list:
    """List `params_<ISO>.json` timestamps from a pipeline_info/ directory
    listing. Each `params_*.json` file corresponds to one Nextflow
    invocation; the earliest is the start of the user's wallclock.

    Returns a list of `datetime.datetime` objects (naive, in the cluster's
    local timezone) sorted ascending. The runs are typically run within the
    same TZ so naive comparison is fine for span computation."""
    import re, requests, datetime as _dt
    listing = requests.get(pipeline_info_url, timeout=30).text
    out: list[_dt.datetime] = []
    seen: set[str] = set()
    for m in _PARAMS_FILENAME_RE.finditer(listing):
        token = m.group(1)
        if token in seen:
            continue
        seen.add(token)
        out.append(_dt.datetime.strptime(token, "%Y-%m-%d_%H-%M-%S"))
    return sorted(out)


def parse_pipeline_report_completion(report_path):
    """Read `pipeline_info/pipeline_report.txt` and return the workflow
    completion datetime (naive, cluster-local). Raises ValueError if the
    line isn't present (e.g. an aborted run)."""
    import datetime as _dt
    text = report_path.read_text(encoding="utf-8")
    m = _PIPELINE_REPORT_COMPLETION_RE.search(text)
    if not m:
        raise ValueError(
            f"'completed at <ISO>' line not found in {report_path}"
        )
    # Parse ISO timestamp (drop microseconds if present, since strptime's
    # %f handles only up to 6 digits and the pipeline report sometimes
    # emits 9-digit nanosecond precision).
    iso = m.group(1)
    if "." in iso:
        head, frac = iso.split(".", 1)
        return _dt.datetime.fromisoformat(head + "." + frac[:6])
    return _dt.datetime.fromisoformat(iso)


def total_wallclock_with_resumes_seconds(
    pipeline_info_url: str,
    report_path,
) -> float:
    """Total wallclock the user actually waited, including all `-resume`
    re-runs: earliest `params_*.json` timestamp to `pipeline_report.txt`'s
    "completed at" datetime. For single-invocation workflows this matches
    the report's reported duration; for multi-invocation workflows it's
    larger because the cluster idle/queue time between resumes counts."""
    starts = list_params_timestamps(pipeline_info_url)
    if not starts:
        raise ValueError(
            f"No params_*.json files listed at {pipeline_info_url}"
        )
    completed = parse_pipeline_report_completion(report_path)
    return (completed - starts[0]).total_seconds()


def count_ms_runs(dataset_id: str, sdrf_path, log_path) -> int:
    """Number of MS runs in the analysis. For cell-line datasets we count
    SDRF rows (one per data file). For benchmarks we count `--f` arguments
    in the DIA-NN command line at the head of `diannsummary.log`. Either
    source is correct for the matching dataset family; we fall back from
    SDRF to the DIA-NN log if the SDRF is unavailable."""
    if sdrf_path.exists():
        n = 0
        with open(sdrf_path, "r", encoding="utf-8") as fh:
            for i, _ in enumerate(fh):
                n = i
        return n  # header excluded
    from analysis.figure_performance_runtime import parse_diann_command
    _, n_files, _ = parse_diann_command(log_path)
    return n_files

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
    """One row per analysis (3 cell-line + 20 benchmark = 23) with:

    - n_runs: number of MS data files in the analysis (SDRF rows for
      cell-lines, --f count from diannsummary.log for benchmarks).
    - total_wallclock_seconds: earliest `params_*.json` timestamp to
      `pipeline_report.txt`'s "completed at" datetime — captures the real
      wall time the user waited, including queue time between `-resume`
      re-runs. For single-invocation workflows this matches the report's
      reported duration; for multi-resume workflows it's larger.
    - n_invocations: how many `params_*.json` invocations the dataset
      went through.
    - n_tasks: number of tasks observed in the report's trace JSON (low
      for `-resume`d cell-line runs since only the post-resume tasks are
      traced).
    """
    rows: list[dict] = []
    for dataset, version, instrument, url, local in iter_analyses():
        if fetch:
            download_if_missing(url, local)
        df = load_report_window_data(local)
        n_tasks = int(len(df))

        # pipeline_info/ directory URL (for params listing + report.txt).
        pipeline_info_url = url.rsplit("/", 1)[0]
        # Cell-line reports already live under pipeline_info/; benchmarks
        # have the report at the version root with pipeline_info/ as a
        # subdir. Compute the URL that hosts params_*.json + pipeline_report.txt.
        if pipeline_info_url.endswith("/pipeline_info"):
            params_listing_url = pipeline_info_url + "/"
            report_txt_url = pipeline_info_url + "/pipeline_report.txt"
            report_txt_local = local.parent / "pipeline_report.txt"
        else:
            params_listing_url = pipeline_info_url + "/pipeline_info/"
            report_txt_url = (pipeline_info_url
                              + "/pipeline_info/pipeline_report.txt")
            report_txt_local = (local.parent / "pipeline_info"
                                / "pipeline_report.txt")
        if fetch:
            download_if_missing(report_txt_url, report_txt_local)

        # Wallclock = the final run's "duration:" line from
        # pipeline_report.txt — well-defined and matches the real compute
        # cost of finishing the workflow (post-resume). The "span from
        # first params_*.json to last completion" metric we tried earlier
        # was confounded by benchmark re-tests done days apart (idle gaps
        # between unrelated invocations dwarfing the actual run time).
        # n_invocations is reported separately so the resume count is
        # visible without polluting Y.
        try:
            wallclock = parse_pipeline_report_duration(report_txt_local)
            starts = list_params_timestamps(params_listing_url)
            n_invocations = len(starts) if starts else 1
        except (ValueError, FileNotFoundError, OSError):
            continue

        # n_runs: SDRF rows for cell-lines, --f count for benchmarks.
        # SDRF caches are dataset-scoped from earlier downloads.
        if dataset in CELL_LINE_ANALYSES:
            sdrf_local = DATA_DIR / dataset / f"{dataset}.sdrf.tsv"
        else:
            sdrf_local = (DATA_DIR / "quantmsdiann_benchmarks" / dataset
                          / f"{dataset}.sdrf.tsv")
        log_local = (DATA_DIR / dataset / "diannsummary.log"
                     if dataset in CELL_LINE_ANALYSES
                     else DATA_DIR / "quantmsdiann_benchmarks" / dataset
                          / version / "quant_tables" / "diannsummary.log")
        try:
            n_runs = count_ms_runs(dataset, sdrf_local, log_local)
        except (FileNotFoundError, ValueError):
            n_runs = 0

        rows.append({
            "dataset": dataset,
            "version": version,
            "instrument": instrument,
            "n_runs": n_runs,
            "n_tasks_observed": n_tasks,
            "n_invocations": n_invocations,
            "wallclock_seconds": wallclock,
            "trace_complete": n_tasks >= MIN_ROWS_COMPLETE,
            "source_file": str(report_txt_local.relative_to(REPO_ROOT)),
        })
    return pd.DataFrame(rows)


def _iter_pxd071075_sweep_trace_paths() -> Iterable[Path]:
    """Yield local nextflow_trace.txt paths for the 5 PXD071075
    cluster-node sweep points. Tolerates both zero-padded (q010) and
    bare (q10) directory layouts. Missing sweep points are skipped
    silently — callers (F2b/F2c aggregators) are happy with whatever
    subset is staged."""
    sweep_dir = DATA_DIR / "queue_size_sweep"
    if not sweep_dir.exists():
        return
    for q in (10, 50, 100, 200, 300):
        for cand in (sweep_dir / f"q{q:03d}", sweep_dir / f"q{q}"):
            trace = cand / "nextflow_trace.txt"
            if trace.exists():
                yield trace
                break


def _load_pxd071075_sweep_traces() -> list[pd.DataFrame]:
    """Return the per-task DataFrames for the PXD071075 sweep traces
    in `load_trace` shape (`step, status, submit, duration_s, finish`).
    Used by `collect_step_runtime_rows` to pool sweep tasks into the
    F2b distribution alongside benchmark tasks."""
    return [load_trace(p) for p in _iter_pxd071075_sweep_trace_paths()]


def _load_pxd071075_sweep_resource_rows() -> list[pd.DataFrame]:
    """Return per-task resource DataFrames for the PXD071075 sweep
    (`step, status, peak_rss_bytes, pct_cpu, duration_s`). Used by
    `collect_step_resource_rows` to pool sweep tasks into the F2c
    distribution alongside benchmark tasks."""
    return [
        load_trace_resources(p)
        for p in _iter_pxd071075_sweep_trace_paths()
    ]


def collect_pxd071075_sweep_rows() -> pd.DataFrame:
    """Add the PXD071075 single-cell `executor.queueSize` sweep points
    to the parallelism scatter. Same schema as
    `collect_parallelism_rows` so the two frames can be concatenated.

    Each sweep point sits at the same n_runs (PXD071075's 2,310 input
    runs) and the same instrument (Orbitrap Eclipse), but varies the
    Nextflow queueSize. The resulting vertical strip in the F2a
    scatter is the "elastic" wallclock dimension — orthogonal to the
    workload-scaling axis the other cohorts demonstrate. The extra
    `queue_size` column lets the renderer annotate each point with
    its queueSize without polluting non-sweep rows.

    Reads `data/queue_size_sweep/q<NNN>/nextflow_trace.txt` for each
    sweep point. Returns an empty DataFrame (right shape) when no
    sweep data is staged."""
    sweep_dir = DATA_DIR / "queue_size_sweep"
    rows: list[dict] = []
    if not sweep_dir.exists():
        return pd.DataFrame(columns=[
            "dataset", "version", "instrument", "n_runs",
            "n_tasks_observed", "n_invocations", "wallclock_seconds",
            "trace_complete", "source_file", "queue_size",
        ])
    # PXD071075's n_runs (count of --f flags) is constant across the
    # sweep. Pull it once from any sweep point's diannsummary.log.
    n_runs = 0
    for q in (10, 50, 100, 200, 300):
        for cand in (sweep_dir / f"q{q:03d}", sweep_dir / f"q{q}"):
            log = cand / "diannsummary.log"
            if log.exists():
                try:
                    with open(log, encoding="utf-8") as fh:
                        text = fh.read()
                    import re as _re
                    n_runs = len(_re.findall(r"--f \S+", text))
                except (FileNotFoundError, OSError):
                    n_runs = 0
                break
        if n_runs:
            break

    for q in (10, 50, 100, 200, 300):
        # Tolerate both zero-padded (q010) and bare (q10) layouts.
        trace_path = None
        for cand in (sweep_dir / f"q{q:03d}", sweep_dir / f"q{q}"):
            if (cand / "nextflow_trace.txt").exists():
                trace_path = cand / "nextflow_trace.txt"
                break
        if trace_path is None:
            continue
        df_trace = load_trace(trace_path)
        if df_trace.empty:
            continue
        wallclock = trace_wallclock_seconds(df_trace)
        n_tasks = int(len(df_trace))
        rows.append({
            "dataset": "PXD071075",
            "version": f"v2_5_0_q{q:03d}",
            "instrument": "Orbitrap Eclipse",
            "n_runs": n_runs,
            "n_tasks_observed": n_tasks,
            "n_invocations": 1,
            "wallclock_seconds": float(wallclock),
            "trace_complete": n_tasks >= MIN_ROWS_COMPLETE,
            "source_file": str(trace_path.relative_to(REPO_ROOT)),
            "queue_size": q,
        })
    return pd.DataFrame(rows)


def collect_step_runtime_rows(*, fetch: bool = True) -> tuple[
    dict[str, list[float]], pd.DataFrame,
]:
    """Aggregate per-step task durations across the **20 benchmark
    analyses + the 5 PXD071075 cluster-node sweep runs**. The 3
    cell-line analyses (PXD003539 / PXD030304 / PXD004701) are still
    excluded — PRIDE publishes their `nextflow_report.html` but not
    `nextflow_trace.txt` and the report's embedded JSON conflates
    different scales (a single 300-file PXD004701 run contributes 300
    PRELIMINARY_ANALYSIS tasks against 20×6=120 benchmark tasks on the
    same step). The PXD071075 sweep is added because per-task duration
    for a given step is roughly invariant to queueSize — only the
    workflow-level parallelism shifts — so pooling the 5 sweep
    traces gives the distribution more weight at the real-cohort
    scale (~2,310 input files) without distorting the per-step shape.

    Returns (durations_by_step, summary_df) where summary_df has columns
    `step, n, p05_seconds, p25_seconds, p50_seconds, p75_seconds,
    p95_seconds, min_seconds, max_seconds`.
    """
    traces: list[pd.DataFrame] = []
    for dataset, _, _, url, local in iter_analyses():
        if dataset in CELL_LINE_ANALYSES:
            continue
        if fetch:
            download_if_missing(url, local)
        traces.append(load_report_window_data(local))
    # Pool the PXD071075 cluster-node sweep traces alongside the
    # benchmarks. Already cached locally; trace.txt has the same
    # step / status / duration columns we need.
    traces.extend(_load_pxd071075_sweep_traces())
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
    svg_path: Path | None = None,
    *,
    ax: plt.Axes | None = None,
    legend_ncol: int = 4,
    legend_bbox_y: float = -0.14,
    composite: bool = False,
    show_legend: bool = True,
) -> None:
    """One-panel scatter: x = number of MS data files (log), y = final-run
    wallclock from pipeline_report.txt (hours), colour = instrument family.
    All dots filled with the same size — the `-resume` invocation count is
    kept in `parallelism_data.tsv` for the audit trail but isn't visually
    encoded because the re-runs reflect SDRF iteration during development,
    not workflow reliability.

    Pass `ax` to draw into an existing axes (composite figures); omit
    `svg_path` in that mode."""
    plot_df = df.copy()
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(6.6, 4.6))

    hours = plot_df["wallclock_seconds"] / 3600.0
    dot_size = 55.0 if composite else 150.0
    label_size = 8 if composite else 10
    tick_size = 7 if composite else 9
    ann_size = 6 if composite else 7

    for instrument in sorted(plot_df["instrument"].unique()):
        mask = plot_df["instrument"] == instrument
        ax.scatter(
            plot_df.loc[mask, "n_runs"],
            hours[mask],
            s=dot_size,
            c=INSTRUMENT_COLOURS.get(instrument, "#9e9e9e"),
            alpha=0.85,
            edgecolors="#222222",
            linewidths=0.6,
            label=instrument,
        )

    # PXD071075 sweep points share an n_runs (~2,310) but span q=10..300
    # in queueSize. Annotate each with its queueSize so the vertical
    # strip reads correctly. For non-sweep cohorts annotate the dataset
    # name once per cohort, on the topmost dot only (so PXD003539 +
    # PXD030304 + PXD004701 each carry one label, not five).
    is_sweep = (
        plot_df.get("queue_size").notna()
        if "queue_size" in plot_df.columns
        else pd.Series([False] * len(plot_df))
    )
    # Per-cohort name annotation — one per non-sweep cohort, attached
    # to the dot with the largest wallclock.
    seen_datasets: set[str] = set()
    for _, row in (
        plot_df[~is_sweep]
        .sort_values("wallclock_seconds", ascending=False)
        .iterrows()
    ):
        ds = row["dataset"]
        if not ds.startswith("PXD") or row["n_runs"] < 100:
            continue
        if ds in seen_datasets:
            continue
        seen_datasets.add(ds)
        ax.annotate(
            ds,
            xy=(row["n_runs"], row["wallclock_seconds"] / 3600.0),
            xytext=(-8, 5), textcoords="offset points",
            fontsize=ann_size, color="#444444", ha="right",
        )
    # PXD071075 sweep annotation.
    # - If exactly one sweep row is shown (the F2a default — q=300
    #   production run), label the dot "PXD071075 (<q> nodes)" once.
    # - If multiple sweep rows are shown, fall back to the original
    #   vertical-strip annotation (per-dot q-label + a single "sweep"
    #   header) — used by callers that pass the full sweep frame.
    sweep_rows = plot_df[is_sweep] if is_sweep.any() else plot_df.iloc[0:0]
    if len(sweep_rows) == 1:
        row = sweep_rows.iloc[0]
        ax.annotate(
            f"{row['dataset']} ({int(row['queue_size'])} nodes)",
            xy=(row["n_runs"], row["wallclock_seconds"] / 3600.0),
            xytext=(-8, 5), textcoords="offset points",
            fontsize=ann_size, color="#1a237e", ha="right",
        )
    elif len(sweep_rows) > 1:
        for _, row in sweep_rows.iterrows():
            ax.annotate(
                f"q={int(row['queue_size'])}",
                xy=(row["n_runs"], row["wallclock_seconds"] / 3600.0),
                xytext=(6, 0), textcoords="offset points",
                fontsize=ann_size, color="#1a237e", ha="left", va="center",
            )
        top_sweep = sweep_rows.sort_values(
            "wallclock_seconds", ascending=False,
        ).iloc[0]
        ax.annotate(
            "PXD071075 (queueSize sweep)",
            xy=(top_sweep["n_runs"], top_sweep["wallclock_seconds"] / 3600.0),
            xytext=(-8, 6), textcoords="offset points",
            fontsize=ann_size, color="#1a237e", ha="right",
        )


    ax.set_xlabel("Number of MS data files", fontsize=label_size)
    ax.set_ylabel(
        "Final-run wall-clock (hours)",
        fontsize=label_size,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=tick_size)
    ax.set_xscale("log")
    xmax = float(plot_df["n_runs"].max()) * 1.6
    ax.set_xlim(max(1, plot_df["n_runs"].min() * 0.7), xmax)
    ax.set_ylim(0, max(hours) * 1.18 if len(hours) else 1.0)
    # Explicit integer tick labels at round powers of 10 spanning the
    # data range (currently n_runs ∈ {6, 120, 300, 2310, 5798}) — avoids
    # matplotlib's default 10⁰ / 10¹ power-notation glyphs which read
    # poorly in vector exports.
    from matplotlib.ticker import FixedLocator, FixedFormatter
    xmin_data = float(plot_df["n_runs"].min())
    candidate_ticks = [10, 100, 1000, 10000]
    xticks = [t for t in candidate_ticks if xmin_data / 1.6 <= t <= xmax]
    if xticks:
        ax.xaxis.set_major_locator(FixedLocator(xticks))
        ax.xaxis.set_major_formatter(
            FixedFormatter([f"{t:,}" for t in xticks])
        )
        ax.xaxis.set_minor_locator(FixedLocator([]))

    # Single instrument-colour legend below the axes — no dot-size
    # dimension to encode, so no second legend needed.
    if show_legend:
        legend_fs = 6 if composite else 8
        legend_title_fs = 7 if composite else 9
        ncol = 2 if composite else legend_ncol
        bbox_y = -0.22 if composite else legend_bbox_y
        ax.legend(
            title="Instrument", loc="upper center",
            bbox_to_anchor=(0.5, bbox_y), fontsize=legend_fs,
            title_fontsize=legend_title_fs,
            frameon=False, borderaxespad=0.0, ncol=ncol,
            columnspacing=1.2 if composite else 1.6,
        )

    if own_fig:
        fig.tight_layout()
        assert svg_path is not None
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)


def render_per_step_boxplot(
    durations: dict[str, list[float]],
    summary: pd.DataFrame,
    svg_path: Path,
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
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# TSV writers
# ---------------------------------------------------------------------------

def write_parallelism_tsv(df: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "dataset", "version", "instrument", "n_runs", "n_tasks_observed",
        "n_invocations", "wallclock_seconds", "trace_complete",
        "source_file",
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
# F2c resources panel — memory + CPU per step
# ---------------------------------------------------------------------------

def cell_line_trace_local_path(dataset: str) -> Path:
    """Local path where a cell-line `nextflow_trace.txt` lands once
    collected from a fresh quantmsdiann rerun (experiment #11). PRIDE
    does not currently publish these files — they're collected with
    `nextflow -with-trace` on the cluster and staged here. If the file
    exists it is included in F2b/F2c/F2a automatically; if not, those
    figures stay benchmarks-only.

    See [docs/superpowers/specs/2026-05-20-experiment-11-cell-line-traces.md]
    for the runbook."""
    return DATA_DIR / dataset / "pipeline_info" / "nextflow_trace.txt"


def has_cell_line_traces() -> bool:
    """Return True iff every cell-line dataset has a local
    `nextflow_trace.txt`. Used by F2c/F2b's collectors to opt in to the
    broader cohort once experiment #11 lands. Until then the
    benchmarks-only scope holds."""
    return all(
        cell_line_trace_local_path(d).exists()
        for d in CELL_LINE_ANALYSES
    )


def collect_step_resource_rows(*, fetch: bool = True) -> tuple[
    dict[str, dict[str, list[float]]], pd.DataFrame,
]:
    """Aggregate per-task resource rows across the benchmark analyses
    + the 5 PXD071075 cluster-node sweep runs.

    Cell-line traces are included automatically when
    `has_cell_line_traces()` returns True — that is, once experiment #11
    has staged `data/<PXD>/pipeline_info/nextflow_trace.txt` for all
    three cell-line cohorts. PXD071075 sweep traces are always
    included when staged under `data/queue_size_sweep/q<N>/` — each
    sweep run executes the same workflow on the same 2,310 input
    files, so pooling all 5 runs increases the effective sample size
    of the per-step resource distribution without distorting it.

    Returns (resources_by_step, summary_df) where summary_df has
    columns `step, n_rss, peak_rss_p50_bytes, peak_rss_p95_bytes,
    pct_cpu_p50, pct_cpu_p95`.
    """
    include_cell_lines = has_cell_line_traces()
    traces: list[pd.DataFrame] = []
    for dataset, version, *_ in iter_analyses():
        if dataset in CELL_LINE_ANALYSES:
            if not include_cell_lines:
                continue
            trace_local = cell_line_trace_local_path(dataset)
            # No URL fetch path — these are local-only artefacts.
        else:
            trace_url = (
                f"{BENCHMARK_BASE}/{dataset}/{version}/nextflow_trace.txt"
            )
            trace_local = (
                DATA_DIR / "quantmsdiann_benchmarks" / dataset / version
                / "nextflow_trace.txt"
            )
            if fetch:
                download_if_missing(trace_url, trace_local)
        traces.append(load_trace_resources(trace_local))
    # Pool the PXD071075 cluster-node sweep traces alongside the
    # benchmarks. Each sweep point is the same workflow at a different
    # queueSize, so per-task resource usage is roughly equivalent
    # across runs — pooling all 5 grows the effective sample without
    # bias.
    traces.extend(_load_pxd071075_sweep_resource_rows())
    resources = aggregate_step_resources(traces)

    summary: list[dict] = []
    for step, metrics in resources.items():
        rss = pd.Series(metrics["peak_rss_bytes"])
        cpu = pd.Series(metrics["pct_cpu"])
        pct_cpu_p50 = float(cpu.median()) if len(cpu) else float("nan")
        pct_cpu_p95 = float(cpu.quantile(0.95)) if len(cpu) else float("nan")
        # Thread efficiency = raw %cpu / 12 (the DIA-NN --threads baseline).
        # Audited alongside raw %cpu so a reader can verify the divisor.
        eff_p50 = (
            pct_cpu_p50 / DIANN_THREADS_BASELINE
            if not pd.isna(pct_cpu_p50) else float("nan")
        )
        eff_p95 = (
            pct_cpu_p95 / DIANN_THREADS_BASELINE
            if not pd.isna(pct_cpu_p95) else float("nan")
        )
        summary.append({
            "step": step,
            "n_rss": int(len(rss)),
            "peak_rss_p50_bytes": (
                float(rss.median()) if len(rss) else float("nan")
            ),
            "peak_rss_p95_bytes": (
                float(rss.quantile(0.95)) if len(rss) else float("nan")
            ),
            "n_cpu": int(len(cpu)),
            "pct_cpu_p50": pct_cpu_p50,
            "pct_cpu_p95": pct_cpu_p95,
            "thread_efficiency_p50": eff_p50,
            "thread_efficiency_p95": eff_p95,
        })
    df = (
        pd.DataFrame(summary)
        .sort_values("peak_rss_p50_bytes", ascending=False)
        .reset_index(drop=True)
    )
    return resources, df


# quantmsdiann always launches DIA-NN with `--threads 12`, set in the
# pipeline's `withLabel: diann` block. The CPU panel uses this as the
# normalisation baseline so the axis reads as "fraction of the
# 12-thread allocation used", not raw %cpu (which is unreadable
# without knowing the Unix convention that 100 % == 1 core).
DIANN_THREADS_BASELINE = 12


def render_resources_boxplot(
    resources: dict[str, dict[str, list[float]]],
    summary: pd.DataFrame,
    svg_path: Path,
) -> None:
    """Two-panel horizontal box plot: (left) `peak_rss` per step in GB,
    (right) **threading efficiency** = `%cpu / (12 * 100 %)` per step.
    Step order is fixed by descending median `peak_rss`, mirrored
    across both panels so a reviewer can read the same step across
    both axes.

    The threading efficiency divisor (12) is `DIANN_THREADS_BASELINE`
    — quantmsdiann's `--threads` setting for the heavy DIA-NN steps.
    DIA-NN steps (INSILICO_LIBRARY_GENERATION / PRELIMINARY_ANALYSIS /
    INDIVIDUAL_ANALYSIS) should sit near 100 % here. Glue steps that
    Nextflow allocates 1-2 cores to (SDRF_PARSING, MSstats, raw
    conversion) sit near 1/12 ≈ 8 % on this scale — this is the
    expected pattern for single-threaded helpers, not a quantmsdiann
    deficiency. The note in the figure footer spells this out."""
    steps = summary["step"].tolist()
    if not steps:
        # No data at all — emit a minimal empty figure so callers don't
        # crash; F2c is benchmarks-only and the data is always present
        # in the published artefacts but defensive.
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "no resource rows", ha="center", va="center",
                transform=ax.transAxes, color="#888888")
        ax.set_axis_off()
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg_path)
        plt.close(fig)
        return
    rss_data = [
        [b / 1024 ** 3 for b in resources[s]["peak_rss_bytes"]]
        for s in steps
    ]
    # Threading efficiency = %cpu / (12 * 100 %), shown as a percentage
    # of the DIA-NN 12-thread allocation. Raw %cpu values are kept in
    # the TSV for auditability.
    cpu_data = [
        [v / DIANN_THREADS_BASELINE for v in resources[s]["pct_cpu"]]
        for s in steps
    ]
    fig, (ax_rss, ax_cpu) = plt.subplots(
        nrows=1, ncols=2,
        figsize=(11.5, max(3.0, 0.36 * len(steps) + 1.0)),
        sharey=True,
    )
    ax_rss.boxplot(
        rss_data, vert=False, widths=0.62,
        labels=steps, showfliers=False,
        medianprops=dict(color="#d62728", linewidth=1.4),
        boxprops=dict(color="#445566"),
        whiskerprops=dict(color="#445566"),
        capprops=dict(color="#445566"),
    )
    ax_rss.set_xlabel("Peak RSS per task (GB)", fontsize=10)
    ax_rss.spines["top"].set_visible(False)
    ax_rss.spines["right"].set_visible(False)
    ax_rss.tick_params(axis="both", labelsize=8)
    ax_rss.invert_yaxis()  # highest-memory step at the top
    # Annotate per-step n_rss to the right of each box.
    for i, step in enumerate(steps, start=1):
        n = len(resources[step]["peak_rss_bytes"])
        if n:
            xpos = max(rss_data[i - 1]) if rss_data[i - 1] else 0
            ax_rss.text(
                xpos * 1.04, i, f"n={n}",
                va="center", ha="left", fontsize=7, color="#666666",
            )

    ax_cpu.boxplot(
        cpu_data, vert=False, widths=0.62,
        labels=steps, showfliers=False,
        medianprops=dict(color="#d62728", linewidth=1.4),
        boxprops=dict(color="#445566"),
        whiskerprops=dict(color="#445566"),
        capprops=dict(color="#445566"),
    )
    ax_cpu.set_xlabel("Threading efficiency (%)", fontsize=10)
    # 100 % reference: full saturation of the 12-thread allocation.
    ax_cpu.axvline(
        100.0, color="#26a69a", linestyle="--", linewidth=0.9, zorder=1,
    )
    # 1/12 reference for single-threaded glue tasks Nextflow allocates
    # 1 core to. Drawn light so it doesn't dominate.
    ax_cpu.axvline(
        100.0 / DIANN_THREADS_BASELINE,
        color="#bdbdbd", linestyle=":", linewidth=0.8, zorder=1,
    )
    ax_cpu.spines["top"].set_visible(False)
    ax_cpu.spines["right"].set_visible(False)
    ax_cpu.tick_params(axis="both", labelsize=8)
    # Explanatory footnote intentionally not rendered into the figure;
    # the divisor (12-thread DIA-NN baseline) and the two reference
    # lines (100 % / 1/12) are spelled out in the manuscript methods
    # text.
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def write_per_step_resources_tsv(
    summary: pd.DataFrame, tsv_path: Path,
) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "step", "n_rss",
        "peak_rss_p50_bytes", "peak_rss_p95_bytes",
        "n_cpu", "pct_cpu_p50", "pct_cpu_p95",
        "thread_efficiency_p50", "thread_efficiency_p95",
    ]
    summary[cols].to_csv(tsv_path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = FIGURES_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    par_df = collect_parallelism_rows(fetch=True)
    sweep_df = collect_pxd071075_sweep_rows()
    if not sweep_df.empty:
        # Audit TSV gets ALL sweep rows (5 points: q=10 through q=300)
        # so the F2d sweep is cross-referencable from F2a's data file.
        full_df = pd.concat([par_df, sweep_df], ignore_index=True, sort=False)
        write_parallelism_tsv(full_df, data_dir / "parallelism_data.tsv")
        # F2a render: keep only the q=300 row (the production-quality
        # cluster run) so PXD071075 contributes one dot, not a vertical
        # strip that crowds the scatter. The five-point sweep curve
        # itself is shown in F2d (queue_size_sweep.svg).
        sweep_top = sweep_df[sweep_df["queue_size"] == 300]
        par_df = pd.concat([par_df, sweep_top], ignore_index=True, sort=False)
    else:
        write_parallelism_tsv(par_df, data_dir / "parallelism_data.tsv")
    render_parallelism_scatter(
        par_df,
        FIGURES_DIR / "parallelism_vs_wallclock.svg",
    )

    durations, summary = collect_step_runtime_rows(fetch=True)
    write_per_step_tsv(summary, data_dir / "runtime_per_step.tsv")
    render_per_step_boxplot(
        durations, summary,
        FIGURES_DIR / "runtime_per_step.svg",
    )

    resources, resources_summary = collect_step_resource_rows(fetch=True)
    write_per_step_resources_tsv(
        resources_summary, data_dir / "resources_per_step.tsv",
    )
    render_resources_boxplot(
        resources, resources_summary,
        FIGURES_DIR / "resources_per_step.svg",
    )

    n_complete = int(par_df["trace_complete"].sum())
    print(
        f"Wrote {len(par_df)} parallelism rows "
        f"({n_complete} complete trace, "
        f"{len(par_df) - n_complete} hollow / pipeline_report.txt fallback) "
        f"and {len(summary)} step-summary rows."
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
