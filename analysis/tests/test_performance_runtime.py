"""Tests for the runtime-vs-threads performance scatter.

Covered: Nextflow `duration: ...` parser, nextflow_trace.txt span aggregator,
DIA-NN command-line parser (threads + --f count + run names), instrument
heuristic, SDRF instrument extraction. No network access (fixtures only)."""
from __future__ import annotations

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Duration parser
# ---------------------------------------------------------------------------

def test_parse_duration_to_seconds_handles_full_hms() -> None:
    """Pipeline reports use compact 'Nh Nm Ns' tokens. 3h 26m 38s should
    round-trip to 12,398 seconds (3*3600 + 26*60 + 38)."""
    from analysis.figure_performance_runtime import parse_duration_to_seconds
    assert parse_duration_to_seconds("3h 26m 38s") == 3 * 3600 + 26 * 60 + 38


def test_parse_duration_to_seconds_handles_subsets() -> None:
    """Tokens may have any subset of h/m/s. nextflow_trace.txt rows often
    show only seconds, only minutes+seconds, or only milliseconds."""
    from analysis.figure_performance_runtime import parse_duration_to_seconds
    assert parse_duration_to_seconds("45s") == 45.0
    assert parse_duration_to_seconds("1m 9s") == 69.0
    assert parse_duration_to_seconds("2h") == 7200.0
    assert parse_duration_to_seconds("") == 0.0
    # Decimal seconds + ms are accepted (nextflow_trace.txt uses them).
    assert parse_duration_to_seconds("2.5s") == 2.5
    assert parse_duration_to_seconds("250ms") == 0.25


def test_parse_duration_rejects_garbage() -> None:
    from analysis.figure_performance_runtime import parse_duration_to_seconds
    with pytest.raises(ValueError):
        parse_duration_to_seconds("nopes")


def test_parse_pipeline_report_duration_finds_line(tmp_path: Path) -> None:
    """The wallclock parser scans the whole report; the line we care about
    isn't always at a known offset (different Nextflow versions emit it at
    different positions). Verify it picks up the canonical phrasing."""
    from analysis.figure_performance_runtime import parse_pipeline_report_duration
    body = (
        "----------------------------------------------------\n"
        "Run completion details\n"
        "----------------------------------------------------\n"
        "Completed: yes\n"
        "Exit status: 0\n"
        "The workflow was completed at 2026-05-17T23:44:10.5+01:00"
        " (duration: 3h 26m 38s)\n"
        "----------------------------------------------------\n"
    )
    p = tmp_path / "pipeline_report.txt"
    p.write_text(body)
    assert parse_pipeline_report_duration(p) == 3 * 3600 + 26 * 60 + 38


def test_parse_pipeline_report_duration_missing_raises(tmp_path: Path) -> None:
    from analysis.figure_performance_runtime import parse_pipeline_report_duration
    p = tmp_path / "broken.txt"
    p.write_text("Run completed; no duration recorded.\n")
    with pytest.raises(ValueError):
        parse_pipeline_report_duration(p)


# ---------------------------------------------------------------------------
# nextflow_trace.txt span aggregator
# ---------------------------------------------------------------------------

def test_aggregate_nextflow_trace_duration_uses_submit_plus_duration(
    tmp_path: Path,
) -> None:
    """For a trace with rows scattered in time, the wallclock is
    `max(submit + duration) - min(submit)`. Build a minimal 3-row trace
    where the last row finishes well after the others."""
    from analysis.figure_performance_runtime import (
        aggregate_nextflow_trace_duration,
    )
    body = (
        "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\t"
        "realtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"
        "1\ta\t1\tFOO\tCOMPLETED\t0\t2026-05-15 09:00:00.000\t30s\t30s\t"
        "100%\t1GB\t1GB\t1MB\t1MB\n"
        "2\tb\t2\tBAR\tCOMPLETED\t0\t2026-05-15 09:00:10.000\t2m\t2m\t"
        "100%\t1GB\t1GB\t1MB\t1MB\n"
        "3\tc\t3\tDIANN_RUN\tCOMPLETED\t0\t2026-05-15 09:05:00.000\t1h 30m\t"
        "1h 30m\t900%\t8GB\t8GB\t10GB\t1GB\n"
    )
    p = tmp_path / "nextflow_trace.txt"
    p.write_text(body)
    # First submit at 09:00:00; last task finishes 09:05:00 + 1h30m = 10:35:00.
    # Span = 1h 35m = 5700 s.
    assert aggregate_nextflow_trace_duration(p) == 5700.0


def test_aggregate_nextflow_trace_duration_empty_trace(tmp_path: Path) -> None:
    """Header-only traces (PXD070049/v2_3_2 publishes this) must not crash;
    we return 0 and let the caller fall back to pipeline_report.txt."""
    from analysis.figure_performance_runtime import (
        aggregate_nextflow_trace_duration,
    )
    p = tmp_path / "trace.txt"
    p.write_text(
        "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\t"
        "realtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"
    )
    assert aggregate_nextflow_trace_duration(p) == 0.0


# ---------------------------------------------------------------------------
# DIA-NN command-line parser
# ---------------------------------------------------------------------------

def test_parse_diann_command_extracts_threads_and_file_count(
    tmp_path: Path,
) -> None:
    """We need both `--threads N` and the count of `--f <path>` arguments,
    plus the basename list (used by the instrument heuristic)."""
    from analysis.figure_performance_runtime import parse_diann_command
    body = (
        "\nDIA-NN 2.5.0 Academia\n"
        "Compiled on Apr 12 2026 10:45:33\n"
        "Current date and time: Sun May 17 23:19:25 2026\n"
        "Logical CPU cores: 128\n"
        "diann --lib lib.parquet --fasta f.fasta "
        "--f run1.mzML --f run2.mzML --f run3.mzML "
        "--threads 24 --verbose 3 --matrices --out diann_report.tsv\n"
        "Thread number set to 24\n"
    )
    p = tmp_path / "log.log"
    p.write_text(body)
    threads, n_files, names = parse_diann_command(p)
    assert threads == 24
    assert n_files == 3
    assert names == ["run1.mzML", "run2.mzML", "run3.mzML"]


def test_parse_diann_command_missing_command_raises(tmp_path: Path) -> None:
    from analysis.figure_performance_runtime import parse_diann_command
    p = tmp_path / "bad.log"
    p.write_text("only header lines, no diann command\n")
    with pytest.raises(ValueError):
        parse_diann_command(p)


# ---------------------------------------------------------------------------
# Instrument inference
# ---------------------------------------------------------------------------

def test_infer_instrument_from_runs_recognises_known_patterns() -> None:
    """The DIA-NN run-name conventions used across the four ProteoBench
    modules are distinctive enough for a substring heuristic to disambiguate
    timsTOF SCP, Astral, and ZenoTOF runs."""
    from analysis.figure_performance_runtime import infer_instrument_from_runs
    assert infer_instrument_from_runs(
        ["LFQ_Astral_DIA_15min_50ng_Condition_A_REP1.raw"]
    ) == "Orbitrap Astral"
    assert infer_instrument_from_runs(
        ["ttSCP_diaPASEF_Condition_A_Sample_Alpha_01_11494.d"]
    ) == "timsTOF SCP"
    assert infer_instrument_from_runs(
        ["LFQ_ZenoTOF8600_ZenoSWATH_85VW_15min_Nano_50ng_A_REP1.mzML"]
    ) == "ZenoTOF 7600"
    # Unrecognised pattern (TripleTOF SWATH file from the Guo dataset).
    assert infer_instrument_from_runs(["guot_L130610_003_SW.mzML"]) is None
    assert infer_instrument_from_runs([]) is None


# ---------------------------------------------------------------------------
# SDRF instrument extraction
# ---------------------------------------------------------------------------

def test_parse_sdrf_instrument_strips_nt_field(tmp_path: Path) -> None:
    """SDRFs store instrument as 'NT=TripleTOF 5600;AC=MS:1000932'; we want
    the bare model name. When multiple instrument values exist across rows,
    return the most common one."""
    from analysis.figure_performance_runtime import parse_sdrf_instrument
    body = (
        "source name\tcomment[instrument]\n"
        "s1\tNT=TripleTOF 5600;AC=MS:1000932\n"
        "s2\tNT=TripleTOF 5600;AC=MS:1000932\n"
        "s3\tNT=TripleTOF 5600;AC=MS:1000932\n"
    )
    p = tmp_path / "x.sdrf.tsv"
    p.write_text(body)
    assert parse_sdrf_instrument(p) == "TripleTOF 5600"


# ---------------------------------------------------------------------------
# Trace-derived analyses: parallelism + per-step durations
# (figure_performance_trace.py)
# ---------------------------------------------------------------------------

def _make_trace(tmp_path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write a minimal nextflow_trace.txt with `(name, submit, duration)`
    tuples; remaining columns are filled with placeholders the parser
    ignores. Returns the path."""
    header = (
        "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\t"
        "realtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"
    )
    lines = [header]
    for i, (name, submit, duration) in enumerate(rows, start=1):
        lines.append(
            f"{i}\thash\t{i}\t{name}\tCOMPLETED\t0\t{submit}\t{duration}\t"
            f"{duration}\t100%\t1GB\t1GB\t1MB\t1MB\n"
        )
    p = tmp_path / "trace.txt"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("".join(lines))
    return p


def test_extract_step_name_handles_qualified_paths() -> None:
    """The trace stores fully-qualified `BIGBIO:...:STEP (input)` strings.
    The step parser must take only the last `:`-segment and strip the
    parenthesised argument. Rows without an argument (SUMMARY_PIPELINE)
    return the bare name."""
    from analysis.figure_performance_trace import extract_step_name
    assert extract_step_name(
        "BIGBIO_QUANTMSDIANN:QUANTMSDIANN:INPUT_CHECK:SAMPLESHEET_CHECK "
        "(PXD049412.sdrf.tsv)"
    ) == "SAMPLESHEET_CHECK"
    assert extract_step_name(
        "BIGBIO_QUANTMSDIANN:QUANTMSDIANN:DIA:PRELIMINARY_ANALYSIS "
        "(run_with_spaces_in_name)"
    ) == "PRELIMINARY_ANALYSIS"
    # No parenthesised argument.
    assert extract_step_name(
        "BIGBIO_QUANTMSDIANN:QUANTMSDIANN:SUMMARY_PIPELINE"
    ) == "SUMMARY_PIPELINE"
    # Defensive: an unqualified bare step name passes through.
    assert extract_step_name("STAND_ALONE") == "STAND_ALONE"
    assert extract_step_name("") == ""


def test_peak_concurrent_tasks_finds_max_overlap(tmp_path: Path) -> None:
    """Peak concurrent count is the largest number of intervals
    `[submit, submit+duration]` that overlap at any instant. Build a
    deliberately staggered 5-task fixture where the peak (3) only occurs in
    a narrow window in the middle of the trace."""
    from analysis.figure_performance_trace import (
        load_trace, peak_concurrent_tasks,
    )
    # Timeline:
    #   A: 09:00:00 -> 09:01:00 (60s)
    #   B: 09:00:30 -> 09:02:00 (90s)   -> A+B overlap 09:00:30-09:01:00
    #   C: 09:01:30 -> 09:03:00 (90s)   -> B+C overlap 09:01:30-09:02:00
    #   D: 09:01:40 -> 09:01:50 (10s)   -> peak=3 in 09:01:40-09:01:50
    #   E: 09:05:00 -> 09:05:10 (10s)   -> isolated
    p = _make_trace(tmp_path, [
        ("A:STEP (x)", "2026-05-15 09:00:00.000", "1m"),
        ("A:STEP (y)", "2026-05-15 09:00:30.000", "1m 30s"),
        ("A:STEP (z)", "2026-05-15 09:01:30.000", "1m 30s"),
        ("A:STEP (w)", "2026-05-15 09:01:40.000", "10s"),
        ("A:STEP (q)", "2026-05-15 09:05:00.000", "10s"),
    ])
    df = load_trace(p)
    peak, med = peak_concurrent_tasks(df)
    assert peak == 3
    # Median is over active-window samples (concurrency > 0); we just
    # require it to be in [1, peak]. A point estimate is brittle to
    # tie-break order, so the range check is the stable assertion.
    assert 1 <= med <= 3


def test_trace_wallclock_seconds_spans_first_submit_to_last_finish(
    tmp_path: Path,
) -> None:
    """The workflow wallclock is `max(submit+duration) - min(submit)`. With
    submits at 09:00:00 / 09:02:00 / 09:04:00 and the last task running for
    1h, the span is 09:00:00 -> 10:04:00 = 1h 4m = 3840 s."""
    from analysis.figure_performance_trace import (
        load_trace, trace_wallclock_seconds,
    )
    p = _make_trace(tmp_path, [
        ("A:FOO (a)", "2026-05-15 09:00:00.000", "30s"),
        ("A:FOO (b)", "2026-05-15 09:02:00.000", "30s"),
        ("A:BAR (c)", "2026-05-15 09:04:00.000", "1h"),
    ])
    df = load_trace(p)
    assert trace_wallclock_seconds(df) == 3840.0


def test_insilico_seconds_sums_only_library_step(tmp_path: Path) -> None:
    """`insilico_seconds` returns the total INSILICO_LIBRARY_GENERATION
    duration, ignoring every other step, and 0.0 when the step is absent."""
    from analysis.figure_performance_trace import insilico_seconds, load_trace
    p = _make_trace(tmp_path, [
        ("A:INDIVIDUAL_ANALYSIS (a)", "2026-05-15 09:00:00.000", "5m"),
        ("A:INSILICO_LIBRARY_GENERATION (lib)", "2026-05-15 09:05:00.000", "12m"),
        ("A:FINAL_QUANTIFICATION (q)", "2026-05-15 09:17:00.000", "8m"),
    ])
    df = load_trace(p)
    assert insilico_seconds(df) == 12 * 60.0

    p2 = _make_trace(tmp_path / "nolib", [
        ("A:INDIVIDUAL_ANALYSIS (a)", "2026-05-15 09:00:00.000", "5m"),
    ])
    assert insilico_seconds(load_trace(p2)) == 0.0


def test_busy_span_ignores_idle_gaps(tmp_path: Path) -> None:
    """`busy_span_seconds` is the union of active intervals, so an idle gap
    between two `-resume` legs is excluded (unlike the raw span). Two 10-min
    tasks an hour apart give 20 min busy, not 70 min span."""
    from analysis.figure_performance_trace import (
        busy_span_seconds, load_trace, trace_wallclock_seconds,
    )
    p = _make_trace(tmp_path, [
        ("A:STEP (leg1)", "2026-05-15 09:00:00.000", "10m"),
        ("A:STEP (leg2)", "2026-05-15 10:00:00.000", "10m"),
    ])
    df = load_trace(p)
    assert busy_span_seconds(df) == 20 * 60.0
    assert trace_wallclock_seconds(df) == 70 * 60.0
    # Overlapping tasks are merged (not double-counted): two 10-min tasks
    # starting 5 min apart cover a single 15-min window.
    p2 = _make_trace(tmp_path / "overlap", [
        ("A:STEP (x)", "2026-05-15 09:00:00.000", "10m"),
        ("A:STEP (y)", "2026-05-15 09:05:00.000", "10m"),
    ])
    assert busy_span_seconds(load_trace(p2)) == 15 * 60.0


def test_aggregate_step_durations_keys_by_step(tmp_path: Path) -> None:
    """`aggregate_step_durations` should collapse tasks across multiple
    traces into a {step -> [duration_s, ...]} map. Empty-step rows must be
    dropped (they correspond to header-only traces)."""
    from analysis.figure_performance_trace import (
        aggregate_step_durations, load_trace,
    )
    t1 = _make_trace(tmp_path / "t1", [
        ("A:DIANN_RUN (x)", "2026-05-15 09:00:00.000", "10m"),
        ("A:SAMPLESHEET_CHECK (s)", "2026-05-15 09:11:00.000", "30s"),
    ])
    t2 = _make_trace(tmp_path / "t2", [
        ("A:DIANN_RUN (y)", "2026-05-15 10:00:00.000", "8m"),
        ("A:DIANN_RUN (z)", "2026-05-15 10:00:00.000", "12m"),
    ])
    by_step = aggregate_step_durations([load_trace(t1), load_trace(t2)])
    assert sorted(by_step.keys()) == ["DIANN_RUN", "SAMPLESHEET_CHECK"]
    assert sorted(by_step["DIANN_RUN"]) == [480.0, 600.0, 720.0]
    assert by_step["SAMPLESHEET_CHECK"] == [30.0]


def test_load_report_window_data_extracts_trace_from_html(tmp_path: Path) -> None:
    """nextflow_report.html embeds the full per-task trace as a JavaScript
    object literal `window.data = {...};`. The parser brace-balances the
    blob, repairs Nextflow's invalid `\\'` escapes, and lifts the trace
    array into a DataFrame with submit/duration in seconds."""
    from analysis.figure_performance_trace import load_report_window_data
    # Minimal HTML mimicking Nextflow's report shape — single trace record,
    # epoch-ms timestamps, the over-escaped apostrophe Nextflow emits in
    # nested shell commands.
    html = (
        "<html><body><script>window.data = {"
        '"trace":[{"task_id":"1",'
        '"process":"BIGBIO:QUANTMSDIANN:DIA:PRELIMINARY_ANALYSIS",'
        '"name":"BIGBIO:QUANTMSDIANN:DIA:PRELIMINARY_ANALYSIS (raw1.mzML)",'
        '"status":"COMPLETED",'
        '"submit":"1700000000000","duration":"125000","realtime":"120000",'
        '"%cpu":"94.5","cmd":"awk -F \\\': \\\' \'{print $1}\'"}'
        "],"
        '"summary":{}};window.payload={};</script></body></html>'
    )
    p = tmp_path / "nextflow_report.html"
    p.write_text(html)
    df = load_report_window_data(p)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["step"] == "PRELIMINARY_ANALYSIS"
    assert row["status"] == "COMPLETED"
    assert row["submit"] == pytest.approx(1_700_000_000.0)
    assert row["duration_s"] == pytest.approx(125.0)
    assert row["finish"] == pytest.approx(1_700_000_125.0)


def test_load_report_window_data_returns_empty_on_missing_or_malformed(
    tmp_path: Path,
) -> None:
    """Parser must not crash on a missing file, a file without the
    `window.data` blob, or a malformed JSON payload."""
    from analysis.figure_performance_trace import load_report_window_data
    missing = tmp_path / "absent.html"
    assert load_report_window_data(missing).empty
    no_blob = tmp_path / "no_blob.html"
    no_blob.write_text("<html>no Nextflow data here</html>")
    assert load_report_window_data(no_blob).empty
    malformed = tmp_path / "bad.html"
    malformed.write_text(
        "<html><script>window.data = {\"trace\":[{not json</script></html>"
    )
    assert load_report_window_data(malformed).empty


def test_load_trace_handles_header_only_trace(tmp_path: Path) -> None:
    """Truncated traces (PXD070049/v2_3_2 ships header-only on PRIDE) must
    not crash the loader or downstream aggregators. `collect_parallelism_rows`
    drops them entirely rather than emitting sentinel rows."""
    from analysis.figure_performance_trace import (
        load_trace, peak_concurrent_tasks, trace_wallclock_seconds,
    )
    p = tmp_path / "trace.txt"
    p.write_text(
        "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\t"
        "realtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"
    )
    df = load_trace(p)
    assert df.empty
    assert peak_concurrent_tasks(df) == (0, 0.0)
    assert trace_wallclock_seconds(df) == 0.0


# ---------------------------------------------------------------------------
# F2c resource-column parsers and aggregator
# ---------------------------------------------------------------------------

def test_parse_size_to_bytes_covers_all_units() -> None:
    """Nextflow stores `peak_rss` as `<float> <unit>` where unit is one of
    B / KB / MB / GB / TB. The parser must round-trip each unit and return
    NaN for empty / `-` cells (Nextflow uses `-` for tasks that exit before
    a measurement is taken)."""
    import math
    from analysis.figure_performance_trace import parse_size_to_bytes
    assert parse_size_to_bytes("165 MB") == 165 * 1024 ** 2
    assert parse_size_to_bytes("4.7 GB") == 4.7 * 1024 ** 3
    assert parse_size_to_bytes("3.2 KB") == 3.2 * 1024
    assert parse_size_to_bytes("1024 B") == 1024.0
    assert math.isnan(parse_size_to_bytes(""))
    assert math.isnan(parse_size_to_bytes("-"))
    assert math.isnan(parse_size_to_bytes("nope"))


def test_parse_pct_cpu_handles_trailing_percent_and_excess() -> None:
    """`%cpu` is `<float>%`. Multi-thread tasks can exceed 100 % — the
    parser must not cap them or it would distort the per-step median for
    DIA-NN steps that run on multiple threads per task."""
    import math
    from analysis.figure_performance_trace import parse_pct_cpu
    assert parse_pct_cpu("94.8%") == 94.8
    assert parse_pct_cpu("1147.0%") == 1147.0
    assert parse_pct_cpu("0%") == 0.0
    assert math.isnan(parse_pct_cpu("-"))
    assert math.isnan(parse_pct_cpu(""))


def test_load_trace_resources_extracts_step_keyed_rows(tmp_path: Path) -> None:
    """Two-row fixture: a SAMPLESHEET_CHECK row and a FINAL_QUANTIFICATION
    row. Resources loader must extract both, attach the step name, and
    parse the human-readable size + percent strings into bytes / floats."""
    from analysis.figure_performance_trace import load_trace_resources
    p = tmp_path / "trace.txt"
    p.write_text(
        "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\t"
        "realtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"
        "1\ta1/b\t1\tPIPE:INPUT_CHECK:SAMPLESHEET_CHECK (foo.tsv)\t"
        "COMPLETED\t0\t2026-05-15 22:17:15.431\t1m 9s\t2s\t94.8%\t"
        "165 MB\t576.2 MB\t41.4 MB\t3.2 KB\n"
        "2\ta1/c\t2\tPIPE:DIANN:FINAL_QUANTIFICATION\tCOMPLETED\t0\t"
        "2026-05-15 22:18:00.000\t2h\t1h 55m\t1147%\t34.5 GB\t40 GB\t"
        "12 GB\t2 GB\n"
    )
    df = load_trace_resources(p)
    assert list(df["step"]) == ["SAMPLESHEET_CHECK", "FINAL_QUANTIFICATION"]
    assert df.iloc[0]["peak_rss_bytes"] == 165 * 1024 ** 2
    assert df.iloc[1]["peak_rss_bytes"] == 34.5 * 1024 ** 3
    assert df.iloc[0]["pct_cpu"] == 94.8
    assert df.iloc[1]["pct_cpu"] == 1147.0


def test_aggregate_step_resources_groups_and_filters_failed() -> None:
    """The aggregator must (1) group rows by step across multiple traces,
    (2) drop FAILED rows under `completed_only=True`, (3) keep COMPLETED
    rows from every trace, (4) carry NaN cells through without polluting
    the metric lists (a row missing RSS but present CPU contributes only
    CPU)."""
    import math
    import pandas as pd
    from analysis.figure_performance_trace import aggregate_step_resources
    trace1 = pd.DataFrame({
        "step": ["FINAL_QUANTIFICATION", "FINAL_QUANTIFICATION"],
        "status": ["COMPLETED", "FAILED"],
        "peak_rss_bytes": [34.5 * 1024 ** 3, 5.0 * 1024 ** 3],
        "pct_cpu": [1147.0, 50.0],
        "duration_s": [7200.0, 30.0],
    })
    trace2 = pd.DataFrame({
        "step": ["FINAL_QUANTIFICATION", "SAMPLESHEET_CHECK"],
        "status": ["COMPLETED", "COMPLETED"],
        "peak_rss_bytes": [40.0 * 1024 ** 3, math.nan],
        "pct_cpu": [1200.0, 90.0],
        "duration_s": [7500.0, 60.0],
    })
    out = aggregate_step_resources([trace1, trace2])
    assert sorted(out.keys()) == ["FINAL_QUANTIFICATION", "SAMPLESHEET_CHECK"]
    assert out["FINAL_QUANTIFICATION"]["peak_rss_bytes"] == [
        34.5 * 1024 ** 3, 40.0 * 1024 ** 3,
    ]
    assert out["FINAL_QUANTIFICATION"]["pct_cpu"] == [1147.0, 1200.0]
    # SAMPLESHEET_CHECK has CPU only; the NaN RSS must not appear.
    assert out["SAMPLESHEET_CHECK"]["peak_rss_bytes"] == []
    assert out["SAMPLESHEET_CHECK"]["pct_cpu"] == [90.0]
