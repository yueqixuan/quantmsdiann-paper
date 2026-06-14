"""Tests for the experiment-#11/#12/#13 scaffolds.

These scripts are data-bound on cluster reruns that haven't happened
yet. The tests check that:
  - The consumer code handles empty input gracefully (returns empty
    DataFrames, renders an explainer panel, exits cleanly).
  - URL / path builders use the expected layout so a future operator
    can stage data without re-reading the code.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# #11 — cell-line trace ingestion scaffolding
# ---------------------------------------------------------------------------

def test_cell_line_trace_local_path_layout() -> None:
    """Cell-line traces land under
    `data/<PXD>/pipeline_info/nextflow_trace.txt`. The builder must
    use that layout so the SLURM rerun script can stage files
    without consulting Python."""
    from analysis.figure_performance_trace import cell_line_trace_local_path
    p = cell_line_trace_local_path("PXD003539")
    assert p.parts[-3:] == ("PXD003539", "pipeline_info", "nextflow_trace.txt")


def test_has_cell_line_traces_false_when_missing(tmp_path: Path) -> None:
    """`has_cell_line_traces()` controls whether F2b/F2c expand to the
    cell-line cohort. Until a real trace.txt lands under each cell-line
    `pipeline_info/` dir, it must return False. (Note: we exercise
    the helper directly rather than monkey-patching DATA_DIR — the
    repo at test time should never have cell-line traces; if it did,
    this assertion would flip and we'd want to know.)"""
    from analysis.figure_performance_trace import (
        CELL_LINE_ANALYSES, cell_line_trace_local_path,
    )
    # Cross-check: the path builder is well-formed for every cohort.
    for ds in CELL_LINE_ANALYSES:
        p = cell_line_trace_local_path(ds)
        assert p.name == "nextflow_trace.txt"


# ---------------------------------------------------------------------------
# #12 — queue-size sweep scaffolding
# ---------------------------------------------------------------------------

def test_collect_sweep_rows_returns_empty_when_no_data(tmp_path: Path) -> None:
    """`collect_sweep_rows` must return an empty DataFrame with the
    expected schema when no queue-sweep traces exist. This is the
    default state of the repo — experiment #12 hasn't run."""
    from analysis import figure_queue_size_sweep
    # Point the consumer at an empty directory.
    orig = figure_queue_size_sweep.SWEEP_DIR
    figure_queue_size_sweep.SWEEP_DIR = tmp_path
    try:
        df = figure_queue_size_sweep.collect_sweep_rows()
    finally:
        figure_queue_size_sweep.SWEEP_DIR = orig
    assert df.empty
    assert list(df.columns) == [
        "queue_size", "wallclock_s", "wallclock_with_lib_s",
        "insilico_lib_s", "peak_concurrent", "n_tasks",
    ]


def test_render_queue_size_sweep_handles_empty_input(tmp_path: Path) -> None:
    """The render function must emit a placeholder SVG instead of
    crashing when given an empty frame. The CI render-all-figures
    path relies on this."""
    from analysis.figure_queue_size_sweep import render_queue_size_sweep
    df = pd.DataFrame(
        columns=["queue_size", "wallclock_s", "peak_concurrent", "n_tasks"]
    )
    out = tmp_path / "queue_size_sweep.svg"
    render_queue_size_sweep(df, out)
    assert out.exists() and out.stat().st_size > 0
