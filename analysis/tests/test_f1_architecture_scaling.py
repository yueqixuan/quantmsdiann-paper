"""Smoke tests for the Fig. 1 composite (architecture + scaling)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
import numpy as np


@pytest.fixture
def minimal_workflow_svg(tmp_path: Path) -> Path:
    svg = tmp_path / "workflow.svg"
    svg.write_text(
        '<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg" '
        'viewBox="0 0 100 100"><rect width="100" height="100" fill="#eee"/>'
        '<text x="10" y="50">workflow</text></svg>',
        encoding="utf-8",
    )
    return svg


def test_render_fig1_writes_svg(tmp_path: Path, minimal_workflow_svg: Path) -> None:
    from analysis.figure_f1_architecture_scaling import render_fig1_architecture_scaling

    sweep_df = pd.DataFrame([
        {"queue_size": 10, "wallclock_s": 3600.0, "peak_concurrent": 10, "n_tasks": 100},
        {"queue_size": 100, "wallclock_s": 900.0, "peak_concurrent": 100, "n_tasks": 100},
    ])
    par_df = pd.DataFrame([
        {
            "dataset": "PXD003539",
            "instrument": "TripleTOF 5600",
            "n_runs": 120,
            "wallclock_seconds": 7200.0,
            "trace_complete": True,
        },
        {
            "dataset": "PXD049412",
            "instrument": "Orbitrap Astral",
            "n_runs": 6,
            "wallclock_seconds": 1800.0,
            "trace_complete": True,
        },
    ])

    out = tmp_path / "fig1.svg"
    dummy = np.full((120, 80, 3), 240, dtype=np.uint8)
    with patch(
        "analysis.figure_f1_architecture_scaling.collect_sweep_rows",
        return_value=sweep_df,
    ), patch(
        "analysis.figure_f1_architecture_scaling._parallelism_plot_frame",
        return_value=par_df,
    ), patch(
        "analysis.figure_f1_architecture_scaling._workflow_raster",
        return_value=dummy,
    ):
        render_fig1_architecture_scaling(
            out,
            workflow_svg=minimal_workflow_svg,
            fetch=False,
        )
    assert out.exists() and out.stat().st_size > 1000
    text = out.read_text(encoding="utf-8")
    assert "<svg" in text[:500]
    # SVG only — no PDF is emitted by the script (the manuscript PDF is made
    # by the paper Makefile via rsvg-convert).
    assert not (tmp_path / "fig1.pdf").exists()
