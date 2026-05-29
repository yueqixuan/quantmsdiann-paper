"""Tests for the Fig 2 accuracy figure (analysis/figure_id_vs_epsilon.py).

Covers the two things that can silently go wrong in the redesigned panel:
1. the per-species expected-ratio reference lines must match ProteoBench's
   authoritative `species_expected_ratio()` (not a stale hand-entered value);
2. `render_accuracy_panels` must write a non-empty SVG even with an empty
   metrics cache (the offline / CI path), so a missing cache degrades to an
   empty panel rather than a crash.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest


def test_per_species_expected_ratios_match_proteobench() -> None:
    """SPECIES_EXPECTED_LOG2_A_vs_B must equal log2 of ProteoBench's
    species_expected_ratio() A_vs_B for every benchmarked module."""
    pytest.importorskip("proteobench")
    from proteobench.io.parsing.parse_settings import ParseSettingsBuilder

    from analysis.figure_id_vs_epsilon import SPECIES_EXPECTED_LOG2_A_vs_B
    from analysis.proteobench_metrics import (
        DATASET_TO_MODULE,
        _proteobench_parse_settings_dir,
    )

    for module_id in set(DATASET_TO_MODULE.values()):
        parser = ParseSettingsBuilder(
            parse_settings_dir=str(_proteobench_parse_settings_dir(module_id)),
            module_id=module_id,
        ).build_parser("DIA-NN")
        authoritative = parser.species_expected_ratio()
        coded = SPECIES_EXPECTED_LOG2_A_vs_B[module_id]
        assert set(coded) == set(authoritative), module_id
        for species, meta in authoritative.items():
            expected_log2 = math.log2(meta["A_vs_B"])
            assert abs(coded[species] - expected_log2) < 1e-6, (
                module_id, species, coded[species], expected_log2,
            )


def test_render_accuracy_panels_writes_svg(tmp_path: Path) -> None:
    """Renders (b)+(c) and writes a non-empty SVG. Resilient to an empty
    metrics cache: per-species panels still draw the expected-ratio lines
    and the community panels simply carry no points."""
    from analysis.figure_id_vs_epsilon import render_accuracy_panels

    svg = tmp_path / "main_accuracy.svg"
    df = render_accuracy_panels(threshold=3, svg_path=svg)
    assert svg.exists() and svg.stat().st_size > 0
    # The audit frame is a DataFrame (possibly empty when the cache is cold).
    assert hasattr(df, "columns")
