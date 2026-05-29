"""F3-version — cell-line cross-version progression.

Scaffold for experiment #13 in [docs/brainstorming.md]. The actual data
collection (3 cell-line cohorts × ≥2 DIA-NN versions = 6+ reruns on the
cluster, with v1.8.1 as the original `quantms` baseline and v2.5.0 as
the current default) happens on the cluster. This file holds the
**consumer** code that loads matrices once they're staged under
`data/<PXD>/<version>/diann_report.pr_matrix.tsv` and
`data/<PXD>/<version>/diann_report.pg_matrix.tsv`, then renders a
side-by-side ID-progression figure per cohort.

When run with the data missing the script prints a "no input yet"
message and exits cleanly so CI / the figure suite never blocks on
experiment #13.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "combined"

# Cell-line cohorts and the DIA-NN versions we plan to sweep. Endpoint
# subset (v1_8_1 + v2_5_0) per §9.2 of the brainstorming doc; the
# intermediate three versions are optional / supplementary.
CELL_LINE_COHORTS = ("PXD003539", "PXD030304", "PXD004701")
VERSION_SWEEP = ("v1_8_1", "v2_1_0", "v2_2_0", "v2_3_2", "v2_5_0")
VERSION_LABELS = {
    "v1_8_1": "1.8.1",
    "v2_1_0": "2.1.0",
    "v2_2_0": "2.2.0",
    "v2_3_2": "2.3.2",
    "v2_5_0": "2.5.0",
}


def cohort_matrix_path(cohort: str, version: str, kind: str) -> Path:
    """Local path to a per-cohort per-version matrix file. `kind` is
    either `"pr"` (precursors) or `"pg"` (protein groups)."""
    return DATA_DIR / cohort / version / f"diann_report.{kind}_matrix.tsv"


def count_matrix_rows(matrix_path: Path) -> int:
    """Row count of a DIA-NN matrix file, minus the header. Returns 0
    when the file is missing."""
    if not matrix_path.exists():
        return 0
    with open(matrix_path, encoding="utf-8") as fh:
        return max(0, sum(1 for _ in fh) - 1)


def iter_present_cells(
    cohorts: Iterable[str] = CELL_LINE_COHORTS,
    versions: Iterable[str] = VERSION_SWEEP,
) -> Iterable[tuple[str, str]]:
    """Yield (cohort, version) for every cohort/version pair whose
    pr_matrix.tsv exists locally. The consumer logic uses this to
    decide what to render — partial sweeps (e.g. v1_8_1 endpoint
    only) are supported."""
    for cohort in cohorts:
        for version in versions:
            if cohort_matrix_path(cohort, version, "pr").exists():
                yield cohort, version


def collect_progression_rows(
    cohorts: Iterable[str] = CELL_LINE_COHORTS,
    versions: Iterable[str] = VERSION_SWEEP,
) -> pd.DataFrame:
    """One row per (cohort, version) pair that has matrix data, with
    `n_precursors` and `n_proteins`. Returns an empty DataFrame with
    that schema when no sweep data is present."""
    rows: list[dict] = []
    for cohort, version in iter_present_cells(cohorts, versions):
        rows.append({
            "cohort": cohort,
            "version": version,
            "version_label": VERSION_LABELS.get(version, version),
            "n_precursors": count_matrix_rows(
                cohort_matrix_path(cohort, version, "pr"),
            ),
            "n_proteins": count_matrix_rows(
                cohort_matrix_path(cohort, version, "pg"),
            ),
        })
    return pd.DataFrame(
        rows,
        columns=[
            "cohort", "version", "version_label",
            "n_precursors", "n_proteins",
        ],
    )


def render_progression(df: pd.DataFrame, svg_path: Path) -> None:
    """Two-panel grouped bar chart (precursors / proteins) with
    DIA-NN version on x-axis, one bar group per cohort.

    When `df` is empty or has data for only one cohort/version pair,
    emits an "experiment #13 data-bound" annotation instead of a
    misleading half-figure."""
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(11.0, 4.2))
    cohorts_present = sorted(df["cohort"].unique()) if len(df) else []
    if not cohorts_present or len(df) < 2:
        for ax in axes:
            ax.text(
                0.5, 0.5,
                "cell-line cross-version sweep is data-bound on "
                "experiment #13.\nRender will populate once "
                "`data/<PXD>/v*_*/diann_report.pr_matrix.tsv` exists.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color="#888888",
            )
            ax.set_axis_off()
        fig.tight_layout()
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
        return

    versions_present = [
        v for v in VERSION_SWEEP if v in df["version"].unique()
    ]
    cohort_colours = {
        "PXD003539": "#1b9e77",
        "PXD030304": "#7570b3",
        "PXD004701": "#d95f02",
    }
    for ax, metric, ylabel in (
        (axes[0], "n_precursors", "Precursors"),
        (axes[1], "n_proteins",  "Protein groups"),
    ):
        n_groups = len(versions_present)
        n_cohorts = len(cohorts_present)
        width = 0.8 / max(1, n_cohorts)
        x_base = list(range(n_groups))
        for i, cohort in enumerate(cohorts_present):
            sub = df[df["cohort"] == cohort].set_index("version")
            heights = [
                sub.loc[v, metric] if v in sub.index else 0
                for v in versions_present
            ]
            offsets = [x + (i - (n_cohorts - 1) / 2) * width for x in x_base]
            ax.bar(
                offsets, heights, width=width,
                color=cohort_colours.get(cohort, "#9e9e9e"),
                label=cohort,
            )
        ax.set_xticks(x_base)
        ax.set_xticklabels(
            [VERSION_LABELS.get(v, v) for v in versions_present],
            fontsize=8,
        )
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_xlabel("DIA-NN version", fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=8)
    axes[0].legend(loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def write_progression_tsv(df: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep="\t", index=False)


def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = FIGURES_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    df = collect_progression_rows()
    svg_path = FIGURES_DIR / "cell_line_version_progression.svg"
    tsv_path = data_dir / "cell_line_version_progression.tsv"
    if df.empty:
        # No data yet — do NOT ship an empty SVG / TSV. Remove any
        # stale artefacts from a previous run and exit cleanly. The
        # consumer is data-bound on experiment #13 cluster reruns.
        if svg_path.exists():
            svg_path.unlink()
        if tsv_path.exists():
            tsv_path.unlink()
        print(
            "F3-version (cell-line cross-version) — no data yet under "
            f"{DATA_DIR}/<PXD>/v*_*/diann_report.pr_matrix.tsv. "
            "No figure rendered; cohort reruns happen on the cluster "
            "(see docs/superpowers/specs/2026-05-20-experiment-13-cell-line-version-sweep.md)."
        )
        return 0
    write_progression_tsv(df, tsv_path)
    render_progression(df, svg_path)
    cohorts = sorted(df["cohort"].unique())
    print(
        f"F3-version rendered from {len(df)} (cohort, version) "
        f"pair(s): {cohorts}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
