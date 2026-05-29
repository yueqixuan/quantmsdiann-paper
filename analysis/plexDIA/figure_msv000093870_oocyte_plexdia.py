"""MSV000093870 plexDIA reanalysis — first plexDIA cohort through quantmsdiann.

Single-cell oocyte plexDIA dataset (Slavov-lab deposition on MassIVE),
reanalysed with quantmsdiann using DIA-NN 2.5.0 in mTRAQ 3-channel mode
(channels 0 / 4 / 8). Each of the 38 Q-Exactive raw files multiplexes three
single cells, one per mTRAQ channel, for up to 114 single-cell proteomes.

This is the first analysis of the plexDIA branch. It characterises the
per-channel, per-single-cell identification depth — the metric that matters
for plexDIA, since the headline run-level numbers conflate three cells.

Inputs are pulled from the public quantmsdiann-benchmarks deposition and
cached on disk:
  https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/
    quantmsdiann-benchmarks/plexDIA/MSV000093870-plexDIA/

Outputs (paper-ready, no titles/footers):
  analysis/figures/plexDIA/MSV000093870/
    main_plexdia_per_cell.{svg,pdf,png}  — proteins & precursors per single cell, by channel
    counts.tsv                           — auditable per-cell and per-channel totals
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from analysis.contaminant_filter import is_target_protein_group

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "plexDIA" / "MSV000093870"
CACHE_DIR = FIGURES_DIR / "data" / "cache"

FTP_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/plexDIA/MSV000093870-plexDIA"
)
REPORT_PARQUET_URL = f"{FTP_BASE}/quant_tables/diann_report.parquet"

# mTRAQ channel display colours (light -> dark, matching the labelling reagent
# delta masses 0 / 4 / 8 Da).
CHANNEL_COLOURS = {"0": "#90caf9", "4": "#1e88e5", "8": "#0d47a1"}
CHANNEL_LABELS = {"0": "mTRAQ-0", "4": "mTRAQ-4", "8": "mTRAQ-8"}

# Columns needed for channel-confident identification counting.
REPORT_COLUMNS = [
    "Run", "Channel", "Protein.Group", "Protein.Ids", "Precursor.Id", "Decoy",
    "Q.Value", "Channel.Q.Value", "Translated.Q.Value", "Precursor.Quantity",
]


def _cached_report() -> Path:
    """Download the DIA-NN report parquet once and cache it on disk."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / "diann_report.parquet"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"Downloading {REPORT_PARQUET_URL} (~245 MB, cached)…", file=sys.stderr)
    with urllib.request.urlopen(REPORT_PARQUET_URL, timeout=600) as resp:
        dest.write_bytes(resp.read())
    return dest


def load_channel_confident(report_path: Path) -> pd.DataFrame:
    """Channel-confident, target-only precursor rows.

    This reproduces how DIA-NN builds its per-channel protein/precursor
    expression matrices in plexDIA mode: a precursor is counted in a channel
    when the run-level ``Q.Value`` and the channel-specific ``Channel.Q.Value``
    both clear 1% FDR, with a positive channel quantity. Target-only under the
    shared conservative contaminant/decoy filter.

    Note: we deliberately do NOT additionally gate on ``Translated.Q.Value``.
    That q-value scores the MBR/translation step specifically; requiring it for
    every precursor (including ones detected directly in the channel) is
    over-conservative for identification counting and depresses per-cell
    counts by ~7% without a corresponding FDR justification. ``Channel.Q.Value``
    is the threshold DIA-NN itself applies when assembling the channel
    matrices (``--matrix-qvalue 0.01``).
    """
    df = pq.read_table(report_path, columns=REPORT_COLUMNS).to_pandas()
    df["Channel"] = df["Channel"].astype(str)
    confident = df[
        (df["Decoy"] == 0)
        & (df["Q.Value"] <= 0.01)
        & (df["Channel.Q.Value"] <= 0.01)
        & (df["Precursor.Quantity"] > 0)
    ].copy()
    confident = confident[confident["Protein.Group"].map(is_target_protein_group)]
    return confident


def per_cell_counts(confident: pd.DataFrame) -> pd.DataFrame:
    """One row per single cell = (Run, Channel): precursor and protein counts."""
    return (
        confident.groupby(["Run", "Channel"])
        .agg(
            precursors=("Precursor.Id", "nunique"),
            proteins=("Protein.Group", "nunique"),
        )
        .reset_index()
        .sort_values(["Channel", "Run"])
    )


def render_per_cell_figure(cells: pd.DataFrame, svg_path: Path) -> None:
    """Two-panel boxplot: proteins and precursors per single cell, by channel."""
    channels = ["0", "4", "8"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4))

    for ax, metric, ylabel in (
        (axes[0], "proteins", "Protein groups per single cell (1% FDR)"),
        (axes[1], "precursors", "Precursors per single cell (1% FDR)"),
    ):
        data = [cells.loc[cells["Channel"] == ch, metric].values for ch in channels]
        bp = ax.boxplot(
            data, widths=0.6, patch_artist=True, showfliers=False,
            medianprops=dict(color="#212121", linewidth=1.4),
        )
        for patch, ch in zip(bp["boxes"], channels):
            patch.set_facecolor(CHANNEL_COLOURS[ch])
            patch.set_alpha(0.85)
            patch.set_edgecolor("#37474f")
        # jittered per-cell points
        rng = np.random.default_rng(0)
        for i, ch in enumerate(channels, start=1):
            y = cells.loc[cells["Channel"] == ch, metric].values
            x = rng.normal(i, 0.05, size=len(y))
            ax.scatter(x, y, s=10, c="#37474f", alpha=0.45, linewidths=0, zorder=3)
        ax.set_xticks(range(1, len(channels) + 1))
        ax.set_xticklabels([CHANNEL_LABELS[ch] for ch in channels])
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)

    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    stem = svg_path.with_suffix("")
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_counts(cells: pd.DataFrame, confident: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("metric\tchannel\tvalue\tsource\n")
        fh.write(f"single_cells_total\tall\t{len(cells)}\tRun x Channel (channel-confident)\n")
        fh.write(f"runs\tall\t{cells['Run'].nunique()}\tdiann_report.parquet\n")
        fh.write(f"unique_protein_groups\tall\t{confident['Protein.Group'].nunique()}\ttarget, channel-confident\n")
        fh.write(f"unique_precursors\tall\t{confident['Precursor.Id'].nunique()}\ttarget, channel-confident\n")
        fh.write(f"median_proteins_per_cell\tall\t{int(cells['proteins'].median())}\tmedian across cells\n")
        fh.write(f"median_precursors_per_cell\tall\t{int(cells['precursors'].median())}\tmedian across cells\n")
        for ch in ["0", "4", "8"]:
            sub = cells[cells["Channel"] == ch]
            fh.write(f"n_cells\t{CHANNEL_LABELS[ch]}\t{len(sub)}\tcells in channel\n")
            fh.write(f"median_proteins_per_cell\t{CHANNEL_LABELS[ch]}\t{int(sub['proteins'].median())}\tmedian\n")
            fh.write(f"median_precursors_per_cell\t{CHANNEL_LABELS[ch]}\t{int(sub['precursors'].median())}\tmedian\n")


def main() -> int:  # pragma: no cover
    report = _cached_report()
    confident = load_channel_confident(report)
    cells = per_cell_counts(confident)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    render_per_cell_figure(cells, FIGURES_DIR / "main_plexdia_per_cell.svg")
    write_counts(cells, confident, FIGURES_DIR / "counts.tsv")

    print(f"Single cells (run x channel): {len(cells)} across {cells['Run'].nunique()} runs")
    print(f"Median proteins/cell: {int(cells['proteins'].median())}, "
          f"median precursors/cell: {int(cells['precursors'].median())}")
    print(f"Unique target protein groups: {confident['Protein.Group'].nunique()}, "
          f"precursors: {confident['Precursor.Id'].nunique()}")
    print(f"Wrote {FIGURES_DIR / 'main_plexdia_per_cell.svg'}")
    print(f"Wrote {FIGURES_DIR / 'counts.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
