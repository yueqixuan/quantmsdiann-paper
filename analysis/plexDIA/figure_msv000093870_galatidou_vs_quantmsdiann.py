"""MSV000093870 plexDIA: quantmsdiann reanalysis vs Galatidou et al. 2024.

The original analysis (Galatidou, Petelski et al., Mol Hum Reprod 2024,
doi:10.1093/molehr/gaae023; code at github.com/SlavovLab/single_cell_oocyte)
quantified single human oocytes by mTRAQ plexDIA (mPOP sample prep, Q
Exactive) and published a post-QC proteins x oocyte matrix
(2022_07_28_Oocyte_ProteinsXcells.csv): relative protein abundance,
mean-normalised across successfully quantified oocytes.

This script compares that published matrix against the quantmsdiann
reanalysis (DIA-NN 2.5.0, mTRAQ 3-channel, channel-confident + conservative
contaminant filter; see figure_msv000093870_oocyte_plexdia.py) on three axes:
protein groups per single cell, total protein groups, and protein-accession
overlap.

Inputs are cached on disk:
  - quantmsdiann report parquet (PRIDE quantmsdiann-benchmarks deposition)
  - original matrix (SlavovLab/single_cell_oocyte GitHub)

Outputs (paper-ready, no titles/footers):
  analysis/figures/plexDIA/MSV000093870/
    main_galatidou_comparison.{svg,pdf,png}
    comparison_counts.tsv
"""
from __future__ import annotations

import re
import sys
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.plexDIA.figure_msv000093870_oocyte_plexdia import (
    FIGURES_DIR,
    CACHE_DIR,
    _cached_report,
    load_channel_confident,
    per_cell_counts,
)

ORIGINAL_MATRIX_URL = (
    "https://raw.githubusercontent.com/SlavovLab/single_cell_oocyte/"
    "main/2022_07_28_Oocyte_ProteinsXcells.csv"
)

ORIG_COLOUR = "#9e9e9e"
QM_COLOUR = "#1e88e5"


def _strip_isoform(accession: str) -> str:
    return str(accession).split("-")[0]


def _cached_original() -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / "Galatidou_2024_ProteinsXcells.csv"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"Downloading original matrix {ORIGINAL_MATRIX_URL} (cached)…", file=sys.stderr)
    with urllib.request.urlopen(ORIGINAL_MATRIX_URL, timeout=120) as resp:
        dest.write_bytes(resp.read())
    return dest


def original_per_cell(matrix_path: Path) -> tuple[pd.Series, set[str]]:
    """Proteins quantified per oocyte and the set of accessions, from the
    published proteins x cell matrix (non-missing entries = quantified)."""
    o = pd.read_csv(matrix_path)
    cells = o.columns[1:]
    per_cell = o[cells].notna().sum(axis=0)
    accessions = {_strip_isoform(a) for a in o["leading.protein"].dropna().astype(str)}
    accessions.discard("")
    return per_cell, accessions


_RUN_PREFIX_RE = re.compile(r"(wAP\d+)")
_ORIG_COL_RE = re.compile(r"(wAP\d+)_d(\d)_")


def _run_prefix(run: str) -> str:
    m = _RUN_PREFIX_RE.match(str(run))
    return m.group(1) if m else str(run)


def qc_matched_table(matrix_path: Path, qm_cells: pd.DataFrame) -> pd.DataFrame:
    """Per-oocyte protein counts for the cells retained by BOTH analyses.

    Each oocyte is keyed by (run prefix, mTRAQ channel). The original matrix
    columns (e.g. ``wAP0021_d0_1A``) give the original count (non-missing
    entries); quantmsdiann gives its count for the same (run, channel) cell.
    Restricting to the intersection removes the cell-QC confound (the original
    dropped low-quality oocytes that quantmsdiann's channel filter retained).
    """
    o = pd.read_csv(matrix_path)
    orig = {}
    for col in o.columns[1:]:
        m = _ORIG_COL_RE.match(col)
        if m:
            orig[(m.group(1), m.group(2))] = int(o[col].notna().sum())

    qm = qm_cells.copy()
    qm["key"] = list(zip(qm["Run"].map(_run_prefix), qm["Channel"].astype(str)))
    rows = []
    for _, r in qm.iterrows():
        if r["key"] in orig:
            rows.append({
                "run": r["key"][0], "channel": r["key"][1],
                "orig_proteins": orig[r["key"]], "qm_proteins": int(r["proteins"]),
            })
    return pd.DataFrame(rows)


def render_qc_matched(matched: pd.DataFrame, svg_path: Path) -> None:
    """Paired per-oocyte comparison on the QC-matched cohort."""
    r = float(np.corrcoef(matched["orig_proteins"], matched["qm_proteins"])[0, 1])
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.3))

    # (a) paired scatter, y = x reference
    ax = axes[0]
    ax.scatter(matched["orig_proteins"], matched["qm_proteins"], s=22,
               c=QM_COLOUR, alpha=0.7, edgecolors="#37474f", linewidths=0.5)
    lo = 0
    hi = max(matched["orig_proteins"].max(), matched["qm_proteins"].max()) * 1.05
    ax.plot([lo, hi], [lo, hi], color="#9e9e9e", linestyle="--", linewidth=1)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Galatidou 2024 proteins / oocyte")
    ax.set_ylabel("quantmsdiann proteins / oocyte")
    ax.text(0.04, 0.96, f"(a)  n={len(matched)} matched oocytes\nPearson r = {r:.2f}",
            transform=ax.transAxes, va="top", fontsize=9)

    # (b) matched boxplots
    ax = axes[1]
    data = [matched["orig_proteins"].values, matched["qm_proteins"].values]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="#212121", linewidth=1.4))
    for patch, c in zip(bp["boxes"], (ORIG_COLOUR, QM_COLOUR)):
        patch.set_facecolor(c); patch.set_alpha(0.85); patch.set_edgecolor("#37474f")
    rng = np.random.default_rng(0)
    for i, vals in enumerate(data, start=1):
        ax.scatter(rng.normal(i, 0.05, size=len(vals)), vals, s=10,
                   c="#37474f", alpha=0.4, linewidths=0, zorder=3)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Galatidou 2024", "quantmsdiann"])
    ax.set_ylabel("Protein groups per oocyte (matched cohort)")
    ax.set_ylim(bottom=0)
    ax.text(0.5, 0.97, "(b)", transform=ax.transAxes, fontweight="bold", va="top")

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)

    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    stem = svg_path.with_suffix("")
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return r


def quantms_accessions(confident: pd.DataFrame) -> set[str]:
    acc: set[str] = set()
    for ids in confident["Protein.Ids"].dropna().unique():
        for a in str(ids).split(";"):
            acc.add(_strip_isoform(a))
    acc.discard("")
    return acc


def render_comparison(
    orig_per_cell: pd.Series,
    qm_cells: pd.DataFrame,
    orig_total: int,
    qm_total: int,
    shared: int,
    orig_only: int,
    qm_only: int,
    svg_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3))

    # (a) proteins per single cell — original vs quantmsdiann
    ax = axes[0]
    data = [orig_per_cell.values, qm_cells["proteins"].values]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="#212121", linewidth=1.4))
    for patch, c in zip(bp["boxes"], (ORIG_COLOUR, QM_COLOUR)):
        patch.set_facecolor(c); patch.set_alpha(0.85); patch.set_edgecolor("#37474f")
    rng = np.random.default_rng(0)
    for i, (vals, c) in enumerate(zip(data, (ORIG_COLOUR, QM_COLOUR)), start=1):
        ax.scatter(rng.normal(i, 0.05, size=len(vals)), vals, s=10,
                   c="#37474f", alpha=0.4, linewidths=0, zorder=3)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"Galatidou 2024\n(n={len(orig_per_cell)})",
                        f"quantmsdiann\n(n={len(qm_cells)})"])
    ax.set_ylabel("Protein groups per single cell")
    ax.set_ylim(bottom=0)
    ax.text(0.5, 0.97, "(a)", transform=ax.transAxes, fontweight="bold", va="top")

    # (b) total protein groups
    ax = axes[1]
    bars = ax.bar([0, 1], [orig_total, qm_total], width=0.6,
                  color=[ORIG_COLOUR, QM_COLOUR], edgecolor="#37474f")
    for b, v in zip(bars, (orig_total, qm_total)):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center",
                va="bottom", fontsize=10)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Galatidou 2024", "quantmsdiann"])
    ax.set_ylabel("Total protein groups (dataset)")
    ax.set_ylim(0, max(orig_total, qm_total) * 1.15)
    ax.text(0.5, 0.97, "(b)", transform=ax.transAxes, fontweight="bold", va="top")

    # (c) protein-group overlap (leading-protein level, like-for-like)
    ax = axes[2]
    cats = ["Shared", "Galatidou\nonly", "quantmsdiann\nonly"]
    vals = [shared, orig_only, qm_only]
    cols = ["#43a047", ORIG_COLOUR, QM_COLOUR]
    bars = ax.bar(range(3), vals, width=0.6, color=cols, edgecolor="#37474f")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center",
                va="bottom", fontsize=10)
    ax.set_xticks(range(3))
    ax.set_xticklabels(cats)
    ax.set_ylabel("Protein groups")
    ax.set_ylim(0, max(vals) * 1.15)
    ax.text(0.5, 0.97, "(c)", transform=ax.transAxes, fontweight="bold", va="top")

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)

    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    stem = svg_path.with_suffix("")
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)


def render_total_pg(orig_total: int, qm_total: int, svg_path: Path) -> None:
    """Compact single-panel total-protein-group comparison for the main
    cell-line-reanalysis figure (Fig 3): original (grey) vs quantmsdiann
    (blue), matching the per-cohort `main_comparison` bar style. Conveys
    the single-cell reanalysis benefit (more protein groups from the same
    raw data) in one column-width panel. Uses the same (7, 5) canvas as the
    per-cohort `main_comparison` panels so all three Fig. 3 sub-panels share
    an aspect ratio and render at equal height."""
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar([0, 1], [orig_total, qm_total], width=0.55,
                  color=[ORIG_COLOUR, QM_COLOUR], edgecolor="#37474f")
    ax.set_xlim(-0.7, 1.7)
    for b, v in zip(bars, (orig_total, qm_total)):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center",
                va="bottom", fontsize=11)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Galatidou 2024\n(original)", "quantmsdiann\n(DIA-NN)"])
    ax.set_ylabel("Protein groups (dataset, 1\\% FDR)".replace("\\%", "%"))
    ax.set_ylim(0, max(orig_total, qm_total) * 1.16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    stem = svg_path.with_suffix("")
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:  # pragma: no cover
    confident = load_channel_confident(_cached_report())
    qm_cells = per_cell_counts(confident)
    qm_total = confident["Protein.Group"].nunique()
    qm_acc = quantms_accessions(confident)
    # Leading protein (group representative) per quantmsdiann protein group, so
    # the overlap is computed against the original's leading.protein set at the
    # SAME granularity (group/leading-protein level) rather than mixing the
    # original's leading proteins with quantmsdiann's fully-expanded accessions.
    qm_lead = {_strip_isoform(str(g).split(";")[0])
               for g in confident["Protein.Group"].dropna().unique()}
    qm_lead.discard("")

    matrix_path = _cached_original()
    orig_per_cell, orig_acc = original_per_cell(matrix_path)
    orig_total = len(orig_acc)

    # Like-for-like overlap at the protein-group / leading-protein level.
    shared = len(orig_acc & qm_lead)
    orig_only = len(orig_acc - qm_lead)
    qm_only = len(qm_lead - orig_acc)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    render_comparison(
        orig_per_cell, qm_cells, orig_total=len(orig_acc), qm_total=qm_total,
        shared=shared, orig_only=orig_only, qm_only=qm_only,
        svg_path=FIGURES_DIR / "main_galatidou_comparison.svg",
    )

    # Compact total-PG panel for the main cell-line-reanalysis figure (Fig 3).
    render_total_pg(
        orig_total=len(orig_acc), qm_total=qm_total,
        svg_path=FIGURES_DIR / "main_galatidou_total_pg.svg",
    )

    # QC-matched, apples-to-apples per-oocyte comparison
    matched = qc_matched_table(matrix_path, qm_cells)
    r = render_qc_matched(matched, FIGURES_DIR / "main_galatidou_qc_matched.svg")

    counts = FIGURES_DIR / "comparison_counts.tsv"
    with counts.open("w") as fh:
        fh.write("metric\tGalatidou_2024\tquantmsdiann\n")
        fh.write(f"cells\t{len(orig_per_cell)}\t{len(qm_cells)}\n")
        fh.write(f"median_proteins_per_cell\t{int(orig_per_cell.median())}\t{int(qm_cells['proteins'].median())}\n")
        fh.write(f"protein_groups_total\t{len(orig_acc)}\t{qm_total}\n")
        fh.write(f"protein_groups_shared\t{shared}\t{shared}\n")
        fh.write(f"protein_groups_unique\t{orig_only}\t{qm_only}\n")
        fh.write(f"quantmsdiann_expanded_accessions\t-\t{len(qm_acc)}\n")
        fh.write(f"qc_matched_cells\t{len(matched)}\t{len(matched)}\n")
        fh.write(f"qc_matched_median_proteins_per_cell\t{int(matched['orig_proteins'].median())}\t{int(matched['qm_proteins'].median())}\n")
        fh.write(f"qc_matched_pearson_r\t{r:.3f}\t{r:.3f}\n")

    print("=== Galatidou 2024 vs quantmsdiann ===")
    print(f"cells:            {len(orig_per_cell)} vs {len(qm_cells)}")
    print(f"median prot/cell: {int(orig_per_cell.median())} vs {int(qm_cells['proteins'].median())}")
    print(f"protein groups:   {len(orig_acc)} vs {qm_total}")
    print(f"protein groups:   shared {shared}, Galatidou-only {orig_only}, quantmsdiann-only {qm_only} "
          f"({100*shared/len(orig_acc):.0f}% of original recovered)")
    print(f"--- QC-matched cohort ({len(matched)} oocytes) ---")
    print(f"median prot/cell: {int(matched['orig_proteins'].median())} (orig) vs "
          f"{int(matched['qm_proteins'].median())} (quantmsdiann), Pearson r={r:.2f}")
    print(f"Wrote {FIGURES_DIR / 'main_galatidou_comparison.svg'}")
    print(f"Wrote {FIGURES_DIR / 'main_galatidou_qc_matched.svg'}")
    print(f"Wrote {counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
