#!/usr/bin/env python
"""PXD064049 (CHP-212 MYCN Deep Visual Proteomics, diaPASEF) reanalysis:
quantmsdiann (DIA-NN 2.5.0, library-free) versus the originally deposited
DIA-NN 1.8.1 analysis on the identical 12 DVP runs.

The original analysis (PRIDE PXD064049) used DIA-NN 1.8.1 library-free with a
plain human FASTA; quantmsdiann re-ran the same raw files with DIA-NN 2.5.0
against an entrapment+contaminant-augmented FASTA, which enforces an
empirically-validated FDR. We therefore compare:

  * main_comparison.svg -- precursors and protein groups at 1% FDR. Precursors
                           (least sensitive to the FASTA choice) are at parity
                           (~99.5%); the quantmsdiann protein-group count is
                           lower because its search used an entrapment-augmented
                           FASTA. counts.tsv also records the entrapment hit
                           rate (fraction of accepted ids mapping to entrapment
                           sequences) as a measure of error control.

Run:  PYTHONPATH=. python -m analysis.figure_pxd064049_spatial_vs_quantmsdiann
"""
from __future__ import annotations

import io
import sys
import urllib.request
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis.contaminant_filter import (
    is_target_protein_group,
    count_target_precursors,
    count_target_protein_groups,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "PXD064049"
CACHE_DIR = FIGURES_DIR / "data" / "cache"

ORIG_COLOUR = "#9e9e9e"
QM_COLOUR = "#1e88e5"

_QB = ("https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
       "quantmsdiann-benchmarks/PXD064049-MYCN-DVP-diaPASEF/quant_tables")
_ORIG_ZIP = ("https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/"
             "PXD064049/DIANN_results.zip")


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or dest.stat().st_size == 0:
        with urllib.request.urlopen(url, timeout=600) as r:
            dest.write_bytes(r.read())
    return dest


def _qm_matrix(kind: str) -> Path:
    """quantmsdiann pr/pg matrix from the benchmarks FTP (cached)."""
    name = f"diann_report.{kind}_matrix.tsv"
    return _download(f"{_QB}/{name}", CACHE_DIR / f"qm_{kind}_matrix.tsv")


def _orig_matrix(kind: str) -> Path:
    """Authors' deposited DIA-NN 1.8.1 pr/pg matrix (cached from the zip)."""
    dest = CACHE_DIR / f"orig_{kind}_matrix.tsv"
    if not dest.exists() or dest.stat().st_size == 0:
        zip_dest = _download(_ORIG_ZIP, CACHE_DIR / "DIANN_results.zip")
        with zipfile.ZipFile(zip_dest) as z:
            member = next(m for m in z.namelist()
                          if m.endswith(f"MYCN_High_Low.{kind}_matrix.tsv"))
            dest.write_bytes(z.read(member))
    return dest


def _entrapment_hit_rate(matrix_path: Path) -> tuple[int, int, float]:
    """(entrapment_passing, target_passing, entrapment_hit_rate_pct): the
    fraction of accepted identifications whose Protein.Group maps to an
    entrapment sequence. This is a direct measure of how many accepted
    groups are entrapment hits; it equals the empirical FDR only when the
    entrapment database is target-sized (1:1 paired entrapment), so we
    report it as an entrapment hit rate rather than a calibrated FDR."""
    pgs = pd.read_csv(matrix_path, sep="\t", usecols=["Protein.Group"],
                      dtype=str)["Protein.Group"].dropna()
    entrap = int(pgs.str.contains("ENTRAP_").sum())
    target = int(pgs.map(is_target_protein_group).sum())
    return entrap, target, (100.0 * entrap / target if target else 0.0)


def render_main_comparison(or_pr: int, qm_pr: int, or_pg: int, qm_pg: int,
                           svg_path: Path) -> None:
    """Main Fig.~3 panel (d): 2-condition x 2-metric grouped bar chart,
    original (DIA-NN 1.8.1, grey) vs quantmsdiann (DIA-NN 2.5.0, blue), for
    precursors and protein groups at 1% FDR. Matches the per-cohort
    `main_comparison` style of the other panels (log y if a metric's
    cross-condition spread exceeds 5x)."""
    conditions = [
        ("Original (DIA-NN 1.8.1)", ORIG_COLOUR, or_pr, or_pg),
        ("quantmsdiann (DIA-NN 2.5.0)", QM_COLOUR, qm_pr, qm_pg),
    ]
    metrics = ["Precursors", "Protein groups"]
    bar_width = 0.27
    x = [0, 1]
    offsets = [bar_width * (i - (len(conditions) - 1) / 2.0)
               for i in range(len(conditions))]

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, (label, color, pr_val, pg_val) in enumerate(conditions):
        values = [pr_val, pg_val]
        bars = ax.bar([xi + offsets[i] for xi in x], values, width=bar_width,
                      color=color, edgecolor="#37474f", label=label)
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height(),
                    f"{val:,}", ha="center", va="bottom", fontsize=9)

    needs_log = any(
        min(v) > 0 and max(v) / min(v) > 5
        for v in ([or_pr, qm_pr], [or_pg, qm_pg])
    )
    ylabel = "Count (1% FDR)"
    if needs_log:
        ax.set_yscale("log")
        ylabel += " (log scale)"
    else:
        ax.set_ylim(0, max(or_pr, qm_pr, or_pg, qm_pg) * 1.18)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    _save(fig, svg_path)


def _save(fig, svg_path: Path) -> None:
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    stem = svg_path.with_suffix("")
    for ext in (".svg", ".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    qm_pr_t = count_target_precursors(_qm_matrix("pr"))[1]
    qm_pg_t = count_target_protein_groups(_qm_matrix("pg"))[1]
    or_pr_t = count_target_precursors(_orig_matrix("pr"))[1]
    or_pg_t = count_target_protein_groups(_orig_matrix("pg"))[1]
    pr_entrap, _, pr_hit = _entrapment_hit_rate(_qm_matrix("pr"))
    pg_entrap, _, pg_hit = _entrapment_hit_rate(_qm_matrix("pg"))

    render_main_comparison(or_pr_t, qm_pr_t, or_pg_t, qm_pg_t,
                           FIGURES_DIR / "main_comparison.svg")

    counts = FIGURES_DIR / "counts.tsv"
    counts.write_text(
        "metric\toriginal_diann181\tquantmsdiann_diann250\t"
        "qm_entrapment_hits\tqm_entrapment_hit_pct\n"
        f"precursors\t{or_pr_t}\t{qm_pr_t}\t{pr_entrap}\t{pr_hit:.3f}\n"
        f"protein_groups\t{or_pg_t}\t{qm_pg_t}\t{pg_entrap}\t{pg_hit:.3f}\n"
    )
    print(f"precursors: original={or_pr_t}  quantmsdiann={qm_pr_t}  "
          f"(entrapment hit rate {pr_hit:.2f}%, {pr_entrap} hits)")
    print(f"protein groups: original={or_pg_t}  quantmsdiann={qm_pg_t}  "
          f"(entrapment hit rate {pg_hit:.2f}%, {pg_entrap} hits)")
    print(f"wrote {FIGURES_DIR}/main_comparison.svg + supp_protein_groups.svg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
