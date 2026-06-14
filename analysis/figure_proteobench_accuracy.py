#!/usr/bin/env python
"""Fig 2 - ProteoBench accuracy: quantmsdiann within the single-machine community.

Restricted to the two datasets that HAVE public ProteoBench community
submissions to compare against (ProteoBench_Module_7 Astral, PXD062685
timsTOF/diaPASEF). The two datasets without predicted-library community
comparators (PXD049412 single-cell, PXD070049 ZenoTOF) are NOT shown here;
they go to the supplement (figure_id_vs_epsilon.py).

Two accuracy styles, two columns (one per dataset):
  Row A  Accuracy per concentration -- the HYE mix is a set of fixed multiples
         (Human 1x -> log2 0, Yeast 2x -> +1, E. coli 1/4x -> -2). Observed
         log2 fold-change per species: community submissions as a grey strip,
         the expected ratio as a dashed line, quantmsdiann versions as coloured
         markers. Shows quantmsdiann lands on the expected ratio within the
         community spread, at every concentration.
  Row B  Overall accuracy -- median |epsilon| (distance to expected, lower is
         better): community as a box+strip, quantmsdiann versions as markers.
         The equivalence headline: distributed quantmsdiann sits inside the
         single-machine community.

All community per-species / global accuracy comes from the cached ProteoBench
submission JSONs (data/quantmsdiann_benchmarks/proteobench/<dataset>.json),
which carry median_log2_empirical_<SP> and median_abs_epsilon_global per
submission. quantmsdiann from the metrics cache via the shared extractors.

Run:  python -m analysis.figure_proteobench_accuracy
Out:  analysis/figures/quantmsdiann_benchmarks/main_accuracy.svg
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()
from analysis.figure_id_vs_epsilon import (
    DATASET_TO_MODULE,
    DIANN_VERSIONS,
    METRICS_CACHE_DIR,
    SPECIES_EXPECTED_LOG2_A_vs_B,
    _COMMUNITY_COMPARATOR_DATASETS,
    _SPECIES_LABEL,
    _VERSION_COLORS,
    _VERSION_LABELS,
    _dataset_display_label,
    extract_qm_per_species_log2,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
COMMUNITY_DIR = REPO_ROOT / "data" / "quantmsdiann_benchmarks" / "proteobench"
OUT = REPO_ROOT / "analysis" / "figures" / "manuscript" / "fig2_proteobench_equivalence.svg"

SPECIES = ("HUMAN", "YEAST", "ECOLI")
COMMUNITY_COLOR = fs.OKABE_ITO["grey"]
# Per Vadim: showcase 1.8.1 vs 2.5.1 Enterprise only (the state-of-the-art
# build); the free standard 2.5.1 is dropped from the comparison.
VERSIONS = ("v1_8_1", "v2_5_1_enterprise")


def _community_entries(dataset: str, threshold: int) -> list[dict]:
    path = COMMUNITY_DIR / f"{dataset}.json"
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as fh:
        payload = json.load(fh)
    out = []
    for entry in payload:
        res = entry.get("results", {}).get(str(threshold))
        if isinstance(res, dict):
            out.append(res)
    return out


def _qm_ids_eps(dataset: str, threshold: int) -> list[tuple[str, float, float]]:
    """(version, nr_prec, median_abs_epsilon_global) per quantms.io version,
    in DIANN_VERSIONS order so the trajectory line reads 1.8.1 -> 2.5.1 -> ent."""
    rows = []
    for ver in VERSIONS:
        p = METRICS_CACHE_DIR / f"{dataset}_{ver}.json"
        if not p.exists():
            continue
        res = json.load(open(p, encoding="utf-8")).get("results", {}).get(str(threshold), {})
        n, e = res.get("nr_prec"), res.get("median_abs_epsilon_global")
        if n is not None and e is not None:
            rows.append((ver, float(n), float(e)))
    return rows


def _per_species(ax, dataset, threshold):
    module = DATASET_TO_MODULE[dataset]
    expected = SPECIES_EXPECTED_LOG2_A_vs_B.get(module, {})
    comm = _community_entries(dataset, threshold)
    qm = extract_qm_per_species_log2(dataset, threshold)
    qm = qm[qm["version"].isin(VERSIONS)]
    present = [s for s in SPECIES if s in expected]
    rng = np.random.default_rng(0)
    for x, sp in enumerate(present):
        ax.hlines(expected[sp], x - 0.34, x + 0.34, color="#444444",
                  linestyle="--", linewidth=1.0, zorder=1)
        cvals = [r.get(f"median_log2_empirical_{sp}") for r in comm]
        cvals = [v for v in cvals if isinstance(v, (int, float)) and not pd.isna(v)]
        if cvals:
            jit = rng.uniform(-0.13, 0.13, len(cvals))
            ax.scatter([x - 0.18 + j for j in jit], cvals, s=22,
                       color=COMMUNITY_COLOR, alpha=0.55, edgecolors="none", zorder=2)
        sub = qm[qm["species"] == sp]
        for _, row in sub.iterrows():
            ax.scatter(x + 0.18, row["mean_log2_empirical"], s=70,
                       color=_VERSION_COLORS.get(row["version"], "#d62728"),
                       edgecolors="white", linewidths=0.7, zorder=3)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([f"{_SPECIES_LABEL[s]}\n(log2={expected[s]:+.2g})" for s in present])
    ax.set_ylabel("Observed log2 fold-change")
    fs.despine(ax)


def _ids_vs_eps(ax, dataset, threshold):
    """Depth vs accuracy: precursors (x) against median |eps| (y). Community
    submissions as a grey cloud; quantms.io versions as a coloured trajectory
    (1.8.1 -> 2.5.1 -> Enterprise) — moving RIGHT (more IDs) at near-constant
    height (accuracy), so the depth-for-accuracy tradeoff reads correctly."""
    comm = _community_entries(dataset, threshold)
    pts = [(r.get("nr_prec"), r.get("median_abs_epsilon_global")) for r in comm]
    pts = [(n, e) for n, e in pts
           if isinstance(n, (int, float)) and isinstance(e, (int, float))
           and not pd.isna(n) and not pd.isna(e)]
    if pts:
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=26,
                   color=COMMUNITY_COLOR, alpha=0.55, edgecolors="none", zorder=2)
    qm = _qm_ids_eps(dataset, threshold)
    if qm:
        ax.plot([q[1] for q in qm], [q[2] for q in qm], color="#999999",
                linewidth=1.0, zorder=2)
        for ver, n, e in qm:
            ax.scatter(n, e, s=85, color=_VERSION_COLORS.get(ver, "#d62728"),
                       edgecolors="white", linewidths=0.8, zorder=3)
    fs.kfmt_axis(ax.xaxis)
    ax.set_xlabel("Precursors quantified (≥3 rep)")
    ax.set_ylabel("Median |ε|  (lower = more accurate)")
    fs.despine(ax)


def render(out: Path, threshold: int = 3) -> Path:
    datasets = list(_COMMUNITY_COMPARATOR_DATASETS)
    fig, axes = plt.subplots(2, len(datasets), figsize=(4.6 * len(datasets), 7.6))
    for j, ds in enumerate(datasets):
        title = _dataset_display_label(ds).splitlines()[0]
        _per_species(axes[0, j], ds, threshold)
        axes[0, j].set_title(title)
        _ids_vs_eps(axes[1, j], ds, threshold)

    # version legend + community
    handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=8,
               markerfacecolor=COMMUNITY_COLOR, markeredgecolor="none",
               label="single-machine community (ProteoBench)")]
    for ver in VERSIONS:
        if any((METRICS_CACHE_DIR / f"{ds}_{ver}.json").exists() for ds in datasets):
            handles.append(Line2D([0], [0], marker="o", linestyle="none", markersize=8,
                           markerfacecolor=_VERSION_COLORS.get(ver, "#d62728"),
                           markeredgecolor="white",
                           label=f"quantms.io DIA-NN {_VERSION_LABELS.get(ver, ver)}"))
    fig.legend(handles=handles, loc="upper center", ncol=min(4, len(handles)),
               bbox_to_anchor=(0.5, 1.01), fontsize=8)

    fig.text(0.01, 0.71, "A  accuracy per concentration", rotation=90,
             va="center", fontsize=11, fontweight="bold")
    fig.text(0.01, 0.28, "B  depth vs accuracy", rotation=90,
             va="center", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0.03, 0, 1, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
