#!/usr/bin/env python
"""Fig 1d - ProteoBench quantification-accuracy concordance vs standalone DIA-NN.

A single concordance panel: for the two ProteoBench DIA modules with public
predicted-from-FASTA **DIA-NN** community submissions (Module 7, Orbitrap
Astral; PXD062685, timsTOF diaPASEF), the measured HYE log2 fold-change (Y) is
plotted against the ProteoBench-expected ratio (X), with the dashed Y=X line.
Standalone DIA-NN community runs (ALL versions) are grey; quantms-diann (1.8.1
and 2.5.1-enterprise) are coloured. Points on the diagonal => quantms-diann
quantifies as accurately as single-machine DIA-NN, independent of release.

Only DIA-NN community submissions are used (other ProteoBench tools excluded),
so "standalone DIA-NN" is the literal comparator. Per-version, per-module
identification DEPTH (precursors / protein groups) is in Supplementary Note 5
(figure_id_vs_epsilon.py); this panel is accuracy only.

Source: data/quantmsdiann_benchmarks/proteobench/<dataset>.json (community) and
the quantms-diann metrics cache (per version).

Run:  python -m analysis.figure_proteobench_accuracy
Out:  analysis/figures/manuscript/fig1d_proteobench_accuracy.svg
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from analysis import figure_style as fs
fs.apply_house_style()
from analysis.figure_id_vs_epsilon import (
    DATASET_TO_MODULE, SPECIES_EXPECTED_LOG2_A_vs_B, _COMMUNITY_COMPARATOR_DATASETS,
    _VERSION_COLORS, _VERSION_LABELS, _SPECIES_LABEL, METRICS_CACHE_DIR,
    extract_qm_per_species_log2,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
COMMUNITY_DIR = REPO_ROOT / "data" / "quantmsdiann_benchmarks" / "proteobench"
OUT = REPO_ROOT / "analysis" / "figures" / "manuscript" / "fig1d_proteobench_accuracy.svg"

SPECIES = ("HUMAN", "YEAST", "ECOLI")
GREY = fs.OKABE_ITO["grey"]
VERSIONS = ("v1_8_1", "v2_5_1_enterprise")
MARKER = {"ProteoBench_Module_7": "o", "PXD062685": "s"}
MARKER_LABEL = {"ProteoBench_Module_7": "Astral (Module 7)", "PXD062685": "timsTOF diaPASEF"}


def diann_community(dataset: str, threshold: int = 3) -> list[dict]:
    """Standalone DIA-NN, predicted-from-FASTA submissions (ALL versions).

    Only DIA-NN entries with DIA-NN-predicted libraries are kept, so other
    ProteoBench tools (AlphaDIA, Spectronaut, ...) and user-supplied-library
    runs are excluded — the grey cloud is literally standalone DIA-NN.
    """
    payload = json.load(open(COMMUNITY_DIR / f"{dataset}.json", encoding="utf-8"))
    out = []
    for e in payload:
        if str(e.get("software_name", "")).strip().upper() != "DIA-NN":
            continue
        pred = e.get("predictors_library")
        if not (isinstance(pred, dict) and str(pred.get("RT", "")).upper() == "DIANN"):
            continue
        r = e.get("results", {}).get(str(threshold))
        if isinstance(r, dict):
            out.append(r)
    return out


def draw(ax, threshold: int = 3, *, with_legend: bool = True, compact: bool = False,
         square: bool = True) -> None:
    """Draw the accuracy concordance into `ax` (reused by the Fig 2 row)."""
    datasets = list(_COMMUNITY_COMPARATOR_DATASETS)
    rng = np.random.default_rng(0)
    lab = 9 if compact else None
    for ds in datasets:
        module = DATASET_TO_MODULE[ds]
        expected = SPECIES_EXPECTED_LOG2_A_vs_B.get(module, {})
        comm = diann_community(ds, threshold)
        qm = extract_qm_per_species_log2(ds, threshold)
        qm = qm[qm["version"].isin(VERSIONS)]
        for sp in SPECIES:
            if sp not in expected:
                continue
            x = expected[sp]
            cv = [c.get(f"median_log2_empirical_{sp}") for c in comm]
            cv = [v for v in cv if isinstance(v, (int, float)) and not np.isnan(v)]
            ax.scatter(x + rng.uniform(-0.05, 0.05, len(cv)), cv, s=22, color=GREY,
                       alpha=0.5, edgecolors="none", zorder=2)
            for _, r in qm[qm["species"] == sp].iterrows():
                ax.scatter(x, r["median_log2_empirical"], s=80, marker=MARKER[ds],
                           color=_VERSION_COLORS.get(r["version"], "#d62728"),
                           edgecolors="white", linewidths=0.8, zorder=3)
    lim = [-2.7, 1.7]
    ax.plot(lim, lim, "--", color="#444444", linewidth=1.1, zorder=1)
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    if square:
        ax.set_aspect("equal")
    ax.set_xlabel("Expected log$_2$ ratio (ProteoBench)", fontsize=lab)
    ax.set_ylabel("Observed log$_2$ ratio", fontsize=lab)
    if compact:
        ax.tick_params(labelsize=8)
    fs.despine(ax)

    if with_legend:
        # short labels, placed in the empty upper-left triangle (points sit on
        # the diagonal), so the legend never overlaps the data.
        handles = [Line2D([0], [0], linestyle="--", color="#444444", label="Y = X"),
                   Line2D([0], [0], marker="o", linestyle="none", ms=7, markerfacecolor=GREY,
                          markeredgecolor="none", label="standalone DIA-NN")]
        for ver in VERSIONS:
            handles.append(Line2D([0], [0], marker="o", linestyle="none", ms=7,
                           markerfacecolor=_VERSION_COLORS.get(ver, "#d62728"), markeredgecolor="white",
                           label=f"quantms-diann {_VERSION_LABELS.get(ver, ver)}"))
        ax.legend(handles=handles, loc="upper left", fontsize=6.5 if compact else 8,
                  frameon=False, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.3)


SHAPE_BY_VERSION = {"v1_8_1": "o", "v2_5_1_enterprise": "^"}


def draw_strip(ax, threshold: int = 3, *, compact: bool = False,
               dataset_colors: dict | None = None) -> None:
    """Per-species accuracy strip+box: one group per HYE species, community runs
    as a jittered grey strip + box (every dot visible), the ProteoBench-expected
    ratio as a dashed line, quantms-diann as large markers. Shows the same
    accuracy/equivalence as the concordance plot but with the community dots
    clearly spread out instead of stacked on the diagonal."""
    datasets = list(_COMMUNITY_COMPARATOR_DATASETS)
    exp: dict[str, float] = {}
    for ds in datasets:
        exp.update(SPECIES_EXPECTED_LOG2_A_vs_B.get(DATASET_TO_MODULE[ds], {}))
    order = sorted([s for s in SPECIES if s in exp], key=lambda s: exp[s])  # ECOLI,-HUMAN,YEAST
    lab = 8 if compact else None
    rng = np.random.default_rng(0)
    for x, sp in enumerate(order):
        ax.hlines(exp[sp], x - 0.42, x + 0.42, color="#444444", linestyle="--",
                  linewidth=1.2, zorder=1)
        cv = []
        for ds in datasets:
            cv += [c.get(f"median_log2_empirical_{sp}") for c in diann_community(ds, threshold)]
        cv = [v for v in cv if isinstance(v, (int, float)) and not np.isnan(v)]
        if cv:
            bp = ax.boxplot([cv], positions=[x], widths=0.62, showfliers=False,
                            patch_artist=True, zorder=2)
            for b in bp["boxes"]:
                b.set(facecolor="#eeeeee", edgecolor="#9e9e9e", linewidth=0.8)
            for w in bp["whiskers"] + bp["caps"]:
                w.set(color="#9e9e9e", linewidth=0.8)
            for med in bp["medians"]:
                med.set(color="#9e9e9e", linewidth=1.0)
            ax.scatter(x + rng.uniform(-0.24, 0.24, len(cv)), cv, s=16, color=GREY,
                       alpha=0.6, edgecolors="white", linewidths=0.3, zorder=3)
        for di, ds in enumerate(datasets):
            qm = extract_qm_per_species_log2(ds, threshold)
            qm = qm[qm["species"] == sp]
            qm = qm[qm["version"].isin(VERSIONS)]
            # colour = instrument/dataset (defined once in the figure's main
            # instrument legend); shape = quantms-diann version.
            colour = (dataset_colors or {}).get(ds)
            ds_dx = -0.05 if di == 0 else 0.05
            for _, r in qm.iterrows():
                ver = r["version"]
                ver_dx = -0.13 if ver == VERSIONS[0] else 0.13
                ax.scatter(x + ver_dx + ds_dx, r["median_log2_empirical"], s=34,
                           marker=SHAPE_BY_VERSION.get(ver, "o"),
                           color=colour or _VERSION_COLORS.get(ver, "#d62728"),
                           edgecolors="black", linewidths=0.6, zorder=4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{_SPECIES_LABEL[s]}\n(exp. {exp[s]:+g})" for s in order], fontsize=lab)
    ax.set_ylabel("Observed log$_2$ ratio", fontsize=lab)
    ax.set_xlim(-0.6, len(order) - 0.4)
    if compact:
        ax.tick_params(labelsize=8)
    fs.despine(ax)
    # Minimal key: colours (= instrument/dataset) are read from the figure's
    # main instrument legend; here we only define the line, the community dots,
    # and the SHAPE = quantms-diann version.
    handles = [Line2D([0], [0], linestyle="--", color="#444444", label="ProteoBench expected"),
               Line2D([0], [0], marker="o", linestyle="none", ms=7, markerfacecolor=GREY,
                      markeredgecolor="white", label="standalone DIA-NN")]
    for ver in VERSIONS:
        handles.append(Line2D([0], [0], marker=SHAPE_BY_VERSION.get(ver, "o"), linestyle="none",
                       ms=7, markerfacecolor="#666666", markeredgecolor="black",
                       label=f"quantms-diann {_VERSION_LABELS.get(ver, ver)}"))
    ax.legend(handles=handles, loc="upper left", fontsize=6.5 if compact else 8,
              frameon=False, handletextpad=0.3, labelspacing=0.3, ncol=1)


def render(out: Path, threshold: int = 3) -> Path:
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    draw(ax, threshold)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
