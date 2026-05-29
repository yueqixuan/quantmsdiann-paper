"""F1b — Identifications vs accuracy scatter.

Per dataset, one panel: x = `nr_prec` at the chosen replicate threshold
(default ≥3), y = `median_abs_epsilon_global` at the same threshold.

- Community ProteoBench submissions: grey cloud (one dot per submission),
  faceted by `predictors_library` palette.
- quantmsdiann (this work): red dots, one per DIA-NN version, joined by a
  thin line in version order so the reader can read the progression.
- Per-dataset metric source: locally-cached `proteobench_metrics/*.json`
  populated by `analysis/proteobench_metrics.py`. If the cache is empty,
  the figure renders an explainer panel rather than crashing.

The accompanying F1c per-species log2 strip plot lives in a sibling
function `render_per_species_log2` so the supp variant ships from the
same script.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
    DATA_DIR as BENCHMARKS_DATA_DIR,
    DIANN_VERSIONS,
    LIBRARY_KIND_EMPIRICAL,
    LIBRARY_KIND_OTHER_TOOL,
    LIBRARY_KIND_PREDICTED,
    LIBRARY_KIND_USER_DEFINED,
    _VERSION_LABELS,
    _dataset_display_label,
    _dataset_sort_key,
    classify_predictors_library,
)
from analysis.proteobench_metrics import (
    DATASET_TO_MODULE,
    METRICS_CACHE_DIR,
    cached_proteobench_metrics,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "quantmsdiann_benchmarks"


LIB_PALETTE = {
    LIBRARY_KIND_EMPIRICAL: "#80cbc4",
    LIBRARY_KIND_PREDICTED: "#1976d2",
    LIBRARY_KIND_USER_DEFINED: "#9575cd",
    LIBRARY_KIND_OTHER_TOOL: "#bdbdbd",
}


# ---------------------------------------------------------------------------
# Community-cohort data extraction
# ---------------------------------------------------------------------------

def extract_community_id_vs_eps(
    dataset: str, threshold: int,
) -> pd.DataFrame:
    """Load the cached ProteoBench submission JSON for a dataset and
    pull `(nr_prec, median_abs_epsilon_global, library_kind, software_name)`
    per submission at the given replicate threshold. Submissions
    missing either metric at this threshold are dropped."""
    cache = BENCHMARKS_DATA_DIR / "proteobench" / f"{dataset}.json"
    if not cache.exists():
        return pd.DataFrame(columns=[
            "nr_prec", "median_abs_epsilon_global",
            "library_kind", "software_name",
        ])
    with open(cache, encoding="utf-8") as fh:
        entries = json.load(fh)
    rows: list[dict] = []
    for entry in entries:
        results = entry.get("results", {})
        thr_key = str(threshold)
        if thr_key not in results:
            continue
        res = results[thr_key]
        nr_prec = res.get("nr_prec")
        eps = res.get("median_abs_epsilon_global")
        if nr_prec is None or eps is None:
            continue
        try:
            nr_prec = int(nr_prec)
            eps = float(eps)
        except (TypeError, ValueError):
            continue
        if pd.isna(eps):
            continue
        rows.append({
            "nr_prec": nr_prec,
            "median_abs_epsilon_global": eps,
            "library_kind": classify_predictors_library(
                entry.get("predictors_library"),
            ),
            "software_name": entry.get("software_name", ""),
        })
    return pd.DataFrame(rows)


def extract_quantmsdiann_id_vs_eps(
    dataset: str, threshold: int,
) -> pd.DataFrame:
    """Pull `(version, nr_prec, median_abs_epsilon_global)` per
    quantmsdiann DIA-NN version from the locally-cached metrics. Rows
    where the metric cache hasn't been computed yet are silently
    skipped — the figure renders only what's available."""
    rows: list[dict] = []
    for version in DIANN_VERSIONS:
        cache_path = METRICS_CACHE_DIR / f"{dataset}_{version}.json"
        if not cache_path.exists():
            continue
        with open(cache_path, encoding="utf-8") as fh:
            payload = json.load(fh)
        res = payload.get("results", {}).get(str(threshold))
        if not res:
            continue
        nr_prec = res.get("nr_prec")
        eps = res.get("median_abs_epsilon_global")
        if nr_prec is None or eps is None or pd.isna(eps):
            continue
        rows.append({
            "version": version,
            "version_label": _VERSION_LABELS.get(version, version),
            "nr_prec": int(nr_prec),
            "median_abs_epsilon_global": float(eps),
        })
    df = pd.DataFrame(rows)
    if len(df):
        df = df.assign(
            order=lambda d: d["version"].map(
                {v: i for i, v in enumerate(DIANN_VERSIONS)}
            )
        ).sort_values("order").drop(columns="order").reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# F1b — main figure
# ---------------------------------------------------------------------------

def render_id_vs_epsilon(
    threshold: int,
    svg_path: Path,
    *,
    datasets: Iterable[str] | None = None,
) -> pd.DataFrame:
    """4-panel scatter (one per dataset). Returns a long-format
    DataFrame of every plotted point for the auditable TSV."""
    datasets = (
        sorted(datasets, key=_dataset_sort_key)
        if datasets is not None
        else sorted(DATASET_TO_MODULE, key=_dataset_sort_key)
    )
    fig, axes = plt.subplots(
        nrows=2, ncols=2, figsize=(11.0, 8.2),
        sharex=False, sharey=False, squeeze=False,
    )
    long_rows: list[dict] = []
    qm_drawn = False
    libs_drawn: set[str] = set()
    for idx, dataset in enumerate(datasets):
        ax = axes[idx // 2][idx % 2]
        community = extract_community_id_vs_eps(dataset, threshold)
        # Apples-to-apples: quantmsdiann predicts its library in-silico from
        # the FASTA (DIA-NN library-free), so compare only against ProteoBench
        # DIA-NN submissions that used the same predicted (DIANN) strategy.
        # Modules with no predicted-library submission (single-cell, ZenoTOF)
        # therefore show the quantmsdiann trajectory alone.
        community = community[
            community["library_kind"] == LIBRARY_KIND_PREDICTED
        ].reset_index(drop=True)
        qm = extract_quantmsdiann_id_vs_eps(dataset, threshold)
        # Community cloud
        for kind in (
            LIBRARY_KIND_EMPIRICAL, LIBRARY_KIND_USER_DEFINED,
            LIBRARY_KIND_OTHER_TOOL, LIBRARY_KIND_PREDICTED,
        ):
            sub = community[community["library_kind"] == kind]
            if not len(sub):
                continue
            ax.scatter(
                sub["nr_prec"], sub["median_abs_epsilon_global"],
                s=40, c=LIB_PALETTE[kind], alpha=0.7,
                edgecolors="#555555", linewidths=0.4,
                label=kind if kind not in libs_drawn else None,
            )
            libs_drawn.add(kind)
        for _, row in community.iterrows():
            long_rows.append({
                "dataset": dataset, "threshold": threshold,
                "source": "proteobench-community",
                "label": row["software_name"],
                "library_kind": row["library_kind"],
                "nr_prec": int(row["nr_prec"]),
                "median_abs_epsilon_global": float(
                    row["median_abs_epsilon_global"]
                ),
            })
        # quantmsdiann trajectory
        if len(qm):
            ax.plot(
                qm["nr_prec"], qm["median_abs_epsilon_global"],
                color="#d62728", linewidth=1.0, alpha=0.7, zorder=2,
            )
            ax.scatter(
                qm["nr_prec"], qm["median_abs_epsilon_global"],
                s=80, c="#d62728", edgecolors="#7f1d1d",
                linewidths=0.8, zorder=3,
                label=(
                    "quantmsdiann (DIA-NN, predicted-from-FASTA lib)"
                    if not qm_drawn else None
                ),
            )
            qm_drawn = True
            for _, row in qm.iterrows():
                ax.annotate(
                    row["version_label"],
                    xy=(row["nr_prec"], row["median_abs_epsilon_global"]),
                    xytext=(6, 6), textcoords="offset points",
                    fontsize=7, color="#7f1d1d", fontweight="bold",
                )
                long_rows.append({
                    "dataset": dataset, "threshold": threshold,
                    "source": "quantmsdiann",
                    "label": f"quantmsdiann {row['version_label']}",
                    "library_kind": LIBRARY_KIND_PREDICTED,
                    "nr_prec": int(row["nr_prec"]),
                    "median_abs_epsilon_global": float(
                        row["median_abs_epsilon_global"]
                    ),
                })
        ax.set_title(
            _dataset_display_label(dataset).splitlines()[0],
            loc="left", fontsize=10, fontweight="bold",
        )
        ax.set_xlabel(f"Precursors quantified (≥{threshold} rep)", fontsize=9)
        ax.set_ylabel("Median |ε| (lower = more accurate)", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=8)
        # Enforce a minimum |ε| span per panel. Modules with a very narrow
        # spread (e.g. ZenoTOF, ~0.007; single-cell, ~0.01) would otherwise be
        # auto-scaled to fill the axis, stretching sub-0.01 version-to-version
        # wobble into a dramatic-looking zig-zag. A fixed floor keeps those
        # differences proportionate to the genuinely larger-spread panels.
        y_all = list(qm["median_abs_epsilon_global"]) + list(
            community["median_abs_epsilon_global"]
        )
        if y_all:
            MIN_EPS_SPAN = 0.05
            y_lo, y_hi = min(y_all), max(y_all)
            if (y_hi - y_lo) < MIN_EPS_SPAN:
                mid = (y_lo + y_hi) / 2.0
                ax.set_ylim(mid - MIN_EPS_SPAN / 2.0, mid + MIN_EPS_SPAN / 2.0)
        # When quantmsdiann lands toward bottom-right (more IDs + lower ε),
        # the panel's headline holds — annotate the cohort sizes.
        if len(qm) and len(community):
            n_lead_x = int((community["nr_prec"] < qm["nr_prec"].median()).sum())
            ax.text(
                0.02, 0.97,
                f"PB community: n={len(community)}    "
                f"qm versions: {len(qm)}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=7, color="#666666",
            )
    # Figure-level legend (data-driven over libraries actually drawn).
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles: list = []
    if qm_drawn:
        handles.append(
            Line2D([0], [0], marker="o", color="#d62728",
                   markerfacecolor="#d62728", markeredgecolor="#7f1d1d",
                   linewidth=1.2, markersize=8,
                   label="quantmsdiann trajectory")
        )
    for k in (
        LIBRARY_KIND_EMPIRICAL, LIBRARY_KIND_PREDICTED,
        LIBRARY_KIND_USER_DEFINED, LIBRARY_KIND_OTHER_TOOL,
    ):
        if k in libs_drawn:
            handles.append(Patch(facecolor=LIB_PALETTE[k], label=k))
    if handles:
        fig.legend(
            handles=handles, loc="upper center",
            bbox_to_anchor=(0.5, 1.02), ncol=min(5, len(handles)),
            fontsize=8, frameon=False,
        )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(long_rows)


# ---------------------------------------------------------------------------
# F1c — per-species log2 fold-change strip plot (supplementary)
# ---------------------------------------------------------------------------

def extract_qm_per_species_log2(
    dataset: str, threshold: int,
) -> pd.DataFrame:
    """Pull `mean_log2_empirical_<SPECIES>` per quantmsdiann version
    for the supplementary panel. Missing-species rows are dropped (the
    singlecell module has no E. coli)."""
    rows: list[dict] = []
    for version in DIANN_VERSIONS:
        cache_path = METRICS_CACHE_DIR / f"{dataset}_{version}.json"
        if not cache_path.exists():
            continue
        with open(cache_path, encoding="utf-8") as fh:
            payload = json.load(fh)
        res = payload.get("results", {}).get(str(threshold), {})
        for species in ("HUMAN", "YEAST", "ECOLI"):
            key = f"mean_log2_empirical_{species}"
            value = res.get(key)
            if value is None:
                continue
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if pd.isna(value):
                continue
            rows.append({
                "dataset": dataset,
                "version": version,
                "version_label": _VERSION_LABELS.get(version, version),
                "species": species,
                "mean_log2_empirical": value,
            })
    return pd.DataFrame(rows)


SPECIES_EXPECTED_LOG2_A_vs_B = {
    # Module-specific expected log2(A_vs_B), validated against ProteoBench's
    # ParseSettingsBuilder(...).build_parser("DIA-NN").species_expected_ratio():
    # the HYE modules (Astral/diaPASEF/ZenoTOF) use A_vs_B ratios Human 1.0,
    # Yeast 2.0, E. coli 0.25 -> log2 0 / +1 / -2; the single-cell module
    # (PXD049412) uses Human 1.2, Yeast 0.2 -> log2 +0.263 / -2.322.
    # test_per_species_expected_ratios asserts this dict matches ProteoBench.
    "quant_lfq_DIA_ion_singlecell": {
        "HUMAN": np.log2(1.2),
        "YEAST": np.log2(0.2),
    },
    "quant_lfq_DIA_ion_diaPASEF": {
        "HUMAN": 0.0,
        "YEAST": 1.0,
        "ECOLI": -2.0,
    },
    "quant_lfq_DIA_ion_ZenoTOF": {
        "HUMAN": 0.0,
        "YEAST": 1.0,
        "ECOLI": -2.0,
    },
    "quant_lfq_DIA_ion_Astral": {
        "HUMAN": 0.0,
        "YEAST": 1.0,
        "ECOLI": -2.0,
    },
}


def render_per_species_log2(
    threshold: int,
    svg_path: Path,
    *,
    datasets: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Strip plot of `mean_log2_empirical_<SPECIES>` per quantmsdiann
    version per dataset, with an expected-ratio reference line per
    species. The supplementary view of F1c."""
    datasets = (
        sorted(datasets, key=_dataset_sort_key)
        if datasets is not None
        else sorted(DATASET_TO_MODULE, key=_dataset_sort_key)
    )
    fig, axes = plt.subplots(
        nrows=len(datasets), ncols=1,
        figsize=(8.2, 2.4 * len(datasets)),
        sharex=False, squeeze=False,
    )
    species_colours = {
        "HUMAN": "#1976d2", "YEAST": "#d62728", "ECOLI": "#388e3c",
    }
    long_rows: list[dict] = []
    for i, dataset in enumerate(datasets):
        ax = axes[i, 0]
        df = extract_qm_per_species_log2(dataset, threshold)
        long_rows.extend(df.to_dict("records"))
        if df.empty:
            ax.text(
                0.5, 0.5, "no cached metrics yet for this dataset",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color="#888888",
            )
            ax.set_axis_off()
            continue
        species_order = list(df["species"].unique())
        for species in species_order:
            sub = df[df["species"] == species]
            ax.scatter(
                sub["version_label"], sub["mean_log2_empirical"],
                s=70, c=species_colours.get(species, "#9e9e9e"),
                edgecolors="#222222", linewidths=0.4,
                label=species, zorder=3,
            )
            expected = SPECIES_EXPECTED_LOG2_A_vs_B.get(
                DATASET_TO_MODULE[dataset], {}
            ).get(species)
            if expected is not None:
                ax.axhline(
                    expected,
                    color=species_colours.get(species, "#9e9e9e"),
                    linewidth=0.8, linestyle="--", alpha=0.6,
                    zorder=1,
                )
        ax.set_title(
            _dataset_display_label(dataset).splitlines()[0],
            loc="left", fontsize=10, fontweight="bold",
        )
        ax.set_ylabel("mean log2 (A/B)", fontsize=9)
        ax.axhline(0, color="#aaaaaa", linewidth=0.6, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=8)
        ax.legend(loc="best", fontsize=7, frameon=False, ncol=3)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(long_rows)


# ---------------------------------------------------------------------------
# Fig 2 accuracy figure — (b) per-species fold-change + (c) vs community
# ---------------------------------------------------------------------------

# Modules with predicted-library DIA-NN community comparators (for panel c).
_COMMUNITY_COMPARATOR_DATASETS = ("ProteoBench_Module_7", "PXD062685")
_SPECIES_X = {"HUMAN": 0, "YEAST": 1, "ECOLI": 2}
_SPECIES_LABEL = {"HUMAN": "Human", "YEAST": "Yeast", "ECOLI": "E. coli"}
_VERSION_BLUES = ["#bbdefb", "#64b5f6", "#1f77b4", "#1565c0", "#0d47a1"]


def render_accuracy_panels(
    threshold: int,
    svg_path: Path,
    *,
    datasets: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Fig 2 accuracy figure. Two stacked message-bearing panels:

      (b) Per-species fold-change accuracy (2x2, one cell per module). x =
          HYE species, y = measured log2 ratio; a dashed line marks the
          ProteoBench-expected ratio per species (accuracy = distance to the
          line) and the five DIA-NN versions are overlaid (light->dark) so
          their tight clustering shows that accuracy is version-invariant.

      (c) quantmsdiann within the predicted-library community (1x2, only the
          modules with predicted-library DIA-NN comparators). Box + strip of
          the community median |eps| with the five quantmsdiann versions
          overlaid as a tight cluster.

    Returns the long-format audit table of every plotted point."""
    datasets = (
        sorted(datasets, key=_dataset_sort_key)
        if datasets is not None
        else sorted(DATASET_TO_MODULE, key=_dataset_sort_key)
    )
    fig = plt.figure(figsize=(9.5, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.95],
                          hspace=0.5, wspace=0.26)
    long_rows: list[dict] = []

    # ---- Panel (b): per-species fold-change accuracy ----
    for idx, dataset in enumerate(datasets):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        df = extract_qm_per_species_log2(dataset, threshold)
        expected = SPECIES_EXPECTED_LOG2_A_vs_B.get(
            DATASET_TO_MODULE.get(dataset, ""), {}
        )
        present = [s for s in ("HUMAN", "YEAST", "ECOLI") if s in expected]
        for species in present:
            x = _SPECIES_X[species]
            ax.hlines(expected[species], x - 0.32, x + 0.32,
                      color="#444444", ls="--", lw=1.1, zorder=1)
            sub = df[df["species"] == species]
            for _, row in sub.iterrows():
                vi = DIANN_VERSIONS.index(row["version"])
                ax.scatter(
                    x - 0.22 + 0.11 * vi, row["mean_log2_empirical"],
                    s=34, color=_VERSION_BLUES[vi],
                    edgecolor="#333333", linewidths=0.35, zorder=3,
                )
                long_rows.append({
                    "panel": "b", "dataset": dataset, "threshold": threshold,
                    "species": species, "version": row["version"],
                    "measured_log2": float(row["mean_log2_empirical"]),
                    "expected_log2": float(expected[species]),
                })
        ax.set_xticks([_SPECIES_X[s] for s in present])
        ax.set_xticklabels([_SPECIES_LABEL[s] for s in present], fontsize=8)
        ax.set_xlim(-0.5, max(_SPECIES_X[s] for s in present) + 0.5
                    if present else 2.5)
        _label_parts = _dataset_display_label(dataset).split("\n")
        ax.set_title(
            _label_parts[1] if len(_label_parts) > 1 else _label_parts[0],
            loc="left", fontsize=9, fontweight="bold",
        )
        ax.set_ylabel("Measured log$_2$ ratio", fontsize=8)
        ax.tick_params(axis="both", labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.text(0.005, 0.99,
             "(b) Per-species fold-change accuracy "
             "(dashed = expected ratio)",
             fontsize=10, fontweight="bold", va="top")
    # Version colour key (light -> dark = oldest -> newest DIA-NN release).
    from matplotlib.lines import Line2D
    version_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=7,
               markerfacecolor=_VERSION_BLUES[i], markeredgecolor="#333333",
               label=_VERSION_LABELS.get(v, v))
        for i, v in enumerate(DIANN_VERSIONS)
    ]
    fig.legend(handles=version_handles, loc="upper right",
               bbox_to_anchor=(0.99, 1.005), ncol=5, fontsize=7.5,
               frameon=False, title="DIA-NN version", title_fontsize=7.5,
               handletextpad=0.2, columnspacing=0.9)

    # ---- Panel (c): quantmsdiann within the predicted-library community ----
    comp = [d for d in datasets if d in _COMMUNITY_COMPARATOR_DATASETS]
    for j, dataset in enumerate(comp):
        ax = fig.add_subplot(gs[2, j])
        community = extract_community_id_vs_eps(dataset, threshold)
        community = community[
            community["library_kind"] == LIBRARY_KIND_PREDICTED
        ]["median_abs_epsilon_global"].astype(float).values
        qm = extract_quantmsdiann_id_vs_eps(dataset, threshold)
        qm_eps = qm["median_abs_epsilon_global"].astype(float).values
        if len(community):
            ax.boxplot([community], positions=[0], widths=0.5,
                       showfliers=False,
                       medianprops=dict(color="#37474f", lw=1.2))
            cx = ([0.0] if len(community) == 1
                  else list(np.linspace(-0.16, 0.16, len(community))))
            ax.scatter(cx, community, s=34, color="#90a4ae",
                       edgecolor="#555555", linewidths=0.3, alpha=0.85,
                       label=f"community (n={len(community)})", zorder=3)
        if len(qm_eps):
            qx = ([1.0] if len(qm_eps) == 1
                  else list(np.linspace(0.84, 1.16, len(qm_eps))))
            ax.scatter(qx, qm_eps, s=60, color="#d62728",
                       edgecolor="#7f1d1d", linewidths=0.6,
                       label="quantmsdiann (5 versions)", zorder=4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["community", "quantmsdiann"], fontsize=8)
        ax.set_xlim(-0.6, 1.6)
        _label_parts = _dataset_display_label(dataset).split("\n")
        ax.set_title(
            _label_parts[1] if len(_label_parts) > 1 else _label_parts[0],
            loc="left", fontsize=9, fontweight="bold",
        )
        ax.set_ylabel("Median |ε|", fontsize=8)
        ax.tick_params(axis="both", labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=6.5, frameon=False, loc="upper right")
        for _, row in qm.iterrows():
            long_rows.append({
                "panel": "c", "dataset": dataset, "threshold": threshold,
                "species": None, "version": row["version"],
                "measured_log2": None,
                "median_abs_epsilon": float(row["median_abs_epsilon_global"]),
            })
    fig.text(0.005, 0.34,
             "(c) quantmsdiann within the predicted-library community "
             "(median |ε|; lower = more accurate)",
             fontsize=10, fontweight="bold", va="top")

    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(long_rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:  # pragma: no cover
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = FIGURES_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    supp_dir = FIGURES_DIR / "supplementary"
    supp_dir.mkdir(parents=True, exist_ok=True)

    # If any cache file is missing the F1b figure still renders with
    # whatever rows exist; callers can re-run after populating the
    # cache via `from analysis.proteobench_metrics import cached_proteobench_metrics`.
    for dataset in DATASET_TO_MODULE:
        for version in DIANN_VERSIONS:
            cache_path = METRICS_CACHE_DIR / f"{dataset}_{version}.json"
            if not cache_path.exists():
                print(
                    f"warn: missing metrics cache "
                    f"{dataset}/{version} — F1b will skip this point. "
                    "Populate via cached_proteobench_metrics(...)."
                )

    # Fig 2 accuracy figure (main): (b) per-species fold-change accuracy +
    # (c) quantmsdiann within the predicted-library community.
    acc_df = render_accuracy_panels(
        threshold=3, svg_path=FIGURES_DIR / "main_accuracy.svg",
    )
    acc_df.to_csv(data_dir / "accuracy_min3.tsv", sep="\t", index=False)
    print(
        f"Fig 2b/c (≥3 rep): {acc_df.shape[0]} points rendered to "
        f"{FIGURES_DIR / 'main_accuracy.svg'}"
    )

    # Identifications-vs-accuracy scatter, demoted to the supplement: it
    # over-resolves the (non-significant) version-to-version |eps| spread, so
    # it is no longer the main accuracy panel — kept for completeness.
    long_df = render_id_vs_epsilon(
        threshold=3, svg_path=supp_dir / "supp_id_vs_epsilon_min3.svg",
    )
    long_df.to_csv(
        data_dir / "id_vs_epsilon_min3.tsv", sep="\t", index=False,
    )
    print(
        f"supp id-vs-ε (≥3 rep): {long_df.shape[0]} points rendered to "
        f"{supp_dir / 'supp_id_vs_epsilon_min3.svg'}"
    )

    # Per-species strip (now subsumed by main panel (b)) retained as a
    # supplementary cross-check.
    supp_df = render_per_species_log2(
        threshold=3, svg_path=supp_dir / "supp_per_species_log2_min3.svg",
    )
    supp_df.to_csv(
        data_dir / "per_species_log2_min3.tsv", sep="\t", index=False,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
