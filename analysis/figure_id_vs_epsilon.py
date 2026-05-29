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
                    "quantmsdiann (DIA-NN, empirical lib)"
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
                    "library_kind": LIBRARY_KIND_EMPIRICAL,
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
    # Module-specific expected log2(A/B). HYE mix on Modules 5/7/10 has
    # the classic 1:2:0.5 design (human=0, yeast=+1, ecoli=-1). Module 9
    # (singlecell) uses 1.2:0.2 (human=+0.26, yeast=-2.32). We display
    # both the expected and empirical values per species per version.
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

    long_df = render_id_vs_epsilon(
        threshold=3, svg_path=FIGURES_DIR / "main_id_vs_epsilon.svg",
    )
    long_df.to_csv(
        data_dir / "id_vs_epsilon_min3.tsv", sep="\t", index=False,
    )
    print(
        f"F1b (≥3 rep): {long_df.shape[0]} points rendered to "
        f"{FIGURES_DIR / 'main_id_vs_epsilon.svg'}"
    )

    supp_df = render_per_species_log2(
        threshold=3, svg_path=supp_dir / "supp_per_species_log2_min3.svg",
    )
    supp_df.to_csv(
        data_dir / "per_species_log2_min3.tsv", sep="\t", index=False,
    )
    print(
        f"F1c (≥3 rep): {supp_df.shape[0]} points rendered to "
        f"{supp_dir / 'supp_per_species_log2_min3.svg'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
