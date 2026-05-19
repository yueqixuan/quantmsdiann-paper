"""quantmsdiann vs ProteoBench benchmarks figure.

Each of the four ProteoBench benchmark datasets has been processed through
quantmsdiann using five DIA-NN versions (1.8.1 -> 2.5.0). This script:

- Pulls each version's `diann_report.{pr,pg}_matrix.tsv` from PRIDE FTP
  (cached on disk).
- Counts data rows = precursors / protein groups quantified at 1% FDR
  with the matrix-level filter that is identical to ProteoBench's
  `nr_prec` definition.
- Fetches public ProteoBench submissions for the matching module from the
  Proteobench/Results_quant_ion_DIA_<module> repo (cached on disk).
- Renders two paper-ready figure sets and an auditable counts.tsv.

We count matrix rows rather than `diannsummary.log` headline numbers because
the log format is not uniform across DIA-NN versions: 1.8.1 doesn't print
the "Protein groups with global q-value <= 0.01" line and 2.1.0-2.3.2 don't
print the "Target precursors at 1% global q-value" line. Matrix row counts
are well-defined across all versions and match how ProteoBench computes
nr_prec from the same pr_matrix.tsv.

ProteoBench ion-level modules report precursors only — there is no protein
count in their public datapoints, so the head-to-head supp is precursors
only.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Iterable, Iterator

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import requests

from analysis.figure_original_vs_quantmsdiann import (
    download_if_missing,
    SUMMARY_LOG_PROTEIN_LINE_RE,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "quantmsdiann_benchmarks"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "quantmsdiann_benchmarks"

PRIDE_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/proteobench/quantmsdiann_results"
)

DIANN_VERSIONS = ("v1_8_1", "v2_1_0", "v2_2_0", "v2_3_2", "v2_5_0")

# Module mapping verified against `docs/available-modules/active-modules/*.md`
# in Proteobench/ProteoBench. Module 7 (DIA Astral 2Th) has no public PXD yet
# (raw files at proteobench.cubimed.rub.de); the others are on ProteomeXchange.
DATASET_TO_MODULE: dict[str, dict[str, str]] = {
    "PXD049412": {
        "label": "Module 9 - DIA single-cell",
        "results_repo": "Proteobench/Results_quant_ion_DIA_singlecell",
    },
    "PXD062685": {
        "label": "Module 5 - DIA diaPASEF",
        "results_repo": "Proteobench/Results_quant_ion_DIA_diaPASEF",
    },
    "PXD070049": {
        "label": "Module 10 - DIA ZenoTOF",
        "results_repo": "Proteobench/Results_quant_ion_DIA_ZenoTOF",
    },
    "ProteoBench_Module_7": {
        "label": "Module 7 - DIA Astral 2Th",
        "results_repo": "Proteobench/Results_quant_ion_DIA_Astral",
    },
}

GITHUB_API_BASE = "https://api.github.com"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com"

SUMMARY_LOG_PRECURSOR_LINE_RE = re.compile(
    r"Target precursors at 1% global q-value:\s*(\d+)"
)


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def parse_diann_summary_log(log_path: Path) -> tuple[int, int]:
    """Return (protein_groups, target_precursors) from a DIA-NN summary log."""
    protein_groups: int | None = None
    precursors: int | None = None
    with open(log_path, encoding="utf-8") as fh:
        for line in fh:
            m = SUMMARY_LOG_PROTEIN_LINE_RE.search(line)
            if m and protein_groups is None:
                protein_groups = int(m.group(1))
                continue
            m = SUMMARY_LOG_PRECURSOR_LINE_RE.search(line)
            if m and precursors is None:
                precursors = int(m.group(1))
    if protein_groups is None:
        raise ValueError(
            "'Protein groups with global q-value <= 0.01: N' not found in log"
        )
    if precursors is None:
        raise ValueError(
            "'Target precursors at 1% global q-value: N' not found in log"
        )
    return protein_groups, precursors


def count_matrix_data_rows(matrix_path: Path) -> int:
    """Count data rows in a TSV (lines minus the header). Used for both
    pr_matrix.tsv (one row per precursor at 1% precursor+pg FDR) and
    pg_matrix.tsv (one row per protein group). Streams the file so we
    never load multi-MB matrices into memory."""
    n = 0
    with open(matrix_path, "rb") as fh:
        for _ in fh:
            n += 1
    # The matrix always has at least a header row; subtract it. Treat an
    # empty file as 0 data rows rather than -1.
    return max(0, n - 1)


LIBRARY_KIND_EMPIRICAL = "empirical"
LIBRARY_KIND_PREDICTED = "predicted (DIANN)"
LIBRARY_KIND_USER_DEFINED = "user-defined speclib"
LIBRARY_KIND_OTHER_TOOL = "other tool"


def classify_predictors_library(predictors_library) -> str:
    """Map a ProteoBench `predictors_library` value to one of three categories
    used to color DIA-NN community submissions in the supp figure:
      - 'empirical' (None / 'None') — DIA-NN built the library from the data
        itself; matches quantmsdiann's two-pass workflow.
      - 'predicted (DIANN)' — DIA-NN's built-in in-silico prediction (the
        FASTA-derived predicted library); larger search space.
      - 'user-defined speclib' — externally uploaded/curated library.
    Any other value returns the original string."""
    if predictors_library is None or predictors_library == "None":
        return LIBRARY_KIND_EMPIRICAL
    if isinstance(predictors_library, dict):
        vals = set(predictors_library.values())
        if vals == {"DIANN"}:
            return LIBRARY_KIND_PREDICTED
        if vals == {"User defined speclib"}:
            return LIBRARY_KIND_USER_DEFINED
        return ",".join(sorted(vals))
    s = str(predictors_library).strip()
    # Some ProteoBench submissions store the predictors_library as a Python
    # repr ("{'RT': 'DIANN', 'IM': 'DIANN', 'MS2_int': 'DIANN'}") instead of
    # a JSON object; try literal-eval before falling through to substring
    # heuristics.
    if s.startswith("{") and s.endswith("}"):
        try:
            import ast
            d = ast.literal_eval(s)
            if isinstance(d, dict):
                return classify_predictors_library(d)
        except (ValueError, SyntaxError):
            pass
    if "User" in s and "speclib" in s.lower():
        return LIBRARY_KIND_USER_DEFINED
    if s.upper() == "DIANN":
        return LIBRARY_KIND_PREDICTED
    return s


def parse_proteobench_datapoints(
    json_path: Path,
) -> Iterator[tuple[str, str, int, str]]:
    """Yield (software_name, software_version, nr_prec, library_kind) per
    ProteoBench submission. The library kind is derived from each entry's
    `predictors_library` field via `classify_predictors_library`; non-DIA-NN
    tools are tagged 'other tool' since their library systems aren't directly
    comparable to DIA-NN's three categories.

    Skips entries that don't carry an integer nr_prec — those are rare
    submission artefacts and rendering them as zero would distort the
    per-module distribution."""
    with open(json_path, encoding="utf-8") as fh:
        payload = json.load(fh)
    for entry in payload:
        nr_prec = entry.get("nr_prec")
        if not isinstance(nr_prec, (int, float)) or pd.isna(nr_prec):
            continue
        software = entry.get("software_name") or ""
        version = entry.get("software_version") or ""
        if normalise_software_name(str(software)) == "dia-nn":
            kind = classify_predictors_library(entry.get("predictors_library"))
        else:
            kind = LIBRARY_KIND_OTHER_TOOL
        yield str(software), str(version), int(nr_prec), kind


def normalise_software_name(name: str) -> str:
    """Collapse case/whitespace variants. DIA-NN appears as 'DIA-NN', 'DIANN',
    'Diann' across ProteoBench submissions; we normalise so highlight
    overlays match every spelling."""
    s = name.strip().lower()
    if s in {"dia-nn", "diann", "dia nn"}:
        return "dia-nn"
    return s


# ---------------------------------------------------------------------------
# ProteoBench fetch + cache
# ---------------------------------------------------------------------------


def fetch_proteobench_module(repo: str, dest: Path) -> Path:
    """Fetch every `<hash>.json` datapoint from a Proteobench results repo
    and write them as a single consolidated list to `dest`. Idempotent:
    skips if `dest` already exists and is non-empty.

    The GitHub contents API is rate-limited (60 req/h unauthenticated) but
    we hit it once per module to list, plus N raw.githubusercontent.com
    fetches; the raw host is not rate-limited."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    listing_url = f"{GITHUB_API_BASE}/repos/{repo}/contents/"
    resp = requests.get(listing_url, timeout=60)
    resp.raise_for_status()
    listing = resp.json()
    if not isinstance(listing, list):
        raise RuntimeError(
            f"GitHub listing for {repo} returned non-list payload: {listing}"
        )
    json_names = sorted(
        item["name"] for item in listing
        if isinstance(item, dict) and item.get("name", "").endswith(".json")
    )
    items: list[dict] = []
    for name in json_names:
        raw_url = f"{GITHUB_RAW_BASE}/{repo}/main/{name}"
        r = requests.get(raw_url, timeout=60)
        r.raise_for_status()
        # ProteoBench JSON files use the bare token `NaN` (not valid JSON);
        # the standard library accepts it because it allows non-finite floats,
        # but `requests.json()` uses the same loader so this Just Works.
        items.append(r.json())
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(items, fh)
    return dest


def consolidate_proteobench_datapoints(
    files: Iterable[Path], dest: Path,
) -> Path:
    """Combine per-submission JSON files into a single list, sorted by file
    name. Used by tests to exercise the same merge logic as
    `fetch_proteobench_module` without hitting the network."""
    items: list[dict] = []
    for p in sorted(files, key=lambda x: x.name):
        with open(p, encoding="utf-8") as fh:
            items.append(json.load(fh))
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(items, fh)
    return dest


# ---------------------------------------------------------------------------
# Long-format aggregator
# ---------------------------------------------------------------------------


def build_long_table(
    quantmsdiann_rows: list[tuple[str, str, int, int]],
    proteobench_rows: dict[str, list[tuple[str, str, int]]],
) -> pd.DataFrame:
    """Combine quantmsdiann rows (dataset, version, precursors, proteins)
    and ProteoBench rows (per dataset: list of (tool, version, precursors))
    into a long-format DataFrame with columns
    [dataset, source, tool, version, precursors, proteins]."""
    rows: list[dict] = []
    for dataset, version, precursors, proteins in quantmsdiann_rows:
        rows.append({
            "dataset": dataset,
            "source": "quantmsdiann",
            "tool": "DIA-NN",
            "version": version,
            "precursors": precursors,
            "proteins": proteins,
            # quantmsdiann's two-pass workflow uses an empirical library
            # built from the data; tag it here so render_vs_proteobench can
            # color-match against community submissions using the same
            # library strategy.
            "library_kind": LIBRARY_KIND_EMPIRICAL,
        })
    for dataset, entries in proteobench_rows.items():
        for tool, version, precursors, lib_kind in entries:
            rows.append({
                "dataset": dataset,
                "source": "proteobench",
                "tool": tool,
                "version": version,
                "precursors": precursors,
                "proteins": None,
                "library_kind": lib_kind,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


_VERSION_LABELS = {
    "v1_8_1": "1.8.1",
    "v2_1_0": "2.1.0",
    "v2_2_0": "2.2.0",
    "v2_3_2": "2.3.2",
    "v2_5_0": "2.5.0",
}


def render_main_overview(
    quantmsdiann_rows: list[tuple[str, str, int, int]],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """4x2 grid: rows = datasets, columns = (precursors, protein groups).
    Each panel is a bar chart of DIA-NN versions. Paper-ready: no title,
    no footer; only axis labels and an in-figure legend."""
    df = pd.DataFrame(
        quantmsdiann_rows,
        columns=["dataset", "version", "precursors", "proteins"],
    )
    datasets = sorted(df["dataset"].unique(), key=_dataset_sort_key)
    fig, axes = plt.subplots(
        nrows=len(datasets), ncols=2,
        figsize=(11.5, 2.4 * len(datasets)),
        squeeze=False,
    )
    for i, dataset in enumerate(datasets):
        sub = df[df["dataset"] == dataset].sort_values(
            "version", key=lambda s: s.map(lambda v: DIANN_VERSIONS.index(v)),
        )
        x_labels = [_VERSION_LABELS.get(v, v) for v in sub["version"]]
        for j, (metric, colour) in enumerate(
            (("precursors", "#1f77b4"), ("proteins", "#d62728"))
        ):
            ax = axes[i, j]
            vals = sub[metric].tolist()
            bars = ax.bar(x_labels, vals, color=colour, width=0.6)
            for bar, v in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f"{v:,}", ha="center", va="bottom", fontsize=8)
            ax.set_ylabel(
                "Precursors (1% FDR)" if metric == "precursors"
                else "Protein groups (1% FDR)",
                fontsize=9,
            )
            if i == len(datasets) - 1:
                ax.set_xlabel("DIA-NN version", fontsize=9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ymax = max(vals) * 1.18 if vals else 1
            ax.set_ylim(0, ymax)
            ax.tick_params(axis="both", labelsize=8)
        # Dataset label as left-most y-axis annotation, so each row's
        # identity reads top-to-bottom alongside the bars. We push it well
        # outside the axis (x=-0.45) and reserve 16% of the figure width
        # for it in tight_layout's rect so it doesn't overlap y-ticks.
        axes[i, 0].text(
            -0.45, 0.5, _dataset_display_label(dataset),
            transform=axes[i, 0].transAxes,
            ha="right", va="center", fontsize=10, fontweight="bold",
        )
    fig.tight_layout(rect=(0.16, 0, 1, 1))
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def render_main_precursors(
    quantmsdiann_rows: list[tuple[str, str, int, int]],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """Precursor-only headline: 4 datasets x 5 DIA-NN versions in a single
    grouped-bar panel. The cross-dataset / cross-version comparison is the
    story we want to lead with (ProteoBench DIA ion modules publish only
    precursor numbers, so it's also the metric with the cleanest external
    benchmark). Paper-ready: no title, no footer."""
    df = pd.DataFrame(
        quantmsdiann_rows,
        columns=["dataset", "version", "precursors", "proteins"],
    )
    datasets = sorted(df["dataset"].unique(), key=_dataset_sort_key)
    version_palette = [
        "#bbdefb", "#64b5f6", "#1f77b4", "#1565c0", "#0d47a1",
    ]
    fig, ax = plt.subplots(figsize=(10.5, 5))
    n_versions = len(DIANN_VERSIONS)
    bar_width = 0.8 / n_versions
    x = list(range(len(datasets)))
    for k, version in enumerate(DIANN_VERSIONS):
        vals = [
            int(df[(df["dataset"] == d) & (df["version"] == version)]
                ["precursors"].iloc[0])
            if not df[(df["dataset"] == d) & (df["version"] == version)].empty
            else 0
            for d in datasets
        ]
        offsets = [xi + (k - (n_versions - 1) / 2) * bar_width for xi in x]
        bars = ax.bar(
            offsets, vals, width=bar_width,
            color=version_palette[k % len(version_palette)],
            label=_VERSION_LABELS.get(version, version),
        )
        for bar, v in zip(bars, vals):
            if v == 0:
                continue
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{v / 1000:.0f}k", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [_dataset_display_label(d).replace("\n", "\n") for d in datasets],
        fontsize=9,
    )
    ax.set_ylabel("Precursors quantified (1% FDR)")
    ymax = df["precursors"].max() * 1.15
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        title="DIA-NN version", loc="upper center",
        bbox_to_anchor=(0.5, -0.18), ncol=n_versions, frameon=False,
        fontsize=8, title_fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def render_vs_proteobench(
    long_df: pd.DataFrame,
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
    *,
    library_kinds: tuple[str, ...] = (LIBRARY_KIND_EMPIRICAL,),
) -> None:
    """One panel per dataset stacked vertically. Horizontal bar chart of
    matching-library ProteoBench DIA-NN submissions plus our five DIA-NN
    versions as red bars below them, sharing the same precursor-count x
    axis. Paper-ready: no title, no footer.

    `library_kinds` filters ProteoBench submissions to only those tagged
    with one of the listed library strategies (default: empirical, matching
    quantmsdiann). Pass e.g. `(LIBRARY_KIND_EMPIRICAL, LIBRARY_KIND_PREDICTED)`
    to broaden the comparison; pass an empty tuple to include all
    submissions."""
    datasets = sorted(long_df["dataset"].unique(), key=_dataset_sort_key)
    if library_kinds:
        long_df = long_df[
            (long_df["source"] == "quantmsdiann")
            | long_df["library_kind"].isin(library_kinds)
        ].copy()
    # Per-panel height proportional to ProteoBench-submission count so a
    # crowded panel doesn't crush a sparse one.
    pb_counts = {
        ds: int(((long_df["dataset"] == ds)
                 & (long_df["source"] == "proteobench")).sum())
        for ds in datasets
    }
    heights = [max(2.5, 0.18 * (pb_counts[ds] + 6)) for ds in datasets]
    fig, axes = plt.subplots(
        nrows=len(datasets), ncols=1,
        figsize=(10, sum(heights) + 1),
        gridspec_kw={"height_ratios": heights},
        squeeze=False,
    )
    for i, dataset in enumerate(datasets):
        ax = axes[i, 0]
        sub = long_df[long_df["dataset"] == dataset].copy()
        sub["norm_tool"] = sub["tool"].map(normalise_software_name)
        # ProteoBench community submissions: sort ascending by precursors so
        # the longest bars are on the right; reduces visual clutter.
        pb = (sub[sub["source"] == "proteobench"]
              .sort_values("precursors", ascending=True)
              .reset_index(drop=True))
        labels = [f"{t} {v}".strip() for t, v in zip(pb["tool"], pb["version"])]
        # Color ProteoBench bars by the library strategy each submission
        # used: matching quantmsdiann's empirical library is the relevant
        # apples-to-apples set; predicted (DIANN in-silico) and user-defined
        # speclibs typically search a larger candidate space. Non-DIA-NN
        # tools stay neutral grey.
        library_palette = {
            LIBRARY_KIND_EMPIRICAL: "#80cbc4",
            LIBRARY_KIND_PREDICTED: "#1976d2",
            LIBRARY_KIND_USER_DEFINED: "#9575cd",
            LIBRARY_KIND_OTHER_TOOL: "#bdbdbd",
        }
        kinds = list(pb["library_kind"]) if "library_kind" in pb.columns else [
            LIBRARY_KIND_OTHER_TOOL for _ in range(len(pb))
        ]
        bar_colours = [library_palette.get(k, "#bdbdbd") for k in kinds]
        ax.barh(
            range(len(pb)), pb["precursors"], color=bar_colours, height=0.7,
        )
        if len(pb) == 0:
            # No comparable ProteoBench submissions for this dataset under
            # the selected library_kinds filter. Anchor a small annotation
            # above the quantmsdiann band so the empty PB region reads as
            # intentional, not a rendering bug.
            ax.text(
                0.02, 0.9,
                "no ProteoBench submission with matching library strategy",
                transform=ax.transAxes, fontsize=8, color="#666666",
                fontstyle="italic",
            )
        # Build proxy artists so the figure-level legend can show each
        # library category once.
        if i == 0:
            from matplotlib.patches import Patch
            ax._library_legend_handles = [
                Patch(facecolor=library_palette[k], label=k)
                for k in (
                    LIBRARY_KIND_EMPIRICAL, LIBRARY_KIND_PREDICTED,
                    LIBRARY_KIND_USER_DEFINED, LIBRARY_KIND_OTHER_TOOL,
                )
                if k in set(kinds) or k == LIBRARY_KIND_EMPIRICAL
            ]
        # Annotate each bar with its tool+version at the right tip; this is
        # much more legible than y-tick labels when N is >15 per panel.
        xmax_pb = float(pb["precursors"].max()) if len(pb) else 0.0
        for k, (precs, lab) in enumerate(zip(pb["precursors"], labels)):
            ax.text(
                precs + xmax_pb * 0.005, k, lab,
                va="center", ha="left", fontsize=6, color="#404040",
            )
        ax.set_yticks([])  # bar-tip annotations replace tick labels
        qm = sub[sub["source"] == "quantmsdiann"].sort_values(
            "version", key=lambda s: s.map(lambda v: DIANN_VERSIONS.index(v)),
        )
        # Plot quantmsdiann versions as red horizontal bars below the
        # ProteoBench block (y = -1, -2, ...), with a gap row at y=-0.5 to
        # visually separate the two groups. A pale red background span +
        # left bracket highlight the quantmsdiann region.
        qm_ys = [-1.0 - k for k in range(len(qm))]
        qm_vals = list(qm["precursors"])
        qm_labels = [
            f"quantmsdiann {_VERSION_LABELS.get(v, v)}"
            for v in qm["version"]
        ]
        if qm_ys:
            top = qm_ys[0] + 0.55
            bottom = qm_ys[-1] - 0.55
            ax.axhspan(bottom, top, color="#ffebee", zorder=0)
            ax.barh(
                qm_ys, qm_vals, color="#d62728", height=0.7,
                edgecolor="#7f1d1d", linewidth=0.6, zorder=2,
                label="quantmsdiann (DIA-NN)" if i == 0 else None,
            )
            xmax_qm = max(qm_vals) if qm_vals else 0.0
            for y, v, lab in zip(qm_ys, qm_vals, qm_labels):
                ax.text(
                    v + max(xmax_pb, xmax_qm) * 0.005, y, lab,
                    va="center", ha="left", fontsize=7,
                    color="#7f1d1d", fontweight="bold",
                )
        else:
            xmax_qm = 0.0
        ax.set_ylim(min(qm_ys) - 1.0 if qm_ys else -1.5, len(pb))
        # Pad x-axis so right-side annotations don't get clipped.
        ax.set_xlim(0, max(xmax_pb, xmax_qm) * 1.25)
        ax.set_xlabel("Precursors quantified", fontsize=9)
        ax.set_title(
            _dataset_display_label(dataset),
            loc="left", fontsize=10, fontweight="bold",
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)
    from matplotlib.patches import Patch
    quantms_handle = Patch(
        facecolor="#d62728", edgecolor="#7f1d1d", linewidth=0.8,
        label="quantmsdiann (DIA-NN, empirical library)",
    )
    lib_handles = getattr(axes[0, 0], "_library_legend_handles", [])
    # Relabel library-kind patches for the legend so they read like
    # "ProteoBench DIA-NN, <kind>" rather than the bare kind string.
    relabeled = []
    for h in lib_handles:
        lab = h.get_label()
        if lab == LIBRARY_KIND_EMPIRICAL:
            relabeled.append(Patch(facecolor=h.get_facecolor(),
                                   label="ProteoBench DIA-NN, empirical library"))
        elif lab == LIBRARY_KIND_PREDICTED:
            relabeled.append(Patch(facecolor=h.get_facecolor(),
                                   label="ProteoBench DIA-NN, DIANN-predicted library"))
        elif lab == LIBRARY_KIND_USER_DEFINED:
            relabeled.append(Patch(facecolor=h.get_facecolor(),
                                   label="ProteoBench DIA-NN, user-defined speclib"))
        elif lab == LIBRARY_KIND_OTHER_TOOL:
            relabeled.append(Patch(facecolor=h.get_facecolor(),
                                   label="ProteoBench other tool"))
        else:
            relabeled.append(h)
    legend_handles = [quantms_handle] + relabeled
    fig.legend(
        handles=legend_handles,
        loc="upper center", bbox_to_anchor=(0.5, 1.0),
        ncol=2, frameon=False, fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def _dataset_sort_key(name: str) -> tuple[int, str]:
    # Module 7 (Astral) first because it's the headline dataset, then PXDs
    # in numeric order — keeps figure rows in a stable order across runs.
    if name == "ProteoBench_Module_7":
        return (0, name)
    return (1, name)


def _dataset_display_label(name: str) -> str:
    info = DATASET_TO_MODULE.get(name)
    if info is None:
        return name
    return f"{name}\n{info['label']}"


# ---------------------------------------------------------------------------
# counts.tsv writer
# ---------------------------------------------------------------------------


def write_counts_tsv(long_df: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["dataset", "source", "tool", "version", "precursors", "proteins"]
    # Stable ordering: quantmsdiann rows first (by dataset, version), then
    # ProteoBench rows (by dataset, descending precursors so the strongest
    # submissions are at the top of each block).
    qm = (long_df[long_df["source"] == "quantmsdiann"]
          .sort_values(["dataset", "version"]))
    pb = (long_df[long_df["source"] == "proteobench"]
          .sort_values(["dataset", "precursors"], ascending=[True, False]))
    out = pd.concat([qm, pb], ignore_index=True)
    # Pandas upgrades int+None columns to float; use pandas' nullable Int64
    # for proteins so ProteoBench rows print as the empty string and
    # quantmsdiann rows print as plain integers (no trailing .0).
    out = out[cols].copy()
    out["precursors"] = out["precursors"].astype("Int64")
    out["proteins"] = out["proteins"].astype("Int64")
    out.to_csv(tsv_path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:  # pragma: no cover
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    quantmsdiann_rows: list[tuple[str, str, int, int]] = []
    for dataset in DATASET_TO_MODULE:
        for version in DIANN_VERSIONS:
            base = f"{PRIDE_BASE}/{dataset}/{version}/quant_tables"
            ddir = DATA_DIR / dataset / version
            try:
                pr_path = download_if_missing(
                    f"{base}/diann_report.pr_matrix.tsv",
                    ddir / "diann_report.pr_matrix.tsv",
                )
                pg_path = download_if_missing(
                    f"{base}/diann_report.pg_matrix.tsv",
                    ddir / "diann_report.pg_matrix.tsv",
                )
            except Exception as exc:
                print(f"WARN: failed to fetch matrices for {dataset}/{version}: "
                      f"{exc}", file=sys.stderr)
                continue
            precursors = count_matrix_data_rows(pr_path)
            proteins = count_matrix_data_rows(pg_path)
            quantmsdiann_rows.append((dataset, version, precursors, proteins))
            print(f"{dataset} {version}: precursors={precursors:,}  "
                  f"proteins={proteins:,}")

    proteobench_rows: dict[str, list[tuple[str, str, int]]] = {}
    for dataset, info in DATASET_TO_MODULE.items():
        cache = DATA_DIR / "proteobench" / f"{dataset}.json"
        try:
            fetch_proteobench_module(info["results_repo"], cache)
        except Exception as exc:
            print(f"WARN: failed to fetch ProteoBench for {dataset}: {exc}",
                  file=sys.stderr)
            proteobench_rows[dataset] = []
            continue
        entries = list(parse_proteobench_datapoints(cache))
        proteobench_rows[dataset] = entries
        print(f"{dataset}: {len(entries)} ProteoBench submissions")

    long_df = build_long_table(quantmsdiann_rows, proteobench_rows)

    print("Rendering precursor-only main panel...")
    render_main_precursors(
        quantmsdiann_rows,
        FIGURES_DIR / "main_benchmarks_precursors.pdf",
        FIGURES_DIR / "main_benchmarks_precursors.png",
        FIGURES_DIR / "main_benchmarks_precursors.svg",
    )

    print("Rendering full main benchmarks overview (precursors + proteins)...")
    render_main_overview(
        quantmsdiann_rows,
        FIGURES_DIR / "main_benchmarks_overview.pdf",
        FIGURES_DIR / "main_benchmarks_overview.png",
        FIGURES_DIR / "main_benchmarks_overview.svg",
    )

    print("Rendering ProteoBench-overlay supp figure...")
    render_vs_proteobench(
        long_df,
        FIGURES_DIR / "supp_vs_proteobench.pdf",
        FIGURES_DIR / "supp_vs_proteobench.png",
        FIGURES_DIR / "supp_vs_proteobench.svg",
    )

    print("Writing auditable counts TSV...")
    write_counts_tsv(long_df, FIGURES_DIR / "counts.tsv")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
