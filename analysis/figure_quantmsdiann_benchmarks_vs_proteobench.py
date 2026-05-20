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


def extract_nr_prec_at_replicate_threshold(
    entry: dict, min_replicates: int,
) -> int | None:
    """Return the precursor count for a single ProteoBench submission
    `entry` at the requested ≥min_replicates threshold.

    Each ProteoBench submission embeds a `results` dict keyed by replicate
    count as a string ('1', '2', ..., '6'). The value at key 'K' is the
    set of metrics computed over precursors quantified in AT LEAST K of the
    6 ProteoBench replicates. So the field we want is:

        entry['results'][str(min_replicates)]['nr_prec']

    Returns None when the field is missing or not numeric. Falls back to
    the top-level `nr_prec` only when `min_replicates == 1` (the top-level
    nr_prec is defined as precursors quantified in ≥1 of 6 samples).

    Robbe Devreese (ProteoBench maintainer) flagged that DIA-NN 1.9.1's
    headline nr_prec is unusually high at the ≥1 threshold but drops back
    in line at ≥3; this function is the canonical extractor for the ≥3
    re-ranking used in the corrected supp figure.
    """
    results = entry.get("results")
    key = str(min_replicates)
    if isinstance(results, dict) and key in results:
        bucket = results[key]
        if isinstance(bucket, dict):
            v = bucket.get("nr_prec")
            if isinstance(v, (int, float)) and not pd.isna(v):
                return int(v)
    if min_replicates == 1:
        v = entry.get("nr_prec")
        if isinstance(v, (int, float)) and not pd.isna(v):
            return int(v)
    return None


def parse_proteobench_datapoints_at_threshold(
    json_path: Path,
    min_replicates: int,
) -> Iterator[tuple[str, str, int, str]]:
    """Same shape as `parse_proteobench_datapoints` but reads precursor
    counts from `entry['results'][str(min_replicates)]['nr_prec']`.

    Skips submissions that have no bucket at the requested threshold (rare;
    every DIA-NN submission inspected during the Slack-driven correction
    review carries all six replicate buckets, but defensive handling here
    keeps the renderer robust against any future submission that doesn't)."""
    with open(json_path, encoding="utf-8") as fh:
        payload = json.load(fh)
    for entry in payload:
        nr_prec = extract_nr_prec_at_replicate_threshold(
            entry, min_replicates,
        )
        if nr_prec is None:
            continue
        software = entry.get("software_name") or ""
        version = entry.get("software_version") or ""
        if normalise_software_name(str(software)) == "dia-nn":
            kind = classify_predictors_library(entry.get("predictors_library"))
        else:
            kind = LIBRARY_KIND_OTHER_TOOL
        yield str(software), str(version), int(nr_prec), kind


def count_pr_matrix_min_replicates(
    matrix_path: Path,
    min_replicates: int,
) -> int:
    """Count precursors in a DIA-NN pr_matrix.tsv that have a non-NA
    intensity in AT LEAST `min_replicates` of the per-run sample columns.

    Sample columns are every column past the standard DIA-NN metadata block
    (everything that is not one of the well-known metadata column names).
    The four ProteoBench DIA pr_matrix.tsv files always have exactly six
    sample columns (Condition_A_REP1..3 + Condition_B_REP1..3), so a
    `min_replicates == 3` count matches Robbe's "≥3 replicate observations"
    bucket for the head-to-head supp figure.
    """
    metadata = {
        "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
        "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
        "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
    }
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    sample_cols = [c for c in df.columns if c not in metadata]
    if not sample_cols:
        return 0
    non_na = df[sample_cols].notna().sum(axis=1)
    return int((non_na >= min_replicates).sum())


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
# Parameter-signature extraction + matching (Step 1 / Step 2)
# ---------------------------------------------------------------------------
#
# The goal of this section is to do parameter-matched comparisons between
# quantmsdiann's DIA-NN invocations (recovered from the `diann ...` command
# line at the top of each `diannsummary.log`) and the parameter set
# published with each ProteoBench submission JSON. The matching surface
# follows the fields the ProteoBench DIA-NN datapoints carry:
#
#   software_name, software_version, predictors_library,
#   quantification_method, protein_inference, enable_match_between_runs,
#   ident_fdr_psm, fixed_mods, variable_mods.
#
# Three buckets are emitted by `param_match_category`:
#   - "exact":  same DIA-NN version + same library + same quant method +
#               same canonicalised fixed/variable mod set
#   - "near":   same DIA-NN version + 1-2 categorical mismatches
#   - "far":    different software OR DIA-NN version-major mismatch OR
#               3+ categorical mismatches


_DIANN_CLI_RE = re.compile(r"^diann\s+(.*)$", re.MULTILINE)
_DIANN_HEADER_VERSION_RE = re.compile(
    r"^DIA-NN\s+([0-9][0-9.]*)", re.MULTILINE
)


def _split_diann_cli(cli: str) -> list[str]:
    """Split a DIA-NN command-line string into tokens. The diannsummary log
    embeds the literal shell-style cmd line with whitespace-separated
    arguments; there are no quoted strings or shell escapes inside it, so a
    plain split is sufficient (and avoids depending on shlex's quoting
    semantics)."""
    return cli.strip().split()


def _parse_diann_mod_arg(arg: str) -> str:
    """Canonicalise a DIA-NN `--fixed-mod` / `--var-mod` argument value into
    'name@site' lowercase form. DIA-NN's cmd-line syntax is
    `name,delta_mass,site` (e.g. `Carbamidomethyl,57.021464,C` or
    `Acetyl,42.010565,*n`). The `*n` site refers to the protein N-term and
    is preserved verbatim."""
    parts = arg.split(",")
    if len(parts) >= 3:
        name = parts[0].strip().lower()
        site = parts[2].strip().lower()
        return f"{name}@{site}"
    return arg.strip().lower()


# Maps UniMod accession -> mod name (lower-case) so we can canonicalise
# mixed ProteoBench mod spellings: `unimod4`, `UniMod:4`, `Carbamidomethyl
# (C)` all collapse to `carbamidomethyl@c`. Only the mods that appear in
# the four ProteoBench DIA modules are listed; an unknown accession falls
# through to a generic `unimod:N@site` token so equality comparisons stay
# stable.
_UNIMOD_NAMES = {
    "1": "acetyl",
    "4": "carbamidomethyl",
    "21": "phospho",
    "35": "oxidation",
    "121": "ggl",
}

# Default DIA-NN site assumption per UniMod accession when ProteoBench
# records a bare `UniMod:N` token with no explicit site (some legacy DIA-NN
# 1.9.2 submissions encode `fixed_mods='UniMod:4'`, `variable_mods='UniMod:35,UniMod:1'`).
# Carbamidomethyl on C and Oxidation on M are the DIA-NN defaults; Acetyl
# is N-terminal by default in DIA-NN.
_UNIMOD_DEFAULT_SITES = {
    "1": "*n",
    "4": "c",
    "35": "m",
    "21": "sty",
    "121": "k",
}


def _canonicalise_proteobench_mod_token(token: str) -> str | None:
    """Canonicalise a single ProteoBench mod token into `name@site`
    lower-case form, or return None if the token is empty / unparseable.

    Handles the four spelling families observed across the four modules:
      - `UniMod:35/15.994915/M`        — accession / delta / site
      - `unimod4`, `UniMod:4`          — accession only (site defaults)
      - `Carbamidomethyl (C)`,
        `Oxidation (M)`                — name + parenthesised site
      - `UniMod:35`, `UniMod:1`        — accession only (defaults)
    """
    s = token.strip()
    if not s:
        return None
    s_lower = s.lower()

    # Form 1: UniMod:N/delta/site  e.g. UniMod:35/15.994915/M
    m = re.match(r"unimod[:_]?(\d+)\s*/[^/]+/\s*([a-z*]+)", s_lower)
    if m:
        acc, site = m.group(1), m.group(2)
        name = _UNIMOD_NAMES.get(acc, f"unimod:{acc}")
        return f"{name}@{site}"

    # Form 2: name (site)  e.g. Carbamidomethyl (C), Oxidation (M)
    m = re.match(r"([a-z]+)\s*\(([a-z\-]+(?:[\s_-]term)?)\)", s_lower)
    if m:
        name, site_raw = m.group(1), m.group(2)
        # Normalise "n-term" / "protein n-term" / "Protein_N-term" → "*n"
        if "n" in site_raw and "term" in site_raw:
            site = "*n"
        elif "c" in site_raw and "term" in site_raw:
            site = "*c"
        else:
            site = site_raw
        return f"{name}@{site}"

    # Form 2b: name@SITE  e.g. Carbamidomethyl@C, Acetyl@Protein_N-term
    m = re.match(r"([a-z]+)\s*@\s*([a-z_\- ]+)", s_lower)
    if m:
        name, site_raw = m.group(1), m.group(2).strip()
        if "n" in site_raw and "term" in site_raw:
            site = "*n"
        elif "c" in site_raw and "term" in site_raw:
            site = "*c"
        else:
            # Collapse e.g. "protein_n-term" with no 'term' marker; pick last
            # alpha run as site letter.
            cleaned = re.sub(r"[^a-z*]+", "", site_raw)
            site = cleaned or site_raw
        return f"{name}@{site}"

    # Form 3: UniMod:N alone, or unimodN
    m = re.match(r"unimod[:_]?(\d+)\s*$", s_lower)
    if m:
        acc = m.group(1)
        name = _UNIMOD_NAMES.get(acc, f"unimod:{acc}")
        site = _UNIMOD_DEFAULT_SITES.get(acc, "?")
        return f"{name}@{site}"

    return s_lower


def _parse_proteobench_mods(blob) -> frozenset[str]:
    """ProteoBench mod fields are comma- or semicolon-separated strings.
    Returns a frozenset of canonical `name@site` tokens; the empty string
    yields an empty set (some submissions record fdr=0.1 + empty mods,
    which we want to compare as 'no mods declared')."""
    if blob is None or (isinstance(blob, float) and pd.isna(blob)):
        return frozenset()
    s = str(blob).strip()
    if not s:
        return frozenset()
    tokens = [t for t in re.split(r"[,;]", s) if t.strip()]
    out: set[str] = set()
    for t in tokens:
        c = _canonicalise_proteobench_mod_token(t)
        if c:
            out.add(c)
    return frozenset(out)


def _normalise_version(version: str) -> str:
    """Strip ProteoBench's '<version> Academia' suffix and surrounding
    whitespace so that '2.5.0 Academia ' compares equal to '2.5.0'."""
    s = str(version).strip()
    s = re.sub(r"\s+Academia\s*$", "", s, flags=re.IGNORECASE)
    return s


def _version_tuple(version: str) -> tuple[int, ...]:
    """Best-effort numeric tuple for major.minor[.patch] comparisons.
    Non-numeric tokens are dropped; an empty tuple compares falsey."""
    s = _normalise_version(version)
    parts = re.findall(r"\d+", s)
    return tuple(int(p) for p in parts)


def extract_quantmsdiann_param_signature(log_path: Path) -> dict:
    """Parse a quantmsdiann `diannsummary.log` into the same parameter
    fingerprint shape ProteoBench publishes. Reads the DIA-NN command line
    (the single `diann ...` line emitted at the top of the log) and the
    `DIA-NN <version>` header line.

    Returns a dict with these canonical keys:
      software_name='DIA-NN', software_version (no 'Academia' suffix),
      predictors_library ('empirical' / 'predicted (DIANN)' / ...),
      quantification_method ('Legacy (direct)' / 'Legacy'),
      protein_inference (string of the --pg-level int),
      enable_match_between_runs (bool, --reanalyse),
      ident_fdr_psm (float, --qvalue),
      fixed_mods (frozenset of canonical tokens),
      variable_mods (frozenset of canonical tokens).
    """
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    cli_match = _DIANN_CLI_RE.search(text)
    if not cli_match:
        raise ValueError(
            f"No 'diann ...' command line found in {log_path}"
        )
    cli = cli_match.group(1)
    tokens = _split_diann_cli(cli)

    version_match = _DIANN_HEADER_VERSION_RE.search(text)
    version = _normalise_version(
        version_match.group(1) if version_match else ""
    )

    fixed_mods: set[str] = set()
    variable_mods: set[str] = set()
    lib_value: str | None = None
    pg_level: str | None = None
    qvalue: float | None = None
    direct_quant = False
    mbr = False
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        nxt = tokens[i + 1] if i + 1 < len(tokens) else None
        if tok == "--lib" and nxt is not None:
            lib_value = nxt
            i += 2
            continue
        if tok == "--fixed-mod" and nxt is not None:
            fixed_mods.add(_parse_diann_mod_arg(nxt))
            i += 2
            continue
        if tok == "--var-mod" and nxt is not None:
            variable_mods.add(_parse_diann_mod_arg(nxt))
            i += 2
            continue
        if tok == "--pg-level" and nxt is not None:
            pg_level = nxt
            i += 2
            continue
        if tok == "--qvalue" and nxt is not None:
            try:
                qvalue = float(nxt)
            except ValueError:
                pass
            i += 2
            continue
        if tok == "--direct-quant":
            direct_quant = True
            i += 1
            continue
        if tok == "--reanalyse":
            mbr = True
            i += 1
            continue
        i += 1

    # quantmsdiann's empirical library is built two-pass from the data;
    # detect the marker substring rather than the exact file extension
    # (1.8.1 uses .speclib, 2.x uses .parquet).
    if lib_value and "empirical_library" in lib_value:
        library_kind = LIBRARY_KIND_EMPIRICAL
    else:
        # Defensive fallback; the four benchmark workflows always pass
        # empirical_library.{speclib,parquet}.
        library_kind = lib_value or LIBRARY_KIND_EMPIRICAL

    # Both `--direct-quant` (2.1.0+) and the pre-2.1.0 absence-of-flag map
    # to DIA-NN's classic / Legacy quantification family, which ProteoBench
    # publishes under the 'Legacy' or 'Legacy (direct)' labels. We use the
    # 'Legacy (direct)' label whenever --direct-quant is present so that
    # comparisons against ProteoBench's 'Legacy (direct)' submissions hit
    # the 'exact' bucket; 1.8.1 (no flag) maps to 'Legacy'.
    quant_method = "Legacy (direct)" if direct_quant else "Legacy"

    return {
        "software_name": "DIA-NN",
        "software_version": version,
        "predictors_library": library_kind,
        "quantification_method": quant_method,
        "protein_inference": pg_level if pg_level is not None else "",
        "enable_match_between_runs": bool(mbr),
        "ident_fdr_psm": qvalue if qvalue is not None else 0.01,
        "fixed_mods": frozenset(fixed_mods),
        "variable_mods": frozenset(variable_mods),
    }


def extract_proteobench_param_signature(entry: dict) -> dict:
    """Project a ProteoBench submission JSON entry into the same canonical
    parameter-signature dict shape as `extract_quantmsdiann_param_signature`.
    The Boolean / numeric coercions are defensive: ProteoBench stores
    `enable_match_between_runs` as a Python bool literal in JSON but some
    older submissions store it as `'True'`/`'False'` strings, and
    `ident_fdr_psm` is occasionally `NaN` for legacy submissions."""
    software_name = str(entry.get("software_name") or "").strip()
    version = _normalise_version(entry.get("software_version") or "")
    library_kind = classify_predictors_library(entry.get("predictors_library"))
    quant_method = str(entry.get("quantification_method") or "").strip()
    prot_inf = entry.get("protein_inference")
    prot_inf_s = "" if prot_inf is None else str(prot_inf).strip()

    mbr_raw = entry.get("enable_match_between_runs")
    if isinstance(mbr_raw, bool):
        mbr = mbr_raw
    elif isinstance(mbr_raw, str):
        mbr = mbr_raw.strip().lower() in {"true", "1", "yes"}
    else:
        mbr = False

    fdr_raw = entry.get("ident_fdr_psm")
    try:
        fdr = float(fdr_raw) if fdr_raw is not None else float("nan")
    except (TypeError, ValueError):
        fdr = float("nan")

    return {
        "software_name": software_name,
        "software_version": version,
        "predictors_library": library_kind,
        "quantification_method": quant_method,
        "protein_inference": prot_inf_s,
        "enable_match_between_runs": mbr,
        "ident_fdr_psm": fdr,
        "fixed_mods": _parse_proteobench_mods(entry.get("fixed_mods")),
        "variable_mods": _parse_proteobench_mods(entry.get("variable_mods")),
    }


def param_match_category(qm_sig: dict, pb_sig: dict) -> str:
    """Categorise a single ProteoBench submission relative to a
    quantmsdiann signature into one of three bins:

      - 'exact': same DIA-NN version + same predictors_library + same
                 quantification_method + same fixed_mods + same
                 variable_mods. (FDR + protein_inference are tolerated
                 if equal modulo numeric form; MBR equality is required
                 because it materially changes the precursor count.)
      - 'near':  same software_name (DIA-NN) + same major.minor version
                 + 1-2 categorical mismatches across the five signature
                 fields above.
      - 'far':   different software (not DIA-NN), OR major-version
                 mismatch, OR 3+ categorical mismatches.

    The function tolerates ProteoBench's 'Academia' version suffix and
    case differences in the categorical strings; mod-set equality uses
    the canonicalised `name@site` token sets produced by
    `_parse_proteobench_mods` / `_parse_diann_mod_arg`.
    """
    if normalise_software_name(pb_sig.get("software_name", "")) != "dia-nn":
        return "far"
    qm_v = _version_tuple(qm_sig["software_version"])
    pb_v = _version_tuple(pb_sig["software_version"])
    # Major-version mismatch is a hard 'far' (e.g. quantmsdiann 1.8.1 vs PB
    # 2.5.0 — even if every other categorical lines up, the binary is a
    # different generation of DIA-NN).
    if qm_v and pb_v and qm_v[0] != pb_v[0]:
        return "far"

    same_version = (
        qm_sig["software_version"] == pb_sig["software_version"]
        and qm_sig["software_version"] != ""
    )

    mismatches = 0
    if qm_sig["predictors_library"] != pb_sig["predictors_library"]:
        mismatches += 1
    if (qm_sig["quantification_method"].lower()
            != pb_sig["quantification_method"].lower()):
        mismatches += 1
    if qm_sig["fixed_mods"] != pb_sig["fixed_mods"]:
        mismatches += 1
    if qm_sig["variable_mods"] != pb_sig["variable_mods"]:
        mismatches += 1
    if (bool(qm_sig["enable_match_between_runs"])
            != bool(pb_sig["enable_match_between_runs"])):
        mismatches += 1

    if same_version and mismatches == 0:
        return "exact"
    # Allow 'near' when the major.minor version aligns (e.g. quantmsdiann
    # 2.5.0 vs PB 2.5.0 Academia → same after suffix strip; quantmsdiann
    # 2.3.2 vs PB 2.3.0 share major.minor 2.3.x and count as 'near').
    qm_majmin = qm_v[:2] if len(qm_v) >= 2 else qm_v
    pb_majmin = pb_v[:2] if len(pb_v) >= 2 else pb_v
    if qm_majmin and pb_majmin and qm_majmin == pb_majmin and mismatches <= 2:
        return "near"
    if same_version and mismatches <= 2:
        return "near"
    return "far"


# ---------------------------------------------------------------------------
# Step 1 — DIA-NN parity panel
# ---------------------------------------------------------------------------


def render_diann_parity(
    quantmsdiann_rows_by_threshold: dict[int, list[tuple[str, str, int, int]]],
    pb_rows_by_threshold: dict[int, dict[str, list[dict]]],
    qm_signatures: dict[tuple[str, str], dict],
    pb_signatures_by_dataset: dict[str, list[tuple[dict, int, dict]]],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Step 1 figure: per dataset, plot quantmsdiann's `nr_prec` at ≥1 and
    ≥3 replicates alongside the matched ProteoBench DIA-NN submissions'
    same metric.

    "Matched" = `param_match_category(qm, pb) == 'exact'`. If no exact
    matches exist for a dataset, the 'near' set is used as a fall-back
    so the panel is non-empty (annotated as such on the figure).

    Returns (parity_long_df, epsilon_df) where:
      - parity_long_df has columns
          [dataset, version, threshold, source, label,
           precursors, match_category]
      - epsilon_df has columns
          [dataset, threshold, n_matched, qm_median, matched_median,
           epsilon_frac, match_level]
    """
    datasets = sorted(
        quantmsdiann_rows_by_threshold[1] and {
            r[0] for r in quantmsdiann_rows_by_threshold[1]
        },
        key=_dataset_sort_key,
    )
    rows: list[dict] = []
    eps_rows: list[dict] = []
    for dataset in datasets:
        for threshold in (1, 3):
            for ds, version, prec, _proteins in quantmsdiann_rows_by_threshold[threshold]:
                if ds != dataset:
                    continue
                rows.append({
                    "dataset": dataset,
                    "threshold": threshold,
                    "source": "quantmsdiann",
                    "version": version,
                    "label": f"quantmsdiann {_VERSION_LABELS.get(version, version)}",
                    "precursors": prec,
                    "match_category": "self",
                })
            # Locate the matched PB submissions for this dataset by
            # cross-referencing each pb entry's full signature against
            # every quantmsdiann signature for the same dataset.
            qm_sigs = [qm_signatures.get((dataset, v)) for v in DIANN_VERSIONS]
            qm_sigs = [s for s in qm_sigs if s is not None]
            matched_entries: list[tuple[dict, int, str]] = []  # (pb_sig, nr_prec, match_level)
            for pb_sig, pb_prec, raw_entry in pb_signatures_by_dataset.get(dataset, []):
                # Skip submissions that don't carry a value for this
                # replicate threshold; otherwise they'd masquerade as zero.
                pb_thr_value = extract_nr_prec_at_replicate_threshold(
                    raw_entry, threshold,
                )
                if pb_thr_value is None:
                    continue
                best_cat = "far"
                for qm_sig in qm_sigs:
                    cat = param_match_category(qm_sig, pb_sig)
                    if cat == "exact":
                        best_cat = "exact"
                        break
                    if cat == "near" and best_cat == "far":
                        best_cat = "near"
                if best_cat in {"exact", "near"}:
                    matched_entries.append((pb_sig, pb_thr_value, best_cat))

            exact_only = [m for m in matched_entries if m[2] == "exact"]
            cohort = exact_only if exact_only else matched_entries
            match_level = (
                "exact" if exact_only
                else ("near" if matched_entries else "none")
            )
            for pb_sig, pb_prec, cat in cohort:
                rows.append({
                    "dataset": dataset,
                    "threshold": threshold,
                    "source": "proteobench-matched",
                    "version": pb_sig["software_version"],
                    "label": f"DIA-NN {pb_sig['software_version']}".strip(),
                    "precursors": pb_prec,
                    "match_category": cat,
                })

            # Per-dataset epsilon: relative gap between quantmsdiann's
            # median (across the five versions at this threshold) and the
            # matched PB cohort's median.
            qm_prec_vals = [
                r[2] for r in quantmsdiann_rows_by_threshold[threshold]
                if r[0] == dataset
            ]
            if cohort and qm_prec_vals:
                qm_med = float(pd.Series(qm_prec_vals).median())
                pb_med = float(pd.Series([m[1] for m in cohort]).median())
                epsilon_frac = (
                    abs(qm_med - pb_med) / pb_med if pb_med else float("nan")
                )
            else:
                qm_med = float(pd.Series(qm_prec_vals).median()) if qm_prec_vals else float("nan")
                pb_med = float("nan")
                epsilon_frac = float("nan")
            eps_rows.append({
                "dataset": dataset,
                "threshold": threshold,
                "n_matched": len(cohort),
                "qm_median": qm_med,
                "matched_median": pb_med,
                "epsilon_frac": epsilon_frac,
                "match_level": match_level,
            })

    parity_long_df = pd.DataFrame(rows)
    epsilon_df = pd.DataFrame(eps_rows)

    # Render: 4 datasets stacked vertically, 2 threshold columns side by side.
    fig, axes = plt.subplots(
        nrows=len(datasets), ncols=2,
        figsize=(11.5, 2.6 * len(datasets)),
        sharex=False, squeeze=False,
    )
    for i, dataset in enumerate(datasets):
        for j, threshold in enumerate((1, 3)):
            ax = axes[i, j]
            sub = parity_long_df[
                (parity_long_df["dataset"] == dataset)
                & (parity_long_df["threshold"] == threshold)
            ].copy()
            qm = sub[sub["source"] == "quantmsdiann"].sort_values(
                "version",
                key=lambda s: s.map(lambda v: DIANN_VERSIONS.index(v)),
            )
            pb = sub[sub["source"] == "proteobench-matched"].sort_values(
                "precursors", ascending=True,
            )
            labels: list[str] = []
            values: list[int] = []
            colours: list[str] = []
            for _, row in pb.iterrows():
                labels.append(row["label"])
                values.append(int(row["precursors"]))
                # Teal for 'exact', light blue for 'near'.
                colours.append(
                    "#26a69a" if row["match_category"] == "exact"
                    else "#90caf9"
                )
            for _, row in qm.iterrows():
                labels.append(row["label"])
                values.append(int(row["precursors"]))
                colours.append("#d62728")
            ys = list(range(len(labels)))
            ax.barh(ys, values, color=colours, height=0.7)
            ax.set_yticks(ys)
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel("Precursors quantified" if i == len(datasets) - 1 else "")
            thr_label = "≥1 replicate" if threshold == 1 else "≥3 replicates"
            ax.set_title(
                f"{_dataset_display_label(dataset).splitlines()[0]} — {thr_label}",
                loc="left", fontsize=9, fontweight="bold",
            )
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(axis="x", labelsize=8)
            xmax = max(values) if values else 1
            ax.set_xlim(0, xmax * 1.18)
            for y, v in zip(ys, values):
                ax.text(v + xmax * 0.005, y, f"{v:,}", va="center",
                        ha="left", fontsize=6, color="#404040")
            # Epsilon annotation in the top-right.
            eps_sub = epsilon_df[
                (epsilon_df["dataset"] == dataset)
                & (epsilon_df["threshold"] == threshold)
            ]
            if len(eps_sub):
                eps_row = eps_sub.iloc[0]
                if pd.notna(eps_row["epsilon_frac"]):
                    ax.text(
                        0.98, 0.05,
                        f"ε = {eps_row['epsilon_frac']*100:.1f}%   "
                        f"n={int(eps_row['n_matched'])} ({eps_row['match_level']})",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=7, color="#444444",
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),
                    )
                else:
                    ax.text(
                        0.98, 0.05,
                        "no parameter-matched DIA-NN submission",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=7, color="#888888", fontstyle="italic",
                    )
    # Single figure-level legend.
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor="#d62728", label="quantmsdiann (DIA-NN, empirical lib)"),
        Patch(facecolor="#26a69a", label="ProteoBench DIA-NN, exact param match"),
        Patch(facecolor="#90caf9", label="ProteoBench DIA-NN, near param match"),
    ]
    fig.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0),
        ncol=3, frameon=False, fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)
    return parity_long_df, epsilon_df


# ---------------------------------------------------------------------------
# Step 2 — match-then-compare supp
# ---------------------------------------------------------------------------


def render_vs_proteobench_matched(
    long_df: pd.DataFrame,
    match_categories: dict[tuple[str, str, str, int], str],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """Same horizontal-bar layout as `render_vs_proteobench`, but bars are
    colour-coded by match category. `match_categories` maps
    `(dataset, tool, version, precursors)` to one of {'exact','near','far'};
    bars that aren't in the map fall back to 'far' (treated as the neutral
    cohort)."""
    datasets = sorted(long_df["dataset"].unique(), key=_dataset_sort_key)
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
    cat_palette = {
        "exact": "#26a69a",
        "near": "#90caf9",
        "far": "#bdbdbd",
    }
    for i, dataset in enumerate(datasets):
        ax = axes[i, 0]
        sub = long_df[long_df["dataset"] == dataset].copy()
        pb = (sub[sub["source"] == "proteobench"]
              .sort_values("precursors", ascending=True)
              .reset_index(drop=True))
        labels = [f"{t} {v}".strip() for t, v in zip(pb["tool"], pb["version"])]
        cats = [
            match_categories.get(
                (dataset, t, v, int(p)), "far"
            )
            for t, v, p in zip(pb["tool"], pb["version"], pb["precursors"])
        ]
        bar_colours = [cat_palette[c] for c in cats]
        ax.barh(range(len(pb)), pb["precursors"], color=bar_colours, height=0.7)
        xmax_pb = float(pb["precursors"].max()) if len(pb) else 0.0
        for k, (precs, lab, cat) in enumerate(zip(pb["precursors"], labels, cats)):
            ax.text(
                precs + xmax_pb * 0.005, k, f"{lab}  [{cat}]",
                va="center", ha="left", fontsize=6, color="#404040",
            )
        ax.set_yticks([])
        qm = sub[sub["source"] == "quantmsdiann"].sort_values(
            "version", key=lambda s: s.map(lambda v: DIANN_VERSIONS.index(v)),
        )
        qm_ys = [-1.0 - k for k in range(len(qm))]
        qm_vals = list(qm["precursors"])
        qm_labels = [
            f"quantmsdiann {_VERSION_LABELS.get(v, v)}"
            for v in qm["version"]
        ]
        xmax_qm = 0.0
        if qm_ys:
            top = qm_ys[0] + 0.55
            bottom = qm_ys[-1] - 0.55
            ax.axhspan(bottom, top, color="#ffebee", zorder=0)
            ax.barh(
                qm_ys, qm_vals, color="#d62728", height=0.7,
                edgecolor="#7f1d1d", linewidth=0.6, zorder=2,
            )
            xmax_qm = max(qm_vals) if qm_vals else 0.0
            for y, v, lab in zip(qm_ys, qm_vals, qm_labels):
                ax.text(
                    v + max(xmax_pb, xmax_qm) * 0.005, y, lab,
                    va="center", ha="left", fontsize=7,
                    color="#7f1d1d", fontweight="bold",
                )
        ax.set_ylim(min(qm_ys) - 1.0 if qm_ys else -1.5, len(pb))
        ax.set_xlim(0, max(xmax_pb, xmax_qm) * 1.3)
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
    handles = [
        Patch(facecolor="#d62728", edgecolor="#7f1d1d", linewidth=0.8,
              label="quantmsdiann (DIA-NN, empirical library)"),
        Patch(facecolor=cat_palette["exact"], label="ProteoBench — exact param match"),
        Patch(facecolor=cat_palette["near"], label="ProteoBench — near param match"),
        Patch(facecolor=cat_palette["far"], label="ProteoBench — far / different stack"),
    ]
    fig.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0),
        ncol=2, frameon=False, fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def median_nr_prec_per_version(
    proteobench_rows_by_threshold: dict[int, dict[str, list[tuple[str, str, int, str]]]],
    quantmsdiann_rows_by_threshold: dict[int, list[tuple[str, str, int, int]]] | None = None,
) -> pd.DataFrame:
    """Build a long-format table of median precursor counts per DIA-NN
    version per module at each replicate threshold. Used to quantify the
    Robbe-predicted 1.9.1 vs 2.3.0 ranking flip across the ≥1 and ≥3
    thresholds.

    `proteobench_rows_by_threshold` is keyed by min_replicates -> dataset ->
    list of (tool, version, nr_prec, library_kind). Optional
    `quantmsdiann_rows_by_threshold` adds quantmsdiann rows for the same
    table (one quantmsdiann analysis per version, so 'median' == the single
    value)."""
    rows = []
    for thr, by_dataset in proteobench_rows_by_threshold.items():
        for dataset, entries in by_dataset.items():
            diann = [(v, n) for tool, v, n, _k in entries
                     if normalise_software_name(tool) == "dia-nn"]
            by_v: dict[str, list[int]] = {}
            for v, n in diann:
                by_v.setdefault(v.strip(), []).append(n)
            for v, ns in sorted(by_v.items()):
                rows.append({
                    "dataset": dataset,
                    "min_replicates": thr,
                    "source": "proteobench",
                    "version": v,
                    "n_submissions": len(ns),
                    "median_nr_prec": int(pd.Series(ns).median()),
                })
    if quantmsdiann_rows_by_threshold is not None:
        for thr, qm_rows in quantmsdiann_rows_by_threshold.items():
            for dataset, version, nr_prec, _ in qm_rows:
                rows.append({
                    "dataset": dataset,
                    "min_replicates": thr,
                    "source": "quantmsdiann",
                    "version": _VERSION_LABELS.get(version, version),
                    "n_submissions": 1,
                    "median_nr_prec": int(nr_prec),
                })
    return pd.DataFrame(rows)


def write_median_table(df: pd.DataFrame, tsv_path: Path) -> None:
    """Auditable TSV of median precursor counts per DIA-NN version per module
    at the ≥1 and ≥3 replicate thresholds. One row per
    (dataset, threshold, source, version)."""
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    out = df.sort_values(
        ["dataset", "min_replicates", "source", "version"]
    )
    out.to_csv(tsv_path, sep="\t", index=False)


def main() -> int:  # pragma: no cover
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # quantmsdiann rows for the headline main panels (≥1-replicate row count,
    # i.e. the raw pr_matrix.tsv row count). Each row is
    # (dataset, version, precursors_min1, proteins).
    quantmsdiann_rows: list[tuple[str, str, int, int]] = []
    # Parallel ≥3-replicate quantmsdiann counts, computed by re-scanning the
    # cached pr_matrix.tsv. (proteins unchanged.)
    quantmsdiann_rows_min3: list[tuple[str, str, int, int]] = []
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
            precursors_min3 = count_pr_matrix_min_replicates(pr_path, 3)
            proteins = count_matrix_data_rows(pg_path)
            quantmsdiann_rows.append((dataset, version, precursors, proteins))
            quantmsdiann_rows_min3.append(
                (dataset, version, precursors_min3, proteins)
            )
            print(f"{dataset} {version}: precursors_min1={precursors:,}  "
                  f"precursors_min3={precursors_min3:,}  "
                  f"proteins={proteins:,}")

    proteobench_rows: dict[str, list[tuple[str, str, int, str]]] = {}
    proteobench_rows_min3: dict[str, list[tuple[str, str, int, str]]] = {}
    for dataset, info in DATASET_TO_MODULE.items():
        cache = DATA_DIR / "proteobench" / f"{dataset}.json"
        try:
            fetch_proteobench_module(info["results_repo"], cache)
        except Exception as exc:
            print(f"WARN: failed to fetch ProteoBench for {dataset}: {exc}",
                  file=sys.stderr)
            proteobench_rows[dataset] = []
            proteobench_rows_min3[dataset] = []
            continue
        entries_min1 = list(parse_proteobench_datapoints_at_threshold(cache, 1))
        entries_min3 = list(parse_proteobench_datapoints_at_threshold(cache, 3))
        proteobench_rows[dataset] = entries_min1
        proteobench_rows_min3[dataset] = entries_min3
        print(f"{dataset}: {len(entries_min1)} ProteoBench submissions "
              f"(min1) / {len(entries_min3)} (min3)")

    long_df = build_long_table(quantmsdiann_rows, proteobench_rows)
    long_df_min3 = build_long_table(
        quantmsdiann_rows_min3, proteobench_rows_min3,
    )

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

    print("Rendering ProteoBench-overlay supp figure (≥1 replicate)...")
    render_vs_proteobench(
        long_df,
        FIGURES_DIR / "supp_vs_proteobench_min1.pdf",
        FIGURES_DIR / "supp_vs_proteobench_min1.png",
        FIGURES_DIR / "supp_vs_proteobench_min1.svg",
    )

    print("Rendering ProteoBench-overlay supp figure (≥3 replicates, "
          "Slack-corrected default)...")
    render_vs_proteobench(
        long_df_min3,
        FIGURES_DIR / "supp_vs_proteobench_min3.pdf",
        FIGURES_DIR / "supp_vs_proteobench_min3.png",
        FIGURES_DIR / "supp_vs_proteobench_min3.svg",
    )

    print("Writing auditable counts TSV (≥1)...")
    write_counts_tsv(long_df, FIGURES_DIR / "counts.tsv")
    print("Writing auditable counts TSV (≥3)...")
    write_counts_tsv(long_df_min3, FIGURES_DIR / "counts_min3.tsv")

    print("Writing per-DIA-NN-version median precursor table (≥1 and ≥3)...")
    median_df = median_nr_prec_per_version(
        {
            1: proteobench_rows,
            3: proteobench_rows_min3,
        },
        {
            1: quantmsdiann_rows,
            3: quantmsdiann_rows_min3,
        },
    )
    write_median_table(
        median_df, FIGURES_DIR / "median_nr_prec_by_version.tsv",
    )

    # ---------------------------------------------------------------------
    # Step 1 / Step 2 — parameter-signature comparison
    # ---------------------------------------------------------------------
    print("Building quantmsdiann parameter signatures from diannsummary.log...")
    qm_signatures: dict[tuple[str, str], dict] = {}
    for dataset in DATASET_TO_MODULE:
        for version in DIANN_VERSIONS:
            log_path = (DATA_DIR / dataset / version / "quant_tables"
                        / "diannsummary.log")
            if not log_path.exists():
                print(f"  skip {dataset}/{version}: no diannsummary.log",
                      file=sys.stderr)
                continue
            try:
                qm_signatures[(dataset, version)] = (
                    extract_quantmsdiann_param_signature(log_path)
                )
            except Exception as exc:
                print(f"  WARN: failed to parse {log_path}: {exc}",
                      file=sys.stderr)

    print("Building ProteoBench parameter signatures from cached JSONs...")
    # Per-dataset list of (pb_signature, top_level_nr_prec, raw_entry); the
    # raw entry is kept so we can pull per-threshold nr_prec on demand.
    pb_signatures_by_dataset: dict[str, list[tuple[dict, int, dict]]] = {}
    for dataset, info in DATASET_TO_MODULE.items():
        cache = DATA_DIR / "proteobench" / f"{dataset}.json"
        if not cache.exists():
            pb_signatures_by_dataset[dataset] = []
            continue
        with open(cache, encoding="utf-8") as fh:
            payload = json.load(fh)
        out_list: list[tuple[dict, int, dict]] = []
        for entry in payload:
            top_nrprec = entry.get("nr_prec")
            if not isinstance(top_nrprec, (int, float)) or pd.isna(top_nrprec):
                top_nrprec = 0
            try:
                sig = extract_proteobench_param_signature(entry)
            except Exception as exc:
                print(f"  WARN: failed to parse pb entry: {exc}",
                      file=sys.stderr)
                continue
            out_list.append((sig, int(top_nrprec), entry))
        pb_signatures_by_dataset[dataset] = out_list

    # Per-(dataset,threshold) match counts for the spec doc.
    match_counts: dict[tuple[str, int], dict[str, int]] = {}
    for dataset, sig_rows in pb_signatures_by_dataset.items():
        qm_sigs_for_dataset = [
            qm_signatures.get((dataset, v)) for v in DIANN_VERSIONS
        ]
        qm_sigs_for_dataset = [s for s in qm_sigs_for_dataset if s is not None]
        counts = {"exact": 0, "near": 0, "far": 0}
        for pb_sig, _prec, _entry in sig_rows:
            best = "far"
            for qm_sig in qm_sigs_for_dataset:
                cat = param_match_category(qm_sig, pb_sig)
                if cat == "exact":
                    best = "exact"
                    break
                if cat == "near" and best == "far":
                    best = "near"
            counts[best] += 1
        match_counts[(dataset, 1)] = counts
        match_counts[(dataset, 3)] = counts  # category is replicate-independent

    print("Rendering Step 1 DIA-NN parity figure...")
    parity_df, epsilon_df = render_diann_parity(
        {1: quantmsdiann_rows, 3: quantmsdiann_rows_min3},
        {1: proteobench_rows, 3: proteobench_rows_min3},
        qm_signatures,
        pb_signatures_by_dataset,
        FIGURES_DIR / "main_diann_quantmsdiann_parity.pdf",
        FIGURES_DIR / "main_diann_quantmsdiann_parity.png",
        FIGURES_DIR / "main_diann_quantmsdiann_parity.svg",
    )
    epsilon_df.to_csv(
        FIGURES_DIR / "diann_quantmsdiann_parity_epsilon.tsv",
        sep="\t", index=False,
    )
    parity_df.to_csv(
        FIGURES_DIR / "diann_quantmsdiann_parity_long.tsv",
        sep="\t", index=False,
    )

    print("Rendering Step 2 match-then-compare supp figure (≥3 default)...")
    # Map every (dataset, tool, version, precursors) -> match category for
    # the bar-colour lookup inside render_vs_proteobench_matched. We build
    # the map against the ≥3 supp; the ≥1 supp uses the same categories.
    match_categories: dict[tuple[str, str, str, int], str] = {}
    for dataset, sig_rows in pb_signatures_by_dataset.items():
        qm_sigs_for_dataset = [
            qm_signatures.get((dataset, v)) for v in DIANN_VERSIONS
        ]
        qm_sigs_for_dataset = [s for s in qm_sigs_for_dataset if s is not None]
        # Iterate the cached entries in their declared order so the lookup
        # remains deterministic even if precursor counts collide.
        for pb_sig, _top_prec, entry in sig_rows:
            for thr in (1, 3):
                val = extract_nr_prec_at_replicate_threshold(entry, thr)
                if val is None:
                    continue
                best = "far"
                for qm_sig in qm_sigs_for_dataset:
                    cat = param_match_category(qm_sig, pb_sig)
                    if cat == "exact":
                        best = "exact"
                        break
                    if cat == "near" and best == "far":
                        best = "near"
                # The supp figure's `tool` column is the verbatim
                # ProteoBench `software_name`; use it as-is.
                tool = entry.get("software_name") or ""
                # Version stored in long_df is the verbatim PB version string
                # ('2.5.0 Academia '), not the normalised form, so we key on
                # that.
                pb_version_raw = entry.get("software_version") or ""
                match_categories[
                    (dataset, str(tool), str(pb_version_raw), int(val))
                ] = best

    render_vs_proteobench_matched(
        long_df_min3,
        match_categories,
        FIGURES_DIR / "supp_vs_proteobench_matched_min3.pdf",
        FIGURES_DIR / "supp_vs_proteobench_matched_min3.png",
        FIGURES_DIR / "supp_vs_proteobench_matched_min3.svg",
    )
    render_vs_proteobench_matched(
        long_df,
        match_categories,
        FIGURES_DIR / "supp_vs_proteobench_matched_min1.pdf",
        FIGURES_DIR / "supp_vs_proteobench_matched_min1.png",
        FIGURES_DIR / "supp_vs_proteobench_matched_min1.svg",
    )

    # Per-dataset match-count summary; printed + saved.
    summary_rows = []
    for dataset in DATASET_TO_MODULE:
        c = match_counts.get((dataset, 1), {"exact": 0, "near": 0, "far": 0})
        print(f"  {dataset}: exact={c['exact']} near={c['near']} far={c['far']}")
        summary_rows.append({
            "dataset": dataset,
            "exact": c["exact"],
            "near": c["near"],
            "far": c["far"],
            "total": c["exact"] + c["near"] + c["far"],
        })
    pd.DataFrame(summary_rows).to_csv(
        FIGURES_DIR / "match_category_counts.tsv", sep="\t", index=False,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
