"""quantmsdiann vs ProteoBench benchmarks figure.

Each of the four ProteoBench benchmark datasets has been processed through
quantmsdiann across DIA-NN versions (1.8.1, 2.5.1, 2.5.1-enterprise). This script:

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
from analysis import figure_style as fs
fs.apply_house_style()
import pandas as pd
import requests

from analysis.contaminant_filter import is_target_protein_group
from analysis.figure_original_vs_quantmsdiann import (
    download_if_missing,
    SUMMARY_LOG_PROTEIN_LINE_RE,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "quantmsdiann_benchmarks"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "quantmsdiann_benchmarks"
SUPP_DIR = FIGURES_DIR / "supplementary"
FIG_DATA_DIR = FIGURES_DIR / "data"

PRIDE_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/proteobench/quantmsdiann_results"
)

# The paper compares the oldest supported DIA-NN release (1.8.1) against the
# current release (2.5.1) in both its academic and enterprise builds. The
# enterprise build is the same 2.5.1 algorithm under a commercial licence, so
# it is rendered as a distinct third series rather than a "newer" version.
DIANN_VERSIONS = ("v1_8_1", "v2_5_1", "v2_5_1_enterprise")

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


def count_matrix_data_rows_split(matrix_path: Path) -> tuple[int, int]:
    """Return `(unfiltered_rows, target_only_rows)` from a DIA-NN
    pr_matrix.tsv / pg_matrix.tsv. Target-only drops rows whose
    `Protein.Group` carries a contaminant / entrapment / decoy prefix
    (see `analysis.contaminant_filter.is_target_protein_group`).

    Used by the benchmark counts: the headline number consumed by
    figures is the target-only count; the unfiltered count is kept in
    the audit TSV for reviewer comparison."""
    unfiltered = count_matrix_data_rows(matrix_path)
    try:
        df = pd.read_csv(
            matrix_path, sep="\t", usecols=["Protein.Group"], dtype=str,
        )
    except (FileNotFoundError, OSError, ValueError):
        return (unfiltered, unfiltered)
    pgs = df["Protein.Group"].dropna()
    target = int(pgs.map(is_target_protein_group).sum())
    return (unfiltered, target)


LIBRARY_KIND_EMPIRICAL = "empirical"
LIBRARY_KIND_PREDICTED = "predicted (DIANN)"
LIBRARY_KIND_USER_DEFINED = "user-defined speclib"
LIBRARY_KIND_OTHER_TOOL = "other tool"

# Library strategy quantmsdiann uses, expressed in ProteoBench's
# `predictors_library` taxonomy. quantmsdiann runs DIA-NN library-free:
# INSILICO_LIBRARY_GENERATION predicts a spectral library in-silico from the
# FASTA with DIA-NN's deep-learning predictor, then ASSEMBLE_EMPIRICAL_LIBRARY
# refines it from the first-pass IDs. There is NO externally supplied
# experimental library. In ProteoBench terms that is `predictors_library =
# {RT,IM,MS2 -> DIANN}` = "predicted (DIANN)" — NOT "empirical" (which marks a
# submission that loaded a pre-existing experimental library without a
# predictor). The `empirical_library.parquet` filename refers to DIA-NN's own
# refined-from-prediction library, not a user experimental library, so the
# apples-to-apples comparison set is the predicted-library DIA-NN submissions.
QUANTMSDIANN_LIBRARY_KIND = LIBRARY_KIND_PREDICTED


def classify_predictors_library(predictors_library) -> str:
    """Map a ProteoBench `predictors_library` value to one of three categories
    used to color DIA-NN community submissions:
      - 'empirical' (None / 'None') — the submission loaded a pre-existing
        experimental/curated spectral library WITHOUT a predictor. This does
        NOT match quantmsdiann, which predicts its library from the FASTA.
      - 'predicted (DIANN)' — DIA-NN's built-in in-silico prediction (the
        FASTA-derived predicted library). This is quantmsdiann's strategy
        (see QUANTMSDIANN_LIBRARY_KIND) and the apples-to-apples set.
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

    Returns the **target-only** count (rows whose Protein.Group passes
    `analysis.contaminant_filter.is_target_protein_group`). Callers that
    need the unfiltered baseline as well should use
    `count_pr_matrix_min_replicates_split`.
    """
    _unf, target = count_pr_matrix_min_replicates_split(
        matrix_path, min_replicates
    )
    return target


def count_pr_matrix_min_replicates_split(
    matrix_path: Path,
    min_replicates: int,
) -> tuple[int, int]:
    """Return `(unfiltered, target_only)` precursor counts at the
    ≥`min_replicates` threshold. Target-only drops rows whose
    Protein.Group carries a contaminant / entrapment / decoy prefix
    (`Cont_`, `CONTAM_`, `ENTRAP_`, `DECOY_`, `decoy_`) per the
    2026-05-21 conservative-filter spec.

    Used by the benchmark headline-count writers so the audit TSV can
    record both numbers and reviewers can verify the ~1 % contamination
    drop predicted by the spec."""
    metadata = {
        "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
        "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
        "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
    }
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    sample_cols = [c for c in df.columns if c not in metadata]
    if not sample_cols:
        return (0, 0)
    non_na = df[sample_cols].notna().sum(axis=1)
    passes = non_na >= min_replicates
    unfiltered = int(passes.sum())
    if "Protein.Group" not in df.columns:
        return (unfiltered, unfiltered)
    target_mask = passes & df["Protein.Group"].map(
        lambda v: is_target_protein_group(v) if isinstance(v, str) else False
    )
    target = int(target_mask.sum())
    return (unfiltered, target)


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
            # quantmsdiann predicts its library in-silico from the FASTA
            # (DIA-NN library-free); in ProteoBench's predictors_library
            # taxonomy that is "predicted (DIANN)". See QUANTMSDIANN_LIBRARY_KIND.
            "library_kind": QUANTMSDIANN_LIBRARY_KIND,
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
    "v2_5_1": "2.5.1",
    "v2_5_1_enterprise": "2.5.1 ent.",
}

# Per-version colour + marker, shared by every ProteoBench panel (precursors
# bars here and the accuracy scatter in figure_id_vs_epsilon). 1.8.1 -> 2.5.1
# read light -> dark blue (older -> newer release); the 2.5.1 enterprise build
# is the same algorithm under a commercial licence, so it gets a distinct
# amber accent instead of continuing the blue ramp.
# Colour-blind-safe DIA-NN version palette (centralised in figure_style):
# light->dark blue = older->newer release; enterprise = distinct accent hue.
_VERSION_COLORS = dict(fs.VERSION_COLORS)
_VERSION_MARKERS = {
    "v1_8_1": "o",
    "v2_5_1": "s",
    "v2_5_1_enterprise": "D",
}


def _render_main_metric(
    quantmsdiann_rows: list[tuple[str, str, int, int]],
    svg_path: Path,
    *,
    metric: str,
    ylabel: str,
    label_fmt,
) -> None:
    """Grouped-bar panel of `metric` (`precursors` or `proteins`) across the
    4 datasets x 3 DIA-NN versions. Paper-ready: no title, no footer."""
    df = pd.DataFrame(
        quantmsdiann_rows,
        columns=["dataset", "version", "precursors", "proteins"],
    )
    datasets = sorted(df["dataset"].unique(), key=_dataset_sort_key)
    fig, ax = plt.subplots(figsize=(10.5, 5))
    n_versions = len(DIANN_VERSIONS)
    bar_width = 0.8 / n_versions
    x = list(range(len(datasets)))
    for k, version in enumerate(DIANN_VERSIONS):
        vals = [
            int(df[(df["dataset"] == d) & (df["version"] == version)]
                [metric].iloc[0])
            if not df[(df["dataset"] == d) & (df["version"] == version)].empty
            else 0
            for d in datasets
        ]
        offsets = [xi + (k - (n_versions - 1) / 2) * bar_width for xi in x]
        bars = ax.bar(
            offsets, vals, width=bar_width,
            color=_VERSION_COLORS.get(version, "#1f77b4"),
            label=_VERSION_LABELS.get(version, version),
        )
        for bar, v in zip(bars, vals):
            if v == 0:
                continue
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    label_fmt(v), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [_dataset_display_label(d).replace("\n", "\n") for d in datasets],
        fontsize=9,
    )
    ax.set_ylabel(ylabel)
    ymax = df[metric].max() * 1.15
    ax.set_ylim(0, ymax)
    fs.kfmt_axis(ax.yaxis)  # match the bar labels' "134k" style on the ticks
    fs.despine(ax)
    ax.legend(
        title="DIA-NN version", loc="upper center",
        bbox_to_anchor=(0.5, -0.18), ncol=n_versions, frameon=False,
        fontsize=8, title_fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


def render_main_precursors(
    quantmsdiann_rows: list[tuple[str, str, int, int]],
    svg_path: Path,
) -> None:
    """Precursor-only headline: 4 datasets x 3 DIA-NN versions, grouped bars.
    Precursor counts come from the DIA-NN report (run-specific q at the
    per-version recommended cut-off + 1% global)."""
    _render_main_metric(
        quantmsdiann_rows, svg_path, metric="precursors",
        ylabel="Precursors quantified (1% global FDR)",
        label_fmt=lambda v: f"{v / 1000:.0f}k",
    )


def render_main_proteins(
    quantmsdiann_rows: list[tuple[str, str, int, int]],
    svg_path: Path,
) -> None:
    """Protein-group headline: 4 datasets x 3 DIA-NN versions, grouped bars.
    Protein groups are the per-run average, counted from the DIA-NN report at
    1% run-specific protein-group q-value (not the pg_matrix row count). This
    per-run depth metric is sensitive to the version improvement; the
    complete-profile (all-runs) companion panel is in the supplementary."""
    _render_main_metric(
        quantmsdiann_rows, svg_path, metric="proteins",
        ylabel="Protein groups per run (1% FDR)",
        label_fmt=lambda v: f"{v / 1000:.1f}k",
    )


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


def write_counts_tsv(
    long_df: pd.DataFrame,
    tsv_path: Path,
    *,
    quantmsdiann_unfiltered_rows: (
        list[tuple[str, str, int, int]] | None
    ) = None,
) -> None:
    """Write the per-dataset counts table for the benchmark figure.

    Columns: dataset, source, tool, version, precursors, proteins,
    filter_policy. `filter_policy` is `"target_only"` on the headline
    rows (rows consumed by the figures) and `"unfiltered"` on the
    companion rows derived from the raw pr_matrix.tsv / pg_matrix.tsv
    line counts. Both rows are written for every (dataset, version)
    pair so a reviewer can audit the conservative
    contaminant/entrapment/decoy filter delta (~1 % per
    ProteoBench DIA module).

    `quantmsdiann_unfiltered_rows` (optional) carries the
    (dataset, version, precursors, proteins) tuples computed BEFORE the
    filter — when present, those rows are appended with
    `filter_policy = "unfiltered"` and `source = "quantmsdiann"`.
    ProteoBench rows always carry `filter_policy = "target_only"`
    because community submissions are filtered upstream by ProteoBench's
    own `contaminant_flag = Cont_` parser.
    """
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "dataset", "source", "tool", "version",
        "precursors", "proteins", "filter_policy",
    ]
    # Stable ordering: quantmsdiann rows first (by dataset, version), then
    # ProteoBench rows (by dataset, descending precursors so the strongest
    # submissions are at the top of each block).
    qm = (long_df[long_df["source"] == "quantmsdiann"]
          .sort_values(["dataset", "version"]))
    pb = (long_df[long_df["source"] == "proteobench"]
          .sort_values(["dataset", "precursors"], ascending=[True, False]))
    out = pd.concat([qm, pb], ignore_index=True).copy()
    out["filter_policy"] = "target_only"
    if quantmsdiann_unfiltered_rows:
        unf_rows = []
        for dataset, version, precursors, proteins in quantmsdiann_unfiltered_rows:
            unf_rows.append({
                "dataset": dataset,
                "source": "quantmsdiann",
                "tool": "DIA-NN",
                "version": version,
                "precursors": precursors,
                "proteins": proteins,
                "library_kind": QUANTMSDIANN_LIBRARY_KIND,
                "filter_policy": "unfiltered",
            })
        out = pd.concat([out, pd.DataFrame(unf_rows)], ignore_index=True)
    # Pandas upgrades int+None columns to float; use pandas' nullable Int64
    # for proteins so ProteoBench rows print as the empty string and
    # quantmsdiann rows print as plain integers (no trailing .0).
    out = out[cols].copy()
    out["precursors"] = out["precursors"].astype("Int64")
    out["proteins"] = out["proteins"].astype("Int64")
    out.to_csv(tsv_path, sep="\t", index=False)


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


def write_median_table(
    df: pd.DataFrame,
    tsv_path: Path,
    *,
    quantmsdiann_rows_by_threshold_unfiltered: (
        dict[int, list[tuple[str, str, int, int]]] | None
    ) = None,
) -> None:
    """Auditable TSV of median precursor counts per DIA-NN version per module
    at the ≥1 and ≥3 replicate thresholds. One row per
    (dataset, threshold, source, version).

    Headline rows carry `filter_policy = "target_only"` (rows whose
    Protein.Group passes the conservative contaminant filter). When
    `quantmsdiann_rows_by_threshold_unfiltered` is provided, this writer
    appends companion `filter_policy = "unfiltered"` rows so a reviewer
    can compare. ProteoBench rows are always `target_only` because
    community submissions are filtered upstream by ProteoBench's own
    `contaminant_flag = Cont_` parser."""
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    if "filter_policy" not in out.columns:
        out["filter_policy"] = "target_only"
    if quantmsdiann_rows_by_threshold_unfiltered:
        extra_rows = []
        for thr, rows_unf in quantmsdiann_rows_by_threshold_unfiltered.items():
            for dataset, version, nr_prec, _ in rows_unf:
                extra_rows.append({
                    "dataset": dataset,
                    "min_replicates": thr,
                    "source": "quantmsdiann",
                    "version": _VERSION_LABELS.get(version, version),
                    "n_submissions": 1,
                    "median_nr_prec": int(nr_prec),
                    "filter_policy": "unfiltered",
                })
        if extra_rows:
            out = pd.concat([out, pd.DataFrame(extra_rows)], ignore_index=True)
    out = out.sort_values(
        ["dataset", "min_replicates", "source", "version", "filter_policy"]
    )
    out.to_csv(tsv_path, sep="\t", index=False)


REPORT_COUNTS_PATH = DATA_DIR / "report_counts.tsv"


def load_report_counts() -> dict[tuple[str, str], dict]:
    """Precursor + protein-group counts read from the DIA-NN *report*
    (`diann_report.parquet` / `.tsv`), NOT the `*_matrix.tsv` files.

    The matrices bake in `--matrix-spec-q` (0.05 run-specific) and, because the
    pipeline sets `--qvalue` to 0.01 for v1.8.1 but 0.05 for v2.5.1/enterprise,
    matrix row counts are filtered at *different* run-specific q-values per
    version and are therefore not comparable. These counts instead apply a
    uniform run AND global q-value <= 0.01 to every version, counting unique
    `Precursor.Id` (min1 / min3 = identified in >=1 / >=3 runs) and unique
    `Protein.Group` at global PG q <= 0.01. Target-only drops
    contaminant/entrapment/decoy groups; the unfiltered companion is kept for
    the audit TSV. Computed once on the cluster (where the multi-GB reports
    live) and staged at data/quantmsdiann_benchmarks/report_counts.tsv. Keyed
    by (dataset, version)."""
    df = pd.read_csv(REPORT_COUNTS_PATH, sep="\t",
                     dtype={"dataset": str, "version": str})
    out: dict[tuple[str, str], dict] = {}
    for _, r in df.iterrows():
        out[(r["dataset"], r["version"])] = {
            k: int(r[k]) for k in (
                "prec_min1_tgt", "prec_min1_unf", "prec_min3_tgt",
                "prec_min3_unf", "proteins_tgt", "proteins_unf",
                "prot_avg_tgt", "prot_complete_tgt",
            )
        }
    return out


def main() -> int:  # pragma: no cover
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    SUPP_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DATA_DIR.mkdir(parents=True, exist_ok=True)

    # quantmsdiann rows for the headline main panels. Counts come from the
    # DIA-NN *report* (not the *_matrix.tsv files) at a uniform run+global
    # q <= 0.01 across all versions (see load_report_counts), target-only.
    # Each row is (dataset, version, precursors_min1, proteins).
    quantmsdiann_rows: list[tuple[str, str, int, int]] = []
    # Parallel ≥3-replicate quantmsdiann counts (target-only).
    quantmsdiann_rows_min3: list[tuple[str, str, int, int]] = []
    # Companion unfiltered rows for the audit TSV (line counts before the
    # contaminant / entrapment / decoy filter). Per 2026-05-21 spec.
    quantmsdiann_rows_unfiltered: list[tuple[str, str, int, int]] = []
    quantmsdiann_rows_min3_unfiltered: list[tuple[str, str, int, int]] = []
    # Per-run-average and complete-profile protein rows (DIA-NN author's metric,
    # run-specific PG.Q <= 1%): these drive the version-comparison protein panels.
    # The union `proteins` above feeds only the auditable counts TSV.
    quantmsdiann_rows_protavg: list[tuple[str, str, int, int]] = []
    quantmsdiann_rows_complete: list[tuple[str, str, int, int]] = []
    report_counts = load_report_counts()
    for dataset in DATASET_TO_MODULE:
        for version in DIANN_VERSIONS:
            c = report_counts.get((dataset, version))
            if c is None:
                print(f"WARN: no report_counts for {dataset}/{version}",
                      file=sys.stderr)
                continue
            precursors, precursors_unf = c["prec_min1_tgt"], c["prec_min1_unf"]
            precursors_min3, precursors_min3_unf = (
                c["prec_min3_tgt"], c["prec_min3_unf"])
            proteins, proteins_unf = c["proteins_tgt"], c["proteins_unf"]
            quantmsdiann_rows.append((dataset, version, precursors, proteins))
            quantmsdiann_rows_min3.append(
                (dataset, version, precursors_min3, proteins)
            )
            quantmsdiann_rows_unfiltered.append(
                (dataset, version, precursors_unf, proteins_unf)
            )
            quantmsdiann_rows_min3_unfiltered.append(
                (dataset, version, precursors_min3_unf, proteins_unf)
            )
            quantmsdiann_rows_protavg.append(
                (dataset, version, precursors, c["prot_avg_tgt"])
            )
            quantmsdiann_rows_complete.append(
                (dataset, version, precursors, c["prot_complete_tgt"])
            )
            print(f"{dataset} {version}: precursors_min1={precursors:,} "
                  f"(unfiltered={precursors_unf:,})  "
                  f"precursors_min3={precursors_min3:,} "
                  f"(unfiltered={precursors_min3_unf:,})  "
                  f"proteins={proteins:,} (unfiltered={proteins_unf:,})  "
                  f"[report, run+global q<=0.01]")

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
        FIGURES_DIR / "main_benchmarks_precursors.svg",
    )

    print("Rendering protein-group main panel (per-run average)...")
    render_main_proteins(
        quantmsdiann_rows_protavg,
        FIGURES_DIR / "main_benchmarks_proteins.svg",
    )

    print("Rendering complete-profile protein panel (supplementary)...")
    _render_main_metric(
        quantmsdiann_rows_complete,
        SUPP_DIR / "supp_benchmarks_proteins_complete.svg",
        metric="proteins",
        ylabel="Complete-profile protein groups (in all runs)",
        label_fmt=lambda v: f"{v / 1000:.1f}k",
    )

    print("Writing auditable counts TSV (≥1)...")
    write_counts_tsv(
        long_df, FIG_DATA_DIR / "counts.tsv",
        quantmsdiann_unfiltered_rows=quantmsdiann_rows_unfiltered,
    )
    print("Writing auditable counts TSV (≥3)...")
    write_counts_tsv(
        long_df_min3, FIG_DATA_DIR / "counts_min3.tsv",
        quantmsdiann_unfiltered_rows=quantmsdiann_rows_min3_unfiltered,
    )

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
        median_df, FIG_DATA_DIR / "median_nr_prec_by_version.tsv",
        quantmsdiann_rows_by_threshold_unfiltered={
            1: quantmsdiann_rows_unfiltered,
            3: quantmsdiann_rows_min3_unfiltered,
        },
    )

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
