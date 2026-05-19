"""Combined cell-lines atlas figure for the quantmsdiann manuscript.

Integrates the three independent cell-line reanalyses already shipped
(PXD003539 NCI-60, PXD030304 ProCan-DepMapSanger 949 lines, PXD004701
BC 76 lines) into one paper-ready multi-panel figure that positions the
quantmsdiann pipeline as a single uniform tool covering broad cancer
cell-line / tissue / proteome space.

Panels:
- A (top-left): dataset-level reproducibility (paper vs quantmsdiann
  headline protein counts per dataset).
- B (top-right): 3-set Venn of normalised cell-line sets.
- C (bottom-left, wider): unified-axis pan-cancer tissue coverage,
  stacked horizontal bars per dataset.
- D (bottom-right): 3-set Venn of UniProt accessions detected by
  quantmsdiann.

Reuses parsers/helpers from the per-dataset scripts via imports — no
duplicated logic, no new downloads. Reads pre-cached JSON caches for
PXD030304 / PXD004701 protein sets (don't re-stream 33 GB parquets) and
the PXD003539 pr_matrix.tsv (already on disk) for its protein set.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis.figure_original_vs_quantmsdiann import (
    PR_METADATA_COLS,
    normalize_cell_line,
)
from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
    parse_procan_mapping,
)
from analysis.figure_pxd004701_sun_vs_quantmsdiann import BC_SUBTYPES  # noqa: F401
from analysis.venn_protein_accessions import extract_accessions_diann


REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = REPO_ROOT / "data"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "combined"

# Per-dataset input paths (all local, no downloads required).
PXD003539_SDRF = DATA_ROOT / "PXD003539" / "PXD003539.sdrf.tsv"
PXD003539_PR_MATRIX = DATA_ROOT / "PXD003539" / "diann_report.pr_matrix.tsv"
PXD003539_COUNTS_TSV = REPO_ROOT / "analysis" / "figures" / "PXD003539" / "counts.tsv"

PXD030304_SDRF = DATA_ROOT / "PXD030304" / "PXD030304.sdrf.tsv"
PXD030304_TISSUE_MAPPING = DATA_ROOT / "PXD030304" / "mapping_file_averaged.txt"
PXD030304_PROTEIN_JSON = DATA_ROOT / "PXD030304" / "diann_per_tissue_procan_filter.json"
PXD030304_COUNTS_TSV = REPO_ROOT / "analysis" / "figures" / "PXD030304" / "counts.tsv"

PXD004701_SDRF = DATA_ROOT / "PXD004701" / "PXD004701.sdrf.tsv"
PXD004701_PROTEIN_JSON = DATA_ROOT / "PXD004701" / "diann_per_subtype_consistency_filter.json"
PXD004701_COUNTS_TSV = REPO_ROOT / "analysis" / "figures" / "PXD004701" / "counts.tsv"

DATASET_COLORS = {
    "PXD003539": "#1b9e77",  # green
    "PXD030304": "#7570b3",  # purple
    "PXD004701": "#d95f02",  # orange
}

DATASET_LABELS = {
    "PXD003539": "PXD003539\n(Guo 2019)",
    "PXD030304": "PXD030304\n(ProCan 2022)",
    "PXD004701": "PXD004701\n(Sun 2023)",
}


# ---------------------------------------------------------------------------
# Headline numbers (Panel A): paper headline + quantmsdiann headline per
# dataset, locked at the same filter family used in each per-dataset spec.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DatasetHeadline:
    paper_count: int      # original-paper apples-to-apples number
    diann_count: int      # quantmsdiann under the same filter
    paper_label: str
    metric: str           # short metric description for tooltip / TSV


DATASET_HEADLINES: dict[str, DatasetHeadline] = {
    # Guo 2019 OpenSWATH deposited matrix (1% FDR) vs quantmsdiann DIA-NN
    # diannsummary.log protein groups at 1% global q-value.
    "PXD003539": DatasetHeadline(
        paper_count=6_556,
        diann_count=6_927,
        paper_label="Guo 2019 (OpenSWATH)",
        metric="Protein groups (1% global FDR)",
    ),
    # ProCan 2022 paper headline (Global.Q.Value <= 0.01 proteotypic) vs
    # quantmsdiann diannsummary.log protein groups at 1% global q-value.
    # The +22% post-ProCan-filter union number lives in the spec/counts.tsv
    # discussion text; the figure bar is the headline.
    "PXD030304": DatasetHeadline(
        paper_count=8_498,
        diann_count=9_370,
        paper_label="ProCan 2022",
        metric="Proteins @ Global.Q.Value <= 0.01",
    ),
    # Sun 2023 paper consistency-filtered (>=10% detection, proteotypic, 1%
    # global Q.Value) vs quantmsdiann under the same two-stage filter.
    "PXD004701": DatasetHeadline(
        paper_count=6_091,
        diann_count=6_296,
        paper_label="Sun 2023",
        metric="Proteins (consistency filter, <=90% missing)",
    ),
}


# ---------------------------------------------------------------------------
# Cell-line set parsing (Panel B and Panel C inputs)
# ---------------------------------------------------------------------------


def cell_lines_from_sdrf(
    sdrf_path: Path,
    cell_line_col: str = "characteristics[cell line]",
) -> set[str]:
    """Return the set of normalised cell-line names from any of the three
    SDRFs. Normalisation via `normalize_cell_line` strips the `NCI-` prefix
    and all non-alphanumeric characters, uppercases — so 'CCRF-CEM',
    'NCI-H226', 'Hs-578-T' collide across datasets with their alternative
    spellings."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    if cell_line_col not in df.columns:
        raise ValueError(
            f"SDRF {sdrf_path} missing required column: {cell_line_col!r}"
        )
    out: set[str] = set()
    for raw in df[cell_line_col].fillna(""):
        norm = normalize_cell_line(raw)
        if norm:
            out.add(norm)
    return out


# ---------------------------------------------------------------------------
# Cell-line -> tissue (per dataset, on the unified ProCan 28-tissue axis)
# ---------------------------------------------------------------------------

# NCI-60 disease / organism-part keyword -> ProCan tissue. Hardcoded mapping
# rules. Keys are lowercase substrings searched (in order) against the
# disease + organism-part text; the first match wins.
PXD003539_TISSUE_RULES: list[tuple[str, str]] = [
    # haematopoietic / lymphoid
    ("leukemia", "Haematopoietic and Lymphoid"),
    ("lymphoma", "Haematopoietic and Lymphoid"),
    ("myeloma", "Haematopoietic and Lymphoid"),
    ("blood", "Haematopoietic and Lymphoid"),
    ("bone marrow", "Haematopoietic and Lymphoid"),
    # central nervous system
    ("glioblastoma", "Central Nervous System"),
    ("gliosarcoma", "Central Nervous System"),
    ("astrocytoma", "Central Nervous System"),
    ("central nervous system", "Central Nervous System"),
    ("brain", "Central Nervous System"),
    # breast
    ("breast", "Breast"),
    # large intestine
    ("colon", "Large Intestine"),
    ("colorectal", "Large Intestine"),
    ("large intestine", "Large Intestine"),
    # lung (incl. pleural mesothelioma which arises in lung-adjacent pleura
    # and is grouped with lung in the NCI-60 categorisation)
    ("lung", "Lung"),
    ("non-small cell", "Lung"),
    ("mesothelioma", "Lung"),
    # skin / melanoma
    ("melanoma", "Skin"),
    ("skin", "Skin"),
    # ovary
    ("ovarian", "Ovary"),
    ("ovary", "Ovary"),
    # prostate
    ("prostate", "Prostate"),
    # kidney / renal
    ("renal", "Kidney"),
    ("kidney", "Kidney"),
]


def harmonise_pxd003539_tissue(
    disease_text: str | None,
    organism_part: str | None = None,
) -> str | None:
    """Map a PXD003539 (disease, organism_part) pair to one of the 28 ProCan
    tissue categories. Returns None if no rule matches (uncategorised lines
    are dropped from Panel C, not silently bucketed into Other tissue).

    The rules cover the 9 NCI-60 cancer types; the matched ProCan-axis
    tissue is the most specific bucket NCI-60 maps onto. Notes:
    - "leukemia / lymphoma / myeloma" maps to ProCan's
      `Haematopoietic and Lymphoid` (NCI-60 has 6 leukemia lines under
      Blood, all of which are present as Haematopoietic & Lymphoid in
      ProCan's mapping_file_averaged.txt).
    - "Pleural epithelioid mesothelioma" (NCI-60 MESO line under Lung)
      maps to ProCan `Lung` — the only mesothelioma category ProCan
      retains.
    - Brain / CNS / glioblastoma / gliosarcoma / astrocytoma all map to
      ProCan `Central Nervous System`."""
    disease = (disease_text or "").strip().lower()
    organism = (organism_part or "").strip().lower()
    haystack = f"{disease} {organism}"
    for needle, tissue in PXD003539_TISSUE_RULES:
        if needle in haystack:
            return tissue
    return None


def cell_line_tissue_pxd003539(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, ProCan_tissue]` for PXD003539.

    Each row of the SDRF carries (cell_line, disease, organism_part);
    rows of the same cell line have identical disease/organism so we
    collapse to a deduplicated mapping. Cell lines whose disease/organism
    text doesn't match any rule are dropped from the result."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    needed = [
        "characteristics[cell line]",
        "characteristics[disease]",
        "characteristics[organism part]",
    ]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD003539 SDRF {sdrf_path} missing columns: {missing}"
        )
    out: dict[str, str] = {}
    for cell, disease, organism in zip(
        df["characteristics[cell line]"],
        df["characteristics[disease]"],
        df["characteristics[organism part]"],
    ):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        tissue = harmonise_pxd003539_tissue(disease, organism)
        if tissue is None:
            continue
        out[norm] = tissue
    return out


def cell_line_tissue_pxd030304(
    sdrf_path: Path,
    mapping_path: Path,
) -> dict[str, str]:
    """Return `dict[normalised_cell_line, ProCan_tissue]` for PXD030304.

    The ProCan figshare `mapping_file_averaged.txt` is the canonical source
    of truth for the 28-tissue axis; we use it directly. The SDRF gives us
    the cell lines actually present in the dataset (the figshare mapping
    covers all 949 cell lines anyway, but going through the SDRF keeps the
    contract identical across the three datasets)."""
    cl_to_tissue = parse_procan_mapping(mapping_path)
    # Build a normalised-name lookup so SDRF spellings (which already match
    # figshare exactly per spec) collide regardless of case/punctuation.
    norm_lookup: dict[str, str] = {}
    for cell, tissue in cl_to_tissue.items():
        n = normalize_cell_line(cell)
        if n:
            norm_lookup[n] = tissue
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    if "characteristics[cell line]" not in df.columns:
        raise ValueError(
            f"PXD030304 SDRF {sdrf_path} missing characteristics[cell line]"
        )
    out: dict[str, str] = {}
    for cell in df["characteristics[cell line]"].fillna(""):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        tissue = norm_lookup.get(norm)
        if tissue is None:
            continue
        out[norm] = tissue
    return out


def cell_line_tissue_pxd004701(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, 'Breast']` for PXD004701.

    All 76 cell lines in PXD004701 are breast-cancer-derived; the BC
    subtype split is internal to that dataset and irrelevant for the
    unified-tissue axis."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    if "characteristics[cell line]" not in df.columns:
        raise ValueError(
            f"PXD004701 SDRF {sdrf_path} missing characteristics[cell line]"
        )
    out: dict[str, str] = {}
    for cell in df["characteristics[cell line]"].fillna(""):
        norm = normalize_cell_line(cell)
        if norm:
            out[norm] = "Breast"
    return out


def combined_tissue_table(
    per_dataset: dict[str, dict[str, str]],
) -> list[tuple[str, dict[str, int]]]:
    """Combine per-dataset cell-line-to-tissue mappings into a single
    `[(tissue, {dataset: cell_line_count})]` list, sorted by total
    cell-line count descending. Tissues with zero contribution from every
    dataset are dropped (always — the input only carries mapped lines)."""
    counts: dict[str, dict[str, int]] = {}
    for dataset, cl_to_tissue in per_dataset.items():
        tissue_counts: dict[str, int] = {}
        for tissue in cl_to_tissue.values():
            tissue_counts[tissue] = tissue_counts.get(tissue, 0) + 1
        for tissue, n in tissue_counts.items():
            counts.setdefault(tissue, {ds: 0 for ds in per_dataset})[dataset] = n
    out = sorted(
        counts.items(),
        key=lambda kv: -sum(kv[1].values()),
    )
    return out


# ---------------------------------------------------------------------------
# Protein-accession sets (Panel D inputs)
# ---------------------------------------------------------------------------


def pxd003539_protein_accessions(pr_matrix_path: Path) -> set[str]:
    """Return the set of UniProt accessions detected by quantmsdiann in
    PXD003539. Reads `diann_report.pr_matrix.tsv` (already cached locally
    at ~67 MB) and collects accessions from `Protein.Group` strings on rows
    with at least one non-NA quantification across the 120 runs.

    Uses `extract_accessions_diann` so the accession-normalisation semantics
    match the PXD030304 / PXD004701 caches (semicolon split, isoform suffix
    stripped, CONTAM_/ENTRAP_/DECOY_ prefix removed)."""
    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD003539 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    if not sample_cols:
        raise ValueError("PXD003539 pr_matrix has no per-run sample columns")
    quantified = df[df[sample_cols].notna().any(axis=1)]
    accessions: set[str] = set()
    for pg in quantified["Protein.Group"]:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions


def _accessions_from_json_cache(json_path: Path) -> set[str]:
    """Helper for the PXD030304 / PXD004701 cached JSONs. Both caches are
    `{group_key: [Protein.Group, ...]}` (group_key is tissue / subtype).
    We union the Protein.Group values across all groups and extract
    accessions."""
    with open(json_path, encoding="utf-8") as fh:
        payload = json.load(fh)
    pg_set: set[str] = set()
    for vs in payload.values():
        pg_set.update(vs)
    accessions: set[str] = set()
    for pg in pg_set:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions


def pxd030304_protein_accessions(json_path: Path) -> set[str]:
    """Read `diann_per_tissue_procan_filter.json` and return the union of
    UniProt accessions across the 28 tissues. The JSON is the cached output
    of `proteins_per_tissue_quantmsdiann_procan_filter` and contains
    Protein.Group strings as deposited in the parquet."""
    return _accessions_from_json_cache(json_path)


def pxd004701_protein_accessions(json_path: Path) -> set[str]:
    """Read `diann_per_subtype_consistency_filter.json` and return the union
    of UniProt accessions across the 3 BC subtypes. The JSON is the cached
    output of `proteins_per_subtype_quantmsdiann_consistency_filter`."""
    return _accessions_from_json_cache(json_path)


# ---------------------------------------------------------------------------
# Panel renderers
# ---------------------------------------------------------------------------


def _annotate_panel_letter(ax, letter: str) -> None:
    ax.text(
        -0.07, 1.02, letter,
        transform=ax.transAxes, fontsize=14, fontweight="bold",
        ha="left", va="bottom",
    )


def _render_panel_a_headlines(ax, headlines: dict[str, DatasetHeadline]) -> None:
    """Grouped bar chart, 3 dataset groups × 2 bars (paper vs quantmsdiann).
    No legend in-panel; bars are colour-coded paper (grey) vs quantmsdiann
    (blue) following the per-dataset figure convention. Paper label sits
    under each group."""
    datasets = list(headlines.keys())
    bar_width = 0.36
    x = list(range(len(datasets)))
    paper_vals = [headlines[d].paper_count for d in datasets]
    diann_vals = [headlines[d].diann_count for d in datasets]
    bars_p = ax.bar(
        [xi - bar_width / 2 for xi in x], paper_vals,
        width=bar_width, color="#9e9e9e", label="Original analysis",
    )
    bars_d = ax.bar(
        [xi + bar_width / 2 for xi in x], diann_vals,
        width=bar_width, color="#1f77b4", label="quantmsdiann (DIA-NN)",
    )
    for bars, vals in ((bars_p, paper_vals), (bars_d, diann_vals)):
        for bar, v in zip(bars, vals):
            ax.text(
                bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{v:,}", ha="center", va="bottom", fontsize=8,
            )
    ax.set_xticks(x)
    ax.set_xticklabels(
        [
            f"{headlines[d].paper_label}\n{d}" for d in datasets
        ],
        fontsize=8,
    )
    ax.set_ylabel("Protein groups")
    ymax = max(max(paper_vals), max(diann_vals)) * 1.18
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper left", frameon=False, fontsize=8)


def _render_panel_b_cellline_venn(
    ax,
    sets: dict[str, set[str]],
    region_sizes: dict[str, int] | None = None,
) -> None:
    """3-set Venn of normalised cell-line names. Uses matplotlib_venn.venn3
    (already a project dependency). Absolute region counts displayed (no
    percentages — see spec)."""
    try:
        from matplotlib_venn import venn3
    except ImportError as exc:  # pragma: no cover - dependency is required
        raise RuntimeError(
            "matplotlib_venn required for Panel B (already in requirements.txt)"
        ) from exc
    ds_order = ["PXD003539", "PXD030304", "PXD004701"]
    s1, s2, s3 = (sets[d] for d in ds_order)
    v = venn3(
        subsets=(s1, s2, s3),
        set_labels=(
            f"PXD003539\n(n={len(s1)})",
            f"PXD030304\n(n={len(s2)})",
            f"PXD004701\n(n={len(s3)})",
        ),
        set_colors=(DATASET_COLORS[d] for d in ds_order),
        alpha=0.55,
        ax=ax,
    )
    for sl in v.set_labels:
        if sl is not None:
            sl.set_fontsize(9)
    for region_id in ("100", "010", "001", "110", "101", "011", "111"):
        lbl = v.get_label_by_id(region_id)
        if lbl is not None:
            lbl.set_fontsize(8)


def _render_panel_c_tissue_coverage(
    ax,
    rows: list[tuple[str, dict[str, int]]],
) -> None:
    """Stacked horizontal bars on a unified tissue axis. Each bar segment
    is a dataset's per-tissue cell-line count; tissues sorted by descending
    total across datasets. Tissues with zero contribution from all three
    datasets never appear (filtered upstream)."""
    if not rows:
        ax.text(0.5, 0.5, "no tissues", ha="center", va="center",
                transform=ax.transAxes)
        return
    tissues = [r[0] for r in rows]
    ds_order = ["PXD003539", "PXD030304", "PXD004701"]
    n = len(tissues)
    y = list(range(n))
    # Plot top tissues at top of chart (largest first).
    y_top = list(reversed(y))
    left = [0] * n
    for ds in ds_order:
        vals = [r[1].get(ds, 0) for r in rows]
        ax.barh(
            y_top, vals, left=left,
            color=DATASET_COLORS[ds],
            label=DATASET_LABELS[ds].replace("\n", " "),
            edgecolor="white", linewidth=0.4,
        )
        left = [l + v for l, v in zip(left, vals)]
    # Annotate total cell-line count at the right edge of each bar.
    totals = [sum(r[1].values()) for r in rows]
    for yi, t in zip(y_top, totals):
        ax.text(t + max(totals) * 0.005, yi, f"{t:,}",
                ha="left", va="center", fontsize=7)
    ax.set_yticks(y_top)
    ax.set_yticklabels(tissues, fontsize=7)
    ax.set_xlabel("Cell lines (sum across datasets)", fontsize=9)
    ax.set_xlim(0, max(totals) * 1.10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=7)


def _render_panel_d_protein_venn(
    ax,
    sets: dict[str, set[str]],
) -> None:
    """3-set Venn of UniProt accessions across the three quantmsdiann
    analyses. Same colour palette as Panel B."""
    try:
        from matplotlib_venn import venn3
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "matplotlib_venn required for Panel D (already in requirements.txt)"
        ) from exc
    ds_order = ["PXD003539", "PXD030304", "PXD004701"]
    s1, s2, s3 = (sets[d] for d in ds_order)
    v = venn3(
        subsets=(s1, s2, s3),
        set_labels=(
            f"PXD003539\n(n={len(s1):,})",
            f"PXD030304\n(n={len(s2):,})",
            f"PXD004701\n(n={len(s3):,})",
        ),
        set_colors=(DATASET_COLORS[d] for d in ds_order),
        alpha=0.55,
        ax=ax,
    )
    for sl in v.set_labels:
        if sl is not None:
            sl.set_fontsize(9)
    for region_id in ("100", "010", "001", "110", "101", "011", "111"):
        lbl = v.get_label_by_id(region_id)
        if lbl is not None:
            lbl.set_fontsize(8)


def render_atlas(
    headlines: dict[str, DatasetHeadline],
    cellline_sets: dict[str, set[str]],
    tissue_rows: list[tuple[str, dict[str, int]]],
    accession_sets: dict[str, set[str]],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """Compose the 4-panel atlas figure with gridspec sized so the
    per-tissue stacked bar (Panel C) gets a wider/taller slot than the
    two Venns. Paper-ready: no figure title, no in-figure footer; panel
    letters A/B/C/D in the top-left corner of each subplot."""
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.05, 1],
        height_ratios=[1, 1.45],
        hspace=0.32, wspace=0.20,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    _render_panel_a_headlines(ax_a, headlines)
    _annotate_panel_letter(ax_a, "A")

    _render_panel_b_cellline_venn(ax_b, cellline_sets)
    _annotate_panel_letter(ax_b, "B")

    _render_panel_c_tissue_coverage(ax_c, tissue_rows)
    _annotate_panel_letter(ax_c, "C")

    _render_panel_d_protein_venn(ax_d, accession_sets)
    _annotate_panel_letter(ax_d, "D")

    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    if svg_path is not None:
        fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# counts.tsv writer
# ---------------------------------------------------------------------------


def _venn_region_sizes_3(
    sets: dict[str, set[str]],
    ds_order: list[str],
) -> dict[str, int]:
    """Return the 7 Venn region sizes for a 3-set Venn:
    `{ds}_only`, `{ds_a}+{ds_b}` for every pair, `all_three`."""
    a, b, c = (sets[d] for d in ds_order)
    return {
        f"{ds_order[0]}_only": len(a - b - c),
        f"{ds_order[1]}_only": len(b - a - c),
        f"{ds_order[2]}_only": len(c - a - b),
        f"{ds_order[0]}+{ds_order[1]}": len(a & b - c),
        f"{ds_order[0]}+{ds_order[2]}": len(a & c - b),
        f"{ds_order[1]}+{ds_order[2]}": len(b & c - a),
        "all_three": len(a & b & c),
    }


def write_combined_counts_tsv(
    tsv_path: Path,
    headlines: dict[str, DatasetHeadline],
    cellline_sets: dict[str, set[str]],
    tissue_rows: list[tuple[str, dict[str, int]]],
    accession_sets: dict[str, set[str]],
) -> None:
    """Auditable TSV with Panel-feeding numbers (Panel A bars, Panel B Venn
    regions, Panel C per-(tissue, dataset) counts, Panel D Venn regions)."""
    ds_order = ["PXD003539", "PXD030304", "PXD004701"]
    rows: list[tuple[str, str, int, str]] = []
    # Panel A: 6 bars (paper + diann per dataset).
    for ds in ds_order:
        h = headlines[ds]
        rows.append((
            f"Panel A | {ds} | original",
            h.paper_label,
            h.paper_count,
            h.metric,
        ))
        rows.append((
            f"Panel A | {ds} | quantmsdiann",
            "quantmsdiann (DIA-NN)",
            h.diann_count,
            h.metric,
        ))
    # Panel B: 7 Venn regions over cell-line sets.
    cl_regions = _venn_region_sizes_3(cellline_sets, ds_order)
    for region, n in cl_regions.items():
        rows.append((
            f"Panel B | cell-line Venn | {region}",
            "normalised cell-line names from SDRF (normalize_cell_line)",
            n,
            "SDRF cell-line set membership at normalised-name level",
        ))
    # Panel C: per (tissue, dataset) cell-line counts.
    for tissue, by_ds in tissue_rows:
        for ds in ds_order:
            rows.append((
                f"Panel C | tissue | {tissue}",
                ds,
                by_ds.get(ds, 0),
                "cell lines mapped to this tissue (unified ProCan 28-tissue axis)",
            ))
    # Panel D: 7 Venn regions over UniProt accessions.
    acc_regions = _venn_region_sizes_3(accession_sets, ds_order)
    for region, n in acc_regions.items():
        rows.append((
            f"Panel D | accession Venn | {region}",
            "UniProt accessions extracted from Protein.Group (extract_accessions_diann)",
            n,
            "PXD003539 from pr_matrix; PXD030304/PXD004701 from cached "
            "per-tissue/per-subtype JSON",
        ))
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", encoding="utf-8") as fh:
        fh.write("metric\tsource\tcount\tnote\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


# ---------------------------------------------------------------------------
# Prerequisite checks (idempotency / fail-loudly)
# ---------------------------------------------------------------------------


PREREQS: list[tuple[Path, str]] = [
    (PXD003539_SDRF, "python -m analysis.figure_original_vs_quantmsdiann"),
    (PXD003539_PR_MATRIX, "python -m analysis.figure_original_vs_quantmsdiann"),
    (PXD030304_SDRF, "python -m analysis.figure_pxd030304_procan_vs_quantmsdiann"),
    (PXD030304_TISSUE_MAPPING, "python -m analysis.figure_pxd030304_procan_vs_quantmsdiann"),
    (PXD030304_PROTEIN_JSON, "python -m analysis.figure_pxd030304_procan_vs_quantmsdiann"),
    (PXD004701_SDRF, "python -m analysis.figure_pxd004701_sun_vs_quantmsdiann"),
    (PXD004701_PROTEIN_JSON, "python -m analysis.figure_pxd004701_sun_vs_quantmsdiann"),
]


def check_prerequisites() -> list[tuple[Path, str]]:
    """Return the list of `(missing_path, instruction)` tuples; empty list
    means all prerequisites are present."""
    return [
        (p, cmd) for p, cmd in PREREQS
        if not (p.exists() and p.stat().st_size > 0)
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:  # pragma: no cover
    missing = check_prerequisites()
    if missing:
        print("Combined atlas requires the per-dataset scripts to run first.",
              file=sys.stderr)
        for path, cmd in missing:
            rel = path.relative_to(REPO_ROOT)
            print(f"  Missing: {rel}", file=sys.stderr)
            print(f"  Produce it via: {cmd}", file=sys.stderr)
        return 1

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading cell-line sets from SDRFs...")
    cellline_sets: dict[str, set[str]] = {
        "PXD003539": cell_lines_from_sdrf(PXD003539_SDRF),
        "PXD030304": cell_lines_from_sdrf(PXD030304_SDRF),
        "PXD004701": cell_lines_from_sdrf(PXD004701_SDRF),
    }
    for d, s in cellline_sets.items():
        print(f"  {d}: {len(s):,} cell lines")

    print("Building per-dataset cell-line -> tissue maps...")
    cl_tissue_3539 = cell_line_tissue_pxd003539(PXD003539_SDRF)
    cl_tissue_30304 = cell_line_tissue_pxd030304(PXD030304_SDRF,
                                                 PXD030304_TISSUE_MAPPING)
    cl_tissue_4701 = cell_line_tissue_pxd004701(PXD004701_SDRF)
    tissue_rows = combined_tissue_table({
        "PXD003539": cl_tissue_3539,
        "PXD030304": cl_tissue_30304,
        "PXD004701": cl_tissue_4701,
    })
    print(f"  unified-axis tissues: {len(tissue_rows)}")
    print("  Top-5 tissues by combined cell-line count:")
    for tissue, by_ds in tissue_rows[:5]:
        total = sum(by_ds.values())
        seg = ", ".join(
            f"{ds}:{by_ds.get(ds, 0)}" for ds in ("PXD003539", "PXD030304", "PXD004701")
        )
        print(f"    {tissue:<32s} total={total:>4d}  ({seg})")

    print("Loading protein-accession sets...")
    acc_3539 = pxd003539_protein_accessions(PXD003539_PR_MATRIX)
    acc_30304 = pxd030304_protein_accessions(PXD030304_PROTEIN_JSON)
    acc_4701 = pxd004701_protein_accessions(PXD004701_PROTEIN_JSON)
    accession_sets = {
        "PXD003539": acc_3539,
        "PXD030304": acc_30304,
        "PXD004701": acc_4701,
    }
    for d, s in accession_sets.items():
        print(f"  {d}: {len(s):,} accessions")

    print("Rendering atlas figure...")
    pdf = FIGURES_DIR / "atlas.pdf"
    png = FIGURES_DIR / "atlas.png"
    svg = FIGURES_DIR / "atlas.svg"
    render_atlas(
        DATASET_HEADLINES, cellline_sets, tissue_rows, accession_sets,
        pdf, png, svg,
    )
    print(f"  saved: {pdf}, {png}, {svg}")

    print("Writing combined counts.tsv...")
    tsv = FIGURES_DIR / "combined_counts.tsv"
    write_combined_counts_tsv(
        tsv, DATASET_HEADLINES, cellline_sets, tissue_rows, accession_sets,
    )
    print(f"  saved: {tsv}")

    ds_order = ["PXD003539", "PXD030304", "PXD004701"]
    print("\nPanel B (cell-line Venn) region sizes:")
    for region, n in _venn_region_sizes_3(cellline_sets, ds_order).items():
        print(f"  {region:<36s} {n:>6,}")
    print("\nPanel D (accession Venn) region sizes:")
    for region, n in _venn_region_sizes_3(accession_sets, ds_order).items():
        print(f"  {region:<36s} {n:>6,}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
