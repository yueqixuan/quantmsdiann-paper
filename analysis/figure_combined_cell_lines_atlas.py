"""Combined cell-lines atlas figure for the quantmsdiann manuscript.

Integrates the five independent cell-line reanalyses already shipped
(PXD003539 NCI-60, PXD030304 ProCan-DepMapSanger 949 lines, PXD004701
BC 76 lines, PXD017199 Tognetti 67 breast lines, PXD041421 Wang 2023)
into two paper-ready figures that position the quantmsdiann pipeline as
a single uniform tool covering broad cancer cell-line / tissue / proteome
space. The original single A–H grid was split in two for readability.

`atlas_overlap.svg` (protein-accession overlap):
- A (top, wide): UpSet plot of UniProt accessions detected by
  quantmsdiann across the five cohorts, with per-cohort headline counts;
  the 5-set UpSet replaces the unreadable 4-set Venn.
- B / C (bottom row): dataset-level reproducibility (paper vs
  quantmsdiann headline counts; PXD017199 / PXD041421 have no paper
  headline so only their quantmsdiann bar is drawn) and the
  detection-count histogram (how many cohorts each protein is seen in,
  giving the pan-cohort core).

`atlas_distribution.svg` (per-tissue coverage):
- A (top, wide): unified-axis pan-cancer tissue coverage, stacked
  horizontal bars per cohort (cell-line counts).
- B (bottom, wide): per-tissue unique-protein stacked bars (parallel to
  A; protein counts instead of cell-line counts).
The breadth-vs-depth scatter (former Panel C/E) was dropped 2026-05-29
and the Expression Atlas overlap (former Panel H) was removed earlier —
it duplicates analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg.

Reuses parsers/helpers from the per-dataset scripts via imports — no
duplicated logic, no new downloads. Reads pre-cached JSON caches for
PXD030304 / PXD004701 protein sets (don't re-stream 33 GB parquets) and
the PXD003539 / PXD017199 pr_matrix.tsv (already on disk) for their
protein sets.
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

from analysis.contaminant_filter import (
    count_target_protein_groups,
    is_target_protein_group,
)
from analysis.figure_original_vs_quantmsdiann import (
    PR_METADATA_COLS,
    count_target_protein_groups_pr_matrix,
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
PXD003539_COUNTS_TSV = REPO_ROOT / "analysis" / "figures" / "PXD003539" / "data" / "counts.tsv"

PXD030304_SDRF = DATA_ROOT / "PXD030304" / "PXD030304.sdrf.tsv"
PXD030304_TISSUE_MAPPING = DATA_ROOT / "PXD030304" / "mapping_file_averaged.txt"
PXD030304_PROTEIN_JSON = DATA_ROOT / "PXD030304" / "diann_per_tissue_procan_filter.json"
PXD030304_PG_MATRIX = DATA_ROOT / "PXD030304" / "diann_report.pg_matrix.tsv"
PXD030304_COUNTS_TSV = REPO_ROOT / "analysis" / "figures" / "PXD030304" / "data" / "counts.tsv"

PXD004701_SDRF = DATA_ROOT / "PXD004701" / "PXD004701.sdrf.tsv"
PXD004701_PROTEIN_JSON = DATA_ROOT / "PXD004701" / "diann_per_subtype_consistency_filter.json"
PXD004701_PG_MATRIX = DATA_ROOT / "PXD004701" / "diann_report.pg_matrix.tsv"
PXD004701_COUNTS_TSV = REPO_ROOT / "analysis" / "figures" / "PXD004701" / "data" / "counts.tsv"

# PXD017199 (Tognetti 2021) — 67 breast (incl. 5 normal-like) lines.
# No JSON cache: per-cell-line accessions are derived live from the pr_matrix.
PXD017199_SDRF = DATA_ROOT / "PXD017199" / "PXD017199.sdrf.tsv"
PXD017199_PR_MATRIX = DATA_ROOT / "PXD017199" / "diann_report.pr_matrix.tsv"
PXD017199_PG_MATRIX = DATA_ROOT / "PXD017199" / "diann_report.pg_matrix.tsv"

# PXD041421 (Wang 2023) — TIQUEST diaPASEF batch-effect testbed, 48 runs
# across A549 (24 reps) + K562 (24 reps). Per 2026-05-21 spec §2 atlas-only
# integration; no per-cohort figure.
PXD041421_SDRF = DATA_ROOT / "PXD041421" / "PXD041421.sdrf.tsv"
PXD041421_PR_MATRIX = DATA_ROOT / "PXD041421" / "diann_report.pr_matrix.tsv"
PXD041421_PG_MATRIX = DATA_ROOT / "PXD041421" / "diann_report.pg_matrix.tsv"

# Expression Atlas E-PROT-73 NCI-60 catalogue (Panel H).
E_PROT_73_TSV = DATA_ROOT / "E-PROT-73-query-results.tsv"

DATASET_COLORS = {
    "PXD003539": "#1b9e77",  # green
    "PXD030304": "#7570b3",  # purple
    "PXD004701": "#d95f02",  # orange
    "PXD017199": "#e6ab02",  # mustard / yellow-amber
    "PXD041421": "#8da0cb",  # lavender (5th cohort, 2026-05-21 spec)
}

DATASET_LABELS = {
    "PXD003539": "PXD003539\n(Guo 2019)",
    "PXD030304": "PXD030304\n(ProCan 2022)",
    "PXD004701": "PXD004701\n(Sun 2023)",
    "PXD017199": "PXD017199\n(Tognetti 2021)",
    "PXD041421": "PXD041421\n(Wang 2023)",
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


# `diann_count` values are lazily computed via
# `count_target_protein_groups(pg_matrix_path)` per the 2026-05-21
# contaminant-filter spec — the headline is the post-filter target-only
# protein-group count on the same `diann_report.pg_matrix.tsv` that the
# per-cohort scripts use. Unfiltered companion numbers go into the audit
# TSV (`write_combined_counts_tsv`).
# PXD003539 has no `pg_matrix.tsv` on disk; we derive the target-only
# count from its `pr_matrix.tsv` unique Protein.Group rows instead.
#
# The lazy computation matters because two of these matrices are large
# (PXD030304 pg_matrix ~200 MB, PXD003539 pr_matrix ~67 MB). The
# `DATASET_HEADLINES` constant therefore initialises with sentinel
# values (`diann_count=0`) and the atlas `main()` calls
# `refresh_dataset_headlines()` to populate them from disk before
# rendering. Tests that don't need the on-disk numbers consume
# `DATASET_HEADLINES` directly with the placeholder values, keeping
# import time fast.


PXD003539_PAPER_COUNT = 6_556
PXD030304_PAPER_COUNT = 8_498
PXD004701_PAPER_COUNT = 6_091


# Default placeholder values used at import time; refresh from disk via
# `refresh_dataset_headlines()` (called in `main()`).
_DEFAULT_DIANN_COUNTS = {
    "PXD003539": 6_927,    # diannsummary.log baseline (unfiltered)
    "PXD030304": 9_370,    # diannsummary.log baseline (unfiltered)
    "PXD004701": 6_296,    # paper-side fallback
    "PXD017199": 10_572,   # diannsummary.log baseline (unfiltered)
    "PXD041421": 9_124,    # diannsummary.log baseline (unfiltered)
}


DATASET_HEADLINES: dict[str, DatasetHeadline] = {
    # Guo 2019 OpenSWATH deposited matrix (1% FDR) vs quantmsdiann.
    "PXD003539": DatasetHeadline(
        paper_count=PXD003539_PAPER_COUNT,
        diann_count=_DEFAULT_DIANN_COUNTS["PXD003539"],
        paper_label="Guo 2019 (OpenSWATH)",
        metric="Protein groups (1% global FDR, target-only)",
    ),
    "PXD030304": DatasetHeadline(
        paper_count=PXD030304_PAPER_COUNT,
        diann_count=_DEFAULT_DIANN_COUNTS["PXD030304"],
        paper_label="ProCan 2022",
        metric="Proteins @ Global.Q.Value <= 0.01 (target-only)",
    ),
    "PXD004701": DatasetHeadline(
        paper_count=PXD004701_PAPER_COUNT,
        diann_count=_DEFAULT_DIANN_COUNTS["PXD004701"],
        paper_label="Sun 2023",
        metric="Proteins (consistency filter, target-only)",
    ),
    # PXD017199 — Tognetti 2021 is a mass-cytometry study; the SWATH-MS
    # data is supplementary and the paper carries no DIA-comparable
    # protein-group headline. Panel A skips the paper bar for this entry
    # (detected via paper_count == 0).
    "PXD017199": DatasetHeadline(
        paper_count=0,
        diann_count=_DEFAULT_DIANN_COUNTS["PXD017199"],
        paper_label="",
        metric="Protein groups (1% global q-value, target-only)",
    ),
    # PXD041421 — Wang 2023 (TIQUEST batch-effect testbed). Atlas-only
    # cohort per 2026-05-21 spec §2; no paper-side DIA headline (deposit
    # is methodological).
    "PXD041421": DatasetHeadline(
        paper_count=0,
        diann_count=_DEFAULT_DIANN_COUNTS["PXD041421"],
        paper_label="",
        metric="Protein groups (1% global q-value, target-only)",
    ),
}


# Companion unfiltered protein-group counts (for the audit TSV).
# Populated by `refresh_dataset_headlines()` alongside the filtered
# headline counts.
DATASET_HEADLINES_UNFILTERED: dict[str, int] = {
    ds: _DEFAULT_DIANN_COUNTS[ds] for ds in DATASET_HEADLINES
}


def _compute_diann_counts(
    pg_matrix_path: Path | None,
    pr_matrix_path: Path | None = None,
) -> tuple[int, int]:
    """Return `(unfiltered, target_only)` Protein.Group counts for a
    cohort's headline. Prefers pg_matrix.tsv; falls back to
    pr_matrix.tsv unique Protein.Group rows when pg_matrix is absent
    (PXD003539 case). Returns `(0, 0)` when neither file exists."""
    if pg_matrix_path is not None and pg_matrix_path.exists():
        return count_target_protein_groups(pg_matrix_path)
    if pr_matrix_path is not None and pr_matrix_path.exists():
        return count_target_protein_groups_pr_matrix(pr_matrix_path)
    return (0, 0)


def refresh_dataset_headlines() -> None:
    """Populate `DATASET_HEADLINES` and `DATASET_HEADLINES_UNFILTERED`
    from the on-disk pg_matrix.tsv / pr_matrix.tsv files. Called by
    `main()` before rendering so the headline bars reflect the
    target-only contaminant-filtered counts (2026-05-21 spec §1.6).

    Tests that don't need the from-disk numbers can either skip this
    call (the module-load defaults are reasonable diannsummary.log
    baselines) or call it themselves to exercise the read path."""
    sources: dict[str, tuple[Path | None, Path | None]] = {
        "PXD003539": (None, PXD003539_PR_MATRIX),
        "PXD030304": (PXD030304_PG_MATRIX, None),
        # PXD004701 headline is the consistency-filtered union (matching the
        # per-cohort Fig. 3b), NOT the raw pg_matrix — handled separately below.
        "PXD017199": (PXD017199_PG_MATRIX, None),
        "PXD041421": (PXD041421_PG_MATRIX, None),
    }
    for ds, (pg, pr) in sources.items():
        unf, target = _compute_diann_counts(pg, pr)
        if target <= 0:
            continue  # keep the import-time default if the file is missing
        h = DATASET_HEADLINES[ds]
        DATASET_HEADLINES[ds] = DatasetHeadline(
            paper_count=h.paper_count,
            diann_count=target,
            paper_label=h.paper_label,
            metric=h.metric,
        )
        DATASET_HEADLINES_UNFILTERED[ds] = unf

    # PXD004701: use the consistency-filtered target-only protein-group union
    # from the per-subtype JSON cache, so the atlas Panel A bar equals the
    # per-cohort comparison (and the metric label "consistency filter,
    # target-only" is truthful) rather than the looser raw pg_matrix count.
    if PXD004701_PROTEIN_JSON.exists() and PXD004701_PROTEIN_JSON.stat().st_size > 0:
        with open(PXD004701_PROTEIN_JSON, encoding="utf-8") as fh:
            payload = json.load(fh)
        pgs: set[str] = set()
        for vs in payload.values():
            pgs.update(v for v in vs if isinstance(v, str))
        target_pgs = {pg for pg in pgs if is_target_protein_group(pg)}
        if target_pgs:
            h = DATASET_HEADLINES["PXD004701"]
            DATASET_HEADLINES["PXD004701"] = DatasetHeadline(
                paper_count=h.paper_count,
                diann_count=len(target_pgs),
                paper_label=h.paper_label,
                metric=h.metric,
            )
            DATASET_HEADLINES_UNFILTERED["PXD004701"] = len(pgs)


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


def _compute_runs_per_cohort() -> dict[str, int]:
    """Return the number of MS runs per cohort from the cached SDRF
    files. Used by Panel E (breadth-vs-depth scatter) to scale the dot
    size for each cohort. Missing SDRFs yield 0 (defensively — the
    panel handles zero by falling back to a constant dot size)."""
    out: dict[str, int] = {}
    for ds, sdrf in (
        ("PXD003539", PXD003539_SDRF),
        ("PXD030304", PXD030304_SDRF),
        ("PXD004701", PXD004701_SDRF),
        ("PXD017199", PXD017199_SDRF),
        ("PXD041421", PXD041421_SDRF),
    ):
        if not Path(sdrf).exists():
            out[ds] = 0
            continue
        try:
            df = pd.read_csv(sdrf, sep="\t", dtype=str)
            out[ds] = int(len(df))
        except (FileNotFoundError, OSError, ValueError, pd.errors.EmptyDataError):
            out[ds] = 0
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


def cell_line_tissue_pxd017199(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, tissue]` for PXD017199.

    PXD017199 is essentially all breast-cancer-derived plus a handful of
    "normal" mammary-epithelial lines (184A1, 184B5, HBL100, MCF10A,
    MCF10F, MCF12A). The disease column distinguishes them:

    - characteristics[disease] == "normal" -> `Healthy (Non-cancer)`
      (matches the tissue category used by `cell_line_tissue_pxd003539`
      for non-tumour rows in the unified ProCan axis).
    - any other disease value -> `Breast`.

    If a single cell line has multiple disease rows (it shouldn't in
    PXD017199 but the SDRF doesn't enforce uniqueness), the first
    matching row wins; ties between Breast and normal are dominated by
    whichever appears first in the SDRF."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    needed = ["characteristics[cell line]", "characteristics[disease]"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD017199 SDRF {sdrf_path} missing columns: {missing}"
        )
    out: dict[str, str] = {}
    for cell, disease in zip(
        df["characteristics[cell line]"], df["characteristics[disease]"],
    ):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        if norm in out:
            continue
        d = (disease or "").strip().lower()
        if d == "normal":
            out[norm] = "Healthy (Non-cancer)"
        else:
            out[norm] = "Breast"
    return out


def cell_line_tissue_pxd041421(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, tissue]` for PXD041421.

    PXD041421 (Wang 2023 batch-effect testbed) carries 2 cell lines:
    A549 (lung adenocarcinoma) and K562 (CML blast-phase). To stay on
    the unified ProCan 28-tissue axis already used by every other
    cohort, the labels match those produced by
    `cell_line_tissue_pxd003539` for the same diseases:

    - A549 -> `Lung`
    - K562 -> `Haematopoietic and Lymphoid`

    The spec calls these "Lung Cancer" and "Leukemia" — the actual
    axis labels are the ProCan headers and remain unchanged so the
    Panel C / F stacked bars merge correctly across cohorts."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    if "characteristics[cell line]" not in df.columns:
        raise ValueError(
            f"PXD041421 SDRF {sdrf_path} missing characteristics[cell line]"
        )
    # Hard-coded mapping; PXD041421 only ever carries A549 and K562.
    LINE_TO_TISSUE = {
        "A549": "Lung",
        "K562": "Haematopoietic and Lymphoid",
    }
    out: dict[str, str] = {}
    for cell in df["characteristics[cell line]"].fillna(""):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        if norm in out:
            continue
        # Look up against both the raw SDRF spelling and its
        # normalised form so an unexpected spelling (e.g. "a549")
        # still hits.
        raw = (cell or "").strip()
        tissue = LINE_TO_TISSUE.get(raw) or LINE_TO_TISSUE.get(
            normalize_cell_line(raw)
        )
        if tissue is None:
            # Defensive: leave unknown lines unmapped rather than guess.
            continue
        out[norm] = tissue
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
# Per-group accession sets (Panels E and F inputs)
# ---------------------------------------------------------------------------


def pxd003539_accessions_per_cell_line(
    pr_matrix_path: Path,
    sdrf_path: Path,
) -> dict[str, set[str]]:
    """Return `dict[normalised_cell_line, set[UniProt accession]]` for
    PXD003539. For each cell line, the value is the union of accessions
    quantified in any of its associated runs (i.e. any non-NA cell in the
    pr_matrix column matching the SDRF `comment[data file]` rewritten to
    `.mzML`).

    Unlike `pxd003539_protein_accessions`, this helper preserves
    per-cell-line granularity so Panel E (rarefaction) and Panel F
    (per-tissue protein counts) can roll up across arbitrary partitions.

    Cell lines with no matching pr_matrix column or no quantified
    precursors are dropped silently."""
    import re

    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD003539 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    sdrf = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    needed = ["characteristics[cell line]", "comment[data file]"]
    sdrf_missing = [c for c in needed if c not in sdrf.columns]
    if sdrf_missing:
        raise ValueError(
            f"PXD003539 SDRF missing required columns: {sdrf_missing}"
        )
    # Build run-column -> normalised cell-line lookup. SDRF stores .wiff;
    # pr_matrix columns are .mzML.
    col_to_cell: dict[str, str] = {}
    for cell, data_file in zip(
        sdrf["characteristics[cell line]"], sdrf["comment[data file]"],
    ):
        if not isinstance(data_file, str) or not data_file:
            continue
        mzml = re.sub(r"\.wiff$", ".mzML", data_file)
        norm = normalize_cell_line(cell)
        if norm:
            col_to_cell[mzml] = norm

    out: dict[str, set[str]] = {}
    for col in sample_cols:
        cell = col_to_cell.get(col)
        if cell is None:
            continue
        mask = df[col].notna()
        if not mask.any():
            continue
        bucket = out.setdefault(cell, set())
        for pg in df.loc[mask, "Protein.Group"]:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
    return out


def pxd017199_protein_accessions(pr_matrix_path: Path) -> set[str]:
    """Return the set of UniProt accessions detected by quantmsdiann in
    PXD017199. Mirrors `pxd003539_protein_accessions`: reads
    `diann_report.pr_matrix.tsv` (~193 MB, already on disk under
    `data/PXD017199/`) and collects accessions from `Protein.Group`
    strings on rows with at least one non-NA quantification across the
    206 runs.

    Uses `extract_accessions_diann` for accession-normalisation parity
    with the other datasets."""
    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD017199 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    if not sample_cols:
        raise ValueError("PXD017199 pr_matrix has no per-run sample columns")
    quantified = df[df[sample_cols].notna().any(axis=1)]
    accessions: set[str] = set()
    for pg in quantified["Protein.Group"]:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions


def pxd017199_accessions_per_cell_line(
    pr_matrix_path: Path,
    sdrf_path: Path,
) -> dict[str, set[str]]:
    """Return `dict[normalised_cell_line, set[UniProt accession]]` for
    PXD017199. Mirrors `pxd003539_accessions_per_cell_line` — but the
    PXD017199 run filenames are `.raw` (not `.wiff`) so no rewrite is
    needed; the SDRF's `comment[data file]` is the exact pr_matrix
    column name.

    Per-cell-line accession sets are the union across that line's runs;
    cell lines with no matching pr_matrix columns or no quantified
    precursors are silently dropped."""
    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD017199 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    sdrf = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    needed = ["characteristics[cell line]", "comment[data file]"]
    sdrf_missing = [c for c in needed if c not in sdrf.columns]
    if sdrf_missing:
        raise ValueError(
            f"PXD017199 SDRF missing required columns: {sdrf_missing}"
        )
    # Build run-column -> normalised cell-line lookup directly from the
    # SDRF; the data-file names already match pr_matrix column headers.
    col_to_cell: dict[str, str] = {}
    for cell, data_file in zip(
        sdrf["characteristics[cell line]"], sdrf["comment[data file]"],
    ):
        if not isinstance(data_file, str) or not data_file:
            continue
        norm = normalize_cell_line(cell)
        if norm:
            col_to_cell[data_file] = norm

    out: dict[str, set[str]] = {}
    for col in sample_cols:
        cell = col_to_cell.get(col)
        if cell is None:
            continue
        mask = df[col].notna()
        if not mask.any():
            continue
        bucket = out.setdefault(cell, set())
        for pg in df.loc[mask, "Protein.Group"]:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
    return out


def pxd041421_protein_accessions(pr_matrix_path: Path) -> set[str]:
    """Return the set of UniProt accessions detected by quantmsdiann in
    PXD041421. Mirrors `pxd017199_protein_accessions` — reads
    `diann_report.pr_matrix.tsv` (48 runs, ~9k protein groups, 64 MB on
    disk) and collects accessions from `Protein.Group` strings on rows
    with at least one non-NA quantification.

    Uses `extract_accessions_diann` so the conservative
    contaminant/entrapment/decoy filter is applied at the row level
    (2026-05-21 spec)."""
    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD041421 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    if not sample_cols:
        raise ValueError("PXD041421 pr_matrix has no per-run sample columns")
    quantified = df[df[sample_cols].notna().any(axis=1)]
    accessions: set[str] = set()
    for pg in quantified["Protein.Group"]:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions


def pxd041421_accessions_per_cell_line(
    pr_matrix_path: Path,
    sdrf_path: Path,
) -> dict[str, set[str]]:
    """Return `dict[normalised_cell_line, set[UniProt accession]]` for
    PXD041421. Mirrors `pxd017199_accessions_per_cell_line`.

    Run-column / data-file alignment quirks for PXD041421:
    - pr_matrix column headers end in `.d` (Bruker timsTOF folder)
    - SDRF `comment[data file]` ends in `.d.zip` (the FTP archive name).
    We strip the trailing `.zip` to map SDRF -> pr_matrix column.

    Per-cell-line accession sets are the union across that line's runs;
    cell lines with no matching pr_matrix columns or no quantified
    precursors are silently dropped."""
    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD041421 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    sdrf = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    needed = ["characteristics[cell line]", "comment[data file]"]
    sdrf_missing = [c for c in needed if c not in sdrf.columns]
    if sdrf_missing:
        raise ValueError(
            f"PXD041421 SDRF missing required columns: {sdrf_missing}"
        )
    import re as _re
    col_to_cell: dict[str, str] = {}
    for cell, data_file in zip(
        sdrf["characteristics[cell line]"], sdrf["comment[data file]"],
    ):
        if not isinstance(data_file, str) or not data_file:
            continue
        # SDRF carries `*.d.zip`; pr_matrix headers carry `*.d`.
        col_name = _re.sub(r"\.zip$", "", data_file, flags=_re.IGNORECASE)
        norm = normalize_cell_line(cell)
        if norm:
            col_to_cell[col_name] = norm

    out: dict[str, set[str]] = {}
    for col in sample_cols:
        cell = col_to_cell.get(col)
        if cell is None:
            continue
        mask = df[col].notna()
        if not mask.any():
            continue
        bucket = out.setdefault(cell, set())
        for pg in df.loc[mask, "Protein.Group"]:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
    return out


def pxd003539_gene_symbols(pr_matrix_path: Path) -> set[str]:
    """Return the set of gene symbols quantified in PXD003539 by
    quantmsdiann. Each `Genes` cell may carry one or more ';'-separated
    symbols; only rows with at least one non-NA quantification across the
    runs are counted (same filter as `pxd003539_protein_accessions`)."""
    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD003539 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    quantified = df[df[sample_cols].notna().any(axis=1)]
    out: set[str] = set()
    for g in quantified["Genes"]:
        if not isinstance(g, str):
            continue
        for piece in g.split(";"):
            sym = piece.strip()
            if sym:
                out.add(sym)
    return out


def _per_group_accessions_from_json(json_path: Path) -> dict[str, set[str]]:
    """Read `{group_key: [Protein.Group, ...]}` JSON cache and return
    `{group_key: set[UniProt accession]}`."""
    with open(json_path, encoding="utf-8") as fh:
        payload = json.load(fh)
    out: dict[str, set[str]] = {}
    for key, pgs in payload.items():
        bucket: set[str] = set()
        for pg in pgs:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
        out[key] = bucket
    return out


def pxd030304_accessions_per_tissue(json_path: Path) -> dict[str, set[str]]:
    """Return `{tissue: set[accession]}` for PXD030304 (28 tissues)."""
    return _per_group_accessions_from_json(json_path)


def pxd004701_accessions_per_subtype(json_path: Path) -> dict[str, set[str]]:
    """Return `{subtype: set[accession]}` for PXD004701 (3 BC subtypes)."""
    return _per_group_accessions_from_json(json_path)


# ---------------------------------------------------------------------------
# Expression Atlas catalogue (Panel H input)
# ---------------------------------------------------------------------------


def expression_atlas_gene_set(tsv_path: Path) -> set[str]:
    """Parse the Expression Atlas E-PROT-73 query-results TSV and return the
    set of unique `Gene Name` values. The file's first lines are
    `#`-prefixed comments and `pd.read_csv(..., comment='#')` skips them.

    Returns an empty set if the file is missing — callers (Panel H) treat
    that as an explainer-only render path."""
    if not tsv_path.exists():
        return set()
    df = pd.read_csv(tsv_path, sep="\t", comment="#", dtype=str)
    if "Gene Name" not in df.columns:
        raise ValueError(
            f"E-PROT-73 TSV {tsv_path} missing required column 'Gene Name'"
        )
    out: set[str] = set()
    for g in df["Gene Name"].fillna(""):
        sym = g.strip()
        if sym:
            out.add(sym)
    return out


# ---------------------------------------------------------------------------
# Panel E / F / G helpers (pure data transforms — separately testable)
# ---------------------------------------------------------------------------


def rarefaction_curve(
    groups: dict[str, set[str]],
    *,
    n_permutations: int = 50,
    seed: int = 42,
) -> list[float]:
    """Return the average cumulative-union size as groups are accumulated
    one-by-one. Order is randomised over `n_permutations` permutations
    using `numpy.random.default_rng(seed)` for reproducibility.

    The returned list has length `len(groups)`; entry `i` is the mean
    size of the union over the first `i+1` groups across permutations.
    For 0 or 1 groups the curve is deterministic (no averaging needed)."""
    import numpy as np

    keys = list(groups.keys())
    n = len(keys)
    if n == 0:
        return []
    if n == 1:
        return [float(len(groups[keys[0]]))]
    rng = np.random.default_rng(seed)
    accum = np.zeros(n, dtype=float)
    for _ in range(n_permutations):
        order = list(rng.permutation(n))
        union: set[str] = set()
        for i, idx in enumerate(order):
            union.update(groups[keys[idx]])
            accum[i] += len(union)
    return list(accum / n_permutations)


def per_tissue_union_accessions(
    tissue_order: list[str],
    cl_tissue_pxd003539: dict[str, str],
    accessions_per_cell_line_pxd003539: dict[str, set[str]],
    accessions_per_tissue_pxd030304: dict[str, set[str]],
    accessions_pxd004701: set[str],
    cl_tissue_pxd017199: dict[str, str] | None = None,
    accessions_per_cell_line_pxd017199: dict[str, set[str]] | None = None,
    cl_tissue_pxd041421: dict[str, str] | None = None,
    accessions_per_cell_line_pxd041421: dict[str, set[str]] | None = None,
) -> list[tuple[str, set[str]]]:
    """Per-tissue UNION of UniProt accessions across all cohorts.

    Returns `[(tissue, union_accession_set)]` in the input order. The
    union answers `how many distinct target-only proteins were
    detected in this tissue across the whole atlas`, which is the
    biologically meaningful number — the previous `per_tissue_protein_counts`
    summed per-cohort counts and so double-counted proteins detected
    in multiple cohorts (e.g., a Breast protein seen in PXD004701,
    PXD017199 and PXD003539 contributed three times to the stack)."""
    # Mirror per_tissue_protein_counts's input-shape contract; only
    # the aggregation step (sum → union) differs.
    tissue_to_cls: dict[str, list[str]] = {}
    for cl, tissue in cl_tissue_pxd003539.items():
        tissue_to_cls.setdefault(tissue, []).append(cl)
    tissue_to_cls_17199: dict[str, list[str]] = {}
    have_17199 = (
        cl_tissue_pxd017199 is not None
        and accessions_per_cell_line_pxd017199 is not None
    )
    if have_17199:
        for cl, tissue in cl_tissue_pxd017199.items():
            tissue_to_cls_17199.setdefault(tissue, []).append(cl)
    tissue_to_cls_41421: dict[str, list[str]] = {}
    have_41421 = (
        cl_tissue_pxd041421 is not None
        and accessions_per_cell_line_pxd041421 is not None
    )
    if have_41421:
        for cl, tissue in cl_tissue_pxd041421.items():
            tissue_to_cls_41421.setdefault(tissue, []).append(cl)

    rows: list[tuple[str, set[str]]] = []
    for tissue in tissue_order:
        union: set[str] = set()
        for cl in tissue_to_cls.get(tissue, []):
            union.update(accessions_per_cell_line_pxd003539.get(cl, set()))
        union.update(accessions_per_tissue_pxd030304.get(tissue, set()))
        if tissue == "Breast":
            union.update(accessions_pxd004701)
        if have_17199:
            for cl in tissue_to_cls_17199.get(tissue, []):
                union.update(
                    accessions_per_cell_line_pxd017199.get(cl, set())
                )
        if have_41421:
            for cl in tissue_to_cls_41421.get(tissue, []):
                union.update(
                    accessions_per_cell_line_pxd041421.get(cl, set())
                )
        rows.append((tissue, union))
    return rows


def per_tissue_protein_counts(
    tissue_order: list[str],
    cl_tissue_pxd003539: dict[str, str],
    accessions_per_cell_line_pxd003539: dict[str, set[str]],
    accessions_per_tissue_pxd030304: dict[str, set[str]],
    accessions_pxd004701: set[str],
    cl_tissue_pxd017199: dict[str, str] | None = None,
    accessions_per_cell_line_pxd017199: dict[str, set[str]] | None = None,
    cl_tissue_pxd041421: dict[str, str] | None = None,
    accessions_per_cell_line_pxd041421: dict[str, set[str]] | None = None,
) -> list[tuple[str, dict[str, int]]]:
    """Roll up per-(tissue, dataset) protein-set sizes onto the unified
    tissue axis (the same `tissue_order` used by Panel C).

    For each tissue:
      - PXD003539: union of per-cell-line accession sets for cell lines
        mapped to this tissue.
      - PXD030304: size of the cached per-tissue accession set.
      - PXD004701: total accession count if and only if the tissue is
        `Breast` (the dataset is breast-only); else 0.
      - PXD017199 (optional, only when both mappings are passed): union
        of per-cell-line accessions for cell lines mapped to this
        tissue (so the 5-6 mammary-normal lines contribute to
        `Healthy (Non-cancer)` and the rest to `Breast`).

    Returns `[(tissue, {dataset: protein_count})]` in the input order.
    The PXD017199 entry is only present in the inner dict when both
    `cl_tissue_pxd017199` and `accessions_per_cell_line_pxd017199` are
    non-None — keeps the function back-compatible with the 3-dataset
    callers in older tests."""
    # Pre-invert PXD003539 cell-line -> tissue to tissue -> cell lines
    tissue_to_cls: dict[str, list[str]] = {}
    for cl, tissue in cl_tissue_pxd003539.items():
        tissue_to_cls.setdefault(tissue, []).append(cl)

    # Same for PXD017199 if provided
    tissue_to_cls_17199: dict[str, list[str]] = {}
    have_17199 = (
        cl_tissue_pxd017199 is not None
        and accessions_per_cell_line_pxd017199 is not None
    )
    if have_17199:
        for cl, tissue in cl_tissue_pxd017199.items():
            tissue_to_cls_17199.setdefault(tissue, []).append(cl)

    # Same for PXD041421 (2026-05-21 spec §2: A549 + K562 -> Lung +
    # Haematopoietic and Lymphoid).
    tissue_to_cls_41421: dict[str, list[str]] = {}
    have_41421 = (
        cl_tissue_pxd041421 is not None
        and accessions_per_cell_line_pxd041421 is not None
    )
    if have_41421:
        for cl, tissue in cl_tissue_pxd041421.items():
            tissue_to_cls_41421.setdefault(tissue, []).append(cl)

    rows: list[tuple[str, dict[str, int]]] = []
    for tissue in tissue_order:
        # PXD003539: union per tissue
        cls = tissue_to_cls.get(tissue, [])
        u_3539: set[str] = set()
        for cl in cls:
            u_3539.update(accessions_per_cell_line_pxd003539.get(cl, set()))
        # PXD030304: direct lookup
        n_30304 = len(accessions_per_tissue_pxd030304.get(tissue, set()))
        # PXD004701: only Breast.
        n_4701 = len(accessions_pxd004701) if tissue == "Breast" else 0
        by_ds = {
            "PXD003539": len(u_3539),
            "PXD030304": n_30304,
            "PXD004701": n_4701,
        }
        if have_17199:
            u_17199: set[str] = set()
            for cl in tissue_to_cls_17199.get(tissue, []):
                u_17199.update(
                    accessions_per_cell_line_pxd017199.get(cl, set())
                )
            by_ds["PXD017199"] = len(u_17199)
        if have_41421:
            u_41421: set[str] = set()
            for cl in tissue_to_cls_41421.get(tissue, []):
                u_41421.update(
                    accessions_per_cell_line_pxd041421.get(cl, set())
                )
            by_ds["PXD041421"] = len(u_41421)
        rows.append((tissue, by_ds))
    return rows


def detection_count_histogram(
    accession_sets: dict[str, set[str]],
) -> dict[int, int]:
    """For every accession in the union of the input sets, count how many
    sets it occurs in, then bucket sizes by that count (1, 2, ..., N).

    Returns `{count: bucket_size}` for `count` in `1..N` (N = number of
    input sets). Missing counts get a 0 entry so the histogram always
    has N bins."""
    n_sets = len(accession_sets)
    union: set[str] = set()
    for s in accession_sets.values():
        union.update(s)
    buckets: dict[int, int] = {k: 0 for k in range(1, n_sets + 1)}
    for acc in union:
        c = sum(1 for s in accession_sets.values() if acc in s)
        if c >= 1:
            buckets[c] = buckets.get(c, 0) + 1
    return buckets


# ---------------------------------------------------------------------------
# Panel renderers
# ---------------------------------------------------------------------------


def _annotate_panel_letter(
    ax, letter: str, *, subtitle: str | None = None,
) -> None:
    """Stamp a bold panel letter in the top-left corner. When
    `subtitle` is provided, render a smaller description next to the
    letter so each panel is self-naming (the bold letter + a short
    description of what the panel shows). Subtitles are paper-ready:
    short noun phrases, no terminal punctuation, fontsize 10.
    """
    ax.text(
        -0.07, 1.02, letter,
        transform=ax.transAxes, fontsize=14, fontweight="bold",
        ha="left", va="bottom",
    )
    if subtitle:
        ax.text(
            -0.02, 1.02, subtitle,
            transform=ax.transAxes, fontsize=10, fontweight="normal",
            ha="left", va="bottom", color="#222222",
        )


def _render_panel_a_headlines(ax, headlines: dict[str, DatasetHeadline]) -> None:
    """Grouped bar chart, N dataset groups × up-to-2 bars (paper vs
    quantmsdiann). Datasets whose `paper_count == 0` get only the
    quantmsdiann bar, centred on the x-tick — no explainer annotation;
    the cohort's x-axis label (PXDxxx + paper-year tag from
    `DATASET_LABELS`) is sufficient context. Bars are colour-coded paper
    (grey) vs quantmsdiann (blue); value labels sit just above each bar
    top with a small fractional pad so they never touch the bar."""
    datasets = list(headlines.keys())
    bar_width = 0.36
    x = list(range(len(datasets)))
    paper_drawn_xs: list[float] = []
    paper_drawn_vals: list[int] = []
    diann_xs: list[float] = []
    diann_vals_drawn: list[int] = []
    for xi, d in zip(x, datasets):
        h = headlines[d]
        if h.paper_count > 0:
            paper_drawn_xs.append(xi - bar_width / 2)
            paper_drawn_vals.append(h.paper_count)
            diann_xs.append(xi + bar_width / 2)
        else:
            # Centre the quantmsdiann bar on the x-tick when no paper bar.
            diann_xs.append(xi)
        diann_vals_drawn.append(h.diann_count)
    bars_p = ax.bar(
        paper_drawn_xs, paper_drawn_vals,
        width=bar_width, color="#9e9e9e", label="Original analysis",
    )
    bars_d = ax.bar(
        diann_xs, diann_vals_drawn,
        width=bar_width, color="#1f77b4", label="quantmsdiann (DIA-NN)",
    )
    all_vals = [h.diann_count for h in headlines.values()] + [
        h.paper_count for h in headlines.values() if h.paper_count > 0
    ]
    ymax = (max(all_vals) if all_vals else 1) * 1.18
    pad = ymax * 0.02
    for bar, v in zip(bars_p, paper_drawn_vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + pad,
            f"{v:,}", ha="center", va="bottom", fontsize=8,
        )
    for bar, v in zip(bars_d, diann_vals_drawn):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + pad,
            f"{v:,}", ha="center", va="bottom", fontsize=8,
        )
    ax.set_xticks(x)
    # X-tick labels: rely on the canonical `DATASET_LABELS` (PXDxxx +
    # paper-year line) so cohorts with no paper bar are still
    # identifiable without any explanatory annotation.
    ax.set_xticklabels(
        [DATASET_LABELS.get(d, d) for d in datasets],
        fontsize=8,
    )
    ax.set_ylabel("Protein groups")
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Legend below the x-axis labels so the two-entry key (paper bar
    # vs quantmsdiann bar) never overlaps a tall cohort bar regardless
    # of which cohort has the highest count. The cohort identities are
    # already on the x-tick labels.
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        frameon=False, fontsize=8, ncol=2,
    )


def _render_upset_in_axes(
    ax,
    sets: dict[str, set[str]],
    ds_order: list[str],
    *,
    bar_color: str = "#445566",
    panel_title: str | None = None,
    panel_letter: str | None = None,
    panel_subtitle: str | None = None,
) -> None:
    """Render an UpSet plot inside an axes slot via a Matplotlib SubFigure.

    `upsetplot.UpSet.plot()` creates its own gridspec on the given
    figure (it does not accept a target Axes). To slot it into a
    pre-existing 4x2 grid cell we:
      1. Hide the host axes `ax`.
      2. Carve a SubFigure out of the host figure at the same
         SubplotSpec as `ax`.
      3. Pass that SubFigure to `upset.plot(fig=subfig)`. UpSet's
         internal gridspec is anchored to the subfig, not the parent
         figure, so it stays inside the panel cell.

    Falls back to a brief text annotation if fewer than 2 datasets have
    content (UpSet needs at least 2 categories)."""
    populated = [d for d in ds_order if sets.get(d)]
    if len(populated) < 2:
        ax.text(
            0.5, 0.5,
            "Insufficient data for UpSet plot\n"
            f"(populated datasets: {len(populated)})",
            ha="center", va="center", transform=ax.transAxes, fontsize=9,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return
    try:
        from upsetplot import UpSet, from_contents
    except ImportError as exc:  # pragma: no cover - dependency required
        raise RuntimeError(
            "upsetplot required for Panel B / D — install via "
            "`pip install upsetplot` (already in analysis/requirements.txt)"
        ) from exc

    fig = ax.figure
    subplotspec = ax.get_subplotspec()

    # Hide the placeholder axes; the SubFigure paints over its cell.
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("none")
    ax.set_visible(False)

    contents = {d: sorted(sets[d]) for d in populated}
    data = from_contents(contents)

    subfig = fig.add_subfigure(subplotspec)

    upset = UpSet(
        data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        sort_categories_by="cardinality",
        facecolor=bar_color,
        # element_size=None makes UpSet's make_grid use the
        # get_window_extent path (which SubFigure supports) instead of
        # the get_figwidth path (which it does not).
        element_size=None,
    )
    upset.plot(fig=subfig)
    if panel_title:
        subfig.suptitle(panel_title, fontsize=10)
    # The host `ax` is hidden so `_annotate_panel_letter(ax, ...)` would
    # be invisible. Stamp the panel letter + subtitle on the SubFigure
    # directly so the panel reads "A — Protein overlap" like its siblings.
    if panel_letter:
        subfig.text(
            0.005, 0.985, panel_letter,
            fontsize=14, fontweight="bold",
            ha="left", va="top",
        )
        if panel_subtitle:
            subfig.text(
                0.045, 0.985, panel_subtitle,
                fontsize=10, ha="left", va="top", color="#222222",
            )


def _render_panel_b_cellline_venn(
    ax,
    sets: dict[str, set[str]],
    region_sizes: dict[str, int] | None = None,
) -> None:
    """UpSet plot of normalised cell-line names across the 4 datasets.

    Replaces the original 3-set venn3 rendering; with PXD017199 added,
    the 4-way cell-line overlap (especially PXD017199-vs-PXD004701
    breast overlap) is unreadable in a 4-set Venn. UpSet renders one
    bar per intersection region, ordered by cardinality."""
    ds_order = [
        d for d in ("PXD003539", "PXD030304", "PXD004701", "PXD017199", "PXD041421")
        if d in sets
    ]
    _render_upset_in_axes(
        ax, sets, ds_order, bar_color="#445566",
        panel_title="Cell-line set intersections (UpSet)",
    )


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
    ds_order = ["PXD003539", "PXD030304", "PXD004701", "PXD017199", "PXD041421"]
    # Only stack datasets that actually contribute something in any row,
    # so panel C remains readable when called with 3-dataset inputs.
    ds_order = [
        ds for ds in ds_order
        if any(r[1].get(ds, 0) > 0 for r in rows)
    ]
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
                ha="left", va="center", fontsize=6)
    ax.set_yticks(y_top)
    ax.set_yticklabels(tissues, fontsize=6)
    ax.set_xlabel("Cell lines (sum across datasets)", fontsize=9)
    ax.set_xlim(0, max(totals) * 1.10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=7)


def _render_panel_d_protein_venn(
    ax,
    sets: dict[str, set[str]],
    *,
    panel_letter: str | None = None,
    panel_subtitle: str | None = None,
) -> None:
    """UpSet plot of UniProt accessions across the 4 quantmsdiann
    analyses. Replaces the original 3-set venn3 rendering for the same
    reason as Panel B — the 4-set Venn is unreadable.

    `panel_letter` / `panel_subtitle` flow through to the SubFigure
    annotation (since the UpSet renderer hides the host axes, normal
    `_annotate_panel_letter(ax, ...)` calls go to an invisible Axes
    and never appear on the rendered SVG)."""
    ds_order = [
        d for d in ("PXD003539", "PXD030304", "PXD004701", "PXD017199", "PXD041421")
        if d in sets
    ]
    _render_upset_in_axes(
        ax, sets, ds_order, bar_color="#445566",
        panel_letter=panel_letter, panel_subtitle=panel_subtitle,
    )


def _render_panel_e_breadth_vs_depth(
    ax,
    cellline_sets: dict[str, set[str]],
    accession_sets: dict[str, set[str]],
    *,
    runs_per_cohort: dict[str, int] | None = None,
) -> None:
    """Per-cohort breadth-vs-depth scatter. One dot per cohort:
      - x = number of distinct cell lines in the cohort SDRF
        (`len(cellline_sets[ds])`), log scale (range 2 → 947)
      - y = target-only union of UniProt accessions
        (`len(accession_sets[ds])`)
      - dot size proportional to the number of MS runs in the cohort
        (per `runs_per_cohort`, if provided)
      - colour = `DATASET_COLORS[ds]`

    Replaces the previous rarefaction curves whose x-axes were
    incompatible across cohorts (cell-lines / tissues / subtypes).
    Single-glance view of where each cohort sits on the
    breadth↔depth tradeoff: PXD030304 is broad (947 lines) but
    matched to a deep proteome; PXD041421 is narrow (2 lines) but
    deep (24 reps each); PXD003539 / PXD017199 / PXD004701 sit in
    the middle. Cohorts missing from either map are silently
    omitted."""
    ds_order = ["PXD003539", "PXD030304", "PXD004701",
                "PXD017199", "PXD041421"]
    any_plotted = False
    xs_all: list[float] = []
    ys_all: list[float] = []
    for ds in ds_order:
        n_cl = len(cellline_sets.get(ds, set()))
        n_acc = len(accession_sets.get(ds, set()))
        if n_cl == 0 or n_acc == 0:
            continue
        n_runs = (runs_per_cohort or {}).get(ds, 0) or 0
        # Dot area scales with runs; clamp to a readable range.
        # No runs metadata → constant medium size.
        if n_runs > 0:
            # min runs ~48 (PXD041421), max ~5800 (PXD030304)
            size = 60 + 240 * (n_runs ** 0.5) / (6000 ** 0.5)
        else:
            size = 120
        ax.scatter(
            [n_cl], [n_acc],
            s=size, c=DATASET_COLORS[ds],
            edgecolors="#222222", linewidths=0.6,
            zorder=3,
            label=(
                f"{ds} ({n_runs:,} runs)" if n_runs > 0
                else ds
            ),
        )
        # Per-cohort annotation above the dot.
        ax.annotate(
            DATASET_LABELS[ds].replace("\n", " "),
            xy=(n_cl, n_acc),
            xytext=(6, 6), textcoords="offset points",
            fontsize=8, color="#222222",
        )
        xs_all.append(n_cl)
        ys_all.append(n_acc)
        any_plotted = True

    if not any_plotted:
        ax.text(0.5, 0.5, "no cohort inputs",
                ha="center", va="center", transform=ax.transAxes)
        return

    ax.set_xscale("log")
    ax.set_xlabel("Cell lines per cohort", fontsize=10)
    ax.set_ylabel("Proteins (union per cohort)", fontsize=10)
    # Pad axes so the per-cohort labels don't fall off the edges.
    if xs_all and ys_all:
        ax.set_xlim(min(xs_all) * 0.4, max(xs_all) * 3.0)
        ax.set_ylim(0, max(ys_all) * 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=8)
    # Dot size encodes √(MS runs) when run-count metadata is provided;
    # the encoding is described in the manuscript methods, not on the
    # figure itself.
    ax.legend(loc="lower right", frameon=False, fontsize=7)


def _render_panel_f_tissue_protein_counts(
    ax,
    rows: list[tuple[str, set[str]]] | list[tuple[str, dict[str, int]]],
) -> None:
    """Per-tissue UNIQUE-protein bars. Single bar per tissue showing
    the union of UniProt accessions across all contributing cohorts —
    the biologically meaningful number (`how many distinct proteins did
    the whole atlas observe in this tissue?`).

    Replaces the previous stacked-bar that summed per-cohort counts
    and so double-counted proteins seen in multiple cohorts. The
    cohort-breakdown is still preserved in the audit TSV
    (`atlas_distribution | Panel B | tissue | <cohort> = <count>`).

    `rows` accepts two shapes for backward compatibility:
      - `[(tissue, set[str])]` — the new union-set shape (preferred).
      - `[(tissue, {dataset: count})]` — the legacy per-cohort dict
        (which we sum, replicating the old visual). Callers should
        migrate to the set-of-accessions form.
    Same tissue ordering as Panel A (per-tissue cell lines) — caller
    passes the already-ordered list."""
    if not rows:
        ax.text(0.5, 0.5, "no tissue/protein rows",
                ha="center", va="center", transform=ax.transAxes)
        return
    tissues = [r[0] for r in rows]
    # Detect the input shape.
    sample_value = rows[0][1]
    if isinstance(sample_value, set):
        counts = [len(r[1]) for r in rows]
        x_label = "Unique proteins per tissue"
    else:
        # Legacy: dict of per-cohort counts → sum.
        counts = [sum(r[1].values()) for r in rows]
        x_label = "Proteins per tissue (legacy sum)"
    n = len(tissues)
    y_top = list(reversed(range(n)))
    max_total = max(counts) if counts else 1
    bar_color = "#7570b3"  # match the ProCan accent used in the atlas
    ax.barh(
        y_top, counts, color=bar_color,
        edgecolor="white", linewidth=0.4,
    )
    for yi, t in zip(y_top, counts):
        ax.text(t + max_total * 0.005, yi, f"{t:,}",
                ha="left", va="center", fontsize=6)
    ax.set_yticks(y_top)
    ax.set_yticklabels(tissues, fontsize=6)
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_xlim(0, max_total * 1.12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _render_panel_g_detection_histogram(
    ax,
    accession_sets: dict[str, set[str]],
) -> None:
    """3-bar (or N-bar) chart of how many accessions are seen in 1, 2 or
    3 datasets. Bars annotated with absolute count and percentage of the
    union. Title-style headline above the bars."""
    buckets = detection_count_histogram(accession_sets)
    n_sets = len(accession_sets)
    if n_sets == 0 or not buckets:
        ax.text(0.5, 0.5, "no detection inputs",
                ha="center", va="center", transform=ax.transAxes)
        return
    xs = sorted(buckets.keys())
    vals = [buckets[k] for k in xs]
    total = sum(vals)
    palette = ["#9e9e9e", "#1f77b4", "#1b7a3a", "#d95f02", "#8e44ad"]
    bars = ax.bar(
        xs, vals,
        color=palette[: len(xs)],
        edgecolor="black", linewidth=0.4,
    )
    ymax = max(vals) * 1.22 if max(vals) > 0 else 1
    pad = ymax * 0.02
    for bar, v in zip(bars, vals):
        pct = (100.0 * v / total) if total > 0 else 0.0
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + pad,
            f"{v:,}\n({pct:.1f}%)",
            ha="center", va="bottom", fontsize=7,
        )
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{k}" for k in xs])
    ax.set_xlabel("Datasets a protein is detected in", fontsize=9)
    ax.set_ylabel("UniProt accessions", fontsize=9)
    ax.set_title(f"Pan-cohort core ({n_sets} datasets)", fontsize=10)
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _render_panel_h_expression_atlas_overlap(
    ax,
    ea_genes: set[str],
    diann_genes: set[str],
) -> None:
    """Expression Atlas (E-PROT-73 NCI-60) vs PXD003539 quantmsdiann gene
    overlap. Renders 3 bars: EA catalogue, quantmsdiann, intersection;
    annotates intersection coverage as a percentage of EA.

    If `ea_genes` is empty (file missing or empty), renders an explainer
    panel and returns without drawing bars."""
    if not ea_genes:
        ax.text(
            0.5, 0.5,
            "Expression Atlas catalogue (E-PROT-73) unavailable\n"
            "Place data/E-PROT-73-query-results.tsv to populate this panel.",
            ha="center", va="center", transform=ax.transAxes, fontsize=9,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return
    inter = ea_genes & diann_genes
    labels = [
        "Expression Atlas\n(NCI-60 catalogue)",
        "quantmsdiann\n(PXD003539)",
        "Intersection",
    ]
    vals = [len(ea_genes), len(diann_genes), len(inter)]
    colors = ["#9e9e9e", DATASET_COLORS["PXD003539"], "#1f77b4"]
    xs = list(range(3))
    bars = ax.bar(xs, vals, color=colors, edgecolor="black", linewidth=0.4)
    for bar, v in zip(bars, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height(),
            f"{v:,}", ha="center", va="bottom", fontsize=8,
        )
    coverage_pct = (100.0 * len(inter) / len(ea_genes)) if ea_genes else 0.0
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Gene symbols", fontsize=9)
    ax.set_title(
        f"PXD003539 covers {coverage_pct:.1f}% of E-PROT-73 gene catalogue",
        fontsize=10,
    )
    ymax = max(vals) * 1.18 if max(vals) > 0 else 1
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def render_atlas_overlap(
    headlines: dict[str, DatasetHeadline],
    cellline_sets: dict[str, set[str]],
    accession_sets: dict[str, set[str]],
    svg_path: Path,
) -> None:
    """Compose the cohort-overlap atlas figure: 3 panels relabelled
    A/B/C after the 2026-05-21 split:

      Row 1 (full width): A — protein-accession UpSet (5 cohorts)
      Row 2 left:         B — per-cohort headline counts
      Row 2 right:        C — pan-cohort detection histogram

    Panel A (the UpSet) gets the wide top row because the 5-set UpSet
    needs horizontal space for both the matrix and the intersection-bar
    chart. The two summary panels (B, C) sit alongside each other below
    where their narrower aspect ratio works.

    `cellline_sets` is retained in the signature for backwards
    compatibility — Panel B-as-cell-line-UpSet was removed in the
    2026-05-21 cleanup because PXD030304's 947 lines dominated the
    inter-cohort intersections."""
    del cellline_sets  # kept in signature for backwards-compatibility
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[1.4, 1.0],
        hspace=0.30, wspace=0.30,
    )
    ax_a = fig.add_subplot(gs[0, :])     # top row spans both columns
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # Panel A: UpSet renders on a SubFigure that hides the host `ax_a`.
    # Stamp the letter/subtitle inside the renderer (via the SubFigure)
    # rather than calling `_annotate_panel_letter(ax_a, ...)` — that
    # would write into the hidden host Axes and disappear from the SVG.
    _render_panel_d_protein_venn(
        ax_a, accession_sets,
        panel_letter="A", panel_subtitle="Protein overlap",
    )

    _render_panel_a_headlines(ax_b, headlines)
    _annotate_panel_letter(ax_b, "B", subtitle="Headline counts")

    _render_panel_g_detection_histogram(ax_c, accession_sets)
    _annotate_panel_letter(ax_c, "C", subtitle="Detection counts")

    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def render_atlas_distribution(
    tissue_rows: list[tuple[str, dict[str, int]]],
    cellline_sets: dict[str, set[str]],
    accession_sets: dict[str, set[str]],
    svg_path: Path,
    *,
    tissue_protein_rows: list[tuple[str, dict[str, int]]] | None = None,
    runs_per_cohort: dict[str, int] | None = None,
) -> None:
    """Compose the per-tissue distribution atlas figure: two stacked
    per-tissue bar panels in a 2x1 grid (figsize 14x14). Layout:

        Row 1: A (wide, per-tissue cell-line stacked bar)
        Row 2: B (wide, per-tissue unique-protein stacked bar)

    Panel H (Expression Atlas vs PXD003539 gene overlap) was removed
    from the atlas because it duplicates `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg`
    — the same NCI-60 gene-set comparison already lives in the
    PXD003539 per-cohort figure where it makes more sense.

    The breadth-vs-depth scatter (former Panel C) was dropped 2026-05-29:
    it crowded the figure and the breadth↔depth tradeoff is already
    legible from the two per-tissue stacked bars (PXD030304's broad,
    shallow per-tissue spread vs the narrow, deep single-line cohorts).
    `cellline_sets`, `accession_sets`, and `runs_per_cohort` are retained
    in the signature for API stability with `render_atlas` but are no
    longer used here.

    All extended inputs default to empty containers so the renderer is
    resilient on partial-data runs."""
    del cellline_sets, accession_sets, runs_per_cohort  # unused since 2026-05-29
    tissue_protein_rows = tissue_protein_rows or []

    fig = plt.figure(figsize=(14, 14))
    gs = fig.add_gridspec(
        2, 1,
        hspace=0.40,
    )
    # Re-lettered 2026-05-29: atlas_distribution now carries two panels,
    # A (per-tissue cell-lines, top) and B (per-tissue unique proteins,
    # bottom). The former breadth-vs-depth scatter (Panel C) was removed.
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])

    _render_panel_c_tissue_coverage(ax_a, tissue_rows)
    _annotate_panel_letter(ax_a, "A", subtitle="Cell lines per tissue")

    _render_panel_f_tissue_protein_counts(ax_b, tissue_protein_rows)
    _annotate_panel_letter(ax_b, "B", subtitle="Proteins per tissue")

    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def render_atlas_main(
    headlines: dict[str, DatasetHeadline],
    tissue_rows: list[tuple[str, dict[str, int]]],
    svg_path: Path,
    *,
    tissue_protein_rows: list[tuple[str, dict[str, int]]] | None = None,
) -> None:
    """Compose the main pan-cohort figure as three stacked, full-width
    panels (drops the protein-accession UpSet overlap and the
    detection-count histogram, which were too small to read in the
    multi-panel layout):

        A (top):    per-cohort headline counts (paper vs quantmsdiann)
        B (middle): per-tissue cell-line coverage (stacked bars)
        C (bottom): per-tissue unique-protein coverage (stacked bars)

    Removing the UpSet panel also removes the upsetplot/numpy-2 render
    dependency, so this figure renders on any supported NumPy."""
    tissue_protein_rows = tissue_protein_rows or []
    fig = plt.figure(figsize=(12, 15))
    gs = fig.add_gridspec(
        3, 1, height_ratios=[0.7, 1.15, 1.15], hspace=0.34,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[2, 0])

    _render_panel_a_headlines(ax_a, headlines)
    _annotate_panel_letter(ax_a, "A", subtitle="Headline counts (paper vs quantmsdiann)")

    _render_panel_c_tissue_coverage(ax_b, tissue_rows)
    _annotate_panel_letter(ax_b, "B", subtitle="Cell lines per tissue")

    _render_panel_f_tissue_protein_counts(ax_c, tissue_protein_rows)
    _annotate_panel_letter(ax_c, "C", subtitle="Proteins per tissue")

    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def render_atlas(
    headlines: dict[str, DatasetHeadline],
    cellline_sets: dict[str, set[str]],
    tissue_rows: list[tuple[str, dict[str, int]]],
    accession_sets: dict[str, set[str]],
    svg_path: Path,
    *,
    tissue_protein_rows: list[tuple[str, dict[str, int]]] | None = None,
    runs_per_cohort: dict[str, int] | None = None,
) -> None:
    """Thin wrapper that renders both the overlap (`atlas_overlap.svg`,
    panels A/B/D/G) and the distribution (`atlas_distribution.svg`,
    panels C/E/F) atlas figures.

    The single `svg_path` argument is interpreted as a hint: the actual
    outputs are written next to it as `<stem>_overlap.svg` and
    `<stem>_distribution.svg`. Panel H (Expression Atlas overlap) was
    removed from the atlas — it lives in
    `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg`."""
    out_dir = svg_path.parent
    stem = svg_path.stem
    overlap_path = out_dir / f"{stem}_overlap.svg"
    distribution_path = out_dir / f"{stem}_distribution.svg"
    render_atlas_overlap(
        headlines, cellline_sets, accession_sets, overlap_path,
    )
    render_atlas_distribution(
        tissue_rows, cellline_sets, accession_sets, distribution_path,
        tissue_protein_rows=tissue_protein_rows,
        runs_per_cohort=runs_per_cohort,
    )


# ---------------------------------------------------------------------------
# counts.tsv writer
# ---------------------------------------------------------------------------


def _set_region_sizes(
    sets: dict[str, set[str]],
    ds_order: list[str],
) -> dict[str, int]:
    """Return all 2^n - 1 non-empty-membership-pattern region sizes
    across n sets in `ds_order`.

    Region-key shapes (matched by the counts.tsv writer and the tests):
      - Singletons: `f"{ds},only"`
      - Pairs: `f"{ds_a},{ds_b}"` (alphabetical pair within `ds_order`)
      - Triples (n>=3): `f"{ds_a},{ds_b},{ds_c}"` for the unique triple,
        OR the special key `"all_three"` when n == 3 (kept for backwards
        compatibility with the 3-set counts.tsv schema).
      - n-tuples for n == 4: the special key `"all_four"`.
      - General m-tuples (3 <= m < n): comma-separated list of dataset
        names in `ds_order`.

    Sizes are the count of accessions that belong to exactly that
    subset of datasets (a partition of the union). The 3-set helper
    `_venn_region_sizes_3` and the test fixture call this function and
    map the keys back through `_legacy_3_set_region_keys` below."""
    from itertools import combinations

    n = len(ds_order)
    if n == 0:
        return {}
    union: set[str] = set()
    for d in ds_order:
        union |= sets[d]
    out: dict[str, int] = {}
    for r in range(1, n + 1):
        for combo in combinations(ds_order, r):
            members = set(combo)
            non_members = [d for d in ds_order if d not in members]
            region = set(union)
            for d in combo:
                region &= sets[d]
            for d in non_members:
                region -= sets[d]
            if r == 1:
                key = f"{combo[0]},only"
            elif r == n and n == 3:
                key = "all_three"
            elif r == n and n == 4:
                key = "all_four"
            elif r == n and n == 5:
                key = "all_five"
            else:
                key = ",".join(combo)
            out[key] = len(region)
    return out


def _venn_region_sizes_3(
    sets: dict[str, set[str]],
    ds_order: list[str],
) -> dict[str, int]:
    """Backwards-compatible 3-set wrapper around `_set_region_sizes`.

    Returns the 7 Venn region sizes using the legacy key shape:
    `{ds}_only`, `{ds_a}+{ds_b}` for every pair, `all_three`.

    Used only by the counts.tsv writer's 3-set rows and the test fixture
    `test_venn_region_sizes_3_partitions_correctly`. New 4-set callers
    should use `_set_region_sizes` directly."""
    if len(ds_order) != 3:
        raise ValueError(
            f"_venn_region_sizes_3 expects 3 datasets, got {len(ds_order)}"
        )
    generic = _set_region_sizes(sets, ds_order)
    a, b, c = ds_order
    return {
        f"{a}_only": generic[f"{a},only"],
        f"{b}_only": generic[f"{b},only"],
        f"{c}_only": generic[f"{c},only"],
        f"{a}+{b}": generic[f"{a},{b}"],
        f"{a}+{c}": generic[f"{a},{c}"],
        f"{b}+{c}": generic[f"{b},{c}"],
        "all_three": generic["all_three"],
    }


def write_combined_counts_tsv(
    tsv_path: Path,
    headlines: dict[str, DatasetHeadline],
    cellline_sets: dict[str, set[str]],
    tissue_rows: list[tuple[str, dict[str, int]]],
    accession_sets: dict[str, set[str]],
    *,
    tissue_protein_rows: (
        list[tuple[str, set[str]]]
        | list[tuple[str, dict[str, int]]]
        | None
    ) = None,
    tissue_protein_rows_per_cohort: (
        list[tuple[str, dict[str, int]]] | None
    ) = None,
    runs_per_cohort: dict[str, int] | None = None,
) -> None:
    """Auditable TSV with Panel-feeding numbers (Panel A bars, Panel B
    set-intersection regions, Panel C per-(tissue, dataset) counts,
    Panel D set-intersection regions, Panels E/F/G rows). Panel H
    (Expression Atlas overlap) was removed because it duplicates
    `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg`.
    The dataset ordering tracks the keys of `cellline_sets` so callers
    can pass 3 or 4 datasets transparently."""
    ds_order = [
        d for d in ("PXD003539", "PXD030304", "PXD004701", "PXD017199", "PXD041421")
        if d in cellline_sets
    ]
    rows: list[tuple[str, str, int, str]] = []
    # Panel A: up to 3 rows per dataset (paper + diann target-only +
    # diann unfiltered). Datasets with paper_count == 0 (PXD017199 /
    # PXD041421) write the "no paper DIA bar" row but still emit the
    # quantmsdiann target-only + unfiltered audit rows.
    for ds in ds_order:
        if ds not in headlines:
            continue
        h = headlines[ds]
        if h.paper_count > 0:
            rows.append((
                f"atlas_overlap | Panel B | {ds} | original",
                h.paper_label,
                h.paper_count,
                h.metric,
            ))
        else:
            rows.append((
                f"atlas_overlap | Panel B | {ds} | original",
                "no paper DIA headline available",
                0,
                h.metric,
            ))
        # Headline post-filter row: matches the diann bar drawn in Panel A.
        rows.append((
            f"atlas_overlap | Panel B | {ds} | quantmsdiann",
            "quantmsdiann (DIA-NN, target-only post-filter)",
            h.diann_count,
            h.metric,
        ))
        # Companion audit row: unfiltered pg_matrix / pr_matrix row count.
        unfiltered = DATASET_HEADLINES_UNFILTERED.get(ds)
        if unfiltered is not None and unfiltered != h.diann_count:
            rows.append((
                f"atlas_overlap | Panel B | {ds} | quantmsdiann",
                "quantmsdiann (DIA-NN, unfiltered pg_matrix/pr_matrix)",
                unfiltered,
                f"{h.metric} — pre conservative contaminant filter "
                "(2026-05-21 spec audit)",
            ))
    # Cell-line UpSet (former atlas_overlap Panel B) removed: PXD030304's
    # 947 lines crowded out the inter-cohort intersections. Per-cohort
    # cell-line counts retained for audit only.
    for ds in ds_order:
        rows.append((
            f"atlas_overlap | dropped (cell-line UpSet) | {ds} | count",
            "normalised cell-line names from SDRF (normalize_cell_line)",
            len(cellline_sets.get(ds, set())),
            "cell-line UpSet removed from atlas_overlap; raw count "
            "retained for audit only",
        ))
    # atlas_distribution | Panel A: per (tissue, dataset) cell-line counts.
    for tissue, by_ds in tissue_rows:
        for ds in ds_order:
            rows.append((
                f"atlas_distribution | Panel A | tissue cell lines | {tissue}",
                ds,
                by_ds.get(ds, 0),
                "cell lines mapped to this tissue (unified ProCan 28-tissue axis)",
            ))
    # atlas_overlap | Panel A: protein-accession UpSet regions (2^n − 1).
    acc_ds_order = [d for d in ds_order if d in accession_sets]
    acc_regions = _set_region_sizes(accession_sets, acc_ds_order)
    for region, n in acc_regions.items():
        rows.append((
            f"atlas_overlap | Panel A | accession intersections | {region}",
            "UniProt accessions extracted from Protein.Group (extract_accessions_diann)",
            n,
            "PXD003539/PXD017199 from pr_matrix; PXD030304/PXD004701 "
            "from cached per-tissue/per-subtype JSON",
        ))
    # Off-figure breadth-vs-depth audit — one row per cohort. The
    # scatter panel that consumed these was dropped from
    # atlas_distribution on 2026-05-29; the numbers are kept here so
    # reviewers can still inspect the breadth↔depth tradeoff.
    runs = runs_per_cohort or {}
    for ds in ds_order:
        n_cl = len(cellline_sets.get(ds, set()))
        n_acc = len(accession_sets.get(ds, set()))
        n_runs = runs.get(ds, 0)
        if n_cl == 0 or n_acc == 0:
            continue
        rows.append((
            f"atlas_distribution | off-figure | breadth-vs-depth | {ds}",
            f"cell_lines={n_cl} runs={n_runs} accessions={n_acc}",
            n_acc,
            f"{ds}: {n_cl} cell line(s) × {n_runs} MS runs → "
            f"{n_acc:,} target-only UniProt accessions",
        ))
    # atlas_distribution | Panel B: unique target-only accessions per
    # tissue (union across cohorts) — the headline number on the bar.
    for tissue, value in (tissue_protein_rows or []):
        union_size = len(value) if isinstance(value, set) else sum(value.values())
        rows.append((
            f"atlas_distribution | Panel B | tissue proteins | {tissue}",
            "union of UniProt accessions across all contributing cohorts",
            union_size,
            "unique target-only proteins detected in this tissue across "
            "the atlas (post-contaminant filter)",
        ))
    # Per-cohort breakdown of the per-tissue protein counts (kept
    # off-figure in this audit row so reviewers can decompose the
    # union bars by source cohort).
    for tissue, by_ds in (tissue_protein_rows_per_cohort or []):
        for ds in ds_order:
            rows.append((
                f"atlas_distribution | Panel B | tissue proteins (per-cohort breakdown) | {tissue}",
                ds,
                by_ds.get(ds, 0),
                "per-cohort target-only protein count for this tissue "
                "(sum across cohorts double-counts proteins seen in "
                "multiple cohorts — figure shows union)",
            ))
    # Panel G: detection-count histogram.
    g_buckets = detection_count_histogram(accession_sets)
    g_total = sum(g_buckets.values())
    n_datasets = len(accession_sets)
    for k, n in sorted(g_buckets.items()):
        pct = (100.0 * n / g_total) if g_total > 0 else 0.0
        rows.append((
            f"atlas_overlap | Panel C | detections-in-n-datasets | n={k}",
            f"UniProt accessions across the {n_datasets} quantmsdiann analyses",
            n,
            f"{pct:.2f}% of pan-cohort union ({g_total:,})",
        ))
    # Panel H removed: the Expression Atlas (E-PROT-73) gene-set
    # overlap was redundant with the per-cohort PXD003539 figure at
    # `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg`.
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
    # PXD017199 is atlas-only — no per-cohort figure; inputs ship directly.
    (PXD017199_SDRF, "stage data/PXD017199/PXD017199.sdrf.tsv (Tognetti 2021)"),
    (PXD017199_PR_MATRIX, "stage data/PXD017199/diann_report.pr_matrix.tsv (Tognetti 2021)"),
    # PXD041421 is atlas-only per 2026-05-21 spec §2; inputs ship directly.
    (PXD041421_SDRF, "stage data/PXD041421/PXD041421.sdrf.tsv (Wang 2023)"),
    (PXD041421_PR_MATRIX, "stage data/PXD041421/diann_report.pr_matrix.tsv (Wang 2023)"),
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

    print("Refreshing DATASET_HEADLINES diann_count from pg_matrix files "
          "(target-only contaminant filter)...")
    refresh_dataset_headlines()
    for ds, h in DATASET_HEADLINES.items():
        unf = DATASET_HEADLINES_UNFILTERED.get(ds)
        print(f"  {ds}: target_only={h.diann_count:,}  unfiltered={unf:,}  "
              f"(delta {unf - h.diann_count:,})")

    print("Loading cell-line sets from SDRFs...")
    cellline_sets: dict[str, set[str]] = {
        "PXD003539": cell_lines_from_sdrf(PXD003539_SDRF),
        "PXD030304": cell_lines_from_sdrf(PXD030304_SDRF),
        "PXD004701": cell_lines_from_sdrf(PXD004701_SDRF),
        "PXD017199": cell_lines_from_sdrf(PXD017199_SDRF),
        "PXD041421": cell_lines_from_sdrf(PXD041421_SDRF),
    }
    for d, s in cellline_sets.items():
        print(f"  {d}: {len(s):,} cell lines")

    print("Building per-dataset cell-line -> tissue maps...")
    cl_tissue_3539 = cell_line_tissue_pxd003539(PXD003539_SDRF)
    cl_tissue_30304 = cell_line_tissue_pxd030304(PXD030304_SDRF,
                                                 PXD030304_TISSUE_MAPPING)
    cl_tissue_4701 = cell_line_tissue_pxd004701(PXD004701_SDRF)
    cl_tissue_17199 = cell_line_tissue_pxd017199(PXD017199_SDRF)
    cl_tissue_41421 = cell_line_tissue_pxd041421(PXD041421_SDRF)
    tissue_rows = combined_tissue_table({
        "PXD003539": cl_tissue_3539,
        "PXD030304": cl_tissue_30304,
        "PXD004701": cl_tissue_4701,
        "PXD017199": cl_tissue_17199,
        "PXD041421": cl_tissue_41421,
    })
    print(f"  unified-axis tissues: {len(tissue_rows)}")
    print("  Top-5 tissues by combined cell-line count:")
    for tissue, by_ds in tissue_rows[:5]:
        total = sum(by_ds.values())
        seg = ", ".join(
            f"{ds}:{by_ds.get(ds, 0)}"
            for ds in ("PXD003539", "PXD030304", "PXD004701", "PXD017199", "PXD041421")
        )
        print(f"    {tissue:<32s} total={total:>4d}  ({seg})")

    print("Loading per-group accession sets (for Panels E and F)...")
    per_cl_3539 = pxd003539_accessions_per_cell_line(
        PXD003539_PR_MATRIX, PXD003539_SDRF,
    )
    per_tissue_30304 = pxd030304_accessions_per_tissue(PXD030304_PROTEIN_JSON)
    per_subtype_4701 = pxd004701_accessions_per_subtype(PXD004701_PROTEIN_JSON)
    per_cl_17199 = pxd017199_accessions_per_cell_line(
        PXD017199_PR_MATRIX, PXD017199_SDRF,
    )
    per_cl_41421 = pxd041421_accessions_per_cell_line(
        PXD041421_PR_MATRIX, PXD041421_SDRF,
    )
    print(f"  PXD003539 per-cell-line groups: {len(per_cl_3539)}")
    print(f"  PXD030304 per-tissue groups: {len(per_tissue_30304)}")
    print(f"  PXD004701 per-subtype groups: {len(per_subtype_4701)}")
    print(f"  PXD017199 per-cell-line groups: {len(per_cl_17199)}")
    print(f"  PXD041421 per-cell-line groups: {len(per_cl_41421)}")

    print("Loading protein-accession sets...")
    # Reuse per-group sets when possible to avoid re-reading the pr_matrix.
    acc_3539: set[str] = set()
    for s in per_cl_3539.values():
        acc_3539.update(s)
    if not acc_3539:
        acc_3539 = pxd003539_protein_accessions(PXD003539_PR_MATRIX)
    acc_30304 = pxd030304_protein_accessions(PXD030304_PROTEIN_JSON)
    acc_4701 = pxd004701_protein_accessions(PXD004701_PROTEIN_JSON)
    acc_17199: set[str] = set()
    for s in per_cl_17199.values():
        acc_17199.update(s)
    if not acc_17199:
        acc_17199 = pxd017199_protein_accessions(PXD017199_PR_MATRIX)
    acc_41421: set[str] = set()
    for s in per_cl_41421.values():
        acc_41421.update(s)
    if not acc_41421:
        acc_41421 = pxd041421_protein_accessions(PXD041421_PR_MATRIX)
    accession_sets = {
        "PXD003539": acc_3539,
        "PXD030304": acc_30304,
        "PXD004701": acc_4701,
        "PXD017199": acc_17199,
        "PXD041421": acc_41421,
    }
    for d, s in accession_sets.items():
        print(f"  {d}: {len(s):,} accessions")

    print("Building per-tissue UNIQUE-accession rows (atlas_distribution Panel B)...")
    tissue_order = [t for t, _ in tissue_rows]
    # Union across cohorts per tissue — the new figure shape.
    tissue_protein_rows = per_tissue_union_accessions(
        tissue_order,
        cl_tissue_3539,
        per_cl_3539,
        per_tissue_30304,
        acc_4701,
        cl_tissue_pxd017199=cl_tissue_17199,
        accessions_per_cell_line_pxd017199=per_cl_17199,
        cl_tissue_pxd041421=cl_tissue_41421,
        accessions_per_cell_line_pxd041421=per_cl_41421,
    )
    # Also compute the per-cohort breakdown for the audit TSV (so the
    # cohort decomposition is preserved off-figure).
    tissue_protein_rows_per_cohort = per_tissue_protein_counts(
        tissue_order,
        cl_tissue_3539,
        per_cl_3539,
        per_tissue_30304,
        acc_4701,
        cl_tissue_pxd017199=cl_tissue_17199,
        accessions_per_cell_line_pxd017199=per_cl_17199,
        cl_tissue_pxd041421=cl_tissue_41421,
        accessions_per_cell_line_pxd041421=per_cl_41421,
    )

    # Panel E inputs (breadth-vs-depth scatter): one point per cohort
    # where dot size scales with the number of MS runs. Run counts come
    # from SDRF row counts where available — quick and matches what
    # the per-cohort design docs report.
    print("Computing per-cohort MS-run counts (Panel E)...")
    runs_per_cohort = _compute_runs_per_cohort()
    for ds in ("PXD003539", "PXD030304", "PXD004701",
               "PXD017199", "PXD041421"):
        print(f"  {ds}: {runs_per_cohort.get(ds, 0):,} runs")

    print("Rendering main pan-cohort figure (headline + per-tissue panels)...")
    main_svg = FIGURES_DIR / "atlas_main.svg"
    render_atlas_main(
        DATASET_HEADLINES, tissue_rows, main_svg,
        tissue_protein_rows=tissue_protein_rows,
    )
    print(f"  saved: {main_svg}")

    print("Writing combined counts.tsv...")
    data_dir = FIGURES_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    tsv = data_dir / "combined_counts.tsv"
    write_combined_counts_tsv(
        tsv, DATASET_HEADLINES, cellline_sets, tissue_rows, accession_sets,
        tissue_protein_rows=tissue_protein_rows,
        tissue_protein_rows_per_cohort=tissue_protein_rows_per_cohort,
        runs_per_cohort=runs_per_cohort,
    )
    print(f"  saved: {tsv}")

    ds_order = ["PXD003539", "PXD030304", "PXD004701", "PXD017199", "PXD041421"]
    print("\nPanel B (cell-line set intersections) region sizes:")
    for region, n in _set_region_sizes(cellline_sets, ds_order).items():
        print(f"  {region:<48s} {n:>6,}")
    print("\nPanel D (accession set intersections) region sizes:")
    for region, n in _set_region_sizes(accession_sets, ds_order).items():
        print(f"  {region:<48s} {n:>6,}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
