"""PXD004701 reanalysis comparison figure (Sun et al. 2023 vs quantmsdiann).

The published 2023 Sun et al. PCT-SWATH analysis of 76 breast cancer cell
lines (Mol Cell Proteomics, doi:10.1016/j.mcpro.2023.100602, PMC10392136)
reports 6,091 SwissProt proteins consistently identified across all samples
after applying a proteotypic-peptide + Global.Q.Value <= 0.01 filter, then
dropping proteins with >90 % missing rate. The PMC supplement is behind a
proof-of-work CAPTCHA so we cannot retrieve the per-cell TNBC / non-TNBC
mapping; the `BC_SUBTYPES` dict below is our reconstruction from the
breast-cancer literature (Heiser 2012, Neve 2006, Lehmann 2011, Cellosaurus).

Outputs (paper-ready, no titles/footers):
- main_comparison.{pdf,png,svg}: 2 conditions x 3 metrics
- supp_proteins_per_subtype.{pdf,png,svg}: per-subtype protein counts
- supp_missing_values_per_run.{pdf,png,svg}: per-run completeness
- counts.tsv: auditable totals + per-subtype numbers
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis.figure_original_vs_quantmsdiann import (
    download_if_missing,
)
from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
    PG_METADATA_COLS,
    PR_METADATA_COLS,
    parse_diann_summary_log,
    per_run_completeness_quantmsdiann,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "PXD004701"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "PXD004701"

PRIDE_QUANT_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines/PXD004701"
)

DIANN_SUMMARY_LOG_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diannsummary.log"
DIANN_PG_MATRIX_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diann_report.pg_matrix.tsv"
DIANN_PR_MATRIX_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diann_report.pr_matrix.tsv"
DIANN_PARQUET_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diann_report.parquet"
QUANTMS_SDRF_URL = f"{PRIDE_QUANT_BASE}/sdrf/PXD004701.sdrf.tsv"

# Headline constants from Sun et al. 2023 (MCP).
SUN_PROTEINS = 6091           # global 1% Global.Q.Value + <=90% missing
SUN_PROTEINS_RAW = 8952       # pre-consistency-filter total (paper Methods)
SUN_PEPTIDES = 90762          # proteotypic peptides under same filter
SUN_LIBRARY_PRECURSORS = 194899   # pan-human CAL library precursors (NOT IDs)
SUN_LIBRARY_PROTEINS = 10323      # pan-human CAL library proteins
SUN_TNBC = 39                 # paper-reported subtype split
SUN_NON_TNBC = 37
SUN_CELL_LINES = 76

# Hardcoded BC subtype classification for the 76 PXD004701 SDRF cell lines.
# Source: standard breast-cancer cell-line literature — Heiser et al. 2012
# (PNAS, doi:10.1073/pnas.1018854108), Neve et al. 2006 (Cancer Cell,
# doi:10.1016/j.ccr.2006.10.008), Lehmann et al. 2011 (JCI,
# doi:10.1172/JCI45014), and Cellosaurus receptor-status annotations.
# Borderline calls documented in the spec:
# - mdamb453: kept TNBC despite AR+ behaviour (canonical assignment in
#   Lehmann 2011's LAR subtype, which is still ER-/PR-/HER2-).
# - skbr7: kept non-TNBC; HER2 weak-positive per Cellosaurus, no canonical
#   TNBC assignment in the cited literature.
# - hcc2185: classified TNBC per Lehmann 2011 (metaplastic carcinoma).
# - 184a1/184b5/hbl100/mcf10a/mcf12a: non-tumorigenic mammary lines kept
#   as `normal-like` because they are not breast cancers; Sun et al. fold
#   them into one of the two cancer subtypes (which 3-line gap explains
#   most of the 39 vs 36 TNBC delta).
BC_SUBTYPES: dict[str, str] = {
    # normal-like (5)
    "184a1": "normal-like",
    "184b5": "normal-like",
    "hbl100": "normal-like",
    "mcf10a": "normal-like",
    "mcf12a": "normal-like",

    # TNBC (ER-/PR-/HER2-): 34 lines
    "bt20": "TNBC",
    "bt549": "TNBC",
    "cal51": "TNBC",
    "cal120": "TNBC",
    "cal148": "TNBC",
    "du4475": "TNBC",
    "evsat": "TNBC",
    "hcc1143": "TNBC",
    "hcc1187": "TNBC",
    "hcc1395": "TNBC",
    "hcc1599": "TNBC",
    "hcc1806": "TNBC",
    "hcc1937": "TNBC",
    "hcc2185": "TNBC",
    "hcc3153": "TNBC",
    "hcc38": "TNBC",
    "hcc70": "TNBC",
    "hdqp1": "TNBC",
    "hs578t": "TNBC",
    "mb157": "TNBC",
    "mdamb157": "TNBC",
    "mdamb231": "TNBC",
    "mdamb436": "TNBC",
    "mdamb453": "TNBC",
    "mdamb468": "TNBC",
    "mfm223": "TNBC",
    "mx1": "TNBC",
    "ocubm": "TNBC",
    "sum102": "TNBC",
    "sum149": "TNBC",
    "sum159": "TNBC",
    "sum190": "TNBC",
    "sum229": "TNBC",
    "macls2": "TNBC",

    # non-TNBC (receptor-positive or HER2+): 35 lines
    "au565": "non-TNBC",
    "bt474": "non-TNBC",
    "bt483": "non-TNBC",
    "cama1": "non-TNBC",
    "efm19": "non-TNBC",
    "efm192a": "non-TNBC",
    "hcc1419": "non-TNBC",
    "hcc1428": "non-TNBC",
    "hcc1569": "non-TNBC",
    "hcc1954": "non-TNBC",
    "hcc202": "non-TNBC",
    "hcc2218": "non-TNBC",
    "hcc2688": "non-TNBC",
    "jimt1": "non-TNBC",
    "kpl1": "non-TNBC",
    "mcf7": "non-TNBC",
    "mdamb134vi": "non-TNBC",
    "mdamb175vii": "non-TNBC",
    "mdamb330": "non-TNBC",
    "mdamb361": "non-TNBC",
    "mdamb415": "non-TNBC",
    "skbr3": "non-TNBC",
    "skbr5": "non-TNBC",
    "skbr7": "non-TNBC",
    "sum185": "non-TNBC",
    "sum225": "non-TNBC",
    "sum44": "non-TNBC",
    "sum52": "non-TNBC",
    "sw527": "non-TNBC",
    "t47d": "non-TNBC",
    "uacc3199": "non-TNBC",
    "uacc893": "non-TNBC",
    "zr751": "non-TNBC",
    "zr7530": "non-TNBC",
    "zr75b": "non-TNBC",
    "600mpe": "non-TNBC",
    "ly2": "non-TNBC",
}


@dataclass(frozen=True)
class Counts:
    sun_proteins: int                       # 6,091 (paper consistency-filtered)
    sun_proteins_raw: int                   # 8,952 (paper pre-filter)
    sun_peptides: int                       # 90,762 (paper)
    sun_tnbc: int                           # 39 (paper)
    sun_non_tnbc: int                       # 37 (paper)
    quantmsdiann_proteins_strict: int       # 7,746 (diannsummary.log)
    quantmsdiann_proteins_consistent: int   # post-consistency-filter union
    quantmsdiann_peptides: int              # unique Stripped.Sequence
    quantmsdiann_precursors: int            # 100,499 (diannsummary.log)


# ---------------------------------------------------------------------------
# SDRF parsing
# ---------------------------------------------------------------------------


def parse_sdrf_data_file_to_cell_line(sdrf_path: Path) -> dict[str, str]:
    """Parse `comment[data file]` -> `characteristics[cell line]`, rewriting
    `.wiff`/`.WIFF` -> `.mzML` so DIA-NN matrix column names match.

    Duplicated from the PXD030304 helper (rather than imported) so PXD004701
    stays self-contained — PXD030304 may evolve its column expectations
    independently."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    needed = ["characteristics[cell line]", "comment[data file]"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"SDRF missing required columns: {missing}")
    out: dict[str, str] = {}
    for cell, data_file in zip(
        df["characteristics[cell line]"], df["comment[data file]"],
    ):
        if not isinstance(data_file, str) or not data_file:
            continue
        mzml = re.sub(r"\.wiff$", ".mzML", data_file, flags=re.IGNORECASE)
        out[mzml] = cell.strip() if isinstance(cell, str) else cell
    return out


# ---------------------------------------------------------------------------
# Sun-style consistency filter on quantmsdiann parquet (two-stage)
# ---------------------------------------------------------------------------


def _compute_or_load_diann_subtype_consistency_filter(
    cache_path: Path,
    parquet_source: str | Path,
    sdrf_path: Path,
    subtype_dict: dict[str, str],
    *,
    qvalue_cutoff: float = 0.01,
    min_detection_fraction: float = 0.10,
) -> dict[str, set[str]]:
    """Side-cache wrapper around
    `proteins_per_subtype_quantmsdiann_consistency_filter`. Streaming 33 GB
    over HTTP takes ~10-20 minutes, so we persist the per-subtype result as
    a small JSON. Delete the JSON to force a fresh stream."""
    import json
    if cache_path.exists() and cache_path.stat().st_size > 0:
        with open(cache_path, encoding="utf-8") as fh:
            payload = json.load(fh)
        return {s: set(vs) for s, vs in payload.items()}
    result = proteins_per_subtype_quantmsdiann_consistency_filter(
        parquet_source, sdrf_path, subtype_dict,
        qvalue_cutoff=qvalue_cutoff,
        min_detection_fraction=min_detection_fraction,
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump({s: sorted(v) for s, v in result.items()}, fh)
    return result


def proteins_per_subtype_quantmsdiann_consistency_filter(
    parquet_source: str | Path,
    sdrf_path: Path,
    subtype_dict: dict[str, str],
    *,
    qvalue_cutoff: float = 0.01,
    min_detection_fraction: float = 0.10,
    batch_size: int = 1_000_000,
) -> dict[str, set[str]]:
    """Sun-style two-stage filter applied to quantmsdiann's long-format
    report.

    Stage 1 (FDR): keep precursor rows where `Proteotypic == 1` AND
    `Global.Q.Value <= qvalue_cutoff` AND whose `Run` maps to one of the
    SDRF cell lines (those rows form the "mapped run space"; runs outside
    this set are ignored entirely).

    Stage 2 (consistency): for each Protein.Group surviving stage 1,
    compute its mapped-run detection fraction (# distinct mapped runs in
    which it appears / total mapped runs). Drop any Protein.Group with
    detection fraction < `min_detection_fraction` (Sun et al. use
    >90 % missing == <10 % detection).

    Stage 3 (subtype aggregation): for the surviving (run, Protein.Group)
    pairs, group by `subtype_dict[cell_line]` and emit per-subtype unions.
    Cell lines mapped to `'unknown'` are excluded from the subtype aggregation
    (they still count toward the consistency-filter denominator).

    Streams the parquet via fsspec's HTTPFileSystem with column projection
    on (`Run`, `Protein.Group`, `Global.Q.Value`, `Proteotypic`) when given
    an HTTPS URL; only those columns' chunks transit the wire."""
    import pyarrow.parquet as pq

    sdrf_run_to_cell = parse_sdrf_data_file_to_cell_line(sdrf_path)
    # SDRF data files end `.mzML` after the .wiff rewrite; the parquet's
    # `Run` column carries the bare basename. Strip the extension.
    sdrf_no_ext: dict[str, str] = {}
    for k, v in sdrf_run_to_cell.items():
        stem = re.sub(r"\.(mzML|wiff)$", "", k, flags=re.IGNORECASE)
        sdrf_no_ext[stem] = v

    mapped_runs = set(sdrf_no_ext.keys())
    total_mapped_runs = len(mapped_runs)
    if total_mapped_runs == 0:
        return {}

    # Stage 1: collect (run, pg) pairs passing FDR + Proteotypic + run is
    # mapped. We need per-pg run-set sizes for stage 2 and per-run pg sets
    # for stage 3, so we keep run-set-per-pg as the canonical store.
    pg_to_run_set: dict[str, set[str]] = {}
    cols = ["Run", "Protein.Group", "Global.Q.Value", "Proteotypic"]
    source = str(parquet_source)
    if source.startswith(("http://", "https://")):
        import fsspec
        fs = fsspec.filesystem("https")
        opener = lambda: fs.open(source, "rb")
    else:
        opener = lambda: open(source, "rb")

    with opener() as fh:
        pf = pq.ParquetFile(fh)
        for batch in pf.iter_batches(batch_size=batch_size, columns=cols):
            runs = batch.column("Run").to_pylist()
            pgs = batch.column("Protein.Group").to_pylist()
            gqv = batch.column("Global.Q.Value").to_pylist()
            prot = batch.column("Proteotypic").to_pylist()
            for r, pg, g, p in zip(runs, pgs, gqv, prot):
                if p != 1 or g is None or g > qvalue_cutoff:
                    continue
                if r not in sdrf_no_ext:
                    continue
                pg_to_run_set.setdefault(pg, set()).add(r)

    # Stage 2: drop pg's with detection fraction below the threshold.
    threshold = min_detection_fraction * total_mapped_runs
    consistent = {
        pg: runs for pg, runs in pg_to_run_set.items()
        if len(runs) >= threshold
    }

    # Stage 3: aggregate to per-subtype unions.
    out: dict[str, set[str]] = {}
    for pg, runs in consistent.items():
        for r in runs:
            cell = sdrf_no_ext.get(r)
            if cell is None:
                continue
            subtype = subtype_dict.get(cell)
            if subtype is None or subtype == "unknown":
                continue
            out.setdefault(subtype, set()).add(pg)
    return out


# ---------------------------------------------------------------------------
# Peptide count from pr_matrix
# ---------------------------------------------------------------------------


def unique_peptides_quantified(pr_matrix_path: Path) -> int:
    """Unique `Stripped.Sequence` count in `pr_matrix.tsv` among rows with at
    least one non-NA per-run value. Mirrors the PXD003539/PXD030304 peptide
    definition.

    pr_matrix.tsv for PXD004701 is ~2 GB; we read in chunks and union the
    set of stripped sequences across chunks."""
    seqs: set[str] = set()
    sample_cols: list[str] | None = None
    for chunk in pd.read_csv(
        pr_matrix_path, sep="\t", dtype=str, chunksize=50_000,
    ):
        if sample_cols is None:
            missing = [c for c in PR_METADATA_COLS if c not in chunk.columns]
            if missing:
                raise ValueError(
                    f"pr_matrix missing metadata columns: {missing}"
                )
            sample_cols = [c for c in chunk.columns if c not in PR_METADATA_COLS]
        any_quant = chunk[sample_cols].notna().any(axis=1)
        for seq in chunk.loc[any_quant, "Stripped.Sequence"]:
            if isinstance(seq, str) and seq:
                seqs.add(seq)
    return len(seqs)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def render_main_figure(
    counts: Counts,
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """Grouped bar chart: Sun et al. 2023 vs quantmsdiann across 3 metrics.

    Metric 1: Protein groups (Sun-style consistency filter) — paper 6,091 vs
    quantmsdiann post-filter union.
    Metric 2: Proteotypic peptides — paper 90,762 vs quantmsdiann unique
    Stripped.Sequence in pr_matrix.tsv.
    Metric 3: Protein groups (strict 1% FDR, no consistency) — paper 8,952
    pre-filter vs quantmsdiann 7,746 from diannsummary.log.

    Paper-ready: no title, no footer."""
    metrics = [
        "Protein groups\n(consistency filter,\n$\\leq$90% missing)",
        "Proteotypic peptides",
        "Protein groups\n(1% FDR, no\nconsistency filter)",
    ]
    sun_vals = [counts.sun_proteins, counts.sun_peptides, counts.sun_proteins_raw]
    diann_vals = [
        counts.quantmsdiann_proteins_consistent,
        counts.quantmsdiann_peptides,
        counts.quantmsdiann_proteins_strict,
    ]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    bar_width = 0.35
    x = list(range(len(metrics)))
    bars_s = ax.bar([xi - bar_width / 2 for xi in x], sun_vals,
                    width=bar_width, color="#9e9e9e",
                    label="Sun et al. 2023")
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals,
                    width=bar_width, color="#1f77b4",
                    label="quantmsdiann (DIA-NN)")
    for bars, vals in ((bars_s, sun_vals), (bars_d, diann_vals)):
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{v:,}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=9)
    ax.set_ylabel("Count")
    ymax = max(max(sun_vals), max(diann_vals)) * 1.18
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", frameon=False)

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def render_proteins_per_subtype(
    diann_per_subtype: dict[str, set[str]],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
    *,
    sun_reference: int = SUN_PROTEINS,
) -> None:
    """Bar chart: quantmsdiann per-subtype protein-group union (TNBC,
    non-TNBC, normal-like). Horizontal reference line at Sun et al.'s 6,091
    global consistency-filtered total — they do not publish per-subtype
    protein counts so a side-by-side per-subtype comparison is not possible.
    Paper-ready: no title, no footer."""
    subtype_order = ["TNBC", "non-TNBC", "normal-like"]
    labels = [s for s in subtype_order if s in diann_per_subtype]
    vals = [len(diann_per_subtype[s]) for s in labels]

    fig, ax = plt.subplots(figsize=(7, 5))
    x = list(range(len(labels)))
    colors = {"TNBC": "#d62728", "non-TNBC": "#1f77b4", "normal-like": "#2ca02c"}
    bars = ax.bar(x, vals, color=[colors.get(s, "#777777") for s in labels],
                  width=0.55)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{v:,}", ha="center", va="bottom", fontsize=10)
    ax.axhline(sun_reference, color="#9e9e9e", linestyle="--", linewidth=1,
               label=f"Sun et al. 2023 global ({sun_reference:,})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Distinct protein groups\n(consistency filter)")
    ymax = max(max(vals, default=0), sun_reference) * 1.15
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", frameon=False)

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def render_per_run_completeness(
    diann_per_run: dict[str, float],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
) -> None:
    """Line plot: per-run fraction of detected protein groups in
    quantmsdiann's pg_matrix. Sun et al. publish no per-run completeness
    data so this is single-condition. Paper-ready: no title, no footer."""
    diann_vals = sorted(diann_per_run.values())
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(range(len(diann_vals)), diann_vals,
            color="#1f77b4", linewidth=0.9,
            label=f"quantmsdiann (DIA-NN) (n={len(diann_vals)} runs)")
    ax.set_xlabel("Run rank (sorted ascending)")
    ax.set_ylabel("Fraction of protein groups\ndetected per run")
    ax.set_ylim(0, 1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False)

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def write_counts_tsv(
    counts: Counts,
    tsv_path: Path,
    *,
    diann_per_subtype: dict[str, set[str]] | None = None,
) -> None:
    """Auditable counts table with inline methodology notes."""
    rows = [
        ("Protein groups (consistency filter)", "Sun et al. 2023 (paper headline)",
         counts.sun_proteins,
         "6,091 SwissProt proteins; proteotypic, Global.Q.Value<=0.01, "
         ">=10% detection across samples (paper, Methods)"),
        ("Protein groups (consistency filter)", "quantmsdiann (DIA-NN)",
         counts.quantmsdiann_proteins_consistent,
         "post-filter union across all subtypes; same two-stage filter "
         "applied to diann_report.parquet"),
        ("Protein groups (1% FDR, no consistency)", "Sun et al. 2023 (paper pre-filter)",
         counts.sun_proteins_raw,
         "8,952 proteins identified pre-consistency-filter (paper, Methods)"),
        ("Protein groups (1% FDR, no consistency)", "quantmsdiann (DIA-NN)",
         counts.quantmsdiann_proteins_strict,
         "from diannsummary.log (Protein groups with global q-value <= 0.01)"),
        ("Proteotypic peptides", "Sun et al. 2023 (paper headline)",
         counts.sun_peptides,
         "90,762 proteotypic peptides under the same consistency filter "
         "(paper, Methods)"),
        ("Proteotypic peptides", "quantmsdiann (DIA-NN)",
         counts.quantmsdiann_peptides,
         "unique Stripped.Sequence in pr_matrix.tsv among rows with "
         ">=1 non-NA sample"),
        ("Precursors", "quantmsdiann (DIA-NN, 1% FDR)",
         counts.quantmsdiann_precursors,
         "from diannsummary.log (Target precursors at 1% global q-value)"),
        ("Spectral library precursors", "Sun et al. 2023 (pan-human CAL)",
         SUN_LIBRARY_PRECURSORS,
         "194,899 library precursors; NOT identified precursors (PXD009597)"),
        ("Spectral library proteins", "Sun et al. 2023 (pan-human CAL)",
         SUN_LIBRARY_PROTEINS,
         "10,323 SwissProt library proteins (PXD009597)"),
        ("Cell lines: TNBC", "Sun et al. 2023 (paper-reported split)",
         counts.sun_tnbc,
         "39 TNBC lines in paper's split of the 76 cell lines"),
        ("Cell lines: non-TNBC", "Sun et al. 2023 (paper-reported split)",
         counts.sun_non_tnbc,
         "37 non-TNBC lines in paper's split of the 76 cell lines"),
    ]
    # BC_SUBTYPES partition rows
    subtype_counts: dict[str, int] = {}
    for cl, st in BC_SUBTYPES.items():
        subtype_counts[st] = subtype_counts.get(st, 0) + 1
    for st in ("TNBC", "non-TNBC", "normal-like", "unknown"):
        if st in subtype_counts:
            rows.append((
                f"Cell lines: {st} (this work mapping)",
                "quantmsdiann (BC_SUBTYPES from Heiser/Neve/Lehmann + Cellosaurus)",
                subtype_counts[st],
                "hardcoded classification; PMC supplement inaccessible "
                "behind a CAPTCHA so per-cell paper assignment not retrievable",
            ))
    # Per-subtype protein counts (quantmsdiann only; paper has no per-subtype
    # number to compare against).
    if diann_per_subtype is not None:
        note = (
            "per-subtype union of Protein.Group from diann_report.parquet "
            "filtered to Proteotypic==1 AND Global.Q.Value<=0.01 AND "
            "consistency filter (>=10% detection across mapped runs)"
        )
        for st in ("TNBC", "non-TNBC", "normal-like"):
            if st in diann_per_subtype:
                rows.append((
                    f"Per-subtype proteins | {st}",
                    "quantmsdiann (DIA-NN)",
                    len(diann_per_subtype[st]),
                    note,
                ))
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", encoding="utf-8") as fh:
        fh.write("metric\tsource\tcount\tnote\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:  # pragma: no cover
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    log_path = download_if_missing(DIANN_SUMMARY_LOG_URL,
                                   DATA_DIR / "diannsummary.log")
    pg_path = download_if_missing(DIANN_PG_MATRIX_URL,
                                  DATA_DIR / "diann_report.pg_matrix.tsv")
    pr_path = download_if_missing(DIANN_PR_MATRIX_URL,
                                  DATA_DIR / "diann_report.pr_matrix.tsv")
    sdrf_path = download_if_missing(QUANTMS_SDRF_URL,
                                    DATA_DIR / "PXD004701.sdrf.tsv")

    print("Parsing DIA-NN summary log...")
    pg, prec = parse_diann_summary_log(log_path)
    print(f"  protein groups: {pg:,}  precursors: {prec:,}")

    print("Counting unique proteotypic peptides in pr_matrix.tsv...")
    pep = unique_peptides_quantified(pr_path)
    print(f"  unique Stripped.Sequence: {pep:,}")

    print("Computing per-subtype consistency-filtered protein sets "
          "(streaming parquet)...")
    # Streams 33 GB of diann_report.parquet over HTTP with column projection;
    # cached to a small JSON.
    diann_per_subtype = _compute_or_load_diann_subtype_consistency_filter(
        DATA_DIR / "diann_per_subtype_consistency_filter.json",
        DIANN_PARQUET_URL,
        sdrf_path,
        BC_SUBTYPES,
    )
    diann_proteins_consistent = set()
    for s, pgs in diann_per_subtype.items():
        diann_proteins_consistent.update(pgs)
        print(f"  {s:<14s} {len(pgs):>6,}")
    print(f"  union across subtypes: {len(diann_proteins_consistent):,}")

    counts = Counts(
        sun_proteins=SUN_PROTEINS,
        sun_proteins_raw=SUN_PROTEINS_RAW,
        sun_peptides=SUN_PEPTIDES,
        sun_tnbc=SUN_TNBC,
        sun_non_tnbc=SUN_NON_TNBC,
        quantmsdiann_proteins_strict=pg,
        quantmsdiann_proteins_consistent=len(diann_proteins_consistent),
        quantmsdiann_peptides=pep,
        quantmsdiann_precursors=prec,
    )

    print("Rendering main figure...")
    render_main_figure(
        counts,
        FIGURES_DIR / "main_comparison.pdf",
        FIGURES_DIR / "main_comparison.png",
        FIGURES_DIR / "main_comparison.svg",
    )

    print("Rendering per-subtype supp figure...")
    render_proteins_per_subtype(
        diann_per_subtype,
        FIGURES_DIR / "supp_proteins_per_subtype.pdf",
        FIGURES_DIR / "supp_proteins_per_subtype.png",
        FIGURES_DIR / "supp_proteins_per_subtype.svg",
    )

    print("Computing per-run completeness (quantmsdiann)...")
    diann_per_run = per_run_completeness_quantmsdiann(pg_path)
    print(f"  {len(diann_per_run)} runs")

    print("Rendering per-run completeness supp figure...")
    render_per_run_completeness(
        diann_per_run,
        FIGURES_DIR / "supp_missing_values_per_run.pdf",
        FIGURES_DIR / "supp_missing_values_per_run.png",
        FIGURES_DIR / "supp_missing_values_per_run.svg",
    )

    print("Writing auditable counts TSV...")
    write_counts_tsv(
        counts, FIGURES_DIR / "counts.tsv",
        diann_per_subtype=diann_per_subtype,
    )

    # Cross-checks (non-gating)
    if pg != 7746:
        print(f"WARN: quantmsdiann protein groups {pg} != expected 7,746",
              file=sys.stderr)
    if prec != 100499:
        print(f"WARN: quantmsdiann precursors {prec} != expected 100,499",
              file=sys.stderr)

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
