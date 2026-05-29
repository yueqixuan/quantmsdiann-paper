"""PXD030304 reanalysis comparison figure (ProCan-DepMapSanger vs quantmsdiann).

The published 2022 ProCan-DepMapSanger analysis (Gonçalves et al., Cancer Cell,
doi:10.1016/j.ccell.2022.06.010) provides the original 949-cell-line proteomic
map. The deposited DIA-NN long-format report on PRIDE is 237 GB; we instead
use the authors' figshare deposit (doi:10.6084/m9.figshare.19345397) which has
small post-processed per-sample protein matrices.

Outputs (paper-ready, no titles/footers):
- main_comparison.{pdf,png,svg}: 2 conditions x 2 metrics
- supp_proteins_per_tissue.{pdf,png,svg}: per-tissue protein counts
- supp_missing_values_per_run.{pdf,png,svg}: per-run completeness
- supp_venn_protein_accessions.{pdf,png,svg}: protein-set Venn
- counts.tsv: auditable totals + per-tissue numbers
"""
from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis.contaminant_filter import count_target_protein_groups
from analysis.figure_original_vs_quantmsdiann import (
    download_if_missing,
    unique_peptides_per_protein_diann,
    proteins_with_min_peptides,
    SUMMARY_LOG_PROTEIN_LINE_RE,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "PXD030304"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "PXD030304"

PRIDE_QUANT_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines/PXD030304"
)

DIANN_SUMMARY_LOG_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diannsummary.log"
DIANN_PG_MATRIX_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diann_report.pg_matrix.tsv"
DIANN_PR_MATRIX_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diann_report.pr_matrix.tsv"
DIANN_PARQUET_URL = f"{PRIDE_QUANT_BASE}/quant_tables/diann_report.parquet"
QUANTMS_SDRF_URL = f"{PRIDE_QUANT_BASE}/sdrf/PXD030304.sdrf.tsv"

# Figshare deposit 19345397 (Gonçalves et al. 2022 supporting data). The
# figshare ndownloader URLs are stable and content-addressed by file id.
# We use the per-replicate matrix (519 MB) — not the averaged matrix — so
# the per-tissue comparison applies identical per-MS-run union semantics on
# both sides; the averaged matrix collapses 6 replicates per cell line and
# would be more permissive than quantmsdiann's per-cell-level filter.
FIGSHARE_FILES = {
    "protein_matrix_8498_replicates.txt": 34411235,
    "peptide_counts_per_protein_per_sample.txt": 34411148,
    "mapping_file_averaged.txt": 34411133,
    "mapping_file_replicates.txt": 34411136,
}
FIGSHARE_BASE = "https://ndownloader.figshare.com/files"

# Headline constants from the Cancer Cell 2022 paper (Gonçalves et al.).
PROCAN_PROTEINS = 8498                # quantified at proteotypic Global.Q.Value <= 0.01
PROCAN_PROTEINS_STRINGENT = 6692      # >= 2 supporting peptides
PROCAN_LIBRARY_PRECURSORS = 144578    # spectral library size (NOT identifications)
PROCAN_LIBRARY_PROTEINS = 12487       # spectral library protein count
PROCAN_MS_RUNS = 6864                 # final dataset (paper); PRIDE lists 6,981

# pg_matrix.tsv metadata columns (everything before per-run sample columns).
# Matches the DIA-NN 2.5 output deposited for PXD030304 — older DIA-NN
# versions add a Protein.Ids column; we only require the columns actually
# present so the schema check stays accurate.
PG_METADATA_COLS = [
    "Protein.Group", "Protein.Names", "Genes",
    "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences",
]


@dataclass(frozen=True)
class Counts:
    procan_proteins: int                       # 8,498 (paper)
    procan_proteins_stringent: int             # 6,692 (paper)
    quantmsdiann_proteins: int                 # post-filter pg_matrix headline
    quantmsdiann_proteins_unfiltered: int      # diannsummary.log (audit only)
    quantmsdiann_proteins_pg_matrix_unfiltered: int  # raw pg_matrix row count
    quantmsdiann_proteins_stringent: int       # >=2 unique peptides, computed
    quantmsdiann_precursors: int               # from diannsummary.log


SUMMARY_LOG_PRECURSOR_LINE_RE = re.compile(
    r"Target precursors at 1% global q-value:\s*(\d+)"
)


def parse_diann_summary_log(log_path: Path) -> tuple[int, int]:
    """Return (protein_groups, target_precursors) from a DIA-NN summary log."""
    protein_groups = precursors = None
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


# ---------------------------------------------------------------------------
# Figshare mapping + protein matrix
# ---------------------------------------------------------------------------


def parse_procan_mapping(mapping_path: Path) -> dict[str, str]:
    """Parse the ProCan-DepMapSanger averaged-sample mapping file and return
    `Cell_line -> Tissue_type`. The 28 tissue categories are the paper's
    canonical condition axis."""
    df = pd.read_csv(mapping_path, sep="\t", dtype=str)
    needed = ["Cell_line", "Tissue_type"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"ProCan mapping missing columns: {missing}")
    out: dict[str, str] = {}
    for cell, tissue in zip(df["Cell_line"], df["Tissue_type"]):
        if isinstance(cell, str) and isinstance(tissue, str):
            out[cell.strip()] = tissue.strip()
    return out


def parse_procan_replicates_mapping(
    mapping_path: Path,
    *,
    exclude_hek293t: bool = True,
) -> dict[str, str]:
    """Return `Automatic_MS_filename -> Tissue_type` from the per-replicate
    ProCan mapping file. By default the 1,064 HEK293T QC runs (Tissue_type
    == 'Control_HEK293T') are excluded so they don't contaminate per-tissue
    counts."""
    df = pd.read_csv(mapping_path, sep="\t", dtype=str)
    needed = ["Automatic_MS_filename", "Tissue_type"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"ProCan replicates mapping missing: {missing}")
    out: dict[str, str] = {}
    for run, tissue in zip(df["Automatic_MS_filename"], df["Tissue_type"]):
        if not (isinstance(run, str) and isinstance(tissue, str)):
            continue
        tissue = tissue.strip()
        if exclude_hek293t and tissue.startswith("Control_HEK"):
            continue
        out[run.strip()] = tissue
    return out


def proteins_per_tissue_procan(
    replicates_matrix_path: Path,
    replicates_mapping_path: Path,
    *,
    chunksize: int = 500,
) -> dict[str, set[str]]:
    """For each ProCan tissue, the union of detected protein IDs across all
    MS runs mapped to that tissue. Applies identical per-MS-run union
    semantics to the quantmsdiann side, so the per-tissue comparison is
    apples-to-apples.

    Detection: cell value in the per-replicate matrix is non-NA. The matrix
    has one row per MS run (`Automatic_MS_filename` in the first column —
    header cell is blank in the file) and one column per protein
    (`<accession>;<name>`). HEK293T QC runs are dropped via the mapping
    file. Chunked-read because the file is ~519 MB."""
    run_to_tissue = parse_procan_replicates_mapping(replicates_mapping_path)
    out: dict[str, set[str]] = {}
    # Header row's first cell is empty; tell pandas the first column has no
    # name and rename it to 'run' on read.
    for chunk in pd.read_csv(
        replicates_matrix_path, sep="\t", dtype=str, chunksize=chunksize,
    ):
        first_col = chunk.columns[0]
        chunk = chunk.rename(columns={first_col: "run"})
        protein_cols = [c for c in chunk.columns if c != "run"]
        chunk["tissue"] = chunk["run"].map(
            lambda r: run_to_tissue.get(r.strip()) if isinstance(r, str) else None
        )
        chunk = chunk[chunk["tissue"].notna()]
        if chunk.empty:
            continue
        # For each tissue, union of proteins with non-NA value in any row
        # of this chunk; vectorised so we don't iterate ~58M cells.
        for tissue, group in chunk.groupby("tissue"):
            any_detected = group[protein_cols].notna().any(axis=0)
            detected = [c for c in protein_cols if bool(any_detected[c])]
            out.setdefault(tissue, set()).update(detected)
    return out


# ---------------------------------------------------------------------------
# SDRF + DIA-NN pg_matrix -> per-tissue protein sets
# ---------------------------------------------------------------------------


def load_sdrf_data_file_to_cell_line(sdrf_path: Path) -> dict[str, str]:
    """Same logic as the PXD003539 SDRF loader: parse `comment[data file]` ->
    `characteristics[cell line]`, rewriting `.wiff` -> `.mzML` so DIA-NN
    matrix column names match.

    Duplicated here (rather than imported) because the PXD003539 helper sits
    inside that module's namespace and we want PXD030304 to stay isolated."""
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
        mzml = re.sub(r"\.wiff$", ".mzML", data_file)
        out[mzml] = cell
    return out


PR_METADATA_COLS = [
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
    "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
]


def _compute_or_load_diann_procan_filter(
    cache_path: Path,
    parquet_source: str,
    sdrf_path: Path,
    procan_mapping_path: Path,
) -> dict[str, set[str]]:
    """Side-cache wrapper around `proteins_per_tissue_quantmsdiann_procan_filter`:
    streaming 33 GB over HTTP takes ~15 minutes, so we persist the per-tissue
    result as a small JSON (tissue -> sorted list of Protein.Group). Delete
    the JSON to force a fresh stream."""
    import json
    if cache_path.exists() and cache_path.stat().st_size > 0:
        with open(cache_path, encoding="utf-8") as fh:
            payload = json.load(fh)
        return {t: set(vs) for t, vs in payload.items()}
    result = proteins_per_tissue_quantmsdiann_procan_filter(
        parquet_source, sdrf_path, procan_mapping_path,
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump({t: sorted(s) for t, s in result.items()}, fh)
    return result


def proteins_per_tissue_quantmsdiann_procan_filter(
    parquet_source: str | Path,
    sdrf_path: Path,
    procan_mapping_path: Path,
    *,
    qvalue_cutoff: float = 0.01,
    batch_size: int = 1_000_000,
) -> dict[str, set[str]]:
    """ProCan-style filter applied to quantmsdiann's long-format report.

    For each ProCan tissue, the set of Protein.Group values that have any
    proteotypic precursor passing Global.Q.Value <= `qvalue_cutoff` in at
    least one MS run mapped to that tissue. Matches Gonçalves et al. 2022
    Methods ("filtered to retain only precursors from proteotypic peptides
    with Global.Q.Value <= 0.01") — global-FDR-only, no per-cell quant FDR.

    `parquet_source` is either a local Path or an HTTPS URL. For the URL
    case we stream the parquet via fsspec's HTTPFileSystem with column
    projection on (`Run`, `Protein.Group`, `Global.Q.Value`, `Proteotypic`);
    only those columns' chunks transit the wire so the 33 GB file never
    needs to be staged locally."""
    import pyarrow.parquet as pq

    sdrf_run_to_cell = load_sdrf_data_file_to_cell_line(sdrf_path)
    # SDRF data files end in `.mzML` (after the .wiff -> .mzML rewrite the
    # loader does); the parquet's `Run` column carries the bare basename
    # without extension. Strip both possible suffixes.
    sdrf_no_ext = {}
    for k, v in sdrf_run_to_cell.items():
        stem = re.sub(r"\.(mzML|wiff)$", "", k, flags=re.IGNORECASE)
        sdrf_no_ext[stem] = v
    cl_to_tissue = parse_procan_mapping(procan_mapping_path)

    out: dict[str, set[str]] = {}
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
                cell = sdrf_no_ext.get(r)
                if cell is None:
                    continue
                tissue = cl_to_tissue.get(cell)
                if tissue is None:
                    continue
                out.setdefault(tissue, set()).add(pg)
    return out


def proteins_per_tissue_quantmsdiann(
    pr_matrix_path: Path,
    sdrf_path: Path,
    procan_mapping_path: Path,
    *,
    chunksize: int = 50_000,
) -> dict[str, set[str]]:
    """For each ProCan tissue, the union of `Protein.Group` values whose
    precursors are quantified in at least one run belonging to that tissue.

    Uses `diann_report.pr_matrix.tsv` (1% precursor + 1% protein-group FDR per
    cell) rather than `pg_matrix.tsv` (5% PG q-value per cell) so we exercise
    DIA-NN's strictest per-cell filter, matching ProCan's per-replicate "any
    precursor identified → protein detected" semantics more closely. Mapping:
    pr_matrix run column -> SDRF cell line -> ProCan Tissue_type.

    Chunked-read because pr_matrix.tsv is ~2 GB for PXD030304."""
    sdrf_run_to_cell = load_sdrf_data_file_to_cell_line(sdrf_path)
    cl_to_tissue = parse_procan_mapping(procan_mapping_path)

    reader = pd.read_csv(
        pr_matrix_path, sep="\t", dtype=str, chunksize=chunksize,
    )
    out: dict[str, set[str]] = {}
    col_to_tissue: dict[str, str] | None = None
    for chunk in reader:
        if col_to_tissue is None:
            missing = [c for c in PR_METADATA_COLS if c not in chunk.columns]
            if missing:
                raise ValueError(f"pr_matrix missing metadata columns: {missing}")
            sample_cols = [c for c in chunk.columns if c not in PR_METADATA_COLS]
            col_to_tissue = {}
            for col in sample_cols:
                cell = sdrf_run_to_cell.get(col)
                if not cell:
                    continue
                tissue = cl_to_tissue.get(cell)
                if tissue is None:
                    continue
                col_to_tissue[col] = tissue
        protein_groups = chunk["Protein.Group"].tolist()
        for col, tissue in col_to_tissue.items():
            mask = chunk[col].notna()
            bucket = out.setdefault(tissue, set())
            for pg, ok in zip(protein_groups, mask):
                if ok and isinstance(pg, str):
                    bucket.add(pg)
    return out


# ---------------------------------------------------------------------------
# Per-run completeness
# ---------------------------------------------------------------------------


def per_run_completeness_procan(peptide_counts_path: Path) -> dict[str, float]:
    """For each MS run row in
    `ProCan-DepMapSanger_peptide_counts_per_protein_per_sample.txt`, the
    fraction of proteins with peptide count > 0. Denominator = total proteins
    in the matrix (typically 8,498).

    We stream the file row-by-row because it's 143 MB and we only need a
    count per row (not the full matrix)."""
    n_proteins: int | None = None
    out: dict[str, float] = {}
    with open(peptide_counts_path, encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "Run":
            raise ValueError(
                f"peptide counts file first column should be 'Run', got "
                f"{header[0]!r}"
            )
        n_proteins = len(header) - 1
        for row in reader:
            if not row:
                continue
            run = row[0]
            n_detected = sum(
                1 for v in row[1:]
                if v and v.strip() and _is_positive_count(v)
            )
            out[run] = n_detected / n_proteins if n_proteins else 0.0
    return out


def _is_positive_count(value: str) -> bool:
    """Return True if `value` parses to a number greater than zero."""
    try:
        return float(value) > 0
    except ValueError:
        return False


def per_run_completeness_quantmsdiann(pg_matrix_path: Path) -> dict[str, float]:
    """For each per-run column in pg_matrix.tsv, the fraction of protein-group
    rows that are non-NA. Denominator = total protein groups in the matrix
    (the quantmsdiann pipeline's identified set)."""
    df = pd.read_csv(pg_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PG_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"pg_matrix missing metadata columns: {missing}")
    sample_cols = [c for c in df.columns if c not in PG_METADATA_COLS]
    n_total = len(df)
    if n_total == 0:
        return {c: 0.0 for c in sample_cols}
    out: dict[str, float] = {}
    for col in sample_cols:
        out[col] = int(df[col].notna().sum()) / n_total
    return out


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def render_main_figure(
    counts: Counts,
    svg_path: Path,
) -> None:
    """Grouped bar chart: 2 conditions x 2 metrics (proteins, ≥2-peptide
    proteins). Paper-ready: no title, no footer."""
    metrics = ["Protein groups", "Protein groups\n($\\geq$2 unique peptides)"]
    procan_vals = [counts.procan_proteins, counts.procan_proteins_stringent]
    diann_vals = [counts.quantmsdiann_proteins,
                  counts.quantmsdiann_proteins_stringent]

    fig, ax = plt.subplots(figsize=(7, 5))
    bar_width = 0.35
    x = list(range(len(metrics)))
    bars_p = ax.bar([xi - bar_width / 2 for xi in x], procan_vals,
                    width=bar_width, color="#9e9e9e",
                    label="ProCan-DepMapSanger 2022")
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals,
                    width=bar_width, color="#1f77b4",
                    label="quantmsdiann (DIA-NN)")
    for bars, vals in ((bars_p, procan_vals), (bars_d, diann_vals)):
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{v:,}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("Protein groups (1% FDR)")
    ymax = max(max(procan_vals), max(diann_vals)) * 1.18
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", frameon=False)

    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


def render_proteins_per_tissue(
    procan_per_tissue: dict[str, set[str]],
    diann_per_tissue: dict[str, set[str]],
    svg_path: Path,
) -> None:
    """Grouped bar chart: 28 tissues x 2 conditions (ProCan vs quantmsdiann).
    Tissues sorted by descending ProCan cell-line count is not directly
    available here; we sort by descending union size across the two
    pipelines (largest tissues first). Paper-ready: no title, no footer."""
    tissues = sorted(
        set(procan_per_tissue) | set(diann_per_tissue),
        key=lambda t: -(
            len(procan_per_tissue.get(t, set()))
            + len(diann_per_tissue.get(t, set()))
        ),
    )
    procan_vals = [len(procan_per_tissue.get(t, set())) for t in tissues]
    diann_vals = [len(diann_per_tissue.get(t, set())) for t in tissues]

    fig, ax = plt.subplots(figsize=(13, 6.5))
    x = list(range(len(tissues)))
    bar_width = 0.4
    bars_p = ax.bar([xi - bar_width / 2 for xi in x], procan_vals,
                    width=bar_width, color="#9e9e9e",
                    label="ProCan-DepMapSanger 2022")
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals,
                    width=bar_width, color="#1f77b4",
                    label="quantmsdiann (DIA-NN)")
    # Value labels intentionally omitted: with 28 tissues × 2 bars and values
    # often within 1-5% of each other, the numeric annotations overlap and
    # become illegible. Exact numbers live in `counts.tsv`.

    ax.set_xticks(x)
    ax.set_xticklabels(tissues, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Distinct protein groups detected")
    ymax = max(max(procan_vals, default=0), max(diann_vals, default=0))
    ax.set_ylim(0, ymax * 1.15 if ymax else 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28),
              ncol=2, frameon=False)

    fig.tight_layout(rect=(0, 0.18, 1, 1))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


def render_per_run_completeness(
    procan_per_run: dict[str, float],
    diann_per_run: dict[str, float],
    svg_path: Path,
) -> None:
    """Line plot: per-run fraction of detected proteins, both pipelines.
    Runs aligned by sorting each pipeline's set independently — they don't
    share a common run-naming convention (ProCan: `180822_e0022_p02_*_s_m04_1`;
    quantmsdiann: `180822_E0022_P02_*_S_M04_1.mzML`). Paper-ready: no title,
    no footer."""
    procan_vals = sorted(procan_per_run.values())
    diann_vals = sorted(diann_per_run.values())

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(range(len(procan_vals)), procan_vals,
            color="#9e9e9e", linewidth=0.8,
            label=f"ProCan-DepMapSanger 2022 (n={len(procan_vals)} runs)")
    ax.plot(range(len(diann_vals)), diann_vals,
            color="#1f77b4", linewidth=0.8,
            label=f"quantmsdiann (DIA-NN) (n={len(diann_vals)} runs)")
    ax.set_xlabel("Run rank (sorted ascending within each pipeline)")
    ax.set_ylabel("Fraction of protein groups\ndetected per run")
    ax.set_ylim(0, 1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False)

    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


def write_counts_tsv(
    counts: Counts,
    tsv_path: Path,
    *,
    procan_per_tissue: dict[str, set[str]] | None = None,
    diann_per_tissue: dict[str, set[str]] | None = None,
) -> None:
    """Auditable counts table. Optionally appends per-tissue protein counts
    (one row per tissue per pipeline) so the supp figure's numbers are also
    machine-readable."""
    rows = [
        ("Protein groups", "ProCan-DepMapSanger 2022 (paper headline)",
         counts.procan_proteins,
         "8,498 proteins quantified at Global.Q.Value <= 0.01 (paper, Results)"),
        ("Protein groups (>=2 peptides)", "ProCan-DepMapSanger 2022 (stringent)",
         counts.procan_proteins_stringent,
         "6,692 proteins with >=2 supporting peptides (paper, Results)"),
        ("Protein groups", "quantmsdiann (DIA-NN, 1% FDR, target-only)",
         counts.quantmsdiann_proteins,
         "post-filter: pg_matrix.tsv rows whose Protein.Group has no "
         "CONTAM_/Cont_/ENTRAP_/DECOY_/decoy_ token (conservative filter, "
         "2026-05-21 spec)"),
        ("Protein groups", "quantmsdiann (DIA-NN, 1% FDR, unfiltered pg_matrix)",
         counts.quantmsdiann_proteins_pg_matrix_unfiltered,
         "raw row count of diann_report.pg_matrix.tsv, pre-filter; includes "
         "CONTAM_/ENTRAP_/DECOY_ rows"),
        ("Protein groups", "quantmsdiann (DIA-NN, 1% FDR, diannsummary.log)",
         counts.quantmsdiann_proteins_unfiltered,
         "audit baseline: diannsummary.log 'Protein groups with global "
         "q-value <= 0.01' line (unfiltered)"),
        ("Protein groups (>=2 peptides)", "quantmsdiann (DIA-NN, 1% FDR)",
         counts.quantmsdiann_proteins_stringent,
         ">=2 unique Stripped.Sequence per Protein.Group (proteotypic) in pr_matrix.tsv"),
        ("Precursors", "quantmsdiann (DIA-NN, 1% FDR)",
         counts.quantmsdiann_precursors,
         "from diannsummary.log (Target precursors at 1% global q-value)"),
        ("Spectral library precursors", "ProCan-DepMapSanger 2022",
         PROCAN_LIBRARY_PRECURSORS,
         "library size (paper, STAR Methods); NOT identified precursors"),
        ("Spectral library proteins", "ProCan-DepMapSanger 2022",
         PROCAN_LIBRARY_PROTEINS,
         "library protein count (paper, STAR Methods)"),
        ("MS runs", "ProCan-DepMapSanger 2022",
         PROCAN_MS_RUNS,
         "paper headline (PRIDE archive lists 6,981)"),
    ]
    if procan_per_tissue is not None or diann_per_tissue is not None:
        procan = procan_per_tissue or {}
        diann = diann_per_tissue or {}
        tissues = sorted(
            set(procan) | set(diann),
            key=lambda t: -(len(procan.get(t, set())) + len(diann.get(t, set()))),
        )
        note_procan = (
            "per-MS-run union over protein_matrix_8498_replicates.txt "
            "(HEK293T QC runs excluded); reflects ProCan's global 1% "
            "Global.Q.Value filtering, NOT per-cell strict identification"
        )
        note_diann = (
            "per-tissue union of Protein.Group from diann_report.parquet "
            "filtered to Proteotypic == 1 AND Global.Q.Value <= 0.01 "
            "(ProCan's filter applied to the long-format report; no "
            "per-cell quant FDR)"
        )
        for t in tissues:
            rows.append((f"Per-tissue proteins | {t}",
                         "ProCan-DepMapSanger 2022 (per-replicate)",
                         len(procan.get(t, set())), note_procan))
            rows.append((f"Per-tissue proteins | {t}",
                         "quantmsdiann (DIA-NN, 1% FDR)",
                         len(diann.get(t, set())), note_diann))
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

    # Reanalysis inputs
    log_path = download_if_missing(DIANN_SUMMARY_LOG_URL,
                                   DATA_DIR / "diannsummary.log")
    pg_path = download_if_missing(DIANN_PG_MATRIX_URL,
                                  DATA_DIR / "diann_report.pg_matrix.tsv")
    pr_path = download_if_missing(DIANN_PR_MATRIX_URL,
                                  DATA_DIR / "diann_report.pr_matrix.tsv")
    sdrf_path = download_if_missing(QUANTMS_SDRF_URL,
                                    DATA_DIR / "PXD030304.sdrf.tsv")

    # ProCan figshare inputs
    fs_paths: dict[str, Path] = {}
    for name, fid in FIGSHARE_FILES.items():
        fs_paths[name] = download_if_missing(
            f"{FIGSHARE_BASE}/{fid}", DATA_DIR / name,
        )

    # Headline counts
    print("Parsing DIA-NN summary log...")
    pg_log, prec = parse_diann_summary_log(log_path)
    print(f"  protein groups (log, unfiltered): {pg_log:,}  precursors: {prec:,}")

    print("Counting pg_matrix.tsv protein-group rows (unfiltered + target-only)...")
    pg_unf, pg_target = count_target_protein_groups(pg_path)
    print(f"  pg_matrix rows: unfiltered={pg_unf:,}  target_only={pg_target:,} "
          f"(delta {pg_unf - pg_target:,})")

    print("Computing quantmsdiann >=2-peptide protein-group count...")
    pep_per_pg = unique_peptides_per_protein_diann(pr_path)
    diann_stringent = proteins_with_min_peptides(pep_per_pg, 2)
    print(f"  >=2 unique peptides: {diann_stringent:,}")

    counts = Counts(
        procan_proteins=PROCAN_PROTEINS,
        procan_proteins_stringent=PROCAN_PROTEINS_STRINGENT,
        # Headline = post-filter pg_matrix count (per 2026-05-21 spec §1.6).
        quantmsdiann_proteins=pg_target,
        quantmsdiann_proteins_unfiltered=pg_log,
        quantmsdiann_proteins_pg_matrix_unfiltered=pg_unf,
        quantmsdiann_proteins_stringent=diann_stringent,
        quantmsdiann_precursors=prec,
    )

    print("Rendering main figure...")
    render_main_figure(
        counts,
        FIGURES_DIR / "main_comparison.svg",
    )

    print("Computing per-tissue protein sets (ProCan, per-replicate, HEK293T excluded)...")
    procan_per_tissue = proteins_per_tissue_procan(
        fs_paths["protein_matrix_8498_replicates.txt"],
        fs_paths["mapping_file_replicates.txt"],
    )
    print(f"  {len(procan_per_tissue)} tissues")

    print("Computing per-tissue protein sets (quantmsdiann, ProCan-style filter)...")
    # Apply Gonçalves et al. 2022's filter (Proteotypic == 1 AND
    # Global.Q.Value <= 0.01, no per-cell quant FDR) so the per-tissue
    # comparison is methodologically equivalent on both sides. Streams the
    # 33 GB diann_report.parquet over HTTP with column projection on the
    # 4 columns we need; the result is cached to a small JSON so subsequent
    # runs are instant.
    diann_per_tissue = _compute_or_load_diann_procan_filter(
        DATA_DIR / "diann_per_tissue_procan_filter.json",
        DIANN_PARQUET_URL,
        sdrf_path,
        fs_paths["mapping_file_averaged.txt"],
    )
    print(f"  {len(diann_per_tissue)} tissues")

    print("Rendering per-tissue supp figure...")
    render_proteins_per_tissue(
        procan_per_tissue, diann_per_tissue,
        FIGURES_DIR / "supp_proteins_per_tissue.svg",
    )

    print("Computing per-run completeness (ProCan)...")
    procan_per_run = per_run_completeness_procan(
        fs_paths["peptide_counts_per_protein_per_sample.txt"],
    )
    print(f"  {len(procan_per_run)} runs")

    print("Computing per-run completeness (quantmsdiann)...")
    diann_per_run = per_run_completeness_quantmsdiann(pg_path)
    print(f"  {len(diann_per_run)} runs")

    print("Rendering per-run completeness supp figure...")
    render_per_run_completeness(
        procan_per_run, diann_per_run,
        FIGURES_DIR / "supp_missing_values_per_run.svg",
    )

    print("Writing auditable counts TSV (with per-tissue rows)...")
    data_dir = FIGURES_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    write_counts_tsv(
        counts, data_dir / "counts.tsv",
        procan_per_tissue=procan_per_tissue,
        diann_per_tissue=diann_per_tissue,
    )

    print("Per-tissue protein counts (ProCan | quantmsdiann):")
    all_t = sorted(
        set(procan_per_tissue) | set(diann_per_tissue),
        key=lambda t: -(
            len(procan_per_tissue.get(t, set()))
            + len(diann_per_tissue.get(t, set()))
        ),
    )
    for t in all_t:
        p = len(procan_per_tissue.get(t, set()))
        d = len(diann_per_tissue.get(t, set()))
        print(f"  {t:32s} {p:>6,} | {d:>6,}")

    # Cross-check (the diannsummary.log unfiltered number is what the
    # historical baseline pinned; the post-filter target count is lower).
    if pg_log != 9370:
        print(f"WARN: quantmsdiann protein groups (log) {pg_log} != "
              f"expected 9,370",
              file=sys.stderr)
    if prec != 153644:
        print(f"WARN: quantmsdiann precursors {prec} != expected 153,644",
              file=sys.stderr)

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
