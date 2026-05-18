from __future__ import annotations

import csv
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "PXD003539"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures"

PRIDE_QUANT_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/quant_tables"
)
PRIDE_SUBMISSION_BASE = (
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD003539"
)

PR_MATRIX_URL = f"{PRIDE_QUANT_BASE}/diann_report.pr_matrix.tsv"
SUMMARY_LOG_URL = f"{PRIDE_QUANT_BASE}/diannsummary.log"
OPENSWATH_MATRIX_URL = f"{PRIDE_SUBMISSION_BASE}/feature_alignment_requant_matrix.tsv"

PR_METADATA_COLS = [
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
    "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
]

# Walzer et al. 2022 (doi:10.1038/s41597-022-01380-9) Supplementary Table S2,
# row PXD003539, 1% FDR, 'top3' protein inference, unfiltered.
WALZER_PEPTIDES = 77014
WALZER_PROTEINS = 7097
# Auxiliary (TSV only): the same row's '50% per group' consistency-filtered
# protein count is 6,867; with 'all' inference at 1% FDR the protein count is
# 5,412. Walzer 2022 does not directly report a precursor-level count.
WALZER_PROTEINS_50PCT_FILTER = 6867

# Expression Atlas E-PROT-73 — Walzer 2022 reanalysis of PXD003539 after their
# full post-processing pipeline (Ensembl gene mapping + 50%-per-group
# consistency filter + decoy removal). Smaller than the raw 7,097 because
# protein-to-gene mapping is lossy.
EPROT73_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/"
    "experiments/E-PROT-73/E-PROT-73.tsv"
)

# Guo et al. 2019 (doi:10.1016/j.isci.2019.10.059) Results, after DIA-expert
# manual curation. Reported in the paper text and matched in Walzer 2022
# Supplementary Table S2, "Original - after filter" column for PXD003539.
# Used for TSV context, NOT plotted as a headline bar.
GUO_CURATED_PEPTIDES = 22554
GUO_CURATED_PROTEINS = 3171

SUMMARY_LOG_PROTEIN_LINE_RE = re.compile(
    r"Protein groups with global q-value <= 0\.01:\s*(\d+)"
)


@dataclass(frozen=True)
class Counts:
    guo_peptides: int
    guo_proteins: int
    guo_precursors: int
    walzer_peptides: int
    walzer_proteins: int
    walzer_ea_genes: int  # new
    diann_peptides: int
    diann_proteins: int
    diann_precursors: int


def download_if_missing(url: str, dest: Path, *, retries: int = 2) -> Path:
    """Download url to dest, skipping if dest already exists and is non-empty."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    last_exc: Exception | None = None
    try:
        for attempt in range(retries + 1):
            try:
                with requests.get(url, stream=True, timeout=120) as resp:
                    resp.raise_for_status()
                    with part.open("wb") as fh:
                        for chunk in resp.iter_content(chunk_size=1 << 20):
                            fh.write(chunk)
                os.replace(part, dest)
                return dest
            except requests.RequestException as exc:
                last_exc = exc
                if attempt < retries:
                    time.sleep(2 ** attempt)
        raise RuntimeError(
            f"Failed to download {url} after {retries + 1} attempts: {last_exc}"
        )
    finally:
        part.unlink(missing_ok=True)


def count_quantified_rows(
    matrix_path: Path,
    metadata_cols: list[str],
    unique_by: str | None = None,
) -> int:
    """Count quantified rows in a DIA-NN-style TSV matrix."""
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    missing = [c for c in metadata_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Matrix is missing expected metadata column(s): {', '.join(missing)}"
        )
    sample_cols = [c for c in df.columns if c not in metadata_cols]
    sample_df = df[sample_cols]
    # read_csv treats "" and "NA" as NaN by default
    quantified_mask = sample_df.notna().any(axis=1)
    quantified_df = df[quantified_mask]
    if unique_by is not None:
        return int(quantified_df[unique_by].nunique())
    return int(quantified_mask.sum())


PEPTIDE_ID_RE = re.compile(r"^(?:DECOY_)?\d+_(?P<modseq>.+)_(?P<charge>\d+)_run0$")
UNIMOD_RE = re.compile(r"\(UniMod:\d+\)")


def _stripped_peptide(peptide_id: str) -> str | None:
    """Extract the unmodified peptide sequence from an OpenSWATH Peptide ID.

    Returns None if the ID does not match the expected format."""
    m = PEPTIDE_ID_RE.match(peptide_id)
    if not m:
        return None
    return UNIMOD_RE.sub("", m.group("modseq"))


def count_openswath_quantified(matrix_path: Path) -> tuple[int, int, int]:
    """Count quantified target precursors, peptides, and protein groups in an OpenSWATH matrix."""
    header_df = pd.read_csv(matrix_path, sep="\t", nrows=0)
    cols = list(header_df.columns)
    if "Peptide" not in cols:
        raise ValueError("Matrix is missing required column: Peptide")
    if "Protein" not in cols:
        raise ValueError("Matrix is missing required column: Protein")
    intensity_cols = [c for c in cols if c.startswith("Intensity_")]
    if not intensity_cols:
        raise ValueError(
            "Matrix has zero Intensity_* columns; expected at least one"
        )
    usecols = ["Peptide", "Protein"] + intensity_cols
    total_precursors = 0
    peptides: set[str] = set()
    protein_groups: set[str] = set()
    for chunk in pd.read_csv(
        matrix_path,
        sep="\t",
        usecols=usecols,
        dtype=str,
        chunksize=50000,
    ):
        # Exclude decoy rows
        is_decoy = (
            chunk["Peptide"].str.startswith("DECOY_", na=False)
            | chunk["Protein"].str.contains("DECOY", case=False, na=False)
        )
        target = chunk[~is_decoy]
        # A row is quantified if at least one Intensity_ column is non-NA / non-empty
        # read_csv treats "" and "NA" as NaN by default
        quantified_mask = target[intensity_cols].notna().any(axis=1)
        quantified = target[quantified_mask]
        total_precursors += len(quantified)
        for pep_id in quantified["Peptide"]:
            stripped = _stripped_peptide(pep_id)
            if stripped is not None:
                peptides.add(stripped)
        protein_groups.update(quantified["Protein"].tolist())
    return total_precursors, len(peptides), len(protein_groups)


def count_eprot73_genes(tsv_path: Path) -> int:
    """Count unique Ensembl gene IDs in the Expression Atlas E-PROT-73 file.

    The file has no preamble: line 1 is the header (`Gene ID`, `Gene Name`, …)
    and data rows follow immediately from line 2.  We read all rows with
    `pd.read_csv` (no `skiprows`) and keep only those whose first column value
    starts with `ENSG`.  The header value `"Gene ID"` is naturally excluded by
    that prefix filter, so no extra row-dropping is needed."""
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, usecols=[0])
    df.columns = ["gene_id"]
    df = df[df["gene_id"].fillna("").str.startswith("ENSG")]
    return int(df["gene_id"].nunique())


def per_run_non_na_fraction_diann(matrix_path: Path) -> dict[str, float]:
    """For each sample column in a DIA-NN matrix, return the fraction of rows
    with a non-NA value in that column. Sample columns are anything not in
    PR_METADATA_COLS."""
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"DIA-NN matrix header missing metadata columns: {missing}"
        )
    samples = [c for c in df.columns if c not in PR_METADATA_COLS]
    total = len(df)
    if total == 0:
        return {s: 0.0 for s in samples}
    return {s: float(df[s].notna().sum()) / total for s in samples}


def per_run_non_na_fraction_openswath(matrix_path: Path) -> dict[str, float]:
    """For each Intensity_<run> column in the OpenSWATH matrix, return the
    fraction of target (non-decoy) rows with a non-NA value. Uses chunked
    reading. Total denominator is the count of target rows (the same for
    every run)."""
    # First pass: count target rows.
    target_total = 0
    intensity_cols: list[str] = []
    header = pd.read_csv(matrix_path, sep="\t", nrows=0)
    cols = list(header.columns)
    if "Peptide" not in cols or "Protein" not in cols:
        raise ValueError("OpenSWATH matrix missing Peptide/Protein columns")
    intensity_cols = [c for c in cols if c.startswith("Intensity_")]
    if not intensity_cols:
        raise ValueError("OpenSWATH matrix has no Intensity_* columns")

    per_run_non_na: dict[str, int] = {c: 0 for c in intensity_cols}
    usecols = ["Peptide", "Protein"] + intensity_cols
    for chunk in pd.read_csv(matrix_path, sep="\t", dtype=str,
                             usecols=usecols, chunksize=50000):
        is_decoy = (
            chunk["Peptide"].str.startswith("DECOY_", na=False)
            | chunk["Protein"].str.contains("DECOY", case=False, na=False)
        )
        targets = chunk[~is_decoy]
        target_total += len(targets)
        for col in intensity_cols:
            per_run_non_na[col] += int(targets[col].notna().sum())

    if target_total == 0:
        return {c: 0.0 for c in intensity_cols}
    return {c: per_run_non_na[c] / target_total for c in intensity_cols}


def parse_summary_log(log_path: Path) -> int:
    """Return protein group count at 1% global FDR from a DIA-NN summary log."""
    with open(log_path, encoding="utf-8") as fh:
        for line in fh:
            m = SUMMARY_LOG_PROTEIN_LINE_RE.search(line)
            if m:
                return int(m.group(1))
    raise ValueError(
        "Line matching 'Protein groups with global q-value <= 0.01:' not found in log"
    )


def render_figure(counts: Counts, pdf_path: Path, png_path: Path) -> None:
    """Render a 2-condition x 2-metric grouped bar chart and save as PDF and PNG."""
    conditions = [
        ("Guo 2019\n(OpenSWATH)", "#9e9e9e",
         counts.guo_peptides, counts.guo_proteins),
        ("quantmsdiann\n(DIA-NN)", "#1f77b4",
         counts.diann_peptides, counts.diann_proteins),
    ]

    metrics = ["Peptides", "Protein groups"]
    bar_width = 0.27
    n_conditions = len(conditions)
    x = [0, 1]  # x positions for the two metric groups

    fig, ax = plt.subplots(figsize=(7, 5))

    offsets = [bar_width * (i - (n_conditions - 1) / 2.0)
               for i in range(n_conditions)]

    for i, (label, color, peptide_val, protein_val) in enumerate(conditions):
        values = [peptide_val, protein_val]
        bars = ax.bar(
            [xi + offsets[i] for xi in x],
            values,
            width=bar_width,
            color=color,
            label=label,
        )
        for bar, val in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height(),
                f"{val:,}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    # Check if log scale is needed for either metric
    peptide_vals = [c[2] for c in conditions]
    protein_vals = [c[3] for c in conditions]
    needs_log = False
    for metric_vals in [peptide_vals, protein_vals]:
        mn, mx = min(metric_vals), max(metric_vals)
        if mn > 0 and mx / mn > 5:
            needs_log = True
            break

    ylabel = "Count (1% FDR)"
    if needs_log:
        ax.set_yscale("log")
        ylabel += " (log scale)"
    else:
        # Headroom so the legend doesn't overlap the tallest bar's value label.
        top = max(max(peptide_vals), max(protein_vals)) * 1.18
        ax.set_ylim(0, top)

    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(loc="upper left", frameon=False)

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(
        0.5, 0.02,
        "All counts at 1% FDR (Guo 2019: OpenSWATH + pyprophet; this work: DIA-NN).",
        ha="center",
        va="bottom",
        fontstyle="italic",
        fontsize=8,
    )

    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def render_missing_values_per_run(
    guo_per_run: dict[str, float],
    diann_per_run: dict[str, float],
    pdf_path: Path,
    png_path: Path,
) -> None:
    """Plot per-run non-NA fraction for both pipelines.

    Runs on the x axis are aligned by the L<date>_<n>_SW stem extracted from
    the column name. Bars sorted by the stem (which sorts chronologically by
    acquisition date)."""
    import re as _re

    def stem(col: str) -> str:
        # diann columns look like 'guot_L130610_003_SW.mzML'
        # openswath columns look like 'Intensity_guot_L130610_003_SW_with_dscore.csv_0_14'
        m = _re.search(r"L\d+_\d+_SW", col)
        return m.group(0) if m else col

    diann_by_stem = {stem(c): v for c, v in diann_per_run.items()}
    guo_by_stem = {stem(c): v for c, v in guo_per_run.items()}
    common = sorted(set(diann_by_stem) & set(guo_by_stem))

    fig, ax = plt.subplots(figsize=(10, 4))
    x = range(len(common))
    ax.plot(list(x), [guo_by_stem[s] for s in common],
            label="Guo 2019 (OpenSWATH)", color="#9e9e9e", linewidth=1.0)
    ax.plot(list(x), [diann_by_stem[s] for s in common],
            label="quantmsdiann (DIA-NN)", color="#1f77b4", linewidth=1.0)
    ax.set_xlabel(f"MS run index ({len(common)} runs, ordered by acquisition date)")
    ax.set_ylabel("Fraction of precursors quantified per run")
    ax.set_ylim(0, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False)

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def write_counts_tsv(counts: Counts, tsv_path: Path) -> None:
    """Write an auditable TSV with metric/source/count/note rows."""
    rows = [
        ("Peptides", "Guo 2019 (OpenSWATH, deposited)",
         counts.guo_peptides,
         "from feature_alignment_requant_matrix.tsv (target rows, >=1 quant)"),
        ("Peptides", "Walzer 2022 (CAL + OpenSWATH, top3)",
         counts.walzer_peptides,
         "Supplementary Table S2 (PXD003539, 1% FDR, top3)"),
        ("Peptides", "quantmsdiann (DIA-NN, 1% FDR)",
         counts.diann_peptides,
         "unique Stripped.Sequence in pr_matrix with >=1 non-NA"),
        ("Protein groups", "Guo 2019 (OpenSWATH, deposited)",
         counts.guo_proteins,
         "unique Protein column values in feature_alignment_requant_matrix.tsv"),
        ("Protein groups", "Walzer 2022 (CAL + OpenSWATH, top3)",
         counts.walzer_proteins,
         "Supplementary Table S2 (PXD003539, 1% FDR, top3, unfiltered)"),
        ("Protein groups", "quantmsdiann (DIA-NN, 1% FDR)",
         counts.diann_proteins,
         "from diannsummary.log (Protein groups with global q-value <= 0.01)"),
        ("Precursors aux", "Guo 2019 (OpenSWATH, deposited)",
         counts.guo_precursors,
         "target rows with >=1 non-NA Intensity in feature_alignment matrix"),
        ("Precursors aux", "quantmsdiann (DIA-NN, 1% FDR)",
         counts.diann_precursors,
         "rows in pr_matrix with >=1 non-NA"),
        ("Curated context", "Guo 2019 (DIA-expert curated, peptides)",
         22554,
         "paper text, not used as a headline bar"),
        ("Curated context", "Guo 2019 (DIA-expert curated, proteins)",
         3171,
         "paper text, not used as a headline bar"),
        ("Filter context", "Walzer 2022 (50% per group filter)",
         6867,
         "Supplementary Table S2 - '50% per group' consistency filter, proteins"),
        ("EA context", "Walzer 2022 (E-PROT-73, Expression Atlas)",
         counts.walzer_ea_genes,
         "unique Ensembl gene IDs in E-PROT-73.tsv (post-processed: gene mapping + "
         "50% per-group filter + decoy removal)"),
    ]

    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "source", "count", "note"])
        for row in rows:
            writer.writerow(row)


def main() -> int:  # pragma: no cover
    # 1. Ensure output directories exist
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # 2. Download input files (idempotent)
    pr_path = DATA_DIR / "diann_report.pr_matrix.tsv"
    log_path = DATA_DIR / "diannsummary.log"
    opensw_path = DATA_DIR / "feature_alignment_requant_matrix.tsv"

    print("Downloading files (skipped if already cached)...")
    download_if_missing(PR_MATRIX_URL, pr_path)
    download_if_missing(SUMMARY_LOG_URL, log_path)
    download_if_missing(OPENSWATH_MATRIX_URL, opensw_path)
    eprot73_path = download_if_missing(EPROT73_URL, DATA_DIR / "E-PROT-73.tsv")

    # 3. Compute counts
    print("Computing quantmsdiann counts...")
    diann_peptides = count_quantified_rows(
        pr_path, PR_METADATA_COLS, unique_by="Stripped.Sequence"
    )
    diann_precursors = count_quantified_rows(pr_path, PR_METADATA_COLS)
    diann_proteins = parse_summary_log(log_path)

    print("Computing Guo 2019 (OpenSWATH) counts...")
    guo_precursors, guo_peptides, guo_proteins = count_openswath_quantified(opensw_path)

    print("Computing Expression Atlas E-PROT-73 gene count...")
    walzer_ea_genes = count_eprot73_genes(eprot73_path)

    # Build the immutable Counts record
    counts = Counts(
        guo_peptides=guo_peptides,
        guo_proteins=guo_proteins,
        guo_precursors=guo_precursors,
        walzer_peptides=WALZER_PEPTIDES,
        walzer_proteins=WALZER_PROTEINS,
        walzer_ea_genes=walzer_ea_genes,
        diann_peptides=diann_peptides,
        diann_proteins=diann_proteins,
        diann_precursors=diann_precursors,
    )

    # 4. Render figure
    pdf_path = FIGURES_DIR / "PXD003539_main_comparison.pdf"
    png_path = FIGURES_DIR / "PXD003539_main_comparison.png"
    render_figure(counts, pdf_path, png_path)
    print(f"Figure saved to {pdf_path} and {png_path}")

    # 5. Write auditable TSV
    tsv_path = FIGURES_DIR / "PXD003539_counts.tsv"
    write_counts_tsv(counts, tsv_path)
    print(f"Counts TSV saved to {tsv_path}")

    # 6. Supplementary figure B: per-run completeness
    print("Computing per-run completeness for supplementary figure B...")
    guo_per_run = per_run_non_na_fraction_openswath(opensw_path)
    diann_per_run = per_run_non_na_fraction_diann(pr_path)
    suppB_pdf = FIGURES_DIR / "PXD003539_supp_missing_values_per_run.pdf"
    suppB_png = FIGURES_DIR / "PXD003539_supp_missing_values_per_run.png"
    render_missing_values_per_run(guo_per_run, diann_per_run, suppB_pdf, suppB_png)
    print(f"Supplementary B figure saved to {suppB_pdf} and {suppB_png}")

    # Also print median per-run fraction for each pipeline.
    import statistics as _stats
    print(
        f"Median per-run completeness: "
        f"Guo={_stats.median(guo_per_run.values()):.2%} "
        f"quantmsdiann={_stats.median(diann_per_run.values()):.2%}"
    )

    # 7. Print 3-line summary
    print(
        f"Peptides:        Guo={guo_peptides:,}  "
        f"Walzer={WALZER_PEPTIDES:,}  "
        f"quantmsdiann={diann_peptides:,}"
    )
    print(
        f"Protein groups:  Guo={guo_proteins:,}  "
        f"Walzer={WALZER_PROTEINS:,}  "
        f"quantmsdiann={diann_proteins:,}"
    )
    print(
        f"Auxiliary:       Guo precursors={guo_precursors:,}  "
        f"quantmsdiann precursors={diann_precursors:,}"
    )
    print(
        f"EA context:     Walzer (E-PROT-73 genes)={counts.walzer_ea_genes:,}"
    )

    # 8. Sanity warnings
    if abs(diann_precursors - 117720) / 117720 > 0.01:
        print(
            f"WARNING: diann_precursors={diann_precursors:,} deviates >1% from expected 117,720",
            file=sys.stderr,
        )
    if diann_proteins != 6927:
        print(
            f"WARNING: diann_proteins={diann_proteins:,} != expected 6,927",
            file=sys.stderr,
        )
    if abs(guo_precursors - 48374) / 48374 > 0.01:
        print(
            f"WARNING: guo_precursors={guo_precursors:,} deviates >1% from expected 48,374",
            file=sys.stderr,
        )
    if abs(guo_proteins - 6556) / 6556 > 0.01:
        print(
            f"WARNING: guo_proteins={guo_proteins:,} deviates >1% from expected 6,556",
            file=sys.stderr,
        )
    if abs(guo_peptides - 40592) / 40592 > 0.05:
        print(
            f"WARNING: guo_peptides={guo_peptides:,} deviates >5% from expected 40,592",
            file=sys.stderr,
        )
    if abs(walzer_ea_genes - 2199) / 2199 > 0.02:
        print(
            f"WARN: Walzer EA gene count {walzer_ea_genes} differs from "
            f"expected 2,199 by >2%",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
