from __future__ import annotations

import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

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

# Guo et al. 2019 (doi:10.1016/j.isci.2019.10.059) Results, after DIA-expert
# manual curation. Reported in the paper text and matched in Walzer 2022
# Supplementary Table S2, "Original - after filter" column for PXD003539.
# Used for TSV context, NOT plotted as a headline bar.
GUO_CURATED_PEPTIDES = 22554
GUO_CURATED_PROTEINS = 3171

SUMMARY_LOG_PROTEIN_LINE_RE = re.compile(
    r"Protein groups with global q-value <= 0\.01:\s*(\d+)"
)


def download_if_missing(url: str, dest: Path, *, retries: int = 2) -> Path:
    """Download url to dest, skipping if dest already exists and is non-empty."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    last_exc: Exception | None = None
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


# Do NOT add the `main()` body in this task — Task 5 fills it in.
def main() -> int:  # pragma: no cover - filled in Task 5
    raise NotImplementedError("main is implemented in Task 5")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
