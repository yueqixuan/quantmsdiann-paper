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
    sample_df = df[sample_cols].replace("NA", pd.NA)
    sample_df = sample_df.replace("", pd.NA)
    quantified_mask = sample_df.notna().any(axis=1)
    quantified_df = df[quantified_mask]
    if unique_by is not None:
        return int(quantified_df[unique_by].nunique())
    return int(quantified_mask.sum())


def count_openswath_quantified(matrix_path: Path) -> tuple[int, int]:
    """Count quantified target precursors and protein groups in an OpenSWATH matrix."""
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
    protein_groups: set[str] = set()
    for chunk in pd.read_csv(
        matrix_path,
        sep="\t",
        usecols=usecols,
        dtype=str,
        chunksize=50000,
    ):
        # Exclude decoy rows
        is_decoy = chunk["Peptide"].str.startswith("DECOY_") | chunk["Protein"].str.contains(
            "DECOY", case=False, na=False
        )
        target = chunk[~is_decoy]
        # A row is quantified if at least one Intensity_ column is non-NA / non-empty
        intensity_df = target[intensity_cols].replace("NA", pd.NA).replace("", pd.NA)
        quantified_mask = intensity_df.notna().any(axis=1)
        quantified = target[quantified_mask]
        total_precursors += len(quantified)
        protein_groups.update(quantified["Protein"].tolist())
    return total_precursors, len(protein_groups)


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
