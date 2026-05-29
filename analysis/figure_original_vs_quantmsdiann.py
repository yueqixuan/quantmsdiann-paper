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

from analysis.contaminant_filter import is_target_protein_group

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "PXD003539"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "PXD003539"

PRIDE_QUANT_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/quant_tables"
)
PRIDE_SUBMISSION_BASE = (
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD003539"
)

PR_MATRIX_URL = f"{PRIDE_QUANT_BASE}/diann_report.pr_matrix.tsv"
SUMMARY_LOG_URL = f"{PRIDE_QUANT_BASE}/diannsummary.log"
DIANN_REPORT_PARQUET_URL = f"{PRIDE_QUANT_BASE}/diann_report.parquet"
DIANN_UNIQUE_GENES_MATRIX_URL = f"{PRIDE_QUANT_BASE}/diann_report.unique_genes_matrix.tsv"
OPENSWATH_MATRIX_URL = f"{PRIDE_SUBMISSION_BASE}/feature_alignment_requant_matrix.tsv"

HGNC_COMPLETE_SET_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/"
    "tsv/tsv/hgnc_complete_set.txt"
)

# quantms-collections SDRF for PXD003539 — runs the cell-line/disease metadata
# for the quantmsdiann side of the per-condition comparison.
QUANTMS_SDRF_URL = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/sdrf/"
    "PXD003539.sdrf.tsv"
)
# Expression Atlas E-PROT-73 experiment-design TSV gives Run -> cell line ->
# disease for the Walzer 2022 reanalysis side.
EA_EXPERIMENT_DESIGN_URL = (
    "https://www.ebi.ac.uk/gxa/experiments-content/E-PROT-73/resources/"
    "experiment-design"
)
# The Downloads HTML page embeds a JSON groupings blob with the
# g<N> -> cell-line and g<N> -> disease maps that don't appear anywhere else
# (E-PROT-73.tsv columns are 'g6.WithInSampleAbundance' style and the matrix
# alone gives no clue which g<N> is which cell line).
EA_DOWNLOADS_URL = "https://www.ebi.ac.uk/gxa/experiments/E-PROT-73/Downloads"

UNIQUE_GENES_METADATA_COLS = [
    "Genes", "N.Sequences", "N.Proteotypic.Sequences",
]

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
    # Headline target-only protein groups (post conservative contaminant
    # filter applied to pr_matrix.tsv unique Protein.Group rows).
    diann_proteins: int
    diann_precursors: int
    # Audit baseline: diannsummary.log "Protein groups with global
    # q-value <= 0.01" line (unfiltered). Default keeps backwards-
    # compatibility with callers that don't yet populate it.
    diann_proteins_unfiltered: int = 0


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


def per_run_real_detection_fraction_diann_parquet(
    parquet_path: Path,
    *,
    qvalue_cutoff: float = 0.01,
    global_qvalue_cutoff: float = 0.01,
) -> dict[str, float]:
    """Per-run fraction of distinct precursors confidently identified, read
    from the DIA-NN long-format report (`diann_report.parquet`).

    For each run, counts distinct Precursor.Id values whose per-run `Q.Value`
    is <= `qvalue_cutoff` AND whose `Global.Q.Value` is <= `global_qvalue_cutoff`.
    The denominator is the count of distinct Precursor.Ids that pass
    `Global.Q.Value <= global_qvalue_cutoff` anywhere in the report — the global
    1% FDR precursor pool, which is the apples-to-apples analogue of OpenSWATH's
    "all target rows in the requant matrix" denominator.

    Why use the parquet instead of the pr_matrix? The matrix is filtered at the
    cell level by --matrix-spec-q 0.05 (spectrum-level quant FDR), not by the
    per-run identification Q.Value. To match OpenSWATH's `score <= 0.01` per-run
    criterion strictly, we have to read Q.Value from the long-format report."""
    import pyarrow.parquet as pq
    cols = ["Run", "Precursor.Id", "Q.Value", "Global.Q.Value"]
    pf = pq.ParquetFile(str(parquet_path))
    schema_names = pf.schema_arrow.names
    missing = [c for c in cols if c not in schema_names]
    if missing:
        raise ValueError(
            f"DIA-NN parquet missing expected columns: {missing}"
        )
    global_precursors: set[str] = set()
    per_run_detected: dict[str, set[str]] = {}
    for batch in pf.iter_batches(columns=cols, batch_size=200_000):
        runs = batch.column("Run").to_pylist()
        pids = batch.column("Precursor.Id").to_pylist()
        qvs = batch.column("Q.Value").to_pylist()
        gqvs = batch.column("Global.Q.Value").to_pylist()
        for r, p, q, g in zip(runs, pids, qvs, gqvs):
            per_run_detected.setdefault(r, set())  # register every run seen
            if g is None or g > global_qvalue_cutoff:
                continue
            global_precursors.add(p)
            if q is None or q > qvalue_cutoff:
                continue
            per_run_detected[r].add(p)
    denom = len(global_precursors)
    if denom == 0:
        return {r: 0.0 for r in per_run_detected}
    return {r: len(s) / denom for r, s in per_run_detected.items()}


def per_run_real_detection_fraction_openswath(
    matrix_path: Path,
    *,
    qvalue_cutoff: float = 0.01,
) -> dict[str, float]:
    """For each `score_<run>` column in the OpenSWATH requantification matrix,
    return the fraction of target (non-decoy) rows whose score for that run
    is below `qvalue_cutoff`.

    The OpenSWATH `feature_alignment_requant_matrix.tsv` is a *requantified*
    matrix: Intensity values are filled in across all runs even when a
    precursor was not directly identified in that run. The corresponding
    `score_<run>` column (pyprophet m_score / q-value) is the truth signal —
    real detections have score <= 0.01, requantified placeholders are tagged
    score = 2.0. Counting non-NA Intensity values is therefore not a real
    measure of completeness; we use the score column instead."""
    header_df = pd.read_csv(matrix_path, sep="\t", nrows=0)
    cols = list(header_df.columns)
    if "Peptide" not in cols or "Protein" not in cols:
        raise ValueError("OpenSWATH matrix missing Peptide/Protein columns")
    score_cols = [c for c in cols if c.startswith("score_")]
    if not score_cols:
        raise ValueError("OpenSWATH matrix has no score_* columns")

    target_total = 0
    per_run_hits: dict[str, int] = {c: 0 for c in score_cols}
    usecols = ["Peptide", "Protein"] + score_cols
    for chunk in pd.read_csv(matrix_path, sep="\t", dtype=str,
                             usecols=usecols, chunksize=50000):
        is_decoy = (
            chunk["Peptide"].str.startswith("DECOY_", na=False)
            | chunk["Protein"].str.contains("DECOY", case=False, na=False)
        )
        targets = chunk[~is_decoy].copy()
        target_total += len(targets)
        for col in score_cols:
            scores = pd.to_numeric(targets[col], errors="coerce")
            per_run_hits[col] += int((scores <= qvalue_cutoff).sum())
    if target_total == 0:
        return {c: 0.0 for c in score_cols}
    return {c: per_run_hits[c] / target_total for c in score_cols}


def unique_peptides_per_protein_diann(matrix_path: Path) -> dict[str, int]:
    """Per Protein.Group, count distinct Stripped.Sequence values from
    proteotypic precursor rows that have at least one non-NA sample value.

    Restricting to Proteotypic == 1 ensures we only count peptides that uniquely
    identify the protein (the natural definition of 'unique peptides per
    protein'). Multiple charge states / modforms of the same peptide collapse
    to a single Stripped.Sequence entry."""
    # Stream in chunks: the ProCan pr_matrix is ~2 GB and loading it whole
    # exhausts memory. We accumulate, per Protein.Group, the set of distinct
    # proteotypic Stripped.Sequence values seen in any quantified row.
    header = pd.read_csv(matrix_path, sep="\t", nrows=0)
    missing = [c for c in PR_METADATA_COLS if c not in header.columns]
    if missing:
        raise ValueError(
            f"DIA-NN matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in header.columns if c not in PR_METADATA_COLS]
    pg_to_peps: dict[str, set[str]] = {}
    for chunk in pd.read_csv(matrix_path, sep="\t", dtype=str, chunksize=100_000):
        proteotypic = chunk[chunk["Proteotypic"] == "1"]
        if proteotypic.empty:
            continue
        quantified = proteotypic[proteotypic[sample_cols].notna().any(axis=1)]
        for pg, seq in zip(quantified["Protein.Group"], quantified["Stripped.Sequence"]):
            if is_target_protein_group(pg):
                pg_to_peps.setdefault(pg, set()).add(seq)
    return {pg: len(s) for pg, s in pg_to_peps.items()}


def unique_peptides_per_protein_openswath(
    matrix_path: Path,
    *,
    qvalue_cutoff: float = 0.01,
) -> dict[str, int]:
    """Per OpenSWATH Protein, count distinct stripped peptide sequences that
    were confidently detected (score <= `qvalue_cutoff`) in at least one run.

    Restricted to proteotypic rows — those whose Protein column starts with
    `1/` (peptide maps to exactly one protein). Decoys are excluded. The
    confidence filter prevents requant-only placeholder rows (where every score
    is 2.0) from inflating the per-protein peptide count."""
    header_df = pd.read_csv(matrix_path, sep="\t", nrows=0)
    cols = list(header_df.columns)
    if "Peptide" not in cols or "Protein" not in cols:
        raise ValueError("OpenSWATH matrix missing Peptide/Protein columns")
    score_cols = [c for c in cols if c.startswith("score_")]
    if not score_cols:
        raise ValueError("OpenSWATH matrix has no score_* columns")
    usecols = ["Peptide", "Protein"] + score_cols
    protein_to_peptides: dict[str, set[str]] = {}
    for chunk in pd.read_csv(matrix_path, sep="\t", dtype=str,
                             usecols=usecols, chunksize=50000):
        is_decoy = (
            chunk["Peptide"].str.startswith("DECOY_", na=False)
            | chunk["Protein"].str.contains("DECOY", case=False, na=False)
        )
        targets = chunk[~is_decoy]
        proteotypic = targets[targets["Protein"].str.startswith("1/", na=False)]
        if proteotypic.empty:
            continue
        score_block = proteotypic[score_cols].apply(
            pd.to_numeric, errors="coerce"
        )
        any_detection = (score_block <= qvalue_cutoff).any(axis=1)
        detected = proteotypic[any_detection]
        for prot, pep_id in zip(detected["Protein"], detected["Peptide"]):
            stripped = _stripped_peptide(pep_id)
            if stripped is None:
                continue
            protein_to_peptides.setdefault(prot, set()).add(stripped)
    return {p: len(s) for p, s in protein_to_peptides.items()}


def proteins_with_min_peptides(counts: dict[str, int], min_k: int) -> int:
    """Number of proteins with at least `min_k` unique peptides."""
    return sum(1 for n in counts.values() if n >= min_k)


def count_quantified_genes_diann(matrix_path: Path) -> int:
    """Count rows in DIA-NN's unique_genes_matrix.tsv with >=1 non-NA sample
    cell. Metadata columns are UNIQUE_GENES_METADATA_COLS (Genes, N.Sequences,
    N.Proteotypic.Sequences); the remaining columns are per-run intensities."""
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    missing = [c for c in UNIQUE_GENES_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"unique_genes_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in UNIQUE_GENES_METADATA_COLS]
    if not sample_cols:
        raise ValueError("unique_genes_matrix has no sample columns")
    return int(df[sample_cols].notna().any(axis=1).sum())


def load_hgnc_symbol_to_ensembl(hgnc_tsv_path: Path) -> dict[str, str]:
    """Parse HGNC's `hgnc_complete_set.txt` and return a mapping from any known
    gene-symbol form (current `symbol`, pipe-separated `alias_symbol` and
    `prev_symbol`) to its `ensembl_gene_id`. Rows without an ensembl_gene_id
    are skipped. Conflicts are resolved last-writer-wins (which only matters
    when a single symbol historically pointed to multiple genes — rare)."""
    df = pd.read_csv(
        hgnc_tsv_path,
        sep="\t",
        dtype=str,
        usecols=["symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id"],
    )
    df = df[df["ensembl_gene_id"].fillna("").str.startswith("ENSG")]
    mapping: dict[str, str] = {}
    for sym, aliases, prevs, ensg in zip(
        df["symbol"], df["alias_symbol"], df["prev_symbol"], df["ensembl_gene_id"],
    ):
        if isinstance(sym, str) and sym:
            mapping[sym] = ensg
        if isinstance(aliases, str) and aliases:
            for a in aliases.split("|"):
                a = a.strip()
                if a:
                    mapping[a] = ensg
        if isinstance(prevs, str) and prevs:
            for p in prevs.split("|"):
                p = p.strip()
                if p:
                    mapping[p] = ensg
    return mapping


def load_walzer_genes_ensembl(eprot73_path: Path) -> set[str]:
    """Return the set of unique Ensembl gene IDs from the Walzer E-PROT-73
    Expression Atlas file. Same filter as count_eprot73_genes (rows whose
    `Gene ID` column starts with `ENSG`)."""
    df = pd.read_csv(eprot73_path, sep="\t", dtype=str, usecols=[0])
    df.columns = ["gene_id"]
    df = df[df["gene_id"].fillna("").str.startswith("ENSG")]
    return set(df["gene_id"].unique())


def quantmsdiann_genes_as_ensembl(
    matrix_path: Path,
    symbol_to_ensembl: dict[str, str],
    *,
    min_detection_fraction: float = 0.0,
) -> tuple[set[str], int]:
    """For each quantified row in the DIA-NN unique_genes_matrix, look up the
    `Genes` column in `symbol_to_ensembl` and collect the Ensembl gene IDs.
    Returns (mapped_ensg_set, unmapped_count).

    `min_detection_fraction` controls how many sample cells must be non-NA for
    a row to count as 'detected':
      - 0.0 (default): any row with >=1 non-NA cell (the bare identification set).
      - 0.5: row must be non-NA in >=ceil(0.5 * n_samples) runs (mimics Walzer's
        50%-per-group consistency filter applied globally across all runs).

    The Genes column can hold multiple symbols separated by `;` (gene-group
    case). We split on `;` and try each; counts as unmapped if none of the
    components are in the mapping."""
    import math
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    missing = [c for c in UNIQUE_GENES_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"unique_genes_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in UNIQUE_GENES_METADATA_COLS]
    non_na_count = df[sample_cols].notna().sum(axis=1)
    min_required = max(1, math.ceil(min_detection_fraction * len(sample_cols)))
    detected = df[non_na_count >= min_required]
    mapped: set[str] = set()
    unmapped = 0
    for raw in detected["Genes"].fillna(""):
        if not raw:
            unmapped += 1
            continue
        ensgs = {
            symbol_to_ensembl[s]
            for s in (sym.strip() for sym in raw.split(";"))
            if s and s in symbol_to_ensembl
        }
        if ensgs:
            mapped.update(ensgs)
        else:
            unmapped += 1
    return mapped, unmapped


_CELL_LINE_NCI_PREFIX_RE = re.compile(r"^NCI[-_/]", flags=re.IGNORECASE)
_CELL_LINE_NONALNUM_RE = re.compile(r"[^A-Za-z0-9]")


def normalize_cell_line(name: str | None) -> str:
    """Normalise a cell-line name so the quantms SDRF and E-PROT-73 spellings
    collide. Strips a leading 'NCI-' / 'NCI_' / 'NCI/' prefix, removes all
    other non-alphanumeric characters, and uppercases.

    Example: 'CCRF-CEM' -> 'CCRFCEM', 'NCI-H226' -> 'H226',
    'Hs-578-T' -> 'HS578T'."""
    if not name:
        return ""
    name = _CELL_LINE_NCI_PREFIX_RE.sub("", name)
    name = _CELL_LINE_NONALNUM_RE.sub("", name)
    return name.upper()


def parse_eprot73_groupings(
    downloads_html_path: Path,
) -> tuple[dict[str, str], dict[str, str]]:
    """Parse the Expression Atlas E-PROT-73 Downloads HTML page and return
    (g_to_cell_line, g_to_disease) mappings.

    The Downloads page embeds an inline JSON object (`content: {...}`) with
    `tabs[0].props.groups` listing primary groupings. The CELL_LINE grouping
    is `[[cell_line_name, [g1, g2, ...]], ...]`; DISEASE follows the same
    shape. We pull both out with a small balanced-bracket scanner so we don't
    depend on the surrounding HTML structure."""
    import json
    html = downloads_html_path.read_text()

    def extract(name: str) -> list[list]:
        m = re.search(
            rf'"name":"{re.escape(name)}",[^{{]*?"groupings":', html,
        )
        if m is None:
            raise ValueError(
                f"E-PROT-73 Downloads HTML missing '{name}' grouping"
            )
        start = m.end()
        if html[start] != "[":
            raise ValueError(
                f"Unexpected character after '{name}' groupings: {html[start]!r}"
            )
        depth = 0
        for i in range(start, len(html)):
            if html[i] == "[":
                depth += 1
            elif html[i] == "]":
                depth -= 1
                if depth == 0:
                    return json.loads(html[start : i + 1])
        raise ValueError(f"Unbalanced brackets in '{name}' grouping")

    cl_groupings = extract("CELL_LINE")
    ds_groupings = extract("DISEASE")
    g_to_cl: dict[str, str] = {}
    for name, gs in cl_groupings:
        for g in gs:
            g_to_cl[g] = name
    g_to_ds: dict[str, str] = {}
    for name, gs in ds_groupings:
        for g in gs:
            g_to_ds[g] = name
    return g_to_cl, g_to_ds


def load_sdrf_data_file_to_cell_line(sdrf_path: Path) -> dict[str, str]:
    """Parse the PXD003539 quantms-collections SDRF and return a mapping from
    DIA-NN matrix column name (i.e. the `comment[data file]` field rewritten
    to a `.mzML` extension) to the row's `characteristics[cell line]`.

    The SDRF stores `*.wiff` filenames but DIA-NN's pr_matrix /
    unique_genes_matrix columns are `*.mzML`. We rewrite the extension so
    downstream callers can look up DIA-NN columns directly."""
    df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    needed = ["characteristics[cell line]", "comment[data file]"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"SDRF missing required columns: {missing}")
    out: dict[str, str] = {}
    for cell_line, data_file in zip(
        df["characteristics[cell line]"], df["comment[data file]"],
    ):
        if not isinstance(data_file, str) or not data_file:
            continue
        mzml = re.sub(r"\.wiff$", ".mzML", data_file)
        out[mzml] = cell_line
    return out


def load_ea_cell_line_to_disease(experiment_design_path: Path) -> dict[str, str]:
    """Parse the E-PROT-73 experiment-design TSV and return a mapping from
    normalised cell-line name (see `normalize_cell_line`) to disease label.

    The TSV has one row per MS Run; every run for a given cell line lists the
    same disease, so we collapse to the unique mapping. Conflicting diseases
    for the same cell line would be unexpected and we raise on collision."""
    df = pd.read_csv(experiment_design_path, sep="\t", dtype=str)
    cell_col = "Sample Characteristic[cell line]"
    disease_col = "Sample Characteristic[disease]"
    missing = [c for c in (cell_col, disease_col) if c not in df.columns]
    if missing:
        raise ValueError(
            f"E-PROT-73 experiment-design missing columns: {missing}"
        )
    out: dict[str, str] = {}
    for cell, disease in zip(df[cell_col], df[disease_col]):
        if not isinstance(cell, str) or not isinstance(disease, str):
            continue
        key = normalize_cell_line(cell)
        if not key:
            continue
        if key in out and out[key] != disease:
            raise ValueError(
                f"Conflicting disease for cell line {cell!r}: "
                f"{out[key]!r} vs {disease!r}"
            )
        out[key] = disease
    return out


def walzer_genes_per_condition(
    eprot73_path: Path,
    downloads_html_path: Path,
) -> dict[str, set[str]]:
    """For each E-PROT-73 disease, return the set of Ensembl gene IDs detected
    in at least one g<N> column belonging to that disease.

    "Detected" means abundance > 0 in the `g<N>.WithInSampleAbundance` column.
    Rows whose `Gene ID` doesn't start with `ENSG` are ignored (Expression
    Atlas occasionally emits non-ENSG summary rows like 'totalGenes')."""
    _, g_to_disease = parse_eprot73_groupings(downloads_html_path)
    df = pd.read_csv(eprot73_path, sep="\t", dtype=str)
    if df.columns[0] != "Gene ID":
        raise ValueError(
            f"E-PROT-73 TSV first column should be 'Gene ID', got "
            f"{df.columns[0]!r}"
        )
    df = df[df["Gene ID"].fillna("").str.startswith("ENSG")].copy()

    # Map each sample column ('g6.WithInSampleAbundance') to its g<N>.
    sample_cols = [c for c in df.columns if ".WithInSampleAbundance" in c]
    col_to_g = {c: c.split(".", 1)[0] for c in sample_cols}

    out: dict[str, set[str]] = {}
    for col in sample_cols:
        g = col_to_g[col]
        disease = g_to_disease.get(g)
        if disease is None:
            continue
        vals = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        detected = df.loc[vals > 0, "Gene ID"]
        out.setdefault(disease, set()).update(detected)
    return out


def quantmsdiann_genes_per_condition(
    unique_genes_matrix_path: Path,
    sdrf_path: Path,
    ea_design_path: Path,
    symbol_to_ensembl: dict[str, str],
    *,
    min_detection_fraction_per_cell_line: float = 0.0,
    min_global_detection_fraction: float = 0.0,
) -> dict[str, set[str]]:
    """For each disease (using the E-PROT-73 9-cancer-type axis), return the
    set of Ensembl gene IDs detected by quantmsdiann in at least one cell
    line belonging to that disease.

    Filtering knobs (apply jointly — a gene must pass both):
      - `min_global_detection_fraction`: gene must be non-NA in
        >= ceil(f * n_runs_total) of all 120 DIA-NN runs. This is the same
        global filter used in `quantmsdiann_genes_as_ensembl` and the
        `supp_walzer_vs_quantms_genes_ensembl` Venn (set to 0.5 there).
      - `min_detection_fraction_per_cell_line`: gene must be non-NA in
        >= ceil(f * n_runs_of_cell_line) of a given cell line's runs to
        count as "detected for that cell line". This mirrors Walzer's
        '50% per group' consistency filter at the replicate-group level
        (Walzer 2022 §Methods), which is what produces the
        g<N>.WithInSampleAbundance values in E-PROT-73. With only 2 runs
        per cell line in this dataset, 0.5 here is essentially a no-op
        (ceil(0.5*2)=1) — the user-facing filter is the global one above.

    The SDRF maps DIA-NN columns to cell line (dashed quantms spelling);
    `load_ea_cell_line_to_disease` maps the normalised cell-line name to one
    of the 9 E-PROT-73 disease labels. Genes column entries are split on `;`
    and each component looked up in the HGNC symbol -> Ensembl mapping;
    unmappable rows are skipped (consistent with `quantmsdiann_genes_as_ensembl`)."""
    import math
    df = pd.read_csv(unique_genes_matrix_path, sep="\t", dtype=str)
    missing = [c for c in UNIQUE_GENES_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"unique_genes_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in UNIQUE_GENES_METADATA_COLS]

    sdrf_run_to_cell = load_sdrf_data_file_to_cell_line(sdrf_path)
    ea_cell_to_disease = load_ea_cell_line_to_disease(ea_design_path)

    # Group sample columns by (cell line, disease).
    cell_to_cols: dict[str, list[str]] = {}
    cell_to_disease: dict[str, str] = {}
    for col in sample_cols:
        cell = sdrf_run_to_cell.get(col)
        if not cell:
            continue
        disease = ea_cell_to_disease.get(normalize_cell_line(cell))
        if disease is None:
            continue
        cell_to_cols.setdefault(cell, []).append(col)
        cell_to_disease[cell] = disease

    # Precompute per-row ENSG list once.
    row_ensgs: list[set[str]] = []
    for raw in df["Genes"].fillna(""):
        if not raw:
            row_ensgs.append(set())
            continue
        ensgs = {
            symbol_to_ensembl[s]
            for s in (sym.strip() for sym in raw.split(";"))
            if s and s in symbol_to_ensembl
        }
        row_ensgs.append(ensgs)

    # Global filter: gene must pass it to be eligible in any condition.
    global_min = max(
        1, math.ceil(min_global_detection_fraction * len(sample_cols))
    )
    global_non_na = df[sample_cols].notna().sum(axis=1)
    global_pass = global_non_na >= global_min

    out: dict[str, set[str]] = {}
    for cell, cols in cell_to_cols.items():
        non_na = df[cols].notna().sum(axis=1)
        min_required = max(
            1,
            math.ceil(min_detection_fraction_per_cell_line * len(cols)),
        )
        detected_mask = (non_na >= min_required) & global_pass
        disease = cell_to_disease[cell]
        bucket = out.setdefault(disease, set())
        for ensgs, ok in zip(row_ensgs, detected_mask):
            if ok:
                bucket.update(ensgs)
    return out


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


def count_target_protein_groups_pr_matrix(pr_matrix_path: Path) -> tuple[int, int]:
    """Return `(unfiltered_unique_pg, target_unique_pg)` from a DIA-NN
    `pr_matrix.tsv`. PXD003539 has no `pg_matrix.tsv` on disk, so the
    target-only protein-group headline must be derived from the
    pr_matrix's unique `Protein.Group` strings.

    Unfiltered = the count of distinct Protein.Group strings among rows
    quantified in at least one run (matches the historical baseline).
    Target = the same set, restricted to Protein.Group strings whose
    every semicolon-separated token has no
    CONTAM_/Cont_/ENTRAP_/DECOY_/decoy_ prefix per the 2026-05-21
    conservative-filter spec.
    """
    df = pd.read_csv(pr_matrix_path, sep="\t", dtype=str)
    missing = [c for c in PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"PXD003539 pr_matrix missing metadata columns: {missing}"
        )
    sample_cols = [c for c in df.columns if c not in PR_METADATA_COLS]
    quantified = df[df[sample_cols].notna().any(axis=1)]
    pg_unique = set(quantified["Protein.Group"].dropna().unique().tolist())
    target_unique = {pg for pg in pg_unique if is_target_protein_group(pg)}
    return (len(pg_unique), len(target_unique))


def render_figure(
    counts: Counts,
    svg_path: Path,
) -> None:
    """Render a 2-condition x 2-metric grouped bar chart. Paper-ready: only
    bars, value labels, axis labels, and legend — no title, no footer."""
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

    ax.legend(loc="upper right", frameon=False)

    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


def render_missing_values_per_run(
    guo_per_run: dict[str, float],
    diann_per_run: dict[str, float],
    svg_path: Path,
) -> None:
    """Plot per-run non-NA fraction for both pipelines. Paper-ready: no title,
    no footer; only lines, axis labels, and legend.

    Runs on the x axis are aligned by the L<date>_<n>_SW stem extracted from
    the column name. Bars sorted by the stem (which sorts chronologically by
    acquisition date)."""
    import re as _re

    def stem(col: str) -> str:
        # diann columns look like 'guot_L130610_003_SW.mzML'
        # openswath score columns look like 'score_guot_L130610_003_SW_with_dscore.csv_0_42'
        # The regex r"L\d+_\d+_SW" matches the shared acquisition-date stem in
        # both formats, so alignment between the two pipelines is unambiguous.
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
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


_DISEASE_LABEL_ORDER = [
    "leukemia", "central nervous system cancer", "breast cancer",
    "colorectal cancer", "lung cancer", "melanoma", "ovarian cancer",
    "prostate cancer", "renal cancer",
]


def render_genes_per_condition(
    walzer_per_cond: dict[str, set[str]],
    diann_per_cond: dict[str, set[str]],
    svg_path: Path,
) -> None:
    """Grouped bar chart: per disease condition, two bars (Walzer 2022
    E-PROT-73 vs quantmsdiann) showing the number of distinct Ensembl gene IDs
    detected. The disease axis is the 9 NCI-60 cancer types defined by
    E-PROT-73's primary DISEASE grouping; both pipelines use this same axis
    after cell-line normalisation. Paper-ready: no title, no footer."""
    conditions = sorted(
        set(walzer_per_cond) | set(diann_per_cond),
        key=lambda d: (
            _DISEASE_LABEL_ORDER.index(d)
            if d in _DISEASE_LABEL_ORDER
            else len(_DISEASE_LABEL_ORDER)
        ),
    )

    walzer_vals = [len(walzer_per_cond.get(c, set())) for c in conditions]
    diann_vals = [len(diann_per_cond.get(c, set())) for c in conditions]

    fig, ax = plt.subplots(figsize=(11, 6))
    x = list(range(len(conditions)))
    bar_width = 0.4
    bars_w = ax.bar([xi - bar_width / 2 for xi in x], walzer_vals,
                    width=bar_width, color="#90caf9",
                    label="Walzer 2022 (E-PROT-73)")
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals,
                    width=bar_width, color="#1f77b4",
                    label="quantmsdiann (DIA-NN)")
    for bars, vals in ((bars_w, walzer_vals), (bars_d, diann_vals)):
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{v:,}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    # Shorten 'central nervous system cancer' for display only.
    display = [c.replace("central nervous system cancer", "CNS cancer")
               for c in conditions]
    ax.set_xticklabels(display, rotation=25, ha="right")
    ax.set_ylabel("Distinct Ensembl gene IDs detected")
    ymax = max(max(walzer_vals, default=0), max(diann_vals, default=0))
    ax.set_ylim(0, ymax * 1.15 if ymax else 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=2, frameon=False)

    fig.tight_layout(rect=(0, 0.12, 1, 1))

    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)


def render_peptides_per_protein(
    guo_peptide_counts: dict[str, int],
    diann_peptide_counts: dict[str, int],
    svg_path: Path,
    thresholds: tuple[int, ...] = (1, 2, 3, 5, 10),
) -> None:
    """Grouped bar chart: number of protein groups with >=k unique peptides,
    for each k in `thresholds`, comparing Guo (OpenSWATH) vs quantmsdiann
    (DIA-NN). Paper-ready: no title, no footer."""
    guo_values = [proteins_with_min_peptides(guo_peptide_counts, k)
                  for k in thresholds]
    diann_values = [proteins_with_min_peptides(diann_peptide_counts, k)
                    for k in thresholds]

    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.38
    x = list(range(len(thresholds)))
    bars_guo = ax.bar([xi - width / 2 for xi in x], guo_values, width=width,
                      color="#9e9e9e", label="Guo 2019 (OpenSWATH)")
    bars_diann = ax.bar([xi + width / 2 for xi in x], diann_values, width=width,
                        color="#1f77b4", label="quantmsdiann (DIA-NN)")
    for bars, vals in [(bars_guo, guo_values), (bars_diann, diann_values)]:
        for bar, v in zip(bars, vals):
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height(),
                f"{v:,}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([f"≥ {k}" for k in thresholds])
    ax.set_xlabel("Minimum unique peptides per protein group")
    ax.set_ylabel("Protein groups (1% FDR)")
    ymax = max(max(guo_values), max(diann_values)) * 1.14
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", frameon=False)

    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
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
        ("Protein groups", "quantmsdiann (DIA-NN, 1% FDR, target-only)",
         counts.diann_proteins,
         "pr_matrix.tsv unique Protein.Group rows (quantified, target-only); "
         "post conservative contaminant filter "
         "(CONTAM_/Cont_/ENTRAP_/DECOY_/decoy_ token) per 2026-05-21 spec"),
        ("Protein groups", "quantmsdiann (DIA-NN, 1% FDR, diannsummary.log)",
         counts.diann_proteins_unfiltered,
         "audit baseline: diannsummary.log 'Protein groups with "
         "global q-value <= 0.01' line (unfiltered)"),
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
    parquet_path = DATA_DIR / "diann_report.parquet"
    diann_genes_path = DATA_DIR / "diann_report.unique_genes_matrix.tsv"

    print("Downloading files (skipped if already cached)...")
    download_if_missing(PR_MATRIX_URL, pr_path)
    download_if_missing(SUMMARY_LOG_URL, log_path)
    download_if_missing(OPENSWATH_MATRIX_URL, opensw_path)
    download_if_missing(DIANN_REPORT_PARQUET_URL, parquet_path)
    download_if_missing(DIANN_UNIQUE_GENES_MATRIX_URL, diann_genes_path)
    eprot73_path = download_if_missing(EPROT73_URL, DATA_DIR / "E-PROT-73.tsv")

    # 3. Compute counts
    print("Computing quantmsdiann counts...")
    diann_peptides = count_quantified_rows(
        pr_path, PR_METADATA_COLS, unique_by="Stripped.Sequence"
    )
    diann_precursors = count_quantified_rows(pr_path, PR_METADATA_COLS)
    # Apply conservative contaminant filter at the Protein.Group level:
    # headline = unique target-only PGs in pr_matrix; the log value is
    # kept as the unfiltered audit baseline.
    diann_proteins_log = parse_summary_log(log_path)
    pg_unf, pg_target = count_target_protein_groups_pr_matrix(pr_path)
    print(f"  pr_matrix unique Protein.Group: unfiltered={pg_unf:,}  "
          f"target_only={pg_target:,} (delta {pg_unf - pg_target:,})")
    diann_proteins = pg_target

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
        diann_proteins_unfiltered=diann_proteins_log,
    )

    # 4. Render figure
    svg_path = FIGURES_DIR / "main_comparison.svg"
    render_figure(counts, svg_path)
    print(f"Figure saved to {svg_path}")

    # 5. Write auditable TSV
    data_dir = FIGURES_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = data_dir / "counts.tsv"
    write_counts_tsv(counts, tsv_path)
    print(f"Counts TSV saved to {tsv_path}")

    # 6. Supplementary figure B: per-run completeness (strict 1% per-run FDR)
    print("Computing per-run completeness for supplementary figure B...")
    guo_per_run = per_run_real_detection_fraction_openswath(opensw_path)
    diann_per_run = per_run_real_detection_fraction_diann_parquet(parquet_path)
    suppB_svg = FIGURES_DIR / "supp_missing_values_per_run.svg"
    render_missing_values_per_run(
        guo_per_run, diann_per_run, suppB_svg,
    )
    print(
        f"Supplementary B figure saved to {suppB_svg}"
    )

    # Also print median per-run fraction for each pipeline.
    import statistics as _stats
    print(
        f"Median per-run real-detection fraction (1% per-run FDR): "
        f"Guo={_stats.median(guo_per_run.values()):.2%} "
        f"quantmsdiann={_stats.median(diann_per_run.values()):.2%}"
    )

    # 6b. Supplementary figure C: peptides-per-protein depth comparison
    print("Computing unique peptides per protein for supplementary figure C...")
    guo_pep_per_prot = unique_peptides_per_protein_openswath(opensw_path)
    diann_pep_per_prot = unique_peptides_per_protein_diann(pr_path)
    suppC_svg = FIGURES_DIR / "supp_peptides_per_protein.svg"
    render_peptides_per_protein(
        guo_pep_per_prot, diann_pep_per_prot, suppC_svg,
    )
    print(
        f"Supplementary C figure saved to {suppC_svg}"
    )
    print(
        f"Proteins with >=2 unique peptides: "
        f"Guo={proteins_with_min_peptides(guo_pep_per_prot, 2):,}  "
        f"quantmsdiann={proteins_with_min_peptides(diann_pep_per_prot, 2):,}"
    )

    # 6c. Supplementary figure D: gene-level Venn at the Ensembl gene ID level
    # (Walzer 2022 vs quantmsdiann). Walzer publishes ENSGs in E-PROT-73;
    # DIA-NN's unique_genes_matrix has HGNC symbols. We map symbols->ENSG
    # via HGNC's symbol + alias_symbol + prev_symbol -> ensembl_gene_id table.
    print("Loading HGNC symbol->Ensembl mapping...")
    hgnc_path = download_if_missing(
        HGNC_COMPLETE_SET_URL, DATA_DIR / "hgnc_complete_set.txt",
    )
    symbol_to_ensg = load_hgnc_symbol_to_ensembl(hgnc_path)
    walzer_ensg = load_walzer_genes_ensembl(eprot73_path)
    # Apply a 50%-of-runs detection filter to DIA-NN, mimicking Walzer's
    # post-processing (Walzer's '50%-per-group filter' is per cell-line
    # replicate group; without the SDRF group map we apply it globally — the
    # closest comparable strictness we can produce).
    diann_ensg, diann_unmapped = quantmsdiann_genes_as_ensembl(
        diann_genes_path, symbol_to_ensg, min_detection_fraction=0.5,
    )
    suppD_svg = FIGURES_DIR / "supp_walzer_vs_quantms_genes_ensembl.svg"
    from analysis.venn_protein_accessions import render_venn_diagram
    render_venn_diagram(
        walzer_ensg, diann_ensg, suppD_svg,
        left_label="Walzer 2022\n(E-PROT-73)",
        right_label="quantmsdiann\n(DIA-NN, $\\geq$50% of runs)",
        left_color="#90caf9",
        right_color="#1f77b4",
    )
    inter_ensg = walzer_ensg & diann_ensg
    print(
        f"Supplementary D figure saved to {suppD_svg}"
    )
    print(
        f"Ensembl gene IDs (DIA-NN >=50% of runs): Walzer={len(walzer_ensg):,}  "
        f"quantmsdiann={len(diann_ensg):,}  "
        f"intersection={len(inter_ensg):,}  Walzer-only="
        f"{len(walzer_ensg - diann_ensg):,}  quantmsdiann-only="
        f"{len(diann_ensg - walzer_ensg):,}  (unmapped HGNC={diann_unmapped:,})"
    )

    # 6c-bis. Supplementary figure F: gene-level detections per cancer
    # condition (Walzer 2022 E-PROT-73 vs quantmsdiann). The 9 NCI-60 cancer
    # types come from E-PROT-73's primary DISEASE grouping; both pipelines
    # are mapped to the same axis via cell-line normalisation.
    print("Computing per-condition gene detections...")
    sdrf_path = download_if_missing(
        QUANTMS_SDRF_URL, DATA_DIR / "PXD003539.sdrf.tsv",
    )
    ea_design_path = download_if_missing(
        EA_EXPERIMENT_DESIGN_URL,
        DATA_DIR / "E-PROT-73-experiment-design.tsv",
    )
    ea_downloads_path = download_if_missing(
        EA_DOWNLOADS_URL, DATA_DIR / "E-PROT-73-downloads.html",
    )
    walzer_per_cond = walzer_genes_per_condition(eprot73_path, ea_downloads_path)
    # Pre-filter with the same global 50%-of-runs threshold as the global
    # Venn (`supp_walzer_vs_quantms_genes_ensembl`). With only 2 runs per
    # cell line in this dataset, a per-cell-line 50% filter is a no-op
    # (ceil(0.5*2)=1), so the global threshold is what actually constrains
    # the set to genes that Walzer's per-group consistency regime would also
    # keep at a comparable strictness.
    diann_per_cond = quantmsdiann_genes_per_condition(
        diann_genes_path, sdrf_path, ea_design_path, symbol_to_ensg,
        min_global_detection_fraction=0.5,
    )
    suppF_svg = FIGURES_DIR / "supp_genes_per_condition.svg"
    render_genes_per_condition(
        walzer_per_cond, diann_per_cond, suppF_svg,
    )
    print(
        f"Supplementary F figure saved to {suppF_svg}"
    )
    print("Per-condition gene detections (Walzer | quantmsdiann):")
    for cond in sorted(set(walzer_per_cond) | set(diann_per_cond)):
        w = len(walzer_per_cond.get(cond, set()))
        d = len(diann_per_cond.get(cond, set()))
        print(f"  {cond:32s} {w:>6,} | {d:>6,}")

    # 6d. Supplementary figure E: Venn of UniProt accessions with >=2 unique
    # peptides (Guo vs quantmsdiann).
    print("Computing Venn of protein accessions with >=2 unique peptides...")
    from analysis.venn_protein_accessions import (
        accessions_with_min_peptides_diann,
        accessions_with_min_peptides_openswath,
    )
    guo_acc = accessions_with_min_peptides_openswath(opensw_path, min_peptides=2)
    diann_acc = accessions_with_min_peptides_diann(pr_path, min_peptides=2)
    suppE_svg = FIGURES_DIR / "supp_venn_protein_accessions.svg"
    render_venn_diagram(guo_acc, diann_acc, suppE_svg)
    inter = guo_acc & diann_acc
    print(
        f"Supplementary E figure saved to {suppE_svg}"
    )
    print(
        f"Accessions (>=2 unique peptides): "
        f"Guo={len(guo_acc):,}  quantmsdiann={len(diann_acc):,}  "
        f"intersection={len(inter):,}  Guo-only={len(guo_acc - diann_acc):,}  "
        f"quantmsdiann-only={len(diann_acc - guo_acc):,}"
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
    # diannsummary.log unfiltered baseline (pinned by the historical
    # spec) — the post-filter target headline is lower.
    if diann_proteins_log != 6927:
        print(
            f"WARNING: diann_proteins (log, unfiltered)={diann_proteins_log:,} "
            f"!= expected 6,927",
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
