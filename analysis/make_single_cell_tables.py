#!/usr/bin/env python
"""Generate the Fig. 3 single-cell input tables from the PUBLIC DIA-NN reports.

This is the reproducibility generator for ``figure_single_cell_combined.py``:
every number in Fig. 3 is derived here from the deposited DIA-NN reports, so
nothing in the figure is hand-entered.

PROVENANCE
==========
Our reanalyses (this pipeline) are published on the PRIDE FTP:

  https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/single-cell/<PXD>/<version>/quant_tables/diann_report.{parquet,tsv}

    PXD046357  HeLa Astral single-cell (Orbitrap Astral) -> "HeLa Astral SC"
    PXD044991  HeLa One-Tip                               -> "HeLa One-Tip"

  DIA-NN versions: ``v1_8_1`` (diann_report.tsv) and ``v2_5_1_enterprise``
  (diann_report.parquet).

Pipeline:        https://github.com/bigbio/quantmsdiann
SDRF tooling:    https://github.com/bigbio/sdrf-pipelines (convert-diann)
Counting logic:  analysis/count_report_ids.py (canonical, reused here).

Counting convention (identical to count_report_ids, the manuscript-wide method)
-------------------------------------------------------------------------------
* protein groups: unique ``Protein.Group`` with ``Global.PG.Q.Value <= 0.01``,
  target-only (drop any group with a CONTAM_/Cont_/ENTRAP_/DECOY_/decoy_ token);
* precursors: distinct ``Precursor.Id`` at run-specific ``Q.Value`` (0.01 for
  1.8.1, 0.05 for >= 2.5.x, DIA-NN's per-version ``--qvalue``) AND
  ``Global.Q.Value <= 0.01``, target-only.

The per-cell and completeness panels count protein groups per run with the SAME
``Global.PG.Q.Value <= 0.01`` target filter, so per-cell counts are a subset of
the cross-run union (= the ``proteins`` total). Dynamic range / CV use
``PG.MaxLFQ`` for the target groups.

NOTE: this counts from the *report* (not the ``*_matrix.tsv`` files); the
matrices bake in ``--matrix-spec-q`` at a version-dependent run q-value and are
not comparable across versions (see count_report_ids.py docstring).

Outputs (data/single_cell/): mv_per_cell.tsv, mv_completeness.tsv,
mv_rank_abundance.tsv, mv_cv.tsv, sc_totals.tsv.

Run:  python -m analysis.make_single_cell_tables
      (downloads ~0.5 GB of reports once; cached under data/single_cell/cache/)
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.count_report_ids import (
    PRECURSOR_Q, DEFAULT_PRECURSOR_Q, Q_THRESHOLD, _CONTAM_RE, count_report,
)

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "data" / "single_cell"
CACHE_DIR = OUT_DIR / "cache"

FTP_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/single-cell"
)
# dataset display name -> PRIDE accession (for labels) and FTP sub-directory.
# One-Tip lives under PXD044991_one-tip/ (PXD044991 also has a mouse-zygote run).
ACC = {"HeLa Astral SC": "PXD046357", "HeLa One-Tip": "PXD044991"}
FTP_DIR = {"HeLa Astral SC": "PXD046357", "HeLa One-Tip": "PXD044991_one-tip"}
VERSIONS = ["1_8_1", "2_5_1_enterprise"]
FLAG = "HeLa Astral SC"  # the dataset carrying the depth/completeness panels

_COLS = [
    "Run", "Precursor.Id", "Protein.Group", "Q.Value", "Global.Q.Value",
    "Global.PG.Q.Value", "PG.Q.Value", "PG.MaxLFQ", "Decoy",
]


def _report_url(ftp_dir: str, version: str) -> str:
    ext = "tsv" if version == "1_8_1" else "parquet"
    return f"{FTP_BASE}/{ftp_dir}/v{version}/quant_tables/diann_report.{ext}"


def _cached_report(ftp_dir: str, version: str) -> Path:
    """Download the deposited DIA-NN report once and cache it on disk."""
    url = _report_url(ftp_dir, version)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / f"{ftp_dir}_v{version}_{Path(url).name}"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"Downloading {url} (cached) ...", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=900) as resp, open(tmp, "wb") as fh:
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
    tmp.replace(dest)
    return dest


def _load(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        have = set(pq.ParquetFile(path).schema_arrow.names)
        return pq.read_table(path, columns=[c for c in _COLS if c in have]).to_pandas()
    return pd.read_csv(path, sep="\t", usecols=lambda c: c in _COLS, low_memory=False)


def _target_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """Rows whose protein group passes Global.PG.Q <= 1% and the target filter."""
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]
    is_target = ~df["Protein.Group"].fillna("").str.contains(_CONTAM_RE)
    return df[(df["Global.PG.Q.Value"] <= Q_THRESHOLD) & is_target]


def build() -> dict[str, pd.DataFrame]:
    per_cell, completeness, rank, cv, totals = [], [], [], [], []
    for ds, acc in ACC.items():
        for version in VERSIONS:
            df = _load(_cached_report(FTP_DIR[ds], version))
            prot = _target_proteins(df)
            pgrun = prot.drop_duplicates(["Run", "Protein.Group"])
            # precursors: canonical count_report_ids prec_min1_tgt (run-specific
            # precursor-q gated AND Global.Q <= 1%, target-only).
            # proteins: cross-run UNION of target protein groups at
            # Global.PG.Q <= 1% -- the protein 1% FDR count, identical to the
            # right endpoint of the completeness panel. (count_report's
            # proteins_tgt additionally gates on precursor-passing rows, which
            # is a slightly smaller set; the union is the value the manuscript
            # reports and that the completeness curve shows.)
            c = count_report(df, precursor_q=PRECURSOR_Q.get(f"v{version}", DEFAULT_PRECURSOR_Q))
            totals.append((ds, version, c["prec_min1_tgt"], int(pgrun["Protein.Group"].nunique())))
            for run, g in pgrun.groupby("Run"):
                per_cell.append((ds, version, int(g["Protein.Group"].nunique())))
            if ds != FLAG:
                continue
            n_runs = pgrun["Run"].nunique()
            seen = pgrun.groupby("Protein.Group")["Run"].nunique()
            for mc in range(1, n_runs + 1):
                completeness.append((ds, version, mc, int((seen >= mc).sum())))
            q = pgrun.dropna(subset=["PG.MaxLFQ"]).copy()
            q["PG.MaxLFQ"] = pd.to_numeric(q["PG.MaxLFQ"], errors="coerce")
            q = q[q["PG.MaxLFQ"] > 0]
            mean_int = q.groupby("Protein.Group")["PG.MaxLFQ"].mean().sort_values(ascending=False)
            for i, val in enumerate(mean_int.values, start=1):
                if i == 1 or i % 10 == 0:
                    rank.append((version, i, float(np.log10(val))))
            agg = q.groupby("Protein.Group")["PG.MaxLFQ"].agg(["mean", "std", "count"])
            agg = agg[agg["count"] >= 3]
            for _, r in agg.iterrows():
                if r["mean"] > 0 and not np.isnan(r["std"]):
                    cv.append((version, float(r["std"] / r["mean"])))
    return {
        "mv_per_cell.tsv": pd.DataFrame(per_cell, columns=["dataset", "version", "pg_count"]),
        "mv_completeness.tsv": pd.DataFrame(completeness, columns=["dataset", "version", "min_cells", "n_proteins"]),
        "mv_rank_abundance.tsv": pd.DataFrame(rank, columns=["version", "rank", "log10_intensity"]),
        "mv_cv.tsv": pd.DataFrame(cv, columns=["version", "cv"]),
        "sc_totals.tsv": pd.DataFrame(totals, columns=["dataset", "version", "precursors", "proteins"]),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, frame in build().items():
        frame.to_csv(OUT_DIR / name, sep="\t", index=False)
        print(f"wrote {OUT_DIR / name} ({len(frame)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
