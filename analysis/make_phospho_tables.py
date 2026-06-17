#!/usr/bin/env python
"""Generate the phospho supplementary figure's input table from public reports.

Reproducibility generator for ``figure_phospho.py`` (Panels A/B): downloads the
deposited DIA-NN phospho reports from the public PRIDE FTP and computes
phosphopeptide and phosphosite counts, so nothing is hand-entered.

PROVENANCE
==========
Our reanalyses are on the PRIDE FTP under
``quantmsdiann-benchmarks/phospho/<dataset>/v<version>/quant_tables/``:
  PXD034128-biological-study, PXD034128-highspeed-DIA, PXD049692, PXD034623
each with ``diann_report.parquet`` and ``diann_report.site_report.parquet``,
for DIA-NN ``v2_5_1`` and ``v2_5_1_enterprise``.

Pipeline: https://github.com/bigbio/quantmsdiann

Definitions (1% FDR, target-only)
---------------------------------
* phosphopeptides: distinct ``Modified.Sequence`` carrying ``UniMod:21`` with
  ``Q.Value <= 0.01`` and ``Global.Q.Value <= 0.01`` (target protein groups);
* sites_all: phospho sites from the DIA-NN site report (``Modification``
  contains ``UniMod:21``), unique by ``(Protein, Site)``;
* sites_classI: the same restricted to localization ``Probability >= 0.99``.

Output: data/phospho/phospho_counts.tsv (consumed by figure_phospho.py).

Run:  python -m analysis.make_phospho_tables
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import pandas as pd

from analysis.count_report_ids import _CONTAM_RE, Q_THRESHOLD

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "data" / "phospho"
CACHE_DIR = OUT_DIR / "cache"
FTP_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/phospho"
)
PHOSPHO = "UniMod:21"
VERSIONS = ["2_5_1", "2_5_1_enterprise"]
# display name -> FTP dataset directory
DATASETS = {
    "PXD034128 biological-study": "PXD034128-biological-study",
    "PXD034128 highspeed-DIA": "PXD034128-highspeed-DIA",
    "PXD049692 NK-phospho": "PXD049692",
    "PXD034623 Galectin1": "PXD034623",
}


def _cached(ftp_dir: str, version: str, fname: str) -> Path:
    url = f"{FTP_BASE}/{ftp_dir}/v{version}/quant_tables/{fname}"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / f"{ftp_dir}_v{version}_{fname}"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"Downloading {url} (cached) ...", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=900) as resp, open(tmp, "wb") as fh:
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
    tmp.replace(dest)
    return dest


def _count(ftp_dir: str, version: str) -> tuple[int, int, int]:
    import pyarrow.parquet as pq
    rep = _cached(ftp_dir, version, "diann_report.parquet")
    r = pq.read_table(rep, columns=[
        "Modified.Sequence", "Q.Value", "Global.Q.Value", "Protein.Group",
    ]).to_pandas()
    is_target = ~r["Protein.Group"].fillna("").str.contains(_CONTAM_RE)
    r = r[(r["Q.Value"] <= Q_THRESHOLD) & (r["Global.Q.Value"] <= Q_THRESHOLD) & is_target]
    phosphopeptides = r.loc[
        r["Modified.Sequence"].str.contains(PHOSPHO, na=False), "Modified.Sequence"
    ].nunique()

    site = _cached(ftp_dir, version, "diann_report.site_report.parquet")
    s = pq.read_table(site, columns=["Protein", "Site", "Modification", "Probability"]).to_pandas()
    ph = s[s["Modification"].astype(str).str.contains(PHOSPHO, na=False)]
    sites_all = ph.drop_duplicates(["Protein", "Site"]).shape[0]
    sites_classI = ph[ph["Probability"] >= 0.99].drop_duplicates(["Protein", "Site"]).shape[0]
    return phosphopeptides, sites_classI, sites_all


def main() -> int:
    rows = []
    for name, ftp_dir in DATASETS.items():
        for version in VERSIONS:
            pp, c1, sa = _count(ftp_dir, version)
            rows.append((name, version, pp, c1, sa))
            print(f"{name} {version}: {pp} phosphopeptides, {c1} class-I, {sa} sites")
    df = pd.DataFrame(rows, columns=["dataset", "version", "phosphopeptides", "sites_classI", "sites_all"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DIR / "phospho_counts.tsv", sep="\t", index=False)
    print(f"wrote {OUT_DIR / 'phospho_counts.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
