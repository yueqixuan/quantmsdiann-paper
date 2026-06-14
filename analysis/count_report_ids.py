"""Count precursors and protein groups from the DIA-NN *report*, not the matrices.

The `*_matrix.tsv` outputs bake in `--matrix-spec-q` (0.05 run-specific) and,
because the quantmsdiann pipeline sets `--qvalue` to 0.01 for DIA-NN 1.8.1 but
0.05 for 2.5.1/enterprise, matrix row counts are filtered at *different*
run-specific q-values per version and are NOT comparable across versions.

This module counts identifications directly from the per-precursor report
(`diann_report.parquet` for DIA-NN >= 2.x, `diann_report.tsv` for 1.8.1). The
run-specific precursor q-value cut-off is DIA-NN's *recommended* per-version
operating point (1% for 1.8.1, 5% for >= 2.5.0 -- this is the `--qvalue` the
pipeline passes, by design, not a parameter slip); the global precursor q and
the protein-group q are held at 1% for every version:

  * a precursor (`Precursor.Id`) is identified in a run when
    `Q.Value <= precursor_q` (run-specific; 0.01 for 1.8.1, 0.05 for >= 2.5.0)
    AND `Global.Q.Value <= 0.01` (global);
  * `prec_min1` / `prec_min3` = precursors identified in >= 1 / >= 3 runs;
  * `proteins` = unique `Protein.Group` with `Global.PG.Q.Value <= 0.01`
    (the cross-run union; used for the vs-deposited cohort figures).
  * `prot_avg` / `prot_complete` = average proteins per run, and proteins
    quantified in *every* run, each at run-specific `PG.Q.Value <= 0.01`
    (target-only). These are sensitive to per-run depth/grouping and are the
    version-comparison metric for the benchmark figure (DIA-NN author's method).

We still count from the report rather than the `*_matrix.tsv` files because the
matrices additionally apply `--matrix-spec-q`/`--matrix-qvalue` and, for the
protein groups, a row count that does not equal the global-q protein count.

Target-only drops contaminant / entrapment / decoy protein groups (the repo's
conservative filter); the unfiltered companion is kept for the audit TSV.

The reports are multi-GB, so this is run ONCE on the cluster (where they live)
and its small output, `report_counts.tsv`, is staged into
`data/quantmsdiann_benchmarks/` and consumed by
`figure_quantmsdiann_benchmarks_vs_proteobench.py`.

Usage (on the cluster):
    python -m analysis.count_report_ids \
        --results-root /hps/.../quantmsdiann_results \
        --out report_counts.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

Q_THRESHOLD = 0.01           # global precursor q and protein-group q (all versions)
# DIA-NN's recommended run-specific precursor q-value (`--qvalue`) per release.
PRECURSOR_Q = {
    "v1_8_1": 0.01,
    "v2_5_1": 0.05,
    "v2_5_1_enterprise": 0.05,
}
DEFAULT_PRECURSOR_Q = 0.01
# Contaminant / entrapment / decoy prefixes (any semicolon-separated token).
_CONTAM_RE = re.compile(
    r"(?:^|;)(?:Cont_|CONTAM_|ENTRAP_|DECOY_|decoy_)", re.IGNORECASE
)
_NEEDED_COLS = [
    "Run", "Precursor.Id", "Protein.Group", "Stripped.Sequence",
    "Q.Value", "Global.Q.Value", "Global.PG.Q.Value", "PG.Q.Value", "Decoy",
]

DATASET_MODULES = (
    "ProteoBench_Module_7", "PXD049412", "PXD062685", "PXD070049",
)
VERSIONS = ("v1_8_1", "v2_5_1", "v2_5_1_enterprise")


def _load_report(report_dir: Path) -> pd.DataFrame:
    """Read the DIA-NN report from `report_dir`, preferring the parquet
    (DIA-NN >= 2.x) and falling back to the classic `diann_report.tsv`
    (1.8.1). Only the columns we need are read."""
    parquet = report_dir / "diann_report.parquet"
    tsv = report_dir / "diann_report.tsv"
    if parquet.exists():
        import pyarrow.parquet as pq
        have = set(pq.ParquetFile(parquet).schema_arrow.names)
        cols = [c for c in _NEEDED_COLS if c in have]
        return pq.read_table(parquet, columns=cols).to_pandas()
    if tsv.exists():
        return pd.read_csv(
            tsv, sep="\t", usecols=lambda c: c in _NEEDED_COLS,
        )
    raise FileNotFoundError(
        f"no diann_report.parquet or .tsv in {report_dir}"
    )


def count_report(df: pd.DataFrame, precursor_q: float = DEFAULT_PRECURSOR_Q) -> dict[str, int]:
    """Compute precursor (min1/min3) and protein-group counts, both
    target-only and unfiltered, from a DIA-NN report frame. `precursor_q` is
    the run-specific precursor q-value cut-off (DIA-NN's per-version
    `--qvalue`); the global precursor q and protein-group q stay at 1%."""
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]
    is_target = ~df["Protein.Group"].fillna("").str.contains(_CONTAM_RE)
    passes = (df["Q.Value"] <= precursor_q) & (df["Global.Q.Value"] <= Q_THRESHOLD)

    def counts(mask: pd.Series) -> tuple[int, int, int]:
        d = df[passes & mask]
        n_runs = d.groupby("Precursor.Id")["Run"].nunique()
        proteins = d.loc[
            d["Global.PG.Q.Value"] <= Q_THRESHOLD, "Protein.Group"
        ].nunique()
        return int((n_runs >= 1).sum()), int((n_runs >= 3).sum()), int(proteins)

    all_mask = pd.Series(True, index=df.index)
    u1, u3, up = counts(all_mask)
    t1, t3, tp = counts(is_target)

    # Peptide-level metrics (target-only), using the same passing-precursor
    # filter. peptides = distinct stripped backbones; proteins_2pep = protein
    # groups (Global.PG.Q <= 1%) supported by >= 2 distinct stripped peptides.
    pep_tgt = prot_2pep_tgt = 0
    if "Stripped.Sequence" in df.columns:
        dt = df[passes & is_target]
        pep_tgt = int(dt["Stripped.Sequence"].nunique())
        prot_pep = dt[dt["Global.PG.Q.Value"] <= Q_THRESHOLD]
        per_prot = prot_pep.groupby("Protein.Group")["Stripped.Sequence"].nunique()
        prot_2pep_tgt = int((per_prot >= 2).sum())

    # Per-run-average and complete-profile protein counts (run-specific
    # PG.Q.Value <= 1%, target-only). Unlike the Global.PG.Q union above these
    # are sensitive to per-run depth and grouping, so they are the metric used
    # for the version-comparison benchmark figure (the DIA-NN author's method:
    # average proteins per run, and proteins quantified in every run).
    prot_avg = prot_complete = 0
    if "PG.Q.Value" in df.columns:
        n_runs = df["Run"].nunique()
        pg = df[(df["PG.Q.Value"] <= Q_THRESHOLD) & is_target][
            ["Run", "Protein.Group"]
        ].drop_duplicates()
        per_run = pg.groupby("Run")["Protein.Group"].nunique()
        prot_avg = int(round(per_run.mean())) if len(per_run) else 0
        in_n_runs = pg.groupby("Protein.Group")["Run"].nunique()
        prot_complete = int((in_n_runs == n_runs).sum()) if n_runs else 0

    return {
        "prec_min1_tgt": t1, "prec_min1_unf": u1,
        "prec_min3_tgt": t3, "prec_min3_unf": u3,
        "proteins_tgt": tp, "proteins_unf": up,
        "proteins_2pep_tgt": prot_2pep_tgt, "peptides_tgt": pep_tgt,
        "prot_avg_tgt": prot_avg, "prot_complete_tgt": prot_complete,
    }


def main(argv: list[str] | None = None) -> int:  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-root", required=True, type=Path,
                    help="quantmsdiann_results dir holding <module>/<version>/quant_tables")
    ap.add_argument("--out", required=True, type=Path,
                    help="output report_counts.tsv path")
    ap.add_argument("--results-suffix", default="",
                    help="suffix for the per-version dir, e.g. '_relaxed' reads "
                         "<module>/<version-without-v><suffix>/quant_tables (the "
                         "--relaxed-prot-inf re-run dirs '1_8_1_relaxed' etc.)")
    args = ap.parse_args(argv)

    rows: list[dict] = []
    for module in DATASET_MODULES:
        for version in VERSIONS:
            if args.results_suffix:
                base = version[1:] if version.startswith("v") else version
                ver_dir = f"{base}{args.results_suffix}"
            else:
                ver_dir = version
            rdir = args.results_root / module / ver_dir / "quant_tables"
            try:
                df = _load_report(rdir)
            except FileNotFoundError as exc:
                print(f"WARN: {exc}", file=sys.stderr)
                continue
            c = count_report(df, precursor_q=PRECURSOR_Q.get(version, DEFAULT_PRECURSOR_Q))
            c.update(dataset=module, version=version)
            rows.append(c)
            print(f"{module} {version}: {c['prec_min1_tgt']:,} prec / "
                  f"{c['proteins_tgt']:,} proteins (target, q<=0.01)")
    cols = ["dataset", "version", "prec_min1_tgt", "prec_min1_unf",
            "prec_min3_tgt", "prec_min3_unf", "proteins_tgt", "proteins_unf",
            "proteins_2pep_tgt", "peptides_tgt",
            "prot_avg_tgt", "prot_complete_tgt"]
    pd.DataFrame(rows)[cols].to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
