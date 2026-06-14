r"""Per-cohort identification aggregates from a DIA-NN report (run ON THE CLUSTER).

The per-cohort reanalysis figures (Sun/PXD004701, ProCan/PXD030304, ...) need a
few summaries derived from the multi-GB `diann_report.parquet` (1.3-40 GB):

  * per-group protein sets -- the union of `Protein.Group` with a proteotypic
    precursor at `Global.Q.Value <= q` in >= 1 run mapped to that group, where
    a "group" is the cohort's biological grouping (subtype, tissue, ...);
    optionally a consistency threshold (Sun: detected in >= 10 %% of a group's
    runs) is applied;
  * `stringent` -- protein groups with >= 2 unique proteotypic peptides
    (`Stripped.Sequence`) at `Global.Q.Value <= q`;
  * `proteins_global` -- unique `Protein.Group` at `Global.PG.Q.Value <= q`
    (the report-based protein count, replacing the pg_matrix row count);
  * `per_run` -- protein count per run (for completeness panels).

Both target-only (contaminant/entrapment/decoy stripped) and unfiltered counts
are emitted. The reports are too large to stage locally, so this is run on the
cluster (one streaming pass, column projection) and its small JSON output is
staged into the paper repo under data/<COHORT>/ for the figure to read.

Usage (cluster, ideally via sbatch):
    python compute_cohort_consistency.py \
        --parquet  .../<COHORT>/results/quant_tables/diann_report.parquet \
        --run-group run_to_group.json \   # {run_basename: group_label}
        --out      cohort_agg.json \
        [--min-detection-fraction 0.10]    # 0 disables the consistency filter
"""
from __future__ import annotations

import argparse
import json
import re
import sys

_CONTAM_RE = re.compile(
    r"(?:^|;)(?:Cont_|CONTAM_|ENTRAP_|DECOY_|decoy_)", re.IGNORECASE
)


def _is_target(pg: str | None) -> bool:
    return not _CONTAM_RE.search(pg or "")


def compute(parquet_path, run_to_group, q=0.01, min_detection_fraction=0.0,
            batch_size=2_000_000):
    """Stream the parquet once and return the aggregate dict. `run_to_group`
    maps run basename -> group label (runs absent from the map are ignored for
    the per-group/per-run sets but still counted for the global protein set)."""
    import pyarrow.parquet as pq

    pf = pq.ParquetFile(str(parquet_path))
    cols = ["Run", "Protein.Group", "Global.Q.Value", "Proteotypic",
            "Stripped.Sequence", "Global.PG.Q.Value"]
    have = set(pf.schema_arrow.names)
    cols = [c for c in cols if c in have]

    pg_runs: dict[str, set[str]] = {}          # pg -> mapped runs passing FDR
    pg_pep: dict[str, set[str]] = {}           # pg -> proteotypic peptides
    pg_global: set[str] = set()                # pg at Global.PG.Q.Value <= q
    per_run: dict[str, set[str]] = {}          # run -> pgs (for completeness)
    group_runs: dict[str, set[str]] = {}       # group -> its mapped runs seen

    for b in pf.iter_batches(batch_size=batch_size, columns=cols):
        R = b.column("Run").to_pylist()
        P = b.column("Protein.Group").to_pylist()
        G = b.column("Global.Q.Value").to_pylist()
        T = b.column("Proteotypic").to_pylist()
        S = b.column("Stripped.Sequence").to_pylist() if "Stripped.Sequence" in cols else [None] * len(R)
        Q = b.column("Global.PG.Q.Value").to_pylist() if "Global.PG.Q.Value" in cols else [None] * len(R)
        for r, pg, g, t, s, q_pg in zip(R, P, G, T, S, Q):
            if q_pg is not None and q_pg <= q:
                pg_global.add(pg)
            if t != 1 or g is None or g > q:
                continue
            pg_pep.setdefault(pg, set()).add(s)
            grp = run_to_group.get(r)
            if grp is None:
                continue
            pg_runs.setdefault(pg, set()).add(r)
            per_run.setdefault(r, set()).add(pg)
            group_runs.setdefault(grp, set()).add(r)

    # Optional consistency filter (Sun-style): keep pg detected in >= fraction
    # of ALL mapped runs.
    total_mapped = len({r for runs in group_runs.values() for r in runs})
    keep = pg_runs
    if min_detection_fraction > 0 and total_mapped:
        thr = min_detection_fraction * total_mapped
        keep = {pg: runs for pg, runs in pg_runs.items() if len(runs) >= thr}

    group_pg: dict[str, set[str]] = {}
    for pg, runs in keep.items():
        for r in runs:
            grp = run_to_group.get(r)
            if grp is None or grp == "unknown":
                continue
            group_pg.setdefault(grp, set()).add(pg)

    return {
        "group_pg": {k: sorted(v) for k, v in group_pg.items()},
        "per_run_protein_count": {r: len(pgs) for r, pgs in per_run.items()},
        "stringent_tgt": sum(1 for pg, ps in pg_pep.items()
                             if len(ps) >= 2 and _is_target(pg)),
        "stringent_unf": sum(1 for ps in pg_pep.values() if len(ps) >= 2),
        "proteins_global_tgt": sum(1 for pg in pg_global if _is_target(pg)),
        "proteins_global_unf": len(pg_global),
        "total_mapped_runs": total_mapped,
    }


def main(argv=None) -> int:  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--run-group", required=True,
                    help="JSON: {run_basename: group_label}")
    ap.add_argument("--out", required=True)
    ap.add_argument("--q", type=float, default=0.01)
    ap.add_argument("--min-detection-fraction", type=float, default=0.0)
    args = ap.parse_args(argv)
    run_to_group = json.load(open(args.run_group))
    agg = compute(args.parquet, run_to_group, q=args.q,
                  min_detection_fraction=args.min_detection_fraction)
    json.dump(agg, open(args.out, "w"))
    print(f"wrote {args.out}: {len(agg['group_pg'])} groups, "
          f"proteins_global_tgt={agg['proteins_global_tgt']}, "
          f"stringent_tgt={agg['stringent_tgt']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
