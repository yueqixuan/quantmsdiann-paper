"""Tests for analysis.compute_cohort_consistency.compute (tiny temp parquet)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


def _write_parquet(tmp_path: Path, rows) -> Path:
    pa = pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq
    df = pd.DataFrame(rows, columns=[
        "Run", "Protein.Group", "Global.Q.Value", "Proteotypic",
        "Stripped.Sequence", "Global.PG.Q.Value",
    ])
    p = tmp_path / "diann_report.parquet"
    pq.write_table(pa.Table.from_pandas(df), p)
    return p


def test_compute_groups_stringent_and_global(tmp_path):
    from analysis.compute_cohort_consistency import compute
    rows = [
        # PG1: proteotypic, 2 peptides, passes FDR in r1(g=A) and r2(g=B)
        ("r1", "PG1", 0.005, 1, "PEPONE", 0.004),
        ("r1", "PG1", 0.005, 1, "PEPTWO", 0.004),
        ("r2", "PG1", 0.006, 1, "PEPONE", 0.004),
        # PG2: only 1 peptide, passes in r1 -> not stringent
        ("r1", "PG2", 0.002, 1, "PEPX", 0.006),
        # PG3: fails proteotypic -> excluded from per-group/stringent; but its
        # global PG q passes so it counts toward proteins_global.
        ("r2", "PG3", 0.001, 0, "PEPY", 0.001),
        # Cont_: contaminant -> unfiltered only
        ("r1", "Cont_ALB", 0.001, 1, "PEPC1", 0.001),
        ("r2", "Cont_ALB", 0.001, 1, "PEPC2", 0.001),
    ]
    p = _write_parquet(tmp_path, rows)
    run2group = {"r1": "A", "r2": "B"}
    agg = compute(p, run2group)
    # group A has PG1, PG2 (+ Cont in unfiltered sense but group_pg keeps all pg)
    assert set(agg["group_pg"].keys()) == {"A", "B"}
    assert "PG1" in agg["group_pg"]["A"] and "PG1" in agg["group_pg"]["B"]
    # stringent target-only: only PG1 has >=2 peptides and is target
    assert agg["stringent_tgt"] == 1
    # stringent unfiltered: PG1 and Cont_ALB (2 peptides across runs)
    assert agg["stringent_unf"] == 2
    # proteins_global target-only: PG1, PG2, PG3 (Cont excluded)
    assert agg["proteins_global_tgt"] == 3
    assert agg["proteins_global_unf"] == 4  # + Cont_ALB
    # per-run protein counts (proteotypic, FDR-passing, mapped)
    assert agg["per_run_protein_count"]["r1"] >= 2


def test_compute_consistency_filter_drops_rare(tmp_path):
    from analysis.compute_cohort_consistency import compute
    # PG1 in 1 of 10 runs (10% detection); PG2 in 1 run only.
    rows = []
    for i in range(10):
        rows.append((f"r{i}", "PG1" if i == 0 else "PGother",
                     0.001, 1, "PEP", 0.001))
    rows.append(("r0", "PG2", 0.001, 1, "PEP2", 0.001))
    p = _write_parquet(tmp_path, rows)
    run2group = {f"r{i}": "A" for i in range(10)}
    # threshold 0.10 * 10 = 1.0 -> PG1 (1 run) kept, PG2 (1 run) kept
    agg = compute(p, run2group, min_detection_fraction=0.10)
    assert "PG1" in agg["group_pg"]["A"]
    # threshold 0.30 * 10 = 3 -> PG1 and PG2 dropped, PGother kept (9 runs)
    agg2 = compute(p, run2group, min_detection_fraction=0.30)
    assert "PG1" not in agg2["group_pg"].get("A", [])
    assert "PGother" in agg2["group_pg"]["A"]
