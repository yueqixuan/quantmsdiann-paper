"""Tests for analysis.count_report_ids.count_report (inline fixtures, no I/O)."""
from __future__ import annotations

import pandas as pd


def _row(run, pid, pg, q, gq, pgq, decoy=0, pgq_run=None):
    d = {
        "Run": run, "Precursor.Id": pid, "Protein.Group": pg,
        "Q.Value": q, "Global.Q.Value": gq, "Global.PG.Q.Value": pgq,
        "Decoy": decoy,
    }
    if pgq_run is not None:
        d["PG.Q.Value"] = pgq_run  # run-specific protein-group q-value
    return d


def test_count_report_per_run_average_and_complete_proteins():
    """prot_avg_tgt = mean proteins/run, prot_complete_tgt = proteins in every
    run, both at run-specific PG.Q.Value <= 0.01, target-only."""
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        # PGa: in all 3 runs -> complete
        _row("r1", "x1", "PGa", 0.001, 0.001, 0.001, pgq_run=0.001),
        _row("r2", "x1", "PGa", 0.001, 0.001, 0.001, pgq_run=0.001),
        _row("r3", "x1", "PGa", 0.001, 0.001, 0.001, pgq_run=0.001),
        # PGb: in 2 runs
        _row("r1", "x2", "PGb", 0.001, 0.001, 0.001, pgq_run=0.001),
        _row("r2", "x2", "PGb", 0.001, 0.001, 0.001, pgq_run=0.001),
        # PGc: in 1 run
        _row("r3", "x3", "PGc", 0.001, 0.001, 0.001, pgq_run=0.001),
        # contaminant -> excluded from per-run/complete
        _row("r1", "x4", "Cont_ALBU", 0.001, 0.001, 0.001, pgq_run=0.001),
        # run-specific PG.Q fails -> excluded
        _row("r2", "x5", "PGd", 0.001, 0.001, 0.001, pgq_run=0.20),
    ])
    c = count_report(df)
    # per run: r1={PGa,PGb}=2, r2={PGa,PGb}=2, r3={PGa,PGc}=2 -> avg 2
    assert c["prot_avg_tgt"] == 2
    # complete (in all 3 runs): only PGa
    assert c["prot_complete_tgt"] == 1


def test_count_report_per_run_metric_absent_without_pgq_column():
    """When the report has no PG.Q.Value column, the per-run/complete metrics
    fall back to 0 (older reports / fixtures) instead of erroring."""
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([_row("r1", "x1", "PGa", 0.001, 0.001, 0.001)])
    c = count_report(df)
    assert c["prot_avg_tgt"] == 0
    assert c["prot_complete_tgt"] == 0


def test_count_report_applies_run_and_global_q_and_replicate_threshold():
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        # P1: passes (run+global<=0.01) in 3 runs -> min1 and min3
        _row("r1", "P1", "PG1", 0.005, 0.005, 0.005),
        _row("r2", "P1", "PG1", 0.004, 0.005, 0.005),
        _row("r3", "P1", "PG1", 0.006, 0.005, 0.005),
        # P2: passes in only 1 run (other run fails run-q) -> min1 only
        _row("r1", "P2", "PG2", 0.002, 0.008, 0.009),
        _row("r2", "P2", "PG2", 0.050, 0.008, 0.009),
        # P3: fails global q everywhere -> not counted
        _row("r1", "P3", "PG3", 0.002, 0.050, 0.050),
        # P4: contaminant -> in unfiltered, out of target-only
        _row("r1", "P4", "Cont_ALBU", 0.001, 0.001, 0.001),
        _row("r2", "P4", "Cont_ALBU", 0.001, 0.001, 0.001),
        _row("r3", "P4", "Cont_ALBU", 0.001, 0.001, 0.001),
        # P5: a decoy row -> dropped entirely
        _row("r1", "P5", "PG5", 0.001, 0.001, 0.001, decoy=1),
    ])
    c = count_report(df)
    # target precursors min1: P1, P2 (P3 fails global, P4 contaminant, P5 decoy)
    assert c["prec_min1_tgt"] == 2
    # target precursors min3: only P1
    assert c["prec_min3_tgt"] == 1
    # unfiltered min1 adds the contaminant P4 -> 3
    assert c["prec_min1_unf"] == 3
    # unfiltered min3 adds P4 (3 runs) -> P1 + P4 = 2
    assert c["prec_min3_unf"] == 2
    # target proteins: PG1, PG2 (PG3 fails, Cont_ excluded, decoy dropped)
    assert c["proteins_tgt"] == 2
    assert c["proteins_unf"] == 3  # + Cont_ALBU


def test_count_report_precursor_q_affects_replicate_counting():
    from analysis.count_report_ids import count_report
    # P1 passes run-q at 1% in run r1 and at 3% (>1%, <=5%) in r2; global ok.
    df = pd.DataFrame([
        _row("r1", "P1", "PG1", 0.005, 0.005, 0.005),
        _row("r2", "P1", "PG1", 0.030, 0.005, 0.005),
    ])
    # At precursor_q=0.01 P1 is identified in only r1 -> min1 yes, not >=2 runs.
    c01 = count_report(df, precursor_q=0.01)
    assert c01["prec_min1_tgt"] == 1
    # At precursor_q=0.05 P1 is identified in both runs.
    c05 = count_report(df, precursor_q=0.05)
    assert c05["prec_min1_tgt"] == 1
    # Difference shows up in run multiplicity: extend to 3 runs to probe min3.
    df3 = pd.concat([df, pd.DataFrame([_row("r3", "P1", "PG1", 0.030, 0.005, 0.005)])],
                    ignore_index=True)
    assert count_report(df3, precursor_q=0.01)["prec_min3_tgt"] == 0  # only r1 counts
    assert count_report(df3, precursor_q=0.05)["prec_min3_tgt"] == 1  # r1+r2+r3


def test_count_report_protein_requires_global_pg_q():
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        # precursor passes but its protein group fails global PG q -> protein
        # not counted, precursor still counted.
        _row("r1", "P1", "PG1", 0.001, 0.001, 0.20),
        _row("r1", "P2", "PG2", 0.001, 0.001, 0.005),
    ])
    c = count_report(df)
    assert c["prec_min1_tgt"] == 2
    assert c["proteins_tgt"] == 1  # only PG2
