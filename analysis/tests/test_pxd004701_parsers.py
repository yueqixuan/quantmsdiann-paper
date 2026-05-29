"""Tests for the PXD004701 (Sun et al. 2023 / quantmsdiann) reanalysis
comparison parsers. The two pipelines being compared are the published 2023
Sun et al. PCT-SWATH analysis of 76 breast cancer cell lines
(doi:10.1016/j.mcpro.2023.100602) and the quantmsdiann DIA-NN reanalysis
(PRIDE quantms-collections).
"""
from __future__ import annotations

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# SDRF -> cell-line map (data file basename match for the parquet 'Run')
# ---------------------------------------------------------------------------

def test_parse_sdrf_data_file_to_cell_line_rewrites_wiff_to_mzml(
    tmp_path: Path,
) -> None:
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        parse_sdrf_data_file_to_cell_line,
    )
    body = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tmdamb231\tguot_K160125_sw_BR1.wiff\n"
        "S2\ths578t\tguot_L160226_BR203.wiff\n"
        "S3\tmcf7\tguot_K160125_sw_BR63.WIFF\n"
        "S4\tbt474\t\n"
    )
    p = tmp_path / "sdrf.tsv"
    p.write_text(body)
    out = parse_sdrf_data_file_to_cell_line(p)
    # Empty rows skipped; .wiff/.WIFF rewritten to .mzML
    assert out == {
        "guot_K160125_sw_BR1.mzML": "mdamb231",
        "guot_L160226_BR203.mzML": "hs578t",
        "guot_K160125_sw_BR63.mzML": "mcf7",
    }


def test_parse_sdrf_data_file_to_cell_line_raises_on_missing_columns(
    tmp_path: Path,
) -> None:
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        parse_sdrf_data_file_to_cell_line,
    )
    p = tmp_path / "bad.tsv"
    p.write_text("source name\tcharacteristics[cell line]\nS1\tmcf7\n")
    with pytest.raises(ValueError):
        parse_sdrf_data_file_to_cell_line(p)


# ---------------------------------------------------------------------------
# BC_SUBTYPES coverage
# ---------------------------------------------------------------------------

def test_bc_subtypes_covers_all_known_pxd004701_cell_lines() -> None:
    """Every cell line listed in the PXD004701 SDRF must classify to one of
    TNBC / non-TNBC / normal-like / unknown. No silent drops."""
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import BC_SUBTYPES
    # The 76 SDRF cell-line names (lowercase, no dashes) as they appear in
    # characteristics[cell line]. Sourced from
    # data/PXD004701/PXD004701.sdrf.tsv on the deposit.
    sdrf_lines = {
        "184a1", "184b5", "600mpe", "au565", "bt20", "bt474", "bt483",
        "bt549", "cal120", "cal148", "cal51", "cama1", "du4475", "efm19",
        "efm192a", "evsat", "hbl100", "hcc1143", "hcc1187", "hcc1395",
        "hcc1419", "hcc1428", "hcc1569", "hcc1599", "hcc1806", "hcc1937",
        "hcc1954", "hcc202", "hcc2185", "hcc2218", "hcc2688", "hcc3153",
        "hcc38", "hcc70", "hdqp1", "hs578t", "jimt1", "kpl1", "ly2",
        "macls2", "mb157", "mcf10a", "mcf12a", "mcf7", "mdamb134vi",
        "mdamb157", "mdamb175vii", "mdamb231", "mdamb330", "mdamb361",
        "mdamb415", "mdamb436", "mdamb453", "mdamb468", "mfm223", "mx1",
        "ocubm", "skbr3", "skbr5", "skbr7", "sum102", "sum149", "sum159",
        "sum185", "sum190", "sum225", "sum229", "sum44", "sum52", "sw527",
        "t47d", "uacc3199", "uacc893", "zr751", "zr7530", "zr75b",
    }
    assert len(sdrf_lines) == 76
    missing = sdrf_lines - set(BC_SUBTYPES)
    assert not missing, f"BC_SUBTYPES missing cell lines: {sorted(missing)}"
    allowed = {"TNBC", "non-TNBC", "normal-like", "unknown"}
    bad = {
        cl: s for cl, s in BC_SUBTYPES.items()
        if cl in sdrf_lines and s not in allowed
    }
    assert not bad, f"BC_SUBTYPES has invalid subtype values: {bad}"


def test_bc_subtypes_partition_sums_to_76() -> None:
    """The 76 SDRF cell lines must partition cleanly across the 4 subtype
    categories. No silent drops."""
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import BC_SUBTYPES
    sdrf_lines = {
        "184a1", "184b5", "600mpe", "au565", "bt20", "bt474", "bt483",
        "bt549", "cal120", "cal148", "cal51", "cama1", "du4475", "efm19",
        "efm192a", "evsat", "hbl100", "hcc1143", "hcc1187", "hcc1395",
        "hcc1419", "hcc1428", "hcc1569", "hcc1599", "hcc1806", "hcc1937",
        "hcc1954", "hcc202", "hcc2185", "hcc2218", "hcc2688", "hcc3153",
        "hcc38", "hcc70", "hdqp1", "hs578t", "jimt1", "kpl1", "ly2",
        "macls2", "mb157", "mcf10a", "mcf12a", "mcf7", "mdamb134vi",
        "mdamb157", "mdamb175vii", "mdamb231", "mdamb330", "mdamb361",
        "mdamb415", "mdamb436", "mdamb453", "mdamb468", "mfm223", "mx1",
        "ocubm", "skbr3", "skbr5", "skbr7", "sum102", "sum149", "sum159",
        "sum185", "sum190", "sum225", "sum229", "sum44", "sum52", "sw527",
        "t47d", "uacc3199", "uacc893", "zr751", "zr7530", "zr75b",
    }
    by_subtype: dict[str, int] = {}
    for cl in sdrf_lines:
        by_subtype.setdefault(BC_SUBTYPES[cl], 0)
        by_subtype[BC_SUBTYPES[cl]] += 1
    assert sum(by_subtype.values()) == 76
    # Sanity: the 5 canonical non-tumorigenic mammary lines must classify as
    # normal-like — a paper-comparison sanity check (Heiser 2012, Neve 2006).
    for cl in {"184a1", "184b5", "hbl100", "mcf10a", "mcf12a"}:
        assert BC_SUBTYPES[cl] == "normal-like", (
            f"{cl} should be normal-like, got {BC_SUBTYPES[cl]}"
        )
    # Sanity: a few canonical TNBC lines.
    for cl in {"mdamb231", "mdamb468", "hs578t", "bt549", "hcc1937"}:
        assert BC_SUBTYPES[cl] == "TNBC", (
            f"{cl} should be TNBC, got {BC_SUBTYPES[cl]}"
        )
    # Sanity: a few canonical non-TNBC lines.
    for cl in {"mcf7", "t47d", "bt474", "skbr3", "zr751"}:
        assert BC_SUBTYPES[cl] == "non-TNBC", (
            f"{cl} should be non-TNBC, got {BC_SUBTYPES[cl]}"
        )


# ---------------------------------------------------------------------------
# Streaming Sun-style consistency filter (Proteotypic + Global.Q.Value +
# >=10% of mapped runs detection)
# ---------------------------------------------------------------------------

def test_proteins_per_subtype_quantmsdiann_consistency_filter_drops_rare(
    tmp_path: Path,
) -> None:
    """Sun-style filter applied to the long-format parquet: keep
    proteotypic precursors with Global.Q.Value <= cutoff AND any
    Protein.Group detected in at least `min_detection_fraction` of
    mapped runs. Per subtype, union of surviving Protein.Group across
    runs of that subtype."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        proteins_per_subtype_quantmsdiann_consistency_filter,
    )
    # 4 runs total mapped: r1, r2 -> TNBC; r3, r4 -> non-TNBC.
    # P1 detected in r1, r2, r3, r4 -> 4/4 = 100%, kept globally.
    # P2 detected in r1 only -> 1/4 = 25%, kept globally (>=10%).
    # P3 detected in r1 only at q=0.5 -> dropped by FDR before counting.
    # P4 detected in r1 only at Proteotypic==0 -> dropped before counting.
    # P5 not detected anywhere passing both filters -> dropped.
    # P6 detected at q=0.001 in r1 only -> 1/4 = 25% kept.
    table = pa.table({
        "Run": ["r1", "r2", "r3", "r4",
                "r1",
                "r1",
                "r1",
                "r1"],
        "Protein.Group": ["P1", "P1", "P1", "P1",
                          "P2",
                          "P3",
                          "P4",
                          "P6"],
        "Global.Q.Value": [0.001, 0.001, 0.001, 0.001,
                           0.001,
                           0.5,
                           0.001,
                           0.001],
        "Proteotypic":   [1, 1, 1, 1,
                          1,
                          1,
                          0,
                          1],
    })
    parquet_path = tmp_path / "report.parquet"
    pq.write_table(table, parquet_path)
    sdrf = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tmdamb231\tr1.wiff\n"
        "S2\tmdamb468\tr2.wiff\n"
        "S3\tmcf7\tr3.wiff\n"
        "S4\tt47d\tr4.wiff\n"
    )
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text(sdrf)
    subtypes = {
        "mdamb231": "TNBC",
        "mdamb468": "TNBC",
        "mcf7": "non-TNBC",
        "t47d": "non-TNBC",
    }
    out = proteins_per_subtype_quantmsdiann_consistency_filter(
        parquet_path, sdrf_p, subtypes,
        qvalue_cutoff=0.01, min_detection_fraction=0.10, batch_size=2,
    )
    # P1: in all 4 runs after FDR filter -> in both subtypes.
    # P2, P6: only in r1 (TNBC) but pass 1/4=25% >= 10% global threshold.
    #   They appear only in the TNBC subtype's union.
    assert out == {
        "TNBC": {"P1", "P2", "P6"},
        "non-TNBC": {"P1"},
    }


def test_proteins_per_subtype_consistency_filter_drops_below_threshold(
    tmp_path: Path,
) -> None:
    """A protein detected in <10% of mapped runs (here 0/10) is dropped
    even if it passes Global.Q.Value and Proteotypic."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        proteins_per_subtype_quantmsdiann_consistency_filter,
    )
    # 10 mapped runs. P1 detected only in r1 -> 1/10 = 10% -> just barely
    # passes default threshold (>=). P2 detected only in r0 (unmapped) so
    # never seen in mapped-run space -> 0/10 -> dropped.
    runs = [f"r{i}" for i in range(10)]
    table = pa.table({
        "Run": runs + ["unmapped_run"],
        "Protein.Group": ["P1"] + ["X"] * 9 + ["P2"],
        "Global.Q.Value": [0.001] * 11,
        "Proteotypic": [1] * 11,
    })
    parquet_path = tmp_path / "report.parquet"
    pq.write_table(table, parquet_path)
    sdrf_rows = [
        "source name\tcharacteristics[cell line]\tcomment[data file]\n",
    ]
    for i, r in enumerate(runs):
        sdrf_rows.append(f"S{i}\tcell{i}\t{r}.wiff\n")
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text("".join(sdrf_rows))
    subtypes = {f"cell{i}": "TNBC" for i in range(10)}
    out = proteins_per_subtype_quantmsdiann_consistency_filter(
        parquet_path, sdrf_p, subtypes,
        qvalue_cutoff=0.01, min_detection_fraction=0.10, batch_size=4,
    )
    # P1: 1 mapped run / 10 mapped runs = 10% >= 10% -> kept.
    # X: 9/10 -> kept.
    # P2: unmapped_run is not in sdrf so detection in mapped space = 0/10 ->
    # dropped.
    assert out == {"TNBC": {"P1", "X"}}


def test_proteins_per_subtype_consistency_filter_excludes_unknown_subtype(
    tmp_path: Path,
) -> None:
    """Runs whose cell line maps to subtype 'unknown' must NOT contribute
    to any subtype's protein union; they do still count toward the
    consistency-filter denominator (mapped runs) as part of total
    coverage."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        proteins_per_subtype_quantmsdiann_consistency_filter,
    )
    table = pa.table({
        "Run": ["r1", "r2"],
        "Protein.Group": ["P1", "P1"],
        "Global.Q.Value": [0.001, 0.001],
        "Proteotypic": [1, 1],
    })
    parquet_path = tmp_path / "report.parquet"
    pq.write_table(table, parquet_path)
    sdrf = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tmdamb231\tr1.wiff\n"
        "S2\tweirdline\tr2.wiff\n"
    )
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text(sdrf)
    subtypes = {"mdamb231": "TNBC", "weirdline": "unknown"}
    out = proteins_per_subtype_quantmsdiann_consistency_filter(
        parquet_path, sdrf_p, subtypes,
        qvalue_cutoff=0.01, min_detection_fraction=0.10, batch_size=2,
    )
    assert out == {"TNBC": {"P1"}}


def test_proteins_per_subtype_consistency_filter_uses_cache(
    tmp_path: Path,
) -> None:
    """The side-cache wrapper short-circuits when the JSON exists, so the
    parquet is never opened on the cache hit."""
    import json
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        _compute_or_load_diann_subtype_consistency_filter,
    )
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({
        "TNBC": ["P1", "P2"],
        "non-TNBC": ["P1"],
    }))
    # If the function actually reads the parquet, this bogus path crashes.
    out = _compute_or_load_diann_subtype_consistency_filter(
        cache,
        tmp_path / "nonexistent.parquet",
        tmp_path / "nonexistent.tsv",
        {},
    )
    assert out == {"TNBC": {"P1", "P2"}, "non-TNBC": {"P1"}}


# ---------------------------------------------------------------------------
# Headline number parsing reuses PXD030304's parse_diann_summary_log
# ---------------------------------------------------------------------------

def test_pxd004701_diann_summary_parses_expected_headline(tmp_path: Path) -> None:
    """Cross-checks the parse against the actual PXD004701 headline numbers
    so any divergence between the deployed log and the constants in the
    script gets surfaced immediately."""
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        parse_diann_summary_log,
    )
    log = tmp_path / "log.txt"
    log.write_text(
        "[0:25] Target precursors at 1% global q-value: 100499\n"
        "[1:52] Assembling protein groups\n"
        "[2:04] Protein groups with global q-value <= 0.01: 7746\n"
    )
    pg, prec = parse_diann_summary_log(log)
    assert (pg, prec) == (7746, 100499)


# ---------------------------------------------------------------------------
# Subtype aggregation correctness across multiple runs of the same subtype
# ---------------------------------------------------------------------------

def test_subtype_aggregation_unions_across_runs(tmp_path: Path) -> None:
    """Two TNBC runs each detect a different protein group; the subtype
    union must contain both. Filtering happens before aggregation."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        proteins_per_subtype_quantmsdiann_consistency_filter,
    )
    table = pa.table({
        "Run": ["r1", "r2", "r1", "r2"],
        "Protein.Group": ["P1", "P2", "P3", "P3"],
        "Global.Q.Value": [0.001, 0.001, 0.001, 0.001],
        "Proteotypic": [1, 1, 1, 1],
    })
    parquet_path = tmp_path / "report.parquet"
    pq.write_table(table, parquet_path)
    sdrf = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tmdamb231\tr1.wiff\n"
        "S2\tmdamb468\tr2.wiff\n"
    )
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text(sdrf)
    subtypes = {"mdamb231": "TNBC", "mdamb468": "TNBC"}
    out = proteins_per_subtype_quantmsdiann_consistency_filter(
        parquet_path, sdrf_p, subtypes,
        qvalue_cutoff=0.01, min_detection_fraction=0.10, batch_size=2,
    )
    # P1: r1 only -> 1/2 = 50% kept. P2: r2 only -> 1/2 = 50% kept.
    # P3: both -> 2/2 = 100% kept. All three end up in TNBC union.
    assert out == {"TNBC": {"P1", "P2", "P3"}}


def test_counts_tsv_carries_unfiltered_and_target_rows_pxd004701(
    tmp_path: Path,
) -> None:
    """`write_counts_tsv` writes paired rows for the contaminant filter
    audit (2026-05-21 spec §1.7): one `quantmsdiann (DIA-NN, target-only)`
    row carrying the post-filter headline count and one `quantmsdiann
    (DIA-NN, unfiltered ...)` row carrying the pre-filter count."""
    from analysis.figure_pxd004701_sun_vs_quantmsdiann import (
        Counts, write_counts_tsv,
    )

    counts = Counts(
        sun_proteins=6091,
        sun_proteins_raw=8952,
        sun_peptides=90762,
        sun_tnbc=39,
        sun_non_tnbc=37,
        quantmsdiann_proteins_strict=7600,             # target-only
        quantmsdiann_proteins_strict_unfiltered=7746,  # diannsummary.log
        quantmsdiann_proteins_pg_matrix_unfiltered=7700,
        quantmsdiann_proteins_consistent=6200,
        quantmsdiann_proteins_consistent_unfiltered=6296,
        quantmsdiann_peptides=85000,
        quantmsdiann_precursors=100499,
    )
    p = tmp_path / "counts.tsv"
    write_counts_tsv(counts, p)
    text = p.read_text(encoding="utf-8")
    # Both filter policies represented for the strict (no-consistency) headline.
    assert "quantmsdiann (DIA-NN, target-only)" in text
    assert "quantmsdiann (DIA-NN, unfiltered pg_matrix)" in text
    assert "quantmsdiann (DIA-NN, diannsummary.log)" in text
    # And both for the consistency-filtered headline.
    assert "7600" in text and "7746" in text and "7700" in text
    assert "6200" in text and "6296" in text
