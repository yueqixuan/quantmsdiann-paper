"""Tests for the quantmsdiann vs ProteoBench benchmarks comparison.

Covered: diannsummary.log parser, ProteoBench datapoint JSON parser,
software-name normalisation, dataset->module mapper. No network."""
from __future__ import annotations

import json
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# diannsummary.log parser
# ---------------------------------------------------------------------------

def test_parse_diann_summary_log_extracts_precursors_and_proteins(
    tmp_path: Path,
) -> None:
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        parse_diann_summary_log,
    )
    body = (
        "[0:00] Loading FASTA\n"
        "[0:02] Target precursors at 1% global q-value: 45770\n"
        "[0:07] Protein groups with global q-value <= 0.01: 8227\n"
        "[0:12] Finished\n"
    )
    p = tmp_path / "diannsummary.log"
    p.write_text(body)
    proteins, precursors = parse_diann_summary_log(p)
    assert proteins == 8227
    assert precursors == 45770


def test_parse_diann_summary_log_missing_protein_line_raises(
    tmp_path: Path,
) -> None:
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        parse_diann_summary_log,
    )
    p = tmp_path / "bad.log"
    p.write_text("[0:02] Target precursors at 1% global q-value: 100\n")
    with pytest.raises(ValueError):
        parse_diann_summary_log(p)


def test_parse_diann_summary_log_missing_precursor_line_raises(
    tmp_path: Path,
) -> None:
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        parse_diann_summary_log,
    )
    p = tmp_path / "bad.log"
    p.write_text("[0:07] Protein groups with global q-value <= 0.01: 50\n")
    with pytest.raises(ValueError):
        parse_diann_summary_log(p)


# ---------------------------------------------------------------------------
# Matrix row counter (precursor + protein-group counts across DIA-NN versions)
# ---------------------------------------------------------------------------

def test_count_matrix_data_rows_subtracts_header(tmp_path: Path) -> None:
    """The pr/pg matrix row count is the headline number used uniformly
    across DIA-NN 1.8.1 ... 2.5.0 because the summary-log format is not
    consistent between versions."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        count_matrix_data_rows,
    )
    p = tmp_path / "pr.tsv"
    p.write_text(
        "Protein.Group\tStripped.Sequence\trun1.raw\n"
        "Q1\tAAAR\t1.0\n"
        "Q2\tBBBR\t2.0\n"
        "Q3\tCCCR\t\n"
    )
    assert count_matrix_data_rows(p) == 3


def test_count_matrix_data_rows_empty_file(tmp_path: Path) -> None:
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        count_matrix_data_rows,
    )
    p = tmp_path / "empty.tsv"
    p.write_text("")
    assert count_matrix_data_rows(p) == 0


def test_count_matrix_data_rows_header_only(tmp_path: Path) -> None:
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        count_matrix_data_rows,
    )
    p = tmp_path / "header.tsv"
    p.write_text("Protein.Group\trun1.raw\n")
    assert count_matrix_data_rows(p) == 0


# ---------------------------------------------------------------------------
# Dataset -> ProteoBench module mapping
# ---------------------------------------------------------------------------

def test_dataset_to_module_mapping_covers_all_four_folders() -> None:
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        DATASET_TO_MODULE,
    )
    assert set(DATASET_TO_MODULE) == {
        "PXD049412", "PXD062685", "PXD070049", "ProteoBench_Module_7",
    }
    # Each module entry must carry both a human-readable label and the
    # results-repo name so the fetcher can target the right GitHub repo.
    for folder, info in DATASET_TO_MODULE.items():
        assert "label" in info
        assert "results_repo" in info
        assert info["results_repo"].startswith("Proteobench/Results_quant_")


# ---------------------------------------------------------------------------
# ProteoBench datapoint JSON parser
# ---------------------------------------------------------------------------

def test_parse_proteobench_datapoints_yields_tool_version_precursors(
    tmp_path: Path,
) -> None:
    """ProteoBench results repos store one JSON per submission. We cache the
    full list in a single `<module>.json` file (list-of-dicts). The parser
    yields the trio `(software_name, software_version, nr_prec)`."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        parse_proteobench_datapoints,
    )
    payload = [
        {
            "software_name": "DIA-NN",
            "software_version": "1.8.1",
            "nr_prec": 12345,
            "predictors_library": None,  # empirical
            "id": "DIANN_1",
        },
        {
            "software_name": "AlphaDIA",
            "software_version": "1.10.4",
            "nr_prec": 67890,
            "id": "AlphaDIA_1",
        },
        {
            "software_name": "FragPipe",
            "software_version": "22.0",
            "nr_prec": 54321,
            "id": "FragPipe_1",
        },
    ]
    p = tmp_path / "module.json"
    p.write_text(json.dumps(payload))
    out = list(parse_proteobench_datapoints(p))
    assert out == [
        ("DIA-NN", "1.8.1", 12345, "empirical"),
        ("AlphaDIA", "1.10.4", 67890, "other tool"),
        ("FragPipe", "22.0", 54321, "other tool"),
    ]


def test_classify_predictors_library_categories() -> None:
    """The library classifier groups DIA-NN submissions into the three
    strategies that drive most of the precursor-count variation: empirical
    (None / 'None'), DIANN in-silico predicted, and user-defined speclibs."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        classify_predictors_library,
    )
    assert classify_predictors_library(None) == "empirical"
    assert classify_predictors_library("None") == "empirical"
    assert classify_predictors_library(
        {"RT": "DIANN", "IM": "DIANN", "MS2_int": "DIANN"}
    ) == "predicted (DIANN)"
    assert classify_predictors_library(
        {"RT": "User defined speclib", "IM": "User defined speclib",
         "MS2_int": "User defined speclib"}
    ) == "user-defined speclib"
    # Bare string fallthroughs.
    assert classify_predictors_library("DIANN") == "predicted (DIANN)"
    assert classify_predictors_library("User defined speclib") == "user-defined speclib"
    # Some submissions store predictors_library as a Python repr; the
    # classifier should literal-eval and recurse before falling through.
    assert classify_predictors_library(
        "{'RT': 'DIANN', 'IM': 'DIANN', 'MS2_int': 'DIANN'}"
    ) == "predicted (DIANN)"


def test_parse_proteobench_datapoints_skips_rows_without_nr_prec(
    tmp_path: Path,
) -> None:
    """Datapoints missing nr_prec (rare submission artefact) are skipped, not
    rendered as zero, so they don't drag the y-axis floor downwards."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        parse_proteobench_datapoints,
    )
    payload = [
        {"software_name": "DIA-NN", "software_version": "1.8.1",
         "nr_prec": 100, "id": "ok"},
        {"software_name": "DIA-NN", "software_version": "broken",
         "nr_prec": None, "id": "broken"},
        {"software_name": "DIA-NN", "software_version": "missing",
         "id": "missing-key"},
    ]
    p = tmp_path / "module.json"
    p.write_text(json.dumps(payload))
    out = list(parse_proteobench_datapoints(p))
    assert out == [("DIA-NN", "1.8.1", 100, "empirical")]


# ---------------------------------------------------------------------------
# Software-name normalisation
# ---------------------------------------------------------------------------

def test_normalise_software_name_collapses_dia_nn_variants() -> None:
    """ProteoBench submissions spell DIA-NN as 'DIA-NN', 'DIANN', 'Diann',
    sometimes with trailing spaces. We collapse to 'dia-nn' so the
    highlight-our-tool overlay in the supp figure matches every spelling."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        normalise_software_name,
    )
    assert normalise_software_name("DIA-NN") == "dia-nn"
    assert normalise_software_name("DIANN") == "dia-nn"
    assert normalise_software_name("diann") == "dia-nn"
    assert normalise_software_name("Diann ") == "dia-nn"
    assert normalise_software_name(" DIA-NN ") == "dia-nn"
    # Non-DIA-NN tools are lowercased and stripped but otherwise preserved.
    assert normalise_software_name("AlphaDIA") == "alphadia"
    assert normalise_software_name("FragPipe ") == "fragpipe"
    assert normalise_software_name("Spectronaut") == "spectronaut"


# ---------------------------------------------------------------------------
# Long-format aggregator
# ---------------------------------------------------------------------------

def test_build_long_table_combines_quantmsdiann_and_proteobench(
    tmp_path: Path,
) -> None:
    """The long-format builder should emit one row per (dataset, source,
    tool, version) tuple. quantmsdiann rows carry both precursors and
    proteins; ProteoBench rows carry precursors only (proteins=None)."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        build_long_table,
    )
    quantmsdiann_rows = [
        ("PXD049412", "v1_8_1", 40000, 7500),
        ("PXD049412", "v2_5_0", 45770, 8227),
    ]
    proteobench_rows = {
        "PXD049412": [
            ("DIA-NN", "1.8.1", 41000, "empirical"),
            ("AlphaDIA", "1.10.4", 65666, "other tool"),
        ],
    }
    df = build_long_table(quantmsdiann_rows, proteobench_rows)
    # 2 quantmsdiann rows + 2 ProteoBench rows for this single dataset.
    assert len(df) == 4
    qm = df[df["source"] == "quantmsdiann"]
    pb = df[df["source"] == "proteobench"]
    assert set(qm["version"]) == {"v1_8_1", "v2_5_0"}
    assert qm["proteins"].notna().all()
    assert pb["proteins"].isna().all()
    assert set(pb["tool"]) == {"DIA-NN", "AlphaDIA"}


# ---------------------------------------------------------------------------
# Cache writer (idempotent)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Replicate-threshold extraction (Slack-driven correction)
# ---------------------------------------------------------------------------

def test_extract_nr_prec_at_replicate_threshold_uses_results_bucket() -> None:
    """ProteoBench submissions embed per-replicate-threshold counts under
    `results[str(min_replicates)]['nr_prec']`. The extractor must read from
    that nested bucket, NOT the top-level `nr_prec` (which is always the
    ≥1-replicate count and overstates DIA-NN 1.9.1 per Robbe's note)."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        extract_nr_prec_at_replicate_threshold,
    )
    entry = {
        "software_name": "DIA-NN",
        "software_version": "1.9.1",
        "nr_prec": 105002,  # headline ≥1 count
        "results": {
            "1": {"nr_prec": 105002},
            "2": {"nr_prec": 105002},
            "3": {"nr_prec": 99707},
            "4": {"nr_prec": 91622},
            "5": {"nr_prec": 81223},
            "6": {"nr_prec": 65666},
        },
    }
    assert extract_nr_prec_at_replicate_threshold(entry, 1) == 105002
    assert extract_nr_prec_at_replicate_threshold(entry, 3) == 99707
    assert extract_nr_prec_at_replicate_threshold(entry, 6) == 65666


def test_extract_nr_prec_falls_back_to_top_level_only_for_min1() -> None:
    """If `results` is missing entirely (rare older submission), the
    top-level `nr_prec` may still be present; that's only equivalent to the
    ≥1-replicate threshold, so we must NOT silently fall back for ≥3."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        extract_nr_prec_at_replicate_threshold,
    )
    entry = {"software_name": "DIA-NN", "nr_prec": 50000}
    assert extract_nr_prec_at_replicate_threshold(entry, 1) == 50000
    assert extract_nr_prec_at_replicate_threshold(entry, 3) is None


def test_parse_proteobench_datapoints_at_threshold_skips_missing_buckets(
    tmp_path: Path,
) -> None:
    """The threshold-aware parser should drop any submission whose results
    block does not carry the requested replicate-threshold bucket. The
    library-kind classification still applies."""
    import json as _json
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        parse_proteobench_datapoints_at_threshold,
    )
    payload = [
        {
            "software_name": "DIA-NN", "software_version": "2.3.0",
            "predictors_library": {"RT": "DIANN", "IM": "DIANN",
                                   "MS2_int": "DIANN"},
            "results": {"1": {"nr_prec": 110000}, "3": {"nr_prec": 95000}},
            "id": "A",
        },
        {
            "software_name": "DIA-NN", "software_version": "1.9.1",
            "predictors_library": None,
            "results": {"1": {"nr_prec": 130000}},  # no 3-bucket
            "id": "B",
        },
        {
            "software_name": "AlphaDIA", "software_version": "1.10",
            "results": {"1": {"nr_prec": 80000}, "3": {"nr_prec": 70000}},
            "id": "C",
        },
    ]
    p = tmp_path / "module.json"
    p.write_text(_json.dumps(payload))
    # ≥1 yields all three rows.
    rows1 = list(parse_proteobench_datapoints_at_threshold(p, 1))
    assert [(t, v, n) for t, v, n, _ in rows1] == [
        ("DIA-NN", "2.3.0", 110000),
        ("DIA-NN", "1.9.1", 130000),
        ("AlphaDIA", "1.10", 80000),
    ]
    # ≥3 drops the DIA-NN 1.9.1 submission (no `3` bucket).
    rows3 = list(parse_proteobench_datapoints_at_threshold(p, 3))
    assert [(t, v, n) for t, v, n, _ in rows3] == [
        ("DIA-NN", "2.3.0", 95000),
        ("AlphaDIA", "1.10", 70000),
    ]


# ---------------------------------------------------------------------------
# pr_matrix ≥N-replicate counter
# ---------------------------------------------------------------------------

def test_count_pr_matrix_min_replicates_counts_non_na_threshold(
    tmp_path: Path,
) -> None:
    """The ≥3-replicate count is precursors with non-NA intensity in at least
    3 of the 6 ProteoBench sample columns. Mirrors how ProteoBench's per-
    bucket `nr_prec` is defined."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        count_pr_matrix_min_replicates,
    )
    p = tmp_path / "pr_matrix.tsv"
    # 4 precursors, 6 sample columns ("Condition_A_REP1..3", "B_REP1..3").
    # P1: 6/6  (passes ≥3 and ≥6)
    # P2: 3/6  (passes ≥3, not ≥6)
    # P3: 2/6  (fails ≥3)
    # P4: 0/6  (fails everything; row exists but is unquantified)
    p.write_text(
        "Protein.Group\tStripped.Sequence\tPrecursor.Id\t"
        "A_REP1\tA_REP2\tA_REP3\tB_REP1\tB_REP2\tB_REP3\n"
        "Q1\tAAAR\tAAAR_2\t1.0\t1.0\t1.0\t1.0\t1.0\t1.0\n"
        "Q2\tBBBR\tBBBR_2\t1.0\t\t1.0\t\t1.0\t\n"
        "Q3\tCCCR\tCCCR_2\t1.0\t1.0\t\t\t\t\n"
        "Q4\tDDDR\tDDDR_2\t\t\t\t\t\t\n"
    )
    assert count_pr_matrix_min_replicates(p, 1) == 3  # P1, P2, P3
    assert count_pr_matrix_min_replicates(p, 3) == 2  # P1, P2
    assert count_pr_matrix_min_replicates(p, 6) == 1  # P1


def test_count_pr_matrix_min_replicates_strips_metadata_columns(
    tmp_path: Path,
) -> None:
    """The matrix carries up to 10 DIA-NN metadata columns (Protein.Group,
    Genes, ...). The counter must only look at sample columns, otherwise
    a string-valued 'Genes' column would always count as 'non-NA' and
    inflate the ≥N tally."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        count_pr_matrix_min_replicates,
    )
    p = tmp_path / "pr_matrix.tsv"
    p.write_text(
        "Protein.Group\tProtein.Ids\tProtein.Names\tGenes\t"
        "First.Protein.Description\tProteotypic\tStripped.Sequence\t"
        "Modified.Sequence\tPrecursor.Charge\tPrecursor.Id\t"
        "A_REP1\tA_REP2\tA_REP3\tB_REP1\tB_REP2\tB_REP3\n"
        # Metadata fully populated, samples empty -> should NOT count.
        "Q1\tQ1\tName\tGENE1\tdesc\tTrue\tAAAR\tAAAR\t2\tAAAR_2\t\t\t\t\t\t\n"
    )
    assert count_pr_matrix_min_replicates(p, 1) == 0


def test_consolidate_proteobench_datapoints_writes_list(
    tmp_path: Path,
) -> None:
    """Several per-submission files get consolidated into one list-of-dicts
    JSON file. Re-running over the same inputs must not reorder or
    duplicate entries (idempotence)."""
    from analysis.figure_quantmsdiann_benchmarks_vs_proteobench import (
        consolidate_proteobench_datapoints,
    )
    files = []
    for i, hash_ in enumerate(["aaa", "bbb"]):
        p = tmp_path / f"{hash_}.json"
        # Sort key is the file name, so 'aaa' must come before 'bbb' in the
        # consolidated output regardless of which order files were created.
        p.write_text(json.dumps({"software_name": f"tool_{i}",
                                 "nr_prec": 1000 * (i + 1)}))
        files.append(p)
    dest = tmp_path / "module.json"
    consolidate_proteobench_datapoints(files, dest)
    out = json.loads(dest.read_text())
    assert isinstance(out, list)
    assert len(out) == 2
    # 'aaa' sorts before 'bbb', so tool_0 must be first.
    assert out[0]["software_name"] == "tool_0"
    assert out[1]["software_name"] == "tool_1"
