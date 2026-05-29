"""Tests for `analysis.proteobench_metrics` — the pr_matrix → ProteoBench
metric pipeline.

No-network tests only. The heavy end-to-end path
(`compute_proteobench_metrics`) is exercised by the cache-population
step at script time; here we test the unit helpers in isolation."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


def test_parse_accessions_from_fasta_returns_primary_only(tmp_path: Path) -> None:
    """SwissProt headers are `>sp|ACC|NAME_SPECIES ...`. The accession
    parser must return only the primary accession token, ignoring the
    species suffix on the entry name."""
    from analysis.proteobench_metrics import _parse_accessions_from_fasta
    f = tmp_path / "fixture.fasta"
    f.write_text(
        ">sp|Q96P70|FYV1_HUMAN cytochrome OS=Homo sapiens OX=9606 GN=FYV1\n"
        "MAAAA\n"
        ">sp|P00330|ADH1_YEAST Alcohol dehydrogenase 1 OS=S288c OX=559292\n"
        "MSAAA\n"
        "AGGGG\n"
        ">tr|A0A0B7|EXAMPLE_OTHER trembl line that should also parse\n"
        "MQQQQ\n"
    )
    acc = _parse_accessions_from_fasta(f)
    assert acc == {"Q96P70", "P00330", "A0A0B7"}


def test_annotate_species_suffix_handles_multi_accession_and_missing() -> None:
    """quantmsdiann's `Protein.Ids` cell can be `;`-separated and may
    include accessions not in the species map (contaminants, decoys,
    proteins from FASTAs we don't carry). The annotator must:
      1. add the suffix to every known accession,
      2. pass through unknown accessions unchanged,
      3. preserve the ordering and `;` separator.
    """
    from analysis.proteobench_metrics import annotate_species_suffix
    species_map = {"Q96P70": "HUMAN", "P00330": "YEAST", "P0A6F5": "ECOLI"}
    assert annotate_species_suffix(
        "Q96P70;P09417;P00330", species_map,
    ) == "Q96P70_HUMAN;P09417;P00330_YEAST"
    assert annotate_species_suffix(
        "P0A6F5", species_map,
    ) == "P0A6F5_ECOLI"
    # Pass-through for unknown / empty.
    assert annotate_species_suffix("FOO_CONTAM", species_map) == "FOO_CONTAM"
    assert annotate_species_suffix("", species_map) == ""
    assert annotate_species_suffix("nan", species_map) == "nan"


def test_melt_pr_matrix_drops_nan_intensities_and_strips_extension(
    tmp_path: Path,
) -> None:
    """The wide-to-long melt must:
      1. produce one row per (precursor, run) pair where the intensity
         is non-NaN,
      2. drop precursors that weren't quantified in a given run,
      3. strip `.raw` / `.d` / `.mzML` suffixes from the `Run` column so
         it matches ProteoBench's condition_mapper keys.
    """
    from analysis.proteobench_metrics import melt_pr_matrix
    p = tmp_path / "pr_matrix.tsv"
    # 2 precursors × 3 runs, one cell missing.
    p.write_text(
        "Protein.Group\tProtein.Ids\tProtein.Names\tGenes\t"
        "First.Protein.Description\tProteotypic\tStripped.Sequence\t"
        "Modified.Sequence\tPrecursor.Charge\tPrecursor.Id\t"
        "RunA.raw\tRunB.d\tRunC.mzML\n"
        "Q96P70\tQ96P70\tFYV1\tFYV1\tdesc\t1\tAAAAR\t(UniMod:1)AAAAR\t2\t"
        "(UniMod:1)AAAAR2\t1234.5\t\t2345.6\n"
        "P09417\tP09417\tQR1\tQR1\tdesc\t1\tBBBBR\tBBBBR\t3\tBBBBR3\t"
        "5678.9\t6789.0\t7890.1\n"
    )
    long_df = melt_pr_matrix(p)
    assert set(long_df.columns) == {
        "Modified.Sequence", "Protein.Ids", "Precursor.Charge",
        "Run", "Precursor.Normalised",
    }
    # 2 precursors × 3 runs = 6, minus 1 NaN cell = 5 rows.
    assert len(long_df) == 5
    # All run extensions stripped.
    assert set(long_df["Run"]) == {"RunA", "RunB", "RunC"}
    # Intensities are floats.
    assert long_df["Precursor.Normalised"].dtype.kind == "f"


def test_iter_metric_rows_emits_only_present_metrics() -> None:
    """Long-format extractor must skip metrics that are None and skip
    `results` keys that don't parse as integer thresholds (defensive
    against any extra keys in the cached payload)."""
    from analysis.proteobench_metrics import iter_metric_rows
    payload = {
        "results": {
            "1": {
                "nr_prec": 1000,
                "median_abs_epsilon_global": 0.3,
                "roc_auc": None,  # should be skipped
            },
            "3": {
                "nr_prec": 800,
                "median_abs_epsilon_global": 0.25,
            },
            "junk": {"nr_prec": 1},  # should be skipped (non-int key)
        },
    }
    rows = iter_metric_rows([("PXD049412", "v2_5_0", payload)])
    triples = {(r["threshold"], r["metric"]) for r in rows}
    assert (1, "nr_prec") in triples
    assert (1, "median_abs_epsilon_global") in triples
    assert (1, "roc_auc") not in triples  # was None
    assert (3, "nr_prec") in triples
    assert (3, "median_abs_epsilon_global") in triples
    # Non-int threshold rejected.
    assert not any(r["threshold"] == "junk" for r in rows)


def test_metrics_cache_path_isolates_dataset_and_version(tmp_path: Path) -> None:
    """Cache layout is one JSON per (dataset, version). The path
    builder must keep them distinct so a v1.8.1 run doesn't clobber a
    v2.5.0 cache for the same dataset."""
    from analysis.proteobench_metrics import metrics_cache_path
    p1 = metrics_cache_path("PXD049412", "v1_8_1")
    p2 = metrics_cache_path("PXD049412", "v2_5_0")
    assert p1 != p2
    assert "PXD049412_v1_8_1" in p1.name
    assert "PXD049412_v2_5_0" in p2.name
