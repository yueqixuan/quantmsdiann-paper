"""Tests for the PXD030304 (ProCan-DepMapSanger) reanalysis comparison
parsers. The two pipelines being compared are the published 2022 ProCan
analysis (Gonçalves et al. Cancer Cell, figshare deposit 19345397) and the
quantmsdiann DIA-NN reanalysis (PRIDE quantms-collections).
"""
from __future__ import annotations

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Figshare mapping file
# ---------------------------------------------------------------------------

def test_parse_procan_mapping_returns_cell_line_to_tissue(tmp_path: Path) -> None:
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        parse_procan_mapping,
    )
    body = (
        "Cell_line\tSIDM\tProject_Identifier\tTissue_type\tCancer_type\tCancer_subtype\n"
        "BC-1\tSIDM00896\tSIDM00896;BC-1\tHaematopoietic and Lymphoid\tBCL\t\n"
        "L-363\tSIDM00312\tSIDM00312;L-363\tHaematopoietic and Lymphoid\tPCM\t\n"
        "HCT-116\tSIDM00783\tSIDM00783;HCT-116\tLarge Intestine\tColon cancer\t\n"
    )
    p = tmp_path / "mapping.txt"
    p.write_text(body)
    out = parse_procan_mapping(p)
    assert out == {
        "BC-1": "Haematopoietic and Lymphoid",
        "L-363": "Haematopoietic and Lymphoid",
        "HCT-116": "Large Intestine",
    }


def test_parse_procan_mapping_raises_on_missing_columns(tmp_path: Path) -> None:
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        parse_procan_mapping,
    )
    p = tmp_path / "bad.txt"
    p.write_text("Cell_line\tFoo\nBC-1\tbar\n")
    with pytest.raises(ValueError):
        parse_procan_mapping(p)


# ---------------------------------------------------------------------------
# ProCan protein matrix -> per-tissue protein sets
# ---------------------------------------------------------------------------

def test_parse_procan_replicates_mapping_excludes_hek293t(tmp_path: Path) -> None:
    """`mapping_file_replicates.txt` carries the per-MS-run cell-line/tissue
    mapping including HEK293T QC runs. The default loader drops them so
    per-tissue counters never see QC controls."""
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        parse_procan_replicates_mapping,
    )
    body = (
        "Automatic_MS_filename\tBatch\tDate\tInstrument\tCell_line\tSIDM"
        "\tProject_Identifier\tTissue_type\tCancer_type\tCancer_subtype\n"
        "run_a\tP01\t1/1/18\tM01\tBC-1\tSIDM00896\tSIDM00896;BC-1"
        "\tHaematopoietic and Lymphoid\tBCL\t\n"
        "run_b\tP01\t1/1/18\tM01\tControl_HEK293T_std_H002\tControl_HEK293T"
        "\tControl_HEK293T;Control_HEK293T_std_H002\tControl_HEK293T\tControl_HEK293T\t\n"
        "run_c\tP01\t1/1/18\tM01\tHCT-116\tSIDM00783\tSIDM00783;HCT-116"
        "\tLarge Intestine\tColon\t\n"
    )
    p = tmp_path / "map_rep.txt"
    p.write_text(body)
    # Default: HEK293T excluded.
    out = parse_procan_replicates_mapping(p)
    assert out == {
        "run_a": "Haematopoietic and Lymphoid",
        "run_c": "Large Intestine",
    }
    # Explicit opt-in keeps the QC row.
    out_with_qc = parse_procan_replicates_mapping(p, exclude_hek293t=False)
    assert out_with_qc == {
        "run_a": "Haematopoietic and Lymphoid",
        "run_b": "Control_HEK293T",
        "run_c": "Large Intestine",
    }


def test_proteins_per_tissue_procan_from_replicates_matrix(tmp_path: Path) -> None:
    """Each row of `protein_matrix_8498_replicates.txt` is one MS run; columns
    are protein IDs (`<accession>;<name>`). The header's first cell is
    empty in the real file. Per tissue we take the union across all MS runs
    mapped to that tissue; HEK293T runs are excluded via the mapping
    file."""
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        proteins_per_tissue_procan,
    )
    # Header: empty first cell, then 3 protein cols (matches the real file's
    # unnamed first-column header).
    matrix = (
        "\tP1;A_HUMAN\tP2;B_HUMAN\tP3;C_HUMAN\n"
        # BC-1 run, detects P1, P2.
        "run_a\t1.0\t2.0\t\n"
        # HCT-116 run #1, detects P1 only.
        "run_c1\t0.5\t\t\n"
        # HCT-116 run #2, detects P3 only -> union with run_c1 -> {P1, P3}.
        "run_c2\t\t\t9.9\n"
        # HEK293T run, must be ignored.
        "run_b\t1.0\t1.0\t1.0\n"
    )
    matrix_p = tmp_path / "matrix.txt"
    matrix_p.write_text(matrix)
    mapping = (
        "Automatic_MS_filename\tBatch\tDate\tInstrument\tCell_line\tSIDM"
        "\tProject_Identifier\tTissue_type\tCancer_type\tCancer_subtype\n"
        "run_a\tP01\t1/1/18\tM01\tBC-1\tSIDM00896\tSIDM00896;BC-1"
        "\tHaematopoietic and Lymphoid\tBCL\t\n"
        "run_c1\tP01\t1/1/18\tM01\tHCT-116\tSIDM00783\tSIDM00783;HCT-116"
        "\tLarge Intestine\tColon\t\n"
        "run_c2\tP01\t1/1/18\tM01\tHCT-116\tSIDM00783\tSIDM00783;HCT-116"
        "\tLarge Intestine\tColon\t\n"
        "run_b\tP01\t1/1/18\tM01\tControl_HEK293T_std_H002\tControl_HEK293T"
        "\tControl_HEK293T;Control_HEK293T_std_H002\tControl_HEK293T\tControl_HEK293T\t\n"
    )
    map_p = tmp_path / "mapping.txt"
    map_p.write_text(mapping)
    out = proteins_per_tissue_procan(matrix_p, map_p, chunksize=2)
    assert out == {
        "Haematopoietic and Lymphoid": {"P1;A_HUMAN", "P2;B_HUMAN"},
        "Large Intestine": {"P1;A_HUMAN", "P3;C_HUMAN"},
    }


# ---------------------------------------------------------------------------
# quantmsdiann pg_matrix -> per-tissue protein sets
# ---------------------------------------------------------------------------

def test_proteins_per_tissue_quantmsdiann_uses_sdrf_and_procan_mapping(
    tmp_path: Path,
) -> None:
    """DIA-NN's pr_matrix has 10 metadata cols then per-run quant cols. SDRF
    maps each run (`comment[data file]` rewritten .wiff -> .mzML) to a cell
    line; ProCan mapping gives cell line -> tissue. A protein group counts
    for a tissue if any of its precursor rows is non-NA in any run belonging
    to that tissue."""
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        proteins_per_tissue_quantmsdiann,
    )
    pr_matrix = (
        "Protein.Group\tProtein.Ids\tProtein.Names\tGenes"
        "\tFirst.Protein.Description\tProteotypic\tStripped.Sequence"
        "\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id"
        "\trun1.mzML\trun2.mzML\trun3.mzML\trun4.mzML\n"
        # Q1 precursor in run1 (BC-1) and run3 (HCT-116).
        "Q1\tQ1\tA\tA\td\t1\tAAAR\tAAAR\t2\tAAAR2\t1.0\t\t2.0\t\n"
        # Q2 precursor only in run2 (BC-1).
        "Q2\tQ2\tB\tB\td\t1\tBBBR\tBBBR\t2\tBBBR2\t\t1.0\t\t\n"
        # Q3 precursor only in run4 (HCT-116).
        "Q3\tQ3\tC\tC\td\t1\tCCCR\tCCCR\t2\tCCCR2\t\t\t\t3.0\n"
        # Q4 has TWO precursor rows; first in run1 (BC-1), second in run3
        # (HCT-116). Tests that union across rows of same Protein.Group works.
        "Q4\tQ4\tD\tD\td\t1\tDDDR\tDDDR\t2\tDDDR2\t1.0\t\t\t\n"
        "Q4\tQ4\tD\tD\td\t1\tEEER\tEEER\t2\tEEER2\t\t\t1.0\t\n"
    )
    pg_p = tmp_path / "pr.tsv"
    pg_p.write_text(pr_matrix)
    sdrf = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tBC-1\trun1.wiff\n"
        "S2\tBC-1\trun2.wiff\n"
        "S3\tHCT-116\trun3.wiff\n"
        "S4\tHCT-116\trun4.wiff\n"
    )
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text(sdrf)
    mapping = (
        "Cell_line\tSIDM\tProject_Identifier\tTissue_type\tCancer_type\tCancer_subtype\n"
        "BC-1\tSIDM00896\tSIDM00896;BC-1\tHaematopoietic and Lymphoid\t\t\n"
        "HCT-116\tSIDM00783\tSIDM00783;HCT-116\tLarge Intestine\t\t\n"
    )
    map_p = tmp_path / "mapping.txt"
    map_p.write_text(mapping)
    out = proteins_per_tissue_quantmsdiann(pg_p, sdrf_p, map_p, chunksize=2)
    assert out == {
        "Haematopoietic and Lymphoid": {"Q1", "Q2", "Q4"},
        "Large Intestine": {"Q1", "Q3", "Q4"},
    }


# ---------------------------------------------------------------------------
# Per-run completeness (ProCan side)
# ---------------------------------------------------------------------------

def test_procan_per_run_completeness_from_peptide_counts(tmp_path: Path) -> None:
    """`peptide_counts_per_protein_per_sample.txt`: 1 row per MS run, 1 col
    per protein, cell value = number of peptides identified for that protein
    in that run. Per-run completeness = fraction of proteins with peptide
    count > 0 (denominator = total proteins in the matrix)."""
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        per_run_completeness_procan,
    )
    body = (
        "Run\tP1\tP2\tP3\tP4\n"
        # run1: 3/4 proteins have at least one peptide.
        "run1\t1\t2\t0\t5\n"
        # run2: 1/4 proteins.
        "run2\t0\t\t0\t1\n"
        # run3: 4/4 proteins.
        "run3\t3\t2\t1\t1\n"
    )
    p = tmp_path / "peptide_counts.txt"
    p.write_text(body)
    out = per_run_completeness_procan(p)
    assert out == {
        "run1": pytest.approx(3 / 4),
        "run2": pytest.approx(1 / 4),
        "run3": pytest.approx(4 / 4),
    }


# ---------------------------------------------------------------------------
# Per-run completeness (quantmsdiann side, from pg_matrix)
# ---------------------------------------------------------------------------

def test_quantmsdiann_per_run_completeness_from_pg_matrix(tmp_path: Path) -> None:
    """For each per-run column in pg_matrix.tsv, the fraction of protein-group
    rows that are non-NA. Mirrors ProCan's per-run completeness denominator
    (total proteins in the pipeline's matrix)."""
    from analysis.figure_pxd030304_procan_vs_quantmsdiann import (
        per_run_completeness_quantmsdiann,
    )
    pg_matrix = (
        "Protein.Group\tProtein.Names\tGenes"
        "\tFirst.Protein.Description\tN.Sequences\tN.Proteotypic.Sequences"
        "\trun1.mzML\trun2.mzML\trun3.mzML\n"
        "Q1\tA\tA\td\t1\t1\t1.0\t\t2.0\n"
        "Q2\tB\tB\td\t1\t1\t\t\t1.0\n"
        "Q3\tC\tC\td\t1\t1\t1.0\t1.0\t1.0\n"
        "Q4\tD\tD\td\t1\t1\t\t\t\n"
    )
    p = tmp_path / "pg.tsv"
    p.write_text(pg_matrix)
    out = per_run_completeness_quantmsdiann(p)
    # 4 protein groups total. run1: 2/4, run2: 1/4, run3: 3/4.
    assert out == {
        "run1.mzML": pytest.approx(2 / 4),
        "run2.mzML": pytest.approx(1 / 4),
        "run3.mzML": pytest.approx(3 / 4),
    }
