from pathlib import Path
import textwrap

import pytest

from analysis.figure_original_vs_quantmsdiann import count_quantified_rows


PR_METADATA = [
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
    "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
]
PG_METADATA = [
    "Protein.Group", "Protein.Names", "Genes",
    "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences",
]


def write_matrix(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "matrix.tsv"
    p.write_text(textwrap.dedent(body).lstrip("\n"))
    return p


def test_counts_pg_matrix_rows_with_at_least_one_non_na(tmp_path: Path) -> None:
    # pg_matrix-style layout: 6 metadata columns then 3 sample columns.
    matrix = write_matrix(
        tmp_path,
        """
        Protein.Group\tProtein.Names\tGenes\tFirst.Protein.Description\tN.Sequences\tN.Proteotypic.Sequences\tRun_A\tRun_B\tRun_C
        P1\tA\tGENE1\tdesc1\t3\t3\t10\t\t20
        P2\tB\tGENE2\tdesc2\t1\t1\t\t\t
        P3\tC\tGENE3\tdesc3\t2\t2\t\tNA\t
        P4\tD\tGENE4\tdesc4\t5\t5\t\t5\t
        """,
    )

    # P1 and P4 each have one real value; P2 is all empty; P3 has only "NA".
    assert count_quantified_rows(matrix, PG_METADATA) == 2


def test_counts_unique_peptides_in_pr_matrix(tmp_path: Path) -> None:
    # pr_matrix-style layout: 10 metadata columns then 3 sample columns.
    # Two precursor rows share Stripped.Sequence "AAAR" (different charges);
    # one of them has a quant, the other doesn't. Expected unique-peptide
    # count = 2 (AAAR, BBBR).
    matrix = write_matrix(
        tmp_path,
        """
        Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tProteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id\tRun_A\tRun_B\tRun_C
        P1\tP1\tA\tGENE1\tdesc1\t1\tAAAR\tAAAR\t2\tAAAR2\t100\t\t
        P1\tP1\tA\tGENE1\tdesc1\t1\tAAAR\tAAAR\t3\tAAAR3\t\t\t
        P2\tP2\tB\tGENE2\tdesc2\t1\tBBBR\tBBBR\t2\tBBBR2\t\t200\t
        P3\tP3\tC\tGENE3\tdesc3\t1\tCCCR\tCCCR\t2\tCCCR2\t\t\t
        """,
    )

    # Quantified precursor rows: AAAR2, BBBR2. Unique peptides: {AAAR, BBBR}.
    assert count_quantified_rows(matrix, PR_METADATA,
                                 unique_by="Stripped.Sequence") == 2

    # And without unique_by, we should count rows: 2.
    assert count_quantified_rows(matrix, PR_METADATA) == 2


def test_raises_if_metadata_column_missing(tmp_path: Path) -> None:
    matrix = write_matrix(
        tmp_path,
        """
        Protein.Group\tProtein.Names\tGenes\tN.Sequences\tN.Proteotypic.Sequences\tRun_A
        P1\tA\tGENE1\t3\t3\t10
        """,
    )

    with pytest.raises(ValueError, match="First.Protein.Description"):
        count_quantified_rows(matrix, PG_METADATA)


def test_count_openswath_quantified_excludes_decoys(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import count_openswath_quantified
    body = """\
Peptide\tProtein\tIntensity_run0\tRT_run0\tscore_run0\tIntensity_run1\tRT_run1\tscore_run1
1_AAAR_2_run0\t1/sp|P1|HUMAN\t100\t10.0\t0.5\t\t\t
2_BBBR_2_run0\t1/sp|P2|HUMAN\t\t\t\t200\t12.0\t0.4
3_CCCR_2_run0\t1/sp|P3|HUMAN\t\t\t\t\t\t
DECOY_4_DDDR_2_run0\t1/sp|P4|HUMAN\t500\t11.0\t0.9\t\t\t
5_EEER_2_run0\tDECOY_1/sp|P5|HUMAN\t600\t11.5\t0.95\t\t\t
6_FFFR_2_run0\t2/sp|P6|HUMAN/sp|P6alt|HUMAN\t700\t9.5\t0.3\t\t\t
"""
    p = tmp_path / "fa.tsv"
    p.write_text(body)
    precursors, peptides, proteins = count_openswath_quantified(p)
    # Quantified target rows: AAAR, BBBR, FFFR. CCCR has no intensity. DECOYs excluded.
    assert precursors == 3
    assert peptides == 3
    # Unique Protein values among quantified targets:
    # "1/sp|P1|HUMAN", "1/sp|P2|HUMAN", "2/sp|P6|HUMAN/sp|P6alt|HUMAN" -> 3
    assert proteins == 3


def test_count_openswath_quantified_raises_on_missing_columns(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import count_openswath_quantified
    p = tmp_path / "bad.tsv"
    p.write_text("Peptide\tProtein\tRT_run0\tscore_run0\nx\ty\t1\t2\n")
    with pytest.raises(ValueError, match="Intensity"):
        count_openswath_quantified(p)  # return value not unpacked; just check it raises


def test_parse_summary_log_finds_protein_total(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import parse_summary_log
    p = tmp_path / "log.txt"
    p.write_text(
        "[0:01] Spectral library loaded\n"
        "[1:31] Protein groups with global q-value <= 0.01: 6927\n"
        "[2:10] Compressed report saved\n"
    )
    assert parse_summary_log(p) == 6927


def test_parse_summary_log_raises_if_line_missing(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import parse_summary_log
    p = tmp_path / "log.txt"
    p.write_text("nothing interesting here\n")
    with pytest.raises(ValueError, match="global q-value"):
        parse_summary_log(p)


def test_count_eprot73_genes_counts_unique_ensembl_ids(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import count_eprot73_genes
    # Mimic the real EA file layout: header on line 1, data from line 2.
    # No preamble — the real file starts directly with the "Gene ID" header.
    body = (
        "Gene ID\tGene Name\tg1.WithInSampleAbundance\n"
        "ENSG00000000001\tFOO\t1.0\n"
        "ENSG00000000002\tBAR\t\n"
        "ENSG00000000001\tFOO\t2.0\n"  # duplicate gene id, distinct row
        "\tEmptyId\t0.5\n"
        "NOTANENSG\tX\t1.0\n"
    )
    p = tmp_path / "e.tsv"
    p.write_text(body)
    # Expect 2 unique Ensembl gene IDs (00000001 and 00000002).
    assert count_eprot73_genes(p) == 2


def test_count_openswath_quantified_handles_missing_peptide_with_decoy_protein(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import count_openswath_quantified
    # An empty Peptide cell with a DECOY Protein must still be excluded.
    body = (
        "Peptide\tProtein\tIntensity_run0\n"
        "1_AAAR_2_run0\t1/sp|P1|HUMAN\t100\n"
        "\tDECOY_1/sp|P2|HUMAN\t999\n"
    )
    p = tmp_path / "fa.tsv"
    p.write_text(body)
    out = count_openswath_quantified(p)
    # Just one target precursor (AAAR), one target protein.
    assert out[0] == 1
    assert out[-1] == 1  # use [-1] so this test still works once a peptide field is added


def test_per_run_non_na_fraction_diann(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import per_run_non_na_fraction_diann
    matrix = write_matrix(
        tmp_path,
        """
        Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tProteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id\tRun_A\tRun_B
        P1\tP1\tA\tG1\td\t1\tAAAR\tAAAR\t2\tAAAR2\t10\t
        P2\tP2\tB\tG2\td\t1\tBBBR\tBBBR\t2\tBBBR2\t20\t30
        P3\tP3\tC\tG3\td\t1\tCCCR\tCCCR\t2\tCCCR2\t\t40
        """,
    )
    out = per_run_non_na_fraction_diann(matrix)
    # 3 rows total; Run_A has 2 non-NA, Run_B has 2 non-NA.
    assert out == {"Run_A": pytest.approx(2 / 3), "Run_B": pytest.approx(2 / 3)}


def test_per_run_non_na_fraction_openswath_excludes_decoys(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import per_run_non_na_fraction_openswath
    body = (
        "Peptide\tProtein\tIntensity_run0\tRT_run0\tscore_run0\tIntensity_run1\tRT_run1\tscore_run1\n"
        "1_AAAR_2_run0\t1/sp|P1|HUMAN\t100\t10\t.5\t\t\t\n"
        "2_BBBR_2_run0\t1/sp|P2|HUMAN\t\t\t\t200\t12\t.4\n"
        "DECOY_3_CCCR_2_run0\t1/sp|P3|HUMAN\t999\t1\t.9\t999\t1\t.9\n"
    )
    p = tmp_path / "fa.tsv"
    p.write_text(body)
    out = per_run_non_na_fraction_openswath(p)
    # 2 target rows; run0 has 1 non-NA, run1 has 1 non-NA among targets.
    assert out == {"Intensity_run0": pytest.approx(1 / 2), "Intensity_run1": pytest.approx(1 / 2)}
