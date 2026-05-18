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
