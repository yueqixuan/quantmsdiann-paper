from pathlib import Path
import textwrap

import pytest

from analysis.venn_protein_accessions import (
    accessions_with_min_peptides_diann,
    accessions_with_min_peptides_openswath,
    extract_accessions_diann,
    extract_accessions_openswath,
)


def test_extract_accessions_diann_single_bare() -> None:
    assert extract_accessions_diann("P12345") == {"P12345"}


def test_extract_accessions_diann_semicolon_list() -> None:
    assert extract_accessions_diann("P12345;P67890") == {"P12345", "P67890"}


def test_extract_accessions_diann_sp_header_form() -> None:
    assert extract_accessions_diann("sp|P12345|HUMAN") == {"P12345"}


def test_extract_accessions_diann_isoform_suffix_stripped() -> None:
    assert extract_accessions_diann("P12345-2") == {"P12345"}
    assert extract_accessions_diann("P12345-2;Q67890-1") == {"P12345", "Q67890"}


def test_extract_accessions_diann_keeps_contam_and_entrap() -> None:
    # Strip CONTAM_/ENTRAP_ prefix so the underlying UniProt accession matches
    # across pipelines; downstream callers can filter if they want.
    assert extract_accessions_diann("CONTAM_P02768-1;P02768") == {"P02768"}
    assert extract_accessions_diann("ENTRAP_Q8NFD5;Q8NFD5") == {"Q8NFD5"}


def test_extract_accessions_diann_empty_or_none() -> None:
    assert extract_accessions_diann("") == set()
    assert extract_accessions_diann(None) == set()  # type: ignore[arg-type]


def test_extract_accessions_openswath_single_proteotypic() -> None:
    assert extract_accessions_openswath("1/sp|P12345|HUMAN") == {"P12345"}


def test_extract_accessions_openswath_multimap() -> None:
    assert extract_accessions_openswath(
        "2/sp|P12345|HUMAN/sp|P67890|HUMAN"
    ) == {"P12345", "P67890"}


def test_extract_accessions_openswath_empty_or_none() -> None:
    assert extract_accessions_openswath("") == set()
    assert extract_accessions_openswath(None) == set()  # type: ignore[arg-type]


def test_extract_accessions_openswath_decoy_yields_accession() -> None:
    # The decoy filter happens upstream; the parser just extracts accessions.
    assert extract_accessions_openswath(
        "DECOY_1/sp|P50990|TCPQ_HUMAN"
    ) == {"P50990"}


def test_accessions_with_min_peptides_diann_filters_and_expands(
    tmp_path: Path,
) -> None:
    body = textwrap.dedent(
        """
        Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tProteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id\tRun_A\tRun_B
        P1\tP1\tA\tG1\td\t1\tAAAR\tAAAR\t2\tAAAR2\t10\t
        P1\tP1\tA\tG1\td\t1\tBBBR\tBBBR\t2\tBBBR2\t\t20
        P2;P2alt\tP2\tB\tG2\td\t1\tDDDR\tDDDR\t2\tDDDR2\t50\t
        P2;P2alt\tP2\tB\tG2\td\t1\tEEER\tEEER\t2\tEEER2\t60\t
        P3\tP3\tC\tG3\td\t1\tFFFR\tFFFR\t2\tFFFR2\t70\t
        CONTAM_P4-1;P4\tP4\tD\tG4\td\t1\tGGGR\tGGGR\t2\tGGGR2\t80\t
        CONTAM_P4-1;P4\tP4\tD\tG4\td\t1\tHHHR\tHHHR\t2\tHHHR2\t90\t
        """
    ).lstrip("\n")
    p = tmp_path / "pr.tsv"
    p.write_text(body)
    out = accessions_with_min_peptides_diann(p, min_peptides=2)
    # P1 has {AAAR, BBBR} -> kept; expands to {P1}
    # P2;P2alt has {DDDR, EEER} -> kept; expands to {P2, P2alt}
    # P3 has {FFFR} -> dropped
    # CONTAM_P4-1;P4 has {GGGR, HHHR} -> kept; expands to {P4}
    assert out == {"P1", "P2", "P2alt", "P4"}


def test_accessions_with_min_peptides_openswath_filters_and_expands(
    tmp_path: Path,
) -> None:
    body = (
        "Peptide\tProtein\tIntensity_run0\tRT_run0\tscore_run0\tIntensity_run1\tRT_run1\tscore_run1\n"
        "1_AAAR_2_run0\t1/sp|P1|HUMAN\t100\t10\t0.005\t150\t10.2\t2.0\n"
        "2_BBBR_2_run0\t1/sp|P1|HUMAN\t110\t10\t0.005\t160\t10.3\t2.0\n"
        "3_CCCR_2_run0\t1/sp|P2|HUMAN\t200\t12\t0.001\t250\t12.1\t0.008\n"
        "4_DDDR_2_run0\t1/sp|P3|HUMAN\t300\t9\t0.001\t300\t9.1\t0.001\n"
        "5_EEER_2_run0\t1/sp|P3|HUMAN\t400\t9\t0.001\t400\t9.1\t0.001\n"
        "DECOY_6_FFFR_2_run0\t1/sp|P4|HUMAN\t500\t9\t0.001\t500\t9.1\t0.001\n"
    )
    p = tmp_path / "fa.tsv"
    p.write_text(body)
    out = accessions_with_min_peptides_openswath(p, min_peptides=2)
    # P1: 2 peptides -> kept; P2: 1 peptide -> dropped; P3: 2 peptides -> kept; decoy P4 excluded
    assert out == {"P1", "P3"}
