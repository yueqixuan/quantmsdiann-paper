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


def test_per_run_real_detection_fraction_openswath_uses_score(tmp_path: Path) -> None:
    from analysis.figure_original_vs_quantmsdiann import per_run_real_detection_fraction_openswath
    body = (
        "Peptide\tProtein\tIntensity_run0\tRT_run0\tscore_run0\tIntensity_run1\tRT_run1\tscore_run1\n"
        # Row 1: real detection in run0 (score 0.005), requant in run1 (score 2.0).
        "1_AAAR_2_run0\t1/sp|P1|HUMAN\t100\t10\t0.005\t150\t10.2\t2.0\n"
        # Row 2: real detection in both runs.
        "2_BBBR_2_run0\t1/sp|P2|HUMAN\t200\t12\t0.001\t250\t12.1\t0.008\n"
        # Row 3: requant in both runs (score 2.0).
        "3_CCCR_2_run0\t1/sp|P3|HUMAN\t300\t9\t2.0\t300\t9.1\t2.0\n"
        # Row 4: decoy (Peptide prefix); should be excluded entirely.
        "DECOY_4_DDDR_2_run0\t1/sp|P4|HUMAN\t999\t1\t0.001\t999\t1\t0.001\n"
    )
    p = tmp_path / "fa.tsv"
    p.write_text(body)
    out = per_run_real_detection_fraction_openswath(p)
    # 3 target rows. score_run0 has 2 detections (rows 1,2); score_run1 has 1 (row 2).
    assert out == {"score_run0": pytest.approx(2 / 3), "score_run1": pytest.approx(1 / 3)}


def test_per_run_real_detection_fraction_diann_parquet(tmp_path: Path) -> None:
    """Per-run completeness from the DIA-NN long-format report.

    Numerator per run: distinct Precursor.Id rows with Q.Value <= 0.01 in that
    run AND Global.Q.Value <= 0.01 (so the per-cell filter is strict 1% per-run
    precursor FDR, matching the OpenSWATH score <= 0.01 criterion).
    Denominator: distinct Precursor.Id with Global.Q.Value <= 0.01 anywhere."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from analysis.figure_original_vs_quantmsdiann import (
        per_run_real_detection_fraction_diann_parquet,
    )
    table = pa.table({
        "Run":            ["A",   "A",   "A",   "B",   "B",   "B",   "C"],
        "Precursor.Id":   ["X2",  "Y2",  "Z2",  "X2",  "Y2",  "Z2",  "X2"],
        "Q.Value":        [0.005, 0.02,  0.001, 0.001, 0.5,   0.005, 0.005],
        "Global.Q.Value": [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.5],
    })
    p = tmp_path / "report.parquet"
    pq.write_table(table, p)
    # Denominator: distinct Precursor.Id with Global.Q.Value <= 0.01 = {X2, Y2, Z2} = 3.
    # Run A: X2 (0.005 OK), Z2 (0.001 OK); Y2 (0.02) fails per-run. -> 2/3
    # Run B: X2 OK, Z2 OK; Y2 fails per-run. -> 2/3
    # Run C: only X2 row, but its Global.Q.Value=0.5 means X2 is NOT in the denominator's pool
    #        from this run's perspective... no, denominator is global. X2 IS in the pool
    #        (it appears in runs A,B with Global=0.001). But this particular C-row has
    #        Global=0.5 -> filtered out, so run C has 0 detections. -> 0/3
    out = per_run_real_detection_fraction_diann_parquet(p)
    assert out == {
        "A": pytest.approx(2 / 3),
        "B": pytest.approx(2 / 3),
        "C": pytest.approx(0 / 3),
    }


def test_unique_peptides_per_protein_diann_counts_proteotypic_quantified(
    tmp_path: Path,
) -> None:
    """Per Protein.Group, count distinct Stripped.Sequence among rows that are
    proteotypic AND have at least one non-NA sample value."""
    from analysis.figure_original_vs_quantmsdiann import (
        unique_peptides_per_protein_diann,
    )
    matrix = write_matrix(
        tmp_path,
        """
        Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tProteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id\tRun_A\tRun_B
        P1\tP1\tA\tG1\td\t1\tAAAR\tAAAR\t2\tAAAR2\t10\t
        P1\tP1\tA\tG1\td\t1\tAAAR\tAAAR\t3\tAAAR3\t\t20
        P1\tP1\tA\tG1\td\t1\tBBBR\tBBBR\t2\tBBBR2\t\t30
        P1\tP1\tA\tG1\td\t0\tCCCR\tCCCR\t2\tCCCR2\t40\t
        P2\tP2\tB\tG2\td\t1\tDDDR\tDDDR\t2\tDDDR2\t50\t
        P3\tP3\tC\tG3\td\t1\tEEER\tEEER\t2\tEEER2\t\t
        """,
    )
    # P1: AAAR (proteotypic, quantified twice as 2 charge states -> 1 peptide),
    #     BBBR (proteotypic, quantified), CCCR (NOT proteotypic, skip) -> 2 unique
    # P2: DDDR (proteotypic, quantified) -> 1 unique
    # P3: EEER (proteotypic but never quantified) -> 0 -> absent from output
    out = unique_peptides_per_protein_diann(matrix)
    assert out == {"P1": 2, "P2": 1}


def test_unique_peptides_per_protein_openswath_filters_decoys_and_multimap(
    tmp_path: Path,
) -> None:
    """Per Protein (proteotypic '1/...' only), count distinct stripped peptides
    that have >=1 score <= 0.01 in any run (= confidently detected, not just
    requant-filled)."""
    from analysis.figure_original_vs_quantmsdiann import (
        unique_peptides_per_protein_openswath,
    )
    body = (
        "Peptide\tProtein\tIntensity_run0\tRT_run0\tscore_run0\tIntensity_run1\tRT_run1\tscore_run1\n"
        # P1 proteotypic, AAAR detected in run0
        "1_AAAR_2_run0\t1/sp|P1|HUMAN\t100\t10\t0.005\t150\t10.2\t2.0\n"
        # P1 proteotypic, same peptide AAAR different charge - should not double-count
        "2_AAAR_3_run0\t1/sp|P1|HUMAN\t110\t10\t0.005\t150\t10.2\t2.0\n"
        # P1 proteotypic, BBBR detected in both runs
        "3_BBBR_2_run0\t1/sp|P1|HUMAN\t200\t12\t0.001\t250\t12.1\t0.008\n"
        # P1 proteotypic, CCCR never confidently detected (both scores 2.0)
        "4_CCCR_2_run0\t1/sp|P1|HUMAN\t300\t9\t2.0\t300\t9.1\t2.0\n"
        # P2 multi-protein peptide - not unique to a single protein
        "5_DDDR_2_run0\t2/sp|P2|HUMAN/sp|P2alt|HUMAN\t400\t9\t0.001\t400\t9.1\t0.001\n"
        # DECOY peptide
        "DECOY_6_EEER_2_run0\t1/sp|P3|HUMAN\t500\t9\t0.001\t500\t9.1\t0.001\n"
        # P4 proteotypic, FFFR detected
        "7_FFFR_2_run0\t1/sp|P4|HUMAN\t600\t9\t0.001\t600\t9.1\t0.001\n"
    )
    p = tmp_path / "fa.tsv"
    p.write_text(body)
    out = unique_peptides_per_protein_openswath(p)
    # P1: {AAAR, BBBR} = 2; P4: {FFFR} = 1; P2/P3 excluded.
    assert out == {"1/sp|P1|HUMAN": 2, "1/sp|P4|HUMAN": 1}


def test_proteins_with_min_peptides_counts_above_threshold() -> None:
    from analysis.figure_original_vs_quantmsdiann import proteins_with_min_peptides
    counts = {"P1": 5, "P2": 2, "P3": 1, "P4": 10, "P5": 1}
    assert proteins_with_min_peptides(counts, 1) == 5
    assert proteins_with_min_peptides(counts, 2) == 3
    assert proteins_with_min_peptides(counts, 3) == 2
    assert proteins_with_min_peptides(counts, 10) == 1
    assert proteins_with_min_peptides(counts, 11) == 0


def test_load_hgnc_symbol_to_ensembl_includes_aliases_and_prev(tmp_path: Path) -> None:
    """HGNC complete-set TSV has: symbol, alias_symbol (pipe-separated),
    prev_symbol (pipe-separated), ensembl_gene_id. The mapping should accept
    the current symbol AND any alias / previous symbol as a key. Empty
    ensembl_gene_id rows are skipped. Last writer wins on conflicts."""
    from analysis.figure_original_vs_quantmsdiann import load_hgnc_symbol_to_ensembl
    body = (
        "hgnc_id\tsymbol\tname\tlocus_group\tlocus_type\tstatus\tlocation\tlocation_sortable\talias_symbol\talias_name\tprev_symbol\tprev_name\tgene_group\tgene_group_id\tdate_approved_reserved\tdate_symbol_changed\tdate_name_changed\tdate_modified\tentrez_id\tensembl_gene_id\n"
        "HGNC:1\tA1\tname1\tprotein\tgene\tApproved\t1q23\t1q23\tA1ALIAS\talias_name\tA1OLD\told_name\t\t\t\t\t\t\t1\tENSG00000001\n"
        "HGNC:2\tB1\tname2\tprotein\tgene\tApproved\t2q23\t2q23\tB1ALIAS1|B1ALIAS2\talias_name\t\t\t\t\t\t\t\t\t2\tENSG00000002\n"
        "HGNC:3\tC1\tname3\tprotein\tgene\tApproved\t3q23\t3q23\t\t\t\t\t\t\t\t\t\t\t3\t\n"
        "HGNC:4\tD1\tname4\tprotein\tgene\tApproved\t4q23\t4q23\t\t\tD1OLD\told_name\t\t\t\t\t\t\t4\tENSG00000004\n"
    )
    p = tmp_path / "hgnc.tsv"
    p.write_text(body)
    out = load_hgnc_symbol_to_ensembl(p)
    # A1 -> ENSG00000001 (with alias A1ALIAS, prev A1OLD)
    assert out["A1"] == "ENSG00000001"
    assert out["A1ALIAS"] == "ENSG00000001"
    assert out["A1OLD"] == "ENSG00000001"
    # B1 -> ENSG00000002, two aliases
    assert out["B1"] == "ENSG00000002"
    assert out["B1ALIAS1"] == "ENSG00000002"
    assert out["B1ALIAS2"] == "ENSG00000002"
    # C1 has no ensembl_gene_id, must NOT be in the map
    assert "C1" not in out
    # D1 with prev symbol
    assert out["D1"] == "ENSG00000004"
    assert out["D1OLD"] == "ENSG00000004"


def test_load_walzer_genes_ensembl_returns_set(tmp_path: Path) -> None:
    """E-PROT-73 has Gene ID (Ensembl) as column 1; return that as a set of
    ENSG IDs, skipping any row whose first column does not start with ENSG."""
    from analysis.figure_original_vs_quantmsdiann import load_walzer_genes_ensembl
    body = (
        "Gene ID\tGene Name\tg1.WithInSampleAbundance\n"
        "ENSG00000000001\tFOO\t1.0\n"
        "ENSG00000000002\tBAR\t\n"
        "ENSG00000000001\tFOO\t2.0\n"  # duplicate gene id -> collapsed
        "\tEmptyId\t0.5\n"
        "NOTANENSG\tX\t1.0\n"
    )
    p = tmp_path / "e.tsv"
    p.write_text(body)
    out = load_walzer_genes_ensembl(p)
    assert out == {"ENSG00000000001", "ENSG00000000002"}


def test_quantmsdiann_genes_as_ensembl_maps_and_reports_unmapped(tmp_path: Path) -> None:
    """Given the DIA-NN unique_genes_matrix and an HGNC symbol->ENSG mapping,
    return (mapped_ensg_set, unmapped_symbol_count). Only quantified rows
    (>=1 non-NA cell) count toward both numbers."""
    from analysis.figure_original_vs_quantmsdiann import quantmsdiann_genes_as_ensembl
    matrix = write_matrix(
        tmp_path,
        """
        Genes\tN.Sequences\tN.Proteotypic.Sequences\tguot_R1.mzML\tguot_R2.mzML
        A1\t2\t2\t1.0\t
        B1\t3\t3\t\t2.0
        UNKNOWN\t5\t5\t1.0\t
        EMPTYROW\t1\t1\t\t
        """,
    )
    mapping = {"A1": "ENSG00000001", "B1": "ENSG00000002"}
    mapped, unmapped_count = quantmsdiann_genes_as_ensembl(matrix, mapping)
    # A1 and B1 are quantified and map -> {ENSG1, ENSG2}; UNKNOWN is quantified
    # but doesn't map -> +1 unmapped; EMPTYROW is unquantified -> ignored.
    assert mapped == {"ENSG00000001", "ENSG00000002"}
    assert unmapped_count == 1


def test_quantmsdiann_genes_as_ensembl_applies_detection_filter(tmp_path: Path) -> None:
    """With min_detection_fraction=0.5, only rows whose non-NA-cell count is
    >= ceil(0.5 * n_sample_cols) survive. Mimics Walzer's '50% per group'
    consistency filter applied globally across all runs."""
    from analysis.figure_original_vs_quantmsdiann import quantmsdiann_genes_as_ensembl
    matrix = write_matrix(
        tmp_path,
        """
        Genes\tN.Sequences\tN.Proteotypic.Sequences\tR1\tR2\tR3\tR4
        A1\t2\t2\t1.0\t2.0\t3.0\t4.0
        B1\t3\t3\t1.0\t\t\t
        C1\t5\t5\t1.0\t2.0\t\t
        D1\t5\t5\t\t\t\t
        """,
    )
    mapping = {"A1": "ENSG1", "B1": "ENSG2", "C1": "ENSG3", "D1": "ENSG4"}
    # 4 sample cols; threshold = ceil(0.5 * 4) = 2.
    # A1 has 4 non-NA -> passes. C1 has 2 -> passes. B1 has 1 -> fails. D1: 0 -> fails.
    mapped, unmapped_count = quantmsdiann_genes_as_ensembl(
        matrix, mapping, min_detection_fraction=0.5,
    )
    assert mapped == {"ENSG1", "ENSG3"}
    assert unmapped_count == 0


def test_count_quantified_genes_diann_excludes_empty_rows(tmp_path: Path) -> None:
    """DIA-NN unique_genes_matrix.tsv layout: 3 metadata cols (Genes,
    N.Sequences, N.Proteotypic.Sequences) followed by one column per run.
    A gene is 'quantified' iff at least one sample cell is non-NA."""
    from analysis.figure_original_vs_quantmsdiann import count_quantified_genes_diann
    matrix = write_matrix(
        tmp_path,
        """
        Genes\tN.Sequences\tN.Proteotypic.Sequences\tguot_R1.mzML\tguot_R2.mzML\tguot_R3.mzML
        A1CF\t2\t2\t\t\t
        A2M\t42\t42\t\t68.3\t
        AAAS\t5\t5\t1.0\t2.0\tNA
        AACS\t3\t3\t\t\t
        """,
    )
    # A2M and AAAS have >=1 non-NA cell. A1CF and AACS are empty.
    assert count_quantified_genes_diann(matrix) == 2
