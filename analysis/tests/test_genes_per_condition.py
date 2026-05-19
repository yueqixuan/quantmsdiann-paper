"""Tests for per-condition (disease/tissue) gene-detection parsing.

The two pipelines we compare label samples differently:

- quantmsdiann (PRIDE quantms-collections SDRF) uses dashed cell-line names
  ('CCRF-CEM', 'NCI-H226') and granular `characteristics[disease]` labels
  ('T-cell childhood acute lymphocytic leukemia').
- Walzer 2022 / Expression Atlas E-PROT-73 uses no-dash cell-line names
  ('CCRFCEM', 'H226') and the 9 NCI-60 cancer types as `disease`
  ('leukemia', 'lung cancer', ...).

To compare gene detections by condition, we normalise cell-line names and use
the E-PROT-73 cell-line -> disease map as the canonical condition axis.
"""
from __future__ import annotations

from pathlib import Path
import textwrap

import pytest


def test_normalize_cell_line_strips_dashes_and_uppercases() -> None:
    from analysis.figure_original_vs_quantmsdiann import normalize_cell_line
    assert normalize_cell_line("CCRF-CEM") == "CCRFCEM"
    assert normalize_cell_line("HCT-116") == "HCT116"
    assert normalize_cell_line("Hs-578-T") == "HS578T"
    assert normalize_cell_line("Malme-3M") == "MALME3M"


def test_normalize_cell_line_strips_nci_prefix() -> None:
    from analysis.figure_original_vs_quantmsdiann import normalize_cell_line
    # quantms SDRF uses 'NCI-H226'; E-PROT-73 uses 'H226'.
    assert normalize_cell_line("NCI-H226") == "H226"
    assert normalize_cell_line("NCI-H522") == "H522"
    # Already-normalised input is idempotent.
    assert normalize_cell_line("H226") == "H226"


def test_normalize_cell_line_handles_none_and_empty() -> None:
    from analysis.figure_original_vs_quantmsdiann import normalize_cell_line
    assert normalize_cell_line("") == ""
    assert normalize_cell_line(None) == ""  # type: ignore[arg-type]


def test_parse_eprot73_groupings_extracts_cell_line_and_disease(
    tmp_path: Path,
) -> None:
    """The Downloads HTML embeds a JSON `content:{...}` blob with primary
    groupings under CELL_LINE and DISEASE. We extract both as
    g<N> -> name dicts."""
    from analysis.figure_original_vs_quantmsdiann import (
        parse_eprot73_groupings,
    )
    # Minimal HTML that mimics the structure of the real Downloads page.
    fake_html = (
        '<html><body><script>x={'
        '"experimentAccession":"E-PROT-73","tabs":[{"props":{"groups":['
        '{"name":"CELL_LINE","primary":true,"selected":"all","groupings":'
        '[["A498",["g1"]],["CCRFCEM",["g6"]],["HCT116",["g16"]]]},'
        '{"name":"DISEASE","primary":true,"selected":"lung cancer","groupings":'
        '[["renal cancer",["g1"]],["leukemia",["g6"]],'
        '["colorectal cancer",["g16"]]]}]}}]}</script></body></html>'
    )
    p = tmp_path / "downloads.html"
    p.write_text(fake_html)
    g_to_cl, g_to_ds = parse_eprot73_groupings(p)
    assert g_to_cl == {"g1": "A498", "g6": "CCRFCEM", "g16": "HCT116"}
    assert g_to_ds == {"g1": "renal cancer", "g6": "leukemia",
                       "g16": "colorectal cancer"}


def test_load_sdrf_data_file_to_cell_line_normalises_to_mzml(
    tmp_path: Path,
) -> None:
    """SDRF `comment[data file]` is `*.wiff` but DIA-NN matrix columns are
    `*.mzML`. The loader rewrites `.wiff` -> `.mzML` so downstream code can
    look up DIA-NN columns by SDRF file name."""
    from analysis.figure_original_vs_quantmsdiann import (
        load_sdrf_data_file_to_cell_line,
    )
    body = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tCCRF-CEM\tguot_L130610_003_SW.wiff\n"
        "S2\tHCT-116\tguot_L130610_005_SW.wiff\n"
        "S3\tNCI-H226\tguot_L130618_004_SW.wiff\n"
    )
    p = tmp_path / "sdrf.tsv"
    p.write_text(body)
    out = load_sdrf_data_file_to_cell_line(p)
    assert out == {
        "guot_L130610_003_SW.mzML": "CCRF-CEM",
        "guot_L130610_005_SW.mzML": "HCT-116",
        "guot_L130618_004_SW.mzML": "NCI-H226",
    }


def test_load_sdrf_data_file_to_cell_line_raises_on_missing_columns(
    tmp_path: Path,
) -> None:
    from analysis.figure_original_vs_quantmsdiann import (
        load_sdrf_data_file_to_cell_line,
    )
    body = "source name\tcomment[data file]\nS1\tx.wiff\n"
    p = tmp_path / "sdrf.tsv"
    p.write_text(body)
    with pytest.raises(ValueError):
        load_sdrf_data_file_to_cell_line(p)


def test_load_ea_cell_line_to_disease(tmp_path: Path) -> None:
    """E-PROT-73 experiment-design TSV gives each MS Run a cell line and a
    disease. We collapse to (normalised cell line) -> disease."""
    from analysis.figure_original_vs_quantmsdiann import (
        load_ea_cell_line_to_disease,
    )
    body = (
        "Run\tSample Characteristic[organism]\tSample Characteristic Ontology Term[organism]"
        "\tSample Characteristic[cell line]\tSample Characteristic Ontology Term[cell line]"
        "\tSample Characteristic[organism part]\tSample Characteristic Ontology Term[organism part]"
        "\tSample Characteristic[disease]\tSample Characteristic Ontology Term[disease]"
        "\tFactor Value[disease]\tFactor Value Ontology Term[disease]"
        "\tFactor Value[cell line]\tFactor Value Ontology Term[cell line]\tAnalysed\n"
        "guot_X1_SW\tHomo sapiens\tNCBITaxon_9606\tCCRFCEM\tEFO_0002128\tblood\tUBERON_0000178"
        "\tleukemia\tEFO_0000565\tleukemia\tEFO_0000565\tCCRFCEM\tEFO_0002128\tYes\n"
        "guot_X2_SW\tHomo sapiens\tNCBITaxon_9606\tCCRFCEM\tEFO_0002128\tblood\tUBERON_0000178"
        "\tleukemia\tEFO_0000565\tleukemia\tEFO_0000565\tCCRFCEM\tEFO_0002128\tYes\n"
        "guot_X3_SW\tHomo sapiens\tNCBITaxon_9606\tHCT116\tEFO_0002824\tcolon\tUBERON_0001155"
        "\tcolorectal cancer\tEFO_0005842\tcolorectal cancer\tEFO_0005842\tHCT116\tEFO_0002824\tYes\n"
    )
    p = tmp_path / "design.tsv"
    p.write_text(body)
    out = load_ea_cell_line_to_disease(p)
    # Both runs of CCRFCEM collapse; key is the normalised cell-line name.
    assert out == {"CCRFCEM": "leukemia", "HCT116": "colorectal cancer"}


def test_walzer_genes_per_condition_from_eprot73(tmp_path: Path) -> None:
    """For each disease, the detected gene set is the union over g<N> columns
    belonging to that disease of genes with abundance > 0. Genes that are 0
    in all g<N> for a disease are not detected for that disease."""
    from analysis.figure_original_vs_quantmsdiann import (
        walzer_genes_per_condition,
    )
    fake_html = (
        '<html><script>x={"tabs":[{"props":{"groups":['
        '{"name":"CELL_LINE","primary":true,"selected":"all","groupings":'
        '[["CCRFCEM",["g6"]],["HL60",["g18"]],["HCT116",["g16"]]]},'
        '{"name":"DISEASE","primary":true,"selected":"x","groupings":'
        '[["leukemia",["g6","g18"]],["colorectal cancer",["g16"]]]}]}}]}'
        '</script></html>'
    )
    html_p = tmp_path / "downloads.html"
    html_p.write_text(fake_html)
    eprot = (
        "Gene ID\tGene Name\tg6.WithInSampleAbundance\tg16.WithInSampleAbundance"
        "\tg18.WithInSampleAbundance\n"
        # ENSG1: detected in CCRFCEM (g6) -> leukemia.
        "ENSG00000000001\tGENE1\t10\t0\t0\n"
        # ENSG2: detected in HCT116 (g16) -> colorectal cancer.
        "ENSG00000000002\tGENE2\t0\t5\t0\n"
        # ENSG3: detected in HL60 (g18) -> also leukemia.
        "ENSG00000000003\tGENE3\t0\t0\t2\n"
        # ENSG4: detected everywhere -> both diseases.
        "ENSG00000000004\tGENE4\t1\t1\t1\n"
        # ENSG5: zero everywhere -> neither.
        "ENSG00000000005\tGENE5\t0\t0\t0\n"
        # Header artifact row that doesn't start with ENSG (sometimes appears
        # in real Expression Atlas dumps) must be skipped.
        "totalGenes\t\t100\t100\t100\n"
    )
    eprot_p = tmp_path / "eprot.tsv"
    eprot_p.write_text(eprot)
    out = walzer_genes_per_condition(eprot_p, html_p)
    assert out == {
        "leukemia": {"ENSG00000000001", "ENSG00000000003", "ENSG00000000004"},
        "colorectal cancer": {"ENSG00000000002", "ENSG00000000004"},
    }


def test_quantmsdiann_genes_per_condition_global_50pct_filter(
    tmp_path: Path,
) -> None:
    """min_global_detection_fraction=0.5 requires non-NA in >=ceil(0.5 * n_total_runs)
    of all DIA-NN runs (across cell lines). Same filter as the global Venn."""
    from analysis.figure_original_vs_quantmsdiann import (
        quantmsdiann_genes_per_condition,
    )
    # 6 runs total -> ceil(0.5*6)=3 required globally.
    matrix = (
        "Genes\tN.Sequences\tN.Proteotypic.Sequences"
        "\tguot_R1_SW.mzML\tguot_R2_SW.mzML"
        "\tguot_R3_SW.mzML\tguot_R4_SW.mzML\tguot_R5_SW.mzML\tguot_R6_SW.mzML\n"
        # GENE1: 3 non-NA globally (pass); in CCRFCEM (R1,R2) and HCT116 (R3)
        # so both diseases.
        "GENE1\t1\t1\t1\t1\t1\t\t\t\n"
        # GENE2: 2 non-NA globally (fail global filter); dropped entirely.
        "GENE2\t1\t1\t1\t1\t\t\t\t\n"
        # GENE3: 6 non-NA globally (pass); both conditions.
        "GENE3\t1\t1\t1\t1\t1\t1\t1\t1\n"
        # GENE4: 4 non-NA globally (pass), but only HCT116 runs (R3..R6)
        # -> only colorectal.
        "GENE4\t1\t1\t\t\t1\t1\t1\t1\n"
    )
    matrix_p = tmp_path / "ugmatrix.tsv"
    matrix_p.write_text(matrix)
    sdrf = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tCCRF-CEM\tguot_R1_SW.wiff\n"
        "S2\tCCRF-CEM\tguot_R2_SW.wiff\n"
        "S3\tHCT-116\tguot_R3_SW.wiff\n"
        "S4\tHCT-116\tguot_R4_SW.wiff\n"
        "S5\tHCT-116\tguot_R5_SW.wiff\n"
        "S6\tHCT-116\tguot_R6_SW.wiff\n"
    )
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text(sdrf)
    design = (
        "Run\tSample Characteristic[organism]\tSample Characteristic Ontology Term[organism]"
        "\tSample Characteristic[cell line]\tSample Characteristic Ontology Term[cell line]"
        "\tSample Characteristic[organism part]\tSample Characteristic Ontology Term[organism part]"
        "\tSample Characteristic[disease]\tSample Characteristic Ontology Term[disease]"
        "\tFactor Value[disease]\tFactor Value Ontology Term[disease]"
        "\tFactor Value[cell line]\tFactor Value Ontology Term[cell line]\tAnalysed\n"
        "guot_R1_SW\t\t\tCCRFCEM\t\tblood\t\tleukemia\t\tleukemia\t\tCCRFCEM\t\tYes\n"
        "guot_R3_SW\t\t\tHCT116\t\tcolon\t\tcolorectal cancer\t\tcolorectal cancer\t\tHCT116\t\tYes\n"
    )
    design_p = tmp_path / "design.tsv"
    design_p.write_text(design)
    symbol_to_ensg = {
        "GENE1": "ENSG00000000001",
        "GENE2": "ENSG00000000002",
        "GENE3": "ENSG00000000003",
        "GENE4": "ENSG00000000004",
    }
    out = quantmsdiann_genes_per_condition(
        matrix_p, sdrf_p, design_p, symbol_to_ensg,
        min_global_detection_fraction=0.5,
    )
    # GENE2 dropped (2/6 < 3 global). GENE1 in both diseases (R1,R2,R3 span
    # CCRFCEM+HCT116). GENE3 in both. GENE4 only HCT116.
    assert out == {
        "leukemia": {"ENSG00000000001", "ENSG00000000003"},
        "colorectal cancer": {
            "ENSG00000000001", "ENSG00000000003", "ENSG00000000004",
        },
    }


def test_quantmsdiann_genes_per_condition_50pct_per_cell_line_filter(
    tmp_path: Path,
) -> None:
    """With min_detection_fraction_per_cell_line=0.5, a gene must be non-NA
    in at least ceil(0.5 * n_runs_of_cell_line) runs of a given cell line to
    count as detected for that cell line — matching Walzer's 50%-per-group
    consistency filter."""
    from analysis.figure_original_vs_quantmsdiann import (
        quantmsdiann_genes_per_condition,
    )
    # CCRF-CEM has 2 runs (R1, R2 -> ceil(0.5*2)=1 required: same as no filter).
    # HCT-116 has 4 runs (R3..R6 -> ceil(0.5*4)=2 required).
    matrix = (
        "Genes\tN.Sequences\tN.Proteotypic.Sequences"
        "\tguot_R1_SW.mzML\tguot_R2_SW.mzML"
        "\tguot_R3_SW.mzML\tguot_R4_SW.mzML\tguot_R5_SW.mzML\tguot_R6_SW.mzML\n"
        # GENE1: 1/2 CCRFCEM (pass), 1/4 HCT116 (fail) -> only leukemia.
        "GENE1\t1\t1\t10\t\t5\t\t\t\n"
        # GENE2: 0/2 CCRFCEM (fail), 2/4 HCT116 (pass) -> only colorectal.
        "GENE2\t1\t1\t\t\t1\t2\t\t\n"
        # GENE3: 2/2 CCRFCEM (pass), 4/4 HCT116 (pass) -> both.
        "GENE3\t1\t1\t1\t1\t1\t1\t1\t1\n"
        # GENE4: 1/2 CCRFCEM (pass), 1/4 HCT116 (fail) -> only leukemia.
        "GENE4\t1\t1\t\t9\t\t3\t\t\n"
    )
    matrix_p = tmp_path / "ugmatrix.tsv"
    matrix_p.write_text(matrix)
    sdrf = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tCCRF-CEM\tguot_R1_SW.wiff\n"
        "S2\tCCRF-CEM\tguot_R2_SW.wiff\n"
        "S3\tHCT-116\tguot_R3_SW.wiff\n"
        "S4\tHCT-116\tguot_R4_SW.wiff\n"
        "S5\tHCT-116\tguot_R5_SW.wiff\n"
        "S6\tHCT-116\tguot_R6_SW.wiff\n"
    )
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text(sdrf)
    design = (
        "Run\tSample Characteristic[organism]\tSample Characteristic Ontology Term[organism]"
        "\tSample Characteristic[cell line]\tSample Characteristic Ontology Term[cell line]"
        "\tSample Characteristic[organism part]\tSample Characteristic Ontology Term[organism part]"
        "\tSample Characteristic[disease]\tSample Characteristic Ontology Term[disease]"
        "\tFactor Value[disease]\tFactor Value Ontology Term[disease]"
        "\tFactor Value[cell line]\tFactor Value Ontology Term[cell line]\tAnalysed\n"
        "guot_R1_SW\t\t\tCCRFCEM\t\tblood\t\tleukemia\t\tleukemia\t\tCCRFCEM\t\tYes\n"
        "guot_R3_SW\t\t\tHCT116\t\tcolon\t\tcolorectal cancer\t\tcolorectal cancer\t\tHCT116\t\tYes\n"
    )
    design_p = tmp_path / "design.tsv"
    design_p.write_text(design)
    symbol_to_ensg = {
        "GENE1": "ENSG00000000001",
        "GENE2": "ENSG00000000002",
        "GENE3": "ENSG00000000003",
        "GENE4": "ENSG00000000004",
    }
    out = quantmsdiann_genes_per_condition(
        matrix_p, sdrf_p, design_p, symbol_to_ensg,
        min_detection_fraction_per_cell_line=0.5,
    )
    assert out == {
        "leukemia": {"ENSG00000000001", "ENSG00000000003", "ENSG00000000004"},
        "colorectal cancer": {"ENSG00000000002", "ENSG00000000003"},
    }


def test_quantmsdiann_genes_per_condition_uses_sdrf_and_ea_design(
    tmp_path: Path,
) -> None:
    """quantmsdiann condition mapping: each DIA-NN run column -> SDRF cell line
    -> (via E-PROT-73 design) -> disease. Detected genes per disease are the
    union of rows with non-NA in any column belonging to that disease, with
    Genes mapped to Ensembl through the HGNC symbol map."""
    from analysis.figure_original_vs_quantmsdiann import (
        quantmsdiann_genes_per_condition,
    )
    matrix = (
        "Genes\tN.Sequences\tN.Proteotypic.Sequences"
        "\tguot_R1_SW.mzML\tguot_R2_SW.mzML\tguot_R3_SW.mzML\n"
        # GENE1 detected only in R1 (CCRFCEM -> leukemia)
        "GENE1\t1\t1\t100\t\t\n"
        # GENE2 detected only in R3 (HCT116 -> colorectal cancer)
        "GENE2\t1\t1\t\t\t50\n"
        # GENE3 detected in R1 and R2; both are CCRFCEM-related so leukemia
        "GENE3\t2\t2\t20\t10\t\n"
        # GENE4 has no HGNC mapping -> dropped from the per-condition sets
        "GENE4\t1\t1\t5\t\t\n"
        # GENE5 detected nowhere -> not in any condition
        "GENE5\t1\t1\t\t\t\n"
    )
    matrix_p = tmp_path / "ugmatrix.tsv"
    matrix_p.write_text(matrix)
    sdrf = (
        "source name\tcharacteristics[cell line]\tcomment[data file]\n"
        "S1\tCCRF-CEM\tguot_R1_SW.wiff\n"
        "S2\tCCRF-CEM\tguot_R2_SW.wiff\n"
        "S3\tHCT-116\tguot_R3_SW.wiff\n"
    )
    sdrf_p = tmp_path / "sdrf.tsv"
    sdrf_p.write_text(sdrf)
    design = (
        "Run\tSample Characteristic[organism]\tSample Characteristic Ontology Term[organism]"
        "\tSample Characteristic[cell line]\tSample Characteristic Ontology Term[cell line]"
        "\tSample Characteristic[organism part]\tSample Characteristic Ontology Term[organism part]"
        "\tSample Characteristic[disease]\tSample Characteristic Ontology Term[disease]"
        "\tFactor Value[disease]\tFactor Value Ontology Term[disease]"
        "\tFactor Value[cell line]\tFactor Value Ontology Term[cell line]\tAnalysed\n"
        "guot_R1_SW\t\t\tCCRFCEM\t\tblood\t\tleukemia\t\tleukemia\t\tCCRFCEM\t\tYes\n"
        "guot_R3_SW\t\t\tHCT116\t\tcolon\t\tcolorectal cancer\t\tcolorectal cancer\t\tHCT116\t\tYes\n"
    )
    design_p = tmp_path / "design.tsv"
    design_p.write_text(design)
    symbol_to_ensg = {
        "GENE1": "ENSG00000000001",
        "GENE2": "ENSG00000000002",
        "GENE3": "ENSG00000000003",
        # GENE4 deliberately missing.
    }
    out = quantmsdiann_genes_per_condition(
        matrix_p, sdrf_p, design_p, symbol_to_ensg,
    )
    assert out == {
        "leukemia": {"ENSG00000000001", "ENSG00000000003"},
        "colorectal cancer": {"ENSG00000000002"},
    }
