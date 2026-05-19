"""Tests for the combined cell-lines atlas figure.

These tests cover the new helpers introduced by the atlas script:
SDRF cell-line parsing + normalisation, PXD003539 disease -> ProCan
tissue harmonisation, the cell-line -> tissue loaders for all three
datasets, the tissue-table combiner, and the JSON-cache accession
extractor. All fixtures are inline; nothing hits the network.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Cell-line normalisation / SDRF parsing
# ---------------------------------------------------------------------------


def _write_sdrf(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a minimal SDRF TSV with arbitrary characteristic columns."""
    cols: list[str] = []
    seen: set[str] = set()
    for r in rows:
        for c in r:
            if c not in seen:
                seen.add(c)
                cols.append(c)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "") for c in cols) + "\n")


def test_cell_lines_from_sdrf_normalises_nci60_lines(tmp_path: Path) -> None:
    """SDRF cell-line spellings collapse via normalize_cell_line so the
    same line is the same set member across datasets (NCI-60's 'NCI-H226'
    is the same as ProCan's 'NCI-H226' and any 'H226' variant)."""
    from analysis.figure_combined_cell_lines_atlas import cell_lines_from_sdrf

    p = tmp_path / "sdrf.tsv"
    _write_sdrf(p, [
        {"characteristics[cell line]": "CCRF-CEM"},
        {"characteristics[cell line]": "NCI-H226"},
        {"characteristics[cell line]": "Hs-578-T"},
        {"characteristics[cell line]": ""},  # silently dropped
    ])
    out = cell_lines_from_sdrf(p)
    assert "CCRFCEM" in out
    assert "H226" in out  # NCI- prefix stripped by normalize_cell_line
    assert "HS578T" in out
    assert len(out) == 3


def test_cell_lines_from_sdrf_overlaps_between_pxd003539_and_pxd030304(
    tmp_path: Path,
) -> None:
    """A known NCI-60 line (CCRF-CEM, NCI-H226) must match across the
    PXD003539 (NCI-60) and PXD030304 (ProCan) SDRFs after normalisation —
    this is the foundation for Panel B's overlap."""
    from analysis.figure_combined_cell_lines_atlas import cell_lines_from_sdrf

    nci60 = tmp_path / "nci60.tsv"
    _write_sdrf(nci60, [
        {"characteristics[cell line]": "CCRF-CEM"},
        {"characteristics[cell line]": "NCI-H226"},
        {"characteristics[cell line]": "MCF7"},
    ])
    procan = tmp_path / "procan.tsv"
    _write_sdrf(procan, [
        # ProCan uses the same dashed spelling as NCI-60 in the SDRF.
        {"characteristics[cell line]": "CCRF-CEM"},
        {"characteristics[cell line]": "NCI-H226"},
        {"characteristics[cell line]": "MCF7"},
        # Lines unique to ProCan
        {"characteristics[cell line]": "HCC1419"},
        {"characteristics[cell line]": "BC-1"},
    ])
    set_a = cell_lines_from_sdrf(nci60)
    set_b = cell_lines_from_sdrf(procan)
    assert set_a <= set_b  # NCI-60 ⊂ ProCan in this fixture
    assert "CCRFCEM" in (set_a & set_b)
    assert "H226" in (set_a & set_b)


def test_cell_lines_from_sdrf_raises_on_missing_column(tmp_path: Path) -> None:
    from analysis.figure_combined_cell_lines_atlas import cell_lines_from_sdrf

    p = tmp_path / "bad.tsv"
    p.write_text("source name\tcharacteristics[organism]\nS1\tHomo sapiens\n")
    with pytest.raises(ValueError):
        cell_lines_from_sdrf(p)


# ---------------------------------------------------------------------------
# PXD003539 disease -> ProCan tissue harmonisation
# ---------------------------------------------------------------------------


def test_harmonise_pxd003539_tissue_haematopoietic(tmp_path: Path) -> None:
    """leukemia / lymphoma / myeloma all collapse to
    `Haematopoietic and Lymphoid`."""
    from analysis.figure_combined_cell_lines_atlas import harmonise_pxd003539_tissue

    assert harmonise_pxd003539_tissue("T-cell childhood acute lymphocytic leukemia",
                                      "Blood") == "Haematopoietic and Lymphoid"
    assert harmonise_pxd003539_tissue("Acute myeloid leukemia",
                                      "Bone marrow") == "Haematopoietic and Lymphoid"
    assert harmonise_pxd003539_tissue("Plasma cell myeloma",
                                      "bone marrow") == "Haematopoietic and Lymphoid"
    assert harmonise_pxd003539_tissue("Burkitt lymphoma",
                                      "Blood") == "Haematopoietic and Lymphoid"


def test_harmonise_pxd003539_tissue_lung_and_mesothelioma(
    tmp_path: Path,
) -> None:
    """Lung adenocarcinoma / non-small cell / large cell and pleural
    mesothelioma all map to ProCan `Lung`."""
    from analysis.figure_combined_cell_lines_atlas import harmonise_pxd003539_tissue

    assert harmonise_pxd003539_tissue("Lung adenocarcinoma", "Lung") == "Lung"
    assert harmonise_pxd003539_tissue("Non-small cell lung carcinoma",
                                      "Lung") == "Lung"
    assert harmonise_pxd003539_tissue("Lung large cell carcinoma",
                                      "Lung") == "Lung"
    assert harmonise_pxd003539_tissue("Pleural epithelioid mesothelioma",
                                      "pleural effusion") == "Lung"


def test_harmonise_pxd003539_tissue_cns(tmp_path: Path) -> None:
    """central nervous system / brain / glioblastoma / gliosarcoma /
    astrocytoma -> `Central Nervous System`."""
    from analysis.figure_combined_cell_lines_atlas import harmonise_pxd003539_tissue

    assert harmonise_pxd003539_tissue("central nervous system cancer",
                                      "Brain") == "Central Nervous System"
    assert harmonise_pxd003539_tissue("Glioblastoma",
                                      "Brain") == "Central Nervous System"
    assert harmonise_pxd003539_tissue("Gliosarcoma",
                                      "Brain") == "Central Nervous System"
    assert harmonise_pxd003539_tissue("Brain astrocytoma",
                                      "Brain") == "Central Nervous System"


def test_harmonise_pxd003539_tissue_colon_renal_breast_skin_ovary(
    tmp_path: Path,
) -> None:
    from analysis.figure_combined_cell_lines_atlas import harmonise_pxd003539_tissue

    assert harmonise_pxd003539_tissue("Colon adenocarcinoma",
                                      "Colon") == "Large Intestine"
    assert harmonise_pxd003539_tissue("colorectal cancer",
                                      "Large Intestine") == "Large Intestine"
    assert harmonise_pxd003539_tissue("Clear cell renal cell carcinoma",
                                      "Kidney") == "Kidney"
    assert harmonise_pxd003539_tissue("papillary renal cell carcinoma",
                                      "Kidney") == "Kidney"
    assert harmonise_pxd003539_tissue("Invasive breast carcinoma",
                                      "Breast") == "Breast"
    assert harmonise_pxd003539_tissue("Cutaneous melanoma",
                                      "Skin") == "Skin"
    assert harmonise_pxd003539_tissue("Amelanotic melanoma",
                                      "Skin") == "Skin"
    assert harmonise_pxd003539_tissue("High grade ovarian serous adenocarcinoma",
                                      "Ovary") == "Ovary"
    assert harmonise_pxd003539_tissue("Prostate carcinoma",
                                      "prostate gland") == "Prostate"


def test_harmonise_pxd003539_tissue_unmatched_returns_none() -> None:
    from analysis.figure_combined_cell_lines_atlas import harmonise_pxd003539_tissue

    assert harmonise_pxd003539_tissue("", "") is None
    assert harmonise_pxd003539_tissue("unknown weird cancer",
                                      "weird tissue") is None


# ---------------------------------------------------------------------------
# combined_tissue_table
# ---------------------------------------------------------------------------


def test_combined_tissue_table_sorts_by_total_descending() -> None:
    """combined_tissue_table aggregates per-dataset cell-line dicts into
    a per-tissue triple of counts, sorted by total descending. Tissues
    with zero contribution from all datasets cannot appear because they
    were never in the input."""
    from analysis.figure_combined_cell_lines_atlas import combined_tissue_table

    per_dataset = {
        "PXD003539": {"clA": "Lung", "clB": "Breast", "clC": "Lung"},
        "PXD030304": {"clD": "Lung", "clE": "Breast", "clF": "Breast",
                      "clG": "Bone"},
        "PXD004701": {"clH": "Breast", "clI": "Breast"},
    }
    out = combined_tissue_table(per_dataset)
    tissues = [r[0] for r in out]
    assert tissues == ["Breast", "Lung", "Bone"]  # 5, 3, 1
    by_tissue = dict(out)
    assert by_tissue["Breast"] == {
        "PXD003539": 1, "PXD030304": 2, "PXD004701": 2,
    }
    assert by_tissue["Lung"] == {
        "PXD003539": 2, "PXD030304": 1, "PXD004701": 0,
    }
    assert by_tissue["Bone"] == {
        "PXD003539": 0, "PXD030304": 1, "PXD004701": 0,
    }


# ---------------------------------------------------------------------------
# PXD030304 / PXD004701 JSON cache -> accession set
# ---------------------------------------------------------------------------


def test_pxd030304_protein_accessions_extracts_from_cached_json(
    tmp_path: Path,
) -> None:
    """The cached JSON is `{tissue: [Protein.Group, ...]}`. We union the
    Protein.Group values and pass each through extract_accessions_diann,
    which splits on ';' and strips isoform suffixes / CONTAM_/ENTRAP_
    prefixes."""
    from analysis.figure_combined_cell_lines_atlas import (
        pxd030304_protein_accessions,
    )

    payload = {
        "Lung": ["P12345", "P12345;Q67890", "P00000-2"],
        "Breast": ["sp|P11111|HUMAN", "Q00000", "P12345"],
        "Bone": [],
    }
    p = tmp_path / "cache.json"
    p.write_text(json.dumps(payload))
    out = pxd030304_protein_accessions(p)
    # P12345 dedup across tissues; P00000-2 -> P00000; sp|P11111|HUMAN -> P11111
    assert out == {"P12345", "Q67890", "P00000", "P11111", "Q00000"}


def test_pxd004701_protein_accessions_extracts_from_cached_json(
    tmp_path: Path,
) -> None:
    from analysis.figure_combined_cell_lines_atlas import (
        pxd004701_protein_accessions,
    )

    payload = {
        "TNBC": ["P02768;Q9NQ29", "CONTAM_P02768;P02768"],
        "non-TNBC": ["P02768"],
        "normal-like": ["O15156"],
    }
    p = tmp_path / "cache.json"
    p.write_text(json.dumps(payload))
    out = pxd004701_protein_accessions(p)
    # CONTAM_ prefix stripped, P02768 + Q9NQ29 + O15156.
    assert out == {"P02768", "Q9NQ29", "O15156"}


# ---------------------------------------------------------------------------
# Venn region helper used by the counts.tsv writer
# ---------------------------------------------------------------------------


def test_venn_region_sizes_3_partitions_correctly() -> None:
    """The 7 region sizes must sum to the size of the 3-way union."""
    from analysis.figure_combined_cell_lines_atlas import _venn_region_sizes_3

    sets = {
        "A": {"x1", "x2", "x3", "x4"},      # only-A: x1; A+B: x2; A+C: x3; ABC: x4
        "B": {"x2", "x4", "x5", "x6"},      # only-B: x5; B+C: x6
        "C": {"x3", "x4", "x6", "x7"},      # only-C: x7
    }
    out = _venn_region_sizes_3(sets, ["A", "B", "C"])
    assert out == {
        "A_only": 1, "B_only": 1, "C_only": 1,
        "A+B": 1, "A+C": 1, "B+C": 1,
        "all_three": 1,
    }
    total = sum(out.values())
    assert total == len(sets["A"] | sets["B"] | sets["C"]) == 7
