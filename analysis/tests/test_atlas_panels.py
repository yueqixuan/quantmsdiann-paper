"""Tests for the new (E/F/G/H) panels of the combined cell-lines atlas.

These exercise the pure data-transform helpers (rarefaction, detection
histogram, per-tissue protein-count rollup) and the Expression Atlas
loader, plus a smoke test that `render_atlas` writes a non-empty SVG
when fed synthetic inputs for all 8 panels.
"""
from __future__ import annotations

from pathlib import Path


# ---------------------------------------------------------------------------
# Panel G — detection-count histogram
# ---------------------------------------------------------------------------

def test_detection_count_histogram_buckets_three_sets() -> None:
    """Each accession's bucket is the number of sets it appears in. The
    histogram has one entry per possible count (1..N)."""
    from analysis.figure_combined_cell_lines_atlas import detection_count_histogram

    sets = {
        "A": {"x1", "x2", "x3", "x4"},  # x1 only-A
        "B": {"x2", "x4", "x5"},        # x5 only-B
        "C": {"x3", "x4", "x6"},        # x6 only-C
    }
    # x1 (A), x5 (B), x6 (C)            -> count=1, bucket size 3
    # x2 (A,B), x3 (A,C)                -> count=2, bucket size 2
    # x4 (A,B,C)                        -> count=3, bucket size 1
    out = detection_count_histogram(sets)
    assert out == {1: 3, 2: 2, 3: 1}


def test_detection_count_histogram_zero_for_missing_bucket() -> None:
    """Buckets with no accessions still appear with size 0 (the panel
    needs N bars even when one of them is empty)."""
    from analysis.figure_combined_cell_lines_atlas import detection_count_histogram

    sets = {
        "A": {"x1"},
        "B": {"x1"},  # all in both -> count=2 only
        "C": {"x1"},  # actually all in all three -> count=3 only
    }
    out = detection_count_histogram(sets)
    assert out == {1: 0, 2: 0, 3: 1}


# ---------------------------------------------------------------------------
# Panel E — rarefaction curve
# ---------------------------------------------------------------------------

def test_rarefaction_curve_monotonic_and_caps_at_union_size() -> None:
    """The curve must be non-decreasing as more groups are accumulated,
    and its final value must equal the size of the union across all
    groups (no matter the permutation)."""
    from analysis.figure_combined_cell_lines_atlas import rarefaction_curve

    groups = {
        "g1": {"a", "b"},
        "g2": {"b", "c"},
        "g3": {"c", "d", "e"},
        "g4": {"a", "e", "f"},
    }
    curve = rarefaction_curve(groups, n_permutations=20, seed=123)
    # Non-decreasing
    for prev, nxt in zip(curve, curve[1:]):
        assert nxt + 1e-9 >= prev, (prev, nxt, curve)
    # Final value equals union size (deterministic regardless of order)
    union = set()
    for s in groups.values():
        union.update(s)
    assert curve[-1] == float(len(union))
    # Curve has one point per group
    assert len(curve) == len(groups)


def test_rarefaction_curve_handles_empty_and_singleton() -> None:
    """0-group input returns an empty curve; 1-group input returns the
    group's own size without averaging."""
    from analysis.figure_combined_cell_lines_atlas import rarefaction_curve

    assert rarefaction_curve({}) == []
    assert rarefaction_curve({"only": {"a", "b", "c"}}) == [3.0]


# ---------------------------------------------------------------------------
# Panel F — per-tissue protein counts
# ---------------------------------------------------------------------------

def test_per_tissue_proteins_zero_for_missing_dataset_contribution() -> None:
    """When a dataset has no cell lines mapped to a tissue (or, for
    PXD004701, the tissue is not Breast), the rollup returns 0 — never
    NaN, never None."""
    from analysis.figure_combined_cell_lines_atlas import per_tissue_protein_counts

    tissue_order = ["Breast", "Lung", "Bone"]
    # PXD003539 has one cell line in Breast, none elsewhere
    cl_tissue_pxd003539 = {"clA": "Breast"}
    per_cl_pxd003539 = {"clA": {"P00001", "P00002"}}
    # PXD030304 has lung-only proteins; Bone empty in its cache
    per_tissue_pxd030304 = {
        "Lung": {"P00003", "P00004", "P00005"},
        "Breast": {"P00001"},
    }
    # PXD004701 is breast-only — all proteins roll into Breast
    pxd004701_union = {"P00006", "P00007"}

    rows = per_tissue_protein_counts(
        tissue_order, cl_tissue_pxd003539, per_cl_pxd003539,
        per_tissue_pxd030304, pxd004701_union,
    )
    by_tissue = dict(rows)
    assert by_tissue["Breast"] == {
        "PXD003539": 2, "PXD030304": 1, "PXD004701": 2,
    }
    # Lung: PXD003539 contributes nothing (no cell line mapped); PXD004701 == 0
    assert by_tissue["Lung"] == {
        "PXD003539": 0, "PXD030304": 3, "PXD004701": 0,
    }
    # Bone: every dataset is 0
    assert by_tissue["Bone"] == {
        "PXD003539": 0, "PXD030304": 0, "PXD004701": 0,
    }
    # No NaN slipped in
    for _, by_ds in rows:
        for v in by_ds.values():
            assert isinstance(v, int)


# ---------------------------------------------------------------------------
# Panel H — Expression Atlas gene set loader
# ---------------------------------------------------------------------------

def test_expression_atlas_gene_set_skips_comment_lines(tmp_path: Path) -> None:
    """Loader must skip `#`-prefixed comment lines (the actual EA file
    starts with 4 such lines) and return the unique set of `Gene Name`
    entries — ignoring blanks and de-duplicating."""
    from analysis.figure_combined_cell_lines_atlas import expression_atlas_gene_set

    p = tmp_path / "ea.tsv"
    p.write_text(
        "# Expression Atlas\n"
        "# Query: synthetic\n"
        "# Timestamp: now\n"
        "Gene ID\tGene Name\tA549\tMCF7\n"
        "ENSG0001\tTP53\t0.1\t0.2\n"
        "ENSG0002\tBRCA1\t\t0.5\n"
        "ENSG0003\t\t0.4\t\n"          # blank gene name dropped
        "ENSG0004\tTP53\t0.7\t0.3\n"   # duplicate gene name collapsed
    )
    out = expression_atlas_gene_set(p)
    assert out == {"TP53", "BRCA1"}


def test_expression_atlas_gene_set_returns_empty_on_missing_file(
    tmp_path: Path,
) -> None:
    """If the EA file is absent, the loader returns an empty set
    (Panel H renders an explainer panel)."""
    from analysis.figure_combined_cell_lines_atlas import expression_atlas_gene_set

    assert expression_atlas_gene_set(tmp_path / "does_not_exist.tsv") == set()


# ---------------------------------------------------------------------------
# render_atlas — smoke test on synthetic inputs (all 8 panels)
# ---------------------------------------------------------------------------

def test_render_atlas_runs_with_extended_inputs(tmp_path: Path) -> None:
    """End-to-end: feed minimal synthetic inputs for the remaining
    panels (Panel H removed; the breadth-vs-depth scatter dropped from
    atlas_distribution on 2026-05-29) and confirm `render_atlas` writes
    two non-empty SVGs."""
    from analysis.figure_combined_cell_lines_atlas import (
        DATASET_HEADLINES, render_atlas,
    )

    cellline_sets = {
        "PXD003539": {"clA", "clB"},
        "PXD030304": {"clB", "clC"},
        "PXD004701": {"clA", "clD"},
    }
    tissue_rows = [
        ("Breast", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 1}),
        ("Lung", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 0}),
    ]
    accession_sets = {
        "PXD003539": {"P0001", "P0002"},
        "PXD030304": {"P0002", "P0003"},
        "PXD004701": {"P0001", "P0003"},
    }
    tissue_protein_rows = [
        ("Breast", {"PXD003539": 2, "PXD030304": 1, "PXD004701": 2}),
        ("Lung", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 0}),
    ]
    runs_per_cohort = {"PXD003539": 100, "PXD030304": 5000, "PXD004701": 300}

    svg = tmp_path / "atlas.svg"
    render_atlas(
        DATASET_HEADLINES,
        cellline_sets,
        tissue_rows,
        accession_sets,
        svg,
        tissue_protein_rows=tissue_protein_rows,
        runs_per_cohort=runs_per_cohort,
    )
    overlap = tmp_path / "atlas_overlap.svg"
    distribution = tmp_path / "atlas_distribution.svg"
    assert overlap.exists() and overlap.stat().st_size > 0
    assert distribution.exists() and distribution.stat().st_size > 0
    for p in (overlap, distribution):
        head = p.read_text(encoding="utf-8")[:256]
        assert "<svg" in head or "<?xml" in head


# ---------------------------------------------------------------------------
# Panel A — paper-bar skip for PXD017199-style entries
# ---------------------------------------------------------------------------

def test_panel_a_skips_paper_bar_for_pxd017199() -> None:
    """Panel A's renderer must skip the grey paper bar when a dataset's
    `paper_count == 0` (Tognetti 2021 / PXD017199 case). With one
    dataset carrying paper_count=0, the resulting Axes should have one
    fewer Rectangle child than `2 * n_datasets` would suggest."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    from analysis.figure_combined_cell_lines_atlas import (
        DatasetHeadline, _render_panel_a_headlines,
    )

    headlines = {
        "DS1": DatasetHeadline(
            paper_count=1000, diann_count=1100,
            paper_label="DS1 paper", metric="proteins",
        ),
        "DS2": DatasetHeadline(
            paper_count=2000, diann_count=2100,
            paper_label="DS2 paper", metric="proteins",
        ),
        # No paper headline — Panel A must skip the grey bar for this one
        "DS017199": DatasetHeadline(
            paper_count=0, diann_count=10572,
            paper_label="", metric="protein groups",
        ),
    }
    fig, ax = plt.subplots(figsize=(6, 4))
    _render_panel_a_headlines(ax, headlines)
    # Count bars: matplotlib bar() adds Rectangle patches. Two
    # additional Rectangles are always added by matplotlib for the
    # axes background; filter to bars whose height matches our values.
    bar_heights = sorted({
        int(round(r.get_height()))
        for r in ax.patches
        if isinstance(r, Rectangle) and r.get_height() > 0.5
    })
    # 2 paper bars + 3 diann bars = 5 distinct heights
    assert 1000 in bar_heights
    assert 2000 in bar_heights
    assert 1100 in bar_heights
    assert 2100 in bar_heights
    assert 10572 in bar_heights
    # The no-paper-bar dataset must NOT contribute a 0-height paper bar.
    n_total_bars = sum(
        1 for r in ax.patches
        if isinstance(r, Rectangle) and r.get_height() > 0.5
    )
    assert n_total_bars == 5  # 2 paper + 3 diann (NOT 6 = 2*3)
    plt.close(fig)


# ---------------------------------------------------------------------------
# render_atlas — 4-dataset synthetic inputs
# ---------------------------------------------------------------------------

def test_render_atlas_runs_with_four_datasets(tmp_path: Path) -> None:
    """End-to-end with 4 datasets including a PXD017199-style entry
    (paper_count=0). Must produce a non-empty SVG without raising."""
    from analysis.figure_combined_cell_lines_atlas import (
        DatasetHeadline, render_atlas,
    )

    headlines = {
        "PXD003539": DatasetHeadline(6556, 6927, "Guo 2019", "Protein groups"),
        "PXD030304": DatasetHeadline(8498, 9370, "ProCan 2022", "Proteins"),
        "PXD004701": DatasetHeadline(6091, 6296, "Sun 2023", "Proteins"),
        "PXD017199": DatasetHeadline(
            0, 10572, "", "Protein groups (1% global q-value)"
        ),
    }
    cellline_sets = {
        "PXD003539": {"clA", "clB"},
        "PXD030304": {"clB", "clC"},
        "PXD004701": {"clA", "clD"},
        "PXD017199": {"clA", "clB", "clE", "clF"},
    }
    tissue_rows = [
        ("Breast", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 1,
                    "PXD017199": 3}),
        ("Lung", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 0,
                  "PXD017199": 0}),
        ("Healthy (Non-cancer)", {
            "PXD003539": 0, "PXD030304": 0, "PXD004701": 0,
            "PXD017199": 1,
        }),
    ]
    accession_sets = {
        "PXD003539": {"P0001", "P0002"},
        "PXD030304": {"P0002", "P0003"},
        "PXD004701": {"P0001", "P0003"},
        "PXD017199": {"P0001", "P0002", "P0003", "P0004"},
    }
    tissue_protein_rows = [
        ("Breast", {"PXD003539": 2, "PXD030304": 1, "PXD004701": 2,
                    "PXD017199": 4}),
        ("Lung", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 0,
                  "PXD017199": 0}),
    ]
    runs_per_cohort = {
        "PXD003539": 120, "PXD030304": 5798, "PXD004701": 300,
        "PXD017199": 206,
    }

    svg = tmp_path / "atlas.svg"
    render_atlas(
        headlines,
        cellline_sets,
        tissue_rows,
        accession_sets,
        svg,
        tissue_protein_rows=tissue_protein_rows,
        runs_per_cohort=runs_per_cohort,
    )
    overlap = tmp_path / "atlas_overlap.svg"
    distribution = tmp_path / "atlas_distribution.svg"
    assert overlap.exists() and overlap.stat().st_size > 0
    assert distribution.exists() and distribution.stat().st_size > 0


def test_render_atlas_runs_with_five_datasets(tmp_path: Path) -> None:
    """End-to-end with all 5 cohorts (PXD003539 / PXD030304 / PXD004701 /
    PXD017199 / PXD041421). Must produce a non-empty SVG without raising
    — exercises the n=5 path through `_set_region_sizes` and the
    extended ds_order in every renderer."""
    from analysis.figure_combined_cell_lines_atlas import (
        DatasetHeadline, render_atlas,
    )

    headlines = {
        "PXD003539": DatasetHeadline(6556, 6840, "Guo 2019", "Protein groups"),
        "PXD030304": DatasetHeadline(8498, 9050, "ProCan 2022", "Proteins"),
        "PXD004701": DatasetHeadline(6091, 6296, "Sun 2023", "Proteins"),
        "PXD017199": DatasetHeadline(
            0, 10380, "", "Protein groups (1% global q-value)",
        ),
        "PXD041421": DatasetHeadline(
            0, 8950, "", "Protein groups (1% global q-value)",
        ),
    }
    cellline_sets = {
        "PXD003539": {"clA", "clB"},
        "PXD030304": {"clB", "clC"},
        "PXD004701": {"clA", "clD"},
        "PXD017199": {"clA", "clB", "clE", "clF"},
        "PXD041421": {"clA", "clC"},
    }
    tissue_rows = [
        ("Breast", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 1,
                    "PXD017199": 3, "PXD041421": 0}),
        ("Lung", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 0,
                  "PXD017199": 0, "PXD041421": 1}),
        ("Haematopoietic and Lymphoid", {
            "PXD003539": 0, "PXD030304": 0, "PXD004701": 0,
            "PXD017199": 0, "PXD041421": 1,
        }),
    ]
    accession_sets = {
        "PXD003539": {"P0001", "P0002"},
        "PXD030304": {"P0002", "P0003"},
        "PXD004701": {"P0001", "P0003"},
        "PXD017199": {"P0001", "P0002", "P0003", "P0004"},
        "PXD041421": {"P0005", "P0006"},  # genuinely new accessions
    }
    tissue_protein_rows = [
        ("Breast", {"PXD003539": 2, "PXD030304": 1, "PXD004701": 2,
                    "PXD017199": 4, "PXD041421": 0}),
        ("Lung", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 0,
                  "PXD017199": 0, "PXD041421": 1}),
    ]
    runs_per_cohort = {
        "PXD003539": 120, "PXD030304": 5798, "PXD004701": 300,
        "PXD017199": 206, "PXD041421": 48,
    }

    svg = tmp_path / "atlas.svg"
    render_atlas(
        headlines,
        cellline_sets,
        tissue_rows,
        accession_sets,
        svg,
        tissue_protein_rows=tissue_protein_rows,
        runs_per_cohort=runs_per_cohort,
    )
    overlap = tmp_path / "atlas_overlap.svg"
    distribution = tmp_path / "atlas_distribution.svg"
    assert overlap.exists() and overlap.stat().st_size > 0
    assert distribution.exists() and distribution.stat().st_size > 0


def test_set_region_sizes_n5_uses_all_five_key() -> None:
    """5-set partition: 2^5 - 1 = 31 region keys. The all-set key shape
    is `"all_five"`, mirroring the all_three / all_four conventions."""
    from analysis.figure_combined_cell_lines_atlas import _set_region_sizes

    sets = {
        "A": {"a", "x", "all"},
        "B": {"b", "x", "all"},
        "C": {"c", "all"},
        "D": {"d", "all"},
        "E": {"e", "all"},
    }
    out = _set_region_sizes(sets, ["A", "B", "C", "D", "E"])
    assert len(out) == 31
    assert "all_five" in out
    assert out["all_five"] == 1   # the "all" accession
    # Partition property: regions sum to the union size.
    union = set()
    for s in sets.values():
        union |= s
    assert sum(out.values()) == len(union)


def test_dataset_headlines_diann_counts_now_target_only(
    tmp_path: Path, monkeypatch
) -> None:
    """`refresh_dataset_headlines()` populates the `DATASET_HEADLINES`
    diann_count with the target-only protein-group count, computed via
    `count_target_protein_groups` on the on-disk pg_matrix.tsv. Synthetic
    pg_matrix fixture: 2 target rows + 1 contam row -> headline = 2."""
    from analysis import figure_combined_cell_lines_atlas as atlas

    # Synthetic pg_matrix file (1 sample column required for shape parity).
    pg = tmp_path / "pg_matrix.tsv"
    pg.write_text(
        "Protein.Group\tProtein.Names\tGenes\tFirst.Protein.Description\t"
        "N.Sequences\tN.Proteotypic.Sequences\tRun_A.d\n"
        "P00001\tn1\tG1\td\t1\t1\t100\n"
        "P00002\tn2\tG2\td\t1\t1\t200\n"
        "CONTAM_P00003\tn3\tG3\td\t1\t1\t300\n"
    )
    # Snapshot the originals so we can restore.
    orig_headlines = atlas.DATASET_HEADLINES["PXD041421"]
    orig_unf = atlas.DATASET_HEADLINES_UNFILTERED.get("PXD041421")
    try:
        # Point only the PXD041421 pg_matrix at the synthetic file; the
        # other 4 cohorts continue to resolve to their on-disk paths and
        # the existing module-load defaults remain accurate for them.
        monkeypatch.setattr(atlas, "PXD041421_PG_MATRIX", pg)
        atlas.refresh_dataset_headlines()
        h = atlas.DATASET_HEADLINES["PXD041421"]
        assert h.diann_count == 2  # 2 target rows, CONTAM dropped
        assert atlas.DATASET_HEADLINES_UNFILTERED["PXD041421"] == 3
    finally:
        # Restore so subsequent tests see the original headline.
        atlas.DATASET_HEADLINES["PXD041421"] = orig_headlines
        if orig_unf is not None:
            atlas.DATASET_HEADLINES_UNFILTERED["PXD041421"] = orig_unf


def test_render_atlas_falls_back_when_runs_metadata_missing(tmp_path: Path) -> None:
    """`render_atlas` must stay resilient when `runs_per_cohort` is
    unset. The breadth-vs-depth scatter that once consumed this metadata
    was dropped from atlas_distribution on 2026-05-29, but the argument
    is retained for API stability — passing None must still write both
    SVGs without error."""
    from analysis.figure_combined_cell_lines_atlas import (
        DATASET_HEADLINES, render_atlas,
    )

    cellline_sets = {
        "PXD003539": {"clA"},
        "PXD030304": {"clB"},
        "PXD004701": {"clC"},
    }
    tissue_rows = [
        ("Breast", {"PXD003539": 1, "PXD030304": 1, "PXD004701": 1}),
    ]
    accession_sets = {
        "PXD003539": {"P0001"},
        "PXD030304": {"P0002"},
        "PXD004701": {"P0003"},
    }
    svg = tmp_path / "atlas.svg"
    render_atlas(
        DATASET_HEADLINES,
        cellline_sets,
        tissue_rows,
        accession_sets,
        svg,
        tissue_protein_rows=[],
        runs_per_cohort=None,  # no metadata — Panel E uses constant dot size
    )
    overlap = tmp_path / "atlas_overlap.svg"
    distribution = tmp_path / "atlas_distribution.svg"
    assert overlap.exists() and overlap.stat().st_size > 0
    assert distribution.exists() and distribution.stat().st_size > 0


# ---------------------------------------------------------------------------
# Panel A — no verbose annotation under no-paper cohorts (post-split cleanup)
# ---------------------------------------------------------------------------

def test_panel_a_no_paper_annotation() -> None:
    """After the atlas-split cleanup, `_render_panel_a_headlines` must
    not emit any explainer text under cohorts whose `paper_count == 0`.
    The x-tick label is the canonical PXDxxx (paper-year) tag from
    `DATASET_LABELS`; the previous "no paper DIA bar — mass-cytometry"
    annotation has been removed."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from analysis.figure_combined_cell_lines_atlas import (
        DatasetHeadline, _render_panel_a_headlines,
    )

    headlines = {
        "PXD003539": DatasetHeadline(
            paper_count=6556, diann_count=6927,
            paper_label="Guo 2019", metric="Protein groups",
        ),
        "PXD017199": DatasetHeadline(
            paper_count=0, diann_count=10572,
            paper_label="", metric="Protein groups",
        ),
        "PXD041421": DatasetHeadline(
            paper_count=0, diann_count=9124,
            paper_label="", metric="Protein groups",
        ),
    }
    fig, ax = plt.subplots(figsize=(6, 4))
    _render_panel_a_headlines(ax, headlines)

    forbidden = ("mass-cytometry", "no paper", "methodological")
    # Inspect every Text child (axis labels, legend, x-tick labels,
    # bar-value annotations) for the deleted explainer phrases.
    text_blob = "\n".join(t.get_text() for t in ax.texts)
    text_blob += "\n" + "\n".join(
        tl.get_text() for tl in ax.get_xticklabels()
    )
    text_blob += "\n" + (ax.get_xlabel() or "")
    text_blob += "\n" + (ax.get_ylabel() or "")
    text_blob_lower = text_blob.lower()
    for phrase in forbidden:
        assert phrase.lower() not in text_blob_lower, (
            f"Panel A still emits {phrase!r}; text was:\n{text_blob}"
        )
    plt.close(fig)
