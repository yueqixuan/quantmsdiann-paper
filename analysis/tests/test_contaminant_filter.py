"""Tests for `analysis.contaminant_filter` — the conservative
contaminant / entrapment / decoy filter applied to every count site.

Conservative policy: a Protein.Group passes the filter iff every
semicolon-separated accession token is a target. Mixed groups
(target + prefixed) are dropped. See
[docs/superpowers/specs/2026-05-21-contaminant-filter-and-pxd041421-design.md]
§1.3 for the rationale."""
from __future__ import annotations

from pathlib import Path


def test_pure_contaminant_returns_false() -> None:
    """A pure-contaminant row must be rejected. The bug we are fixing
    used to strip the prefix and return the bare accession, polluting
    our target catalog with phantom accessions."""
    from analysis.contaminant_filter import is_target_protein_group
    assert not is_target_protein_group("CONTAM_HORSE_ALB")
    assert not is_target_protein_group("CONTAM_P02768-1")


def test_pure_entrap_returns_false() -> None:
    """Entrapment proteins are randomised sequences used to estimate
    FDR — they should never appear in the reported target catalog."""
    from analysis.contaminant_filter import is_target_protein_group
    assert not is_target_protein_group("ENTRAP_Q3KR16")
    assert not is_target_protein_group("ENTRAP_RANDOMSEQ")


def test_pure_decoy_returns_false() -> None:
    """DIA-NN's q-value filter normally excludes decoys, but
    defensively the consumer must drop any row whose Protein.Group
    has a DECOY_ token."""
    from analysis.contaminant_filter import is_target_protein_group
    assert not is_target_protein_group("DECOY_P02768")
    assert not is_target_protein_group("decoy_P02768")


def test_lowercase_cont_recognised() -> None:
    """ProteoBench FASTAs use `Cont_` (capitalised, single underscore)
    while cell-line FASTAs use `CONTAM_` (all-caps). The filter must
    recognise both — the previous regex was upper-case only and let
    ~1,200 `Cont_` precursor rows per benchmark dataset leak into our
    F1a precursor counts."""
    from analysis.contaminant_filter import is_target_protein_group
    assert not is_target_protein_group("Cont_P60712")
    assert not is_target_protein_group("Cont_Q3T052")


def test_mixed_prefix_target_dropped_under_conservative() -> None:
    """`CONTAM_P02768;P02768` — a contaminant entry sharing peptides
    with real human albumin — must be dropped under the conservative
    policy chosen on 2026-05-21. The previous behaviour stripped the
    prefix and kept P02768, which is the LIBERAL interpretation.
    Switching to CONSERVATIVE drops the row entirely."""
    from analysis.contaminant_filter import is_target_protein_group
    assert not is_target_protein_group("CONTAM_P02768;P02768")
    assert not is_target_protein_group("P02768;CONTAM_P02768")
    assert not is_target_protein_group("ENTRAP_Q1;Q2")


def test_pure_target_passes() -> None:
    """Target-only rows must pass the filter unchanged."""
    from analysis.contaminant_filter import is_target_protein_group
    assert is_target_protein_group("P02768")
    assert is_target_protein_group("P02768;Q9NQ29")
    assert is_target_protein_group("P02768;Q9NQ29;O15156")


def test_swissprot_style_token_with_prefix_in_entry_name() -> None:
    """Some Protein.Group cells carry the SwissProt
    `sp|ACCESSION|ENTRY_NAME` shape. The contaminant flag may live on
    the entry-name field, not on the leading token. The filter must
    catch this case too."""
    from analysis.contaminant_filter import is_target_protein_group
    assert not is_target_protein_group("sp|P12345|CONTAM_HORSE_ALB")
    # A clean SwissProt entry passes through.
    assert is_target_protein_group("sp|P12345|ALBU_HUMAN")


def test_empty_or_whitespace_returns_false() -> None:
    """Defensive: empty / whitespace / None inputs return False.
    Callers normally exclude these via dropna() upstream; this
    behaviour just keeps the predicate well-defined."""
    from analysis.contaminant_filter import is_target_protein_group
    assert not is_target_protein_group(None)
    assert not is_target_protein_group("")
    assert not is_target_protein_group("   ")
    assert not is_target_protein_group(";")


def test_strip_known_prefix_is_idempotent() -> None:
    """`strip_known_prefix` is the per-token normaliser used AFTER
    a row has passed the conservative filter — it cleans up
    accessions that came through. Idempotent: applying twice gives
    the same result as once."""
    from analysis.contaminant_filter import strip_known_prefix
    assert strip_known_prefix("P02768") == "P02768"
    assert strip_known_prefix("CONTAM_P02768") == "P02768"
    assert strip_known_prefix("CONTAM_CONTAM_P02768") == "P02768"
    assert strip_known_prefix("Cont_P60712") == "P60712"
    assert strip_known_prefix("ENTRAP_Q1") == "Q1"
    # Already-clean accession is unchanged after two passes.
    one = strip_known_prefix("CONTAM_P02768")
    two = strip_known_prefix(one)
    assert one == two == "P02768"


# ---------------------------------------------------------------------------
# Matrix-level counters
# ---------------------------------------------------------------------------

def test_count_target_protein_groups_writes_pre_and_post(tmp_path: Path) -> None:
    """The matrix-level counter must return both unfiltered and
    filtered counts so the audit TSV can carry both for transparency.
    Fixture: 5 PG rows of which 2 are pure-contaminant, 1 mixed (also
    dropped under conservative), 2 pure-target."""
    from analysis.contaminant_filter import count_target_protein_groups
    p = tmp_path / "pg_matrix.tsv"
    p.write_text(
        "Protein.Group\tProtein.Ids\tProtein.Names\n"
        "P02768\tP02768\tALBU_HUMAN\n"
        "Q9NQ29\tQ9NQ29\tLUC7L_HUMAN\n"
        "CONTAM_HORSE_ALB\tCONTAM_HORSE_ALB\tALBU_HORSE\n"
        "CONTAM_P02768;P02768\tCONTAM_P02768;P02768\tALBU_HUMAN\n"
        "ENTRAP_Q3KR16\tENTRAP_Q3KR16\tQ3KR16\n"
    )
    unfiltered, target = count_target_protein_groups(p)
    assert unfiltered == 5
    assert target == 2


def test_count_target_protein_groups_handles_missing_file(tmp_path: Path) -> None:
    """A missing pg_matrix returns (0, 0) defensively rather than
    raising — used by audit-writers that may be invoked on a partial
    dataset."""
    from analysis.contaminant_filter import count_target_protein_groups
    unfiltered, target = count_target_protein_groups(
        tmp_path / "does_not_exist.tsv",
    )
    assert unfiltered == 0
    assert target == 0


def test_count_target_precursors_filters_at_pg_level(tmp_path: Path) -> None:
    """Same predicate applied to a pr_matrix-shaped fixture; verifies
    the matrix-level wrapper works on the precursor-row schema (just
    a different fixture, same Protein.Group column)."""
    from analysis.contaminant_filter import count_target_precursors
    p = tmp_path / "pr_matrix.tsv"
    p.write_text(
        "Protein.Group\tStripped.Sequence\n"
        "P02768\tAAAR\n"
        "Cont_P60712\tBBBR\n"
        "Q9NQ29\tCCCR\n"
        "ENTRAP_Q1;Q2\tDDDR\n"
    )
    unfiltered, target = count_target_precursors(p)
    assert unfiltered == 4
    assert target == 2
