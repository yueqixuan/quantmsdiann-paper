"""Conservative contaminant / entrapment / decoy filter.

Single canonical predicate used by every count site in the analysis
codebase. A Protein.Group string passes the filter iff **every**
semicolon-separated accession token is a target — no prefix from the
recognised contaminant / entrapment / decoy set.

The conservative policy means rows with even a single prefixed token
are dropped, including mixed groups like `CONTAM_P02768;P02768` (a
contaminant entry sharing peptides with real human albumin). Rationale
in [docs/superpowers/specs/2026-05-21-contaminant-filter-and-pxd041421-design.md]
§1.3: when DIA-NN's protein-grouping has placed a target inside a
contaminant-named group, the inference is ambiguous; conservatively
excluding such rows guards against contamination of the target catalog.

Recognised prefixes are case-sensitive on the exact token below. Two
prefix conventions co-exist in our data:

- `CONTAM_`, `ENTRAP_`, `DECOY_` — cell-line FASTA
  `Homo-sapiens-uniprot-reviewed-entrap-contaminants-202605.fasta`.
- `Cont_` — ProteoBench `ProteoBenchFASTA_*` series.
- `decoy_` — lowercase variant kept for defensive parity (not currently
  observed in our data; low cost to include).
"""
from __future__ import annotations

import re
from pathlib import Path

# Order matters only for the regex alternation — performance is not
# sensitive at our scale (millions of tokens at most).
_FILTER_PREFIXES: tuple[str, ...] = (
    "CONTAM_",
    "Cont_",
    "ENTRAP_",
    "DECOY_",
    "decoy_",
)

_HAS_PREFIX_RE = re.compile(
    "^(?:" + "|".join(re.escape(p) for p in _FILTER_PREFIXES) + ")"
)

# Strips the prefix from the start of a token. Reused by `_clean_token`
# in venn_protein_accessions when the row passes the filter (and the
# accession needs to be normalised, e.g. `sp|P12345|HUMAN` → `P12345`).
_STRIP_PREFIX_RE = re.compile(
    "^(?:" + "|".join(re.escape(p) for p in _FILTER_PREFIXES) + ")+"
)


def _token_has_prefix(token: str) -> bool:
    """Return True iff `token` starts with a recognised
    contaminant / entrapment / decoy prefix."""
    token = token.strip()
    if not token:
        return False
    # When the token is in `sp|ACC|NAME_SPECIES` form, the prefix is
    # on the SwissProt entry name (last `|` field), not the leading
    # accession. We check both the raw token and the post-`|`-strip
    # form so the filter catches `sp|P12345|CONTAM_HORSE_ALB` too.
    if _HAS_PREFIX_RE.match(token):
        return True
    if "|" in token:
        parts = token.split("|")
        # Last field after `|` is the entry name (NAME_SPECIES); check
        # that against the prefix list as well as the middle accession
        # field (in case the prefix decorates either).
        for part in parts[1:]:
            if _HAS_PREFIX_RE.match(part.strip()):
                return True
    return False


def is_target_protein_group(pg_string: str | None) -> bool:
    """Return True iff every semicolon-separated accession in
    `pg_string` is a target (carries no contaminant / entrapment /
    decoy prefix).

    Empty / whitespace-only / None inputs return False — defensive,
    because callers normally use this predicate to filter
    `notna`-checked rows. If a downstream caller needs different
    behaviour for missing values it should handle them explicitly.
    """
    if not pg_string:
        return False
    s = str(pg_string).strip()
    if not s:
        return False
    pieces = [p for p in s.split(";") if p.strip()]
    if not pieces:
        # Degenerate input like `";"` or `";;"` — no actual tokens.
        return False
    for piece in pieces:
        if _token_has_prefix(piece):
            return False
    return True


def strip_known_prefix(token: str) -> str:
    """Strip any recognised contaminant prefix from the start of
    `token`. Used for normalising target accessions ONCE the row has
    passed the filter. Idempotent."""
    return _STRIP_PREFIX_RE.sub("", token)


# ---------------------------------------------------------------------------
# Convenience filters for matrix files
# ---------------------------------------------------------------------------

def count_target_protein_groups(pg_matrix_path: Path) -> tuple[int, int]:
    """Return `(unfiltered_count, target_count)` from a DIA-NN
    `pg_matrix.tsv`-style file. The unfiltered count is the total
    distinct `Protein.Group` rows (one row per group). The target
    count drops rows whose Protein.Group fails `is_target_protein_group`.

    Used by the per-cohort headline writers to record both numbers in
    the audit TSV. Defensively handles missing `Protein.Group`
    columns (returns (0, 0))."""
    import pandas as pd
    try:
        df = pd.read_csv(
            pg_matrix_path, sep="\t", usecols=["Protein.Group"], dtype=str,
        )
    except (FileNotFoundError, OSError, ValueError):
        return (0, 0)
    pgs = df["Protein.Group"].dropna()
    unfiltered = int(len(pgs))
    target = int(pgs.map(is_target_protein_group).sum())
    return (unfiltered, target)


def count_target_precursors(pr_matrix_path: Path) -> tuple[int, int]:
    """Return `(unfiltered_count, target_count)` from a DIA-NN
    `pr_matrix.tsv`-style file. Counts precursor rows; the target
    count drops rows whose Protein.Group fails `is_target_protein_group`.
    """
    import pandas as pd
    try:
        df = pd.read_csv(
            pr_matrix_path, sep="\t", usecols=["Protein.Group"], dtype=str,
        )
    except (FileNotFoundError, OSError, ValueError):
        return (0, 0)
    pgs = df["Protein.Group"].dropna()
    unfiltered = int(len(pgs))
    target = int(pgs.map(is_target_protein_group).sum())
    return (unfiltered, target)
