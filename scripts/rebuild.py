#!/usr/bin/env python3
"""rebuild — reproduce every manuscript figure (and the PDFs) from original
sources, in one command.

This is a single self-contained module: every former ``analysis/*.py`` and
``analysis/plexDIA/*.py`` helper/figure/data-prep module has been inlined here
in a flat namespace. Colliding top-level names from different source modules
were disambiguated with a ``_<module>__<name>`` prefix; the canonical copy of
each shared helper/constant keeps its bare name. The two module-as-namespace
aliases (``fs`` for the old figure_style, ``acc`` for figure_proteobench_accuracy)
are re-pointed at this module via ``sys.modules[__name__]`` after all defs.

All data is pulled from the public PRIDE FTP (deposited DIA-NN reports and the
original published matrices); nothing here re-runs DIA-NN from raw vendor files
(that is the upstream HPC reanalysis, out of scope). Every figure number obeys
the filter rule documented in ``methods.md`` §1.

Usage:
    python -m scripts.rebuild --all            # data prep -> figures -> PDFs
    python -m scripts.rebuild --figures-only    # skip data prep
    python -m scripts.rebuild --data-only       # only the data-prep stages
    python -m scripts.rebuild --only NAME [...]  # run specific stage(s)
    python -m scripts.rebuild --no-pdf          # skip the LaTeX/PDF build
    python -m scripts.rebuild --list            # print all stages + provenance
    python -m scripts.rebuild --keep-going      # don't stop on first failure (default)
    python -m scripts.rebuild --fail-fast       # stop at the first failing stage

Every stage runs in-process (no ``python -m`` subprocess). The atlas (Fig S13)
needs a numpy<2 environment (see methods.md / environment.yml); run rebuild
inside the ``quantmsdiann`` conda env to reproduce it.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import argparse
import collections
import csv
import io
import json
import os
import re
import subprocess
import sys
import time
import types
import urllib.request
import zipfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Iterator

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import requests
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

# Repo root: scripts/ is one level under the repo root (same depth analysis/ was).
REPO = Path(__file__).resolve().parents[1]

# Module-as-namespace aliases. The originals did `from analysis import
# figure_style as fs` and `from analysis import figure_proteobench_accuracy as
# acc`; with everything merged here, those namespaces ARE this module. These are
# bound up front (figure_style/figure_proteobench_accuracy bodies are inlined
# first) so module-level `fs.apply_house_style()` etc. resolve at import time.
fs = sys.modules[__name__]
acc = sys.modules[__name__]
figure_proteobench_accuracy = sys.modules[__name__]


# ======================================================================
# inlined from analysis/figure_style.py
# ======================================================================

"""Shared figure house-style for the quantmsdiann manuscript.

One import unifies typography, spines, palette, and number formatting across
every ``figure_*.py`` script so the paper reads as a single coherent set.

Usage (top of each figure script, after the matplotlib import)::

    from analysis import figure_style as fs
    fs.apply_house_style()

Then use the shared palette/formatters instead of inline hex + ad-hoc labels::

    ax.bar(x, y, color=fs.COMPARISON["quantmsdiann"])
    fs.kfmt_axis(ax.yaxis)          # 140000 -> "140k" on the tick labels
    ax.set_ylabel("Protein groups (1% FDR)")
    fs.despine(ax)

Design notes
------------
* Palette is **Okabe-Ito** (colour-blind-safe, prints legibly in grayscale).
* The three DIA-NN versions keep the original light->dark "age" semantics but
  on colour-blind-safe hues; the enterprise build is a distinct accent hue.
* Two-way "original vs quantmsdiann" comparisons use neutral grey vs one accent
  (the scheme the cohort figures already used well — now centralised).
"""
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
OKABE_ITO = {'black': '#000000', 'orange': '#E69F00', 'sky_blue': '#56B4E9', 'bluish_green': '#009E73', 'yellow': '#F0E442', 'blue': '#0072B2', 'vermillion': '#D55E00', 'reddish_purple': '#CC79A7', 'grey': '#9E9E9E'}
VERSION_COLORS = {'v1_8_1': OKABE_ITO['sky_blue'], '1_8_1': OKABE_ITO['sky_blue'], 'v2_5_1': OKABE_ITO['blue'], '2_5_1': OKABE_ITO['blue'], 'v2_5_1_enterprise': OKABE_ITO['vermillion'], '2_5_1_enterprise': OKABE_ITO['vermillion']}
COMPARISON = {'original': OKABE_ITO['grey'], 'quantmsdiann': OKABE_ITO['blue']}
CATEGORICAL_CYCLE = [OKABE_ITO['blue'], OKABE_ITO['vermillion'], OKABE_ITO['bluish_green'], OKABE_ITO['orange'], OKABE_ITO['sky_blue'], OKABE_ITO['reddish_purple'], '#882255', '#117733', '#888888', '#332288', '#999933', '#44AA99']
INSTRUMENT_COLORS = {'Orbitrap Astral': OKABE_ITO['blue'], 'Orbitrap Eclipse': '#332288', 'Orbitrap Exploris 480': OKABE_ITO['sky_blue'], 'Q Exactive': '#44AA99', 'TripleTOF 5600': OKABE_ITO['orange'], 'TripleTOF 6600': OKABE_ITO['vermillion'], 'ZenoTOF 7600': OKABE_ITO['bluish_green'], 'timsTOF SCP': OKABE_ITO['reddish_purple'], 'timsTOF Pro': '#882255', 'timsTOF HT': '#999933', 'unknown': OKABE_ITO['grey']}

def categorical_colors(n):
    """Return *n* colours from the colour-blind-leaning categorical cycle."""
    if n <= len(CATEGORICAL_CYCLE):
        return CATEGORICAL_CYCLE[:n]
    out = list(CATEGORICAL_CYCLE)
    while len(out) < n:
        out.append(CATEGORICAL_CYCLE[len(out) % len(CATEGORICAL_CYCLE)])
    return out[:n]

def apply_house_style():
    """Apply the manuscript-wide matplotlib rcParams. Idempotent."""
    mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'], 'font.size': 16, 'axes.titlesize': 18, 'axes.titleweight': 'bold', 'axes.labelsize': 17, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 14, 'legend.title_fontsize': 14, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8, 'xtick.direction': 'out', 'ytick.direction': 'out', 'legend.frameon': False, 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300, 'svg.fonttype': 'none'})

def kfmt(value, _pos=None):
    """Abbreviate large counts: 134000 -> '134k', 1300 -> '1.3k', 950 -> '950'.

    Safe as a matplotlib ``FuncFormatter`` (ignores the position arg).
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return str(value)
    a = abs(v)
    if a >= 1000:
        s = f'{v / 1000:.1f}'.rstrip('0').rstrip('.')
        return f'{s}k'
    return f'{v:.0f}'

def kfmt_axis(axis):
    """Apply :func:`kfmt` to a matplotlib axis (e.g. ``ax.yaxis``)."""
    axis.set_major_formatter(FuncFormatter(kfmt))

def commafmt(value, _pos=None):
    """Thousands-separated integer: 9542 -> '9,542'. Safe as a FuncFormatter."""
    try:
        return f'{int(round(float(value))):,}'
    except (TypeError, ValueError):
        return str(value)

def despine(ax, top=True, right=True, left=False, bottom=False):
    """Hide the requested spines (top/right by default)."""
    for name, off in (('top', top), ('right', right), ('left', left), ('bottom', bottom)):
        if off:
            ax.spines[name].set_visible(False)

def style_boxplot(bp, color=None, median_color=None):
    """Restyle a matplotlib ``ax.boxplot`` result to the house style.

    Kills the default cyan boxes / red medians / loud fliers in favour of a
    neutral box with an accent median and discreet outliers.
    """
    color = color or OKABE_ITO['blue']
    median_color = median_color or OKABE_ITO['vermillion']
    for box in bp.get('boxes', []):
        box.set(color=color, linewidth=1.0)
    for whisk in bp.get('whiskers', []):
        whisk.set(color=color, linewidth=0.9)
    for cap in bp.get('caps', []):
        cap.set(color=color, linewidth=0.9)
    for med in bp.get('medians', []):
        med.set(color=median_color, linewidth=1.4)
    for fl in bp.get('fliers', []):
        fl.set(marker='o', markersize=2.5, markerfacecolor='none', markeredgecolor='#999999', alpha=0.6)


# ======================================================================
# inlined from analysis/collect_paper_numbers.py
# ======================================================================

"""collect_paper_numbers — aggregate every number cited in the manuscript into
ONE auditable file, generated (never hand-typed) from the figure-data TSVs that
the figure scripts emit under the filter rule (methods.md §1).

Outputs:
  * data/paper_numbers.tsv          -- key, value, source, note (human audit)
  * paper/generated_numbers.tex     -- \\newcommand macros the prose can \\input

Run after the figures (their data TSVs are the inputs):
    python -m analysis.collect_paper_numbers
or via:  python -m scripts.rebuild --all   (runs as the final data stage)

Each source is read defensively: a missing input warns and is skipped, so the
file always reflects whatever has been regenerated so far.
"""
import sys
from pathlib import Path
import pandas as pd
_collect_paper_numbers__REPO = Path(__file__).resolve().parents[1]
OUT_TSV = _collect_paper_numbers__REPO / 'data' / 'paper_numbers.tsv'
OUT_TEX = _collect_paper_numbers__REPO / 'paper' / 'generated_numbers.tex'
MODULE_LABEL = {'ProteoBench_Module_7': 'M7', 'PXD049412': 'M9', 'PXD062685': 'M5', 'PXD070049': 'M10'}

def _pct(old: float, new: float) -> float:
    return 100.0 * (new - old) / old if old else 0.0

class Numbers:
    """Accumulator: (key, value, source, note). `key` (if given) becomes a LaTeX
    macro \\num<Key>; rows without a key are audit-only (TSV)."""

    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str, str]] = []

    def add(self, key: str, value, source: str, note: str='') -> None:
        self.rows.append((key, str(value), source, note))

    def to_tsv(self, path: Path) -> None:
        pd.DataFrame(self.rows, columns=['key', 'value', 'source', 'note']).to_csv(path, sep='\t', index=False)

    def to_tex(self, path: Path) -> None:
        lines = ['% AUTO-GENERATED by analysis/collect_paper_numbers.py — do not edit.', '% \\input{generated_numbers} then use e.g. \\numMSevenProtGainEnt.']
        for key, value, source, _ in self.rows:
            if not key:
                continue
            lines.append(f'\\newcommand{{\\num{key}}}{{{value}}}')
        path.write_text('\n'.join(lines) + '\n')

def _read(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        print(f'WARN: missing {path.relative_to(_collect_paper_numbers__REPO)}', file=sys.stderr)
        return None
    return pd.read_csv(path, sep='\t')

def collect_benchmarks(N: Numbers) -> None:
    df = _read(_collect_paper_numbers__REPO / 'data' / 'quantmsdiann_benchmarks' / 'report_counts.tsv')
    if df is None:
        return
    df = df.set_index(['dataset', 'version'])
    for module, label in MODULE_LABEL.items():
        try:
            v18, v25, vent = (df.loc[module, v] for v in ('v1_8_1', 'v2_5_1', 'v2_5_1_enterprise'))
        except KeyError:
            continue
        src = 'report_counts.tsv'
        N.add(f'{label}PrecGlobalV18', int(v18['prec_global']), src)
        N.add(f'{label}PrecGlobalEnt', int(vent['prec_global']), src)
        N.add(f'{label}ProtGlobalV18', int(v18['prot_global']), src)
        N.add(f'{label}ProtGlobalEnt', int(vent['prot_global']), src)
        N.add(f'{label}PrecGainEnt', f"{_pct(v18['prec_global'], vent['prec_global']):.0f}", src, 'precursor %gain v1.8.1 -> 2.5.1-ent (Lib.Q.Value)')
        N.add(f'{label}ProtGainEnt', f"{_pct(v18['prot_global'], vent['prot_global']):.0f}", src, 'protein-group %gain v1.8.1 -> 2.5.1-ent (Lib.PG.Q.Value)')

def collect_reanalysis(N: Numbers) -> None:
    df = _read(_collect_paper_numbers__REPO / 'analysis' / 'figures' / 'reanalysis' / 'data' / 'reanalysis_improvement.tsv')
    if df is None:
        return
    for _, r in df.iterrows():
        ds = str(r['dataset'])
        gain = _pct(float(r['original']), float(r['reanalysis']))
        N.add('', int(r['original']), 'reanalysis_improvement.tsv', f"{ds} {r['label']} original ({r.get('metric', '')})")
        N.add('', int(r['reanalysis']), 'reanalysis_improvement.tsv', f"{ds} {r['label']} reanalysis (+{gain:.0f}%)")

def collect_single_cell(N: Numbers) -> None:
    df = _read(_collect_paper_numbers__REPO / 'data' / 'single_cell' / 'sc_totals.tsv')
    if df is not None:
        for _, r in df.iterrows():
            tag = f"{r['dataset']} {r['version']}"
            N.add('', int(r['precursors']), 'sc_totals.tsv', f'{tag} precursors (global)')
            N.add('', int(r['proteins']), 'sc_totals.tsv', f'{tag} protein groups (global)')
    pc = _read(_collect_paper_numbers__REPO / 'data' / 'single_cell' / 'mv_per_cell.tsv')
    if pc is not None and {'dataset', 'version', 'n_proteins'}.issubset(pc.columns):
        med = pc.groupby(['dataset', 'version'])['n_proteins'].median()
        for (ds, ver), v in med.items():
            N.add('', int(v), 'mv_per_cell.tsv', f'{ds} {ver} median protein groups/cell')

def collect_plexdia(N: Numbers) -> None:
    df = _read(_collect_paper_numbers__REPO / 'analysis' / 'figures' / 'plexDIA' / 'MSV000093870' / 'comparison_counts.tsv')
    if df is None:
        return
    df = df.set_index('metric')
    for metric in df.index:
        for col in df.columns:
            N.add('', df.loc[metric, col], 'plexDIA/comparison_counts.tsv', f'plexDIA {metric} ({col})')

def collect_phospho(N: Numbers) -> None:
    df = _read(_collect_paper_numbers__REPO / 'data' / 'phospho' / 'phospho_counts.tsv')
    if df is None:
        return
    for _, r in df.iterrows():
        tag = f"{r['dataset']} {r['version']}"
        for col in ('phosphopeptides', 'sites_classI', 'sites_all'):
            if col in r:
                N.add('', int(r[col]), 'phospho_counts.tsv', f'{tag} {col}')

def collect_atlas(N: Numbers) -> None:
    df = _read(_collect_paper_numbers__REPO / 'analysis' / 'figures' / 'combined' / 'data' / 'combined_counts.tsv')
    if df is None:
        return
    for _, r in df.iterrows():
        N.add('', r.get('count', ''), 'combined_counts.tsv', f"{r.get('metric', '')} | {r.get('source', '')}")

def collect_paper_numbers_main() -> int:
    N = Numbers()
    collect_benchmarks(N)
    collect_reanalysis(N)
    collect_single_cell(N)
    collect_plexdia(N)
    collect_phospho(N)
    collect_atlas(N)
    N.to_tsv(OUT_TSV)
    N.to_tex(OUT_TEX)
    n_macros = sum((1 for k, *_ in N.rows if k))
    print(f'wrote {OUT_TSV.relative_to(_collect_paper_numbers__REPO)} ({len(N.rows)} numbers) and {OUT_TEX.relative_to(_collect_paper_numbers__REPO)} ({n_macros} macros)')
    return 0


# ======================================================================
# inlined from analysis/contaminant_filter.py
# ======================================================================

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
import re
from pathlib import Path
_FILTER_PREFIXES: tuple[str, ...] = ('CONTAM_', 'Cont_', 'ENTRAP_', 'DECOY_', 'decoy_')
_HAS_PREFIX_RE = re.compile('^(?:' + '|'.join((re.escape(p) for p in _FILTER_PREFIXES)) + ')')
_STRIP_PREFIX_RE = re.compile('^(?:' + '|'.join((re.escape(p) for p in _FILTER_PREFIXES)) + ')+')

def _token_has_prefix(token: str) -> bool:
    """Return True iff `token` starts with a recognised
    contaminant / entrapment / decoy prefix."""
    token = token.strip()
    if not token:
        return False
    if _HAS_PREFIX_RE.match(token):
        return True
    if '|' in token:
        parts = token.split('|')
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
    pieces = [p for p in s.split(';') if p.strip()]
    if not pieces:
        return False
    for piece in pieces:
        if _token_has_prefix(piece):
            return False
    return True

def strip_known_prefix(token: str) -> str:
    """Strip any recognised contaminant prefix from the start of
    `token`. Used for normalising target accessions ONCE the row has
    passed the filter. Idempotent."""
    return _STRIP_PREFIX_RE.sub('', token)


# ======================================================================
# inlined from analysis/count_matrix.py
# ======================================================================

"""Reproducible row counting for DIA-NN ``*_pg_matrix.tsv`` / ``*_pr_matrix.tsv``.

Per methods.md §1 ("Original / deposited side of reanalysis comparisons"),
DIA-NN matrix files are ALREADY q-filtered count matrices and carry no
q-value columns, so there is no q-filter left to apply. The single
reproducible number is the count of matrix ROWS that are quantified:

  * for a ``*_pg_matrix.tsv`` a row is one protein group;
  * for a ``*_pr_matrix.tsv`` a row is one precursor.

Counting rule (no other filter is admissible):
  * NO contaminant/target filter -- ``CONTAM_`` / ``Cont_`` / ``ENTRAP_``
    prefixed rows are counted like any other;
  * NO positive-quantity filter -- a literal ``0`` quantity is "quantified"
    and is counted; only an empty / ``NA`` cell means "not measured";
  * a row counts iff it has at least one non-empty quantity across the
    sample columns (matrices can contain fully-empty rows);
  * decoys are not present in matrices, so there is nothing to drop.

The sample columns are every column that is not one of the fixed DIA-NN
metadata columns (``PG_METADATA`` / ``PR_METADATA`` below).
"""
from pathlib import Path
import pandas as pd
PG_METADATA = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description', 'N.Sequences', 'N.Proteotypic.Sequences']
PR_METADATA = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Precursor.Id']
_PG_REQUIRED = ['Protein.Group', 'First.Protein.Description']
_PR_REQUIRED = ['Protein.Group', 'Precursor.Id', 'First.Protein.Description']

def count_matrix_rows(matrix_path: Path, metadata_cols: list[str]) -> int:
    """Number of quantified rows in a DIA-NN pg/pr count matrix.

    A row is "quantified" iff at least one of its sample columns (every column
    not in ``metadata_cols`` that is actually present in the matrix) holds a
    non-empty, non-``NA`` value. Zeros count; contaminant/entrapment rows count.

    Raises ``ValueError`` if a required anchor column for the chosen layout is
    missing from the header (guards against mis-parsed / wrong-kind matrices).
    """
    required = _PR_REQUIRED if 'Precursor.Id' in metadata_cols else _PG_REQUIRED
    df = pd.read_csv(matrix_path, sep='\t', dtype=str)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"matrix {matrix_path} is missing metadata column(s): {', '.join(missing)}")
    sample_cols = [c for c in df.columns if c not in metadata_cols]
    if not sample_cols:
        return 0
    samples = df[sample_cols].replace({'NA': pd.NA, '': pd.NA})
    quantified = samples.notna().any(axis=1)
    return int(quantified.sum())


# ======================================================================
# inlined from analysis/count_report_ids.py
# ======================================================================

"""Count precursors and protein groups from the DIA-NN *report*, not the matrices.

The `*_matrix.tsv` outputs bake in `--matrix-spec-q` (0.05 run-specific) and,
because the quantmsdiann pipeline sets `--qvalue` to 0.01 for DIA-NN 1.8.1 but
0.05 for 2.5.1/enterprise, matrix row counts are filtered at *different*
run-specific q-values per version and are NOT comparable across versions.

This module counts identifications directly from the per-precursor report
(`diann_report.parquet` for DIA-NN >= 2.x, `diann_report.tsv` for 1.8.1).

Filter rule (Vadim review, 2026-06-21 -- scientific-correctness requirement for
submission). Reported quantities fall into two classes, each with EXACTLY one
admissible filter and nothing else (no contaminant/target filter, no
positive-quantity filter, zeros counted):

  * Per-run numbers (within a single run / cell):
      - protein groups: `PG.Q.Value <= 0.01` only.
      - precursors:     `Q.Value <= 0.01` only (flat, all versions, per
        methods.md §1; the matrices' baked per-version `--matrix-spec-q` does
        not apply to the report-side count, which is gated at 0.01 throughout).
  * Global numbers (dataset union / totals):
      - protein groups: `Lib.PG.Q.Value <= 0.01` only.
      - precursors:     `Lib.Q.Value <= 0.01` only.

Emitted metrics:
  * `prec_min1` / `prec_min3` -- precursors identified in >= 1 / >= 3 runs at the
    per-run `Q.Value` cut-off (replicate-reproducibility view).
  * `prec_global` -- distinct precursors with `Lib.Q.Value <= 0.01` (global total).
  * `prot_global` -- distinct protein groups with `Lib.PG.Q.Value <= 0.01`.
  * `prot_perrun_avg` / `prot_complete` -- average protein groups per run, and
    protein groups quantified in *every* run, at run-specific `PG.Q.Value <= 0.01`.
  * `peptides` -- distinct stripped peptides among `Lib.Q.Value`-passing precursors.
  * `prot_2pep` -- global protein groups (`Lib.PG.Q.Value`) with >= 2 distinct
    stripped peptides among `Lib.Q.Value`-passing precursors.

Decoys (`Decoy == 1`) are dropped everywhere: this removes the FDR null model and
is not one of the forbidden "filters".

The reports are multi-GB; this is run on the downloaded reports (from the public
PRIDE FTP) and its small output, `report_counts.tsv`, is staged into
`data/quantmsdiann_benchmarks/` and consumed by
`figure_quantmsdiann_benchmarks_vs_proteobench.py`.

Usage:
    python -m analysis.count_report_ids         --results-root /path/to/quantmsdiann_results         --out report_counts.tsv
"""
import argparse
import sys
from pathlib import Path
import pandas as pd
Q_THRESHOLD = 0.01
PRECURSOR_Q = {'v1_8_1': 0.01, 'v2_5_1': 0.01, 'v2_5_1_enterprise': 0.01}
DEFAULT_PRECURSOR_Q = 0.01
_NEEDED_COLS = ['Run', 'Precursor.Id', 'Protein.Group', 'Stripped.Sequence', 'Q.Value', 'PG.Q.Value', 'Lib.Q.Value', 'Lib.PG.Q.Value', 'Decoy']
DATASET_MODULES = ('ProteoBench_Module_7', 'PXD049412', 'PXD062685', 'PXD070049')
_count_report_ids__VERSIONS = ('v1_8_1', 'v2_5_1', 'v2_5_1_enterprise')

def _load_report(report_dir: Path) -> pd.DataFrame:
    """Read the DIA-NN report from `report_dir`, preferring the parquet
    (DIA-NN >= 2.x) and falling back to the classic `diann_report.tsv`
    (1.8.1). Only the columns we need are read."""
    parquet = report_dir / 'diann_report.parquet'
    tsv = report_dir / 'diann_report.tsv'
    if parquet.exists():
        import pyarrow.parquet as pq
        have = set(pq.ParquetFile(parquet).schema_arrow.names)
        cols = [c for c in _NEEDED_COLS if c in have]
        return pq.read_table(parquet, columns=cols).to_pandas()
    if tsv.exists():
        return pd.read_csv(tsv, sep='\t', usecols=lambda c: c in _NEEDED_COLS)
    raise FileNotFoundError(f'no diann_report.parquet or .tsv in {report_dir}')

def count_report(df: pd.DataFrame, precursor_q: float=DEFAULT_PRECURSOR_Q) -> dict[str, int]:
    """Compute precursor and protein-group counts from a DIA-NN report frame
    under the Vadim filter rule (see module docstring): per-run protein groups
    on `PG.Q.Value` only; per-run precursors on `Q.Value` only; global protein
    groups on `Lib.PG.Q.Value` only; global precursors on `Lib.Q.Value` only.
    No contaminant/target filter, no positive-quantity filter (zeros counted)."""
    if 'Decoy' in df.columns:
        df = df[df['Decoy'] == 0]
    passing = df[df['Q.Value'] <= precursor_q]
    n_runs_per_prec = passing.groupby('Precursor.Id')['Run'].nunique()
    prec_min1 = int((n_runs_per_prec >= 1).sum())
    prec_min3 = int((n_runs_per_prec >= 3).sum())
    prec_global = 0
    if 'Lib.Q.Value' in df.columns:
        prec_global = int(df.loc[df['Lib.Q.Value'] <= Q_THRESHOLD, 'Precursor.Id'].nunique())
    prot_global = 0
    if 'Lib.PG.Q.Value' in df.columns:
        global_pgs = df.loc[df['Lib.PG.Q.Value'] <= Q_THRESHOLD, 'Protein.Group'].dropna().unique()
        prot_global = int(len(global_pgs))
        # methods.md §4 protein-inference sanity check: among unique global
        # protein groups at 1% q-value, no more than 1.2% may carry multiple
        # accessions (';' in Protein.Group). Exceeding it signals inflated
        # inference (e.g. DIA-NN 1.8.1 without --relaxed-prot-inf) and is flagged.
        if len(global_pgs):
            multi_frac = 100.0 * sum(';' in str(pg) for pg in global_pgs) / len(global_pgs)
            if multi_frac > 1.2:
                print(f"WARNING [methods.md §4]: {multi_frac:.2f}% of global protein "
                      f"groups are multi-accession (>1.2% cap) — possible inference "
                      f"inflation; counts may not be comparable.", file=sys.stderr)
    prot_perrun_avg = prot_complete = 0
    if 'PG.Q.Value' in df.columns:
        n_runs = df['Run'].nunique()
        pg = df[df['PG.Q.Value'] <= Q_THRESHOLD][['Run', 'Protein.Group']].drop_duplicates()
        per_run = pg.groupby('Run')['Protein.Group'].nunique()
        prot_perrun_avg = int(round(per_run.mean())) if len(per_run) else 0
        in_n_runs = pg.groupby('Protein.Group')['Run'].nunique()
        prot_complete = int((in_n_runs == n_runs).sum()) if n_runs else 0
    peptides = prot_2pep = 0
    if 'Stripped.Sequence' in df.columns and 'Lib.Q.Value' in df.columns:
        dq = df[df['Lib.Q.Value'] <= Q_THRESHOLD]
        peptides = int(dq['Stripped.Sequence'].nunique())
        if 'Lib.PG.Q.Value' in df.columns:
            dp = dq[dq['Lib.PG.Q.Value'] <= Q_THRESHOLD]
            per_prot = dp.groupby('Protein.Group')['Stripped.Sequence'].nunique()
            prot_2pep = int((per_prot >= 2).sum())
    return {'prec_min1': prec_min1, 'prec_min3': prec_min3, 'prec_global': prec_global, 'prot_global': prot_global, 'prot_perrun_avg': prot_perrun_avg, 'prot_complete': prot_complete, 'peptides': peptides, 'prot_2pep': prot_2pep}

def count_report_ids_main(argv: list[str] | None=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results-root', required=True, type=Path, help='quantmsdiann_results dir holding <module>/<version>/quant_tables')
    ap.add_argument('--out', required=True, type=Path, help='output report_counts.tsv path')
    ap.add_argument('--results-suffix', default='', help="suffix for the per-version dir, e.g. '_relaxed' reads <module>/<version-without-v><suffix>/quant_tables (the --relaxed-prot-inf re-run dirs '1_8_1_relaxed' etc.)")
    args = ap.parse_args(argv)
    rows: list[dict] = []
    for module in DATASET_MODULES:
        for version in _count_report_ids__VERSIONS:
            if args.results_suffix:
                base = version[1:] if version.startswith('v') else version
                ver_dir = f'{base}{args.results_suffix}'
            else:
                ver_dir = version
            rdir = args.results_root / module / ver_dir / 'quant_tables'
            try:
                df = _load_report(rdir)
            except FileNotFoundError as exc:
                print(f'WARN: {exc}', file=sys.stderr)
                continue
            c = count_report(df, precursor_q=PRECURSOR_Q.get(version, DEFAULT_PRECURSOR_Q))
            c.update(dataset=module, version=version)
            rows.append(c)
            print(f"{module} {version}: {c['prec_global']:,} prec (global) / {c['prot_global']:,} proteins (global, Lib.q<=0.01)")
    cols = ['dataset', 'version', 'prec_min1', 'prec_min3', 'prec_global', 'prot_global', 'prot_perrun_avg', 'prot_complete', 'peptides', 'prot_2pep']
    pd.DataFrame(rows)[cols].to_csv(args.out, sep='\t', index=False)
    print(f'wrote {args.out}')
    return 0


# ======================================================================
# inlined from analysis/figure_combined_cell_lines_atlas.py
# ======================================================================

"""Combined cell-lines atlas figure for the quantmsdiann manuscript.

Integrates the five independent cell-line reanalyses already shipped
(PXD003539 NCI-60, PXD030304 ProCan-DepMapSanger 949 lines, PXD004701
BC 76 lines, PXD017199 Tognetti 67 breast lines, PXD041421 Wang 2023)
into two paper-ready figures that position the quantmsdiann pipeline as
a single uniform tool covering broad cancer cell-line / tissue / proteome
space. The original single A–H grid was split in two for readability.

`atlas_overlap.svg` (protein-accession overlap):
- A (top, wide): UpSet plot of UniProt accessions detected by
  quantmsdiann across the five cohorts, with per-cohort headline counts;
  the 5-set UpSet replaces the unreadable 4-set Venn.
- B / C (bottom row): dataset-level reproducibility (paper vs
  quantmsdiann headline counts; PXD017199 / PXD041421 have no paper
  headline so only their quantmsdiann bar is drawn) and the
  detection-count histogram (how many cohorts each protein is seen in,
  giving the pan-cohort core).

`atlas_distribution.svg` (per-tissue coverage):
- A (top, wide): unified-axis pan-cancer tissue coverage, stacked
  horizontal bars per cohort (cell-line counts).
- B (bottom, wide): per-tissue unique-protein stacked bars (parallel to
  A; protein counts instead of cell-line counts).
The breadth-vs-depth scatter (former Panel C/E) was dropped 2026-05-29
and the Expression Atlas overlap (former Panel H) was removed earlier —
it duplicates analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg.

Reuses parsers/helpers from the per-dataset scripts via imports — no
duplicated logic, no new downloads. Reads pre-cached JSON caches for
PXD030304 / PXD004701 protein sets (don't re-stream 33 GB parquets) and
the PXD003539 / PXD017199 pr_matrix.tsv (already on disk) for their
protein sets.
"""
import json
import sys
from dataclasses import dataclass
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
_figure_combined_cell_lines_atlas__REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = _figure_combined_cell_lines_atlas__REPO_ROOT / 'data'
_figure_combined_cell_lines_atlas__FIGURES_DIR = _figure_combined_cell_lines_atlas__REPO_ROOT / 'analysis' / 'figures' / 'combined'
PXD003539_SDRF = DATA_ROOT / 'PXD003539' / 'PXD003539.sdrf.tsv'
PXD003539_PR_MATRIX = DATA_ROOT / 'PXD003539' / 'diann_report.pr_matrix.tsv'
PXD003539_REPORT_PARQUET = DATA_ROOT / 'PXD003539' / 'diann_report.parquet'
PXD003539_COUNTS_TSV = _figure_combined_cell_lines_atlas__REPO_ROOT / 'analysis' / 'figures' / 'PXD003539' / 'data' / 'counts.tsv'
PXD030304_SDRF = DATA_ROOT / 'PXD030304' / 'PXD030304.sdrf.tsv'
PXD030304_TISSUE_MAPPING = DATA_ROOT / 'PXD030304' / 'mapping_file_averaged.txt'
PXD030304_PROTEIN_JSON = DATA_ROOT / 'PXD030304' / 'diann_per_tissue_procan_filter.json'
PXD030304_PG_MATRIX = DATA_ROOT / 'PXD030304' / 'diann_report.pg_matrix.tsv'
PXD030304_COUNTS_TSV = _figure_combined_cell_lines_atlas__REPO_ROOT / 'analysis' / 'figures' / 'PXD030304' / 'data' / 'counts.tsv'
PXD004701_SDRF = DATA_ROOT / 'PXD004701' / 'PXD004701.sdrf.tsv'
PXD004701_PROTEIN_JSON = DATA_ROOT / 'PXD004701' / 'diann_per_subtype_consistency_filter.json'
PXD003539_PROTEIN_JSON = DATA_ROOT / 'PXD003539' / 'diann_report_protein_counts.json'
PXD004701_PG_MATRIX = DATA_ROOT / 'PXD004701' / 'diann_report.pg_matrix.tsv'
PXD004701_COUNTS_TSV = _figure_combined_cell_lines_atlas__REPO_ROOT / 'analysis' / 'figures' / 'PXD004701' / 'data' / 'counts.tsv'
PXD017199_SDRF = DATA_ROOT / 'PXD017199' / 'PXD017199.sdrf.tsv'
PXD017199_PR_MATRIX = DATA_ROOT / 'PXD017199' / 'diann_report.pr_matrix.tsv'
PXD017199_PG_MATRIX = DATA_ROOT / 'PXD017199' / 'diann_report.pg_matrix.tsv'
PXD041421_SDRF = DATA_ROOT / 'PXD041421' / 'PXD041421.sdrf.tsv'
PXD041421_PR_MATRIX = DATA_ROOT / 'PXD041421' / 'diann_report.pr_matrix.tsv'
PXD041421_PG_MATRIX = DATA_ROOT / 'PXD041421' / 'diann_report.pg_matrix.tsv'
E_PROT_73_TSV = DATA_ROOT / 'E-PROT-73-query-results.tsv'

# Public-FTP sources for the bulk cell-line DIA-NN count matrices the reanalysis,
# venn and atlas figures read. Three live under quantmsdiann-benchmarks/cell-lines
# (with a version dir), two under the absolute-expression-2.0 collection (no
# version dir). ensure_cell_line_matrices() auto-fetches any that are missing so
# the rebuild is self-sufficient (no manual staging).
_BENCH_CL_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/cell-lines'
_ABSEXP_CL_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines'
# Per-cohort quant_tables base URLs (three on quantmsdiann-benchmarks with a
# version dir; two on the absolute-expression-2.0 collection, no version dir).
_CELL_LINE_QT_BASE: dict[str, str] = {
    'PXD003539': f'{_BENCH_CL_BASE}/PXD003539/v2_5_1_enterprise/quant_tables',
    'PXD030304': f'{_BENCH_CL_BASE}/PXD030304/v2_5_1/quant_tables',
    'PXD004701': f'{_BENCH_CL_BASE}/PXD004701/v2_5_1/quant_tables',
    'PXD017199': f'{_ABSEXP_CL_BASE}/PXD017199/quant_tables',
    'PXD041421': f'{_ABSEXP_CL_BASE}/PXD041421/quant_tables',
}
# DIA-NN count matrices each cohort's figures actually read (shared across
# stages, NOT purged). Per-cohort to avoid fetching the unused 2 GB ProCan
# pr_matrix: pxd030304 reads only pg, and atlas reads ProCan/Sun from their
# cached JSONs. The full diann_report.parquet is fetched only on demand
# (with_report) and IS purged after the stage that needs it.
_PR, _PG, _UG = 'diann_report.pr_matrix.tsv', 'diann_report.pg_matrix.tsv', 'diann_report.unique_genes_matrix.tsv'
_CELL_LINE_MATRIX_FILES = {
    'PXD003539': (_PR, _PG, _UG),
    'PXD030304': (_PG,),
    'PXD004701': (_PG,),
    'PXD017199': (_PR, _PG),
    'PXD041421': (_PR, _PG),
}

def ensure_cell_line_matrices(*accessions: str, with_report: bool = False) -> None:
    """Download missing cell-line DIA-NN quant_tables files from the public FTP
    so the rebuild is self-sufficient (no manual staging). Fetches the count
    matrices each given cohort's figures read (all cohorts if none given);
    with_report=True also pulls the large diann_report.parquet. No-op for files
    already present; files absent on the FTP are skipped."""
    for acc in (accessions or tuple(_CELL_LINE_QT_BASE)):
        base = _CELL_LINE_QT_BASE.get(acc)
        if not base:
            continue
        wanted = list(_CELL_LINE_MATRIX_FILES.get(acc, (_PR, _PG, _UG)))
        if with_report:
            wanted.append('diann_report.parquet')
        for fn in wanted:
            dest = DATA_ROOT / acc / fn
            if dest.exists() and dest.stat().st_size > 0:
                continue
            try:
                print(f'Fetching {acc}/{fn} ...', file=sys.stderr)
                download_if_missing(f'{base}/{fn}', dest)
            except Exception as exc:
                print(f'  (skipping {acc}/{fn}: {exc})', file=sys.stderr)
DATASET_COLORS = {'PXD003539': '#1b9e77', 'PXD030304': '#7570b3', 'PXD004701': '#d95f02', 'PXD017199': '#e6ab02', 'PXD041421': '#8da0cb'}
DATASET_LABELS = {'PXD003539': 'PXD003539\n(Guo 2019)', 'PXD030304': 'PXD030304\n(ProCan 2022)', 'PXD004701': 'PXD004701\n(Sun 2023)', 'PXD017199': 'PXD017199\n(Tognetti 2021)', 'PXD041421': 'PXD041421\n(Wang 2023)'}

@dataclass(frozen=True)
class DatasetHeadline:
    paper_count: int
    diann_count: int
    paper_label: str
    metric: str
PXD003539_PAPER_COUNT = 6556
PXD030304_PAPER_COUNT = 8498
PXD004701_PAPER_COUNT = 6091
_DEFAULT_DIANN_COUNTS = {'PXD003539': 7018, 'PXD030304': 9606, 'PXD004701': 7961, 'PXD017199': 10713, 'PXD041421': 9183}
_GLOBAL_PG_METRIC = 'Protein groups (global rule: matrix rows / Lib.PG.Q.Value <= 0.01; no filters)'
DATASET_HEADLINES: dict[str, DatasetHeadline] = {'PXD003539': DatasetHeadline(paper_count=PXD003539_PAPER_COUNT, diann_count=_DEFAULT_DIANN_COUNTS['PXD003539'], paper_label='Guo 2019 (OpenSWATH)', metric=_GLOBAL_PG_METRIC), 'PXD030304': DatasetHeadline(paper_count=PXD030304_PAPER_COUNT, diann_count=_DEFAULT_DIANN_COUNTS['PXD030304'], paper_label='ProCan 2022', metric=_GLOBAL_PG_METRIC), 'PXD004701': DatasetHeadline(paper_count=PXD004701_PAPER_COUNT, diann_count=_DEFAULT_DIANN_COUNTS['PXD004701'], paper_label='Sun 2023', metric=_GLOBAL_PG_METRIC), 'PXD017199': DatasetHeadline(paper_count=0, diann_count=_DEFAULT_DIANN_COUNTS['PXD017199'], paper_label='', metric=_GLOBAL_PG_METRIC), 'PXD041421': DatasetHeadline(paper_count=0, diann_count=_DEFAULT_DIANN_COUNTS['PXD041421'], paper_label='', metric=_GLOBAL_PG_METRIC)}

def _report_global_protein_groups(report_parquet: Path) -> int:
    """Global protein-group count from a DIA-NN report parquet under the
    methods.md §1 rule: distinct `Protein.Group` with
    `Lib.PG.Q.Value <= 0.01`, decoys (`Decoy == 1`) dropped. No
    contaminant/target filter, no positive-quantity filter. Returns 0 if
    the file or required columns are absent."""
    if not report_parquet.exists():
        return 0
    import pyarrow.parquet as pq
    have = set(pq.ParquetFile(report_parquet).schema_arrow.names)
    needed = {'Protein.Group', 'Lib.PG.Q.Value'}
    if not needed <= have:
        return 0
    cols = [c for c in ('Protein.Group', 'Lib.PG.Q.Value', 'Decoy') if c in have]
    df = pq.read_table(report_parquet, columns=cols).to_pandas()
    if 'Decoy' in df.columns:
        df = df[df['Decoy'] == 0]
    return int(df.loc[df['Lib.PG.Q.Value'] <= 0.01, 'Protein.Group'].nunique())

def refresh_dataset_headlines() -> None:
    """Populate `DATASET_HEADLINES` from the on-disk inputs under the
    methods.md §1 global protein-group rule (2026-06-21 Vadim review):

      * deposited DIA-NN count matrices (PXD030304 / PXD004701 /
        PXD017199 / PXD041421) -> `count_matrix_rows()`: a row counts iff
        it has >=1 non-empty sample value (zeros counted), with NO
        contaminant/target filter;
      * the PXD003539 report parquet -> distinct `Protein.Group` at
        `Lib.PG.Q.Value <= 0.01`, decoys dropped.

    There is a single, filter-free number per cohort — the previous
    "unfiltered vs target-only" duality has been removed. Cohorts whose
    input file is missing keep the import-time default (which is the same
    rule, precomputed).

    Tests that don't need the from-disk numbers can either skip this
    call (the module-load defaults are accurate) or call it themselves to
    exercise the read path."""
    matrix_sources: dict[str, Path] = {'PXD030304': PXD030304_PG_MATRIX, 'PXD004701': PXD004701_PG_MATRIX, 'PXD017199': PXD017199_PG_MATRIX, 'PXD041421': PXD041421_PG_MATRIX}
    for ds, pg in matrix_sources.items():
        if not (pg.exists() and pg.stat().st_size > 0):
            continue
        n = count_matrix_rows(pg, PG_METADATA)
        if n <= 0:
            continue
        h = DATASET_HEADLINES[ds]
        DATASET_HEADLINES[ds] = DatasetHeadline(paper_count=h.paper_count, diann_count=n, paper_label=h.paper_label, metric=h.metric)
    n3539 = _report_global_protein_groups(PXD003539_REPORT_PARQUET)
    if n3539 > 0:
        h = DATASET_HEADLINES['PXD003539']
        DATASET_HEADLINES['PXD003539'] = DatasetHeadline(paper_count=h.paper_count, diann_count=n3539, paper_label=h.paper_label, metric=h.metric)

def cell_lines_from_sdrf(sdrf_path: Path, cell_line_col: str='characteristics[cell line]') -> set[str]:
    """Return the set of normalised cell-line names from any of the three
    SDRFs. Normalisation via `normalize_cell_line` strips the `NCI-` prefix
    and all non-alphanumeric characters, uppercases — so 'CCRF-CEM',
    'NCI-H226', 'Hs-578-T' collide across datasets with their alternative
    spellings."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    if cell_line_col not in df.columns:
        raise ValueError(f'SDRF {sdrf_path} missing required column: {cell_line_col!r}')
    out: set[str] = set()
    for raw in df[cell_line_col].fillna(''):
        norm = normalize_cell_line(raw)
        if norm:
            out.add(norm)
    return out

def _compute_runs_per_cohort() -> dict[str, int]:
    """Return the number of MS runs per cohort from the cached SDRF
    files. Used by Panel E (breadth-vs-depth scatter) to scale the dot
    size for each cohort. Missing SDRFs yield 0 (defensively — the
    panel handles zero by falling back to a constant dot size)."""
    out: dict[str, int] = {}
    for ds, sdrf in (('PXD003539', PXD003539_SDRF), ('PXD030304', PXD030304_SDRF), ('PXD004701', PXD004701_SDRF), ('PXD017199', PXD017199_SDRF), ('PXD041421', PXD041421_SDRF)):
        if not Path(sdrf).exists():
            out[ds] = 0
            continue
        try:
            df = pd.read_csv(sdrf, sep='\t', dtype=str)
            out[ds] = int(len(df))
        except (FileNotFoundError, OSError, ValueError, pd.errors.EmptyDataError):
            out[ds] = 0
    return out
PXD003539_TISSUE_RULES: list[tuple[str, str]] = [('leukemia', 'Haematopoietic and Lymphoid'), ('lymphoma', 'Haematopoietic and Lymphoid'), ('myeloma', 'Haematopoietic and Lymphoid'), ('blood', 'Haematopoietic and Lymphoid'), ('bone marrow', 'Haematopoietic and Lymphoid'), ('glioblastoma', 'Central Nervous System'), ('gliosarcoma', 'Central Nervous System'), ('astrocytoma', 'Central Nervous System'), ('central nervous system', 'Central Nervous System'), ('brain', 'Central Nervous System'), ('breast', 'Breast'), ('colon', 'Large Intestine'), ('colorectal', 'Large Intestine'), ('large intestine', 'Large Intestine'), ('lung', 'Lung'), ('non-small cell', 'Lung'), ('mesothelioma', 'Lung'), ('melanoma', 'Skin'), ('skin', 'Skin'), ('ovarian', 'Ovary'), ('ovary', 'Ovary'), ('prostate', 'Prostate'), ('renal', 'Kidney'), ('kidney', 'Kidney')]

def harmonise_pxd003539_tissue(disease_text: str | None, organism_part: str | None=None) -> str | None:
    """Map a PXD003539 (disease, organism_part) pair to one of the 28 ProCan
    tissue categories. Returns None if no rule matches (uncategorised lines
    are dropped from Panel C, not silently bucketed into Other tissue).

    The rules cover the 9 NCI-60 cancer types; the matched ProCan-axis
    tissue is the most specific bucket NCI-60 maps onto. Notes:
    - "leukemia / lymphoma / myeloma" maps to ProCan's
      `Haematopoietic and Lymphoid` (NCI-60 has 6 leukemia lines under
      Blood, all of which are present as Haematopoietic & Lymphoid in
      ProCan's mapping_file_averaged.txt).
    - "Pleural epithelioid mesothelioma" (NCI-60 MESO line under Lung)
      maps to ProCan `Lung` — the only mesothelioma category ProCan
      retains.
    - Brain / CNS / glioblastoma / gliosarcoma / astrocytoma all map to
      ProCan `Central Nervous System`."""
    disease = (disease_text or '').strip().lower()
    organism = (organism_part or '').strip().lower()
    haystack = f'{disease} {organism}'
    for needle, tissue in PXD003539_TISSUE_RULES:
        if needle in haystack:
            return tissue
    return None

def cell_line_tissue_pxd003539(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, ProCan_tissue]` for PXD003539.

    Each row of the SDRF carries (cell_line, disease, organism_part);
    rows of the same cell line have identical disease/organism so we
    collapse to a deduplicated mapping. Cell lines whose disease/organism
    text doesn't match any rule are dropped from the result."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'characteristics[disease]', 'characteristics[organism part]']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f'PXD003539 SDRF {sdrf_path} missing columns: {missing}')
    out: dict[str, str] = {}
    for cell, disease, organism in zip(df['characteristics[cell line]'], df['characteristics[disease]'], df['characteristics[organism part]']):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        tissue = harmonise_pxd003539_tissue(disease, organism)
        if tissue is None:
            continue
        out[norm] = tissue
    return out

def cell_line_tissue_pxd030304(sdrf_path: Path, mapping_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, ProCan_tissue]` for PXD030304.

    The ProCan figshare `mapping_file_averaged.txt` is the canonical source
    of truth for the 28-tissue axis; we use it directly. The SDRF gives us
    the cell lines actually present in the dataset (the figshare mapping
    covers all 949 cell lines anyway, but going through the SDRF keeps the
    contract identical across the three datasets)."""
    cl_to_tissue = parse_procan_mapping(mapping_path)
    norm_lookup: dict[str, str] = {}
    for cell, tissue in cl_to_tissue.items():
        n = normalize_cell_line(cell)
        if n:
            norm_lookup[n] = tissue
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    if 'characteristics[cell line]' not in df.columns:
        raise ValueError(f'PXD030304 SDRF {sdrf_path} missing characteristics[cell line]')
    out: dict[str, str] = {}
    for cell in df['characteristics[cell line]'].fillna(''):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        tissue = norm_lookup.get(norm)
        if tissue is None:
            continue
        out[norm] = tissue
    return out

def cell_line_tissue_pxd004701(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, 'Breast']` for PXD004701.

    All 76 cell lines in PXD004701 are breast-cancer-derived; the BC
    subtype split is internal to that dataset and irrelevant for the
    unified-tissue axis."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    if 'characteristics[cell line]' not in df.columns:
        raise ValueError(f'PXD004701 SDRF {sdrf_path} missing characteristics[cell line]')
    out: dict[str, str] = {}
    for cell in df['characteristics[cell line]'].fillna(''):
        norm = normalize_cell_line(cell)
        if norm:
            out[norm] = 'Breast'
    return out

def cell_line_tissue_pxd017199(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, tissue]` for PXD017199.

    PXD017199 is essentially all breast-cancer-derived plus a handful of
    "normal" mammary-epithelial lines (184A1, 184B5, HBL100, MCF10A,
    MCF10F, MCF12A). The disease column distinguishes them:

    - characteristics[disease] == "normal" -> `Healthy (Non-cancer)`
      (matches the tissue category used by `cell_line_tissue_pxd003539`
      for non-tumour rows in the unified ProCan axis).
    - any other disease value -> `Breast`.

    If a single cell line has multiple disease rows (it shouldn't in
    PXD017199 but the SDRF doesn't enforce uniqueness), the first
    matching row wins; ties between Breast and normal are dominated by
    whichever appears first in the SDRF."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'characteristics[disease]']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f'PXD017199 SDRF {sdrf_path} missing columns: {missing}')
    out: dict[str, str] = {}
    for cell, disease in zip(df['characteristics[cell line]'], df['characteristics[disease]']):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        if norm in out:
            continue
        d = (disease or '').strip().lower()
        if d == 'normal':
            out[norm] = 'Healthy (Non-cancer)'
        else:
            out[norm] = 'Breast'
    return out

def cell_line_tissue_pxd041421(sdrf_path: Path) -> dict[str, str]:
    """Return `dict[normalised_cell_line, tissue]` for PXD041421.

    PXD041421 (Wang 2023 batch-effect testbed) carries 2 cell lines:
    A549 (lung adenocarcinoma) and K562 (CML blast-phase). To stay on
    the unified ProCan 28-tissue axis already used by every other
    cohort, the labels match those produced by
    `cell_line_tissue_pxd003539` for the same diseases:

    - A549 -> `Lung`
    - K562 -> `Haematopoietic and Lymphoid`

    The spec calls these "Lung Cancer" and "Leukemia" — the actual
    axis labels are the ProCan headers and remain unchanged so the
    Panel C / F stacked bars merge correctly across cohorts."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    if 'characteristics[cell line]' not in df.columns:
        raise ValueError(f'PXD041421 SDRF {sdrf_path} missing characteristics[cell line]')
    LINE_TO_TISSUE = {'A549': 'Lung', 'K562': 'Haematopoietic and Lymphoid'}
    out: dict[str, str] = {}
    for cell in df['characteristics[cell line]'].fillna(''):
        norm = normalize_cell_line(cell)
        if not norm:
            continue
        if norm in out:
            continue
        raw = (cell or '').strip()
        tissue = LINE_TO_TISSUE.get(raw) or LINE_TO_TISSUE.get(normalize_cell_line(raw))
        if tissue is None:
            continue
        out[norm] = tissue
    return out

def combined_tissue_table(per_dataset: dict[str, dict[str, str]]) -> list[tuple[str, dict[str, int]]]:
    """Combine per-dataset cell-line-to-tissue mappings into a single
    `[(tissue, {dataset: cell_line_count})]` list, sorted by total
    cell-line count descending. Tissues with zero contribution from every
    dataset are dropped (always — the input only carries mapped lines)."""
    counts: dict[str, dict[str, int]] = {}
    for dataset, cl_to_tissue in per_dataset.items():
        tissue_counts: dict[str, int] = {}
        for tissue in cl_to_tissue.values():
            tissue_counts[tissue] = tissue_counts.get(tissue, 0) + 1
        for tissue, n in tissue_counts.items():
            counts.setdefault(tissue, {ds: 0 for ds in per_dataset})[dataset] = n
    out = sorted(counts.items(), key=lambda kv: -sum(kv[1].values()))
    return out

def pxd003539_protein_accessions(pr_matrix_path: Path) -> set[str]:
    """Return the set of UniProt accessions detected by quantmsdiann in
    PXD003539. Reads `diann_report.pr_matrix.tsv` (already cached locally
    at ~67 MB) and collects accessions from `Protein.Group` strings on rows
    with at least one non-NA quantification across the 120 runs.

    Uses `extract_accessions_diann` so the accession-normalisation semantics
    match the PXD030304 / PXD004701 caches (semicolon split, isoform suffix
    stripped, CONTAM_/ENTRAP_/DECOY_ prefix removed)."""
    df = pd.read_csv(pr_matrix_path, sep='\t', dtype=str)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'PXD003539 pr_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    if not sample_cols:
        raise ValueError('PXD003539 pr_matrix has no per-run sample columns')
    quantified = df[df[sample_cols].notna().any(axis=1)]
    accessions: set[str] = set()
    for pg in quantified['Protein.Group']:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions

def _accessions_from_json_cache(json_path: Path) -> set[str]:
    """Helper for the PXD030304 / PXD004701 cached JSONs. Both caches are
    `{group_key: [Protein.Group, ...]}` (group_key is tissue / subtype).
    We union the Protein.Group values across all groups and extract
    accessions."""
    with open(json_path, encoding='utf-8') as fh:
        payload = json.load(fh)
    pg_set: set[str] = set()
    for vs in payload.values():
        pg_set.update(vs)
    accessions: set[str] = set()
    for pg in pg_set:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions

def pxd030304_protein_accessions(json_path: Path) -> set[str]:
    """Read `diann_per_tissue_procan_filter.json` and return the union of
    UniProt accessions across the 28 tissues. The JSON is the cached output
    of `proteins_per_tissue_quantmsdiann_procan_filter` and contains
    Protein.Group strings as deposited in the parquet."""
    return _accessions_from_json_cache(json_path)

def pxd004701_protein_accessions(json_path: Path) -> set[str]:
    """Read `diann_per_subtype_consistency_filter.json` and return the union
    of UniProt accessions across the 3 BC subtypes. The JSON is the cached
    output of `proteins_per_subtype_quantmsdiann_consistency_filter`."""
    return _accessions_from_json_cache(json_path)

def pxd003539_accessions_per_cell_line(pr_matrix_path: Path, sdrf_path: Path) -> dict[str, set[str]]:
    """Return `dict[normalised_cell_line, set[UniProt accession]]` for
    PXD003539. For each cell line, the value is the union of accessions
    quantified in any of its associated runs (i.e. any non-NA cell in the
    pr_matrix column matching the SDRF `comment[data file]` rewritten to
    `.mzML`).

    Unlike `pxd003539_protein_accessions`, this helper preserves
    per-cell-line granularity so Panel E (rarefaction) and Panel F
    (per-tissue protein counts) can roll up across arbitrary partitions.

    Cell lines with no matching pr_matrix column or no quantified
    precursors are dropped silently."""
    import re
    df = pd.read_csv(pr_matrix_path, sep='\t', dtype=str)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'PXD003539 pr_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    sdrf = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'comment[data file]']
    sdrf_missing = [c for c in needed if c not in sdrf.columns]
    if sdrf_missing:
        raise ValueError(f'PXD003539 SDRF missing required columns: {sdrf_missing}')
    col_to_cell: dict[str, str] = {}
    for cell, data_file in zip(sdrf['characteristics[cell line]'], sdrf['comment[data file]']):
        if not isinstance(data_file, str) or not data_file:
            continue
        mzml = re.sub('\\.wiff$', '.mzML', data_file)
        norm = normalize_cell_line(cell)
        if norm:
            col_to_cell[mzml] = norm
    out: dict[str, set[str]] = {}
    for col in sample_cols:
        cell = col_to_cell.get(col)
        if cell is None:
            continue
        mask = df[col].notna()
        if not mask.any():
            continue
        bucket = out.setdefault(cell, set())
        for pg in df.loc[mask, 'Protein.Group']:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
    return out

def pxd017199_protein_accessions(pr_matrix_path: Path) -> set[str]:
    """Return the set of UniProt accessions detected by quantmsdiann in
    PXD017199. Mirrors `pxd003539_protein_accessions`: reads
    `diann_report.pr_matrix.tsv` (~193 MB, already on disk under
    `data/PXD017199/`) and collects accessions from `Protein.Group`
    strings on rows with at least one non-NA quantification across the
    206 runs.

    Uses `extract_accessions_diann` for accession-normalisation parity
    with the other datasets."""
    df = pd.read_csv(pr_matrix_path, sep='\t', dtype=str)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'PXD017199 pr_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    if not sample_cols:
        raise ValueError('PXD017199 pr_matrix has no per-run sample columns')
    quantified = df[df[sample_cols].notna().any(axis=1)]
    accessions: set[str] = set()
    for pg in quantified['Protein.Group']:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions

def pxd017199_accessions_per_cell_line(pr_matrix_path: Path, sdrf_path: Path) -> dict[str, set[str]]:
    """Return `dict[normalised_cell_line, set[UniProt accession]]` for
    PXD017199. Mirrors `pxd003539_accessions_per_cell_line` — but the
    PXD017199 run filenames are `.raw` (not `.wiff`) so no rewrite is
    needed; the SDRF's `comment[data file]` is the exact pr_matrix
    column name.

    Per-cell-line accession sets are the union across that line's runs;
    cell lines with no matching pr_matrix columns or no quantified
    precursors are silently dropped."""
    df = pd.read_csv(pr_matrix_path, sep='\t', dtype=str)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'PXD017199 pr_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    sdrf = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'comment[data file]']
    sdrf_missing = [c for c in needed if c not in sdrf.columns]
    if sdrf_missing:
        raise ValueError(f'PXD017199 SDRF missing required columns: {sdrf_missing}')
    col_to_cell: dict[str, str] = {}
    for cell, data_file in zip(sdrf['characteristics[cell line]'], sdrf['comment[data file]']):
        if not isinstance(data_file, str) or not data_file:
            continue
        norm = normalize_cell_line(cell)
        if norm:
            col_to_cell[data_file] = norm
    out: dict[str, set[str]] = {}
    for col in sample_cols:
        cell = col_to_cell.get(col)
        if cell is None:
            continue
        mask = df[col].notna()
        if not mask.any():
            continue
        bucket = out.setdefault(cell, set())
        for pg in df.loc[mask, 'Protein.Group']:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
    return out

def pxd041421_protein_accessions(pr_matrix_path: Path) -> set[str]:
    """Return the set of UniProt accessions detected by quantmsdiann in
    PXD041421. Mirrors `pxd017199_protein_accessions` — reads
    `diann_report.pr_matrix.tsv` (48 runs, ~9k protein groups, 64 MB on
    disk) and collects accessions from `Protein.Group` strings on rows
    with at least one non-NA quantification.

    Uses `extract_accessions_diann` so the conservative
    contaminant/entrapment/decoy filter is applied at the row level
    (2026-05-21 spec)."""
    df = pd.read_csv(pr_matrix_path, sep='\t', dtype=str)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'PXD041421 pr_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    if not sample_cols:
        raise ValueError('PXD041421 pr_matrix has no per-run sample columns')
    quantified = df[df[sample_cols].notna().any(axis=1)]
    accessions: set[str] = set()
    for pg in quantified['Protein.Group']:
        if isinstance(pg, str):
            accessions.update(extract_accessions_diann(pg))
    return accessions

def pxd041421_accessions_per_cell_line(pr_matrix_path: Path, sdrf_path: Path) -> dict[str, set[str]]:
    """Return `dict[normalised_cell_line, set[UniProt accession]]` for
    PXD041421. Mirrors `pxd017199_accessions_per_cell_line`.

    Run-column / data-file alignment quirks for PXD041421:
    - pr_matrix column headers end in `.d` (Bruker timsTOF folder)
    - SDRF `comment[data file]` ends in `.d.zip` (the FTP archive name).
    We strip the trailing `.zip` to map SDRF -> pr_matrix column.

    Per-cell-line accession sets are the union across that line's runs;
    cell lines with no matching pr_matrix columns or no quantified
    precursors are silently dropped."""
    df = pd.read_csv(pr_matrix_path, sep='\t', dtype=str)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'PXD041421 pr_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    sdrf = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'comment[data file]']
    sdrf_missing = [c for c in needed if c not in sdrf.columns]
    if sdrf_missing:
        raise ValueError(f'PXD041421 SDRF missing required columns: {sdrf_missing}')
    import re as _re
    col_to_cell: dict[str, str] = {}
    for cell, data_file in zip(sdrf['characteristics[cell line]'], sdrf['comment[data file]']):
        if not isinstance(data_file, str) or not data_file:
            continue
        col_name = _re.sub('\\.zip$', '', data_file, flags=_re.IGNORECASE)
        norm = normalize_cell_line(cell)
        if norm:
            col_to_cell[col_name] = norm
    out: dict[str, set[str]] = {}
    for col in sample_cols:
        cell = col_to_cell.get(col)
        if cell is None:
            continue
        mask = df[col].notna()
        if not mask.any():
            continue
        bucket = out.setdefault(cell, set())
        for pg in df.loc[mask, 'Protein.Group']:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
    return out

def pxd003539_gene_symbols(pr_matrix_path: Path) -> set[str]:
    """Return the set of gene symbols quantified in PXD003539 by
    quantmsdiann. Each `Genes` cell may carry one or more ';'-separated
    symbols; only rows with at least one non-NA quantification across the
    runs are counted (same filter as `pxd003539_protein_accessions`)."""
    df = pd.read_csv(pr_matrix_path, sep='\t', dtype=str)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'PXD003539 pr_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    quantified = df[df[sample_cols].notna().any(axis=1)]
    out: set[str] = set()
    for g in quantified['Genes']:
        if not isinstance(g, str):
            continue
        for piece in g.split(';'):
            sym = piece.strip()
            if sym:
                out.add(sym)
    return out

def _per_group_accessions_from_json(json_path: Path) -> dict[str, set[str]]:
    """Read `{group_key: [Protein.Group, ...]}` JSON cache and return
    `{group_key: set[UniProt accession]}`."""
    with open(json_path, encoding='utf-8') as fh:
        payload = json.load(fh)
    out: dict[str, set[str]] = {}
    for key, pgs in payload.items():
        bucket: set[str] = set()
        for pg in pgs:
            if isinstance(pg, str):
                bucket.update(extract_accessions_diann(pg))
        out[key] = bucket
    return out

def pxd030304_accessions_per_tissue(json_path: Path) -> dict[str, set[str]]:
    """Return `{tissue: set[accession]}` for PXD030304 (28 tissues)."""
    return _per_group_accessions_from_json(json_path)

def pxd004701_accessions_per_subtype(json_path: Path) -> dict[str, set[str]]:
    """Return `{subtype: set[accession]}` for PXD004701 (3 BC subtypes)."""
    return _per_group_accessions_from_json(json_path)

def expression_atlas_gene_set(tsv_path: Path) -> set[str]:
    """Parse the Expression Atlas E-PROT-73 query-results TSV and return the
    set of unique `Gene Name` values. The file's first lines are
    `#`-prefixed comments and `pd.read_csv(..., comment='#')` skips them.

    Returns an empty set if the file is missing — callers (Panel H) treat
    that as an explainer-only render path."""
    if not tsv_path.exists():
        return set()
    df = pd.read_csv(tsv_path, sep='\t', comment='#', dtype=str)
    if 'Gene Name' not in df.columns:
        raise ValueError(f"E-PROT-73 TSV {tsv_path} missing required column 'Gene Name'")
    out: set[str] = set()
    for g in df['Gene Name'].fillna(''):
        sym = g.strip()
        if sym:
            out.add(sym)
    return out

def rarefaction_curve(groups: dict[str, set[str]], *, n_permutations: int=50, seed: int=42) -> list[float]:
    """Return the average cumulative-union size as groups are accumulated
    one-by-one. Order is randomised over `n_permutations` permutations
    using `numpy.random.default_rng(seed)` for reproducibility.

    The returned list has length `len(groups)`; entry `i` is the mean
    size of the union over the first `i+1` groups across permutations.
    For 0 or 1 groups the curve is deterministic (no averaging needed)."""
    import numpy as np
    keys = list(groups.keys())
    n = len(keys)
    if n == 0:
        return []
    if n == 1:
        return [float(len(groups[keys[0]]))]
    rng = np.random.default_rng(seed)
    accum = np.zeros(n, dtype=float)
    for _ in range(n_permutations):
        order = list(rng.permutation(n))
        union: set[str] = set()
        for i, idx in enumerate(order):
            union.update(groups[keys[idx]])
            accum[i] += len(union)
    return list(accum / n_permutations)

def per_tissue_union_accessions(tissue_order: list[str], cl_tissue_pxd003539: dict[str, str], accessions_per_cell_line_pxd003539: dict[str, set[str]], accessions_per_tissue_pxd030304: dict[str, set[str]], accessions_pxd004701: set[str], cl_tissue_pxd017199: dict[str, str] | None=None, accessions_per_cell_line_pxd017199: dict[str, set[str]] | None=None, cl_tissue_pxd041421: dict[str, str] | None=None, accessions_per_cell_line_pxd041421: dict[str, set[str]] | None=None) -> list[tuple[str, set[str]]]:
    """Per-tissue UNION of UniProt accessions across all cohorts.

    Returns `[(tissue, union_accession_set)]` in the input order. The
    union answers `how many distinct target-only proteins were
    detected in this tissue across the whole atlas`, which is the
    biologically meaningful number — the previous `per_tissue_protein_counts`
    summed per-cohort counts and so double-counted proteins detected
    in multiple cohorts (e.g., a Breast protein seen in PXD004701,
    PXD017199 and PXD003539 contributed three times to the stack)."""
    tissue_to_cls: dict[str, list[str]] = {}
    for cl, tissue in cl_tissue_pxd003539.items():
        tissue_to_cls.setdefault(tissue, []).append(cl)
    tissue_to_cls_17199: dict[str, list[str]] = {}
    have_17199 = cl_tissue_pxd017199 is not None and accessions_per_cell_line_pxd017199 is not None
    if have_17199:
        for cl, tissue in cl_tissue_pxd017199.items():
            tissue_to_cls_17199.setdefault(tissue, []).append(cl)
    tissue_to_cls_41421: dict[str, list[str]] = {}
    have_41421 = cl_tissue_pxd041421 is not None and accessions_per_cell_line_pxd041421 is not None
    if have_41421:
        for cl, tissue in cl_tissue_pxd041421.items():
            tissue_to_cls_41421.setdefault(tissue, []).append(cl)
    rows: list[tuple[str, set[str]]] = []
    for tissue in tissue_order:
        union: set[str] = set()
        for cl in tissue_to_cls.get(tissue, []):
            union.update(accessions_per_cell_line_pxd003539.get(cl, set()))
        union.update(accessions_per_tissue_pxd030304.get(tissue, set()))
        if tissue == 'Breast':
            union.update(accessions_pxd004701)
        if have_17199:
            for cl in tissue_to_cls_17199.get(tissue, []):
                union.update(accessions_per_cell_line_pxd017199.get(cl, set()))
        if have_41421:
            for cl in tissue_to_cls_41421.get(tissue, []):
                union.update(accessions_per_cell_line_pxd041421.get(cl, set()))
        rows.append((tissue, union))
    return rows

def per_tissue_protein_counts(tissue_order: list[str], cl_tissue_pxd003539: dict[str, str], accessions_per_cell_line_pxd003539: dict[str, set[str]], accessions_per_tissue_pxd030304: dict[str, set[str]], accessions_pxd004701: set[str], cl_tissue_pxd017199: dict[str, str] | None=None, accessions_per_cell_line_pxd017199: dict[str, set[str]] | None=None, cl_tissue_pxd041421: dict[str, str] | None=None, accessions_per_cell_line_pxd041421: dict[str, set[str]] | None=None) -> list[tuple[str, dict[str, int]]]:
    """Roll up per-(tissue, dataset) protein-set sizes onto the unified
    tissue axis (the same `tissue_order` used by Panel C).

    For each tissue:
      - PXD003539: union of per-cell-line accession sets for cell lines
        mapped to this tissue.
      - PXD030304: size of the cached per-tissue accession set.
      - PXD004701: total accession count if and only if the tissue is
        `Breast` (the dataset is breast-only); else 0.
      - PXD017199 (optional, only when both mappings are passed): union
        of per-cell-line accessions for cell lines mapped to this
        tissue (so the 5-6 mammary-normal lines contribute to
        `Healthy (Non-cancer)` and the rest to `Breast`).

    Returns `[(tissue, {dataset: protein_count})]` in the input order.
    The PXD017199 entry is only present in the inner dict when both
    `cl_tissue_pxd017199` and `accessions_per_cell_line_pxd017199` are
    non-None — keeps the function back-compatible with the 3-dataset
    callers in older tests."""
    tissue_to_cls: dict[str, list[str]] = {}
    for cl, tissue in cl_tissue_pxd003539.items():
        tissue_to_cls.setdefault(tissue, []).append(cl)
    tissue_to_cls_17199: dict[str, list[str]] = {}
    have_17199 = cl_tissue_pxd017199 is not None and accessions_per_cell_line_pxd017199 is not None
    if have_17199:
        for cl, tissue in cl_tissue_pxd017199.items():
            tissue_to_cls_17199.setdefault(tissue, []).append(cl)
    tissue_to_cls_41421: dict[str, list[str]] = {}
    have_41421 = cl_tissue_pxd041421 is not None and accessions_per_cell_line_pxd041421 is not None
    if have_41421:
        for cl, tissue in cl_tissue_pxd041421.items():
            tissue_to_cls_41421.setdefault(tissue, []).append(cl)
    rows: list[tuple[str, dict[str, int]]] = []
    for tissue in tissue_order:
        cls = tissue_to_cls.get(tissue, [])
        u_3539: set[str] = set()
        for cl in cls:
            u_3539.update(accessions_per_cell_line_pxd003539.get(cl, set()))
        n_30304 = len(accessions_per_tissue_pxd030304.get(tissue, set()))
        n_4701 = len(accessions_pxd004701) if tissue == 'Breast' else 0
        by_ds = {'PXD003539': len(u_3539), 'PXD030304': n_30304, 'PXD004701': n_4701}
        if have_17199:
            u_17199: set[str] = set()
            for cl in tissue_to_cls_17199.get(tissue, []):
                u_17199.update(accessions_per_cell_line_pxd017199.get(cl, set()))
            by_ds['PXD017199'] = len(u_17199)
        if have_41421:
            u_41421: set[str] = set()
            for cl in tissue_to_cls_41421.get(tissue, []):
                u_41421.update(accessions_per_cell_line_pxd041421.get(cl, set()))
            by_ds['PXD041421'] = len(u_41421)
        rows.append((tissue, by_ds))
    return rows

def detection_count_histogram(accession_sets: dict[str, set[str]]) -> dict[int, int]:
    """For every accession in the union of the input sets, count how many
    sets it occurs in, then bucket sizes by that count (1, 2, ..., N).

    Returns `{count: bucket_size}` for `count` in `1..N` (N = number of
    input sets). Missing counts get a 0 entry so the histogram always
    has N bins."""
    n_sets = len(accession_sets)
    union: set[str] = set()
    for s in accession_sets.values():
        union.update(s)
    buckets: dict[int, int] = {k: 0 for k in range(1, n_sets + 1)}
    for acc in union:
        c = sum((1 for s in accession_sets.values() if acc in s))
        if c >= 1:
            buckets[c] = buckets.get(c, 0) + 1
    return buckets

def _annotate_panel_letter(ax, letter: str, *, subtitle: str | None=None) -> None:
    """Stamp a bold panel letter in the top-left corner. When
    `subtitle` is provided, render a smaller description next to the
    letter so each panel is self-naming (the bold letter + a short
    description of what the panel shows). Subtitles are paper-ready:
    short noun phrases, no terminal punctuation, fontsize 10.
    """
    ax.text(-0.07, 1.02, letter, transform=ax.transAxes, fontsize=14, fontweight='bold', ha='left', va='bottom')
    if subtitle:
        ax.text(-0.02, 1.02, subtitle, transform=ax.transAxes, fontsize=10, fontweight='normal', ha='left', va='bottom', color='#222222')

def _render_panel_a_headlines(ax, headlines: dict[str, DatasetHeadline]) -> None:
    """Grouped bar chart, N dataset groups × up-to-2 bars (paper vs
    quantmsdiann). Datasets whose `paper_count == 0` get only the
    quantmsdiann bar, centred on the x-tick — no explainer annotation;
    the cohort's x-axis label (PXDxxx + paper-year tag from
    `DATASET_LABELS`) is sufficient context. Bars are colour-coded paper
    (grey) vs quantmsdiann (blue); value labels sit just above each bar
    top with a small fractional pad so they never touch the bar."""
    datasets = list(headlines.keys())
    bar_width = 0.34
    x = list(range(len(datasets)))
    paper_drawn_xs: list[float] = []
    paper_drawn_vals: list[int] = []
    diann_xs: list[float] = []
    diann_vals_drawn: list[int] = []
    for xi, d in zip(x, datasets):
        h = headlines[d]
        if h.paper_count > 0:
            paper_drawn_xs.append(xi - bar_width / 2)
            paper_drawn_vals.append(h.paper_count)
            diann_xs.append(xi + bar_width / 2)
        else:
            diann_xs.append(xi)
        diann_vals_drawn.append(h.diann_count)
    bars_p = ax.bar(paper_drawn_xs, paper_drawn_vals, width=bar_width, color=fs.COMPARISON['original'], label='Original analysis')
    bars_d = ax.bar(diann_xs, diann_vals_drawn, width=bar_width, color=fs.COMPARISON['quantmsdiann'], label='quantmsdiann (DIA-NN)')
    all_vals = [h.diann_count for h in headlines.values()] + [h.paper_count for h in headlines.values() if h.paper_count > 0]
    ymax = (max(all_vals) if all_vals else 1) * 1.3
    pad = ymax * 0.02
    for bar, v in zip(bars_p, paper_drawn_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + pad, f'{v:,}', ha='center', va='bottom', fontsize=13, fontweight='bold')
    for bar, v, d in zip(bars_d, diann_vals_drawn, datasets):
        pc = headlines[d].paper_count
        label = f'{v:,}\n(+{round(100 * (v - pc) / pc)}%)' if pc > 0 else f'{v:,}'
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + pad, label, ha='center', va='bottom', fontsize=12, fontweight='bold', color=fs.COMPARISON['quantmsdiann'])
    ax.set_xticks(x)
    ax.set_xticklabels([DATASET_LABELS.get(d, d) for d in datasets], fontsize=12)
    ax.set_xlim(-1.2, len(datasets) + 0.2)
    ax.set_ylabel('Protein groups', fontsize=13)
    ax.set_ylim(0, ymax)
    ax.tick_params(axis='y', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper left', frameon=False, fontsize=12)

def _render_upset_in_axes(ax, sets: dict[str, set[str]], ds_order: list[str], *, bar_color: str=fs.OKABE_ITO['blue'], panel_title: str | None=None, panel_letter: str | None=None, panel_subtitle: str | None=None) -> None:
    """Render an UpSet plot inside an axes slot via a Matplotlib SubFigure.

    `upsetplot.UpSet.plot()` creates its own gridspec on the given
    figure (it does not accept a target Axes). To slot it into a
    pre-existing 4x2 grid cell we:
      1. Hide the host axes `ax`.
      2. Carve a SubFigure out of the host figure at the same
         SubplotSpec as `ax`.
      3. Pass that SubFigure to `upset.plot(fig=subfig)`. UpSet's
         internal gridspec is anchored to the subfig, not the parent
         figure, so it stays inside the panel cell.

    Falls back to a brief text annotation if fewer than 2 datasets have
    content (UpSet needs at least 2 categories)."""
    populated = [d for d in ds_order if sets.get(d)]
    if len(populated) < 2:
        ax.text(0.5, 0.5, f'Insufficient data for UpSet plot\n(populated datasets: {len(populated)})', ha='center', va='center', transform=ax.transAxes, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return
    try:
        from upsetplot import UpSet, from_contents
    except ImportError as exc:
        raise RuntimeError('upsetplot required for Panel B / D — install via `pip install upsetplot` (already in analysis/requirements.txt)') from exc
    fig = ax.figure
    subplotspec = ax.get_subplotspec()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor('none')
    ax.set_visible(False)
    contents = {d: sorted(sets[d]) for d in populated}
    data = from_contents(contents)
    subfig = fig.add_subfigure(subplotspec)
    upset = UpSet(data, subset_size='count', show_counts=True, sort_by='cardinality', sort_categories_by='cardinality', facecolor=bar_color, element_size=None)
    upset.plot(fig=subfig)
    if panel_title:
        subfig.suptitle(panel_title, fontsize=10)
    if panel_letter:
        subfig.text(0.005, 0.985, panel_letter, fontsize=14, fontweight='bold', ha='left', va='top')
        if panel_subtitle:
            subfig.text(0.045, 0.985, panel_subtitle, fontsize=10, ha='left', va='top', color='#222222')

def _render_panel_b_cellline_venn(ax, sets: dict[str, set[str]], region_sizes: dict[str, int] | None=None) -> None:
    """UpSet plot of normalised cell-line names across the 4 datasets.

    Replaces the original 3-set venn3 rendering; with PXD017199 added,
    the 4-way cell-line overlap (especially PXD017199-vs-PXD004701
    breast overlap) is unreadable in a 4-set Venn. UpSet renders one
    bar per intersection region, ordered by cardinality."""
    ds_order = [d for d in ('PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421') if d in sets]
    _render_upset_in_axes(ax, sets, ds_order, bar_color=fs.OKABE_ITO['blue'], panel_title='Cell-line set intersections (UpSet)')

def _render_panel_c_tissue_coverage(ax, rows: list[tuple[str, dict[str, int]]]) -> None:
    """Stacked horizontal bars on a unified tissue axis. Each bar segment
    is a dataset's per-tissue cell-line count; tissues sorted by descending
    total across datasets. Tissues with zero contribution from all three
    datasets never appear (filtered upstream)."""
    if not rows:
        ax.text(0.5, 0.5, 'no tissues', ha='center', va='center', transform=ax.transAxes)
        return
    tissues = [r[0] for r in rows]
    ds_order = ['PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421']
    ds_order = [ds for ds in ds_order if any((r[1].get(ds, 0) > 0 for r in rows))]
    n = len(tissues)
    y = list(range(n))
    y_top = list(reversed(y))
    left = [0] * n
    for ds in ds_order:
        vals = [r[1].get(ds, 0) for r in rows]
        ax.barh(y_top, vals, left=left, color=DATASET_COLORS[ds], label=DATASET_LABELS[ds].replace('\n', ' '), edgecolor='white', linewidth=0.4)
        left = [l + v for l, v in zip(left, vals)]
    totals = [sum(r[1].values()) for r in rows]
    for yi, t in zip(y_top, totals):
        ax.text(t + max(totals) * 0.005, yi, f'{t:,}', ha='left', va='center', fontsize=11)
    ax.set_yticks(y_top)
    ax.set_yticklabels(tissues, fontsize=11)
    ax.tick_params(axis='x', labelsize=11)
    ax.set_xlabel('Cell lines (sum across datasets)', fontsize=12)
    ax.set_xlim(0, max(totals) * 1.1)
    fs.despine(ax)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=min(3, len(ds_order)), frameon=False, fontsize=10, columnspacing=1.2, handletextpad=0.5)

def _render_panel_d_protein_venn(ax, sets: dict[str, set[str]], *, panel_letter: str | None=None, panel_subtitle: str | None=None) -> None:
    """UpSet plot of UniProt accessions across the 4 quantmsdiann
    analyses. Replaces the original 3-set venn3 rendering for the same
    reason as Panel B — the 4-set Venn is unreadable.

    `panel_letter` / `panel_subtitle` flow through to the SubFigure
    annotation (since the UpSet renderer hides the host axes, normal
    `_annotate_panel_letter(ax, ...)` calls go to an invisible Axes
    and never appear on the rendered SVG)."""
    ds_order = [d for d in ('PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421') if d in sets]
    _render_upset_in_axes(ax, sets, ds_order, bar_color=fs.OKABE_ITO['blue'], panel_letter=panel_letter, panel_subtitle=panel_subtitle)

def _render_panel_e_breadth_vs_depth(ax, cellline_sets: dict[str, set[str]], accession_sets: dict[str, set[str]], *, runs_per_cohort: dict[str, int] | None=None) -> None:
    """Per-cohort breadth-vs-depth scatter. One dot per cohort:
      - x = number of distinct cell lines in the cohort SDRF
        (`len(cellline_sets[ds])`), log scale (range 2 → 947)
      - y = target-only union of UniProt accessions
        (`len(accession_sets[ds])`)
      - dot size proportional to the number of MS runs in the cohort
        (per `runs_per_cohort`, if provided)
      - colour = `DATASET_COLORS[ds]`

    Replaces the previous rarefaction curves whose x-axes were
    incompatible across cohorts (cell-lines / tissues / subtypes).
    Single-glance view of where each cohort sits on the
    breadth↔depth tradeoff: PXD030304 is broad (947 lines) but
    matched to a deep proteome; PXD041421 is narrow (2 lines) but
    deep (24 reps each); PXD003539 / PXD017199 / PXD004701 sit in
    the middle. Cohorts missing from either map are silently
    omitted."""
    ds_order = ['PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421']
    any_plotted = False
    xs_all: list[float] = []
    ys_all: list[float] = []
    for ds in ds_order:
        n_cl = len(cellline_sets.get(ds, set()))
        n_acc = len(accession_sets.get(ds, set()))
        if n_cl == 0 or n_acc == 0:
            continue
        n_runs = (runs_per_cohort or {}).get(ds, 0) or 0
        if n_runs > 0:
            size = 60 + 240 * n_runs ** 0.5 / 6000 ** 0.5
        else:
            size = 120
        ax.scatter([n_cl], [n_acc], s=size, c=DATASET_COLORS[ds], edgecolors='#222222', linewidths=0.6, zorder=3, label=f'{ds} ({n_runs:,} runs)' if n_runs > 0 else ds)
        ax.annotate(DATASET_LABELS[ds].replace('\n', ' '), xy=(n_cl, n_acc), xytext=(6, 6), textcoords='offset points', fontsize=8, color='#222222')
        xs_all.append(n_cl)
        ys_all.append(n_acc)
        any_plotted = True
    if not any_plotted:
        ax.text(0.5, 0.5, 'no cohort inputs', ha='center', va='center', transform=ax.transAxes)
        return
    ax.set_xscale('log')
    ax.set_xlabel('Cell lines per cohort', fontsize=10)
    ax.set_ylabel('Proteins (union per cohort)', fontsize=10)
    if xs_all and ys_all:
        ax.set_xlim(min(xs_all) * 0.4, max(xs_all) * 3.0)
        ax.set_ylim(0, max(ys_all) * 1.15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=8)
    ax.legend(loc='lower right', frameon=False, fontsize=10)

def _render_panel_f_tissue_protein_counts(ax, rows: list[tuple[str, set[str]]] | list[tuple[str, dict[str, int]]]) -> None:
    """Per-tissue UNIQUE-protein bars. Single bar per tissue showing
    the union of UniProt accessions across all contributing cohorts —
    the biologically meaningful number (`how many distinct proteins did
    the whole atlas observe in this tissue?`).

    Replaces the previous stacked-bar that summed per-cohort counts
    and so double-counted proteins seen in multiple cohorts. The
    cohort-breakdown is still preserved in the audit TSV
    (`atlas_distribution | Panel B | tissue | <cohort> = <count>`).

    `rows` accepts two shapes for backward compatibility:
      - `[(tissue, set[str])]` — the new union-set shape (preferred).
      - `[(tissue, {dataset: count})]` — the legacy per-cohort dict
        (which we sum, replicating the old visual). Callers should
        migrate to the set-of-accessions form.
    Same tissue ordering as Panel A (per-tissue cell lines) — caller
    passes the already-ordered list."""
    if not rows:
        ax.text(0.5, 0.5, 'no tissue/protein rows', ha='center', va='center', transform=ax.transAxes)
        return
    tissues = [r[0] for r in rows]
    sample_value = rows[0][1]
    if isinstance(sample_value, set):
        counts = [len(r[1]) for r in rows]
        x_label = 'Unique proteins per tissue'
    else:
        counts = [sum(r[1].values()) for r in rows]
        x_label = 'Proteins per tissue (legacy sum)'
    n = len(tissues)
    y_top = list(reversed(range(n)))
    max_total = max(counts) if counts else 1
    bar_color = '#7570b3'
    ax.barh(y_top, counts, color=bar_color, edgecolor='white', linewidth=0.4)
    for yi, t in zip(y_top, counts):
        ax.text(t + max_total * 0.005, yi, f'{t:,}', ha='left', va='center', fontsize=11)
    ax.set_yticks(y_top)
    ax.set_yticklabels(tissues, fontsize=11)
    ax.tick_params(axis='x', labelsize=11)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_xlim(0, max_total * 1.12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def _render_panel_g_detection_histogram(ax, accession_sets: dict[str, set[str]]) -> None:
    """3-bar (or N-bar) chart of how many accessions are seen in 1, 2 or
    3 datasets. Bars annotated with absolute count and percentage of the
    union. Title-style headline above the bars."""
    buckets = detection_count_histogram(accession_sets)
    n_sets = len(accession_sets)
    if n_sets == 0 or not buckets:
        ax.text(0.5, 0.5, 'no detection inputs', ha='center', va='center', transform=ax.transAxes)
        return
    xs = sorted(buckets.keys())
    vals = [buckets[k] for k in xs]
    total = sum(vals)
    palette = ['#9e9e9e', '#1f77b4', '#1b7a3a', '#d95f02', '#8e44ad']
    bars = ax.bar(xs, vals, color=palette[:len(xs)], edgecolor='black', linewidth=0.4)
    ymax = max(vals) * 1.22 if max(vals) > 0 else 1
    pad = ymax * 0.02
    for bar, v in zip(bars, vals):
        pct = 100.0 * v / total if total > 0 else 0.0
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + pad, f'{v:,}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=7)
    ax.set_xticks(xs)
    ax.set_xticklabels([f'{k}' for k in xs])
    ax.set_xlabel('Datasets a protein is detected in', fontsize=12)
    ax.set_ylabel('UniProt accessions', fontsize=12)
    ax.set_title(f'Pan-cohort core ({n_sets} datasets)', fontsize=10)
    ax.set_ylim(0, ymax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def _render_panel_h_expression_atlas_overlap(ax, ea_genes: set[str], diann_genes: set[str]) -> None:
    """Expression Atlas (E-PROT-73 NCI-60) vs PXD003539 quantmsdiann gene
    overlap. Renders 3 bars: EA catalogue, quantmsdiann, intersection;
    annotates intersection coverage as a percentage of EA.

    If `ea_genes` is empty (file missing or empty), renders an explainer
    panel and returns without drawing bars."""
    if not ea_genes:
        ax.text(0.5, 0.5, 'Expression Atlas catalogue (E-PROT-73) unavailable\nPlace data/E-PROT-73-query-results.tsv to populate this panel.', ha='center', va='center', transform=ax.transAxes, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return
    inter = ea_genes & diann_genes
    labels = ['Expression Atlas\n(NCI-60 catalogue)', 'quantmsdiann\n(PXD003539)', 'Intersection']
    vals = [len(ea_genes), len(diann_genes), len(inter)]
    colors = ['#9e9e9e', DATASET_COLORS['PXD003539'], '#1f77b4']
    xs = list(range(3))
    bars = ax.bar(xs, vals, color=colors, edgecolor='black', linewidth=0.4)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{v:,}', ha='center', va='bottom', fontsize=13, fontweight='bold')
    coverage_pct = 100.0 * len(inter) / len(ea_genes) if ea_genes else 0.0
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Gene symbols', fontsize=12)
    ax.set_title(f'PXD003539 covers {coverage_pct:.1f}% of E-PROT-73 gene catalogue', fontsize=10)
    ymax = max(vals) * 1.18 if max(vals) > 0 else 1
    ax.set_ylim(0, ymax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def render_atlas_overlap(headlines: dict[str, DatasetHeadline], cellline_sets: dict[str, set[str]], accession_sets: dict[str, set[str]], svg_path: Path) -> None:
    """Compose the cohort-overlap atlas figure: 3 panels relabelled
    A/B/C after the 2026-05-21 split:

      Row 1 (full width): A — protein-accession UpSet (5 cohorts)
      Row 2 left:         B — per-cohort headline counts
      Row 2 right:        C — pan-cohort detection histogram

    Panel A (the UpSet) gets the wide top row because the 5-set UpSet
    needs horizontal space for both the matrix and the intersection-bar
    chart. The two summary panels (B, C) sit alongside each other below
    where their narrower aspect ratio works.

    `cellline_sets` is retained in the signature for backwards
    compatibility — Panel B-as-cell-line-UpSet was removed in the
    2026-05-21 cleanup because PXD030304's 947 lines dominated the
    inter-cohort intersections."""
    del cellline_sets
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.4, 1.0], hspace=0.3, wspace=0.3)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    _render_panel_d_protein_venn(ax_a, accession_sets, panel_letter='A', panel_subtitle='Protein overlap')
    _render_panel_a_headlines(ax_b, headlines)
    _annotate_panel_letter(ax_b, 'B', subtitle='Headline counts')
    _render_panel_g_detection_histogram(ax_c, accession_sets)
    _annotate_panel_letter(ax_c, 'C', subtitle='Detection counts')
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def render_atlas_distribution(tissue_rows: list[tuple[str, dict[str, int]]], cellline_sets: dict[str, set[str]], accession_sets: dict[str, set[str]], svg_path: Path, *, tissue_protein_rows: list[tuple[str, dict[str, int]]] | None=None, runs_per_cohort: dict[str, int] | None=None) -> None:
    """Compose the per-tissue distribution atlas figure: two stacked
    per-tissue bar panels in a 2x1 grid (figsize 14x14). Layout:

        Row 1: A (wide, per-tissue cell-line stacked bar)
        Row 2: B (wide, per-tissue unique-protein stacked bar)

    Panel H (Expression Atlas vs PXD003539 gene overlap) was removed
    from the atlas because it duplicates `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg`
    — the same NCI-60 gene-set comparison already lives in the
    PXD003539 per-cohort figure where it makes more sense.

    The breadth-vs-depth scatter (former Panel C) was dropped 2026-05-29:
    it crowded the figure and the breadth↔depth tradeoff is already
    legible from the two per-tissue stacked bars (PXD030304's broad,
    shallow per-tissue spread vs the narrow, deep single-line cohorts).
    `cellline_sets`, `accession_sets`, and `runs_per_cohort` are retained
    in the signature for API stability with `render_atlas` but are no
    longer used here.

    All extended inputs default to empty containers so the renderer is
    resilient on partial-data runs."""
    del cellline_sets, accession_sets, runs_per_cohort
    tissue_protein_rows = tissue_protein_rows or []
    fig = plt.figure(figsize=(14, 14))
    gs = fig.add_gridspec(2, 1, hspace=0.4)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    _render_panel_c_tissue_coverage(ax_a, tissue_rows)
    _annotate_panel_letter(ax_a, 'A', subtitle='Cell lines per tissue')
    _render_panel_f_tissue_protein_counts(ax_b, tissue_protein_rows)
    _annotate_panel_letter(ax_b, 'B', subtitle='Proteins per tissue')
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def _render_tissue_combined(ax, tissue_rows: list[tuple[str, dict[str, int]]], protein_rows: list[tuple[str, set[str]]] | list[tuple[str, dict[str, int]]]) -> None:
    """Combined per-tissue panel merging old panels b and c onto one shared
    tissue axis, drawn back-to-back (population-pyramid style):
      - LEFT  of centre -> cell lines, stacked by contributing cohort
      - RIGHT of centre -> unique UniProt proteins, single colour
    The two metrics differ ~100x in scale, so each side is normalised to its
    own maximum (bar *length* is a fraction of that metric's max); real counts
    are annotated at every bar end and tissue labels appear once on the left."""
    if not tissue_rows:
        ax.text(0.5, 0.5, 'no tissue rows', ha='center', va='center', transform=ax.transAxes)
        return
    from matplotlib.patches import Patch
    ds_order = ['PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421']
    tissues = [t for t, _ in tissue_rows]
    prot_map = {t: len(v) if isinstance(v, set) else sum(v.values()) for t, v in protein_rows}
    n = len(tissues)
    y = list(reversed(range(n)))
    bh = 0.62
    cell_tot = [sum(d.values()) for _, d in tissue_rows]
    max_c = max(cell_tot) if cell_tot else 1
    max_p = max((prot_map.get(t, 0) for t in tissues)) or 1
    prot_color = '#34495e'
    drawn: set[str] = set()
    for yi, (t, d) in zip(y, tissue_rows):
        left = 0.0
        for ds in ds_order:
            c = d.get(ds, 0)
            if not c:
                continue
            w = c / max_c
            ax.barh(yi, -w, left=-left, height=bh, color=DATASET_COLORS[ds], edgecolor='white', linewidth=0.4, label=DATASET_LABELS[ds].replace('\n', ' ') if ds not in drawn else None)
            drawn.add(ds)
            left += w
        ax.text(-left - 0.015, yi, f'{int(round(left * max_c)):,}', ha='right', va='center', fontsize=9.5)
    for yi, t in zip(y, tissues):
        p = prot_map.get(t, 0)
        w = p / max_p
        ax.barh(yi, w, height=bh, color=prot_color, edgecolor='white', linewidth=0.4)
        ax.text(w + 0.015, yi, f'{p:,}', ha='left', va='center', fontsize=9.5, color=prot_color)
    ax.axvline(0, color='#888888', linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(tissues, fontsize=11)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(-1.32, 1.3)
    ax.set_xticks([])
    for sp in ('top', 'right', 'bottom'):
        ax.spines[sp].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.text(-0.5, n - 0.15, f'← Cell lines per tissue (max {max_c:,})', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.text(0.5, n - 0.15, f'Unique proteins per tissue (max {max_p:,}) →', ha='center', va='bottom', fontsize=12, fontweight='bold', color=prot_color)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(facecolor=prot_color, edgecolor='white'))
    labels.append('Unique proteins')
    ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.04), frameon=False, fontsize=10, ncol=len(labels))

def render_atlas_main(headlines: dict[str, DatasetHeadline], tissue_rows: list[tuple[str, dict[str, int]]], svg_path: Path, *, tissue_protein_rows: list[tuple[str, dict[str, int]]] | None=None) -> None:
    """Compose the main pan-cohort figure as three stacked, full-width
    panels (drops the protein-accession UpSet overlap and the
    detection-count histogram, which were too small to read in the
    multi-panel layout):

        A (top):    per-cohort headline counts (paper vs quantmsdiann)
        B (middle): per-tissue cell-line coverage (stacked bars)
        C (bottom): per-tissue unique-protein coverage (stacked bars)

    Removing the UpSet panel also removes the upsetplot/numpy-2 render
    dependency, so this figure renders on any supported NumPy."""
    tissue_protein_rows = tissue_protein_rows or []
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.5, 1.0], hspace=0.3)
    ax_a = fig.add_subplot(gs[0])
    ax_bc = fig.add_subplot(gs[1])
    _render_panel_a_headlines(ax_a, headlines)
    _annotate_panel_letter(ax_a, 'a', subtitle='Protein groups (≥2 peptides): deposited vs quantms.io')
    _render_tissue_combined(ax_bc, tissue_rows, tissue_protein_rows)
    _annotate_panel_letter(ax_bc, 'b', subtitle='Per-tissue cell lines and unique proteins')
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def render_atlas(headlines: dict[str, DatasetHeadline], cellline_sets: dict[str, set[str]], tissue_rows: list[tuple[str, dict[str, int]]], accession_sets: dict[str, set[str]], svg_path: Path, *, tissue_protein_rows: list[tuple[str, dict[str, int]]] | None=None, runs_per_cohort: dict[str, int] | None=None) -> None:
    """Thin wrapper that renders both the overlap (`atlas_overlap.svg`,
    panels A/B/D/G) and the distribution (`atlas_distribution.svg`,
    panels C/E/F) atlas figures.

    The single `svg_path` argument is interpreted as a hint: the actual
    outputs are written next to it as `<stem>_overlap.svg` and
    `<stem>_distribution.svg`. Panel H (Expression Atlas overlap) was
    removed from the atlas — it lives in
    `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg`."""
    out_dir = svg_path.parent
    stem = svg_path.stem
    overlap_path = out_dir / f'{stem}_overlap.svg'
    distribution_path = out_dir / f'{stem}_distribution.svg'
    render_atlas_overlap(headlines, cellline_sets, accession_sets, overlap_path)
    render_atlas_distribution(tissue_rows, cellline_sets, accession_sets, distribution_path, tissue_protein_rows=tissue_protein_rows, runs_per_cohort=runs_per_cohort)

def _set_region_sizes(sets: dict[str, set[str]], ds_order: list[str]) -> dict[str, int]:
    """Return all 2^n - 1 non-empty-membership-pattern region sizes
    across n sets in `ds_order`.

    Region-key shapes (matched by the counts.tsv writer and the tests):
      - Singletons: `f"{ds},only"`
      - Pairs: `f"{ds_a},{ds_b}"` (alphabetical pair within `ds_order`)
      - Triples (n>=3): `f"{ds_a},{ds_b},{ds_c}"` for the unique triple,
        OR the special key `"all_three"` when n == 3 (kept for backwards
        compatibility with the 3-set counts.tsv schema).
      - n-tuples for n == 4: the special key `"all_four"`.
      - General m-tuples (3 <= m < n): comma-separated list of dataset
        names in `ds_order`.

    Sizes are the count of accessions that belong to exactly that
    subset of datasets (a partition of the union). The 3-set helper
    `_venn_region_sizes_3` and the test fixture call this function and
    map the keys back through `_legacy_3_set_region_keys` below."""
    from itertools import combinations
    n = len(ds_order)
    if n == 0:
        return {}
    union: set[str] = set()
    for d in ds_order:
        union |= sets[d]
    out: dict[str, int] = {}
    for r in range(1, n + 1):
        for combo in combinations(ds_order, r):
            members = set(combo)
            non_members = [d for d in ds_order if d not in members]
            region = set(union)
            for d in combo:
                region &= sets[d]
            for d in non_members:
                region -= sets[d]
            if r == 1:
                key = f'{combo[0]},only'
            elif r == n and n == 3:
                key = 'all_three'
            elif r == n and n == 4:
                key = 'all_four'
            elif r == n and n == 5:
                key = 'all_five'
            else:
                key = ','.join(combo)
            out[key] = len(region)
    return out

def _venn_region_sizes_3(sets: dict[str, set[str]], ds_order: list[str]) -> dict[str, int]:
    """Backwards-compatible 3-set wrapper around `_set_region_sizes`.

    Returns the 7 Venn region sizes using the legacy key shape:
    `{ds}_only`, `{ds_a}+{ds_b}` for every pair, `all_three`.

    Used only by the counts.tsv writer's 3-set rows and the test fixture
    `test_venn_region_sizes_3_partitions_correctly`. New 4-set callers
    should use `_set_region_sizes` directly."""
    if len(ds_order) != 3:
        raise ValueError(f'_venn_region_sizes_3 expects 3 datasets, got {len(ds_order)}')
    generic = _set_region_sizes(sets, ds_order)
    a, b, c = ds_order
    return {f'{a}_only': generic[f'{a},only'], f'{b}_only': generic[f'{b},only'], f'{c}_only': generic[f'{c},only'], f'{a}+{b}': generic[f'{a},{b}'], f'{a}+{c}': generic[f'{a},{c}'], f'{b}+{c}': generic[f'{b},{c}'], 'all_three': generic['all_three']}

def write_combined_counts_tsv(tsv_path: Path, headlines: dict[str, DatasetHeadline], cellline_sets: dict[str, set[str]], tissue_rows: list[tuple[str, dict[str, int]]], accession_sets: dict[str, set[str]], *, tissue_protein_rows: list[tuple[str, set[str]]] | list[tuple[str, dict[str, int]]] | None=None, tissue_protein_rows_per_cohort: list[tuple[str, dict[str, int]]] | None=None, runs_per_cohort: dict[str, int] | None=None) -> None:
    """Auditable TSV with Panel-feeding numbers (Panel A bars, Panel B
    set-intersection regions, Panel C per-(tissue, dataset) counts,
    Panel D set-intersection regions, Panels E/F/G rows). Panel H
    (Expression Atlas overlap) was removed because it duplicates
    `analysis/figures/PXD003539/supp_walzer_vs_quantms_genes_ensembl.svg`.
    The dataset ordering tracks the keys of `cellline_sets` so callers
    can pass 3 or 4 datasets transparently."""
    ds_order = [d for d in ('PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421') if d in cellline_sets]
    rows: list[tuple[str, str, int, str]] = []
    for ds in ds_order:
        if ds not in headlines:
            continue
        h = headlines[ds]
        if h.paper_count > 0:
            rows.append((f'atlas_overlap | Panel B | {ds} | original', h.paper_label, h.paper_count, h.metric))
        else:
            rows.append((f'atlas_overlap | Panel B | {ds} | original', 'no paper DIA headline available', 0, h.metric))
        rows.append((f'atlas_overlap | Panel B | {ds} | quantmsdiann', 'quantmsdiann (DIA-NN, methods.md §1 global rule, no filters)', h.diann_count, h.metric))
    for ds in ds_order:
        rows.append((f'atlas_overlap | dropped (cell-line UpSet) | {ds} | count', 'normalised cell-line names from SDRF (normalize_cell_line)', len(cellline_sets.get(ds, set())), 'cell-line UpSet removed from atlas_overlap; raw count retained for audit only'))
    for tissue, by_ds in tissue_rows:
        for ds in ds_order:
            rows.append((f'atlas_distribution | Panel A | tissue cell lines | {tissue}', ds, by_ds.get(ds, 0), 'cell lines mapped to this tissue (unified ProCan 28-tissue axis)'))
    acc_ds_order = [d for d in ds_order if d in accession_sets]
    acc_regions = _set_region_sizes(accession_sets, acc_ds_order)
    for region, n in acc_regions.items():
        rows.append((f'atlas_overlap | Panel A | accession intersections | {region}', 'UniProt accessions extracted from Protein.Group (extract_accessions_diann)', n, 'PXD003539/PXD017199 from pr_matrix; PXD030304/PXD004701 from cached per-tissue/per-subtype JSON'))
    runs = runs_per_cohort or {}
    for ds in ds_order:
        n_cl = len(cellline_sets.get(ds, set()))
        n_acc = len(accession_sets.get(ds, set()))
        n_runs = runs.get(ds, 0)
        if n_cl == 0 or n_acc == 0:
            continue
        rows.append((f'atlas_distribution | off-figure | breadth-vs-depth | {ds}', f'cell_lines={n_cl} runs={n_runs} accessions={n_acc}', n_acc, f'{ds}: {n_cl} cell line(s) × {n_runs} MS runs → {n_acc:,} UniProt accessions'))
    for tissue, value in tissue_protein_rows or []:
        union_size = len(value) if isinstance(value, set) else sum(value.values())
        rows.append((f'atlas_distribution | Panel B | tissue proteins | {tissue}', 'union of UniProt accessions across all contributing cohorts', union_size, 'unique proteins detected in this tissue across the atlas (accession-overlap geometry)'))
    for tissue, by_ds in tissue_protein_rows_per_cohort or []:
        for ds in ds_order:
            rows.append((f'atlas_distribution | Panel B | tissue proteins (per-cohort breakdown) | {tissue}', ds, by_ds.get(ds, 0), 'per-cohort protein count for this tissue (sum across cohorts double-counts proteins seen in multiple cohorts — figure shows union)'))
    g_buckets = detection_count_histogram(accession_sets)
    g_total = sum(g_buckets.values())
    n_datasets = len(accession_sets)
    for k, n in sorted(g_buckets.items()):
        pct = 100.0 * n / g_total if g_total > 0 else 0.0
        rows.append((f'atlas_overlap | Panel C | detections-in-n-datasets | n={k}', f'UniProt accessions across the {n_datasets} quantmsdiann analyses', n, f'{pct:.2f}% of pan-cohort union ({g_total:,})'))
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, 'w', encoding='utf-8') as fh:
        fh.write('metric\tsource\tcount\tnote\n')
        for r in rows:
            fh.write('\t'.join((str(x) for x in r)) + '\n')
PREREQS: list[tuple[Path, str]] = [(PXD003539_SDRF, 'python -m analysis.figure_original_vs_quantmsdiann'), (PXD003539_PR_MATRIX, 'python -m analysis.figure_original_vs_quantmsdiann'), (PXD030304_SDRF, 'python -m analysis.figure_pxd030304_procan_vs_quantmsdiann'), (PXD030304_TISSUE_MAPPING, 'python -m analysis.figure_pxd030304_procan_vs_quantmsdiann'), (PXD030304_PROTEIN_JSON, 'python -m analysis.figure_pxd030304_procan_vs_quantmsdiann'), (PXD004701_SDRF, 'python -m analysis.figure_pxd004701_sun_vs_quantmsdiann'), (PXD004701_PROTEIN_JSON, 'python -m analysis.figure_pxd004701_sun_vs_quantmsdiann'), (PXD017199_SDRF, 'stage data/PXD017199/PXD017199.sdrf.tsv (Tognetti 2021)'), (PXD017199_PR_MATRIX, 'stage data/PXD017199/diann_report.pr_matrix.tsv (Tognetti 2021)'), (PXD041421_SDRF, 'stage data/PXD041421/PXD041421.sdrf.tsv (Wang 2023)'), (PXD041421_PR_MATRIX, 'stage data/PXD041421/diann_report.pr_matrix.tsv (Wang 2023)')]

def check_prerequisites() -> list[tuple[Path, str]]:
    """Return the list of `(missing_path, instruction)` tuples; empty list
    means all prerequisites are present."""
    return [(p, cmd) for p, cmd in PREREQS if not (p.exists() and p.stat().st_size > 0)]

def figure_combined_cell_lines_atlas_main() -> int:
    # atlas reads accessions from pr_matrix for PXD003539/017199/041421 only;
    # ProCan/Sun come from cached JSONs, so skip their (large) matrices.
    # PXD003539 with its report parquet so refresh_dataset_headlines() actually
    # recomputes the headline from the FTP instead of falling back to the
    # precomputed constant; 017199/041421 need only their pr/pg matrices.
    ensure_cell_line_matrices('PXD003539', with_report=True)
    ensure_cell_line_matrices('PXD017199', 'PXD041421')
    missing = check_prerequisites()
    if missing:
        print('Combined atlas requires the per-dataset scripts to run first.', file=sys.stderr)
        for path, cmd in missing:
            rel = path.relative_to(_figure_combined_cell_lines_atlas__REPO_ROOT)
            print(f'  Missing: {rel}', file=sys.stderr)
            print(f'  Produce it via: {cmd}', file=sys.stderr)
        return 1
    _figure_combined_cell_lines_atlas__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    print('Refreshing DATASET_HEADLINES diann_count from on-disk inputs (methods.md §1 global rule: matrix rows / Lib.PG.Q.Value, no filters)...')
    refresh_dataset_headlines()
    for ds, h in DATASET_HEADLINES.items():
        print(f'  {ds}: protein_groups={h.diann_count:,}')
    print('Loading cell-line sets from SDRFs...')
    cellline_sets: dict[str, set[str]] = {'PXD003539': cell_lines_from_sdrf(PXD003539_SDRF), 'PXD030304': cell_lines_from_sdrf(PXD030304_SDRF), 'PXD004701': cell_lines_from_sdrf(PXD004701_SDRF), 'PXD017199': cell_lines_from_sdrf(PXD017199_SDRF), 'PXD041421': cell_lines_from_sdrf(PXD041421_SDRF)}
    for d, s in cellline_sets.items():
        print(f'  {d}: {len(s):,} cell lines')
    print('Building per-dataset cell-line -> tissue maps...')
    cl_tissue_3539 = cell_line_tissue_pxd003539(PXD003539_SDRF)
    cl_tissue_30304 = cell_line_tissue_pxd030304(PXD030304_SDRF, PXD030304_TISSUE_MAPPING)
    cl_tissue_4701 = cell_line_tissue_pxd004701(PXD004701_SDRF)
    cl_tissue_17199 = cell_line_tissue_pxd017199(PXD017199_SDRF)
    cl_tissue_41421 = cell_line_tissue_pxd041421(PXD041421_SDRF)
    tissue_rows = combined_tissue_table({'PXD003539': cl_tissue_3539, 'PXD030304': cl_tissue_30304, 'PXD004701': cl_tissue_4701, 'PXD017199': cl_tissue_17199, 'PXD041421': cl_tissue_41421})
    print(f'  unified-axis tissues: {len(tissue_rows)}')
    print('  Top-5 tissues by combined cell-line count:')
    for tissue, by_ds in tissue_rows[:5]:
        total = sum(by_ds.values())
        seg = ', '.join((f'{ds}:{by_ds.get(ds, 0)}' for ds in ('PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421')))
        print(f'    {tissue:<32s} total={total:>4d}  ({seg})')
    print('Loading per-group accession sets (for Panels E and F)...')
    per_cl_3539 = pxd003539_accessions_per_cell_line(PXD003539_PR_MATRIX, PXD003539_SDRF)
    per_tissue_30304 = pxd030304_accessions_per_tissue(PXD030304_PROTEIN_JSON)
    per_subtype_4701 = pxd004701_accessions_per_subtype(PXD004701_PROTEIN_JSON)
    per_cl_17199 = pxd017199_accessions_per_cell_line(PXD017199_PR_MATRIX, PXD017199_SDRF)
    per_cl_41421 = pxd041421_accessions_per_cell_line(PXD041421_PR_MATRIX, PXD041421_SDRF)
    print(f'  PXD003539 per-cell-line groups: {len(per_cl_3539)}')
    print(f'  PXD030304 per-tissue groups: {len(per_tissue_30304)}')
    print(f'  PXD004701 per-subtype groups: {len(per_subtype_4701)}')
    print(f'  PXD017199 per-cell-line groups: {len(per_cl_17199)}')
    print(f'  PXD041421 per-cell-line groups: {len(per_cl_41421)}')
    print('Loading protein-accession sets...')
    acc_3539: set[str] = set()
    for s in per_cl_3539.values():
        acc_3539.update(s)
    if not acc_3539:
        acc_3539 = pxd003539_protein_accessions(PXD003539_PR_MATRIX)
    acc_30304 = pxd030304_protein_accessions(PXD030304_PROTEIN_JSON)
    acc_4701 = pxd004701_protein_accessions(PXD004701_PROTEIN_JSON)
    acc_17199: set[str] = set()
    for s in per_cl_17199.values():
        acc_17199.update(s)
    if not acc_17199:
        acc_17199 = pxd017199_protein_accessions(PXD017199_PR_MATRIX)
    acc_41421: set[str] = set()
    for s in per_cl_41421.values():
        acc_41421.update(s)
    if not acc_41421:
        acc_41421 = pxd041421_protein_accessions(PXD041421_PR_MATRIX)
    accession_sets = {'PXD003539': acc_3539, 'PXD030304': acc_30304, 'PXD004701': acc_4701, 'PXD017199': acc_17199, 'PXD041421': acc_41421}
    for d, s in accession_sets.items():
        print(f'  {d}: {len(s):,} accessions')
    print('Building per-tissue UNIQUE-accession rows (atlas_distribution Panel B)...')
    tissue_order = [t for t, _ in tissue_rows]
    tissue_protein_rows = per_tissue_union_accessions(tissue_order, cl_tissue_3539, per_cl_3539, per_tissue_30304, acc_4701, cl_tissue_pxd017199=cl_tissue_17199, accessions_per_cell_line_pxd017199=per_cl_17199, cl_tissue_pxd041421=cl_tissue_41421, accessions_per_cell_line_pxd041421=per_cl_41421)
    tissue_protein_rows_per_cohort = per_tissue_protein_counts(tissue_order, cl_tissue_3539, per_cl_3539, per_tissue_30304, acc_4701, cl_tissue_pxd017199=cl_tissue_17199, accessions_per_cell_line_pxd017199=per_cl_17199, cl_tissue_pxd041421=cl_tissue_41421, accessions_per_cell_line_pxd041421=per_cl_41421)
    print('Computing per-cohort MS-run counts (Panel E)...')
    runs_per_cohort = _compute_runs_per_cohort()
    for ds in ('PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421'):
        print(f'  {ds}: {runs_per_cohort.get(ds, 0):,} runs')
    print('Rendering main pan-cohort figure (headline + per-tissue panels)...')
    main_svg = _figure_combined_cell_lines_atlas__FIGURES_DIR / 'atlas_main.svg'
    _ri = pd.read_csv(_figure_combined_cell_lines_atlas__REPO_ROOT / 'analysis' / 'figures' / 'reanalysis' / 'data' / 'reanalysis_improvement.tsv', sep='\t').set_index('dataset')
    _panel_a_labels = {'PXD003539': 'Guo 2019 (OpenSWATH)', 'PXD030304': 'ProCan 2022'}
    panel_a_headlines = {acc: DatasetHeadline(int(_ri.loc[acc, 'original']), int(_ri.loc[acc, 'reanalysis']), _panel_a_labels[acc], 'Protein groups (>=2 peptides)') for acc in ('PXD003539', 'PXD030304')}
    TOP_TISSUES = 12
    tissue_rows_top = tissue_rows[:TOP_TISSUES]
    _top_set = {t for t, _ in tissue_rows_top}
    tissue_protein_rows_top = [(t, d) for t, d in tissue_protein_rows if t in _top_set]
    render_atlas_main(panel_a_headlines, tissue_rows_top, main_svg, tissue_protein_rows=tissue_protein_rows_top)
    print(f'  saved: {main_svg}')
    import shutil
    manuscript_svg = _figure_combined_cell_lines_atlas__REPO_ROOT / 'analysis' / 'figures' / 'manuscript' / 'fig4_cellline_atlas.svg'
    manuscript_svg.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(main_svg, manuscript_svg)
    print(f'  saved: {manuscript_svg}')
    print('Writing combined counts.tsv...')
    data_dir = _figure_combined_cell_lines_atlas__FIGURES_DIR / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    tsv = data_dir / 'combined_counts.tsv'
    write_combined_counts_tsv(tsv, DATASET_HEADLINES, cellline_sets, tissue_rows, accession_sets, tissue_protein_rows=tissue_protein_rows, tissue_protein_rows_per_cohort=tissue_protein_rows_per_cohort, runs_per_cohort=runs_per_cohort)
    print(f'  saved: {tsv}')
    ds_order = ['PXD003539', 'PXD030304', 'PXD004701', 'PXD017199', 'PXD041421']
    print('\nPanel B (cell-line set intersections) region sizes:')
    for region, n in _set_region_sizes(cellline_sets, ds_order).items():
        print(f'  {region:<48s} {n:>6,}')
    print('\nPanel D (accession set intersections) region sizes:')
    for region, n in _set_region_sizes(accession_sets, ds_order).items():
        print(f'  {region:<48s} {n:>6,}')
    return 0


# ======================================================================
# inlined from analysis/figure_fig2_validation.py
# ======================================================================

"""Fig 2 - quantmsdiann validation row (scaling + ProteoBench accuracy).

One row of three panels (the former Fig 1 b/c/d, now a standalone figure so the
workflow gets Fig 1 to itself):
  (a) wall-clock versus cluster nodes (PXD071075 single-cell sweep)
  (b) wall-clock to finish each reanalysis, one bar per dataset
  (c) ProteoBench quantification-accuracy concordance vs standalone DIA-NN

Reuses the existing per-panel renderers (composite/ax mode) so the numbers stay
identical to the standalone figures.

Out: analysis/figures/manuscript/fig2_validation.svg
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
fs.apply_house_style()
_figure_fig2_validation__REPO = Path(__file__).resolve().parents[1]
PERF = _figure_fig2_validation__REPO / 'analysis' / 'figures' / 'performance' / 'data'
_figure_fig2_validation__OUT = _figure_fig2_validation__REPO / 'analysis' / 'figures' / 'manuscript' / 'fig2_validation.svg'

def _figure_fig2_validation__render(out: Path) -> Path:
    from matplotlib.patches import Patch
    dq = pd.read_csv(PERF / 'queue_size_sweep.tsv', sep='\t')
    dp = pd.read_csv(PERF / 'parallelism_data.tsv', sep='\t')
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 4.6))
    render_queue_size_sweep(dq, ax=ax[0], composite=True)
    render_parallelism_scatter(dp, ax=ax[1], composite=True, show_legend=False, short_labels=True)
    ds_colors = {}
    for ds in _COMMUNITY_COMPARATOR_DATASETS:
        inst = dp.loc[dp['dataset'] == ds, 'instrument']
        if len(inst):
            ds_colors[ds] = INSTRUMENT_COLOURS.get(inst.iloc[0], '#9e9e9e')
    acc.draw_strip(ax[2], compact=True, dataset_colors=ds_colors)
    for a, lab in zip(ax, 'abc'):
        a.text(-0.06, 1.05, f'({lab})', transform=a.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
    insts = [i for i in dict.fromkeys(dp['instrument']) if isinstance(i, str)]
    handles = [Patch(facecolor=INSTRUMENT_COLOURS.get(i, '#9e9e9e'), edgecolor='#222222', label=i) for i in insts]
    fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=7, frameon=False, title='Instrument / dataset (panels b, c)', title_fontsize=7.5, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.13, 1, 1), w_pad=0.4)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out

def figure_fig2_validation_main() -> int:
    print(f'wrote {_figure_fig2_validation__render(_figure_fig2_validation__OUT)}')
    return 0


# ======================================================================
# inlined from analysis/figure_quantmsdiann_benchmarks_vs_proteobench.py
# ======================================================================

"""quantmsdiann vs ProteoBench benchmarks figure.

Each of the four ProteoBench benchmark datasets has been processed through
quantmsdiann across DIA-NN versions (1.8.1, 2.5.1, 2.5.1-enterprise). This script:

- Reads the quantmsdiann headline counts from `report_counts.tsv` (produced by
  the `report_counts` stage of scripts/rebuild.py under the methods.md §1 rule:
  global totals at Lib.PG.Q.Value / Lib.Q.Value, per-run at PG.Q.Value /
  Q.Value; no contaminant/target filter, zeros counted).
- Fetches public ProteoBench submissions for the matching module from the
  Proteobench/Results_quant_ion_DIA_<module> repo (cached on disk).
- Renders the paper-ready figure sets and an auditable counts.tsv.

ProteoBench ion-level modules report precursors only — there is no protein
count in their public datapoints, so the head-to-head supp is precursors
only.
"""
import json
import re
import sys
from pathlib import Path
from typing import Iterable, Iterator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
import requests
_figure_quantmsdiann_benchmarks_vs_proteobench__REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = _figure_quantmsdiann_benchmarks_vs_proteobench__REPO_ROOT / 'data' / 'quantmsdiann_benchmarks'
_figure_quantmsdiann_benchmarks_vs_proteobench__FIGURES_DIR = _figure_quantmsdiann_benchmarks_vs_proteobench__REPO_ROOT / 'analysis' / 'figures' / 'quantmsdiann_benchmarks'
SUPP_DIR = _figure_quantmsdiann_benchmarks_vs_proteobench__FIGURES_DIR / 'supplementary'
FIG_DATA_DIR = _figure_quantmsdiann_benchmarks_vs_proteobench__FIGURES_DIR / 'data'
PRIDE_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/proteobench/quantmsdiann_results'
DIANN_VERSIONS = ('v1_8_1', 'v2_5_1', 'v2_5_1_enterprise')
_figure_quantmsdiann_benchmarks_vs_proteobench__DATASET_TO_MODULE: dict[str, dict[str, str]] = {'PXD049412': {'label': 'Module 9 - DIA single-cell', 'results_repo': 'Proteobench/Results_quant_ion_DIA_singlecell'}, 'PXD062685': {'label': 'Module 5 - DIA diaPASEF', 'results_repo': 'Proteobench/Results_quant_ion_DIA_diaPASEF'}, 'PXD070049': {'label': 'Module 10 - DIA ZenoTOF', 'results_repo': 'Proteobench/Results_quant_ion_DIA_ZenoTOF'}, 'ProteoBench_Module_7': {'label': 'Module 7 - DIA Astral 2Th', 'results_repo': 'Proteobench/Results_quant_ion_DIA_Astral'}}
GITHUB_API_BASE = 'https://api.github.com'
GITHUB_RAW_BASE = 'https://raw.githubusercontent.com'
_figure_quantmsdiann_benchmarks_vs_proteobench__SUMMARY_LOG_PRECURSOR_LINE_RE = re.compile('Target precursors at 1% global q-value:\\s*(\\d+)')

def _figure_quantmsdiann_benchmarks_vs_proteobench__parse_diann_summary_log(log_path: Path) -> tuple[int, int]:
    """Return (protein_groups, target_precursors) from a DIA-NN summary log."""
    protein_groups: int | None = None
    precursors: int | None = None
    with open(log_path, encoding='utf-8') as fh:
        for line in fh:
            m = SUMMARY_LOG_PROTEIN_LINE_RE.search(line)
            if m and protein_groups is None:
                protein_groups = int(m.group(1))
                continue
            m = _figure_quantmsdiann_benchmarks_vs_proteobench__SUMMARY_LOG_PRECURSOR_LINE_RE.search(line)
            if m and precursors is None:
                precursors = int(m.group(1))
    if protein_groups is None:
        raise ValueError("'Protein groups with global q-value <= 0.01: N' not found in log")
    if precursors is None:
        raise ValueError("'Target precursors at 1% global q-value: N' not found in log")
    return (protein_groups, precursors)
LIBRARY_KIND_EMPIRICAL = 'empirical'
LIBRARY_KIND_PREDICTED = 'predicted (DIANN)'
LIBRARY_KIND_USER_DEFINED = 'user-defined speclib'
LIBRARY_KIND_OTHER_TOOL = 'other tool'
QUANTMSDIANN_LIBRARY_KIND = LIBRARY_KIND_PREDICTED

def classify_predictors_library(predictors_library) -> str:
    """Map a ProteoBench `predictors_library` value to one of three categories
    used to color DIA-NN community submissions:
      - 'empirical' (None / 'None') — the submission loaded a pre-existing
        experimental/curated spectral library WITHOUT a predictor. This does
        NOT match quantmsdiann, which predicts its library from the FASTA.
      - 'predicted (DIANN)' — DIA-NN's built-in in-silico prediction (the
        FASTA-derived predicted library). This is quantmsdiann's strategy
        (see QUANTMSDIANN_LIBRARY_KIND) and the apples-to-apples set.
      - 'user-defined speclib' — externally uploaded/curated library.
    Any other value returns the original string."""
    if predictors_library is None or predictors_library == 'None':
        return LIBRARY_KIND_EMPIRICAL
    if isinstance(predictors_library, dict):
        vals = set(predictors_library.values())
        if vals == {'DIANN'}:
            return LIBRARY_KIND_PREDICTED
        if vals == {'User defined speclib'}:
            return LIBRARY_KIND_USER_DEFINED
        return ','.join(sorted(vals))
    s = str(predictors_library).strip()
    if s.startswith('{') and s.endswith('}'):
        try:
            import ast
            d = ast.literal_eval(s)
            if isinstance(d, dict):
                return classify_predictors_library(d)
        except (ValueError, SyntaxError):
            pass
    if 'User' in s and 'speclib' in s.lower():
        return LIBRARY_KIND_USER_DEFINED
    if s.upper() == 'DIANN':
        return LIBRARY_KIND_PREDICTED
    return s

def parse_proteobench_datapoints(json_path: Path) -> Iterator[tuple[str, str, int, str]]:
    """Yield (software_name, software_version, nr_prec, library_kind) per
    ProteoBench submission. The library kind is derived from each entry's
    `predictors_library` field via `classify_predictors_library`; non-DIA-NN
    tools are tagged 'other tool' since their library systems aren't directly
    comparable to DIA-NN's three categories.

    Skips entries that don't carry an integer nr_prec — those are rare
    submission artefacts and rendering them as zero would distort the
    per-module distribution."""
    with open(json_path, encoding='utf-8') as fh:
        payload = json.load(fh)
    for entry in payload:
        nr_prec = entry.get('nr_prec')
        if not isinstance(nr_prec, (int, float)) or pd.isna(nr_prec):
            continue
        software = entry.get('software_name') or ''
        version = entry.get('software_version') or ''
        if normalise_software_name(str(software)) == 'dia-nn':
            kind = classify_predictors_library(entry.get('predictors_library'))
        else:
            kind = LIBRARY_KIND_OTHER_TOOL
        yield (str(software), str(version), int(nr_prec), kind)

def extract_nr_prec_at_replicate_threshold(entry: dict, min_replicates: int) -> int | None:
    """Return the precursor count for a single ProteoBench submission
    `entry` at the requested ≥min_replicates threshold.

    Each ProteoBench submission embeds a `results` dict keyed by replicate
    count as a string ('1', '2', ..., '6'). The value at key 'K' is the
    set of metrics computed over precursors quantified in AT LEAST K of the
    6 ProteoBench replicates. So the field we want is:

        entry['results'][str(min_replicates)]['nr_prec']

    Returns None when the field is missing or not numeric. Falls back to
    the top-level `nr_prec` only when `min_replicates == 1` (the top-level
    nr_prec is defined as precursors quantified in ≥1 of 6 samples).

    Robbe Devreese (ProteoBench maintainer) flagged that DIA-NN 1.9.1's
    headline nr_prec is unusually high at the ≥1 threshold but drops back
    in line at ≥3; this function is the canonical extractor for the ≥3
    re-ranking used in the corrected supp figure.
    """
    results = entry.get('results')
    key = str(min_replicates)
    if isinstance(results, dict) and key in results:
        bucket = results[key]
        if isinstance(bucket, dict):
            v = bucket.get('nr_prec')
            if isinstance(v, (int, float)) and (not pd.isna(v)):
                return int(v)
    if min_replicates == 1:
        v = entry.get('nr_prec')
        if isinstance(v, (int, float)) and (not pd.isna(v)):
            return int(v)
    return None

def parse_proteobench_datapoints_at_threshold(json_path: Path, min_replicates: int) -> Iterator[tuple[str, str, int, str]]:
    """Same shape as `parse_proteobench_datapoints` but reads precursor
    counts from `entry['results'][str(min_replicates)]['nr_prec']`.

    Skips submissions that have no bucket at the requested threshold (rare;
    every DIA-NN submission inspected during the Slack-driven correction
    review carries all six replicate buckets, but defensive handling here
    keeps the renderer robust against any future submission that doesn't)."""
    with open(json_path, encoding='utf-8') as fh:
        payload = json.load(fh)
    for entry in payload:
        nr_prec = extract_nr_prec_at_replicate_threshold(entry, min_replicates)
        if nr_prec is None:
            continue
        software = entry.get('software_name') or ''
        version = entry.get('software_version') or ''
        if normalise_software_name(str(software)) == 'dia-nn':
            kind = classify_predictors_library(entry.get('predictors_library'))
        else:
            kind = LIBRARY_KIND_OTHER_TOOL
        yield (str(software), str(version), int(nr_prec), kind)

def normalise_software_name(name: str) -> str:
    """Collapse case/whitespace variants. DIA-NN appears as 'DIA-NN', 'DIANN',
    'Diann' across ProteoBench submissions; we normalise so highlight
    overlays match every spelling."""
    s = name.strip().lower()
    if s in {'dia-nn', 'diann', 'dia nn'}:
        return 'dia-nn'
    return s

def fetch_proteobench_module(repo: str, dest: Path) -> Path:
    """Fetch every `<hash>.json` datapoint from a Proteobench results repo
    and write them as a single consolidated list to `dest`. Idempotent:
    skips if `dest` already exists and is non-empty.

    The GitHub contents API is rate-limited (60 req/h unauthenticated) but
    we hit it once per module to list, plus N raw.githubusercontent.com
    fetches; the raw host is not rate-limited."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    listing_url = f'{GITHUB_API_BASE}/repos/{repo}/contents/'
    resp = requests.get(listing_url, timeout=60)
    resp.raise_for_status()
    listing = resp.json()
    if not isinstance(listing, list):
        raise RuntimeError(f'GitHub listing for {repo} returned non-list payload: {listing}')
    json_names = sorted((item['name'] for item in listing if isinstance(item, dict) and item.get('name', '').endswith('.json')))
    items: list[dict] = []
    for name in json_names:
        raw_url = f'{GITHUB_RAW_BASE}/{repo}/main/{name}'
        r = requests.get(raw_url, timeout=60)
        r.raise_for_status()
        items.append(r.json())
    with open(dest, 'w', encoding='utf-8') as fh:
        json.dump(items, fh)
    return dest

def consolidate_proteobench_datapoints(files: Iterable[Path], dest: Path) -> Path:
    """Combine per-submission JSON files into a single list, sorted by file
    name. Used by tests to exercise the same merge logic as
    `fetch_proteobench_module` without hitting the network."""
    items: list[dict] = []
    for p in sorted(files, key=lambda x: x.name):
        with open(p, encoding='utf-8') as fh:
            items.append(json.load(fh))
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, 'w', encoding='utf-8') as fh:
        json.dump(items, fh)
    return dest

def build_long_table(quantmsdiann_rows: list[tuple[str, str, int, int]], proteobench_rows: dict[str, list[tuple[str, str, int]]]) -> pd.DataFrame:
    """Combine quantmsdiann rows (dataset, version, precursors, proteins)
    and ProteoBench rows (per dataset: list of (tool, version, precursors))
    into a long-format DataFrame with columns
    [dataset, source, tool, version, precursors, proteins]."""
    rows: list[dict] = []
    for dataset, version, precursors, proteins in quantmsdiann_rows:
        rows.append({'dataset': dataset, 'source': 'quantmsdiann', 'tool': 'DIA-NN', 'version': version, 'precursors': precursors, 'proteins': proteins, 'library_kind': QUANTMSDIANN_LIBRARY_KIND})
    for dataset, entries in proteobench_rows.items():
        for tool, version, precursors, lib_kind in entries:
            rows.append({'dataset': dataset, 'source': 'proteobench', 'tool': tool, 'version': version, 'precursors': precursors, 'proteins': None, 'library_kind': lib_kind})
    return pd.DataFrame(rows)
_VERSION_LABELS = {'v1_8_1': '1.8.1', 'v2_1_0': '2.1.0', 'v2_2_0': '2.2.0', 'v2_3_2': '2.3.2', 'v2_5_0': '2.5.0', 'v2_5_1': '2.5.1', 'v2_5_1_enterprise': '2.5.1 ent.'}
_VERSION_COLORS = dict(fs.VERSION_COLORS)
_VERSION_MARKERS = {'v1_8_1': 'o', 'v2_5_1': 's', 'v2_5_1_enterprise': 'D'}

def _render_main_metric(quantmsdiann_rows: list[tuple[str, str, int, int]], svg_path: Path, *, metric: str, ylabel: str, label_fmt) -> None:
    """Grouped-bar panel of `metric` (`precursors` or `proteins`) across the
    4 datasets x 3 DIA-NN versions. Paper-ready: no title, no footer."""
    df = pd.DataFrame(quantmsdiann_rows, columns=['dataset', 'version', 'precursors', 'proteins'])
    datasets = sorted(df['dataset'].unique(), key=_dataset_sort_key)
    fig, ax = plt.subplots(figsize=(10.5, 5))
    n_versions = len(DIANN_VERSIONS)
    bar_width = 0.8 / n_versions
    x = list(range(len(datasets)))
    for k, version in enumerate(DIANN_VERSIONS):
        vals = [int(df[(df['dataset'] == d) & (df['version'] == version)][metric].iloc[0]) if not df[(df['dataset'] == d) & (df['version'] == version)].empty else 0 for d in datasets]
        offsets = [xi + (k - (n_versions - 1) / 2) * bar_width for xi in x]
        bars = ax.bar(offsets, vals, width=bar_width, color=_VERSION_COLORS.get(version, '#1f77b4'), label=_VERSION_LABELS.get(version, version))
        for bar, v in zip(bars, vals):
            if v == 0:
                continue
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), label_fmt(v), ha='center', va='bottom', fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([_dataset_display_label(d).replace('\n', '\n') for d in datasets], fontsize=9)
    ax.set_ylabel(ylabel)
    ymax = df[metric].max() * 1.15
    ax.set_ylim(0, ymax)
    fs.kfmt_axis(ax.yaxis)
    fs.despine(ax)
    ax.legend(title='DIA-NN version', loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=n_versions, frameon=False, fontsize=8, title_fontsize=9)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def render_main_precursors(quantmsdiann_rows: list[tuple[str, str, int, int]], svg_path: Path) -> None:
    """Precursor-only headline: 4 datasets x 3 DIA-NN versions, grouped bars.
    Precursor counts come from the DIA-NN report (run-specific q at the
    per-version recommended cut-off + 1% global)."""
    _render_main_metric(quantmsdiann_rows, svg_path, metric='precursors', ylabel='Precursors quantified (1% global FDR)', label_fmt=lambda v: f'{v / 1000:.0f}k')

def render_main_proteins(quantmsdiann_rows: list[tuple[str, str, int, int]], svg_path: Path) -> None:
    """Protein-group headline: 4 datasets x 3 DIA-NN versions, grouped bars.
    Protein groups are the per-run average, counted from the DIA-NN report at
    1% run-specific protein-group q-value (not the pg_matrix row count). This
    per-run depth metric is sensitive to the version improvement; the
    complete-profile (all-runs) companion panel is in the supplementary."""
    _render_main_metric(quantmsdiann_rows, svg_path, metric='proteins', ylabel='Protein groups per run (1% FDR)', label_fmt=lambda v: f'{v / 1000:.1f}k')

def _dataset_sort_key(name: str) -> tuple[int, str]:
    if name == 'ProteoBench_Module_7':
        return (0, name)
    return (1, name)

def _dataset_display_label(name: str) -> str:
    info = _figure_quantmsdiann_benchmarks_vs_proteobench__DATASET_TO_MODULE.get(name)
    if info is None:
        return name
    return f"{name}\n{info['label']}"

def _figure_quantmsdiann_benchmarks_vs_proteobench__write_counts_tsv(long_df: pd.DataFrame, tsv_path: Path, *, quantmsdiann_unfiltered_rows: list[tuple[str, str, int, int]] | None=None) -> None:
    """Write the per-dataset counts table for the benchmark figure.

    Columns: dataset, source, tool, version, precursors, proteins,
    filter_policy. `filter_policy` is `"target_only"` on the headline
    rows (rows consumed by the figures) and `"unfiltered"` on the
    companion rows derived from the raw pr_matrix.tsv / pg_matrix.tsv
    line counts. Both rows are written for every (dataset, version)
    pair so a reviewer can audit the conservative
    contaminant/entrapment/decoy filter delta (~1 % per
    ProteoBench DIA module).

    `quantmsdiann_unfiltered_rows` (optional) carries the
    (dataset, version, precursors, proteins) tuples computed BEFORE the
    filter — when present, those rows are appended with
    `filter_policy = "unfiltered"` and `source = "quantmsdiann"`.
    ProteoBench rows always carry `filter_policy = "target_only"`
    because community submissions are filtered upstream by ProteoBench's
    own `contaminant_flag = Cont_` parser.
    """
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ['dataset', 'source', 'tool', 'version', 'precursors', 'proteins', 'filter_policy']
    qm = long_df[long_df['source'] == 'quantmsdiann'].sort_values(['dataset', 'version'])
    pb = long_df[long_df['source'] == 'proteobench'].sort_values(['dataset', 'precursors'], ascending=[True, False])
    out = pd.concat([qm, pb], ignore_index=True).copy()
    out['filter_policy'] = 'target_only'
    if quantmsdiann_unfiltered_rows:
        unf_rows = []
        for dataset, version, precursors, proteins in quantmsdiann_unfiltered_rows:
            unf_rows.append({'dataset': dataset, 'source': 'quantmsdiann', 'tool': 'DIA-NN', 'version': version, 'precursors': precursors, 'proteins': proteins, 'library_kind': QUANTMSDIANN_LIBRARY_KIND, 'filter_policy': 'unfiltered'})
        out = pd.concat([out, pd.DataFrame(unf_rows)], ignore_index=True)
    out = out[cols].copy()
    out['precursors'] = out['precursors'].astype('Int64')
    out['proteins'] = out['proteins'].astype('Int64')
    out.to_csv(tsv_path, sep='\t', index=False)

def median_nr_prec_per_version(proteobench_rows_by_threshold: dict[int, dict[str, list[tuple[str, str, int, str]]]], quantmsdiann_rows_by_threshold: dict[int, list[tuple[str, str, int, int]]] | None=None) -> pd.DataFrame:
    """Build a long-format table of median precursor counts per DIA-NN
    version per module at each replicate threshold. Used to quantify the
    Robbe-predicted 1.9.1 vs 2.3.0 ranking flip across the ≥1 and ≥3
    thresholds.

    `proteobench_rows_by_threshold` is keyed by min_replicates -> dataset ->
    list of (tool, version, nr_prec, library_kind). Optional
    `quantmsdiann_rows_by_threshold` adds quantmsdiann rows for the same
    table (one quantmsdiann analysis per version, so 'median' == the single
    value)."""
    rows = []
    for thr, by_dataset in proteobench_rows_by_threshold.items():
        for dataset, entries in by_dataset.items():
            diann = [(v, n) for tool, v, n, _k in entries if normalise_software_name(tool) == 'dia-nn']
            by_v: dict[str, list[int]] = {}
            for v, n in diann:
                by_v.setdefault(v.strip(), []).append(n)
            for v, ns in sorted(by_v.items()):
                rows.append({'dataset': dataset, 'min_replicates': thr, 'source': 'proteobench', 'version': v, 'n_submissions': len(ns), 'median_nr_prec': int(pd.Series(ns).median())})
    if quantmsdiann_rows_by_threshold is not None:
        for thr, qm_rows in quantmsdiann_rows_by_threshold.items():
            for dataset, version, nr_prec, _ in qm_rows:
                rows.append({'dataset': dataset, 'min_replicates': thr, 'source': 'quantmsdiann', 'version': _VERSION_LABELS.get(version, version), 'n_submissions': 1, 'median_nr_prec': int(nr_prec)})
    return pd.DataFrame(rows)

def write_median_table(df: pd.DataFrame, tsv_path: Path, *, quantmsdiann_rows_by_threshold_unfiltered: dict[int, list[tuple[str, str, int, int]]] | None=None) -> None:
    """Auditable TSV of median precursor counts per DIA-NN version per module
    at the ≥1 and ≥3 replicate thresholds. One row per
    (dataset, threshold, source, version).

    Headline rows carry `filter_policy = "target_only"` (rows whose
    Protein.Group passes the conservative contaminant filter). When
    `quantmsdiann_rows_by_threshold_unfiltered` is provided, this writer
    appends companion `filter_policy = "unfiltered"` rows so a reviewer
    can compare. ProteoBench rows are always `target_only` because
    community submissions are filtered upstream by ProteoBench's own
    `contaminant_flag = Cont_` parser."""
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    if 'filter_policy' not in out.columns:
        out['filter_policy'] = 'target_only'
    if quantmsdiann_rows_by_threshold_unfiltered:
        extra_rows = []
        for thr, rows_unf in quantmsdiann_rows_by_threshold_unfiltered.items():
            for dataset, version, nr_prec, _ in rows_unf:
                extra_rows.append({'dataset': dataset, 'min_replicates': thr, 'source': 'quantmsdiann', 'version': _VERSION_LABELS.get(version, version), 'n_submissions': 1, 'median_nr_prec': int(nr_prec), 'filter_policy': 'unfiltered'})
        if extra_rows:
            out = pd.concat([out, pd.DataFrame(extra_rows)], ignore_index=True)
    out = out.sort_values(['dataset', 'min_replicates', 'source', 'version', 'filter_policy'])
    out.to_csv(tsv_path, sep='\t', index=False)
REPORT_COUNTS_PATH = DATA_DIR / 'report_counts.tsv'

def load_report_counts() -> dict[tuple[str, str], dict]:
    """Precursor + protein-group counts read from the DIA-NN *report*
    (`diann_report.parquet` / `.tsv`), NOT the `*_matrix.tsv` files.

    The matrices bake in `--matrix-spec-q` (0.05 run-specific) and, because the
    pipeline sets `--qvalue` to 0.01 for v1.8.1 but 0.05 for v2.5.1/enterprise,
    matrix row counts are filtered at *different* run-specific q-values per
    version and are therefore not comparable. These counts instead follow the
    Vadim filter rule (see analysis.count_report_ids): global totals on
    `Lib.Q.Value` / `Lib.PG.Q.Value` (`prec_global` / `prot_global`),
    replicate precursors on `Q.Value` (`prec_min1` / `prec_min3`), and per-run
    protein groups on `PG.Q.Value` (`prot_perrun_avg` / `prot_complete`). No
    contaminant/target filter. Recomputed from the public FTP reports by
    the rebuild `report_counts` stage and staged at
    data/quantmsdiann_benchmarks/report_counts.tsv. Keyed by (dataset, version)."""
    df = pd.read_csv(REPORT_COUNTS_PATH, sep='\t', dtype={'dataset': str, 'version': str})
    out: dict[tuple[str, str], dict] = {}
    for _, r in df.iterrows():
        out[r['dataset'], r['version']] = {k: int(r[k]) for k in ('prec_min1', 'prec_min3', 'prec_global', 'prot_global', 'prot_perrun_avg', 'prot_complete')}
    return out

def figure_quantmsdiann_benchmarks_vs_proteobench_main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    _figure_quantmsdiann_benchmarks_vs_proteobench__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    SUPP_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DATA_DIR.mkdir(parents=True, exist_ok=True)
    quantmsdiann_rows: list[tuple[str, str, int, int]] = []
    quantmsdiann_rows_min3: list[tuple[str, str, int, int]] = []
    quantmsdiann_rows_unfiltered: list[tuple[str, str, int, int]] = []
    quantmsdiann_rows_min3_unfiltered: list[tuple[str, str, int, int]] = []
    quantmsdiann_rows_protavg: list[tuple[str, str, int, int]] = []
    quantmsdiann_rows_complete: list[tuple[str, str, int, int]] = []
    report_counts = load_report_counts()
    for dataset in _figure_quantmsdiann_benchmarks_vs_proteobench__DATASET_TO_MODULE:
        for version in DIANN_VERSIONS:
            c = report_counts.get((dataset, version))
            if c is None:
                print(f'WARN: no report_counts for {dataset}/{version}', file=sys.stderr)
                continue
            precursors = c['prec_global']
            precursors_min3 = c['prec_min3']
            proteins = c['prot_global']
            quantmsdiann_rows.append((dataset, version, precursors, proteins))
            quantmsdiann_rows_min3.append((dataset, version, precursors_min3, proteins))
            quantmsdiann_rows_unfiltered.append((dataset, version, precursors, proteins))
            quantmsdiann_rows_min3_unfiltered.append((dataset, version, precursors_min3, proteins))
            quantmsdiann_rows_protavg.append((dataset, version, precursors, c['prot_perrun_avg']))
            quantmsdiann_rows_complete.append((dataset, version, precursors, c['prot_complete']))
            print(f"{dataset} {version}: precursors_global={precursors:,}  precursors_min3={precursors_min3:,}  proteins_global={proteins:,}  prot_perrun_avg={c['prot_perrun_avg']:,}  [report, Vadim rule]")
    proteobench_rows: dict[str, list[tuple[str, str, int, str]]] = {}
    proteobench_rows_min3: dict[str, list[tuple[str, str, int, str]]] = {}
    for dataset, info in _figure_quantmsdiann_benchmarks_vs_proteobench__DATASET_TO_MODULE.items():
        cache = DATA_DIR / 'proteobench' / f'{dataset}.json'
        try:
            fetch_proteobench_module(info['results_repo'], cache)
        except Exception as exc:
            print(f'WARN: failed to fetch ProteoBench for {dataset}: {exc}', file=sys.stderr)
            proteobench_rows[dataset] = []
            proteobench_rows_min3[dataset] = []
            continue
        entries_min1 = list(parse_proteobench_datapoints_at_threshold(cache, 1))
        entries_min3 = list(parse_proteobench_datapoints_at_threshold(cache, 3))
        proteobench_rows[dataset] = entries_min1
        proteobench_rows_min3[dataset] = entries_min3
        print(f'{dataset}: {len(entries_min1)} ProteoBench submissions (min1) / {len(entries_min3)} (min3)')
    long_df = build_long_table(quantmsdiann_rows, proteobench_rows)
    long_df_min3 = build_long_table(quantmsdiann_rows_min3, proteobench_rows_min3)
    print('Rendering precursor-only main panel...')
    render_main_precursors(quantmsdiann_rows, _figure_quantmsdiann_benchmarks_vs_proteobench__FIGURES_DIR / 'main_benchmarks_precursors.svg')
    print('Rendering protein-group main panel (per-run average)...')
    render_main_proteins(quantmsdiann_rows_protavg, _figure_quantmsdiann_benchmarks_vs_proteobench__FIGURES_DIR / 'main_benchmarks_proteins.svg')
    print('Rendering complete-profile protein panel (supplementary)...')
    _render_main_metric(quantmsdiann_rows_complete, SUPP_DIR / 'supp_benchmarks_proteins_complete.svg', metric='proteins', ylabel='Complete-profile protein groups (in all runs)', label_fmt=lambda v: f'{v / 1000:.1f}k')
    print('Writing auditable counts TSV (≥1)...')
    _figure_quantmsdiann_benchmarks_vs_proteobench__write_counts_tsv(long_df, FIG_DATA_DIR / 'counts.tsv', quantmsdiann_unfiltered_rows=quantmsdiann_rows_unfiltered)
    print('Writing auditable counts TSV (≥3)...')
    _figure_quantmsdiann_benchmarks_vs_proteobench__write_counts_tsv(long_df_min3, FIG_DATA_DIR / 'counts_min3.tsv', quantmsdiann_unfiltered_rows=quantmsdiann_rows_min3_unfiltered)
    print('Writing per-DIA-NN-version median precursor table (≥1 and ≥3)...')
    median_df = median_nr_prec_per_version({1: proteobench_rows, 3: proteobench_rows_min3}, {1: quantmsdiann_rows, 3: quantmsdiann_rows_min3})
    write_median_table(median_df, FIG_DATA_DIR / 'median_nr_prec_by_version.tsv', quantmsdiann_rows_by_threshold_unfiltered={1: quantmsdiann_rows_unfiltered, 3: quantmsdiann_rows_min3_unfiltered})
    return 0


# ======================================================================
# inlined from analysis/figure_id_vs_epsilon.py
# ======================================================================

"""F1b — Identifications vs accuracy scatter.

Per dataset, one panel: x = `nr_prec` at the chosen replicate threshold
(default ≥3), y = `median_abs_epsilon_global` at the same threshold.

- Community ProteoBench submissions: grey cloud (one dot per submission),
  faceted by `predictors_library` palette.
- quantmsdiann (this work): red dots, one per DIA-NN version, joined by a
  thin line in version order so the reader can read the progression.
- Per-dataset metric source: locally-cached `proteobench_metrics/*.json`
  populated by `analysis/proteobench_metrics.py`. If the cache is empty,
  the figure renders an explainer panel rather than crashing.

The accompanying F1c per-species log2 strip plot lives in a sibling
function `render_per_species_log2` so the supp variant ships from the
same script.
"""
import json
import sys
from pathlib import Path
from typing import Iterable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import numpy as np
import pandas as pd
_figure_id_vs_epsilon__REPO_ROOT = Path(__file__).resolve().parent.parent
_figure_id_vs_epsilon__FIGURES_DIR = _figure_id_vs_epsilon__REPO_ROOT / 'analysis' / 'figures' / 'quantmsdiann_benchmarks'
LIB_PALETTE = {LIBRARY_KIND_EMPIRICAL: '#80cbc4', LIBRARY_KIND_PREDICTED: '#1976d2', LIBRARY_KIND_USER_DEFINED: '#9575cd', LIBRARY_KIND_OTHER_TOOL: '#bdbdbd'}

def extract_community_id_vs_eps(dataset: str, threshold: int) -> pd.DataFrame:
    """Load the cached ProteoBench submission JSON for a dataset and
    pull `(nr_prec, median_abs_epsilon_global, library_kind, software_name)`
    per submission at the given replicate threshold. Submissions
    missing either metric at this threshold are dropped."""
    cache = DATA_DIR / 'proteobench' / f'{dataset}.json'
    if not cache.exists():
        return pd.DataFrame(columns=['nr_prec', 'median_abs_epsilon_global', 'library_kind', 'software_name'])
    with open(cache, encoding='utf-8') as fh:
        entries = json.load(fh)
    rows: list[dict] = []
    for entry in entries:
        results = entry.get('results', {})
        thr_key = str(threshold)
        if thr_key not in results:
            continue
        res = results[thr_key]
        nr_prec = res.get('nr_prec')
        eps = res.get('median_abs_epsilon_global')
        if nr_prec is None or eps is None:
            continue
        try:
            nr_prec = int(nr_prec)
            eps = float(eps)
        except (TypeError, ValueError):
            continue
        if pd.isna(eps):
            continue
        rows.append({'nr_prec': nr_prec, 'median_abs_epsilon_global': eps, 'library_kind': classify_predictors_library(entry.get('predictors_library')), 'software_name': entry.get('software_name', '')})
    return pd.DataFrame(rows)

def extract_quantmsdiann_id_vs_eps(dataset: str, threshold: int) -> pd.DataFrame:
    """Pull `(version, nr_prec, median_abs_epsilon_global)` per
    quantmsdiann DIA-NN version from the locally-cached metrics. Rows
    where the metric cache hasn't been computed yet are silently
    skipped — the figure renders only what's available."""
    rows: list[dict] = []
    for version in DIANN_VERSIONS:
        cache_path = METRICS_CACHE_DIR / f'{dataset}_{version}.json'
        if not cache_path.exists():
            continue
        with open(cache_path, encoding='utf-8') as fh:
            payload = json.load(fh)
        res = payload.get('results', {}).get(str(threshold))
        if not res:
            continue
        nr_prec = res.get('nr_prec')
        eps = res.get('median_abs_epsilon_global')
        if nr_prec is None or eps is None or pd.isna(eps):
            continue
        rows.append({'version': version, 'version_label': _VERSION_LABELS.get(version, version), 'nr_prec': int(nr_prec), 'median_abs_epsilon_global': float(eps)})
    df = pd.DataFrame(rows)
    if len(df):
        df = df.assign(order=lambda d: d['version'].map({v: i for i, v in enumerate(DIANN_VERSIONS)})).sort_values('order').drop(columns='order').reset_index(drop=True)
    return df

def render_id_vs_epsilon(threshold: int, svg_path: Path, *, datasets: Iterable[str] | None=None) -> pd.DataFrame:
    """4-panel scatter (one per dataset). Returns a long-format
    DataFrame of every plotted point for the auditable TSV."""
    datasets = sorted(datasets, key=_dataset_sort_key) if datasets is not None else sorted(DATASET_TO_MODULE, key=_dataset_sort_key)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11.0, 8.2), sharex=False, sharey=False, squeeze=False)
    long_rows: list[dict] = []
    qm_drawn = False
    libs_drawn: set[str] = set()
    for idx, dataset in enumerate(datasets):
        ax = axes[idx // 2][idx % 2]
        community = extract_community_id_vs_eps(dataset, threshold)
        community = community[(community['library_kind'] == LIBRARY_KIND_PREDICTED) & (community['software_name'] == 'DIA-NN')].reset_index(drop=True)
        qm = extract_quantmsdiann_id_vs_eps(dataset, threshold)
        for kind in (LIBRARY_KIND_EMPIRICAL, LIBRARY_KIND_USER_DEFINED, LIBRARY_KIND_OTHER_TOOL, LIBRARY_KIND_PREDICTED):
            sub = community[community['library_kind'] == kind]
            if not len(sub):
                continue
            ax.scatter(sub['nr_prec'], sub['median_abs_epsilon_global'], s=40, c=LIB_PALETTE[kind], alpha=0.7, edgecolors='#555555', linewidths=0.4, label=kind if kind not in libs_drawn else None)
            libs_drawn.add(kind)
        for _, row in community.iterrows():
            long_rows.append({'dataset': dataset, 'threshold': threshold, 'source': 'proteobench-community', 'label': row['software_name'], 'library_kind': row['library_kind'], 'nr_prec': int(row['nr_prec']), 'median_abs_epsilon_global': float(row['median_abs_epsilon_global'])})
        if len(qm):
            ax.plot(qm['nr_prec'], qm['median_abs_epsilon_global'], color='#d62728', linewidth=1.0, alpha=0.7, zorder=2)
            ax.scatter(qm['nr_prec'], qm['median_abs_epsilon_global'], s=80, c='#d62728', edgecolors='#7f1d1d', linewidths=0.8, zorder=3, label='quantmsdiann (DIA-NN, predicted-from-FASTA lib)' if not qm_drawn else None)
            qm_drawn = True
            for _, row in qm.iterrows():
                ax.annotate(row['version_label'], xy=(row['nr_prec'], row['median_abs_epsilon_global']), xytext=(6, 6), textcoords='offset points', fontsize=7, color='#7f1d1d', fontweight='bold')
                long_rows.append({'dataset': dataset, 'threshold': threshold, 'source': 'quantmsdiann', 'label': f"quantmsdiann {row['version_label']}", 'library_kind': LIBRARY_KIND_PREDICTED, 'nr_prec': int(row['nr_prec']), 'median_abs_epsilon_global': float(row['median_abs_epsilon_global'])})
        ax.set_title(_dataset_display_label(dataset).splitlines()[0], loc='left', fontsize=10, fontweight='bold')
        ax.set_xlabel(f'Precursors quantified (≥{threshold} rep)', fontsize=9)
        ax.set_ylabel('Median |ε| (lower = more accurate)', fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', labelsize=8)
        y_all = list(qm['median_abs_epsilon_global']) + list(community['median_abs_epsilon_global'])
        if y_all:
            MIN_EPS_SPAN = 0.05
            y_lo, y_hi = (min(y_all), max(y_all))
            if y_hi - y_lo < MIN_EPS_SPAN:
                mid = (y_lo + y_hi) / 2.0
                ax.set_ylim(mid - MIN_EPS_SPAN / 2.0, mid + MIN_EPS_SPAN / 2.0)
        if len(qm) and len(community):
            n_lead_x = int((community['nr_prec'] < qm['nr_prec'].median()).sum())
            ax.text(0.02, 0.97, f'PB community: n={len(community)}    qm versions: {len(qm)}', transform=ax.transAxes, ha='left', va='top', fontsize=7, color='#666666')
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles: list = []
    if qm_drawn:
        handles.append(Line2D([0], [0], marker='o', color='#d62728', markerfacecolor='#d62728', markeredgecolor='#7f1d1d', linewidth=1.2, markersize=8, label='quantmsdiann (distributed, this work)'))
    LEGEND_LABEL = {LIBRARY_KIND_PREDICTED: 'single-machine community (ProteoBench, predicted lib)'}
    for k in (LIBRARY_KIND_EMPIRICAL, LIBRARY_KIND_PREDICTED, LIBRARY_KIND_USER_DEFINED, LIBRARY_KIND_OTHER_TOOL):
        if k in libs_drawn:
            handles.append(Patch(facecolor=LIB_PALETTE[k], label=LEGEND_LABEL.get(k, k)))
    if handles:
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=min(5, len(handles)), fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    return pd.DataFrame(long_rows)

def extract_qm_per_species_log2(dataset: str, threshold: int) -> pd.DataFrame:
    """Pull `mean_log2_empirical_<SPECIES>` per quantmsdiann version
    for the supplementary panel. Missing-species rows are dropped (the
    singlecell module has no E. coli)."""
    rows: list[dict] = []
    for version in DIANN_VERSIONS:
        cache_path = METRICS_CACHE_DIR / f'{dataset}_{version}.json'
        if not cache_path.exists():
            continue
        with open(cache_path, encoding='utf-8') as fh:
            payload = json.load(fh)
        res = payload.get('results', {}).get(str(threshold), {})
        for species in ('HUMAN', 'YEAST', 'ECOLI'):
            key = f'mean_log2_empirical_{species}'
            value = res.get(key)
            if value is None:
                continue
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if pd.isna(value):
                continue
            median = res.get(f'median_log2_empirical_{species}')
            try:
                median = float(median)
            except (TypeError, ValueError):
                median = value
            if pd.isna(median):
                median = value
            rows.append({'dataset': dataset, 'version': version, 'version_label': _VERSION_LABELS.get(version, version), 'species': species, 'mean_log2_empirical': value, 'median_log2_empirical': median})
    return pd.DataFrame(rows)
SPECIES_EXPECTED_LOG2_A_vs_B = {'quant_lfq_DIA_ion_singlecell': {'HUMAN': np.log2(1.2), 'YEAST': np.log2(0.2)}, 'quant_lfq_DIA_ion_diaPASEF': {'HUMAN': 0.0, 'YEAST': 1.0, 'ECOLI': -2.0}, 'quant_lfq_DIA_ion_ZenoTOF': {'HUMAN': 0.0, 'YEAST': 1.0, 'ECOLI': -2.0}, 'quant_lfq_DIA_ion_Astral': {'HUMAN': 0.0, 'YEAST': 1.0, 'ECOLI': -2.0}}

def render_per_species_log2(threshold: int, svg_path: Path, *, datasets: Iterable[str] | None=None) -> pd.DataFrame:
    """Strip plot of `mean_log2_empirical_<SPECIES>` per quantmsdiann
    version per dataset, with an expected-ratio reference line per
    species. The supplementary view of F1c."""
    datasets = sorted(datasets, key=_dataset_sort_key) if datasets is not None else sorted(DATASET_TO_MODULE, key=_dataset_sort_key)
    fig, axes = plt.subplots(nrows=len(datasets), ncols=1, figsize=(8.2, 2.4 * len(datasets)), sharex=False, squeeze=False)
    species_colours = {'HUMAN': '#1976d2', 'YEAST': '#d62728', 'ECOLI': '#388e3c'}
    long_rows: list[dict] = []
    for i, dataset in enumerate(datasets):
        ax = axes[i, 0]
        df = extract_qm_per_species_log2(dataset, threshold)
        long_rows.extend(df.to_dict('records'))
        if df.empty:
            ax.text(0.5, 0.5, 'no cached metrics yet for this dataset', transform=ax.transAxes, ha='center', va='center', fontsize=9, color='#888888')
            ax.set_axis_off()
            continue
        species_order = list(df['species'].unique())
        for species in species_order:
            sub = df[df['species'] == species]
            ax.scatter(sub['version_label'], sub['mean_log2_empirical'], s=70, c=species_colours.get(species, '#9e9e9e'), edgecolors='#222222', linewidths=0.4, label=species, zorder=3)
            expected = SPECIES_EXPECTED_LOG2_A_vs_B.get(DATASET_TO_MODULE[dataset], {}).get(species)
            if expected is not None:
                ax.axhline(expected, color=species_colours.get(species, '#9e9e9e'), linewidth=0.8, linestyle='--', alpha=0.6, zorder=1)
        ax.set_title(_dataset_display_label(dataset).splitlines()[0], loc='left', fontsize=10, fontweight='bold')
        ax.set_ylabel('mean log2 (A/B)', fontsize=9)
        ax.axhline(0, color='#aaaaaa', linewidth=0.6, zorder=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', labelsize=8)
        ax.legend(loc='best', fontsize=7, frameon=False, ncol=3)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    return pd.DataFrame(long_rows)
_COMMUNITY_COMPARATOR_DATASETS = ('ProteoBench_Module_7', 'PXD062685')
_SPECIES_X = {'HUMAN': 0, 'YEAST': 1, 'ECOLI': 2}
_SPECIES_LABEL = {'HUMAN': 'Human', 'YEAST': 'Yeast', 'ECOLI': 'E. coli'}

def render_accuracy_panels(threshold: int, svg_path: Path, *, datasets: Iterable[str] | None=None) -> pd.DataFrame:
    """Fig 2 accuracy figure. Two stacked message-bearing panels:

      (b) Per-species fold-change accuracy (2x2, one cell per module). x =
          HYE species, y = measured log2 ratio; a dashed line marks the
          ProteoBench-expected ratio per species (accuracy = distance to the
          line) and the three DIA-NN configurations are overlaid (1.8.1 and
          2.5.1 light->dark blue, 2.5.1-enterprise amber) so their tight
          clustering shows that accuracy is version-invariant.

      (c) quantmsdiann within the predicted-library community (1x2, only the
          modules with predicted-library DIA-NN comparators). Box + strip of
          the community median |eps| with the three quantmsdiann
          configurations overlaid as a tight cluster.

    Returns the long-format audit table of every plotted point."""
    datasets = sorted(datasets, key=_dataset_sort_key) if datasets is not None else sorted(DATASET_TO_MODULE, key=_dataset_sort_key)
    fig = plt.figure(figsize=(9.5, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.95], hspace=0.62, wspace=0.26)
    long_rows: list[dict] = []
    for idx, dataset in enumerate(datasets):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        df = extract_qm_per_species_log2(dataset, threshold)
        expected = SPECIES_EXPECTED_LOG2_A_vs_B.get(DATASET_TO_MODULE.get(dataset, ''), {})
        present = [s for s in ('HUMAN', 'YEAST', 'ECOLI') if s in expected]
        for species in present:
            x = _SPECIES_X[species]
            ax.hlines(expected[species], x - 0.32, x + 0.32, color='#444444', ls='--', lw=1.1, zorder=1)
            sub = df[df['species'] == species]
            n_ver = len(DIANN_VERSIONS)
            for _, row in sub.iterrows():
                vi = DIANN_VERSIONS.index(row['version'])
                dx = (vi - (n_ver - 1) / 2) * 0.14
                ax.scatter(x + dx, row['mean_log2_empirical'], s=20, marker=_VERSION_MARKERS.get(row['version'], 'o'), color=_VERSION_COLORS.get(row['version'], '#1f77b4'), edgecolor='#333333', linewidths=0.3, zorder=3)
                long_rows.append({'panel': 'b', 'dataset': dataset, 'threshold': threshold, 'species': species, 'version': row['version'], 'measured_log2': float(row['mean_log2_empirical']), 'expected_log2': float(expected[species])})
        ax.set_xticks([_SPECIES_X[s] for s in present])
        ax.set_xticklabels([_SPECIES_LABEL[s] for s in present], fontsize=8)
        ax.set_xlim(-0.5, max((_SPECIES_X[s] for s in present)) + 0.5 if present else 2.5)
        _label_parts = _dataset_display_label(dataset).split('\n')
        ax.set_title(_label_parts[1] if len(_label_parts) > 1 else _label_parts[0], loc='left', fontsize=9, fontweight='bold')
        ax.set_ylabel('Measured log$_2$ ratio', fontsize=8)
        ax.tick_params(axis='both', labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    fig.text(0.005, 0.99, '(b) Per-species fold-change accuracy (dashed = expected ratio)', fontsize=10, fontweight='bold', va='top')
    from matplotlib.lines import Line2D
    version_handles = [Line2D([0], [0], marker=_VERSION_MARKERS.get(v, 'o'), linestyle='none', markersize=6, markerfacecolor=_VERSION_COLORS.get(v, '#1f77b4'), markeredgecolor='#333333', label=_VERSION_LABELS.get(v, v)) for v in DIANN_VERSIONS]
    fig.legend(handles=version_handles, loc='upper right', bbox_to_anchor=(0.99, 1.005), ncol=len(DIANN_VERSIONS), fontsize=7.5, frameon=False, title='DIA-NN version', title_fontsize=7.5, handletextpad=0.2, columnspacing=0.9)
    comp = [d for d in datasets if d in _COMMUNITY_COMPARATOR_DATASETS]
    panel_c_top = 0.0
    for j, dataset in enumerate(comp):
        ax = fig.add_subplot(gs[2, j])
        panel_c_top = max(panel_c_top, ax.get_position().y1)
        community_df = extract_community_id_vs_eps(dataset, threshold)
        community_df = community_df[(community_df['library_kind'] == LIBRARY_KIND_PREDICTED) & (community_df['software_name'] == 'DIA-NN')]
        for _, _crow in community_df.iterrows():
            long_rows.append({'panel': 'c', 'dataset': dataset, 'threshold': threshold, 'source': 'community', 'software_name': _crow['software_name'], 'median_abs_epsilon': float(_crow['median_abs_epsilon_global'])})
        community = community_df['median_abs_epsilon_global'].astype(float).values
        qm = extract_quantmsdiann_id_vs_eps(dataset, threshold)
        qm_eps = qm['median_abs_epsilon_global'].astype(float).values
        if len(community):
            ax.boxplot([community], positions=[0], widths=0.5, showfliers=False, medianprops=dict(color='#37474f', lw=1.2))
            cx = [0.0] if len(community) == 1 else list(np.linspace(-0.16, 0.16, len(community)))
            ax.scatter(cx, community, s=16, color='#90a4ae', edgecolor='#555555', linewidths=0.3, alpha=0.85, label=f'community (n={len(community)})', zorder=3)
        if len(qm_eps):
            ax.boxplot([qm_eps], positions=[1], widths=0.5, showfliers=False, boxprops=dict(color='#7f1d1d'), whiskerprops=dict(color='#7f1d1d'), capprops=dict(color='#7f1d1d'), medianprops=dict(color='#d62728', lw=1.2))
            qx = list(np.linspace(0.84, 1.16, len(qm_eps))) if len(qm_eps) > 1 else [1.0]
            qm_versions = list(qm['version'])
            for k, (xk, yk) in enumerate(zip(qx, qm_eps)):
                ver = qm_versions[k] if k < len(qm_versions) else None
                ax.scatter(xk, yk, s=22, marker=_VERSION_MARKERS.get(ver, 'o'), color='#d62728', edgecolor='#7f1d1d', linewidths=0.5, zorder=4, label=f'quantmsdiann ({len(DIANN_VERSIONS)} versions)' if k == 0 else None)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['community', 'quantmsdiann'], fontsize=8)
        ax.set_xlim(-0.6, 1.6)
        _label_parts = _dataset_display_label(dataset).split('\n')
        ax.set_title(_label_parts[1] if len(_label_parts) > 1 else _label_parts[0], loc='left', fontsize=9, fontweight='bold')
        ax.set_ylabel('Median |ε|', fontsize=8)
        ax.tick_params(axis='both', labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(fontsize=6.5, frameon=False, loc='upper right')
        for _, row in qm.iterrows():
            long_rows.append({'panel': 'c', 'dataset': dataset, 'threshold': threshold, 'species': None, 'version': row['version'], 'measured_log2': None, 'median_abs_epsilon': float(row['median_abs_epsilon_global'])})
    fig.text(0.005, panel_c_top + 0.025, '(c) quantmsdiann within the predicted-library community', fontsize=9, fontweight='bold', va='bottom')
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    return pd.DataFrame(long_rows)

def figure_id_vs_epsilon_main() -> int:
    _figure_id_vs_epsilon__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = _figure_id_vs_epsilon__FIGURES_DIR / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    supp_dir = _figure_id_vs_epsilon__FIGURES_DIR / 'supplementary'
    supp_dir.mkdir(parents=True, exist_ok=True)
    for dataset in DATASET_TO_MODULE:
        for version in DIANN_VERSIONS:
            cache_path = METRICS_CACHE_DIR / f'{dataset}_{version}.json'
            if not cache_path.exists():
                print(f'warn: missing metrics cache {dataset}/{version} — F1b will skip this point. Populate via cached_proteobench_metrics(...).')
    acc_df = render_accuracy_panels(threshold=3, svg_path=supp_dir / 'supp_accuracy_all_modules.svg')
    acc_df.to_csv(data_dir / 'accuracy_min3.tsv', sep='\t', index=False)
    print(f"Supp accuracy (≥3 rep): {acc_df.shape[0]} points rendered to {supp_dir / 'supp_accuracy_all_modules.svg'}")
    long_df = render_id_vs_epsilon(threshold=3, svg_path=supp_dir / 'supp_id_vs_epsilon_min3.svg')
    long_df.to_csv(data_dir / 'id_vs_epsilon_min3.tsv', sep='\t', index=False)
    print(f"supp id-vs-ε (≥3 rep): {long_df.shape[0]} points rendered to {supp_dir / 'supp_id_vs_epsilon_min3.svg'}")
    supp_df = render_per_species_log2(threshold=3, svg_path=supp_dir / 'supp_per_species_log2_min3.svg')
    supp_df.to_csv(data_dir / 'per_species_log2_min3.tsv', sep='\t', index=False)
    return 0


# ======================================================================
# inlined from analysis/figure_mdc_cluster_runtime.py
# ======================================================================

"""Per-step wall-clock of one quantmsdiann run on an independent (non-EBI) cluster.

Companion to the EBI per-step figure (`runtime_per_step.svg`): same horizontal
per-step box-plot style, rendered by the *same* function
(`figure_performance_trace.render_per_step_boxplot`) so panels (a) EBI and
(b) non-EBI in Supplementary Fig. S1 are directly comparable. Shows that
quantmsdiann runs unchanged on a third-party SLURM cluster — a 279-raw-file
cohort processed by a different group at the Max Delbruck Center (MDC Berlin)
under the SLURM+Singularity profile (DIA-NN 2.5.0, quantmsdiann v2.1.0; total
wall-clock 7h49m). The per-file stages run in parallel across the cluster, so
the summed task time collapses to a few hours of wall-clock; the run-level
facts live in the figure caption, matching the title-free house style of the
EBI panel.

Privacy: built ONLY from an anonymised per-step duration table (step name +
seconds + memory + %CPU); raw file names, sample identifiers, internal paths
and contact details from the source run are never imported or stored.

Source: analysis/figures/performance/data/mdc_step_durations.tsv
Out:    analysis/figures/performance/mdc_cluster_runtime.svg
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import pandas as pd
_figure_mdc_cluster_runtime__REPO = Path(__file__).resolve().parents[1]
_figure_mdc_cluster_runtime__DATA = _figure_mdc_cluster_runtime__REPO / 'analysis' / 'figures' / 'performance' / 'data' / 'mdc_step_durations.tsv'
EBI_TSV = _figure_mdc_cluster_runtime__REPO / 'analysis' / 'figures' / 'performance' / 'data' / 'runtime_per_step.tsv'
_figure_mdc_cluster_runtime__OUT = _figure_mdc_cluster_runtime__REPO / 'analysis' / 'figures' / 'performance' / 'mdc_cluster_runtime.svg'

def _ebi_panel_height() -> float:
    """Figure height of the EBI per-step panel (S1a) so this companion panel
    (S1b) renders at the same height when laid out side by side. Mirrors the
    auto-height rule in render_per_step_boxplot, keyed on the EBI step count."""
    n_ebi = len(pd.read_csv(EBI_TSV, sep='\t'))
    return max(3.5, 0.45 * n_ebi + 1.5)

def _figure_mdc_cluster_runtime__render(out: Path) -> Path:
    df = pd.read_csv(_figure_mdc_cluster_runtime__DATA, sep='\t')
    df = df[df['status'] == 'COMPLETED']
    durations = {step: sub['duration_s'].tolist() for step, sub in df.groupby('step')}
    order = sorted(durations, key=lambda s: pd.Series(durations[s]).median(), reverse=True)
    summary = pd.DataFrame({'step': order})
    render_per_step_boxplot(durations, summary, out, fig_h=_ebi_panel_height())
    return out

def figure_mdc_cluster_runtime_main() -> int:
    print(f'wrote {_figure_mdc_cluster_runtime__render(_figure_mdc_cluster_runtime__OUT)}')
    return 0


# ======================================================================
# inlined from analysis/figure_original_vs_quantmsdiann.py
# ======================================================================

import csv
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
import requests
_figure_original_vs_quantmsdiann__REPO_ROOT = Path(__file__).resolve().parent.parent
_figure_original_vs_quantmsdiann__DATA_DIR = _figure_original_vs_quantmsdiann__REPO_ROOT / 'data' / 'PXD003539'
_figure_original_vs_quantmsdiann__FIGURES_DIR = _figure_original_vs_quantmsdiann__REPO_ROOT / 'analysis' / 'figures' / 'PXD003539'
_figure_original_vs_quantmsdiann__PRIDE_QUANT_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/cell-lines/PXD003539/v2_5_1_enterprise/quant_tables'
PRIDE_SUBMISSION_BASE = 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD003539'
PR_MATRIX_URL = f'{_figure_original_vs_quantmsdiann__PRIDE_QUANT_BASE}/diann_report.pr_matrix.tsv'
SUMMARY_LOG_URL = f'{_figure_original_vs_quantmsdiann__PRIDE_QUANT_BASE}/diannsummary.log'
DIANN_REPORT_PARQUET_URL = f'{_figure_original_vs_quantmsdiann__PRIDE_QUANT_BASE}/diann_report.parquet'
DIANN_UNIQUE_GENES_MATRIX_URL = f'{_figure_original_vs_quantmsdiann__PRIDE_QUANT_BASE}/diann_report.unique_genes_matrix.tsv'
OPENSWATH_MATRIX_URL = f'{PRIDE_SUBMISSION_BASE}/feature_alignment_requant_matrix.tsv'
HGNC_COMPLETE_SET_URL = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
_figure_original_vs_quantmsdiann__QUANTMS_SDRF_URL = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/sdrf/PXD003539.sdrf.tsv'
EA_EXPERIMENT_DESIGN_URL = 'https://www.ebi.ac.uk/gxa/experiments-content/E-PROT-73/resources/experiment-design'
EA_DOWNLOADS_URL = 'https://www.ebi.ac.uk/gxa/experiments/E-PROT-73/Downloads'
UNIQUE_GENES_METADATA_COLS = ['Genes', 'N.Sequences', 'N.Proteotypic.Sequences']
_figure_original_vs_quantmsdiann__PR_METADATA_COLS = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Precursor.Id']
WALZER_PEPTIDES = 77014
WALZER_PROTEINS = 7097
WALZER_PROTEINS_50PCT_FILTER = 6867
EPROT73_URL = 'https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-PROT-73/E-PROT-73.tsv'
GUO_CURATED_PEPTIDES = 22554
GUO_CURATED_PROTEINS = 3171
SUMMARY_LOG_PROTEIN_LINE_RE = re.compile('Protein groups with global q-value <= 0\\.01:\\s*(\\d+)')

@dataclass(frozen=True)
class _figure_original_vs_quantmsdiann__Counts:
    guo_peptides: int
    guo_proteins: int
    guo_precursors: int
    walzer_peptides: int
    walzer_proteins: int
    walzer_ea_genes: int
    diann_peptides: int
    diann_proteins: int
    diann_precursors: int
    diann_proteins_unfiltered: int = 0

def download_if_missing(url: str, dest: Path, *, retries: int=2) -> Path:
    """Download url to dest, skipping if dest already exists and is non-empty."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + '.part')
    last_exc: Exception | None = None
    try:
        for attempt in range(retries + 1):
            try:
                with requests.get(url, stream=True, timeout=120) as resp:
                    resp.raise_for_status()
                    with part.open('wb') as fh:
                        for chunk in resp.iter_content(chunk_size=1 << 20):
                            fh.write(chunk)
                os.replace(part, dest)
                return dest
            except requests.RequestException as exc:
                last_exc = exc
                if attempt < retries:
                    time.sleep(2 ** attempt)
        raise RuntimeError(f'Failed to download {url} after {retries + 1} attempts: {last_exc}')
    finally:
        part.unlink(missing_ok=True)

# When set, raw downloads (DIA-NN reports / matrices) are kept on disk for fast
# repeated runs. Default OFF: the rebuild streams the public-FTP data one file at
# a time and discards each after it has been counted/rendered, so it never needs
# the whole FTP dataset on disk at once. Only the small derived TSV/JSON outputs
# persist, and those are what the figure stages consume.
KEEP_DOWNLOADS = os.environ.get('REBUILD_KEEP_DOWNLOADS', '').lower() not in ('', '0', 'false', 'no')

def discard_download(*paths: Path) -> None:
    """Delete large downloaded reports/matrices once consumed (unless
    REBUILD_KEEP_DOWNLOADS is set), bounding peak disk to a single file."""
    if KEEP_DOWNLOADS:
        return
    for p in paths:
        try:
            Path(p).unlink(missing_ok=True)
        except OSError:
            pass

# Filename patterns of the large raw per-precursor FTP downloads (full DIA-NN
# reports, deposited zips, site reports) that are consumed once and re-derivable.
# NOTE: the cell-line `*_matrix.tsv` count matrices are deliberately NOT purged:
# they are small (~1.5 GB total) and shared across stages (a per-cohort figure
# downloads them; `atlas`/`venn` read the same files later), so deleting them
# between stages breaks those downstream figures.
_RAW_DOWNLOAD_GLOBS = ('diann_report.parquet', 'diann_report.tsv',
                       '*.zip', 'site_report*.parquet', 'report.parquet')

def purge_raw_downloads() -> int:
    """Delete the large raw FTP downloads under data/ and analysis/figures/,
    keeping the derived outputs. No-op when REBUILD_KEEP_DOWNLOADS is set. Called
    after each stage by the orchestrator so the rebuild never holds the whole
    public-FTP dataset on disk at once. Returns the number of files removed."""
    if KEEP_DOWNLOADS:
        return 0
    n = 0
    for root in (REPO / 'data', REPO / 'analysis' / 'figures'):
        if not root.exists():
            continue
        for pat in _RAW_DOWNLOAD_GLOBS:
            for p in root.rglob(pat):
                try:
                    p.unlink()
                    n += 1
                except OSError:
                    pass
    return n

def count_quantified_rows(matrix_path: Path, metadata_cols: list[str], unique_by: str | None=None) -> int:
    """Count quantified rows in a DIA-NN-style TSV matrix."""
    df = pd.read_csv(matrix_path, sep='\t', dtype=str)
    missing = [c for c in metadata_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Matrix is missing expected metadata column(s): {', '.join(missing)}")
    sample_cols = [c for c in df.columns if c not in metadata_cols]
    sample_df = df[sample_cols]
    quantified_mask = sample_df.notna().any(axis=1)
    quantified_df = df[quantified_mask]
    if unique_by is not None:
        return int(quantified_df[unique_by].nunique())
    return int(quantified_mask.sum())
PEPTIDE_ID_RE = re.compile('^(?:DECOY_)?\\d+_(?P<modseq>.+)_(?P<charge>\\d+)_run0$')
UNIMOD_RE = re.compile('\\(UniMod:\\d+\\)')

def _stripped_peptide(peptide_id: str) -> str | None:
    """Extract the unmodified peptide sequence from an OpenSWATH Peptide ID.

    Returns None if the ID does not match the expected format."""
    m = PEPTIDE_ID_RE.match(peptide_id)
    if not m:
        return None
    return UNIMOD_RE.sub('', m.group('modseq'))

def count_openswath_quantified(matrix_path: Path) -> tuple[int, int, int]:
    """Count quantified target precursors, peptides, and protein groups in an OpenSWATH matrix."""
    header_df = pd.read_csv(matrix_path, sep='\t', nrows=0)
    cols = list(header_df.columns)
    if 'Peptide' not in cols:
        raise ValueError('Matrix is missing required column: Peptide')
    if 'Protein' not in cols:
        raise ValueError('Matrix is missing required column: Protein')
    intensity_cols = [c for c in cols if c.startswith('Intensity_')]
    if not intensity_cols:
        raise ValueError('Matrix has zero Intensity_* columns; expected at least one')
    usecols = ['Peptide', 'Protein'] + intensity_cols
    total_precursors = 0
    peptides: set[str] = set()
    protein_groups: set[str] = set()
    for chunk in pd.read_csv(matrix_path, sep='\t', usecols=usecols, dtype=str, chunksize=50000):
        is_decoy = chunk['Peptide'].str.startswith('DECOY_', na=False) | chunk['Protein'].str.contains('DECOY', case=False, na=False)
        target = chunk[~is_decoy]
        quantified_mask = target[intensity_cols].notna().any(axis=1)
        quantified = target[quantified_mask]
        total_precursors += len(quantified)
        for pep_id in quantified['Peptide']:
            stripped = _stripped_peptide(pep_id)
            if stripped is not None:
                peptides.add(stripped)
        protein_groups.update(quantified['Protein'].tolist())
    return (total_precursors, len(peptides), len(protein_groups))

def count_eprot73_genes(tsv_path: Path) -> int:
    """Count unique Ensembl gene IDs in the Expression Atlas E-PROT-73 file.

    The file has no preamble: line 1 is the header (`Gene ID`, `Gene Name`, …)
    and data rows follow immediately from line 2.  We read all rows with
    `pd.read_csv` (no `skiprows`) and keep only those whose first column value
    starts with `ENSG`.  The header value `"Gene ID"` is naturally excluded by
    that prefix filter, so no extra row-dropping is needed."""
    df = pd.read_csv(tsv_path, sep='\t', dtype=str, usecols=[0])
    df.columns = ['gene_id']
    df = df[df['gene_id'].fillna('').str.startswith('ENSG')]
    return int(df['gene_id'].nunique())

def per_run_real_detection_fraction_diann_parquet(parquet_path: Path, *, qvalue_cutoff: float=0.01, global_qvalue_cutoff: float=0.01) -> dict[str, float]:
    """Per-run fraction of distinct precursors confidently identified, read
    from the DIA-NN long-format report (`diann_report.parquet`).

    Per methods.md §1: the per-run numerator uses the per-run precursor rule
    (`Q.Value <= qvalue_cutoff`) ONLY; the denominator is the global precursor
    pool (`Lib.Q.Value <= global_qvalue_cutoff` anywhere in the report). No
    `Global.Q.Value` gate, no contaminant filter, zeros counted. The global
    pool is the apples-to-apples analogue of OpenSWATH's "all target rows in
    the requant matrix" denominator.

    Why use the parquet instead of the pr_matrix? The matrix is filtered at the
    cell level by --matrix-spec-q 0.05 (spectrum-level quant FDR), not by the
    per-run identification Q.Value. To match OpenSWATH's `score <= 0.01` per-run
    criterion strictly, we have to read Q.Value from the long-format report."""
    import pyarrow.parquet as pq
    cols = ['Run', 'Precursor.Id', 'Q.Value', 'Lib.Q.Value']
    pf = pq.ParquetFile(str(parquet_path))
    schema_names = pf.schema_arrow.names
    missing = [c for c in cols if c not in schema_names]
    if missing:
        raise ValueError(f'DIA-NN parquet missing expected columns: {missing}')
    global_precursors: set[str] = set()
    per_run_detected: dict[str, set[str]] = {}
    for batch in pf.iter_batches(columns=cols, batch_size=200000):
        runs = batch.column('Run').to_pylist()
        pids = batch.column('Precursor.Id').to_pylist()
        qvs = batch.column('Q.Value').to_pylist()
        lqvs = batch.column('Lib.Q.Value').to_pylist()
        for r, p, q, lq in zip(runs, pids, qvs, lqvs):
            per_run_detected.setdefault(r, set())
            if lq is not None and lq <= global_qvalue_cutoff:
                global_precursors.add(p)
            if q is None or q > qvalue_cutoff:
                continue
            per_run_detected[r].add(p)
    denom = len(global_precursors)
    if denom == 0:
        return {r: 0.0 for r in per_run_detected}
    return {r: len(s & global_precursors) / denom for r, s in per_run_detected.items()}

def per_run_real_detection_fraction_openswath(matrix_path: Path, *, qvalue_cutoff: float=0.01) -> dict[str, float]:
    """For each `score_<run>` column in the OpenSWATH requantification matrix,
    return the fraction of target (non-decoy) rows whose score for that run
    is below `qvalue_cutoff`.

    The OpenSWATH `feature_alignment_requant_matrix.tsv` is a *requantified*
    matrix: Intensity values are filled in across all runs even when a
    precursor was not directly identified in that run. The corresponding
    `score_<run>` column (pyprophet m_score / q-value) is the truth signal —
    real detections have score <= 0.01, requantified placeholders are tagged
    score = 2.0. Counting non-NA Intensity values is therefore not a real
    measure of completeness; we use the score column instead."""
    header_df = pd.read_csv(matrix_path, sep='\t', nrows=0)
    cols = list(header_df.columns)
    if 'Peptide' not in cols or 'Protein' not in cols:
        raise ValueError('OpenSWATH matrix missing Peptide/Protein columns')
    score_cols = [c for c in cols if c.startswith('score_')]
    if not score_cols:
        raise ValueError('OpenSWATH matrix has no score_* columns')
    target_total = 0
    per_run_hits: dict[str, int] = {c: 0 for c in score_cols}
    usecols = ['Peptide', 'Protein'] + score_cols
    for chunk in pd.read_csv(matrix_path, sep='\t', dtype=str, usecols=usecols, chunksize=50000):
        is_decoy = chunk['Peptide'].str.startswith('DECOY_', na=False) | chunk['Protein'].str.contains('DECOY', case=False, na=False)
        targets = chunk[~is_decoy].copy()
        target_total += len(targets)
        for col in score_cols:
            scores = pd.to_numeric(targets[col], errors='coerce')
            per_run_hits[col] += int((scores <= qvalue_cutoff).sum())
    if target_total == 0:
        return {c: 0.0 for c in score_cols}
    return {c: per_run_hits[c] / target_total for c in score_cols}

def unique_peptides_per_protein_diann(matrix_path: Path) -> dict[str, int]:
    """Per Protein.Group, count distinct Stripped.Sequence values from
    proteotypic precursor rows that have at least one non-NA sample value.

    Restricting to Proteotypic == 1 ensures we only count peptides that uniquely
    identify the protein (the natural definition of 'unique peptides per
    protein'). Multiple charge states / modforms of the same peptide collapse
    to a single Stripped.Sequence entry.

    Per methods.md §1 there is NO contaminant/target filter: every quantified
    proteotypic row counts (the matrix carries no q-value column to re-apply).
    The >=k-unique-peptides threshold is a study-defined comparison criterion,
    not one of the forbidden q/contaminant filters, and is kept."""
    header = pd.read_csv(matrix_path, sep='\t', nrows=0)
    missing = [c for c in _figure_original_vs_quantmsdiann__PR_METADATA_COLS if c not in header.columns]
    if missing:
        raise ValueError(f'DIA-NN matrix missing metadata columns: {missing}')
    sample_cols = [c for c in header.columns if c not in _figure_original_vs_quantmsdiann__PR_METADATA_COLS]
    pg_to_peps: dict[str, set[str]] = {}
    for chunk in pd.read_csv(matrix_path, sep='\t', dtype=str, chunksize=100000):
        proteotypic = chunk[chunk['Proteotypic'] == '1']
        if proteotypic.empty:
            continue
        quantified = proteotypic[proteotypic[sample_cols].notna().any(axis=1)]
        for pg, seq in zip(quantified['Protein.Group'], quantified['Stripped.Sequence']):
            pg_to_peps.setdefault(pg, set()).add(seq)
    return {pg: len(s) for pg, s in pg_to_peps.items()}

def unique_peptides_per_protein_openswath(matrix_path: Path, *, qvalue_cutoff: float=0.01) -> dict[str, int]:
    """Per OpenSWATH Protein, count distinct stripped peptide sequences that
    were confidently detected (score <= `qvalue_cutoff`) in at least one run.

    Restricted to proteotypic rows — those whose Protein column starts with
    `1/` (peptide maps to exactly one protein). Decoys are excluded. The
    confidence filter prevents requant-only placeholder rows (where every score
    is 2.0) from inflating the per-protein peptide count."""
    header_df = pd.read_csv(matrix_path, sep='\t', nrows=0)
    cols = list(header_df.columns)
    if 'Peptide' not in cols or 'Protein' not in cols:
        raise ValueError('OpenSWATH matrix missing Peptide/Protein columns')
    score_cols = [c for c in cols if c.startswith('score_')]
    if not score_cols:
        raise ValueError('OpenSWATH matrix has no score_* columns')
    usecols = ['Peptide', 'Protein'] + score_cols
    protein_to_peptides: dict[str, set[str]] = {}
    for chunk in pd.read_csv(matrix_path, sep='\t', dtype=str, usecols=usecols, chunksize=50000):
        is_decoy = chunk['Peptide'].str.startswith('DECOY_', na=False) | chunk['Protein'].str.contains('DECOY', case=False, na=False)
        targets = chunk[~is_decoy]
        proteotypic = targets[targets['Protein'].str.startswith('1/', na=False)]
        if proteotypic.empty:
            continue
        score_block = proteotypic[score_cols].apply(pd.to_numeric, errors='coerce')
        any_detection = (score_block <= qvalue_cutoff).any(axis=1)
        detected = proteotypic[any_detection]
        for prot, pep_id in zip(detected['Protein'], detected['Peptide']):
            stripped = _stripped_peptide(pep_id)
            if stripped is None:
                continue
            protein_to_peptides.setdefault(prot, set()).add(stripped)
    return {p: len(s) for p, s in protein_to_peptides.items()}

def proteins_with_min_peptides(counts: dict[str, int], min_k: int) -> int:
    """Number of proteins with at least `min_k` unique peptides."""
    return sum((1 for n in counts.values() if n >= min_k))

def count_quantified_genes_diann(matrix_path: Path) -> int:
    """Count rows in DIA-NN's unique_genes_matrix.tsv with >=1 non-NA sample
    cell. Metadata columns are UNIQUE_GENES_METADATA_COLS (Genes, N.Sequences,
    N.Proteotypic.Sequences); the remaining columns are per-run intensities."""
    df = pd.read_csv(matrix_path, sep='\t', dtype=str)
    missing = [c for c in UNIQUE_GENES_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'unique_genes_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in UNIQUE_GENES_METADATA_COLS]
    if not sample_cols:
        raise ValueError('unique_genes_matrix has no sample columns')
    return int(df[sample_cols].notna().any(axis=1).sum())

def load_hgnc_symbol_to_ensembl(hgnc_tsv_path: Path) -> dict[str, str]:
    """Parse HGNC's `hgnc_complete_set.txt` and return a mapping from any known
    gene-symbol form (current `symbol`, pipe-separated `alias_symbol` and
    `prev_symbol`) to its `ensembl_gene_id`. Rows without an ensembl_gene_id
    are skipped. Conflicts are resolved last-writer-wins (which only matters
    when a single symbol historically pointed to multiple genes — rare)."""
    df = pd.read_csv(hgnc_tsv_path, sep='\t', dtype=str, usecols=['symbol', 'alias_symbol', 'prev_symbol', 'ensembl_gene_id'])
    df = df[df['ensembl_gene_id'].fillna('').str.startswith('ENSG')]
    mapping: dict[str, str] = {}
    for sym, aliases, prevs, ensg in zip(df['symbol'], df['alias_symbol'], df['prev_symbol'], df['ensembl_gene_id']):
        if isinstance(sym, str) and sym:
            mapping[sym] = ensg
        if isinstance(aliases, str) and aliases:
            for a in aliases.split('|'):
                a = a.strip()
                if a:
                    mapping[a] = ensg
        if isinstance(prevs, str) and prevs:
            for p in prevs.split('|'):
                p = p.strip()
                if p:
                    mapping[p] = ensg
    return mapping

def load_walzer_genes_ensembl(eprot73_path: Path) -> set[str]:
    """Return the set of unique Ensembl gene IDs from the Walzer E-PROT-73
    Expression Atlas file. Same filter as count_eprot73_genes (rows whose
    `Gene ID` column starts with `ENSG`)."""
    df = pd.read_csv(eprot73_path, sep='\t', dtype=str, usecols=[0])
    df.columns = ['gene_id']
    df = df[df['gene_id'].fillna('').str.startswith('ENSG')]
    return set(df['gene_id'].unique())

def quantmsdiann_genes_as_ensembl(matrix_path: Path, symbol_to_ensembl: dict[str, str], *, min_detection_fraction: float=0.0) -> tuple[set[str], int]:
    """For each quantified row in the DIA-NN unique_genes_matrix, look up the
    `Genes` column in `symbol_to_ensembl` and collect the Ensembl gene IDs.
    Returns (mapped_ensg_set, unmapped_count).

    `min_detection_fraction` controls how many sample cells must be non-NA for
    a row to count as 'detected':
      - 0.0 (default): any row with >=1 non-NA cell (the bare identification set).
      - 0.5: row must be non-NA in >=ceil(0.5 * n_samples) runs (mimics Walzer's
        50%-per-group consistency filter applied globally across all runs).

    The Genes column can hold multiple symbols separated by `;` (gene-group
    case). We split on `;` and try each; counts as unmapped if none of the
    components are in the mapping."""
    import math
    df = pd.read_csv(matrix_path, sep='\t', dtype=str)
    missing = [c for c in UNIQUE_GENES_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'unique_genes_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in UNIQUE_GENES_METADATA_COLS]
    non_na_count = df[sample_cols].notna().sum(axis=1)
    min_required = max(1, math.ceil(min_detection_fraction * len(sample_cols)))
    detected = df[non_na_count >= min_required]
    mapped: set[str] = set()
    unmapped = 0
    for raw in detected['Genes'].fillna(''):
        if not raw:
            unmapped += 1
            continue
        ensgs = {symbol_to_ensembl[s] for s in (sym.strip() for sym in raw.split(';')) if s and s in symbol_to_ensembl}
        if ensgs:
            mapped.update(ensgs)
        else:
            unmapped += 1
    return (mapped, unmapped)
_CELL_LINE_NCI_PREFIX_RE = re.compile('^NCI[-_/]', flags=re.IGNORECASE)
_CELL_LINE_NONALNUM_RE = re.compile('[^A-Za-z0-9]')

def normalize_cell_line(name: str | None) -> str:
    """Normalise a cell-line name so the quantms SDRF and E-PROT-73 spellings
    collide. Strips a leading 'NCI-' / 'NCI_' / 'NCI/' prefix, removes all
    other non-alphanumeric characters, and uppercases.

    Example: 'CCRF-CEM' -> 'CCRFCEM', 'NCI-H226' -> 'H226',
    'Hs-578-T' -> 'HS578T'."""
    if not name:
        return ''
    name = _CELL_LINE_NCI_PREFIX_RE.sub('', name)
    name = _CELL_LINE_NONALNUM_RE.sub('', name)
    return name.upper()

def parse_eprot73_groupings(downloads_html_path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Parse the Expression Atlas E-PROT-73 Downloads HTML page and return
    (g_to_cell_line, g_to_disease) mappings.

    The Downloads page embeds an inline JSON object (`content: {...}`) with
    `tabs[0].props.groups` listing primary groupings. The CELL_LINE grouping
    is `[[cell_line_name, [g1, g2, ...]], ...]`; DISEASE follows the same
    shape. We pull both out with a small balanced-bracket scanner so we don't
    depend on the surrounding HTML structure."""
    import json
    html = downloads_html_path.read_text()

    def extract(name: str) -> list[list]:
        m = re.search(f'"name":"{re.escape(name)}",[^{{]*?"groupings":', html)
        if m is None:
            raise ValueError(f"E-PROT-73 Downloads HTML missing '{name}' grouping")
        start = m.end()
        if html[start] != '[':
            raise ValueError(f"Unexpected character after '{name}' groupings: {html[start]!r}")
        depth = 0
        for i in range(start, len(html)):
            if html[i] == '[':
                depth += 1
            elif html[i] == ']':
                depth -= 1
                if depth == 0:
                    return json.loads(html[start:i + 1])
        raise ValueError(f"Unbalanced brackets in '{name}' grouping")
    cl_groupings = extract('CELL_LINE')
    ds_groupings = extract('DISEASE')
    g_to_cl: dict[str, str] = {}
    for name, gs in cl_groupings:
        for g in gs:
            g_to_cl[g] = name
    g_to_ds: dict[str, str] = {}
    for name, gs in ds_groupings:
        for g in gs:
            g_to_ds[g] = name
    return (g_to_cl, g_to_ds)

def _figure_original_vs_quantmsdiann__load_sdrf_data_file_to_cell_line(sdrf_path: Path) -> dict[str, str]:
    """Parse the PXD003539 quantms-collections SDRF and return a mapping from
    DIA-NN matrix column name (i.e. the `comment[data file]` field rewritten
    to a `.mzML` extension) to the row's `characteristics[cell line]`.

    The SDRF stores `*.wiff` filenames but DIA-NN's pr_matrix /
    unique_genes_matrix columns are `*.mzML`. We rewrite the extension so
    downstream callers can look up DIA-NN columns directly."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'comment[data file]']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f'SDRF missing required columns: {missing}')
    out: dict[str, str] = {}
    for cell_line, data_file in zip(df['characteristics[cell line]'], df['comment[data file]']):
        if not isinstance(data_file, str) or not data_file:
            continue
        mzml = re.sub('\\.wiff$', '.mzML', data_file)
        out[mzml] = cell_line
    return out

def load_ea_cell_line_to_disease(experiment_design_path: Path) -> dict[str, str]:
    """Parse the E-PROT-73 experiment-design TSV and return a mapping from
    normalised cell-line name (see `normalize_cell_line`) to disease label.

    The TSV has one row per MS Run; every run for a given cell line lists the
    same disease, so we collapse to the unique mapping. Conflicting diseases
    for the same cell line would be unexpected and we raise on collision."""
    df = pd.read_csv(experiment_design_path, sep='\t', dtype=str)
    cell_col = 'Sample Characteristic[cell line]'
    disease_col = 'Sample Characteristic[disease]'
    missing = [c for c in (cell_col, disease_col) if c not in df.columns]
    if missing:
        raise ValueError(f'E-PROT-73 experiment-design missing columns: {missing}')
    out: dict[str, str] = {}
    for cell, disease in zip(df[cell_col], df[disease_col]):
        if not isinstance(cell, str) or not isinstance(disease, str):
            continue
        key = normalize_cell_line(cell)
        if not key:
            continue
        if key in out and out[key] != disease:
            raise ValueError(f'Conflicting disease for cell line {cell!r}: {out[key]!r} vs {disease!r}')
        out[key] = disease
    return out

def walzer_genes_per_condition(eprot73_path: Path, downloads_html_path: Path) -> dict[str, set[str]]:
    """For each E-PROT-73 disease, return the set of Ensembl gene IDs detected
    in at least one g<N> column belonging to that disease.

    "Detected" means abundance > 0 in the `g<N>.WithInSampleAbundance` column.
    Rows whose `Gene ID` doesn't start with `ENSG` are ignored (Expression
    Atlas occasionally emits non-ENSG summary rows like 'totalGenes')."""
    _, g_to_disease = parse_eprot73_groupings(downloads_html_path)
    df = pd.read_csv(eprot73_path, sep='\t', dtype=str)
    if df.columns[0] != 'Gene ID':
        raise ValueError(f"E-PROT-73 TSV first column should be 'Gene ID', got {df.columns[0]!r}")
    df = df[df['Gene ID'].fillna('').str.startswith('ENSG')].copy()
    sample_cols = [c for c in df.columns if '.WithInSampleAbundance' in c]
    col_to_g = {c: c.split('.', 1)[0] for c in sample_cols}
    out: dict[str, set[str]] = {}
    for col in sample_cols:
        g = col_to_g[col]
        disease = g_to_disease.get(g)
        if disease is None:
            continue
        vals = pd.to_numeric(df[col], errors='coerce').fillna(0.0)
        detected = df.loc[vals > 0, 'Gene ID']
        out.setdefault(disease, set()).update(detected)
    return out

def quantmsdiann_genes_per_condition(unique_genes_matrix_path: Path, sdrf_path: Path, ea_design_path: Path, symbol_to_ensembl: dict[str, str], *, min_detection_fraction_per_cell_line: float=0.0, min_global_detection_fraction: float=0.0) -> dict[str, set[str]]:
    """For each disease (using the E-PROT-73 9-cancer-type axis), return the
    set of Ensembl gene IDs detected by quantmsdiann in at least one cell
    line belonging to that disease.

    Filtering knobs (apply jointly — a gene must pass both):
      - `min_global_detection_fraction`: gene must be non-NA in
        >= ceil(f * n_runs_total) of all 120 DIA-NN runs. This is the same
        global filter used in `quantmsdiann_genes_as_ensembl` and the
        `supp_walzer_vs_quantms_genes_ensembl` Venn (set to 0.5 there).
      - `min_detection_fraction_per_cell_line`: gene must be non-NA in
        >= ceil(f * n_runs_of_cell_line) of a given cell line's runs to
        count as "detected for that cell line". This mirrors Walzer's
        '50% per group' consistency filter at the replicate-group level
        (Walzer 2022 §Methods), which is what produces the
        g<N>.WithInSampleAbundance values in E-PROT-73. With only 2 runs
        per cell line in this dataset, 0.5 here is essentially a no-op
        (ceil(0.5*2)=1) — the user-facing filter is the global one above.

    The SDRF maps DIA-NN columns to cell line (dashed quantms spelling);
    `load_ea_cell_line_to_disease` maps the normalised cell-line name to one
    of the 9 E-PROT-73 disease labels. Genes column entries are split on `;`
    and each component looked up in the HGNC symbol -> Ensembl mapping;
    unmappable rows are skipped (consistent with `quantmsdiann_genes_as_ensembl`)."""
    import math
    df = pd.read_csv(unique_genes_matrix_path, sep='\t', dtype=str)
    missing = [c for c in UNIQUE_GENES_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'unique_genes_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in UNIQUE_GENES_METADATA_COLS]
    sdrf_run_to_cell = _figure_original_vs_quantmsdiann__load_sdrf_data_file_to_cell_line(sdrf_path)
    ea_cell_to_disease = load_ea_cell_line_to_disease(ea_design_path)
    cell_to_cols: dict[str, list[str]] = {}
    cell_to_disease: dict[str, str] = {}
    for col in sample_cols:
        cell = sdrf_run_to_cell.get(col)
        if not cell:
            continue
        disease = ea_cell_to_disease.get(normalize_cell_line(cell))
        if disease is None:
            continue
        cell_to_cols.setdefault(cell, []).append(col)
        cell_to_disease[cell] = disease
    row_ensgs: list[set[str]] = []
    for raw in df['Genes'].fillna(''):
        if not raw:
            row_ensgs.append(set())
            continue
        ensgs = {symbol_to_ensembl[s] for s in (sym.strip() for sym in raw.split(';')) if s and s in symbol_to_ensembl}
        row_ensgs.append(ensgs)
    global_min = max(1, math.ceil(min_global_detection_fraction * len(sample_cols)))
    global_non_na = df[sample_cols].notna().sum(axis=1)
    global_pass = global_non_na >= global_min
    out: dict[str, set[str]] = {}
    for cell, cols in cell_to_cols.items():
        non_na = df[cols].notna().sum(axis=1)
        min_required = max(1, math.ceil(min_detection_fraction_per_cell_line * len(cols)))
        detected_mask = (non_na >= min_required) & global_pass
        disease = cell_to_disease[cell]
        bucket = out.setdefault(disease, set())
        for ensgs, ok in zip(row_ensgs, detected_mask):
            if ok:
                bucket.update(ensgs)
    return out

def parse_summary_log(log_path: Path) -> int:
    """Return protein group count at 1% global FDR from a DIA-NN summary log."""
    with open(log_path, encoding='utf-8') as fh:
        for line in fh:
            m = SUMMARY_LOG_PROTEIN_LINE_RE.search(line)
            if m:
                return int(m.group(1))
    raise ValueError("Line matching 'Protein groups with global q-value <= 0.01:' not found in log")

def report_global_counts_diann(parquet_path: Path) -> dict[str, int]:
    """Global (dataset-total) counts from the DIA-NN long-format report under
    the methods.md §1 rule: protein groups on ``Lib.PG.Q.Value <= 0.01`` only,
    precursors on ``Lib.Q.Value <= 0.01`` only, and ``>= 2`` distinct stripped
    peptides per global protein group. No contaminant/target filter, no
    positive-quantity filter; decoys (``Decoy == 1``) are dropped.

    Reads only the columns ``count_report`` needs and returns its ``prot_global``
    (global protein groups), ``prec_global`` (global precursors), ``peptides``
    (distinct stripped peptides), and ``prot_2pep`` (>= 2-peptide protein
    groups) keys."""
    import pyarrow.parquet as pq
    have = set(pq.ParquetFile(parquet_path).schema_arrow.names)
    cols = [c for c in _NEEDED_COLS if c in have]
    df = pq.read_table(parquet_path, columns=cols).to_pandas()
    return count_report(df, precursor_q=DEFAULT_PRECURSOR_Q)

def render_figure(counts: Counts, svg_path: Path) -> None:
    """Render a 2-condition x 2-metric grouped bar chart. Paper-ready: only
    bars, value labels, axis labels, and legend — no title, no footer."""
    conditions = [('Guo 2019\n(OpenSWATH)', '#9e9e9e', counts.guo_peptides, counts.guo_proteins), ('quantmsdiann\n(DIA-NN)', '#1f77b4', counts.diann_peptides, counts.diann_proteins)]
    metrics = ['Peptides', 'Protein groups']
    bar_width = 0.27
    n_conditions = len(conditions)
    x = [0, 1]
    fig, ax = plt.subplots(figsize=(7, 5))
    offsets = [bar_width * (i - (n_conditions - 1) / 2.0) for i in range(n_conditions)]
    for i, (label, color, peptide_val, protein_val) in enumerate(conditions):
        values = [peptide_val, protein_val]
        bars = ax.bar([xi + offsets[i] for xi in x], values, width=bar_width, color=color, label=label)
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height(), f'{val:,}', ha='center', va='bottom', fontsize=9)
    peptide_vals = [c[2] for c in conditions]
    protein_vals = [c[3] for c in conditions]
    needs_log = False
    for metric_vals in [peptide_vals, protein_vals]:
        mn, mx = (min(metric_vals), max(metric_vals))
        if mn > 0 and mx / mn > 5:
            needs_log = True
            break
    ylabel = 'Count (1% FDR)'
    if needs_log:
        ax.set_yscale('log')
        ylabel += ' (log scale)'
    else:
        top = max(max(peptide_vals), max(protein_vals)) * 1.18
        ax.set_ylim(0, top)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def render_missing_values_per_run(guo_per_run: dict[str, float], diann_per_run: dict[str, float], svg_path: Path) -> None:
    """Plot per-run non-NA fraction for both pipelines. Paper-ready: no title,
    no footer; only lines, axis labels, and legend.

    Runs on the x axis are aligned by the L<date>_<n>_SW stem extracted from
    the column name. Bars sorted by the stem (which sorts chronologically by
    acquisition date)."""
    import re as _re

    def stem(col: str) -> str:
        m = _re.search('L\\d+_\\d+_SW', col)
        return m.group(0) if m else col
    diann_by_stem = {stem(c): v for c, v in diann_per_run.items()}
    guo_by_stem = {stem(c): v for c, v in guo_per_run.items()}
    common = sorted(set(diann_by_stem) & set(guo_by_stem))
    fig, ax = plt.subplots(figsize=(10, 4))
    x = range(len(common))
    ax.plot(list(x), [guo_by_stem[s] for s in common], label='Guo 2019 (OpenSWATH)', color=fs.COMPARISON['original'], linewidth=1.0)
    ax.plot(list(x), [diann_by_stem[s] for s in common], label='quantmsdiann (DIA-NN)', color=fs.COMPARISON['quantmsdiann'], linewidth=1.0)
    ax.set_xlabel(f'MS run index ({len(common)} runs, ordered by acquisition date)')
    ax.set_ylabel('Fraction of precursors quantified per run')
    ax.set_ylim(0, 1.05)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='lower right', frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)
_DISEASE_LABEL_ORDER = ['leukemia', 'central nervous system cancer', 'breast cancer', 'colorectal cancer', 'lung cancer', 'melanoma', 'ovarian cancer', 'prostate cancer', 'renal cancer']

def render_genes_per_condition(walzer_per_cond: dict[str, set[str]], diann_per_cond: dict[str, set[str]], svg_path: Path) -> None:
    """Grouped bar chart: per disease condition, two bars (Walzer 2022
    E-PROT-73 vs quantmsdiann) showing the number of distinct Ensembl gene IDs
    detected. The disease axis is the 9 NCI-60 cancer types defined by
    E-PROT-73's primary DISEASE grouping; both pipelines use this same axis
    after cell-line normalisation. Paper-ready: no title, no footer."""
    conditions = sorted(set(walzer_per_cond) | set(diann_per_cond), key=lambda d: _DISEASE_LABEL_ORDER.index(d) if d in _DISEASE_LABEL_ORDER else len(_DISEASE_LABEL_ORDER))
    walzer_vals = [len(walzer_per_cond.get(c, set())) for c in conditions]
    diann_vals = [len(diann_per_cond.get(c, set())) for c in conditions]
    fig, ax = plt.subplots(figsize=(11, 6))
    x = list(range(len(conditions)))
    bar_width = 0.4
    bars_w = ax.bar([xi - bar_width / 2 for xi in x], walzer_vals, width=bar_width, color='#90caf9', label='Walzer 2022 (E-PROT-73)')
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals, width=bar_width, color=fs.COMPARISON['quantmsdiann'], label='quantmsdiann (DIA-NN)')
    for bars, vals in ((bars_w, walzer_vals), (bars_d, diann_vals)):
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{v:,}', ha='center', va='bottom', fontsize=8)
    ax.set_xticks(x)
    display = [c.replace('central nervous system cancer', 'CNS cancer') for c in conditions]
    ax.set_xticklabels(display, rotation=25, ha='right')
    ax.set_ylabel('Distinct Ensembl gene IDs detected')
    ymax = max(max(walzer_vals, default=0), max(diann_vals, default=0))
    ax.set_ylim(0, ymax * 1.15 if ymax else 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def render_peptides_per_protein(guo_peptide_counts: dict[str, int], diann_peptide_counts: dict[str, int], svg_path: Path, thresholds: tuple[int, ...]=(1, 2, 3, 5, 10)) -> None:
    """Grouped bar chart: number of protein groups with >=k unique peptides,
    for each k in `thresholds`, comparing Guo (OpenSWATH) vs quantmsdiann
    (DIA-NN). Paper-ready: no title, no footer."""
    guo_values = [proteins_with_min_peptides(guo_peptide_counts, k) for k in thresholds]
    diann_values = [proteins_with_min_peptides(diann_peptide_counts, k) for k in thresholds]
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.38
    x = list(range(len(thresholds)))
    bars_guo = ax.bar([xi - width / 2 for xi in x], guo_values, width=width, color=fs.COMPARISON['original'], label='Guo 2019 (OpenSWATH)')
    bars_diann = ax.bar([xi + width / 2 for xi in x], diann_values, width=width, color=fs.COMPARISON['quantmsdiann'], label='quantmsdiann (DIA-NN)')
    for bars, vals in [(bars_guo, guo_values), (bars_diann, diann_values)]:
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height(), f'{v:,}', ha='center', va='bottom', fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([f'≥ {k}' for k in thresholds])
    ax.set_xlabel('Minimum unique peptides per protein group')
    ax.set_ylabel('Protein groups (1% FDR)')
    ymax = max(max(guo_values), max(diann_values)) * 1.14
    ax.set_ylim(0, ymax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def _figure_original_vs_quantmsdiann__write_counts_tsv(counts: Counts, tsv_path: Path) -> None:
    """Write an auditable TSV with metric/source/count/note rows."""
    rows = [('Peptides', 'Guo 2019 (OpenSWATH, deposited)', counts.guo_peptides, 'from feature_alignment_requant_matrix.tsv (target rows, >=1 quant)'), ('Peptides', 'Walzer 2022 (CAL + OpenSWATH, top3)', counts.walzer_peptides, 'Supplementary Table S2 (PXD003539, 1% FDR, top3)'), ('Peptides', 'quantmsdiann (DIA-NN, 1% FDR)', counts.diann_peptides, 'unique Stripped.Sequence in pr_matrix with >=1 non-NA'), ('Protein groups', 'Guo 2019 (OpenSWATH, deposited)', counts.guo_proteins, 'unique Protein column values in feature_alignment_requant_matrix.tsv'), ('Protein groups', 'Walzer 2022 (CAL + OpenSWATH, top3)', counts.walzer_proteins, 'Supplementary Table S2 (PXD003539, 1% FDR, top3, unfiltered)'), ('Protein groups', 'quantmsdiann (DIA-NN, Lib.PG.Q.Value)', counts.diann_proteins, 'global rule: distinct Protein.Group at Lib.PG.Q.Value<=0.01 in diann_report.parquet per methods.md §1; no contaminant/target filter, no Global.PG.Q.Value gate'), ('Protein groups', 'quantmsdiann (DIA-NN, diannsummary.log)', counts.diann_proteins_unfiltered, "audit baseline: diannsummary.log 'Protein groups with global q-value <= 0.01' line"), ('Precursors aux', 'Guo 2019 (OpenSWATH, deposited)', counts.guo_precursors, 'target rows with >=1 non-NA Intensity in feature_alignment matrix'), ('Precursors aux', 'quantmsdiann (DIA-NN, 1% FDR)', counts.diann_precursors, 'rows in pr_matrix with >=1 non-NA'), ('Curated context', 'Guo 2019 (DIA-expert curated, peptides)', 22554, 'paper text, not used as a headline bar'), ('Curated context', 'Guo 2019 (DIA-expert curated, proteins)', 3171, 'paper text, not used as a headline bar'), ('Filter context', 'Walzer 2022 (50% per group filter)', 6867, "Supplementary Table S2 - '50% per group' consistency filter, proteins"), ('EA context', 'Walzer 2022 (E-PROT-73, Expression Atlas)', counts.walzer_ea_genes, 'unique Ensembl gene IDs in E-PROT-73.tsv (post-processed: gene mapping + 50% per-group filter + decoy removal)')]
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['metric', 'source', 'count', 'note'])
        for row in rows:
            writer.writerow(row)

def figure_original_vs_quantmsdiann_main() -> int:
    ensure_cell_line_matrices('PXD003539', with_report=True)
    _figure_original_vs_quantmsdiann__DATA_DIR.mkdir(parents=True, exist_ok=True)
    _figure_original_vs_quantmsdiann__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    pr_path = _figure_original_vs_quantmsdiann__DATA_DIR / 'diann_report.pr_matrix.tsv'
    log_path = _figure_original_vs_quantmsdiann__DATA_DIR / 'diannsummary.log'
    opensw_path = _figure_original_vs_quantmsdiann__DATA_DIR / 'feature_alignment_requant_matrix.tsv'
    parquet_path = _figure_original_vs_quantmsdiann__DATA_DIR / 'diann_report.parquet'
    diann_genes_path = _figure_original_vs_quantmsdiann__DATA_DIR / 'diann_report.unique_genes_matrix.tsv'
    print('Resolving inputs (filesystem; external baselines cached)...')
    download_if_missing(OPENSWATH_MATRIX_URL, opensw_path)
    eprot73_path = download_if_missing(EPROT73_URL, _figure_original_vs_quantmsdiann__DATA_DIR / 'E-PROT-73.tsv')
    print('Computing quantmsdiann counts...')
    diann_peptides = count_quantified_rows(pr_path, _figure_original_vs_quantmsdiann__PR_METADATA_COLS, unique_by='Stripped.Sequence')
    diann_precursors = count_quantified_rows(pr_path, _figure_original_vs_quantmsdiann__PR_METADATA_COLS)
    diann_proteins_log = parse_summary_log(log_path)
    report_counts = report_global_counts_diann(parquet_path)
    diann_proteins = int(report_counts['prot_global'])
    diann_prot_2pep = int(report_counts['prot_2pep'])
    print(f'  report protein groups (Lib.PG.Q.Value<=0.01): {diann_proteins:,}  (>=2 peptides: {diann_prot_2pep:,})')
    print('Computing Guo 2019 (OpenSWATH) counts...')
    guo_precursors, guo_peptides, guo_proteins = count_openswath_quantified(opensw_path)
    print('Computing Expression Atlas E-PROT-73 gene count...')
    walzer_ea_genes = count_eprot73_genes(eprot73_path)
    counts = _figure_original_vs_quantmsdiann__Counts(guo_peptides=guo_peptides, guo_proteins=guo_proteins, guo_precursors=guo_precursors, walzer_peptides=WALZER_PEPTIDES, walzer_proteins=WALZER_PROTEINS, walzer_ea_genes=walzer_ea_genes, diann_peptides=diann_peptides, diann_proteins=diann_proteins, diann_precursors=diann_precursors, diann_proteins_unfiltered=diann_proteins_log)
    svg_path = _figure_original_vs_quantmsdiann__FIGURES_DIR / 'main_comparison.svg'
    render_figure(counts, svg_path)
    print(f'Figure saved to {svg_path}')
    data_dir = _figure_original_vs_quantmsdiann__FIGURES_DIR / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = data_dir / 'counts.tsv'
    _figure_original_vs_quantmsdiann__write_counts_tsv(counts, tsv_path)
    print(f'Counts TSV saved to {tsv_path}')
    print('Computing per-run completeness for supplementary figure B...')
    guo_per_run = per_run_real_detection_fraction_openswath(opensw_path)
    diann_per_run = per_run_real_detection_fraction_diann_parquet(parquet_path)
    suppB_svg = _figure_original_vs_quantmsdiann__FIGURES_DIR / 'supp_missing_values_per_run.svg'
    render_missing_values_per_run(guo_per_run, diann_per_run, suppB_svg)
    print(f'Supplementary B figure saved to {suppB_svg}')
    import statistics as _stats
    print(f'Median per-run real-detection fraction (1% per-run FDR): Guo={_stats.median(guo_per_run.values()):.2%} quantmsdiann={_stats.median(diann_per_run.values()):.2%}')
    print('Computing unique peptides per protein for supplementary figure C...')
    guo_pep_per_prot = unique_peptides_per_protein_openswath(opensw_path)
    diann_pep_per_prot = unique_peptides_per_protein_diann(pr_path)
    suppC_svg = _figure_original_vs_quantmsdiann__FIGURES_DIR / 'supp_peptides_per_protein.svg'
    render_peptides_per_protein(guo_pep_per_prot, diann_pep_per_prot, suppC_svg)
    print(f'Supplementary C figure saved to {suppC_svg}')
    print(f'Proteins with >=2 unique peptides: Guo={proteins_with_min_peptides(guo_pep_per_prot, 2):,}  quantmsdiann={proteins_with_min_peptides(diann_pep_per_prot, 2):,}')
    print('Loading HGNC symbol->Ensembl mapping...')
    hgnc_path = download_if_missing(HGNC_COMPLETE_SET_URL, _figure_original_vs_quantmsdiann__DATA_DIR / 'hgnc_complete_set.txt')
    symbol_to_ensg = load_hgnc_symbol_to_ensembl(hgnc_path)
    walzer_ensg = load_walzer_genes_ensembl(eprot73_path)
    diann_ensg, diann_unmapped = quantmsdiann_genes_as_ensembl(diann_genes_path, symbol_to_ensg, min_detection_fraction=0.5)
    suppD_svg = _figure_original_vs_quantmsdiann__FIGURES_DIR / 'supp_walzer_vs_quantms_genes_ensembl.svg'
    render_venn_diagram(walzer_ensg, diann_ensg, suppD_svg, left_label='Walzer 2022\n(E-PROT-73)', right_label='quantmsdiann\n(DIA-NN, $\\geq$50% of runs)', left_color='#90caf9', right_color=fs.COMPARISON['quantmsdiann'])
    inter_ensg = walzer_ensg & diann_ensg
    print(f'Supplementary D figure saved to {suppD_svg}')
    print(f'Ensembl gene IDs (DIA-NN >=50% of runs): Walzer={len(walzer_ensg):,}  quantmsdiann={len(diann_ensg):,}  intersection={len(inter_ensg):,}  Walzer-only={len(walzer_ensg - diann_ensg):,}  quantmsdiann-only={len(diann_ensg - walzer_ensg):,}  (unmapped HGNC={diann_unmapped:,})')
    print('Computing per-condition gene detections...')
    sdrf_path = download_if_missing(_figure_original_vs_quantmsdiann__QUANTMS_SDRF_URL, _figure_original_vs_quantmsdiann__DATA_DIR / 'PXD003539.sdrf.tsv')
    ea_design_path = download_if_missing(EA_EXPERIMENT_DESIGN_URL, _figure_original_vs_quantmsdiann__DATA_DIR / 'E-PROT-73-experiment-design.tsv')
    ea_downloads_path = download_if_missing(EA_DOWNLOADS_URL, _figure_original_vs_quantmsdiann__DATA_DIR / 'E-PROT-73-downloads.html')
    walzer_per_cond = walzer_genes_per_condition(eprot73_path, ea_downloads_path)
    diann_per_cond = quantmsdiann_genes_per_condition(diann_genes_path, sdrf_path, ea_design_path, symbol_to_ensg, min_global_detection_fraction=0.5)
    suppF_svg = _figure_original_vs_quantmsdiann__FIGURES_DIR / 'supp_genes_per_condition.svg'
    render_genes_per_condition(walzer_per_cond, diann_per_cond, suppF_svg)
    print(f'Supplementary F figure saved to {suppF_svg}')
    print('Per-condition gene detections (Walzer | quantmsdiann):')
    for cond in sorted(set(walzer_per_cond) | set(diann_per_cond)):
        w = len(walzer_per_cond.get(cond, set()))
        d = len(diann_per_cond.get(cond, set()))
        print(f'  {cond:32s} {w:>6,} | {d:>6,}')
    print('Computing Venn of protein accessions with >=2 unique peptides...')
    guo_acc = accessions_with_min_peptides_openswath(opensw_path, min_peptides=2)
    diann_acc = accessions_with_min_peptides_diann(pr_path, min_peptides=2)
    suppE_svg = _figure_original_vs_quantmsdiann__FIGURES_DIR / 'supp_venn_protein_accessions.svg'
    render_venn_diagram(guo_acc, diann_acc, suppE_svg)
    inter = guo_acc & diann_acc
    print(f'Supplementary E figure saved to {suppE_svg}')
    print(f'Accessions (>=2 unique peptides): Guo={len(guo_acc):,}  quantmsdiann={len(diann_acc):,}  intersection={len(inter):,}  Guo-only={len(guo_acc - diann_acc):,}  quantmsdiann-only={len(diann_acc - guo_acc):,}')
    print(f'Peptides:        Guo={guo_peptides:,}  Walzer={WALZER_PEPTIDES:,}  quantmsdiann={diann_peptides:,}')
    print(f'Protein groups:  Guo={guo_proteins:,}  Walzer={WALZER_PROTEINS:,}  quantmsdiann={diann_proteins:,}')
    print(f'Auxiliary:       Guo precursors={guo_precursors:,}  quantmsdiann precursors={diann_precursors:,}')
    print(f'EA context:     Walzer (E-PROT-73 genes)={counts.walzer_ea_genes:,}')
    if abs(diann_precursors - 117720) / 117720 > 0.01:
        print(f'WARNING: diann_precursors={diann_precursors:,} deviates >1% from expected 117,720', file=sys.stderr)
    if diann_proteins_log != 6927:
        print(f'WARNING: diann_proteins (log, unfiltered)={diann_proteins_log:,} != expected 6,927', file=sys.stderr)
    if abs(guo_precursors - 48374) / 48374 > 0.01:
        print(f'WARNING: guo_precursors={guo_precursors:,} deviates >1% from expected 48,374', file=sys.stderr)
    if abs(guo_proteins - 6556) / 6556 > 0.01:
        print(f'WARNING: guo_proteins={guo_proteins:,} deviates >1% from expected 6,556', file=sys.stderr)
    if abs(guo_peptides - 40592) / 40592 > 0.05:
        print(f'WARNING: guo_peptides={guo_peptides:,} deviates >5% from expected 40,592', file=sys.stderr)
    if abs(walzer_ea_genes - 2199) / 2199 > 0.02:
        print(f'WARN: Walzer EA gene count {walzer_ea_genes} differs from expected 2,199 by >2%', file=sys.stderr)
    return 0


# ======================================================================
# inlined from analysis/figure_performance_runtime.py
# ======================================================================

"""Performance scatter: workflow runtime vs DIA-NN threads.

One dot per analysis we have run end-to-end through the quantmsdiann nf-core
workflow:

  - 3 cell-line datasets (PXD003539, PXD030304, PXD004701) processed once
    with DIA-NN 2.5.0.
  - 4 ProteoBench benchmark datasets (PXD049412, PXD062685, PXD070049,
    ProteoBench_Module_7) processed once per DIA-NN release
    (1.8.1, 2.1.0, 2.2.0, 2.3.2, 2.5.0).

For each analysis we plot:

  x = DIA-NN `--threads N` from the quant_tables/diannsummary.log command line
  y = total Nextflow wallclock duration in hours, parsed from
      pipeline_info/pipeline_report.txt (the "duration: 3h 26m 38s" line)
  size = number of MS runs (count of `--f` arguments in the DIA-NN command
      line; matches SDRF row count for cell-line datasets)
  color = MS instrument family (TripleTOF 5600 / 6600, timsTOF SCP,
      Orbitrap Astral, ZenoTOF 7600)

Paper-ready output (no title, no footer) is written to
  analysis/figures/performance/runtime_vs_threads.{pdf,png,svg}

An auditable TSV (one row per dot) is written to
  analysis/figures/performance/runtime_data.tsv

We use pipeline_report.txt rather than aggregating nextflow_trace.txt because
some uploaded trace files are truncated (PXD062685 traces only contain the
SAMPLESHEET_CHECK + SDRF_PARSING rows; PXD070049/v2_3_2 is header-only),
whereas pipeline_report.txt always carries the run-level total. Both come
from the same Nextflow run, so the wallclock is identical when both are
intact.
"""
import re
import sys
from pathlib import Path
from typing import Iterable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
_figure_performance_runtime__REPO_ROOT = Path(__file__).resolve().parent.parent
_figure_performance_runtime__DATA_DIR = _figure_performance_runtime__REPO_ROOT / 'data'
_figure_performance_runtime__FIGURES_DIR = _figure_performance_runtime__REPO_ROOT / 'analysis' / 'figures' / 'performance'
CELL_LINE_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/absolute-expression-2.0/cell-lines'
BENCHMARK_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/proteobench/quantmsdiann_results'
CELL_LINE_DATASETS = ('PXD003539', 'PXD004701', 'PXD030304')
BENCHMARK_DATASETS = ('PXD049412', 'PXD062685', 'PXD070049', 'ProteoBench_Module_7')
_figure_performance_runtime__DIANN_VERSIONS = ('v1_8_1', 'v2_1_0', 'v2_2_0', 'v2_3_2', 'v2_5_0')
BENCHMARK_INSTRUMENT = {'PXD049412': 'Orbitrap Astral', 'PXD062685': 'timsTOF SCP', 'PXD070049': 'ZenoTOF 7600', 'ProteoBench_Module_7': 'Orbitrap Astral'}
_DURATION_PART_RE = re.compile('(?P<value>\\d+(?:\\.\\d+)?)\\s*(?P<unit>h|m|s|ms)\\b')
_PIPELINE_DURATION_RE = re.compile('duration:\\s*([0-9hms\\.\\sM]+?)\\s*\\)')
_THREADS_RE = re.compile('--threads\\s+(\\d+)\\b')
_F_ARG_RE = re.compile('--f\\s+(\\S+)')

def parse_duration_to_seconds(text: str) -> float:
    """Parse a Nextflow-style duration token like "3h 26m 38s" or "1m 9s"
    into seconds (float). Returns 0.0 for an empty string. Tokens may carry
    decimal seconds ("2.5s") and milliseconds ("250ms")."""
    s = text.strip()
    if not s:
        return 0.0
    total = 0.0
    matched = False
    for m in _DURATION_PART_RE.finditer(s):
        matched = True
        v = float(m.group('value'))
        unit = m.group('unit')
        if unit == 'h':
            total += v * 3600.0
        elif unit == 'm':
            total += v * 60.0
        elif unit == 's':
            total += v
        elif unit == 'ms':
            total += v / 1000.0
    if not matched:
        raise ValueError(f'Unrecognised duration string: {text!r}')
    return total

def parse_pipeline_report_duration(report_path: Path) -> float:
    """Extract the Nextflow wallclock seconds from a pipeline_report.txt
    file. Raises ValueError if the `duration:` line cannot be found.

    Format (one line near the bottom):
      The workflow was completed at <iso8601> (duration: 3h 26m 38s)
    """
    with open(report_path, encoding='utf-8') as fh:
        for line in fh:
            m = _PIPELINE_DURATION_RE.search(line)
            if m:
                return parse_duration_to_seconds(m.group(1))
    raise ValueError(f'`duration: ... )` not found in {report_path}')

def aggregate_nextflow_trace_duration(trace_path: Path) -> float:
    """Wallclock seconds inferred from a Nextflow trace.txt file by computing
    `max(submit + duration) - min(submit)` across rows.

    This is provided as a fallback / cross-check. We don't use it for the
    figure because several published traces in the benchmark set are
    truncated (header only or first two rows only); pipeline_report.txt is
    the authoritative source. Tested via fixtures.

    Returns 0.0 if the trace has zero data rows. Raises ValueError if a
    required column is missing on a non-empty trace.
    """
    df = pd.read_csv(trace_path, sep='\t', dtype=str)
    if df.empty:
        return 0.0
    for col in ('submit', 'duration'):
        if col not in df.columns:
            raise ValueError(f'Trace {trace_path} missing required column {col!r}')
    submit = pd.to_datetime(df['submit'], errors='coerce')
    dur_s = df['duration'].fillna('').map(parse_duration_to_seconds)
    finish = submit + pd.to_timedelta(dur_s, unit='s')
    if submit.isna().all():
        return 0.0
    span = finish.max() - submit.min()
    return float(span.total_seconds())

def parse_diann_command(log_path: Path) -> tuple[int, int, list[str]]:
    """Return (threads, n_files, run_basenames) parsed from a diannsummary.log
    file. The DIA-NN command is always on a single early line starting with
    `diann --lib ...`.

    Raises ValueError if the command line cannot be located or `--threads`
    is missing.
    """
    cmd: str | None = None
    with open(log_path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('diann ') and '--threads' in line:
                cmd = line
                break
    if cmd is None:
        raise ValueError(f'DIA-NN command line not found in {log_path}')
    m = _THREADS_RE.search(cmd)
    if not m:
        raise ValueError(f'--threads N flag not found in {log_path}')
    threads = int(m.group(1))
    runs = [m.group(1) for m in _F_ARG_RE.finditer(cmd)]
    return (threads, len(runs), runs)

def infer_instrument_from_runs(run_names: Iterable[str]) -> str | None:
    """Heuristic instrument-family classifier from run/raw file basenames.
    Used as a sanity check when SDRF is not available. Returns None if the
    name pattern is not recognised."""
    if not run_names:
        return None
    sample = ' '.join(list(run_names)[:8]).lower()
    if 'timstof' in sample or 'ttscp' in sample or sample.endswith('.d'):
        return 'timsTOF SCP'
    if 'zenotof' in sample or 'zeno' in sample:
        return 'ZenoTOF 7600'
    if 'astral' in sample:
        return 'Orbitrap Astral'
    return None

def parse_sdrf_instrument(sdrf_path: Path) -> str | None:
    """Return the comment[instrument] value (NT= field stripped) from an SDRF.
    Returns None if the column or NT= value is absent. If multiple distinct
    instruments appear, returns the most common one (PXD030304 SDRF mixes
    6 instrument identifiers but all map to the same NT=TripleTOF 6600)."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str, nrows=10000)
    col = None
    for c in df.columns:
        if c.strip().lower() == 'comment[instrument]':
            col = c
            break
    if col is None:
        return None
    vals = df[col].dropna().astype(str)
    if vals.empty:
        return None

    def _strip(v: str) -> str:
        if 'NT=' in v:
            inner = v.split('NT=', 1)[1]
            return inner.split(';', 1)[0].strip()
        return v.strip()
    cleaned = vals.map(_strip)
    return cleaned.value_counts().idxmax()

def collect_runtime_rows(*, fetch: bool=True) -> pd.DataFrame:
    """Build the per-analysis runtime table by reading (and optionally
    fetching) the small set of text files we need.

    Columns:
      dataset, version, instrument, threads, n_runs, wallclock_seconds,
      source_file_for_duration, source_file_for_threads
    """
    rows: list[dict] = []
    for dataset in CELL_LINE_DATASETS:
        dset_dir = _figure_performance_runtime__DATA_DIR / dataset
        report = dset_dir / 'pipeline_info' / 'pipeline_report.txt'
        log = dset_dir / 'diannsummary.log'
        sdrf = dset_dir / f'{dataset}.sdrf.tsv'
        if fetch:
            download_if_missing(f'{CELL_LINE_BASE}/{dataset}/pipeline_info/pipeline_report.txt', report)
            download_if_missing(f'{CELL_LINE_BASE}/{dataset}/quant_tables/diannsummary.log', log)
            download_if_missing(f'{CELL_LINE_BASE}/{dataset}/sdrf/{dataset}.sdrf.tsv', sdrf)
        secs = parse_pipeline_report_duration(report)
        threads, n_files, run_names = parse_diann_command(log)
        instrument = parse_sdrf_instrument(sdrf) if sdrf.exists() else None
        if instrument is None:
            instrument = infer_instrument_from_runs(run_names) or 'unknown'
        rows.append({'dataset': dataset, 'version': 'v2_5_0', 'instrument': instrument, 'threads': threads, 'n_runs': n_files, 'wallclock_seconds': secs, 'source_file_for_duration': str(report.relative_to(_figure_performance_runtime__REPO_ROOT)), 'source_file_for_threads': str(log.relative_to(_figure_performance_runtime__REPO_ROOT))})
    for dataset in BENCHMARK_DATASETS:
        for version in _figure_performance_runtime__DIANN_VERSIONS:
            base = f'{BENCHMARK_BASE}/{dataset}/{version}'
            dset_dir = _figure_performance_runtime__DATA_DIR / 'quantmsdiann_benchmarks' / dataset / version
            report = dset_dir / 'pipeline_info' / 'pipeline_report.txt'
            log = dset_dir / 'quant_tables' / 'diannsummary.log'
            if fetch:
                download_if_missing(f'{base}/pipeline_info/pipeline_report.txt', report)
                download_if_missing(f'{base}/quant_tables/diannsummary.log', log)
            secs = parse_pipeline_report_duration(report)
            threads, n_files, run_names = parse_diann_command(log)
            instrument = infer_instrument_from_runs(run_names) or BENCHMARK_INSTRUMENT.get(dataset, 'unknown')
            rows.append({'dataset': dataset, 'version': version, 'instrument': instrument, 'threads': threads, 'n_runs': n_files, 'wallclock_seconds': secs, 'source_file_for_duration': str(report.relative_to(_figure_performance_runtime__REPO_ROOT)), 'source_file_for_threads': str(log.relative_to(_figure_performance_runtime__REPO_ROOT))})
    return pd.DataFrame(rows)
INSTRUMENT_COLOURS = dict(fs.INSTRUMENT_COLORS)


# ======================================================================
# inlined from analysis/figure_performance_trace.py
# ======================================================================

"""Performance plots derived from per-task nextflow_trace.txt files.

Two paper-ready panels are produced from the same set of 20
quantmsdiann x ProteoBench benchmark runs (4 datasets x 5 DIA-NN versions);
the 3 cell-line datasets are excluded because their nextflow_trace.txt files
are not published on PRIDE.

Plot 1 (`parallelism_vs_wallclock`)
-----------------------------------

One dot per analysis. X-axis = peak concurrent task count observed during
the workflow (max overlap of [submit, submit+duration] intervals across all
trace rows). Y-axis = total workflow wallclock in hours, defined as
`max(submit + duration) - min(submit)` across the same trace. Colour =
instrument family (re-uses `INSTRUMENT_COLOURS` from
`analysis.figure_performance_runtime`). MS run count is constant (6) across
all benchmarks so the dot size is held at a single readable value and the
size dimension is intentionally not used.

Plot 2 (`runtime_per_step`)
---------------------------

Horizontal box plot of per-task `duration` in seconds, one row per workflow
step, ordered by descending median. The step name is the last colon-
separated segment of `name` before the parenthesised input argument
(e.g. `SAMPLESHEET_CHECK`, `DIANN_RUN`-equivalents like `PRELIMINARY_ANALYSIS`
and `INDIVIDUAL_ANALYSIS`). All 20 traces are aggregated. FAILED tasks
(Nextflow retry attempts) are excluded from the per-step duration
distributions; they remain in the parallelism analysis because they did
occupy cluster slots. The x-axis is log scale when the observed range spans
more than two orders of magnitude.

Truncated traces
----------------

PRIDE publishes nextflow_trace.txt files for all 20 analyses but several
copies are truncated (PXD062685 x5 versions ship the first two rows only,
PXD070049/v2_3_2 ships header-only). For Plot 1 we still emit a TSV row for
every analysis but mark `complete = False` for those with <= 2 data rows and
omit them from the scatter (legend annotation explains the count). For
Plot 2 we ingest whatever rows exist from every trace (the partial PXD062685
rows contribute valid SAMPLESHEET_CHECK / SDRF_PARSING durations).
"""
import sys
from pathlib import Path
from typing import Iterable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
_figure_performance_trace__REPO_ROOT = Path(__file__).resolve().parent.parent
_figure_performance_trace__DATA_DIR = _figure_performance_trace__REPO_ROOT / 'data'
_figure_performance_trace__FIGURES_DIR = _figure_performance_trace__REPO_ROOT / 'analysis' / 'figures' / 'performance'
MIN_ROWS_COMPLETE = 5

def extract_step_name(full_name: str) -> str:
    """Return the step identifier from a Nextflow `name` cell.

    The trace stores fully-qualified process paths like
    `BIGBIO_QUANTMSDIANN:QUANTMSDIANN:INPUT_CHECK:SAMPLESHEET_CHECK (PXD049412.sdrf.tsv)`.
    The step is the last `:`-separated segment, with any parenthesised input
    argument stripped. Rows without an argument (e.g. `SUMMARY_PIPELINE`)
    return the bare step name.
    """
    s = (full_name or '').strip()
    if not s:
        return ''
    paren = s.find(' (')
    if paren != -1:
        s = s[:paren]
    if ':' in s:
        s = s.rsplit(':', 1)[1]
    return s.strip()

def load_report_window_data(report_path: Path) -> pd.DataFrame:
    """Read a Nextflow-generated `nextflow_report.html` and return its
    embedded `window.data.trace` table as a normalised DataFrame.

    Nextflow embeds the full per-task trace inside the HTML report as a
    JavaScript object literal. We brace-balance the `window.data = {...}`
    assignment, repair its invalid `\\'` escapes (Nextflow over-escapes
    single quotes inside JSON strings), and lift the `trace` array.

    The report is strictly better than `nextflow_trace.txt` for our purposes:
    it's published for every dataset (under `pipeline_info/` for cell-line
    runs, alongside the dataset for benchmark runs), and it contains the
    full task list even when the sibling `trace.txt` is upstream-truncated.

    Returns columns: `step`, `status`, `submit` (epoch seconds, float),
    `duration_s` (float), `finish` (epoch seconds, float). Returns an empty
    DataFrame with those columns when the HTML cannot be parsed or the
    trace array is empty."""
    import json
    cols = ['step', 'status', 'submit', 'duration_s', 'finish']
    try:
        html = report_path.read_text(encoding='utf-8')
    except (FileNotFoundError, OSError):
        return pd.DataFrame(columns=cols)
    start = html.find('window.data = ')
    if start < 0:
        return pd.DataFrame(columns=cols)
    i = html.find('{', start)
    depth, in_str, esc, end = (0, False, False, -1)
    for j in range(i, len(html)):
        ch = html[j]
        if esc:
            esc = False
            continue
        if ch == '\\':
            esc = True
            continue
        if ch == '"':
            in_str = not in_str
            continue
        if in_str:
            continue
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    if end < 0:
        return pd.DataFrame(columns=cols)
    blob = html[i:end].replace("\\'", "'")
    try:
        data = json.loads(blob)
    except json.JSONDecodeError:
        return pd.DataFrame(columns=cols)
    trace = data.get('trace', [])
    if not trace:
        return pd.DataFrame(columns=cols)
    rows: list[dict] = []
    for rec in trace:
        submit_ms = rec.get('submit')
        duration_ms = rec.get('duration')
        try:
            submit_s = float(submit_ms) / 1000.0 if submit_ms is not None else float('nan')
        except (TypeError, ValueError):
            submit_s = float('nan')
        try:
            duration_s = float(duration_ms) / 1000.0 if duration_ms is not None else 0.0
        except (TypeError, ValueError):
            duration_s = 0.0
        process = rec.get('process') or rec.get('name') or ''
        step = extract_step_name(process)
        rows.append({'step': step, 'status': rec.get('status') or '', 'submit': submit_s, 'duration_s': duration_s, 'finish': submit_s + duration_s})
    return pd.DataFrame(rows, columns=cols)

def load_trace(trace_path: Path) -> pd.DataFrame:
    """Read a nextflow_trace.txt as a normalised DataFrame.

    Returns columns: `step`, `status`, `submit`, `duration_s`, `finish`.
    Drops rows with unparseable submit timestamps. Empty traces (header
    only) return an empty DataFrame with the expected columns. FAILED tasks
    (Nextflow retry attempts) are kept so that callers can decide whether
    to include them: the parallelism analysis counts them (they did occupy
    cluster slots while running), the per-step duration analysis filters to
    COMPLETED.
    """
    df = pd.read_csv(trace_path, sep='\t', dtype=str)
    if df.empty:
        return pd.DataFrame(columns=['step', 'status', 'submit', 'duration_s', 'finish'])
    submit = pd.to_datetime(df['submit'], errors='coerce')
    dur_s = df['duration'].fillna('').map(parse_duration_to_seconds)
    step = df['name'].fillna('').map(extract_step_name)
    status = df.get('status', pd.Series([''] * len(df))).fillna('')
    out = pd.DataFrame({'step': step, 'status': status, 'submit': submit, 'duration_s': dur_s})
    out = out.dropna(subset=['submit']).reset_index(drop=True)
    out['finish'] = out['submit'] + pd.to_timedelta(out['duration_s'], unit='s')
    return out

def peak_concurrent_tasks(df: pd.DataFrame) -> tuple[int, float]:
    """Compute (peak, median) concurrent-task counts from a normalised
    trace DataFrame produced by `load_trace`.

    Peak = max number of intervals [submit, finish] that overlap at any
    instant. Median = median concurrency *over the union of event points*
    (so it counts the concurrency seen at each distinct submit/finish
    instant). This biases the median toward periods with more activity,
    which is what we want for a "how parallel did this workflow get" metric.

    Returns (0, 0.0) for an empty frame.
    """
    if df.empty:
        return (0, 0.0)
    events: list[tuple[pd.Timestamp, int]] = []
    for _, row in df.iterrows():
        events.append((row['submit'], +1))
        events.append((row['finish'], -1))
    events.sort(key=lambda x: (x[0], x[1]))
    concurrency = 0
    peak = 0
    samples: list[int] = []
    for _, delta in events:
        concurrency += delta
        peak = max(peak, concurrency)
        if concurrency > 0:
            samples.append(concurrency)
    if not samples:
        return (peak, 0.0)
    med = float(pd.Series(samples).median())
    return (peak, med)
_SIZE_UNITS = {'': 1, 'B': 1, 'KB': 1024, 'MB': 1024 ** 2, 'GB': 1024 ** 3, 'TB': 1024 ** 4}

def parse_size_to_bytes(text: str) -> float:
    """Parse a Nextflow size string like `"165 MB"` / `"4.7 GB"` / `"3.2 KB"`
    to bytes. Empty / `"-"` / unparseable strings return NaN."""
    s = (text or '').strip()
    if not s or s in {'-', '0'}:
        return float('nan') if not s or s == '-' else 0.0
    parts = s.split()
    if len(parts) == 1:
        try:
            return float(parts[0])
        except ValueError:
            return float('nan')
    try:
        value = float(parts[0])
    except ValueError:
        return float('nan')
    unit = parts[1].strip().upper()
    factor = _SIZE_UNITS.get(unit)
    if factor is None:
        return float('nan')
    return value * factor

def parse_pct_cpu(text: str) -> float:
    """Parse a Nextflow `%cpu` cell like `"94.8%"` to a float (94.8). Empty
    / `"-"` strings return NaN; values are *not* capped at 100 because
    multi-thread tasks can legitimately exceed 100 %."""
    s = (text or '').strip().rstrip('%').strip()
    if not s or s == '-':
        return float('nan')
    try:
        return float(s)
    except ValueError:
        return float('nan')

def load_trace_resources(trace_path: Path) -> pd.DataFrame:
    """Read a `nextflow_trace.txt` and return per-task resource rows.

    Columns: `step`, `status`, `peak_rss_bytes`, `pct_cpu`,
    `duration_s`. Drops rows with empty step or status. FAILED rows
    are kept; callers filter as needed (the F2c box plot uses COMPLETED
    rows only — same convention as `aggregate_step_durations`).

    Empty / header-only traces return an empty DataFrame with the
    expected columns."""
    cols = ['step', 'status', 'peak_rss_bytes', 'pct_cpu', 'duration_s']
    try:
        df = pd.read_csv(trace_path, sep='\t', dtype=str)
    except (FileNotFoundError, OSError, pd.errors.EmptyDataError):
        return pd.DataFrame(columns=cols)
    if df.empty:
        return pd.DataFrame(columns=cols)
    needed = {'name', 'status', '%cpu', 'peak_rss', 'duration'}
    if not needed.issubset(df.columns):
        return pd.DataFrame(columns=cols)
    out = pd.DataFrame({'step': df['name'].fillna('').map(extract_step_name), 'status': df['status'].fillna(''), 'peak_rss_bytes': df['peak_rss'].fillna('').map(parse_size_to_bytes), 'pct_cpu': df['%cpu'].fillna('').map(parse_pct_cpu), 'duration_s': df['duration'].fillna('').map(parse_duration_to_seconds)})
    out = out[out['step'].astype(bool)].reset_index(drop=True)
    return out

def aggregate_step_resources(traces: Iterable[pd.DataFrame], *, completed_only: bool=True) -> dict[str, dict[str, list[float]]]:
    """Aggregate per-task `peak_rss_bytes` and `pct_cpu` across traces,
    keyed first by step, then by metric. Mirrors
    `aggregate_step_durations` but for resource columns. Returns a dict
    `{step: {"peak_rss_bytes": [...], "pct_cpu": [...]}}`. NaN values
    are dropped per-metric so a row missing one column still contributes
    the other."""
    out: dict[str, dict[str, list[float]]] = {}
    for df in traces:
        if df.empty:
            continue
        sub = df
        if completed_only and 'status' in df.columns:
            sub = df[df['status'] == 'COMPLETED']
        for step, rss, cpu in zip(sub['step'], sub['peak_rss_bytes'], sub['pct_cpu']):
            if not step:
                continue
            bucket = out.setdefault(step, {'peak_rss_bytes': [], 'pct_cpu': []})
            if pd.notna(rss):
                bucket['peak_rss_bytes'].append(float(rss))
            if pd.notna(cpu):
                bucket['pct_cpu'].append(float(cpu))
    return out

def trace_wallclock_seconds(df: pd.DataFrame) -> float:
    """`max(submit + duration) - min(submit)` in seconds. Returns 0.0 for an
    empty frame. Handles both pandas Timestamp/Timedelta (from the legacy
    nextflow_trace.txt loader) and plain epoch-second floats (from the new
    nextflow_report.html loader)."""
    if df.empty:
        return 0.0
    span = df['finish'].max() - df['submit'].min()
    if hasattr(span, 'total_seconds'):
        return float(span.total_seconds())
    return float(span)
INSILICO_STEP = 'INSILICO_LIBRARY_GENERATION'

def insilico_seconds(df: pd.DataFrame) -> float:
    """Total wall-clock seconds spent in INSILICO_LIBRARY_GENERATION for a
    single run, summed over its task rows (normally one). Accepts either a
    `load_trace`/`load_report_window_data` frame (`step` + `duration_s`).
    Returns 0.0 when the step is absent (e.g. a user-supplied library run or
    a truncated trace). COMPLETED-only when a status column is present, so a
    failed-then-retried attempt is not double counted."""
    if df is None or df.empty or 'step' not in df.columns:
        return 0.0
    sub = df[df['step'] == INSILICO_STEP]
    if 'status' in sub.columns and sub['status'].astype(bool).any():
        completed = sub[sub['status'] == 'COMPLETED']
        if not completed.empty:
            sub = completed
    if sub.empty or 'duration_s' not in sub.columns:
        return 0.0
    return float(sub['duration_s'].fillna(0.0).sum())

def busy_span_seconds(df: pd.DataFrame) -> float:
    """Total wall-clock seconds during which at least one task was running —
    the union of all [submit, finish] intervals. Unlike `trace_wallclock_seconds`
    (raw span) this ignores idle gaps between `-resume` invocations, so it is a
    resume-robust "active compute time". For a single uninterrupted run it
    equals the raw span (and the pipeline_report.txt duration). Returns 0.0 for
    an empty frame."""
    if df is None or df.empty:
        return 0.0
    iv = sorted(((s, f) for s, f in zip(df['submit'], df['finish']) if pd.notna(s) and pd.notna(f)))
    if not iv:
        return 0.0
    total = cur_s = cur_e = None
    for s, e in iv:
        if cur_s is None:
            cur_s, cur_e, total = (s, e, 0.0)
        elif s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += (cur_e - cur_s).total_seconds()
            cur_s, cur_e = (s, e)
    total += (cur_e - cur_s).total_seconds()
    return float(total)

def aggregate_step_durations(traces: Iterable[pd.DataFrame], *, completed_only: bool=True) -> dict[str, list[float]]:
    """Aggregate per-task `duration_s` values across many traces, keyed by
    step name. Empty / missing-step rows are dropped. When `completed_only`
    is True (default) we drop FAILED rows so the duration distribution
    reflects task success times — important for the per-step box plot
    because PXD070049/v1_8_1 has ~14 failed SDRF_PARSING retry attempts
    that would otherwise dominate the distribution."""
    out: dict[str, list[float]] = {}
    for df in traces:
        if df.empty:
            continue
        sub = df
        if completed_only and 'status' in df.columns:
            sub = df[df['status'] == 'COMPLETED']
        for step, dur in zip(sub['step'], sub['duration_s']):
            if not step:
                continue
            out.setdefault(step, []).append(float(dur))
    return out
_PARAMS_FILENAME_RE = __import__('re').compile('params_(\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2})\\.json')
_PIPELINE_REPORT_COMPLETION_RE = __import__('re').compile('completed at\\s+(\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}(?:\\.\\d+)?)')

def list_params_timestamps(pipeline_info_url: str) -> list:
    """List `params_<ISO>.json` timestamps from a pipeline_info/ directory
    listing. Each `params_*.json` file corresponds to one Nextflow
    invocation; the earliest is the start of the user's wallclock.

    Returns a list of `datetime.datetime` objects (naive, in the cluster's
    local timezone) sorted ascending. The runs are typically run within the
    same TZ so naive comparison is fine for span computation."""
    import re, requests, datetime as _dt
    listing = requests.get(pipeline_info_url, timeout=30).text
    out: list[_dt.datetime] = []
    seen: set[str] = set()
    for m in _PARAMS_FILENAME_RE.finditer(listing):
        token = m.group(1)
        if token in seen:
            continue
        seen.add(token)
        out.append(_dt.datetime.strptime(token, '%Y-%m-%d_%H-%M-%S'))
    return sorted(out)

def parse_pipeline_report_completion(report_path):
    """Read `pipeline_info/pipeline_report.txt` and return the workflow
    completion datetime (naive, cluster-local). Raises ValueError if the
    line isn't present (e.g. an aborted run)."""
    import datetime as _dt
    text = report_path.read_text(encoding='utf-8')
    m = _PIPELINE_REPORT_COMPLETION_RE.search(text)
    if not m:
        raise ValueError(f"'completed at <ISO>' line not found in {report_path}")
    iso = m.group(1)
    if '.' in iso:
        head, frac = iso.split('.', 1)
        return _dt.datetime.fromisoformat(head + '.' + frac[:6])
    return _dt.datetime.fromisoformat(iso)

def total_wallclock_with_resumes_seconds(pipeline_info_url: str, report_path) -> float:
    """Total wallclock the user actually waited, including all `-resume`
    re-runs: earliest `params_*.json` timestamp to `pipeline_report.txt`'s
    "completed at" datetime. For single-invocation workflows this matches
    the report's reported duration; for multi-invocation workflows it's
    larger because the cluster idle/queue time between resumes counts."""
    starts = list_params_timestamps(pipeline_info_url)
    if not starts:
        raise ValueError(f'No params_*.json files listed at {pipeline_info_url}')
    completed = parse_pipeline_report_completion(report_path)
    return (completed - starts[0]).total_seconds()

def count_ms_runs(dataset_id: str, sdrf_path, log_path) -> int:
    """Number of MS runs in the analysis. For cell-line datasets we count
    SDRF rows (one per data file). For benchmarks we count `--f` arguments
    in the DIA-NN command line at the head of `diannsummary.log`. Either
    source is correct for the matching dataset family; we fall back from
    SDRF to the DIA-NN log if the SDRF is unavailable."""
    if sdrf_path.exists():
        n = 0
        with open(sdrf_path, 'r', encoding='utf-8') as fh:
            for i, _ in enumerate(fh):
                n = i
        return n
    _, n_files, _ = parse_diann_command(log_path)
    return n_files
CELL_LINE_ANALYSES: dict[str, dict[str, str]] = {'PXD003539': {'instrument': 'TripleTOF 5600', 'version': 'v2_5_0'}, 'PXD030304': {'instrument': 'TripleTOF 6600', 'version': 'v2_5_0'}, 'PXD004701': {'instrument': 'TripleTOF 5600', 'version': 'v2_5_0'}}
PERF_TRACE_DIR = _figure_performance_trace__DATA_DIR / 'perf_traces'
PERF_DATASETS: tuple[tuple[str, str, int], ...] = (('PXD049412', 'Orbitrap Astral', 6), ('PXD062685', 'timsTOF SCP', 6), ('PXD070049', 'ZenoTOF 7600', 6), ('ProteoBench_Module_7', 'Orbitrap Astral', 6), ('PXD003539', 'TripleTOF 5600', 120), ('PXD004701', 'TripleTOF 5600', 300), ('PXD030304', 'TripleTOF 6600', 5798), ('MSV000093870', 'Q Exactive', 38), ('PXD034128', 'timsTOF Pro', 7), ('PXD034623', 'Orbitrap Exploris 480', 63), ('PXD064049', 'timsTOF SCP', 12), ('PXD049692', 'timsTOF HT', 10), ('PXD046357', 'Orbitrap Astral', 12), ('PXD017199', 'Q Exactive', 206), ('PXD041421', 'timsTOF Pro', 48))

def iter_perf_trace_paths() -> Iterable[Path]:
    """Yield the staged per-dataset `nextflow_trace.txt` paths under
    data/perf_traces/ (one representative fresh run per dataset). Filesystem
    only — these are staged from the cluster, never fetched from FTP. Missing
    datasets are skipped."""
    for dataset, *_ in PERF_DATASETS:
        trace = PERF_TRACE_DIR / dataset / 'nextflow_trace.txt'
        if trace.exists():
            yield trace

def _n_runs_from_trace(df: pd.DataFrame, fallback: int) -> int:
    """Count MS data files = per-file INDIVIDUAL_ANALYSIS tasks in the trace
    (all statuses, since `-resume`d runs mark earlier-leg tasks CACHED rather
    than COMPLETED). Falls back to the known SDRF file count when the trace
    layout doesn't emit that step."""
    if df is None or df.empty or 'step' not in df.columns:
        return fallback
    n = int((df['step'] == 'INDIVIDUAL_ANALYSIS').sum())
    return n if n > 0 else fallback

def collect_parallelism_rows(*, fetch: bool=True) -> pd.DataFrame:
    """One row per dataset for the Fig 1c "wall-clock time to finish" bars,
    built from the fresh per-dataset traces staged under data/perf_traces/.

    Wall-clock = the run's pipeline_report.txt "duration" when present, else the
    resume-robust `busy_span_seconds` of the trace (the two agree for single
    uninterrupted runs). `insilico_lib_seconds` is the INSILICO_LIBRARY_GENERATION
    time from the SAME trace; `wallclock_seconds` is wall-clock minus that, so the
    bar reflects per-file quantification time rather than the one-time,
    cohort-independent library-prediction cost. `wallclock_with_lib_seconds`
    keeps the raw value for the audit TSV. `fetch` is accepted for API
    compatibility but unused — the traces are staged locally.
    """
    rows: list[dict] = []
    for dataset, instrument, fallback_runs in PERF_DATASETS:
        tdir = PERF_TRACE_DIR / dataset
        trace_path = tdir / 'nextflow_trace.txt'
        if not trace_path.exists():
            continue
        df = load_trace(trace_path)
        n_tasks = int(len(df))
        report_txt = tdir / 'pipeline_report.txt'
        wallclock = None
        if report_txt.exists():
            try:
                wallclock = parse_pipeline_report_duration(report_txt)
            except (ValueError, FileNotFoundError, OSError):
                wallclock = None
        if wallclock is None:
            wallclock = busy_span_seconds(df)
        separable = n_tasks >= MIN_ROWS_COMPLETE
        lib_seconds = insilico_seconds(df) if separable else 0.0
        n_runs = _n_runs_from_trace(df, fallback_runs)
        rows.append({'dataset': dataset, 'version': 'v2_5_1', 'instrument': instrument, 'n_runs': n_runs, 'n_tasks_observed': n_tasks, 'n_invocations': 1, 'wallclock_seconds': max(0.0, float(wallclock) - float(lib_seconds)), 'wallclock_with_lib_seconds': float(wallclock), 'insilico_lib_seconds': float(lib_seconds), 'insilico_separable': bool(separable), 'trace_complete': separable, 'source_file': str(trace_path.relative_to(_figure_performance_trace__REPO_ROOT))})
    return pd.DataFrame(rows)

def _iter_pxd071075_sweep_trace_paths() -> Iterable[Path]:
    """Yield local nextflow_trace.txt paths for the 5 PXD071075
    cluster-node sweep points. Tolerates both zero-padded (q010) and
    bare (q10) directory layouts. Missing sweep points are skipped
    silently — callers (F2b/F2c aggregators) are happy with whatever
    subset is staged."""
    sweep_dir = _figure_performance_trace__DATA_DIR / 'queue_size_sweep'
    if not sweep_dir.exists():
        return
    for q in (10, 50, 100, 200, 300):
        for cand in (sweep_dir / f'q{q:03d}', sweep_dir / f'q{q}'):
            trace = cand / 'nextflow_trace.txt'
            if trace.exists():
                yield trace
                break

def _load_pxd071075_sweep_traces() -> list[pd.DataFrame]:
    """Return the per-task DataFrames for the PXD071075 sweep traces
    in `load_trace` shape (`step, status, submit, duration_s, finish`).
    Used by `collect_step_runtime_rows` to pool sweep tasks into the
    F2b distribution alongside benchmark tasks."""
    return [load_trace(p) for p in _iter_pxd071075_sweep_trace_paths()]

def _load_pxd071075_sweep_resource_rows() -> list[pd.DataFrame]:
    """Return per-task resource DataFrames for the PXD071075 sweep
    (`step, status, peak_rss_bytes, pct_cpu, duration_s`). Used by
    `collect_step_resource_rows` to pool sweep tasks into the F2c
    distribution alongside benchmark tasks."""
    return [load_trace_resources(p) for p in _iter_pxd071075_sweep_trace_paths()]

def collect_pxd071075_sweep_rows() -> pd.DataFrame:
    """Add the PXD071075 single-cell `executor.queueSize` sweep points
    to the parallelism scatter. Same schema as
    `collect_parallelism_rows` so the two frames can be concatenated.

    Each sweep point sits at the same n_runs (PXD071075's 2,310 input
    runs) and the same instrument (Orbitrap Eclipse), but varies the
    Nextflow queueSize. The resulting vertical strip in the F2a
    scatter is the "elastic" wallclock dimension — orthogonal to the
    workload-scaling axis the other cohorts demonstrate. The extra
    `queue_size` column lets the renderer annotate each point with
    its queueSize without polluting non-sweep rows.

    Reads `data/queue_size_sweep/q<NNN>/nextflow_trace.txt` for each
    sweep point. Returns an empty DataFrame (right shape) when no
    sweep data is staged."""
    sweep_dir = _figure_performance_trace__DATA_DIR / 'queue_size_sweep'
    rows: list[dict] = []
    if not sweep_dir.exists():
        return pd.DataFrame(columns=['dataset', 'version', 'instrument', 'n_runs', 'n_tasks_observed', 'n_invocations', 'wallclock_seconds', 'trace_complete', 'source_file', 'queue_size'])
    n_runs = 0
    for q in (10, 50, 100, 200, 300):
        for cand in (sweep_dir / f'q{q:03d}', sweep_dir / f'q{q}'):
            log = cand / 'diannsummary.log'
            if log.exists():
                try:
                    with open(log, encoding='utf-8') as fh:
                        text = fh.read()
                    import re as _re
                    n_runs = len(_re.findall('--f \\S+', text))
                except (FileNotFoundError, OSError):
                    n_runs = 0
                break
        if n_runs:
            break
    for q in (10, 50, 100, 200, 300):
        trace_path = None
        for cand in (sweep_dir / f'q{q:03d}', sweep_dir / f'q{q}'):
            if (cand / 'nextflow_trace.txt').exists():
                trace_path = cand / 'nextflow_trace.txt'
                break
        if trace_path is None:
            continue
        df_trace = load_trace(trace_path)
        if df_trace.empty:
            continue
        wallclock = trace_wallclock_seconds(df_trace)
        n_tasks = int(len(df_trace))
        rows.append({'dataset': 'PXD071075', 'version': f'v2_5_0_q{q:03d}', 'instrument': 'Orbitrap Eclipse', 'n_runs': n_runs, 'n_tasks_observed': n_tasks, 'n_invocations': 1, 'wallclock_seconds': float(wallclock), 'trace_complete': n_tasks >= MIN_ROWS_COMPLETE, 'source_file': str(trace_path.relative_to(_figure_performance_trace__REPO_ROOT)), 'queue_size': q})
    return pd.DataFrame(rows)

def collect_step_runtime_rows(*, fetch: bool=True) -> tuple[dict[str, list[float]], pd.DataFrame]:
    """Aggregate per-step task durations across the staged per-dataset fresh
    traces (data/perf_traces/<dataset>/nextflow_trace.txt — one representative
    run per dataset) plus the 5 PXD071075 cluster-node sweep runs. Filesystem
    only; nothing is fetched from FTP. The PXD071075 sweep dominates the
    high-count per-file steps because each sweep point runs the full ~2,310-file
    cohort, so the per-step shape is weighted toward real-cohort scale.

    Returns (durations_by_step, summary_df) where summary_df has columns
    `step, n, p05_seconds, p25_seconds, p50_seconds, p75_seconds,
    p95_seconds, min_seconds, max_seconds`.
    """
    traces: list[pd.DataFrame] = [load_trace(p) for p in iter_perf_trace_paths()]
    traces.extend(_load_pxd071075_sweep_traces())
    durations = aggregate_step_durations(traces)
    summary: list[dict] = []
    for step, vals in durations.items():
        s = pd.Series(vals)
        summary.append({'step': step, 'n': int(len(s)), 'p05_seconds': float(s.quantile(0.05)), 'p25_seconds': float(s.quantile(0.25)), 'p50_seconds': float(s.quantile(0.5)), 'p75_seconds': float(s.quantile(0.75)), 'p95_seconds': float(s.quantile(0.95)), 'min_seconds': float(s.min()), 'max_seconds': float(s.max())})
    df = pd.DataFrame(summary).sort_values('p50_seconds', ascending=False)
    return (durations, df.reset_index(drop=True))

def render_parallelism_scatter(df: pd.DataFrame, svg_path: Path | None=None, *, ax: plt.Axes | None=None, legend_ncol: int=4, legend_bbox_y: float=-0.14, composite: bool=False, show_legend: bool=True, short_labels: bool=False) -> None:
    """Horizontal bar chart of the wall-clock time each reanalysis took to
    finish (hours, from pipeline_report.txt), one bar per analysis, ordered
    by cohort size (number of MS data files) and coloured by instrument
    family. Time-to-finish stays within a few hours from a handful to
    thousands of files because per-file work is parallelised across the
    cluster, so it does not grow with cohort size.

    Pass `ax` to draw into an existing axes (composite figures); omit
    `svg_path` in that mode."""
    plot_df = df.copy()
    own_fig = ax is None
    reps = []
    for _ds, g in plot_df.groupby('dataset'):
        if 'version' in g.columns and (g['version'] == 'v2_5_0').any():
            r = g[g['version'] == 'v2_5_0'].iloc[0]
        else:
            gs = g.sort_values('wallclock_seconds')
            r = gs.iloc[len(gs) // 2]
        reps.append(r)
    rep = pd.DataFrame(reps).copy()
    rep['hours'] = rep['wallclock_seconds'] / 3600.0
    rep = rep.sort_values('n_runs').reset_index(drop=True)

    def _bar_label(row) -> str:
        if short_labels:
            return f"{row['dataset']} ({int(row['n_runs']):,})"
        base = f"{row['dataset']}  ({int(row['n_runs']):,} files"
        q = row.get('queue_size')
        if q is not None and pd.notna(q):
            base += f', {int(q)} nodes'
        return base + ')'
    rep['label'] = rep.apply(_bar_label, axis=1)
    if own_fig:
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
    label_size = 8 if composite else 12
    tick_size = 6.5 if composite else 10
    ann_size = 6 if composite else 9
    y = list(range(len(rep)))
    colours = [INSTRUMENT_COLOURS.get(i, '#9e9e9e') for i in rep['instrument']]
    ax.barh(y, rep['hours'], color=colours, edgecolor='#222222', linewidth=0.6, height=0.7)
    hmax = float(rep['hours'].max()) if len(rep) else 1.0
    for yi, h in zip(y, rep['hours']):
        ax.text(h + hmax * 0.012, yi, f'{h:.1f} h', va='center', ha='left', fontsize=ann_size, color='#333333')
    ax.set_yticks(y)
    ax.set_yticklabels(rep['label'], fontsize=tick_size)
    ax.invert_yaxis()
    ax.set_xlabel('Wall-clock time to finish (hours)', fontsize=label_size)
    ax.set_xlim(0, hmax * 1.16)
    ax.tick_params(axis='x', labelsize=tick_size)
    fs.despine(ax)
    if show_legend:
        from matplotlib.patches import Patch
        seen: dict[str, str] = {}
        for inst in rep['instrument']:
            seen.setdefault(inst, INSTRUMENT_COLOURS.get(inst, '#9e9e9e'))
        handles = [Patch(facecolor=c, edgecolor='#222222', label=i) for i, c in seen.items()]
        legend_fs = 6 if composite else 10
        legend_title_fs = 7 if composite else 10
        ncol = 2 if composite else legend_ncol
        bbox_y = -0.22 if composite else legend_bbox_y
        ax.legend(handles=handles, title='Instrument', loc='upper center', bbox_to_anchor=(0.5, bbox_y), fontsize=legend_fs, title_fontsize=legend_title_fs, frameon=False, borderaxespad=0.0, ncol=ncol, columnspacing=1.2 if composite else 1.6)
    if own_fig:
        fig.tight_layout()
        assert svg_path is not None
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg_path, bbox_inches='tight')
        plt.close(fig)

def render_per_step_boxplot(durations: dict[str, list[float]], summary: pd.DataFrame, svg_path: Path, fig_h: float | None=None) -> None:
    """Horizontal box plot of per-task durations, one row per step, ordered
    by descending median. `fig_h` overrides the auto-scaled figure height
    (used to match a companion panel's height when placed side by side)."""
    steps = summary['step'].tolist()
    data = [durations[s] for s in steps]
    if fig_h is None:
        fig_h = max(3.5, 0.45 * len(steps) + 1.5)
    fig, ax = plt.subplots(figsize=(8.0, fig_h))
    bp = ax.boxplot(data, vert=False, whis=(5, 95), showfliers=True, flierprops=dict(marker='o', markerfacecolor='#888888', markeredgecolor='none', markersize=3.5, alpha=0.6), medianprops=dict(color=fs.OKABE_ITO['vermillion'], linewidth=1.6), boxprops=dict(color=fs.OKABE_ITO['blue'], linewidth=1.2), whiskerprops=dict(color=fs.OKABE_ITO['blue'], linewidth=1.0), capprops=dict(color=fs.OKABE_ITO['blue'], linewidth=1.0), patch_artist=False)
    del bp
    ax.set_yticks(range(1, len(steps) + 1))
    ax.set_yticklabels(steps, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Task duration (s)', fontsize=10)
    fs.despine(ax)
    ax.tick_params(axis='both', labelsize=9)
    flat = [v for vals in data for v in vals if v > 0]
    if flat:
        ratio = max(flat) / max(min(flat), 1e-06)
        if ratio > 100:
            ax.set_xscale('log')
            ax.set_xlim(left=max(0.5, min(flat) * 0.7))
    ax.grid(True, axis='x', which='both', alpha=0.25, linestyle=':')
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def write_parallelism_tsv(df: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ['dataset', 'version', 'instrument', 'n_runs', 'n_tasks_observed', 'n_invocations', 'wallclock_seconds', 'wallclock_with_lib_seconds', 'insilico_lib_seconds', 'insilico_separable', 'trace_complete', 'source_file']
    cols = [c for c in cols if c in df.columns]
    out = df[cols].copy().sort_values(['dataset', 'version'])
    out.to_csv(tsv_path, sep='\t', index=False)

def write_per_step_tsv(summary: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ['step', 'n', 'p05_seconds', 'p25_seconds', 'p50_seconds', 'p75_seconds', 'p95_seconds', 'min_seconds', 'max_seconds']
    summary[cols].to_csv(tsv_path, sep='\t', index=False)

def cell_line_trace_local_path(dataset: str) -> Path:
    """Local path where a cell-line `nextflow_trace.txt` lands once
    collected from a fresh quantmsdiann rerun (experiment #11). PRIDE
    does not currently publish these files — they're collected with
    `nextflow -with-trace` on the cluster and staged here. If the file
    exists it is included in F2b/F2c/F2a automatically; if not, those
    figures stay benchmarks-only.

    See [docs/superpowers/specs/2026-05-20-experiment-11-cell-line-traces.md]
    for the runbook."""
    return _figure_performance_trace__DATA_DIR / dataset / 'pipeline_info' / 'nextflow_trace.txt'

def has_cell_line_traces() -> bool:
    """Return True iff every cell-line dataset has a local
    `nextflow_trace.txt`. Used by F2c/F2b's collectors to opt in to the
    broader cohort once experiment #11 lands. Until then the
    benchmarks-only scope holds."""
    return all((cell_line_trace_local_path(d).exists() for d in CELL_LINE_ANALYSES))

def collect_step_resource_rows(*, fetch: bool=True) -> tuple[dict[str, dict[str, list[float]]], pd.DataFrame]:
    """Aggregate per-task resource rows across the staged per-dataset fresh
    traces (data/perf_traces/) plus the 5 PXD071075 cluster-node sweep runs.
    Filesystem only; nothing is fetched from FTP. Each sweep point runs the
    same workflow on the same ~2,310 input files, so pooling all 5 grows the
    effective sample size of the per-step resource distribution without bias.

    Returns (resources_by_step, summary_df) where summary_df has
    columns `step, n_rss, peak_rss_p50_bytes, peak_rss_p95_bytes,
    pct_cpu_p50, pct_cpu_p95`.
    """
    traces: list[pd.DataFrame] = [load_trace_resources(p) for p in iter_perf_trace_paths()]
    traces.extend(_load_pxd071075_sweep_resource_rows())
    resources = aggregate_step_resources(traces)
    summary: list[dict] = []
    for step, metrics in resources.items():
        rss = pd.Series(metrics['peak_rss_bytes'])
        cpu = pd.Series(metrics['pct_cpu'])
        pct_cpu_p50 = float(cpu.median()) if len(cpu) else float('nan')
        pct_cpu_p95 = float(cpu.quantile(0.95)) if len(cpu) else float('nan')
        eff_p50 = pct_cpu_p50 / DIANN_THREADS_BASELINE if not pd.isna(pct_cpu_p50) else float('nan')
        eff_p95 = pct_cpu_p95 / DIANN_THREADS_BASELINE if not pd.isna(pct_cpu_p95) else float('nan')
        summary.append({'step': step, 'n_rss': int(len(rss)), 'peak_rss_p50_bytes': float(rss.median()) if len(rss) else float('nan'), 'peak_rss_p95_bytes': float(rss.quantile(0.95)) if len(rss) else float('nan'), 'n_cpu': int(len(cpu)), 'pct_cpu_p50': pct_cpu_p50, 'pct_cpu_p95': pct_cpu_p95, 'thread_efficiency_p50': eff_p50, 'thread_efficiency_p95': eff_p95})
    df = pd.DataFrame(summary).sort_values('peak_rss_p50_bytes', ascending=False).reset_index(drop=True)
    return (resources, df)
DIANN_THREADS_BASELINE = 12

def render_resources_boxplot(resources: dict[str, dict[str, list[float]]], summary: pd.DataFrame, svg_path: Path) -> None:
    """Two-panel horizontal box plot: (left) `peak_rss` per step in GB,
    (right) **threading efficiency** = `%cpu / (12 * 100 %)` per step.
    Step order is fixed by descending median `peak_rss`, mirrored
    across both panels so a reviewer can read the same step across
    both axes.

    The threading efficiency divisor (12) is `DIANN_THREADS_BASELINE`
    — quantmsdiann's `--threads` setting for the heavy DIA-NN steps.
    DIA-NN steps (INSILICO_LIBRARY_GENERATION / PRELIMINARY_ANALYSIS /
    INDIVIDUAL_ANALYSIS) should sit near 100 % here. Glue steps that
    Nextflow allocates 1-2 cores to (SDRF_PARSING, MSstats, raw
    conversion) sit near 1/12 ≈ 8 % on this scale — this is the
    expected pattern for single-threaded helpers, not a quantmsdiann
    deficiency. The note in the figure footer spells this out."""
    steps = summary['step'].tolist()
    if not steps:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, 'no resource rows', ha='center', va='center', transform=ax.transAxes, color='#888888')
        ax.set_axis_off()
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg_path)
        plt.close(fig)
        return
    rss_data = [[b / 1024 ** 3 for b in resources[s]['peak_rss_bytes']] for s in steps]
    cpu_data = [[v / DIANN_THREADS_BASELINE for v in resources[s]['pct_cpu']] for s in steps]
    fig, (ax_rss, ax_cpu) = plt.subplots(nrows=1, ncols=2, figsize=(9.5, max(4.5, 0.55 * len(steps) + 1.5)), sharey=True)
    _box_kw = dict(medianprops=dict(color=fs.OKABE_ITO['vermillion'], linewidth=1.4), boxprops=dict(color=fs.OKABE_ITO['blue']), whiskerprops=dict(color=fs.OKABE_ITO['blue']), capprops=dict(color=fs.OKABE_ITO['blue']))
    ax_rss.boxplot(rss_data, vert=False, widths=0.62, tick_labels=steps, showfliers=False, **_box_kw)
    ax_rss.set_xlabel('Peak RSS per task (GB)', fontsize=10)
    all_rss = [v for vals in rss_data for v in vals if v > 0]
    rss_hi = max(all_rss) if all_rss else 1.0
    rss_lo = min(all_rss) if all_rss else 0.05
    ax_rss.set_xscale('log')
    ax_rss.set_xlim(rss_lo * 0.6, rss_hi * 2.2)
    fs.despine(ax_rss)
    ax_rss.tick_params(axis='both', labelsize=8)
    ax_rss.invert_yaxis()
    for i, step in enumerate(steps, start=1):
        n = len(resources[step]['peak_rss_bytes'])
        if n:
            xpos = max(rss_data[i - 1]) if rss_data[i - 1] else rss_lo
            ax_rss.text(xpos * 1.15, i, f'n={n}', va='center', ha='left', fontsize=7, color='#666666', clip_on=True)
    ax_cpu.boxplot(cpu_data, vert=False, widths=0.62, tick_labels=steps, showfliers=False, **_box_kw)
    ax_cpu.set_xlabel('Threading efficiency (%)', fontsize=10)
    ax_cpu.axvline(100.0, color=fs.OKABE_ITO['bluish_green'], linestyle='--', linewidth=0.9, zorder=1)
    ax_cpu.axvline(100.0 / DIANN_THREADS_BASELINE, color='#bdbdbd', linestyle=':', linewidth=0.8, zorder=1)
    fs.despine(ax_cpu)
    ax_cpu.tick_params(axis='both', labelsize=8)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def write_per_step_resources_tsv(summary: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ['step', 'n_rss', 'peak_rss_p50_bytes', 'peak_rss_p95_bytes', 'n_cpu', 'pct_cpu_p50', 'pct_cpu_p95', 'thread_efficiency_p50', 'thread_efficiency_p95']
    summary[cols].to_csv(tsv_path, sep='\t', index=False)

def figure_performance_trace_main() -> int:
    _figure_performance_trace__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = _figure_performance_trace__FIGURES_DIR / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    par_df = collect_parallelism_rows(fetch=True)
    sweep_df = collect_pxd071075_sweep_rows()
    if not sweep_df.empty:
        full_df = pd.concat([par_df, sweep_df], ignore_index=True, sort=False)
        write_parallelism_tsv(full_df, data_dir / 'parallelism_data.tsv')
        sweep_top = sweep_df[sweep_df['queue_size'] == 300]
        par_df = pd.concat([par_df, sweep_top], ignore_index=True, sort=False)
    else:
        write_parallelism_tsv(par_df, data_dir / 'parallelism_data.tsv')
    render_parallelism_scatter(par_df, _figure_performance_trace__FIGURES_DIR / 'parallelism_vs_wallclock.svg')
    durations, summary = collect_step_runtime_rows(fetch=True)
    write_per_step_tsv(summary, data_dir / 'runtime_per_step.tsv')
    render_per_step_boxplot(durations, summary, _figure_performance_trace__FIGURES_DIR / 'runtime_per_step.svg')
    resources, resources_summary = collect_step_resource_rows(fetch=True)
    write_per_step_resources_tsv(resources_summary, data_dir / 'resources_per_step.tsv')
    render_resources_boxplot(resources, resources_summary, _figure_performance_trace__FIGURES_DIR / 'resources_per_step.svg')
    n_complete = int(par_df['trace_complete'].sum())
    print(f'Wrote {len(par_df)} parallelism rows ({n_complete} complete trace, {len(par_df) - n_complete} hollow / pipeline_report.txt fallback) and {len(summary)} step-summary rows.')
    return 0


# ======================================================================
# inlined from analysis/figure_phospho.py
# ======================================================================

"""Phospho panel - quantms.io DIA-NN reanalysis, standard 2.5.1 vs 2.5.1 Enterprise.

Three panels (all counts at 1% FDR, max_mods=2; see data/phospho/phospho_counts.tsv):
  A  Phosphopeptides  -- distinct phospho-bearing modified precursors
     (UniMod:21 in Modified.Sequence; Lib.Q.Value <= 0.01, no target filter).
  B  Class-I phosphosites -- localized sites from the DIA-NN site report
     (Modification == UniMod:21, localization Probability >= 0.99), unique by
     protein + residue/site.
  C  PXD049692 deposited-vs-reanalysis -- distinct stripped phosphopeptide
     backbones on the identical 10 diaPASEF runs: the originally deposited
     Spectronaut directDIA report vs the quantms.io DIA-NN 2.5.1 Enterprise
     reanalysis. Stripped backbones are the engine-independent metric (the
     Spectronaut report has no localization-probability column for a site-level
     head-to-head). This supersedes the stale serine-only 2.5.0 comparison.

Datasets in A/B: PXD034128 (phospho-enriched, two acquisitions), PXD049692
(NK-cell Fe-NTA diaPASEF). PXD034623 (Galectin-1 DIA) is NOT phospho-enriched
(~280 phosphopeptides vs >20k for PXD034128) so it is recorded in the TSV but
left out of the bar panels, where it would be an unreadable sliver.

Data provenance
---------------
Panels A/B read data/phospho/phospho_counts.tsv, regenerated by
``analysis/make_phospho_tables.py`` from the deposited DIA-NN reports on the
public PRIDE FTP (quantmsdiann-benchmarks/phospho/<dataset>/v<version>/
quant_tables/diann_report{,.site_report}.parquet). Pipeline:
https://github.com/bigbio/quantmsdiann
Panel C: the reanalysis backbone count is from the same reports; the deposited
Spectronaut directDIA count is from the originally deposited *_PH_Report.tsv
(PRIDE archive PXD049692), cached in analysis/figures/PXD049692/counts.tsv.

Run:  python -m analysis.figure_phospho
Out:  analysis/figures/supplementary/supp_phospho.svg
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
fs.apply_house_style()
_figure_phospho__REPO_ROOT = Path(__file__).resolve().parents[1]
COUNTS = _figure_phospho__REPO_ROOT / 'data' / 'phospho' / 'phospho_counts.tsv'
_figure_phospho__OUT = _figure_phospho__REPO_ROOT / 'analysis' / 'figures' / 'supplementary' / 'supp_phospho.svg'
_figure_phospho__VERSIONS = ['2_5_1', '2_5_1_enterprise']
VLABEL = {'2_5_1': '2.5.1', '2_5_1_enterprise': '2.5.1 Enterprise'}
_figure_phospho__VCOL = {v: fs.VERSION_COLORS[v] for v in _figure_phospho__VERSIONS}
BAR_DATASETS = ['PXD034128 biological-study', 'PXD034128 highspeed-DIA', 'PXD049692 NK-phospho']

def _phospho_backbones():
    import pandas as pd
    df = pd.read_csv(FIGURES_DIR_049692 / 'counts.tsv', sep='\t').set_index('metric')
    r = df.loc['phosphopeptides_stripped']
    return (int(r['original_spectronaut']), int(r[[c for c in df.columns if c.startswith('quantmsdiann')][0]]))
FIGURES_DIR_049692 = _figure_phospho__REPO_ROOT / 'analysis' / 'figures' / 'PXD049692'
PXD049692_DEPOSITED, PXD049692_QUANTMSDIANN = _phospho_backbones()

def _short(label: str) -> str:
    pxd, _, rest = label.partition(' ')
    return f'{pxd}\n{rest}'

def _val(df, ds, ver, col):
    s = df[(df['dataset'] == ds) & (df['version'] == ver)]
    return int(s[col].iloc[0]) if len(s) else 0

def _version_panel(ax, df, col, title, ylabel):
    bw = 0.38
    xs = range(len(BAR_DATASETS))
    for k, ver in enumerate(_figure_phospho__VERSIONS):
        vals = [_val(df, d, ver, col) for d in BAR_DATASETS]
        off = (k - 0.5) * bw
        bars = ax.bar([x + off for x in xs], vals, bw, color=_figure_phospho__VCOL[ver], edgecolor='white', linewidth=0.6, label=VLABEL[ver])
        if ver == _figure_phospho__VERSIONS[-1]:
            for x, d, bar in zip(xs, BAR_DATASETS, bars):
                lo, hi = (_val(df, d, _figure_phospho__VERSIONS[0], col), _val(df, d, ver, col))
                if lo:
                    ax.annotate(f'+{round(100 * (hi - lo) / lo)}%', (bar.get_x() + bar.get_width() / 2, hi), textcoords='offset points', xytext=(0, 2), ha='center', va='bottom', fontsize=8, fontweight='bold', color=_figure_phospho__VCOL[ver])
    ax.set_xticks(list(xs))
    ax.set_xticklabels([_short(d) for d in BAR_DATASETS], fontsize=8.5)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fs.kfmt_axis(ax.yaxis)
    fs.despine(ax)

def _deposited_panel(ax):
    vals = [PXD049692_DEPOSITED, PXD049692_QUANTMSDIANN]
    cols = [fs.COMPARISON['original'], fs.COMPARISON['quantmsdiann']]
    bars = ax.bar([0, 1], vals, 0.55, color=cols, edgecolor='white', linewidth=0.6)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f'{v:,}', ha='center', va='bottom', fontsize=9.5)
    pct = round(100 * (vals[1] - vals[0]) / vals[0])
    ax.annotate(f'+{pct}%', (bars[1].get_x() + bars[1].get_width() / 2, vals[1]), textcoords='offset points', xytext=(0, 16), ha='center', va='bottom', fontsize=9, fontweight='bold', color=fs.COMPARISON['quantmsdiann'])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Original\n(Spectronaut\ndirectDIA)', 'quantms.io\n(DIA-NN 2.5.1\nEnterprise)'], fontsize=8.5)
    ax.set_xlim(-0.7, 1.7)
    ax.set_ylim(0, max(vals) * 1.22)
    ax.set_title('PXD049692 — deposited vs reanalysis')
    ax.set_ylabel('phosphopeptide backbones (1% FDR)')
    fs.kfmt_axis(ax.yaxis)
    fs.despine(ax)

def _figure_phospho__render(out: Path) -> Path:
    df = pd.read_csv(COUNTS, sep='\t')
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.3))
    _version_panel(axes[0], df, 'phosphopeptides', 'Phosphopeptides', 'phosphopeptides')
    _version_panel(axes[1], df, 'sites_classI', 'Class-I phosphosites (loc ≥ 0.99)', 'localized sites')
    _deposited_panel(axes[2])
    handles = [Line2D([0], [0], marker='s', linestyle='none', markersize=9, markerfacecolor=_figure_phospho__VCOL[v], markeredgecolor='white', label=f'DIA-NN {VLABEL[v]}') for v in _figure_phospho__VERSIONS]
    handles += [Line2D([0], [0], marker='s', linestyle='none', markersize=9, markerfacecolor=fs.COMPARISON['original'], markeredgecolor='white', label='deposited (original)'), Line2D([0], [0], marker='s', linestyle='none', markersize=9, markerfacecolor=fs.COMPARISON['quantmsdiann'], markeredgecolor='white', label='quantms.io reanalysis')]
    fig.legend(handles=handles, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.03), fontsize=9)
    for a, lab in zip(axes, 'ABC'):
        a.text(-0.16, 1.06, lab, transform=a.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='right')
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out

def figure_phospho_main() -> int:
    print(f'wrote {_figure_phospho__render(_figure_phospho__OUT)}')
    return 0


# ======================================================================
# inlined from analysis/figure_proteobench_accuracy.py
# ======================================================================

"""Fig 1d - ProteoBench quantification-accuracy concordance vs standalone DIA-NN.

A single concordance panel: for the two ProteoBench DIA modules with public
predicted-from-FASTA **DIA-NN** community submissions (Module 7, Orbitrap
Astral; PXD062685, timsTOF diaPASEF), the measured HYE log2 fold-change (Y) is
plotted against the ProteoBench-expected ratio (X), with the dashed Y=X line.
Standalone DIA-NN community runs (ALL versions) are grey; quantms-diann (1.8.1
and 2.5.1-enterprise) are coloured. Points on the diagonal => quantms-diann
quantifies as accurately as single-machine DIA-NN, independent of release.

Only DIA-NN community submissions are used (other ProteoBench tools excluded),
so "standalone DIA-NN" is the literal comparator. Per-version, per-module
identification DEPTH (precursors / protein groups) is in Supplementary Note 5
(figure_id_vs_epsilon.py); this panel is accuracy only.

Source: data/quantmsdiann_benchmarks/proteobench/<dataset>.json (community) and
the quantms-diann metrics cache (per version).

Run:  python -m analysis.figure_proteobench_accuracy
Out:  analysis/figures/manuscript/fig1d_proteobench_accuracy.svg
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
fs.apply_house_style()
REPO_ROOT = Path(__file__).resolve().parents[1]
COMMUNITY_DIR = REPO_ROOT / 'data' / 'quantmsdiann_benchmarks' / 'proteobench'
OUT = REPO_ROOT / 'analysis' / 'figures' / 'manuscript' / 'fig1d_proteobench_accuracy.svg'
SPECIES = ('HUMAN', 'YEAST', 'ECOLI')
GREY = fs.OKABE_ITO['grey']
VERSIONS = ('v1_8_1', 'v2_5_1_enterprise')
MARKER = {'ProteoBench_Module_7': 'o', 'PXD062685': 's'}
MARKER_LABEL = {'ProteoBench_Module_7': 'Astral (Module 7)', 'PXD062685': 'timsTOF diaPASEF'}

def diann_community(dataset: str, threshold: int=3) -> list[dict]:
    """Standalone DIA-NN, predicted-from-FASTA submissions (ALL versions).

    Only DIA-NN entries with DIA-NN-predicted libraries are kept, so other
    ProteoBench tools (AlphaDIA, Spectronaut, ...) and user-supplied-library
    runs are excluded — the grey cloud is literally standalone DIA-NN.
    """
    payload = json.load(open(COMMUNITY_DIR / f'{dataset}.json', encoding='utf-8'))
    out = []
    for e in payload:
        if str(e.get('software_name', '')).strip().upper() != 'DIA-NN':
            continue
        pred = e.get('predictors_library')
        if not (isinstance(pred, dict) and str(pred.get('RT', '')).upper() == 'DIANN'):
            continue
        r = e.get('results', {}).get(str(threshold))
        if isinstance(r, dict):
            out.append(r)
    return out

def draw(ax, threshold: int=3, *, with_legend: bool=True, compact: bool=False, square: bool=True) -> None:
    """Draw the accuracy concordance into `ax` (reused by the Fig 2 row)."""
    datasets = list(_COMMUNITY_COMPARATOR_DATASETS)
    rng = np.random.default_rng(0)
    lab = 9 if compact else None
    for ds in datasets:
        module = DATASET_TO_MODULE[ds]
        expected = SPECIES_EXPECTED_LOG2_A_vs_B.get(module, {})
        comm = diann_community(ds, threshold)
        qm = extract_qm_per_species_log2(ds, threshold)
        qm = qm[qm['version'].isin(VERSIONS)]
        for sp in SPECIES:
            if sp not in expected:
                continue
            x = expected[sp]
            cv = [c.get(f'median_log2_empirical_{sp}') for c in comm]
            cv = [v for v in cv if isinstance(v, (int, float)) and (not np.isnan(v))]
            ax.scatter(x + rng.uniform(-0.05, 0.05, len(cv)), cv, s=22, color=GREY, alpha=0.5, edgecolors='none', zorder=2)
            for _, r in qm[qm['species'] == sp].iterrows():
                ax.scatter(x, r['median_log2_empirical'], s=80, marker=MARKER[ds], color=_VERSION_COLORS.get(r['version'], '#d62728'), edgecolors='white', linewidths=0.8, zorder=3)
    lim = [-2.7, 1.7]
    ax.plot(lim, lim, '--', color='#444444', linewidth=1.1, zorder=1)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    if square:
        ax.set_aspect('equal')
    ax.set_xlabel('Expected log$_2$ ratio (ProteoBench)', fontsize=lab)
    ax.set_ylabel('Observed log$_2$ ratio', fontsize=lab)
    if compact:
        ax.tick_params(labelsize=8)
    fs.despine(ax)
    if with_legend:
        handles = [Line2D([0], [0], linestyle='--', color='#444444', label='Y = X'), Line2D([0], [0], marker='o', linestyle='none', ms=7, markerfacecolor=GREY, markeredgecolor='none', label='standalone DIA-NN')]
        for ver in VERSIONS:
            handles.append(Line2D([0], [0], marker='o', linestyle='none', ms=7, markerfacecolor=_VERSION_COLORS.get(ver, '#d62728'), markeredgecolor='white', label=f'quantms-diann {_VERSION_LABELS.get(ver, ver)}'))
        ax.legend(handles=handles, loc='upper left', fontsize=6.5 if compact else 8, frameon=False, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.3)
SHAPE_BY_VERSION = {'v1_8_1': 'o', 'v2_5_1_enterprise': '^'}

def draw_strip(ax, threshold: int=3, *, compact: bool=False, dataset_colors: dict | None=None) -> None:
    """Per-species accuracy strip+box: one group per HYE species, community runs
    as a jittered grey strip + box (every dot visible), the ProteoBench-expected
    ratio as a dashed line, quantms-diann as large markers. Shows the same
    accuracy/equivalence as the concordance plot but with the community dots
    clearly spread out instead of stacked on the diagonal."""
    datasets = list(_COMMUNITY_COMPARATOR_DATASETS)
    exp: dict[str, float] = {}
    for ds in datasets:
        exp.update(SPECIES_EXPECTED_LOG2_A_vs_B.get(DATASET_TO_MODULE[ds], {}))
    order = sorted([s for s in SPECIES if s in exp], key=lambda s: exp[s])
    lab = 8 if compact else None
    rng = np.random.default_rng(0)
    for x, sp in enumerate(order):
        ax.hlines(exp[sp], x - 0.42, x + 0.42, color='#444444', linestyle='--', linewidth=1.2, zorder=1)
        cv = []
        for ds in datasets:
            cv += [c.get(f'median_log2_empirical_{sp}') for c in diann_community(ds, threshold)]
        cv = [v for v in cv if isinstance(v, (int, float)) and (not np.isnan(v))]
        if cv:
            bp = ax.boxplot([cv], positions=[x], widths=0.62, showfliers=False, patch_artist=True, zorder=2)
            for b in bp['boxes']:
                b.set(facecolor='#eeeeee', edgecolor='#9e9e9e', linewidth=0.8)
            for w in bp['whiskers'] + bp['caps']:
                w.set(color='#9e9e9e', linewidth=0.8)
            for med in bp['medians']:
                med.set(color='#9e9e9e', linewidth=1.0)
            ax.scatter(x + rng.uniform(-0.24, 0.24, len(cv)), cv, s=16, color=GREY, alpha=0.6, edgecolors='white', linewidths=0.3, zorder=3)
        for di, ds in enumerate(datasets):
            qm = extract_qm_per_species_log2(ds, threshold)
            qm = qm[qm['species'] == sp]
            qm = qm[qm['version'].isin(VERSIONS)]
            colour = (dataset_colors or {}).get(ds)
            ds_dx = -0.05 if di == 0 else 0.05
            for _, r in qm.iterrows():
                ver = r['version']
                ver_dx = -0.13 if ver == VERSIONS[0] else 0.13
                ax.scatter(x + ver_dx + ds_dx, r['median_log2_empirical'], s=34, marker=SHAPE_BY_VERSION.get(ver, 'o'), color=colour or _VERSION_COLORS.get(ver, '#d62728'), edgecolors='black', linewidths=0.6, zorder=4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f'{_SPECIES_LABEL[s]}\n(exp. {exp[s]:+g})' for s in order], fontsize=lab)
    ax.set_ylabel('Observed log$_2$ ratio', fontsize=lab)
    ax.set_xlim(-0.6, len(order) - 0.4)
    if compact:
        ax.tick_params(labelsize=8)
    fs.despine(ax)
    handles = [Line2D([0], [0], linestyle='--', color='#444444', label='ProteoBench expected'), Line2D([0], [0], marker='o', linestyle='none', ms=7, markerfacecolor=GREY, markeredgecolor='white', label='standalone DIA-NN')]
    for ver in VERSIONS:
        handles.append(Line2D([0], [0], marker=SHAPE_BY_VERSION.get(ver, 'o'), linestyle='none', ms=7, markerfacecolor='#666666', markeredgecolor='black', label=f'quantms-diann {_VERSION_LABELS.get(ver, ver)}'))
    ax.legend(handles=handles, loc='upper left', fontsize=6.5 if compact else 8, frameon=False, handletextpad=0.3, labelspacing=0.3, ncol=1)

def render(out: Path, threshold: int=3) -> Path:
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    draw(ax, threshold)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out

def figure_proteobench_accuracy_main() -> int:
    print(f'wrote {render(OUT)}')
    return 0


# ======================================================================
# inlined from analysis/figure_pxd004701_sun_vs_quantmsdiann.py
# ======================================================================

"""PXD004701 reanalysis comparison figure (Sun et al. 2023 vs quantmsdiann).

The published 2023 Sun et al. PCT-SWATH analysis of 76 breast cancer cell
lines (Mol Cell Proteomics, doi:10.1016/j.mcpro.2023.100602, PMC10392136)
reports 6,091 SwissProt proteins consistently identified across all samples
after applying a proteotypic-peptide + Global.Q.Value <= 0.01 filter, then
dropping proteins with >90 % missing rate. The PMC supplement is behind a
proof-of-work CAPTCHA so we cannot retrieve the per-cell TNBC / non-TNBC
mapping; the `BC_SUBTYPES` dict below is our reconstruction from the
breast-cancer literature (Heiser 2012, Neve 2006, Lehmann 2011, Cellosaurus).

Outputs (paper-ready, no titles/footers):
- main_comparison.{pdf,png,svg}: 2 conditions x 3 metrics
- supp_proteins_per_subtype.{pdf,png,svg}: per-subtype protein counts
- supp_missing_values_per_run.{pdf,png,svg}: per-run completeness
- counts.tsv: auditable totals + per-subtype numbers
"""
import re
import sys
from dataclasses import dataclass
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
_figure_pxd004701_sun_vs_quantmsdiann__REPO_ROOT = Path(__file__).resolve().parent.parent
_figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR = _figure_pxd004701_sun_vs_quantmsdiann__REPO_ROOT / 'data' / 'PXD004701'
_figure_pxd004701_sun_vs_quantmsdiann__FIGURES_DIR = _figure_pxd004701_sun_vs_quantmsdiann__REPO_ROOT / 'analysis' / 'figures' / 'PXD004701'
FTP_QUANT_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/cell-lines/PXD004701/v2_5_1/quant_tables'

def _stage_from_ftp() -> None:
    """Download the small core reanalysis files from the PRIDE FTP into
    DATA_DIR if not already present (so a fresh checkout can render the
    headline counts without cluster access)."""
    import urllib.request
    _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name in ('diannsummary.log', 'diann_report.pg_matrix.tsv', 'diann_report.pr_matrix.tsv'):
        dest = _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / name
        if dest.exists() and dest.stat().st_size > 0:
            continue
        print(f'  staging {name} from FTP ...', file=sys.stderr)
        with urllib.request.urlopen(f'{FTP_QUANT_BASE}/{name}', timeout=900) as r:
            dest.write_bytes(r.read())
SUN_PROTEINS = 6091
SUN_PROTEINS_RAW = 8952
SUN_PEPTIDES = 90762
SUN_LIBRARY_PRECURSORS = 194899
SUN_LIBRARY_PROTEINS = 10323
SUN_TNBC = 39
SUN_NON_TNBC = 37
SUN_CELL_LINES = 76
BC_SUBTYPES: dict[str, str] = {'184a1': 'normal-like', '184b5': 'normal-like', 'hbl100': 'normal-like', 'mcf10a': 'normal-like', 'mcf12a': 'normal-like', 'bt20': 'TNBC', 'bt549': 'TNBC', 'cal51': 'TNBC', 'cal120': 'TNBC', 'cal148': 'TNBC', 'du4475': 'TNBC', 'evsat': 'TNBC', 'hcc1143': 'TNBC', 'hcc1187': 'TNBC', 'hcc1395': 'TNBC', 'hcc1599': 'TNBC', 'hcc1806': 'TNBC', 'hcc1937': 'TNBC', 'hcc2185': 'TNBC', 'hcc3153': 'TNBC', 'hcc38': 'TNBC', 'hcc70': 'TNBC', 'hdqp1': 'TNBC', 'hs578t': 'TNBC', 'mb157': 'TNBC', 'mdamb157': 'TNBC', 'mdamb231': 'TNBC', 'mdamb436': 'TNBC', 'mdamb453': 'TNBC', 'mdamb468': 'TNBC', 'mfm223': 'TNBC', 'mx1': 'TNBC', 'ocubm': 'TNBC', 'sum102': 'TNBC', 'sum149': 'TNBC', 'sum159': 'TNBC', 'sum190': 'TNBC', 'sum229': 'TNBC', 'macls2': 'TNBC', 'au565': 'non-TNBC', 'bt474': 'non-TNBC', 'bt483': 'non-TNBC', 'cama1': 'non-TNBC', 'efm19': 'non-TNBC', 'efm192a': 'non-TNBC', 'hcc1419': 'non-TNBC', 'hcc1428': 'non-TNBC', 'hcc1569': 'non-TNBC', 'hcc1954': 'non-TNBC', 'hcc202': 'non-TNBC', 'hcc2218': 'non-TNBC', 'hcc2688': 'non-TNBC', 'jimt1': 'non-TNBC', 'kpl1': 'non-TNBC', 'mcf7': 'non-TNBC', 'mdamb134vi': 'non-TNBC', 'mdamb175vii': 'non-TNBC', 'mdamb330': 'non-TNBC', 'mdamb361': 'non-TNBC', 'mdamb415': 'non-TNBC', 'skbr3': 'non-TNBC', 'skbr5': 'non-TNBC', 'skbr7': 'non-TNBC', 'sum185': 'non-TNBC', 'sum225': 'non-TNBC', 'sum44': 'non-TNBC', 'sum52': 'non-TNBC', 'sw527': 'non-TNBC', 't47d': 'non-TNBC', 'uacc3199': 'non-TNBC', 'uacc893': 'non-TNBC', 'zr751': 'non-TNBC', 'zr7530': 'non-TNBC', 'zr75b': 'non-TNBC', '600mpe': 'non-TNBC', 'ly2': 'non-TNBC'}

@dataclass(frozen=True)
class _figure_pxd004701_sun_vs_quantmsdiann__Counts:
    sun_proteins: int
    sun_proteins_raw: int
    sun_peptides: int
    sun_tnbc: int
    sun_non_tnbc: int
    quantmsdiann_proteins_strict: int
    quantmsdiann_proteins_strict_unfiltered: int
    quantmsdiann_proteins_pg_matrix: int
    quantmsdiann_proteins_consistent: int
    quantmsdiann_peptides: int
    quantmsdiann_precursors: int

def parse_sdrf_data_file_to_cell_line(sdrf_path: Path) -> dict[str, str]:
    """Parse `comment[data file]` -> `characteristics[cell line]`, rewriting
    `.wiff`/`.WIFF` -> `.mzML` so DIA-NN matrix column names match.

    Duplicated from the PXD030304 helper (rather than imported) so PXD004701
    stays self-contained — PXD030304 may evolve its column expectations
    independently."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'comment[data file]']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f'SDRF missing required columns: {missing}')
    out: dict[str, str] = {}
    for cell, data_file in zip(df['characteristics[cell line]'], df['comment[data file]']):
        if not isinstance(data_file, str) or not data_file:
            continue
        mzml = re.sub('\\.wiff$', '.mzML', data_file, flags=re.IGNORECASE)
        out[mzml] = cell.strip() if isinstance(cell, str) else cell
    return out

def _compute_or_load_diann_subtype_consistency_filter(cache_path: Path, parquet_source: str | Path, sdrf_path: Path, subtype_dict: dict[str, str], *, qvalue_cutoff: float=0.01, min_detection_fraction: float=0.1) -> dict[str, set[str]]:
    """Side-cache wrapper around
    `proteins_per_subtype_quantmsdiann_consistency_filter`. Streaming 33 GB
    over HTTP takes ~10-20 minutes, so we persist the per-subtype result as
    a small JSON. Delete the JSON to force a fresh stream."""
    import json
    if cache_path.exists() and cache_path.stat().st_size > 0:
        with open(cache_path, encoding='utf-8') as fh:
            payload = json.load(fh)
        return {s: set(vs) for s, vs in payload.items()}
    result = proteins_per_subtype_quantmsdiann_consistency_filter(parquet_source, sdrf_path, subtype_dict, qvalue_cutoff=qvalue_cutoff, min_detection_fraction=min_detection_fraction)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'w', encoding='utf-8') as fh:
        json.dump({s: sorted(v) for s, v in result.items()}, fh)
    return result

def proteins_per_subtype_quantmsdiann_consistency_filter(parquet_source: str | Path, sdrf_path: Path, subtype_dict: dict[str, str], *, qvalue_cutoff: float=0.01, min_detection_fraction: float=0.1, batch_size: int=1000000) -> dict[str, set[str]]:
    """Sun-style two-stage filter applied to quantmsdiann's long-format
    report.

    Stage 1 (FDR): keep precursor rows where `Proteotypic == 1` AND
    `Lib.Q.Value <= qvalue_cutoff` (the methods.md §1 global precursor rule;
    NO Global.Q.Value gate, NO contaminant/target filter) AND whose `Run` maps
    to one of the SDRF cell lines (those rows form the "mapped run space"; runs
    outside this set are ignored entirely).

    Stage 2 (consistency): for each Protein.Group surviving stage 1,
    compute its mapped-run detection fraction (# distinct mapped runs in
    which it appears / total mapped runs). Drop any Protein.Group with
    detection fraction < `min_detection_fraction` (Sun et al. use
    >90 % missing == <10 % detection).

    Stage 3 (subtype aggregation): for the surviving (run, Protein.Group)
    pairs, group by `subtype_dict[cell_line]` and emit per-subtype unions.
    Cell lines mapped to `'unknown'` are excluded from the subtype aggregation
    (they still count toward the consistency-filter denominator).

    Streams the parquet via fsspec's HTTPFileSystem with column projection
    on (`Run`, `Protein.Group`, `Lib.Q.Value`, `Proteotypic`) when given
    an HTTPS URL; only those columns' chunks transit the wire."""
    import pyarrow.parquet as pq
    sdrf_run_to_cell = parse_sdrf_data_file_to_cell_line(sdrf_path)
    sdrf_no_ext: dict[str, str] = {}
    for k, v in sdrf_run_to_cell.items():
        stem = re.sub('\\.(mzML|wiff)$', '', k, flags=re.IGNORECASE)
        sdrf_no_ext[stem] = v
    mapped_runs = set(sdrf_no_ext.keys())
    total_mapped_runs = len(mapped_runs)
    if total_mapped_runs == 0:
        return {}
    pg_to_run_set: dict[str, set[str]] = {}
    cols = ['Run', 'Protein.Group', 'Lib.Q.Value', 'Proteotypic']
    source = str(parquet_source)
    if source.startswith(('http://', 'https://')):
        import fsspec
        fs = fsspec.filesystem('https')
        opener = lambda: fs.open(source, 'rb')
    else:
        opener = lambda: open(source, 'rb')
    with opener() as fh:
        pf = pq.ParquetFile(fh)
        for batch in pf.iter_batches(batch_size=batch_size, columns=cols):
            runs = batch.column('Run').to_pylist()
            pgs = batch.column('Protein.Group').to_pylist()
            lqv = batch.column('Lib.Q.Value').to_pylist()
            prot = batch.column('Proteotypic').to_pylist()
            for r, pg, q, p in zip(runs, pgs, lqv, prot):
                if p != 1 or q is None or q > qvalue_cutoff:
                    continue
                if r not in sdrf_no_ext:
                    continue
                pg_to_run_set.setdefault(pg, set()).add(r)
    threshold = min_detection_fraction * total_mapped_runs
    consistent = {pg: runs for pg, runs in pg_to_run_set.items() if len(runs) >= threshold}
    out: dict[str, set[str]] = {}
    for pg, runs in consistent.items():
        for r in runs:
            cell = sdrf_no_ext.get(r)
            if cell is None:
                continue
            subtype = subtype_dict.get(cell)
            if subtype is None or subtype == 'unknown':
                continue
            out.setdefault(subtype, set()).add(pg)
    return out

def unique_peptides_quantified(pr_matrix_path: Path) -> int:
    """Unique `Stripped.Sequence` count in `pr_matrix.tsv` among rows with at
    least one non-NA per-run value. Mirrors the PXD003539/PXD030304 peptide
    definition.

    pr_matrix.tsv for PXD004701 is ~2 GB; we read in chunks and union the
    set of stripped sequences across chunks."""
    seqs: set[str] = set()
    sample_cols: list[str] | None = None
    for chunk in pd.read_csv(pr_matrix_path, sep='\t', dtype=str, chunksize=50000):
        if sample_cols is None:
            missing = [c for c in _figure_pxd030304_procan_vs_quantmsdiann__PR_METADATA_COLS if c not in chunk.columns]
            if missing:
                raise ValueError(f'pr_matrix missing metadata columns: {missing}')
            sample_cols = [c for c in chunk.columns if c not in _figure_pxd030304_procan_vs_quantmsdiann__PR_METADATA_COLS]
        any_quant = chunk[sample_cols].notna().any(axis=1)
        for seq in chunk.loc[any_quant, 'Stripped.Sequence']:
            if isinstance(seq, str) and seq:
                seqs.add(seq)
    return len(seqs)

def _figure_pxd004701_sun_vs_quantmsdiann__render_main_figure(counts: Counts, svg_path: Path) -> None:
    """Grouped bar chart: Sun et al. 2023 vs quantmsdiann across 3 metrics.

    Metric 1: Protein groups (Sun-style consistency filter) — paper 6,091 vs
    quantmsdiann post-filter union.
    Metric 2: Proteotypic peptides — paper 90,762 vs quantmsdiann unique
    Stripped.Sequence in pr_matrix.tsv.
    Metric 3: Protein groups (strict 1% global protein-group FDR, no
    consistency filter) — Sun 8,952 vs quantmsdiann 7,886, counted from the
    DIA-NN report (Global.PG.Q.Value <= 0.01), not the pg_matrix. Sun's
    library-based search against a curated pan-human spectral library
    (10,323 proteins) explains the residual gap, whereas at the consistency
    filter (Metric 1) quantmsdiann's library-free analysis exceeds Sun.

    Paper-ready: no title, no footer."""
    metrics = ['Protein groups\n(consistency filter,\n$\\leq$90% missing)', 'Proteotypic peptides', 'Protein groups\n(1% FDR, no\nconsistency filter)']
    sun_vals = [counts.sun_proteins, counts.sun_peptides, counts.sun_proteins_raw]
    diann_vals = [counts.quantmsdiann_proteins_consistent, counts.quantmsdiann_peptides, counts.quantmsdiann_proteins_strict]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    bar_width = 0.35
    x = list(range(len(metrics)))
    bars_s = ax.bar([xi - bar_width / 2 for xi in x], sun_vals, width=bar_width, color=fs.COMPARISON['original'], label='Sun et al. 2023')
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals, width=bar_width, color=fs.COMPARISON['quantmsdiann'], label='quantmsdiann (DIA-NN)')
    for bars, vals in ((bars_s, sun_vals), (bars_d, diann_vals)):
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{v:,}', ha='center', va='bottom', fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=9)
    ax.set_ylabel('Count')
    ymax = max(max(sun_vals), max(diann_vals)) * 1.18
    ax.set_ylim(0, ymax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def render_proteins_per_subtype(diann_per_subtype: dict[str, set[str]], svg_path: Path, *, sun_reference: int=SUN_PROTEINS) -> None:
    """Bar chart: quantmsdiann per-subtype protein-group union (TNBC,
    non-TNBC, normal-like). Horizontal reference line at Sun et al.'s 6,091
    global consistency-filtered total — they do not publish per-subtype
    protein counts so a side-by-side per-subtype comparison is not possible.
    Paper-ready: no title, no footer."""
    subtype_order = ['TNBC', 'non-TNBC', 'normal-like']
    labels = [s for s in subtype_order if s in diann_per_subtype]
    vals = [len(diann_per_subtype[s]) for s in labels]
    fig, ax = plt.subplots(figsize=(7, 5))
    x = list(range(len(labels)))
    colors = {'TNBC': '#d62728', 'non-TNBC': '#1f77b4', 'normal-like': '#2ca02c'}
    bars = ax.bar(x, vals, color=[colors.get(s, '#777777') for s in labels], width=0.55)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{v:,}', ha='center', va='bottom', fontsize=10)
    ax.axhline(sun_reference, color=fs.COMPARISON['original'], linestyle='--', linewidth=1, label=f'Sun et al. 2023 global ({sun_reference:,})')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Distinct protein groups\n(consistency filter)')
    ymax = max(max(vals, default=0), sun_reference) * 1.15
    ax.set_ylim(0, ymax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def _figure_pxd004701_sun_vs_quantmsdiann__render_per_run_completeness(diann_per_run: dict[str, float], svg_path: Path, *, run_to_subtype: dict[str, str] | None=None) -> None:
    """Strip plot: one dot per MS run, y = fraction of the
    quantmsdiann union catalog detected in that run, x-jittered within
    a per-subtype lane. Sun et al. publish no per-run completeness
    data, so we cannot overlay an original-paper curve; instead the
    figure carries information by **splitting runs by BC subtype**
    (TNBC / non-TNBC / normal-like), with per-subtype median bars
    and an overall median reference.

    `run_to_subtype` maps the bare run identifier (matrix column name
    with `.mzML` / `.raw` / `.d` stripped) to one of the BC_SUBTYPES
    labels. Runs whose subtype can't be resolved fall into an
    `unmapped` lane so they remain visible on the figure rather than
    being silently dropped. Paper-ready: no title, no footer."""
    import re as _re

    def _bare_run(run_label: str) -> str:
        return _re.sub('\\.(mzML|raw|d)$', '', run_label, flags=_re.IGNORECASE)
    fig, ax = plt.subplots(figsize=(10, 4.8))
    by_subtype: dict[str, list[tuple[str, float]]] = {}
    for run, v in diann_per_run.items():
        subtype = 'unmapped'
        if run_to_subtype:
            subtype = run_to_subtype.get(_bare_run(run), 'unmapped')
        by_subtype.setdefault(subtype, []).append((run, v))
    SUBTYPE_ORDER = ['TNBC', 'non-TNBC', 'normal-like', 'unmapped']
    subtype_palette = {'TNBC': '#d62728', 'non-TNBC': '#1f77b4', 'normal-like': '#2ca02c', 'unmapped': '#9e9e9e'}
    subtypes_present = [s for s in SUBTYPE_ORDER if s in by_subtype and by_subtype[s]]
    n_total = sum((len(v) for v in by_subtype.values()))
    median_overall = pd.Series([v for vals in by_subtype.values() for _, v in vals]).median()
    import numpy as _np
    rng = _np.random.default_rng(0)
    for i, subtype in enumerate(subtypes_present):
        ys = [v for _, v in by_subtype[subtype]]
        xs = rng.uniform(i - 0.28, i + 0.28, size=len(ys))
        ax.scatter(xs, ys, s=22, color=subtype_palette[subtype], edgecolors='#222222', linewidths=0.3, alpha=0.75, label=f'{subtype} (n={len(ys)})')
        med = pd.Series(ys).median()
        ax.hlines(med, i - 0.34, i + 0.34, color=subtype_palette[subtype], linewidth=2.2, zorder=4)
        ax.text(i, med, f'  median={med:.2f}', va='center', ha='left', fontsize=8, color='#222222')
    ax.axhline(median_overall, color='#444444', linestyle='--', linewidth=0.8, label=f'overall median ({median_overall:.2f})', zorder=2)
    ax.set_xticks(range(len(subtypes_present)))
    ax.set_xticklabels(subtypes_present, fontsize=10)
    ax.set_xlabel(f'Breast-cancer subtype (n_total = {n_total} runs)', fontsize=10)
    ax.set_ylabel('Fraction of protein groups\ndetected per run\n(quantmsdiann DIA-NN)', fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(-0.6, len(subtypes_present) - 0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='lower right', frameon=False, fontsize=8)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def _figure_pxd004701_sun_vs_quantmsdiann__write_counts_tsv(counts: Counts, tsv_path: Path, *, diann_per_subtype: dict[str, set[str]] | None=None) -> None:
    """Auditable counts table with inline methodology notes."""
    rows = [('Protein groups (consistency filter)', 'Sun et al. 2023 (paper headline)', counts.sun_proteins, '6,091 SwissProt proteins; proteotypic, Global.Q.Value<=0.01, >=10% detection across samples (paper, Methods)'), ('Protein groups (consistency filter)', 'quantmsdiann (DIA-NN)', counts.quantmsdiann_proteins_consistent, 'consistency-filter union (Proteotypic==1 AND Lib.Q.Value<=0.01, >=10% detection across mapped runs) per methods.md §1; NO contaminant/target filter, NO Global.Q.Value gate'), ('Protein groups (1% FDR, no consistency)', 'Sun et al. 2023 (paper pre-filter)', counts.sun_proteins_raw, '8,952 proteins identified pre-consistency-filter (paper, Methods)'), ('Protein groups (1% FDR, no consistency)', 'quantmsdiann (DIA-NN, Lib.PG.Q.Value)', counts.quantmsdiann_proteins_strict, "report-based: distinct Protein.Group at Lib.PG.Q.Value<=0.01 in diann_report.parquet (DIA-NN 2.5.1) per methods.md §1; no contaminant filter; comparable to Sun's 8,952"), ('Protein groups (1% FDR, no consistency)', 'quantmsdiann (DIA-NN, pg_matrix rows)', counts.quantmsdiann_proteins_pg_matrix, 'count_matrix_rows on diann_report.pg_matrix.tsv (>=1 non-empty sample; zeros counted; no filter) per methods.md §1'), ('Protein groups (1% FDR, no consistency)', 'quantmsdiann (DIA-NN, diannsummary.log)', counts.quantmsdiann_proteins_strict_unfiltered, 'audit baseline: diannsummary.log headline'), ('Proteotypic peptides', 'Sun et al. 2023 (paper headline)', counts.sun_peptides, '90,762 proteotypic peptides under the same consistency filter (paper, Methods)'), ('Proteotypic peptides', 'quantmsdiann (DIA-NN)', counts.quantmsdiann_peptides, 'unique Stripped.Sequence in pr_matrix.tsv among rows with >=1 non-NA sample'), ('Precursors', 'quantmsdiann (DIA-NN, 1% FDR)', counts.quantmsdiann_precursors, 'from diannsummary.log (Target precursors at 1% global q-value)'), ('Spectral library precursors', 'Sun et al. 2023 (pan-human CAL)', SUN_LIBRARY_PRECURSORS, '194,899 library precursors; NOT identified precursors (PXD009597)'), ('Spectral library proteins', 'Sun et al. 2023 (pan-human CAL)', SUN_LIBRARY_PROTEINS, '10,323 SwissProt library proteins (PXD009597)'), ('Cell lines: TNBC', 'Sun et al. 2023 (paper-reported split)', counts.sun_tnbc, "39 TNBC lines in paper's split of the 76 cell lines"), ('Cell lines: non-TNBC', 'Sun et al. 2023 (paper-reported split)', counts.sun_non_tnbc, "37 non-TNBC lines in paper's split of the 76 cell lines")]
    subtype_counts: dict[str, int] = {}
    for cl, st in BC_SUBTYPES.items():
        subtype_counts[st] = subtype_counts.get(st, 0) + 1
    for st in ('TNBC', 'non-TNBC', 'normal-like', 'unknown'):
        if st in subtype_counts:
            rows.append((f'Cell lines: {st} (this work mapping)', 'quantmsdiann (BC_SUBTYPES from Heiser/Neve/Lehmann + Cellosaurus)', subtype_counts[st], 'hardcoded classification; PMC supplement inaccessible behind a CAPTCHA so per-cell paper assignment not retrievable'))
    if diann_per_subtype is not None:
        note = 'per-subtype union of Protein.Group from diann_report.parquet filtered to Proteotypic==1 AND Global.Q.Value<=0.01 AND consistency filter (>=10% detection across mapped runs)'
        for st in ('TNBC', 'non-TNBC', 'normal-like'):
            if st in diann_per_subtype:
                rows.append((f'Per-subtype proteins | {st}', 'quantmsdiann (DIA-NN)', len(diann_per_subtype[st]), note))
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, 'w', encoding='utf-8') as fh:
        fh.write('metric\tsource\tcount\tnote\n')
        for r in rows:
            fh.write('\t'.join((str(x) for x in r)) + '\n')

def figure_pxd004701_sun_vs_quantmsdiann_main() -> int:
    ensure_cell_line_matrices('PXD004701')
    _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR.mkdir(parents=True, exist_ok=True)
    _figure_pxd004701_sun_vs_quantmsdiann__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    _stage_from_ftp()
    log_path = _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / 'diannsummary.log'
    pg_path = _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / 'diann_report.pg_matrix.tsv'
    pr_path = _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / 'diann_report.pr_matrix.tsv'
    sdrf_path = _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / 'PXD004701.sdrf.tsv'
    print('Parsing DIA-NN summary log...')
    pg_log, prec = parse_diann_summary_log(log_path)
    print(f'  protein groups (log, unfiltered): {pg_log:,}  precursors: {prec:,}')
    print('Counting pg_matrix.tsv quantified rows (count_matrix_rows, no filter)...')
    pg_matrix_rows = count_matrix_rows(pg_path, PG_METADATA_COLS)
    print(f'  pg_matrix quantified rows: {pg_matrix_rows:,}')
    print('Counting unique proteotypic peptides in pr_matrix.tsv...')
    pep = unique_peptides_quantified(pr_path)
    print(f'  unique Stripped.Sequence: {pep:,}')
    import json as _json
    with open(_figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / 'diann_report_protein_counts.json', encoding='utf-8') as _fh:
        _rep_prot = _json.load(_fh)
    report_proteins = int(_rep_prot['prot_global'])
    print(f'  report protein groups (Lib.PG.Q.Value<=0.01): {report_proteins:,}')
    print('Computing per-subtype consistency-filtered protein sets (streaming parquet)...')
    diann_per_subtype = _compute_or_load_diann_subtype_consistency_filter(_figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / 'diann_per_subtype_consistency_filter.json', _figure_pxd004701_sun_vs_quantmsdiann__DATA_DIR / 'diann_report.parquet', sdrf_path, BC_SUBTYPES)
    diann_proteins_consistent = set()
    for s, pgs in diann_per_subtype.items():
        diann_proteins_consistent.update(pgs)
        print(f'  {s:<14s} {len(pgs):>6,}')
    print(f'  union across subtypes (Lib.Q.Value rule, no contaminant filter): {len(diann_proteins_consistent):,}')
    counts = _figure_pxd004701_sun_vs_quantmsdiann__Counts(sun_proteins=SUN_PROTEINS, sun_proteins_raw=SUN_PROTEINS_RAW, sun_peptides=SUN_PEPTIDES, sun_tnbc=SUN_TNBC, sun_non_tnbc=SUN_NON_TNBC, quantmsdiann_proteins_strict=report_proteins, quantmsdiann_proteins_strict_unfiltered=pg_log, quantmsdiann_proteins_pg_matrix=pg_matrix_rows, quantmsdiann_proteins_consistent=len(diann_proteins_consistent), quantmsdiann_peptides=pep, quantmsdiann_precursors=prec)
    print('Rendering per-subtype supp figure...')
    render_proteins_per_subtype(diann_per_subtype, _figure_pxd004701_sun_vs_quantmsdiann__FIGURES_DIR / 'supp_proteins_per_subtype.svg')
    print('Computing per-run completeness (quantmsdiann)...')
    diann_per_run = per_run_completeness_quantmsdiann(pg_path)
    print(f'  {len(diann_per_run)} runs')
    sdrf_run_to_cell = parse_sdrf_data_file_to_cell_line(sdrf_path)
    run_to_subtype: dict[str, str] = {}
    for data_file, cell in sdrf_run_to_cell.items():
        bare = re.sub('\\.(mzML|raw|d)$', '', data_file, flags=re.IGNORECASE)
        cell_norm = (cell or '').strip().lower()
        subtype = BC_SUBTYPES.get(cell_norm)
        if subtype:
            run_to_subtype[bare] = subtype
    print('Rendering per-run completeness supp figure (by subtype)...')
    _figure_pxd004701_sun_vs_quantmsdiann__render_per_run_completeness(diann_per_run, _figure_pxd004701_sun_vs_quantmsdiann__FIGURES_DIR / 'supp_missing_values_per_run.svg', run_to_subtype=run_to_subtype)
    print('Writing auditable counts TSV...')
    data_dir = _figure_pxd004701_sun_vs_quantmsdiann__FIGURES_DIR / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    _figure_pxd004701_sun_vs_quantmsdiann__write_counts_tsv(counts, data_dir / 'counts.tsv', diann_per_subtype=diann_per_subtype)
    if pg_log != 8008:
        print(f'WARN: quantmsdiann protein groups (log) {pg_log} != expected 8,008', file=sys.stderr)
    if prec != 104238:
        print(f'WARN: quantmsdiann precursors {prec} != expected 104,238', file=sys.stderr)
    return 0


# ======================================================================
# inlined from analysis/figure_pxd030304_procan_vs_quantmsdiann.py
# ======================================================================

"""PXD030304 reanalysis comparison figure (ProCan-DepMapSanger vs quantmsdiann).

The published 2022 ProCan-DepMapSanger analysis (Gonçalves et al., Cancer Cell,
doi:10.1016/j.ccell.2022.06.010) provides the original 949-cell-line proteomic
map. The deposited DIA-NN long-format report on PRIDE is 237 GB; we instead
use the authors' figshare deposit (doi:10.6084/m9.figshare.19345397) which has
small post-processed per-sample protein matrices.

Outputs (paper-ready, no titles/footers):
- main_comparison.{pdf,png,svg}: 2 conditions x 2 metrics
- supp_proteins_per_tissue.{pdf,png,svg}: per-tissue protein counts
- supp_missing_values_per_run.{pdf,png,svg}: per-run completeness
- supp_venn_protein_accessions.{pdf,png,svg}: protein-set Venn
- counts.tsv: auditable totals + per-tissue numbers
"""
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import numpy as np
import pandas as pd
_figure_pxd030304_procan_vs_quantmsdiann__REPO_ROOT = Path(__file__).resolve().parent.parent
_figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR = _figure_pxd030304_procan_vs_quantmsdiann__REPO_ROOT / 'data' / 'PXD030304'
_figure_pxd030304_procan_vs_quantmsdiann__FIGURES_DIR = _figure_pxd030304_procan_vs_quantmsdiann__REPO_ROOT / 'analysis' / 'figures' / 'PXD030304'
_figure_pxd030304_procan_vs_quantmsdiann__PRIDE_QUANT_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/cell-lines/PXD030304/v2_5_1'
DIANN_SUMMARY_LOG_URL = f'{_figure_pxd030304_procan_vs_quantmsdiann__PRIDE_QUANT_BASE}/quant_tables/diannsummary.log'
DIANN_PG_MATRIX_URL = f'{_figure_pxd030304_procan_vs_quantmsdiann__PRIDE_QUANT_BASE}/quant_tables/diann_report.pg_matrix.tsv'
DIANN_PR_MATRIX_URL = f'{_figure_pxd030304_procan_vs_quantmsdiann__PRIDE_QUANT_BASE}/quant_tables/diann_report.pr_matrix.tsv'
DIANN_PARQUET_URL = f'{_figure_pxd030304_procan_vs_quantmsdiann__PRIDE_QUANT_BASE}/quant_tables/diann_report.parquet'
_figure_pxd030304_procan_vs_quantmsdiann__QUANTMS_SDRF_URL = f'{_figure_pxd030304_procan_vs_quantmsdiann__PRIDE_QUANT_BASE}/sdrf/PXD030304.sdrf.tsv'
FIGSHARE_FILES = {'protein_matrix_8498_replicates.txt': 34411235, 'peptide_counts_per_protein_per_sample.txt': 34411148, 'mapping_file_averaged.txt': 34411133, 'mapping_file_replicates.txt': 34411136}
FIGSHARE_BASE = 'https://ndownloader.figshare.com/files'
PROCAN_PROTEINS = 8498
PROCAN_PROTEINS_STRINGENT = 6692
PROCAN_LIBRARY_PRECURSORS = 144578
PROCAN_LIBRARY_PROTEINS = 12487
PROCAN_MS_RUNS = 6864
PG_METADATA_COLS = ['Protein.Group', 'Protein.Names', 'Genes', 'First.Protein.Description', 'N.Sequences', 'N.Proteotypic.Sequences']

@dataclass(frozen=True)
class _figure_pxd030304_procan_vs_quantmsdiann__Counts:
    procan_proteins: int
    procan_proteins_stringent: int
    quantmsdiann_proteins: int
    quantmsdiann_proteins_unfiltered: int
    quantmsdiann_proteins_pg_matrix: int
    quantmsdiann_proteins_stringent: int
    quantmsdiann_precursors: int
_figure_pxd030304_procan_vs_quantmsdiann__SUMMARY_LOG_PRECURSOR_LINE_RE = re.compile('Target precursors at 1% global q-value:\\s*(\\d+)')

def parse_diann_summary_log(log_path: Path) -> tuple[int, int]:
    """Return (protein_groups, target_precursors) from a DIA-NN summary log."""
    protein_groups = precursors = None
    with open(log_path, encoding='utf-8') as fh:
        for line in fh:
            m = SUMMARY_LOG_PROTEIN_LINE_RE.search(line)
            if m and protein_groups is None:
                protein_groups = int(m.group(1))
                continue
            m = _figure_pxd030304_procan_vs_quantmsdiann__SUMMARY_LOG_PRECURSOR_LINE_RE.search(line)
            if m and precursors is None:
                precursors = int(m.group(1))
    if protein_groups is None:
        raise ValueError("'Protein groups with global q-value <= 0.01: N' not found in log")
    if precursors is None:
        raise ValueError("'Target precursors at 1% global q-value: N' not found in log")
    return (protein_groups, precursors)

def report_global_pg_lib(parquet_source: str | Path, *, qvalue_cutoff: float=0.01, batch_size: int=2000000) -> int:
    """Distinct global protein groups under methods.md §1: ``Protein.Group``
    values with ``Lib.PG.Q.Value <= qvalue_cutoff``. Decoys (``Decoy == 1``)
    dropped; NO contaminant/target filter, NO ``Global.PG.Q.Value`` gate.

    Streams the DIA-NN long-format report with column projection on
    (``Protein.Group``, ``Lib.PG.Q.Value``, ``Decoy``). Accepts a local Path
    or an HTTPS URL (streamed via fsspec)."""
    import pyarrow.parquet as pq
    source = str(parquet_source)
    if source.startswith(('http://', 'https://')):
        import fsspec
        opener = lambda: fsspec.filesystem('https').open(source, 'rb')
    else:
        opener = lambda: open(source, 'rb')
    cols = ['Protein.Group', 'Lib.PG.Q.Value', 'Decoy']
    pgs: set[str] = set()
    with opener() as fh:
        pf = pq.ParquetFile(fh)
        have = set(pf.schema_arrow.names)
        use = [c for c in cols if c in have]
        if 'Lib.PG.Q.Value' not in use or 'Protein.Group' not in use:
            raise ValueError(f'report {source} missing Protein.Group/Lib.PG.Q.Value')
        for batch in pf.iter_batches(batch_size=batch_size, columns=use):
            pg = batch.column('Protein.Group').to_pylist()
            q = batch.column('Lib.PG.Q.Value').to_pylist()
            dec = batch.column('Decoy').to_pylist() if 'Decoy' in use else [0] * len(pg)
            for g, v, d in zip(pg, q, dec):
                if d == 1 or v is None or v > qvalue_cutoff:
                    continue
                pgs.add(g)
    return len(pgs)

def parse_procan_mapping(mapping_path: Path) -> dict[str, str]:
    """Parse the ProCan-DepMapSanger averaged-sample mapping file and return
    `Cell_line -> Tissue_type`. The 28 tissue categories are the paper's
    canonical condition axis."""
    df = pd.read_csv(mapping_path, sep='\t', dtype=str)
    needed = ['Cell_line', 'Tissue_type']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f'ProCan mapping missing columns: {missing}')
    out: dict[str, str] = {}
    for cell, tissue in zip(df['Cell_line'], df['Tissue_type']):
        if isinstance(cell, str) and isinstance(tissue, str):
            out[cell.strip()] = tissue.strip()
    return out

def parse_procan_replicates_mapping(mapping_path: Path, *, exclude_hek293t: bool=True) -> dict[str, str]:
    """Return `Automatic_MS_filename -> Tissue_type` from the per-replicate
    ProCan mapping file. By default the 1,064 HEK293T QC runs (Tissue_type
    == 'Control_HEK293T') are excluded so they don't contaminate per-tissue
    counts."""
    df = pd.read_csv(mapping_path, sep='\t', dtype=str)
    needed = ['Automatic_MS_filename', 'Tissue_type']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f'ProCan replicates mapping missing: {missing}')
    out: dict[str, str] = {}
    for run, tissue in zip(df['Automatic_MS_filename'], df['Tissue_type']):
        if not (isinstance(run, str) and isinstance(tissue, str)):
            continue
        tissue = tissue.strip()
        if exclude_hek293t and tissue.startswith('Control_HEK'):
            continue
        out[run.strip()] = tissue
    return out

def proteins_per_tissue_procan(replicates_matrix_path: Path, replicates_mapping_path: Path, *, chunksize: int=500) -> dict[str, set[str]]:
    """For each ProCan tissue, the union of detected protein IDs across all
    MS runs mapped to that tissue. Applies identical per-MS-run union
    semantics to the quantmsdiann side, so the per-tissue comparison is
    apples-to-apples.

    Detection: cell value in the per-replicate matrix is non-NA. The matrix
    has one row per MS run (`Automatic_MS_filename` in the first column —
    header cell is blank in the file) and one column per protein
    (`<accession>;<name>`). HEK293T QC runs are dropped via the mapping
    file. Chunked-read because the file is ~519 MB."""
    run_to_tissue = parse_procan_replicates_mapping(replicates_mapping_path)
    out: dict[str, set[str]] = {}
    for chunk in pd.read_csv(replicates_matrix_path, sep='\t', dtype=str, chunksize=chunksize):
        first_col = chunk.columns[0]
        chunk = chunk.rename(columns={first_col: 'run'})
        protein_cols = [c for c in chunk.columns if c != 'run']
        chunk['tissue'] = chunk['run'].map(lambda r: run_to_tissue.get(r.strip()) if isinstance(r, str) else None)
        chunk = chunk[chunk['tissue'].notna()]
        if chunk.empty:
            continue
        for tissue, group in chunk.groupby('tissue'):
            any_detected = group[protein_cols].notna().any(axis=0)
            detected = [c for c in protein_cols if bool(any_detected[c])]
            out.setdefault(tissue, set()).update(detected)
    return out

def _figure_pxd030304_procan_vs_quantmsdiann__load_sdrf_data_file_to_cell_line(sdrf_path: Path) -> dict[str, str]:
    """Same logic as the PXD003539 SDRF loader: parse `comment[data file]` ->
    `characteristics[cell line]`, rewriting `.wiff` -> `.mzML` so DIA-NN
    matrix column names match.

    Duplicated here (rather than imported) because the PXD003539 helper sits
    inside that module's namespace and we want PXD030304 to stay isolated."""
    df = pd.read_csv(sdrf_path, sep='\t', dtype=str)
    needed = ['characteristics[cell line]', 'comment[data file]']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f'SDRF missing required columns: {missing}')
    out: dict[str, str] = {}
    for cell, data_file in zip(df['characteristics[cell line]'], df['comment[data file]']):
        if not isinstance(data_file, str) or not data_file:
            continue
        mzml = re.sub('\\.wiff$', '.mzML', data_file)
        out[mzml] = cell
    return out
_figure_pxd030304_procan_vs_quantmsdiann__PR_METADATA_COLS = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Precursor.Id']

def _compute_or_load_diann_procan_filter(cache_path: Path, parquet_source: str, sdrf_path: Path, procan_mapping_path: Path) -> dict[str, set[str]]:
    """Side-cache wrapper around `proteins_per_tissue_quantmsdiann_procan_filter`:
    streaming 33 GB over HTTP takes ~15 minutes, so we persist the per-tissue
    result as a small JSON (tissue -> sorted list of Protein.Group). Delete
    the JSON to force a fresh stream."""
    import json
    if cache_path.exists() and cache_path.stat().st_size > 0:
        with open(cache_path, encoding='utf-8') as fh:
            payload = json.load(fh)
        return {t: set(vs) for t, vs in payload.items()}
    result = proteins_per_tissue_quantmsdiann_procan_filter(parquet_source, sdrf_path, procan_mapping_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'w', encoding='utf-8') as fh:
        json.dump({t: sorted(s) for t, s in result.items()}, fh)
    return result

def proteins_per_tissue_quantmsdiann_procan_filter(parquet_source: str | Path, sdrf_path: Path, procan_mapping_path: Path, *, qvalue_cutoff: float=0.01, batch_size: int=1000000) -> dict[str, set[str]]:
    """ProCan-style filter applied to quantmsdiann's long-format report.

    For each ProCan tissue, the set of Protein.Group values that have any
    proteotypic precursor passing Lib.Q.Value <= `qvalue_cutoff` (the
    methods.md §1 global precursor rule; NO Global.Q.Value gate, NO
    contaminant/target filter) in at least one MS run mapped to that tissue.
    This keeps Gonçalves et al. 2022's study-defined comparison criterion
    (proteotypic peptides, per-tissue union) while using the admissible
    Lib.Q.Value cut-off.

    `parquet_source` is either a local Path or an HTTPS URL. For the URL
    case we stream the parquet via fsspec's HTTPFileSystem with column
    projection on (`Run`, `Protein.Group`, `Lib.Q.Value`, `Proteotypic`);
    only those columns' chunks transit the wire so the 33 GB file never
    needs to be staged locally."""
    import pyarrow.parquet as pq
    sdrf_run_to_cell = _figure_pxd030304_procan_vs_quantmsdiann__load_sdrf_data_file_to_cell_line(sdrf_path)
    sdrf_no_ext = {}
    for k, v in sdrf_run_to_cell.items():
        stem = re.sub('\\.(mzML|wiff)$', '', k, flags=re.IGNORECASE)
        sdrf_no_ext[stem] = v
    cl_to_tissue = parse_procan_mapping(procan_mapping_path)
    out: dict[str, set[str]] = {}
    cols = ['Run', 'Protein.Group', 'Lib.Q.Value', 'Proteotypic']
    source = str(parquet_source)
    if source.startswith(('http://', 'https://')):
        import fsspec
        fs = fsspec.filesystem('https')
        opener = lambda: fs.open(source, 'rb')
    else:
        opener = lambda: open(source, 'rb')
    with opener() as fh:
        pf = pq.ParquetFile(fh)
        for batch in pf.iter_batches(batch_size=batch_size, columns=cols):
            runs = batch.column('Run').to_pylist()
            pgs = batch.column('Protein.Group').to_pylist()
            lqv = batch.column('Lib.Q.Value').to_pylist()
            prot = batch.column('Proteotypic').to_pylist()
            for r, pg, q, p in zip(runs, pgs, lqv, prot):
                if p != 1 or q is None or q > qvalue_cutoff:
                    continue
                cell = sdrf_no_ext.get(r)
                if cell is None:
                    continue
                tissue = cl_to_tissue.get(cell)
                if tissue is None:
                    continue
                out.setdefault(tissue, set()).add(pg)
    return out

def proteins_per_tissue_quantmsdiann(pr_matrix_path: Path, sdrf_path: Path, procan_mapping_path: Path, *, chunksize: int=50000) -> dict[str, set[str]]:
    """For each ProCan tissue, the union of `Protein.Group` values whose
    precursors are quantified in at least one run belonging to that tissue.

    Uses `diann_report.pr_matrix.tsv` (1% precursor + 1% protein-group FDR per
    cell) rather than `pg_matrix.tsv` (5% PG q-value per cell) so we exercise
    DIA-NN's strictest per-cell filter, matching ProCan's per-replicate "any
    precursor identified → protein detected" semantics more closely. Mapping:
    pr_matrix run column -> SDRF cell line -> ProCan Tissue_type.

    Chunked-read because pr_matrix.tsv is ~2 GB for PXD030304."""
    sdrf_run_to_cell = _figure_pxd030304_procan_vs_quantmsdiann__load_sdrf_data_file_to_cell_line(sdrf_path)
    cl_to_tissue = parse_procan_mapping(procan_mapping_path)
    reader = pd.read_csv(pr_matrix_path, sep='\t', dtype=str, chunksize=chunksize)
    out: dict[str, set[str]] = {}
    col_to_tissue: dict[str, str] | None = None
    for chunk in reader:
        if col_to_tissue is None:
            missing = [c for c in _figure_pxd030304_procan_vs_quantmsdiann__PR_METADATA_COLS if c not in chunk.columns]
            if missing:
                raise ValueError(f'pr_matrix missing metadata columns: {missing}')
            sample_cols = [c for c in chunk.columns if c not in _figure_pxd030304_procan_vs_quantmsdiann__PR_METADATA_COLS]
            col_to_tissue = {}
            for col in sample_cols:
                cell = sdrf_run_to_cell.get(col)
                if not cell:
                    continue
                tissue = cl_to_tissue.get(cell)
                if tissue is None:
                    continue
                col_to_tissue[col] = tissue
        protein_groups = chunk['Protein.Group'].tolist()
        for col, tissue in col_to_tissue.items():
            mask = chunk[col].notna()
            bucket = out.setdefault(tissue, set())
            for pg, ok in zip(protein_groups, mask):
                if ok and isinstance(pg, str):
                    bucket.add(pg)
    return out

def per_run_completeness_procan(peptide_counts_path: Path) -> dict[str, float]:
    """For each MS run row in
    `ProCan-DepMapSanger_peptide_counts_per_protein_per_sample.txt`, the
    fraction of proteins with peptide count > 0. Denominator = total proteins
    in the matrix (typically 8,498).

    We stream the file row-by-row because it's 143 MB and we only need a
    count per row (not the full matrix)."""
    n_proteins: int | None = None
    out: dict[str, float] = {}
    with open(peptide_counts_path, encoding='utf-8', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        if not header or header[0] != 'Run':
            raise ValueError(f"peptide counts file first column should be 'Run', got {header[0]!r}")
        n_proteins = len(header) - 1
        for row in reader:
            if not row:
                continue
            run = row[0]
            n_detected = sum((1 for v in row[1:] if v and v.strip() and _is_positive_count(v)))
            out[run] = n_detected / n_proteins if n_proteins else 0.0
    return out

def _is_positive_count(value: str) -> bool:
    """Return True if `value` parses to a number greater than zero."""
    try:
        return float(value) > 0
    except ValueError:
        return False

def per_run_completeness_quantmsdiann(pg_matrix_path: Path) -> dict[str, float]:
    """For each per-run column in pg_matrix.tsv, the fraction of protein-group
    rows that are non-NA. Denominator = total protein groups in the matrix
    (the quantmsdiann pipeline's identified set)."""
    df = pd.read_csv(pg_matrix_path, sep='\t', dtype=str)
    missing = [c for c in PG_METADATA_COLS if c not in df.columns]
    if missing:
        raise ValueError(f'pg_matrix missing metadata columns: {missing}')
    sample_cols = [c for c in df.columns if c not in PG_METADATA_COLS]
    n_total = len(df)
    if n_total == 0:
        return {c: 0.0 for c in sample_cols}
    out: dict[str, float] = {}
    for col in sample_cols:
        out[col] = int(df[col].notna().sum()) / n_total
    return out

def protein_detection_freq_procan(peptide_counts_path: Path) -> list[float]:
    """Per-protein detection frequency for ProCan = (# runs with peptide
    count > 0) / (total runs), one value per protein column. Streamed because
    the matrix is 143 MB; we accumulate a per-protein counter across rows."""
    counts: list[int] | None = None
    n_runs = 0
    with open(peptide_counts_path, encoding='utf-8', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        if not header or header[0] != 'Run':
            raise ValueError(f"peptide counts file first column should be 'Run', got {header[0]!r}")
        counts = [0] * (len(header) - 1)
        for row in reader:
            if not row:
                continue
            n_runs += 1
            for j, v in enumerate(row[1:]):
                if v and v.strip() and _is_positive_count(v):
                    counts[j] += 1
    if not n_runs or counts is None:
        return []
    return [c / n_runs for c in counts]

def protein_detection_freq_quantmsdiann(pg_matrix_path: Path) -> list[float]:
    """Per-protein detection frequency for quantmsdiann = (# non-NA run
    columns) / (total runs), one value per protein-group row of pg_matrix.tsv."""
    df = pd.read_csv(pg_matrix_path, sep='\t', dtype=str)
    sample_cols = [c for c in df.columns if c not in PG_METADATA_COLS]
    n_runs = len(sample_cols)
    if n_runs == 0 or len(df) == 0:
        return []
    detected = df[sample_cols].notna().sum(axis=1)
    return (detected / n_runs).tolist()

def _figure_pxd030304_procan_vs_quantmsdiann__render_main_figure(counts: Counts, svg_path: Path) -> None:
    """Grouped bar chart: 2 conditions x 2 metrics (proteins, ≥2-peptide
    proteins). Paper-ready: no title, no footer."""
    metrics = ['Protein groups', 'Protein groups\n($\\geq$2 unique peptides)']
    procan_vals = [counts.procan_proteins, counts.procan_proteins_stringent]
    diann_vals = [counts.quantmsdiann_proteins, counts.quantmsdiann_proteins_stringent]
    fig, ax = plt.subplots(figsize=(7, 5))
    bar_width = 0.35
    x = list(range(len(metrics)))
    bars_p = ax.bar([xi - bar_width / 2 for xi in x], procan_vals, width=bar_width, color=fs.COMPARISON['original'], label='ProCan-DepMapSanger 2022')
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals, width=bar_width, color=fs.COMPARISON['quantmsdiann'], label='quantmsdiann (DIA-NN)')
    for bars, vals in ((bars_p, procan_vals), (bars_d, diann_vals)):
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{v:,}', ha='center', va='bottom', fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel('Protein groups (1% FDR)')
    ymax = max(max(procan_vals), max(diann_vals)) * 1.18
    ax.set_ylim(0, ymax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def render_proteins_per_tissue(procan_per_tissue: dict[str, set[str]], diann_per_tissue: dict[str, set[str]], svg_path: Path) -> None:
    """Grouped bar chart: 28 tissues x 2 conditions (ProCan vs quantmsdiann).
    Tissues sorted by descending ProCan cell-line count is not directly
    available here; we sort by descending union size across the two
    pipelines (largest tissues first). Paper-ready: no title, no footer."""
    tissues = sorted(set(procan_per_tissue) | set(diann_per_tissue), key=lambda t: -(len(procan_per_tissue.get(t, set())) + len(diann_per_tissue.get(t, set()))))
    procan_vals = [len(procan_per_tissue.get(t, set())) for t in tissues]
    diann_vals = [len(diann_per_tissue.get(t, set())) for t in tissues]
    fig, ax = plt.subplots(figsize=(13, 6.5))
    x = list(range(len(tissues)))
    bar_width = 0.4
    bars_p = ax.bar([xi - bar_width / 2 for xi in x], procan_vals, width=bar_width, color=fs.COMPARISON['original'], label='ProCan-DepMapSanger 2022')
    bars_d = ax.bar([xi + bar_width / 2 for xi in x], diann_vals, width=bar_width, color=fs.COMPARISON['quantmsdiann'], label='quantmsdiann (DIA-NN)')
    ax.set_xticks(x)
    ax.set_xticklabels(tissues, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Distinct protein groups detected')
    ymax = max(max(procan_vals, default=0), max(diann_vals, default=0))
    ax.set_ylim(0, ymax * 1.15 if ymax else 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.18, 1, 1))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def _figure_pxd030304_procan_vs_quantmsdiann__render_per_run_completeness(procan_freqs: list[float], diann_freqs: list[float], svg_path: Path) -> None:
    """Data-completeness curve: number of protein groups detected in at least
    X% of runs, for each pipeline. This is the fair way to show missing values
    — unlike per-run fraction it is NOT confounded by total proteome size, so a
    pipeline that identifies *more* proteins is no longer artificially pushed
    down. Left edge (>=0%) = total identified; right edge (>=100%) = the fully
    complete core present in every run. Paper-ready: no title, no footer."""
    procan = np.asarray(procan_freqs, dtype=float)
    diann = np.asarray(diann_freqs, dtype=float)
    thresholds = np.linspace(0.0, 1.0, 101)
    procan_curve = [int((procan >= t).sum()) for t in thresholds]
    diann_curve = [int((diann >= t).sum()) for t in thresholds]
    x = thresholds * 100.0
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(x, procan_curve, color=fs.COMPARISON['original'], linewidth=1.8, label=f'ProCan-DepMapSanger 2022 ({len(procan):,} protein groups)')
    ax.plot(x, diann_curve, color=fs.COMPARISON['quantmsdiann'], linewidth=1.8, label=f'quantmsdiann (DIA-NN) ({len(diann):,} protein groups)')
    ax.set_xlabel('Detected in ≥ X% of runs (data completeness)')
    ax.set_ylabel('Protein groups')
    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0)
    fs.kfmt_axis(ax.yaxis)
    fs.despine(ax)
    ax.legend(loc='upper right', frameon=False)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def _figure_pxd030304_procan_vs_quantmsdiann__write_counts_tsv(counts: Counts, tsv_path: Path, *, procan_per_tissue: dict[str, set[str]] | None=None, diann_per_tissue: dict[str, set[str]] | None=None) -> None:
    """Auditable counts table. Optionally appends per-tissue protein counts
    (one row per tissue per pipeline) so the supp figure's numbers are also
    machine-readable."""
    rows = [('Protein groups', 'ProCan-DepMapSanger 2022 (paper headline)', counts.procan_proteins, '8,498 proteins quantified at Global.Q.Value <= 0.01 (paper, Results)'), ('Protein groups (>=2 peptides)', 'ProCan-DepMapSanger 2022 (stringent)', counts.procan_proteins_stringent, '6,692 proteins with >=2 supporting peptides (paper, Results)'), ('Protein groups', 'quantmsdiann (DIA-NN, Lib.PG.Q.Value)', counts.quantmsdiann_proteins, 'global rule: distinct Protein.Group at Lib.PG.Q.Value<=0.01 in diann_report.parquet per methods.md §1; no contaminant/target filter, no Global.PG.Q.Value gate'), ('Protein groups', 'quantmsdiann (DIA-NN, pg_matrix rows)', counts.quantmsdiann_proteins_pg_matrix, 'count_matrix_rows on diann_report.pg_matrix.tsv (>=1 non-empty sample; zeros counted; no filter) per methods.md §1'), ('Protein groups', 'quantmsdiann (DIA-NN, diannsummary.log)', counts.quantmsdiann_proteins_unfiltered, "audit baseline: diannsummary.log 'Protein groups with global q-value <= 0.01' line"), ('Protein groups (>=2 peptides)', 'quantmsdiann (DIA-NN, Lib rule)', counts.quantmsdiann_proteins_stringent, '>=2 distinct Stripped.Sequence per global Protein.Group (Lib.PG.Q.Value<=0.01 / Lib.Q.Value<=0.01), no contaminant filter'), ('Precursors', 'quantmsdiann (DIA-NN, 1% FDR)', counts.quantmsdiann_precursors, 'from diannsummary.log (Target precursors at 1% global q-value)'), ('Spectral library precursors', 'ProCan-DepMapSanger 2022', PROCAN_LIBRARY_PRECURSORS, 'library size (paper, STAR Methods); NOT identified precursors'), ('Spectral library proteins', 'ProCan-DepMapSanger 2022', PROCAN_LIBRARY_PROTEINS, 'library protein count (paper, STAR Methods)'), ('MS runs', 'ProCan-DepMapSanger 2022', PROCAN_MS_RUNS, 'paper headline (PRIDE archive lists 6,981)')]
    if procan_per_tissue is not None or diann_per_tissue is not None:
        procan = procan_per_tissue or {}
        diann = diann_per_tissue or {}
        tissues = sorted(set(procan) | set(diann), key=lambda t: -(len(procan.get(t, set())) + len(diann.get(t, set()))))
        note_procan = "per-MS-run union over protein_matrix_8498_replicates.txt (HEK293T QC runs excluded); reflects ProCan's global 1% Global.Q.Value filtering, NOT per-cell strict identification"
        note_diann = 'per-tissue union of Protein.Group from diann_report.parquet filtered to Proteotypic == 1 AND Lib.Q.Value <= 0.01 (methods.md §1 global precursor rule; no Global.Q.Value gate, no contaminant filter; no per-cell quant FDR)'
        for t in tissues:
            rows.append((f'Per-tissue proteins | {t}', 'ProCan-DepMapSanger 2022 (per-replicate)', len(procan.get(t, set())), note_procan))
            rows.append((f'Per-tissue proteins | {t}', 'quantmsdiann (DIA-NN, 1% FDR)', len(diann.get(t, set())), note_diann))
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, 'w', encoding='utf-8') as fh:
        fh.write('metric\tsource\tcount\tnote\n')
        for r in rows:
            fh.write('\t'.join((str(x) for x in r)) + '\n')

def figure_pxd030304_procan_vs_quantmsdiann_main() -> int:
    ensure_cell_line_matrices('PXD030304')
    _figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR.mkdir(parents=True, exist_ok=True)
    _figure_pxd030304_procan_vs_quantmsdiann__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    log_path = _figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR / 'diannsummary.log'
    pg_path = _figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR / 'diann_report.pg_matrix.tsv'
    sdrf_path = _figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR / 'PXD030304.sdrf.tsv'
    import json as _json
    with open(_figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR / 'diann_report_protein_counts.json', encoding='utf-8') as _fh:
        _rep = _json.load(_fh)
    fs_paths = {name: _figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR / name for name in FIGSHARE_FILES}
    print('Parsing DIA-NN summary log...')
    pg_log, prec = parse_diann_summary_log(log_path)
    print(f'  protein groups (log, unfiltered): {pg_log:,}  precursors: {prec:,}')
    print('Counting pg_matrix.tsv quantified rows (count_matrix_rows, no filter)...')
    pg_matrix_rows = count_matrix_rows(pg_path, PG_METADATA_COLS)
    print(f'  pg_matrix quantified rows: {pg_matrix_rows:,}')
    report_proteins = int(_rep['prot_global'])
    diann_stringent = int(_rep['prot_2pep'])
    print(f'  report protein groups (Lib.PG.Q.Value<=0.01): {report_proteins:,}')
    print(f'  report >=2-peptide protein groups (Lib rule): {diann_stringent:,}')
    counts = _figure_pxd030304_procan_vs_quantmsdiann__Counts(procan_proteins=PROCAN_PROTEINS, procan_proteins_stringent=PROCAN_PROTEINS_STRINGENT, quantmsdiann_proteins=report_proteins, quantmsdiann_proteins_unfiltered=pg_log, quantmsdiann_proteins_pg_matrix=pg_matrix_rows, quantmsdiann_proteins_stringent=diann_stringent, quantmsdiann_precursors=prec)
    print('Rendering main figure...')
    _figure_pxd030304_procan_vs_quantmsdiann__render_main_figure(counts, _figure_pxd030304_procan_vs_quantmsdiann__FIGURES_DIR / 'main_comparison.svg')
    print('Computing per-tissue protein sets (ProCan, per-replicate, HEK293T excluded)...')
    procan_per_tissue = proteins_per_tissue_procan(fs_paths['protein_matrix_8498_replicates.txt'], fs_paths['mapping_file_replicates.txt'])
    print(f'  {len(procan_per_tissue)} tissues')
    print('Loading per-tissue protein sets (quantmsdiann, ProCan-style filter)...')
    diann_per_tissue = _compute_or_load_diann_procan_filter(_figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR / 'diann_per_tissue_procan_filter.json', _figure_pxd030304_procan_vs_quantmsdiann__DATA_DIR / 'diann_report.parquet', sdrf_path, fs_paths['mapping_file_averaged.txt'])
    print(f'  {len(diann_per_tissue)} tissues')
    print('Rendering per-tissue supp figure...')
    render_proteins_per_tissue(procan_per_tissue, diann_per_tissue, _figure_pxd030304_procan_vs_quantmsdiann__FIGURES_DIR / 'supp_proteins_per_tissue.svg')
    print('Computing per-protein detection frequency (ProCan)...')
    procan_freqs = protein_detection_freq_procan(fs_paths['peptide_counts_per_protein_per_sample.txt'])
    print(f'  {len(procan_freqs):,} protein groups')
    print('Computing per-protein detection frequency (quantmsdiann)...')
    diann_freqs = protein_detection_freq_quantmsdiann(pg_path)
    print(f'  {len(diann_freqs):,} protein groups')
    print('Rendering data-completeness supp figure...')
    _figure_pxd030304_procan_vs_quantmsdiann__render_per_run_completeness(procan_freqs, diann_freqs, _figure_pxd030304_procan_vs_quantmsdiann__FIGURES_DIR / 'supp_missing_values_per_run.svg')
    print('Writing auditable counts TSV (with per-tissue rows)...')
    data_dir = _figure_pxd030304_procan_vs_quantmsdiann__FIGURES_DIR / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    _figure_pxd030304_procan_vs_quantmsdiann__write_counts_tsv(counts, data_dir / 'counts.tsv', procan_per_tissue=procan_per_tissue, diann_per_tissue=diann_per_tissue)
    print('Per-tissue protein counts (ProCan | quantmsdiann):')
    all_t = sorted(set(procan_per_tissue) | set(diann_per_tissue), key=lambda t: -(len(procan_per_tissue.get(t, set())) + len(diann_per_tissue.get(t, set()))))
    for t in all_t:
        p = len(procan_per_tissue.get(t, set()))
        d = len(diann_per_tissue.get(t, set()))
        print(f'  {t:32s} {p:>6,} | {d:>6,}')
    if pg_log != 9680:
        print(f'WARN: quantmsdiann protein groups (log) {pg_log} != expected 9,680', file=sys.stderr)
    if prec != 156411:
        print(f'WARN: quantmsdiann precursors {prec} != expected 156,411', file=sys.stderr)
    return 0


# ======================================================================
# inlined from analysis/figure_pxd064049_spatial_vs_quantmsdiann.py
# ======================================================================

"""PXD064049 (CHP-212 MYCN Deep Visual Proteomics, diaPASEF) reanalysis:
quantmsdiann (DIA-NN 2.5.1-enterprise, library-free, plain FASTA) versus the originally deposited
DIA-NN 1.8.1 analysis on the identical 12 DVP runs.

The original analysis (PRIDE PXD064049) used DIA-NN 1.8.1 library-free with a
plain human FASTA; quantmsdiann re-ran the same raw files with DIA-NN
2.5.1-enterprise against the same plain (contaminant-only, NO entrapment)
human FASTA, so both sides search the same space. We therefore compare:

  * main_comparison.svg -- precursors and protein groups. Both sides ship only
                           DIA-NN ``*_pr_matrix.tsv`` / ``*_pg_matrix.tsv``
                           (already q-filtered count matrices, no q-value
                           columns), so each number is the reproducible count
                           of quantified matrix ROWS: a row counts if it has
                           >= 1 non-empty quantity, with NO contaminant/target
                           filter and zeros counted (methods.md §1). The newer
                           build recovers more of both (precursors 17,287 ->
                           20,705; protein groups 2,947 -> 3,099). counts.tsv
                           also records the entrapment hit rate (now 0, since
                           the plain FASTA has no entrapment sequences) for
                           audit parity with earlier runs.

Run:  PYTHONPATH=. python -m analysis.figure_pxd064049_spatial_vs_quantmsdiann
"""
import io
import sys
import urllib.request
import zipfile
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
_figure_pxd064049_spatial_vs_quantmsdiann__REPO_ROOT = Path(__file__).resolve().parent.parent
_figure_pxd064049_spatial_vs_quantmsdiann__FIGURES_DIR = _figure_pxd064049_spatial_vs_quantmsdiann__REPO_ROOT / 'analysis' / 'figures' / 'PXD064049'
_figure_pxd064049_spatial_vs_quantmsdiann__CACHE_DIR = _figure_pxd064049_spatial_vs_quantmsdiann__FIGURES_DIR / 'data' / 'cache'
_figure_pxd064049_spatial_vs_quantmsdiann__ORIG_COLOUR = '#9e9e9e'
_figure_pxd064049_spatial_vs_quantmsdiann__QM_COLOUR = '#1e88e5'
_QB = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/spatial/PXD064049/v2_5_1_enterprise/quant_tables'
_ORIG_ZIP = 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD064049/DIANN_results.zip'

def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or dest.stat().st_size == 0:
        with urllib.request.urlopen(url, timeout=600) as r:
            dest.write_bytes(r.read())
    return dest

def _qm_matrix(kind: str) -> Path:
    """quantmsdiann pr/pg matrix from the benchmarks FTP (cached)."""
    name = f'diann_report.{kind}_matrix.tsv'
    return _download(f'{_QB}/{name}', _figure_pxd064049_spatial_vs_quantmsdiann__CACHE_DIR / f'qm_{kind}_matrix.tsv')

def _orig_matrix(kind: str) -> Path:
    """Authors' deposited DIA-NN 1.8.1 pr/pg matrix (cached from the zip)."""
    dest = _figure_pxd064049_spatial_vs_quantmsdiann__CACHE_DIR / f'orig_{kind}_matrix.tsv'
    if not dest.exists() or dest.stat().st_size == 0:
        zip_dest = _download(_ORIG_ZIP, _figure_pxd064049_spatial_vs_quantmsdiann__CACHE_DIR / 'DIANN_results.zip')
        with zipfile.ZipFile(zip_dest) as z:
            member = next((m for m in z.namelist() if m.endswith(f'MYCN_High_Low.{kind}_matrix.tsv')))
            dest.write_bytes(z.read(member))
    return dest

def _entrapment_hit_rate(matrix_path: Path) -> tuple[int, int, float]:
    """(entrapment_passing, target_passing, entrapment_hit_rate_pct): the
    fraction of accepted identifications whose Protein.Group maps to an
    entrapment sequence. This is a direct measure of how many accepted
    groups are entrapment hits; it equals the empirical FDR only when the
    entrapment database is target-sized (1:1 paired entrapment), so we
    report it as an entrapment hit rate rather than a calibrated FDR."""
    pgs = pd.read_csv(matrix_path, sep='\t', usecols=['Protein.Group'], dtype=str)['Protein.Group'].dropna()
    entrap = int(pgs.str.contains('ENTRAP_').sum())
    target = int(pgs.map(is_target_protein_group).sum())
    return (entrap, target, 100.0 * entrap / target if target else 0.0)

def render_main_comparison(or_pr: int, qm_pr: int, or_pg: int, qm_pg: int, svg_path: Path) -> None:
    """Main Fig.~3 panel (d): 2-condition x 2-metric grouped bar chart,
    original (DIA-NN 1.8.1, grey) vs quantmsdiann (DIA-NN 2.5.1-enterprise, blue), for
    precursors and protein groups at 1% FDR. Matches the per-cohort
    `main_comparison` style of the other panels (log y if a metric's
    cross-condition spread exceeds 5x)."""
    conditions = [('Original (DIA-NN 1.8.1)', _figure_pxd064049_spatial_vs_quantmsdiann__ORIG_COLOUR, or_pr, or_pg), ('quantmsdiann (DIA-NN 2.5.1-enterprise)', _figure_pxd064049_spatial_vs_quantmsdiann__QM_COLOUR, qm_pr, qm_pg)]
    metrics = ['Precursors', 'Protein groups']
    bar_width = 0.27
    x = [0, 1]
    offsets = [bar_width * (i - (len(conditions) - 1) / 2.0) for i in range(len(conditions))]
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, (label, color, pr_val, pg_val) in enumerate(conditions):
        values = [pr_val, pg_val]
        bars = ax.bar([xi + offsets[i] for xi in x], values, width=bar_width, color=color, edgecolor='#37474f', label=label)
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height(), f'{val:,}', ha='center', va='bottom', fontsize=9)
    needs_log = any((min(v) > 0 and max(v) / min(v) > 5 for v in ([or_pr, qm_pr], [or_pg, qm_pg])))
    ylabel = 'Count (1% FDR)'
    if needs_log:
        ax.set_yscale('log')
        ylabel += ' (log scale)'
    else:
        ax.set_ylim(0, max(or_pr, qm_pr, or_pg, qm_pg) * 1.18)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', frameon=False, fontsize=9)
    fig.tight_layout()
    _save(fig, svg_path)

def _save(fig, svg_path: Path) -> None:
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def figure_pxd064049_spatial_vs_quantmsdiann_main() -> int:
    _figure_pxd064049_spatial_vs_quantmsdiann__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    qm_pr_t = count_matrix_rows(_qm_matrix('pr'), PR_METADATA)
    qm_pg_t = count_matrix_rows(_qm_matrix('pg'), PG_METADATA)
    or_pr_t = count_matrix_rows(_orig_matrix('pr'), PR_METADATA)
    or_pg_t = count_matrix_rows(_orig_matrix('pg'), PG_METADATA)
    pr_entrap, _, pr_hit = _entrapment_hit_rate(_qm_matrix('pr'))
    pg_entrap, _, pg_hit = _entrapment_hit_rate(_qm_matrix('pg'))
    render_main_comparison(or_pr_t, qm_pr_t, or_pg_t, qm_pg_t, _figure_pxd064049_spatial_vs_quantmsdiann__FIGURES_DIR / 'main_comparison.svg')
    counts = _figure_pxd064049_spatial_vs_quantmsdiann__FIGURES_DIR / 'counts.tsv'
    counts.write_text(f'metric\toriginal_diann181\tquantmsdiann_diann251_enterprise\tqm_entrapment_hits\tqm_entrapment_hit_pct\nprecursors\t{or_pr_t}\t{qm_pr_t}\t{pr_entrap}\t{pr_hit:.3f}\nprotein_groups\t{or_pg_t}\t{qm_pg_t}\t{pg_entrap}\t{pg_hit:.3f}\n')
    print(f'precursors: original={or_pr_t}  quantmsdiann={qm_pr_t}  (entrapment hit rate {pr_hit:.2f}%, {pr_entrap} hits)')
    print(f'protein groups: original={or_pg_t}  quantmsdiann={qm_pg_t}  (entrapment hit rate {pg_hit:.2f}%, {pg_entrap} hits)')
    print(f'wrote {_figure_pxd064049_spatial_vs_quantmsdiann__FIGURES_DIR}/main_comparison.svg + supp_protein_groups.svg')
    return 0


# ======================================================================
# inlined from analysis/figure_queue_size_sweep.py
# ======================================================================

"""F2d — Nextflow queueSize scaling.

Per-run sweep on a single cohort (PXD071075 single-cell early
results, Wang 2025). The x-axis is **Nextflow `executor.queueSize`**
— the maximum number of concurrent Nextflow tasks. Each task is one
SLURM job (single-core process in this configuration), so queueSize
also bounds the total CPU cores in use; `run_metadata.json` in each
sweep dir reports `queue_size == sweep_cores`.

5 sweep points: 10 / 50 / 100 / 200 / 300. The earlier 4-point sweep
(10 / 20 / 100 / 200) had a q20 outlier from cluster contention; the
re-run replaces q20 with q50 and adds q300 for the high end.

When run with no data present, `main()` exits cleanly with a
"no input yet" message; the suite never fails just because the
sweep hasn't been collected yet.

Spec: docs/superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md
"""
import sys
from pathlib import Path
from typing import Iterable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import pandas as pd
_figure_queue_size_sweep__REPO_ROOT = Path(__file__).resolve().parent.parent
SWEEP_DIR = _figure_queue_size_sweep__REPO_ROOT / 'data' / 'queue_size_sweep'
_figure_queue_size_sweep__FIGURES_DIR = _figure_queue_size_sweep__REPO_ROOT / 'analysis' / 'figures' / 'performance'
DEFAULT_QUEUE_SIZES = (10, 50, 100, 200, 300)

def iter_sweep_traces(queue_sizes: Iterable[int]=DEFAULT_QUEUE_SIZES) -> Iterable[tuple[int, Path]]:
    """Yield (queue_size, trace_path) for every sweep point whose
    `nextflow_trace.txt` exists locally. Looks for both `q<N>/` and
    zero-padded `q<NNN>/` directory layouts (the PRIDE deposit uses
    zero-padded `v2_5_0_sweep_010cores` etc.; the staged local layout
    follows that convention). Points without a trace are silently
    skipped — `main()` decides whether the remaining set is sufficient
    to render the figure."""
    for q in queue_sizes:
        candidates = (SWEEP_DIR / f'q{q:03d}' / 'nextflow_trace.txt', SWEEP_DIR / f'q{q}' / 'nextflow_trace.txt')
        for path in candidates:
            if path.exists():
                yield (q, path)
                break

def collect_sweep_rows(queue_sizes: Iterable[int]=DEFAULT_QUEUE_SIZES) -> pd.DataFrame:
    """Per-queueSize row: `(queue_size, wallclock_s, peak_concurrent,
    n_tasks)`. Returns an empty DataFrame with that schema when no
    sweep traces are present.

    The wallclock is computed from
    `max(submit+duration) − min(submit)` across all rows in the trace
    (FAILED retries included — they did occupy slots), matching the
    existing benchmark trace-wallclock semantics, then the one-time
    INSILICO_LIBRARY_GENERATION step duration is subtracted (see
    `insilico_seconds`) so the panel reflects quantification scaling
    rather than the cohort-independent library-prediction cost. The raw
    span is retained in `wallclock_with_lib_s` for the audit TSV."""
    rows: list[dict] = []
    for q, path in iter_sweep_traces(queue_sizes):
        df = load_trace(path)
        wallclock_s = trace_wallclock_seconds(df)
        lib_s = insilico_seconds(df)
        peak, _med = peak_concurrent_tasks(df)
        rows.append({'queue_size': q, 'wallclock_s': max(0.0, float(wallclock_s) - float(lib_s)), 'wallclock_with_lib_s': float(wallclock_s), 'insilico_lib_s': float(lib_s), 'peak_concurrent': int(peak), 'n_tasks': int(len(df))})
    return pd.DataFrame(rows, columns=['queue_size', 'wallclock_s', 'wallclock_with_lib_s', 'insilico_lib_s', 'peak_concurrent', 'n_tasks']).sort_values('queue_size').reset_index(drop=True)

def render_queue_size_sweep(df: pd.DataFrame, svg_path: Path | None=None, *, ax: plt.Axes | None=None, composite: bool=False) -> None:
    """Wallclock vs cluster-node count for the PXD071075 single-cell
    queueSize sweep. Log-log axes with explicit decimal tick labels
    (10/50/100/200/300 — no 10^1 / 10^2 power notation). Dot fill is
    proportional to queueSize so the eye reads the sweep direction at
    a glance. The wallclock value (hours) is printed next to each dot
    for direct readability.

    When `df` has fewer than 2 points the function still emits an SVG
    but with a "needs more sweep points" placeholder.

    Pass `ax` to draw into an existing axes (composite figures); omit
    `svg_path` in that mode."""
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(5.6, 4.4))
    if len(df) < 2:
        ax.text(0.5, 0.5, f'sweep is data-bound on experiment #12 — need at least 2 sweep points, have {len(df)}', transform=ax.transAxes, ha='center', va='center', fontsize=9, color='#888888')
        ax.set_axis_off()
        if own_fig:
            fig.tight_layout()
            assert svg_path is not None
            svg_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(svg_path, bbox_inches='tight')
            plt.close(fig)
        return
    df = df.sort_values('queue_size').reset_index(drop=True)
    xs = df['queue_size'].astype(float).values
    ys = (df['wallclock_s'] / 3600.0).values
    dot_color = '#1976d2'
    dot_size = 90 if composite else 140
    line_width = 1.4 if composite else 1.8
    ann_size = 6 if composite else 8.5
    label_size = 8 if composite else 11
    tick_size = 6.5 if composite else 10
    ax.plot(xs, ys, color=dot_color, linewidth=line_width, alpha=0.55, zorder=2)
    ax.scatter(xs, ys, s=dot_size, c=dot_color, edgecolors='#0d47a1', linewidths=0.8, zorder=3)
    for xi, yi in zip(xs, ys):
        ax.annotate(f'{yi:.1f} h', xy=(xi, yi), xytext=(6, 6), textcoords='offset points', fontsize=ann_size, color='#1a237e', fontweight='bold', ha='left', va='bottom')
    ax.set_xscale('log')
    ax.set_yscale('log')
    from matplotlib.ticker import FixedLocator, FixedFormatter
    ax.xaxis.set_major_locator(FixedLocator(xs))
    ax.xaxis.set_major_formatter(FixedFormatter([str(int(x)) for x in xs]))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    y_ticks = []
    for cand in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0):
        if ys.min() / 1.4 <= cand <= ys.max() * 1.4:
            y_ticks.append(cand)
    if y_ticks:
        ax.yaxis.set_major_locator(FixedLocator(y_ticks))
        ax.yaxis.set_major_formatter(FixedFormatter([f'{t:g}' for t in y_ticks]))
    ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel('Cluster nodes', fontsize=label_size)
    ax.set_ylabel('Workflow wallclock (hours)', fontsize=label_size)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=tick_size)
    ax.set_xlim(xs.min() / 1.4, xs.max() * 1.8)
    ax.set_ylim(ys.min() / 1.6, ys.max() * 1.6)
    if own_fig:
        fig.tight_layout()
        assert svg_path is not None
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg_path, bbox_inches='tight')
        plt.close(fig)

def write_sweep_tsv(df: pd.DataFrame, tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep='\t', index=False)

def figure_queue_size_sweep_main() -> int:
    _figure_queue_size_sweep__FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data_dir = _figure_queue_size_sweep__FIGURES_DIR / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    df = collect_sweep_rows()
    if df.empty:
        print(f'F2d (queue-size sweep) — no trace data under {SWEEP_DIR}/q<N>/nextflow_trace.txt yet; skipping render. See docs/superpowers/specs/2026-05-20-experiment-12-queue-size-scaling.md.')
        return 0
    write_sweep_tsv(df, data_dir / 'queue_size_sweep.tsv')
    render_queue_size_sweep(df, _figure_queue_size_sweep__FIGURES_DIR / 'queue_size_sweep.svg')
    print(f"F2d (queue-size sweep) rendered from {len(df)} sweep point(s) ({list(df['queue_size'])})")
    return 0


# ======================================================================
# inlined from analysis/figure_reanalysis_improvement.py
# ======================================================================

"""Reanalysis-improvement figure: original analysis vs quantmsdiann reanalysis.

Two panels of paired horizontal bars, one row per public DIA deposit:
  (a) Protein groups  — all seven deposits.
  (b) Precursors      — only the deposits whose ORIGINAL analysis was DIA-NN
      (HeLa Astral, spatial DVP); "precursor" is a DIA-NN concept, so
      OpenSWATH/PCT-SWATH/protein-only originals have no comparable count.

Original (grey) = deposit matrix / published headline; reanalysis (coloured by
the DIA-NN version used) = counted from the precursor report (parquet), per the
counting rule Vadim specified. Single-cell PG keeps the pg_matrix union (so the
match-between-runs gain stays visible); NCI-60/ProCan use >=2-peptide; Sun its
consistency filter; plexDIA channel-confident mTRAQ. Provenance per row is in
the data TSV.

Source: analysis/figures/reanalysis/data/reanalysis_improvement.tsv
Out:    analysis/figures/manuscript/fig_reanalysis_improvement.svg
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
fs.apply_house_style()
_figure_reanalysis_improvement__REPO = Path(__file__).resolve().parents[1]
_figure_reanalysis_improvement__DATA = _figure_reanalysis_improvement__REPO / 'analysis' / 'figures' / 'reanalysis' / 'data' / 'reanalysis_improvement.tsv'
_figure_reanalysis_improvement__OUT = _figure_reanalysis_improvement__REPO / 'analysis' / 'figures' / 'manuscript' / 'fig_reanalysis_improvement.svg'
_figure_reanalysis_improvement__ORIG_COLOUR = '#9e9e9e'
VERSION_COLOUR = {'2.5.0': '#90caf9', '2.5.1': '#1976d2', '2.5.1-enterprise': '#ff8f00'}

def _panel(ax, df, orig_col, new_col, xlabel):
    """Paired horizontal bars (grey original / version-coloured reanalysis),
    sorted by gain (largest at top). `df` already filtered to rows with data."""
    df = df.copy()
    df['gain'] = df[new_col] / df[orig_col] - 1.0
    df = df.sort_values('gain', ascending=True).reset_index(drop=True)
    xmax = df[new_col].max()
    bar_h = 0.36
    for i, row in df.iterrows():
        vcol = VERSION_COLOUR.get(row['diann_version'], '#1976d2')
        ax.barh(i + bar_h / 2 + 0.02, row[orig_col], height=bar_h, color=_figure_reanalysis_improvement__ORIG_COLOUR, edgecolor='#37474f', linewidth=0.6, zorder=2)
        ax.barh(i - bar_h / 2 - 0.02, row[new_col], height=bar_h, color=vcol, edgecolor='#37474f', linewidth=0.6, zorder=2)
        ax.text(row[orig_col] + xmax * 0.008, i + bar_h / 2 + 0.02, f'{int(row[orig_col]):,}', va='center', ha='left', fontsize=8, color='#555555')
        ax.text(row[new_col] + xmax * 0.008, i - bar_h / 2 - 0.02, f'{int(row[new_col]):,}', va='center', ha='left', fontsize=8, fontweight='bold', color='#222222')
        ax.text(xmax * 1.17, i, f"+{row['gain'] * 100:.0f}%", va='center', ha='right', fontsize=10, fontweight='bold', color=vcol)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{r['dataset']} ({r['label']})\n{r['original_engine']} → DIA-NN {r['diann_version']}" for _, r in df.iterrows()], fontsize=8.0)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_xlim(0, xmax * 1.2)
    ax.set_ylim(-0.7, len(df) - 0.3)
    fs.despine(ax)
    ax.tick_params(axis='y', length=0)
    return df

def _figure_reanalysis_improvement__render(out: Path) -> Path:
    df = pd.read_csv(_figure_reanalysis_improvement__DATA, sep='\t')
    pg = df.dropna(subset=['reanalysis'])
    prec = df.dropna(subset=['orig_precursors', 'new_precursors']).copy()
    prec['orig_precursors'] = prec['orig_precursors'].astype(int)
    prec['new_precursors'] = prec['new_precursors'].astype(int)
    n_pg, n_pr = (len(pg), len(prec))
    fig = plt.figure(figsize=(8.4, 0.66 * (n_pg + n_pr) + 2.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[n_pg, n_pr], hspace=0.32)
    ax_pg = fig.add_subplot(gs[0])
    ax_pr = fig.add_subplot(gs[1])
    _panel(ax_pg, pg, 'original', 'reanalysis', 'Protein groups')
    _panel(ax_pr, prec, 'orig_precursors', 'new_precursors', 'Precursors')
    ax_pg.set_title('(a) Protein groups — all reanalysed deposits', loc='left', fontsize=11, fontweight='bold')
    ax_pr.set_title('(b) Precursors — deposits with a DIA-NN original', loc='left', fontsize=11, fontweight='bold')
    handles = [Patch(facecolor=_figure_reanalysis_improvement__ORIG_COLOUR, edgecolor='#37474f', label='Original analysis')]
    for v, c in VERSION_COLOUR.items():
        if (df['diann_version'] == v).any():
            handles.append(Patch(facecolor=c, edgecolor='#37474f', label=f'quantmsdiann (DIA-NN {v})'))
    ax_pg.legend(handles=handles, loc='lower right', frameon=False, fontsize=8, title='reanalysis', title_fontsize=8.5)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out

def figure_reanalysis_improvement_main() -> int:
    print(f'wrote {_figure_reanalysis_improvement__render(_figure_reanalysis_improvement__OUT)}')
    return 0


# ======================================================================
# inlined from analysis/figure_single_cell_combined.py
# ======================================================================

"""Single-cell reanalysis figure - DIA-NN 1.8.1 vs 2.5.1 Enterprise, two
label-free single-cell cohorts: HeLa Astral (PXD046357) and A549/H460
(PXD049412; the 20x/40x A549 library runs are excluded from per-cell counts).
2x2 layout:
  A  (top, spanning) total precursors + total protein groups + per-cell protein
     groups (box + jitter), both cohorts x build, with per-build %change.
  B  Data-completeness curve (HeLa Astral flagship; >= N cells), y from 0.
  C  CV across cells -- quantitative precision (both cohorts; Astral solid,
     A549/H460 dashed).

Data provenance
---------------
All numbers are derived from the deposited DIA-NN reports by
``analysis/make_single_cell_tables.py`` (run it to (re)generate the inputs):
  * mv_{per_cell,completeness,rank_abundance,cv}.tsv and sc_totals.tsv
    <- PRIDE FTP quantmsdiann-benchmarks/single-cell/{PXD046357,PXD049412}/
       v{1_8_1,2_5_1_enterprise}/quant_tables/diann_report.{tsv,parquet}
       (our reanalysis; counting via analysis/count_report_ids.py under the
       methods.md filter rule: per-cell PG.Q.Value, totals Lib.*).
Pipeline: https://github.com/bigbio/quantmsdiann

Data: data/single_cell/mv_*.tsv, data/single_cell/sc_totals.tsv.

Run:  python -m analysis.figure_single_cell_combined
Out:  analysis/figures/manuscript/fig3_single_cell_combined.svg
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
fs.apply_house_style()
_figure_single_cell_combined__REPO = Path(__file__).resolve().parents[1]
D = _figure_single_cell_combined__REPO / 'data' / 'single_cell'
_figure_single_cell_combined__OUT = _figure_single_cell_combined__REPO / 'analysis' / 'figures' / 'manuscript' / 'fig3_single_cell_combined.svg'
VERS = ['1_8_1', '2_5_1_enterprise']
VLAB = {'1_8_1': '1.8.1', '2_5_1_enterprise': '2.5.1 Enterprise'}
_figure_single_cell_combined__VCOL = {v: fs.VERSION_COLORS[v] for v in VERS}
_figure_single_cell_combined__ACC = {'HeLa Astral SC': 'PXD046357', 'A549/H460 SC': 'PXD049412'}
ACC_SHORT = {'Astral': 'PXD046357', 'A549/H460': 'PXD049412'}
DS_STYLE = {'HeLa Astral SC': '-', 'A549/H460 SC': '--'}
_figure_single_cell_combined__FLAG = 'HeLa Astral SC'

def _completeness(ax):
    df = pd.read_csv(D / 'mv_completeness.tsv', sep='\t')
    df = df[df['dataset'] == _figure_single_cell_combined__FLAG]
    for v in VERS:
        s = df[df['version'] == v].sort_values('min_cells')
        ax.plot(s['min_cells'], s['n_proteins'], linestyle='-', marker='o', ms=3.5, lw=1.8, color=_figure_single_cell_combined__VCOL[v])
    ax.set_xlabel('quantified in ≥ N cells')
    ax.set_ylabel('protein groups')
    ax.set_ylim(bottom=0)
    ax.set_title('Data completeness')
    fs.kfmt_axis(ax.yaxis)
    fs.despine(ax)

def _cv(ax):
    df = pd.read_csv(D / 'mv_cv.tsv', sep='\t')
    bins = np.linspace(0, 1.5, 46)
    for ds in df['dataset'].unique():
        for v in VERS:
            # Show the distribution on [0, 1.5]; the long high-CV tail of the
            # A549/H460 cohort is dropped from the view (NOT clipped into the
            # last bin, which would create a spurious spike at 1.5). density
            # re-normalises over the in-range values, so the curve shapes stay
            # comparable across cohorts.
            cv = df[(df['dataset'] == ds) & (df['version'] == v)]['cv']
            cv = cv[(cv >= 0) & (cv <= 1.5)]
            ax.hist(cv, bins=bins, density=True, histtype='step', linewidth=1.5, color=_figure_single_cell_combined__VCOL[v], linestyle=DS_STYLE.get(ds, '-'))
    ax.set_xlabel('CV across cells')
    ax.set_ylabel('density')
    ax.set_xlim(0, 1.5)
    ax.set_title('Quantitative precision')
    fs.despine(ax)
_SHORT = {'HeLa Astral SC': 'Astral', 'A549/H460 SC': 'A549/H460'}

def _load_totals() -> dict:
    df = pd.read_csv(D / 'sc_totals.tsv', sep='\t')
    out: dict = {}
    for ds, g in df.groupby('dataset'):
        gv = g.set_index('version')
        out[_SHORT.get(ds, ds)] = {m: (int(gv.loc['1_8_1', m]), int(gv.loc['2_5_1_enterprise', m])) for m in ('precursors', 'proteins')}
    return {k: out[k] for k in ('Astral', 'A549/H460') if k in out}
_FULL = {'Astral': 'HeLa Astral SC', 'A549/H460': 'A549/H460 SC'}

def _merged(ax):
    """Merged A+B panel: three x-sections sharing the figure, 1.8.1 vs 2.5.1
    Enterprise. (i) total precursors (left axis), (ii) total protein groups and
    (iii) per-cell protein groups (box+jitter) both on the right axis (same
    scale). Replaces the separate totals + per-cell panels."""
    bw = 0.36
    ax2 = ax.twinx()
    TOTALS = _load_totals()
    dsx = list(TOTALS)
    percell = pd.read_csv(D / 'mv_per_cell.tsv', sep='\t')
    prec_x = {d: i for i, d in enumerate(dsx)}
    prot_x = {d: i + len(dsx) + 0.6 for i, d in enumerate(dsx)}
    cell_x = {d: i + 2 * len(dsx) + 1.2 for i, d in enumerate(dsx)}
    rng = np.random.default_rng(0)
    for k, v in enumerate(VERS):
        idx = 0 if v == '1_8_1' else 1
        for d in dsx:
            xp = prec_x[d] + (k - 0.5) * bw
            hp = TOTALS[d]['precursors'][idx]
            ax.bar(xp, hp, bw, color=_figure_single_cell_combined__VCOL[v], edgecolor='white', linewidth=0.6)
            if v != '1_8_1':
                lo = TOTALS[d]['precursors'][0]
                ax.annotate(f'{round(100 * (hp - lo) / lo):+d}%', (xp, hp), textcoords='offset points', xytext=(0, 3), ha='center', va='bottom', fontsize=9, fontweight='bold', color=_figure_single_cell_combined__VCOL[v])
            xt = prot_x[d] + (k - 0.5) * bw
            hg = TOTALS[d]['proteins'][idx]
            ax2.bar(xt, hg, bw, color=_figure_single_cell_combined__VCOL[v], edgecolor='white', linewidth=0.6)
            if v != '1_8_1':
                lo = TOTALS[d]['proteins'][0]
                ax2.annotate(f'{round(100 * (hg - lo) / lo):+d}%', (xt, hg), textcoords='offset points', xytext=(0, 3), ha='center', va='bottom', fontsize=9, fontweight='bold', color=_figure_single_cell_combined__VCOL[v])
            xc = cell_x[d] + (k - 0.5) * bw
            vals = percell[(percell['dataset'] == _FULL[d]) & (percell['version'] == v)]['pg_count'].values
            bp = ax2.boxplot([vals], positions=[xc], widths=bw * 0.85, patch_artist=True, showfliers=False)
            fs.style_boxplot(bp, color=_figure_single_cell_combined__VCOL[v])
            ax2.scatter(xc + rng.uniform(-0.07, 0.07, len(vals)), vals, s=9, color=_figure_single_cell_combined__VCOL[v], alpha=0.6, edgecolors='none', zorder=3)
    ax.axvline(len(dsx) - 0.2, color='#cccccc', linewidth=0.8)
    ax.axvline(prot_x[dsx[-1]] + 0.7, color='#cccccc', linewidth=0.8)
    ticks = [prec_x[d] for d in dsx] + [prot_x[d] for d in dsx] + [cell_x[d] for d in dsx]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{ACC_SHORT.get(d, d)}\n({d})' for d in dsx] * 3, fontsize=8.5)
    ax.set_xlim(-0.7, cell_x[dsx[-1]] + 0.7)
    ax.set_ylabel('precursors')
    ax2.set_ylabel('protein groups')
    prec_max = max((TOTALS[d]['precursors'][i] for d in dsx for i in (0, 1)))
    prot_max = max((TOTALS[d]['proteins'][i] for d in dsx for i in (0, 1)))
    cell_max = float(percell[percell['dataset'].isin([_FULL[d] for d in dsx])]['pg_count'].max())
    ax.set_ylim(0, prec_max * 1.15)
    ax2.set_ylim(0, max(prot_max, cell_max) * 1.12)
    fs.kfmt_axis(ax.yaxis)
    fs.kfmt_axis(ax2.yaxis)
    ax.set_title('Total identifications and per-cell protein groups')
    for xs, lab in ((prec_x, 'precursors\n(total)'), (prot_x, 'protein groups\n(total)'), (cell_x, 'protein groups\n(per cell)')):
        ax.text(np.mean([xs[d] for d in dsx]), -0.16, lab, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=10, fontweight='bold')
    for sp in ('top',):
        ax.spines[sp].set_visible(False)
        ax2.spines[sp].set_visible(False)

def _figure_single_cell_combined__render(out: Path) -> Path:
    fig = plt.figure(figsize=(10.5, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.42, wspace=0.28)
    ax_merged = fig.add_subplot(gs[0, :])
    ax_comp = fig.add_subplot(gs[1, 0])
    ax_cv = fig.add_subplot(gs[1, 1])
    _merged(ax_merged)
    _completeness(ax_comp)
    _cv(ax_cv)
    handles = [Line2D([0], [0], color=_figure_single_cell_combined__VCOL[v], marker='o', linewidth=2, markersize=8, label=f'DIA-NN {VLAB[v]}') for v in VERS]
    handles += [Line2D([0], [0], color='#555555', linestyle='-', linewidth=2, label='PXD046357 (HeLa Astral)'), Line2D([0], [0], color='#555555', linestyle='--', linewidth=2, label='PXD049412 (A549/H460)')]
    fig.legend(handles=handles, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.01), fontsize=11)
    for a, lab in zip([ax_merged, ax_comp, ax_cv], 'ABC'):
        a.text(-0.06, 1.05, lab, transform=a.transAxes, fontsize=17, fontweight='bold', va='bottom', ha='right')
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=2.6, w_pad=2.2)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out

def figure_single_cell_combined_main() -> int:
    print(f'wrote {_figure_single_cell_combined__render(_figure_single_cell_combined__OUT)}')
    return 0


# ======================================================================
# inlined from analysis/make_phospho_tables.py
# ======================================================================

"""Generate the phospho supplementary figure's input table from public reports.

Reproducibility generator for ``figure_phospho.py`` (Panels A/B): downloads the
deposited DIA-NN phospho reports from the public PRIDE FTP and computes
phosphopeptide and phosphosite counts, so nothing is hand-entered.

PROVENANCE
==========
Our reanalyses are on the PRIDE FTP under
``quantmsdiann-benchmarks/phospho/<dataset>/v<version>/quant_tables/``:
  PXD034128-biological-study, PXD034128-highspeed-DIA, PXD049692, PXD034623
each with ``diann_report.parquet`` and ``diann_report.site_report.parquet``,
for DIA-NN ``v2_5_1`` and ``v2_5_1_enterprise``.

Pipeline: https://github.com/bigbio/quantmsdiann

Definitions (Vadim filter rule; see methods.md §1)
--------------------------------------------------
* phosphopeptides: distinct ``Modified.Sequence`` carrying ``UniMod:21`` at the
  global precursor rule ``Lib.Q.Value <= 0.01`` (no contaminant/target filter);
* sites_all: phospho sites from the DIA-NN site report (``Modification``
  contains ``UniMod:21``), unique by ``(Protein, Site)``;
* sites_classI: the same restricted to localization ``Probability >= 0.99``.

Output: data/phospho/phospho_counts.tsv (consumed by figure_phospho.py).

Run:  python -m analysis.make_phospho_tables
"""
import sys
import urllib.request
from pathlib import Path
import pandas as pd
_make_phospho_tables__REPO = Path(__file__).resolve().parents[1]
_make_phospho_tables__OUT_DIR = _make_phospho_tables__REPO / 'data' / 'phospho'
_make_phospho_tables__CACHE_DIR = _make_phospho_tables__OUT_DIR / 'cache'
_make_phospho_tables__FTP_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/phospho'
PHOSPHO = 'UniMod:21'
_make_phospho_tables__VERSIONS = ['2_5_1', '2_5_1_enterprise']
DATASETS = {'PXD034128 biological-study': 'PXD034128-biological-study', 'PXD034128 highspeed-DIA': 'PXD034128-highspeed-DIA', 'PXD049692 NK-phospho': 'PXD049692', 'PXD034623 Galectin1': 'PXD034623'}

def _cached(ftp_dir: str, version: str, fname: str) -> Path:
    url = f'{_make_phospho_tables__FTP_BASE}/{ftp_dir}/v{version}/quant_tables/{fname}'
    _make_phospho_tables__CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = _make_phospho_tables__CACHE_DIR / f'{ftp_dir}_v{version}_{fname}'
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f'Downloading {url} (cached) ...', file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + '.part')
    with urllib.request.urlopen(url, timeout=900) as resp, open(tmp, 'wb') as fh:
        while (chunk := resp.read(1 << 20)):
            fh.write(chunk)
    tmp.replace(dest)
    return dest

def _count(ftp_dir: str, version: str) -> tuple[int, int, int]:
    import pyarrow.parquet as pq
    rep = _cached(ftp_dir, version, 'diann_report.parquet')
    r = pq.read_table(rep, columns=['Modified.Sequence', 'Lib.Q.Value', 'Protein.Group']).to_pandas()
    r = r[r['Lib.Q.Value'] <= Q_THRESHOLD]
    phosphopeptides = r.loc[r['Modified.Sequence'].str.contains(PHOSPHO, na=False), 'Modified.Sequence'].nunique()
    site = _cached(ftp_dir, version, 'diann_report.site_report.parquet')
    s = pq.read_table(site, columns=['Protein', 'Site', 'Modification', 'Probability']).to_pandas()
    ph = s[s['Modification'].astype(str).str.contains(PHOSPHO, na=False)]
    sites_all = ph.drop_duplicates(['Protein', 'Site']).shape[0]
    sites_classI = ph[ph['Probability'] >= 0.99].drop_duplicates(['Protein', 'Site']).shape[0]
    return (phosphopeptides, sites_classI, sites_all)

def make_phospho_tables_main() -> int:
    rows = []
    for name, ftp_dir in DATASETS.items():
        for version in _make_phospho_tables__VERSIONS:
            pp, c1, sa = _count(ftp_dir, version)
            rows.append((name, version, pp, c1, sa))
            print(f'{name} {version}: {pp} phosphopeptides, {c1} class-I, {sa} sites')
    df = pd.DataFrame(rows, columns=['dataset', 'version', 'phosphopeptides', 'sites_classI', 'sites_all'])
    _make_phospho_tables__OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(_make_phospho_tables__OUT_DIR / 'phospho_counts.tsv', sep='\t', index=False)
    print(f"wrote {_make_phospho_tables__OUT_DIR / 'phospho_counts.tsv'}")
    return 0


# ======================================================================
# inlined from analysis/make_single_cell_tables.py
# ======================================================================

"""Generate the Fig. 3 single-cell input tables from the PUBLIC DIA-NN reports.

This is the reproducibility generator for ``figure_single_cell_combined.py``:
every number in Fig. 3 is derived here from the deposited DIA-NN reports, so
nothing in the figure is hand-entered.

PROVENANCE
==========
Our reanalyses (this pipeline) are published on the PRIDE FTP:

  https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/single-cell/<PXD>/<version>/quant_tables/diann_report.{parquet,tsv}

    PXD046357  HeLa Astral single-cell (Orbitrap Astral) -> "HeLa Astral SC"

  DIA-NN versions: ``v1_8_1`` (diann_report.tsv) and ``v2_5_1_enterprise``
  (diann_report.parquet).

Pipeline:        https://github.com/bigbio/quantmsdiann
SDRF tooling:    https://github.com/bigbio/sdrf-pipelines (convert-diann)
Counting logic:  analysis/count_report_ids.py (canonical, reused here).

Counting convention (Vadim filter rule; see methods.md §1)
----------------------------------------------------------
* totals are GLOBAL numbers: protein groups = distinct ``Protein.Group`` at
  ``Lib.PG.Q.Value <= 0.01`` (``prot_global``); precursors = distinct
  ``Precursor.Id`` at ``Lib.Q.Value <= 0.01`` (``prec_global``). No
  contaminant/target filter; zeros counted.
* the per-cell and completeness panels are PER-RUN numbers: protein groups per
  run at ``PG.Q.Value <= 0.01`` only (no target/global filter).

Dynamic range / CV use ``PG.MaxLFQ`` of the per-run protein groups (a
quantitative metric, not an identification count).

NOTE: this counts from the *report* (not the ``*_matrix.tsv`` files); the
matrices bake in ``--matrix-spec-q`` at a version-dependent run q-value and are
not comparable across versions (see count_report_ids.py docstring).

Outputs (data/single_cell/): mv_per_cell.tsv, mv_completeness.tsv,
mv_rank_abundance.tsv, mv_cv.tsv, sc_totals.tsv.

Run:  python -m analysis.make_single_cell_tables
      (downloads ~0.5 GB of reports once; cached under data/single_cell/cache/)
"""
import sys
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
_make_single_cell_tables__REPO = Path(__file__).resolve().parents[1]
_make_single_cell_tables__OUT_DIR = _make_single_cell_tables__REPO / 'data' / 'single_cell'
_make_single_cell_tables__CACHE_DIR = _make_single_cell_tables__OUT_DIR / 'cache'
_make_single_cell_tables__FTP_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/single-cell'
_make_single_cell_tables__ACC = {'HeLa Astral SC': 'PXD046357', 'A549/H460 SC': 'PXD049412'}
FTP_DIR = {'HeLa Astral SC': 'PXD046357', 'A549/H460 SC': 'PXD049412'}
_make_single_cell_tables__VERSIONS = ['1_8_1', '2_5_1_enterprise']
_make_single_cell_tables__FLAG = 'HeLa Astral SC'
_COLS = ['Run', 'Precursor.Id', 'Protein.Group', 'Q.Value', 'PG.Q.Value', 'Lib.Q.Value', 'Lib.PG.Q.Value', 'PG.MaxLFQ', 'Decoy']

def _report_url(ftp_dir: str, version: str) -> str:
    ext = 'tsv' if version == '1_8_1' else 'parquet'
    # Both single-cell cohorts live under the single-cell/ subtree on the FTP.
    return f'{_make_single_cell_tables__FTP_BASE}/{ftp_dir}/v{version}/quant_tables/diann_report.{ext}'

def _make_single_cell_tables___cached_report(ftp_dir: str, version: str) -> Path:
    """Download the deposited DIA-NN report once and cache it on disk."""
    url = _report_url(ftp_dir, version)
    _make_single_cell_tables__CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = _make_single_cell_tables__CACHE_DIR / f'{ftp_dir}_v{version}_{Path(url).name}'
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f'Downloading {url} (cached) ...', file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + '.part')
    with urllib.request.urlopen(url, timeout=900) as resp, open(tmp, 'wb') as fh:
        while (chunk := resp.read(1 << 20)):
            fh.write(chunk)
    tmp.replace(dest)
    return dest

def _load(path: Path) -> pd.DataFrame:
    if path.suffix == '.parquet':
        import pyarrow.parquet as pq
        have = set(pq.ParquetFile(path).schema_arrow.names)
        return pq.read_table(path, columns=[c for c in _COLS if c in have]).to_pandas()
    return pd.read_csv(path, sep='\t', usecols=lambda c: c in _COLS, low_memory=False)

def _perrun_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """Per-run protein-group rows under the Vadim rule: PG.Q.Value <= 1% only
    (no contaminant/target filter, no global filter). Decoys dropped."""
    if 'Decoy' in df.columns:
        df = df[df['Decoy'] == 0]
    return df[df['PG.Q.Value'] <= Q_THRESHOLD]

def build() -> dict[str, pd.DataFrame]:
    per_cell, completeness, rank, cv, totals = ([], [], [], [], [])
    for ds, acc in _make_single_cell_tables__ACC.items():
        for version in _make_single_cell_tables__VERSIONS:
            _sc_report = _make_single_cell_tables___cached_report(FTP_DIR[ds], version)
            df = _load(_sc_report)
            df = df[~df['Run'].astype(str).str.contains('20xSC|40xSC', case=False, regex=True)]
            prot = _perrun_proteins(df)
            pgrun = prot.drop_duplicates(['Run', 'Protein.Group'])
            c = count_report(df, precursor_q=PRECURSOR_Q.get(f'v{version}', DEFAULT_PRECURSOR_Q))
            totals.append((ds, version, c['prec_global'], c['prot_global']))
            for run, g in pgrun.groupby('Run'):
                per_cell.append((ds, version, int(g['Protein.Group'].nunique())))
            n_runs = pgrun['Run'].nunique()
            seen = pgrun.groupby('Protein.Group')['Run'].nunique()
            for mc in range(1, n_runs + 1):
                completeness.append((ds, version, mc, int((seen >= mc).sum())))
            q = pgrun.dropna(subset=['PG.MaxLFQ']).copy()
            q['PG.MaxLFQ'] = pd.to_numeric(q['PG.MaxLFQ'], errors='coerce')
            q = q[q['PG.MaxLFQ'] > 0]
            agg = q.groupby('Protein.Group')['PG.MaxLFQ'].agg(['mean', 'std', 'count'])
            agg = agg[agg['count'] >= 3]
            for _, r in agg.iterrows():
                if r['mean'] > 0 and (not np.isnan(r['std'])):
                    cv.append((ds, version, float(r['std'] / r['mean'])))
            if ds == _make_single_cell_tables__FLAG:
                mean_int = q.groupby('Protein.Group')['PG.MaxLFQ'].mean().sort_values(ascending=False)
                for i, val in enumerate(mean_int.values, start=1):
                    if i == 1 or i % 10 == 0:
                        rank.append((version, i, float(np.log10(val))))
            discard_download(_sc_report)  # bound disk: drop the report once tabulated
    return {'mv_per_cell.tsv': pd.DataFrame(per_cell, columns=['dataset', 'version', 'pg_count']), 'mv_completeness.tsv': pd.DataFrame(completeness, columns=['dataset', 'version', 'min_cells', 'n_proteins']), 'mv_rank_abundance.tsv': pd.DataFrame(rank, columns=['version', 'rank', 'log10_intensity']), 'mv_cv.tsv': pd.DataFrame(cv, columns=['dataset', 'version', 'cv']), 'sc_totals.tsv': pd.DataFrame(totals, columns=['dataset', 'version', 'precursors', 'proteins'])}

def make_single_cell_tables_main() -> int:
    _make_single_cell_tables__OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, frame in build().items():
        frame.to_csv(_make_single_cell_tables__OUT_DIR / name, sep='\t', index=False)
        print(f'wrote {_make_single_cell_tables__OUT_DIR / name} ({len(frame)} rows)')
    return 0


# ======================================================================
# inlined from analysis/plexDIA/figure_msv000093870_galatidou_vs_quantmsdiann.py
# ======================================================================

"""MSV000093870 plexDIA: quantmsdiann reanalysis vs Galatidou et al. 2024.

The original analysis (Galatidou, Petelski et al., Mol Hum Reprod 2024,
doi:10.1093/molehr/gaae023; code at github.com/SlavovLab/single_cell_oocyte)
quantified single human oocytes by mTRAQ plexDIA (mPOP sample prep, Q
Exactive) and published a post-QC proteins x oocyte matrix
(2022_07_28_Oocyte_ProteinsXcells.csv): relative protein abundance,
mean-normalised across successfully quantified oocytes.

This script compares that published matrix against the quantmsdiann
reanalysis (DIA-NN 2.5.0, mTRAQ 3-channel, channel-confident + conservative
contaminant filter; see figure_msv000093870_oocyte_plexdia.py) on three axes:
protein groups per single cell, total protein groups, and protein-accession
overlap.

Inputs are cached on disk:
  - quantmsdiann report parquet (PRIDE quantmsdiann-benchmarks deposition)
  - original matrix (SlavovLab/single_cell_oocyte GitHub)

Outputs (paper-ready, no titles/footers):
  analysis/figures/plexDIA/MSV000093870/
    main_galatidou_comparison.{svg,pdf,png}
    comparison_counts.tsv
"""
import re
import sys
import urllib.request
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import numpy as np
import pandas as pd
ORIGINAL_MATRIX_URL = 'https://raw.githubusercontent.com/SlavovLab/single_cell_oocyte/main/2022_07_28_Oocyte_ProteinsXcells.csv'
_figure_msv000093870_galatidou_vs_quantmsdiann__ORIG_COLOUR = '#9e9e9e'
_figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR = '#1e88e5'

def _strip_isoform(accession: str) -> str:
    return str(accession).split('-')[0]

def _cached_original() -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / 'Galatidou_2024_ProteinsXcells.csv'
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f'Downloading original matrix {ORIGINAL_MATRIX_URL} (cached)…', file=sys.stderr)
    with urllib.request.urlopen(ORIGINAL_MATRIX_URL, timeout=120) as resp:
        dest.write_bytes(resp.read())
    return dest

def original_per_cell(matrix_path: Path) -> tuple[pd.Series, set[str]]:
    """Proteins quantified per oocyte and the set of accessions, from the
    published proteins x cell matrix (non-missing entries = quantified)."""
    o = pd.read_csv(matrix_path)
    cells = o.columns[1:]
    per_cell = o[cells].notna().sum(axis=0)
    accessions = {_strip_isoform(a) for a in o['leading.protein'].dropna().astype(str)}
    accessions.discard('')
    return (per_cell, accessions)
_RUN_PREFIX_RE = re.compile('(wAP\\d+)')
_ORIG_COL_RE = re.compile('(wAP\\d+)_d(\\d)_')

def _run_prefix(run: str) -> str:
    m = _RUN_PREFIX_RE.match(str(run))
    return m.group(1) if m else str(run)

def qc_matched_table(matrix_path: Path, qm_cells: pd.DataFrame) -> pd.DataFrame:
    """Per-oocyte protein counts for the cells retained by BOTH analyses.

    Each oocyte is keyed by (run prefix, mTRAQ channel). The original matrix
    columns (e.g. ``wAP0021_d0_1A``) give the original count (non-missing
    entries); quantmsdiann gives its count for the same (run, channel) cell.
    Restricting to the intersection removes the cell-QC confound (the original
    dropped low-quality oocytes that quantmsdiann's channel filter retained).
    """
    o = pd.read_csv(matrix_path)
    orig = {}
    for col in o.columns[1:]:
        m = _ORIG_COL_RE.match(col)
        if m:
            orig[m.group(1), m.group(2)] = int(o[col].notna().sum())
    qm = qm_cells.copy()
    qm['key'] = list(zip(qm['Run'].map(_run_prefix), qm['Channel'].astype(str)))
    rows = []
    for _, r in qm.iterrows():
        if r['key'] in orig:
            rows.append({'run': r['key'][0], 'channel': r['key'][1], 'orig_proteins': orig[r['key']], 'qm_proteins': int(r['proteins'])})
    return pd.DataFrame(rows)

def render_qc_matched(matched: pd.DataFrame, svg_path: Path) -> None:
    """Paired per-oocyte comparison on the QC-matched cohort."""
    r = float(np.corrcoef(matched['orig_proteins'], matched['qm_proteins'])[0, 1])
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.3))
    ax = axes[0]
    ax.scatter(matched['orig_proteins'], matched['qm_proteins'], s=22, c=_figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR, alpha=0.7, edgecolors='#37474f', linewidths=0.5)
    lo = 0
    hi = max(matched['orig_proteins'].max(), matched['qm_proteins'].max()) * 1.05
    ax.plot([lo, hi], [lo, hi], color=fs.COMPARISON['original'], linestyle='--', linewidth=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel('Galatidou 2024 proteins / oocyte')
    ax.set_ylabel('quantmsdiann proteins / oocyte')
    ax.text(0.04, 0.96, f'(a)  n={len(matched)} matched oocytes\nPearson r = {r:.2f}', transform=ax.transAxes, va='top', fontsize=9)
    ax = axes[1]
    data = [matched['orig_proteins'].values, matched['qm_proteins'].values]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color='#212121', linewidth=1.4))
    for patch, c in zip(bp['boxes'], (_figure_msv000093870_galatidou_vs_quantmsdiann__ORIG_COLOUR, _figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR)):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
        patch.set_edgecolor('#37474f')
    rng = np.random.default_rng(0)
    for i, vals in enumerate(data, start=1):
        ax.scatter(rng.normal(i, 0.05, size=len(vals)), vals, s=10, c='#37474f', alpha=0.4, linewidths=0, zorder=3)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Galatidou 2024', 'quantmsdiann'])
    ax.set_ylabel('Protein groups per oocyte (matched cohort)')
    ax.set_ylim(bottom=0)
    ax.text(0.5, 0.97, '(b)', transform=ax.transAxes, fontweight='bold', va='top')
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    return r

def quantms_accessions(confident: pd.DataFrame) -> set[str]:
    acc: set[str] = set()
    for ids in confident['Protein.Ids'].dropna().unique():
        for a in str(ids).split(';'):
            acc.add(_strip_isoform(a))
    acc.discard('')
    return acc

def render_comparison(orig_per_cell: pd.Series, qm_cells: pd.DataFrame, orig_total: int, qm_total: int, shared: int, orig_only: int, qm_only: int, svg_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3))
    ax = axes[0]
    data = [orig_per_cell.values, qm_cells['proteins'].values]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color='#212121', linewidth=1.4))
    for patch, c in zip(bp['boxes'], (_figure_msv000093870_galatidou_vs_quantmsdiann__ORIG_COLOUR, _figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR)):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
        patch.set_edgecolor('#37474f')
    rng = np.random.default_rng(0)
    for i, (vals, c) in enumerate(zip(data, (_figure_msv000093870_galatidou_vs_quantmsdiann__ORIG_COLOUR, _figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR)), start=1):
        ax.scatter(rng.normal(i, 0.05, size=len(vals)), vals, s=10, c='#37474f', alpha=0.4, linewidths=0, zorder=3)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'Galatidou 2024\n(n={len(orig_per_cell)})', f'quantmsdiann\n(n={len(qm_cells)})'])
    ax.set_ylabel('Protein groups per single cell')
    ax.set_ylim(bottom=0)
    ax.text(0.5, 0.97, '(a)', transform=ax.transAxes, fontweight='bold', va='top')
    ax = axes[1]
    bars = ax.bar([0, 1], [orig_total, qm_total], width=0.6, color=[_figure_msv000093870_galatidou_vs_quantmsdiann__ORIG_COLOUR, _figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR], edgecolor='#37474f')
    for b, v in zip(bars, (orig_total, qm_total)):
        ax.text(b.get_x() + b.get_width() / 2, v, f'{v:,}', ha='center', va='bottom', fontsize=10)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Galatidou 2024', 'quantmsdiann'])
    ax.set_ylabel('Total protein groups (dataset)')
    ax.set_ylim(0, max(orig_total, qm_total) * 1.15)
    ax.text(0.5, 0.97, '(b)', transform=ax.transAxes, fontweight='bold', va='top')
    ax = axes[2]
    cats = ['Shared', 'Galatidou\nonly', 'quantmsdiann\nonly']
    vals = [shared, orig_only, qm_only]
    cols = ['#43a047', _figure_msv000093870_galatidou_vs_quantmsdiann__ORIG_COLOUR, _figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR]
    bars = ax.bar(range(3), vals, width=0.6, color=cols, edgecolor='#37474f')
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f'{v:,}', ha='center', va='bottom', fontsize=10)
    ax.set_xticks(range(3))
    ax.set_xticklabels(cats)
    ax.set_ylabel('Protein groups')
    ax.set_ylim(0, max(vals) * 1.15)
    ax.text(0.5, 0.97, '(c)', transform=ax.transAxes, fontweight='bold', va='top')
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def render_total_pg(orig_total: int, qm_total: int, svg_path: Path) -> None:
    """Compact single-panel total-protein-group comparison for the main
    cell-line-reanalysis figure (Fig 3): original (grey) vs quantmsdiann
    (blue), matching the per-cohort `main_comparison` bar style. Conveys
    the single-cell reanalysis benefit (more protein groups from the same
    raw data) in one column-width panel. Uses the same (7, 5) canvas as the
    per-cohort `main_comparison` panels so all three Fig. 3 sub-panels share
    an aspect ratio and render at equal height."""
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar([0, 1], [orig_total, qm_total], width=0.55, color=[_figure_msv000093870_galatidou_vs_quantmsdiann__ORIG_COLOUR, _figure_msv000093870_galatidou_vs_quantmsdiann__QM_COLOUR], edgecolor='#37474f')
    ax.set_xlim(-0.7, 1.7)
    for b, v in zip(bars, (orig_total, qm_total)):
        ax.text(b.get_x() + b.get_width() / 2, v, f'{v:,}', ha='center', va='bottom', fontsize=11)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Galatidou 2024\n(original)', 'quantmsdiann\n(DIA-NN)'])
    ax.set_ylabel('Protein groups (dataset, 1\\% FDR)'.replace('\\%', '%'))
    ax.set_ylim(0, max(orig_total, qm_total) * 1.16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def figure_msv000093870_galatidou_vs_quantmsdiann_main() -> int:
    confident = load_channel_confident(_cached_report())
    qm_cells = per_cell_counts(confident)
    qm_total = confident['Protein.Group'].nunique()
    qm_acc = quantms_accessions(confident)
    qm_lead = {_strip_isoform(str(g).split(';')[0]) for g in confident['Protein.Group'].dropna().unique()}
    qm_lead.discard('')
    matrix_path = _cached_original()
    orig_per_cell, orig_acc = original_per_cell(matrix_path)
    orig_total = len(orig_acc)
    shared = len(orig_acc & qm_lead)
    orig_only = len(orig_acc - qm_lead)
    qm_only = len(qm_lead - orig_acc)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    render_comparison(orig_per_cell, qm_cells, orig_total=len(orig_acc), qm_total=qm_total, shared=shared, orig_only=orig_only, qm_only=qm_only, svg_path=FIGURES_DIR / 'main_galatidou_comparison.svg')
    render_total_pg(orig_total=len(orig_acc), qm_total=qm_total, svg_path=FIGURES_DIR / 'main_galatidou_total_pg.svg')
    matched = qc_matched_table(matrix_path, qm_cells)
    r = render_qc_matched(matched, FIGURES_DIR / 'main_galatidou_qc_matched.svg')
    counts = FIGURES_DIR / 'comparison_counts.tsv'
    with counts.open('w') as fh:
        fh.write('metric\tGalatidou_2024\tquantmsdiann\n')
        fh.write(f'cells\t{len(orig_per_cell)}\t{len(qm_cells)}\n')
        fh.write(f"median_proteins_per_cell\t{int(orig_per_cell.median())}\t{int(qm_cells['proteins'].median())}\n")
        fh.write(f'protein_groups_total\t{len(orig_acc)}\t{qm_total}\n')
        fh.write(f'protein_groups_shared\t{shared}\t{shared}\n')
        fh.write(f'protein_groups_unique\t{orig_only}\t{qm_only}\n')
        fh.write(f'quantmsdiann_expanded_accessions\t-\t{len(qm_acc)}\n')
        fh.write(f'qc_matched_cells\t{len(matched)}\t{len(matched)}\n')
        fh.write(f"qc_matched_median_proteins_per_cell\t{int(matched['orig_proteins'].median())}\t{int(matched['qm_proteins'].median())}\n")
        fh.write(f'qc_matched_pearson_r\t{r:.3f}\t{r:.3f}\n')
    print('=== Galatidou 2024 vs quantmsdiann ===')
    print(f'cells:            {len(orig_per_cell)} vs {len(qm_cells)}')
    print(f"median prot/cell: {int(orig_per_cell.median())} vs {int(qm_cells['proteins'].median())}")
    print(f'protein groups:   {len(orig_acc)} vs {qm_total}')
    print(f'protein groups:   shared {shared}, Galatidou-only {orig_only}, quantmsdiann-only {qm_only} ({100 * shared / len(orig_acc):.0f}% of original recovered)')
    print(f'--- QC-matched cohort ({len(matched)} oocytes) ---')
    print(f"median prot/cell: {int(matched['orig_proteins'].median())} (orig) vs {int(matched['qm_proteins'].median())} (quantmsdiann), Pearson r={r:.2f}")
    print(f"Wrote {FIGURES_DIR / 'main_galatidou_comparison.svg'}")
    print(f"Wrote {FIGURES_DIR / 'main_galatidou_qc_matched.svg'}")
    print(f'Wrote {counts}')
    return 0


# ======================================================================
# inlined from analysis/plexDIA/figure_msv000093870_oocyte_plexdia.py
# ======================================================================

"""MSV000093870 plexDIA reanalysis — first plexDIA cohort through quantmsdiann.

Single-cell oocyte plexDIA dataset (Slavov-lab deposition on MassIVE),
reanalysed with quantmsdiann using DIA-NN 2.5.0 in mTRAQ 3-channel mode
(channels 0 / 4 / 8). Each of the 38 Q-Exactive raw files multiplexes three
single cells, one per mTRAQ channel, for up to 114 single-cell proteomes.

This is the first analysis of the plexDIA branch. It characterises the
per-channel, per-single-cell identification depth — the metric that matters
for plexDIA, since the headline run-level numbers conflate three cells.

Inputs are pulled from the public quantmsdiann-benchmarks deposition and
cached on disk:
  https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/
    quantmsdiann-benchmarks/plexDIA/MSV000093870-plexDIA/

Outputs (paper-ready, no titles/footers):
  analysis/figures/plexDIA/MSV000093870/
    main_plexdia_per_cell.{svg,pdf,png}  — proteins & precursors per single cell, by channel
    counts.tsv                           — auditable per-cell and per-channel totals
"""
import sys
import urllib.request
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
_figure_msv000093870_oocyte_plexdia__REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = _figure_msv000093870_oocyte_plexdia__REPO_ROOT / 'analysis' / 'figures' / 'plexDIA' / 'MSV000093870'
CACHE_DIR = FIGURES_DIR / 'data' / 'cache'
_figure_msv000093870_oocyte_plexdia__FTP_BASE = 'https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/single-cell/MSV000093870/v2_5_1'
REPORT_PARQUET_URL = f'{_figure_msv000093870_oocyte_plexdia__FTP_BASE}/quant_tables/diann_report.parquet'
CHANNEL_COLOURS = {'0': '#90caf9', '4': '#1e88e5', '8': '#0d47a1'}
CHANNEL_LABELS = {'0': 'mTRAQ-0', '4': 'mTRAQ-4', '8': 'mTRAQ-8'}
REPORT_COLUMNS = ['Run', 'Channel', 'Protein.Group', 'Protein.Ids', 'Precursor.Id', 'Decoy', 'Channel.Q.Value', 'PG.Q.Value']

def _cached_report() -> Path:
    """Download the DIA-NN report parquet once and cache it on disk. The cache
    filename is deliberately NOT 'diann_report.parquet' so purge_raw_downloads()
    leaves it: both plexDIA stages (plexdia_per_cell, plexdia_vs_galatidou) share
    this report, and purging it between them would force a wasteful ~245 MB
    re-download."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / 'msv093870_plexdia_report.parquet'
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f'Downloading {REPORT_PARQUET_URL} (~245 MB, cached)…', file=sys.stderr)
    with urllib.request.urlopen(REPORT_PARQUET_URL, timeout=600) as resp:
        dest.write_bytes(resp.read())
    return dest

def load_channel_confident(report_path: Path) -> pd.DataFrame:
    """Decoy-dropped per-channel precursor rows for per-cell counting.

    Filter rule (Vadim review, 2026-06-21): per-cell (= per Run x Channel)
    numbers admit only the per-run q-values, and nothing else. In plexDIA the
    *channel* is the per-run unit (a single cell = one channel of one run), so
    the per-cell precursor q-value is ``Channel.Q.Value`` (the channel-level
    analog of ``Q.Value``; run-level ``Q.Value`` passes a precursor in all
    channels and is not a per-cell number). Per-cell protein groups use
    ``PG.Q.Value``. We drop the previous positive-quantity filter (zero
    quantities are counted) and the contaminant/target filter. Decoys
    (``Decoy == 1``) are dropped. The cut-offs are applied in
    :func:`per_cell_counts`.
    """
    df = pq.read_table(report_path, columns=REPORT_COLUMNS).to_pandas()
    df['Channel'] = df['Channel'].astype(str)
    if 'Decoy' in df.columns:
        df = df[df['Decoy'] == 0]
    return df.copy()

def per_cell_counts(confident: pd.DataFrame) -> pd.DataFrame:
    """One row per single cell = (Run, Channel): precursors (Channel.Q.Value
    <= 1%) and protein groups (PG.Q.Value <= 1%), per the Vadim per-run rule
    (the channel is the per-run unit in plexDIA)."""
    prec = confident[confident['Channel.Q.Value'] <= 0.01].groupby(['Run', 'Channel'])['Precursor.Id'].nunique().rename('precursors')
    prot = confident[confident['PG.Q.Value'] <= 0.01].groupby(['Run', 'Channel'])['Protein.Group'].nunique().rename('proteins')
    cells = pd.concat([prec, prot], axis=1).fillna(0).astype(int).reset_index()
    return cells.sort_values(['Channel', 'Run'])

def render_per_cell_figure(cells: pd.DataFrame, svg_path: Path) -> None:
    """Two-panel boxplot: proteins and precursors per single cell, by channel."""
    channels = ['0', '4', '8']
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4))
    for ax, metric, ylabel in ((axes[0], 'proteins', 'Protein groups per single cell (1% FDR)'), (axes[1], 'precursors', 'Precursors per single cell (1% FDR)')):
        data = [cells.loc[cells['Channel'] == ch, metric].values for ch in channels]
        bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color='#212121', linewidth=1.4))
        for patch, ch in zip(bp['boxes'], channels):
            patch.set_facecolor(CHANNEL_COLOURS[ch])
            patch.set_alpha(0.85)
            patch.set_edgecolor('#37474f')
        rng = np.random.default_rng(0)
        for i, ch in enumerate(channels, start=1):
            y = cells.loc[cells['Channel'] == ch, metric].values
            x = rng.normal(i, 0.05, size=len(y))
            ax.scatter(x, y, s=10, c='#37474f', alpha=0.45, linewidths=0, zorder=3)
        ax.set_xticks(range(1, len(channels) + 1))
        ax.set_xticklabels([CHANNEL_LABELS[ch] for ch in channels])
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_ylim(bottom=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

def write_counts(cells: pd.DataFrame, confident: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as fh:
        fh.write('metric\tchannel\tvalue\tsource\n')
        fh.write(f'single_cells_total\tall\t{len(cells)}\tRun x Channel (channel-confident)\n')
        fh.write(f"runs\tall\t{cells['Run'].nunique()}\tdiann_report.parquet\n")
        fh.write(f"unique_protein_groups\tall\t{confident['Protein.Group'].nunique()}\ttarget, channel-confident\n")
        fh.write(f"unique_precursors\tall\t{confident['Precursor.Id'].nunique()}\ttarget, channel-confident\n")
        fh.write(f"median_proteins_per_cell\tall\t{int(cells['proteins'].median())}\tmedian across cells\n")
        fh.write(f"median_precursors_per_cell\tall\t{int(cells['precursors'].median())}\tmedian across cells\n")
        for ch in ['0', '4', '8']:
            sub = cells[cells['Channel'] == ch]
            fh.write(f'n_cells\t{CHANNEL_LABELS[ch]}\t{len(sub)}\tcells in channel\n')
            fh.write(f"median_proteins_per_cell\t{CHANNEL_LABELS[ch]}\t{int(sub['proteins'].median())}\tmedian\n")
            fh.write(f"median_precursors_per_cell\t{CHANNEL_LABELS[ch]}\t{int(sub['precursors'].median())}\tmedian\n")

def figure_msv000093870_oocyte_plexdia_main() -> int:
    report = _cached_report()
    confident = load_channel_confident(report)
    cells = per_cell_counts(confident)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    render_per_cell_figure(cells, FIGURES_DIR / 'main_plexdia_per_cell.svg')
    write_counts(cells, confident, FIGURES_DIR / 'counts.tsv')
    print(f"Single cells (run x channel): {len(cells)} across {cells['Run'].nunique()} runs")
    print(f"Median proteins/cell: {int(cells['proteins'].median())}, median precursors/cell: {int(cells['precursors'].median())}")
    print(f"Unique target protein groups: {confident['Protein.Group'].nunique()}, precursors: {confident['Precursor.Id'].nunique()}")
    print(f"Wrote {FIGURES_DIR / 'main_plexdia_per_cell.svg'}")
    print(f"Wrote {FIGURES_DIR / 'counts.tsv'}")
    return 0


# ======================================================================
# inlined from analysis/proteobench_metrics.py
# ======================================================================

"""ProteoBench quant-accuracy metrics for our 20 benchmark analyses.

This module wraps ProteoBench's `QuantScoresHYE` + `QuantDatapointHYE`
pipeline so we can compute the same per-replicate-threshold metrics
ProteoBench publishes (median_abs_epsilon_global, ROC AUC, per-species
log2 fold-change, CV_median, ...) on **our** `diann_report.pr_matrix.tsv`
files, without having to either resubmit to ProteoBench or recompute the
metric formulas ourselves.

Two adapters are needed because quantmsdiann's pr_matrix.tsv format
differs from the DIA-NN main report ProteoBench's parsers expect:

1. **Wide → long melt.** pr_matrix.tsv stores one column per run; we
   melt to a long-format DataFrame whose columns match what
   `ParseSettingsBuilder` expects (`Modified.Sequence`, `Protein.Ids`,
   `Precursor.Charge`, `Run`, `Precursor.Normalised`).
2. **Species re-annotation.** quantmsdiann strips species suffixes from
   the Protein.Ids column (bare `Q96P70` instead of `Q96P70;Q96P70_HUMAN`).
   ProteoBench's species detection works by substring match on the
   Protein.Ids string for `_HUMAN` / `_YEAST` / `_ECOLI`. We rebuild a
   UniProt-accession → species map once from three SwissProt FASTA
   streams cached under `data/quantmsdiann_benchmarks/uniprot/` and use
   it to add the species suffix to each accession before passing to
   ProteoBench.

Per-analysis results are cached as a single JSON under
`data/quantmsdiann_benchmarks/proteobench_metrics/<dataset>_<version>.json`
so the figure rebuild is offline thereafter. 20 invocations total when
the cache is cold.
"""
import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Iterable
import pandas as pd
_proteobench_metrics__REPO_ROOT = Path(__file__).resolve().parent.parent
_proteobench_metrics__DATA_DIR = _proteobench_metrics__REPO_ROOT / 'data' / 'quantmsdiann_benchmarks'
UNIPROT_DIR = _proteobench_metrics__DATA_DIR / 'uniprot'
METRICS_CACHE_DIR = _proteobench_metrics__DATA_DIR / 'proteobench_metrics'
DATASET_TO_MODULE = {'PXD049412': 'quant_lfq_DIA_ion_singlecell', 'PXD062685': 'quant_lfq_DIA_ion_diaPASEF', 'PXD070049': 'quant_lfq_DIA_ion_ZenoTOF', 'ProteoBench_Module_7': 'quant_lfq_DIA_ion_Astral'}
MODULE_TO_PARSE_DIR = {'quant_lfq_DIA_ion_singlecell': ('Quant', 'lfq', 'DIA', 'ion', 'singlecell'), 'quant_lfq_DIA_ion_diaPASEF': ('Quant', 'lfq', 'DIA', 'ion', 'diaPASEF'), 'quant_lfq_DIA_ion_ZenoTOF': ('Quant', 'lfq', 'DIA', 'ion', 'ZenoTOF'), 'quant_lfq_DIA_ion_Astral': ('Quant', 'lfq', 'DIA', 'ion', 'Astral')}
PR_MATRIX_META_COLS = {'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Precursor.Id'}
_RUN_EXT_RE = re.compile('\\.(raw|d|mzML)$', re.IGNORECASE)
_SPECIES_FASTAS = {'HUMAN': UNIPROT_DIR / 'human_reviewed.fasta', 'YEAST': UNIPROT_DIR / 'yeast_reviewed.fasta', 'ECOLI': UNIPROT_DIR / 'ecoli_reviewed.fasta'}

def _parse_accessions_from_fasta(fasta_path: Path) -> set[str]:
    """Return the set of UniProt primary accessions from a SwissProt FASTA.
    Each `>sp|ACC|NAME_SPECIES ...` header contributes the ACC token."""
    out: set[str] = set()
    with open(fasta_path, encoding='utf-8') as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            parts = line[1:].split('|', 2)
            if len(parts) >= 2:
                out.add(parts[1].strip())
    return out

@lru_cache(maxsize=1)
def build_species_map() -> dict[str, str]:
    """Build a UniProt accession → species string map from the three
    cached SwissProt FASTAs. Returns a dict keyed by accession. Raises
    FileNotFoundError if a FASTA is missing — those files must be
    fetched once via the README instructions."""
    out: dict[str, str] = {}
    for species, path in _SPECIES_FASTAS.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing SwissProt cache: {path}. Fetch via\n  curl -s -o {path} 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:<id>&format=fasta&compressed=false'")
        for acc in _parse_accessions_from_fasta(path):
            out[acc] = species
    return out

def annotate_species_suffix(protein_ids: str, species_map: dict[str, str]) -> str:
    """Re-annotate a semicolon-separated `Protein.Ids` cell from
    quantmsdiann's bare-accession output (`Q96P70;P09417`) to the
    species-suffixed form ProteoBench expects (`Q96P70_HUMAN;P09417_HUMAN`).

    Accessions not in the species map are passed through unchanged
    (e.g. contaminants, decoys). The exact suffix style matches
    ProteoBench's species_mapper rules — a simple `_HUMAN` / `_YEAST` /
    `_ECOLI` token anywhere in the cell makes `str.contains("_HUMAN")`
    return True."""
    if not protein_ids or protein_ids == 'nan':
        return protein_ids
    out_tokens: list[str] = []
    for acc in str(protein_ids).split(';'):
        acc_clean = acc.strip()
        species = species_map.get(acc_clean)
        if species:
            out_tokens.append(f'{acc_clean}_{species}')
        else:
            out_tokens.append(acc_clean)
    return ';'.join(out_tokens)

def melt_pr_matrix(matrix_path: Path) -> pd.DataFrame:
    """Read a DIA-NN `diann_report.pr_matrix.tsv` and melt it to the
    long-format DataFrame ProteoBench's `convert_to_standard_format`
    consumes.

    Returns a DataFrame with columns
    `[Modified.Sequence, Protein.Ids, Precursor.Charge, Run,
    Precursor.Normalised]`. NaN intensities are dropped (a precursor
    that wasn't quantified in a given run contributes no row). The
    `Run` column has the `.raw` / `.d` / `.mzML` extension stripped to
    match ProteoBench's bare-filename condition_mapper keys.
    """
    wide = pd.read_csv(matrix_path, sep='\t')
    run_cols = [c for c in wide.columns if c not in PR_MATRIX_META_COLS]
    long_df = wide.melt(id_vars=['Modified.Sequence', 'Protein.Ids', 'Precursor.Charge'], value_vars=run_cols, var_name='Run', value_name='Precursor.Normalised').dropna(subset=['Precursor.Normalised'])
    long_df['Run'] = long_df['Run'].str.replace(_RUN_EXT_RE, '', regex=True)
    return long_df.reset_index(drop=True)

def _proteobench_parse_settings_dir(module_id: str) -> Path:
    """Return the absolute path to the per-module parse-settings dir
    bundled with the installed `proteobench` package."""
    import proteobench
    return Path(proteobench.__file__).parent / 'io' / 'parsing' / 'io_parse_settings' / Path(*MODULE_TO_PARSE_DIR[module_id])

def compute_proteobench_metrics(matrix_path: Path, module_id: str, *, software_version: str='', proteobench_version_pin: str | None=None) -> dict:
    """Run ProteoBench's metric computation on one quantmsdiann
    `pr_matrix.tsv`. Returns the per-replicate-threshold `results`
    dict (keys `"1"..."6"`) plus headline aliases used by the figures.

    The function is intentionally side-effect-free — it does NOT write
    a cache file; the caller (`cached_proteobench_metrics`) handles
    caching so unit tests can exercise this path in isolation."""
    from proteobench.io.parsing.parse_settings import ParseSettingsBuilder
    from proteobench.score.quantscoresHYE import QuantScoresHYE
    from proteobench.datapoint.quant_datapoint import QuantDatapointHYE
    species_map = build_species_map()
    long_df = melt_pr_matrix(matrix_path)
    long_df['Protein.Ids'] = long_df['Protein.Ids'].astype(str).apply(lambda s: annotate_species_suffix(s, species_map))
    parse_settings_dir = _proteobench_parse_settings_dir(module_id)
    parser = ParseSettingsBuilder(parse_settings_dir=str(parse_settings_dir), module_id=module_id).build_parser('DIA-NN')
    standard_format, replicate_to_raw = parser.convert_to_standard_format(long_df)
    score = QuantScoresHYE('precursor ion', parser.species_expected_ratio(), parser.species_dict())
    intermediate = score.generate_intermediate(standard_format, replicate_to_raw)
    user_input = {'software_version': software_version, 'search_engine': 'DIA-NN', 'search_engine_version': software_version, 'ident_fdr_psm': 0.01, 'ident_fdr_peptide': None, 'ident_fdr_protein': 0.01, 'enable_match_between_runs': False, 'precursor_mass_tolerance': '', 'fragment_mass_tolerance': '', 'enzyme': 'Trypsin/P', 'allowed_miscleavages': 1, 'min_peptide_length': 7, 'max_peptide_length': 50}
    datapoint = QuantDatapointHYE.generate_datapoint(intermediate, 'DIA-NN', user_input)
    results_str_keys = {str(k): {kk: _jsonable(vv) for kk, vv in v.items()} for k, v in datapoint.results.items()}
    return {'matrix_path': str(matrix_path), 'module_id': module_id, 'software_version': software_version, 'proteobench_version': proteobench_version_pin or getattr(__import__('proteobench'), '__version__', ''), 'n_long_rows': int(len(long_df)), 'n_standard_rows': int(len(standard_format)), 'n_intermediate_rows': int(len(intermediate)), 'replicate_to_raw': {k: list(v) for k, v in replicate_to_raw.items()}, 'results': results_str_keys}

def _jsonable(value):
    """Coerce a metric value to a JSON-serialisable scalar. NumPy
    scalars get cast to Python types; lists pass through; everything
    else is stringified as a fallback."""
    if value is None:
        return None
    if hasattr(value, 'item'):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, (int, float, str, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return str(value)

def metrics_cache_path(dataset: str, version: str) -> Path:
    return METRICS_CACHE_DIR / f'{dataset}_{version}.json'

def cached_proteobench_metrics(dataset: str, version: str, *, fetch: bool=True) -> dict:
    """Cache-or-compute wrapper. Returns the metrics dict either from
    the on-disk JSON cache or by computing it fresh.

    `fetch=False` skips computation entirely — the cache file must
    already exist (raises FileNotFoundError otherwise). Used by the
    figure rebuild path so it never silently invokes ProteoBench's
    parser on a CI runner that lacks the SwissProt FASTAs."""
    cache_path = metrics_cache_path(dataset, version)
    if cache_path.exists():
        with open(cache_path, encoding='utf-8') as fh:
            return json.load(fh)
    if not fetch:
        raise FileNotFoundError(f'Cache miss for {dataset}/{version} at {cache_path} and fetch=False. Run with fetch=True (default) to populate.')
    module_id = DATASET_TO_MODULE[dataset]
    matrix = _proteobench_metrics__DATA_DIR / dataset / version / 'diann_report.pr_matrix.tsv'
    payload = compute_proteobench_metrics(matrix, module_id, software_version=version)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2)
    return payload

def iter_metric_rows(payloads: Iterable[tuple[str, str, dict]]) -> list[dict]:
    """Yield long-format rows from a sequence of (dataset, version,
    payload) tuples. Each row is one (dataset, version, threshold,
    metric) combination. Headline metrics + per-species variants are
    emitted; the figure consumer selects what it needs."""
    rows: list[dict] = []
    headline_metrics = ['nr_prec', 'nr_prec_HUMAN', 'nr_prec_YEAST', 'nr_prec_ECOLI', 'median_abs_epsilon_global', 'mean_abs_epsilon_global', 'median_abs_epsilon_eq_species', 'mean_abs_epsilon_eq_species', 'median_abs_epsilon_HUMAN', 'median_abs_epsilon_YEAST', 'median_abs_epsilon_ECOLI', 'mean_log2_empirical_HUMAN', 'mean_log2_empirical_YEAST', 'mean_log2_empirical_ECOLI', 'median_log2_empirical_HUMAN', 'median_log2_empirical_YEAST', 'median_log2_empirical_ECOLI', 'CV_median', 'CV_q75', 'CV_q90', 'CV_q95', 'roc_auc', 'roc_auc_directional', 'variance_epsilon_global']
    for dataset, version, payload in payloads:
        results = payload.get('results', {})
        for threshold_str, metric_dict in results.items():
            try:
                threshold = int(threshold_str)
            except ValueError:
                continue
            for metric in headline_metrics:
                value = metric_dict.get(metric)
                if value is None:
                    continue
                rows.append({'dataset': dataset, 'version': version, 'threshold': threshold, 'metric': metric, 'value': value})
    return rows


# ======================================================================
# inlined from analysis/venn_protein_accessions.py
# ======================================================================

import re
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fs.apply_house_style()
_venn_protein_accessions__REPO_ROOT = Path(__file__).resolve().parent.parent
_venn_protein_accessions__DATA_DIR = _venn_protein_accessions__REPO_ROOT / 'data' / 'PXD003539'
_venn_protein_accessions__FIGURES_DIR = _venn_protein_accessions__REPO_ROOT / 'analysis' / 'figures' / 'PXD003539'
_ACCESSION_RE = re.compile('^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:[A-Z0-9]{1,5})?$')

def _clean_token(tok: str) -> str | None:
    """Normalise an accession token. The caller has already verified
    via `is_target_protein_group` that the surrounding Protein.Group
    has no contaminant/entrapment/decoy prefix, so `strip_known_prefix`
    here is defensive — covers any future asymmetry between the row-level
    filter and per-token normalisation."""
    tok = tok.strip()
    if not tok:
        return None
    if '|' in tok:
        parts = tok.split('|')
        if len(parts) >= 2 and parts[1]:
            tok = parts[1]
        else:
            return None
    tok = strip_known_prefix(tok)
    if '-' in tok:
        tok = tok.split('-', 1)[0]
    if not tok:
        return None
    return tok

def extract_accessions_diann(protein_group: str | None) -> set[str]:
    """Return the set of target UniProt accessions in `protein_group`.

    Conservative policy: returns the empty set if **any**
    semicolon-separated accession token carries a contaminant /
    entrapment / decoy prefix. Mixed groups like
    `CONTAM_P02768;P02768` are dropped entirely — the conservative
    interpretation per the 2026-05-21 contaminant-filter spec.

    Pure-target rows pass through unchanged. Empty / None inputs
    return the empty set."""
    if not protein_group:
        return set()
    if not is_target_protein_group(protein_group):
        return set()
    out: set[str] = set()
    for piece in protein_group.split(';'):
        cleaned = _clean_token(piece)
        if cleaned:
            out.add(cleaned)
    return out

def extract_accessions_openswath(protein_str: str | None) -> set[str]:
    """OpenSWATH protein-string extractor. Same conservative filter
    as `extract_accessions_diann`: drop the row if any token after
    the leading count carries a prefix. The leading count field
    (`1`, `2`, `DECOY_1`, ...) is *not* itself an accession and is
    always skipped — but if it carries a `DECOY_` prefix we treat the
    whole row as a decoy and drop it.
    """
    if not protein_str:
        return set()
    parts = protein_str.split('/')
    if parts and parts[0].startswith(('DECOY_', 'decoy_')):
        return set()
    accession_tokens = [p for p in parts[1:] if p.strip()]
    if not is_target_protein_group(';'.join(accession_tokens)):
        return set()
    out: set[str] = set()
    for piece in accession_tokens:
        cleaned = _clean_token(piece)
        if cleaned:
            out.add(cleaned)
    return out

def accessions_with_min_peptides_diann(pr_matrix_path: Path, *, min_peptides: int=2) -> set[str]:
    counts = unique_peptides_per_protein_diann(pr_matrix_path)
    out: set[str] = set()
    for pg, n in counts.items():
        if n >= min_peptides:
            out.update(extract_accessions_diann(pg))
    return out

def accessions_with_min_peptides_openswath(matrix_path: Path, *, min_peptides: int=2) -> set[str]:
    counts = unique_peptides_per_protein_openswath(matrix_path)
    out: set[str] = set()
    for prot, n in counts.items():
        if n >= min_peptides:
            out.update(extract_accessions_openswath(prot))
    return out

def render_venn_diagram(guo_acc: set[str], diann_acc: set[str], svg_path: Path, *, left_label: str='Guo 2019 (OpenSWATH)', right_label: str='quantmsdiann (DIA-NN)', left_color: str='#9e9e9e', right_color: str='#1f77b4', title: str | None=None, footer: str | None=None) -> None:
    full_left = f'{left_label}\n(n={len(guo_acc):,})'
    full_right = f'{right_label}\n(n={len(diann_acc):,})'
    inter = guo_acc & diann_acc
    guo_only = guo_acc - diann_acc
    diann_only = diann_acc - guo_acc
    total = len(guo_acc | diann_acc) or 1
    fig, ax = plt.subplots(figsize=(7, 5.5))
    try:
        from matplotlib_venn import venn2
        v = venn2(subsets=(len(guo_only), len(diann_only), len(inter)), set_labels=(full_left, full_right), set_colors=(left_color, right_color), alpha=0.55, ax=ax)
        labels = {'10': (len(guo_only), guo_only), '01': (len(diann_only), diann_only), '11': (len(inter), inter)}
        for region_id, (count, _) in labels.items():
            lbl = v.get_label_by_id(region_id)
            if lbl is None:
                continue
            pct = 100.0 * count / total
            lbl.set_text(f'{count:,}\n({pct:.1f}%)')
            lbl.set_fontsize(10)
        for sl in v.set_labels:
            if sl is not None:
                sl.set_fontsize(11)
    except ImportError:
        from matplotlib.patches import Circle
        ax.set_xlim(-2, 4)
        ax.set_ylim(-2, 2)
        ax.set_aspect('equal')
        ax.axis('off')
        c_left = Circle((0, 0), 1.3, color=left_color, alpha=0.55, linewidth=0)
        c_right = Circle((1.4, 0), 1.3, color=right_color, alpha=0.55, linewidth=0)
        ax.add_patch(c_left)
        ax.add_patch(c_right)
        ax.text(-0.9, 0, f'{len(guo_only):,}\n({100 * len(guo_only) / total:.1f}%)', ha='center', va='center', fontsize=10)
        ax.text(2.3, 0, f'{len(diann_only):,}\n({100 * len(diann_only) / total:.1f}%)', ha='center', va='center', fontsize=10)
        ax.text(0.7, 0, f'{len(inter):,}\n({100 * len(inter) / total:.1f}%)', ha='center', va='center', fontsize=10)
        ax.text(-0.6, 1.5, full_left, ha='center', va='bottom', fontsize=11)
        ax.text(2.0, 1.5, full_right, ha='center', va='bottom', fontsize=11)
    fig.tight_layout()
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg_path)
    plt.close(fig)

def venn_protein_accessions_main() -> int:
    ensure_cell_line_matrices('PXD003539', with_report=True)  # 2-set venn (Guo vs quantmsdiann) is PXD003539 only
    pr_path = _venn_protein_accessions__DATA_DIR / 'diann_report.pr_matrix.tsv'
    opensw_path = _venn_protein_accessions__DATA_DIR / 'feature_alignment_requant_matrix.tsv'
    if not pr_path.exists():
        print(f'Missing input: {pr_path}', file=sys.stderr)
        return 1
    if not opensw_path.exists():
        print(f'Missing input: {opensw_path}', file=sys.stderr)
        return 1
    print('Computing DIA-NN accession set (>=2 unique peptides)...')
    diann_acc = accessions_with_min_peptides_diann(pr_path, min_peptides=2)
    print('Computing Guo/OpenSWATH accession set (>=2 unique peptides)...')
    guo_acc = accessions_with_min_peptides_openswath(opensw_path, min_peptides=2)
    inter = guo_acc & diann_acc
    guo_only = guo_acc - diann_acc
    diann_only = diann_acc - guo_acc
    print(f'Guo total:        {len(guo_acc):,}')
    print(f'DIA-NN total:     {len(diann_acc):,}')
    print(f'Intersection:     {len(inter):,}')
    print(f'Guo only:         {len(guo_only):,}')
    print(f'DIA-NN only:      {len(diann_only):,}')
    svg_path = _venn_protein_accessions__FIGURES_DIR / 'supp_venn_protein_accessions.svg'
    render_venn_diagram(guo_acc, diann_acc, svg_path)
    print(f'Figure saved to {svg_path}')
    return 0


# ======================================================================
# Orchestrator (CLI, STAGES registry).
# ======================================================================

_BENCH_FTP = ("https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/"
              "proteomes/quantmsdiann-benchmarks/proteobench")
_CELLLINES_FTP = ("https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/"
                  "proteomes/quantmsdiann-benchmarks/cell-lines")


def recompute_report_counts() -> int:
    """Download each ProteoBench cohort's DIA-NN report from the public PRIDE
    FTP and recount under the methods.md §1 rule -> report_counts.tsv. The
    counting primitive is the inlined ``count_report`` (was
    analysis.count_report_ids.count_report)."""
    root = REPO / "data" / "quantmsdiann_benchmarks"
    cols = ["dataset", "version", "prec_min1", "prec_min3", "prec_global",
            "prot_global", "prot_perrun_avg", "prot_complete", "peptides", "prot_2pep"]
    rows = []
    for dataset in DATASET_MODULES:
        for version in _count_report_ids__VERSIONS:
            qt = root / dataset / version / "quant_tables"
            parquet, tsv = qt / "diann_report.parquet", qt / "diann_report.tsv"
            if not ((parquet.exists() and parquet.stat().st_size)
                    or (tsv.exists() and tsv.stat().st_size)):
                fn = "diann_report.tsv" if version == "v1_8_1" else "diann_report.parquet"
                download_if_missing(f"{_BENCH_FTP}/{dataset}/{version}/quant_tables/{fn}", qt / fn)
            c = count_report(_load_report(qt),
                             precursor_q=PRECURSOR_Q.get(version, DEFAULT_PRECURSOR_Q))
            c.update(dataset=dataset, version=version)
            rows.append(c)
            print(f"  {dataset} {version}: prec_global={c['prec_global']:,} "
                  f"prot_global={c['prot_global']:,}", flush=True)
            discard_download(parquet, tsv)  # bound disk: drop the report once counted
    pd.DataFrame(rows)[cols].to_csv(root / "report_counts.tsv", sep="\t", index=False)
    print(f"  wrote {root / 'report_counts.tsv'}", flush=True)
    return 0


def recompute_reanalysis_pg_counts() -> int:
    """Stream each bulk cell-line cohort's report from the FTP and write the
    per-cohort diann_report_protein_counts.json (prot_global at Lib.PG.Q.Value,
    prot_2pep at Lib.Q.Value) under methods.md §1. Decoys dropped, nothing else
    filtered. Consumed by the PXD004701 / PXD030304 figure code."""
    import fsspec
    reports = {
        "PXD004701": f"{_CELLLINES_FTP}/PXD004701/v2_5_1/quant_tables/diann_report.parquet",
        "PXD030304": f"{_CELLLINES_FTP}/PXD030304/v2_5_1/quant_tables/diann_report.parquet",
    }
    cols = ["Protein.Group", "Stripped.Sequence", "Lib.Q.Value", "Lib.PG.Q.Value", "Decoy"]
    for acc_id, url in reports.items():
        lib_pg: set = set()
        peps_pg: dict = collections.defaultdict(set)
        with fsspec.filesystem("https").open(url, "rb") as fh:
            pf = pq.ParquetFile(fh)
            use = [c for c in cols if c in set(pf.schema_arrow.names)]
            for batch in pf.iter_batches(batch_size=2_000_000, columns=use):
                d = batch.to_pydict()
                n = len(d["Protein.Group"])
                dec = d.get("Decoy", [0] * n)
                for i in range(n):
                    if dec[i] == 1:
                        continue
                    lpg, pg = d["Lib.PG.Q.Value"][i], d["Protein.Group"][i]
                    if lpg is not None and lpg <= 0.01:
                        lib_pg.add(pg)
                        lq = d["Lib.Q.Value"][i]
                        if lq is not None and lq <= 0.01:
                            peps_pg[pg].add(d["Stripped.Sequence"][i])
        counts = {"prot_global": len(lib_pg),
                  "prot_2pep": sum(1 for s in peps_pg.values() if len(s) >= 2)}
        out = REPO / "data" / acc_id / "diann_report_protein_counts.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(counts))
        print(f"  {acc_id}: prot_global={counts['prot_global']:,} "
              f"prot_2pep={counts['prot_2pep']:,} -> {out}", flush=True)
    return 0


# Stage = (name, callable, description). Every target is an in-process callable
# returning an int return code. Order matters: data prep before figures.
DATA_PREP = [
    ("report_counts", recompute_report_counts,
     "Download FTP reports, recount under the filter rule -> report_counts.tsv"),
    ("reanalysis_pg_counts", recompute_reanalysis_pg_counts,
     "Recount bulk-cohort report protein groups (Lib rule) -> per-cohort JSON caches"),
    ("single_cell_tables", make_single_cell_tables_main,
     "Per-cell / completeness / CV tables for the single-cell figure"),
    ("phospho_tables", make_phospho_tables_main,
     "Phosphopeptide / phosphosite tables (Lib.Q.Value; site Prob>=0.99)"),
]

FIGURES = [
    ("benchmarks", figure_quantmsdiann_benchmarks_vs_proteobench_main,
     "Fig 2 — ProteoBench benchmark panels + counts.tsv"),
    ("queue_sweep", figure_queue_size_sweep_main,
     "queue_size_sweep.tsv (consumed by the Fig 2 validation composite)"),
    ("fig2_validation", figure_fig2_validation_main,
     "Fig 2 validation composite"),
    ("id_vs_epsilon", figure_id_vs_epsilon_main,
     "ProteoBench id-vs-epsilon panel"),
    ("proteobench_accuracy", figure_proteobench_accuracy_main,
     "ProteoBench accuracy panels"),
    ("reanalysis_improvement", figure_reanalysis_improvement_main,
     "Reanalysis-recovery figure (original vs quantmsdiann)"),
    ("single_cell_combined", figure_single_cell_combined_main,
     "Single-cell figure (per-cell depth, completeness, CV)"),
    ("plexdia_per_cell", figure_msv000093870_oocyte_plexdia_main,
     "plexDIA per-cell depth (MSV000093870)"),
    ("plexdia_vs_galatidou", figure_msv000093870_galatidou_vs_quantmsdiann_main,
     "plexDIA deposited vs quantmsdiann"),
    ("pxd003539", figure_original_vs_quantmsdiann_main,
     "PXD003539 (NCI-60) panels"),
    ("pxd004701", figure_pxd004701_sun_vs_quantmsdiann_main,
     "PXD004701 (Sun) panels"),
    ("pxd030304", figure_pxd030304_procan_vs_quantmsdiann_main,
     "PXD030304 (ProCan) panels — streams a 2GB matrix"),
    ("pxd064049_spatial", figure_pxd064049_spatial_vs_quantmsdiann_main,
     "PXD064049 spatial DVP panels"),
    ("atlas", figure_combined_cell_lines_atlas_main,
     "Pan-cohort of DIA reanalyses (Fig S13) — needs numpy<2"),
    ("phospho", figure_phospho_main,
     "Phosphoproteomics supplementary figure"),
    ("venn", venn_protein_accessions_main,
     "Protein-accession overlap (supplementary)"),
    ("performance_trace", figure_performance_trace_main,
     "Per-step runtime + resources (runtime_per_step.svg, resources_per_step.svg)"),
    ("mdc_cluster_runtime", figure_mdc_cluster_runtime_main,
     "MDC cluster runtime"),
]

# Runs LAST: aggregates every cited number from the figure-data TSVs into one
# file (data/paper_numbers.tsv + paper/generated_numbers.tex).
COLLECT = [
    ("paper_numbers", collect_paper_numbers_main,
     "Aggregate ALL manuscript numbers -> data/paper_numbers.tsv + paper/generated_numbers.tex"),
]

ALL_STAGES = DATA_PREP + FIGURES + COLLECT
BY_NAME = {name: (name, fn, desc) for name, fn, desc in ALL_STAGES}


def print_list() -> None:
    print("DATA PREP:")
    for name, fn, desc in DATA_PREP:
        print(f"  {name:22s} {desc}")
    print("\nFIGURES:")
    for name, fn, desc in FIGURES:
        print(f"  {name:22s} {desc}")
    print("\nNUMBERS:")
    for name, fn, desc in COLLECT:
        print(f"  {name:22s} {desc}")
    print("\nPDFs: paper/Makefile -> figures, pdf, supplementary")


def run_stage(name: str, target) -> tuple[str, bool, float]:
    """Run one stage in-process. `target` is a callable returning an int rc."""
    t0 = time.time()
    print(f"\n=== [{name}] (in-process) ===", flush=True)
    try:
        rc = int(target() or 0)
    except Exception as exc:  # keep the run going; report as a failed stage
        import traceback
        traceback.print_exc()
        print(f"  ERROR: {exc}", flush=True)
        rc = 1
    dt = time.time() - t0
    ok = rc == 0
    print(f"=== [{name}] {'OK' if ok else 'FAILED (rc=%d)' % rc} in {dt:.0f}s ===",
          flush=True)
    return name, ok, dt


def run_make(target: str) -> tuple[str, bool, float]:
    print(f"\n=== [make {target}] ===", flush=True)
    t0 = time.time()
    rc = subprocess.run(["make", target], cwd=REPO / "paper").returncode
    dt = time.time() - t0
    return f"make {target}", rc == 0, dt


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="data prep + figures + PDFs")
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--data-only", action="store_true")
    ap.add_argument("--only", nargs="+", metavar="NAME", help="run specific stage(s)")
    ap.add_argument("--no-pdf", action="store_true", help="skip the PDF build")
    ap.add_argument("--list", action="store_true", help="list stages and exit")
    ap.add_argument("--fail-fast", action="store_true", help="stop at first failure")
    ap.add_argument("--keep-going", action="store_true",
                    help="continue past failures (default)")
    args = ap.parse_args(argv)

    if args.list:
        print_list()
        return 0

    if args.only:
        unknown = [n for n in args.only if n not in BY_NAME]
        if unknown:
            print(f"unknown stage(s): {', '.join(unknown)}", file=sys.stderr)
            print("use --list to see valid names", file=sys.stderr)
            return 2
        stages = [BY_NAME[n] for n in args.only]
        build_pdf = False
    elif args.data_only:
        stages, build_pdf = DATA_PREP, False
    elif args.figures_only:
        stages, build_pdf = FIGURES + COLLECT, not args.no_pdf
    elif args.all:
        stages, build_pdf = ALL_STAGES, not args.no_pdf
    else:
        ap.print_help()
        return 0

    fail_fast = args.fail_fast and not args.keep_going
    results: list[tuple[str, bool, float]] = []
    for name, fn, _ in stages:
        r = run_stage(name, fn)
        results.append(r)
        purge_raw_downloads()  # drop raw FTP downloads so disk stays bounded across stages
        if not r[1] and fail_fast:
            break

    if build_pdf and (not fail_fast or all(ok for _, ok, _ in results)):
        for target in ("figures", "pdf", "supplementary"):
            results.append(run_make(target))

    print("\n================ SUMMARY ================")
    failed = [n for n, ok, _ in results if not ok]
    for name, ok, dt in results:
        print(f"  {'OK  ' if ok else 'FAIL'}  {name:30s} {dt:6.0f}s")
    if failed:
        print(f"\n{len(failed)} stage(s) failed: {', '.join(failed)}")
        return 1
    print("\nAll stages succeeded.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
