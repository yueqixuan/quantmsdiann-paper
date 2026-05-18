# PXD003539 reanalysis figure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a paper-ready figure comparing precursor + protein counts from Guo et al. 2019 (original SWATH analysis of NCI-60, PXD003539) vs the quantmsdiann reanalysis, demonstrating that quantmsdiann recovers substantially more identifications.

**Architecture:** One self-contained Python script under `analysis/` downloads the DIA-NN report matrices from the PRIDE FTP (cached locally under `data/PXD003539/`), counts entries quantified in ≥ 1 of the 120 runs, and renders a grouped bar chart (PDF + PNG) plus a TSV with the underlying numbers. Original-paper baselines live as hardcoded constants with inline citations to Guo 2019.

**Tech Stack:** Python 3, `pandas`, `matplotlib`, `requests`, `pytest`.

**Spec:** `docs/superpowers/specs/2026-05-18-pxd003539-reanalysis-figure-design.md`

---

## File Structure

Files this plan creates or modifies:

- `analysis/figure_original_vs_quantmsdiann.py` — the one analysis script (downloader, counter, plotter, main).
- `analysis/requirements.txt` — Python dependencies pinned to permissive minimums.
- `analysis/tests/test_count_quantified_rows.py` — unit test for the counting function (the only logic worth covering).
- `analysis/figures/.gitkeep` — keep the empty output directory tracked.
- `data/.gitignore` — ignore downloaded PRIDE artefacts.
- `.gitignore` — extend to ignore Python build/cache artefacts and `data/`.
- `README.md` — minimal note pointing to `analysis/figure_original_vs_quantmsdiann.py` and how to run it.

The script is intentionally a single file. There are no shared modules to extract until we have a second dataset that proves the abstraction.

---

### Task 1: Scaffold the analysis directory

**Files:**
- Create: `analysis/requirements.txt`
- Create: `analysis/figures/.gitkeep`
- Create: `data/.gitignore`
- Modify: `.gitignore`
- Create: `README.md`

- [ ] **Step 1: Initialize git if needed**

Run from repo root `/Users/yperez/work/articles/quantmsdiann`:
```bash
git init 2>/dev/null || true
git status
```
Expected: clean status (or untracked `.github/`, `.gitignore`, `analysis/`).

- [ ] **Step 2: Write `analysis/requirements.txt`**

```
pandas>=2.0
matplotlib>=3.7
requests>=2.31
pytest>=8.0
```

- [ ] **Step 3: Create the figures output directory**

Create empty file `analysis/figures/.gitkeep` (zero bytes). This keeps the directory tracked before any PDFs are generated.

- [ ] **Step 4: Write `data/.gitignore`**

```
*
!.gitignore
```
This commits the cache directory marker but ignores its contents (the downloaded TSVs).

- [ ] **Step 5: Extend the root `.gitignore`**

Append these lines to `/Users/yperez/work/articles/quantmsdiann/.gitignore`:
```
# Python
__pycache__/
*.pyc
.pytest_cache/
.venv/
venv/

# Analysis cache
data/
```

- [ ] **Step 6: Write a minimal `README.md`**

```markdown
# quantmsdiann reanalysis manuscript

Source repository for the quantmsdiann reanalysis paper.

## Analyses

| Dataset    | Script                                              | Description                                                              |
| ---------- | --------------------------------------------------- | ------------------------------------------------------------------------ |
| PXD003539  | `analysis/figure_original_vs_quantmsdiann.py`       | NCI-60 PCT-SWATH (Guo 2019) — original vs quantmsdiann ID counts plot.   |

## Running an analysis

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r analysis/requirements.txt
python analysis/figure_original_vs_quantmsdiann.py
```

Outputs land in `analysis/figures/`. Downloaded inputs are cached under
`data/<PXD-id>/` and are git-ignored.
```

- [ ] **Step 7: Commit**

```bash
git add analysis/requirements.txt analysis/figures/.gitkeep data/.gitignore .gitignore README.md
git commit -m "chore: scaffold analysis directory for quantmsdiann reanalysis paper"
```

---

### Task 2: Confirm Guo 2019 baseline numbers

The script needs two constants: `ORIGINAL_PRECURSORS` and `ORIGINAL_PROTEINS`. These come from Guo et al. 2019, iScience (NCI-60 PCT-SWATH paper). Confirm them before writing the script so we don't ship a figure with wrong baselines.

**Files:**
- Modify: `docs/superpowers/specs/2026-05-18-pxd003539-reanalysis-figure-design.md` (record the confirmed values + source location)

- [ ] **Step 1: Locate the Guo 2019 paper text**

The article PubMed ID is 31733513 and the DOI is `10.1016/j.isci.2019.10.059`. Try fetching the full text via PubMed Central / Europe PMC. The PubMed MCP tool `get_full_text_article` accepts a PMID. If that fails, the open-access XML is at:
```
https://europepmc.org/article/MED/31733513
```

- [ ] **Step 2: Extract the two numbers**

Look in the Results and Methods sections for sentences of the form "we quantified N proteins / M precursors across the NCI-60 panel". The supplementary materials index is at the bottom of the iScience article page; the protein-quantification matrix is typically `mmc2.xlsx` or similar.

Concretely:
- **Precursors:** total unique peptide precursors quantified across all 120 SWATH runs in the published analysis. If the paper only reports a per-cell-line number, sum the union (preferred) or fall back to the largest single-run total reported.
- **Proteins:** total unique SwissProt protein groups quantified across all 120 runs in the published analysis.

Record both numbers, the exact sentence/table they came from, and the file path / page number of the source.

- [ ] **Step 3: Update the spec with the confirmed numbers**

In `docs/superpowers/specs/2026-05-18-pxd003539-reanalysis-figure-design.md`, replace the "Open questions" entry about confirming the numbers with a "Baseline values" subsection under "Inputs → Original analysis", e.g.:

```markdown
**Baseline values (Guo et al. 2019):**

- `ORIGINAL_PRECURSORS = <N>` — from <exact source, e.g. "Table S2 column `peptides_quantified`, unique values across cell lines" or "Results, paragraph 2: 'a total of N peptide precursors…'">.
- `ORIGINAL_PROTEINS = <M>` — from <exact source>.
```

- [ ] **Step 4: Commit the spec update**

```bash
git add docs/superpowers/specs/2026-05-18-pxd003539-reanalysis-figure-design.md
git commit -m "docs(spec): confirm Guo 2019 baseline counts for PXD003539 figure"
```

---

### Task 3: Write the failing test for `count_quantified_rows`

This is the only function in the script with non-trivial logic; everything else is glue. TDD it.

**Files:**
- Create: `analysis/tests/__init__.py` (empty)
- Create: `analysis/tests/test_count_quantified_rows.py`

- [ ] **Step 1: Create the test file**

```python
# analysis/tests/test_count_quantified_rows.py
from pathlib import Path
import textwrap
import pytest

from analysis.figure_original_vs_quantmsdiann import count_quantified_rows


def write_matrix(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "matrix.tsv"
    p.write_text(textwrap.dedent(body).lstrip("\n"))
    return p


def test_counts_rows_with_at_least_one_non_na(tmp_path: Path) -> None:
    # Minimal DIA-NN-style pr_matrix layout: metadata columns end at
    # First.Protein.Description; everything after is per-run quantification.
    matrix = write_matrix(
        tmp_path,
        """
        Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tRun_A\tRun_B\tRun_C
        P1\tP1\tA\tGENE1\tdesc1\t10\t\t20
        P2\tP2\tB\tGENE2\tdesc2\t\t\t
        P3\tP3\tC\tGENE3\tdesc3\t\t5\t
        """,
    )

    n = count_quantified_rows(matrix)

    # P1 and P3 each have one non-empty quant value; P2 has none.
    assert n == 2


def test_raises_if_schema_metadata_anchor_missing(tmp_path: Path) -> None:
    matrix = write_matrix(
        tmp_path,
        """
        Protein.Group\tWeird.Column\tRun_A\tRun_B
        P1\tx\t10\t\t
        """,
    )

    with pytest.raises(ValueError, match="First.Protein.Description"):
        count_quantified_rows(matrix)
```

- [ ] **Step 2: Run the tests and watch them fail**

```bash
cd /Users/yperez/work/articles/quantmsdiann
pytest analysis/tests -v
```
Expected: import error — `analysis.figure_original_vs_quantmsdiann` does not exist yet. That's the failing state we want before Task 4.

---

### Task 4: Implement `count_quantified_rows` and supporting helpers

**Files:**
- Create: `analysis/__init__.py` (empty, to make it a package for the test import)
- Create: `analysis/figure_original_vs_quantmsdiann.py`

- [ ] **Step 1: Create `analysis/__init__.py`**

Empty file (zero bytes).

- [ ] **Step 2: Create the script skeleton with constants, downloader, and counter**

```python
# analysis/figure_original_vs_quantmsdiann.py
"""Compare Guo et al. 2019 (PXD003539, SWATH/OpenSWATH) vs the quantmsdiann
reanalysis on precursor and protein-group counts. Renders a grouped bar chart
to `analysis/figures/`."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "PXD003539"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures"

PRIDE_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantms-collections/absolute-expression-2.0/cell-lines/PXD003539/quant_tables"
)
PR_MATRIX_URL = f"{PRIDE_BASE}/diann_report.pr_matrix.tsv"
PG_MATRIX_URL = f"{PRIDE_BASE}/diann_report.pg_matrix.tsv"
SUMMARY_LOG_URL = f"{PRIDE_BASE}/diannsummary.log"

# Baseline values from Guo et al. 2019 (iScience, PMID 31733513).
# These MUST be filled in from the paper before the script can render.
# See docs/superpowers/specs/2026-05-18-pxd003539-reanalysis-figure-design.md
# for the exact source citation per number.
ORIGINAL_PRECURSORS: int | None = None  # TODO(task-2): set from Guo 2019.
ORIGINAL_PROTEINS: int | None = None    # TODO(task-2): set from Guo 2019.

METADATA_ANCHOR = "First.Protein.Description"


# ---------------------------------------------------------------------------
# Downloader
# ---------------------------------------------------------------------------

def download_if_missing(url: str, dest: Path, *, retries: int = 2) -> Path:
    """Fetch `url` to `dest` exactly once. Idempotent; retries `retries` times
    on network failure before re-raising."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest

    last_err: Exception | None = None
    for attempt in range(retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        if chunk:
                            f.write(chunk)
                tmp.replace(dest)
            return dest
        except requests.RequestException as e:
            last_err = e
            if attempt == retries:
                break
    raise RuntimeError(f"failed to download {url}: {last_err}") from last_err


# ---------------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------------

def _sample_columns(header: Iterable[str]) -> list[str]:
    cols = list(header)
    if METADATA_ANCHOR not in cols:
        raise ValueError(
            f"expected metadata anchor column {METADATA_ANCHOR!r} in matrix "
            f"header; got columns: {cols!r}"
        )
    return cols[cols.index(METADATA_ANCHOR) + 1 :]


def count_quantified_rows(matrix_path: Path) -> int:
    """Count rows with at least one non-NA quantification across the per-run
    sample columns. Sample columns are everything after the DIA-NN metadata
    anchor (`First.Protein.Description`)."""
    df = pd.read_csv(matrix_path, sep="\t", dtype=str, keep_default_na=False)
    samples = _sample_columns(df.columns)
    # Treat empty strings and the literal "NA" as missing.
    quant = df[samples].replace({"": pd.NA, "NA": pd.NA})
    has_any = quant.notna().any(axis=1)
    return int(has_any.sum())


# ---------------------------------------------------------------------------
# main (filled in Task 5)
# ---------------------------------------------------------------------------

def main() -> int:  # pragma: no cover - exercised end-to-end in Task 5
    raise NotImplementedError("main is implemented in Task 5")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
```

- [ ] **Step 3: Run the tests and verify they pass**

```bash
cd /Users/yperez/work/articles/quantmsdiann
python -m pip install -r analysis/requirements.txt
pytest analysis/tests -v
```
Expected: both tests pass.

- [ ] **Step 4: Commit**

```bash
git add analysis/__init__.py analysis/figure_original_vs_quantmsdiann.py analysis/tests/__init__.py analysis/tests/test_count_quantified_rows.py
git commit -m "feat(analysis): add count_quantified_rows with TDD coverage"
```

---

### Task 5: Implement `main` — download, count, render, write outputs

**Files:**
- Modify: `analysis/figure_original_vs_quantmsdiann.py`

- [ ] **Step 1: Replace the placeholder `main` and add the renderer**

Replace the `main` placeholder block (everything from the `# main (filled in Task 5)` comment down to the end of file, exclusive of nothing — replace through end-of-file) with:

```python
# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Counts:
    original_precursors: int
    original_proteins: int
    quantmsdiann_precursors: int
    quantmsdiann_proteins: int


def render_figure(counts: Counts, out_pdf: Path, out_png: Path) -> None:
    import matplotlib.pyplot as plt  # local import keeps test runs lighter

    categories = ["Precursors", "Protein groups"]
    original = [counts.original_precursors, counts.original_proteins]
    reanalysis = [counts.quantmsdiann_precursors, counts.quantmsdiann_proteins]

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    x = range(len(categories))
    width = 0.38
    bars_a = ax.bar([i - width / 2 for i in x], original, width,
                    label="Original (Guo 2019)", color="#9e9e9e",
                    edgecolor="black", linewidth=0.4)
    bars_b = ax.bar([i + width / 2 for i in x], reanalysis, width,
                    label="quantmsdiann", color="#1f77b4",
                    edgecolor="black", linewidth=0.4)

    # Switch to log y if any quantmsdiann/original ratio exceeds 5x.
    ratios = [
        counts.quantmsdiann_precursors / max(counts.original_precursors, 1),
        counts.quantmsdiann_proteins / max(counts.original_proteins, 1),
    ]
    if max(ratios) > 5:
        ax.set_yscale("log")

    ax.set_xticks(list(x))
    ax.set_xticklabels(categories)
    ax.set_ylabel("Count (quantified in ≥ 1 run)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left")

    for bars in (bars_a, bars_b):
        for b in bars:
            h = b.get_height()
            ax.annotate(f"{int(h):,}",
                        xy=(b.get_x() + b.get_width() / 2, h),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def write_counts_tsv(counts: Counts, out_tsv: Path,
                     log_precursor_total: int | None,
                     log_protein_total: int | None) -> None:
    rows = [
        ("Precursors", "Original (Guo 2019)", counts.original_precursors, ""),
        ("Precursors", "quantmsdiann (matrix, ≥1 run)",
         counts.quantmsdiann_precursors,
         f"diannsummary.log total: {log_precursor_total}"
         if log_precursor_total is not None else ""),
        ("Protein groups", "Original (Guo 2019)", counts.original_proteins, ""),
        ("Protein groups", "quantmsdiann (matrix, ≥1 run)",
         counts.quantmsdiann_proteins,
         f"diannsummary.log total: {log_protein_total}"
         if log_protein_total is not None else ""),
    ]
    df = pd.DataFrame(rows, columns=["metric", "source", "count", "note"])
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Summary-log parsing (best-effort sanity check)
# ---------------------------------------------------------------------------

def parse_summary_log(log_path: Path) -> tuple[int | None, int | None]:
    """Return (precursors_at_1pct_global_q, protein_groups_at_1pct_global_q).
    Either may be None if not found."""
    pr, pg = None, None
    for line in log_path.read_text().splitlines():
        if "Target precursors at 1% global q-value" in line:
            pr = int(line.rsplit(":", 1)[1].strip())
        elif "Protein groups with global q-value <= 0.01" in line:
            pg = int(line.rsplit(":", 1)[1].strip())
    return pr, pg


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    missing = [name for name, val in [
        ("ORIGINAL_PRECURSORS", ORIGINAL_PRECURSORS),
        ("ORIGINAL_PROTEINS", ORIGINAL_PROTEINS),
    ] if val is None]
    if missing:
        print("ERROR: baseline constants not set: " + ", ".join(missing),
              file=sys.stderr)
        print("Set them at the top of this script (see Task 2 in the plan).",
              file=sys.stderr)
        return 2

    pr_path = download_if_missing(PR_MATRIX_URL, DATA_DIR / "diann_report.pr_matrix.tsv")
    pg_path = download_if_missing(PG_MATRIX_URL, DATA_DIR / "diann_report.pg_matrix.tsv")
    log_path = download_if_missing(SUMMARY_LOG_URL, DATA_DIR / "diannsummary.log")

    n_precursors = count_quantified_rows(pr_path)
    n_proteins = count_quantified_rows(pg_path)
    log_pr, log_pg = parse_summary_log(log_path)

    counts = Counts(
        original_precursors=ORIGINAL_PRECURSORS,  # type: ignore[arg-type]
        original_proteins=ORIGINAL_PROTEINS,      # type: ignore[arg-type]
        quantmsdiann_precursors=n_precursors,
        quantmsdiann_proteins=n_proteins,
    )

    out_pdf = FIGURES_DIR / "PXD003539_original_vs_quantmsdiann.pdf"
    out_png = FIGURES_DIR / "PXD003539_original_vs_quantmsdiann.png"
    out_tsv = FIGURES_DIR / "PXD003539_counts.tsv"

    render_figure(counts, out_pdf, out_png)
    write_counts_tsv(counts, out_tsv, log_pr, log_pg)

    print(f"precursors: original={counts.original_precursors:,} "
          f"quantmsdiann={counts.quantmsdiann_precursors:,}")
    print(f"proteins:   original={counts.original_proteins:,} "
          f"quantmsdiann={counts.quantmsdiann_proteins:,}")
    if log_pr is not None and abs(log_pr - n_precursors) / log_pr > 0.01:
        print(f"WARN: matrix precursor count {n_precursors} differs from "
              f"diannsummary.log total {log_pr} by >1%", file=sys.stderr)
    if log_pg is not None and n_proteins < log_pg:
        print(f"WARN: matrix protein count {n_proteins} is below the "
              f"diannsummary.log total {log_pg}", file=sys.stderr)

    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")
    print(f"wrote {out_tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

(Delete the old `main` placeholder and the existing `if __name__ == "__main__"` block; the replacement above ends with the same guard.)

- [ ] **Step 2: Re-run the existing unit tests to make sure nothing broke**

```bash
pytest analysis/tests -v
```
Expected: both `count_quantified_rows` tests still pass.

- [ ] **Step 3: Commit**

```bash
git add analysis/figure_original_vs_quantmsdiann.py
git commit -m "feat(analysis): wire download, render, and TSV output into main"
```

---

### Task 6: End-to-end run and verification

**Files:** none modified — this task runs the script and inspects outputs.

- [ ] **Step 1: Run the script**

```bash
cd /Users/yperez/work/articles/quantmsdiann
python analysis/figure_original_vs_quantmsdiann.py
```

Expected stdout (numbers in the second column will match Task 2's confirmed values; the third column should be in the same ballpark as the `diannsummary.log` totals — 117,720 precursors and 6,927 protein groups, with proteins likely a bit higher because the pg_matrix is at 5 % q-value):

```
precursors: original=<N> quantmsdiann=~117,000
proteins:   original=<M> quantmsdiann=~7,000–9,000
wrote /Users/yperez/work/articles/quantmsdiann/analysis/figures/PXD003539_original_vs_quantmsdiann.pdf
wrote /Users/yperez/work/articles/quantmsdiann/analysis/figures/PXD003539_original_vs_quantmsdiann.png
wrote /Users/yperez/work/articles/quantmsdiann/analysis/figures/PXD003539_counts.tsv
```

If the `ORIGINAL_*` constants are still `None`, the script exits 2 with a clear message — finish Task 2 first.

- [ ] **Step 2: Verify outputs exist and look sane**

```bash
ls -la analysis/figures/
cat analysis/figures/PXD003539_counts.tsv
```
Expected: three files (`.pdf`, `.png`, `.tsv`) and a 4-row TSV with the original and quantmsdiann counts plus the cross-check note for each row.

- [ ] **Step 3: Eyeball the PNG**

Open `analysis/figures/PXD003539_original_vs_quantmsdiann.png` and confirm:
1. Two grouped bars per metric, original on the left (grey), quantmsdiann on the right (blue).
2. Numeric labels on top of each bar match the TSV.
3. quantmsdiann > original for both metrics. If not, stop — either the baseline numbers or the counting rule is wrong, and the comparison is not yet meaningful.

- [ ] **Step 4: Commit the rendered figure**

```bash
git add analysis/figures/PXD003539_original_vs_quantmsdiann.pdf analysis/figures/PXD003539_original_vs_quantmsdiann.png analysis/figures/PXD003539_counts.tsv
git commit -m "feat(analysis): add first reanalysis figure (PXD003539 original vs quantmsdiann)"
```

---

## Self-review notes

- Spec coverage: every section of the design doc maps to a task (scaffolding → 1, baselines → 2, counting rule → 3+4, downloader → 4, figure + TSV → 5, end-to-end validation → 6).
- No placeholders: each step contains either a concrete command, a concrete code block, or a concrete file edit. The two `ORIGINAL_*` constants start at `None` deliberately (the script enforces they must be set) and are filled by Task 2 with explicit citations.
- Type consistency: `Counts` dataclass field names are reused unchanged in `render_figure`, `write_counts_tsv`, and `main`; `count_quantified_rows`, `_sample_columns`, `download_if_missing`, and `parse_summary_log` signatures match how they're called.
