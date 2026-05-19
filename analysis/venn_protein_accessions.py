from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis.figure_original_vs_quantmsdiann import (
    unique_peptides_per_protein_diann,
    unique_peptides_per_protein_openswath,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "PXD003539"
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "PXD003539"

_ACCESSION_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:[A-Z0-9]{1,5})?$")
_PREFIX_RE = re.compile(r"^(?:CONTAM_|ENTRAP_|DECOY_)+")


def _clean_token(tok: str) -> str | None:
    tok = tok.strip()
    if not tok:
        return None
    if "|" in tok:
        # sp|P12345|HUMAN -> P12345 (middle field)
        parts = tok.split("|")
        if len(parts) >= 2 and parts[1]:
            tok = parts[1]
        else:
            return None
    tok = _PREFIX_RE.sub("", tok)
    if "-" in tok:
        tok = tok.split("-", 1)[0]
    if not tok:
        return None
    return tok


def extract_accessions_diann(protein_group: str | None) -> set[str]:
    if not protein_group:
        return set()
    out: set[str] = set()
    for piece in protein_group.split(";"):
        cleaned = _clean_token(piece)
        if cleaned:
            out.add(cleaned)
    return out


def extract_accessions_openswath(protein_str: str | None) -> set[str]:
    if not protein_str:
        return set()
    parts = protein_str.split("/")
    # First field is the count (e.g. "1", "2", "DECOY_1") — skip it.
    out: set[str] = set()
    for piece in parts[1:]:
        cleaned = _clean_token(piece)
        if cleaned:
            out.add(cleaned)
    return out


def accessions_with_min_peptides_diann(
    pr_matrix_path: Path, *, min_peptides: int = 2
) -> set[str]:
    counts = unique_peptides_per_protein_diann(pr_matrix_path)
    out: set[str] = set()
    for pg, n in counts.items():
        if n >= min_peptides:
            out.update(extract_accessions_diann(pg))
    return out


def accessions_with_min_peptides_openswath(
    matrix_path: Path, *, min_peptides: int = 2
) -> set[str]:
    counts = unique_peptides_per_protein_openswath(matrix_path)
    out: set[str] = set()
    for prot, n in counts.items():
        if n >= min_peptides:
            out.update(extract_accessions_openswath(prot))
    return out


def render_venn_diagram(
    guo_acc: set[str],
    diann_acc: set[str],
    pdf_path: Path,
    png_path: Path,
    svg_path: Path | None = None,
    *,
    left_label: str = "Guo 2019 (OpenSWATH)",
    right_label: str = "quantmsdiann (DIA-NN)",
    left_color: str = "#9e9e9e",
    right_color: str = "#1f77b4",
    # title and footer kept as parameters for backwards-compatible callers,
    # but ignored: paper-ready figures carry no title or footnote.
    title: str | None = None,
    footer: str | None = None,
) -> None:
    full_left = f"{left_label}\n(n={len(guo_acc):,})"
    full_right = f"{right_label}\n(n={len(diann_acc):,})"
    inter = guo_acc & diann_acc
    guo_only = guo_acc - diann_acc
    diann_only = diann_acc - guo_acc
    total = len(guo_acc | diann_acc) or 1

    fig, ax = plt.subplots(figsize=(7, 5.5))

    try:
        from matplotlib_venn import venn2
        v = venn2(
            subsets=(len(guo_only), len(diann_only), len(inter)),
            set_labels=(full_left, full_right),
            set_colors=(left_color, right_color),
            alpha=0.55,
            ax=ax,
        )
        labels = {
            "10": (len(guo_only), guo_only),
            "01": (len(diann_only), diann_only),
            "11": (len(inter), inter),
        }
        for region_id, (count, _) in labels.items():
            lbl = v.get_label_by_id(region_id)
            if lbl is None:
                continue
            pct = 100.0 * count / total
            lbl.set_text(f"{count:,}\n({pct:.1f}%)")
            lbl.set_fontsize(10)
        for sl in v.set_labels:
            if sl is not None:
                sl.set_fontsize(11)
    except ImportError:
        from matplotlib.patches import Circle
        ax.set_xlim(-2, 4)
        ax.set_ylim(-2, 2)
        ax.set_aspect("equal")
        ax.axis("off")
        c_left = Circle((0, 0), 1.3, color=left_color, alpha=0.55, linewidth=0)
        c_right = Circle((1.4, 0), 1.3, color=right_color, alpha=0.55, linewidth=0)
        ax.add_patch(c_left)
        ax.add_patch(c_right)
        ax.text(-0.9, 0, f"{len(guo_only):,}\n({100*len(guo_only)/total:.1f}%)",
                ha="center", va="center", fontsize=10)
        ax.text(2.3, 0, f"{len(diann_only):,}\n({100*len(diann_only)/total:.1f}%)",
                ha="center", va="center", fontsize=10)
        ax.text(0.7, 0, f"{len(inter):,}\n({100*len(inter)/total:.1f}%)",
                ha="center", va="center", fontsize=10)
        ax.text(-0.6, 1.5, full_left, ha="center", va="bottom", fontsize=11)
        ax.text(2.0, 1.5, full_right, ha="center", va="bottom", fontsize=11)

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    if svg_path is not None:
        fig.savefig(svg_path)
    plt.close(fig)


def main() -> int:
    pr_path = DATA_DIR / "diann_report.pr_matrix.tsv"
    opensw_path = DATA_DIR / "feature_alignment_requant_matrix.tsv"
    if not pr_path.exists():
        print(f"Missing input: {pr_path}", file=sys.stderr)
        return 1
    if not opensw_path.exists():
        print(f"Missing input: {opensw_path}", file=sys.stderr)
        return 1

    print("Computing DIA-NN accession set (>=2 unique peptides)...")
    diann_acc = accessions_with_min_peptides_diann(pr_path, min_peptides=2)
    print("Computing Guo/OpenSWATH accession set (>=2 unique peptides)...")
    guo_acc = accessions_with_min_peptides_openswath(opensw_path, min_peptides=2)

    inter = guo_acc & diann_acc
    guo_only = guo_acc - diann_acc
    diann_only = diann_acc - guo_acc

    print(f"Guo total:        {len(guo_acc):,}")
    print(f"DIA-NN total:     {len(diann_acc):,}")
    print(f"Intersection:     {len(inter):,}")
    print(f"Guo only:         {len(guo_only):,}")
    print(f"DIA-NN only:      {len(diann_only):,}")

    pdf_path = FIGURES_DIR / "supp_venn_protein_accessions.pdf"
    png_path = FIGURES_DIR / "supp_venn_protein_accessions.png"
    svg_path = FIGURES_DIR / "supp_venn_protein_accessions.svg"
    render_venn_diagram(guo_acc, diann_acc, pdf_path, png_path, svg_path)
    print(f"Figure saved to {pdf_path}, {png_path}, {svg_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
