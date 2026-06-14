"""Fig. 1 — quantmsdiann pipeline architecture and runtime scaling.

Composite MCP-ready panel (7.2 x 6.0 in) built with matplotlib GridSpec:

  (a) left column  — subway workflow diagram (rasterised from source SVG,
                     full canvas, no crop, natural aspect ratio preserved)
  (b) right top    — cluster queueSize sweep (PXD071075)
  (c) right bottom — wallclock vs MS run count across cohorts

Output: analysis/figures/manuscript/fig1_architecture_scaling.svg
"""
from __future__ import annotations

import io
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analysis import figure_style as fs
fs.apply_house_style()
import numpy as np
import pandas as pd
from PIL import Image

from analysis.figure_performance_trace import (
    collect_parallelism_rows,
    render_parallelism_scatter,
)
from analysis.figure_queue_size_sweep import (
    collect_sweep_rows,
    render_queue_size_sweep,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = REPO_ROOT / "analysis" / "figures" / "manuscript"
WORKFLOW_SRC = REPO_ROOT / "paper" / "figures" / "source" / "quantmsdiann_workflow.svg"
WORKFLOW_FALLBACK_URL = (
    "https://raw.githubusercontent.com/bigbio/quantmsdiann/main/"
    "docs/images/quantmsdiann_workflow.svg"
)

# MCP full-page figure size (inches).
# Enlarged so the embedded workflow schematic (panel a) is legible at print
# size; the prior 7.2x6.0 shrank its node labels below readable.
FIG_W = 11.5
FIG_H = 7.4

# DPI used when ImageMagick rasterises the workflow SVG. No -crop is
# applied — the prior hard-coded "535x782+85+38" slice was tuned for an
# earlier 880x830-pt source canvas at 72 DPI and silently truncated the
# subway diagram once the script switched to -density 300.
WORKFLOW_RASTER_DPI = 300


def ensure_workflow_svg() -> Path:
    if WORKFLOW_SRC.exists() and WORKFLOW_SRC.stat().st_size > 0:
        return WORKFLOW_SRC
    WORKFLOW_SRC.parent.mkdir(parents=True, exist_ok=True)
    try:
        import urllib.request
        with urllib.request.urlopen(WORKFLOW_FALLBACK_URL, timeout=30) as resp:
            WORKFLOW_SRC.write_bytes(resp.read())
    except OSError as exc:
        raise FileNotFoundError(
            f"Workflow SVG not found at {WORKFLOW_SRC} and download failed: {exc}"
        ) from exc
    return WORKFLOW_SRC


def _annotate_panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.06, 1.03, letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def _workflow_raster(workflow_svg: Path, *, dpi: int = WORKFLOW_RASTER_DPI) -> np.ndarray:
    """Rasterise the full workflow SVG to an RGB numpy array. No crop is
    applied — keeping the diagram intact is the whole point of the panel.

    Uses rsvg-convert, NOT ImageMagick: the source SVG styles its connecting
    lines via a CSS `<style>` block (class selectors), which ImageMagick's SVG
    renderer ignores -- so `magick` dropped every edge and only the nodes
    survived. rsvg-convert honours the CSS, so the lines render."""
    cmd = [
        "rsvg-convert",
        "--background-color", "white",
        "--dpi-x", str(dpi), "--dpi-y", str(dpi),
        str(workflow_svg),
    ]
    result = subprocess.run(cmd, capture_output=True, check=True)
    with Image.open(io.BytesIO(result.stdout)) as img:
        return np.asarray(img.convert("RGB"))


def _draw_workflow_panel(ax: plt.Axes, workflow_svg: Path) -> None:
    img = _workflow_raster(workflow_svg)
    # aspect="equal" preserves the SVG's native aspect ratio inside the
    # panel; the prior aspect="auto" stretched the diagram to whatever
    # the GridSpec column happened to be, distorting line spacing and
    # blowing up text when combined with the legacy crop.
    ax.imshow(img, interpolation="lanczos", aspect="equal", zorder=0)
    ax.set_axis_off()


def _parallelism_plot_frame(fetch: bool = True) -> pd.DataFrame:
    par_df = collect_parallelism_rows(fetch=fetch)
    from analysis.figure_performance_trace import collect_pxd071075_sweep_rows
    sweep_df = collect_pxd071075_sweep_rows()
    if not sweep_df.empty:
        sweep_top = sweep_df[sweep_df["queue_size"] == 300]
        par_df = pd.concat([par_df, sweep_top], ignore_index=True, sort=False)
    return par_df


def compose_fig1(
    *,
    workflow_svg: Path,
    sweep_df: pd.DataFrame,
    par_df: pd.DataFrame,
) -> plt.Figure:
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white", layout="constrained")
    # Tighten the constrained-layout padding so b/c sit close to a.
    fig.get_layout_engine().set(w_pad=0.01, h_pad=0.03, wspace=0.0, hspace=0.04)
    # Left column wide enough to give the subway diagram its native
    # 880x830 aspect ratio without crowding the right-column scatters.
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.95, 1.0],  # panel a (schematic) wider/bigger
        height_ratios=[1.0, 1.0],
        wspace=0.0,                # pull b/c in close to a
    )

    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    _draw_workflow_panel(ax_a, workflow_svg)
    render_queue_size_sweep(sweep_df, ax=ax_b, composite=True)
    render_parallelism_scatter(
        par_df,
        ax=ax_c,
        legend_ncol=2,
        legend_bbox_y=-0.32,
        composite=True,
        show_legend=True,
    )

    _annotate_panel_letter(ax_a, "a")
    _annotate_panel_letter(ax_b, "b")
    _annotate_panel_letter(ax_c, "c")

    return fig


def render_fig1_architecture_scaling(
    out_path: Path,
    *,
    workflow_svg: Path | None = None,
    fetch: bool = True,
) -> None:
    fig = compose_fig1(
        workflow_svg=workflow_svg or ensure_workflow_svg(),
        sweep_df=collect_sweep_rows(),
        par_df=_parallelism_plot_frame(fetch=fetch),
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    stem = out_path.with_suffix("")
    save_kw = {"facecolor": "white", "dpi": 300}
    fig.savefig(stem.with_suffix(".svg"), **save_kw)
    plt.close(fig)


def main() -> int:  # pragma: no cover
    # SVG only — the manuscript PDF is produced by the paper Makefile
    # (rsvg-convert of the SVG into paper/figures/pdf/), not here.
    out = FIGURES_DIR / "fig1_architecture_scaling.svg"
    render_fig1_architecture_scaling(out)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
