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

from __future__ import annotations

import matplotlib as mpl
from matplotlib.ticker import FuncFormatter

# --- Okabe-Ito colour-blind-safe palette -----------------------------------
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "grey": "#9E9E9E",
}

# DIA-NN version palette (oldest -> newest = light blue -> dark blue; the
# enterprise build is a distinct vermillion accent). Keyed by the version
# tokens the scripts already use (with and without the ``v`` prefix).
VERSION_COLORS = {
    "v1_8_1": OKABE_ITO["sky_blue"],
    "1_8_1": OKABE_ITO["sky_blue"],
    "v2_5_1": OKABE_ITO["blue"],
    "2_5_1": OKABE_ITO["blue"],
    "v2_5_1_enterprise": OKABE_ITO["vermillion"],
    "2_5_1_enterprise": OKABE_ITO["vermillion"],
}

# Two-way comparison: published baseline vs our reanalysis.
COMPARISON = {
    "original": OKABE_ITO["grey"],
    "quantmsdiann": OKABE_ITO["blue"],
}

# Ordered colour-blind-leaning cycle for categorical encodings that need many
# colours (e.g. the per-instrument runtime bars). First 8 are Okabe-Ito; the
# tail extends with a few muted Tol hues. >12 categories will always be hard to
# read — prefer grouping by vendor over adding colours.
CATEGORICAL_CYCLE = [
    OKABE_ITO["blue"],
    OKABE_ITO["vermillion"],
    OKABE_ITO["bluish_green"],
    OKABE_ITO["orange"],
    OKABE_ITO["sky_blue"],
    OKABE_ITO["reddish_purple"],
    "#882255",  # wine
    "#117733",  # forest
    "#888888",  # mid grey
    "#332288",  # indigo
    "#999933",  # olive
    "#44AA99",  # teal
]


# Stable, colour-blind-leaning instrument palette, grouped by vendor family
# (Orbitrap = blues, SCIEX = warm/green, Bruker timsTOF = purple/wine/olive) so
# the same instrument is the same colour in every figure that encodes it.
INSTRUMENT_COLORS = {
    "Orbitrap Astral": OKABE_ITO["blue"],
    "Orbitrap Eclipse": "#332288",            # indigo
    "Orbitrap Exploris 480": OKABE_ITO["sky_blue"],
    "Q Exactive": "#44AA99",                  # teal
    "TripleTOF 5600": OKABE_ITO["orange"],
    "TripleTOF 6600": OKABE_ITO["vermillion"],
    "ZenoTOF 7600": OKABE_ITO["bluish_green"],
    "timsTOF SCP": OKABE_ITO["reddish_purple"],
    "timsTOF Pro": "#882255",                 # wine
    "timsTOF HT": "#999933",                  # olive
    "unknown": OKABE_ITO["grey"],
}


def categorical_colors(n):
    """Return *n* colours from the colour-blind-leaning categorical cycle."""
    if n <= len(CATEGORICAL_CYCLE):
        return CATEGORICAL_CYCLE[:n]
    # cycle with repetition only if forced; caller should rethink the encoding
    out = list(CATEGORICAL_CYCLE)
    while len(out) < n:
        out.append(CATEGORICAL_CYCLE[len(out) % len(CATEGORICAL_CYCLE)])
    return out[:n]


# --- global rcParams (typography, spines, output) --------------------------
def apply_house_style():
    """Apply the manuscript-wide matplotlib rcParams. Idempotent."""
    mpl.rcParams.update({
        # typography — one clean sans across every figure
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "legend.title_fontsize": 10,
        # spines — drop top/right everywhere (matches the cohort template)
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        # ticks
        "xtick.direction": "out",
        "ytick.direction": "out",
        # legend — frameless, sits cleanly on white
        "legend.frameon": False,
        # figure/output
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "svg.fonttype": "none",  # keep text as text in the SVG
    })


# --- number formatting -----------------------------------------------------
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
        s = f"{v / 1000:.1f}".rstrip("0").rstrip(".")
        return f"{s}k"
    return f"{v:.0f}"


def kfmt_axis(axis):
    """Apply :func:`kfmt` to a matplotlib axis (e.g. ``ax.yaxis``)."""
    axis.set_major_formatter(FuncFormatter(kfmt))


def commafmt(value, _pos=None):
    """Thousands-separated integer: 9542 -> '9,542'. Safe as a FuncFormatter."""
    try:
        return f"{int(round(float(value))):,}"
    except (TypeError, ValueError):
        return str(value)


def despine(ax, top=True, right=True, left=False, bottom=False):
    """Hide the requested spines (top/right by default)."""
    for name, off in (("top", top), ("right", right),
                      ("left", left), ("bottom", bottom)):
        if off:
            ax.spines[name].set_visible(False)


def style_boxplot(bp, color=None, median_color=None):
    """Restyle a matplotlib ``ax.boxplot`` result to the house style.

    Kills the default cyan boxes / red medians / loud fliers in favour of a
    neutral box with an accent median and discreet outliers.
    """
    color = color or OKABE_ITO["blue"]
    median_color = median_color or OKABE_ITO["vermillion"]
    for box in bp.get("boxes", []):
        box.set(color=color, linewidth=1.0)
    for whisk in bp.get("whiskers", []):
        whisk.set(color=color, linewidth=0.9)
    for cap in bp.get("caps", []):
        cap.set(color=color, linewidth=0.9)
    for med in bp.get("medians", []):
        med.set(color=median_color, linewidth=1.4)
    for fl in bp.get("fliers", []):
        fl.set(marker="o", markersize=2.5, markerfacecolor="none",
               markeredgecolor="#999999", alpha=0.6)
