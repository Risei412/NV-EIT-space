"""prl_style.py -- shared matplotlib style for the PRL main-text figures.

Reuses no physics; a pure presentation module. Colorblind-safe categorical
palette (Okabe-Ito), validated with the dataviz skill's validate_palette.js:
all six checks PASS (worst adjacent CVD deltaE 9.6, normal-vision floor 16.4).
Contrast-vs-surface is WARN for the lighter entries, so every mark that uses
one of those colors is paired with a direct label or marker-shape encoding
(never color alone) per the skill's non-negotiables.

Column widths follow APS/PRL page geometry (single column 3.375in, double
column 6.75in). No LaTeX dependency (mathtext only) so figures render in any
environment.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe categorical palette (validated, see docstring).
COLORS = {
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "blue": "#0072B2",
    "black": "#000000",
    "gray": "#7F7F7F",
}

# fixed categorical assignment order (never cycled/reassigned per-plot)
CAT_ORDER = ["blue", "vermillion", "green", "orange", "purple", "sky"]

MARKERS = ["o", "s", "^", "D", "v", "P"]

COL_1 = 3.375  # inches, PRL single column
COL_2 = 6.75   # inches, PRL double column

FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")


def apply_style():
    plt.rcParams.update({
        "font.size": 8.5,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7,
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "axes.edgecolor": "#333333",
        "axes.labelcolor": "#111111",
        "text.color": "#111111",
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "axes.grid": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "mathtext.fontset": "dejavuserif",
        "font.family": "serif",
        "legend.frameon": False,
        "legend.handlelength": 1.6,
    })


def panel_label(ax, letter, x=-0.16, y=1.04):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="right")


def save(fig, name):
    os.makedirs(FIGURES_DIR, exist_ok=True)
    pdf = os.path.join(FIGURES_DIR, f"{name}.pdf")
    png = os.path.join(FIGURES_DIR, f"{name}.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    return pdf, png
