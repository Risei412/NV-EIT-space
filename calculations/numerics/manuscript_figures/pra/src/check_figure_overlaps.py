"""Detect text collisions inside the manuscript figures.

Legends, axis labels, panel labels, colorbar labels and titles are all drawn
independently, so nothing stops matplotlib from putting one on top of another.
This walks every figure, asks the renderer for the actual drawn extent of each
text-bearing artist, and reports overlapping pairs -- the same check the eye
would do, but exhaustive and repeatable.

Usage: python check_figure_overlaps.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matplotlib                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402

import make_figures as mf              # noqa: E402

TOL = 1.0          # pixels of overlap tolerated


def artists_of(fig):
    """(label, bbox) for every text-bearing artist actually drawn."""
    r = fig.canvas.get_renderer()
    out = []
    for ax in fig.axes:
        for name, art in (("xlabel", ax.xaxis.label),
                          ("ylabel", ax.yaxis.label),
                          ("title", ax.title)):
            if art.get_text().strip():
                out.append((f"{_axid(fig, ax)}.{name}"
                            f" {art.get_text()[:24]!r}",
                            art.get_window_extent(r)))
        leg = ax.get_legend()
        if leg is not None:
            out.append((f"{_axid(fig, ax)}.legend",
                        leg.get_window_extent(r)))
        for t in ax.texts:
            if t.get_text().strip():
                out.append((f"{_axid(fig, ax)}.text {t.get_text()[:22]!r}",
                            t.get_window_extent(r)))
        for nm, off in (("xoffset", ax.xaxis.get_offset_text()),
                        ("yoffset", ax.yaxis.get_offset_text())):
            if off.get_visible() and off.get_text().strip():
                out.append((f"{_axid(fig, ax)}.{nm} {off.get_text()[:12]!r}",
                            off.get_window_extent(r)))
        for nm, ticks in (("xtick", ax.get_xticklabels()),
                          ("ytick", ax.get_yticklabels())):
            for t in ticks:
                if t.get_text().strip():
                    out.append((f"{_axid(fig, ax)}.{nm} {t.get_text()[:8]!r}",
                                t.get_window_extent(r)))
    for t in fig.texts:
        if t.get_text().strip():
            out.append((f"fig.text {t.get_text()[:22]!r}",
                        t.get_window_extent(r)))
    if fig.legends:
        for i, leg in enumerate(fig.legends):
            out.append((f"fig.legend{i}", leg.get_window_extent(r)))
    return out


def _axid(fig, ax):
    return f"ax{fig.axes.index(ax)}"


def graphics_of(fig):
    """(label, bbox) for arrows and rules drawn on *schematic* axes.

    Only axes with the frame switched off are inspected.  On a data plot,
    an annotation sitting on a curve is normal; on a level diagram, a label
    sitting on an arrow is the defect this pass exists to find.
    """
    r = fig.canvas.get_renderer()
    out = []
    for ax in fig.axes:
        if ax.axison:
            continue
        for i, ln in enumerate(ax.lines):
            out.append((f"{_axid(fig, ax)}.line{i}", ln.get_window_extent(r)))
        for i, p in enumerate(ax.patches):
            out.append((f"{_axid(fig, ax)}.patch{i}",
                        p.get_window_extent(r)))
        for i, a in enumerate(ax.texts):
            # annotate() arrows live on the Annotation object
            arr = getattr(a, "arrow_patch", None)
            if arr is not None:
                out.append((f"{_axid(fig, ax)}.arrow{i}",
                            arr.get_window_extent(r)))
    return out


PANEL = ("(a)", "(b)", "(c)")


def _benign(n1, n2):
    """Pairs whose proximity is by design, not a defect."""
    pair = (n1, n2)
    # panel letters are deliberately placed in the margin over the tick strip
    if any(f"text {p!r}" in n for p in PANEL for n in pair):
        if any(".xtick" in n or ".ytick" in n for n in pair):
            return True
    # an axis label necessarily sits just outside its own tick labels
    for a, b in (pair, pair[::-1]):
        if a.endswith("label") or ".xlabel" in a or ".ylabel" in a:
            if ".xtick" in b or ".ytick" in b:
                if a.split(".")[0] == b.split(".")[0]:
                    return True
    return False


def overlap(a, b):
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return dx, dy


def check(name, builder):
    """Build the figure and inspect it BEFORE the builder's _save closes it."""
    plt.close("all")
    captured = []
    real_save = mf._save

    def spy(fig, nm):
        fig.canvas.draw()
        captured.append(fig)
        return None          # skip the write; we only want the geometry

    mf._save = spy
    try:
        builder()
    finally:
        mf._save = real_save

    issues = 0
    for fig in captured:
        items = artists_of(fig)
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                (n1, b1), (n2, b2) = items[i], items[j]
                if _benign(n1, n2):
                    continue
                dx, dy = overlap(b1, b2)
                if dx > TOL and dy > TOL:
                    print(f"  {name}: {n1}  X  {n2}   "
                          f"({dx:.0f} x {dy:.0f} px)")
                    issues += 1
        # text over arrows/rules on schematic axes
        graphics = graphics_of(fig)
        for n1, b1 in items:
            if ".text " not in n1:
                continue
            for n2, b2 in graphics:
                if n1.split(".")[0] != n2.split(".")[0]:
                    continue
                dx, dy = overlap(b1, b2)
                if dx > TOL and dy > TOL:
                    print(f"  {name}: {n1}  X  {n2}   "
                          f"({dx:.0f} x {dy:.0f} px)")
                    issues += 1
    plt.close("all")
    return issues


def main():
    mf.st.apply_style()
    checks = [
        ("fig1", mf.fig1_level_scheme),
        ("fig2", lambda: mf.fig2_spectra(False)),
        ("fig3", mf.fig3_phase_diagram),
        ("fig4", mf.fig4_contrast_vs_T),
        ("fig5", mf.fig5_bperp_scaling),
        ("fig6", mf.fig6_observables),
    ]
    total = 0
    for name, fn in checks:
        print(f"--- {name}")
        try:
            total += check(name, fn)
        except Exception as exc:
            print(f"  {name}: FAILED to build ({type(exc).__name__}: {exc})")
            total += 1
    print(f"\n{total} text-collision(s)")
    return total


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
