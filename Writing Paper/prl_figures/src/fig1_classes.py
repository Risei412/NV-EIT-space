"""fig1_classes.py -- PRL Fig. 1: general structure + the three integer
suppression classes.

(a) Schematic: input c -> fast dissipative sector D -> internal dynamics A0
    -> readout p, with the selection-rule-cancellation -> first-nonzero-
    moment -> integer-order chain (strategy doc Sec.9, boxed equation).
(b) |K(Gamma)| log-log for the three physical classes n=1,2,3 (Gate C data):
    group-IV SiV (M0 != 0), NV ms=0<->-1 (d=1), NV ms=-1<->+1 (d=2), with
    Gamma^-n reference lines.

No new physics: reuses nv_reduced_kernel and group_iv_full (Gate C) unchanged.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path
from matplotlib.patches import PathPatch

HERE = os.path.dirname(os.path.abspath(__file__))
NOGO_SRC = os.path.join(HERE, "..", "..", "..", "No-go theorem", "src")
GATEC_SRC = os.path.join(HERE, "..", "..", "..", "New no-go theory",
                         "GateC_material_independence", "src")
PHASE_SRC = os.path.join(HERE, "..", "..", "..", "New no-go theory", "src")
for _p in (HERE, NOGO_SRC, GATEC_SRC, PHASE_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import prl_style as sty  # noqa: E402
import nv_reduced_kernel as nvk  # noqa: E402
import group_iv_full as gf  # noqa: E402
import core  # noqa: E402

_NV_H = nvk.H_3E()


def _panel_a(ax):
    """Schematic input -> D -> A0 -> readout, with the selection-rule chain."""
    ax.set_xlim(0, 10); ax.set_ylim(0, 5)
    ax.axis("off")

    boxes = [
        (0.1, 1.6, 1.5, 1.8, "input\n" + r"$c$", sty.COLORS["blue"]),
        (2.15, 0.9, 2.55, 3.2, "fast dissipative\nsector\n" + r"$D$", sty.COLORS["vermillion"]),
        (5.3, 0.9, 2.55, 3.2, "internal\ndynamics\n" + r"$A_0(z)$", sty.COLORS["green"]),
        (8.4, 1.6, 1.5, 1.8, "readout\n" + r"$p$", sty.COLORS["blue"]),
    ]
    for x, y, w, h, label, color in boxes:
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.12",
                             linewidth=1.1, edgecolor=color, facecolor="white", zorder=2)
        ax.add_patch(box)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=6.8, zorder=3)

    arrow_y = 2.5
    pairs = [(1.6, 2.15), (4.7, 5.3), (7.85, 8.4)]
    for x0, x1 in pairs:
        ax.add_patch(FancyArrowPatch((x0, arrow_y), (x1, arrow_y),
                                     arrowstyle="-|>", mutation_scale=10,
                                     color="#333333", linewidth=1.0, zorder=1))

    ax.text(5.0, 0.35,
            r"selection rule $\Rightarrow$ moment cancellation $\Rightarrow$ "
            r"integer order $\Rightarrow$ response class",
            ha="center", va="center", fontsize=7.2, style="italic")
    ax.text(5.0, 4.55,
            r"$K_{pc}(\Gamma,z)=p^\dagger[\Gamma D+A_0(z)]^{-1}c \;\sim\; M_m(z)\,\Gamma^{-(m+1)}$",
            ha="center", va="center", fontsize=8)


def _panel_b(ax):
    """Three integer classes: |K(Gamma)| with Gamma^-n reference lines."""
    ks_giv = np.logspace(4, 8, 40)
    K_giv = np.array([gf.full_response(k, "SiV") for k in ks_giv])

    ks_nv = np.logspace(2, 5, 40)
    K_nv2 = nvk.kernel(_NV_H, (0, -1), ks_nv)
    K_nv3 = nvk.kernel(_NV_H, (-1, 1), ks_nv)

    series = [
        (ks_giv, K_giv, 1, "group-IV (SiV), $M_0\\neq0$", sty.COLORS["blue"], "o"),
        (ks_nv, K_nv2, 2, r"NV $m_s{=}0{\leftrightarrow}{-}1$, $d{=}1$", sty.COLORS["vermillion"], "s"),
        (ks_nv, K_nv3, 3, r"NV $m_s{=}{-}1{\leftrightarrow}{+}1$, $d{=}2$", sty.COLORS["green"], "^"),
    ]
    for ks, K, n, label, color, marker in series:
        y = np.abs(K)
        sl = core.fit_nu_loglog(ks, K)["nu_global"]
        ax.loglog(ks, y, marker, ms=3, color=color, mfc="none", mew=1.0,
                  label=f"{label}: slope $={sl:.2f}$")
        ref = y[len(y) // 2] * (ks[len(ks) // 2] / ks) ** n
        ax.loglog(ks, ref, "--", color=color, linewidth=0.8, alpha=0.6)

    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"$|K(\Gamma)|$")
    ax.set_ylim(1e-16, 1e-1)
    ax.legend(loc="lower left", fontsize=6.3, handlelength=1.3)


def build(quick=False):
    sty.apply_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(sty.COL_2, 2.7))
    _panel_a(ax1)
    _panel_b(ax2)
    sty.panel_label(ax1, "a", x=-0.02, y=1.02)
    sty.panel_label(ax2, "b")
    fig.tight_layout(w_pad=2.0)
    return sty.save(fig, "fig1_classes")


if __name__ == "__main__":
    print(build())
