"""fig3_material_independence.py -- PRL Fig. 3: material independence +
blind prediction.

(a) Three-class scaling collapse. Every system is rescaled by its own
    reference point, which removes the material-specific prefactor and Gamma
    scale but not the exponent: diamond and non-diamond realizations of a
    class fall on one master curve, and the three classes separate into three
    slopes. Every class now has both hosts -- group-IV / SC generic (n=1),
    NV 0<->-1 / SC protected (n=2), NV -1<->+1 / three-mode chain (n=3).
    The inset compensates by Gamma^(n+delta) for delta = 0, +-1: flat only at
    the correct integer, so the panel asserts THIS integer, not merely a power
    law (Gate C data).
(b) Superconducting blind-prediction witness: transfer efficiency |K|^2 for
    generic (kappa^-2) and protected (kappa^-4), with the blind-predicted
    integer power laws overlaid (Gate B data).

No new physics: reuses group_iv_full.py (Gate C) and model_sc_transfer.py
(Gate B) unchanged.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
import _paths  # noqa: E402
_paths.add_sources("nogo", "gate_b", "gate_c", "phase")

import prl_style as sty  # noqa: E402
import nv_reduced_kernel as nvk  # noqa: E402
import group_iv_full as gf  # noqa: E402
import model_sc_transfer as sc  # noqa: E402
import chain3_witness as c3  # noqa: E402

_NV_H = nvk.H_3E()


def _systems():
    """(class n, label, host, Gamma grid, |K|) for every realization.

    Class 3 carries a non-diamond realization too (the three-mode chain of
    Gate C); it used to be NV only, which left the steepest class resting on a
    single material.
    """
    ks_giv = np.logspace(4, 8, 40)
    ks_sc = np.logspace(5, 9, 40)
    ks_nv = np.logspace(2, 5, 40)
    ks_ch = np.logspace(2, 5, 40)
    return [
        (1, "group-IV SiV", "diamond", ks_giv,
         np.abs([gf.full_response(k, "SiV") for k in ks_giv])),
        (1, "SC generic", "non-diamond", ks_sc,
         np.abs([sc.transfer_kernel(k, tuning="generic") for k in ks_sc])),
        (2, "NV $0\\!\\leftrightarrow\\!-1$", "diamond", ks_nv,
         np.abs(nvk.kernel(_NV_H, (0, -1), ks_nv))),
        (2, "SC protected", "non-diamond", ks_sc,
         np.abs([sc.transfer_kernel(k, tuning="protected") for k in ks_sc])),
        (3, "NV $-1\\!\\leftrightarrow\\!+1$", "diamond", ks_nv,
         np.abs(nvk.kernel(_NV_H, (-1, 1), ks_nv))),
        (3, "3-mode chain", "non-diamond", ks_ch,
         np.abs([c3.kernel(k) for k in ks_ch])),
    ]


_CLASS_COLOR = {1: "blue", 2: "vermillion", 3: "green"}


def _panel_a(ax):
    """Scaling collapse onto three master curves, one per integer class.

    Each system is rescaled by its own reference point (Gamma_ref at the
    geometric centre of its window, and |K(Gamma_ref)|), which removes the
    material-specific prefactor and the material-specific Gamma scale but NOT
    the exponent. Diamond and non-diamond realizations of the same class then
    fall on one line, and the three classes separate into three slopes -- which
    is the material-independence claim.

    The earlier version divided each curve by its own plateau of Gamma^n |K|,
    so every curve was flattened to 1 by construction and the panel carried no
    visual information; the exponent-specific content is kept in the inset.
    """
    for n, label, host, ks, K in _systems():
        color = sty.COLORS[_CLASS_COLOR[n]]
        j = len(ks) // 2
        x = ks / ks[j]
        y = K / K[j]
        marker = "o" if host == "diamond" else "s"
        ax.loglog(x, y, marker, ms=2.6, color=color, mfc="none", mew=0.8,
                  label=f"{label} ({host}), $n{{=}}{n}$")

    xs = np.array([1e-2, 1e2])
    for n in (1, 2, 3):
        ax.loglog(xs, xs ** (-float(n)), "-", lw=0.8,
                  color=sty.COLORS[_CLASS_COLOR[n]], alpha=0.55, zorder=0)
        ax.annotate(rf"$\Gamma^{{-{n}}}$", xy=(xs[1], xs[1] ** (-float(n))),
                    xytext=(-2, 2), textcoords="offset points",
                    fontsize=6.0, color=sty.COLORS[_CLASS_COLOR[n]],
                    ha="right", va="bottom")

    ax.set_xlabel(r"$\Gamma/\Gamma_{\rm ref}$")
    ax.set_ylabel(r"$|K(\Gamma)|\,/\,|K(\Gamma_{\rm ref})|$")
    ax.set_xlim(2e-2, 5e1)
    ax.legend(loc="lower left", fontsize=5.6, handlelength=1.1, ncol=1)
    _inset_exponent_specificity(ax)


def _inset_exponent_specificity(ax):
    """Compensated curves at n and at n+-1.

    Gamma^n |K| is flat only for the correct integer; Gamma^(n+-1) |K| tilts
    with slope +-1. So the panel asserts not just "a power law" but "this
    integer and not its neighbours".
    """
    ins = ax.inset_axes([0.60, 0.62, 0.38, 0.35])
    for n, label, host, ks, K in _systems():
        color = sty.COLORS[_CLASS_COLOR[n]]
        x = ks / ks[len(ks) // 2]
        for dn, alpha, lw in ((0, 1.0, 0.9), (-1, 0.28, 0.6), (1, 0.28, 0.6)):
            comp = ks ** (n + dn) * K
            comp = comp / np.mean(comp[len(comp) // 3: 2 * len(comp) // 3])
            ins.loglog(x, comp, "-", lw=lw, color=color, alpha=alpha)
    ins.set_ylim(1e-2, 1e2)
    ins.set_xlabel(r"$\Gamma/\Gamma_{\rm ref}$", fontsize=5.2, labelpad=1)
    ins.set_ylabel(r"$\Gamma^{n+\delta}|K|$", fontsize=5.2, labelpad=1)
    ins.tick_params(labelsize=4.6, pad=1)
    ins.text(0.03, 0.06, r"solid $\delta{=}0$; faint $\delta{=}\pm1$",
             transform=ins.transAxes, fontsize=4.6)


def _panel_b(ax):
    """SC blind-prediction: |K|^2 generic (kappa^-2) vs protected (kappa^-4)."""
    ks = np.logspace(3, 9, 60)
    eff_g = np.abs([sc.transfer_kernel(k, tuning="generic") for k in ks]) ** 2
    eff_p = np.abs([sc.transfer_kernel(k, tuning="protected") for k in ks]) ** 2

    ax.loglog(ks, eff_g, "o", ms=2.5, color=sty.COLORS["blue"], mfc="none",
             mew=0.8, label=r"generic, blind pred. $\nu_{\rm obs}=2$")
    ax.loglog(ks, eff_p, "s", ms=2.5, color=sty.COLORS["vermillion"], mfc="none",
             mew=0.8, label=r"protected, blind pred. $\nu_{\rm obs}=4$")

    ref_g = eff_g[-1] * (ks[-1] / ks) ** 2
    ref_p = eff_p[-1] * (ks[-1] / ks) ** 4
    ax.loglog(ks, ref_g, "--", color=sty.COLORS["blue"], linewidth=0.9, alpha=0.7)
    ax.loglog(ks, ref_p, "--", color=sty.COLORS["vermillion"], linewidth=0.9, alpha=0.7)

    ax.set_xlabel(r"$\kappa$ (bus decay)")
    ax.set_ylabel(r"transfer efficiency $|K(\kappa)|^2$")
    ax.legend(loc="lower left", fontsize=6.2, handlelength=1.3)


def build(quick=False):
    sty.apply_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(sty.COL_2, 2.7))
    _panel_a(ax1)
    _panel_b(ax2)
    sty.panel_label(ax1, "a")
    sty.panel_label(ax2, "b")
    fig.tight_layout(w_pad=2.0)
    return sty.save(fig, "fig3_material_independence")


if __name__ == "__main__":
    print(build())
