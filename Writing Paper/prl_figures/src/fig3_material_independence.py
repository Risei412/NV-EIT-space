"""fig3_material_independence.py -- PRL Fig. 3: material independence +
blind prediction.

(a) Three-class scaling collapse, Gamma^n |K|, overlaying diamond (group-IV
    n=1, NV n=2,3) and non-diamond (superconducting generic n=1, protected
    n=2) systems -- the same integer exponent across materials (Gate C data).
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
NOGO_SRC = os.path.join(HERE, "..", "..", "..", "No-go theorem", "src")
GATEB_SRC = os.path.join(HERE, "..", "..", "..", "New no-go theory",
                         "GateB_superconducting_witness", "src")
GATEC_SRC = os.path.join(HERE, "..", "..", "..", "New no-go theory",
                         "GateC_material_independence", "src")
for _p in (HERE, NOGO_SRC, GATEB_SRC, GATEC_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import prl_style as sty  # noqa: E402
import nv_reduced_kernel as nvk  # noqa: E402
import group_iv_full as gf  # noqa: E402
import model_sc_transfer as sc  # noqa: E402

_NV_H = nvk.H_3E()


def _panel_a(ax):
    """Gamma^n |K|, each curve normalized by its own plateau value, overlaying
    diamond (o) and non-diamond (s) systems: the coefficient is material-
    specific but the collapse to a common height ~1 shows the shared integer
    exponent n."""
    c1, c2 = sty.COLORS["blue"], sty.COLORS["vermillion"]

    def norm_plateau(ks, y, n):
        comp = ks ** n * y
        return comp / np.mean(comp[-len(comp) // 3:])

    ks_giv = np.logspace(4, 8, 30)
    K_giv = np.abs([gf.full_response(k, "SiV") for k in ks_giv])
    ax.semilogx(ks_giv, norm_plateau(ks_giv, K_giv, 1), "o", ms=3, color=c1,
               mfc="none", mew=0.9, label="group-IV (diamond), $n{=}1$")

    ks_scg = np.logspace(5, 9, 30)
    K_scg = np.abs([sc.transfer_kernel(k, tuning="generic") for k in ks_scg])
    ax.semilogx(ks_scg, norm_plateau(ks_scg, K_scg, 1), "s", ms=3, color=c1,
               mfc="none", mew=0.9, label="SC generic (non-diamond), $n{=}1$")

    ks_nv2 = np.logspace(2, 5, 30)
    K_nv2 = np.abs(nvk.kernel(_NV_H, (0, -1), ks_nv2))
    ax.semilogx(ks_nv2, norm_plateau(ks_nv2, K_nv2, 2), "o", ms=3, color=c2,
               mfc="none", mew=0.9, label="NV (diamond), $n{=}2$")

    K_scp = np.abs([sc.transfer_kernel(k, tuning="protected") for k in ks_scg])
    ax.semilogx(ks_scg, norm_plateau(ks_scg, K_scp, 2), "s", ms=3, color=c2,
               mfc="none", mew=0.9, label="SC protected (non-diamond), $n{=}2$")

    ax.axhline(1.0, color="#999999", linewidth=0.7, zorder=0)
    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"$\Gamma^{n}\,|K(\Gamma)|$ / plateau value")
    ax.set_ylim(0, 2)
    ax.legend(loc="lower left", fontsize=6.2, handlelength=1.2, ncol=1)


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
