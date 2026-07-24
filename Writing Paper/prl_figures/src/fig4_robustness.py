"""fig4_robustness.py -- PRL Fig. 4: robustness/crossover + experimental
discriminability.

(a) Exact (NV, symmetry-protected, unbreakable) vs approximate (superconducting,
    tuned cancellation) class: nu_eff(Gamma) family for several symmetry-
    breaking eps, and the crossover scale Gamma*(eps) ~ 1/eps (inset).
(b) Gamma(T) platform reach: NV phonon rate (~8.5 decades, 4-300 K) and
    group-IV Bose-law rates (narrow), against the required Gamma dynamic
    range to resolve an adjacent class.

No new physics: reuses run_gate_d.py (Gate D) unchanged.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
GATED_SRC = os.path.join(HERE, "..", "..", "..", "New no-go theory",
                         "GateD_robustness_discriminability", "src")
for _p in (HERE, GATED_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import prl_style as sty  # noqa: E402
import run_gate_d as rd  # noqa: E402


def _panel_a(ax, quick):
    sc_res = rd.sc_approximate_class(quick=quick)
    palette = [sty.COLORS["blue"], sty.COLORS["green"], sty.COLORS["vermillion"],
              sty.COLORS["purple"]]
    for (eps, (gm, ne)), color in zip(sc_res["curves"].items(), palette):
        ax.semilogx(gm, ne, "-", color=color, linewidth=1.1,
                   label=fr"$\epsilon={eps:.0e}$")
    ax.axhline(2, color="#555555", ls="--", linewidth=0.8)
    ax.axhline(1, color="#555555", ls=":", linewidth=0.8)
    ax.text(0.02, 0.06, r"NV $m_s{=}{-}1{\leftrightarrow}{+}1$: exact class, "
            r"$\nu{=}3$ unbroken $\forall$ strain, $B_\perp$",
            transform=ax.transAxes, fontsize=6.0, style="italic", va="bottom")
    ax.set_xlabel(r"$\Gamma$ (bus decay $\kappa$)")
    ax.set_ylabel(r"effective order $\nu_{\rm eff}$")
    ax.set_ylim(0.5, 2.5)
    ax.legend(loc="lower left", fontsize=6.0, handlelength=1.2,
             title="SC protected, broken by $\\epsilon$", title_fontsize=6.0,
             bbox_to_anchor=(0.0, 0.16))

    inset = ax.inset_axes([0.56, 0.50, 0.40, 0.40])
    eps = np.array([r["eps"] for r in sc_res["rows"]])
    gstar = np.array([r["gamma_star"] for r in sc_res["rows"]])
    inset.loglog(eps, gstar, "o-", color=sty.COLORS["black"], ms=3, linewidth=0.9)
    inset.loglog(eps, gstar[0] * (eps[0] / eps), "--", color="#999999", linewidth=0.7)
    inset.set_xlabel(r"$\epsilon$", fontsize=6.0, labelpad=1)
    inset.set_ylabel(r"$\Gamma_\ast$", fontsize=6.0, labelpad=0)
    inset.tick_params(labelsize=5.0)
    inset.text(0.5, 1.20, fr"$\Gamma_\ast\propto\epsilon^{{{sc_res['crossover_power']:.2f}}}$",
              transform=inset.transAxes, fontsize=6.0, ha="center", va="bottom")


def _panel_b(ax):
    gT = rd.gamma_T_mapping()
    T = gT["T"]
    ax.semilogy(T, gT["nv_korb_GHz"], "o-", color=sty.COLORS["blue"],
              label=fr"NV $k_{{\rm orb}}\propto T^5$ (${gT['nv_decades']:.1f}$ dec)")
    ax.semilogy(T, gT["siv_GHz"], "s-", color=sty.COLORS["vermillion"],
              label=fr"SiV Bose (${gT['siv_decades']:.1f}$ dec)")
    ax.semilogy(T, gT["snv_GHz"], "^-", color=sty.COLORS["green"],
              label=fr"SnV Bose (${gT['snv_decades']:.1f}$ dec)")
    ax.axhspan(1, 10 ** (1 + 1), color=sty.COLORS["sky"], alpha=0.12, zorder=0)
    ax.text(6, 3, "required range\n($\\sim$1$-$2 decades)", fontsize=5.8,
           color="#333333", va="center")
    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$\Gamma(T)$ (GHz)")
    ax.legend(loc="lower right", fontsize=6.2, handlelength=1.3)


def build(quick=False):
    sty.apply_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(sty.COL_2, 2.9))
    _panel_a(ax1, quick)
    _panel_b(ax2)
    sty.panel_label(ax1, "a")
    sty.panel_label(ax2, "b")
    fig.tight_layout(w_pad=2.0)
    return sty.save(fig, "fig4_robustness")


if __name__ == "__main__":
    print(build())
