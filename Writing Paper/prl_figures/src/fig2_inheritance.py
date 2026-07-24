"""fig2_inheritance.py -- PRL Fig. 2: observable-order inheritance + crossover.

(a) nu_obs = n12 + n21 - nu_den on the synthetic three-class family
    (generic nu_obs=2, protected 4, doubly-protected 6; Gate A data).
(b) NV EIT pre-asymptotic vs asymptotic crossover: nu_eff(Gamma), Gamma_cross,
    and the physical operating point Gamma(300K) marked (Gate A separate_regimes).

No new physics: reuses gate_a_observable.py and model_specs.py (Gate A) unchanged.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
GATEA_SRC = os.path.join(HERE, "..", "..", "..", "New no-go theory",
                         "PhaseO_observable_inheritance", "src")
for _p in (HERE, GATEA_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import prl_style as sty  # noqa: E402
import gate_a_observable as gao  # noqa: E402
import model_specs as ms  # noqa: E402


def _panel_a(ax):
    """nu_obs = n12+n21-nu_den on the synthetic 3-class family."""
    gammas = np.logspace(2, 7, 60)
    specs = ms.synthetic_specs()
    colors = [sty.COLORS["blue"], sty.COLORS["vermillion"], sty.COLORS["green"]]
    markers = ["o", "s", "^"]
    labels = [r"generic: $\nu_{\rm obs}=2$", r"protected: $\nu_{\rm obs}=4$",
             r"doubly-protected: $\nu_{\rm obs}=6$"]
    for spec, color, marker, label in zip(specs, colors, markers, labels):
        v = gao.verify_nu_obs_loglog(spec, gammas, z=0.0)
        y = np.abs(v["vals"])
        ax.loglog(gammas, y, marker, ms=2.5, color=color, mfc="none", mew=0.9,
                  label=f"{label} (fit ${v['nu_tail']:.2f}$)")
    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"$|R_{\rm obs}(\Gamma)|$")
    ax.legend(loc="lower left", fontsize=6.2, handlelength=1.2)


def _panel_b(ax):
    """NV pre-asymptotic vs asymptotic crossover."""
    nv = ms.nv_spec()
    gammas = np.logspace(0, 10, 140)
    v = gao.verify_nu_obs_loglog(nv, gammas, z=0.0)
    ne, gm = v["nu_eff"], v["gamma_mid"]

    reg = gao.separate_regimes(nv, gammas, z=0.0,
                               gamma_phys=nv.meta["gamma_phys_300K"],
                               generic_order=3.0, asymptotic_order=4.0)

    ax.semilogx(gm, ne, "-", color=sty.COLORS["blue"], linewidth=1.3)
    ax.axhline(4, color=sty.COLORS["green"], ls="--", linewidth=0.9,
              label=r"asymptotic $\nu=4$")
    ax.axhline(3, color=sty.COLORS["vermillion"], ls=":", linewidth=0.9,
              label=r"pre-asymptotic $\nu\approx3$")
    ax.axvline(reg["gamma_cross"], color="#555555", ls="-.", linewidth=0.9,
              label=fr"$\Gamma_{{\rm cross}}\approx{reg['gamma_cross']:.1e}$")
    ax.axvline(nv.meta["gamma_phys_300K"], color=sty.COLORS["purple"], ls="-",
              linewidth=1.1, alpha=0.8,
              label=fr"$\Gamma(300\,{{\rm K}})\approx{nv.meta['gamma_phys_300K']:.1e}$")

    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"effective order $\nu_{\rm eff}$")
    ax.set_ylim(0, 4.6)
    ax.legend(loc="lower right", fontsize=6.0, handlelength=1.2)


def build(quick=False):
    sty.apply_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(sty.COL_2, 2.7))
    _panel_a(ax1)
    _panel_b(ax2)
    sty.panel_label(ax1, "a")
    sty.panel_label(ax2, "b")
    fig.tight_layout(w_pad=2.0)
    return sty.save(fig, "fig2_inheritance")


if __name__ == "__main__":
    print(build())
