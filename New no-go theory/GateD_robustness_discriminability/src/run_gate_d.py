"""run_gate_d.py -- Gate D driver: physical significance (PRL Priorities 6 & 8).

Final gate. Two results complete the PRL central claim:

P6  Robustness / crossover fan.
  * EXACT (symmetry-protected) class vs APPROXIMATE (tuned) class.
      - NV ms=-1<->+1 (n=3): M0=M1=0 are structural zeros; strain, transverse
        field and detuning do NOT lift them -> order stays 3, no crossover
        (Gamma* = infinity). An exact class.
      - Superconducting protected (nu=2): a tuned cancellation
        g_A1 g_B1 = -g_A2 g_B2; a small coupling imbalance eps breaks it to
        nu=1 with a crossover Gamma*(eps) ~ 1/eps. An approximate class.
  * Crossover-scale law Gamma*(eps) proportional to 1/eps (numerically).
  * Effective-exponent map nu_eff(Gamma, eps); the high-order class is
    observable for Gamma < Gamma*(eps), a window that opens to all Gamma as
    eps -> 0.
  * Intervention-scaling nu(q) fan (Phase N) as a general corroboration of the
    crossover as a Newton-polygon phenomenon.

P8  Experimental discriminability.
  * Required Gamma dynamic range to resolve an adjacent-class exponent
    difference (Delta nu = 1) to a target precision.
  * Platform reach: superconducting kappa is directly tunable over decades;
    NV phonon rate Gamma(T)=k_orb ~ T^5 spans ~8.5 decades over 4-300 K
    (temperature-tuned exponent measurement); group-IV Bose-law Gamma(T) is
    narrow (SiV ~1.6, SnV ~0.6 decade).
  * Optical read-out feasibility via signal_chain.py (ZPL OD -> SNR): the
    minimum detectable contrast at a feasible density x integration time, so
    the log-log slope (the class) is measurable.

Usage:  python run_gate_d.py [--quick] [--smoke]
Outputs: results/tables/gates_summary_gateD.json, gate_d_robustness.csv,
         gate_d_discriminability.csv, results/figures/*.png (+ .pdf)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE_SRC = os.path.join(HERE, "..", "..", "src")
NOGO_SRC = os.path.join(HERE, "..", "..", "..", "No-go theorem", "src")
GATEB_SRC = os.path.join(HERE, "..", "..", "GateB_superconducting_witness", "src")
PHASEN12_SRC = os.path.join(HERE, "..", "..", "PhaseN", "priority_1_2")
for _p in (HERE, PHASE_SRC, NOGO_SRC, GATEB_SRC, PHASEN12_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import core                              # noqa: E402
import nv_reduced_kernel as nvk          # noqa: E402
import model_sc_transfer as sc           # noqa: E402
import phonon_rates as pr                # noqa: E402
import group_iv_model as giv             # noqa: E402
import nv_model as nv                    # noqa: E402
import signal_chain as sig               # noqa: E402
import phase_n_exact_core as pn          # noqa: E402

RESULTS = os.path.join(HERE, "..", "results")
TABLES = os.path.join(RESULTS, "tables")
FIGS = os.path.join(RESULTS, "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

_NV_H = nvk.H_3E()


def _slope(g, v):
    return core.fit_nu_loglog(np.asarray(g, float), np.asarray(v))["nu_global"]


# =====================================================================
# P6 -- robustness / crossover
# =====================================================================
def nv_exact_class():
    """NV ms=-1<->+1 (n=3): perturbations do NOT lift M0=M1=0 -> exact class."""
    ks = np.logspace(2, 6, 40)
    perturbations = [("unperturbed", {}),
                     ("xi_x=0.05", dict(xi_x=0.05)), ("xi_x=0.2", dict(xi_x=0.2)),
                     ("xi_y=0.1", dict(xi_y=0.1)), ("Bx=0.2", dict(Bx=0.2))]
    rows = []
    z_rows = []
    for label, kw in perturbations:
        H = nvk.H_3E(**kw)
        M = nvk.moments(H, (-1, 1), 3)
        sl = _slope(ks, nvk.kernel(H, (-1, 1), ks))
        rows.append(dict(perturbation=label, abs_M0=float(abs(M[0])),
                         abs_M1=float(abs(M[1])), slope=float(sl),
                         order_unchanged=bool(abs(sl - 3) < 0.02)))
    # also detuning z (probe frequency) leaves the structural zeros intact
    for z in (0.0, 0.5, 1.0):
        M = nvk.moments(_NV_H, (-1, 1), 3, z=z)
        z_rows.append(dict(z=z, abs_M0=float(abs(M[0])), abs_M1=float(abs(M[1]))))
    exact = all(r["order_unchanged"] and r["abs_M0"] < 1e-9 and r["abs_M1"] < 1e-9
                for r in rows) and all(r["abs_M0"] < 1e-9 and r["abs_M1"] < 1e-9 for r in z_rows)
    return dict(rows=rows, z_rows=z_rows, is_exact_class=bool(exact))


def sc_approximate_class(quick=False):
    """Superconducting protected nu=2: eps breaks it to nu=1 with Gamma*(eps)."""
    ks = np.logspace(4, 13, 90 if quick else 140)
    epsilons = [1e-8, 1e-9, 1e-10, 1e-11]
    rows = []
    curves = {}
    for eps in epsilons:
        K = np.array([sc.transfer_kernel(k, tuning="broken", eps=eps) for k in ks])
        fit = core.fit_nu_loglog(ks, K)
        ne, gm = fit["nu_eff"], fit["gamma_mid"]
        below = np.where(ne < 1.5)[0]
        kstar = float(gm[below[0]]) if len(below) else float("inf")
        rows.append(dict(eps=float(eps), gamma_star=kstar))
        curves[eps] = (gm, ne)
    # protected (eps=0) stays nu=2
    K0 = np.array([sc.transfer_kernel(k, tuning="protected") for k in ks])
    nu0 = float(core.fit_nu_loglog(ks, K0)["nu_eff"][-1])
    # crossover law Gamma*(eps) ~ eps^power
    epsv = np.array([r["eps"] for r in rows])
    gstar = np.array([r["gamma_star"] for r in rows])
    power = float(np.polyfit(np.log10(epsv), np.log10(gstar), 1)[0])  # expect -1
    return dict(rows=rows, curves=curves, nu_protected_eps0=nu0,
                crossover_power=power, law_is_inverse=bool(abs(power + 1) < 0.1))


def nu_q_fan():
    """Phase N intervention-scaling fan nu(q) = 4-q, 2+q, 4 (Newton polygon)."""
    num, den = pn.master_polynomials()
    def expected(q):
        q = Fraction(q)
        return 4 - q if q <= 1 else (2 + q if q <= 2 else Fraction(4))
    rows, ok = [], True
    for q in [Fraction(0), Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3)]:
        order, _cn, _cd = pn.path_order(q, num, den)
        m = (order == expected(q))
        ok = ok and m
        rows.append(dict(q=str(q), nu=str(order), expected=str(expected(q)), match=bool(m)))
    return dict(rows=rows, fan_ok=bool(ok))


# =====================================================================
# P8 -- experimental discriminability
# =====================================================================
def required_gamma_range():
    """Decades of clean asymptotic Gamma needed to pin a power-law exponent to
    +-0.1 (so an adjacent class Delta nu = 1 is unambiguously resolved). Uses
    the clean NV Gamma^-2 kernel; fit the slope over windows of increasing
    width and report the smallest width reaching |slope - 2| < 0.1."""
    decades_needed = None
    per = []
    for W in (1, 2, 3, 4):
        ks = np.logspace(3, 3 + W, 12 * W)
        err = abs(_slope(ks, nvk.kernel(_NV_H, (0, -1), ks)) - 2.0)
        per.append(dict(window_decades=W, slope_err=float(err)))
        if decades_needed is None and err < 0.1:
            decades_needed = W
    return dict(per_window=per, decades_needed=decades_needed or 99,
                # Delta nu = 1 separation needs ~ this many decades of clean Gamma
                delta_nu_resolvable=bool(decades_needed is not None))


def gamma_T_mapping():
    """Gamma(T) reach per platform vs the required dynamic range."""
    Ts = [4, 10, 20, 50, 100, 200, 300]
    d = 1.683
    nv_g = [float(nv.korb_GHz(T, d)) for T in Ts]
    siv_g = [float(giv.gamma_orb_GHz("SiV", T)) for T in Ts]
    snv_g = [float(giv.gamma_orb_GHz("SnV", T)) for T in Ts]
    def decades(vals):
        v = [x for x in vals if x > 0]
        return float(np.log10(max(v) / min(v)))
    return dict(
        T=Ts, nv_korb_GHz=nv_g, siv_GHz=siv_g, snv_GHz=snv_g,
        nv_decades=decades(nv_g), siv_decades=decades(siv_g), snv_decades=decades(snv_g),
        sc_decades_available=9.0,  # kappa engineerable over >= 9 decades (Gate B sweep)
        regime_SiV_300K=float(giv.thermal_regime("SiV", 300.0)),
        regime_SnV_300K=float(giv.thermal_regime("SnV", 300.0)),
    )


def optical_snr():
    """Minimum detectable ZPL contrast at a feasible operating point (1 ppm,
    1 h), so the response slope (the class) is measurable via optical OD."""
    lam, n_refr, DW, gamma_inh = 637.0, 2.41, 0.035, 30.0
    eta, power, sigma_tech, target = 0.1, 1e-6, 1e-6, 5.0
    n_nv, L, tau = 1.76e17, 0.05, 3600.0  # 1 ppm, 0.5 mm, 1 hour
    T = 50.0
    gamma_h = float(nv.gamma_oc_GHz(T, 1.683))
    sigma = sig.sigma_zpl_cm2(lam, n_refr, DW, nv.GRAD, gamma_h)
    f_spec = sig.spectral_fraction(gamma_h, gamma_inh)
    alpha = sig.alpha_cm(sigma, n_nv, 0.25, 1 / 3, f_spec)
    od_sector = sig.od(alpha, L)
    od_total = od_sector
    c_min = sig.min_detectable_contrast(target, od_sector, od_total, power, lam, tau,
                                        eta, sigma_tech)
    # a modest sector-EIT contrast (order 1e-2) is well above c_min -> detectable
    representative_contrast = 1e-2
    detectable = bool(np.isfinite(c_min) and representative_contrast > c_min)
    return dict(gamma_h_GHz=gamma_h, sigma_zpl_cm2=sigma, spectral_fraction=f_spec,
                od_sector=od_sector, min_detectable_contrast=float(c_min),
                representative_contrast=representative_contrast, detectable=detectable,
                operating_point=dict(n_nv_cm3=n_nv, L_cm=L, tau_s=tau, T_K=T,
                                     density_ppm=1.0))


# =====================================================================
# figures
# =====================================================================
def make_figures(sc_res, gT):
    # robustness fan: nu_eff(Gamma) for several eps + Gamma*(eps)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.6))
    for eps, (gm, ne) in sc_res["curves"].items():
        ax1.semilogx(gm, ne, "-", lw=1.0, label=fr"$\epsilon={eps:.0e}$")
    ax1.axhline(2, color="gray", ls="--", lw=0.7)
    ax1.axhline(1, color="gray", ls=":", lw=0.7)
    ax1.set_xlabel(r"$\Gamma$"); ax1.set_ylabel(r"$\nu_{\rm eff}$")
    ax1.set_title(r"Approximate class: $\nu:2\to1$ at $\Gamma_\ast(\epsilon)$")
    ax1.set_ylim(0.5, 2.5); ax1.legend(fontsize=7)
    eps = np.array([r["eps"] for r in sc_res["rows"]])
    gstar = np.array([r["gamma_star"] for r in sc_res["rows"]])
    ax2.loglog(eps, gstar, "o-", label="numerical")
    ax2.loglog(eps, gstar[0] * (eps[0] / eps), "k--", lw=0.8, label=r"$\propto 1/\epsilon$")
    ax2.set_xlabel(r"$\epsilon$ (symmetry breaking)"); ax2.set_ylabel(r"$\Gamma_\ast$")
    ax2.set_title(r"Crossover scale $\Gamma_\ast(\epsilon)\propto 1/\epsilon$")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(FIGS, "fig_gateD_robustness.png")
    fig.savefig(p, dpi=200); fig.savefig(p.replace(".png", ".pdf")); plt.close(fig)

    # Gamma(T) mapping
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    T = gT["T"]
    ax.semilogy(T, gT["nv_korb_GHz"], "o-", label=r"NV $k_{\rm orb}\propto T^5$")
    ax.semilogy(T, gT["siv_GHz"], "s-", label="SiV Bose")
    ax.semilogy(T, gT["snv_GHz"], "^-", label="SnV Bose")
    ax.set_xlabel("T (K)"); ax.set_ylabel(r"$\Gamma(T)$ (GHz)")
    ax.set_title(r"Finite-$T$ rate mapping (NV spans ~8.5 decades, group-IV narrow)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(FIGS, "fig_gateD_gammaT.png")
    fig.savefig(p, dpi=200); fig.savefig(p.replace(".png", ".pdf")); plt.close(fig)


# =====================================================================
# main
# =====================================================================
def main():
    ap = argparse.ArgumentParser(description="Gate D: robustness + discriminability")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    quick = args.quick or args.smoke

    t0 = time.time()
    nv_exact = nv_exact_class()
    sc_res = sc_approximate_class(quick=quick)
    fan = nu_q_fan()
    req = required_gamma_range()
    gT = gamma_T_mapping()
    snr = optical_snr()

    # -- gates --------------------------------------------------------
    g_d1 = bool(nv_exact["is_exact_class"] and sc_res["nu_protected_eps0"] > 1.9
                and all(np.isfinite(r["gamma_star"]) for r in sc_res["rows"]))
    g_d2 = bool(sc_res["law_is_inverse"])
    g_d3 = bool(fan["fan_ok"]
                and all(sc_res["rows"][i]["gamma_star"] < sc_res["rows"][i + 1]["gamma_star"]
                        for i in range(len(sc_res["rows"]) - 1)))
    g_d4 = bool(req["delta_nu_resolvable"] and req["decades_needed"] <= 4)
    g_d5 = bool(gT["nv_decades"] >= req["decades_needed"]
                and gT["sc_decades_available"] >= req["decades_needed"]
                and gT["siv_decades"] < gT["nv_decades"])
    g_d6 = bool(snr["detectable"] and np.isfinite(snr["min_detectable_contrast"]))

    gates = dict(
        G_D1_exact_vs_approximate_class=g_d1,
        G_D2_crossover_scale_law=g_d2,
        G_D3_exponent_map_and_fan=g_d3,
        G_D4_required_gamma_range=g_d4,
        G_D5_gammaT_platform_reach=g_d5,
        G_D6_optical_snr_discriminability=g_d6,
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        description="Gate D: physical significance -- robustness/crossover (P6) and "
                    "experimental discriminability (P8). Completes the PRL central claim.",
        p6_nv_exact_class=nv_exact,
        p6_sc_approximate_class={k: v for k, v in sc_res.items() if k != "curves"},
        p6_nu_q_fan=fan,
        p8_required_gamma_range=req,
        p8_gamma_T_mapping=gT,
        p8_optical_snr=snr,
        quick=quick, gates=gates, runtime_s=round(time.time() - t0, 2),
    )

    out_json = os.path.join(TABLES, "gates_summary_gateD.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    with open(os.path.join(TABLES, "gate_d_robustness.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["eps", "gamma_star"])
        for r in sc_res["rows"]:
            w.writerow([r["eps"], r["gamma_star"]])
    with open(os.path.join(TABLES, "gate_d_discriminability.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["T_K", "NV_korb_GHz", "SiV_GHz", "SnV_GHz"])
        for i, T in enumerate(gT["T"]):
            w.writerow([T, gT["nv_korb_GHz"][i], gT["siv_GHz"][i], gT["snv_GHz"][i]])

    make_figures(sc_res, gT)

    print(json.dumps(gates, indent=2))
    print(f"\nP6 exact class (NV -1<->+1): {nv_exact['is_exact_class']} "
          f"(M0=M1=0 under all perturbations, order stays 3)")
    print(f"P6 approximate class (SC): Gamma*(eps) power={sc_res['crossover_power']:.3f} "
          f"(expect -1); protected eps=0 nu={sc_res['nu_protected_eps0']:.2f}")
    print(f"P6 nu(q) fan ok: {fan['fan_ok']}")
    print(f"P8 required Gamma range: {req['decades_needed']} decades for +-0.1 slope")
    print(f"P8 Gamma(T) decades: NV={gT['nv_decades']:.1f} SiV={gT['siv_decades']:.1f} "
          f"SnV={gT['snv_decades']:.1f} SC>={gT['sc_decades_available']:.0f}")
    print(f"P8 optical: C_min={snr['min_detectable_contrast']:.2e} detectable={snr['detectable']}")
    print(f"\nwrote {out_json}")
    print(f"runtime {summary['runtime_s']} s")
    return summary


if __name__ == "__main__":
    main()
