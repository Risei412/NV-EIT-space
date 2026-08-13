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
    minimum detectable contrast at a feasible density x integration time,
    compared against the ENSEMBLE-AVERAGED contrasts measured by the PRA
    campaign's Gate 5 rather than an assumed single-defect value, and at the
    same 70 K candidate temperature as those contrasts.
  * Slope budget (G-D7): whether the exponent is actually measurable, i.e.
    whether the Gamma window that stays above the detection floor is as wide
    as the noisy log-log fit needs. Optical NV fails this once ensemble
    washout is included -- at nu = 4 one decade of Gamma costs four decades of
    signal -- while engineered-dissipation platforms pass comfortably. The
    law is platform-independent, so one realization that resolves the class
    satisfies P8; which platforms do and do not is reported explicitly.

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
GATEC_SRC = os.path.join(HERE, "..", "..", "GateC_material_independence", "src")
# phase_n_exact_core lives in New no-go theory/src (PHASE_SRC). The theory-side
# PhaseN tree moved to the SMRT repository, so there is no PhaseN/priority_1_2
# path here to add.
for _p in (HERE, PHASE_SRC, NOGO_SRC, GATEB_SRC, GATEC_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import core                              # noqa: E402
import experimental_budget as eb         # noqa: E402
import chain3_witness as chain3          # noqa: E402
import nv_reduced_kernel as nvk          # noqa: E402
import model_sc_transfer as sc           # noqa: E402
import phonon_rates as pr                # noqa: E402
import group_iv_model as giv             # noqa: E402
import nv_model as nv                    # noqa: E402
import signal_chain as sig               # noqa: E402
import phase_n_exact_core as pn          # noqa: E402

RESULTS = os.path.join(HERE, "..", "..", "..", "..", "..", "results")
TABLES = os.path.join(RESULTS, "tables")
FIGS = os.path.join(RESULTS, "figures")
CERTS = os.path.join(RESULTS, "certificates")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(CERTS, exist_ok=True)
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
    """Decades of Gamma needed to resolve an adjacent class AT FINITE SNR.

    This used to fit the noiseless analytic kernel and ask for |slope - 2| <
    0.1, which one decade satisfies to 2e-5 -- a statement about floating-point
    arithmetic, not about an experiment. The window width is now decided by the
    spread of the fitted slope under measurement noise, requiring
    |bias| + 2*std < 0.5 so two adjacent integer classes separate at 2 sigma.

    Reported across a range of per-point precisions: sigma_rel = 0.20 is the
    SNR = 5 detection threshold, 0.01 a well-resolved measurement.
    """
    def nv_kernel(gs):
        return nvk.kernel(_NV_H, (0, -1), gs)

    per_sigma = {}
    for sigma_rel in (0.01, 0.03, 0.10, 0.20):
        r = eb.decades_needed(nv_kernel, 1e3, sigma_rel, seed=20260801)
        per_sigma[str(sigma_rel)] = r
    at_threshold = per_sigma["0.2"]["decades_needed"]
    return dict(
        criterion="|bias| + 2*std < %.1f on the fitted log-log slope" % eb.SLOPE_RESOLUTION,
        per_sigma_rel=per_sigma,
        decades_needed=at_threshold if at_threshold is not None else 99,
        decades_needed_well_resolved=per_sigma["0.01"]["decades_needed"],
        delta_nu_resolvable=bool(at_threshold is not None),
    )


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
        # Two different things, previously conflated into a hard-coded 9.0.
        # The sweep range is how far the numerics ran; the engineerable range is
        # how far a bus can actually be tuned. Only the latter is experimental
        # reach, and it is nearly six decades smaller.
        sc_decades_swept=sc.decades(sc.KAPPA_ASYMPTOTIC),
        sc_decades_engineerable=sc.decades(sc.KAPPA_ENGINEERABLE),
        sc_kappa_engineerable_GHz=list(sc.KAPPA_ENGINEERABLE),
        regime_SiV_300K=float(giv.thermal_regime("SiV", 300.0)),
        regime_SnV_300K=float(giv.thermal_regime("SnV", 300.0)),
    )


# The PRA campaign's candidate point is 70 K; c_min must be evaluated there,
# not at some other temperature, or the comparison mixes two operating points.
CANDIDATE_T_K = 70.0
# Gate 5's own pass criterion is that post-selection (one orientation class,
# spectral-hole selection) plus field shimming keeps the feature above 1e-3, so
# that scenario is the stated recommended operating point and the one to gate
# on. Everything else is reported alongside it, not folded into a verdict.
GATING_SCENARIO = "post_selected_shimmed"


def _c_min_at(T, tau=3600.0, n_nv=1.76e17, L=0.05):
    lam, n_refr, DW, gamma_inh = 637.0, 2.41, 0.035, 30.0
    eta, power, sigma_tech, target = 0.1, 1e-6, 1e-6, 5.0
    gamma_h = float(nv.gamma_oc_GHz(T, 1.683))
    sigma = sig.sigma_zpl_cm2(lam, n_refr, DW, nv.GRAD, gamma_h)
    f_spec = sig.spectral_fraction(gamma_h, gamma_inh)
    alpha = sig.alpha_cm(sigma, n_nv, 0.25, 1 / 3, f_spec)
    od_sector = sig.od(alpha, L)
    c_min = sig.min_detectable_contrast(target, od_sector, od_sector, power, lam, tau,
                                        eta, sigma_tech)
    return dict(T_K=T, gamma_h_GHz=gamma_h, sigma_zpl_cm2=sigma,
                spectral_fraction=f_spec, od_sector=od_sector,
                min_detectable_contrast=float(c_min), target_snr=target,
                operating_point=dict(n_nv_cm3=n_nv, L_cm=L, tau_s=tau,
                                     density_ppm=1.0))


def optical_snr():
    """Detectability of the NV EIT contrast against the optical floor.

    The representative contrast used to be the literal 1e-2. That is
    essentially the SINGLE-DEFECT value from the PRA campaign's ensemble study
    (0.0136) -- i.e. the value before any ensemble averaging. Every averaged
    scenario in that same table is one to two orders of magnitude smaller and
    the worst of them sits level with the detection floor, so the hard-coded
    number was quietly assuming away the washout the study had measured.

    All five scenarios are now read from that table and reported; the verdict
    is taken at the recommended operating point (GATING_SCENARIO), and c_min is
    evaluated at the same 70 K candidate temperature as the contrasts rather
    than at 50 K as before.
    """
    scen = eb.load_gate5_contrast()
    floor = _c_min_at(CANDIDATE_T_K)
    floor_50K = _c_min_at(50.0)
    c_min = floor["min_detectable_contrast"]

    per_scenario = {}
    for name, row in scen.items():
        c = row["Cmax"]
        per_scenario[name] = dict(
            contrast=c,
            washout_factor=row["washout_factor"],
            margin_over_floor=float(c / c_min) if c_min > 0 else float("inf"),
            sigma_rel_per_point=eb.sigma_rel_from_contrast(c, c_min),
            detectable=bool(np.isfinite(c_min) and c > c_min),
        )

    gated = per_scenario[GATING_SCENARIO]
    return dict(
        provenance=dict(
            contrast_source=os.path.relpath(eb.GATE5_CSV, HERE),
            note=("contrasts are ensemble-averaged results of the PRA campaign's "
                  "Gate 5; the previous hard-coded 1e-2 corresponded to the "
                  "un-averaged single-defect scenario"),
            gating_scenario=GATING_SCENARIO,
            gating_rationale=("Gate 5's own pass criterion is that the "
                              "post-selected, field-shimmed configuration stays "
                              "above 1e-3"),
        ),
        floor_at_candidate_T=floor,
        floor_at_50K_previous_convention=floor_50K,
        min_detectable_contrast=c_min,
        per_scenario=per_scenario,
        representative_contrast=gated["contrast"],
        detectable=gated["detectable"],
        worst_case_scenario="high_density",
        worst_case_detectable=per_scenario["high_density"]["detectable"],
    )


def slope_budget(snr):
    """Can each platform actually MEASURE its exponent?

    Combines the Gamma window that stays above the detection floor with the
    slope precision achievable there. A platform resolves its class only if the
    usable window is at least as wide as the width the noisy fit needs.

    The optical NV entry is the hard case by construction: its observable order
    is 4, so one decade of Gamma costs four decades of signal and the window
    above the floor closes almost immediately. Platforms where the dissipation
    is an engineered knob read out by a transmission amplitude do not pay that.
    """
    c_min = snr["min_detectable_contrast"]
    rows = []

    for scenario in ("single", GATING_SCENARIO, "high_density"):
        s = snr["per_scenario"][scenario]
        win = eb.usable_gamma_window(nu=4.0, contrast_ref=s["contrast"],
                                     gamma_ref=1e3, c_min=c_min)
        sigma_rel = min(s["sigma_rel_per_point"], 0.2)
        need = eb.decades_needed(
            lambda gs: nvk.kernel(_NV_H, (-1, 1), gs), 1e3, sigma_rel,
            seed=20260801)["decades_needed"]
        rows.append(dict(
            platform=f"NV optical EIT ({scenario})", nu=4.0,
            readout="ZPL absorption contrast vs Gamma(T)",
            window_decades=win["decades"], sigma_rel=sigma_rel,
            decades_needed=need,
            class_resolvable=bool(need is not None and win["decades"] >= need),
        ))

    # Engineered-dissipation platform: kappa is a design knob swept directly and
    # the readout is a transmission amplitude, so the precision does not decay
    # with kappa the way an optical contrast does.
    sc_window = sc.decades(sc.KAPPA_ENGINEERABLE)
    for label, nu, tuning in (("SC transfer (generic)", 1.0, "generic"),
                              ("SC transfer (protected)", 2.0, "protected")):
        need = eb.decades_needed(
            lambda gs, _t=tuning: np.array([sc.transfer_kernel(g, tuning=_t) for g in gs]),
            1e6, 0.01, seed=20260801)["decades_needed"]
        rows.append(dict(
            platform=label, nu=nu, readout="bus transmission amplitude vs kappa",
            window_decades=sc_window, sigma_rel=0.01, decades_needed=need,
            class_resolvable=bool(need is not None and sc_window >= need),
        ))

    # Non-diamond class-3 chain, the Gate C witness, with the same 1% readout.
    chain_need = eb.decades_needed(
        lambda gs: np.array([chain3.kernel(g) for g in gs]), 1e2, 0.01, seed=20260801
    )["decades_needed"]
    rows.append(dict(
        platform="3-mode chain (class 3)", nu=3.0,
        readout="marked-port transmission amplitude vs engineered loss",
        window_decades=3.0, sigma_rel=0.01, decades_needed=chain_need,
        class_resolvable=bool(chain_need is not None and 3.0 >= chain_need),
    ))

    return dict(rows=rows,
                any_platform_resolvable=bool(any(r["class_resolvable"] for r in rows)),
                optical_nv_resolvable=bool(
                    any(r["class_resolvable"] for r in rows
                        if r["platform"].startswith("NV optical"))))


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
    budget = slope_budget(snr)

    # -- gates --------------------------------------------------------
    g_d1 = bool(nv_exact["is_exact_class"] and sc_res["nu_protected_eps0"] > 1.9
                and all(np.isfinite(r["gamma_star"]) for r in sc_res["rows"]))
    g_d2 = bool(sc_res["law_is_inverse"])
    g_d3 = bool(fan["fan_ok"]
                and all(sc_res["rows"][i]["gamma_star"] < sc_res["rows"][i + 1]["gamma_star"]
                        for i in range(len(sc_res["rows"]) - 1)))
    g_d4 = bool(req["delta_nu_resolvable"] and req["decades_needed"] <= 4)
    # Reach is measured against the ENGINEERABLE kappa range, not the numerical
    # sweep width that used to be hard-coded as 9.0 decades.
    g_d5 = bool(gT["nv_decades"] >= req["decades_needed"]
                and gT["sc_decades_engineerable"] >= req["decades_needed"]
                and gT["siv_decades"] < gT["nv_decades"])
    g_d6 = bool(snr["detectable"] and np.isfinite(snr["min_detectable_contrast"]))
    # G-D7: the exponent must be measurable somewhere, combining the Gamma
    # window that stays above the detection floor with the slope precision
    # available in it. "At least one platform" is deliberate: the law is
    # platform-independent, so one realization where the class is readable is
    # what P8 asks for. Which platforms do and do not qualify is reported.
    g_d7 = bool(budget["any_platform_resolvable"])

    gates = dict(
        G_D1_exact_vs_approximate_class=g_d1,
        G_D2_crossover_scale_law=g_d2,
        G_D3_exponent_map_and_fan=g_d3,
        G_D4_required_gamma_range=g_d4,
        G_D5_gammaT_platform_reach=g_d5,
        G_D6_optical_snr_discriminability=g_d6,
        G_D7_slope_measurable=g_d7,
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
        p8_slope_budget=budget,
        quick=quick, gates=gates, runtime_s=round(time.time() - t0, 2),
    )

    out_json = os.path.join(CERTS, "gates_summary_gateD.json")
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
    print(f"P8 required Gamma range: {req['decades_needed']} decades at sigma_rel=0.20 "
          f"({req['criterion']})")
    print(f"P8 Gamma(T) decades: NV={gT['nv_decades']:.1f} SiV={gT['siv_decades']:.1f} "
          f"SnV={gT['snv_decades']:.1f} SC engineerable={gT['sc_decades_engineerable']:.1f} "
          f"(swept {gT['sc_decades_swept']:.0f})")
    print(f"P8 optical: C_min={snr['min_detectable_contrast']:.2e} "
          f"(T={snr['floor_at_candidate_T']['T_K']:.0f} K) "
          f"gating={snr['provenance']['gating_scenario']} "
          f"C={snr['representative_contrast']:.3e} detectable={snr['detectable']}")
    print("P8 slope budget:")
    for r in budget["rows"]:
        print(f"    {r['platform']:<38} nu={r['nu']:.0f} window={r['window_decades']:.2f} dec "
              f"need={r['decades_needed']} resolvable={r['class_resolvable']}")
    print(f"\nwrote {out_json}")
    print(f"runtime {summary['runtime_s']} s")
    return summary


if __name__ == "__main__":
    main()
