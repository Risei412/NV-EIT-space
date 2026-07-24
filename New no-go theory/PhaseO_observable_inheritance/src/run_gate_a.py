"""run_gate_a.py -- Gate A driver: observable-order inheritance (Priority 1).

Runs the model-independent inheritance verification across every observable
class and emits a PASS/FAIL gate summary.

Gates (strategy doc Sec.13 P1 pass conditions):
  G-O1  inheritance law nu_obs = n12 + n21 - nu_den derived & confirmed for
        >= 3 distinct observables (bilinear susceptibility, quadratic QFI,
        frozen-source Raman/difference).
  G-O2  generic vs symmetry-protected cancellation distinguished by the
        EXACT rational-degree certificate (not a float slope).
  G-O3  predicted nu_obs matches the full log-log (deep-tail) slope.
  G-O4  pre-asymptotic effective exponent separated from the true asymptotic
        index, with the crossover scale Gamma_cross located.
  G-O5  predictor / exact certificate / log-log fit agree across all models,
        and the Phase N q-fan and z* resonance corroborate the modifiers.

Usage:
    python run_gate_a.py [--quick] [--smoke]

Outputs:
    results/tables/gates_summary_gateO.json
    results/tables/gate_a_models.csv
    results/figures/fig_gateO_inheritance.png (+ .pdf)
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
ROOMT_SRC = os.path.join(HERE, "..", "..", "RoomT", "src")
PHASEN12_SRC = os.path.join(HERE, "..", "..", "PhaseN", "priority_1_2")
PHASEN3_SRC = os.path.join(HERE, "..", "..", "PhaseN", "priority_3_frequency")
for _p in (HERE, PHASE_SRC, NOGO_SRC, ROOMT_SRC, PHASEN12_SRC, PHASEN3_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import gate_a_observable as gao          # noqa: E402
import model_specs as ms                 # noqa: E402
import core                              # noqa: E402

RESULTS = os.path.join(HERE, "..", "results")
TABLES = os.path.join(RESULTS, "tables")
FIGS = os.path.join(RESULTS, "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def _clean(obj):
    """Recursively make a structure JSON-serializable."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.complexfloating, complex)):
        return {"re": float(np.real(obj)), "im": float(np.imag(obj))}
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    if isinstance(obj, Fraction):
        return str(obj)
    return obj


# ----------------------------------------------------------------------
# 1. bilinear inheritance: synthetic three-class family + physical NV
# ----------------------------------------------------------------------
def run_bilinear(quick=False):
    gammas_syn = np.logspace(2, 7, 40 if quick else 90)
    gammas_nv = np.logspace(0, 10, 60 if quick else 140)

    rows = []
    curves = {}

    for spec in ms.synthetic_specs():
        chk = gao.check_model(spec, gammas_syn, z=0.0, kmax=8)
        rows.append(chk)
        v = gao.verify_nu_obs_loglog(spec, gammas_syn, z=0.0)
        curves[spec.label] = (v["gammas"], np.abs(v["vals"]))

    nv = ms.nv_spec()
    chk_nv = gao.check_model(nv, gammas_nv, z=0.0, kmax=8, slope_tol=0.06)
    rows.append(chk_nv)
    v_nv = gao.verify_nu_obs_loglog(nv, gammas_nv, z=0.0)
    curves[nv.label] = (v_nv["gammas"], np.abs(v_nv["vals"]))

    # pre-asymptotic vs asymptotic separation on the physical NV observable
    reg_nv = gao.separate_regimes(
        nv, gammas_nv, z=0.0, gamma_phys=nv.meta["gamma_phys_300K"],
        generic_order=3.0, asymptotic_order=4.0)
    reg_nv["gamma_phys_300K"] = nv.meta["gamma_phys_300K"]
    reg_nv["gamma_phys_10K"] = nv.meta["gamma_phys_10K"]

    return rows, curves, reg_nv


# ----------------------------------------------------------------------
# 2. NV exact certificate via the already-validated RoomT step 3
# ----------------------------------------------------------------------
def run_nv_exact():
    import step3_merged_manifold_moments as step3
    cert = step3.build_symbolic_certificate()
    analysis = step3.analyze_certificate(cert)
    return dict(
        M0_exact_zero=bool(analysis["M0_is_exact_zero"]),
        nu_K12=analysis["nu_K12"], nu_K21=analysis["nu_K21"],
        nu_S2=analysis["nu_S2"], nu_R=analysis["nu_R"],
        inheritance_ok=bool(analysis["nu_R"] == analysis["nu_K12"] + analysis["nu_K21"]),
        smrt_class=analysis["smrt_class"],
    )


# ----------------------------------------------------------------------
# 3. quadratic (QFI) observable: nu_obs = 2*nu  (Phase M model)
# ----------------------------------------------------------------------
def run_qfi(quick=False):
    import model_metro_lindblad as mm
    gammas = np.logspace(0.5, 4.0, 24 if quick else 48)
    theta, kappa0, lam, phi = 0.0, 1.5, 0.0, 0.0

    xnorm, fqs = [], []
    for g in gammas:
        x_S, _df, _dc, rho_full, _rc = mm.x_S_lindblad(g, theta, lam, kappa0=kappa0, phi=phi)
        xnorm.append(np.linalg.norm(x_S))
        F, _drop, _pmin = mm.qfi(rho_full, mm.unvec(x_S))
        fqs.append(abs(F))
    xnorm = np.array(xnorm)
    fqs = np.array(fqs)

    ntail = max(4, len(gammas) // 2)
    nu_x = core.fit_nu_loglog(gammas[-ntail:], xnorm[-ntail:])["nu_global"]
    nu_F = core.fit_nu_loglog(gammas[-ntail:], fqs[-ntail:])["nu_global"]

    return dict(
        nu_x=float(nu_x), nu_F=float(nu_F),
        nu_obs_pred=float(2 * nu_x),
        ratio_F_over_2x=float(nu_F / (2 * nu_x)) if nu_x else None,
        quadratic_inheritance_ok=bool(abs(nu_F - 2 * nu_x) < 0.15),
        gammas=gammas.tolist(), xnorm=xnorm.tolist(), fqs=fqs.tolist(),
    )


# ----------------------------------------------------------------------
# 4. frozen-source difference (Phase P): cancellation-promoted order
# ----------------------------------------------------------------------
def run_frozen_source():
    import model_physical as mp
    import sympy as sp

    # generic: J45 = 1/2, phi = 0  -> order 3
    nu_generic, _Qg, _Ng = mp.certificate_R_S(0, sp.Rational(1, 2), 0)
    # tuned:   J45 = J45_star, phi = pi -> cancellation promotes order 3 -> 4
    J45s = mp.J45_star()
    nu_tuned, _Qt, _Nt = mp.certificate_R_S(0, J45s, sp.pi)

    return dict(
        nu_generic=int(nu_generic), nu_tuned=int(nu_tuned),
        J45_star=str(J45s),
        promotion_ok=bool(nu_tuned == nu_generic + 1),
    )


# ----------------------------------------------------------------------
# 5. Phase N modifiers: q-fan and z* resonance
# ----------------------------------------------------------------------
def _expected_fan(q):
    q = Fraction(q)
    if q <= 1:
        return 4 - q
    if q <= 2:
        return 2 + q
    return Fraction(4)


def run_qfan():
    import phase_n_exact_core as pn
    num, den = pn.master_polynomials()
    qs = [Fraction(0), Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3)]
    fan = []
    ok = True
    for q in qs:
        order, _cn, _cd = pn.path_order(q, num, den)
        exp = _expected_fan(q)
        match = (order == exp)
        ok = ok and match
        fan.append(dict(q=str(q), nu=str(order), expected=str(exp), match=bool(match)))
    return dict(fan=fan, fan_ok=bool(ok))


def run_zstar():
    import phase_n_frequency_core as fq
    moments = fq.ideal_moment_polynomials(kmax=7)
    z_generic = Fraction(1)
    z_star = Fraction(543, 280)
    order_generic = fq.first_nonzero_moment_at(moments, z_generic)
    order_star = fq.first_nonzero_moment_at(moments, z_star)
    return dict(
        z_generic=str(z_generic), order_generic=order_generic,
        z_star=str(z_star), order_star=order_star,
        resonance_promotion_ok=bool(order_star is not None and order_generic is not None
                                    and order_star > order_generic),
    )


# ----------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------
def make_figure(curves, path_png):
    plt.figure(figsize=(7.2, 5.4))
    for label, (g, y) in curves.items():
        plt.loglog(g, y, "o-", ms=2.5, lw=1.0, label=label)
    plt.xlabel(r"$\Gamma$")
    plt.ylabel(r"$|R_{\mathrm{obs}}(\Gamma)|$")
    plt.title("Gate A: observable-order inheritance "
              r"($\nu_{\mathrm{obs}}=n_{12}+n_{21}-\nu_{\rm den}$)")
    plt.legend(fontsize=7, loc="lower left")
    plt.tight_layout()
    plt.savefig(path_png, dpi=200)
    plt.savefig(path_png.replace(".png", ".pdf"))
    plt.close()


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Gate A: observable-order inheritance")
    ap.add_argument("--quick", action="store_true", help="fast CI-sized sweep")
    ap.add_argument("--smoke", action="store_true", help="alias of --quick")
    args = ap.parse_args()
    quick = args.quick or args.smoke

    t0 = time.time()
    bilinear_rows, curves, reg_nv = run_bilinear(quick=quick)
    nv_exact = run_nv_exact()
    qfi = run_qfi(quick=quick)
    frozen = run_frozen_source()
    qfan = run_qfan()
    zstar = run_zstar()

    # -- gate assembly ------------------------------------------------
    bilinear_ok = all(r["agree"] for r in bilinear_rows)

    # G-O1: inheritance law for >= 3 distinct observables
    #   (bilinear susceptibility, quadratic QFI, frozen-source difference)
    g_o1 = bool(bilinear_ok and nv_exact["inheritance_ok"]
                and qfi["quadratic_inheritance_ok"] and frozen["promotion_ok"])

    # G-O2: generic vs protected distinguished by exact certificate
    syn = {r["label"]: r for r in bilinear_rows}
    generic = next((r for r in bilinear_rows if r["mechanism"] == "generic"), None)
    protected = [r for r in bilinear_rows if r["mechanism"] == "symmetry-protected"]
    g_o2 = bool(
        generic is not None and generic["nu_obs_cert"] == 2
        and any(r["nu_obs_cert"] == 4 for r in protected)
        and any(r["nu_obs_cert"] == 6 for r in protected)
        and nv_exact["M0_exact_zero"]
    )

    # G-O3: predicted == full log-log tail slope (per-model fit_agree)
    g_o3 = bool(all(r["fit_agree"] for r in bilinear_rows))

    # G-O4: pre-asymptotic vs asymptotic separated on physical NV
    g_o4 = bool(
        reg_nv["gamma_cross"] is not None
        and reg_nv["is_preasymptotic"] is True
        and abs(reg_nv["nu_asymptotic"] - 4.0) < 0.1
    )

    # G-O5: three-route consistency + Phase N modifiers
    g_o5 = bool(bilinear_ok and qfan["fan_ok"] and zstar["resonance_promotion_ok"])

    gates = dict(
        G_O1_inheritance_three_observables=g_o1,
        G_O2_generic_vs_protected_exact=g_o2,
        G_O3_predict_matches_loglog=g_o3,
        G_O4_preasymptotic_separated=g_o4,
        G_O5_three_route_and_modifiers=g_o5,
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        description="Gate A: observable-order inheritance nu_obs = n12 + n21 - nu_den",
        quick=quick,
        bilinear_models=bilinear_rows,
        nv_exact_certificate=nv_exact,
        qfi_quadratic=_clean({k: v for k, v in qfi.items()
                              if k not in ("gammas", "xnorm", "fqs")}),
        frozen_source_difference=frozen,
        phase_n_q_fan=qfan,
        phase_n_z_star=zstar,
        preasymptotic_separation_NV=_clean(reg_nv),
        gates=gates,
        runtime_s=round(time.time() - t0, 2),
    )

    out_json = os.path.join(TABLES, "gates_summary_gateO.json")
    with open(out_json, "w") as f:
        json.dump(_clean(summary), f, indent=2)

    # per-model CSV
    out_csv = os.path.join(TABLES, "gate_a_models.csv")
    with open(out_csv, "w", newline="") as f:
        cols = ["label", "readout_mode", "mechanism", "n12", "n21", "nu_den",
                "nu_obs_pred", "nu_obs_cert", "nu_obs_tail_fit",
                "pred_cert_agree", "fit_agree", "agree"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in bilinear_rows:
            w.writerow({k: r.get(k) for k in cols})

    fig_png = os.path.join(FIGS, "fig_gateO_inheritance.png")
    make_figure(curves, fig_png)

    print(json.dumps(gates, indent=2))
    print(f"\nbilinear models:")
    for r in bilinear_rows:
        print(f"  {r['label']:42s} pred={r['nu_obs_pred']} cert={r['nu_obs_cert']} "
              f"tail_fit={r['nu_obs_tail_fit']:.3f} agree={r['agree']}")
    print(f"NV exact: nu_R={nv_exact['nu_R']} = nu_K12+nu_K21="
          f"{nv_exact['nu_K12']}+{nv_exact['nu_K21']}  M0=0:{nv_exact['M0_exact_zero']}")
    print(f"QFI: nu_x={qfi['nu_x']:.3f} nu_F={qfi['nu_F']:.3f} "
          f"(nu_F/2nu_x={qfi['ratio_F_over_2x']:.3f})")
    print(f"frozen-source: generic nu={frozen['nu_generic']} -> tuned nu={frozen['nu_tuned']}")
    print(f"q-fan ok: {qfan['fan_ok']};  z* {zstar['order_generic']}->{zstar['order_star']}")
    print(f"\nwrote {out_json}\nwrote {out_csv}\nwrote {fig_png}")
    print(f"runtime {summary['runtime_s']} s")
    return summary


if __name__ == "__main__":
    main()
