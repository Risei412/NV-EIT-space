"""run_gate_b.py -- Gate B driver: non-EIT / non-diamond superconducting
witness with a pre-registered BLIND prediction (PRL Priorities 2 & 7).

The witness is a superconducting dissipative-state-transfer through two lossy
bus modes (model_sc_transfer.py). The observable is the A -> B transfer
amplitude K(kappa) = p_B^dag [kappa D + A0]^-1 c_A and its efficiency |K|^2 --
a conversion/transport quantity, NOT an EIT dark-state transparency. The same
path-moment law predicts an integer dissipation order set by the leg
selection rule:
    generic  (sum gAi gBi != 0):  nu_K = 1,  efficiency 2
    protected(sum gAi gBi  = 0):  nu_K = 2,  efficiency 4

BLIND PROTOCOL (P7): the exact kappa-degree nu_K is computed from the reduced
pencil (H, jumps, input, readout) ONLY and recorded in `blind_prediction`
BEFORE the full-GKSL sweep; the full Lindblad log-log slope then confirms it.

Gates (strategy doc Sec.13 P2/P7):
  G-B1  non-EIT witness well-posed: transfer observable, full GKSL built,
        reduced pencil == full-Liouvillian coherence sector (~1e-9).
  G-B2  blind prediction fixed (nu_K, efficiency 2*nu_K) before the full sweep.
  G-B3  full-GKSL deep-tail log-log slope matches the blind 2*nu_K.
  G-B4  selection rule changes the order (generic nu_K=1 vs protected nu_K=2),
        distinguished by the EXACT rational-degree certificate.
  G-B5  small symmetry breaking eps -> high->low crossover, kappa*(eps) grows
        as eps -> 0 (protected recovered).

Usage:
    python run_gate_b.py [--quick] [--smoke]

Outputs:
    results/tables/gates_summary_gateB.json
    results/tables/gate_b_scaling.csv
    results/figures/fig_gateB_transfer_scaling.png (+ .pdf)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE_SRC = os.path.join(HERE, "..", "..", "src")
PHASEO_SRC = os.path.join(HERE, "..", "..", "PhaseO_observable_inheritance", "src")
for _p in (HERE, PHASE_SRC, PHASEO_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import core                              # noqa: E402
import gate_a_observable as gao          # noqa: E402
import model_sc_transfer as m            # noqa: E402

RESULTS = os.path.join(HERE, "..", "results")
TABLES = os.path.join(RESULTS, "tables")
FIGS = os.path.join(RESULTS, "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)


def _sectorspec(tuning, eps=0.0):
    """Map the transfer kernel onto a gate_a SectorSpec (quadratic readout =
    transfer efficiency |K|^2). dc = source A, dp = readout B."""
    Dnum = m.D_reduced()
    B_of_z = m.B_of_z_factory(tuning, eps)
    D_sym, B_sym_of_z, c_sym, p_sym = m.symbolic_pencil(tuning, eps)
    return gao.SectorSpec(
        D=Dnum, B_of_z=B_of_z, dp=m.p_readout(), dc=m.c_source(),
        g_eff=0.0, beta=1.0, readout_mode="quadratic",
        label=f"SC-transfer/{tuning}",
        D_sym=D_sym, B_sym_of_z=B_sym_of_z, dp_sym=p_sym, dc_sym=c_sym,
    )


# ----------------------------------------------------------------------
# G-B2 / G-B4: blind prediction from the reduced pencil only
# ----------------------------------------------------------------------
def blind_prediction():
    Gamma = sp.symbols("Gamma")
    out = {}
    for tuning in ("generic", "protected"):
        D_sym, B_sym_of_z, c_sym, p_sym = m.symbolic_pencil(tuning)
        nu_K, _Q, _N = core.certificate_deg_nu(D_sym, B_sym_of_z, c_sym, p_sym,
                                               Gamma, 0, z_sym=None)
        out[tuning] = dict(
            nu_K=int(nu_K),
            nu_efficiency_pred=int(2 * nu_K),
        )
    return out


# ----------------------------------------------------------------------
# G-B1: reduced pencil == full-GKSL coherence sector
# ----------------------------------------------------------------------
def reduced_vs_full(kappas):
    rows = []
    for tuning in ("generic", "protected"):
        errs = []
        for k in kappas:
            Kred = m.transfer_kernel(k, tuning=tuning)
            Kfull = m.full_transfer_amplitude(k, tuning=tuning) / (-1j)  # remove i[V,.] drive factor
            errs.append(abs(Kfull - Kred) / max(abs(Kred), 1e-300))
        rows.append(dict(tuning=tuning, max_rel_err=float(max(errs))))
    resid = float(m.verify_rho0_steady(kappas[len(kappas) // 2]))
    return rows, resid


# ----------------------------------------------------------------------
# G-B3: full-GKSL deep-tail slope vs the blind prediction
# ----------------------------------------------------------------------
def full_model_slopes(kappas_tail):
    rows = {}
    curves = {}
    for tuning in ("generic", "protected"):
        Kf = np.array([m.full_transfer_amplitude(k, tuning=tuning) for k in kappas_tail])
        amp = core.fit_nu_loglog(kappas_tail, Kf)["nu_global"]
        eff = core.fit_nu_loglog(kappas_tail, np.abs(Kf) ** 2)["nu_global"]
        # exact efficiency certificate via the reused Gate A harness (quadratic)
        cert = gao.certify_nu_obs_exact(_sectorspec(tuning), z_val=0.0)
        rows[tuning] = dict(amp_slope=float(amp), eff_slope=float(eff),
                            eff_cert=int(cert["nu_cert"]), nu_K_cert=int(cert["nu_cert"] // 2))
        curves[tuning] = (kappas_tail, np.abs(Kf) ** 2)
    return rows, curves


# ----------------------------------------------------------------------
# G-B5: symmetry-breaking crossover kappa*(eps)
# ----------------------------------------------------------------------
def crossover_fan(epsilons, quick=False):
    kappas = np.logspace(4, 12, 60 if quick else 130)
    rows = []
    curves = {}
    for eps in epsilons:
        vals = np.array([m.transfer_kernel(k, tuning="broken", eps=eps) for k in kappas])
        fit = core.fit_nu_loglog(kappas, vals)
        nu_eff, gmid = fit["nu_eff"], fit["gamma_mid"]
        plateau = float(nu_eff[np.argmin(np.abs(gmid - 3e4))])  # in the protected plateau
        below = np.where(nu_eff < 1.5)[0]
        kstar = float(gmid[below[0]]) if len(below) else float("inf")
        rows.append(dict(eps=float(eps), nu_eff_plateau=plateau, kappa_star=kstar))
        curves[eps] = (gmid, nu_eff)
    # protected (eps=0) stays at 2
    vals0 = np.array([m.transfer_kernel(k, tuning="protected") for k in kappas])
    nu_eff0 = core.fit_nu_loglog(kappas, vals0)["nu_eff"][-1]
    # kappa*(eps) should grow (approx 1/eps) as eps -> 0
    finite = [r for r in rows if np.isfinite(r["kappa_star"])]
    monotonic = all(finite[i]["kappa_star"] < finite[i + 1]["kappa_star"]
                    for i in range(len(finite) - 1)) if len(finite) > 1 else False
    # eps sorted descending -> kappa* ascending
    return rows, curves, float(nu_eff0), bool(monotonic)


# ----------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------
def make_figure(eff_curves, cross_curves, path_png):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for tuning, (k, y) in eff_curves.items():
        ax1.loglog(k, y, "o-", ms=2.5, lw=1.0, label=f"{tuning}")
    k = list(eff_curves.values())[0][0]
    y0 = list(eff_curves["generic"][1])[0]
    ax1.loglog(k, y0 * (k[0] / k) ** 2, "k--", lw=0.8, label=r"$\kappa^{-2}$ (generic)")
    yp = list(eff_curves["protected"][1])[0]
    ax1.loglog(k, yp * (k[0] / k) ** 4, "k:", lw=0.8, label=r"$\kappa^{-4}$ (protected)")
    ax1.set_xlabel(r"$\kappa$ (GHz)"); ax1.set_ylabel(r"transfer efficiency $|K|^2$")
    ax1.set_title("Gate B: superconducting transfer, selection-rule order")
    ax1.legend(fontsize=8)

    for eps, (gm, ne) in cross_curves.items():
        ax2.semilogx(gm, ne, "-", lw=1.0, label=fr"$\epsilon={eps:.0e}$")
    ax2.axhline(2, color="gray", ls="--", lw=0.7)
    ax2.axhline(1, color="gray", ls=":", lw=0.7)
    ax2.set_xlabel(r"$\kappa$ (GHz)"); ax2.set_ylabel(r"effective order $\nu_{\rm eff}$")
    ax2.set_title(r"symmetry-breaking crossover $\nu:2\to1$, $\kappa_\ast(\epsilon)\propto1/\epsilon$")
    ax2.set_ylim(0.5, 2.5); ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(path_png, dpi=200)
    fig.savefig(path_png.replace(".png", ".pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Gate B: superconducting non-EIT witness")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    quick = args.quick or args.smoke

    t0 = time.time()

    # --- BLIND: fix the prediction before any full-model fit ---
    blind = blind_prediction()

    kappas_chk = np.logspace(3, 7, 5)
    rvf_rows, rho0_resid = reduced_vs_full(kappas_chk)

    kappas_tail = np.logspace(5, 9, 25 if quick else 50)
    full_rows, eff_curves = full_model_slopes(kappas_tail)

    epsilons = [1e-8, 1e-9, 1e-10, 1e-11]
    cross_rows, cross_curves, nu_eff0_protected, kstar_monotonic = crossover_fan(epsilons, quick=quick)

    # -- gate assembly ------------------------------------------------
    g_b1 = bool(all(r["max_rel_err"] < 1e-7 for r in rvf_rows) and rho0_resid < 1e-9)

    g_b2 = bool(blind["generic"]["nu_K"] == 1 and blind["protected"]["nu_K"] == 2
                and blind["generic"]["nu_efficiency_pred"] == 2
                and blind["protected"]["nu_efficiency_pred"] == 4)

    g_b3 = bool(
        abs(full_rows["generic"]["eff_slope"] - blind["generic"]["nu_efficiency_pred"]) < 0.1
        and abs(full_rows["protected"]["eff_slope"] - blind["protected"]["nu_efficiency_pred"]) < 0.1
    )

    g_b4 = bool(
        full_rows["generic"]["eff_cert"] == 2 and full_rows["protected"]["eff_cert"] == 4
        and full_rows["generic"]["eff_cert"] != full_rows["protected"]["eff_cert"]
    )

    g_b5 = bool(kstar_monotonic and abs(nu_eff0_protected - 2.0) < 0.05
                and all(r["nu_eff_plateau"] > 1.5 for r in cross_rows))

    gates = dict(
        G_B1_reduced_equals_full=g_b1,
        G_B2_blind_prediction_fixed=g_b2,
        G_B3_full_model_matches_blind=g_b3,
        G_B4_selection_rule_changes_order=g_b4,
        G_B5_symmetry_breaking_crossover=g_b5,
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        description="Gate B: non-EIT non-diamond superconducting dissipative-transfer "
                    "witness; path-moment law predicts integer transfer-efficiency order.",
        witness="superconducting transmon + two lossy bus modes (dissipative state transfer)",
        observable="A->B transfer amplitude K(kappa)=p_B^dag[kappa D+A0]^-1 c_A; efficiency |K|^2",
        blind_prediction=blind,
        reduced_vs_full=dict(rows=rvf_rows, rho0_steady_residual=rho0_resid),
        full_model_slopes=full_rows,
        crossover_fan=dict(rows=cross_rows, nu_eff_protected_eps0=nu_eff0_protected,
                           kappa_star_monotonic_in_1_over_eps=kstar_monotonic),
        quick=quick,
        gates=gates,
        runtime_s=round(time.time() - t0, 2),
    )

    out_json = os.path.join(TABLES, "gates_summary_gateB.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    out_csv = os.path.join(TABLES, "gate_b_scaling.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tuning", "amp_slope", "eff_slope", "eff_cert", "blind_nu_K", "blind_eff"])
        for tuning in ("generic", "protected"):
            fr = full_rows[tuning]
            w.writerow([tuning, f"{fr['amp_slope']:.4f}", f"{fr['eff_slope']:.4f}",
                        fr["eff_cert"], blind[tuning]["nu_K"], blind[tuning]["nu_efficiency_pred"]])

    fig_png = os.path.join(FIGS, "fig_gateB_transfer_scaling.png")
    make_figure(eff_curves, cross_curves, fig_png)

    print(json.dumps(gates, indent=2))
    print("\nBLIND prediction (fixed before full sweep):")
    for t in ("generic", "protected"):
        print(f"  {t:10s} nu_K={blind[t]['nu_K']}  efficiency order={blind[t]['nu_efficiency_pred']}")
    print("full-GKSL deep-tail slopes:")
    for t in ("generic", "protected"):
        fr = full_rows[t]
        print(f"  {t:10s} amp={fr['amp_slope']:.3f} eff={fr['eff_slope']:.3f} (cert eff={fr['eff_cert']})")
    print(f"reduced-vs-full max rel err: "
          f"{max(r['max_rel_err'] for r in rvf_rows):.1e}; rho0 residual {rho0_resid:.1e}")
    print("crossover kappa*(eps):")
    for r in cross_rows:
        print(f"  eps={r['eps']:.0e}  plateau nu_eff={r['nu_eff_plateau']:.2f}  kappa*={r['kappa_star']:.2e}")
    print(f"protected(eps=0) nu_eff->{nu_eff0_protected:.2f};  kappa* monotonic: {kstar_monotonic}")
    print(f"\nwrote {out_json}\nwrote {out_csv}\nwrote {fig_png}")
    print(f"runtime {summary['runtime_s']} s")
    return summary


if __name__ == "__main__":
    main()
