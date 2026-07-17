"""Phase D driver: phenomena invisible to the old no-go framework
(phase_d_beyond_old_theory_plan.md, runs D1-D4).

Outputs into ../results/figures and ../results/gates_summary_phaseD.json,
and appends a section to ../results/summary.md.
"""

import json
import os

import numpy as np
from scipy.optimize import brentq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core import fit_nu_loglog, condition_number_of_projector

import model_oblique as D1
import model_paths as D3
import model_cancel as D2
import model_masquerade as D4

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
FIGDIR = os.path.join(RESULTS, "figures")
os.makedirs(FIGDIR, exist_ok=True)


def _default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    return str(o)


# ----------------------------------------------------------------------
# D1: oblique Riesz protection vs the old orthogonal-projector formula
# ----------------------------------------------------------------------
def run_D1():
    gate = {}
    a_values = np.linspace(0.0, 5.0, 26)
    scan = D1.scan_a(a_values)

    # Hermitian control: Riesz == orthogonal at a=0
    control_err = abs(scan[0]["r0_riesz"] - scan[0]["r0_orth"])

    # false-suppression point: r0_orth crosses zero, r0_riesz stays O(1)
    def r0_orth_of_a(a):
        return D1.scan_a([a])[0]["r0_orth"]

    vals = [s["r0_orth"] for s in scan]
    a_star = None
    for i in range(len(a_values) - 1):
        if vals[i] * vals[i + 1] < 0:
            a_star = brentq(r0_orth_of_a, a_values[i], a_values[i + 1], xtol=1e-12)
            break
    res_star = D1.scan_a([a_star])[0] if a_star is not None else None

    # plateau verification at a=0 (control) and a=2 (oblique) against both formulas
    Gamma_max = 1e7
    plateau_checks = {}
    for a in (0.0, 1.0, 2.0, 5.0):
        _, _, R = D1.R_S_gamma(Gamma_max, a)
        s = D1.scan_a([a])[0]
        plateau_checks[str(a)] = {
            "plateau": abs(R),
            "r0_riesz": abs(s["r0_riesz"]),
            "r0_orth": abs(s["r0_orth"]),
            "err_riesz": abs(abs(R) - abs(s["r0_riesz"])) / max(abs(s["r0_riesz"]), 1e-12),
            "err_orth": abs(abs(R) - abs(s["r0_orth"])) / max(abs(s["r0_orth"]), 1e-12),
        }

    # non-monotonic dissipation dependence (interior maximum in |R_S(Gamma)|)
    a_bump, dQ1_bump, dQ2_bump = 4.0, 0.05, 3.0
    gammas_bump = np.logspace(0.5, 6, 60)
    Rs_bump = np.array([
        abs(D1.R_S_gamma_custom(g, a_bump, dQ1=dQ1_bump, dQ2=dQ2_bump)
            if hasattr(D1, "R_S_gamma_custom") else
            _r_s_custom(g, a_bump, dQ1_bump, dQ2_bump))
        for g in gammas_bump
    ])
    imax = int(np.argmax(Rs_bump))
    interior_bump = 2 < imax < len(gammas_bump) - 3

    gate["control_hermitian_a0_riesz_vs_orth_diff"] = float(control_err)
    gate["control_pass"] = bool(control_err < 1e-8)
    gate["false_suppression_a_star"] = float(a_star) if a_star is not None else None
    gate["false_suppression_r0_riesz"] = abs(res_star["r0_riesz"]) if res_star else None
    gate["false_suppression_r0_orth"] = abs(res_star["r0_orth"]) if res_star else None
    gate["false_suppression_pass"] = bool(
        res_star is not None and abs(res_star["r0_orth"]) < 1e-6 and abs(res_star["r0_riesz"]) > 1.0
    )
    gate["plateau_checks"] = plateau_checks
    gate["plateau_pass"] = bool(all(v["err_riesz"] < 1e-4 for v in plateau_checks.values()))
    gate["nonmonotonic_bump"] = {
        "a": a_bump, "Gamma_at_max": float(gammas_bump[imax]),
        "R_max": float(Rs_bump[imax]), "R_first": float(Rs_bump[0]), "R_last": float(Rs_bump[-1]),
        "interior_max": bool(interior_bump),
    }
    gate["overall_pass"] = bool(gate["control_pass"] and gate["false_suppression_pass"]
                                 and gate["plateau_pass"] and interior_bump)

    # figures
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    a_plot = np.array([s["a"] for s in scan])
    r0_riesz_plot = np.array([s["r0_riesz"] for s in scan])
    r0_orth_plot = np.array([s["r0_orth"] for s in scan])
    axes[0].plot(a_plot, r0_riesz_plot, "o-", label=r"$r_0$ (Riesz, correct)")
    axes[0].plot(a_plot, r0_orth_plot, "s--", label=r"$r_0$ (orthogonal, old)")
    if a_star is not None:
        axes[0].axvline(a_star, color="k", ls=":", lw=1, label=r"$a^*$ (false suppression)")
    axes[0].axhline(0, color="gray", lw=0.5)
    axes[0].set_xlabel("$a$ (non-normality of $D$)")
    axes[0].set_ylabel(r"$r_0$")
    axes[0].set_title("D1a: oblique vs orthogonal protected coefficient")
    axes[0].legend(fontsize=8)

    axes[1].loglog(gammas_bump, Rs_bump, "o-", ms=3, label=f"$a={a_bump}$ (oblique)")
    gammas_ctrl = gammas_bump
    Rs_ctrl = np.array([_r_s_custom(g, 0.0, dQ1_bump, dQ2_bump) for g in gammas_ctrl])
    axes[1].loglog(gammas_ctrl, np.abs(Rs_ctrl), "s--", ms=3, label="$a=0$ (Hermitian control)")
    axes[1].set_xlabel(r"$\Gamma$")
    axes[1].set_ylabel(r"$|R_S(\Gamma)|$")
    axes[1].set_title("D1b: dissipation-enhanced sector response")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figD1_oblique_protection.png"), dpi=150)
    plt.close(fig)

    return gate


def _r_s_custom(Gamma, a, dQ1, dQ2):
    D = D1.D_of_a(a)
    B_full = D1.make_B(D1.KAPPA_FULL, 1.3, 0.9, 0.55, dQ1, dQ2)
    B_cut = D1.make_B(0.0, 1.3, 0.9, 0.55, dQ1, dQ2)
    A_full = Gamma * D + B_full
    A_cut = Gamma * D + B_cut
    cf = np.conj(D1.P_VEC) @ np.linalg.solve(A_full, D1.C_VEC)
    cc = np.conj(D1.P_VEC) @ np.linalg.solve(A_cut, D1.C_VEC)
    return cf - cc


# ----------------------------------------------------------------------
# D3: scaling-path dependence of the class
# ----------------------------------------------------------------------
def run_D3():
    gate = {}
    gate["calibration_max_diff"] = D3.calibration_check()
    cfT, cfU = D3.chi_full_at_calibration()
    gate["chi_full_agreement_at_Gamma0"] = abs(cfT - cfU)

    gammas = np.logspace(1, 6, 25)
    nus = {}
    for path in ("T", "U"):
        Rs = [D3.R_S_gamma(g, path)[2] for g in gammas]
        fit = fit_nu_loglog(gammas, Rs)
        nus[path] = fit["nu_global"]
    gate["nu_path_T"] = nus["T"]
    gate["nu_path_U"] = nus["U"]
    gate["paths_disagree"] = bool(abs(nus["T"] - nus["U"]) > 0.5)

    thetas = np.linspace(0.0, np.pi / 2, 31)
    nu_fan = []
    for th in thetas:
        Rs = [D3.R_S_gamma_fan(g, th) for g in gammas]
        fit = fit_nu_loglog(gammas, Rs)
        nu_fan.append(fit["nu_global"])
    nu_fan = np.array(nu_fan)
    gate["fan_theta0_nu"] = float(nu_fan[0])
    gate["fan_interior_nu_mean"] = float(np.mean(nu_fan[5:-5]))
    gate["fan_shows_transition"] = bool(nu_fan[0] < 0.3 and np.mean(nu_fan[5:-5]) > 1.5)

    gate["overall_pass"] = bool(
        gate["calibration_max_diff"] < 1e-9 and gate["chi_full_agreement_at_Gamma0"] < 1e-9
        and gate["paths_disagree"] and gate["fan_shows_transition"]
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for path, marker in [("T", "o-"), ("U", "s-")]:
        Rs = [abs(D3.R_S_gamma(g, path)[2]) for g in gammas]
        axes[0].loglog(gammas, Rs, marker, label=f"path {path} ($\\nu={nus[path]:.2f}$)")
    axes[0].set_xlabel(r"$\Gamma$")
    axes[0].set_ylabel(r"$|R_S(\Gamma)|$")
    axes[0].set_title("D3a: two scaling paths, one physical point")
    axes[0].legend(fontsize=8)

    axes[1].plot(thetas, nu_fan, "o-")
    axes[1].set_xlabel(r"$\theta$ (scaling-path direction)")
    axes[1].set_ylabel(r"$\nu[R_S]$")
    axes[1].set_title("D3b: class fan vs path direction")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figD3_scaling_path.png"), dpi=150)
    plt.close(fig)

    return gate


# ----------------------------------------------------------------------
# D2: cancellation-promoted suppression / symmetry-free exact zero
# ----------------------------------------------------------------------
def run_D2():
    gate = {}
    v23_star = D2.tune_v23_for_small_m2()
    D, B_of_z, c, p = D2.build_promotion_model(v23=v23_star, v12=0.05)
    mu = D2.moment_ladder(D, B_of_z, c, p, z=0.0, kmax=4)

    gammas = np.logspace(1, 7, 30)
    Rs = []
    for g in gammas:
        A = g * D + B_of_z(0.0)
        Rs.append(np.conj(p) @ np.linalg.solve(A, c))
    fit = fit_nu_loglog(gammas, Rs)

    gate["v23_star"] = float(v23_star)
    gate["moments_v12_0p05"] = [complex(m) for m in mu]
    gate["nu_direct_fit"] = fit["nu_global"]
    gate["naive_graph_bound_2hop"] = 3
    gate["order_promoted"] = bool(fit["nu_global"] > 3.5)

    # control: detune v23 away from cancellation -> reopens at generic order
    D_ctrl, B_ctrl, c_ctrl, p_ctrl = D2.build_promotion_model(v23=v23_star + 0.1, v12=0.05)
    Rs_ctrl = []
    for g in gammas:
        A = g * D_ctrl + B_ctrl(0.0)
        Rs_ctrl.append(np.conj(p_ctrl) @ np.linalg.solve(A, c_ctrl))
    fit_ctrl = fit_nu_loglog(gammas, Rs_ctrl)
    gate["control_detuned_nu"] = fit_ctrl["nu_global"]
    gate["control_demotes_to_naive_order"] = bool(abs(fit_ctrl["nu_global"] - 3.0) < 0.3)

    # exact symbolic Class I certificate, no reducing symmetry
    cert = D2.build_exact_zero_certificate()
    moments_zero = all(m == 0 for m in cert["moments"])
    import sympy as sp
    Mfull_t1 = cert["M_full_symbolic"].subs(cert["t"], 1)
    K01 = D2.cross_kernel(Mfull_t1, 0, 1, sp.Rational(1, 2))
    K01_cut = D2.cross_kernel(cert["M_cut"], 0, 1, sp.Rational(1, 2))

    gate["exact_class_I_moments_zero_for_all_t"] = bool(moments_zero)
    gate["exact_class_I_K01_full"] = float(K01)
    gate["exact_class_I_K01_cut"] = float(K01_cut)
    gate["exact_class_I_kernel_nonzero"] = bool(abs(float(K01)) > 1e-6)

    gate["overall_pass"] = bool(
        gate["order_promoted"] and gate["control_demotes_to_naive_order"]
        and gate["exact_class_I_moments_zero_for_all_t"] and gate["exact_class_I_kernel_nonzero"]
    )

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(gammas, np.abs(Rs), "o-", label=f"near-cancellation ($\\nu\\approx{fit['nu_global']:.1f}$, naive bound 3)")
    ax.loglog(gammas, np.abs(Rs_ctrl), "s--", label=f"detuned control ($\\nu\\approx{fit_ctrl['nu_global']:.1f}$)")
    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"$|R_S(\Gamma)|$")
    ax.set_title("D2: cancellation-promoted suppression order")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figD2_cancellation_promotion.png"), dpi=150)
    plt.close(fig)

    return gate


# ----------------------------------------------------------------------
# D4: preasymptotic masquerade vs finite polynomial certificate
# ----------------------------------------------------------------------
def run_D4():
    gate = {}
    exact_moments = D2  # placeholder not used
    moments_at_star = D4.symbolic_cancellation_check()
    gate["exact_moments_at_v02_star"] = [str(m) for m in moments_at_star]
    gate["m2_exactly_zero_at_v02_star"] = bool(moments_at_star[2] == 0)

    delta = 1e-5
    v02 = D4.V02_STAR + delta
    gammas = np.logspace(0, 8, 80)
    Rs = np.array([D4.R_gamma(g, v02) for g in gammas])
    loglog = fit_nu_loglog(gammas, Rs)
    nu_eff = loglog["nu_eff"]
    gamma_mid = loglog["gamma_mid"]

    plateau_mask = np.abs(nu_eff - 4.0) < 0.05
    plateau_decades = 0.0
    if np.any(plateau_mask):
        idx = np.where(plateau_mask)[0]
        plateau_decades = float(np.log10(gamma_mid[idx[-1]]) - np.log10(gamma_mid[idx[0]]))
    true_nu_tail = float(nu_eff[-1])

    cert_nu = D4.certificate_true_nu(delta)

    gate["delta"] = delta
    gate["false_plateau_decades_at_nu4"] = plateau_decades
    gate["true_nu_tail_measured"] = true_nu_tail
    gate["certificate_nu_exact"] = cert_nu
    gate["certificate_matches_tail"] = bool(abs(cert_nu - true_nu_tail) < 0.1)
    gate["overall_pass"] = bool(
        plateau_decades > 1.0 and gate["certificate_matches_tail"] and cert_nu == 3
    )

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.semilogx(gamma_mid, nu_eff, "o-", ms=3, label=r"$\nu_{\rm eff}(\Gamma)$ (finite-window fit)")
    ax.axhline(4, color="r", ls="--", lw=1, label="false plateau ($\\nu=4$)")
    ax.axhline(cert_nu, color="g", ls="--", lw=1, label=f"exact certificate ($\\nu={cert_nu}$)")
    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"$\nu_{\rm eff}$")
    ax.set_title("D4: preasymptotic masquerade")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figD4_masquerade.png"), dpi=150)
    plt.close(fig)

    return gate


def main():
    print("Running D1 (oblique Riesz protection)...")
    g1 = run_D1()
    print("Running D3 (scaling-path dependence)...")
    g3 = run_D3()
    print("Running D2 (cancellation-promoted suppression)...")
    g2 = run_D2()
    print("Running D4 (preasymptotic masquerade)...")
    g4 = run_D4()

    summary = {
        "D1_oblique_protection": g1,
        "D3_scaling_path": g3,
        "D2_cancellation_promotion": g2,
        "D4_masquerade": g4,
        "phase_d_overall_pass": bool(g1["overall_pass"] and g2["overall_pass"]
                                      and g3["overall_pass"] and g4["overall_pass"]),
    }

    with open(os.path.join(RESULTS, "gates_summary_phaseD.json"), "w") as f:
        json.dump(summary, f, indent=2, default=_default)

    print("\n=== Phase D summary ===")
    for k in ("D1_oblique_protection", "D3_scaling_path", "D2_cancellation_promotion", "D4_masquerade"):
        print(f"{k}: overall_pass = {summary[k]['overall_pass']}")
    print("phase_d_overall_pass =", summary["phase_d_overall_pass"])


if __name__ == "__main__":
    main()
