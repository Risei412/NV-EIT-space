"""Phase B driver: singular-D models, Figures 1/3/4, Gates 1/2/4.

Outputs into ../results/figures and ../results/tables (relative to this
file), plus ../results/gates_summary.json.
"""

import json
import os

import numpy as np
from scipy.optimize import brentq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core import (
    transfer, krylov_class1_certificate, nu_from_moments,
    fit_nu_loglog, protected_F0, riesz_projection_diag,
)
from model_protected import (
    build_class1_example, build_class2_example, build_class3_example,
    transition_matrices, r0_of_lambda, R_S_gamma,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
FIGDIR = os.path.join(RESULTS, "figures")
TABDIR = os.path.join(RESULTS, "tables")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TABDIR, exist_ok=True)

GAMMAS_FIG1 = np.logspace(1, 6, 25)


# ----------------------------------------------------------------------
# Figure 1: theory unit test (Class I, II nu=1, II nu=2, III)
# ----------------------------------------------------------------------
def figure1_and_gate1():
    gate1 = {}

    # -- Class I: Krylov-orthogonal, embedded in a Gamma-family with D=I so
    #    that R_S(Gamma) = F(Gamma) is EXACTLY zero for all Gamma (Theorem I
    #    (iv)), not merely asymptotically small.
    M, r, l = build_class1_example()
    n = M.shape[0]
    max_mu, moments = krylov_class1_certificate(M, r, l)
    D_I = np.eye(n)
    vals_I = []
    for g in GAMMAS_FIG1:
        A = g * D_I + M  # A_Gamma = Gamma*I + M; class unaffected by shift
        vals_I.append(transfer(A, l, r))
    vals_I = np.abs(vals_I)
    gate1["class_I"] = {
        "max_abs_moment": max_mu,
        "certificate": "Class I (exact, machine precision)" if max_mu < 1e-10 else "FAILED",
        "max_abs_F_over_gammas": float(np.max(vals_I)),
    }

    # -- Class II, nu=1 (d=0: c and p on the same site) and nu=2 (d=1)
    class2_results = {}
    for d, label in [(0, "nu1"), (1, "nu2"), (3, "nu4")]:
        D, B_of_z, c, p, dim = build_class2_example(d)
        nu_moment, mu = nu_from_moments(D, B_of_z, c, p, z=0.3 + 0.1j, kmax=dim + 1)
        vals = []
        for g in GAMMAS_FIG1:
            A = g * D + B_of_z(0.3 + 0.1j)
            vals.append(transfer(A, p, c))
        vals = np.array(vals)
        fit = fit_nu_loglog(GAMMAS_FIG1, vals)
        class2_results[label] = {
            "d": d,
            "nu_moment_method": nu_moment,
            "nu_direct_fit": fit["nu_global"],
            "values": vals,
        }
        gate1[f"class_II_{label}"] = {
            "nu_moment_method": nu_moment,
            "nu_direct_fit": fit["nu_global"],
            "match": abs(nu_moment - fit["nu_global"]) < 0.05,
        }

    # -- Class III: protected O(1) channel
    D3, B_of_z3, c3, p3, S, k = build_class3_example()
    B0 = B_of_z3(0.0)
    mask = np.array([True] * k + [False] * (D3.shape[0] - k))
    # work in the S-diagonalizing basis where D is exactly diagonal
    Dinv_S = np.linalg.inv(S)
    D_diag = Dinv_S @ D3 @ S
    B_diag = Dinv_S @ B0 @ S
    c_diag = Dinv_S @ c3
    p_diag = S.T @ p3  # p^dag P == (S^-T p)^dag restricted, using p^dag = p^T (real p)
    F0 = protected_F0(B_diag, mask, c_diag, p_diag)
    vals_III = []
    for g in GAMMAS_FIG1:
        A = g * D3 + B0
        vals_III.append(transfer(A, p3, c3))
    vals_III = np.array(vals_III)
    plateau_err = abs(vals_III[-1] - F0) / abs(F0)
    gate1["class_III"] = {
        "F0_protected_coefficient": [F0.real, F0.imag],
        "F_at_Gamma_max": [vals_III[-1].real, vals_III[-1].imag],
        "relative_plateau_error": float(plateau_err),
        "match": bool(plateau_err < 1e-4),
    }

    # -- plot
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(GAMMAS_FIG1, np.maximum(vals_I, 1e-16), "o-", label="Class I ($\\nu=\\infty$)")
    ax.loglog(GAMMAS_FIG1, np.abs(class2_results["nu1"]["values"]), "s-", label="Class II, $\\nu=1$")
    ax.loglog(GAMMAS_FIG1, np.abs(class2_results["nu2"]["values"]), "^-", label="Class II, $\\nu=2$")
    ax.loglog(GAMMAS_FIG1, np.abs(class2_results["nu4"]["values"]), "d-", label="Class II, $\\nu=4$")
    ax.loglog(GAMMAS_FIG1, np.abs(vals_III), "v-", label="Class III ($\\nu=0$)")
    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"$|F_\Gamma|$ (or $|R_{S,\Gamma}|$)")
    ax.set_title("Figure 1: theory unit test")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig1_theory_unit_test.png"), dpi=150)
    plt.close(fig)

    return gate1


# ----------------------------------------------------------------------
# Priority 3: hidden class transition + Figure 3
# ----------------------------------------------------------------------
def priority3_transition():
    lam_c = brentq(r0_of_lambda, 0.5, 0.9, xtol=1e-13)

    lams = np.linspace(0.0, 1.4, 141)
    r0_vals = np.array([r0_of_lambda(l) for l in lams])

    Gamma_probe = (1e3, 1e5)  # bracketing pair for local nu_eff(lambda)
    chi_full_vals = []
    R_S_vals_lo = []
    R_S_vals_hi = []
    for l in lams:
        cf_lo, cc_lo, R_lo = R_S_gamma(Gamma_probe[0], l)
        cf_hi, cc_hi, R_hi = R_S_gamma(Gamma_probe[1], l)
        chi_full_vals.append(cf_hi)
        R_S_vals_lo.append(R_lo)
        R_S_vals_hi.append(R_hi)
    chi_full_vals = np.array(chi_full_vals)
    R_S_vals_lo = np.array(R_S_vals_lo)
    R_S_vals_hi = np.array(R_S_vals_hi)

    nu_eff_lambda = -(np.log(np.abs(R_S_vals_hi)) - np.log(np.abs(R_S_vals_lo))) / \
        (np.log(Gamma_probe[1]) - np.log(Gamma_probe[0]))

    fig, axes = plt.subplots(4, 1, figsize=(6.5, 9), sharex=True)
    axes[0].plot(lams, nu_eff_lambda)
    axes[0].axvline(lam_c, color="k", ls="--", lw=1)
    axes[0].set_ylabel(r"$\nu_{\rm eff}(\lambda)$")
    axes[1].semilogy(lams, np.abs(r0_vals))
    axes[1].axvline(lam_c, color="k", ls="--", lw=1)
    axes[1].set_ylabel(r"$|r_0(\lambda)|$")
    axes[2].plot(lams, np.abs(chi_full_vals))
    axes[2].axvline(lam_c, color="k", ls="--", lw=1)
    axes[2].set_ylabel(r"$|\chi_{\rm full}|$")
    axes[3].semilogy(lams, np.abs(R_S_vals_hi))
    axes[3].axvline(lam_c, color="k", ls="--", lw=1)
    axes[3].set_ylabel(r"$|R_S|$ ($\Gamma=10^5$)")
    axes[3].set_xlabel(r"$\lambda$")
    fig.suptitle(f"Figure 3: hidden class transition ($\\lambda_c={lam_c:.4f}$)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig3_hidden_class_transition.png"), dpi=150)
    plt.close(fig)

    gate2 = {
        "lambda_c": float(lam_c),
        "r0_at_lambda_c": [r0_of_lambda(lam_c).real if hasattr(r0_of_lambda(lam_c), "real") else float(r0_of_lambda(lam_c)), 0.0],
    }

    gammas = np.logspace(2, 6, 15)
    for lam, tag in [(lam_c, "at_lambda_c"), (lam_c + 0.3, "off_lambda_c")]:
        Rs, Cfs = [], []
        for g in gammas:
            cf, cc, R = R_S_gamma(g, lam)
            Rs.append(R)
            Cfs.append(cf)
        fitR = fit_nu_loglog(gammas, Rs)
        fitF = fit_nu_loglog(gammas, Cfs)
        gate2[tag] = {
            "lambda": float(lam),
            "nu_R_S_direct_fit": fitR["nu_global"],
            "nu_chi_full_direct_fit": fitF["nu_global"],
            "chi_full_plateau": abs(Cfs[-1]),
        }

    gate2["passes"] = bool(
        abs(gate2["off_lambda_c"]["nu_chi_full_direct_fit"]) < 0.05
        and abs(gate2["at_lambda_c"]["nu_chi_full_direct_fit"]) < 0.05
        and gate2["off_lambda_c"]["nu_R_S_direct_fit"] < 0.1
        and abs(gate2["at_lambda_c"]["nu_R_S_direct_fit"] - 1.0) < 0.1
    )

    return lam_c, gate2


# ----------------------------------------------------------------------
# Priority 4: scaling collapse + Figure 4
# ----------------------------------------------------------------------
def priority4_collapse(lam_c):
    gammas_collapse = [1e3, 3e3, 1e4, 3e4, 1e5]
    window = 3.0  # X = Gamma*(lambda - lambda_c) window half-width

    fig, ax = plt.subplots(figsize=(6, 4.5))
    curves = {}
    for g in gammas_collapse:
        half_width_lambda = window / g
        lams_loc = lam_c + np.linspace(-half_width_lambda, half_width_lambda, 121)
        Y = []
        for l in lams_loc:
            _, _, R = R_S_gamma(g, l)
            Y.append(g * R)
        X = g * (lams_loc - lam_c)
        Y = np.array(Y)
        ax.plot(X, Y.real, label=f"$\\Gamma={g:.0e}$")
        curves[str(g)] = {"X": X, "Y_real": Y.real}

    ax.set_xlabel(r"$X=\Gamma(\lambda-\lambda_c)$")
    ax.set_ylabel(r"$Y=\Gamma\,\mathrm{Re}\,R_{S,\Gamma}$")
    ax.set_title("Figure 4: scaling collapse")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig4_scaling_collapse.png"), dpi=150)
    plt.close(fig)

    # collapse quality: RMS spread of Y across curves on a common X grid,
    # relative to the RMS amplitude of Y, checked for two different pairs
    # of Gamma spanning >= 1 decade
    common_X = np.linspace(-window * 0.8, window * 0.8, 81)
    Y_interp = []
    for g in gammas_collapse:
        X = curves[str(g)]["X"]
        Y = curves[str(g)]["Y_real"]
        Y_interp.append(np.interp(common_X, X, Y))
    Y_interp = np.array(Y_interp)
    spread = np.std(Y_interp, axis=0)
    amplitude = np.mean(np.abs(Y_interp), axis=0) + 1e-12
    rel_spread = float(np.mean(spread / amplitude))

    gate4 = {
        "gammas_used": gammas_collapse,
        "relative_collapse_spread": rel_spread,
        "passes": bool(rel_spread < 0.05),
    }
    return gate4


def main():
    gate1 = figure1_and_gate1()
    lam_c, gate2 = priority3_transition()
    gate4 = priority4_collapse(lam_c)

    summary = {"gate1_classifier_consistency": gate1,
               "gate2_hidden_transition": gate2,
               "gate4_scaling_collapse": gate4}

    def _default(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        if isinstance(o, complex):
            return [o.real, o.imag]
        raise TypeError(str(type(o)))

    with open(os.path.join(RESULTS, "gates_summary_phaseB.json"), "w") as f:
        json.dump(summary, f, indent=2, default=_default)

    print("Gate 1 (classifier consistency):")
    for k, v in gate1.items():
        print(" ", k, "->", {kk: vv for kk, vv in v.items() if kk != "values"})
    print("\nGate 2 (hidden class transition):")
    print(" lambda_c =", lam_c, " passes =", gate2["passes"])
    print("\nGate 4 (scaling collapse):")
    print(" relative_collapse_spread =", gate4["relative_collapse_spread"], " passes =", gate4["passes"])


if __name__ == "__main__":
    main()
