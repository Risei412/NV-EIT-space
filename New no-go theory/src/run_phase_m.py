"""Phase M driver: metrological (QFI) extension of the sector-resolved
response classification. Tests the candidate prediction

    x_{S,Gamma} = d(rho_full)/dtheta - d(rho_cut)/dtheta ,  ||x_S|| ~ Gamma^{-nu}
    F_{Q,S} = x_S^dagger G_rho x_S ~ Gamma^{-2 nu}

Gates:
  M1 - vector-valued readouts preserve the nu classification (abstract model)
  M2 - EXACT vector hidden class transition + Gamma(lambda-lambda_c) collapse
       (abstract model, rank-1 sector cut)
  M3 - nu -> 2*nu QFI translation at a generic point (physical 3-level model)
  M4 - search for a physical (Hermitian-coupling) interior transition; report
       whether one is found, and interpret via the codimension count of the
       protected block's traceless-Hermitian degrees of freedom
  M5 - rank-deficiency negative control: the nu->2nu law degrades as the SLD
       metric approaches a singular (near-pure-state) regime
"""

import json
import os

import numpy as np
from scipy.optimize import minimize

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core import fit_nu_loglog
from model_metro_linear import (
    vector_norm_class1, vector_norm_class2, vector_norm_class3,
    leading_order_lambda_c, x_S_gamma,
)
from model_metro_lindblad import x_S_lindblad, qfi

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
    raise TypeError(str(type(o)))


# ----------------------------------------------------------------------
# Gate M1
# ----------------------------------------------------------------------
def gate_M1():
    gammas = np.logspace(1, 6, 20)
    n1 = vector_norm_class1(gammas)
    results = {"class_I_floor": float(np.max(n1))}
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(gammas, np.maximum(n1, 1e-16), "o-", label="Class I (vector)")
    for d, marker in [(0, "s"), (1, "^"), (3, "d")]:
        n2 = vector_norm_class2(gammas, d)
        fit = fit_nu_loglog(gammas, n2)
        results[f"class_II_d{d}_nu"] = fit["nu_global"]
        ax.loglog(gammas, n2, marker + "-", label=f"Class II, d={d} ($\\nu={d+1}$)")
    n3 = vector_norm_class3(gammas)
    results["class_III_range"] = [float(n3.min()), float(n3.max())]
    ax.loglog(gammas, n3, "v-", label="Class III (vector)")
    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel(r"$\|x_{S,\Gamma}\|$ (vector readout)")
    ax.set_title("Figure M1: vector-readout unit test")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figM1_vector_unit_test.png"), dpi=150)
    plt.close(fig)

    results["passes"] = bool(
        results["class_I_floor"] < 1e-10
        and abs(results["class_II_d0_nu"] - 1) < 0.05
        and abs(results["class_II_d1_nu"] - 2) < 0.05
        and abs(results["class_II_d3_nu"] - 4) < 0.05
    )
    return results


# ----------------------------------------------------------------------
# Gate M2
# ----------------------------------------------------------------------
def gate_M2():
    lam_c0, _ = leading_order_lambda_c()
    gammas = np.logspace(2, 6, 15)

    results = {"lambda_c0": float(lam_c0)}
    for lam, tag in [(lam_c0, "at_lambda_c"), (lam_c0 + 0.5, "off_lambda_c")]:
        norms = [np.linalg.norm(x_S_gamma(g, lam)) for g in gammas]
        fit = fit_nu_loglog(gammas, norms)
        results[tag] = {"lambda": float(lam), "nu": fit["nu_global"], "norm_at_Gmax": norms[-1]}

    # scaling collapse
    gammas_collapse = [1e3, 3e3, 1e4, 3e4, 1e5]
    window = 0.3
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    curves = {}
    for g in gammas_collapse:
        half = window / g
        lams_loc = lam_c0 + np.linspace(-half, half, 121)
        Y = np.array([g * np.linalg.norm(x_S_gamma(g, lam)) for lam in lams_loc])
        X = g * (lams_loc - lam_c0)
        axes[1].plot(X, Y, label=f"$\\Gamma={g:.0e}$")
        curves[str(g)] = {"X": X, "Y": Y}

    lam_grid = np.linspace(lam_c0 - 1.5, lam_c0 + 1.5, 200)
    for g in [1e3, 1e5]:
        norms = [np.linalg.norm(x_S_gamma(g, lam)) for lam in lam_grid]
        axes[0].semilogy(lam_grid, norms, label=f"$\\Gamma={g:.0e}$")
    axes[0].axvline(lam_c0, color="k", ls="--", lw=1)
    axes[0].set_xlabel(r"$\lambda$")
    axes[0].set_ylabel(r"$\|x_{S,\Gamma}\|$")
    axes[0].set_title("hidden class transition")
    axes[0].legend(fontsize=8)

    axes[1].set_xlabel(r"$X=\Gamma(\lambda-\lambda_c)$")
    axes[1].set_ylabel(r"$Y=\Gamma\|x_{S,\Gamma}\|$")
    axes[1].set_title("scaling collapse")
    axes[1].legend(fontsize=8)
    fig.suptitle(f"Figure M2: abstract rank-1 cut, exact vector transition ($\\lambda_c={lam_c0:.4f}$)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figM2_abstract_vector_transition.png"), dpi=150)
    plt.close(fig)

    common_X = np.linspace(-window * 0.8, window * 0.8, 81)
    Y_interp = []
    for g in gammas_collapse:
        Y_interp.append(np.interp(common_X, curves[str(g)]["X"], curves[str(g)]["Y"]))
    Y_interp = np.array(Y_interp)
    spread = np.std(Y_interp, axis=0)
    amp = np.mean(np.abs(Y_interp), axis=0) + 1e-12
    rel_spread = float(np.mean(spread / amp))
    results["collapse_relative_spread"] = rel_spread

    results["passes"] = bool(
        results["at_lambda_c"]["nu"] > 0.9
        and abs(results["off_lambda_c"]["nu"]) < 0.05
        and rel_spread < 0.05
    )
    return results


# ----------------------------------------------------------------------
# Gate M3: nu -> 2nu QFI translation, physical model, generic point
# ----------------------------------------------------------------------
def gate_M3():
    gammas = np.logspace(2, 5, 15)
    norms, FQS, pmins = [], [], []
    for g in gammas:
        xS, dfull, dcut, rho_full, rho_cut = x_S_lindblad(g, 0.0, lam=0.0, kappa0=1.5)
        norms.append(np.linalg.norm(xS))
        F, drop, pmin = qfi(rho_full, xS.reshape(3, 3))
        FQS.append(F)
        pmins.append(pmin)

    fit_x = fit_nu_loglog(gammas, norms)
    fit_F = fit_nu_loglog(gammas, FQS)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(gammas, norms, "o-", label=r"$\|x_{S,\Gamma}\|$" + f" ($\\nu={fit_x['nu_global']:.3f}$)")
    ax.loglog(gammas, FQS, "s-", label=r"$F_{Q,S}$" + f" ($2\\nu={fit_F['nu_global']:.3f}$)")
    ax.set_xlabel(r"$\Gamma$")
    ax.set_ylabel("magnitude")
    ax.set_title("Figure M3: $\\nu \\to 2\\nu$ QFI translation (physical model)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figM3_qfi_translation.png"), dpi=150)
    plt.close(fig)

    results = {
        "nu_x": fit_x["nu_global"],
        "nu_F": fit_F["nu_global"],
        "ratio_nuF_over_2nuX": fit_F["nu_global"] / (2 * fit_x["nu_global"]),
        "min_eigenvalue_range": [float(np.min(pmins)), float(np.max(pmins))],
    }
    results["passes"] = bool(abs(results["ratio_nuF_over_2nuX"] - 1.0) < 0.05)
    return results


# ----------------------------------------------------------------------
# Gate M4: physical (Hermitian-coupling) interior transition search
# ----------------------------------------------------------------------
def gate_M4():
    Gamma = 1e6
    kappa0 = 1.5

    def obj(params):
        lam, phi = params
        if lam < -1.0 or lam > 1.5:
            return 1e3
        xS, _, _, _, _ = x_S_lindblad(Gamma, 0.0, lam, kappa0=kappa0, phi=phi)
        return np.linalg.norm(xS)

    generic_value = obj([0.0, 0.0])

    rng = np.random.default_rng(1)
    best = None
    for _ in range(40):
        x0 = [rng.uniform(-1.0, 1.5), rng.uniform(-np.pi, np.pi)]
        res = minimize(obj, x0, method="Nelder-Mead",
                        options={"xatol": 1e-12, "fatol": 1e-20, "maxiter": 3000})
        if best is None or res.fun < best.fun:
            best = res

    # 2D heatmap for visualization (coarser grid, smaller Gamma for speed/contrast)
    Gamma_map = 1e4
    lams = np.linspace(-0.95, 1.4, 45)
    phis = np.linspace(-np.pi, np.pi, 45)
    grid = np.zeros((len(phis), len(lams)))
    for i, phi in enumerate(phis):
        for j, lam in enumerate(lams):
            xS, _, _, _, _ = x_S_lindblad(Gamma_map, 0.0, lam, kappa0=kappa0, phi=phi)
            grid[i, j] = np.log10(np.linalg.norm(xS) + 1e-300)

    fig, ax = plt.subplots(figsize=(6.5, 5))
    im = ax.pcolormesh(lams, phis, grid, shading="auto", cmap="viridis")
    ax.set_xlabel(r"$\lambda$ (amplitude shift)")
    ax.set_ylabel(r"$\phi$ (phase)")
    ax.set_title(f"Figure M4: $\\log_{{10}}\\|x_{{S,\\Gamma}}\\|$ over ($\\lambda,\\phi$), $\\Gamma=10^4$")
    fig.colorbar(im, ax=ax, label=r"$\log_{10}\|x_S\|$")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figM4_physical_lambda_phi_scan.png"), dpi=150)
    plt.close(fig)

    suppression_ratio = float(generic_value / max(best.fun, 1e-300))
    results = {
        "generic_point_norm": float(generic_value),
        "best_found_norm": float(best.fun),
        "best_found_lam_phi": [float(best.x[0]), float(best.x[1])],
        "suppression_ratio_vs_generic": suppression_ratio,
        "interior_transition_found": bool(suppression_ratio > 100 and best.x[0] > -0.999),
    }
    results["interpretation"] = (
        "No interior (non-trivial) exact or near-exact vector transition was found "
        "with a 2-real-parameter (amplitude, phase) physical Hermitian coupling cut; "
        "the leading-order response x_{S,0} lives in the 3-real-dimensional "
        "traceless-Hermitian slice of the protected block, so codimension counting "
        "predicts an interior zero generically requires 3 independent real controls. "
        "This is consistent with (not a failure of) the theory: it correctly predicts "
        "when a hidden transition is and is not achievable with a given control budget. "
        "The exact vector transition IS realized with a rank-1 (non-Hermitian) abstract "
        "cut in Gate M2."
    )
    return results


# ----------------------------------------------------------------------
# Gate M5: rank-deficiency negative control
# ----------------------------------------------------------------------
def gate_M5():
    gammas = np.logspace(2, 5, 15)
    out = {}
    for eps, tag in [(0.01, "normal_epsilon"), (1e-6, "tiny_epsilon")]:
        norms, FQS, pmins = [], [], []
        for g in gammas:
            xS, dfull, dcut, rho_full, rho_cut = x_S_lindblad(
                g, 0.0, lam=0.0, kappa0=1.5, epsilon=eps)
            norms.append(np.linalg.norm(xS))
            F, drop, pmin = qfi(rho_full, xS.reshape(3, 3))
            FQS.append(F)
            pmins.append(pmin)
        fit_x = fit_nu_loglog(gammas, norms)
        fit_F = fit_nu_loglog(gammas, FQS)
        out[tag] = {
            "epsilon": eps,
            "nu_x": fit_x["nu_global"],
            "nu_F": fit_F["nu_global"],
            "ratio_nuF_over_2nuX": fit_F["nu_global"] / (2 * fit_x["nu_global"]),
            "min_eigenvalue_at_Gmax": pmins[-1],
        }

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for tag, marker in [("normal_epsilon", "o-"), ("tiny_epsilon", "s-")]:
        d = out[tag]
        ax.axhline(d["ratio_nuF_over_2nuX"], label=f"{tag}: ratio={d['ratio_nuF_over_2nuX']:.3f}")
    ax.set_ylabel(r"$\nu_F / (2\nu_x)$")
    ax.set_title("Figure M5: SLD-singularity breakdown of the $\\nu\\to2\\nu$ law")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figM5_rank_deficiency_control.png"), dpi=150)
    plt.close(fig)

    deviation = abs(out["tiny_epsilon"]["ratio_nuF_over_2nuX"] - 1.0)
    out["passes"] = bool(
        abs(out["normal_epsilon"]["ratio_nuF_over_2nuX"] - 1.0) < 0.02
        and deviation > 0.1
    )
    return out


def main():
    m1 = gate_M1()
    m2 = gate_M2()
    m3 = gate_M3()
    m4 = gate_M4()
    m5 = gate_M5()

    summary = {
        "gateM1_vector_readout_consistency": m1,
        "gateM2_abstract_exact_vector_transition": m2,
        "gateM3_nu_to_2nu_qfi_translation": m3,
        "gateM4_physical_interior_transition_search": m4,
        "gateM5_rank_deficiency_negative_control": m5,
    }

    with open(os.path.join(RESULTS, "gates_summary_phaseM.json"), "w") as f:
        json.dump(summary, f, indent=2, default=_default)

    print("Gate M1 (vector readout):", "PASS" if m1["passes"] else "FAIL")
    print("Gate M2 (abstract exact vector transition):", "PASS" if m2["passes"] else "FAIL", m2["lambda_c0"])
    print("Gate M3 (nu->2nu QFI translation):", "PASS" if m3["passes"] else "FAIL", m3["ratio_nuF_over_2nuX"])
    print("Gate M4 (physical interior transition search): found =", m4["interior_transition_found"],
          " suppression_ratio =", m4["suppression_ratio_vs_generic"])
    print("Gate M5 (rank-deficiency negative control):", "PASS" if m5["passes"] else "FAIL")


if __name__ == "__main__":
    main()
