"""Phase H: physical GKSL realization of the scalar hidden-class transition.

This phase deliberately stays inside the SMRT claim boundary: one sector,
one fixed source, one fixed scalar readout, and one full-cut master response.
It does not scan observables, compare functional norms, or construct
irreducible multi-sector signals reserved for RISEI.
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from core import fit_nu_loglog
import model_hidden_physical as M


HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
FIGDIR = os.path.join(RESULTS, "figures")
os.makedirs(FIGDIR, exist_ok=True)


def _pair(z):
    return [float(np.real(z)), float(np.imag(z))]


def finite_jacobian(fun, x, steps=(1e-6, 1e-6)):
    cols = []
    for k, h in enumerate(steps):
        xp = np.array(x, dtype=float)
        xm = np.array(x, dtype=float)
        xp[k] += h
        xm[k] -= h
        cols.append((np.asarray(fun(xp)) - np.asarray(fun(xm))) / (2 * h))
    return np.column_stack(cols)


def gate_H0_physicality():
    jstar = M.j12_star()
    q = M.params()
    D = M.reduced_D()
    rates = [q["gamma1"], q["gamma2"], q["gamma3_fixed"], q["d3_fast"]]
    residuals = []
    hermiticity = []
    for cut in [False, True]:
        H, _, L = M.full_components(17.0, cut=cut)
        hermiticity.append(float(np.linalg.norm(H - H.conj().T)))
        residuals.append(M.trace_preservation_residual(L))
    gate = {
        "claim_boundary": {
            "single_sector": True,
            "fixed_source": True,
            "fixed_readout_rho10": True,
            "observable_scan": False,
            "multi_sector_composition": False,
        },
        "J12_star": _pair(jstar),
        "J12_star_abs": float(abs(jstar)),
        "rates": rates,
        "all_rates_nonnegative": bool(min(rates) >= 0),
        "reduced_D_rank": int(np.linalg.matrix_rank(D)),
        "reduced_D_kernel_dimension": int(D.shape[0] - np.linalg.matrix_rank(D)),
        "hamiltonian_hermiticity_residual_max": max(hermiticity),
        "trace_preservation_residual_max": max(residuals),
    }
    gate["overall_pass"] = bool(
        gate["all_rates_nonnegative"]
        and gate["reduced_D_rank"] == 1
        and gate["reduced_D_kernel_dimension"] == 2
        and gate["hamiltonian_hermiticity_residual_max"] < 1e-12
        and gate["trace_preservation_residual_max"] < 1e-12
    )
    return gate


def gate_H1_analytic_hidden_point():
    jstar = M.j12_star()
    gstar, phistar = abs(jstar), np.angle(jstar)

    def r0_controls(v):
        j = v[0] * np.exp(1j * v[1])
        r = M.protected_r0(j)
        return np.array([r.real, r.imag])

    jac = finite_jacobian(r0_controls, [gstar, phistar])
    svals = np.linalg.svd(jac, compute_uv=False)
    r0 = M.protected_r0(jstar)

    # The mechanism is hidden only in the fixed rho10 readout: the full and
    # cut protected states remain different internally at the cancellation.
    cP = M.source_vector()[:2]
    Bf = M.reduced_B(jstar)[:2, :2]
    Bc = M.reduced_B(jstar, cut=True)[:2, :2]
    xf = np.linalg.solve(Bf, cP)
    xc = np.linalg.solve(Bc, cP)
    min_sv = min(
        np.linalg.svd(Bf, compute_uv=False)[-1],
        np.linalg.svd(Bc, compute_uv=False)[-1],
    )
    gate = {
        "analytic_condition": "J12_star = conjugate(-i*a*c2/c1)",
        "J12_star": _pair(jstar),
        "amplitude": float(gstar),
        "phase": float(phistar),
        "r0_at_star": _pair(r0),
        "abs_r0_at_star": float(abs(r0)),
        "transverse_jacobian": jac.tolist(),
        "jacobian_singular_values": svals.tolist(),
        "jacobian_rank": int(np.linalg.matrix_rank(jac, tol=1e-8)),
        "protected_internal_state_difference": float(np.linalg.norm(xf - xc)),
        "protected_block_min_singular_value": float(min_sv),
        "chi_full_cut_common_plateau": _pair(xf[0]),
    }
    gate["overall_pass"] = bool(
        abs(r0) < 1e-12
        and gate["jacobian_rank"] == 2
        and gate["protected_internal_state_difference"] > 1e-2
        and min_sv > 0.1
    )
    return gate


def gate_H2_exponents():
    gammas = np.logspace(2, 7, 31)
    cases = {"at_hidden_point": 0.0, "off_hidden_point": 0.2}
    table = {}
    curves = {}
    for name, lam in cases.items():
        full, cut, response = [], [], []
        for g in gammas:
            cf, cc, rs, _, _ = M.reduced_responses(g, lambda_amp=lam)
            full.append(cf)
            cut.append(cc)
            response.append(rs)
        full = np.asarray(full)
        cut = np.asarray(cut)
        response = np.asarray(response)
        tail = gammas >= 1e4
        table[name] = {
            "lambda_amp": lam,
            "nu_R_S": fit_nu_loglog(gammas[tail], response[tail])["nu_global"],
            "nu_chi_full": fit_nu_loglog(gammas[tail], full[tail])["nu_global"],
            "nu_chi_cut": fit_nu_loglog(gammas[tail], cut[tail])["nu_global"],
            "chi_full_at_Gamma_max": _pair(full[-1]),
            "R_S_at_Gamma_max": _pair(response[-1]),
            "Gamma_R_S_at_Gamma_max": _pair(gammas[-1] * response[-1]),
        }
        curves[name] = (full, response)

    at, off = table["at_hidden_point"], table["off_hidden_point"]
    gate = {"gammas": gammas.tolist(), "table": table}
    gate["overall_pass"] = bool(
        abs(at["nu_R_S"] - 1.0) < 0.02
        and abs(off["nu_R_S"]) < 0.02
        and abs(at["nu_chi_full"]) < 0.02
        and abs(off["nu_chi_full"]) < 0.02
        and abs(complex(*at["Gamma_R_S_at_Gamma_max"])) > 1e-3
    )
    return gate, gammas, curves


def gate_H3_three_control_manifold():
    """Continue the codimension-2 zero through the third control Delta2."""
    jstar = M.j12_star()
    points = []
    exponent_checks = []
    for delta2 in np.linspace(-0.3, 0.3, 13):
        r0 = M.protected_r0(jstar, delta2_control=delta2)
        Bf = M.reduced_B(jstar, delta2_control=delta2)[:2, :2]
        Bc = M.reduced_B(jstar, delta2_control=delta2, cut=True)[:2, :2]
        points.append({
            "delta2_control": float(delta2),
            "J12_abs": float(abs(jstar)),
            "J12_phase": float(np.angle(jstar)),
            "abs_r0": float(abs(r0)),
            "min_block_singular_value": float(min(
                np.linalg.svd(Bf, compute_uv=False)[-1],
                np.linalg.svd(Bc, compute_uv=False)[-1],
            )),
        })
    for delta2 in [-0.3, 0.0, 0.3]:
        gammas = np.logspace(3, 7, 21)
        vals = [M.reduced_responses(g, delta2_control=delta2)[2] for g in gammas]
        exponent_checks.append({
            "delta2_control": delta2,
            "nu_R_S": fit_nu_loglog(gammas[-13:], np.asarray(vals)[-13:])["nu_global"],
        })
    max_r0 = max(p["abs_r0"] for p in points)
    min_sv = min(p["min_block_singular_value"] for p in points)
    gate = {
        "controls": ["abs(J12)", "arg(J12)", "Delta2"],
        "zero_manifold_dimension": 1,
        "points": points,
        "max_abs_r0_on_line": max_r0,
        "min_protected_block_singular_value_on_line": min_sv,
        "exponent_checks": exponent_checks,
    }
    gate["overall_pass"] = bool(
        max_r0 < 1e-12
        and min_sv > 0.1
        and all(abs(d["nu_R_S"] - 1.0) < 0.02 for d in exponent_checks)
    )
    return gate


def gate_H4_scaling_collapse():
    gammas = np.array([1e3, 3e3, 1e4, 3e4, 1e5])
    common_X = np.linspace(-3.0, 3.0, 121)
    curves = []
    for g in gammas:
        vals = []
        for x in common_X:
            lam = x / g
            vals.append(g * M.reduced_responses(g, lambda_amp=lam)[2])
        curves.append(np.asarray(vals))
    curves = np.asarray(curves)
    mean = np.mean(curves, axis=0)
    spread = np.sqrt(np.mean(np.abs(curves - mean) ** 2, axis=0))
    relative = spread / (np.abs(mean) + 1e-12)
    gate = {
        "gammas": gammas.tolist(),
        "X_window": [float(common_X[0]), float(common_X[-1])],
        "mean_relative_spread": float(np.mean(relative)),
        "max_relative_spread": float(np.max(relative)),
    }
    gate["overall_pass"] = bool(gate["max_relative_spread"] < 0.01)
    return gate, common_X, curves


def gate_H5_full_liouvillian():
    gammas = [10.0, 100.0, 1e3, 1e4, 1e5]
    cases = {"at_hidden_point": 0.0, "off_hidden_point": 0.2}
    rows = []
    max_component_error = 0.0
    max_response_error = 0.0
    for name, lam in cases.items():
        for g in gammas:
            cf, cc, rs, xf, xc = M.reduced_responses(g, lambda_amp=lam)
            cff, xff, _, _, _, _ = M.full_liouvillian_response(g, lambda_amp=lam, cut=False)
            ccf, xcf, _, _, _, _ = M.full_liouvillian_response(g, lambda_amp=lam, cut=True)
            rsf = cff - ccf
            component_error = max(np.max(np.abs(xf - xff)), np.max(np.abs(xc - xcf)))
            response_error = abs(rs - rsf)
            max_component_error = max(max_component_error, component_error)
            max_response_error = max(max_response_error, response_error)
            rows.append({
                "case": name,
                "Gamma": g,
                "component_abs_error": float(component_error),
                "R_S_abs_error": float(response_error),
            })
    gate = {
        "reduced_dimension": 3,
        "full_liouvillian_dimension": 16,
        "rows": rows,
        "max_component_abs_error": float(max_component_error),
        "max_R_S_abs_error": float(max_response_error),
    }
    gate["overall_pass"] = bool(
        max_component_error < 1e-12 and max_response_error < 1e-12
    )
    return gate


def make_figure(gammas, exponent_curves, X, collapse_curves):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for name, (full, response) in exponent_curves.items():
        label = "hidden point" if name == "at_hidden_point" else "off point"
        axes[0, 0].loglog(gammas, np.abs(response), label=label)
        axes[0, 1].semilogx(gammas, np.abs(full), label=label)
    axes[0, 0].set_xlabel(r"$\Gamma$")
    axes[0, 0].set_ylabel(r"$|R_S|$")
    axes[0, 0].set_title(r"Hidden transition: $\nu[R_S]:0\to1$")
    axes[0, 0].legend()
    axes[0, 1].set_xlabel(r"$\Gamma$")
    axes[0, 1].set_ylabel(r"$|\chi_{\rm full}|$")
    axes[0, 1].set_title(r"Full response remains Class III ($\nu=0$)")
    axes[0, 1].legend()
    for g, curve in zip([1e3, 3e3, 1e4, 3e4, 1e5], collapse_curves):
        axes[1, 0].plot(X, curve.real, label=f"{g:.0e}")
        axes[1, 1].plot(X, curve.imag, label=f"{g:.0e}")
    axes[1, 0].set_xlabel(r"$X=\Gamma\lambda$")
    axes[1, 0].set_ylabel(r"$\mathrm{Re}(\Gamma R_S)$")
    axes[1, 0].set_title("Universal collapse (real)")
    axes[1, 1].set_xlabel(r"$X=\Gamma\lambda$")
    axes[1, 1].set_ylabel(r"$\mathrm{Im}(\Gamma R_S)$")
    axes[1, 1].set_title("Universal collapse (imaginary)")
    axes[1, 1].legend(title=r"$\Gamma$", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figH1_physical_hidden_transition.png"), dpi=170)
    plt.close(fig)


def main():
    gates = {}
    gates["H0_physicality_and_claim_boundary"] = gate_H0_physicality()
    gates["H1_analytic_hidden_point"] = gate_H1_analytic_hidden_point()
    gates["H2_hidden_transition_exponents"], gammas, exponent_curves = gate_H2_exponents()
    gates["H3_three_control_zero_manifold"] = gate_H3_three_control_manifold()
    gates["H4_scaling_collapse"], X, collapse_curves = gate_H4_scaling_collapse()
    gates["H5_reduced_vs_full_liouvillian"] = gate_H5_full_liouvillian()
    gates["phase_h_overall_pass"] = bool(all(
        gate["overall_pass"] for gate in gates.values()
        if isinstance(gate, dict) and "overall_pass" in gate
    ))
    make_figure(gammas, exponent_curves, X, collapse_curves)
    out = os.path.join(RESULTS, "gates_summary_phaseH.json")
    with open(out, "w") as f:
        json.dump(gates, f, indent=2)
    for name, gate in gates.items():
        if isinstance(gate, dict) and "overall_pass" in gate:
            print(f"{name}: {gate['overall_pass']}")
    print(f"phase_h_overall_pass: {gates['phase_h_overall_pass']}")
    print(out)


if __name__ == "__main__":
    main()
