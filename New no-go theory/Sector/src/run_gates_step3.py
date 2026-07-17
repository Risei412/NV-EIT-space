"""Runs Gates U1, U2, U5 (Step 3 of the operational-cut revision plan) on
the minimal Lambda / ground-coherence models, and writes a JSON summary to
Sector/results/gates_summary_sector.json.
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import operational_cut as oc
import model_lambda_operational as mlo
import model_ground_coherence_lindblad as mgc

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def gate_u1():
    """Algebraic-operational equivalence: kappa -> infinity convergence of
    chi_op(kappa) to chi_cut^(S), rate O(kappa^-1), coefficient match to
    the closed form of the Corollary."""
    Delta, delta, Omega_c, gamma31, gamma21 = 2.0, -1.5, 3.0 + 0.5j, 1.0, 0.05
    M = mlo.M_matrix(Delta, delta, Omega_c, gamma31, gamma21)
    D_S = mlo.D_S()
    c, p = mlo.source_readout()

    v = oc.validate_cut_generator(D_S)
    assert v["semisimple"], "D_S zero eigenvalue not semisimple"
    assert v["is_hermitian"], "D_S expected Hermitian for Lemma 0B'"

    # Fit only the asymptotic tail (large kappa): at small kappa the
    # O(kappa^-2) term contaminates the log-log slope/coefficient fit.
    kappas = np.logspace(3, 7, 20)
    res = oc.compare_cut_equivalence(M, D_S, c, p, kappas)

    chi_alg = 1.0 / (gamma31 - 1j * Delta)  # Theorem 2B algebraic cut
    coeff_closed = mlo.closed_form_kappa_coefficient(Delta, Omega_c, gamma31)

    ideal_vs_alg_err = abs(res["chi_ideal"] - chi_alg)
    slope = res["fit_slope"]
    fitted_coeff = res["fitted_coeff_mag"]
    closed_coeff_mag = abs(coeff_closed)
    coeff_rel_err = abs(fitted_coeff - closed_coeff_mag) / closed_coeff_mag

    passed = (
        ideal_vs_alg_err < 1e-9
        and abs(slope - (-1.0)) < 0.02
        and coeff_rel_err < 0.02
    )
    return {
        "name": "U1_algebraic_operational_equivalence",
        "ideal_vs_algebraic_cut_error": ideal_vs_alg_err,
        "fit_slope": slope,
        "fitted_kappa_coefficient_mag": fitted_coeff,
        "closed_form_coefficient_mag": closed_coeff_mag,
        "coefficient_relative_error": coeff_rel_err,
        "pass": bool(passed),
    }


def gate_u2():
    """Kernel universality: two admissible cut generators with the same
    P_S (here, same trivial 1-d kernel by construction, but structurally
    different -- one real rate, one complex rate with positive real part,
    i.e. dephasing plus an extra detuning-like anti-Hermitian addition on
    the sector) agree in the kappa -> infinity limit but differ at finite
    kappa."""
    Delta, delta, Omega_c, gamma31, gamma21 = 2.0, -1.5, 3.0 + 0.5j, 1.0, 0.05
    M = mlo.M_matrix(Delta, delta, Omega_c, gamma31, gamma21)
    c, p = mlo.source_readout()

    D_S1 = np.array([[0.0, 0.0], [0.0, 1.0]], dtype=complex)
    D_S2 = np.array([[0.0, 0.0], [0.0, 3.0 + 0.7j]], dtype=complex)  # Re>0

    v1 = oc.validate_cut_generator(D_S1)
    v2 = oc.validate_cut_generator(D_S2)
    same_PS = np.allclose(v1["P_S"], v2["P_S"], atol=1e-10)

    ideal1 = oc.ideal_cut_response(M, D_S1, c, p)
    ideal2 = oc.ideal_cut_response(M, D_S2, c, p)
    ideal_diff = abs(ideal1 - ideal2)

    kappas = np.logspace(1, 4, 10)
    op1 = np.array([oc.operational_cut_response(M, D_S1, c, p, k) for k in kappas])
    op2 = np.array([oc.operational_cut_response(M, D_S2, c, p, k) for k in kappas])
    finite_diff = np.abs(op1 - op2)

    passed = same_PS and ideal_diff < 1e-9 and finite_diff[0] > 1e-6
    return {
        "name": "U2_kernel_universality",
        "same_P_S": bool(same_PS),
        "ideal_limit_difference": ideal_diff,
        "finite_kappa_difference_at_kappa10": float(finite_diff[0]),
        "finite_kappa_difference_at_kappa_max": float(finite_diff[-1]),
        "pass": bool(passed),
    }


def gate_u5():
    """Frozen vs self-consistent cut: verify C_S(rho0) = 0 (condition C4)
    and that the steady state of L0 + kappa*C_S equals rho0 for every
    kappa (Lemma 0D'), hence R_S^reprep = 0."""
    Gamma, Delta3 = 5.0, 1.3
    L0 = mgc.L0(Gamma, Delta3)
    C_S = mgc.C_S_super()
    rho0 = mgc.rho0(Gamma, Delta3)

    def C_S_apply(rho):
        return mgc.unvec(C_S @ mgc.vec(rho))

    noninvasive, C_S_rho0 = oc.check_noninvasive(C_S_apply, rho0)

    kappas = [0.0, 1.0, 10.0, 100.0, 1e4, 1e6]
    max_dev = 0.0
    for k in kappas:
        L_cut = L0 + k * C_S
        rho_cut = mgc.steady_state(L_cut)
        max_dev = max(max_dev, float(np.linalg.norm(rho_cut - rho0)))

    passed = noninvasive and max_dev < 1e-8
    return {
        "name": "U5_frozen_vs_self_consistent",
        "noninvasive_C_S_rho0_norm": float(np.linalg.norm(C_S_rho0)),
        "steady_state_deviation_across_kappa": max_dev,
        "pass": bool(passed),
    }


if __name__ == "__main__":
    results = {"U1": gate_u1(), "U2": gate_u2(), "U5": gate_u5()}
    for key, r in results.items():
        status = "PASS" if r["pass"] else "FAIL"
        print(f"[{key}] {r['name']}: {status}")
        for k, v in r.items():
            if k in ("name", "pass"):
                continue
            print(f"    {k}: {v}")
    out_path = os.path.join(RESULTS_DIR, "gates_summary_sector.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(f"\nWritten: {out_path}")
