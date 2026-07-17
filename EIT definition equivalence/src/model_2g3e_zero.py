"""
E4: Search for an exact dark-state-free EIT transfer zero in a 2g+3e model.

Generalizes the 2g+2e no-go boundary (Theorem IV, `2g2e_package/docs/`):
at ideal two-photon resonance (gamma_g = 0) the ground-coherence
dephasing beta cancels out of the Schur complement, giving the
amplitude-independent structural numerator

    N_full = S_p * S_c - K_pc * K_cp = det(D^dagger G D)

with D = [d_p, d_c] the (N_e x 2) dipole matrix and G = diag(1/a_1,...,1/a_Ne)
the excited-manifold resolvent. For N_e = 2 this is |det D|^2 * det(G): a
single product that cannot vanish while rank D = 2 and G is regular
(Theorem IV, the 2g+2e no-go boundary).

For N_e = 3, Cauchy-Binet expands det(D^dagger G D) into a SUM over the
three 2-element subsets of {1,2,3} of excited-manifold indices:

    det(D^dagger G D) = sum_{|S|=2} det(D_S)^* * det(G_{S,S}) * det(D_S)
                       = sum_{|S|=2} |det(D_S)|^2 * det(G_{S,S})

(since G is diagonal, G_{S,S} = diag over S, det(G_{S,S}) = prod_{j in S} 1/a_j,
and det(D_S)^* * det(D_S) = |det(D_S)|^2 for each 2x2 minor D_S of D formed by
rows S). This is a sum of terms with fixed non-negative real coefficients
|det(D_S)|^2 but COMPLEX phases 1/a_j(delta), so, unlike the N_e=2 case, the
terms can destructively interfere and produce an exact zero at rank D = 2 (no
optical dark vector) and no closed-form obstruction stops it.

This script searches for such a zero: fix theta (probe/control angle in the
e1,e2 plane, inherited from the 2g+2e baseline) and the third-state coupling
amplitude, then solve the 2 real equations Re(N_full)=Im(N_full)=0 for two
real unknowns (probe detuning delta, the third-state closed-loop phase Phi)
at a chosen third-excited-state detuning Delta3.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, fsolve

# ---- baseline parameters, inherited from the 2g+2e package -----------------
gamma_e = 1.0
Delta_e2 = 8.0          # detuning of e2 (same as "Delta_e" in the 2g+2e model)
theta = np.pi / 4        # probe/control angle spanning e1,e2 (rank-2 baseline)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def dipoles(dp3_amp: float, dc3_amp: float, phi: float):
    """Return d_p, d_c in C^3. dp3 fixed real >= 0 by phase convention;
    dc3 carries the free closed-loop phase Phi relative to d_p3."""
    c, s = np.cos(theta), np.sin(theta)
    d_p = np.array([1.0, 0.0, dp3_amp], dtype=complex)
    d_c = np.array([c, s, dc3_amp * np.exp(1j * phi)], dtype=complex)
    return d_p, d_c


def excited_resolvent(delta, Delta3, Gamma1=0.0, Gamma2=0.0, Gamma3=0.0):
    a1 = gamma_e + Gamma1 - 1j * delta
    a2 = gamma_e + Gamma2 + 1j * Delta_e2 - 1j * delta
    a3 = gamma_e + Gamma3 + 1j * Delta3 - 1j * delta
    return np.array([a1, a2, a3])


def N_full(delta, Delta3, dp3_amp, dc3_amp, phi, **Gammas):
    """Structural numerator det(D^dagger G D) = S_p S_c - K_pc K_cp."""
    d_p, d_c = dipoles(dp3_amp, dc3_amp, phi)
    a = excited_resolvent(delta, Delta3, **Gammas)
    G_diag = 1.0 / a
    S_p = np.sum(np.abs(d_p) ** 2 * G_diag)
    S_c = np.sum(np.abs(d_c) ** 2 * G_diag)
    # G = diag(1/a_j) is complex (non-Hermitian for complex a_j), so
    # K_cp = d_c^dagger G d_p is NOT conj(K_pc) in general -- compute both.
    K_pc = np.sum(np.conj(d_p) * d_c * G_diag)
    K_cp = np.sum(np.conj(d_c) * d_p * G_diag)
    return S_p * S_c - K_pc * K_cp, S_c


def rank_D(dp3_amp, dc3_amp, phi):
    d_p, d_c = dipoles(dp3_amp, dc3_amp, phi)
    D = np.column_stack([d_p, d_c])
    return np.linalg.matrix_rank(D, tol=1e-10)


def cauchy_binet_terms(delta, Delta3, dp3_amp, dc3_amp, phi, **Gammas):
    """Explicit 3-term Cauchy-Binet decomposition, for verification."""
    d_p, d_c = dipoles(dp3_amp, dc3_amp, phi)
    D = np.column_stack([d_p, d_c])
    a = excited_resolvent(delta, Delta3, **Gammas)
    terms = {}
    for S in [(0, 1), (0, 2), (1, 2)]:
        D_S = D[list(S), :]
        detD_S = np.linalg.det(D_S)
        detG_S = 1.0 / (a[S[0]] * a[S[1]])
        terms[S] = (np.abs(detD_S) ** 2) * detG_S
    return terms


def residual(x, Delta3, dp3_amp, dc3_amp):
    delta, phi = x
    val, _ = N_full(delta, Delta3, dp3_amp, dc3_amp, phi)
    return [val.real, val.imag]


def search(Delta3_list, dp3_amp=0.5, dc3_amp=0.5, n_delta0=25, n_phi0=25,
           delta_range=(-3.0, 3.0)):
    """Multi-start Newton search for a real (delta, phi) solving N_full=0."""
    solutions = []
    delta0_grid = np.linspace(delta_range[0], delta_range[1], n_delta0)
    phi0_grid = np.linspace(0.0, 2 * np.pi, n_phi0, endpoint=False)
    for Delta3 in Delta3_list:
        found_for_this_Delta3 = []
        for d0 in delta0_grid:
            for p0 in phi0_grid:
                sol, info, ier, msg = fsolve(
                    residual, x0=[d0, p0], args=(Delta3, dp3_amp, dc3_amp),
                    full_output=True, xtol=1e-13,
                )
                if ier != 1:
                    continue
                delta_sol, phi_sol = sol
                res = residual([delta_sol, phi_sol], Delta3, dp3_amp, dc3_amp)
                if max(abs(res[0]), abs(res[1])) > 1e-9:
                    continue
                r = rank_D(dp3_amp, dc3_amp, phi_sol)
                if r < 2:
                    continue  # would be a dark state, not the result we want
                _, S_c = N_full(delta_sol, Delta3, dp3_amp, dc3_amp, phi_sol)
                if abs(S_c) < 1e-8:
                    continue  # pole, not a regular zero
                # dedupe
                dup = False
                for prev in found_for_this_Delta3:
                    if abs(prev[0] - delta_sol) < 1e-6 and abs(
                        ((prev[1] - phi_sol + np.pi) % (2 * np.pi)) - np.pi
                    ) < 1e-6:
                        dup = True
                        break
                if not dup:
                    found_for_this_Delta3.append((delta_sol, phi_sol))
        for delta_sol, phi_sol in found_for_this_Delta3:
            solutions.append(
                {
                    "Delta3": Delta3,
                    "delta": delta_sol,
                    "phi": phi_sol,
                    "rank_D": rank_D(dp3_amp, dc3_amp, phi_sol),
                }
            )
    return solutions


def wide_search():
    """Broader search: also vary dp3_amp, dc3_amp, and restrict delta0 grid
    to bracket each Delta3 so that a3(delta) can pass near resonance
    (needed for the e1-e3/e2-e3 Cauchy-Binet terms to be large enough to
    cancel the e1-e2 term, whose magnitude |det D_01|=sin(theta) is fixed)."""
    all_solutions = []
    amp_grid = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
    Delta3_grid = np.concatenate([
        np.linspace(-4, 4, 17),
    ])
    for dp3_amp in amp_grid:
        for dc3_amp in amp_grid:
            for Delta3 in Delta3_grid:
                sols = search(
                    [Delta3], dp3_amp=dp3_amp, dc3_amp=dc3_amp,
                    n_delta0=15, n_phi0=15,
                    delta_range=(Delta3 - 2.0, Delta3 + 2.0),
                )
                for s in sols:
                    s["dp3_amp"] = dp3_amp
                    s["dc3_amp"] = dc3_amp
                all_solutions.extend(sols)
    return all_solutions


def main():
    dp3_amp, dc3_amp = 0.5, 0.5
    Delta3_list = np.linspace(-6, 6, 25)
    solutions = search(Delta3_list, dp3_amp=dp3_amp, dc3_amp=dc3_amp)

    report = {
        "baseline": {
            "gamma_e": gamma_e, "Delta_e2": Delta_e2, "theta": theta,
            "dp3_amp": dp3_amp, "dc3_amp": dc3_amp,
        },
        "n_Delta3_scanned": len(Delta3_list),
        "n_exact_zeros_found": len(solutions),
        "solutions": solutions[:20],
    }

    if solutions:
        # verify the first solution in full detail, incl. Cauchy-Binet terms
        s0 = solutions[0]
        terms = cauchy_binet_terms(
            s0["delta"], s0["Delta3"], dp3_amp, dc3_amp, s0["phi"]
        )
        total = sum(terms.values())
        report["verification_example"] = {
            "params": s0,
            "cauchy_binet_terms": {str(k): [v.real, v.imag] for k, v in terms.items()},
            "sum_of_terms": [total.real, total.imag],
            "N_full_direct": [
                N_full(s0["delta"], s0["Delta3"], dp3_amp, dc3_amp, s0["phi"])[0].real,
                N_full(s0["delta"], s0["Delta3"], dp3_amp, dc3_amp, s0["phi"])[0].imag,
            ],
        }
        print(f"FOUND {len(solutions)} exact dark-state-free real-axis zeros.")
        print("Example:", s0)
        print("Cauchy-Binet terms (should sum to N_full, all nonzero individually):")
        for k, v in terms.items():
            print(f"  det(G_{{{k}}}) term = {v:.6f}")
        print("Sum of terms:", total, " vs N_full direct:",
              report["verification_example"]["N_full_direct"])
    else:
        print("NO exact real-axis dark-state-free zero found in this scan.")
        print("=> candidate extended no-go result; widen scan or change ansatz before concluding.")

    def _default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        raise TypeError(f"Object of type {o.__class__.__name__} is not JSON serializable")

    out_path = RESULTS_DIR / "gate_E4_2g3e_search.json"
    out_path.write_text(json.dumps(report, indent=2, default=_default))
    print("Report written to", out_path)


if __name__ == "__main__":
    main()
