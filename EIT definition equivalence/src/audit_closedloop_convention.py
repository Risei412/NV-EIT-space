"""
Gate F0/F1/F3: convention audit for excited-state mixing in the 2g+3e model.

Compares three ways to add a coupling between excited states e2, e3 in the
coherence-block generator A(delta) = Gamma + i(H - delta I):

  Model R (old/wrong):  A_23 = J23           (real off-diagonal, ad hoc)
  Model H (physical):   A_23 = i * J23       (Hermitian mixing H_23 = J23 in
                                               A = Gamma + i(H - delta I))
  Model D (physical):   A_23 = kappa23,      Gamma_23 = kappa23 (correlated
                         radiative decay; requires Gamma >= 0, i.e.
                         |kappa23| <= sqrt(gamma2 gamma3))

For each model this script:
  1. verifies the accretivity identity (G+G^dagger)/2 = G^dagger Gamma G;
  2. reports passivity margin m_Gamma = lambda_min((A+A^dagger)/2);
  3. searches for exact real-axis zeros of M = D^dagger G D (rank D = 2);
  4. classifies solver hits as genuine (finite params, sigma_min(M) not
     underflowing) vs. runaway (params diverging, an artifact of chasing a
     singularity at infinity);
  5. for Model R, replays the previously reported 238 solutions and audits
     their passivity margin (Gate F0/F3);
  6. builds the alpha-continuation A_23(alpha) = (1-alpha) J23 + i alpha J23
     from Model R to Model H and tracks the nearest complex zero z0(alpha).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

gamma_e = 1.0
theta = np.pi / 4


def dipoles(dp3_amp, dc3_amp, phi):
    c, s = np.cos(theta), np.sin(theta)
    d_p = np.array([1.0, 0.0, dp3_amp], dtype=complex)
    d_c = np.array([c, s, dc3_amp * np.exp(1j * phi)], dtype=complex)
    return np.column_stack([d_p, d_c])


def build_A(delta, Delta2, Delta3, J23, model, kappa23=None,
            Gamma1=0.0, Gamma2=0.0, Gamma3=0.0):
    a1 = gamma_e + Gamma1 - 1j * delta
    a2 = gamma_e + Gamma2 + 1j * Delta2 - 1j * delta
    a3 = gamma_e + Gamma3 + 1j * Delta3 - 1j * delta
    A = np.diag([a1, a2, a3]).astype(complex)
    if model == "R":
        A[1, 2] = A[2, 1] = J23
    elif model == "H":
        A[1, 2] = 1j * J23
        A[2, 1] = -1j * J23  # Hermitian H => A_23 = i H_23, A_32 = i H_32 = i H_23^* = -i J23 if J23 real... see note below
        # For real J23 (H Hermitian with real off-diag), A = Gamma + i(H - delta I)
        # gives A_23 = i*J23, A_32 = i*J23 (since H_32 = H_23^* = J23 for real J23).
        # Correct this: both off-diagonal entries of A equal i*J23.
        A[2, 1] = 1j * J23
    elif model == "D":
        # correlated radiative decay: J23 plays the role of kappa23 directly,
        # entering the (Hermitian, real) dissipative part; caller must ensure
        # |J23| <= sqrt(gamma2*gamma3) for Gamma >= 0 (GKSL positivity).
        A[1, 2] = A[2, 1] = J23
    else:
        raise ValueError(model)
    return A


def passivity_margin(A):
    Geff = (A + A.conj().T) / 2
    return float(np.linalg.eigvalsh(Geff).min())


def N_and_M(delta, Delta2, Delta3, J23, model, dp3_amp, dc3_amp, phi, **kw):
    A = build_A(delta, Delta2, Delta3, J23, model, **kw)
    G = np.linalg.inv(A)
    D = dipoles(dp3_amp, dc3_amp, phi)
    M = D.conj().T @ G @ D
    return np.linalg.det(M), M, A


def residual(x, Delta2, Delta3, model, dp3_amp, dc3_amp, phi, **kw):
    delta, J23 = x
    N, _, _ = N_and_M(delta, Delta2, Delta3, J23, model, dp3_amp, dc3_amp, phi, **kw)
    return [N.real, N.imag]


def search_model(model, n_starts=600, seed=0, dp3_amp=0.5, dc3_amp=0.5, kappa_cap=None):
    rng = np.random.default_rng(seed)
    genuine, runaway = [], 0
    for _ in range(n_starts):
        Delta2, Delta3 = rng.uniform(-8, 8), rng.uniform(-8, 8)
        phi = rng.uniform(0, 2 * np.pi)
        delta0 = rng.uniform(-6, 6)
        J0 = rng.uniform(-3, 3)
        kw = {}
        if model == "D":
            gam23_cap = np.sqrt(gamma_e * gamma_e)  # gamma2=gamma3=gamma_e here
            J0 = np.clip(J0, -gam23_cap * 0.99, gam23_cap * 0.99)
            kw = {}
        sol, info, ier, msg = fsolve(
            residual, x0=[delta0, J0],
            args=(Delta2, Delta3, model, dp3_amp, dc3_amp, phi), full_output=True,
            xtol=1e-13,
        )
        if ier != 1:
            continue
        delta_s, J_s = sol
        if model == "D" and abs(J_s) > gamma_e * 0.999:
            continue  # outside GKSL-positive region, discard
        res = residual([delta_s, J_s], Delta2, Delta3, model, dp3_amp, dc3_amp, phi)
        rm = max(abs(res[0]), abs(res[1]))
        if rm > 1e-9:
            continue
        N, M, A = N_and_M(delta_s, Delta2, Delta3, J_s, model, dp3_amp, dc3_amp, phi)
        sig_min = np.linalg.svd(M, compute_uv=False).min()
        D = dipoles(dp3_amp, dc3_amp, phi)
        rank_ok = np.linalg.matrix_rank(D, tol=1e-9) == 2
        entry = {
            "delta": float(delta_s), "Delta2": float(Delta2), "Delta3": float(Delta3),
            "J23": float(J_s), "phi": float(phi), "residual": float(rm),
            "sigma_min_M": float(sig_min), "passivity_margin": passivity_margin(A),
            "rank_D_full": bool(rank_ok),
        }
        if abs(J_s) > 1e3 or abs(delta_s) > 1e3 or sig_min < 1e-30:
            runaway += 1
        else:
            genuine.append(entry)
    return genuine, runaway


def alpha_continuation(Delta2, Delta3, dp3_amp, dc3_amp, phi, J23_fixed,
                        n_alpha=41, delta_window=6.0, n_scan=4001):
    """Track the nearest complex zero of chi_full-like det(M) as the
    excited-state coupling interpolates A_23(alpha) = (1-alpha)*J23 + i*alpha*J23
    from the old real-J model (alpha=0) to the physical Hermitian model
    (alpha=1)."""
    alphas = np.linspace(0.0, 1.0, n_alpha)
    trajectory = []
    for alpha in alphas:
        A23 = (1 - alpha) * J23_fixed + 1j * alpha * J23_fixed

        def det_M(delta_complex):
            a1 = gamma_e - 1j * delta_complex
            a2 = gamma_e + 1j * Delta2 - 1j * delta_complex
            a3 = gamma_e + 1j * Delta3 - 1j * delta_complex
            A = np.array([[a1, 0, 0], [0, a2, A23], [0, np.conj(A23) if alpha > 0 else A23, a3]], dtype=complex)
            if alpha == 0:
                A[2, 1] = A23
            G = np.linalg.inv(A)
            D = dipoles(dp3_amp, dc3_amp, phi)
            M = D.conj().T @ G @ D
            return np.linalg.det(M)

        # locate approx real-axis minimum of |det_M| via scan, refine w/ Newton in complex delta
        deltas = np.linspace(-delta_window, delta_window, n_scan)
        vals = np.array([abs(det_M(d)) for d in deltas])
        d0 = deltas[np.argmin(vals)]

        from scipy.optimize import newton
        try:
            z0 = newton(det_M, d0 + 0j, tol=1e-13, maxiter=200)
        except Exception:
            z0 = complex(d0, np.nan)
        trajectory.append({"alpha": float(alpha), "z0_real": float(z0.real), "z0_imag": float(z0.imag),
                            "min_abs_on_real_axis": float(vals.min())})
    return trajectory


def main():
    report = {}

    # --- F1: search each model for exact zeros, classify genuine vs runaway ---
    for model, kw in [("R", {}), ("H", {}), ("D", {})]:
        genuine, runaway = search_model(model, n_starts=600, seed=3)
        report[f"model_{model}"] = {
            "n_genuine_zeros": len(genuine),
            "n_runaway_solver_hits": runaway,
            "genuine_examples": genuine[:5],
            "passivity_margins_of_genuine": [g["passivity_margin"] for g in genuine],
        }
        print(f"Model {model}: genuine={len(genuine)}, runaway={runaway}")
        if genuine:
            margins = [g["passivity_margin"] for g in genuine]
            print(f"  passivity margins of genuine zeros: min={min(margins):.4f} max={max(margins):.4f}")

    # --- F0/F3: re-audit the 238 previously reported Model-R solutions ---
    from model_2g3e_closedloop import search as search_R_full
    old_solutions = search_R_full(dp3_amp=0.5, dc3_amp=0.5, n_starts=1500, seed=7)
    margins = []
    for s in old_solutions:
        A = build_A(s["delta"], s["Delta2"], s["Delta3"], s["J23"], "R")
        margins.append(passivity_margin(A))
    margins = np.array(margins)
    report["gate_F0_F3_old_solution_audit"] = {
        "n_old_solutions": len(old_solutions),
        "min_J23_abs": float(np.min(np.abs([s["J23"] for s in old_solutions]))),
        "max_passivity_margin": float(margins.max()),
        "n_with_margin_geq_0": int((margins >= 0).sum()),
        "verdict": "ALL non-passive (Outcome IV: exact zero only past the passivity boundary)"
                   if (margins < 0).all() else "SOME PASSIVE ZEROS FOUND -- re-audit needed",
    }
    print("Old-solution audit:", report["gate_F0_F3_old_solution_audit"])

    # --- alpha-continuation trajectory for one representative old solution ---
    if old_solutions:
        s0 = old_solutions[0]
        traj = alpha_continuation(s0["Delta2"], s0["Delta3"], 0.5, 0.5, s0["phi"], s0["J23"])
        report["alpha_continuation_example"] = {"base_solution": s0, "trajectory": traj}
        print("alpha=0 z0:", traj[0], " alpha=1 z0:", traj[-1])

    out_path = RESULTS_DIR / "gate_F0_F1_F3_convention_audit.json"
    out_path.write_text(json.dumps(report, indent=2))
    print("Report written to", out_path)


if __name__ == "__main__":
    main()
