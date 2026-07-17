"""Fast vectorized scan for E4, companion to model_2g3e_zero.py.

Instead of brute-force multistart Newton over (delta, phi) for every
(Delta3, dp3_amp, dc3_amp) combination (too slow), evaluate N_full on a
vectorized (delta, phi) grid for each outer parameter combo, locate the
grid cell with the smallest |N_full|, and refine only that candidate with
fsolve. Much cheaper, still exhaustive over the outer parameters.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve

gamma_e = 1.0
Delta_e2 = 8.0
theta = np.pi / 4

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def N_full_grid(delta, phi, Delta3, dp3_amp, dc3_amp):
    """Vectorized over delta (1D array) and phi (1D array) via broadcasting;
    returns N_full on the outer-product grid, shape (len(delta), len(phi))."""
    d = delta[:, None]
    p = phi[None, :]
    c, s = np.cos(theta), np.sin(theta)

    a1 = gamma_e - 1j * d
    a2 = gamma_e + 1j * Delta_e2 - 1j * d
    a3 = gamma_e + 1j * Delta3 - 1j * d
    g1, g2, g3 = 1.0 / a1, 1.0 / a2, 1.0 / a3

    dp = np.array([1.0, 0.0, dp3_amp], dtype=complex)
    # d_c depends on phi: dc3_amp * exp(i phi), broadcast over phi axis
    dc0 = c * np.ones_like(p)
    dc1 = s * np.ones_like(p)
    dc2 = dc3_amp * np.exp(1j * p) * np.ones_like(d)

    S_p = np.abs(dp[0]) ** 2 * g1 + np.abs(dp[1]) ** 2 * g2 + np.abs(dp[2]) ** 2 * g3
    S_c = np.abs(dc0) ** 2 * g1 + np.abs(dc1) ** 2 * g2 + np.abs(dc2) ** 2 * g3
    K_pc = np.conj(dp[0]) * dc0 * g1 + np.conj(dp[1]) * dc1 * g2 + np.conj(dp[2]) * dc2 * g3
    K_cp = np.conj(dc0) * dp[0] * g1 + np.conj(dc1) * dp[1] * g2 + np.conj(dc2) * dp[2] * g3
    return S_p * S_c - K_pc * K_cp, S_c


def N_full_scalar(x, Delta3, dp3_amp, dc3_amp):
    delta, phi = x
    val, _ = N_full_grid(np.array([delta]), np.array([phi]), Delta3, dp3_amp, dc3_amp)
    v = val[0, 0]
    return [v.real, v.imag]


def rank_D(dp3_amp, dc3_amp, phi):
    c, s = np.cos(theta), np.sin(theta)
    d_p = np.array([1.0, 0.0, dp3_amp])
    d_c = np.array([c, s, dc3_amp * np.exp(1j * phi)])
    D = np.column_stack([d_p, d_c])
    return np.linalg.matrix_rank(D, tol=1e-9)


def scan_one(Delta3, dp3_amp, dc3_amp, n_delta=200, n_phi=200,
             delta_span=4.0):
    delta = np.linspace(Delta3 - delta_span, Delta3 + delta_span, n_delta)
    # also always include a window around delta=0 (probe resonance)
    delta = np.union1d(delta, np.linspace(-delta_span, delta_span, n_delta))
    phi = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)

    val, S_c = N_full_grid(delta, phi, Delta3, dp3_amp, dc3_amp)
    mag = np.abs(val)
    # avoid poles (S_c ~ 0)
    mag = np.where(np.abs(S_c) < 1e-6, np.inf, mag)
    idx = np.unravel_index(np.argmin(mag), mag.shape)
    d0, p0 = delta[idx[0]], phi[idx[1]]
    best_mag = mag[idx]

    sol, info, ier, msg = fsolve(
        N_full_scalar, x0=[d0, p0], args=(Delta3, dp3_amp, dc3_amp),
        full_output=True, xtol=1e-13,
    )
    if ier != 1:
        return None
    delta_sol, phi_sol = sol
    res = N_full_scalar([delta_sol, phi_sol], Delta3, dp3_amp, dc3_amp)
    resmag = max(abs(res[0]), abs(res[1]))
    if resmag > 1e-9:
        return None
    r = rank_D(dp3_amp, dc3_amp, phi_sol)
    if r < 2:
        return None
    _, S_c_sol = N_full_grid(np.array([delta_sol]), np.array([phi_sol]), Delta3, dp3_amp, dc3_amp)
    if abs(S_c_sol[0, 0]) < 1e-8:
        return None
    return {
        "Delta3": float(Delta3), "dp3_amp": float(dp3_amp), "dc3_amp": float(dc3_amp),
        "delta": float(delta_sol), "phi": float(phi_sol % (2 * np.pi)),
        "rank_D": int(r), "grid_best_mag": float(best_mag), "residual": float(resmag),
    }


def main():
    amp_grid = [0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0]
    Delta3_grid = np.linspace(-8, 8, 33)

    solutions = []
    tested = 0
    for dp3_amp in amp_grid:
        for dc3_amp in amp_grid:
            for Delta3 in Delta3_grid:
                tested += 1
                s = scan_one(Delta3, dp3_amp, dc3_amp)
                if s is not None:
                    solutions.append(s)

    report = {
        "baseline": {"gamma_e": gamma_e, "Delta_e2": Delta_e2, "theta": theta},
        "n_configs_tested": tested,
        "n_exact_zeros_found": len(solutions),
        "solutions": solutions[:30],
    }
    print(f"tested {tested} (Delta3, dp3_amp, dc3_amp) configs")
    print(f"found {len(solutions)} exact dark-state-free real-axis zeros")
    for s in solutions[:10]:
        print(s)

    out_path = RESULTS_DIR / "gate_E4_2g3e_fast_scan.json"
    out_path.write_text(json.dumps(report, indent=2))
    print("Report written to", out_path)


if __name__ == "__main__":
    main()
