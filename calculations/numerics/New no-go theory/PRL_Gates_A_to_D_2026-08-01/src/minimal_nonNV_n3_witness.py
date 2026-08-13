"""Reproduce the minimal non-NV n=3 path-moment witness.

Model: ground plus three lossy single-excitation states in a chain 1-2-3.
The source drives state 1 and the marked readout is state 3.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parent
J12, J23 = 1.0, 0.8
DET = (0.22, -0.31, 0.13)
FLOOR = 1e-3
N = 4


def hamiltonian(eps13: float = 0.0) -> np.ndarray:
    h = np.zeros((N, N), complex)
    h[1, 1], h[2, 2], h[3, 3] = DET
    h[1, 2] = h[2, 1] = J12
    h[2, 3] = h[3, 2] = J23
    h[1, 3] = h[3, 1] = eps13
    return h


def jumps(gamma: float) -> list[np.ndarray]:
    ans = []
    for j in (1, 2, 3):
        c = np.zeros((N, N), complex)
        c[0, j] = np.sqrt(2.0 * (gamma + FLOOR))
        ans.append(c)
    return ans


def liouvillian(h: np.ndarray, cs: list[np.ndarray]) -> np.ndarray:
    dim = N * N
    l = np.zeros((dim, dim), complex)
    for idx in range(dim):
        rho = np.zeros(dim, complex)
        rho[idx] = 1.0
        rho = rho.reshape(N, N)
        drho = -1j * (h @ rho - rho @ h)
        for c in cs:
            cd = c.conj().T
            drho += c @ rho @ cd - 0.5 * (cd @ c @ rho + rho @ cd @ c)
        l[:, idx] = drho.reshape(-1)
    return l


def solve_trace(l: np.ndarray, rhs: np.ndarray, tr_value: float) -> np.ndarray:
    m = l.copy()
    b = rhs.reshape(-1).copy()
    tr = np.zeros(N * N, complex)
    for i in range(N):
        tr[i * N + i] = 1.0
    m[0, :] = tr
    b[0] = tr_value
    return np.linalg.solve(m, b).reshape(N, N)


def response(gamma: float, eps13: float = 0.0, eta: float = 0.0) -> tuple[complex, float]:
    h = hamiltonian(eps13)
    l0 = liouvillian(h, jumps(gamma))
    rho0 = np.zeros((N, N), complex)
    rho0[0, 0] = 1.0
    v = np.zeros((N, N), complex)
    v[0, 1] = v[1, 0] = 1.0
    rho1 = solve_trace(l0, 1j * (v @ rho0 - rho0 @ v), 0.0)
    rho2 = solve_trace(l0, 1j * (v @ rho1 - rho1 @ v), 0.0)
    amp = rho1[3, 0] + eta * rho1[1, 0]
    return complex(amp), float(rho2[3, 3].real)


def fit_order(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> float:
    return float(-np.polyfit(np.log(x[mask]), np.log(np.abs(y[mask])), 1)[0])


def symbolic_certificate() -> dict[str, str]:
    gamma, j1, j2, d1, d2, d3, g0 = sp.symbols(
        "Gamma J12 J23 d1 d2 d3 gamma0", real=True
    )
    h = sp.Matrix([[d1, j1, 0], [j1, d2, j2], [0, j2, d3]])
    a0 = g0 * sp.eye(3) + sp.I * h
    c = sp.Matrix([1, 0, 0])
    p = sp.Matrix([0, 0, 1])
    x = -a0
    moments = [sp.simplify((p.T * x**n * c)[0]) for n in range(3)]
    k = sp.simplify((p.T * (gamma * sp.eye(3) + a0).inv() * c)[0])
    kres = sp.factor(k.subs({d1: 0, d2: 0, d3: 0, g0: 0}))
    return {"M0": str(moments[0]), "M1": str(moments[1]), "M2": str(moments[2]), "K_resonant": str(kres)}


def main() -> None:
    gammas = np.logspace(0.8, 4.2, 96)
    tail = gammas > 80
    clean = np.array([response(g)[0] for g in gammas])
    pop = np.array([response(g)[1] for g in gammas])
    bypass = np.array([response(g, eps13=1e-3)[0] for g in gammas])
    leak = np.array([response(g, eta=1e-4)[0] for g in gammas])
    result = {
        "certificate": symbolic_certificate(),
        "clean_amplitude_order": fit_order(gammas, clean, tail),
        "clean_population_order": fit_order(gammas, pop, tail),
        "bypass_order": fit_order(gammas, bypass, gammas > 2e3),
        "readout_leakage_order": fit_order(gammas, leak, gammas > 2e3),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
