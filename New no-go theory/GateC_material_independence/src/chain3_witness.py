"""chain3_witness.py -- minimal NON-DIAMOND class-3 (n=3) path-moment witness.

Gate C claims material independence: the integer suppression class is set by
the input-dynamics-readout path structure, not by the host material. Classes
1 and 2 each had a diamond and a non-diamond realization (group-IV / SC
generic, NV 0<->-1 / SC protected), but class 3 had only NV -1<->+1. A class
demonstrated in exactly one material is a property of that material until
shown otherwise, so this module supplies the missing non-diamond n=3 witness.

Model. A three-mode chain in the single-excitation manifold, sharing a common
ground state |g>:

    source -> mode 1 --J12-- mode 2 --J23-- mode 3 -> marked readout

with engineered loss on all three active modes and no direct 1-3 coupling.
With c = |1>, p = |3> the path moments are

    M0 = 0            (p and c are different modes: no direct overlap)
    M1 = 0            (no direct 1-3 link: the graph distance is 2)
    M2 = -J12*J23     (the only path is 1-2-3)

so the amplitude scales as Gamma^-3 and a fixed-readout population as
Gamma^-6. Nothing here is diamond-specific: the class follows from the chain
topology alone, and breaking it (adding a 1-3 bypass) drops the order, which
is the point of the material-independence claim.

This is a coherent-transport witness, not an EIT dark state and not a defect
center, so it also serves as a second non-EIT witness alongside Gate B.

Adapted from the frozen independent-audit bundle
(New no-go theory/PRL_Gates_A_to_D_2026-08-01/src/minimal_nonNV_n3_witness.py),
restructured here to expose the same kernel/full-GKSL interface the other Gate
C systems use.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

# Dimensionless design point: J's set the coherent scale, the detunings are
# deliberately generic (no accidental degeneracy that could create an extra
# cancellation), FLOOR is the Gamma-independent residual loss.
J12, J23 = 1.0, 0.8
DET = (0.22, -0.31, 0.13)
FLOOR = 1e-3
N = 4  # |g> + three single-excitation mode states


def hamiltonian(eps13: float = 0.0) -> np.ndarray:
    """Chain Hamiltonian; eps13 opens the class-breaking direct 1-3 bypass."""
    h = np.zeros((N, N), complex)
    h[1, 1], h[2, 2], h[3, 3] = DET
    h[1, 2] = h[2, 1] = J12
    h[2, 3] = h[3, 2] = J23
    h[1, 3] = h[3, 1] = eps13
    return h


# ----------------------------------------------------------------------
# reduced pencil: A(Gamma) = Gamma*D + A0, D = identity on the active modes
# ----------------------------------------------------------------------
def D_reduced() -> np.ndarray:
    return np.eye(3, dtype=complex)


def A0_reduced(eps13: float = 0.0) -> np.ndarray:
    h = hamiltonian(eps13)[1:, 1:]
    return FLOOR * np.eye(3, dtype=complex) + 1j * h


def c_source() -> np.ndarray:
    return np.array([1.0, 0.0, 0.0], dtype=complex)


def p_readout(leakage: complex = 0.0) -> np.ndarray:
    """Marked readout on mode 3; `leakage` admits a spurious mode-1 component."""
    return np.array([leakage, 0.0, 1.0], dtype=complex)


def kernel(gamma: float, eps13: float = 0.0, leakage: complex = 0.0) -> complex:
    """K(Gamma) = p^dag [Gamma*D + A0]^-1 c -- the reduced transport amplitude."""
    A = gamma * D_reduced() + A0_reduced(eps13)
    return complex(np.conj(p_readout(leakage)) @ np.linalg.solve(A, c_source()))


def moments(nmax: int = 3, eps13: float = 0.0) -> np.ndarray:
    """Path moments M_k = p^dag (-A0)^k c; the first nonzero index + 1 is n."""
    A0 = A0_reduced(eps13)
    p, c = p_readout(), c_source()
    v = c.copy()
    out = []
    for _ in range(nmax + 1):
        out.append(complex(np.conj(p) @ v))
        v = (-A0) @ v
    return np.asarray(out)


def exact_certificate() -> dict:
    """Exact symbolic moments and the resonant kernel, in closed form.

    Kept symbolic so the M2 = -J12*J23 statement is a proof about the chain
    topology rather than a numerical coincidence at this design point.
    """
    Gamma, j1, j2, d1, d2, d3, g0 = sp.symbols(
        "Gamma J12 J23 d1 d2 d3 gamma0", real=True)
    h = sp.Matrix([[d1, j1, 0], [j1, d2, j2], [0, j2, d3]])
    a0 = g0 * sp.eye(3) + sp.I * h
    c = sp.Matrix([1, 0, 0])
    p = sp.Matrix([0, 0, 1])
    mom = [sp.simplify((p.T * (-a0) ** n * c)[0]) for n in range(3)]
    k = sp.simplify((p.T * (Gamma * sp.eye(3) + a0).inv() * c)[0])
    k_res = sp.factor(k.subs({d1: 0, d2: 0, d3: 0, g0: 0}))
    return {"M0": str(mom[0]), "M1": str(mom[1]), "M2": str(mom[2]),
            "K_resonant": str(k_res)}


# ----------------------------------------------------------------------
# full GKSL, for the reduced-equals-full audit the other Gate C systems get
# ----------------------------------------------------------------------
def jump_operators(gamma: float) -> list:
    ops = []
    for j in (1, 2, 3):
        c = np.zeros((N, N), complex)
        c[0, j] = np.sqrt(2.0 * (gamma + FLOOR))
        ops.append(c)
    return ops


def liouvillian(h: np.ndarray, cs: list) -> np.ndarray:
    dim = N * N
    lio = np.zeros((dim, dim), complex)
    for idx in range(dim):
        rho = np.zeros(dim, complex)
        rho[idx] = 1.0
        rho = rho.reshape(N, N)
        drho = -1j * (h @ rho - rho @ h)
        for c in cs:
            cd = c.conj().T
            drho += c @ rho @ cd - 0.5 * (cd @ c @ rho + rho @ cd @ c)
        lio[:, idx] = drho.reshape(-1)
    return lio


def _solve_trace(lio: np.ndarray, rhs: np.ndarray, tr_value: float) -> np.ndarray:
    m = lio.copy()
    b = rhs.reshape(-1).copy()
    tr = np.zeros(N * N, complex)
    for i in range(N):
        tr[i * N + i] = 1.0
    m[0, :] = tr
    b[0] = tr_value
    return np.linalg.solve(m, b).reshape(N, N)


def full_response(gamma: float, eps13: float = 0.0) -> tuple:
    """First- and second-order weak-drive response of the full 16-dim GKSL.

    Returns (amplitude rho1[3,0], mode-3 population rho2[3,3]).
    """
    h = hamiltonian(eps13)
    lio = liouvillian(h, jump_operators(gamma))
    rho0 = np.zeros((N, N), complex)
    rho0[0, 0] = 1.0
    v = np.zeros((N, N), complex)
    v[0, 1] = v[1, 0] = 1.0
    rho1 = _solve_trace(lio, 1j * (v @ rho0 - rho0 @ v), 0.0)
    rho2 = _solve_trace(lio, 1j * (v @ rho1 - rho1 @ v), 0.0)
    return complex(rho1[3, 0]), float(rho2[3, 3].real)


def full_vs_reduced_max_rel_err(gammas) -> float:
    """The full GKSL amplitude and the reduced kernel must be the same object."""
    worst = 0.0
    for g in gammas:
        amp, _ = full_response(g)
        red = kernel(g)
        # the weak-drive amplitude carries a fixed -i from the commutator source
        worst = max(worst, abs(abs(amp) - abs(red)) / max(abs(red), 1e-300))
    return float(worst)
