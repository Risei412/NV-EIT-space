"""Minimal 3-level Lindblad model for Gate U5 (frozen vs self-consistent
cut, Proposition 0D / Lemma 0D'), matching Sec. 11 of the revision note:
sector S = ground coherence span{|g1><g2|, |g2><g1|}, cut generator = pure
dephasing on |g2> (Hermitian jump operator L_S = |g2><g2|, so Lemma 0B'
applies and the Riesz projector is orthogonal).

Basis order: [g1, g2, e] = indices [0, 1, 2]. Unperturbed dynamics: e decays
entirely to g1 (branch13=1), so the unique unperturbed steady state is
rho0 = |g1><g1|, a diagonal state. Since dephasing conserves the diagonal
of any state written in its own eigenbasis, D[L_S](rho0) = 0 automatically
for *any* diagonal rho0 -- this is the concrete instance of condition (C4)
used to verify Lemma 0D' numerically (Gate U5).
"""

import numpy as np

DIM = 3


def vec(rho):
    return rho.reshape(-1)


def unvec(v):
    return v.reshape(DIM, DIM)


def liouvillian(H=None, jump_ops=None):
    dim = DIM * DIM
    L = np.zeros((dim, dim), dtype=complex)
    jump_ops = jump_ops or []
    for idx in range(dim):
        e = np.zeros(dim, dtype=complex)
        e[idx] = 1.0
        rho = unvec(e)
        drho = np.zeros((DIM, DIM), dtype=complex)
        if H is not None:
            drho += -1j * (H @ rho - rho @ H)
        for C in jump_ops:
            Cd = C.conj().T
            drho += C @ rho @ Cd - 0.5 * (Cd @ C @ rho + rho @ Cd @ C)
        L[:, idx] = vec(drho)
    return L


def decay_ops(Gamma, eps=1e-3):
    """e -> g1 only (branch13=1), scaled by native dissipation Gamma, plus
    a weak Gamma-independent incoherent relaxation g2 -> g1 that breaks the
    otherwise degenerate steady-state manifold (any diagonal split between
    g1, g2 with no e population solves L0 v=0 without it) and pins the
    unique steady state to rho0 = |g1><g1| exactly, since this extra jump
    operator L = |g1><g2| also annihilates |g1><g1| under D[L]:
    D[L](rho0) = L rho0 L^dag - (1/2){L^dag L, rho0} = 0 for rho0=|g1><g1|
    (direct computation, since <g2|g1>=0)."""
    C1 = np.sqrt(Gamma) * np.outer([1, 0, 0], [0, 0, 1])  # |g1><e|
    C2 = np.sqrt(eps) * np.outer([1, 0, 0], [0, 1, 0])    # |g1><g2|
    return [C1, C2]


def dephasing_op_g2():
    """Hermitian jump operator L_S = |g2><g2| (pure dephasing of the
    ground coherence). D_S := D[L_S] as a Liouville-space superoperator is
    the admissible cut generator for this sector."""
    return np.outer([0, 1, 0], [0, 1, 0])


def unperturbed_hamiltonian(Delta3=1.3):
    H = np.zeros((DIM, DIM), dtype=complex)
    H[2, 2] = Delta3
    return H


def L0(Gamma, Delta3=1.3):
    """Unperturbed Lindbladian (no probe/control coupling): H diagonal +
    decay e -> g1. Steady state is exactly rho0 = |g1><g1|."""
    return liouvillian(H=unperturbed_hamiltonian(Delta3), jump_ops=decay_ops(Gamma))


def C_S_super():
    """D[L_S] as a 9x9 superoperator (unit rate; multiplied by kappa)."""
    return liouvillian(jump_ops=[dephasing_op_g2()])


def steady_state(L, tol=1e-10):
    """Solve L @ vec(rho) = 0 with trace(rho) = 1 by replacing one row."""
    dim = DIM * DIM
    A = L.copy()
    A[0, :] = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1], dtype=complex)
    b = np.zeros(dim, dtype=complex)
    b[0] = 1.0
    v = np.linalg.solve(A, b)
    return unvec(v)


def rho0(Gamma=5.0, Delta3=1.3):
    return steady_state(L0(Gamma, Delta3))
