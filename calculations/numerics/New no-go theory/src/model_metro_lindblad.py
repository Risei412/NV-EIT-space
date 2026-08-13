"""Phase M3-M5: physical 3-level Lambda Lindblad model, steady-state QFI,
and the numerical test of the central metrological prediction

    F_{Q,S} = x_S^dagger G_rho x_S ~ ||x_{S,Gamma}||^2 ~ Gamma^{-2 nu} .

Level scheme: |1>, |2> ground manifold, |3> excited. H couples |2>-|3>
with strength kappa (the "sector S" control field); theta = two-photon
detuning (coefficient of |2><2|) is the estimated parameter. Fast
radiative decay |3> -> |1>, |3> -> |2> is scaled by Gamma (this is the
singular D of the abstract theorems: its kernel is the "protected"
4-dimensional block of density-matrix elements not involving level 3).
Fixed-rate (Gamma-independent) ground-state dephasing and a weak
all-to-all incoherent "reset" channel (rate epsilon) keep the steady
state unique and full rank in both the full and sector-cut generators;
epsilon is the same in both, consistent with the frozen-source cut
convention (D, and all fixed-rate dissipators, are shared).

The sector cut severs only the kappa coupling (H_cut has kappa=0),
leaving D, the fixed dissipators, and the readout basis untouched.
"""

import numpy as np

DIM = 3  # Hilbert space dimension


def vec(rho):
    return rho.reshape(-1)


def unvec(v):
    return v.reshape(DIM, DIM)


def liouvillian_block(H=None, jump_ops=None):
    """Build the 9x9 superoperator for -i[H,.] plus sum_C D[C], using an
    explicit basis loop (robust, convention-free: vec = row-major flatten,
    used consistently everywhere in this module)."""
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


def fast_decay_ops(branch13=0.5, branch23=0.5):
    """Unit-rate (Gamma will multiply this block) radiative decay from |3>
    to |1> and |2>. D = sum of these dissipators, singular (kernel =
    Liouville elements with no level-3 population/coherence)."""
    C1 = np.sqrt(branch13) * np.outer([1, 0, 0], [0, 0, 1])  # |1><3|
    C2 = np.sqrt(branch23) * np.outer([0, 1, 0], [0, 0, 1])  # |2><3|
    return [C1, C2]


def fixed_dissipators(gamma21=0.02, epsilon=0.01):
    """Gamma-independent, fixed-rate dissipators shared by full and cut:
    ground-coherence dephasing (rate gamma21) plus a weak all-to-all
    incoherent-hopping reset channel (rate epsilon) that guarantees a
    unique, full-rank steady state in both models."""
    dephasing = np.sqrt(gamma21) * (np.diag([1.0, -1.0, 0.0]))
    ops = [dephasing]
    basis = np.eye(DIM)
    for i in range(DIM):
        for j in range(DIM):
            if i != j:
                ops.append(np.sqrt(epsilon / (2 * (DIM - 1))) * np.outer(basis[i], basis[j]))
    return ops


def hamiltonian(theta, kappa, Delta3=1.3):
    H = np.zeros((DIM, DIM), dtype=complex)
    H[2, 2] = Delta3
    H[1, 1] = theta
    H[1, 2] = kappa / 2
    H[2, 1] = np.conj(kappa) / 2
    return H


def dH_dtheta():
    dH = np.zeros((DIM, DIM), dtype=complex)
    dH[1, 1] = 1.0
    return dH


TRACE_ROW = 0  # replace this row of the Liouvillian with the trace constraint


def build_A(Gamma, theta, kappa, gamma21=0.02, epsilon=0.01, Delta3=1.3,
            branch13=0.5, branch23=0.5):
    H = hamiltonian(theta, kappa, Delta3)
    L_H = liouvillian_block(H=H)
    L_fast = liouvillian_block(jump_ops=fast_decay_ops(branch13, branch23))
    L_fixed = liouvillian_block(jump_ops=fixed_dissipators(gamma21, epsilon))
    A = L_H + Gamma * L_fast + L_fixed
    A[TRACE_ROW, :] = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1], dtype=complex)
    return A


def build_dA_dtheta(Gamma, theta, kappa):
    dH = dH_dtheta()
    dL = liouvillian_block(H=dH)
    dL[TRACE_ROW, :] = 0.0
    return dL


def steady_state(Gamma, theta, kappa, **kwargs):
    A = build_A(Gamma, theta, kappa, **kwargs)
    c = np.zeros(DIM * DIM, dtype=complex)
    c[TRACE_ROW] = 1.0
    v = np.linalg.solve(A, c)
    return unvec(v), A, c


def steady_state_and_dtheta(Gamma, theta, kappa, **kwargs):
    """rho(theta) and d(rho)/dtheta by implicit differentiation of
    A_Gamma(theta) vec(rho) = c (c is theta-independent: the trace target)."""
    rho, A, c = steady_state(Gamma, theta, kappa, **kwargs)
    dA = build_dA_dtheta(Gamma, theta, kappa)
    rhs = -dA @ vec(rho)
    drho_vec = np.linalg.solve(A, rhs)
    return rho, unvec(drho_vec)


def x_S_lindblad(Gamma, theta, lam, kappa0=1.5, phi=0.0, **kwargs):
    """x_{S,Gamma} = d(rho_full)/dtheta - d(rho_cut)/dtheta, vec'd (9,).

    kappa_full = (kappa0 + lam) * exp(i*phi): a physically natural
    2-real-parameter control (amplitude via lam, phase via phi) on the
    single Hermitian coupling being cut. A single real parameter (phi=0
    fixed) cannot generically zero the full vector x_{S,0} (its leading
    order lives in a >1-dimensional slice of the protected block, unlike
    the rank-1 abstract construction in model_metro_linear.py); both
    parameters are swept in Gate M4 to search for a deep interior
    suppression point.
    """
    kappa_full = (kappa0 + lam) * np.exp(1j * phi)
    rho_full, drho_full = steady_state_and_dtheta(Gamma, theta, kappa_full, **kwargs)
    rho_cut, drho_cut = steady_state_and_dtheta(Gamma, theta, 0.0, **kwargs)
    return vec(drho_full) - vec(drho_cut), vec(drho_full), vec(drho_cut), rho_full, rho_cut


def qfi(rho, drho, eig_tol=1e-12):
    """Quantum Fisher information (SLD) F = 2 sum_ij |<i|drho|j>|^2/(p_i+p_j),
    dropping near-zero denominator terms (near rank-deficient rho)."""
    p, V = np.linalg.eigh(rho)
    drho_eig = V.conj().T @ drho @ V
    F = 0.0
    n = len(p)
    dropped_weight = 0.0
    for i in range(n):
        for j in range(n):
            denom = p[i] + p[j]
            if denom > eig_tol:
                F += 2 * abs(drho_eig[i, j]) ** 2 / denom
            else:
                dropped_weight += abs(drho_eig[i, j]) ** 2
    return float(F.real), float(dropped_weight), float(np.min(p))


def qfi_cross_term(rho, x_S_mat, drho_cut_mat, eig_tol=1e-12):
    """Cross term 2 Re <x_S, G_rho drho_cut> so that
    F_Q,full - F_Q,cut = F_Q,S + cross_term is available as a diagnostic
    (F_Q,S alone is the quantity the prediction concerns)."""
    p, V = np.linalg.eigh(rho)
    xS_eig = V.conj().T @ x_S_mat @ V
    dc_eig = V.conj().T @ drho_cut_mat @ V
    n = len(p)
    C = 0.0
    for i in range(n):
        for j in range(n):
            denom = p[i] + p[j]
            if denom > eig_tol:
                C += 2 * (xS_eig[i, j] * np.conj(dc_eig[i, j])).real / denom
    return float(C)
