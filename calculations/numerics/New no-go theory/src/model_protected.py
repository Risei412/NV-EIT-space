"""Phase B: singular-dissipation models for Theorems I-III (priorities
document Sec. 2.2 and Priorities 2-4).

Four constructions:

  build_class1_example  -- Krylov-orthogonal block model (Theorem I, nu=inf)
  build_class2_example   -- full-rank D, sector chain of length d
                             (Theorem II / Corollary graph selection rule,
                             nu_diss = d+1)
  build_class3_example   -- singular D, semisimple kernel, non-Hermitian,
                             O(1) protected channel (Theorem III, nu=0)
  build_transition_model -- singular D, 2-dim protected block P (indices
                             0,1) + 2-dim damped block Q (indices 2,3).
                             The sector cut S severs the *off-diagonal
                             coupling inside the protected block itself*
                             (B_P[0,1], B_P[1,0]), i.e. the coupling
                             between the "optical" (index 0) and "ground
                             coherence" (index 1) components of the
                             protected subspace. D, c, p are identical
                             between full and cut (frozen-source cut).
                             This is the generic case of the EIT
                             sector-cut corollary in which B_P,full !=
                             B_P,cut, so delta_chi_{S,0} is generically
                             nonzero and can be tuned through zero by the
                             control parameter lambda -- the hidden class
                             transition of Priority 3.
"""

import numpy as np


# ----------------------------------------------------------------------
# Class I: engineered Krylov-orthogonal model
# ----------------------------------------------------------------------
def build_class1_example(seed=7):
    rng = np.random.default_rng(seed)
    M1 = rng.standard_normal((3, 3))
    M2 = rng.standard_normal((3, 3))
    M = np.block([[M1, np.zeros((3, 3))], [np.zeros((3, 3)), M2]])
    r = np.concatenate([rng.standard_normal(3), np.zeros(3)])
    l = np.concatenate([np.zeros(3), rng.standard_normal(3)])
    return M, r, l


# ----------------------------------------------------------------------
# Class II: sector chain of length d (graph selection rule, Corollary)
# ----------------------------------------------------------------------
def build_class2_example(d, seed=11):
    """n = d+1 site chain, c on site 0, p on site d -> nu_diss = d+1."""
    n = d + 1
    rng = np.random.default_rng(seed)
    D = np.diag(rng.uniform(1.0, 2.0, n))
    V = np.zeros((n, n))
    for i in range(n - 1):
        V[i, i + 1] = V[i + 1, i] = rng.uniform(0.5, 1.5)

    diag_coeff = rng.uniform(0.5, 1.0, n)

    def B_of_z(z, _diag=diag_coeff, _V=V):
        return np.diag(_diag * z) + _V

    c = np.eye(n)[0]
    p = np.eye(n)[d]
    return D, B_of_z, c, p, n


# ----------------------------------------------------------------------
# Class III: non-Hermitian D, semisimple kernel, protected O(1) channel
# ----------------------------------------------------------------------
def build_class3_example(seed=3, n=5, k=2):
    rng = np.random.default_rng(seed)
    S = rng.standard_normal((n, n)) + 0.1 * np.eye(n)
    dvals = np.concatenate([np.zeros(k), rng.uniform(1.0, 3.0, n - k)])
    D = S @ np.diag(dvals) @ np.linalg.inv(S)
    B0 = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    c = rng.standard_normal(n)
    p = rng.standard_normal(n)

    def B_of_z(z):
        return B0  # z-independent for this unit-test model

    return D, B_of_z, c, p, S, k


# ----------------------------------------------------------------------
# Hidden-class-transition model (Priority 3-4)
# ----------------------------------------------------------------------
def transition_matrices(lam, e1=3.0, kappa0=0.3, e2=2.0, d3=1.0, d4=1.7,
                         c_vec=None, p_vec=None):
    """Returns D (4x4, singular, diagonal, ker = {0,1}), B_full(lam),
    B_cut(lam), c, p.

    Protected block (indices 0,1):
        B_P,full(lam) = [[e1, kappa0 + lam], [kappa0 + lam, e2]]
        B_P,cut(lam)  = diag(e1, e2)            (coupling severed)

    Damped block (indices 2,3): D_Q = diag(d3, d4), fixed detuning-like
    diagonal B_QQ, plus a small fixed leakage between P and Q that is
    identical in full and cut (only the P-internal coupling is cut).

    delta_chi_{S,0}(kappa) = F0_full(kappa) - F0_cut always has the
    trivial root kappa=0 (cutting a zero coupling changes nothing); with
    the default (e1,e2,c_vec,p_vec) below it also has a *nontrivial* root
    near kappa* ~= 1.0, well separated from both kappa=0 and the pole of
    B_P,full at kappa = sqrt(e1*e2) ~= 2.449 -- this is the hidden class
    transition point used in Priority 3/4.
    """
    if c_vec is None:
        c_vec = np.array([1.0, 0.1743, 0.3, -0.2])
    if p_vec is None:
        p_vec = np.array([1.0, 0.1743, -0.1, 0.4])

    D = np.diag([0.0, 0.0, d3, d4])

    leak = 0.15  # fixed P<->Q leakage, identical in full and cut

    def make_B(kappa):
        B = np.zeros((4, 4))
        B[0, 0], B[1, 1] = e1, e2
        B[0, 1] = B[1, 0] = kappa
        B[2, 2], B[3, 3] = 1.1, 0.9  # fixed damped-sector detuning offsets
        B[0, 2] = B[2, 0] = leak
        B[1, 3] = B[3, 1] = leak
        return B

    B_full = make_B(kappa0 + lam)
    B_cut = make_B(0.0)
    B_cut[0, 0], B_cut[1, 1] = e1, e2  # keep diagonal, coupling severed
    return D, B_full, B_cut, c_vec, p_vec


def r0_of_lambda(lam, **kwargs):
    """delta_chi_{S,0}(lambda) = p^dag P [B_P,full^-1 - B_P,cut^-1] P c,
    evaluated directly from the 2x2 protected blocks (EIT sector-cut
    corollary, leading order)."""
    D, B_full, B_cut, c, p = transition_matrices(lam, **kwargs)
    mask = np.array([True, True, False, False])
    B_P_full = B_full[np.ix_(mask, mask)]
    B_P_cut = B_cut[np.ix_(mask, mask)]
    Pc = c[mask]
    pP = p[mask]
    F0_full = pP @ np.linalg.solve(B_P_full, Pc)
    F0_cut = pP @ np.linalg.solve(B_P_cut, Pc)
    return F0_full - F0_cut


def R_S_gamma(Gamma, lam, **kwargs):
    """R_S(Gamma, lambda) = chi_full - chi_cut computed by direct linear
    solve of (Gamma D + B)x = c (exact, no asymptotic expansion)."""
    D, B_full, B_cut, c, p = transition_matrices(lam, **kwargs)
    A_full = Gamma * D + B_full
    A_cut = Gamma * D + B_cut
    chi_full = np.conj(p) @ np.linalg.solve(A_full, c)
    chi_cut = np.conj(p) @ np.linalg.solve(A_cut, c)
    return chi_full, chi_cut, chi_full - chi_cut
