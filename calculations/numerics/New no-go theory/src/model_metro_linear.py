"""Phase M1-M2 (abstract arena): vector-valued tangent response and its
EXACT hidden class transition via a rank-1 sector cut.

Motivates the metrological prediction:
  x_{S,Gamma} = d(rho_full)/dtheta - d(rho_cut)/dtheta
in the abstract master-response arena, where theta enters only through the
source term c(theta) (equivalently, via implicit differentiation of
A_Gamma(theta) x = c). Norm scaling ||x_{S,Gamma}|| ~ Gamma^{-nu} is tested
against the same three estimators used in the scalar case (direct fit,
moment method, protected-coefficient method), generalized to vector
readouts.

Rank-1 cut and exact vector transition
---------------------------------------
For a Hermitian (physically realizable) coupling being severed, the
protected-block difference B_P,full^-1 - B_P,cut^-1 is generically rank >=2
in the 2-dim protected block, so a single real control parameter lambda
cannot zero the full *vector* x_{S,0}(lambda) except at the trivial point
(coupling = 0) -- confirmed analytically before coding this module. A
*rank-1* cut, B_full(lambda) = B_base + u w(lambda)^dagger with w(lambda) =
w0 + lambda*w1 linear in lambda, instead gives (Sherman-Morrison)

  x_{S,Gamma}(lambda) = s(Gamma,lambda) * [A_cut,Gamma^{-1} u]

i.e. a FIXED direction times a scalar. The whole vector vanishes together
whenever the scalar vanishes. Taking the Gamma -> infinity (Riesz
projection) limit of the scalar's numerator gives a closed-form,
Gamma-independent root lambda_c^0 -- no root search needed.

This rank-1 structure is not meant to model a literal Hamiltonian coupling
(which is Hermitian, i.e. rank-2); it is a legitimate example within the
theorem's abstract scope (Theorem III places no Hermiticity assumption on
B), used here to establish the *exact* existence of a vector hidden class
transition. The physical (Hermitian) Lindblad model in
model_metro_lindblad.py is tested separately for an approximate analogue.
"""

import numpy as np

from model_protected import build_class1_example, build_class2_example, build_class3_example
from core import transfer, fit_nu_loglog


# ----------------------------------------------------------------------
# Gate M1: vector-valued readout preserves nu (Class I/II/III unit tests)
# ----------------------------------------------------------------------
def vector_norm_class1(gammas, seed=7):
    M, r, l = build_class1_example(seed)
    n = M.shape[0]
    rng = np.random.default_rng(seed + 100)
    # l2 must also respect the unobservable subspace (sector 2 only, indices
    # 3:6) for the Class-I exact-zero certificate to survive a second,
    # independent readout direction.
    l2 = np.concatenate([np.zeros(n // 2), rng.standard_normal(n // 2)])
    p_list = [l, l2]
    norms = []
    for g in gammas:
        A = g * np.eye(n) + M
        x = np.linalg.solve(A, r)
        vals = [np.conj(p) @ x for p in p_list]
        norms.append(np.linalg.norm(vals))
    return np.array(norms)


def vector_norm_class2(gammas, d, seed=11):
    """Two independent length-(d+1) chains (different random couplings,
    seed and seed+1), block-combined so a single source feeds both and
    two independent readouts (one per chain, both at graph distance d)
    give a genuinely 2-component vector -- avoiding the pitfall of a
    second readout accidentally sitting closer to the source (which
    would carry its own smaller effective nu and swamp the signal)."""
    D1, B1, c1, p1, n1 = build_class2_example(d, seed)
    D2, B2, c2, p2, n2 = build_class2_example(d, seed + 1)
    D = np.block([[D1, np.zeros((n1, n2))], [np.zeros((n2, n1)), D2]])
    c = np.concatenate([c1, c2])
    p_full1 = np.concatenate([p1, np.zeros(n2)])
    p_full2 = np.concatenate([np.zeros(n1), p2])
    z0 = 0.3 + 0.1j

    def B_of_z(z):
        return np.block([[B1(z), np.zeros((n1, n2))], [np.zeros((n2, n1)), B2(z)]])

    norms = []
    for g in gammas:
        A = g * D + B_of_z(z0)
        x = np.linalg.solve(A, c)
        vals = [np.conj(p_full1) @ x, np.conj(p_full2) @ x]
        norms.append(np.linalg.norm(vals))
    return np.array(norms)


def vector_norm_class3(gammas, seed=3):
    D, B_of_z, c, p, S, k = build_class3_example(seed)
    n = D.shape[0]
    rng = np.random.default_rng(seed + 100)
    p2 = rng.standard_normal(n)
    p_list = [p, p2]
    norms = []
    for g in gammas:
        A = g * D + B_of_z(0.0)
        x = np.linalg.solve(A, c)
        vals = [np.conj(pp) @ x for pp in p_list]
        norms.append(np.linalg.norm(vals))
    return np.array(norms)


# ----------------------------------------------------------------------
# Gate M2: rank-1 cut, exact vector hidden class transition + collapse
# ----------------------------------------------------------------------
def rank1_transition_setup():
    """Fixed, deterministic model data (n=4, P={0,1}, Q={2,3})."""
    D = np.diag([0.0, 0.0, 1.0, 1.7])
    B_base = np.array([
        [2.1, 0.4, 0.15, 0.05],
        [0.4, 1.6, 0.05, 0.10],
        [0.15, 0.05, 1.9, 0.20],
        [0.05, 0.10, 0.20, 2.3],
    ])
    u = np.array([0.8, -0.5, 0.3, 0.1])
    w0 = np.array([0.6, 0.9, 0.0, 0.0])
    w1 = np.array([1.0, -0.4, 0.0, 0.0])
    c1 = np.array([1.0, 0.6, -0.3, 0.2])
    mask = np.array([True, True, False, False])
    return D, B_base, u, w0, w1, c1, mask


def leading_order_lambda_c():
    """Closed-form root of the leading-order (Gamma->infinity) numerator
    w(lambda)^dagger [P B_P,cut^-1 P] c1 = 0, using only the protected
    (P-P) block of B_base -- no numerical root search required."""
    D, B_base, u, w0, w1, c1, mask = rank1_transition_setup()
    B_P = B_base[np.ix_(mask, mask)]
    y0 = np.zeros_like(c1)
    y0[mask] = np.linalg.solve(B_P, c1[mask])  # y0 = P B_P,cut^{-1} P c1
    num0 = w0 @ y0
    num1 = w1 @ y0
    lam_c0 = -num0 / num1
    return lam_c0, y0


def x_S_gamma(Gamma, lam):
    """Exact x_{S,Gamma}(lambda) = A_full,Gamma(lambda)^{-1} c1 - A_cut,Gamma^{-1} c1
    computed by direct linear solves (no asymptotic approximation)."""
    D, B_base, u, w0, w1, c1, mask = rank1_transition_setup()
    w = w0 + lam * w1
    B_full = B_base + np.outer(u, w)
    A_full = Gamma * D + B_full
    A_cut = Gamma * D + B_base
    x_full = np.linalg.solve(A_full, c1)
    x_cut = np.linalg.solve(A_cut, c1)
    return x_full - x_cut


def scan_lambda(gammas_probe, lam_grid):
    """||x_S,Gamma(lambda)|| on a grid, for local nu_eff(lambda) estimation."""
    out = {}
    for g in gammas_probe:
        out[g] = np.array([np.linalg.norm(x_S_gamma(g, lam)) for lam in lam_grid])
    return out
