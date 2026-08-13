"""Phase D, run D2: cancellation-promoted suppression and the
symmetry-free exact zero.

Part 1 (order promotion, numeric): a 4-site sector graph with two parallel
2-hop paths (0-1-3 and 0-2-3) whose contributions to the first nonzero
master moment can be tuned to nearly cancel, while a genuine 3-hop path
(0-1-2-3, opened by a small cross-coupling) sets the next moment. Compared
naively (per branch) the graph-distance bound predicts nu = 3 (2 hops + 1);
the actual first nonzero *master* moment is pushed to k=3 (nu = 4) by the
interference the old per-kernel analysis does not see.

Part 2 (exact, symbolic, augmented Krylov certificate, Theorem I / 4.2):
M_full = M_cut + t * u v^T, a rank-1 sector-opening perturbation. Choosing v
orthogonal (over the rationals) to {c, M_cut c} makes the augmented Krylov
moments l^dag M^k r (k=0,1,2, on the doubled space M=diag(M_full,M_cut),
r=(c,c), l=(p,-p)) vanish identically for *every* t -- an exact, all-order
Class I master response -- while the per-branch cross kernel
p_0^dag (zI-M_full)^-1 p_1 (the direct analogue of K12*K21) is nonzero
pointwise and no reducing symmetry projector commutes with M_full.
"""

import numpy as np
import sympy as sp


# ----------------------------------------------------------------------
# Part 1: order promotion via near-cancelling parallel paths
# ----------------------------------------------------------------------
def build_promotion_model(v01=1.0, v02=1.0, v13=0.6, v23=-0.6, v12=0.05,
                           d=(1.0, 1.3, 1.7, 2.1)):
    """4-site graph 0,1,2,3; c on 0, p on 3.
    Edges: 0-1 (v01), 0-2 (v02), 1-3 (v13), 2-3 (v23), 1-2 (v12, small).
    D = diag(d). Returns D, B_of_z (z-independent), c, p."""
    D = np.diag(d)
    B = np.zeros((4, 4))
    B[0, 1] = B[1, 0] = v01
    B[0, 2] = B[2, 0] = v02
    B[1, 3] = B[3, 1] = v13
    B[2, 3] = B[3, 2] = v23
    B[1, 2] = B[2, 1] = v12

    def B_of_z(z, _B=B):
        return _B

    c = np.array([1.0, 0.0, 0.0, 0.0])
    p = np.array([0.0, 0.0, 0.0, 1.0])
    return D, B_of_z, c, p


def moment_ladder(D, B_of_z, c, p, z, kmax=4):
    from core import moments_invertible_D
    return moments_invertible_D(D, B_of_z, c, p, z, kmax)


def tune_v23_for_small_m2(v01=1.0, v02=1.0, v13=0.6, d=(1.0, 1.3, 1.7, 2.1),
                           target=1e-8):
    """Root-find v23 so that the 2-hop master moment m_2 is anomalously
    small (near-cancellation between the two parallel 2-hop paths)."""
    from scipy.optimize import brentq

    def m2_of_v23(v23):
        D, B_of_z, c, p = build_promotion_model(v01, v02, v13, v23, v12=0.0, d=d)
        mu = moment_ladder(D, B_of_z, c, p, z=0.0, kmax=3)
        return mu[2].real

    lo, hi = -0.9, -0.3
    v23_star = brentq(m2_of_v23, lo, hi, xtol=1e-14)
    return v23_star


# ----------------------------------------------------------------------
# Part 2: exact symbolic Class I with no reducing symmetry, K != 0
# ----------------------------------------------------------------------
def build_exact_zero_certificate():
    """Rational M_cut, u; solve for v (rational) orthogonal to {c, M_cut c}.
    Returns (M_cut, u, v, c, p) with all entries sp.Rational, plus a
    verification that the augmented Krylov moments vanish for symbolic t."""
    M_cut = sp.Matrix([
        [sp.Rational(2, 1), sp.Rational(1, 3), sp.Rational(-1, 5)],
        [sp.Rational(1, 4), sp.Rational(3, 2), sp.Rational(1, 7)],
        [sp.Rational(-1, 6), sp.Rational(2, 9), sp.Rational(5, 4)],
    ])
    c = sp.Matrix([1, 2, 1])
    p = sp.Matrix([1, -1, 3])
    u = sp.Matrix([1, 1, 1])

    Mc = M_cut * c
    # v must satisfy v^T c = 0 and v^T (M_cut c) = 0 -> nullspace of the
    # 2x3 matrix whose rows are c^T and (M_cut c)^T
    Rows = sp.Matrix.vstack(c.T, Mc.T)
    ns = Rows.nullspace()
    assert len(ns) == 1, "expected a 1-dim nullspace in this generic 3-dim example"
    v = ns[0]
    v = v / sp.gcd(list(v))  # clear common factor for a clean rational vector

    t = sp.symbols("t")
    M_full = M_cut + t * (u * v.T)

    # augmented Krylov moments mu_k = p^T M_full^k c - p^T M_cut^k c, k=0,1,2
    moments = []
    for k in range(3):
        mu = (p.T * M_full**k * c)[0, 0] - (p.T * M_cut**k * c)[0, 0]
        moments.append(sp.simplify(mu))

    return {
        "M_cut": M_cut, "u": u, "v": v, "c": c, "p": p, "t": t,
        "M_full_symbolic": M_full,
        "moments": moments,
    }


def cross_kernel(M, i, j, z):
    """K_ij(z) = e_i^T (zI - M)^{-1} e_j, the abstract analogue of the
    coherent-Raman kernel K12 (lecture Eq. 11)."""
    n = M.shape[0]
    R = (z * sp.eye(n) - M).inv()
    return R[i, j]
