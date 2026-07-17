"""Phase D, run D4: preasymptotic masquerade vs the finite polynomial
certificate.

4-site graph 0-1-2-3 (chain v01,v12,v23) plus two competing 2-hop
shortcuts, 0-1-3 (v13) and 0-2-3 (v02). At v02 = v02_star = -0.35 the two
2-hop contributions to the master moment m_2 cancel exactly (m_2 = 0,
verified independently by the exact symbolic certificate below); backing
off by a small detuning delta leaves m_2 anomalously small but nonzero
while the genuine 3-hop chain moment m_3 stays O(1). The result is a
multi-decade false plateau at nu_eff = 4 (dominated by the O(1) m_3 term)
before crossing over, at Gamma ~ m_3/m_2, to the true asymptotic order
nu = 3 (set by the small but nonzero m_2).
"""

import numpy as np
import sympy as sp

V01, V12, V23, V13 = 1.0, 1.0, 1.0, 0.3
D_DIAG = (1.0, 1.2, 1.4, 1.6)
V02_STAR = -0.35  # exact cancellation point (verified symbolically below)


def build(v02, v01=V01, v12=V12, v23=V23, v13=V13, d=D_DIAG):
    D = np.diag(d)
    B = np.zeros((4, 4))
    B[0, 1] = B[1, 0] = v01
    B[1, 2] = B[2, 1] = v12
    B[2, 3] = B[3, 2] = v23
    B[1, 3] = B[3, 1] = v13
    B[0, 2] = B[2, 0] = v02

    def B_of_z(z, _B=B):
        return _B

    c = np.array([1.0, 0.0, 0.0, 0.0])
    p = np.array([0.0, 0.0, 0.0, 1.0])
    return D, B_of_z, c, p


def R_gamma(Gamma, v02):
    D, B_of_z, c, p = build(v02)
    A = Gamma * D + B_of_z(0.0)
    return np.conj(p) @ np.linalg.solve(A, c)


def symbolic_cancellation_check():
    """Exact rational verification that m_2 = 0 identically at v02_star
    (Theorem 8.1 style finite exact test, not a floating-point coincidence)."""
    v02 = sp.Rational(-35, 100)
    D = sp.diag(*[sp.nsimplify(x) for x in D_DIAG])
    B = sp.zeros(4, 4)
    B[0, 1] = B[1, 0] = sp.nsimplify(V01)
    B[1, 2] = B[2, 1] = sp.nsimplify(V12)
    B[2, 3] = B[3, 2] = sp.nsimplify(V23)
    B[1, 3] = B[3, 1] = sp.nsimplify(V13)
    B[0, 2] = B[2, 0] = v02
    Dinv = D.inv()
    X = Dinv * B
    v = Dinv * sp.Matrix([1, 0, 0, 0])
    p = sp.Matrix([0, 0, 0, 1])
    m0 = (p.T * v)[0, 0]
    v = X * v
    m1 = (p.T * v)[0, 0]
    v = X * v
    m2 = (p.T * v)[0, 0]
    v = X * v
    m3 = (p.T * v)[0, 0]
    return [sp.simplify(m) for m in (m0, m1, m2, m3)]


def certificate_true_nu(delta):
    """Finite polynomial certificate (Theorem 8.1 Step 5): exact degree
    difference deg_Gamma Q - deg_Gamma N, computed once, valid for every
    Gamma (no asymptotic fit, no false-plateau ambiguity)."""
    from core import certificate_deg_nu

    Gamma = sp.symbols("Gamma")
    v02 = sp.Rational(-35, 100) + sp.nsimplify(delta)
    D = sp.diag(*[sp.nsimplify(x) for x in D_DIAG])
    B = sp.zeros(4, 4)
    B[0, 1] = B[1, 0] = sp.nsimplify(V01)
    B[1, 2] = B[2, 1] = sp.nsimplify(V12)
    B[2, 3] = B[3, 2] = sp.nsimplify(V23)
    B[1, 3] = B[3, 1] = sp.nsimplify(V13)
    B[0, 2] = B[2, 0] = v02
    c = sp.Matrix([1, 0, 0, 0])
    p = sp.Matrix([0, 0, 0, 1])

    nu, Qp, Np = certificate_deg_nu(D, lambda z: B, c, p, Gamma, z_val=0)
    return nu
