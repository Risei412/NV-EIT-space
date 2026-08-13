"""nv_reduced_kernel.py -- shared, importable copy of the validated NV excited-
manifold reduced-kernel machinery (A_Gamma = Gamma*D + A0, D=I, X=-A0,
M_n = p^dag X^n c, K(Gamma) = p^dag (Gamma*D+A0)^-1 c), extracted verbatim
from verify_nv_3E_graph_distance_PRL.py (v6.2 conventions, already validated
there: T1-T8 all PASS, slopes -2.000/-3.000). Kept as a separate module (not
imported from the verify_* script, which runs its full test battery as a
side effect on import) so moment_order_common_pipeline.py can reuse the
identical NV kernel alongside the new group-IV kernel (group_iv_model.py).
"""
import numpy as np

sq2 = 1.0/np.sqrt(2.0)
Sz = np.diag([-1.0, 0.0, 1.0]).astype(complex)
Sx = sq2*np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], complex)
Sy = sq2*np.array([[0, 1j, 0], [-1j, 0, 1j], [0, -1j, 0]], complex)
I3 = np.eye(3, dtype=complex)
sz_o = np.array([[1, 0], [0, -1]], complex)
sx_o = np.array([[0, 1], [1, 0]], complex)
Lz_o = np.array([[0, -1j], [1j, 0]], complex)
I2 = np.eye(2, dtype=complex)

LAM_Z, D_PAR, DELTA1, DELTA2, LAM_PERP = 5.33, 1.42/3, 1.55/2, 0.20, 0.20

def H_3E(lam_perp=LAM_PERP, delta2=DELTA2, delta1=DELTA1,
         xi_x=0.0, xi_y=0.0, Bx=0.0, lam_z=LAM_Z, d_par=D_PAR):
    H = lam_z*np.kron(Lz_o, Sz)
    H += d_par*np.kron(I2, Sz@Sz - (2.0/3.0)*I3)
    H += delta1*np.kron(sz_o, (Sy@Sy - Sx@Sx))
    H += delta2*np.kron(sx_o, (Sx@Sz + Sz@Sx))
    H += lam_perp*(np.kron(sx_o, Sx) + np.kron(Lz_o, Sy))
    H += xi_x*np.kron(sz_o, I3) + xi_y*np.kron(sx_o, I3)
    H += Bx*np.kron(I2, Sx)
    return H

D_SHAPE = np.eye(6, dtype=complex)
GAMMA_RAD = 0.0157
E_SPIN = {-1: np.array([1, 0, 0], complex), 0: np.array([0, 1, 0], complex),
          +1: np.array([0, 0, 1], complex)}
ORB_X = np.array([1, 0], complex); ORB_Y = np.array([0, 1], complex)

def legs(pair):
    ms_src, ms_det = pair
    return np.kron(ORB_Y, E_SPIN[ms_det]), np.kron(ORB_X, E_SPIN[ms_src])

def A0_of(H, z=0.0):
    return GAMMA_RAD*np.eye(6) + 1j*(H - z*np.eye(6))

def moments(H, pair, nmax, z=0.0):
    p, c = legs(pair)
    X = -np.linalg.solve(D_SHAPE, A0_of(H, z))
    nu = np.linalg.solve(D_SHAPE, c)
    out, v = [], nu.copy()
    for _ in range(nmax+1):
        out.append(p.conj()@v); v = X@v
    return np.array(out)

def kernel(H, pair, Gammas, z=0.0):
    p, c = legs(pair)
    A0m = A0_of(H, z)
    return np.array([p.conj()@np.linalg.solve(G*D_SHAPE + A0m, c) for G in Gammas])
