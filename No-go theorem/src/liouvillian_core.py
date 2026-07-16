"""liouvillian_core.py -- generic dense Lindblad superoperator utilities.
Shared full-Liouvillian backbone for NV (nv_system.py) and group-IV
(group_iv_system.py). Row-stacking vectorization convention:
  vec(rho)[i*n+j] = rho[i,j],  L = -i(I(x)H - H^T(x)I) + sum_k D[L_k]
  D[L] = L(x)L^* - 0.5*(L^dag L)^T(x)I - 0.5*I(x)(L^dag L)
This matches d/dt vec(rho) = L @ vec(rho) for drho/dt = -i[H,rho]
  + sum_k (L_k rho L_k^dag - 0.5{L_k^dag L_k, rho}).
"""
import numpy as np

def vec(rho):
    return rho.reshape(-1)

def unvec(v, n):
    return v.reshape(n, n)

def liouvillian(H, Ls):
    n = H.shape[0]; I = np.eye(n, dtype=complex)
    L = -1j*(np.kron(I, H) - np.kron(H.T, I))
    for Lk in Ls:
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk) - 0.5*np.kron(I, LdL) - 0.5*np.kron(LdL.T, I)
    return L

def steady_state(L, n=None):
    """Null vector of L (smallest |eigenvalue|), normalized to trace 1.
    Returns (rho0_vec, residual_norm, spectral_gap)."""
    if n is None: n = int(round(np.sqrt(L.shape[0])))
    w, V = np.linalg.eig(L)
    idx = np.argsort(np.abs(w))
    w0 = w[idx[0]]; v0 = V[:, idx[0]]
    rho0 = unvec(v0, n)
    rho0 = 0.5*(rho0 + rho0.conj().T)  # hermitize (steady state must be)
    tr = np.trace(rho0)
    if abs(tr) < 1e-300: raise ValueError("degenerate/traceless steady state")
    rho0 = rho0 / tr
    residual = np.linalg.norm(L @ vec(rho0))
    gap = float(np.abs(w[idx[1]])) if len(idx) > 1 else np.inf
    return vec(rho0), float(residual), gap

def first_order(L, V, rho0_vec):
    """Solve L @ x = -V @ rho0_vec for the first-order (linear-response)
    correction, restricted to the traceless subspace (x has zero trace
    since rho0 already carries the full population). Uses least-squares
    on the (singular, one zero mode) L."""
    n = int(round(np.sqrt(L.shape[0])))
    b = -(V @ rho0_vec)
    x, *_ = np.linalg.lstsq(L, b, rcond=None)
    rho1 = unvec(x, n)
    tr = np.trace(rho1)
    rho1 = rho1 - (tr/n)*np.eye(n, dtype=complex)  # project out trace (gauge choice)
    x = vec(rho1)
    residual = np.linalg.norm(L @ x - b)
    return x, float(residual)

def commutator_super(Op):
    """Superoperator for -i[Op,.] acting alone (used to build probe drive V)."""
    n = Op.shape[0]; I = np.eye(n, dtype=complex)
    return -1j*(np.kron(I, Op) - np.kron(Op.T, I))
