"""Phase P, physical model (paper3_smrt_numerical_plan.md, P0).

A 5-level "diamond" Lindblad model realizing the SMRT interference-controlled
exponent promotion nu: 3 -> 4 in a genuine CPTP open quantum system.

States: |1> ground, |2> probe-excited, |3>,|4> intermediates, |5> readout.
Coherent paths from |2> to |5>: 2-3-5 (couplings J23, J35) and 2-4-5
(couplings J24, J45 e^{i phi}), plus a weak direct route J25 and a weak
cross-coupling J34 that opens the 3-hop path.

Sector S = {J35, J45 e^{i phi}} -- the two couplings that close *both*
coherent branches into the readout state |5>. Cutting S (J35 = J45 = 0)
removes both routes simultaneously, so the leading moment of R_S =
chi_full - chi_cut is the coherent SUM of both branch amplitudes -- the
interference structure the SMRT strategy requires (see plan P0.1).

Two representations are provided:

  1. The weak-probe REDUCED linear system on the optical coherences
     x = (rho21, rho31, rho41, rho51), which is exactly of SMRT form
     A(z) x = c,  A(z) = Gamma*D + B(z),
     and is graph-isomorphic to the already-validated abstract D2 model
     (model_cancel.build_promotion_model): site 0 <-> rho21,
     sites {1,2} <-> {rho31, rho41}, site 3 <-> rho51, v12 <-> J34.

  2. The FULL vectorized Liouvillian (5x5 Hilbert space, 25x25 vectorized,
     24x24 on the traceless subspace) with an implicit-differentiation
     linear-response solver (no finite-epsilon subtraction), used for
     Calculation H (reduced vs. full) and the second-order readout rho55
     in Calculation G.
"""

import numpy as np
import sympy as sp

# ----------------------------------------------------------------------
# Symbols (exact-rational parameters throughout; see plan P0.5)
# ----------------------------------------------------------------------
Gamma, z = sp.symbols("Gamma z", real=True)
d2, d3, d4, d5 = sp.symbols("d2 d3 d4 d5", positive=True)
D2_, D3_, D4_, D5_ = sp.symbols("Delta2 Delta3 Delta4 Delta5", real=True)  # detunings
J23, J24, J25, J34, J35, J45 = sp.symbols("J23 J24 J25 J34 J35 J45", real=True)
phi = sp.symbols("phi", real=True)
I = sp.I

DEFAULT_PARAMS = dict(
    d2=sp.Rational(1, 1), d3=sp.Rational(13, 10), d4=sp.Rational(17, 10), d5=sp.Rational(21, 10),
    Delta2=0, Delta3=0, Delta4=0, Delta5=0,
    J23=sp.Rational(1, 1), J24=sp.Rational(1, 1), J25=sp.Rational(1, 10),
    J34=sp.Rational(1, 20), J35=sp.Rational(3, 5),
)


# ----------------------------------------------------------------------
# Reduced 4x4 system: A(z) x = c,  A = Gamma*D + B(z)
# Basis order: (rho21, rho31, rho41, rho51)
# ----------------------------------------------------------------------
def D_matrix(params=None):
    p = dict(DEFAULT_PARAMS); p.update(params or {})
    return sp.diag(p["d2"], p["d3"], p["d4"], p["d5"]) / 2


def B_matrix(z_val, J45_val, phi_val, cut=False, params=None):
    """B(z) for the reduced system. If cut=True, J35=J45=0 (sector removed)."""
    p = dict(DEFAULT_PARAMS); p.update(params or {})
    J35v = 0 if cut else p["J35"]
    J45v = 0 if cut else J45_val
    phiv = phi_val
    M = sp.zeros(4, 4)
    # diagonal: i*(Delta_j - z)
    M[0, 0] = I * (p["Delta2"] - z_val)
    M[1, 1] = I * (p["Delta3"] - z_val)
    M[2, 2] = I * (p["Delta4"] - z_val)
    M[3, 3] = I * (p["Delta5"] - z_val)
    # off-diagonal couplings (Hermitian H -> i*J symmetric structure)
    M[0, 1] = I * p["J23"]; M[1, 0] = I * p["J23"]
    M[0, 2] = I * p["J24"]; M[2, 0] = I * p["J24"]
    M[0, 3] = I * p["J25"]; M[3, 0] = I * p["J25"]
    M[1, 2] = I * p["J34"]; M[2, 1] = I * p["J34"]
    M[1, 3] = I * J35v; M[3, 1] = I * J35v
    M[2, 3] = I * J45v * sp.exp(I * phiv); M[3, 2] = I * J45v * sp.exp(-I * phiv)
    return M


def c_vector(eps=sp.Rational(1, 1)):
    v = sp.zeros(4, 1)
    v[0, 0] = -I * eps
    return v


def p_vector(which="rho51"):
    idx = {"rho21": 0, "rho31": 1, "rho41": 2, "rho51": 3}[which]
    v = sp.zeros(4, 1)
    v[idx, 0] = 1
    return v


def A_matrix(Gamma_val, z_val, J45_val, phi_val, cut=False, params=None):
    return Gamma_val * D_matrix(params) + B_matrix(z_val, J45_val, phi_val, cut=cut, params=params)


def R_S_symbolic(Gamma_val, z_val, J45_val, phi_val, params=None, readout="rho51"):
    """R_S(Gamma,z) = p^dag [A_full^-1 - A_cut^-1] c, exact symbolic."""
    Af = A_matrix(Gamma_val, z_val, J45_val, phi_val, cut=False, params=params)
    Ac = A_matrix(Gamma_val, z_val, J45_val, phi_val, cut=True, params=params)
    c = c_vector()
    p = p_vector(readout)
    xf = Af.solve(c)
    xc = Ac.solve(c)
    return (p.T * (xf - xc))[0, 0]


def chi_full_symbolic(Gamma_val, z_val, J45_val, phi_val, params=None, readout="rho51"):
    Af = A_matrix(Gamma_val, z_val, J45_val, phi_val, cut=False, params=params)
    xf = Af.solve(c_vector())
    return (p_vector(readout).T * xf)[0, 0]


def chi_cut_symbolic(Gamma_val, z_val, J45_val, phi_val, params=None, readout="rho51"):
    Ac = A_matrix(Gamma_val, z_val, J45_val, phi_val, cut=True, params=params)
    xc = Ac.solve(c_vector())
    return (p_vector(readout).T * xc)[0, 0]


# ----------------------------------------------------------------------
# Cancellation condition (plan P0.4): eta * e^{i phi} = -1
# eta = (J23 J35 d4) / (J24 J45 d3)
# Closed form: phi_star = pi, eta_star = 1  =>  J45_star = -(J23 J35 d4)/(J24 d3)
# (real J45, phi=pi realizes the e^{i phi} = -1 branch)
# ----------------------------------------------------------------------
def eta_of(J45_val, params=None):
    p = dict(DEFAULT_PARAMS); p.update(params or {})
    return (p["J23"] * p["J35"] * p["d4"]) / (p["J24"] * J45_val * p["d3"])


def J45_star(params=None, phi_val=sp.pi):
    """Solve eta * e^{i phi} = -1 for J45 at fixed phi (real J45 branch)."""
    p = dict(DEFAULT_PARAMS); p.update(params or {})
    # eta * exp(i phi) = -1  =>  eta = -exp(-i phi)
    eta_target = -sp.exp(-I * phi_val)
    assert sp.im(sp.simplify(eta_target)) == 0 or phi_val in (0, sp.pi), \
        "phi must give real eta for a real J45 solution"
    eta_target = sp.re(sp.simplify(eta_target))
    J45v = (p["J23"] * p["J35"] * p["d4"]) / (p["J24"] * eta_target * p["d3"])
    return sp.nsimplify(J45v)


# ----------------------------------------------------------------------
# Moment ladder (Theorem II): mu_k = p^dag (D^-1 B(z))^k D^-1 c
# on the DIFFERENCE (full - cut), evaluated exactly symbolically.
# ----------------------------------------------------------------------
def moments_difference(z_val, J45_val, phi_val, kmax=4, params=None, readout="rho51"):
    Dm = D_matrix(params)
    Dinv = Dm.inv()
    Bf = B_matrix(z_val, J45_val, phi_val, cut=False, params=params)
    Bc = B_matrix(z_val, J45_val, phi_val, cut=True, params=params)
    c = c_vector()
    p = p_vector(readout)

    Xf = Dinv * Bf
    Xc = Dinv * Bc
    vf = Dinv * c
    vc = Dinv * c
    out = []
    for _ in range(kmax):
        mf = (p.T * vf)[0, 0]
        mc = (p.T * vc)[0, 0]
        out.append(sp.simplify(mf - mc))
        vf = Xf * vf
        vc = Xc * vc
    return out


# ----------------------------------------------------------------------
# Polynomial certificate (Theorem 8.1): nu_cert = deg_Gamma Q - deg_Gamma N
# for R_S = N(Gamma,z)/Q(Gamma,z), assembled via a common denominator of the
# two rational functions chi_full and chi_cut.
# ----------------------------------------------------------------------
def certificate_R_S(z_val, J45_val, phi_val, params=None, readout="rho51"):
    Af = A_matrix(Gamma, z_val, J45_val, phi_val, cut=False, params=params)
    Ac = A_matrix(Gamma, z_val, J45_val, phi_val, cut=True, params=params)
    c = c_vector()
    p = p_vector(readout)
    Qf = sp.together(Af.det())
    Qc = sp.together(Ac.det())
    Nf = sp.together((p.T * Af.adjugate() * c)[0, 0])
    Nc = sp.together((p.T * Ac.adjugate() * c)[0, 0])
    # R_S = Nf/Qf - Nc/Qc = (Nf*Qc - Nc*Qf) / (Qf*Qc)
    N = sp.expand(Nf * Qc - Nc * Qf)
    Q = sp.expand(Qf * Qc)
    Qp = sp.Poly(Q, Gamma)
    if N == 0:
        return None, Qp, None
    Np = sp.Poly(N, Gamma)
    nu = Qp.degree() - Np.degree()
    return int(nu), Qp, Np


# ----------------------------------------------------------------------
# Fast numeric (float / mpmath) evaluation for Gamma sweeps
# ----------------------------------------------------------------------
def lambdify_R_S(z_val, J45_val, phi_val, params=None, readout="rho51"):
    """Return callables R_S(Gamma) in float64 and in mpmath, built once from
    the symbolic 4x4 solve (fast repeated evaluation over a Gamma grid)."""
    import mpmath as mp

    Bf = B_matrix(z_val, J45_val, phi_val, cut=False, params=params)
    Bc = B_matrix(z_val, J45_val, phi_val, cut=True, params=params)
    Dm = D_matrix(params)
    c = c_vector()
    p = p_vector(readout)

    Bf_n = np.array(Bf.evalf(), dtype=complex)
    Bc_n = np.array(Bc.evalf(), dtype=complex)
    D_n = np.array(Dm.evalf(), dtype=float)
    c_n = np.array(c.evalf(), dtype=complex).flatten()
    p_n = np.array(p.evalf(), dtype=complex).flatten()

    def R_S_f64(Gval):
        Af = Gval * D_n + Bf_n
        Ac = Gval * D_n + Bc_n
        xf = np.linalg.solve(Af, c_n)
        xc = np.linalg.solve(Ac, c_n)
        cond = max(np.linalg.cond(Af), np.linalg.cond(Ac))
        return complex(p_n @ xf - p_n @ xc), float(cond)

    def R_S_mp(Gval, dps=50):
        mp.mp.dps = dps
        Af = [[mp.mpc(Gval) * complex(D_n[i, j]) + complex(Bf_n[i, j]) for j in range(4)] for i in range(4)]
        Ac = [[mp.mpc(Gval) * complex(D_n[i, j]) + complex(Bc_n[i, j]) for j in range(4)] for i in range(4)]
        Af_m = mp.matrix(Af); Ac_m = mp.matrix(Ac)
        c_m = mp.matrix([complex(x) for x in c_n])
        p_m = [complex(x) for x in p_n]
        xf = mp.lu_solve(Af_m, c_m)
        xc = mp.lu_solve(Ac_m, c_m)
        rf = sum(p_m[i] * xf[i] for i in range(4))
        rc = sum(p_m[i] * xc[i] for i in range(4))
        return rf - rc

    return R_S_f64, R_S_mp


# ----------------------------------------------------------------------
# CPTP check
# ----------------------------------------------------------------------
def cptp_check(params=None):
    """All decay rates Gamma*d_j/2 >= 0 for Gamma >= 0 (jump operators
    L_j = sqrt(Gamma d_j) |1><j|); simply requires d_j > 0."""
    p = dict(DEFAULT_PARAMS); p.update(params or {})
    ds = [p["d2"], p["d3"], p["d4"], p["d5"]]
    return all(float(d) > 0 for d in ds)


# ----------------------------------------------------------------------
# Full vectorized Liouvillian (5-level Hilbert space)
# ----------------------------------------------------------------------
N_LEVELS = 5


def hamiltonian_numeric(J45_val, phi_val, cut=False, params=None, eps=0.0, z_val=0.0):
    """5x5 Hamiltonian in the rotating frame, including the weak probe on
    |1><2| with amplitude eps (real). z_val enters via a frame shift Delta_j
    -> Delta_j - z applied uniformly (probe detuning)."""
    p = dict(DEFAULT_PARAMS); p.update(params or {})
    J35v = 0.0 if cut else float(p["J35"])
    J45v = 0.0 if cut else float(J45_val)
    H = np.zeros((5, 5), dtype=complex)
    Deltas = [0.0, float(p["Delta2"]) - z_val, float(p["Delta3"]) - z_val,
              float(p["Delta4"]) - z_val, float(p["Delta5"]) - z_val]
    for j in range(1, 5):
        H[j, j] = Deltas[j]
    H[0, 1] = eps; H[1, 0] = eps  # probe |1><2| + h.c.
    H[1, 2] = float(p["J23"]); H[2, 1] = float(p["J23"])
    H[1, 3] = float(p["J24"]); H[3, 1] = float(p["J24"])
    H[1, 4] = float(p["J25"]); H[4, 1] = float(p["J25"])
    H[2, 3] = float(p["J34"]); H[3, 2] = float(p["J34"])
    H[2, 4] = J35v; H[4, 2] = J35v
    phi_f = float(phi_val)
    H[3, 4] = J45v * np.exp(1j * phi_f); H[4, 3] = J45v * np.exp(-1j * phi_f)
    return H


def jump_operators(Gamma_val, params=None):
    p = dict(DEFAULT_PARAMS); p.update(params or {})
    ds = [float(p["d2"]), float(p["d3"]), float(p["d4"]), float(p["d5"])]
    Ls = []
    for j, d in zip(range(1, 5), ds):
        L = np.zeros((5, 5), dtype=complex)
        L[0, j] = np.sqrt(Gamma_val * d)
        Ls.append(L)
    return Ls


def liouvillian(H, Ls):
    """Vectorized Liouvillian L such that vec(drho/dt) = L @ vec(rho),
    row-major vec (vec(rho)_ij -> index i*n+j), using the standard
    Lindblad vectorization: L = -i(H_kron_I - I_kron_H^T) + sum_k D[L_k]."""
    n = H.shape[0]
    Ivec = np.eye(n)
    Lsup = -1j * (np.kron(H, Ivec) - np.kron(Ivec, H.T))
    for Lk in Ls:
        Lkd = Lk.conj().T
        Lsup += np.kron(Lk, Lk.conj())
        Lsup -= 0.5 * (np.kron(Lkd @ Lk, Ivec) + np.kron(Ivec, (Lkd @ Lk).T))
    return Lsup


def steady_state(Gamma_val, J45_val, phi_val, cut=False, params=None, eps=0.0, z_val=0.0):
    """Solve L @ vec(rho) = 0 with trace(rho)=1 by replacing one row."""
    H = hamiltonian_numeric(J45_val, phi_val, cut=cut, params=params, eps=eps, z_val=z_val)
    Ls = jump_operators(Gamma_val, params=params)
    L = liouvillian(H, Ls)
    n = H.shape[0]
    M = L.copy()
    b = np.zeros(n * n, dtype=complex)
    # replace the equation for element (0,0) with the trace condition
    trace_row = np.zeros(n * n, dtype=complex)
    for i in range(n):
        trace_row[i * n + i] = 1.0
    M[0, :] = trace_row
    b[0] = 1.0
    vecrho = np.linalg.solve(M, b)
    return vecrho.reshape(n, n)


def implicit_linear_response(Gamma_val, J45_val, phi_val, cut=False, params=None, z_val=0.0):
    """O(eps) linear response without finite-eps subtraction: solve
    L0 @ vec(rho1) = -i vec([V_probe, rho0]) on the traceless subspace,
    where rho0 is the eps=0 steady state (= |1><1| here) and V_probe =
    |1><2| + h.c. Returns rho1 (5x5, complex, traceless)."""
    H0 = hamiltonian_numeric(J45_val, phi_val, cut=cut, params=params, eps=0.0, z_val=z_val)
    Ls = jump_operators(Gamma_val, params=params)
    L0 = liouvillian(H0, Ls)
    n = H0.shape[0]

    rho0 = np.zeros((n, n), dtype=complex)
    rho0[0, 0] = 1.0  # ground-state steady state at eps=0 (verified below)

    Vp = np.zeros((n, n), dtype=complex)
    Vp[0, 1] = 1.0; Vp[1, 0] = 1.0
    # steady state (rho0+eps*rho1) of dot{rho}=L0(rho)-i*eps[V,rho] to O(eps):
    # L0 rho1 = i[V,rho0]  (not -i[V,rho0])
    rhs = 1j * (Vp @ rho0 - rho0 @ Vp)
    rhs_vec = rhs.reshape(-1)

    # project onto traceless subspace: solve L0 x = rhs with an extra
    # trace(x)=0 constraint replacing a redundant row (L0 is singular,
    # rank n^2-1, with kernel spanned by steady states)
    M = L0.copy()
    b = rhs_vec.copy()
    trace_row = np.zeros(n * n, dtype=complex)
    for i in range(n):
        trace_row[i * n + i] = 1.0
    M[0, :] = trace_row
    b[0] = 0.0
    x = np.linalg.solve(M, b)
    return x.reshape(n, n)


def verify_rho0_is_steady(Gamma_val, J45_val, phi_val, cut=False, params=None):
    """Check |1><1| is indeed the eps=0 steady state (all decay returns to
    ground, no coherent coupling touches |1> except via the probe)."""
    H0 = hamiltonian_numeric(J45_val, phi_val, cut=cut, params=params, eps=0.0)
    Ls = jump_operators(Gamma_val, params=params)
    L0 = liouvillian(H0, Ls)
    rho0 = np.zeros((5, 5), dtype=complex); rho0[0, 0] = 1.0
    resid = L0 @ rho0.reshape(-1)
    return float(np.max(np.abs(resid)))
