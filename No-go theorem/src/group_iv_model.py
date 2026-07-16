"""group_iv_model.py -- SiV-/SnV- excited-manifold Hamiltonian, in the SAME
reduced-kernel formalism as nv_reduced_kernel.py (A_Gamma = Gamma*D + A0,
D=I, X=-A0, M_n = p^dag X^n c, K(Gamma) = p^dag (Gamma*D+A0)^-1 c), so that
NV and group-IV moment orders can be compared in one common computational
pipeline (moment_order_common_pipeline.py).

Physical model (D3d symmetry, orbital doublet (e_x,e_y) x spin-1/2):
  H = H_SO + H_strain + H_Z
  H_SO      = (Delta_e/2) * kron(sz_orb, sz_spin)      (excited spin-orbit)
  H_strain  = xi_x*kron(sz_orb,I_spin) + xi_y*kron(sx_orb,I_spin)
              (xi_x,xi_y map to d(exx-eyy)+f*ezx , -2d*exy+f*eyz)
  H_Z       = Bx*kron(I_orb, sx_spin)                  (transverse spin Zeeman)

Parameters (fixed, from SiV_SnV_phonon_AIC_parameters.md, GHz, ordinary freq):
  SiV-: Delta_e/2pi = 255 GHz,  d_e=1.8 PHz/strain, f_e=-0.720 PHz/strain
  SnV-: Delta_e/2pi = 3000 GHz, d_e=0.921 PHz/strain, f_e=-2.592 PHz/strain

Dipole legs: "same-spin orbital-Lambda" -- unlike NV (probe/control differ in
ORBITAL branch AND spin, giving p^dag c = 0 identically because <OX|OY>=0),
group-IV probe/control share the SAME excited ket at zeroth order (D3d E(x)E
dipole tensor allows both ground orbital branches to reach a common bright
excited sublevel with the same spin projection), giving M0 = p^dag c != 0.
theta is a leg-mixing angle representing the (unmeasured) ratio of cross-
branch to same-branch dipole matrix elements; theta=0 is the maximal-overlap
case (pure same-ket Lambda) used as the default/representative claim.
"""
import numpy as np

PARAMS = {
    'SiV': dict(Delta_e=255.0, d_e=1.8, f_e=-0.720),      # GHz, PHz/strain
    'SnV': dict(Delta_e=3000.0, d_e=0.921, f_e=-2.592),
}

sz_orb = np.array([[1, 0], [0, -1]], complex)
sx_orb = np.array([[0, 1], [1, 0]], complex)
sz_spin = np.array([[1, 0], [0, -1]], complex)
sx_spin = np.array([[0, 1], [1, 0]], complex)
I2 = np.eye(2, dtype=complex)
ORB_X = np.array([1, 0], complex); ORB_Y = np.array([0, 1], complex)
SPIN_UP = np.array([1, 0], complex); SPIN_DN = np.array([0, 1], complex)

def H_groupIV(material, xi_x=0.0, xi_y=0.0, Bx=0.0):
    Delta_e = PARAMS[material]['Delta_e']
    H = (Delta_e/2.0)*np.kron(sz_orb, sz_spin)
    H += xi_x*np.kron(sz_orb, I2) + xi_y*np.kron(sx_orb, I2)
    H += Bx*np.kron(I2, sx_spin)
    return H

D_SHAPE = np.eye(4, dtype=complex)
GAMMA_RAD = 0.0157  # same order-of-magnitude radiative HWHM convention as NV

def legs(theta=0.0, spin=SPIN_UP):
    """probe -> e_x branch; control -> cos(theta) e_x + sin(theta) e_y branch,
    both with the SAME spin projection (spin-conserving optical dipole)."""
    c = np.kron(ORB_X, spin)
    p = np.cos(theta)*np.kron(ORB_X, spin) + np.sin(theta)*np.kron(ORB_Y, spin)
    return p, c

def A0_of(H, z=0.0):
    return GAMMA_RAD*np.eye(4) + 1j*(H - z*np.eye(4))

def moments(H, nmax, theta=0.0, z=0.0):
    p, c = legs(theta)
    X = -np.linalg.solve(D_SHAPE, A0_of(H, z))
    nu = np.linalg.solve(D_SHAPE, c)
    out, v = [], nu.copy()
    for _ in range(nmax+1):
        out.append(p.conj()@v); v = X@v
    return np.array(out)

def kernel(H, Gammas, theta=0.0, z=0.0):
    p, c = legs(theta)
    A0m = A0_of(H, z)
    return np.array([p.conj()@np.linalg.solve(G*D_SHAPE + A0m, c) for G in Gammas])

# ---- phonon-limited orbital relaxation, Bose-law form (Part I/II of
# SiV_SnV_phonon_AIC_parameters.md; NOT a single T^n over the full range).
# Normalization constants (NORM_GHZ) are schematic/representative -- fitted
# per-sample constants (A1,A5,C) are not given in the source material -- and
# are chosen only so that Gamma_orb(300K) lands at the order of magnitude of
# reported room-temperature group-IV orbital relaxation rates (GHz for SiV,
# slower for SnV due to its much larger splitting).
H_PLANCK, KB = 6.62607015e-34, 1.380649e-23
NORM_GHZ = {'SiV': 3.0e-3, 'SnV': 3.0e-3}  # A_1-equivalent, GHz per (rad/ns)^3-ish unit

def n_bose(omega_GHz, T):
    if T <= 0: return 0.0
    x = H_PLANCK*(omega_GHz*1e9)/(KB*T)
    return 1.0/np.expm1(x) if x > 1e-12 else 1e12

def gamma_orb_GHz(material, T):
    """One-phonon Bose-law orbital relaxation rate, Gamma_orb(T) =
    A*omega_Delta^3*[2 n_B(omega_Delta,T)+1] (Part I Sec 4.1 / Part II Sec 7.1)."""
    Delta_e = PARAMS[material]['Delta_e']  # GHz, ordinary frequency
    nB = n_bose(Delta_e, T)
    return NORM_GHZ[material]*Delta_e**3*(2*nB + 1)*1e-6  # scaled to GHz-order at 300K

def thermal_regime(material, T):
    """k_B T / (h * Delta_e): >>1 high-T (power-law) limit, <<1 low-T (activated)."""
    Delta_e = PARAMS[material]['Delta_e']
    return (KB*T)/(H_PLANCK*Delta_e*1e9)
