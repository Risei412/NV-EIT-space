"""Phase D, run D3: the class belongs to the scaling path, not the material.

One physical 4-dimensional generator, two inequivalent decompositions
A(z,Gamma) = Gamma*D_path + B_path(z), calibrated to agree exactly at a
single reference point Gamma0:

  Path T ("temperature-like"): only the phonon/orbital-hopping channel is
    scaled; D_T = diag(0,0,d3,d4) is singular (protected block = indices
    0,1, the ground-coherence pair). The radiative rate is *frozen* into
    B_T at its Gamma0 value.

  Path U ("uniform broadening"): every dissipation channel (radiative and
    phonon) is scaled together; D_U = diag(gU,gU,d3,d4) is full rank.

Calibration identity (exact by construction, not fitted):
    Gamma0*D_T + B_T(z) == Gamma0*D_U + B_U(z)   for every z.

Sector cut: identical protected-block internal coupling kappa (indices
0,1) severed in the cut model, exactly as in model_oblique / model_protected.
"""

import numpy as np

GAMMA0 = 10.0
GU = 0.7        # radiative rate absorbed uniformly along path U
D3, D4 = 1.1, 0.9
E1, E2 = 1.3, 0.9
LEAK = 0.5
KAPPA_FULL = 0.6


def _B0(kappa):
    B = np.zeros((4, 4))
    B[0, 0], B[1, 1] = E1, E2
    B[0, 1] = B[1, 0] = kappa
    B[2, 2], B[3, 3] = D3, D4
    B[0, 2] = B[2, 0] = LEAK
    B[1, 3] = B[3, 1] = LEAK
    return B


C_VEC = np.array([1.0, 0.4, 0.18, -0.11])
P_VEC = np.array([1.0, -0.25, 0.13, 0.22])


def D_T():
    return np.diag([0.0, 0.0, D3, D4])


def D_U():
    return np.diag([GU, GU, D3, D4])


def B_T(kappa):
    """Radiative term frozen at its Gamma0 value (B0 already carries the
    fixed D3,D4 diagonal terms untouched -- only the missing Gamma0*D_U
    radiative contribution, absent from D_T, is added here so that
    Gamma0*D_T + B_T == Gamma0*D_U + B_U exactly)."""
    B = _B0(kappa)
    B[0, 0] += GAMMA0 * GU
    B[1, 1] += GAMMA0 * GU
    return B


def B_U(kappa):
    return _B0(kappa)


def calibration_check():
    A_T = GAMMA0 * D_T() + B_T(KAPPA_FULL)
    A_U = GAMMA0 * D_U() + B_U(KAPPA_FULL)
    return float(np.max(np.abs(A_T - A_U)))


def R_S_gamma(Gamma, path, kappa_full=KAPPA_FULL):
    if path == "T":
        D, Bfun = D_T(), B_T
    elif path == "U":
        D, Bfun = D_U(), B_U
    else:
        raise ValueError(path)
    A_full = Gamma * D + Bfun(kappa_full)
    A_cut = Gamma * D + Bfun(0.0)
    cf = np.conj(P_VEC) @ np.linalg.solve(A_full, C_VEC)
    cc = np.conj(P_VEC) @ np.linalg.solve(A_cut, C_VEC)
    return cf, cc, cf - cc


def chi_full_at_calibration():
    """chi_full evaluated at Gamma0 on both paths -- must agree exactly,
    confirming this is one physical system probed two ways."""
    cfT, _, _ = R_S_gamma(GAMMA0, "T")
    cfU, _, _ = R_S_gamma(GAMMA0, "U")
    return cfT, cfU


def theta_fan_D(theta):
    """D(theta) = cos(theta)*D_phonon_only + sin(theta)*D_rad_only, a
    one-parameter family of scaling directions (B held fixed at B_U, i.e.
    not recalibrated per theta -- this is a separate demonstration, class
    as a function of scaling-path *direction*, not a physical-point match).
    D(theta) = diag(sin(theta)*gU, sin(theta)*gU, cos(theta)*d3, cos(theta)*d4)
    is singular whenever sin(theta)=0 (path T, protected 0,1 block) or
    cos(theta)=0 (protected 2,3 block instead); full rank in between."""
    D_phonon = np.diag([0.0, 0.0, D3, D4])
    D_rad = np.diag([GU, GU, 0.0, 0.0])
    return np.cos(theta) * D_phonon + np.sin(theta) * D_rad


def R_S_gamma_fan(Gamma, theta, kappa_full=KAPPA_FULL):
    D = theta_fan_D(theta)
    A_full = Gamma * D + B_U(kappa_full)
    A_cut = Gamma * D + B_U(0.0)
    cf = np.conj(P_VEC) @ np.linalg.solve(A_full, C_VEC)
    cc = np.conj(P_VEC) @ np.linalg.solve(A_cut, C_VEC)
    return cf - cc
