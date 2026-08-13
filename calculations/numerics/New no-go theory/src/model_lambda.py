"""Phase A: standard Lambda model (priorities document, Sec. 2.1).

M(Delta, delta) = [[gamma31 - i*Delta,    i*Omega_c/2      ],
                    [i*conj(Omega_c)/2,   gamma21 - i*delta]]

chi_full = p^dag M^-1 c with c = p = e1 (probe injection/readout on the
optical coherence rho31), i.e. chi_full = (M^-1)_{11}.

chi_cut is the frozen-source cut: same D (gamma31, gamma21), same c, p,
only the Omega_c coupling between the optical and ground-coherence sectors
is severed (M_cut is M with the off-diagonal set to zero). This is *not*
Omega_c -> 0 in the full spectrum (that would also change the diagonal
control-induced light shift/dephasing contributions were they present);
here, since Omega_c enters only off-diagonally, cutting the coupling and
setting Omega_c=0 in this minimal model coincide numerically, but the two
are computed via structurally different routes (direct formula vs.
generic M-zeroing) as a self-consistency check.
"""

import numpy as np


def M_matrix(Delta, delta, Omega_c, gamma31, gamma21):
    return np.array([
        [gamma31 - 1j * Delta, 1j * Omega_c / 2],
        [1j * np.conj(Omega_c) / 2, gamma21 - 1j * delta],
    ], dtype=complex)


def chi_full_analytic(Delta, delta, Omega_c, gamma31, gamma21):
    num = gamma21 - 1j * delta
    den = (gamma31 - 1j * Delta) * (gamma21 - 1j * delta) + abs(Omega_c) ** 2 / 4
    return num / den


def chi_cut_analytic(Delta, delta, Omega_c, gamma31, gamma21):
    return 1.0 / (gamma31 - 1j * Delta)


def chi_full_numeric(Delta, delta, Omega_c, gamma31, gamma21):
    M = M_matrix(Delta, delta, Omega_c, gamma31, gamma21)
    c = np.array([1.0, 0.0], dtype=complex)
    p = np.array([1.0, 0.0], dtype=complex)
    return np.conj(p) @ np.linalg.solve(M, c)


def chi_cut_numeric(Delta, delta, Omega_c, gamma31, gamma21):
    """Sector cut: sever the Omega_c coupling that links the ground
    coherence sector to the optical response, keep D, c, p fixed."""
    M = M_matrix(Delta, delta, Omega_c, gamma31, gamma21)
    M_cut = M.copy()
    M_cut[0, 1] = 0.0
    M_cut[1, 0] = 0.0
    c = np.array([1.0, 0.0], dtype=complex)
    p = np.array([1.0, 0.0], dtype=complex)
    return np.conj(p) @ np.linalg.solve(M_cut, c)


def full_cut_R(Delta, delta, Omega_c, gamma31, gamma21):
    cf = chi_full_analytic(Delta, delta, Omega_c, gamma31, gamma21)
    cc = chi_cut_analytic(Delta, delta, Omega_c, gamma31, gamma21)
    return cf, cc, cf - cc


def self_consistency_check(n=200, seed=0):
    """Analytic vs. generic numeric (M-zeroing) cross-check, Phase A step 2."""
    rng = np.random.default_rng(seed)
    max_err_full = 0.0
    max_err_cut = 0.0
    for _ in range(n):
        Delta = rng.uniform(-10, 10)
        delta = rng.uniform(-10, 10)
        Omega_c = rng.uniform(0.01, 10) * np.exp(1j * rng.uniform(0, 2 * np.pi))
        gamma31 = rng.uniform(0.5, 2.0)
        gamma21 = rng.uniform(1e-4, 1.0)
        cf_a = chi_full_analytic(Delta, delta, Omega_c, gamma31, gamma21)
        cf_n = chi_full_numeric(Delta, delta, Omega_c, gamma31, gamma21)
        cc_a = chi_cut_analytic(Delta, delta, Omega_c, gamma31, gamma21)
        cc_n = chi_cut_numeric(Delta, delta, Omega_c, gamma31, gamma21)
        max_err_full = max(max_err_full, abs(cf_a - cf_n))
        max_err_cut = max(max_err_cut, abs(cc_a - cc_n))
    return max_err_full, max_err_cut


if __name__ == "__main__":
    ef, ec = self_consistency_check()
    print(f"[Phase A self-check] max|chi_full analytic-numeric| = {ef:.2e}, "
          f"max|chi_cut analytic-numeric| = {ec:.2e} -> "
          f"{'PASS' if max(ef, ec) < 1e-10 else 'FAIL'}")
