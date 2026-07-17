"""
Gate F5-B: full GKSL (density-matrix) treatment of a 2g+3e model.

Motivation: every result so far (Proposition A, the F4 floor theorem, the
2g+2e/2g+3e no-go, F5-A) was derived from a REDUCED coherence-block
generator A(delta) acting only on the coherence variables (optical
coherences + ground coherence), not from a genuine full Liouvillian that
also carries POPULATIONS. The reduced-block derivation implicitly assumes
that eliminating populations adiabatically cannot inject anything that
breaks passivity. This script checks that assumption directly by solving
the full N=5 (2 ground + 3 excited) density-matrix master equation:

  1. build the full GKSL Liouvillian at a given control strength Omega_c
     (probe off), find the unique control-dressed steady state rho_0;
  2. linearize in a weak probe field: rho(t) = rho_0 + delta_rho e^{-i
     delta t} + h.c., solve (L_0 - i*delta) delta_rho = -i [H_p, rho_0]
     for the coherence sourced by the probe drive on the DRESSED steady
     state (not the bare ground state -- this is the "population channel"
     the reduced-block calculation eliminates);
  3. evaluate the MATCHED response chi_full(delta) = Tr[mu_p^dagger
     delta_rho] (same dipole operator for source and readout -- the
     standard transmission/absorption observable) and check whether its
     real part can ever go negative or hit an exact zero;
  4. as a control, also evaluate two UNMATCHED observables (a Raman
     ground-coherence readout Tr[|g2><g1| delta_rho], and a
     cross-excited-state readout Tr[|e2><e1| delta_rho]) to see whether a
     natural (non-rotated-detector) unmatched escape exists inside the
     full physical model.

No approximation beyond weak-probe linearization and the standard
rotating-frame/RWA GKSL setup is made; populations, coherences, and their
coupling are all kept.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ---- basis: 0,1 = g1,g2 ; 2,3,4 = e1,e2,e3 -----------------------------
NG, NE = 2, 3
DIM = NG + NE


def basis_op(i, j):
    op = np.zeros((DIM, DIM), dtype=complex)
    op[i, j] = 1.0
    return op


def build_hamiltonian(Delta, Omega_p, Omega_c, dp, dc, Jmix):
    """Delta: excited detunings (len 3). Omega_p: probe Rabi (couples g1-e_j
    via dp, complex len-3). Omega_c: control Rabi (couples g2-e_j via dc).
    Jmix: Hermitian 3x3 real mixing among excited states."""
    H = np.zeros((DIM, DIM), dtype=complex)
    for j in range(NE):
        H[NG + j, NG + j] = Delta[j]
    for j in range(NE):
        H[NG + j, 0] += Omega_p * dp[j] / 2
        H[0, NG + j] += np.conj(Omega_p * dp[j] / 2)
        H[NG + j, 1] += Omega_c * dc[j] / 2
        H[1, NG + j] += np.conj(Omega_c * dc[j] / 2)
    for i in range(NE):
        for k in range(i + 1, NE):
            H[NG + i, NG + k] += Jmix[i, k]
            H[NG + k, NG + i] += Jmix[i, k]
    return H


def build_jumps(gamma_rad, gamma_g, gamma_dephase_e, gamma_relax=0.0):
    """gamma_rad[j][g]: decay rate e_j -> g_g. gamma_g: ground-coherence
    dephasing rate (pure dephasing between g1,g2). gamma_dephase_e: pure
    dephasing rate for each excited state (diagonal, models e.g. phonon
    dephasing). gamma_relax: symmetric ground-state population relaxation
    (T1-like spin-lattice relaxation g1<->g2) -- without this, Omega_p=0
    steady state funnels ALL population into g1 (a trivial control-dark
    trap with zero coherence), leaving nothing for the control field to
    dress; a small relaxation keeps both ground populations (and hence
    control-induced coherences) nonzero in steady state."""
    jumps = []
    for j in range(NE):
        for g in range(NG):
            rate = gamma_rad[j][g]
            if rate > 0:
                jumps.append(np.sqrt(rate) * basis_op(g, NG + j))
    if gamma_g > 0:
        Lg = np.sqrt(gamma_g) * (basis_op(0, 0) - basis_op(1, 1)) / np.sqrt(2)
        jumps.append(Lg)
    for j in range(NE):
        if gamma_dephase_e[j] > 0:
            Lj = np.sqrt(gamma_dephase_e[j]) * basis_op(NG + j, NG + j)
            jumps.append(Lj)
    if gamma_relax > 0:
        jumps.append(np.sqrt(gamma_relax) * basis_op(1, 0))  # g1 -> g2
        jumps.append(np.sqrt(gamma_relax) * basis_op(0, 1))  # g2 -> g1
    return jumps


def liouvillian(H, jumps):
    """Vectorized Liouvillian acting on vec(rho) (column-major/Fortran vec
    convention, vec(A X B) = (B^T ⊗ A) vec(X)):
    L = -i(I⊗H - H^T⊗I) + sum_k [conj(L_k)⊗L_k
    - 1/2 (I⊗L_k^dagger L_k + (L_k^dagger L_k)^T⊗I)]."""
    I = np.eye(DIM, dtype=complex)
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for Lk in jumps:
        LkdLk = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk) - 0.5 * (np.kron(I, LkdLk) + np.kron(LkdLk.T, I))
    return L


def steady_state(L):
    """Null vector of L (via smallest-|eigenvalue| eigenvector), normalized
    to trace 1 and Hermitian."""
    w, V = np.linalg.eig(L)
    idx = np.argmin(np.abs(w))
    vec = V[:, idx]
    rho = vec.reshape(DIM, DIM, order="F")
    rho = (rho + rho.conj().T) / 2
    rho = rho / np.trace(rho).real
    return rho, w[idx]


def solve_linear_response(L0, H_p_op, rho0, delta):
    """Solve for the harmonic response rho1(t) = drho * exp(-i*delta*t) of
    d(rho1)/dt = L0[rho1] + S*exp(-i*delta*t), S = -i[H_p_op, rho0]
    (the source term from a weak probe of envelope exp(-i*delta*t) acting
    on the control-dressed steady state rho0). Substituting the ansatz
    gives -i*delta*drho = L0[drho] + S, i.e.
        (L0 + i*delta*I) drho = -S = i*[H_p_op, rho0].
    (Verified against the exact two-level Lorentzian chi=1/a1 in the
    isolated e1-g1 sanity check.)"""
    I = np.eye(DIM * DIM, dtype=complex)
    comm = H_p_op @ rho0 - rho0 @ H_p_op
    rhs = 1j * comm.reshape(-1, order="F")
    # L0 alone is singular (its null space is the steady state); at delta=0
    # exactly, Aop is singular too, but rhs is guaranteed in its range
    # (traceless commutator source), so use the least-squares/minimum-norm
    # solution instead of a direct solve.
    Aop = L0 + 1j * delta * I
    vec_drho, *_ = np.linalg.lstsq(Aop, rhs, rcond=1e-12)
    return vec_drho.reshape(DIM, DIM, order="F")


def run_point(Delta, Omega_c, dp, dc, Jmix, gamma_rad, gamma_g, gamma_dephase_e,
              delta_scan, gamma_relax=0.0):
    H0 = build_hamiltonian(Delta, 0.0, Omega_c, dp, dc, Jmix)  # control only
    jumps = build_jumps(gamma_rad, gamma_g, gamma_dephase_e, gamma_relax)
    L0 = liouvillian(H0, jumps)
    rho0, eig0 = steady_state(L0)

    # H_p(t) = Omega_p * H_p_unit * cos(omega t): only the resonant
    # (rotating-wave) coupling is retained, H_p_unit = (dp[0]/2)|e1><g1| + h.c.
    H_p_unit = basis_op(NG + 0, 0) * (dp[0] / 2) + basis_op(0, NG + 0) * (np.conj(dp[0]) / 2)

    # chi is a SINGLE coherence matrix element (the standard optical-Bloch
    # convention used throughout this project, e.g. chi_full = <e1|rho|g1>
    # in the 2g+2e closed form), NOT Tr[mu_full @ drho] with the full
    # Hermitian dipole -- the latter double-counts co-/counter-rotating
    # contributions of drho and corrupts the sign of Re(chi).
    matched_idx = (NG + 0, 0)          # <e1|drho|g1>: same leg as the probe
    raman_idx = (1, 0)                 # <g2|drho|g1>: ground coherence
    cross_idx = (NG + 1, NG + 0)       # <e2|drho|e1>: cross-excited coherence

    # Normalization: solving (L0-i*delta)drho = -i*[H_p_unit,rho0] introduces
    # a factor of -i*(dp0/2) on the source side relative to the "chi = 1/a"
    # convention used throughout this project (e.g. chi_cut=1/a1 in the
    # 2g+2e closed form has Re>=0 for a passive Lorentzian). Multiplying by
    # i/(dp0/2) removes that bookkeeping factor and restores the standard
    # convention where Re(chi) >= 0 for a passive two-level absorber.
    norm = 1j / (dp[0] / 2)

    results = []
    for delta in delta_scan:
        drho = solve_linear_response(L0, H_p_unit, rho0, delta)
        chi_matched = norm * drho[matched_idx]
        chi_raman = norm * drho[raman_idx]
        chi_cross = norm * drho[cross_idx]
        results.append({
            "delta": float(delta),
            "chi_matched": [float(chi_matched.real), float(chi_matched.imag)],
            "chi_raman_unmatched": [float(chi_raman.real), float(chi_raman.imag)],
            "chi_cross_unmatched": [float(chi_cross.real), float(chi_cross.imag)],
        })
    return results, rho0, eig0


def main():
    rng = np.random.default_rng(0)
    Delta = np.array([0.0, 3.0, -2.0])
    Omega_c = 0.6
    dp = np.array([1.0, 0.0, 0.3], dtype=complex)
    dc = np.array([0.7, 1.0, 0.4 * np.exp(1j * 0.9)], dtype=complex)
    Jmix = np.zeros((3, 3))
    Jmix[0, 1] = Jmix[1, 0] = 0.5
    Jmix[1, 2] = Jmix[2, 1] = 0.3
    gamma_rad = [[0.6, 0.4], [0.5, 0.5], [0.45, 0.55]]  # e_j -> [g1, g2] rates
    gamma_g = 0.02
    gamma_dephase_e = [0.05, 0.05, 0.05]

    delta_scan = np.linspace(-6, 6, 1201)
    results, rho0, eig0 = run_point(Delta, Omega_c, dp, dc, Jmix, gamma_rad,
                                     gamma_g, gamma_dephase_e, delta_scan,
                                     gamma_relax=0.01)

    re_matched = np.array([r["chi_matched"][0] for r in results])
    min_idx = np.argmin(re_matched)
    print("steady-state eigenvalue (should be ~0):", eig0)
    print("rho0 populations (g1,g2,e1,e2,e3):", np.real(np.diag(rho0)))
    print(f"min Re(chi_matched) = {re_matched.min():.6e} at delta={delta_scan[min_idx]:.4f}")
    print(f"any Re(chi_matched) < 0 (gain-like)? {bool((re_matched < -1e-10).any())}")
    print(f"n points with |chi_matched| < 1e-6: "
          f"{sum(1 for r in results if abs(complex(*r['chi_matched'])) < 1e-6)}")

    re_raman = np.array([r["chi_raman_unmatched"][0] for r in results])
    print(f"Raman-unmatched: min/max Re = {re_raman.min():.4e}/{re_raman.max():.4e} "
          f"(sign change: {bool(re_raman.min() < 0 < re_raman.max())})")

    report = {
        "params": {
            "Delta": Delta.tolist(), "Omega_c": Omega_c,
            "gamma_rad": gamma_rad, "gamma_g": gamma_g,
            "gamma_dephase_e": gamma_dephase_e,
        },
        "rho0_populations": np.real(np.diag(rho0)).tolist(),
        "steady_state_eigenvalue": [float(eig0.real), float(eig0.imag)],
        "min_Re_chi_matched": float(re_matched.min()),
        "min_Re_chi_matched_at_delta": float(delta_scan[min_idx]),
        "any_gain_in_matched": bool((re_matched < -1e-10).any()),
        "raman_unmatched_sign_change": bool(re_raman.min() < 0 < re_raman.max()),
        "raman_unmatched_min_max": [float(re_raman.min()), float(re_raman.max())],
        "scan_sample": results[::100],
    }
    out_path = RESULTS_DIR / "gate_F5B_full_gksl.json"
    out_path.write_text(json.dumps(report, indent=2))
    print("Report written to", out_path)


if __name__ == "__main__":
    main()
