"""Numerical verification of the three central theorems.

Checks, with explicit tolerances:
  1. Theorem I  (Class I):  engineered Krylov-orthogonal model -> all moments
     mu_k = l^dag M^k r vanish and delta_chi_S(z) = 0 at machine precision
     for every test frequency.
  2. Theorem II (Class II): full-rank-D model with sector-graph distance d ->
     first nonzero moment at order m = d, and the fitted slope of
     log|F_Gamma| vs log Gamma equals -(m+1) = -nu_diss.
  3. Theorem III (Class III): singular D with semisimple kernel ->
     F_Gamma plateaus at F0 = p^dag P B_P^{-1} P c, and the subleading
     correction decays as 1/Gamma.
  4. Counterexample lemma: Pc != 0 and p^dag P != 0 with F0 = 0.

Run:  python three_theorem_verification.py
"""

import numpy as np

rng = np.random.default_rng(7)
TOL = 1e-10


def transfer(A, p, c):
    return p.conj() @ np.linalg.solve(A, c)


# ----------------------------------------------------------------------
# 1. Theorem I: exact structural zero via Krylov orthogonality
# ----------------------------------------------------------------------
def check_theorem_I():
    n = 6
    # Block-diagonal M: two invariant sectors that never mix.
    M1 = rng.standard_normal((3, 3))
    M2 = rng.standard_normal((3, 3))
    M = np.block([[M1, np.zeros((3, 3))], [np.zeros((3, 3)), M2]])
    r = np.concatenate([rng.standard_normal(3), np.zeros(3)])   # source in sector 1
    l = np.concatenate([np.zeros(3), rng.standard_normal(3)])   # readout in sector 2

    moments = []
    v = r.copy()
    for _ in range(n):
        moments.append(l.conj() @ v)
        v = M @ v
    max_moment = max(abs(m) for m in moments)

    zs = rng.standard_normal(5) + 1j * rng.standard_normal(5) + 10.0
    max_chi = max(abs(l.conj() @ np.linalg.solve(z * np.eye(n) - M, r)) for z in zs)

    ok = max_moment < TOL and max_chi < TOL
    print(f"[Theorem I ] max |mu_k| = {max_moment:.2e}, "
          f"max |delta_chi(z)| = {max_chi:.2e}  -> {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------
# 2. Theorem II: suppression index nu_diss = d + 1 (sector chain)
# ----------------------------------------------------------------------
def check_theorem_II():
    # 4-site chain: c lives on site 0, p on site 3, nearest-neighbour V only.
    # Graph distance d = 3 -> first nonzero moment m = 3, nu_diss = 4.
    n, d_expected = 4, 3
    D = np.diag(rng.uniform(1.0, 2.0, n))          # full rank
    V = np.zeros((n, n))
    for i in range(n - 1):
        V[i, i + 1] = V[i + 1, i] = rng.uniform(0.5, 1.5)
    z = 0.3 + 0.1j
    B = np.diag(rng.uniform(0.5, 1.0, n) * z) + V
    c = np.eye(n)[0]
    p = np.eye(n)[3]

    # Moments mu_k = p^dag (D^-1 B)^k D^-1 c
    Dinv = np.linalg.inv(D)
    X = Dinv @ B
    v = Dinv @ c
    moments = []
    for _ in range(n + 2):
        moments.append(p.conj() @ v)
        v = X @ v
    m = next(k for k, mu in enumerate(moments) if abs(mu) > TOL)
    assert all(abs(mu) < TOL for mu in moments[:m])

    # Fit slope of log|F_Gamma| vs log Gamma
    gammas = np.logspace(3, 5, 12)
    F = np.array([transfer(g * D + B, p, c) for g in gammas])
    slope = np.polyfit(np.log(gammas), np.log(np.abs(F)), 1)[0]

    nu = m + 1
    ok = (m == d_expected) and abs(slope + nu) < 0.01
    print(f"[Theorem II] first nonzero moment m = {m} (expected {d_expected}), "
          f"fitted slope = {slope:.4f} (expected {-nu})  -> {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------
# 3. Theorem III: protected O(1) channel under singular dissipation
# ----------------------------------------------------------------------
def riesz_projection_kernel(D, tol=1e-9):
    """Riesz projection onto the zero eigenvalue via eigendecomposition.

    Asserts semisimplicity (eigenvector matrix well conditioned and
    algebraic = geometric multiplicity of 0).
    """
    w, V = np.linalg.eig(D)
    zero = np.abs(w) < tol
    Vinv = np.linalg.inv(V)
    # semisimplicity check: D restricted to the zero eigenspace vanishes
    P = V[:, zero] @ Vinv[zero, :]
    assert np.linalg.norm(D @ P) < 1e-8, "zero eigenvalue not semisimple"
    return P


def check_theorem_III():
    n, k = 5, 2                        # 2-dimensional protected subspace
    # Non-Hermitian D with semisimple kernel: D = S diag(0,0,d3,d4,d5) S^-1
    S = rng.standard_normal((n, n)) + 0.1 * np.eye(n)
    dvals = np.concatenate([np.zeros(k), rng.uniform(1.0, 3.0, n - k)])
    D = S @ np.diag(dvals) @ np.linalg.inv(S)
    B = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    c = rng.standard_normal(n)
    p = rng.standard_normal(n)

    P = riesz_projection_kernel(D)

    # F0 = p^dag P (P B P)^{-1} P c evaluated on the protected subspace basis
    Pbasis = S[:, :k]                                   # basis of ker D
    Pdual = np.linalg.inv(S)[:k, :]                     # dual basis rows
    B_P = Pdual @ B @ Pbasis
    F0 = (p.conj() @ Pbasis) @ np.linalg.solve(B_P, Pdual @ c)

    gammas = np.logspace(3, 6, 10)
    F = np.array([transfer(g * D + B, p, c) for g in gammas])
    plateau_err = abs(F[-1] - F0) / abs(F0)
    # subleading correction should scale as 1/Gamma
    corr_slope = np.polyfit(np.log(gammas), np.log(np.abs(F - F0)), 1)[0]

    ok = plateau_err < 1e-5 and abs(corr_slope + 1) < 0.05
    print(f"[Theorem III] |F(Gamma_max) - F0|/|F0| = {plateau_err:.2e}, "
          f"correction slope = {corr_slope:.4f} (expected -1)  -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------
# 4. Counterexample: Pc != 0, p^dag P != 0, yet F0 = 0
# ----------------------------------------------------------------------
def check_counterexample():
    # Protected block: B_P = [[0,1],[1,0]], endpoints on component 1.
    n, k = 4, 2
    D = np.diag([0.0, 0.0, 1.0, 2.0])
    B = np.zeros((n, n), dtype=complex)
    B[0, 1] = B[1, 0] = 1.0
    B[2:, 2:] = rng.standard_normal((2, 2)) + 3 * np.eye(2)
    B[0, 2] = B[2, 0] = 0.4          # leakage P<->Q allowed
    c = np.array([1.0, 0.0, 0.3, 0.0])
    p = np.array([1.0, 0.0, 0.0, 0.2])

    P = np.diag([1.0, 1.0, 0.0, 0.0])
    Pc_norm = np.linalg.norm(P @ c)
    pP_norm = np.linalg.norm(P.T @ p)
    B_P = B[:k, :k]
    F0 = p[:k].conj() @ np.linalg.solve(B_P, c[:k])

    gammas = np.logspace(3, 6, 8)
    F = np.array([transfer(g * D + B, p, c) for g in gammas])
    decays = np.abs(F[-1]) < 1e-4    # response dies despite nonzero projections

    ok = Pc_norm > 0.5 and pP_norm > 0.5 and abs(F0) < TOL and decays
    print(f"[Counterexample] |Pc| = {Pc_norm:.2f}, |p^dag P| = {pP_norm:.2f}, "
          f"F0 = {abs(F0):.2e}, |F(Gamma_max)| = {abs(F[-1]):.2e}  -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    results = [check_theorem_I(), check_theorem_II(),
               check_theorem_III(), check_counterexample()]
    print("\nAll checks passed." if all(results) else "\nSOME CHECKS FAILED.")
    raise SystemExit(0 if all(results) else 1)
