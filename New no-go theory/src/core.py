"""Common numerical infrastructure for the sector-resolved response
classification (Theorems I-III, three_theorems_proofs.tex).

Implements, independently of any specific model:
  - the scalar transfer function F_Gamma(z) = p^dag A_Gamma(z)^-1 c;
  - three independent estimators of the suppression index nu:
      (1) direct log-log fit of ||R||_K vs Gamma (local effective index),
      (2) moment method (Theorem II, mu_k = p^dag (D^-1 B)^k D^-1 c),
      (3) protected-coefficient method (Theorem III, Riesz projection,
          F0 = p^dag P B_P^-1 P c, F1 from the Schur-complement expansion);
  - the Class-I Krylov certificate (Theorem I): mu_k = l^dag M^k r = 0 for
    k < n, exact at machine precision -- never inferred from a small float.

All quantities are kept explicit (full response, cut response, and the
difference R_S) as required by the priorities document: chi_full alone is
never computed without chi_cut and R_S alongside it.
"""

import numpy as np

TOL_CLASS_I = 1e-10


def transfer(A, p, c):
    """F = p^dag A^-1 c."""
    return np.conj(p) @ np.linalg.solve(A, c)


def krylov_class1_certificate(M, r, l):
    """Exact (machine-precision) Class-I test: mu_k = l^dag M^k r, k<n.

    Returns (max|mu_k|, moments). Class I iff max|mu_k| < TOL_CLASS_I; this
    is a finite certificate (Theorem I), not an asymptotic inference.
    """
    n = M.shape[0]
    moments = []
    v = r.copy()
    for _ in range(n):
        moments.append(np.conj(l) @ v)
        v = M @ v
    moments = np.array(moments)
    return float(np.max(np.abs(moments))), moments


def moments_invertible_D(D, B_of_z, c, p, z, kmax):
    """mu_k(z) = p^dag (D^-1 B(z))^k D^-1 c, k = 0..kmax-1 (Theorem II)."""
    Dinv = np.linalg.inv(D)
    X = Dinv @ B_of_z(z)
    v = Dinv @ c
    out = []
    for _ in range(kmax):
        out.append(np.conj(p) @ v)
        v = X @ v
    return np.array(out)


def nu_from_moments(D, B_of_z, c, p, z, kmax, tol=TOL_CLASS_I):
    """First nonzero moment index m -> nu_diss = m + 1 (Theorem II), or
    None if all kmax moments vanish (candidate Class I at this z)."""
    mu = moments_invertible_D(D, B_of_z, c, p, z, kmax)
    nz = np.nonzero(np.abs(mu) > tol)[0]
    if len(nz) == 0:
        return None, mu
    m = int(nz[0])
    return m + 1, mu


def riesz_projection_diag(D, tol=1e-9):
    """Riesz projection onto ker D for a D that is exactly block-diagonal
    (or diagonal) in the working basis, i.e. D = blockdiag(0_k, D_Q) with
    D_Q invertible. Returns boolean mask of the kernel (protected) indices.

    Verifies semisimplicity (D @ P == 0) exactly, matching Lemma "Riesz
    splitting of a semisimple kernel".
    """
    diag_like = np.allclose(D, np.diag(np.diag(D)), atol=1e-12)
    if diag_like:
        mask = np.abs(np.diag(D)) < tol
    else:
        w, V = np.linalg.eig(D)
        mask = np.abs(w) < tol
        Vinv = np.linalg.inv(V)
        P = (V[:, mask] @ Vinv[mask, :]).real if np.allclose(V.imag, 0) else V[:, mask] @ Vinv[mask, :]
        assert np.linalg.norm(D @ P) < 1e-7, "zero eigenvalue of D is not semisimple"
    return mask


def protected_F0(B_z, mask, c, p):
    """F0(z) = p^dag P B_P(z)^-1 P c, with P the coordinate projection onto
    `mask` (valid when D is already block-diagonal in the working basis,
    i.e. P is the exact coordinate projector -- Theorem III leading order).
    """
    B_P = B_z[np.ix_(mask, mask)]
    Pc = c[mask]
    pP = np.conj(p[mask])
    return pP @ np.linalg.solve(B_P, Pc)


def riesz_projection_exact(D, tol=1e-9):
    """Riesz spectral projector onto the (assumed semisimple) zero eigenspace
    of a general (possibly non-normal) matrix D, via eigendecomposition:

        P = V[:, mask] @ Vinv[mask, :],   D @ P == 0,  P @ D == 0.

    Unlike `riesz_projection_diag`, this does not assume D is diagonal or
    Hermitian; it is the genuine oblique (non-orthogonal) Riesz projector of
    Assumption 6.1 in the master-response package, used to contrast against
    the orthogonal projector `proj_orth_kernel` (Phase D, run D1).
    """
    w, V = np.linalg.eig(D)
    mask = np.abs(w) < tol
    Vinv = np.linalg.inv(V)
    P = V[:, mask] @ Vinv[mask, :]
    assert np.linalg.norm(D @ P) < 1e-7, "zero eigenvalue of D is not semisimple"
    assert np.linalg.norm(P @ D) < 1e-7, "zero eigenvalue of D is not semisimple"
    return P


def proj_orth_kernel(D, tol=1e-9):
    """Orthogonal projector onto ker D (valid formula only for D = D^dag >= 0;
    used in Phase D run D1 as the 'old-theory' projector to contrast against
    the correct oblique Riesz projector for non-Hermitian D)."""
    U, s, _ = np.linalg.svd(D)
    mask = s < tol
    if not np.any(mask):
        # ker D trivial in the numerical sense; still return the projector
        # onto the smallest-singular-value direction for comparison plots
        mask = s == s.min()
    Uk = U[:, mask]
    return Uk @ Uk.conj().T


def condition_number_of_projector(P):
    """kappa(P): a nonzero, non-orthogonal projector has singular values
    {1 (multiplicity r), 0 (multiplicity n-r)} only if orthogonal; for an
    oblique projector the nonzero singular values differ from 1, and
    kappa(P) = sigma_max(P)/sigma_min_nonzero(P) measures the obliqueness
    (kappa=1 iff P is an orthogonal projector)."""
    s = np.linalg.svd(P, compute_uv=False)
    s_nonzero = s[s > 1e-12]
    if len(s_nonzero) == 0:
        return 1.0
    return float(s_nonzero.max() / s_nonzero.min())


def certificate_deg_nu(D_sym, B_sym_of_z, c_sym, p_sym, Gamma, z_val, z_sym=None):
    """Exact finite polynomial certificate (Theorem 8.1, Step 5): construct
    N(Gamma,z) = p^dag adj(Gamma D + B(z)) c and Q(Gamma,z) = det(Gamma D + B(z))
    as exact SymPy polynomials in Gamma (at a fixed numeric or symbolic z),
    and return nu = deg_Gamma Q - deg_Gamma N together with the polynomials.

    All arithmetic is exact (sympy Rational/complex-rational), never floating
    point, per Assumption 1.3 / Remark 1.4 of the master-response package.
    """
    import sympy as sp

    Bz = B_sym_of_z(z_val)
    A = Gamma * D_sym + Bz
    Q = sp.together(A.det())
    Qp = sp.Poly(sp.expand(Q), Gamma)
    adjA = A.adjugate()
    N = sp.together((p_sym.conjugate().T * adjA * c_sym)[0, 0])
    Nexp = sp.expand(N)
    if Nexp == 0:
        return None, Qp, None  # Class I candidate: numerator identically zero
    Np = sp.Poly(Nexp, Gamma)
    nu = Qp.degree() - Np.degree()
    return int(nu), Qp, Np


def fit_nu_loglog(gammas, values):
    """Global log-log slope (least squares) and local effective index
    nu_eff(Gamma) from adjacent points (finite-difference slope), per
    Priority 2 of the priorities document.

    nu = -d log|R|/d log Gamma.
    """
    gammas = np.asarray(gammas, dtype=float)
    vals = np.asarray(np.abs(values), dtype=float)
    logG = np.log(gammas)
    logV = np.log(np.maximum(vals, 1e-300))
    slope_global = -np.polyfit(logG, logV, 1)[0]
    nu_eff = -np.diff(logV) / np.diff(logG)
    gamma_mid = np.sqrt(gammas[:-1] * gammas[1:])
    return {
        "nu_global": float(slope_global),
        "nu_eff": nu_eff,
        "gamma_mid": gamma_mid,
    }


def max_abs_on_grid(func_of_z, zgrid):
    """||f||_K = max_{z in K} |f(z)| over a sampled compact window K."""
    vals = np.array([func_of_z(z) for z in zgrid])
    return float(np.max(np.abs(vals))), vals
