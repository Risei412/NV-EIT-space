"""Numerical infrastructure for the operational sector cut (Theorem 0A-0E,
Sector/section0_operational_foundations.tex) and its algebraic-cut special
case (Sector/theorem2B_operational_realization.tex).

This module is intentionally model-agnostic and does not modify
`New no-go theory/src/core.py`; it reuses the Riesz-projection machinery
already implemented there instead of duplicating it, per Definition 0.2
(admissible cut generator) and Theorem 0A.

    D_S           : the (Gamma-independent) selective cut generator, with
                     nonzero eigenvalues taken with positive real part
                     (damping convention of Sec. 4.1 of the revision note).
    A(z)          : the fixed-Gamma, fixed-z generator A_Gamma(z) (may be a
                     plain matrix for the algebraic/quasi-static case used
                     in the Lambda-model block realization).
    kappa         : cut-intervention strength; the ideal cut is
                     kappa -> infinity (Theorem 0A).
"""

import sys
import os
import importlib.util

import numpy as np

_CORE_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "src", "core.py"
)
_spec = importlib.util.spec_from_file_location("nogo_core", _CORE_PATH)
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)


def validate_cut_generator(D_S, tol=1e-9):
    """Check admissibility conditions (C2)-(C3) of Definition 0.2 that are
    checkable from D_S alone (C1, C4, C5 depend on the full Lindbladian and
    reference state and are checked separately, e.g. by
    `check_noninvasive`).

    Returns a dict with:
      - 'semisimple'      : bool, zero eigenvalue of D_S is semisimple (C3)
      - 'P_S', 'Q_S'       : the Riesz projectors (oblique in general)
      - 'is_hermitian'     : bool, D_S is Hilbert-Schmidt self-adjoint on
                              this matrix representation (triggers Lemma 0B')
      - 'nonzero_real_parts': real parts of the nonzero spectrum (must all
                              be > 0 under the damping convention, i.e.
                              selective damping, condition (C2))
    """
    w = np.linalg.eigvals(D_S)
    nz = w[np.abs(w) > tol]
    is_hermitian = np.allclose(D_S, D_S.conj().T, atol=1e-10)
    P_S = core.riesz_projection_exact(D_S, tol=tol)
    Q_S = np.eye(D_S.shape[0]) - P_S
    semisimple = np.linalg.norm(D_S @ P_S) < 1e-7 and np.linalg.norm(P_S @ D_S) < 1e-7
    selective_damping_ok = bool(np.all(nz.real > -tol)) if len(nz) else True
    return {
        "semisimple": semisimple,
        "P_S": P_S,
        "Q_S": Q_S,
        "is_hermitian": is_hermitian,
        "nonzero_real_parts": nz.real,
        "selective_damping_ok": selective_damping_ok,
    }


def check_noninvasive(C_S_apply, rho0, tol=1e-10):
    """Condition (C4): C_S(rho0) == 0. `C_S_apply` is a callable
    rho -> C_S(rho) (superoperator action), `rho0` a density matrix."""
    out = C_S_apply(rho0)
    return float(np.linalg.norm(out)) < tol, out


def operational_cut_response(A, D_S, c, p, kappa):
    """chi_{Gamma,kappa} = p^dag (A + kappa D_S)^-1 c, a single finite-kappa
    evaluation (pre-limit quantity of Theorem 0A)."""
    return core.transfer(A + kappa * D_S, p, c)


def ideal_cut_response(A, D_S, c, p, tol=1e-9):
    """chi_cut^(S) = p^dag P_S [P_S A P_S]^-1 P_S c, formula (0.2) of
    Theorem 0A, computed directly from the Riesz projector (not as a
    kappa -> infinity numerical limit)."""
    P_S = core.riesz_projection_exact(D_S, tol=tol)
    mask = np.abs(np.diag(P_S)) > 0.5  # valid when P_S is a coordinate
    # projector (block-diagonal cut generator in the working basis, as in
    # the Theorem 2B / Lambda-model realization); for a genuinely oblique
    # P_S use the generic projector contraction below instead.
    if np.allclose(P_S, np.diag(np.diag(P_S)), atol=1e-10):
        A_P = A[np.ix_(mask, mask)]
        Pc = c[mask]
        pP = np.conj(p[mask])
        return pP @ np.linalg.solve(A_P, Pc)
    # generic oblique projector: solve on the range of P_S via its rank
    Pc = P_S @ c
    A_P = P_S @ A @ P_S
    # invert on Ran(P_S): use pseudo-inverse restricted to the range basis
    r = np.linalg.matrix_rank(P_S, tol=1e-8)
    U, s, Vh = np.linalg.svd(P_S)
    basis = U[:, :r]
    A_small = basis.conj().T @ A_P @ basis
    c_small = basis.conj().T @ Pc
    x_small = np.linalg.solve(A_small, c_small)
    x = basis @ x_small
    return np.conj(p) @ x


def compare_cut_equivalence(A, D_S, c, p, kappas):
    """Gate U1/U2 helper: chi_op(kappa) for a range of kappa, the ideal
    limit chi_cut^(S), and the fitted O(kappa^-1) slope/coefficient."""
    chi_ideal = ideal_cut_response(A, D_S, c, p)
    chi_op = np.array([operational_cut_response(A, D_S, c, p, k) for k in kappas])
    diff = chi_op - chi_ideal
    kappas = np.asarray(kappas, dtype=float)
    # fit log|diff| = log|coeff| - 1*log(kappa) over the tail (large kappa)
    logk = np.log(kappas)
    logd = np.log(np.abs(diff))
    slope, intercept = np.polyfit(logk, logd, 1)
    fitted_coeff_mag = np.exp(intercept)
    return {
        "chi_ideal": chi_ideal,
        "chi_op": chi_op,
        "diff": diff,
        "fit_slope": float(slope),        # should be approx -1
        "fitted_coeff_mag": float(fitted_coeff_mag),
    }
