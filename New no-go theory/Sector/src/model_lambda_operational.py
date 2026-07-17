"""Minimal Lambda-model realization of the operational sector cut, matching
Corollary "Algebraic-operational equivalence for the block cut"
(Sector/theorem2B_operational_realization.tex) and Sec. 11 of the revision
note. Reuses the exact block structure of `New no-go theory/src/model_lambda.py`
(A = gamma31 - i*Delta, coupling Omega_c, sector G_S = gamma21 - i*delta),
adding the explicit cut generator D_S = diag(0, 1) so that

    M_kappa = M + kappa * D_S = [[A, B], [C, G_S + kappa]]

is exactly A_{Gamma,kappa} of Theorem 0A specialized to this block.
"""

import numpy as np


def M_matrix(Delta, delta, Omega_c, gamma31, gamma21):
    return np.array([
        [gamma31 - 1j * Delta, 1j * Omega_c / 2],
        [1j * np.conj(Omega_c) / 2, gamma21 - 1j * delta],
    ], dtype=complex)


def D_S():
    """Selective cut generator on the ground-coherence sector (index 1)."""
    return np.array([[0.0, 0.0], [0.0, 1.0]], dtype=complex)


def source_readout():
    c = np.array([1.0, 0.0], dtype=complex)
    p = np.array([1.0, 0.0], dtype=complex)
    return c, p


def closed_form_kappa_coefficient(Delta, Omega_c, gamma31):
    """chi_0 d^dag A^-1 B C A^-1 b_p, the explicit O(kappa^-1) coefficient
    derived in the Corollary by direct block elimination:
        coeff = B*C/A^2 = -|Omega_c|^2 / (4 A^2),  A = gamma31 - i*Delta.
    """
    A = gamma31 - 1j * Delta
    return -(abs(Omega_c) ** 2 / 4.0) / (A ** 2)
