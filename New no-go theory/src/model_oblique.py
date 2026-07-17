"""Phase D, run D1: oblique (Riesz) protection vs the old orthogonal-
projector formula, and dissipation-enhanced (non-monotonic) sector response.

Model: D is block upper-triangular,
    D(a) = [[0, 0, a, 0],
            [0, 0, 0, a],
            [0, 0, 1, 0],
            [0, 0, 0, 2]]
The right null space of D is span{e1,e2} for every a (rows 3,4 of Dv=0 force
v3=v4=0 whenever a != 0, independent of a), so the *orthogonal* projector onto
ker D is P_orth = diag(1,1,0,0) for every a. But the correct Riesz (spectral)
projector -- the one Theorem 6.2 actually requires for a non-Hermitian D --
picks up an a-dependent oblique component from the eigenvectors of the
nonzero eigenvalues (1 and 2), and differs from P_orth for a != 0.

Sector cut: the protected-block internal coupling kappa (between indices
0,1 -- the "EIT ground-coherence" pair) is present in B_full and severed in
B_cut (frozen-source cut, exactly as in model_protected.transition_matrices).
Fixed leakage terms couple the protected block to the damped block
(B[0,2]=B[2,0]=leak, B[1,3]=B[3,1]=leak), identical in full and cut.
"""

import numpy as np

from core import riesz_projection_exact, proj_orth_kernel, condition_number_of_projector


def D_of_a(a):
    return np.array([
        [0.0, 0.0, a,   0.0],
        [0.0, 0.0, 0.0, a],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 2.0],
    ])


def make_B(kappa, e1=1.3, e2=0.9, leak=0.55, dQ1=1.1, dQ2=0.8):
    B = np.zeros((4, 4))
    B[0, 0], B[1, 1] = e1, e2
    B[0, 1] = B[1, 0] = kappa
    B[2, 2], B[3, 3] = dQ1, dQ2
    B[0, 2] = B[2, 0] = leak
    B[1, 3] = B[3, 1] = leak
    return B


C_VEC = np.array([1.0, 0.35, 0.22, -0.13])
P_VEC = np.array([1.0, -0.22, 0.16, 0.27])
KAPPA_FULL = 0.6  # protected-block internal coupling, present in full model


def full_cut_B(kappa_full=KAPPA_FULL):
    return make_B(kappa_full), make_B(0.0)


def _restricted_operator(B, W_basis, W_coord):
    """B_P matrix (r x r) representing P B P restricted to Im(P) = span(W_basis
    columns), where P = W_basis @ W_coord and W_coord @ W_basis = I_r."""
    return W_coord @ B @ W_basis


def r0_via_projector(D, B_full, B_cut, c, p, projector="riesz", tol=1e-9):
    """delta_chi_{S,0} = p^dag P [B_P,full^-1 - B_P,cut^-1] P c using either
    the correct oblique Riesz projector or the (generically wrong for
    non-Hermitian D) orthogonal projector."""
    if projector == "riesz":
        w, V = np.linalg.eig(D)
        mask = np.abs(w) < tol
        Vinv = np.linalg.inv(V)
        W_basis = V[:, mask]
        W_coord = Vinv[mask, :]
        P = W_basis @ W_coord
        assert np.linalg.norm(D @ P) < 1e-7
    elif projector == "orth":
        U, s, _ = np.linalg.svd(D)
        mask = s < tol
        W_basis = U[:, mask]
        W_coord = W_basis.conj().T
        P = W_basis @ W_coord
    else:
        raise ValueError(projector)

    p_coords = np.conj(p) @ W_basis
    c_full_coords = W_coord @ c
    c_cut_coords = c_full_coords  # c is identical between full/cut (frozen source)

    BP_full = _restricted_operator(B_full, W_basis, W_coord)
    BP_cut = _restricted_operator(B_cut, W_basis, W_coord)

    F0_full = p_coords @ np.linalg.solve(BP_full, c_full_coords)
    F0_cut = p_coords @ np.linalg.solve(BP_cut, c_cut_coords)
    return F0_full - F0_cut, P


def R_S_gamma(Gamma, a, kappa_full=KAPPA_FULL):
    D = D_of_a(a)
    B_full, B_cut = full_cut_B(kappa_full)
    A_full = Gamma * D + B_full
    A_cut = Gamma * D + B_cut
    chi_full = np.conj(P_VEC) @ np.linalg.solve(A_full, C_VEC)
    chi_cut = np.conj(P_VEC) @ np.linalg.solve(A_cut, C_VEC)
    return chi_full, chi_cut, chi_full - chi_cut


def scan_a(a_values, kappa_full=KAPPA_FULL):
    out = []
    for a in a_values:
        D = D_of_a(a)
        B_full, B_cut = full_cut_B(kappa_full)
        r0_riesz, P_riesz = r0_via_projector(D, B_full, B_cut, C_VEC, P_VEC, "riesz")
        r0_orth, P_orth = r0_via_projector(D, B_full, B_cut, C_VEC, P_VEC, "orth")
        kappa_P = condition_number_of_projector(P_riesz)
        out.append({
            "a": a,
            "r0_riesz": r0_riesz,
            "r0_orth": r0_orth,
            "kappa_P": kappa_P,
        })
    return out
