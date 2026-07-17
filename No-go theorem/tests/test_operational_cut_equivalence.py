"""Regression test: gate2_candidate_full_vs_reduced.py's ground-coherence
cut (an ad hoc "zero the S<->X blocks of L" construction) is exactly the
algebraic-cut special case D_S = diag(0, I_S) of the operational sector cut
(Theorem 0A / Corollary "Algebraic-operational equivalence for the block
cut", New no-go theory/Sector/theorem2B_operational_realization.tex),
implemented independently in New no-go theory/Sector/src/operational_cut.py.

This closes the room-temperature no-go plan's preflight item: the existing
NV gate2 cut and the new operational-cut machinery must be shown to agree,
not merely asserted to agree, before Step 2 (operational cut audit) reuses
either one on the full NV Liouvillian.

Checks:
  1. chi_cut from gate2's own block-zeroing construction equals
     operational_cut.ideal_cut_response(A_full, D_S, c, p) with D_S the
     coordinate cut generator for the same two-index subspace S.
  2. The finite-kappa intervention chi_op(kappa) = p^dag (A+kappa D_S)^-1 c
     converges to that same ideal value as O(kappa^-1) (Gate 2 of the plan).
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
SECTOR_SRC = os.path.join(HERE, "..", "..", "New no-go theory", "Sector", "src")
sys.path.insert(0, SRC)
sys.path.insert(0, SECTOR_SRC)

import gate2_candidate_full_vs_reduced as g2
import operational_cut as oc


def _candidate_full_generator_and_cut_subspace(d2=0.0):
    """Rebuild A(z) = -L (steady-state generator convention: L rho = 0 at
    z=0 defines the reference; the weak-probe first-order system solves
    (L - i*z) drho = source, i.e. A(z) = i*z*I - L in gate2's `first_order`
    convention) and the ground-coherence index pair S, reusing gate2's own
    `build_full` so there is exactly one source of truth for the candidate
    Liouvillian."""
    H, Ls, Vp, dp, meta = g2.build_full(70.0, g2.rp.BX0, g2.rp.BZ0, d2, isc=False)
    N, p_idx, c_idx = meta['N'], meta['p_idx'], meta['c_idx']
    from liouvillian_core import liouvillian
    L = liouvillian(H, Ls)
    S = [c_idx * N + p_idx, p_idx * N + c_idx]
    return L, N, S, meta


def test_gate2_cut_equals_operational_ideal_cut():
    L, N, S, meta = _candidate_full_generator_and_cut_subspace()
    dim = N * N

    # gate2's own construction: zero the S<->complement blocks of L.
    X = [k for k in range(dim) if k not in S]
    Lc_gate2 = L.copy()
    Lc_gate2[np.ix_(S, X)] = 0
    Lc_gate2[np.ix_(X, S)] = 0

    # operational_cut's algebraic D_S = diag(0, I_S) construction: retained
    # subspace is the complement of S (P_EIT keeps X, i.e. suppresses S).
    D_S = np.zeros((dim, dim), complex)
    for s in S:
        D_S[s, s] = 1.0

    # Use a generic complex source/readout pair sourced/read OFF the cut
    # subspace S (c_S = p_S = 0) -- matching the physical setup (probe
    # polarization / injected coherence never lives on the ground-coherence
    # indices) and condition (C4) of Definition 0.2 (the cut generator's
    # sector is non-invasive w.r.t. the declared source/readout). Both use
    # the SAME A, c, p so the comparison isolates the cut construction.
    rng = np.random.default_rng(0)
    c = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    p = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    c[S] = 0.0
    p[S] = 0.0
    z = 0.37  # generic off-resonance probe detuning, arbitrary units of L

    A = 1j * z * np.eye(dim) - L
    A_cut_gate2 = 1j * z * np.eye(dim) - Lc_gate2

    chi_full = oc.core.transfer(A, p, c)
    chi_cut_gate2 = np.conj(p) @ np.linalg.solve(A_cut_gate2, c)
    chi_cut_operational = oc.ideal_cut_response(A, D_S, c, p)

    assert abs(chi_cut_gate2 - chi_cut_operational) / abs(chi_cut_operational) < 1e-9, (
        "gate2's ad hoc block-zeroed cut and the operational_cut.py "
        "algebraic-cut D_S construction disagree on the same candidate "
        "Liouvillian -- they must be numerically identical (Theorem 2B)."
    )
    return A, D_S, c, p, chi_cut_operational


def test_finite_kappa_converges_O_kappa_inverse():
    A, D_S, c, p, chi_ideal = test_gate2_cut_equals_operational_ideal_cut()
    kappas = np.array([1e2, 1e3, 1e4, 1e5, 1e6])
    result = oc.compare_cut_equivalence(A, D_S, c, p, kappas)
    assert abs(result["fit_slope"] + 1.0) < 0.05, (
        f"finite-kappa intervention should converge as O(kappa^-1); "
        f"fitted slope={result['fit_slope']:.4f}"
    )
    assert abs(result["chi_op"][-1] - chi_ideal) / abs(chi_ideal) < 1e-2


if __name__ == "__main__":
    test_gate2_cut_equals_operational_ideal_cut()
    test_finite_kappa_converges_O_kappa_inverse()
    print("PASS: gate2 cut == operational_cut.py ideal cut; O(kappa^-1) convergence confirmed.")
