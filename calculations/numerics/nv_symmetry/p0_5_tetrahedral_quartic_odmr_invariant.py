"""P0-5: basis-invariant tetrahedral quartic ODMR invariant.

Clean replacement for using a one-sided ms=+/- channel at B=0.
For each of the four NV orientations, diagonalize the ground-state spin-1
Hamiltonian

    H_a = D Sz^2 + gamma B_a . S,

and identify only the nondegenerate eigenvalue continuously connected to
ms=0. The sum of the two ODMR transition frequencies for orientation a is
basis invariant inside the degenerate ms=+/- doublet:

    Omega_sum,a = Tr(H_a) - 3 lambda_0,a = 2D - 3 lambda_0,a.

Average this over all four tetrahedral orientations. The mean of all eight
ODMR transition frequencies is one half of this quantity and therefore does
not require orientation labeling at the level of the mathematical observable.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

D = 2.877
GAMMA = 28.02495164
S2 = 1 / np.sqrt(2)
SZ = np.diag([-1.0, 0.0, 1.0]).astype(complex)
SX = S2 * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], complex)
SY = S2 * np.array([[0, 1j, 0], [-1j, 0, 1j], [0, -1j, 0]], complex)
AXES = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) / np.sqrt(3)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def omega_sum_one(Bglob, n):
    Bglob = np.asarray(Bglob, float)
    n = np.asarray(n, float)
    u = float(Bglob @ n)
    bp2 = max(float(Bglob @ Bglob - u * u), 0.0)
    H = D * (SZ @ SZ) + GAMMA * (np.sqrt(bp2) * SX + u * SZ)
    w = np.linalg.eigvalsh(H)
    lam0 = float(w[int(np.argmin(np.abs(w)))])
    return float(np.sum(w) - 3 * lam0)


def omega_bar(Bglob):
    return float(np.mean([omega_sum_one(Bglob, n) for n in AXES]))


def centroid8(Bglob):
    return 0.5 * omega_bar(Bglob)


def analytic_omega_bar(Bglob):
    B = np.asarray(Bglob, float)
    b2 = float(B @ B)
    s4 = float(np.sum(B**4))
    return float(
        2 * D
        + 2 * GAMMA**2 / D * b2
        - 2 * GAMMA**4 / D**3 * b2**2
        + 4 * GAMMA**4 / (3 * D**3) * s4
    )


def direction(name, r):
    if name == "100":
        return np.array([r, 0.0, 0.0])
    if name == "110":
        return np.array([r, r, 0.0]) / np.sqrt(2)
    if name == "111":
        return np.array([r, r, r]) / np.sqrt(3)
    raise ValueError(name)


def main():
    coef_s4 = 4 * GAMMA**4 / (3 * D**3)
    coef_100_minus_110 = 2 * GAMMA**4 / (3 * D**3)
    rows = []
    for r in (1e-4, 5e-4, 1e-3, 2e-3, 5e-3):
        b100 = direction("100", r)
        b110 = direction("110", r)
        b111 = direction("111", r)
        o100 = omega_bar(b100)
        o110 = omega_bar(b110)
        o111 = omega_bar(b111)
        rows.append(
            dict(
                B_T=r,
                omega100_GHz=o100,
                omega110_GHz=o110,
                omega111_GHz=o111,
                delta100_110_GHz=o100 - o110,
                delta100_110_over_B4=(o100 - o110) / r**4,
                analytic_delta100_110_over_B4=coef_100_minus_110,
                relative_formula_error_100=abs(o100 - analytic_omega_bar(b100)) / max(abs(o100 - 2 * D), 1e-300),
                relative_formula_error_110=abs(o110 - analytic_omega_bar(b110)) / max(abs(o110 - 2 * D), 1e-300),
            )
        )

    M1 = AXES.sum(axis=0)
    M2 = np.einsum("ai,aj->ij", AXES, AXES)
    M3 = np.einsum("ai,aj,ak->ijk", AXES, AXES, AXES)

    out = dict(
        model="NV ground-state spin-1 Hamiltonian; four tetrahedral orientations; both ODMR transitions summed",
        D_GHz=D,
        gamma_GHz_per_T=GAMMA,
        exact_characteristic_polynomial="lambda^3 - 2D lambda^2 + (D^2-p^2-q^2)lambda + D p^2 = 0",
        local_small_field_m0_branch="lambda0 = -p^2/D + (p^4-p^2 q^2)/D^3 + O(B^6)",
        tetrahedral_moments=dict(M1=M1.tolist(), M2=M2.tolist(), M3_xyz=float(M3[0, 1, 2])),
        averaged_transition_sum_formula=(
            "Omega_bar = 2D + (2 gamma^2/D) B^2 - (2 gamma^4/D^3)(B^2)^2 "
            "+ (4 gamma^4/(3D^3))(Bx^4+By^4+Bz^4) + O(B^6)"
        ),
        all_eight_line_centroid_formula=(
            "f_centroid = D + (gamma^2/D) B^2 - (gamma^4/D^3)(B^2)^2 "
            "+ (2 gamma^4/(3D^3))(Bx^4+By^4+Bz^4) + O(B^6)"
        ),
        tensor_selection=dict(
            rank1="forbidden",
            rank2="isotropic only",
            rank3="forbidden in this field-even spectral sum",
            rank4="first direction-dependent tensor",
        ),
        rank4_coefficients=dict(
            S4_in_Omega_bar=coef_s4,
            difference_100_minus_110_per_B4=coef_100_minus_110,
            S4_in_eight_line_centroid=0.5 * coef_s4,
        ),
        numerical_validation=rows,
        matched_pair=dict(
            fields="equal-magnitude [100] and [110]",
            lower_aggregate_invariants="same B^2; both have BxByBz=0; orientation labels not used by centroid",
            distinguishing_invariant="S4=Bx^4+By^4+Bz^4: r^4 versus r^4/2",
            leading_centroid_difference="(gamma^4/(3D^3)) r^4 + O(r^6)",
        ),
        status=dict(
            mathematical_witness="verified analytically and by direct diagonalization",
            PRL_novelty="HOLD",
            reason="quartic invariant is a natural consequence of standard spin-1 Hamiltonian plus tetrahedral averaging; operational unresolved-spectrum capability still must be established",
        ),
    )
    path = REPO / "results" / "tables" / "p0_5_tetrahedral_quartic_odmr_invariant.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
