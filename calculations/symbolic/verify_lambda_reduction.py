"""Check the response kernel against the EIT literature, symbolically.

The manuscript builds its no-go criterion on a Schur complement of the
first-order response,

    ( A    B  ) ( x )   ( b_p )
    ( C   G_g ) ( s ) = (  0  ),

with x the optical-coherence block and s the targeted ground coherence, and
defines the pathway cut by setting B = C = 0.  That construction is only worth
anything if, in the case where the answer is already known, it returns the
known answer.  This script checks that it does, in the ideal three-level
Lambda system, against the standard weak-probe result of Fleischhauer,
Imamoglu and Marangos, Rev. Mod. Phys. 77, 633 (2005).

Four things are checked, and all four are exact identities rather than
numerical agreement:

  1. eliminating x through the Schur complement gives the same rho_e1 as
     solving the 2x2 system directly -- i.e. the "never divide by G_g" route
     costs nothing;
  2. that rho_e1 is the textbook EIT coherence, with the control field
     entering as the continued-fraction term |Omega_c|^2/4 over the two-photon
     denominator;
  3. at perfect two-photon resonance with no ground dephasing, chi_full = 0
     while the sector-resolved correction is maximal, dXi_S = -1;
  4. for N_e resolved excited branches the same elimination returns the
     Mishina-type multilevel susceptibility of the group's theory note
     (EIT_no_go_go_theory_v6_2_English.tex, Sec. 9.7), in which each of
     S_1, S_2, K_12 and K_21 is a coherent sum over branches.

Point 3 is the paper's central conceptual claim -- that chi_full = 0 is what
perfect EIT looks like, so a criterion phrased on the total susceptibility
would report ideal EIT as an absence of EIT -- and here it is a two-line
consequence rather than an assertion.

Point 4 is the case NV actually is: the two ^3E branches make K_12 a coherent
sum, not a single amplitude, so the three-level reduction alone would not
license the criterion used in the manuscript.

Usage: python verify_lambda_reduction.py     (exit 0 iff all four hold)
"""
from __future__ import annotations

import sys

import sympy as sp

I = sp.I


def build():
    """The ideal Lambda system's first-order block system.

    |1> is probe-coupled, |2> control-coupled, |e> the single excited state;
    the probe is weak, so all population sits in |1>.  These are the standard
    optical Bloch equations for that configuration, written in the block form
    of the manuscript's Eq. (2).
    """
    g31, g21 = sp.symbols("gamma_31 gamma_21")
    dp, dc = sp.symbols("delta_p delta_c")
    Oc, Op = sp.symbols("Omega_c Omega_p")

    A = g31 + I * dp                      # optical coherence rho_e1
    B = -I * Oc / 2                       # control couples s into the x row
    C = -I * sp.conjugate(Oc) / 2         # and x into the s row
    Gg = g21 + I * (dp - dc)              # two-photon denominator
    bp = I * Op / 2                       # weak-probe source

    return dict(A=A, B=B, C=C, Gg=Gg, bp=bp,
                g31=g31, g21=g21, dp=dp, dc=dc, Oc=Oc, Op=Op)


def multilevel(Ne):
    """Sec. 9.7's multi-excited-state reduction, for N_e resolved branches.

    The optical block is diagonal in the branch basis, A = diag(a_j), and the
    probe and control reach the branches through dipole vectors d_1 and d_2.
    Nothing here is specialised to NV: the a_j, d_1j and d_2j are free complex
    symbols, so the identity holds for any diagonal optical generator.

    The control enters the optical rows as B and the ground-coherence row as
    C = -B^dagger.  That relative sign is not a convention -- it is what the
    optical Bloch equations give, and it is why S_g comes out as G_g + beta S_2
    with a plus rather than a minus.
    """
    a = sp.symbols(f"a1:{Ne + 1}")
    d1 = sp.symbols(f"d1_1:{Ne + 1}")
    d2 = sp.symbols(f"d2_1:{Ne + 1}")
    beta = sp.Symbol("beta", positive=True)
    Gg = sp.Symbol("G_g")

    A = sp.diag(*a)
    D1 = sp.Matrix(Ne, 1, list(d1))
    D2 = sp.Matrix(Ne, 1, list(d2))
    B = -I * sp.sqrt(beta) * D2
    C = -B.T.applyfunc(sp.conjugate)

    Sg = Gg - (C * A.inv() * B)[0, 0]
    x = A.inv() * D1 + A.inv() * B * ((C * A.inv() * D1)[0, 0]) / Sg
    Xi = (D1.T.applyfunc(sp.conjugate) * x)[0, 0]

    def coherent_sum(u, w):
        return sum(sp.conjugate(u[j]) * w[j] / a[j] for j in range(Ne))

    S1, S2 = coherent_sum(d1, d1), coherent_sum(d2, d2)
    K12, K21 = coherent_sum(d1, d2), coherent_sum(d2, d1)
    note = S1 - beta * K12 * K21 / (Gg + beta * S2)

    return (sp.simplify(sp.together(Xi - note)) == 0
            and sp.simplify(Sg - (Gg + beta * S2)) == 0)


def main():
    v = build()
    A, B, C, Gg, bp = v["A"], v["B"], v["C"], v["Gg"], v["bp"]

    direct = sp.Matrix([[A, B], [C, Gg]]).solve(sp.Matrix([bp, 0]))[0]

    # the manuscript's route: eliminate x, so G_g may go to zero safely
    Sg = Gg - C * A**-1 * B
    schur = A**-1 * bp + A**-1 * B * Sg**-1 * C * A**-1 * bp

    # Fleischhauer, Imamoglu and Marangos, Rev. Mod. Phys. 77, 633 (2005)
    beta = sp.simplify(v["Oc"] * sp.conjugate(v["Oc"])) / 4
    textbook = bp / (A + beta / Gg)

    checks = [
        ("Schur complement == direct solve",
         sp.simplify(direct - schur) == 0),
        ("Schur complement == textbook rho_e1",
         sp.simplify(sp.together(direct - textbook)) == 0),
    ]

    # the ideal limit, where the total response vanishes and the correction
    # does not
    ideal = {v["dp"]: 0, v["dc"]: 0, v["g21"]: 0}
    cut = bp / A                                    # B = C = 0
    chi_full = sp.simplify(direct.subs(ideal))
    dXi = sp.simplify(((direct - cut) / cut).subs(ideal))
    checks.append(("ideal limit: chi_full == 0", chi_full == 0))
    checks.append(("ideal limit: dXi_S == -1", dXi == -1))

    # NV is the N_e = 2 case; the others are there to show nothing about the
    # identity depends on the branch count
    for Ne in (1, 2, 3, 4):
        checks.append((f"multilevel N_e = {Ne} == theory note Sec. 9.7",
                       multilevel(Ne)))

    ok = True
    for name, passed in checks:
        print(f"  {'PASS' if passed else 'FAIL'}  {name}")
        ok = ok and bool(passed)

    print("\nrho_e1 (general) =")
    sp.pprint(sp.simplify(direct))
    print("\nrho_e1 (textbook form) = b_p / (A + beta/G_g),  beta = "
          "|Omega_c|^2/4")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
