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

Three things are checked, and all three are exact identities rather than
numerical agreement:

  1. eliminating x through the Schur complement gives the same rho_e1 as
     solving the 2x2 system directly -- i.e. the "never divide by G_g" route
     costs nothing;
  2. that rho_e1 is the textbook EIT coherence, with the control field
     entering as the continued-fraction term |Omega_c|^2/4 over the two-photon
     denominator;
  3. at perfect two-photon resonance with no ground dephasing, chi_full = 0
     while the sector-resolved correction is maximal, dXi_S = -1.

Point 3 is the paper's central conceptual claim -- that chi_full = 0 is what
perfect EIT looks like, so a criterion phrased on the total susceptibility
would report ideal EIT as an absence of EIT -- and here it is a two-line
consequence rather than an assertion.

Usage: python verify_lambda_reduction.py     (exit 0 iff all three hold)
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
