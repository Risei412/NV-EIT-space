"""Step 3 of the room-temperature SMRT no-go plan (plan Sec. 5, Step 3):
merged-manifold moment analysis of the NV branch-resolved spin-Lambda's
EIT sector, on the SAME model validated in Step 1/2 (No-go theorem/src/
nv_model.py: probe on orbital branch X / ms=0, control on orbital branch
Y / ms=-1, Bx=By=0). Gamma is the native excited-state dissipation scale
(gamma_oc_GHz(T,d) in nv_model.py's convention -- physically the
phonon-driven orbital-branch merging rate that grows with T).

At Bx=By=0, nv_model.Hgs(Bvec) is EXACTLY diagonal (Sz, Sz@Sz are diagonal
and the only other term is Bz*Sz, also diagonal), so the ground "dressed"
eigenvectors are exactly the bare |ms> basis vectors -- no numerical
eigendecomposition is needed for this analysis, dp = e_1 and dc = e_3 in
the 6-dim (orbital x spin) excited-manifold basis are EXACT, closed-form
objects. This makes the whole moment/certificate calculation exact
symbolic algebra rather than a numerical approximation.

Three exact results are certified here, using sympy's adjugate/determinant
method (`certificate_deg_nu`-style, New no-go theory/src/core.py Theorem
8.1), NOT a small-float inference:

  M0 = dp^T dc = 0 EXACTLY (trivial orthogonality: dp, dc live on
       orthogonal orbital branches X, Y at leading order -- the leading
       probe-control interference moment vanishes by a group-theoretic
       selection rule, not a numerical coincidence, and holds regardless
       of the values of any Hamiltonian parameter).

  nu_K = 2: the first NONZERO moment of the kernel K12(Gamma) = dp^dag
       G(Gamma) dc is at order Gamma^-2 (M1 != 0, generic in the physical
       parameters Dperp, Lperp -- kept as free symbols in the certificate
       specifically to show the nonvanishing is structural, not a
       coincidence of the literature-fitted values).

  nu_R = 4: the EIT observable itself, R_EIT ~ dXi(Gamma) = -beta*K12*K21
       / (geff + beta*S2), inherits nu_K12 + nu_K21 = 2 + 2 = 4 (NOT the
       naive nu_K = 2), because R_EIT enters the response as the PRODUCT
       K12*K21 while the denominator's leading term (geff, the ground
       decoherence, degree-0 in Gamma) dominates over its Gamma-dependent
       correction (S2 ~ Gamma^-1) as Gamma -> infinity -- exactly the
       nu_K != nu_R distinction the plan's Sec. 3 (Step 3) warns about.

Gate 3 (plan Sec. 5, Step 3) requires this be shown via an exact identity
or symmetry (not "the float was 1e-14"), the first nonzero moment and its
coefficient identified, and a numeric log-log fit confirming the same
exponent independent of precision -- all four are produced below.

Outputs: RoomT/results/gates_summary_step3.json,
         RoomT/results/figures/fig_step3_moment_scaling.png
"""
from __future__ import annotations
import os
import sys
import json
import time

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOMT_ROOT = os.path.dirname(HERE)
NOGO_SRC = os.path.join(ROOMT_ROOT, "..", "..", "No-go theorem", "src")
PHASE_SRC = os.path.join(ROOMT_ROOT, "..", "src")
sys.path.insert(0, NOGO_SRC)
sys.path.insert(0, PHASE_SRC)

import nv_model as nm  # noqa: E402
import core  # noqa: E402  (New no-go theory/src/core.py: fit_nu_loglog)

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

D_STRAIN = 1.683
REFERENCE_OC = 1.0
REFERENCE_GG = 6.3e-5


def kron(A, B):
    return sp.Matrix(sp.BlockMatrix([[A[i, j] * B for j in range(A.cols)] for i in range(A.rows)]))


def build_symbolic_certificate():
    """Exact adjugate/determinant certificate for K12, K21, S2 on the
    Bx=By=0 excited-manifold Hamiltonian, with Dperp, Lperp kept as FREE
    symbols (all other parameters at their exact literature-fitted
    rational values, matching nv_model.py) so the moment-order result is
    proven generic in exactly the two couplings that determine it (see
    module docstring)."""
    Gamma = sp.symbols("Gamma")
    Dperp, Lperp = sp.symbols("Dperp Lperp", positive=True)

    s2 = 1 / sp.sqrt(2)
    Sz = sp.diag(-1, 0, 1)
    Sx = s2 * sp.Matrix([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    Sy = s2 * sp.Matrix([[0, sp.I, 0], [-sp.I, 0, sp.I], [0, -sp.I, 0]])
    I3 = sp.eye(3)
    sz_o = sp.diag(-1, 1)
    sx_o = sp.Matrix([[0, 1], [1, 0]])
    sy_o = sp.Matrix([[0, sp.I], [-sp.I, 0]])
    I2 = sp.eye(2)

    Dpar = sp.Rational(142, 100)
    Lpar = sp.Rational(533, 100)
    dx = sp.nsimplify(D_STRAIN, rational=True)  # phi=0 => dx=d, dy=0

    H = kron(I2, Dpar * (Sz * Sz - sp.Rational(2, 3) * I3))
    H += -Lpar * kron(sy_o, Sz)
    H += Dperp * (kron(sz_o, Sy * Sy - Sx * Sx) - kron(sx_o, Sx * Sz + Sz * Sx))
    H += Lperp * (kron(sz_o, Sx * Sz + Sz * Sx) - kron(sx_o, Sy * Sz + Sz * Sy))
    H += dx * kron(sz_o, I3)
    H = sp.Matrix(H)

    dp = sp.Matrix([0, 1, 0, 0, 0, 0])  # e_1: probe leg (orbital X, ms=0)
    dc = sp.Matrix([0, 0, 0, 1, 0, 0])  # e_3: control leg (orbital Y, ms=-1)

    M0 = sp.simplify((dp.T * dc)[0, 0])

    A = Gamma * sp.eye(6) + sp.I * 2 * sp.pi * H  # z = 0
    Q = sp.expand(A.det())
    Qp = sp.Poly(Q, Gamma)
    adjA = A.adjugate()
    N12 = sp.Poly(sp.expand((dp.T * adjA * dc)[0, 0]), Gamma)
    N21 = sp.Poly(sp.expand((dc.T * adjA * dp)[0, 0]), Gamma)
    N22 = sp.Poly(sp.expand((dc.T * adjA * dc)[0, 0]), Gamma)

    H13_exact = sp.simplify(H[1, 3])  # exact closed form of the first nonzero moment's generator

    return dict(
        Gamma=Gamma, Dperp=Dperp, Lperp=Lperp, H13_exact=H13_exact,
        M0=M0, Qp=Qp, N12=N12, N21=N21, N22=N22,
    )


def analyze_certificate(cert):
    Qp, N12, N21, N22 = cert["Qp"], cert["N12"], cert["N21"], cert["N22"]
    nu_K12 = Qp.degree() - N12.degree()
    nu_K21 = Qp.degree() - N21.degree()
    nu_S2 = Qp.degree() - N22.degree()

    # R_EIT ~ dXi = -beta*N12*N21 / (Q*(geff*Q + beta*N22)); geff (ground
    # decoherence) is Gamma-independent and nonzero, so deg(Q) > deg(N22)
    # (nu_S2 >= 1, confirmed above) guarantees the geff*Q term dominates
    # the denominator as Gamma -> infinity, independent of the specific
    # positive values of geff and beta -- kept symbolic here for that
    # reason (see module docstring).
    deg_num = N12.degree() + N21.degree()
    deg_den = Qp.degree() + max(Qp.degree(), N22.degree())
    nu_R = deg_den - deg_num

    # leading coefficient of R_EIT ~ coeff * Gamma^-nu_R, exact, as a
    # function of the (symbolic) beta, geff and the leading coefficients
    # of N12, N21, Q (LC(N22) does not enter at leading order since
    # deg(Q) > deg(N22)).
    beta_s, geff_s = sp.symbols("beta geff", positive=True)
    lead_coeff = sp.simplify(
        -beta_s * N12.LC() * N21.LC() / (geff_s * Qp.LC() ** 2)
    )

    class_ = "Class III (protected, REJECTS no-go)" if nu_R == 0 else (
        "Class I (exact zero to all checked orders)" if N12.as_expr() == 0 and N21.as_expr() == 0
        else "Class II (finite suppression order)"
    )

    return dict(
        M0=str(cert["M0"]), M0_is_exact_zero=bool(cert["M0"] == 0),
        H13_exact_closed_form=str(cert["H13_exact"]),
        deg_Q=Qp.degree(), deg_N12=N12.degree(), deg_N21=N21.degree(), deg_N22=N22.degree(),
        nu_K12=nu_K12, nu_K21=nu_K21, nu_S2=nu_S2, nu_R=nu_R,
        leading_coefficient_of_R_EIT_tail=str(lead_coeff),
        smrt_class=class_,
    )


def numeric_cross_check(N=300):
    """Independent numeric log-log fit of |K12|, |K12*K21|, |dXi| vs
    Gamma using the ACTUAL nv_model.py machinery (double precision, the
    same code path Step 1/2 validated), reused as a cross-check on the
    exact symbolic degrees above -- two independent estimators per Gate 3
    (New no-go theorem convention: symbolic exact + numeric asymptotic).

    IMPORTANT: R_EIT's denominator den = geff + beta*S2(Gamma) has a
    genuine CROSSOVER, not just a fit-window artifact. Since nu_S2 = 1
    (S2 ~ Gamma^-1) while geff is Gamma-independent, den is dominated by
    beta*S2 (NOT geff) for Gamma << Gamma_cross ~ beta/geff, giving a
    PRE-ASYMPTOTIC nu_R,eff ~ 3 (den ~ Gamma^-1 cancels one power of the
    Gamma^-4 naive product estimate) that only crosses over to the true
    asymptotic nu_R = 4 for Gamma >> Gamma_cross. At the Step 1 reference
    point (Oc=1 GHz, gg=6.3e-5 GHz), Gamma_cross ~ 7.5e4 GHz -- ABOVE the
    physical Gamma(300 K) ~ 1.3e4 GHz (see `crossover_diagnostics`) -- so
    the numeric fit here deliberately uses a deep tail (Gamma > 1e7,
    Gamma_cross) to confirm the TRUE Gamma -> infinity degree against the
    exact symbolic certificate; the pre-asymptotic ~3 regime and its
    physical relevance at 300 K is reported separately for Step 4."""
    Bvec = (0.0, 0.0, 0.02)
    w, U = nm.dressed_ground(Bvec)
    H = nm.Hes(Bvec, d=D_STRAIN)
    dp, dc = nm.legs(U)
    Ep, fid, _ = nm.probe_line(H, dp)

    gammas = np.logspace(0, 10, N)
    K12s, K21s, dXis = [], [], []
    beta = (nm.TWOPI * REFERENCE_OC) ** 2 / 4
    geff = 2 * REFERENCE_GG + 2e-6
    for g in gammas:
        G = np.linalg.inv(g * np.eye(6) + 1j * nm.TWOPI * (H - 0.0 * np.eye(6)))
        K12 = np.vdot(dp, G @ dc)
        K21 = np.vdot(dc, G @ dp)
        S2 = np.vdot(dc, G @ dc)
        den = geff + beta * S2
        dXi = -beta * K12 * K21 / den
        K12s.append(K12); K21s.append(K21); dXis.append(dXi)
    K12s, K21s, dXis = map(np.array, (K12s, K21s, dXis))

    # K12 alone has no crossover (nu_K12=2 holds from Gamma>~1 already);
    # R_EIT needs the deep tail (Gamma > Gamma_cross) for its TRUE
    # asymptotic degree.
    tail_K = gammas > 1e3
    tail_R = gammas > 1e8
    fit_K12 = core.fit_nu_loglog(gammas[tail_K], K12s[tail_K])
    fit_K12K21 = core.fit_nu_loglog(gammas[tail_K], K12s[tail_K] * K21s[tail_K])
    fit_R = core.fit_nu_loglog(gammas[tail_R], dXis[tail_R])

    # crossover diagnostics: local effective index of R_EIT vs Gamma, and
    # where it crosses the 3<->4 midpoint (3.5), compared to the actual
    # physical Gamma(T) at the reference points used in Step 1.
    full_fit = core.fit_nu_loglog(gammas, dXis)
    nu_eff, gmid = full_fit["nu_eff"], full_fit["gamma_mid"]
    cross_idx = int(np.argmin(np.abs(nu_eff - 3.5)))
    gamma_cross = float(gmid[cross_idx])
    gamma_300K = nm.gamma_oc_GHz(300.0, D_STRAIN)
    gamma_10K = nm.gamma_oc_GHz(10.0, D_STRAIN)

    return dict(
        gammas=gammas, K12s=K12s, K21s=K21s, dXis=dXis,
        nu_K12_fit=fit_K12["nu_global"], nu_K12K21_fit=fit_K12K21["nu_global"],
        nu_R_fit=fit_R["nu_global"],
        crossover_diagnostics=dict(
            nu_eff_at_small_gamma=float(nu_eff[0]),
            nu_eff_at_large_gamma=float(nu_eff[-1]),
            gamma_cross_GHz=gamma_cross,
            gamma_300K_GHz=float(gamma_300K),
            gamma_10K_GHz=float(gamma_10K),
            note_300K_below_crossover=bool(gamma_300K < gamma_cross),
        ),
    )


def precision_stability_check():
    """Same tail fit at float64 vs a higher-precision (long double via
    numpy longdouble is unavailable for complex; use finer/wider Gamma
    sampling plus a residual check instead) to confirm the exponent
    estimate is not a float64 artifact -- Gate 3's "精度を変えても指数が
    変わらない" requirement. Uses the EXACT symbolic Q/N12/N21/N22
    polynomials evaluated at high (50-digit) mpmath precision as the
    higher-precision branch, cross-checked against the float64 numeric
    fit above."""
    import mpmath as mp
    mp.mp.dps = 50

    cert = build_symbolic_certificate()
    Qp, N12, N21, N22 = cert["Qp"], cert["N12"], cert["N21"], cert["N22"]
    Dperp_val = sp.Rational(155, 200)          # 1.55/2
    Lperp_val = sp.sqrt(2) * sp.Rational(1, 10)  # 0.20/sqrt(2) rationalized

    def subs_eval(poly, Gamma_val):
        expr = poly.as_expr().subs({cert["Dperp"]: Dperp_val, cert["Lperp"]: Lperp_val,
                                     cert["Gamma"]: Gamma_val})
        return complex(sp.N(expr, 50))

    gammas = [mp.mpf(10) ** k for k in range(3, 9)]
    K12_hp = []
    for g in gammas:
        q = subs_eval(Qp, g)
        n12 = subs_eval(N12, g)
        K12_hp.append(n12 / q)
    K12_hp = np.array(K12_hp, dtype=complex)
    gammas_f = np.array([float(g) for g in gammas])
    fit_hp = core.fit_nu_loglog(gammas_f, K12_hp)
    return dict(nu_K12_highprecision_fit=fit_hp["nu_global"])


def run():
    t0 = time.time()
    cert = build_symbolic_certificate()
    analysis = analyze_certificate(cert)
    print(f"symbolic certificate built in {time.time()-t0:.1f} s")

    numeric = numeric_cross_check()
    precision = precision_stability_check()

    gates = dict(
        M0_exact_zero=bool(analysis["M0_is_exact_zero"]),
        nu_K12_symbolic_equals_2=bool(analysis["nu_K12"] == 2),
        nu_R_symbolic_equals_4=bool(analysis["nu_R"] == 4),
        nu_K_neq_nu_R=bool(analysis["nu_K12"] != analysis["nu_R"]),
        numeric_fit_matches_symbolic_K12=bool(abs(numeric["nu_K12_fit"] - analysis["nu_K12"]) < 0.05),
        numeric_fit_matches_symbolic_R=bool(abs(numeric["nu_R_fit"] - analysis["nu_R"]) < 0.05),
        precision_stable=bool(abs(precision["nu_K12_highprecision_fit"] - analysis["nu_K12"]) < 0.05),
        not_class_III=bool("Class III" not in analysis["smrt_class"]),
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        model="No-go theorem/src/nv_model.py, Bx=By=0 exact-diagonal ground manifold",
        exact_symbolic_certificate=analysis,
        numeric_cross_check=dict(
            nu_K12_fit=numeric["nu_K12_fit"],
            nu_K12K21_fit=numeric["nu_K12K21_fit"],
            nu_R_fit=numeric["nu_R_fit"],
        ),
        crossover_diagnostics_flag_for_step4=numeric["crossover_diagnostics"],
        precision_stability=precision,
        gates=gates,
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step3.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure ---
    fig, ax = plt.subplots(figsize=(7, 5.5))
    g = numeric["gammas"]
    ax.loglog(g, np.abs(numeric["K12s"]), label="|K12(Gamma)|")
    ax.loglog(g, np.abs(numeric["K12s"] * numeric["K21s"]), label="|K12*K21(Gamma)|")
    ax.loglog(g, np.abs(numeric["dXis"]), label="|R_EIT(Gamma)| = |dXi|")
    tail = g > 1e3
    tail_deep = g > 1e8
    ref2 = np.abs(numeric["K12s"])[tail][0] * (g[tail][0] / g[tail]) ** 2
    ref4 = np.abs(numeric["dXis"])[tail_deep][0] * (g[tail_deep][0] / g[tail_deep]) ** 4
    ax.loglog(g[tail], ref2, "k--", label="Gamma^-2 reference")
    ax.loglog(g[tail_deep], ref4, "k:", label="Gamma^-4 reference (true asymptote)")
    gc = numeric["crossover_diagnostics"]["gamma_cross_GHz"]
    g300 = numeric["crossover_diagnostics"]["gamma_300K_GHz"]
    ax.axvline(gc, color="gray", linestyle="-.", alpha=0.6, label=f"Gamma_cross~{gc:.1e}")
    ax.axvline(g300, color="red", linestyle="-.", alpha=0.6, label=f"Gamma(300K)~{g300:.1e}")
    ax.set_xlabel("Gamma (GHz)"); ax.set_ylabel("magnitude")
    ax.set_title("Step 3: merged-manifold moment scaling (nu_K=2, nu_R=4 asymptotic)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step3_moment_scaling.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps(gates, indent=2))
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
