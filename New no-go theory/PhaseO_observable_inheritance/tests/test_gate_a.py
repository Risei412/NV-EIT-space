"""Regression tests for Gate A observable-order inheritance.

Pins the known suppression orders so a future change to gate_a_observable.py
or the reused pencils is caught. Written pytest-style (bare asserts) but also
runnable standalone: `python tests/test_gate_a.py`.

Only the fast routes are exercised here (moment predictor, exact rational
certificate, Phase N q-fan / z*); the full log-log sweeps live in the runner
(run_gate_a.py) gate JSON.
"""
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
PHASE_SRC = os.path.join(HERE, "..", "..", "src")
NOGO = os.path.join(HERE, "..", "..", "..", "No-go theorem", "src")
ROOMT = os.path.join(HERE, "..", "..", "RoomT", "src")
for _p in (SRC, PHASE_SRC, NOGO, ROOMT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import gate_a_observable as gao
import model_specs as ms


def test_synthetic_three_classes_predict_and_certify():
    """generic / protected / doubly-protected -> nu_obs = 2 / 4 / 6, agreeing
    between the moment predictor and the exact rational-degree certificate."""
    expected = {"synthetic-generic (nu_obs=2)": 2,
                "synthetic-protected (nu_obs=4)": 4,
                "synthetic-doubly-protected (nu_obs=6)": 6}
    for spec in ms.synthetic_specs():
        pred = gao.predict_nu_obs(spec, z=0.0, kmax=8)
        cert = gao.certify_nu_obs_exact(spec, z_val=0.0)
        exp = expected[spec.label]
        assert pred["nu_obs_pred"] == exp, (spec.label, pred["nu_obs_pred"], exp)
        assert cert["nu_cert"] == exp, (spec.label, cert["nu_cert"], exp)
        # generic law nu_obs = n12 + n21 (g_eff != 0 => nu_den = 0)
        assert cert["nu_cert"] == cert["nu_K12"] + cert["nu_K21"]


def test_generic_vs_protected_selection_rule():
    """The generic case has a nonzero leading overlap (M0 != 0, n12 = 1); the
    protected cases have an EXACT M0 = dp^dag dc = 0 and n12 >= 2."""
    specs = {s.label: s for s in ms.synthetic_specs()}
    g = gao.predict_nu_obs(specs["synthetic-generic (nu_obs=2)"])
    assert g["mechanism"] == "generic" and g["n12"] == 1 and not g["protected"]
    for lbl in ("synthetic-protected (nu_obs=4)", "synthetic-doubly-protected (nu_obs=6)"):
        p = gao.predict_nu_obs(specs[lbl])
        assert p["mechanism"] == "symmetry-protected"
        assert p["protected"] and abs(p["m0_overlap"]) < 1e-9
        assert p["n12"] >= 2


def test_nv_predicts_order_four():
    """Physical NV EIT pencil: dp^dag dc = 0 -> n12 = 2 -> nu_obs = 4."""
    nv = ms.nv_spec()
    pred = gao.predict_nu_obs(nv, z=0.0, kmax=8)
    assert pred["n12"] == 2 and pred["n21"] == 2
    assert pred["nu_obs_pred"] == 4
    assert abs(pred["m0_overlap"]) < 1e-9


def test_nv_carries_its_own_exact_certificate():
    """The physical NV spec must certify nu_obs = 4 through Gate A's own
    rational-degree route, not only through the predictor. Constituent orders
    are pinned too, so a change that keeps nu_obs = 4 by cancelling errors is
    still caught."""
    nv = ms.nv_spec()
    cert = gao.certify_nu_obs_exact(nv, z_val=0.0)
    assert cert["nu_cert"] == 4, cert
    assert cert["nu_K12"] == 2 and cert["nu_K21"] == 2, cert
    assert cert["nu_S2"] == 1, cert


def test_nv_numeric_and_symbolic_pencil_agree():
    """The attached symbolic pencil must be the same object as the numeric one.

    nv_model.legs() builds dp/dc from an eigen-decomposition, while the
    certificate uses RoomT step3's exact rational pencil; nothing but this
    check ties the two together, and without it the certificate would be
    decorative."""
    nv = ms.nv_spec()
    chk = nv.meta["pencil_check"]
    assert chk["leg_overlap_ok"]
    assert chk["pencil_max_abs_err"] < 1e-9, chk


def test_check_model_does_not_pass_a_missing_certificate():
    """Regression: check_model used to compute
        pred_cert_agree = (nu_cert is None) or (nu_pred == nu_cert)
    so a spec with no symbolic pencil satisfied the certificate criterion
    vacuously and reported pred_cert_agree=True. A certificate-free spec must
    now report cert_present=False, pred_cert_agree=None, and fail 'agree'."""
    import numpy as np
    bare = ms.nv_spec(symbolic=False)
    assert bare.D_sym is None
    gammas = np.logspace(0, 10, 40)
    chk = gao.check_model(bare, gammas, z=0.0, kmax=8, require_cert=True)
    assert chk["cert_present"] is False
    assert chk["pred_cert_agree"] is None
    assert chk["fit_agree"] is True, "the fit itself should still be fine"
    assert chk["agree"] is False, "no certificate must not count as agreement"


def test_gate_a_cert_matches_roomt_step3():
    """Gate A's certificate and the independent RoomT step3 certificate are two
    routes to the same number; they must not drift apart."""
    import step3_merged_manifold_moments as step3
    analysis = step3.analyze_certificate(step3.build_symbolic_certificate())
    cert = gao.certify_nu_obs_exact(ms.nv_spec(), z_val=0.0)
    assert analysis["nu_R"] == cert["nu_cert"] == 4
    assert analysis["nu_K12"] == cert["nu_K12"]
    assert analysis["nu_S2"] == cert["nu_S2"]


def test_qfan_is_v_shaped():
    """Phase N intervention-scaling fan nu(q) = 4-q, 2+q, 4."""
    import phase_n_exact_core as pn
    num, den = pn.master_polynomials()
    for q, exp in [(Fraction(0), Fraction(4)), (Fraction(1), Fraction(3)),
                   (Fraction(3, 2), Fraction(7, 2)), (Fraction(2), Fraction(4)),
                   (Fraction(3), Fraction(4))]:
        order, _cn, _cd = pn.path_order(q, num, den)
        assert order == exp, (q, order, exp)


def test_zstar_promotes_order():
    """Isolated resonance z* = 543/280 promotes the S34 order 4 -> 5."""
    import phase_n_frequency_core as fq
    moments = fq.ideal_moment_polynomials(kmax=7)
    assert fq.first_nonzero_moment_at(moments, Fraction(1)) == 4
    assert fq.first_nonzero_moment_at(moments, Fraction(543, 280)) == 5


def test_frozen_source_promotion():
    """Phase P frozen-source difference: cancellation promotes order 3 -> 4."""
    import model_physical as mp
    import sympy as sp
    nu_generic, _Q, _N = mp.certificate_R_S(0, sp.Rational(1, 2), 0)
    nu_tuned, _Q2, _N2 = mp.certificate_R_S(0, mp.J45_star(), sp.pi)
    assert nu_generic == 3
    assert nu_tuned == 4


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
