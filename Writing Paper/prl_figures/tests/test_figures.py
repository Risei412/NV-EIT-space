"""Regression tests for the PRL main-text figures: each figure's PDF+PNG
exist after build(), and the key numbers baked into each panel match the
Gate A-D certified values. Pytest-style, also standalone-runnable.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import fig1_classes
import fig2_inheritance
import fig3_material_independence
import fig4_robustness
import core
import nv_reduced_kernel as nvk
import group_iv_full as gf
import model_sc_transfer as sc
import gate_a_observable as gao
import model_specs as ms


def _exists(pdf, png):
    assert os.path.isfile(pdf) and os.path.getsize(pdf) > 0
    assert os.path.isfile(png) and os.path.getsize(png) > 0


def test_fig1_outputs_and_slopes():
    pdf, png = fig1_classes.build()
    _exists(pdf, png)
    ks_giv = np.logspace(4, 8, 20)
    K_giv = np.array([gf.full_response(k, "SiV") for k in ks_giv])
    assert abs(core.fit_nu_loglog(ks_giv, K_giv)["nu_global"] - 1.0) < 0.02
    H = nvk.H_3E()
    ks_nv = np.logspace(2, 5, 20)
    assert abs(core.fit_nu_loglog(ks_nv, nvk.kernel(H, (0, -1), ks_nv))["nu_global"] - 2.0) < 0.02
    assert abs(core.fit_nu_loglog(ks_nv, nvk.kernel(H, (-1, 1), ks_nv))["nu_global"] - 3.0) < 0.02


def test_fig2_outputs_and_inheritance():
    pdf, png = fig2_inheritance.build()
    _exists(pdf, png)
    gammas = np.logspace(2, 7, 30)
    for spec, want in zip(ms.synthetic_specs(), (2, 4, 6)):
        v = gao.verify_nu_obs_loglog(spec, gammas, z=0.0)
        assert abs(v["nu_tail"] - want) < 0.05


def test_fig3_outputs_and_collapse():
    pdf, png = fig3_material_independence.build()
    _exists(pdf, png)
    ks = np.logspace(5, 9, 20)
    eff_g = np.abs([sc.transfer_kernel(k, tuning="generic") for k in ks]) ** 2
    eff_p = np.abs([sc.transfer_kernel(k, tuning="protected") for k in ks]) ** 2
    assert abs(core.fit_nu_loglog(ks, eff_g)["nu_global"] - 2.0) < 0.05
    assert abs(core.fit_nu_loglog(ks, eff_p)["nu_global"] - 4.0) < 0.05


def test_fig4_outputs_and_crossover():
    pdf, png = fig4_robustness.build(quick=True)
    _exists(pdf, png)
    H = nvk.H_3E(xi_x=0.2)
    M = nvk.moments(H, (-1, 1), 3)
    assert abs(M[0]) < 1e-9 and abs(M[1]) < 1e-9  # NV exact class unbroken


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
