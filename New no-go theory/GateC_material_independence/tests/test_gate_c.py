"""Regression tests for Gate C (material independence).

Pins the group-IV full-GKSL Gamma^-1 certification and the three physical
suppression classes n=1,2,3. Pytest-style, also standalone-runnable.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
PHASE_SRC = os.path.join(HERE, "..", "..", "src")
NOGO_SRC = os.path.join(HERE, "..", "..", "..", "No-go theorem", "src")
GATEB_SRC = os.path.join(HERE, "..", "..", "GateB_superconducting_witness", "src")
for _p in (SRC, PHASE_SRC, NOGO_SRC, GATEB_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import core
import group_iv_full as gf
import nv_reduced_kernel as nvk
import model_sc_transfer as sc

_NV_H = nvk.H_3E()


def _slope(g, v):
    return core.fit_nu_loglog(np.asarray(g, float), np.asarray(v))["nu_global"]


def test_group_iv_full_gksl_gamma_inverse():
    """group-IV full GKSL: slope -1 in both dephasing and hopping modes."""
    for material, (lo, hi) in (("SiV", (4, 8)), ("SnV", (5, 10))):
        g = np.logspace(lo, hi, 20)
        for mode in ("dephasing", "hopping"):
            Kf = np.array([gf.full_response(x, material, mode=mode) for x in g])
            assert abs(_slope(g, Kf) - 1.0) < 1e-2, (material, mode)


def test_group_iv_reduced_equals_full():
    """Full-GKSL coherence sector == reduced kernel to machine precision."""
    for material in ("SiV", "SnV"):
        for x in (1e4, 1e6, 1e8):
            Kf = gf.full_response(x, material) / (-1j)
            Kr = gf.reduced_kernel_response(x, material)
            assert abs(Kf - Kr) / abs(Kr) < 1e-7, material


def test_group_iv_M0_nonzero_gamma_times_R():
    """M0 = p^dag c != 0 and Gamma*R(Gamma) -> M0."""
    for material in ("SiV", "SnV"):
        m0 = gf.M0(material)
        assert abs(m0) > 1e-8
        R = gf.full_response(1e8, material) / (-1j)
        assert abs(1e8 * R - m0) / abs(m0) < 1e-3


def test_nv_two_and_three_classes():
    """NV ms=0<->-1 -> Gamma^-2, ms=-1<->+1 -> Gamma^-3."""
    g = np.logspace(2, 5, 31)
    assert abs(_slope(g, nvk.kernel(_NV_H, (0, -1), g)) - 2.0) < 0.05
    assert abs(_slope(g, nvk.kernel(_NV_H, (-1, 1), g)) - 3.0) < 0.05


def test_material_independence_class_1_and_2():
    """The same integer class appears in diamond and non-diamond (SC)."""
    gd = np.logspace(4, 8, 20)
    gsc = np.logspace(5, 9, 20)
    # class 1
    assert abs(_slope(gd, [gf.full_response(x, "SiV") for x in gd]) - 1.0) < 0.05
    assert abs(_slope(gsc, [sc.transfer_kernel(x, tuning="generic") for x in gsc]) - 1.0) < 0.05
    # class 2
    gnv = np.logspace(2, 5, 31)
    assert abs(_slope(gnv, nvk.kernel(_NV_H, (0, -1), gnv)) - 2.0) < 0.05
    assert abs(_slope(gsc, [sc.transfer_kernel(x, tuning="protected") for x in gsc]) - 2.0) < 0.05


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
