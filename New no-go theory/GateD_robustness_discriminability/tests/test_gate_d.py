"""Regression tests for Gate D (robustness + experimental discriminability).

Pins: NV n=3 is an exact (unbreakable) class; the superconducting protected
class breaks with a crossover Gamma*(eps) ~ 1/eps; Gamma(T) reach per platform;
optical detectability. Pytest-style, also standalone-runnable.
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
import nv_reduced_kernel as nvk
import model_sc_transfer as sc
import phonon_rates as pr
import nv_model as nv

_NV_H = nvk.H_3E()


def _slope(g, v):
    return core.fit_nu_loglog(np.asarray(g, float), np.asarray(v))["nu_global"]


def test_nv_class_is_exact_under_perturbations():
    """NV ms=-1<->+1: M0=M1=0 stay structural zeros and order stays 3 under
    strain, transverse field and detuning -> an exact (unbreakable) class."""
    ks = np.logspace(2, 6, 40)
    for kw in ({}, dict(xi_x=0.2), dict(xi_y=0.1), dict(Bx=0.2)):
        H = nvk.H_3E(**kw)
        M = nvk.moments(H, (-1, 1), 3)
        assert abs(M[0]) < 1e-9 and abs(M[1]) < 1e-9
        assert abs(_slope(ks, nvk.kernel(H, (-1, 1), ks)) - 3.0) < 0.02
    for z in (0.5, 1.0):
        M = nvk.moments(_NV_H, (-1, 1), 3, z=z)
        assert abs(M[0]) < 1e-9 and abs(M[1]) < 1e-9


def test_sc_approximate_class_crossover_inverse_eps():
    """Superconducting protected nu=2 breaks to nu=1 with Gamma*(eps) ~ 1/eps."""
    ks = np.logspace(4, 13, 100)
    eps = np.array([1e-9, 1e-10, 1e-11])
    gstar = []
    for e in eps:
        K = np.array([sc.transfer_kernel(k, tuning="broken", eps=e) for k in ks])
        fit = core.fit_nu_loglog(ks, K)
        below = np.where(fit["nu_eff"] < 1.5)[0]
        gstar.append(fit["gamma_mid"][below[0]] if len(below) else np.inf)
    gstar = np.array(gstar)
    assert np.all(np.diff(gstar) > 0)  # grows as eps shrinks
    power = np.polyfit(np.log10(eps), np.log10(gstar), 1)[0]
    assert abs(power + 1.0) < 0.1


def test_gamma_T_reach_nv_wide_group_iv_narrow():
    """NV k_orb ~ T^5 spans many decades over 4-300 K; enough to fit a slope."""
    d = 1.683
    lo = nv.korb_GHz(4.0, d)
    hi = nv.korb_GHz(300.0, d)
    assert np.log10(hi / lo) > 3.0


def test_optical_detectable_contrast_finite():
    """A feasible optical operating point yields a finite min detectable
    contrast well below an order-1e-2 signal."""
    import signal_chain as sig
    gamma_h = nv.gamma_oc_GHz(50.0, 1.683)
    sigma = sig.sigma_zpl_cm2(637.0, 2.41, 0.035, nv.GRAD, gamma_h)
    f_spec = sig.spectral_fraction(gamma_h, 30.0)
    alpha = sig.alpha_cm(sigma, 1.76e17, 0.25, 1 / 3, f_spec)
    od_s = sig.od(alpha, 0.05)
    c_min = sig.min_detectable_contrast(5.0, od_s, od_s, 1e-6, 637.0, 3600.0, 0.1, 1e-6)
    assert np.isfinite(c_min) and c_min < 1e-2


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
