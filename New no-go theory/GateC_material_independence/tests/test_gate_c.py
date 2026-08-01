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
import chain3_witness as c3
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


def test_chain3_is_an_exact_class_three_witness():
    """The non-diamond n=3 witness must be exact, not a good fit: the first two
    path moments vanish identically because the chain has no direct 1-3 link,
    and the only surviving path gives M2 = -J12*J23."""
    cert = c3.exact_certificate()
    assert cert["M0"] == "0" and cert["M1"] == "0"
    assert cert["M2"] == "-J12*J23"
    M = np.abs(c3.moments(3))
    assert M[0] < 1e-12 and M[1] < 1e-12
    assert abs(M[2] - c3.J12 * c3.J23) < 1e-12


def test_chain3_full_gksl_orders():
    """Full 16-dim GKSL: amplitude order 3 and fixed-readout population order 6,
    with the reduced kernel reproducing the full model to machine precision."""
    gammas = np.logspace(2, 4, 20)
    amps, pops = [], []
    for g in gammas:
        a, p = c3.full_response(g)
        amps.append(abs(a))
        pops.append(p)
    amp_order = -np.polyfit(np.log(gammas), np.log(amps), 1)[0]
    pop_order = -np.polyfit(np.log(gammas), np.log(pops), 1)[0]
    assert abs(amp_order - 3.0) < 0.05, amp_order
    assert abs(pop_order - 6.0) < 0.05, pop_order
    assert c3.full_vs_reduced_max_rel_err(np.logspace(2, 4, 5)) < 1e-7


def test_chain3_class_drops_when_the_topology_is_broken():
    """Material independence is a claim about topology, so breaking the topology
    must change the class. Opening a direct 1-3 bypass makes M1 nonzero and
    drops the amplitude order from 3 towards 2."""
    gammas = np.logspace(3, 6, 24)
    broken = np.array([c3.kernel(g, eps13=1e-2) for g in gammas])
    order = -np.polyfit(np.log(gammas), np.log(np.abs(broken)), 1)[0]
    assert abs(order - 2.0) < 0.1, order
    assert abs(c3.moments(3, eps13=1e-2)[1]) > 1e-9


def test_every_class_has_a_diamond_and_a_non_diamond_host():
    """Class 3 previously had only NV. Pin that each of n = 1,2,3 is now
    realized in both host families with the predicted slope."""
    gammas_nv = np.logspace(2, 5, 24)
    gammas_sc = np.logspace(5, 9, 24)
    gammas_gi = np.logspace(4, 8, 24)
    cases = [
        (1, "diamond", np.array([gf.full_response(g, "SiV", mode="dephasing") for g in gammas_gi]), gammas_gi),
        (1, "non-diamond", np.array([sc.transfer_kernel(g, tuning="generic") for g in gammas_sc]), gammas_sc),
        (2, "diamond", nvk.kernel(_NV_H, (0, -1), gammas_nv), gammas_nv),
        (2, "non-diamond", np.array([sc.transfer_kernel(g, tuning="protected") for g in gammas_sc]), gammas_sc),
        (3, "diamond", nvk.kernel(_NV_H, (-1, 1), gammas_nv), gammas_nv),
        (3, "non-diamond", np.array([c3.kernel(g) for g in gammas_nv]), gammas_nv),
    ]
    hosts = {}
    for n, host, vals, gs in cases:
        nu = core.fit_nu_loglog(gs, vals)["nu_global"]
        assert abs(nu - n) < 0.05, (n, host, nu)
        hosts.setdefault(n, set()).add(host)
    for n in (1, 2, 3):
        assert hosts[n] == {"diamond", "non-diamond"}, (n, hosts[n])


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
