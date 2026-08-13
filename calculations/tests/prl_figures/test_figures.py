"""Regression tests for the PRL main-text figures: each figure's PDF+PNG
exist after build(), and the key numbers baked into each panel match the
Gate A-D certified values. Pytest-style, also standalone-runnable.
"""
import os
import sys

import numpy as np

# sys.path is set up by conftest.py, so these imports no longer depend on
# fig1_classes having been imported first for its side effects.
HERE = os.path.dirname(os.path.abspath(__file__))

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
import chain3_witness as c3
import run_gate_d as rd


def test_fig2b_crossover_scales():
    """Fig. 2(b) writes Gamma_cross and Gamma(300 K) into the panel as vertical
    lines, and the whole message is that the physical operating point sits
    BELOW the crossover -- i.e. NV at 300 K is pre-asymptotic. None of those
    numbers was covered by a test."""
    nv = ms.nv_spec()
    reg = gao.separate_regimes(
        nv, np.logspace(0, 10, 140), z=0.0,
        gamma_phys=nv.meta["gamma_phys_300K"],
        generic_order=3.0, asymptotic_order=4.0)
    assert abs(reg["gamma_cross"] / 7.18e4 - 1.0) < 0.03, reg["gamma_cross"]
    assert abs(nv.meta["gamma_phys_300K"] / 1.335e4 - 1.0) < 0.01
    assert reg["is_preasymptotic"] is True
    assert reg["gamma_cross"] > nv.meta["gamma_phys_300K"]
    assert abs(reg["nu_eff_at_phys"] - 3.15) < 0.1, reg["nu_eff_at_phys"]
    assert abs(reg["nu_asymptotic"] - 4.0) < 0.01, reg["nu_asymptotic"]


def test_fig4a_crossover_power_label():
    """Fig. 4(a)'s inset prints the fitted exponent of Gamma*(eps) ~ eps^p.
    Pin p = -1 and the monotonicity the fan is drawn from."""
    res = rd.sc_approximate_class(quick=True)
    assert abs(res["crossover_power"] + 1.0) < 0.05, res["crossover_power"]
    assert res["law_is_inverse"]
    gstars = [r["gamma_star"] for r in res["rows"]]
    assert all(gstars[i] < gstars[i + 1] for i in range(len(gstars) - 1)), gstars
    assert f"{res['crossover_power']:.2f}" in ("-1.00", "-0.99", "-1.01")


def test_fig4b_platform_reach_labels():
    """Fig. 4(b) labels NV's Gamma(T) span as ~8.5 decades and shows group-IV
    as narrower."""
    gT = rd.gamma_T_mapping()
    assert round(gT["nv_decades"], 1) == 8.5, gT["nv_decades"]
    assert gT["siv_decades"] < gT["nv_decades"]
    assert gT["snv_decades"] < gT["siv_decades"]


def test_fig4b_required_range_band_comes_from_gate_d():
    """The required-Gamma-range band must be read from Gate D rather than drawn
    from a literal, so the figure cannot drift away from the gate."""
    req = rd.required_gamma_range()
    assert req["delta_nu_resolvable"]
    assert 1 <= req["decades_needed"] <= 4, req["decades_needed"]


def test_fig3a_collapse_is_exponent_specific():
    """Fig. 3(a) claims a specific integer, not just a power law. Compensating
    by Gamma^n must flatten each system, and by Gamma^(n+-1) must not."""
    for n, label, host, ks, K in fig3_material_independence._systems():
        for dn, flat in ((0, True), (-1, False), (1, False)):
            comp = ks ** (n + dn) * K
            spread = float(np.std(comp) / np.mean(comp))
            if flat:
                assert spread < 2e-2, (label, dn, spread)
            else:
                assert spread > 0.5, (label, dn, spread)


def test_fig3a_every_class_has_both_hosts():
    """Class 3 used to appear in Fig. 3(a) as NV only."""
    hosts = {}
    for n, label, host, ks, K in fig3_material_independence._systems():
        hosts.setdefault(n, set()).add(host)
    for n in (1, 2, 3):
        assert hosts[n] == {"diamond", "non-diamond"}, (n, hosts.get(n))


def test_fig3a_chain_witness_is_class_three():
    """The non-diamond class-3 curve in Fig. 3(a) really is order 3."""
    ks = np.logspace(2, 5, 40)
    K = np.array([c3.kernel(k) for k in ks])
    nu = core.fit_nu_loglog(ks, K)["nu_global"]
    assert abs(nu - 3.0) < 0.05, nu


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
