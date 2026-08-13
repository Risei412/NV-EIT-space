"""Regression tests for Gate D (robustness + experimental discriminability).

Pins: NV n=3 is an exact (unbreakable) class; the superconducting protected
class breaks with a crossover Gamma*(eps) ~ 1/eps; Gamma(T) reach per platform;
optical detectability. Pytest-style, also standalone-runnable.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "numerics", "New no-go theory", "GateD_robustness_discriminability", "src")
PHASE_SRC = os.path.join(HERE, "..", "..", "numerics", "New no-go theory", "src")
NOGO_SRC = os.path.join(HERE, "..", "..", "numerics", "No-go theorem", "src")
GATEB_SRC = os.path.join(HERE, "..", "..", "numerics", "New no-go theory", "GateB_superconducting_witness", "src")
GATEC_SRC = os.path.join(HERE, "..", "..", "numerics", "New no-go theory", "GateC_material_independence", "src")
for _p in (SRC, PHASE_SRC, NOGO_SRC, GATEB_SRC, GATEC_SRC):
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


def test_representative_contrast_is_read_from_gate5_not_hard_coded():
    """Regression: the representative contrast used to be the literal 1e-2,
    which is the un-averaged single-defect value. It must now come from the
    ensemble table, and the gating scenario must be one of the averaged ones."""
    import experimental_budget as eb
    import run_gate_d as rd
    scen = eb.load_gate5_contrast()
    assert set(scen) >= {"single", "low_density", "high_density",
                         "post_selected", "post_selected_shimmed"}
    assert rd.GATING_SCENARIO != "single"
    snr = rd.optical_snr()
    assert snr["representative_contrast"] == scen[rd.GATING_SCENARIO]["Cmax"]
    assert abs(snr["representative_contrast"] - 1e-2) > 1e-3, \
        "gating contrast must not be the old hard-coded single-defect value"
    # every scenario is reported, not just the gating one
    assert set(snr["per_scenario"]) == set(scen)


def test_detection_floor_is_evaluated_at_the_candidate_temperature():
    """The contrasts are computed at the 70 K candidate point, so c_min has to
    be too; it used to be taken at 50 K, comparing two different operating
    points."""
    import run_gate_d as rd
    snr = rd.optical_snr()
    assert rd.CANDIDATE_T_K == 70.0
    assert snr["floor_at_candidate_T"]["T_K"] == 70.0
    assert snr["min_detectable_contrast"] == \
        snr["floor_at_candidate_T"]["min_detectable_contrast"]


def test_slope_precision_matches_the_closed_form():
    """Monte-Carlo slope spread must agree with sigma/(sqrt(N)*std(ln Gamma)).
    Two estimators of the same quantity, so a mistake in either shows up."""
    import experimental_budget as eb
    gs = np.logspace(3, 6, 40)
    vals = gs ** -2.0
    prec = eb.slope_precision(gs, vals, sigma_rel=0.05, n_trials=800, seed=7)
    assert abs(prec["std"] - prec["closed_form_std"]) / prec["closed_form_std"] < 0.15
    assert abs(prec["bias"]) < 0.02


def test_noisy_fit_needs_more_than_the_noiseless_one():
    """The old criterion fitted a noiseless kernel and was satisfied to 2e-5 by
    one decade. With noise the requirement must actually bind: a narrow window
    at the detection-threshold precision must NOT resolve an adjacent class."""
    import experimental_budget as eb
    gs = np.logspace(3, 3.3, 6)
    vals = gs ** -3.0
    prec = eb.slope_precision(gs, vals, sigma_rel=0.20, n_trials=400, seed=11)
    assert not prec["resolves_adjacent_class"], prec


def test_optical_nv_cannot_read_its_own_exponent_after_washout():
    """With nu = 4, one decade of Gamma costs four decades of signal, so the
    window above the detection floor closes before the fit can resolve the
    class. Pinning this keeps the limitation visible: the measurable-design-rule
    claim rests on engineered-dissipation platforms, not on optical NV."""
    import run_gate_d as rd
    budget = rd.slope_budget(rd.optical_snr())
    nv_rows = [r for r in budget["rows"] if r["platform"].startswith("NV optical")]
    gated = [r for r in nv_rows if rd.GATING_SCENARIO in r["platform"]]
    assert gated and not gated[0]["class_resolvable"], gated
    assert not budget["optical_nv_resolvable"] or \
        all(not r["class_resolvable"] for r in nv_rows if "single" not in r["platform"])
    # but the law is still measurable somewhere
    assert budget["any_platform_resolvable"]
    resolvable = [r["platform"] for r in budget["rows"] if r["class_resolvable"]]
    assert any("SC transfer" in p or "chain" in p for p in resolvable), resolvable


def test_sc_reach_uses_engineerable_range_not_sweep_width():
    """sc_decades_available was the literal 9.0, attributed to a sweep. Reach
    must come from the range a bus can be tuned over, which is far smaller."""
    import model_sc_transfer as m
    import run_gate_d as rd
    gT = rd.gamma_T_mapping()
    assert "sc_decades_available" not in gT
    assert gT["sc_decades_engineerable"] == m.decades(m.KAPPA_ENGINEERABLE)
    assert gT["sc_decades_engineerable"] < gT["sc_decades_swept"]
    assert gT["sc_decades_engineerable"] < 4.0


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
