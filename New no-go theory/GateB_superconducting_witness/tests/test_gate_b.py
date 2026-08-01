"""Regression tests for Gate B (superconducting non-EIT witness).

Pins the blind-predicted integer orders and the reduced==full agreement so a
future change to model_sc_transfer.py is caught. Pytest-style but also
standalone-runnable: `python tests/test_gate_b.py`.
"""
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
PHASE_SRC = os.path.join(HERE, "..", "..", "src")
for _p in (SRC, PHASE_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import core
import model_sc_transfer as m


def test_blind_certificate_orders():
    """generic -> nu_K = 1, protected -> nu_K = 2 (exact rational degree,
    valid for the singular bus-dissipation D)."""
    Gamma = sp.symbols("Gamma")
    exp = {"generic": 1, "protected": 2}
    for tuning, want in exp.items():
        D_sym, B_sym_of_z, c_sym, p_sym = m.symbolic_pencil(tuning)
        nu_K, _Q, _N = core.certificate_deg_nu(D_sym, B_sym_of_z, c_sym, p_sym,
                                               Gamma, 0, z_sym=None)
        assert nu_K == want, (tuning, nu_K, want)


def test_full_gksl_deep_tail_slopes():
    """Full Lindblad transfer efficiency scales as kappa^-2 (generic) and
    kappa^-4 (protected) in the deep-kappa tail."""
    kappas = np.logspace(5, 9, 25)
    want = {"generic": 2.0, "protected": 4.0}
    for tuning, w in want.items():
        Kf = np.array([m.full_transfer_amplitude(k, tuning=tuning) for k in kappas])
        eff = core.fit_nu_loglog(kappas, np.abs(Kf) ** 2)["nu_global"]
        assert abs(eff - w) < 0.05, (tuning, eff, w)


def test_reduced_equals_full():
    """The reduced amplitude pencil equals the full-Liouvillian coherence
    sector up to the trivial -i drive factor (~1e-9)."""
    for tuning in ("generic", "protected"):
        for k in (1e3, 1e5, 1e7):
            Kred = m.transfer_kernel(k, tuning=tuning)
            Kfull = m.full_transfer_amplitude(k, tuning=tuning) / (-1j)
            assert abs(Kfull - Kred) / abs(Kred) < 1e-7, (tuning, k)


def test_ground_state_is_steady():
    """|g><g| is the drive-off steady state of the full Lindbladian."""
    assert m.verify_rho0_steady(1e3) < 1e-12


def test_symmetry_breaking_crossover_moves_with_eps():
    """Breaking the interference (eps != 0) restores nu_K = 1 at large kappa,
    with the crossover scale kappa*(eps) growing as eps -> 0."""
    kappas = np.logspace(4, 12, 90)
    kstar = {}
    for eps in (1e-9, 1e-10, 1e-11):
        vals = np.array([m.transfer_kernel(k, tuning="broken", eps=eps) for k in kappas])
        fit = core.fit_nu_loglog(kappas, vals)
        below = np.where(fit["nu_eff"] < 1.5)[0]
        kstar[eps] = fit["gamma_mid"][below[0]] if len(below) else np.inf
    assert kstar[1e-9] < kstar[1e-10] < kstar[1e-11]
    # protected (eps = 0) stays at nu = 2
    vals0 = np.array([m.transfer_kernel(k, tuning="protected") for k in kappas])
    assert abs(core.fit_nu_loglog(kappas, vals0)["nu_eff"][-1] - 2.0) < 0.05


def test_asymptotic_window_is_not_engineerable():
    """The kappa range where the integer order is read is far above the range a
    bus can actually be tuned over. Pinning this keeps the two from being
    conflated: the number of decades a numerical sweep covers is not evidence
    of experimental reach."""
    assert m.KAPPA_ASYMPTOTIC[0] > m.KAPPA_ENGINEERABLE[1] * 1e5
    assert m.decades(m.KAPPA_ENGINEERABLE) < 3.0


def test_physical_window_order_is_not_the_asymptotic_integer():
    """Inside the engineerable range kappa sits far below the couplings, so the
    response has not entered its 1/kappa tail and the measured efficiency order
    is ~0 for both tunings -- nowhere near the asymptotic 2 and 4.

    This is the substance behind Gate B's
    STRUCTURAL_PASS_EXPERIMENTAL_WINDOW_UNRESOLVED verdict; if a future change
    made the window reachable this test should be updated deliberately, not
    silently."""
    lo, hi = m.KAPPA_ENGINEERABLE
    kappas = np.logspace(np.log10(lo), np.log10(hi), 120)
    for tuning, asymptotic in (("generic", 2.0), ("protected", 4.0)):
        K = np.array([m.transfer_kernel(k, tuning=tuning) for k in kappas])
        eff_order = core.fit_nu_loglog(kappas, np.abs(K) ** 2)["nu_global"]
        assert abs(eff_order) < 0.1, (tuning, eff_order)
        assert abs(eff_order - asymptotic) > 0.5, (tuning, eff_order)


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
