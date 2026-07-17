"""Gate U8 (information-metric naming audit): checks whether
Delta F_Q = F_Q[rho_full] - F_Q[rho_cut] and I_sector_metric =
g_rho_bar(x_S, x_S) are the same quantity, justifying (or not) the rename
away from "F_{Q,S}" recommended in Sec. 8 of the revision note and Sec.
4.3 of step4_amendments.md.

Two parts:
  (1) A synthetic qubit counterexample with hand-chosen drho_full,
      drho_cut proving Delta F_Q != I_sector_metric *in general* (a
      nonzero cross term is exhibited explicitly, in closed form).
  (2) An empirical check on the repository's physical 3-level Lindblad
      model (`New no-go theory/src/model_metro_lindblad.py`, unmodified),
      reporting the actual cross term found there -- which turns out to
      be (numerically) negligible at the tested point, a genuine feature
      of that particular kappa_cut=0 construction (drho_cut vanishes
      identically there since a purely diagonal cut Hamiltonian carries
      no theta-dependence into the steady state), consistent with the
      near-unit ratio nu_F/(2 nu_x) already reported for Gate M3 in the
      top-level README. This is not evidence against the naming
      distinction -- part (1) already settles that -- it documents why
      the distinction happens to be numerically invisible in that one
      existing model.
"""

import json
import os
import sys
import importlib.util

import numpy as np

_MODEL_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "src", "model_metro_lindblad.py"
)
_spec = importlib.util.spec_from_file_location("model_metro_lindblad", _MODEL_PATH)
mml = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mml)

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def synthetic_qubit_counterexample():
    """Hand-built qubit example: rho_bar = diag(0.7, 0.3) (a faithful
    reference state), drho_full = [[0,a],[a*,0]] with a=0.3,
    drho_cut = [[0,b],[b*,0]] with b=0.1j. Both off-diagonal, so the SLD
    QFI formula collapses to F(rho,drho) = 4|drho_01|^2 (since
    p0+p1=1 identically for a qubit)."""
    p0, p1 = 0.7, 0.3
    a = 0.3 + 0.0j
    b = 0.0 + 0.1j
    x = a - b
    F_full = 4 * abs(a) ** 2
    F_cut = 4 * abs(b) ** 2
    delta_F_Q = F_full - F_cut
    I_sector_metric = 4 * abs(x) ** 2
    # cross_term as defined by qfi_cross_term's convention (2 Re<x,G drho_cut>
    # summed over both off-diagonal index pairs); the correct bilinear-form
    # identity is F_full = I_sector_metric + F_cut + 2*cross_term, i.e.
    #   Delta_F_Q = F_full - F_cut = I_sector_metric + 2*cross_term
    # (verified numerically against `qfi`/`qfi_cross_term` directly).
    cross_term = 4 * (x * np.conj(b)).real
    reconstruction_residual = delta_F_Q - I_sector_metric - 2 * cross_term
    return {
        "rho_bar_populations": [p0, p1],
        "a_drho_full_01": str(a),
        "b_drho_cut_01": str(b),
        "F_Q_full": F_full,
        "F_Q_cut": F_cut,
        "Delta_F_Q": delta_F_Q,
        "I_sector_metric": I_sector_metric,
        "cross_term": cross_term,
        "reconstruction_residual_should_be_zero": reconstruction_residual,
        "quantities_are_distinct": abs(cross_term) > 1e-9,
    }


def physical_model_check(Gamma=50.0, theta=0.3, lam=0.7, phi=0.4):
    x_S, drho_full, drho_cut, rho_full, rho_cut = mml.x_S_lindblad(
        Gamma, theta, lam, phi=phi
    )
    x_S_mat = mml.unvec(x_S)

    F_full, dropped_full, pmin_full = mml.qfi(rho_full, mml.unvec(drho_full))
    F_cut, dropped_cut, pmin_cut = mml.qfi(rho_cut, mml.unvec(drho_cut))
    delta_F_Q = F_full - F_cut

    # I_sector_metric = g_{rho_bar}(x_S, x_S), rho_bar = rho_full (a valid,
    # faithful common reference state per the Assumption in Sec. 8).
    I_sector_metric, dropped_I, pmin_I = mml.qfi(rho_full, x_S_mat)

    cross_term = mml.qfi_cross_term(rho_full, x_S_mat, mml.unvec(drho_cut))

    # See synthetic_qubit_counterexample for the verified factor of 2; this
    # residual is only exactly zero when rho_full and rho_cut share the same
    # metric (approximately true here since drho_cut = 0 identically).
    reconstruction_residual = delta_F_Q - I_sector_metric - 2 * cross_term

    return {
        "params": {"Gamma": Gamma, "theta": theta, "lam": lam, "phi": phi},
        "F_Q_full": F_full,
        "F_Q_cut": F_cut,
        "Delta_F_Q": delta_F_Q,
        "I_sector_metric": I_sector_metric,
        "cross_term": cross_term,
        "reconstruction_residual_should_be_zero": reconstruction_residual,
        "Delta_F_Q_minus_I_sector_metric": delta_F_Q - I_sector_metric,
        "note": (
            "cross_term ~ 0 here because kappa_cut=0 makes the cut "
            "Hamiltonian purely diagonal, so drho_cut vanishes identically "
            "and x_S = drho_full exactly; see module docstring part (2)."
        ),
    }


if __name__ == "__main__":
    synthetic = synthetic_qubit_counterexample()
    physical = physical_model_check()

    print("[U8-part1] Synthetic qubit counterexample "
          "(proves Delta_F_Q != I_sector_metric in general):")
    for k, v in synthetic.items():
        print(f"    {k}: {v}")
    part1_pass = (
        abs(synthetic["reconstruction_residual_should_be_zero"]) < 1e-9
        and synthetic["quantities_are_distinct"]
    )
    print(f"    -> {'PASS' if part1_pass else 'FAIL'}")

    print("\n[U8-part2] Physical 3-level Lindblad model (informative, not "
          "a pass/fail check):")
    for k, v in physical.items():
        print(f"    {k}: {v}")

    result = {"part1_synthetic_counterexample": synthetic,
              "part1_pass": bool(part1_pass),
              "part2_physical_model": physical,
              "pass": bool(part1_pass)}

    out_path = os.path.join(RESULTS_DIR, "gate_u8_naming_audit.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\nWritten: {out_path}")
