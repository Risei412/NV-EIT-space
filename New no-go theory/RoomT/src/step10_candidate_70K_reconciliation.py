"""Step 10: reconcile the 70 K PRL candidate point with the room-temperature
campaign's Criterion-E1 verdict.

Open question raised by the review draft
(`Writing Paper/drafts/NV_EIT_nogo_review_ja.md`, Sec. 4.5 item 1):

  * `No-go theorem/`'s PRL candidate sits at T=70 K, Bx=BX0=0.2323 T and
    reports an ensemble contrast C ~ 1.4e-2 (Gate E anchor).
  * This campaign's Step 5 found that at T=300 K, Bx=0.2 T fails Criterion
    E1 (C < 0: the control field INCREASES absorption), and Step 5's
    reference-configuration sign crossover is bracketed at [50, 77] K.

Read naively, those two look contradictory: is the 70 K candidate on the
transparency side of the crossover or not?

They are NOT the same configuration. The two geometries differ in every
selection-rule-relevant way:

                          RoomT Steps 1-9        PRL candidate
    probe polarization    orbital X              orbital Y
    control polarization  orbital Y              orbital Y   (SAME branch)
    control spin branch   ms = -1                ms = +1
    Bz                    0.02 T                 0.005 T
    Bx                    0 (0.2 in Step 5's     0.2323 T
                          adversarial sweep)
    Omega_c               1.0 GHz                0.1 GHz
    readout point         two-photon resonance   dressed eigenvalue j=3

So Step 5's rejected "Bx=0.2 T" point is the RoomT geometry evaluated at a
large transverse field -- not the candidate. Nothing in Steps 1-9 has ever
evaluated the candidate geometry, and the sign crossover bracket [50,77] K
is a statement about the RoomT reference configuration only.

This step closes the gap by evaluating the CANDIDATE geometry with the same
full N=9 Lindblad model (Model L) and the same Criteria E1/E2/E4 that Step 6
applied to the RoomT geometry, and by mapping the candidate's own C(T) sign
crossover. `gate2_candidate_full_vs_reduced.py`'s defaults ARE the candidate
geometry, so no new model construction is needed.

Gates:
  G1 candidate_passes_E1_at_70K   -- C > 0 somewhere in the scanned window
  G2 candidate_passes_E2          -- C == 0 exactly with the control off
  G3 candidate_E4_quadratic       -- C ~ Omega_c^2 at small Omega_c
  G4 candidate_sign_crossover_above_70K -- the candidate's own crossover
                                     temperature is above 70 K
  G5 roomT_geometry_differs       -- the two geometries are recorded as
                                     distinct (documentation gate)

Outputs: RoomT/results/gates_summary_step10.json
"""
from __future__ import annotations
import os
import sys
import json

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOMT_ROOT = os.path.dirname(HERE)
NOGO_SRC = os.path.join(ROOMT_ROOT, "..", "..", "No-go theorem", "src")
sys.path.insert(0, NOGO_SRC)

import gate2_candidate_full_vs_reduced as g2  # noqa: E402
import run_prl_prediction as rp  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Candidate geometry (= gate2/run_prl_prediction defaults, restated explicitly
# so this file is self-documenting rather than relying on imported constants).
CAND = dict(Bx=rp.BX0, Bz=rp.BZ0, ctrl="+1", ppol=rp.Y, cpol=rp.Y, d=rp.D)
CAND_OC = rp.OC          # 0.1 GHz
CAND_T = 70.0
HALF_GHZ = 0.005         # +/-5 MHz, the window gate2 itself uses for this dip
N_GRID = 121


def spectrum(T, Oc, half_GHz=HALF_GHZ, n=N_GRID):
    d2s = np.linspace(-half_GHz, half_GHz, n)
    return g2.full_spectrum(d2s, T=T, Oc=Oc, **CAND)


def peak_C(T, Oc, half_GHz=HALF_GHZ, n=N_GRID):
    """Largest correctly-signed (transparency) contrast in the window, and
    the largest-|C| feature for reference."""
    sp = spectrum(T, Oc, half_GHz, n)
    C = np.asarray(sp["C"], float)
    i_pos = int(np.argmax(C))
    i_abs = int(np.argmax(np.abs(C)))
    return dict(C_max_signed=float(C[i_pos]),
                d2_at_max_signed_MHz=float(sp["d2_MHz"][i_pos]),
                C_at_max_abs=float(C[i_abs]),
                d2_at_max_abs_MHz=float(sp["d2_MHz"][i_abs]),
                on_window_edge=bool(i_pos in (0, n - 1)),
                A_cut_at_max=float(sp["A_cut"][i_pos]),
                rho_cp_over_Op=float(sp["diag"]["rho_cp_over_Op"]))


def main():
    out = {
        "purpose": "reconcile the 70 K PRL candidate with Step 5/6's E1 verdict",
        "candidate_geometry": {
            "Bx_T": float(rp.BX0), "Bz_T": float(rp.BZ0),
            "control_spin_branch": "+1",
            "probe_polarization": "orbital Y", "control_polarization": "orbital Y",
            "same_orbital_branch": True,
            "d_strain_GHz": float(rp.D), "Oc_GHz": float(CAND_OC),
        },
        "roomT_geometry": {
            "Bx_T": 0.0, "Bz_T": 0.02,
            "control_spin_branch": "-1",
            "probe_polarization": "orbital X", "control_polarization": "orbital Y",
            "same_orbital_branch": False,
            "d_strain_GHz": 1.683, "Oc_GHz": 1.0,
        },
        "model": "gate2_candidate_full_vs_reduced.full_spectrum (full N=9 Lindblad)",
    }

    # --- E1 at the candidate point -------------------------------------
    e1 = peak_C(CAND_T, CAND_OC)
    out["E1_at_70K"] = e1

    # --- E2: control off ------------------------------------------------
    off = spectrum(CAND_T, 0.0)
    out["E2_control_off"] = {
        "max_abs_C": float(np.max(np.abs(np.asarray(off["C"], float)))),
    }

    # --- E4: quadratic emergence at the E1 peak location ----------------
    # The candidate's own Omega_c = 0.1 GHz already sits in the power-broadened
    # saturation shoulder, so a fit spanning it underestimates the exponent.
    # Fit only the unsaturated decades (|C| below 1% of the reference value).
    d2_star = e1["d2_at_max_signed_MHz"] * 1e-3
    Ocs = np.array([1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
    Cs = []
    for Oc in Ocs:
        sp = g2.full_spectrum(np.array([d2_star]), T=CAND_T, Oc=float(Oc), **CAND)
        Cs.append(float(sp["C"][0]))
    Cs = np.asarray(Cs)
    small = (np.abs(Cs) > 0) & (np.abs(Cs) < 0.01 * abs(e1["C_max_signed"]))
    slope = float(np.polyfit(np.log10(Ocs[small]), np.log10(np.abs(Cs[small])), 1)[0]) \
        if small.sum() >= 2 else float("nan")
    slope_all = float(np.polyfit(np.log10(Ocs[np.abs(Cs) > 0]),
                                 np.log10(np.abs(Cs[np.abs(Cs) > 0])), 1)[0])
    out["E4_control_power_scaling"] = {
        "d2_star_MHz": float(d2_star * 1e3),
        "Oc_GHz": Ocs.tolist(), "C": Cs.tolist(),
        "loglog_slope_unsaturated": slope,
        "n_points_unsaturated": int(small.sum()),
        "loglog_slope_all_points_includes_saturation": slope_all,
    }

    # --- the candidate's OWN temperature dependence ---------------------
    Ts = [4.0, 10.0, 30.0, 50.0, 70.0, 77.0, 100.0, 150.0, 300.0]
    rows = []
    for T in Ts:
        r = peak_C(T, CAND_OC)
        r["T_K"] = float(T)
        rows.append(r)
    out["candidate_temperature_sweep"] = rows

    Cs_T = [r["C_max_signed"] for r in rows]
    crossover = None
    for (Ta, Ca), (Tb, Cb) in zip(zip(Ts, Cs_T), list(zip(Ts, Cs_T))[1:]):
        if Ca > 0 and Cb <= 0:
            crossover = [float(Ta), float(Tb)]
            break
    out["candidate_sign_crossover_bracket_K"] = crossover

    # --- detectability ceiling for the candidate geometry ---------------
    # Step 5/7's illustrative detection threshold, reused unchanged so the two
    # campaigns are compared on one yardstick.
    eps = 1.5344920423293595e-07
    out["detection_threshold"] = eps
    bracket = None
    for (Ta, Ca), (Tb, Cb) in zip(zip(Ts, Cs_T), list(zip(Ts, Cs_T))[1:]):
        if Ca >= eps > Cb:
            bracket = [float(Ta), float(Tb)]
            break
    out["candidate_detectability_ceiling_bracket_K"] = bracket
    if bracket is not None:
        # log-linear interpolation in (T, log10 C) between the bracketing points
        ia = Ts.index(bracket[0])
        Ta, Tb = bracket
        Ca, Cb = Cs_T[ia], Cs_T[ia + 1]
        frac = (np.log10(Ca) - np.log10(eps)) / (np.log10(Ca) - np.log10(Cb))
        out["candidate_detectability_ceiling_K_interp"] = float(Ta + frac * (Tb - Ta))

    gates = {
        "candidate_passes_E1_at_70K": bool(e1["C_max_signed"] > 0
                                           and not e1["on_window_edge"]),
        "candidate_passes_E2": bool(out["E2_control_off"]["max_abs_C"] < 1e-10),
        "candidate_E4_quadratic": bool(abs(slope - 2.0) < 0.3),
        "candidate_sign_crossover_above_70K": bool(
            crossover is None or crossover[0] >= 70.0),
        "roomT_geometry_differs": True,
    }
    gates["overall_pass"] = all(gates.values())
    out["gates"] = gates

    path = os.path.join(RESULTS_DIR, "gates_summary_step10.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in
                      ("E1_at_70K", "E2_control_off", "E4_control_power_scaling",
                       "candidate_sign_crossover_bracket_K", "gates")}, indent=1))
    print("wrote", path)


if __name__ == "__main__":
    main()
