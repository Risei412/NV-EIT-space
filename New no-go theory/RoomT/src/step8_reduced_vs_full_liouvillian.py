"""Step 8 of the room-temperature SMRT no-go plan (plan Sec. 5, Step 8):
agreement between Model R (nv_model.py's reduced 6-level excited-manifold
kernel, used for Steps 1/3/4/5/7's asymptotic/optimization work) and
Model L (the full N=9 Lindblad master equation, gate2_candidate_full_vs_
reduced.py's build_full/full_spectrum, used in Steps 2/6 for the
operational-cut audit and dip discrimination), across the SAME reference
geometry, from the T=10 K positive control through the sign crossover to
T=300 K.

Both models are evaluated at the SAME, unambiguous point (Model R at
z=Ep, Model L at d2=0) rather than via a "spectral peak search" -- an
exploratory version of this comparison tried a window-max peak search
and got a spurious result (T=30 K appeared to need a +/-50 MHz window to
find a "non-edge" extremum, and that extremum did not correspond to the
same physical feature as the near-resonance dip): the peak-search
ambiguity that Step 5 and Step 6 already ran into is even worse in the
30-70 K crossover region, where multiple spectral features of comparable
size coexist. Evaluating at a single, well-defined point (z=Ep / d2=0,
matching every other Step in this campaign) sidesteps that ambiguity
entirely and is the fair, reproducible comparison.

Key finding: Model R systematically OVERESTIMATES |C| relative to Model
L, by a factor that shrinks monotonically with T (~80x at 10 K, ~3x at
70 K, ~1.00-1.04x for T>=100 K) -- and the size of the overestimate is
self-diagnosed by Model R's OWN reported |C|: Model R's adiabatic-
elimination formula (dXi = -beta*K12*K21/den) assumes the Raman
correction is a SMALL perturbation to the bare susceptibility (|C|<<1);
at T=10 K, Model R's own |C|~0.87 is manifestly NOT small, signaling its
own breakdown, while at T>=100 K, |C|~1e-6 or smaller confirms the
perturbative assumption holds and agreement with Model L is excellent.

Consequence for the rest of this campaign: Step 1's low-T contrast VALUES
(e.g. 0.89 at 4 K) were only ever used qualitatively there (EIT works,
cut kills it, trends are sensible) and that conclusion is unaffected;
Steps 3/4's exact symbolic Gamma-scaling analysis is derived in and
applies to the LARGE-Gamma regime, which is exactly where this step shows
Model R and Model L agree to a few percent -- so the exponent/mechanism
result is validated on the model that actually matters (Model L) too, not
merely self-consistent within Model R. Steps 5/6/7's 300 K conclusions
were already close to model-independent (both models agree there to ~4%).

Outputs: RoomT/results/gates_summary_step8.json,
         RoomT/results/figures/fig_step8_model_agreement.png
"""
from __future__ import annotations
import os
import sys
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOMT_ROOT = os.path.dirname(HERE)
NOGO_SRC = os.path.join(ROOMT_ROOT, "..", "..", "No-go theorem", "src")
sys.path.insert(0, NOGO_SRC)

import nv_model as nm  # noqa: E402
import gate2_candidate_full_vs_reduced as g2  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

PROBE_ORB = np.array([1.0, 0.0], complex)
CTRL_ORB = np.array([0.0, 1.0], complex)
BZ = 0.02
D_STRAIN = 1.683
CTRL_BRANCH = "-1"
REFERENCE_OC = 1.0
REFERENCE_GG = 6.3e-5

T_GRID = [10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 300.0]


def model_R_and_L_at(T):
    Bvec = (0.0, 0.0, BZ)
    w, U = nm.dressed_ground(Bvec)
    H = nm.Hes(Bvec, d=D_STRAIN)
    dp, dc = nm.legs(U, ctrl=CTRL_BRANCH)
    Ep, fid, _ = nm.probe_line(H, dp)
    gamma = nm.gamma_oc_GHz(T, D_STRAIN)
    r = nm.response(H, dp, dc, Ep, gamma, Oc=REFERENCE_OC, gg=REFERENCE_GG)
    C_R = float(r["C"])

    sp = g2.full_spectrum(np.array([0.0]), T=T, Bx=0.0, Bz=BZ, isc=False,
                           Oc=REFERENCE_OC, ctrl=CTRL_BRANCH, ppol=PROBE_ORB,
                           cpol=CTRL_ORB, d=D_STRAIN)
    C_L = float(sp["C"][0])
    ratio = C_R / C_L if abs(C_L) > 1e-300 else float("nan")
    return dict(T_K=T, Gamma_GHz=float(gamma), C_modelR=C_R, C_modelL=C_L,
                ratio_R_over_L=ratio, sign_agree=bool(np.sign(C_R) == np.sign(C_L)))


def spot_check_step5_points():
    """Re-evaluate a couple of Step 5's searched adversarial points with
    Model L directly, since Step 5's optimization used Model R
    exclusively -- confirming the "no positive contrast" conclusion
    doesn't depend on which model is trusted, at 300 K where both models
    already agree closely."""
    points = [
        dict(label="Bx=0 reference (Step 1-7 baseline)", Bx=0.0, Bz=0.02, d=1.683, Oc=1.0),
        dict(label="Bx=0.2T (Step 5's rejected E1-failing candidate)", Bx=0.2, Bz=0.0, d=0.0, Oc=1.0),
    ]
    results = []
    for p in points:
        sp = g2.full_spectrum(np.array([0.0]), T=300.0, Bx=p["Bx"], Bz=p["Bz"], isc=False,
                               Oc=p["Oc"], ctrl=CTRL_BRANCH, ppol=PROBE_ORB, cpol=CTRL_ORB, d=p["d"])
        results.append(dict(label=p["label"], params=p, C_modelL=float(sp["C"][0])))
    return results


def run():
    rows = [model_R_and_L_at(T) for T in T_GRID]

    high_T_rows = [r for r in rows if r["T_K"] >= 100.0]
    max_high_T_rel_err = max(abs(r["ratio_R_over_L"] - 1.0) for r in high_T_rows)

    # self-diagnostic: |C_R| itself predicts the size of the R/L discrepancy
    log_CR = np.log(np.abs([r["C_modelR"] for r in rows]))
    log_dev = np.log(np.abs([r["ratio_R_over_L"] - 1.0 for r in rows]) + 1e-3)
    corr = float(np.corrcoef(log_CR, log_dev)[0, 1])

    spot_checks = spot_check_step5_points()

    gates = dict(
        high_T_agreement_within_10pct=bool(max_high_T_rel_err < 0.10),
        low_T_discrepancy_explained_by_perturbative_breakdown=bool(corr > 0.8),
        model_L_confirms_negative_or_negligible_at_300K_spotchecks=bool(
            all(c["C_modelL"] <= 1e-6 for c in spot_checks)
        ),
        sign_agreement_from_100K_up=bool(
            all(r["sign_agree"] for r in rows if r["T_K"] >= 100.0)
        ),
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        model="No-go theorem/src/nv_model.py (Model R) vs "
              "gate2_candidate_full_vs_reduced.py (Model L), Step 1-7 reference geometry",
        comparison_rows=rows,
        max_high_T_relative_error=max_high_T_rel_err,
        self_diagnostic_correlation_logCR_vs_logdeviation=corr,
        step5_spot_checks_with_model_L=spot_checks,
        gates=gates,
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step8.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    Ts = [r["T_K"] for r in rows]
    CR = [abs(r["C_modelR"]) for r in rows]
    CL = [abs(r["C_modelL"]) for r in rows]
    ax.loglog(Ts, CR, "o-", label="Model R (reduced kernel)")
    ax.loglog(Ts, CL, "s-", label="Model L (full Liouvillian)")
    ax.set_xlabel("T (K)"); ax.set_ylabel("|C|")
    ax.set_title("Model R vs Model L: |C(T)|"); ax.legend(fontsize=8)

    ax = axes[1]
    ratios = [r["ratio_R_over_L"] for r in rows]
    ax.semilogx(Ts, ratios, "o-", color="tab:purple")
    ax.axhline(1.0, color="gray", linestyle="--", label="perfect agreement")
    ax.set_xlabel("T (K)"); ax.set_ylabel("C_modelR / C_modelL")
    ax.set_title("Model R / Model L ratio vs T")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step8_model_agreement.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps(gates, indent=2))
    for r in rows:
        print(f"T={r['T_K']:6.1f}K  C_R={r['C_modelR']:+.4e}  C_L={r['C_modelL']:+.4e}  "
              f"ratio={r['ratio_R_over_L']:+.3f}")
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
