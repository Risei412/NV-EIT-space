"""Step 1 of the room-temperature SMRT no-go plan
(SMRT_NV_room_temperature_EIT_no_go_numerical_plan.md, Sec. 5, Step 1):
"positive control" -- show the model and code reproduce the KNOWN
low-temperature NV branch-resolved spin-Lambda EIT before any room-
temperature no-go claim is attempted (plan: "Gate 1をを通過しないcodeで
室温no-goを論じてはならない").

Model: reuses `No-go theorem/src/nv_model.py` verbatim (no new physics) --
this is exactly the ordinary NV ZPL spin-Lambda in question: ground states
|ms=0> (probe leg) and |ms=-1> (control leg), connected through the two
orbital branches (X, Y) of the thermally-split 3E manifold. gamma_oc_GHz(T,d)
is the phonon-driven X<->Y orbital-branch hopping rate (Happacher et al.,
already validated in phonon_rates.py) that grows with T and washes out the
branch distinction the Lambda needs -- exactly the SMRT "thermal manifold
merging" mechanism, at low T where it is still small.

Cut used here (Step 1 only needs the coarse full/no-control comparison;
the operational, GKSL-admissible sector cut with its O(kappa^-1)
convergence check is Step 2, run_step2_operational_cut_audit.py):
setting the control Rabi frequency Oc=0 kills K12*K21's contribution to
the response (dXi ~ Oc^2), which is the same "sever the coupling that
creates the ground Raman coherence" cut already validated as equivalent
to the operational algebraic cut for this Lambda block (see
model_lambda.py docstring, and Sector/model_lambda_operational.py's
D_S = diag(0,1) construction, which is exactly this block).

Outputs: RoomT/results/gates_summary_step1.json,
         RoomT/results/figures/fig_step1_low_T_positive_control.png
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

PHASE_A_SRC = os.path.join(ROOMT_ROOT, "..", "src")
sys.path.insert(0, PHASE_A_SRC)
import model_lambda as ml  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Candidate Lambda geometry: small Bz to make the (ms=-1, ms=+1) branch
# assignment in dressed_ground() unambiguous (matches No-go
# theorem/tests/test_core.py's B tuples); strain d fixed at the value used
# throughout the existing gate scripts (gate2_candidate_full_vs_reduced.py).
BVEC = (0.0, 0.0, 0.02)
D_STRAIN = 1.683

LOW_T_GRID = [4.0, 10.0, 20.0, 30.0]
OC_GRID = [0.2, 0.5, 1.0, 2.0, 4.0, 8.0]
GG_GRID = [1e-6, 6.3e-6, 6.3e-5, 6.3e-4, 6.3e-3]

REFERENCE_T = 10.0
REFERENCE_OC = 1.0
REFERENCE_GG = 6.3e-5


def _build_model():
    w, U = nm.dressed_ground(BVEC)
    H = nm.Hes(BVEC, d=D_STRAIN)
    dp, dc = nm.legs(U)
    Ep, fid, _ = nm.probe_line(H, dp)
    return H, dp, dc, Ep


def spectrum(H, dp, dc, gamma, Oc, gg, z_near=2.0, n=401, cut=False):
    """Full and sector-cut (Oc set to 0 in the response formula) spectra on
    a window around the excited-state resonance Ep. Returns
    (z_grid, Acut_bare, A_full, A_cut_response) where Acut_bare = Im(S1)
    is the model's own control-independent bare-absorption reference and
    A_cut_response is chi_cut^(EIT) = Acut_bare (since dXi = 0 identically
    at Oc = 0) -- kept as a separate array so this Step-1 comparison and
    the Step-2 operational-cut audit report the same quantity under the
    same name."""
    zz, Ep = nm.window(H, dp, "A", z_near=z_near, n=n)
    Oc_eff = 0.0 if cut else Oc
    rows = [nm.response(H, dp, dc, z, gamma, Oc=Oc_eff, gg=gg) for z in zz]
    Acut = np.array([r["Acut"] for r in rows])
    dA = np.array([r["C"] for r in rows]) * Acut
    Afull = Acut - dA
    return zz, Acut, Afull, Ep


def fwhm_of_transparency(zz, Cvals):
    """FWHM of the EIT contrast feature |C(z)| (Step 1's operational proxy
    for 'EIT linewidth increases with control power', plan Sec. 5 Step 1
    required results)."""
    absC = np.abs(Cvals)
    peak = absC.max()
    if peak < 1e-300:
        return 0.0
    half = 0.5 * peak
    above = absC >= half
    if not np.any(above):
        return 0.0
    idx = np.where(above)[0]
    return float(zz[idx[-1]] - zz[idx[0]])


def run():
    H, dp, dc, Ep = _build_model()

    # --- (a) low-T positive control: nonzero R_EIT, killed by the cut ---
    T_results = []
    for T in LOW_T_GRID:
        gamma = nm.gamma_oc_GHz(T, D_STRAIN)
        r_full = nm.response(H, dp, dc, Ep, gamma, Oc=REFERENCE_OC, gg=REFERENCE_GG)
        r_cut = nm.response(H, dp, dc, Ep, gamma, Oc=0.0, gg=REFERENCE_GG)
        T_results.append(dict(T=T, gamma_oc_GHz=gamma, C_full=r_full["C"], C_cut=r_cut["C"]))

    C0_by_T = [abs(row["C_full"]) for row in T_results]
    cut_kills_feature = all(abs(row["C_cut"]) < 1e-12 for row in T_results)
    nonzero_at_low_T = C0_by_T[0] > 1e-9
    # allow small numerical non-monotonicity but require the overall trend
    # to fall as T rises (qualitative reproduction of the known suppression
    # of the branch-resolved Lambda with increasing phonon-driven mixing)
    monotone_decreasing_in_T = all(
        C0_by_T[i + 1] <= C0_by_T[i] * 1.05 for i in range(len(C0_by_T) - 1)
    )

    # --- (b) EIT linewidth grows with control power (fixed low T) ---
    # nv_model.response() ties probe and control to a single shared
    # detuning z (it evaluates the Raman/two-photon-resonant kernel at one
    # frequency, per its archived-test usage in No-go theorem/tests/
    # test_core.py -- there is no independent two-photon-detuning axis to
    # scan). The standard "EIT linewidth broadens with control power"
    # signature needs that independent axis, so this gate uses the clean
    # three-level Lambda model (model_lambda.py, Phase A positive control,
    # plan Sec. 3.2 item 1) instead, with gamma31 set to the same physical
    # NV optical/orbital linewidth used above (gamma_oc_GHz at REFERENCE_T)
    # so the two models share the same dissipation scale.
    gamma_ref = nm.gamma_oc_GHz(REFERENCE_T, D_STRAIN)
    widths = []
    window_edge_clipped = []
    for Oc in OC_GRID:
        # window scales with the expected power-broadened width
        # (~ Oc^2/gamma_ref) plus margin, so the FWHM measurement is never
        # clipped by a fixed window (that would silently plateau the fit).
        half_span = max(2.0, 3.0 * Oc ** 2 / gamma_ref)
        delta_grid = np.linspace(-half_span, half_span, 2001)
        chi_full = np.array([
            ml.chi_full_analytic(0.0, d, Oc, gamma_ref, REFERENCE_GG) for d in delta_grid
        ])
        # absorption quadrature of chi = 1/(gamma - i*Delta) is Re(chi)
        # (peaked, Lorentzian at Delta=0), not Im(chi) (dispersive, zero at
        # Delta=0) -- matches nv_model.response()'s Acut = Re(S1) convention.
        chi_cut = ml.chi_cut_analytic(0.0, 0.0, Oc, gamma_ref, REFERENCE_GG)
        contrast = (chi_cut.real - chi_full.real) / chi_cut.real
        w = fwhm_of_transparency(delta_grid, contrast)
        widths.append(w)
        window_edge_clipped.append(bool(w > 1.8 * half_span))
    width_increasing_with_Oc = all(
        widths[i + 1] >= widths[i] * 0.98 for i in range(len(widths) - 1)
    ) and not any(window_edge_clipped)

    # --- (c) contrast decreases with ground-state decoherence gg ---
    C_by_gg = []
    for gg in GG_GRID:
        r = nm.response(H, dp, dc, Ep, gamma_ref, Oc=REFERENCE_OC, gg=gg)
        C_by_gg.append(abs(r["C"]))
    monotone_decreasing_in_gg = all(
        C_by_gg[i + 1] <= C_by_gg[i] * 1.05 for i in range(len(C_by_gg) - 1)
    )

    gates = dict(
        nonzero_at_low_T=bool(nonzero_at_low_T),
        cut_kills_feature=bool(cut_kills_feature),
        monotone_decreasing_in_T=bool(monotone_decreasing_in_T),
        width_increasing_with_control_power=bool(width_increasing_with_Oc),
        no_window_edge_clipping=bool(not any(window_edge_clipped)),
        monotone_decreasing_in_ground_decoherence=bool(monotone_decreasing_in_gg),
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        model="No-go theorem/src/nv_model.py (branch-resolved NV spin-Lambda)",
        Bvec_T=BVEC,
        d_strain_GHz=D_STRAIN,
        Ep_GHz=Ep,
        low_T_grid_K=LOW_T_GRID,
        low_T_results=T_results,
        Oc_grid=OC_GRID,
        eit_fwhm_GHz_vs_Oc=widths,
        window_edge_clipped=window_edge_clipped,
        gg_grid_GHz=GG_GRID,
        C_vs_gg=C_by_gg,
        gates=gates,
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step1.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure ---
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    zz, Acut, Afull, _ = spectrum(H, dp, dc, gamma_ref, REFERENCE_OC, REFERENCE_GG)
    _, Acut_cut, Afull_cut, _ = spectrum(H, dp, dc, gamma_ref, REFERENCE_OC, REFERENCE_GG, cut=True)
    ax = axes[0, 0]
    ax.plot(zz - Ep, Acut, "k--", label="bare (chi_cut)")
    ax.plot(zz - Ep, Afull, "r-", label="full (control on)")
    ax.plot(zz - Ep, Afull_cut, "b:", label="sector cut (Oc=0)")
    ax.set_xlabel("probe detuning z - Ep (GHz)"); ax.set_ylabel("Im chi (a.u.)")
    ax.set_title(f"T={REFERENCE_T} K spectrum"); ax.legend(fontsize=7)

    ax = axes[0, 1]
    ax.semilogy(LOW_T_GRID, C0_by_T, "o-")
    ax.set_xlabel("T (K)"); ax.set_ylabel("|C_EIT|"); ax.set_title("Contrast vs low T")

    ax = axes[1, 0]
    ax.plot(OC_GRID, widths, "s-")
    ax.set_xlabel("Omega_c (GHz)"); ax.set_ylabel("EIT FWHM (GHz)")
    ax.set_title("Linewidth vs control power")

    ax = axes[1, 1]
    ax.loglog(GG_GRID, C_by_gg, "^-")
    ax.set_xlabel("gamma_ground (GHz)"); ax.set_ylabel("|C_EIT|")
    ax.set_title("Contrast vs ground decoherence")

    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step1_low_T_positive_control.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps(gates, indent=2))
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
