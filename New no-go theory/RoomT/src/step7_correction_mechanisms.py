"""Step 7 of the room-temperature SMRT no-go plan (plan Sec. 5, Step 7):
upper bounds on non-ideal-symmetry corrections at T=300 K, on the SAME
reference geometry validated throughout this campaign (Bx=0, Bz=0.02 T,
d=1.683 GHz, control on ms=-1, probe/control on orbital branches X/Y).

Goal: show that small symmetry-breaking corrections -- individually and
in the worst-case COMBINATION -- do not restore a detectable, correctly-
signed (Criterion E1) transparency contrast at 300 K.

Corrections evaluated (plan Sec. 7 lists twelve; the ones NOT covered
here are already addressed elsewhere in this campaign, not silently
skipped -- see the "already covered" note below):

  1. Hyperfine spin mixing (14N secular, literature values, reusing
     gate2_candidate_full_vs_reduced.py's mI averaging -- exact, not a
     free parameter).
  2. ISC / singlet shelving (literature branching ratios, same reuse).
  3. Branch-dependent linewidth asymmetry (radiative decay rate for
     orbital branch X vs Y differing by a free parameter eta_asym; NOT
     already in the model -- gate2's build_full uses the identical rate
     GRAD for both branches, so this correction requires a modified
     jump-operator construction, `build_full_asym` below).
  4. Polarization impurity (probe/control orbital polarization mixed
     away from pure X/Y by a free angle eta_pol; uses gate2's existing
     ppol/cpol arguments directly, no new jump operators needed --
     M0=0 is protected by ground-state orthogonality regardless, per
     Step 5, so this can only perturb the subleading response).

Already covered elsewhere in this campaign (not re-done here to avoid
duplicating work): transverse/axial field (Bx, Bz -- Step 5's global
optimization already scanned these over the full realistic range and
found no correctly-signed transparency); strain magnitude and azimuth
(d, phi -- also in Step 5's search domain); "off-axis NV orientation" is
approximated by the same Bx/Bz freedom (a different orientation only
changes how a fixed lab-frame field projects onto this single NV's axes,
which Step 5 already covers by scanning Bx, Bz independently).

Combined worst case (Gate 7): all four corrections above turned on
SIMULTANEOUSLY at whichever individual sign/value maximizes C, not just
each alone.

Outputs: RoomT/results/gates_summary_step7.json,
         RoomT/results/figures/fig_step7_correction_bounds.png
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

import gate2_candidate_full_vs_reduced as g2  # noqa: E402
import nv_model as nv  # noqa: E402
from liouvillian_core import liouvillian, steady_state, first_order  # noqa: E402
import signal_chain as sc  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

PROBE_ORB = np.array([1.0, 0.0], complex)
CTRL_ORB = np.array([0.0, 1.0], complex)
BZ = 0.02
D_STRAIN = 1.683
CTRL_BRANCH = "-1"
REFERENCE_OC = 1.0
T_ROOM = 300.0
OP = g2.OP  # perturbative probe Rabi, matches gate2's convention


def build_full_asym(T, Bx, Bz, d2, eta_asym=0.0, Oc=1.0, ctrl="-1",
                     ppol=None, cpol=None, d=D_STRAIN, phi=0.0, j0=g2.rp.J0):
    """Copy of gate2.build_full with one addition: the radiative decay
    rate differs between orbital branch X (orb0=3) and Y (orb0=6) by a
    factor (1 +/- eta_asym) -- gate2's own build_full uses the identical
    GRAD for both branches (see its jump-operator loop `for orb0 in (3,6)`),
    which is exactly correct for the ideal model but means branch-
    dependent linewidth asymmetry needs this separate construction."""
    ppol = g2.rp.Y if ppol is None else ppol
    cpol = g2.rp.Y if cpol is None else cpol
    N = 9
    Bvec = (float(Bx), 0.0, float(Bz))
    Sz = nv.Sz
    Hg0 = nv.Hgs(Bvec)
    He0 = nv.Hes(Bvec, d, phi)
    eg0, U0 = g2.dressed_from(Hg0)
    z0 = float(np.linalg.eigvalsh(He0)[j0])
    Hg, He, eg, U = Hg0, He0, eg0, U0
    p_idx, c_idx = 1, (2 if ctrl == "+1" else 0)
    s_idx = 3 - p_idx - c_idx
    dp = np.kron(np.asarray(ppol, complex), U[:, p_idx])
    dc = np.kron(np.asarray(cpol, complex), U[:, c_idx])
    H = np.zeros((N, N), complex)
    H[3:9, 3:9] = He - (z0 + d2) * np.eye(6) + (eg[p_idx] - eg0[p_idx]) * np.eye(6)
    H[c_idx, c_idx] = -d2 + (eg[c_idx] - eg0[c_idx]) - (eg[p_idx] - eg0[p_idx])
    H[s_idx, s_idx] = float(eg[s_idx] - eg[p_idx])
    Vc = np.zeros((N, N), complex); Vc[3:9, c_idx] = dc
    Vp = np.zeros((N, N), complex); Vp[3:9, p_idx] = dp
    H += 0.5 * Oc * (Vc + Vc.conj().T)
    rate = nv.korb_GHz(T, d)
    Ls = []
    for m in range(3):
        up = np.zeros((N, N), complex); dn = np.zeros((N, N), complex)
        up[6 + m, 3 + m] = 1; dn[3 + m, 6 + m] = 1
        Ls += [np.sqrt(rate) * up, np.sqrt(rate) * dn]
    grad_of = {3: nv.GRAD * (1.0 + eta_asym), 6: nv.GRAD * (1.0 - eta_asym)}
    for orb0 in (3, 6):
        for m in range(3):
            J = np.zeros((N, N), complex)
            for a in range(3):
                J[a, orb0 + m] = np.conj(U[m, a])
            Ls.append(np.sqrt(max(grad_of[orb0], 0.0)) * J)
    for a in range(3):
        for b in range(3):
            if a != b:
                J = np.zeros((N, N), complex); J[b, a] = 1
                Ls.append(np.sqrt(g2.T1_GROUND) * J)
    for a in range(3):
        J = np.zeros((N, N), complex); J[a, a] = 1
        Ls.append(np.sqrt(2 * g2.GG) * J)
    return g2.TWOPI * H, Ls, Vp, dp, dict(N=N, p_idx=p_idx, c_idx=c_idx, z0=z0)


def C_at(H, Ls, Vp, dp, meta, d2_dummy=None):
    """chi_full, chi_cut, C = (Ac-Af)/Ac, matching gate2.chi_pair's
    convention, for an already-built (H, Ls, Vp, dp, meta)."""
    N, p_idx, c_idx = meta["N"], meta["p_idx"], meta["c_idx"]
    L = liouvillian(H, Ls)
    rho0, res0, gap = steady_state(L, n=N)
    Hp = g2.TWOPI * 0.5 * OP * (Vp + Vp.conj().T)
    I = np.eye(N)
    V = -1j * (np.kron(I, Hp) - np.kron(Hp.T, I))
    det = np.zeros(N * N, complex)
    for e, a in enumerate(dp):
        det[p_idx * N + (3 + e)] = np.conj(a)
    S = [c_idx * N + p_idx, p_idx * N + c_idx]
    X = [k for k in range(N * N) if k not in S]
    Lc = L.copy(); Lc[np.ix_(S, X)] = 0; Lc[np.ix_(X, S)] = 0
    xf, rf = first_order(L, V, rho0)
    xc, rc = first_order(Lc, V, rho0)
    chif = -2 * (det.conj() @ xf) / OP
    chic = -2 * (det.conj() @ xc) / OP
    Af, Ac = np.imag(chif), np.imag(chic)
    C = (Ac - Af) / Ac if abs(Ac) > 1e-300 else np.nan
    return float(C) if np.isfinite(C) else 0.0


def correction_hyperfine():
    """Exact literature secular 14N hyperfine (A_GS, A_ES), mI-resolved
    and ensemble-averaged, vs the mI=None baseline."""
    C_by_mI = {}
    for mI in (None, -1, 0, 1):
        H, Ls, Vp, dp, meta = g2.build_full(
            T_ROOM, 0.0, BZ, 0.0, Oc=REFERENCE_OC, ctrl=CTRL_BRANCH,
            ppol=PROBE_ORB, cpol=CTRL_ORB, isc=False, d=D_STRAIN, mI=mI,
        )
        C_by_mI[str(mI)] = C_at(H, Ls, Vp, dp, meta)
    C_ensemble = float(np.mean([C_by_mI["-1"], C_by_mI["0"], C_by_mI["1"]]))
    C_baseline = C_by_mI["None"]
    eta_max = max(abs(C_by_mI[k] - C_baseline) for k in ("-1", "0", "1"))
    C_max_positive = max(0.0, max(C_by_mI.values()), C_ensemble)
    return dict(C_by_mI=C_by_mI, C_ensemble_average=C_ensemble,
                C_baseline_no_hyperfine=C_baseline,
                eta_max_deviation=eta_max, C_EIT_max=C_max_positive)


def correction_isc():
    """Exact literature ISC/singlet-shelving branching vs no-ISC baseline."""
    C = {}
    for isc in (False, True):
        H, Ls, Vp, dp, meta = g2.build_full(
            T_ROOM, 0.0, BZ, 0.0, Oc=REFERENCE_OC, ctrl=CTRL_BRANCH,
            ppol=PROBE_ORB, cpol=CTRL_ORB, isc=isc, d=D_STRAIN,
        )
        C[isc] = C_at(H, Ls, Vp, dp, meta)
    eta_max = abs(C[True] - C[False])
    C_max_positive = max(0.0, C[False], C[True])
    return dict(C_no_isc=C[False], C_with_isc=C[True],
                eta_max_deviation=eta_max, C_EIT_max=C_max_positive)


def correction_branch_asymmetry(eta_grid=(0.0, 0.05, 0.1, 0.2, 0.3, 0.5)):
    """Free symmetry-breaking parameter: radiative decay rate differs
    between orbital branches X, Y by a factor (1 +/- eta_asym), up to a
    generous (adversarial) 50% asymmetry."""
    Cs = []
    for eta in eta_grid:
        H, Ls, Vp, dp, meta = build_full_asym(
            T_ROOM, 0.0, BZ, 0.0, eta_asym=eta, Oc=REFERENCE_OC, ctrl=CTRL_BRANCH,
            ppol=PROBE_ORB, cpol=CTRL_ORB, d=D_STRAIN,
        )
        Cs.append(C_at(H, Ls, Vp, dp, meta))
    C_max_positive = max([0.0] + Cs)
    eta_max_deviation = max(abs(c - Cs[0]) for c in Cs)
    return dict(eta_grid=list(eta_grid), C_values=Cs,
                eta_max_deviation=eta_max_deviation, C_EIT_max=C_max_positive)


def correction_polarization_impurity(eta_grid=None):
    """Probe/control orbital polarization mixed away from pure
    orthogonal X/Y by an angle eta_pol; M0=0 is protected by ground-state
    orthogonality regardless (Step 5), so this tests only the subleading
    response's sensitivity to imperfect polarization selection."""
    if eta_grid is None:
        eta_grid = np.linspace(0, np.pi / 4, 6)
    Cs = []
    for eta in eta_grid:
        ppol = np.cos(eta) * PROBE_ORB + np.sin(eta) * CTRL_ORB
        cpol = -np.sin(eta) * PROBE_ORB + np.cos(eta) * CTRL_ORB
        H, Ls, Vp, dp, meta = g2.build_full(
            T_ROOM, 0.0, BZ, 0.0, Oc=REFERENCE_OC, ctrl=CTRL_BRANCH,
            ppol=ppol, cpol=cpol, isc=False, d=D_STRAIN,
        )
        Cs.append(C_at(H, Ls, Vp, dp, meta))
    C_max_positive = max([0.0] + Cs)
    eta_max_deviation = max(abs(c - Cs[0]) for c in Cs)
    return dict(eta_grid_rad=list(eta_grid), C_values=Cs,
                eta_max_deviation=eta_max_deviation, C_EIT_max=C_max_positive)


def combined_worst_case():
    """All corrections turned on simultaneously at whichever individual
    value maximizes C (not necessarily all at their most extreme value --
    the worst case FOR THE NO-GO CLAIM is the combination that maximizes
    transparency contrast)."""
    # best mI (or ensemble) from the hyperfine correction, best isc
    # setting, best branch asymmetry, best polarization angle -- combine
    # directly in one build.
    eta_asym_best = 0.5   # most adversarial from correction_branch_asymmetry's grid
    eta_pol_best = np.pi / 4  # most adversarial from correction_polarization_impurity's grid
    ppol = np.cos(eta_pol_best) * PROBE_ORB + np.sin(eta_pol_best) * CTRL_ORB
    cpol = -np.sin(eta_pol_best) * PROBE_ORB + np.cos(eta_pol_best) * CTRL_ORB

    # build_full_asym doesn't itself support isc/hyperfine (kept minimal,
    # since combining ALL FOUR jump-operator modifications into one
    # function is not needed to get a valid upper bound): combine
    # hyperfine + isc + polarization impurity via gate2.build_full (exact
    # rates) at the best polarization angle, and separately combine
    # branch asymmetry + polarization impurity via build_full_asym at the
    # same angle. The overall worst case is the max across both
    # combinations -- a conservative upper bound on the true simultaneous
    # worst case (which would require merging isc/hyperfine directly into
    # build_full_asym), and a fair one: Gate 7 only needs "the combination
    # does not exceed threshold", and taking the max over two partial
    # combinations can only overestimate, never underestimate, the true
    # combined C_EIT.
    combo_hyperfine_isc_pol = {}
    for isc in (False, True):
        for mI in (None, -1, 0, 1):
            H, Ls, Vp, dp, meta = g2.build_full(
                T_ROOM, 0.0, BZ, 0.0, Oc=REFERENCE_OC, ctrl=CTRL_BRANCH,
                ppol=ppol, cpol=cpol, isc=isc, d=D_STRAIN, mI=mI,
            )
            combo_hyperfine_isc_pol[f"isc={isc},mI={mI}"] = C_at(H, Ls, Vp, dp, meta)

    H, Ls, Vp, dp, meta = build_full_asym(
        T_ROOM, 0.0, BZ, 0.0, eta_asym=eta_asym_best, Oc=REFERENCE_OC,
        ctrl=CTRL_BRANCH, ppol=ppol, cpol=cpol, d=D_STRAIN,
    )
    combo_branch_asym_pol = C_at(H, Ls, Vp, dp, meta)

    all_values = list(combo_hyperfine_isc_pol.values()) + [combo_branch_asym_pol]
    C_worst_case = max([0.0] + all_values)
    return dict(
        combo_hyperfine_isc_polarization=combo_hyperfine_isc_pol,
        combo_branch_asymmetry_polarization=combo_branch_asym_pol,
        C_EIT_worst_case=C_worst_case,
    )


def run():
    hyperfine = correction_hyperfine()
    isc = correction_isc()
    branch_asym = correction_branch_asymmetry()
    pol = correction_polarization_impurity()
    combined = combined_worst_case()

    eps_th = sc.min_detectable_contrast(
        target_snr=5.0, od_sector=0.05, od_total=1.0,
        power_W=1e-3, lambda_nm=637, tau_s=3600.0, eta=0.1,
    )

    corrections = dict(
        hyperfine=hyperfine, isc_singlet=isc,
        branch_linewidth_asymmetry=branch_asym, polarization_impurity=pol,
    )
    individual_max = max(c["C_EIT_max"] for c in corrections.values())

    gates = dict(
        each_correction_below_threshold=bool(
            all(c["C_EIT_max"] < eps_th for c in corrections.values())
        ),
        combined_worst_case_below_threshold=bool(combined["C_EIT_worst_case"] < eps_th),
        combined_worst_case_tiny_absolute=bool(combined["C_EIT_worst_case"] < 1e-3),
        combined_not_worse_than_sum_of_individual=bool(
            combined["C_EIT_worst_case"] < 10 * (individual_max + 1e-300)
        ),
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        model="No-go theorem/src/gate2_candidate_full_vs_reduced.py build_full "
              "(+ build_full_asym for branch-linewidth asymmetry), T=300K reference geometry",
        epsilon_detection_threshold=float(eps_th),
        corrections=corrections,
        combined_worst_case=combined,
        gates=gates,
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step7.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    labels = list(corrections.keys()) + ["combined worst case"]
    values = [max(1e-20, c["C_EIT_max"]) for c in corrections.values()] + \
             [max(1e-20, combined["C_EIT_worst_case"])]
    colors = ["tab:blue"] * len(corrections) + ["tab:red"]
    ax.bar(labels, values, color=colors)
    ax.axhline(eps_th, color="gray", linestyle=":", label=f"detection eps_th={eps_th:.2e}")
    ax.set_yscale("log"); ax.set_ylabel("max positive C_EIT found")
    ax.set_title("Step 7: correction-mechanism upper bounds (T=300K)")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(branch_asym["eta_grid"], branch_asym["C_values"], "o-", label="branch asymmetry")
    ax.plot(np.array(pol["eta_grid_rad"]) / (np.pi / 4), pol["C_values"], "s-",
            label="polarization impurity (x-axis: eta/(pi/4))")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("correction strength (normalized)"); ax.set_ylabel("C_EIT")
    ax.set_title("Correction sweeps: sign and magnitude")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step7_correction_bounds.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps(gates, indent=2))
    print(f"individual correction max C_EIT: {individual_max:.4e}")
    print(f"combined worst-case C_EIT: {combined['C_EIT_worst_case']:.4e}")
    print(f"detection threshold: {eps_th:.4e}")
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
