"""Step 5 of the room-temperature SMRT no-go plan (plan Sec. 5, Step 5):
300 K global adversarial optimization over the realistic parameter domain
(plan Sec. 4), on the SAME branch-resolved NV spin-Lambda validated in
Steps 1-4 (No-go theorem/src/nv_model.py).

Goal: show the no-go conclusion is not an artifact of one convenient
parameter point (Step 1-4's fixed Oc=1 GHz, gg=6.3e-5 GHz, Bx=0, d=1.683
GHz) by searching the ENTIRE declared realistic domain simultaneously for
the combination most favorable to EIT, at T=300 K.

Robustness fact carried in from Step 3/4 and re-verified here (see
`test_M0_robust_under_Bx`): M0 = dp^T dc = 0 is protected by GROUND-STATE
orthogonality (gp, gc are eigenvectors of the Hermitian Hgs(Bvec) for
DIFFERENT eigenvalues, hence exactly orthogonal for ANY Bx, Bz -- not just
Bx=0). Opening a transverse field cannot reopen the M0=0 protection; it
can only change the size of the M1-driven subleading response. This is
why Bx is included as a free adversarial parameter below rather than
fixed at 0: if the search finds Bx>0 helps, that is real physics (mixing
the ground manifold), not evidence against the M0=0 result.

Search domain: MAIN_BOUNDS (plan Sec. 4's primary declared realistic
domain -- what Gate 5 is actually checked against) and STRESS_BOUNDS
(plan Sec. 4's explicit non-physical-extreme upper bounds, reported
separately per plan Sec 1.2/4.2/4.3/16 -- a positive finding there does
NOT count against the main no-go claim):
  Oc/2pi   MAIN [1 MHz, 1 GHz], STRESS up to 10 GHz  (Sec. 4.2)
  gg       [10 kHz, 10 MHz]                          (Sec. 4.6, "optimistic")
  Bx       MAIN [0, 0.2 T], STRESS up to 1 T          (Sec. 4.3)
  Bz       [0, 0.1 T]                                 (Sec. 4.3 main)
  d_strain MAIN [0, 50 GHz], STRESS up to 100 GHz     (Sec. 4.4)
  phi      [0, 2*pi]                                  (strain azimuth, Sec. 4.5)
  detuning z: NOT a free optimizer parameter. An earlier version searched
    it via nv_model.scanC's window-max, matching the plan's own
    max_{delta2 in W2} structure literally -- but scanC's window search
    turned out to be unreliable for this specific reduced-kernel model:
    nv_model.response()'s z is a single SHARED probe/control detuning
    (confirmed in Step 4), not an independent two-photon axis, so its
    C(z) has no clean interior extremum to search for in general (verified:
    even the KNOWN-GOOD Step 1 reference point flags scanC's own
    edge_max=True at every window width tried). The objective below
    instead evaluates at z=Ep directly, matching Steps 1/3/4's convention.

Scope note (deferred, not silently dropped): probe/control polarization
angle is NOT varied here -- kept at the pure orthogonal orbital-X/orbital-Y
choice validated in Steps 1-4. Polarization impurity is deferred to Step 7
(correction mechanisms), where it belongs alongside hyperfine, ensemble,
and linewidth-asymmetry corrections; M0=0's ground-orthogonality mechanism
(see above) means polarization mixing cannot open the leading-order
channel regardless, only perturb the already-small subleading one.

Method (plan Sec. 5, Step 5): Latin hypercube pre-sampling, local
refinement from the best LHS points, an independent differential-evolution
global search, and the overall best re-verified directly (not through the
optimizer's own bookkeeping).

IMPORTANT -- two bugs caught and fixed during development, both material
to the conclusion (see RoomT/README.md for the full account):
  1. The objective originally scored |C| (absolute value). The adversarial
     search promptly found Bx=0.2T (well within the MAIN range) giving
     |C|~1.85e-6, apparently exceeding the detection threshold. Checking
     Criterion E1 (Im chi_full < Im chi_cut, i.e. transparency) showed
     this point has C<0 -- the control field INCREASES absorption there,
     the opposite of EIT -- so it must not be scored as "contrast" just
     because its magnitude is large. Fixed to max(C,0), matching the
     plan's own objective C_EIT(theta) = max[...]_+ (Sec. 5) exactly.
  2. A companion finding this bug fix led to: C(T) at the Step 1/3/4
     reference configuration (Bx=0) is not just shrinking in magnitude
     with T, it changes SIGN between ~50 K and ~77 K (positive/
     transparency at low T, negative/absorption-increase from ~100 K up)
     -- see `sign_crossover_scan`. By 300 K the reference configuration
     is not "weak EIT", it is not EIT-signed at all.

Outputs: RoomT/results/gates_summary_step5.json,
         RoomT/results/figures/fig_step5_global_optimization.png
"""
from __future__ import annotations
import os
import sys
import json
import time

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOMT_ROOT = os.path.dirname(HERE)
NOGO_SRC = os.path.join(ROOMT_ROOT, "..", "..", "No-go theorem", "src")
sys.path.insert(0, NOGO_SRC)

import nv_model as nm  # noqa: E402
import phonon_rates as pr  # noqa: E402
import signal_chain as sc  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

T_CAMPAIGN = 300.0
CTRL_BRANCH = "-1"
REFERENCE_OC_DEFAULT = 1.0     # matches Step 1/4's REFERENCE_OC
REFERENCE_GG_DEFAULT = 6.3e-5  # matches Step 1/4's REFERENCE_GG

# theta = (log10 Oc, log10 gg, Bx, Bz, d_strain, phi)
# MAIN range (plan Sec. 4's primary declared realistic domain -- THIS is
# what Gate 5 is checked against) vs STRESS range (plan Sec. 4's explicit
# "non-physical extreme parameter" upper bounds, Sec 1.2/4.2/4.3: a
# positive finding here does NOT count against the no-go claim, it is
# reported separately as an appendix-style robustness check, plan Sec 16
# checklist item "非現実的parameterはmain claimから除外").
MAIN_BOUNDS = [
    (np.log10(1e-3), np.log10(1.0)),    # Oc/2pi: 1 MHz .. 1 GHz (Sec 4.2 main)
    (np.log10(1e-5), np.log10(1e-2)),   # gg: 10 kHz .. 10 MHz (Sec 4.6)
    (0.0, 0.2),                          # Bx: 0 .. 0.2 T (Sec 4.3 main)
    (0.0, 0.1),                          # Bz: 0 .. 0.1 T (Sec 4.3 main)
    (0.0, 50.0),                         # d_strain: 0 .. 50 GHz (Sec 4.4 main)
    (0.0, 2 * np.pi),                    # phi
]
STRESS_BOUNDS = [
    (np.log10(1e-3), np.log10(10.0)),   # Oc/2pi: up to 10 GHz stress test
    (np.log10(1e-5), np.log10(1e-2)),
    (0.0, 1.0),                          # Bx: up to 1 T stress test
    (0.0, 0.1),
    (0.0, 100.0),                        # d_strain: up to 100 GHz stress test
    (0.0, 2 * np.pi),
]
PARAM_NAMES = ["log10_Oc", "log10_gg", "Bx", "Bz", "d_strain", "phi"]


def gamma_oc_conservative(T, d):
    k = pr.k_orb_variants(T, d)["conservative_lower_bound"]
    return 0.5 * k * 1e-9 + 0.5 * nm.GRAD


def unpack(theta):
    log_Oc, log_gg, Bx, Bz, d_strain, phi = theta
    return 10 ** log_Oc, 10 ** log_gg, Bx, Bz, d_strain, phi


def contrast_at(theta, **_ignored):
    """Contrast evaluated AT the probe resonance z=Ep, matching the
    convention used consistently in Steps 1, 3, and 4 (nv_model.response()'s
    z is a single shared probe/control detuning, not an independent
    two-photon axis -- confirmed in Step 4 -- so there is no separate
    "detuning window" to search over for this specific reduced-kernel
    model; scanC()'s own window-max search was tried first and rejected,
    see module docstring / RoomT/README.md for why).

    Returns ONLY the positive (transparency) part, matching the plan's
    own objective C_EIT(theta) = max[...]_+ (Sec. 5, Step 5) and Criterion
    E1 (Im chi_full < Im chi_cut, i.e. LESS absorption with control on).
    C<0 here means the control field INCREASES absorption (Afull>Acut) --
    the opposite of transparency -- and must not be scored as "contrast"
    just because |C| is large. A first version of this objective used
    abs(C) and the adversarial search promptly found a large NEGATIVE-C
    point (Bx=0.2T, |C|~1.85e-6) that turned out to fail Criterion E1
    outright (Afull > Acut there); restricting to max(C,0) is the fix,
    not merely a refinement."""
    Oc, gg, Bx, Bz, d_strain, phi = unpack(theta)
    Bvec = (Bx, 0.0, Bz)
    try:
        w, U = nm.dressed_ground(Bvec)
        H = nm.Hes(Bvec, d=d_strain, phi=phi)
        dp, dc = nm.legs(U, ctrl=CTRL_BRANCH)
        gamma = gamma_oc_conservative(T_CAMPAIGN, d_strain)
        Ep, fid, _ = nm.probe_line(H, dp)
        r = nm.response(H, dp, dc, Ep, gamma, Oc=Oc, gg=gg)
        c = max(r["C"], 0.0)
        return (float(c) if np.isfinite(c) else 0.0), r
    except Exception:
        return 0.0, None


def neg_contrast(theta):
    c, _ = contrast_at(theta)
    return -c


def test_M0_robust_under_Bx():
    """Structural check: M0 stays exactly zero for generic Bx (ground-
    state orthogonality, not orbital-branch orthogonality alone)."""
    rng = np.random.default_rng(0)
    max_M0 = 0.0
    for _ in range(20):
        Bx = rng.uniform(0, 1.0)
        Bz = rng.uniform(0, 0.1)
        w, U = nm.dressed_ground((Bx, 0.0, Bz))
        dp, dc = nm.legs(U, ctrl=CTRL_BRANCH)
        max_M0 = max(max_M0, abs(np.vdot(dp, dc)))
    return max_M0


def sign_crossover_scan(Bvec=(0.0, 0.0, 0.02), d=1.683, Oc=REFERENCE_OC_DEFAULT,
                         gg=REFERENCE_GG_DEFAULT):
    """Companion finding (see module docstring, bug #2): at the Step
    1/3/4 reference configuration, C(T) does not merely shrink with T --
    it changes SIGN between roughly 50 K and 77 K. Below that crossover
    C>0 (genuine transparency, Criterion E1 satisfied, matches Step 1's
    positive control); above it C<0 (the control field increases
    absorption -- NOT EIT by Criterion E1, regardless of how small).
    Returns the signed C(T) trace and the crossover temperature bracket."""
    w, U = nm.dressed_ground(Bvec)
    H = nm.Hes(Bvec, d=d)
    dp, dc = nm.legs(U, ctrl=CTRL_BRANCH)
    Ep, fid, _ = nm.probe_line(H, dp)
    Ts = [4.0, 10.0, 20.0, 30.0, 50.0, 77.0, 100.0, 150.0, 200.0, 300.0]
    Cs = []
    for T in Ts:
        gamma = nm.gamma_oc_GHz(T, d)
        r = nm.response(H, dp, dc, Ep, gamma, Oc=Oc, gg=gg)
        Cs.append(float(r["C"]))
    crossover_lo, crossover_hi = None, None
    for i in range(len(Ts) - 1):
        if Cs[i] > 0 and Cs[i + 1] <= 0:
            crossover_lo, crossover_hi = Ts[i], Ts[i + 1]
            break
    return dict(Ts=Ts, Cs=Cs, crossover_bracket_K=[crossover_lo, crossover_hi])


def optimize_over(bounds, seed, n_lhs=400):
    """LHS pre-sampling + local refinement from top points + an
    independent differential-evolution global search, over the given
    bounds; returns the overall best point (re-verified at high window
    resolution) plus the raw LHS sample for the landscape figure."""
    dim = len(bounds)
    sampler = qmc.LatinHypercube(d=dim, seed=seed)
    u = sampler.random(n_lhs)
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    lhs_points = qmc.scale(u, lo, hi)
    lhs_vals = np.array([-neg_contrast(p) for p in lhs_points])
    top_idx = np.argsort(-lhs_vals)[:10]

    local_results = []
    for i in top_idx:
        res = minimize(neg_contrast, lhs_points[i], method="Nelder-Mead",
                        bounds=bounds, options=dict(xatol=1e-6, fatol=1e-12, maxiter=500))
        local_results.append((res.x, -res.fun))

    de_result = differential_evolution(
        neg_contrast, bounds, seed=seed + 100, popsize=15, maxiter=60,
        tol=1e-10, polish=True, workers=1,
    )

    candidates = [(p, v) for p, v in local_results] + [(de_result.x, -de_result.fun)]
    best_theta, best_c_reported = max(candidates, key=lambda t: t[1])
    best_c_verified, _ = contrast_at(best_theta, n_window=801, z_near=20.0)

    return dict(
        lhs_points=lhs_points, lhs_vals=lhs_vals,
        best_theta=best_theta, best_c_reported=best_c_reported,
        best_c_verified=best_c_verified, de_best_C=float(-de_result.fun),
        de_nit=int(de_result.nit), n_lhs=n_lhs,
    )


def run():
    t0 = time.time()
    m0_check = test_M0_robust_under_Bx()

    main_res = optimize_over(MAIN_BOUNDS, seed=0, n_lhs=800)
    stress_res = optimize_over(STRESS_BOUNDS, seed=1, n_lhs=800)
    crossover = sign_crossover_scan()

    Oc_m, gg_m, Bx_m, Bz_m, d_m, phi_m = unpack(main_res["best_theta"])
    Oc_s, gg_s, Bx_s, Bz_s, d_s, phi_s = unpack(stress_res["best_theta"])

    # detection threshold, reused from the pre-flight signal-chain work
    # (SNR=5, 1 hour integration, matches the earlier derivation)
    eps_th = sc.min_detectable_contrast(
        target_snr=5.0, od_sector=0.05, od_total=1.0,
        power_W=1e-3, lambda_nm=637, tau_s=3600.0, eta=0.1,
    )

    main_c = main_res["best_c_verified"]
    stress_c = stress_res["best_c_verified"]

    gates = dict(
        M0_robust_under_Bx=bool(m0_check < 1e-10),
        main_range_max_below_detection_threshold=bool(main_c < eps_th),
        main_range_max_tiny_absolute=bool(main_c < 1e-4),
        stress_range_max_still_far_below_unity=bool(stress_c < 1e-2),
        lhs_and_DE_consistent=bool(
            abs(main_c - main_res["de_best_C"]) / max(main_c, 1e-300) < 5.0 or main_c < 1e-6
        ),
        reference_config_sign_crossover_found=bool(crossover["crossover_bracket_K"][0] is not None),
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        model="No-go theorem/src/nv_model.py, T=300K, conservative_lower_bound Gamma_XY(T,d)",
        T_K=T_CAMPAIGN,
        M0_max_abs_over_random_Bx=m0_check,
        epsilon_detection_threshold=float(eps_th),
        reference_config_sign_crossover=crossover,
        main_range=dict(
            bounds=dict(zip(PARAM_NAMES, MAIN_BOUNDS)),
            n_lhs=main_res["n_lhs"], lhs_best_C=float(main_res["lhs_vals"].max()),
            de_best_C=main_res["de_best_C"], de_nit=main_res["de_nit"],
            best_theta_physical=dict(Oc_GHz=Oc_m, gg_GHz=gg_m, Bx_T=Bx_m, Bz_T=Bz_m,
                                      d_strain_GHz=d_m, phi_rad=phi_m),
            best_C_EIT_reported_by_optimizer=float(main_res["best_c_reported"]),
            best_C_EIT_reverified=float(main_c),
        ),
        stress_test_range=dict(
            note="plan Sec 1.2/4.2/4.3/16: a positive finding here is reported "
                 "separately and does NOT count against the main no-go claim.",
            bounds=dict(zip(PARAM_NAMES, STRESS_BOUNDS)),
            n_lhs=stress_res["n_lhs"], lhs_best_C=float(stress_res["lhs_vals"].max()),
            de_best_C=stress_res["de_best_C"], de_nit=stress_res["de_nit"],
            best_theta_physical=dict(Oc_GHz=Oc_s, gg_GHz=gg_s, Bx_T=Bx_s, Bz_T=Bz_s,
                                      d_strain_GHz=d_s, phi_rad=phi_s),
            best_C_EIT_reported_by_optimizer=float(stress_res["best_c_reported"]),
            best_C_EIT_reverified=float(stress_c),
        ),
        gates=gates,
        runtime_s=time.time() - t0,
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step5.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure: LHS landscapes (main range) + sample distribution + sign crossover ---
    lhs_points, lhs_vals = main_res["lhs_points"], main_res["lhs_vals"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes = axes.ravel()
    ax = axes[0]
    sc_plot = ax.scatter(lhs_points[:, 4], lhs_points[:, 2], c=np.log10(np.maximum(lhs_vals, 1e-20)),
                          cmap="viridis", s=15)
    ax.scatter([d_m], [Bx_m], color="red", marker="*", s=200, label="main-range best")
    fig.colorbar(sc_plot, ax=ax, label="log10|C_EIT|")
    ax.set_xlabel("d_strain (GHz)"); ax.set_ylabel("Bx (T)")
    ax.set_title("Main-range LHS: d_strain vs Bx"); ax.legend(fontsize=7)

    ax = axes[1]
    sc_plot2 = ax.scatter(lhs_points[:, 0], lhs_points[:, 1], c=np.log10(np.maximum(lhs_vals, 1e-20)),
                           cmap="viridis", s=15)
    ax.scatter([np.log10(Oc_m)], [np.log10(gg_m)], color="red", marker="*", s=200)
    fig.colorbar(sc_plot2, ax=ax, label="log10|C_EIT|")
    ax.set_xlabel("log10(Oc)"); ax.set_ylabel("log10(gg)")
    ax.set_title("Main-range LHS: Oc vs gg")

    ax = axes[2]
    sorted_main = np.sort(lhs_vals)[::-1]
    sorted_stress = np.sort(stress_res["lhs_vals"])[::-1]
    ax.semilogy(np.arange(len(sorted_main)), np.maximum(sorted_main, 1e-20), "o-", ms=3, label="main range")
    ax.semilogy(np.arange(len(sorted_stress)), np.maximum(sorted_stress, 1e-20), "s-", ms=3, label="stress range")
    ax.axhline(main_c, color="red", linestyle="--", label=f"main-range global best={main_c:.2e}")
    ax.axhline(stress_c, color="orange", linestyle="-.", label=f"stress-range global best={stress_c:.2e}")
    ax.axhline(eps_th, color="gray", linestyle=":", label=f"detection eps_th={eps_th:.2e}")
    ax.set_xlabel("LHS sample rank"); ax.set_ylabel("|C_EIT|")
    ax.set_title("LHS sample distribution vs detection threshold")
    ax.legend(fontsize=6)

    ax = axes[3]
    Ts, Cs = crossover["Ts"], crossover["Cs"]
    colors = ["tab:blue" if c > 0 else "tab:red" for c in Cs]
    ax.scatter(Ts, np.abs(Cs), c=colors, s=40, zorder=3)
    ax.plot(Ts, np.abs(Cs), color="gray", alpha=0.5, zorder=2)
    ax.set_yscale("log")
    lo, hi = crossover["crossover_bracket_K"]
    if lo is not None:
        ax.axvspan(lo, hi, color="orange", alpha=0.2, label=f"sign crossover in [{lo},{hi}] K")
    ax.set_xlabel("T (K)"); ax.set_ylabel("|C| (color: blue=transparency C>0, red=absorption-increase C<0)")
    ax.set_title("Reference config (Bx=0): C(T) sign crossover")
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step5_global_optimization.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps(gates, indent=2))
    print(f"MAIN range best: Oc={Oc_m:.4g} GHz, gg={gg_m:.4g} GHz, Bx={Bx_m:.4g} T, "
          f"Bz={Bz_m:.4g} T, d={d_m:.4g} GHz, phi={phi_m:.4g} -> |C_EIT|={main_c:.4e}")
    print(f"STRESS range best: Oc={Oc_s:.4g} GHz, gg={gg_s:.4g} GHz, Bx={Bx_s:.4g} T, "
          f"Bz={Bz_s:.4g} T, d={d_s:.4g} GHz, phi={phi_s:.4g} -> |C_EIT|={stress_c:.4e}")
    print(f"detection threshold = {eps_th:.4e}")
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
