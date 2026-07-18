"""Step 9 of the room-temperature SMRT no-go plan (plan Sec. 5, Step 9;
plan Sec. 9.3's logical order): convert the LOCAL response result (Steps
3-8: R_EIT is structurally suppressed and, at the reference point, not
even correctly signed) into an actual DETECTABLE-SIGNAL statement,
keeping the two claims separate as the plan insists -- a large optical
depth can AMPLIFY an existing small EIT response into more photons, but
it cannot manufacture a response that Steps 3-8 already showed is absent.

Reuses `No-go theorem/src/signal_chain.py` (built during this campaign's
pre-flight work, before Step 1) unmodified: `sigma_zpl_cm2`, `alpha_cm`,
`od`, `delta_od`, `transmission`, `delta_T_over_T`, `snr`,
`required_tau_s`, `min_detectable_contrast`.

Two contrast values are carried through the chain, both honestly labeled:

  C_reference = |C_EIT(300K)| at the Step 1-8 reference configuration
      (Bx=0, Bz=0.02T, d=1.683 GHz, Oc=1 GHz, gg=6.3e-5 GHz) -- this is
      the campaign's actual finding, and it is NEGATIVE (fails Criterion
      E1, Step 6) -- there is no genuine transparency signal to detect
      here at all, only an absorption increase of this tiny magnitude.

  C_generous_bound = the largest |C| magnitude found ANYWHERE in this
      campaign, sign-agnostic (Step 5's stress-test-range adversarial
      optimum, which itself FAILED Criterion E1 -- Step 6). Used here
      only as a deliberately over-generous upper bound: "even granting
      this magnitude, AND granting it the correct sign for transparency
      (which Step 6 already showed it does not have), is it detectable?"

Homogeneous linewidth at 300 K uses this campaign's OWN Gamma_XY(300K,d)
value (conservative_lower_bound model, Step 4) -- not an assumed number
-- so the ZPL cross-section itself is already suppressed by the same
phonon broadening driving the whole SMRT mechanism, a physically
DISTINCT and compounding reason detection is hard (independent of how
small C_EIT is).

Outputs: RoomT/results/gates_summary_step9.json,
         RoomT/results/figures/fig_step9_snr_map.png
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
import phonon_rates as pr  # noqa: E402
import signal_chain as sc  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

LAMBDA_NM = 637.0
N_REFR = 2.4
DEBYE_WALLER = 0.03
GAMMA_INH_GHZ = 5.0        # illustrative inhomogeneous width, GHz
D_STRAIN = 1.683
C_REFERENCE = 1.815093868469757e-09     # |C| at reference point (Step 6/8), sign-agnostic magnitude
C_GENEROUS_BOUND = 2.4029e-06           # largest |C| found anywhere (Step 5 stress range), sign-agnostic

POWER_W = 1e-3
ETA_DETECT = 0.1
TARGET_SNR = 5.0
SIGMA_TECH = 1e-6  # relative technical-noise floor (plan Sec 9.2 explicitly
# requires this ingredient). 1e-6 is already an aggressive, optimistic
# assumption -- real referenced/lock-in transmission measurements typically
# sit at 1e-4 to 1e-5 relative stability; state-of-the-art specialized
# setups can approach 1e-7. Omitting this (sigma_tech=0, pure shot-noise
# limit) was an early bug in this script: it made a sub-nanosecond-scale
# relative signal look "achievable in 12.8 hours" by integrating shot
# noise down indefinitely, when in reality no real experiment is shot-
# noise-limited at a 1e-9 relative modulation level over many-hour
# integration -- ANY realistic technical noise floor (checked down to
# 1e-9, already far beyond current capability) makes it undetectable at
# ANY integration time (see `sigma_tech_sensitivity` in the output).

N_NV_GRID = np.logspace(15, 18, 13)   # cm^-3 (MAIN: up to ~5.7 ppm)
L_GRID = np.logspace(-3, -1, 13)      # cm (MAIN: 10 um .. 1 mm)
N_NV_GRID_STRESS = np.logspace(15, 19, 13)  # cm^-3 (STRESS: up to ~57 ppm)
L_GRID_STRESS = np.logspace(-3, 0, 13)      # cm (STRESS: 10 um .. 1 cm)

REALISTIC_TAU_CEILING_S = 24 * 3600.0  # "reasonable experiment duration"


def gamma_h_300K():
    k = pr.k_orb_variants(300.0, D_STRAIN)["conservative_lower_bound"]
    return 0.5 * k * 1e-9 + 0.5 * nm.GRAD  # matches this campaign's Gamma_XY(300K) convention


def build_snr_map(contrast, n_nv_grid, L_grid, sigma_tech=SIGMA_TECH):
    gamma_h = gamma_h_300K()
    gamma_rad = nm.GRAD
    sigma = sc.sigma_zpl_cm2(LAMBDA_NM, N_REFR, DEBYE_WALLER, gamma_rad, gamma_h)
    f_spec = sc.spectral_fraction(gamma_h, GAMMA_INH_GHZ)

    tau_map = np.zeros((len(n_nv_grid), len(L_grid)))
    for i, n_nv in enumerate(n_nv_grid):
        for j, L in enumerate(L_grid):
            alpha = sc.alpha_cm(sigma, n_nv, f_spectral=f_spec)
            od_sector = sc.od(alpha, L)
            d_od = sc.delta_od(od_sector, contrast)
            tau_map[i, j] = sc.required_tau_s(
                TARGET_SNR, d_od, od_sector, POWER_W, LAMBDA_NM, ETA_DETECT, sigma_tech=sigma_tech
            )
    return dict(gamma_h_GHz=gamma_h, sigma_zpl_cm2=sigma, spectral_fraction=f_spec, tau_map=tau_map)


def sigma_tech_sensitivity(contrast, n_nv, L):
    """How undetectability depends on the assumed technical noise floor,
    at one representative (n_NV, L) point -- shows the conclusion holds
    even for optimistic (small) sigma_tech, not just the main value."""
    gamma_h = gamma_h_300K()
    sigma = sc.sigma_zpl_cm2(LAMBDA_NM, N_REFR, DEBYE_WALLER, nm.GRAD, gamma_h)
    f_spec = sc.spectral_fraction(gamma_h, GAMMA_INH_GHZ)
    alpha = sc.alpha_cm(sigma, n_nv, f_spectral=f_spec)
    od_sector = sc.od(alpha, L)
    d_od = sc.delta_od(od_sector, contrast)
    out = {}
    for st in (0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        tau = sc.required_tau_s(TARGET_SNR, d_od, od_sector, POWER_W, LAMBDA_NM, ETA_DETECT, sigma_tech=st)
        out[str(st)] = tau
    return dict(n_nv_cm3=n_nv, L_cm=L, od_sector=od_sector, rel_signal_d_od=d_od, tau_by_sigma_tech=out)


def run():
    ref = build_snr_map(C_REFERENCE, N_NV_GRID, L_GRID)
    gen = build_snr_map(C_GENEROUS_BOUND, N_NV_GRID, L_GRID)
    ref_stress = build_snr_map(C_REFERENCE, N_NV_GRID_STRESS, L_GRID_STRESS)
    gen_stress = build_snr_map(C_GENEROUS_BOUND, N_NV_GRID_STRESS, L_GRID_STRESS)

    min_tau_ref = float(np.min(ref["tau_map"]))
    min_tau_gen = float(np.min(gen["tau_map"]))
    min_tau_ref_stress = float(np.min(ref_stress["tau_map"]))
    min_tau_gen_stress = float(np.min(gen_stress["tau_map"]))
    idx_ref = np.unravel_index(np.argmin(ref["tau_map"]), ref["tau_map"].shape)
    idx_gen = np.unravel_index(np.argmin(gen["tau_map"]), gen["tau_map"].shape)
    idx_gen_stress = np.unravel_index(np.argmin(gen_stress["tau_map"]), gen_stress["tau_map"].shape)

    # converse: what contrast WOULD a "reasonable" ensemble sample
    # (n_NV=1e17 cm^-3, L=0.1 cm, within the MAIN realistic domain) need
    # to hit SNR=5 in 1 hour?
    gamma_h = ref["gamma_h_GHz"]
    sigma = ref["sigma_zpl_cm2"]
    f_spec = ref["spectral_fraction"]
    n_nv_typical, L_typical = 1e17, 0.1
    alpha_typical = sc.alpha_cm(sigma, n_nv_typical, f_spectral=f_spec)
    od_typical = sc.od(alpha_typical, L_typical)
    eps_th_1hr = sc.min_detectable_contrast(
        TARGET_SNR, od_typical, od_typical, POWER_W, LAMBDA_NM, 3600.0, ETA_DETECT, sigma_tech=SIGMA_TECH
    )

    sens_ref = sigma_tech_sensitivity(C_REFERENCE, float(N_NV_GRID[idx_ref[0]]), float(L_GRID[idx_ref[1]]))
    sens_gen_stress = sigma_tech_sensitivity(
        C_GENEROUS_BOUND, float(N_NV_GRID_STRESS[idx_gen_stress[0]]), float(L_GRID_STRESS[idx_gen_stress[1]])
    )

    def is_undetectable(tau):
        return bool(not np.isfinite(tau) or tau > REALISTIC_TAU_CEILING_S)

    gates = dict(
        reference_point_undetectable_main_domain=is_undetectable(min_tau_ref),
        generous_bound_undetectable_main_domain=is_undetectable(min_tau_gen),
        reference_point_undetectable_stress_domain=is_undetectable(min_tau_ref_stress),
        required_contrast_1hr_exceeds_reference=bool(
            not np.isfinite(eps_th_1hr) or eps_th_1hr > 10 * C_REFERENCE
        ),
        sigma_tech_1e9_already_sufficient_to_exclude_reference=bool(
            not np.isfinite(sens_ref["tau_by_sigma_tech"]["1e-09"])
        ),
    )
    gates["overall_pass"] = bool(all(gates.values()))

    def clean(x):
        return None if not np.isfinite(x) else x

    summary = dict(
        model="No-go theorem/src/signal_chain.py applied to this campaign's C_EIT(300K) values",
        physical_parameters=dict(
            lambda_nm=LAMBDA_NM, n_refr=N_REFR, debye_waller=DEBYE_WALLER,
            gamma_rad_GHz=nm.GRAD, gamma_h_300K_GHz=gamma_h,
            gamma_inh_GHz=GAMMA_INH_GHZ, spectral_fraction=f_spec,
            sigma_zpl_cm2=sigma, power_W=POWER_W, eta_detect=ETA_DETECT,
            target_snr=TARGET_SNR, sigma_tech=SIGMA_TECH,
            realistic_tau_ceiling_s=REALISTIC_TAU_CEILING_S,
        ),
        C_reference=C_REFERENCE,
        C_generous_bound=C_GENEROUS_BOUND,
        main_domain=dict(
            n_nv_range_cm3=[float(N_NV_GRID[0]), float(N_NV_GRID[-1])],
            L_range_cm=[float(L_GRID[0]), float(L_GRID[-1])],
            min_required_tau_s_at_reference=clean(min_tau_ref),
            min_required_tau_s_at_generous_bound=clean(min_tau_gen),
        ),
        stress_domain=dict(
            note="plan Sec 1.2/4/16 convention: extreme (n_NV, L) combinations "
                 "reported separately, not counted toward the main verdict.",
            n_nv_range_cm3=[float(N_NV_GRID_STRESS[0]), float(N_NV_GRID_STRESS[-1])],
            L_range_cm=[float(L_GRID_STRESS[0]), float(L_GRID_STRESS[-1])],
            min_required_tau_s_at_reference=clean(min_tau_ref_stress),
            min_required_tau_s_at_generous_bound=clean(min_tau_gen_stress),
        ),
        converse_required_contrast_for_1hr_SNR5_typical_sample=dict(
            n_nv_cm3=n_nv_typical, L_cm=L_typical, od_sector=od_typical,
            required_contrast=clean(eps_th_1hr),
        ),
        sigma_tech_sensitivity_at_reference_best_point=sens_ref,
        sigma_tech_sensitivity_at_generous_stress_extreme=sens_gen_stress,
        gates=gates,
        verdict=(
            "Theoretical no-go (Steps 3-8): R_EIT is structurally suppressed at 300K "
            "and, at the reference point, not even correctly signed for transparency. "
            "Detection no-go (this step): including a realistic technical-noise floor "
            "(plan Sec 9.2), no (NV density, thickness) combination in the realistic "
            "MAIN ensemble domain reaches SNR=5 within a 24-hour integration for either "
            "the reference-point contrast or the generous, sign-agnostic upper bound "
            "found anywhere in this campaign. An early version of this analysis omitted "
            "the technical-noise floor (sigma_tech=0, pure shot-noise limit), which made "
            "the reference point look marginally 'achievable' in 12.8 hours by "
            "integrating shot noise down indefinitely -- adding ANY realistic noise "
            "floor (checked down to 1e-9, already far beyond current capability) makes "
            "it undetectable at any integration time. Larger optical depth changes "
            "where on this map you sit, not whether the underlying (already-absent) "
            "signal exists."
        ),
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step9.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure: required integration time map (STRESS domain, generous
    # bound -- most likely to show finite structure once a realistic
    # technical-noise floor is included; the MAIN domain is mostly
    # undetectable-at-any-time, per the gates above) + sigma_tech sweep ---
    DISPLAY_CAP = 1e15  # stand-in for "inf" on the log color scale
    tau_display = np.where(np.isfinite(gen_stress["tau_map"]), gen_stress["tau_map"], DISPLAY_CAP)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    im = ax.pcolormesh(np.log10(L_GRID_STRESS), np.log10(N_NV_GRID_STRESS),
                        np.log10(np.maximum(tau_display, 1e-3)), shading="auto", cmap="viridis")
    fig.colorbar(im, ax=ax, label="log10(required tau, s) [capped at 1e15 = undetectable]")
    if np.any(tau_display < REALISTIC_TAU_CEILING_S):
        ax.contour(np.log10(L_GRID_STRESS), np.log10(N_NV_GRID_STRESS),
                   np.log10(np.maximum(tau_display, 1e-3)),
                   levels=[np.log10(REALISTIC_TAU_CEILING_S)], colors="red")
    ax.set_xlabel("log10(thickness, cm)"); ax.set_ylabel("log10(n_NV, cm^-3)")
    ax.set_title(f"STRESS domain, generous |C| bound,\nsigma_tech={SIGMA_TECH:.0e}: required tau for SNR=5")

    ax = axes[1]
    sigma_techs = sorted(float(k) for k in sens_ref["tau_by_sigma_tech"].keys())
    tau_ref_vs_st = [sens_ref["tau_by_sigma_tech"][str(st)] for st in sigma_techs]
    tau_gen_vs_st = [sens_gen_stress["tau_by_sigma_tech"][str(st)] for st in sigma_techs]
    tau_ref_disp = [min(t, DISPLAY_CAP) if np.isfinite(t) else DISPLAY_CAP for t in tau_ref_vs_st]
    tau_gen_disp = [min(t, DISPLAY_CAP) if np.isfinite(t) else DISPLAY_CAP for t in tau_gen_vs_st]
    x = [max(st, 1e-12) for st in sigma_techs]
    ax.loglog(x, tau_ref_disp, "o-", label="reference point (best n_NV,L in MAIN domain)")
    ax.loglog(x, tau_gen_disp, "s-", label="generous bound (best n_NV,L in STRESS domain)")
    ax.axhline(REALISTIC_TAU_CEILING_S, color="gray", linestyle=":", label="24h ceiling")
    ax.set_xlabel("assumed sigma_tech (relative)"); ax.set_ylabel("required tau (s), capped at 1e15")
    ax.set_title("Required integration time vs assumed technical-noise floor")
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step9_snr_map.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps(gates, indent=2))
    print(summary["verdict"])
    print(f"min required tau at reference (main domain): {min_tau_ref}")
    print(f"min required tau at generous bound (main domain): {min_tau_gen}")
    print(f"min required tau at reference (stress domain): {min_tau_ref_stress}")
    print(f"min required tau at generous bound (stress domain): {min_tau_gen_stress}")
    print(f"required contrast for 1hr SNR=5 (typical sample): {eps_th_1hr:.3e}")
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
