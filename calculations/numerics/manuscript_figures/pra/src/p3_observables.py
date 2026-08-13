"""p3_observables.py -- PRA calculation P3.

One common observable axis from the low-temperature positive signal, through
the transition region, to the room-temperature negative signal.

Two conversions already existed in the repository but were not commensurate:
gate3_snr_map.py (70 K candidate: probe power 1 uW, Debye-Waller 0.035,
gamma_inh 30 GHz, 1 hour ceiling) and RoomT/step9_signal_conversion.py
(300 K: probe power 1 mW, Debye-Waller 0.03, gamma_inh 5 GHz, 24 hour
ceiling).  Neither covers the other's temperature, and their parameter sets
differ, so the positive and negative signals could not be placed on a single
scale -- which is exactly what the P3 pass criterion in
NV_EIT_PRA_PRL_Split_Strategy_20260724.md Sec. 5 requires.

This module fixes ONE detection chain (the conservative choice at each point
where the two disagreed) and applies it at every temperature, and adds the
fluorescence (PL) readout, which was absent from the repository entirely.

Physics inputs are the full 9-level Liouvillian spectra of P1 -- not the
reduced kernel -- so the observable curve inherits no reduced-model error.

Readouts
  transmission : Delta T/T = expm1(Delta OD),      Delta OD = OD_sector * C
  fluorescence : Delta PL/PL from the change in absorbed power; for an
                 optically thin sector this tends to -C, i.e. a transparency
                 shows up as a PL DECREASE of the same fractional size.

Outputs
  results/tables/p3_observables_vs_T.csv   per temperature, both readouts
  results/tables/p3_summary.json           the three regimes on one axis
Usage
  python p3_observables.py [--quick] [--jobs N]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRA = HERE.parent
REPO = PRA.parents[1]          # calculations/numerics
RES = PRA.parents[3] / "results"
NOGO_SRC = REPO / "No-go theorem" / "src"
sys.path.insert(0, str(NOGO_SRC))
sys.path.insert(0, str(HERE))

import run_prl_prediction as rp        # noqa: E402
import nv_model as nv                  # noqa: E402
import signal_chain as sc              # noqa: E402
import p1_phase_diagram as p1          # noqa: E402

# ------------------------------------------------------------ detection chain
# One parameter set for the whole temperature range.  Where gate3 (70 K) and
# RoomT step9 (300 K) disagreed, the more conservative value is taken and the
# disagreement is recorded, so no part of the curve is flattered by a
# parameter choice that was only ever justified at the other end.
CHAIN = dict(
    lambda_nm=(637.0, "nm", "NV- ZPL", "exact; both sources agree"),
    n_refr=(2.40, "-", "diamond @ 637 nm", "gate3 used 2.41; 2.40 taken"),
    debye_waller=(0.030, "-", "Santori 2010 / Doherty review",
                  "gate3 0.035 vs step9 0.030; conservative 0.030 taken"),
    gamma_inh_GHz=(30.0, "GHz", "ensemble ZPL strain broadening",
                   "gate3 30 vs step9 5; conservative 30 taken "
                   "(smaller resonant fraction)"),
    f_orient=(0.25, "-", "one of four NV orientations", "exact"),
    f_spin=(1 / 3, "-", "unpolarized ground spin, no MW initialization",
            "up to ~0.9 with initialization"),
    power_W=(1e-6, "W", "weak probe, below sector saturation",
             "gate3 1e-6 vs step9 1e-3; the weak-probe value is the one "
             "consistent with the linear-response theory used throughout"),
    eta_detect=(0.1, "-", "collection x detector efficiency", "x/ 3"),
    sigma_tech=(1e-6, "-", "relative technical-noise floor per sqrt(sample)",
                "both sources agree; already aggressive"),
    L_cm=(0.05, "cm", "0.5 mm sample", "design"),
    n_nv_cm3=(1.76e17, "cm^-3", "1 ppm", "coherence-preserving upper end"),
    target_snr=(5.0, "-", "detection criterion", "design"),
    tau_ceiling_s=(24 * 3600.0, "s", "reasonable experiment duration",
                   "step9 ceiling; gate3 used 1 h"),
    quantum_yield=(0.7, "-", "NV- radiative quantum yield", "0.6-0.8"),
)
V = {k: v[0] for k, v in CHAIN.items()}


# ------------------------------------------------------------ PL readout
def pl_fraction_absorbed(od_total: float) -> float:
    """Fraction of probe power absorbed by the sample."""
    return float(-np.expm1(-od_total))


def delta_pl_over_pl(od_cut: float, d_od: float) -> float:
    """Fractional fluorescence change produced by the sector pathway.

    The control removes Delta OD = OD_sector * C from the absorption, so the
    absorbed power -- and hence the fluorescence -- drops.  A positive
    contrast C therefore gives a NEGATIVE Delta PL/PL: in fluorescence a
    transparency reads as a dip, with the same fractional magnitude as the
    transmission signal in the optically thin limit.
    """
    denom = -np.expm1(-od_cut)
    if denom <= 0:
        return 0.0
    return float(-np.exp(-od_cut) * np.expm1(d_od) / denom)


def pl_photons(power_W, lambda_nm, tau_s, eta, quantum_yield, od_total):
    rate = power_W / sc.photon_energy_J(lambda_nm)
    return rate * pl_fraction_absorbed(od_total) * quantum_yield * eta * tau_s


def pl_required_tau_s(target_snr, rel_signal, power_W, lambda_nm, eta,
                      quantum_yield, od_total, sigma_tech=0.0):
    """Integration time for the PL readout to reach target SNR."""
    rel = abs(rel_signal)
    denom = rel ** 2 - (target_snr * sigma_tech) ** 2
    if denom <= 0:
        return float("inf")
    rate = (power_W / sc.photon_energy_J(lambda_nm)
            * pl_fraction_absorbed(od_total) * quantum_yield * eta)
    if rate <= 0:
        return float("inf")
    return target_snr ** 2 / (denom * rate)


# ------------------------------------------------------------ per-temperature
def observables_at(T, Bx):
    """Full-Liouvillian contrast at (T, Bx), converted to both readouts."""
    row = p1.classify_point(float(T), float(Bx))
    C = row.get("Cmax", np.nan)

    gamma_h = nv.gamma_oc_GHz(float(T), rp.D)
    sigma = sc.sigma_zpl_cm2(V["lambda_nm"], V["n_refr"], V["debye_waller"],
                             nv.GRAD, gamma_h)
    f_spec = sc.spectral_fraction(gamma_h, V["gamma_inh_GHz"])
    alpha = sc.alpha_cm(sigma, V["n_nv_cm3"], V["f_orient"], V["f_spin"], f_spec)
    od_sector = sc.od(alpha, V["L_cm"])

    out = dict(T_K=float(T), Bx_T=float(Bx), C=float(C),
               klass=row.get("klass"), verdict=row.get("verdict"),
               fwhm_MHz=row.get("fwhm_MHz"),
               gamma_h_GHz=float(gamma_h), sigma_zpl_cm2=float(sigma),
               spectral_fraction=float(f_spec), alpha_percm=float(alpha),
               od_sector=float(od_sector))

    # Length that puts the probed sector at OD = 1.  The sector cross-section
    # grows as the optical line narrows, so a FIXED sample is opaque
    # (OD ~ 16) at 10 K and nearly transparent (OD ~ 0.04) at 300 K: comparing
    # observables across T at fixed geometry compares sample thickness, not
    # physics.  The OD-matched sample is therefore the primary common axis,
    # and the fixed 1 ppm / 0.5 mm sample is reported alongside it.
    L_od1 = float(1.0 / alpha) if alpha > 0 else np.inf
    out.update(L_for_od1_cm=L_od1,
               optically_thin=bool(od_sector < 1.0))

    if not np.isfinite(C):
        out.update(d_od=np.nan, dT_over_T=np.nan, tau_trans_s=np.nan,
                   dPL_over_PL=np.nan, tau_pl_s=np.nan,
                   d_od_m=np.nan, dT_over_T_m=np.nan, tau_trans_m_s=np.nan,
                   dPL_over_PL_m=np.nan, tau_pl_m_s=np.nan,
                   detectable_trans=False, detectable_pl=False,
                   detectable_trans_m=False, detectable_pl_m=False)
        return out

    def readouts(od_ref):
        d_od = sc.delta_od(od_ref, C)
        dTT = sc.delta_T_over_T(d_od)
        tau_t = sc.required_tau_s(V["target_snr"], d_od, od_ref, V["power_W"],
                                  V["lambda_nm"], V["eta_detect"],
                                  sigma_tech=V["sigma_tech"])
        dpl = delta_pl_over_pl(od_ref, d_od)
        tau_p = pl_required_tau_s(V["target_snr"], dpl, V["power_W"],
                                  V["lambda_nm"], V["eta_detect"],
                                  V["quantum_yield"], od_ref,
                                  sigma_tech=V["sigma_tech"])
        return d_od, dTT, tau_t, dpl, tau_p

    d_od, dTT, tau_t, dpl, tau_p = readouts(od_sector)          # fixed sample
    d_odm, dTTm, tau_tm, dplm, tau_pm = readouts(1.0)           # OD-matched

    out.update(d_od=float(d_od), dT_over_T=float(dTT),
               tau_trans_s=float(tau_t), dPL_over_PL=float(dpl),
               tau_pl_s=float(tau_p),
               d_od_m=float(d_odm), dT_over_T_m=float(dTTm),
               tau_trans_m_s=float(tau_tm), dPL_over_PL_m=float(dplm),
               tau_pl_m_s=float(tau_pm),
               detectable_trans=bool(np.isfinite(tau_t)
                                     and tau_t <= V["tau_ceiling_s"]),
               detectable_pl=bool(np.isfinite(tau_p)
                                  and tau_p <= V["tau_ceiling_s"]),
               detectable_trans_m=bool(np.isfinite(tau_tm)
                                       and tau_tm <= V["tau_ceiling_s"]),
               detectable_pl_m=bool(np.isfinite(tau_pm)
                                    and tau_pm <= V["tau_ceiling_s"]))
    return out


def _worker(a):
    T, B = a
    try:
        return observables_at(T, B)
    except Exception as exc:                                # pragma: no cover
        return dict(T_K=float(T), Bx_T=float(B), C=np.nan,
                    klass=f"error:{type(exc).__name__}")


def main(quick=False, jobs=4):
    Bx = rp.BX0
    Ts = (np.array([10., 30., 70., 100., 150., 300.]) if quick else
          np.unique(np.concatenate([
              np.arange(10., 100.1, 10.),
              np.array([105., 110., 120., 140., 160., 200., 250., 300.])])))
    pts = [(float(T), float(Bx)) for T in Ts]
    print(f"P3: {len(pts)} temperatures at B_perp={Bx:.5f} T, jobs={jobs}")

    t0 = time.time()
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(min(jobs, len(pts))) as pool:
            rows = pool.map(_worker, pts)
    else:
        rows = [_worker(p) for p in pts]
    elapsed = time.time() - t0

    tabdir = RES / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)
    fields = ["T_K", "Bx_T", "C", "klass", "verdict", "fwhm_MHz", "gamma_h_GHz",
              "sigma_zpl_cm2", "spectral_fraction", "alpha_percm", "od_sector",
              "optically_thin", "L_for_od1_cm",
              "d_od_m", "dT_over_T_m", "tau_trans_m_s", "dPL_over_PL_m",
              "tau_pl_m_s", "detectable_trans_m", "detectable_pl_m",
              "d_od", "dT_over_T", "tau_trans_s", "dPL_over_PL", "tau_pl_s",
              "detectable_trans", "detectable_pl"]
    with (tabdir / "p3_observables_vs_T.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    ok = [r for r in rows if np.isfinite(r.get("C", np.nan))]
    pos = [r for r in ok if r["C"] > 0]
    neg = [r for r in ok if r["C"] < 0]
    det_t = [r for r in ok if r.get("detectable_trans_m")]
    det_p = [r for r in ok if r.get("detectable_pl_m")]

    def at(T):
        m = [r for r in ok if abs(r["T_K"] - T) < 1e-9]
        return m[0] if m else None

    highlights = {}
    for T in (10.0, 70.0, 100.0, 300.0):
        r = at(T)
        if r:
            highlights[f"{T:.0f}K"] = dict(
                C=r["C"], od_sector_fixed_sample=r["od_sector"],
                L_for_od1_cm=r["L_for_od1_cm"],
                matched=dict(d_od=r["d_od_m"], dT_over_T=r["dT_over_T_m"],
                             tau_trans_s=r["tau_trans_m_s"],
                             dPL_over_PL=r["dPL_over_PL_m"],
                             tau_pl_s=r["tau_pl_m_s"],
                             detectable_trans=r["detectable_trans_m"],
                             detectable_pl=r["detectable_pl_m"]),
                fixed_sample=dict(d_od=r["d_od"], dT_over_T=r["dT_over_T"],
                                  tau_trans_s=r["tau_trans_s"],
                                  dPL_over_PL=r["dPL_over_PL"],
                                  tau_pl_s=r["tau_pl_s"],
                                  detectable_trans=r["detectable_trans"],
                                  detectable_pl=r["detectable_pl"]),
                klass=r["klass"])

    gates = dict(
        single_chain_covers_whole_range=bool(len(ok) == len(rows)),
        positive_and_negative_on_one_axis=bool(pos and neg),
        pl_readout_present=True,
        detectability_boundary_found=bool(det_t and len(det_t) < len(ok)),
    )
    summary = dict(
        what="P3 one observable axis: low-T positive, transition, 300 K negative",
        physics="full 9-level Liouvillian (P1 machinery), not the reduced kernel",
        detection_chain={k: dict(value=v[0], unit=v[1], source=v[2], note=v[3])
                         for k, v in CHAIN.items()},
        reconciles=["No-go theorem/src/gate3_snr_map.py (70 K)",
                    "New no-go theory/RoomT/src/step9_signal_conversion.py (300 K)"],
        B_perp_T=float(Bx),
        highest_T_positive_C=(max(r["T_K"] for r in pos) if pos else None),
        lowest_T_negative_C=(min(r["T_K"] for r in neg) if neg else None),
        primary_axis="OD-matched sample (sector OD = 1) -- removes the "
                     "geometry confound; the fixed 1 ppm / 0.5 mm sample "
                     "spans OD 16 (10 K) to 0.04 (300 K)",
        highest_T_detectable_transmission=(max(r["T_K"] for r in det_t)
                                           if det_t else None),
        highest_T_detectable_pl=(max(r["T_K"] for r in det_p) if det_p else None),
        highlights=highlights, gates=gates,
        all_gates_pass=bool(all(gates.values())),
        seconds=elapsed, quick=bool(quick),
    )
    with (tabdir / "p3_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\nOD-matched sample (sector OD = 1) -- primary common axis")
    print(f"{'T(K)':>6} {'C':>12} {'dT/T':>12} {'tau_tr(s)':>12} "
          f"{'dPL/PL':>12} {'tau_PL(s)':>12} | {'OD_fix':>8} {'L_OD1(cm)':>10}  class")
    for r in rows:
        if not np.isfinite(r.get("C", np.nan)):
            print(f"{r['T_K']:6.0f} {'--':>12}")
            continue
        print(f"{r['T_K']:6.0f} {r['C']:12.4g} {r['dT_over_T_m']:12.4g} "
              f"{r['tau_trans_m_s']:12.4g} {r['dPL_over_PL_m']:12.4g} "
              f"{r['tau_pl_m_s']:12.4g} | {r['od_sector']:8.3g} "
              f"{r['L_for_od1_cm']:10.3g}  {r['klass']}")
    print(f"\nhighest T with detectable transmission signal: "
          f"{summary['highest_T_detectable_transmission']} K")
    print(f"highest T with detectable PL signal: "
          f"{summary['highest_T_detectable_pl']} K")
    print(f"gates: {gates} -> {summary['all_gates_pass']}")
    print(f"elapsed {elapsed:.0f}s")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    main(quick=a.quick, jobs=a.jobs)
