"""p1_phase_diagram.py -- PRA calculation P1.

Full-Liouvillian T-B_perp phase diagram of the NV(-) ZPL spin-Lambda
channel, with every grid point classified jointly by

  (i)   the SIGN of the peak sector contrast C = (A_cut - A_full)/A_cut,
  (ii)  the LINESHAPE of the control-on absorption A_full(delta_2), and
  (iii) the EIT/ATS model-selection verdict on that lineshape,

into exactly one of

  'transparency'  genuine positive transparency (C>0, robust EIT lineshape)
  'absorption'    control-induced absorption    (C<0, sign reversed)
  'ATS'           Autler-Townes doublet         (robust ATS lineshape)
  'unresolved'    below detection, or lineshape model selection inconclusive

This is the object required by NV_EIT_PRA_PRL_Split_Strategy_20260724.md
Sec. 5 (P1); its pass criterion is that sign, lineshape and the EIT/ATS
verdict are used *simultaneously* to fix the upper temperature limit of NV
EIT.  Nothing here is new physics: the spectra come from the validated
9-level Lindblad pipeline of gate2_candidate_full_vs_reduced.build_full and
the model comparison from gate1_candidate_aic_bootstrap.fit_all, both used
unchanged.

Per grid point the two-photon window is chosen adaptively (the transparency
feature narrows and broadens by orders of magnitude across the grid), then
the fine spectrum is classified.

Outputs
  results/tables/p1_phase_diagram.csv     one row per grid point
  results/raw/p1_phase_diagram.npz     C, dAIC, class arrays on the grid
  results/tables/p1_summary.json          boundaries + validation + gates
Usage
  python p1_phase_diagram.py [--quick] [--jobs N]
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
RAW = RES / "raw"
NOGO_SRC = REPO / "No-go theorem" / "src"
sys.path.insert(0, str(NOGO_SRC))

import run_prl_prediction as rp                      # noqa: E402
import nv_model as nv                                # noqa: E402
import gate2_candidate_full_vs_reduced as g2         # noqa: E402
import gate1_candidate_aic_bootstrap as g1           # noqa: E402

# ----------------------------------------------------------------- config
BZ0 = rp.BZ0            # 0.005 T axial bias, as at the candidate
OC = rp.OC              # 0.1 GHz control Rabi, as at the candidate
J0 = rp.J0              # excited branch index of the candidate
D_STRAIN = rp.D         # 1.683 GHz

# Detection floor on |C|.  RoomT step5 fixes eps = 1.534e-7 as the contrast
# a realistic measurement can still resolve; a point whose |C| falls below it
# is reported as 'unresolved' no matter how clean its lineshape fit is.
C_DETECT = 1.534e-7

# Adaptive-window control
WIN_START_MHz = 5.0     # first attempt: +-5 MHz (the gate2 candidate window)
WIN_MAX_MHz = 4000.0
WIN_MIN_MHz = 1e-4
N_COARSE = 81
N_FINE = 161
MIN_PTS_ACROSS_FWHM = 8
FILLS_WINDOW_FRAC = 0.6     # feature this wide relative to the window -> widen
MAX_ATTEMPTS = 8


def _metrics(d_MHz, C):
    """Signed peak contrast, its position, and the FWHM of |C|."""
    i0 = int(np.argmax(np.abs(C)))
    cmax = float(C[i0])
    half = abs(cmax) / 2.0
    idx = np.where(np.abs(C) >= half)[0]
    if len(idx) > 1:
        fwhm = float(d_MHz[idx[-1]] - d_MHz[idx[0]])
    else:
        fwhm = float("nan")
    at_edge = bool(i0 <= 1 or i0 >= len(C) - 2)
    return cmax, float(d_MHz[i0]), fwhm, at_edge


def adaptive_spectrum(T, Bx, half_MHz=WIN_START_MHz, Oc=OC):
    """Settle a two-photon window that actually contains the feature.

    Widens while the peak sits on the edge, while no half-maximum crossing
    exists, or while the feature spans most of the window; narrows while the
    feature is under-sampled.  A single resolution (N_FINE) is used for every
    attempt so that the containment test and the reported metrics are derived
    from the same grid.  Returns (d2_MHz, A_full, A_cut, C, info).
    """
    info = {"attempts": 0, "widened": 0, "narrowed": 0}
    last = None
    sp = d = Af = Ac = C = None
    for attempt in range(MAX_ATTEMPTS):
        info["attempts"] = attempt + 1
        d2s = np.linspace(-half_MHz * 1e-3, half_MHz * 1e-3, N_FINE)
        sp = g2.full_spectrum(d2s, T=T, Bx=Bx, Bz=BZ0, Oc=Oc)
        d, Af, Ac, C = sp["d2_MHz"], sp["A_full"], sp["A_cut"], sp["C"]
        if not np.all(np.isfinite(C)):
            info.update(half_MHz=float(half_MHz), status="nonfinite")
            return d, Af, Ac, C, info
        cmax, ctr, fwhm, at_edge = _metrics(d, C)
        spacing = float(d[1] - d[0])
        width = float(d[-1] - d[0])

        # Not contained: peak on the edge, no half-maximum crossing at all,
        # or a feature spanning most of the window.  All three mean: widen.
        too_wide = at_edge or (not np.isfinite(fwhm)) or \
            (fwhm > FILLS_WINDOW_FRAC * width)
        too_narrow = np.isfinite(fwhm) and fwhm > 0 and \
            fwhm < MIN_PTS_ACROSS_FWHM * spacing

        if too_wide and half_MHz < WIN_MAX_MHz and last != "narrow":
            half_MHz = min(half_MHz * 6.0, WIN_MAX_MHz)
            info["widened"] += 1
            last = "widen"
            continue
        if too_narrow and half_MHz > WIN_MIN_MHz and last != "widen":
            half_MHz = max(fwhm * 4.0, WIN_MIN_MHz)
            info["narrowed"] += 1
            last = "narrow"
            continue
        info.update(half_MHz=float(half_MHz), status="ok")
        return d, Af, Ac, C, info

    # exhausted the attempt budget without settling
    info.update(half_MHz=float(half_MHz),
                status="no_feature" if half_MHz >= WIN_MAX_MHz else "unsettled")
    return d, Af, Ac, C, info


def classify_point(T, Bx, Oc=OC):
    """Sign + lineshape + EIT/ATS verdict at one (T, B_perp) grid point."""
    t0 = time.time()
    d, Af, Ac, C, info = adaptive_spectrum(T, Bx, Oc=Oc)
    row = dict(T_K=float(T), Bx_T=float(Bx), Bz_T=float(BZ0), Oc_GHz=float(Oc),
               half_MHz=float(info.get("half_MHz", np.nan)),
               attempts=int(info["attempts"]))

    row["status"] = info["status"]
    if info["status"] != "ok" or C is None or not np.all(np.isfinite(C)):
        why = {"nonfinite": "nonfinite response",
               "no_feature": "no two-photon feature up to the widest window",
               "unsettled": "window did not settle within the attempt budget"}
        cm = float(C[int(np.argmax(np.abs(C)))]) if (
            C is not None and np.all(np.isfinite(C))) else np.nan
        row.update(Cmax=cm, center_MHz=np.nan, fwhm_MHz=np.nan,
                   delta_aic=np.nan, verdict="n/a", best_model="n/a",
                   window_ok=False, at_edge=True,
                   klass="unresolved",
                   reason=why.get(info["status"], info["status"]),
                   Acut_peak=np.nan, seconds=time.time() - t0)
        return row

    cmax, ctr, fwhm, at_edge = _metrics(d, C)
    # is the feature actually contained in the settled window?
    window_ok = (not at_edge) and np.isfinite(fwhm) and \
        fwhm <= FILLS_WINDOW_FRAC * float(d[-1] - d[0])
    # C is a RATIO; where the cut absorption is itself tiny (small B_perp
    # closes the probe leg) a large |C| can accompany a negligible absolute
    # change.  Record the absolute numbers so the two are never confused --
    # this is the |C|-maximization trap flagged in the strategy document.
    ipk = int(np.argmax(np.abs(C)))
    row.update(Cmax=float(cmax), center_MHz=float(ctr), fwhm_MHz=float(fwhm),
               Acut_peak=float(np.max(np.abs(Ac))),
               Acut_at_peak=float(Ac[ipk]), Afull_at_peak=float(Af[ipk]),
               dA_at_peak=float(Ac[ipk] - Af[ipk]),
               at_edge=bool(at_edge), window_ok=bool(window_ok))

    # lineshape model selection on the control-on absorption
    try:
        fit = g1.fit_all(d, Af)
        dAIC = float(fit["delta_aic_ats_eit"])
        verdict = str(fit["verdict"])
        best = str(fit["best"])
    except Exception as exc:                                # pragma: no cover
        dAIC, verdict, best = float("nan"), "fit failed", f"error:{type(exc).__name__}"

    row.update(delta_aic=dAIC, verdict=verdict, best_model=best)

    # ---- joint decision: containment, magnitude, sign, then lineshape ----
    if not window_ok:
        klass = "unresolved"
        reason = ("feature not contained in the settled window "
                  f"(edge={at_edge}, fwhm={fwhm:.4g} MHz)")
    elif not np.isfinite(cmax) or abs(cmax) < C_DETECT:
        klass, reason = "unresolved", f"|C|<{C_DETECT:g} (below detection)"
    elif cmax < 0:
        klass, reason = "absorption", "sign reversed: control increases absorption"
    elif verdict == "robust ATS":
        klass, reason = "ATS", f"robust ATS lineshape (dAIC={dAIC:.1f})"
    elif verdict == "robust EIT":
        klass, reason = "transparency", f"robust EIT lineshape (dAIC={dAIC:.1f})"
    else:
        klass, reason = "unresolved", f"lineshape {verdict} (dAIC={dAIC:.1f})"

    row.update(klass=klass, reason=reason, seconds=time.time() - t0)
    return row


def _worker(args):
    T, Bx = args
    try:
        return classify_point(T, Bx)
    except Exception as exc:                                # pragma: no cover
        return dict(T_K=float(T), Bx_T=float(Bx), klass="unresolved",
                    reason=f"exception {type(exc).__name__}: {exc}",
                    Cmax=np.nan, delta_aic=np.nan, verdict="n/a",
                    best_model="n/a", fwhm_MHz=np.nan, center_MHz=np.nan,
                    Acut_peak=np.nan, half_MHz=np.nan, attempts=0, seconds=0.0)


def grids(quick=False):
    if quick:
        Ts = np.array([30.0, 50.0, 70.0, 90.0, 110.0])
        Bs = np.array([0.0, 0.1, rp.BX0, 0.4])
    else:
        Ts = np.arange(20.0, 115.0 + 1e-9, 5.0)                # 20 points
        Bs = np.unique(np.concatenate([np.arange(0.0, 0.5 + 1e-9, 0.05),
                                       [rp.BX0]]))             # 12 points
    return Ts, Bs


CLASS_CODE = {"transparency": 0, "ATS": 1, "absorption": 2, "unresolved": 3}


def main(quick=False, jobs=4):
    Ts, Bs = grids(quick)
    pts = [(float(T), float(B)) for T in Ts for B in Bs]
    print(f"P1: {len(Ts)} temperatures x {len(Bs)} fields = {len(pts)} "
          f"full-Liouvillian grid points, jobs={jobs}")

    t0 = time.time()
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(jobs) as pool:
            rows = []
            for i, r in enumerate(pool.imap(_worker, pts), 1):
                rows.append(r)
                if i % 10 == 0 or i == len(pts):
                    print(f"  {i}/{len(pts)}  ({time.time()-t0:.0f}s)", flush=True)
    else:
        rows = []
        for i, p in enumerate(pts, 1):
            rows.append(_worker(p))
            if i % 5 == 0 or i == len(pts):
                print(f"  {i}/{len(pts)}  ({time.time()-t0:.0f}s)", flush=True)
    elapsed = time.time() - t0

    tabdir = RES / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)

    fields = ["T_K", "Bx_T", "Bz_T", "Oc_GHz", "Cmax", "center_MHz", "fwhm_MHz",
              "Acut_peak", "Acut_at_peak", "Afull_at_peak", "dA_at_peak",
              "delta_aic", "verdict", "best_model", "klass",
              "reason", "half_MHz", "attempts", "at_edge", "window_ok", "status",
              "seconds"]
    with (tabdir / "p1_phase_diagram.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # grid arrays (T major, B minor)
    nt, nb = len(Ts), len(Bs)
    Cg = np.full((nt, nb), np.nan)
    Ag = np.full((nt, nb), np.nan)
    Fg = np.full((nt, nb), np.nan)
    Kg = np.full((nt, nb), CLASS_CODE["unresolved"], dtype=int)
    by = {(r["T_K"], r["Bx_T"]): r for r in rows}
    for i, T in enumerate(Ts):
        for j, B in enumerate(Bs):
            r = by[(float(T), float(B))]
            Cg[i, j] = r.get("Cmax", np.nan)
            Ag[i, j] = r.get("delta_aic", np.nan)
            Fg[i, j] = r.get("fwhm_MHz", np.nan)
            Kg[i, j] = CLASS_CODE.get(r.get("klass", "unresolved"), 3)
    np.savez(RAW / "p1_phase_diagram.npz", T_K=Ts, Bx_T=Bs, C=Cg,
             delta_aic=Ag, fwhm_MHz=Fg, klass=Kg,
             class_code=json.dumps(CLASS_CODE))

    # ---- boundaries: highest T that is still 'transparency', per field ----
    boundary = {}
    for j, B in enumerate(Bs):
        ok = [Ts[i] for i in range(nt) if Kg[i, j] == CLASS_CODE["transparency"]]
        boundary[f"{B:.5f}"] = (float(max(ok)) if ok else None)

    # Upper sign-reversal boundary: the lowest T above which the contrast is
    # negative at EVERY higher grid temperature.  (A plain "first negative T"
    # is the wrong statistic here: at low field the surface is non-monotonic
    # and carries a separate low-T control-induced-absorption lobe.)
    signflip = {}
    for j, B in enumerate(Bs):
        Tr = None
        for i in range(nt):
            tail = Cg[i:, j]
            fin = tail[np.isfinite(tail)]
            if len(fin) and np.all(fin < 0):
                Tr = float(Ts[i])
                break
        signflip[f"{B:.5f}"] = Tr

    # The separate low-temperature negative lobe (present only at small
    # B_perp, where the probe leg is nearly closed and |C| is a ratio of two
    # small numbers): report its extent so it is never read as the main
    # sign reversal.
    lowT_lobe = {}
    for j, B in enumerate(Bs):
        hi = None
        for i in range(nt):
            if np.isfinite(Cg[i, j]) and Cg[i, j] < 0:
                if Tr_j := signflip[f"{B:.5f}"]:
                    if Ts[i] >= Tr_j:
                        break
                hi = float(Ts[i])
            elif np.isfinite(Cg[i, j]) and Cg[i, j] > 0 and hi is not None:
                break
        lowT_lobe[f"{B:.5f}"] = hi

    counts = {k: int(np.sum(Kg == v)) for k, v in CLASS_CODE.items()}

    # ---- validation against the archived 70 K candidate (gate2 base row) ----
    val = {}
    ref = next((r for r in rows
                if abs(r["T_K"] - 70.0) < 1e-9 and abs(r["Bx_T"] - rp.BX0) < 1e-9),
               None)
    if ref is not None:
        C_REF = 0.013836          # gate2_candidate_comparison.json, config 'base'
        val = dict(candidate_C_here=float(ref["Cmax"]), candidate_C_gate2=C_REF,
                   rel_diff=float(abs(ref["Cmax"] - C_REF) / abs(C_REF)),
                   candidate_class=ref["klass"], candidate_verdict=ref["verdict"])
        val["reproduces_gate2"] = bool(val["rel_diff"] < 0.05)

    gates = dict(
        grid_complete=bool(all(np.isfinite(Cg[i, j]) or Kg[i, j] == 3
                               for i in range(nt) for j in range(nb))),
        has_transparency_region=bool(counts["transparency"] > 0),
        has_absorption_region=bool(counts["absorption"] > 0),
        candidate_reproduced=bool(val.get("reproduces_gate2", False)),
    )

    summary = dict(
        what="P1 full-Liouvillian T-B_perp phase diagram, jointly classified",
        model="9-level NV Lindblad (gate2_candidate_full_vs_reduced.build_full)",
        classifier="4-model AIC (gate1_candidate_aic_bootstrap.fit_all), "
                   "robust gate |dAIC|>=6",
        detection_floor_C=C_DETECT,
        fixed=dict(Bz_T=BZ0, Oc_GHz=OC, branch=J0, d_strain_GHz=D_STRAIN,
                   polarization="Y/Y", control_on="ms=+1"),
        T_K=[float(x) for x in Ts], Bx_T=[float(x) for x in Bs],
        class_counts=counts,
        transparency_upper_T_by_field=boundary,
        sign_reversal_T_by_field=signflip,
        lowT_negative_lobe_upper_T_by_field=lowT_lobe,
        validation=val, gates=gates, all_gates_pass=bool(all(gates.values())),
        n_points=len(pts), seconds=elapsed, quick=bool(quick),
    )
    with (tabdir / "p1_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\nclass counts: {counts}")
    print(f"transparency upper T by field: {boundary}")
    print(f"upper sign-reversal T by field: {signflip}")
    print(f"low-T negative lobe upper T by field: {lowT_lobe}")
    if val:
        print(f"candidate check: C={val['candidate_C_here']:.6g} vs gate2 "
              f"{val['candidate_C_gate2']:.6g} (rel {val['rel_diff']:.3%}) "
              f"-> {val['candidate_class']}")
    print(f"gates: {gates} -> {summary['all_gates_pass']}")
    print(f"elapsed {elapsed:.0f}s")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    main(quick=a.quick, jobs=a.jobs)
