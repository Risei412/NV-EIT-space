"""p7_control_power_edges.py -- PRA calculation P7.

Control-power dependence of the two temperature edges of the NV EIT island.

P1 established that at the working field the island is bounded in temperature
on both sides, and that crossing either edge takes the coherence-mediated
response negative: the control field stops inducing transparency and starts
inducing absorption.  The two edges are attributed to different mechanisms --
an Autler-Townes crossover at the cold edge, phonon-driven collapse of the
leading Raman path at the warm edge -- but P1 fixes the control Rabi frequency
at the candidate value (Oc = 0.1 GHz) and so cannot separate them.

This calculation sweeps Oc and locates both edges at each power.  The
attribution above predicts

    the COLD edge moves with Oc,          because it is set by Oc vs gamma_oc,
    the WARM edge does not,               because it is set by the phonon rate,

which is a statement an experiment can test in a single sample without any
absolute calibration: only the *displacement* of each edge under a change of
control power is needed, not the contrast magnitude.  Neither the conventional
dephasing account (which predicts no edge at all, only a monotonic fade) nor a
pure Autler-Townes account (which cannot produce the warm edge) reproduces the
pair.

The located quantity is the SIGN edge: the temperature at which sign(C)
changes.  Its laboratory reading is which way the fringe points, and it is
the observable the room-temperature note and the manuscript both rest on.
The EIT/ATS model-selection verdict is recorded at every evaluated point but
is NOT bisected: at some (Oc, T) the adaptive two-photon window fails to
settle -- e.g. Oc = 0.03 GHz at 70 K widens to the 4 GHz ceiling and returns
'inconclusive' -- so the verdict is contaminated by a window artifact that
the sign is not.  Verdict boundaries are therefore reported as brackets at
whatever resolution the sign bisection happened to sample, and flagged as
such.

Everything comes from p1_phase_diagram.classify_point, used unchanged, so
every number here is directly comparable with the published phase diagram.

CAUTION carried over from P1: at small B_perp the probe leg is nearly closed,
A_cut is tiny, and |C| is a ratio of two small numbers.  At the working field
used here (B_perp = 0.232 T) it is not: P1 gives A_cut = 0.76 and dA = -1.17
at the cold edge, both O(1).  This script stays at the working field and
records A_cut and dA at every point so the reader can check that the ratio
never carries the conclusion.

Outputs
  results/tables/p7_control_power_map.csv   one row per (Oc, T) evaluation
  results/tables/p7_edges_vs_control.csv    one row per Oc: the located edges
  results/tables/p7_summary.json            slopes, gates, registered prediction
Usage
  python p7_control_power_edges.py [--quick] [--jobs N]
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

import run_prl_prediction as rp                      # noqa: E402
import nv_model as nv                                # noqa: E402
import p1_phase_diagram as p1                        # noqa: E402

# ----------------------------------------------------------------- config
BX0 = rp.BX0            # 0.23226 T working field, as at the candidate
D_STRAIN = rp.D         # 1.683 GHz
OC0 = rp.OC             # 0.1 GHz, the P1 / candidate control power

# Control-power grid: geometric, ~1.4 decades, containing the working point.
OC_GRID = [0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.22, 0.32, 0.50]
OC_GRID_QUICK = [0.03, 0.10, 0.32]

# Fixed wide brackets.  T_LO must be below the cold edge and T_HI above the
# warm edge at EVERY control power, and the island interior anchor must be
# positive; all three are checked, not assumed, and a failure is reported as
# "not bracketed" rather than silently bisected.
T_LO = 5.0              # below the cold edge: expect C < 0
T_MID = 60.0            # island interior, cold side: expect C > 0
T_ANCHOR = 70.0         # island interior, the candidate temperature
T_HI = 130.0            # above the warm edge: expect C < 0

TOL_K = 0.5             # bisection tolerance on every located edge
MAX_ROUNDS = 9          # log2((T_MID-T_LO)/TOL) ~ 6.8; 9 is slack


# ------------------------------------------------------------------ probes
def _positive(row):
    """Sign predicate: is the peak sector contrast positive here?"""
    c = row.get("Cmax", np.nan)
    if c is None or not np.isfinite(c):
        return None
    return bool(c > 0.0)


def _transparency(row):
    """Verdict predicate: joint class is genuine transparency."""
    k = row.get("klass")
    if k is None:
        return None
    return bool(k == "transparency")


def _worker(args):
    """One full-Liouvillian classification.  (Oc, T) -> annotated p1 row."""
    Oc, T = args
    try:
        row = p1.classify_point(float(T), BX0, Oc=float(Oc))
    except Exception as exc:                                # pragma: no cover
        row = dict(T_K=float(T), Bx_T=BX0, Cmax=np.nan,
                   klass="unresolved", verdict="n/a",
                   reason=f"exception {type(exc).__name__}: {exc}")
    row["Oc_GHz"] = float(Oc)
    row["gamma_oc_GHz"] = float(nv.gamma_oc_GHz(float(T), D_STRAIN))
    row["gamma_over_Oc"] = row["gamma_oc_GHz"] / float(Oc)
    row["positive"] = _positive(row)
    row["transparency"] = _transparency(row)
    return row


def _key(Oc, T):
    return (round(float(Oc), 9), round(float(T), 6))


def _run(pts, jobs, label, cache):
    """Evaluate the (Oc, T) points not already cached."""
    todo = [p for p in pts if _key(*p) not in cache]
    if not todo:
        return
    t0 = time.time()
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(jobs) as pool:
            out = list(pool.imap(_worker, todo))
    else:
        out = [_worker(p) for p in todo]
    for p, r in zip(todo, out):
        cache[_key(*p)] = r
    print(f"  {label}: {len(todo)} points in {time.time()-t0:.0f}s", flush=True)


# ----------------------------------------------------------------- edges
def _bisect_round(cache, jobs, tasks, rnd):
    """One lockstep bisection round over every open edge; tasks mutated."""
    mids = [(t["Oc"], 0.5 * sum(t["bracket"])) for t in tasks if not t["done"]]
    if not mids:
        return False
    _run(mids, jobs, f"bisect round {rnd} ({len(mids)} open)", cache)
    for t in tasks:
        if t["done"]:
            continue
        lo, hi = t["bracket"]          # (T where False, T where True)
        mid = 0.5 * (lo + hi)
        v = _positive(cache[_key(t["Oc"], mid)])
        if v is None:
            t["done"], t["status"] = True, "undefined at midpoint"
            continue
        t["bracket"] = (lo, mid) if v else (mid, hi)
        if abs(t["bracket"][1] - t["bracket"][0]) <= TOL_K:
            t["done"], t["status"] = True, "converged"
    return True


def main(quick=False, jobs=4):
    oc_grid = OC_GRID_QUICK if quick else OC_GRID
    tol_note = TOL_K
    ends = (T_LO, T_MID, T_ANCHOR, T_HI)

    print(f"P7: {len(oc_grid)} control powers at B_perp = {BX0:.5f} T, "
          f"brackets [{T_LO},{T_MID}] and [{T_ANCHOR},{T_HI}] K, "
          f"tol {tol_note} K, jobs={jobs}")

    t0 = time.time()
    cache = {}
    _run([(o, T) for o in oc_grid for T in ends], jobs, "endpoints", cache)

    # ---- validate every bracket before bisecting into it -----------------
    tasks, brackets = [], {}
    for o in oc_grid:
        p_lo = _positive(cache[_key(o, T_LO)])
        p_mid = _positive(cache[_key(o, T_MID)])
        p_anc = _positive(cache[_key(o, T_ANCHOR)])
        p_hi = _positive(cache[_key(o, T_HI)])
        b = dict(positive_at=dict(T_LO=p_lo, T_MID=p_mid,
                                  T_ANCHOR=p_anc, T_HI=p_hi))
        if p_lo is False and p_mid is True:
            tasks.append(dict(Oc=float(o), side="cold",
                              bracket=(T_LO, T_MID), done=False, status="open"))
            b["cold"] = "bracketed"
        else:
            b["cold"] = f"not bracketed (C>0 at {T_LO} K: {p_lo}, " \
                        f"at {T_MID} K: {p_mid})"
        if p_hi is False and p_anc is True:
            tasks.append(dict(Oc=float(o), side="warm",
                              bracket=(T_HI, T_ANCHOR), done=False,
                              status="open"))
            b["warm"] = "bracketed"
        else:
            b["warm"] = f"not bracketed (C>0 at {T_ANCHOR} K: {p_anc}, " \
                        f"at {T_HI} K: {p_hi})"
        brackets[float(o)] = b

    print(f"bisecting {len(tasks)} edges")
    rnd = 0
    while rnd < MAX_ROUNDS and _bisect_round(cache, jobs, tasks, rnd + 1):
        rnd += 1
    for t in tasks:
        lo, hi = t["bracket"]
        t["T_K"] = 0.5 * (lo + hi)
        t["half_width_K"] = 0.5 * abs(hi - lo)
    elapsed = time.time() - t0

    tabdir = RES / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)

    # ---- full evaluation map -------------------------------------------
    fields = ["Oc_GHz", "T_K", "Bx_T", "Bz_T", "Cmax", "positive",
              "transparency", "center_MHz", "fwhm_MHz", "Acut_at_peak",
              "dA_at_peak", "gamma_oc_GHz", "gamma_over_Oc", "delta_aic",
              "verdict", "best_model", "klass", "reason", "half_MHz",
              "window_ok", "status", "seconds"]
    allrows = sorted(cache.values(), key=lambda r: (r["Oc_GHz"], r["T_K"]))
    with (tabdir / "p7_control_power_map.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in allrows:
            w.writerow(r)

    # ---- verdict boundaries, read off the sampled points ----------------
    # Not bisected (see module docstring): reported as the bracket between
    # the sampled temperatures where klass enters and leaves 'transparency'.
    def verdict_brackets(o):
        pts = sorted(((r["T_K"], _transparency(r)) for r in allrows
                      if abs(r["Oc_GHz"] - o) < 1e-9), key=lambda x: x[0])
        pts = [(t, v) for t, v in pts if v is not None]
        lo = hi = None
        for (ta, va), (tb, vb) in zip(pts, pts[1:]):
            if (not va) and vb and lo is None:
                lo = (ta, tb)
            if va and (not vb):
                hi = (ta, tb)
        return lo, hi

    # ---- one row per control power -------------------------------------
    by = {}
    for t in tasks:
        by.setdefault(t["Oc"], {})[t["side"]] = t
    erows = []
    for o in oc_grid:
        d = by.get(float(o), {})
        vlo, vhi = verdict_brackets(float(o))
        row = dict(Oc_GHz=float(o), Bx_T=BX0,
                   klass_at_70K=cache[_key(o, T_ANCHOR)].get("klass"),
                   Cmax_at_70K=cache[_key(o, T_ANCHOR)].get("Cmax"))
        for side in ("cold", "warm"):
            t = d.get(side)
            row[f"T_{side}_sign_K"] = round(t["T_K"], 3) if t else None
            row[f"T_{side}_sign_halfwidth_K"] = (round(t["half_width_K"], 3)
                                                 if t else None)
            row[f"{side}_status"] = (t["status"] if t
                                     else brackets[float(o)][side])
        row["verdict_enter_transparency_K"] = (f"{vlo[0]:.2f}-{vlo[1]:.2f}"
                                               if vlo else None)
        row["verdict_leave_transparency_K"] = (f"{vhi[0]:.2f}-{vhi[1]:.2f}"
                                               if vhi else None)
        tc = d.get("cold")
        if tc:
            g = float(nv.gamma_oc_GHz(tc["T_K"], D_STRAIN))
            row["gamma_oc_at_cold_edge_GHz"] = round(g, 6)
            row["gamma_over_Oc_at_cold_edge"] = round(g / float(o), 4)
        else:
            row["gamma_oc_at_cold_edge_GHz"] = None
            row["gamma_over_Oc_at_cold_edge"] = None
        tw = d.get("warm")
        row["gamma_oc_at_warm_edge_GHz"] = (
            round(float(nv.gamma_oc_GHz(tw["T_K"], D_STRAIN)), 6)
            if tw else None)
        erows.append(row)

    with (tabdir / "p7_edges_vs_control.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(erows[0].keys()))
        w.writeheader()
        for r in erows:
            w.writerow(r)

    # ---- slopes and the registered prediction ---------------------------
    def _slope(side):
        oc = np.array([r["Oc_GHz"] for r in erows
                       if r[f"T_{side}_sign_K"] is not None], float)
        T = np.array([r[f"T_{side}_sign_K"] for r in erows
                      if r[f"T_{side}_sign_K"] is not None], float)
        if len(T) < 2:
            return dict(n=int(len(T)), slope_K_per_decade=None, span_K=None,
                        T_min=None, T_max=None)
        return dict(n=int(len(T)),
                    slope_K_per_decade=float(np.polyfit(np.log10(oc), T, 1)[0]),
                    span_K=float(T.max() - T.min()),
                    T_min=float(T.min()), T_max=float(T.max()))

    cold, warm = _slope("cold"), _slope("warm")

    # gamma_oc / Oc at the cold edge: a constant means the cold edge is fixed
    # by the Autler-Townes condition and by nothing else.
    ratios = [r["gamma_over_Oc_at_cold_edge"] for r in erows
              if r["gamma_over_Oc_at_cold_edge"] is not None]
    m = float(np.mean(ratios)) if ratios else None
    ats_scaling = dict(n=len(ratios), values=ratios, mean=m,
                       rel_spread=(float((max(ratios) - min(ratios)) / m)
                                   if ratios and m else None))

    ratio = None
    if cold["span_K"] is not None and warm["span_K"] is not None:
        ratio = (cold["span_K"] / warm["span_K"]) if warm["span_K"] > 0 \
            else float("inf")

    gates = dict(
        cold_edge_bracketed=bool(cold["n"] >= 2),
        warm_edge_bracketed=bool(warm["n"] >= 2),
        cold_edge_moves=bool(cold["span_K"] is not None
                             and cold["span_K"] >= 5.0),
        warm_edge_static=bool(warm["span_K"] is not None
                              and warm["span_K"] <= 2.0),
        edges_separable=bool(ratio is not None and ratio >= 3.0),
        island_interior_positive=bool(
            all(_positive(cache[_key(o, T_ANCHOR)]) for o in oc_grid)),
    )

    summary = dict(
        what="P7 control-power dependence of the two temperature edges of the "
             "NV EIT island, at the working transverse field",
        model="9-level NV Lindblad via p1_phase_diagram.classify_point "
              "(unchanged classifier, adaptive window, detection floor)",
        registered_prediction=(
            "The cold edge is an Autler-Townes crossover and moves with the "
            "control Rabi frequency; the warm edge is set by the "
            "phonon-induced orbital hopping rate and does not. An experiment "
            "that sweeps control power at fixed field therefore displaces one "
            "edge and not the other, which no single-mechanism account "
            "reproduces."),
        located_quantity="sign edge: temperature at which sign(C) changes",
        verdict_note=("EIT/ATS verdict boundaries are reported as brackets "
                      "between sampled temperatures, not bisected: the "
                      "adaptive two-photon window fails to settle at some "
                      "(Oc, T) and contaminates the verdict, but not the "
                      "sign."),
        fixed=dict(Bx_T=BX0, Bz_T=p1.BZ0, branch=p1.J0,
                   d_strain_GHz=D_STRAIN, polarization="Y/Y",
                   control_on="ms=+1", detection_floor_C=p1.C_DETECT),
        Oc_GHz=[float(o) for o in oc_grid], Oc_reference_GHz=OC0,
        brackets=dict(T_LO=T_LO, T_MID=T_MID, T_ANCHOR=T_ANCHOR, T_HI=T_HI),
        bracket_validation=brackets,
        bisection_tol_K=TOL_K, rounds=rnd,
        edges=erows, cold_edge=cold, warm_edge=warm,
        cold_over_warm_span_ratio=ratio,
        ats_scaling_gamma_oc_over_Oc_at_cold_edge=ats_scaling,
        caution=("At small B_perp the probe leg is nearly closed and |C| is a "
                 "ratio of two small numbers; this run stays at the working "
                 "field, where P1 gives A_cut = O(1) at the cold edge. A_cut "
                 "and dA are recorded per point in p7_control_power_map.csv."),
        gates=gates, all_gates_pass=bool(all(gates.values())),
        n_evaluations=len(cache), seconds=elapsed, quick=bool(quick),
    )
    with (tabdir / "p7_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    # ---- console report -------------------------------------------------
    print("\n  Oc[GHz]   T_cold[K]    T_warm[K]   gamma_oc/Oc@cold   "
          "enter-transp    leave-transp")
    for r in erows:
        def f(k):
            v = r[f"T_{k}_sign_K"]
            return f"{v:9.2f}" if v is not None else "        -"
        g = r["gamma_over_Oc_at_cold_edge"]
        print(f"  {r['Oc_GHz']:7.3f} {f('cold')} {f('warm')}   "
              f"{(f'{g:.4f}' if g is not None else '-'):>14}   "
              f"{str(r['verdict_enter_transparency_K']):>13} "
              f"{str(r['verdict_leave_transparency_K']):>15}")
    print(f"\ncold edge: span {cold['span_K']} K, "
          f"slope {cold['slope_K_per_decade']} K/decade")
    print(f"warm edge: span {warm['span_K']} K, "
          f"slope {warm['slope_K_per_decade']} K/decade")
    print(f"span ratio cold/warm: {ratio}")
    print(f"gamma_oc/Oc at cold edge: {ats_scaling['values']}")
    print(f"   mean {ats_scaling['mean']}, rel spread {ats_scaling['rel_spread']}")
    print(f"gates: {gates} -> {summary['all_gates_pass']}")
    print(f"{len(cache)} evaluations, elapsed {elapsed:.0f}s")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    main(quick=a.quick, jobs=a.jobs)
