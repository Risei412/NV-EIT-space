"""p5_bperp_scaling.py -- transverse-field opening of the Raman pathway.

At B_perp = 0 the NV spin-Lambda leading Raman vertex is closed: the optical
dipole is orbital, so the two legs of the Lambda reach orthogonal spin states
and the leading path moment vanishes.  A transverse field mixes the spin
states and opens the pathway.  Because it enters once in each leg, the
sector-mediated contrast is expected to open QUADRATICALLY,

    C  ~  B_perp^2 ,

and it is that exponent -- not the mere existence of a feature -- that makes
the transverse-field channel a falsifiable prediction rather than a fitted
knob.

The archived campaign established this at 300 K (reduced kernel exponent
2.026; full Liouvillian 2.178, `bperp_full_vs_reduced_slopes.json`).  Here the
exponent is measured with the FULL Liouvillian at the temperatures the PRA
actually claims, in the small-field regime where the expansion applies, and
its breakdown at larger field is located.

Outputs
  results/tables/p5_bperp_scaling.csv    C(B_perp) per temperature
  results/tables/p5_summary.json         fitted exponents + gates
Usage
  python p5_bperp_scaling.py [--quick] [--jobs N]
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
REPO = PRA.parents[1]
NOGO_SRC = REPO / "No-go theorem" / "src"
sys.path.insert(0, str(NOGO_SRC))
sys.path.insert(0, str(HERE))

import p1_phase_diagram as p1        # noqa: E402

# Small-field window: large enough that the opened pathway dominates the
# residual B=0 background, small enough that the quadratic term still leads.
B_SMALL = np.round(np.logspace(np.log10(1.5e-3), np.log10(0.06), 12), 6)
B_WIDE = np.unique(np.concatenate([
    B_SMALL, np.array([0.09, 0.13, 0.18, 0.25, 0.35, 0.5])]))

# The B_perp = 0 response is NOT zero: spin-orbit and spin-spin terms inside
# the 3E manifold leave a small residual Raman path that survives with no
# transverse field at all.  That residual is a DIFFERENT pathway from the one
# the transverse field opens, and at the bottom of the field window the two
# are comparable, so fitting the raw C(B) mixes them and biases the exponent
# upward.  The field-opened part is C(B) - C(0), and that is what must show
# the quadratic law.  A point is used in the fit only when the opened part
# stands clear of the residual by this factor.
SATURATION_FRACTION = 0.30   # fit only below this fraction of the plateau
MIN_FIT_POINTS = 5


def _worker(a):
    T, B = a
    try:
        r = p1.classify_point(float(T), float(B))
        return dict(T_K=float(T), Bx_T=float(B), C=float(r.get("Cmax", np.nan)),
                    dA=float(r.get("dA_at_peak", np.nan)),
                    klass=r.get("klass"), verdict=r.get("verdict"))
    except Exception as exc:                                # pragma: no cover
        return dict(T_K=float(T), Bx_T=float(B), C=np.nan, dA=np.nan,
                    klass=f"error:{type(exc).__name__}", verdict="n/a")


def fit_power_law(B, C):
    """Fit C(B) = C_res + a * B**n over the pre-saturation range.

    Subtracting the zero-field value and fitting a pure power law does not
    work here: the residual pathway is not the value at exactly B = 0 (the
    zero-field point is structurally degenerate and does not join smoothly
    onto the small-field plateau), and at the temperatures where the island
    exists the plateau and the onset of saturation are less than a decade
    apart, so a two-point log slope sweeps continuously through 2 instead of
    resting on it.  Fitting the residual as a free parameter, over the range
    where the response has not yet saturated, is the well-posed version of
    the same question.
    """
    from scipy.optimize import curve_fit
    B = np.asarray(B, float); C = np.asarray(C, float)
    ok = np.isfinite(B) & np.isfinite(C) & (B > 0)
    B, C = B[ok], C[ok]
    if len(B) < MIN_FIT_POINTS:
        return dict(n=np.nan, n_err=np.nan, C_res=np.nan, r2=np.nan,
                    n_points=int(len(B)), fit_range_T=None)
    C_sat = float(C[np.argmax(B)])
    keep = np.abs(C) < abs(SATURATION_FRACTION * C_sat) if C_sat else \
        np.ones_like(C, bool)
    if keep.sum() < MIN_FIT_POINTS:                    # widen if too strict
        keep = np.argsort(B) < max(MIN_FIT_POINTS, len(B) // 2)
    x, y = B[keep], C[keep]

    def model(b, c_res, a, n):
        return c_res + a * b ** n

    a0 = (y.max() - y.min()) / (x.max() ** 2 - x.min() ** 2 + 1e-30)
    try:
        popt, pcov = curve_fit(model, x, y, p0=(float(y.min()), float(a0), 2.0),
                               maxfev=40000)
        perr = np.sqrt(np.diag(pcov))
        yhat = model(x, *popt)
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        return dict(n=float(popt[2]), n_err=float(perr[2]),
                    C_res=float(popt[0]), a=float(popt[1]), r2=float(r2),
                    n_points=int(keep.sum()),
                    fit_range_T=[float(x.min()), float(x.max())])
    except Exception:
        return dict(n=np.nan, n_err=np.nan, C_res=np.nan, r2=np.nan,
                    n_points=int(keep.sum()), fit_range_T=None)


def main(quick=False, jobs=4):
    # 40 K is excluded: there the zero-field residual is itself large and
    # negative (the low-temperature control-induced-absorption lobe reaches
    # down to zero field), so "the field-opened part" is not separable.
    Ts = [70.0] if quick else [55.0, 70.0, 85.0]
    Bs = B_SMALL if quick else B_WIDE
    pts = [(T, B) for T in Ts for B in Bs] + [(T, 0.0) for T in Ts]
    print(f"P5: {len(Ts)} temperatures x {len(Bs)} fields = {len(pts)} points")

    t0 = time.time()
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(jobs) as pool:
            rows = pool.map(_worker, pts)
    else:
        rows = [_worker(p) for p in pts]
    elapsed = time.time() - t0

    base0 = {r["T_K"]: r["C"] for r in rows if r["Bx_T"] == 0.0}
    rows = [r for r in rows if r["Bx_T"] > 0.0]

    tabdir = PRA / "results" / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)
    with (tabdir / "p5_bperp_scaling.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["T_K", "Bx_T", "C", "dA", "klass",
                                           "verdict"])
        w.writeheader()
        w.writerows(rows)

    fits = {}
    for T in Ts:
        sub = sorted([r for r in rows if r["T_K"] == T], key=lambda r: r["Bx_T"])
        B = np.array([r["Bx_T"] for r in sub])
        C = np.array([r["C"] for r in sub])
        f = fit_power_law(B, C)
        f["C_zero_field_point"] = float(base0.get(T, np.nan))
        f["C_saturated"] = float(C[np.argmax(B)])
        # Window robustness: the three-parameter fit is only meaningful if it
        # survives a change of upper cutoff.  Where the plateau and the onset
        # of saturation are less than a decade apart there is not enough
        # dynamic range and the exponent is not determined; that has to be
        # reported rather than hidden by a single lucky window.
        wr = {}
        for cut in (0.0157, 0.0307, 0.0605):
            m = B <= cut
            if m.sum() >= MIN_FIT_POINTS:
                g = fit_power_law(B[m], C[m])
                if np.isfinite(g["n"]) and np.isfinite(g.get("n_err", np.nan)) \
                        and g["r2"] > 0.99:
                    wr[f"{cut:g}"] = dict(n=g["n"], n_err=g["n_err"], r2=g["r2"])
        f["n_by_window"] = wr
        ns = [v["n"] for v in wr.values()]
        f["window_spread"] = float(max(ns) - min(ns)) if len(ns) > 1 else None
        f["window_robust"] = bool(len(ns) > 1 and (max(ns) - min(ns)) < 0.5)
        fits[f"{T:.0f}"] = f
        rs = (f"[{f['fit_range_T'][0]:.4g},{f['fit_range_T'][1]:.4g}] T"
              if f["fit_range_T"] else "no usable range")
        wr = ", ".join(f"{k}T:{v['n']:.2f}" for k, v in f["n_by_window"].items())
        print(f"  T={T:5.1f} K : n = {f['n']:.3f} +- {f['n_err']:.3f}  "
              f"(R^2={f['r2']:.4f}, {f['n_points']} pts over {rs})  "
              f"C_sat={f['C_saturated']:.4g}")
        print(f"           window check [{wr}]  spread={f['window_spread']}  "
              f"robust={f['window_robust']}")

    exps = [v["n"] for v in fits.values() if np.isfinite(v["n"])]
    errs = [v["n_err"] for v in fits.values() if np.isfinite(v.get("n_err", np.nan))]
    robust = {k: v for k, v in fits.items() if v.get("window_robust")}
    gates = dict(
        fits_are_clean=bool(all(v["r2"] > 0.98 for v in fits.values()
                                if np.isfinite(v["r2"]))),
        opening_is_superlinear=bool(exps and all(e > 1.5 for e in exps)),
        # The quadratic law is claimed ONLY where the exponent survives a
        # change of fitting window; elsewhere the dynamic range between the
        # residual plateau and saturation is too small to determine it.
        quadratic_where_window_robust=bool(
            robust and all(abs(v["n"] - 2.0) <= max(3 * v["n_err"], 0.3)
                           for v in robust.values())),
        window_robust_somewhere=bool(len(robust) > 0),
    )
    summary = dict(
        what="P5 transverse-field opening exponent of the sector contrast",
        model="full 9-level Liouvillian (P1 machinery)",
        expectation="C = C_res + a*B_perp^2: one symmetry-breaking insertion "
                    "per Lambda leg, on top of a field-independent residual "
                    "pathway through 3E spin mixing",
        method="three-parameter fit C = C_res + a*B^n below saturation",
        archived_reference=dict(
            source="No-go theorem/results/tables/bperp_full_vs_reduced_slopes.json",
            reduced_300K=2.1615, full_300K=2.1780),
        saturation_fraction=SATURATION_FRACTION,
        fits=fits,
        caveat="At the temperatures where the transparency island exists the "
               "field range between the residual plateau and saturation is "
               "less than a decade, so the exponent is window-dependent and "
               "the asymptotic quadratic law is NOT verified there; it is "
               "quoted only where the window check passes.",
        window_robust_temperatures=sorted(robust),
        gates=gates, all_gates_pass=bool(all(gates.values())),
        seconds=elapsed, quick=bool(quick),
    )
    with (tabdir / "p5_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"gates: {gates} -> {summary['all_gates_pass']}")
    print(f"elapsed {elapsed:.0f}s")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    main(quick=a.quick, jobs=a.jobs)
