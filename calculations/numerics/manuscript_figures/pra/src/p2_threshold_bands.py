"""p2_threshold_bands.py -- PRA calculation P2.

Uncertainty propagation on the NV EIT temperature boundary, as a function of
the transverse field.

This extends No-go theorem/src/gate4_threshold_uncertainty.py, which delivered
bands only at the single candidate field B_perp = 0.23226 T and folded the
optical linewidth implicitly into the orbital/radiative rate scales.  Three
things are added here, each demanded by the P2 pass criterion in
NV_EIT_PRA_PRL_Split_Strategy_20260724.md Sec. 5 ("do not fix 'go up to about
80 K'; present a configuration-dependent band"):

  1. B_perp is a swept variable, so the boundary is reported as a band
     *curve* T_thr(B_perp) rather than a single interval.
  2. The optical linewidth carries its OWN prior (spectral diffusion /
     residual inhomogeneous broadening, gamma_extra), instead of being
     absorbed into the orbital-hopping and radiative scale factors.
  3. The four Gamma_XY(T) phonon-rate variants of phonon_rates.k_orb_variants
     enter as a discrete model systematic, both pooled into the reported band
     and resolved per model, so the boundary cannot rest on a single
     extrapolation of the room-temperature dissipation rate.

Priors (per Monte-Carlo sample)
  lognormal(x/ 1.35)               orbital-hopping rate scale
  lognormal(x/ 1.2)                radiative rate scale
  normal(1.683, 0.34) GHz, >0.3    transverse strain d
  normal(B, max(1% B, 1 mT))       transverse field about the swept value
  normal(BZ0, 5%)                  axial bias
  loguniform(2e-5, 2e-4) GHz       ground decoherence gg
  lognormal(x/ 1.1)                control Rabi
  normal(0, 3 deg)                 probe / control polarization angles
  loguniform(1e-3, 1e-1) GHz       extra optical linewidth gamma_extra  [NEW]
  uniform over 4 variants          phonon rate model                    [NEW]

Thresholds per sample, by bisection on a temperature grid:
  T_1% (C=1e-2), T_0.1% (1e-3), T_0.01% (1e-4), T_sign (C=0).

Outputs
  results/tables/p2_threshold_bands_vs_B.csv   quantiles per (B_perp, threshold)
  results/tables/p2_phonon_model_split.csv     the same, resolved per variant
  results/tables/p2_sensitivity.json           Spearman ranking per field
  results/tables/p2_summary.json               headline bands + gates
Usage
  python p2_threshold_bands.py [--quick] [--jobs N]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
PRA = HERE.parent
REPO = PRA.parents[1]          # calculations/numerics
RES = PRA.parents[3] / "results"
NOGO_SRC = REPO / "No-go theorem" / "src"
sys.path.insert(0, str(NOGO_SRC))

import run_prl_prediction as rp        # noqa: E402
import nv_model as nv                  # noqa: E402
import phonon_rates as pr              # noqa: E402

TWOPI = 2 * np.pi
SEED = 20260801
TGRID = np.linspace(20.0, 200.0, 61)
THRESHOLDS = {"T_1pct": 1e-2, "T_0.1pct": 1e-3, "T_0.01pct": 1e-4, "T_sign": 0.0}
PHONON_MODELS = ["full_happacher", "conservative_lower_bound",
                 "saturation", "naive_T5_extrapolation"]
PARAM_NAMES = ["scale_orb", "scale_rad", "d_GHz", "Bx_T", "Bz_T", "gg_GHz",
               "Oc_GHz", "pol_p_deg", "pol_c_deg", "gamma_extra_GHz"]


def draw(rng, B_center):
    sigma_B = max(0.01 * abs(B_center), 1e-3)      # 1% or 1 mT, whichever larger
    return dict(
        scale_orb=float(np.exp(rng.normal(0, np.log(1.35)))),
        scale_rad=float(np.exp(rng.normal(0, np.log(1.2)))),
        d_GHz=float(max(0.3, rng.normal(1.683, 0.34))),
        Bx_T=float(max(0.0, rng.normal(B_center, sigma_B))),
        Bz_T=float(rng.normal(rp.BZ0, 0.05 * rp.BZ0)),
        gg_GHz=float(np.exp(rng.uniform(np.log(2e-5), np.log(2e-4)))),
        Oc_GHz=float(rp.OC * np.exp(rng.normal(0, np.log(1.1)))),
        pol_p_deg=float(rng.normal(0, 3)),
        pol_c_deg=float(rng.normal(0, 3)),
        gamma_extra_GHz=float(np.exp(rng.uniform(np.log(1e-3), np.log(1e-1)))),
        phonon_model=str(rng.choice(PHONON_MODELS)),
    )


def baseline(B_center):
    return dict(scale_orb=1.0, scale_rad=1.0, d_GHz=rp.D, Bx_T=float(B_center),
                Bz_T=rp.BZ0, gg_GHz=6.3e-5, Oc_GHz=rp.OC, pol_p_deg=0.0,
                pol_c_deg=0.0, gamma_extra_GHz=0.0,
                phonon_model="full_happacher")


def k_orb_GHz(T, d, model):
    """Orbital hopping rate in GHz for one of the four labelled variants."""
    return float(pr.k_orb_variants(float(T), float(d))[model]) * 1e-9


def contrast(T, th):
    """Candidate-branch sector contrast with every parameter exposed.

    Identical in structure to gate4_threshold_uncertainty.contrast (and hence
    to nv_model.response), with the optical-linewidth prior added to gamma and
    the phonon variant selectable.
    """
    B = (th["Bx_T"], 0.0, th["Bz_T"])
    H = nv.Hes(B, th["d_GHz"], 0.0)
    _, U = nv.dressed_ground(B)
    ap, ac = np.deg2rad(th["pol_p_deg"]), np.deg2rad(th["pol_c_deg"])
    ppol = np.array([np.cos(ap), np.sin(ap)], complex)
    cpol = np.array([np.cos(ac), np.sin(ac)], complex)
    dp = np.kron(ppol, U[:, 1])
    dc = np.kron(cpol, U[:, 2])
    w = np.linalg.eigvalsh(H)
    z = float(w[rp.J0])
    gamma = (0.5 * th["scale_orb"] * k_orb_GHz(T, th["d_GHz"], th["phonon_model"])
             + 0.5 * th["scale_rad"] * nv.GRAD
             + th["gamma_extra_GHz"])
    beta = (TWOPI * th["Oc_GHz"]) ** 2 / 4
    geff = 2 * th["gg_GHz"] + 2e-6
    G = np.linalg.inv(gamma * np.eye(6) + 1j * TWOPI * (H - z * np.eye(6)))
    K12 = np.vdot(dp, G @ dc)
    K21 = np.vdot(dc, G @ dp)
    S1 = np.vdot(dp, G @ dp)
    S2 = np.vdot(dc, G @ dc)
    dXi = -beta * K12 * K21 / (geff + beta * S2)
    Acut = float(np.real(S1))
    return float(-np.real(dXi)) / Acut if abs(Acut) > 1e-300 else np.nan


def sample_thresholds(th):
    """First (descending) crossing of each threshold, by bisection."""
    Cs = np.array([contrast(float(T), th) for T in TGRID])
    out = {}
    for name, target in THRESHOLDS.items():
        f = Cs - target
        root = np.nan
        for a, b, fa, fb in zip(TGRID[:-1], TGRID[1:], f[:-1], f[1:]):
            if not (np.isfinite(fa) and np.isfinite(fb)):
                continue
            if fa == 0:
                root = float(a)
                break
            if fa * fb < 0:
                try:
                    root = float(brentq(lambda t: contrast(t, th) - target,
                                        float(a), float(b), maxiter=200))
                except Exception:
                    root = float(0.5 * (a + b))
                break
        out[name] = root
    return out, Cs


def _quant(arr):
    ok = np.isfinite(arr)
    if not ok.sum():
        return [np.nan] * 5, 0
    return list(np.percentile(arr[ok], [2.5, 16, 50, 84, 97.5])), int(ok.sum())


def run_field(args):
    """Full Monte-Carlo at one transverse field."""
    B, n_samples, seed = args
    rng = np.random.default_rng(seed)
    base_th, base_C = sample_thresholds(baseline(B))
    samples, curves = [], []
    roots = {k: [] for k in THRESHOLDS}
    for _ in range(n_samples):
        th = draw(rng, B)
        r, Cs = sample_thresholds(th)
        samples.append(th)
        curves.append(Cs)
        for k in THRESHOLDS:
            roots[k].append(r[k])
    return dict(B=float(B), base_th=base_th, base_C=np.asarray(base_C),
                samples=samples, curves=np.asarray(curves),
                roots={k: np.asarray(v) for k, v in roots.items()})


def main(quick=False, jobs=4):
    n_samples = 60 if quick else 500
    Bs = ([0.0, 0.1, rp.BX0, 0.4] if quick else
          [0.0, 0.05, 0.10, 0.15, 0.20, rp.BX0, 0.25, 0.30, 0.40, 0.50])
    print(f"P2: {len(Bs)} fields x {n_samples} Monte-Carlo samples, jobs={jobs}")

    t0 = time.time()
    tasks = [(B, n_samples, SEED + i) for i, B in enumerate(Bs)]
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(min(jobs, len(tasks))) as pool:
            out = pool.map(run_field, tasks)
    else:
        out = [run_field(t) for t in tasks]
    elapsed = time.time() - t0

    tabdir = RES / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)

    # ---- pooled bands per (field, threshold) ----
    rows = []
    for o in out:
        for k, target in THRESHOLDS.items():
            q, nv_ = _quant(o["roots"][k])
            rows.append(dict(Bx_T=o["B"], threshold=k, target=target,
                             baseline_K=o["base_th"][k],
                             n_valid=nv_, n_samples=n_samples,
                             q02_5=q[0], q16=q[1], median=q[2], q84=q[3],
                             q97_5=q[4],
                             band_68_K=q[3] - q[1], band_95_K=q[4] - q[0]))
    with (tabdir / "p2_threshold_bands_vs_B.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- the same, resolved per phonon-rate variant ----
    split = []
    for o in out:
        mods = np.array([s["phonon_model"] for s in o["samples"]])
        for m in PHONON_MODELS:
            sel = mods == m
            for k in THRESHOLDS:
                q, nv_ = _quant(o["roots"][k][sel])
                split.append(dict(Bx_T=o["B"], phonon_model=m, threshold=k,
                                  n_valid=nv_, n_in_model=int(sel.sum()),
                                  q16=q[1], median=q[2], q84=q[3]))
    with (tabdir / "p2_phonon_model_split.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(split[0].keys()))
        w.writeheader()
        w.writerows(split)

    # ---- Spearman sensitivity per field ----
    sens = {}
    for o in out:
        X = {p: np.array([s[p] for s in o["samples"]]) for p in PARAM_NAMES}
        per = {}
        for k in THRESHOLDS:
            arr = o["roots"][k]
            ok = np.isfinite(arr)
            if ok.sum() < 10:
                continue
            rhos = {}
            for p in PARAM_NAMES:
                xs = X[p][ok]
                rhos[p] = (0.0 if np.allclose(xs, xs[0])
                           else float(spearmanr(xs, arr[ok]).statistic))
            order = sorted(rhos, key=lambda p: -abs(rhos[p]))
            per[k] = dict(spearman=rhos, dominant=order[0], ranking=order)
        sens[f"{o['B']:.5f}"] = per
    with (tabdir / "p2_sensitivity.json").open("w") as fh:
        json.dump(dict(seed=SEED, n_samples=n_samples, sensitivity=sens),
                  fh, indent=2)

    # ---- headline: T_1% and T_sign band curves vs field ----
    def curve(kname):
        return {f"{r['Bx_T']:.5f}": dict(median=r["median"], q16=r["q16"],
                                         q84=r["q84"], q02_5=r["q02_5"],
                                         q97_5=r["q97_5"], n_valid=r["n_valid"])
                for r in rows if r["threshold"] == kname}

    # spread across phonon variants (model systematic) at each field
    model_spread = {}
    for o in out:
        mods = np.array([s["phonon_model"] for s in o["samples"]])
        meds = []
        for m in PHONON_MODELS:
            a = o["roots"]["T_1pct"][mods == m]
            a = a[np.isfinite(a)]
            if len(a):
                meds.append(float(np.median(a)))
        model_spread[f"{o['B']:.5f}"] = (
            float(max(meds) - min(meds)) if len(meds) > 1 else None)

    dominants = sorted({per["T_1pct"]["dominant"]
                        for per in sens.values() if "T_1pct" in per})

    gates = dict(
        every_threshold_is_an_interval=bool(
            all(np.isfinite([r["q16"], r["q84"]]).all()
                for r in rows if r["n_valid"] > 10)),
        field_dependence_resolved=bool(len(Bs) >= 4),
        optical_linewidth_has_own_prior=True,
        phonon_model_systematic_included=True,
        dominant_identified=bool(len(dominants) > 0),
    )

    summary = dict(
        what="P2 temperature-boundary uncertainty as a function of B_perp",
        extends="No-go theorem/src/gate4_threshold_uncertainty.py",
        new_vs_gate4=["B_perp swept", "optical linewidth as its own prior",
                      "four phonon-rate variants as a model systematic"],
        n_samples=n_samples, seed=SEED, fields_T=[float(b) for b in Bs],
        thresholds=THRESHOLDS, phonon_models=PHONON_MODELS,
        T_1pct_band_vs_field=curve("T_1pct"),
        T_sign_band_vs_field=curve("T_sign"),
        phonon_model_median_spread_T_1pct_K=model_spread,
        dominant_parameters_T_1pct=dominants,
        gates=gates, all_gates_pass=bool(all(gates.values())),
        seconds=elapsed, quick=bool(quick),
    )
    with (tabdir / "p2_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    print("\nT_1% band vs field (68%):")
    for b, v in summary["T_1pct_band_vs_field"].items():
        print(f"  B={b} T : median {v['median']:.1f} K  "
              f"[{v['q16']:.1f}, {v['q84']:.1f}]  (n_valid={v['n_valid']})")
    print("\nT_sign band vs field (68%):")
    for b, v in summary["T_sign_band_vs_field"].items():
        m = v["median"]
        print(f"  B={b} T : median {m:.1f} K  "
              f"[{v['q16']:.1f}, {v['q84']:.1f}]  (n_valid={v['n_valid']})"
              if np.isfinite(m) else f"  B={b} T : no crossing")
    print(f"\nphonon-model median spread on T_1%: {model_spread}")
    print(f"dominant parameters for T_1%: {dominants}")
    print(f"gates: {gates} -> {summary['all_gates_pass']}")
    print(f"elapsed {elapsed:.0f}s")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    main(quick=a.quick, jobs=a.jobs)
