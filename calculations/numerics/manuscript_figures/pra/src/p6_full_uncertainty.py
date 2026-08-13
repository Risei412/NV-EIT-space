"""p6_full_uncertainty.py -- full-Liouvillian check on the P2 bands.

P4 shows that the reduced kernel does NOT reproduce the full Liouvillian above
about 90 K at any field (0/6 fields within 10% at 100, 105 and 110 K).  P2
propagated its uncertainty through the reduced kernel, so its sign-reversal
band -- which sits at ~100 K -- is quoted from a model that is not valid
there.  The band therefore has to be re-derived with the full model before it
can enter the paper.

Doing the whole of P2 with the full Liouvillian is not affordable (500 samples
x 61 temperatures x a spectrum scan each).  What is affordable, and what is
actually required, is a check at the field the paper quotes: a smaller
Monte-Carlo with the SAME priors, evaluated with the full 9-level Liouvillian
on a temperature grid that brackets both thresholds.  If the full-model band
overlaps the reduced-model band the P2 numbers stand (with the validity
statement); if it does not, the full-model band is the one to publish.

`build_full_p` below is a parameter-exposed copy of
gate2_candidate_full_vs_reduced.build_full -- the same construction, with the
orbital/radiative scales, the ground decoherence, an extra optical dephasing
and the phonon-rate variant lifted into arguments.  This mirrors what
gate4_threshold_uncertainty.contrast does for the reduced model, whose
docstring likewise records that it "mirrors rp.branch_value / nv_model.response,
but gamma and gg are exposed".

Outputs
  results/tables/p6_full_uncertainty.csv   per-sample thresholds
  results/tables/p6_summary.json           full vs reduced band comparison
Usage
  python p6_full_uncertainty.py [--quick] [--jobs N] [--samples N]
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
import phonon_rates as pr              # noqa: E402
import gate2_candidate_full_vs_reduced as g2   # noqa: E402
from liouvillian_core import liouvillian, steady_state, first_order  # noqa: E402
import p2_threshold_bands as p2        # noqa: E402

TWOPI = 2 * np.pi
SEED = 20260801
T1_GROUND = 1e-6
OP = 1e-5
THRESHOLDS = {"T_1pct": 1e-2, "T_sign": 0.0}


def build_full_p(T, Bx, Bz, d2, th):
    """gate2.build_full with the uncertainty parameters exposed.

    Differences from the original, all of them parameters rather than
    structure: the orbital hopping rate carries `scale_orb` and a selectable
    phonon variant, the radiative rate carries `scale_rad`, the ground
    dephasing is `gg` rather than the module constant, and an extra optical
    dephasing `gamma_extra` is applied to the excited manifold (the
    spectral-diffusion / residual-inhomogeneity prior of P2).
    """
    N = 9
    Bvec = (float(Bx), 0.0, float(Bz))
    Hg0 = nv.Hgs(Bvec)
    He0 = nv.Hes(Bvec, th["d_GHz"], 0.0)
    eg0, U = g2.dressed_from(Hg0)
    z0 = float(np.linalg.eigvalsh(He0)[rp.J0])

    p_idx, c_idx = 1, 2
    s_idx = 3 - p_idx - c_idx
    ap, ac = np.deg2rad(th["pol_p_deg"]), np.deg2rad(th["pol_c_deg"])
    ppol = np.array([np.cos(ap), np.sin(ap)], complex)
    cpol = np.array([np.cos(ac), np.sin(ac)], complex)
    dp = np.kron(ppol, U[:, p_idx])
    dc = np.kron(cpol, U[:, c_idx])

    H = np.zeros((N, N), complex)
    H[3:9, 3:9] = He0 - (z0 + d2) * np.eye(6)
    H[c_idx, c_idx] = -d2
    H[s_idx, s_idx] = float(eg0[s_idx] - eg0[p_idx])
    Vc = np.zeros((N, N), complex); Vc[3:9, c_idx] = dc
    Vp = np.zeros((N, N), complex); Vp[3:9, p_idx] = dp
    H += 0.5 * th["Oc_GHz"] * (Vc + Vc.conj().T)

    rate = th["scale_orb"] * float(
        pr.k_orb_variants(float(T), th["d_GHz"])[th["phonon_model"]]) * 1e-9
    grad = th["scale_rad"] * nv.GRAD
    Ls = []
    for m in range(3):                       # orbital hopping X<->Y
        up = np.zeros((N, N), complex); dn = np.zeros((N, N), complex)
        up[6 + m, 3 + m] = 1; dn[3 + m, 6 + m] = 1
        Ls += [np.sqrt(rate) * up, np.sqrt(rate) * dn]
    for orb0 in (3, 6):                      # radiative, spin conserving
        for m in range(3):
            J = np.zeros((N, N), complex)
            for a in range(3):
                J[a, orb0 + m] = np.conj(U[m, a])
            Ls.append(np.sqrt(grad) * J)
    for a in range(3):                       # ground T1
        for b in range(3):
            if a != b:
                J = np.zeros((N, N), complex); J[b, a] = 1
                Ls.append(np.sqrt(T1_GROUND) * J)
    for a in range(3):                       # ground dephasing
        J = np.zeros((N, N), complex); J[a, a] = 1
        Ls.append(np.sqrt(2 * th["gg_GHz"]) * J)
    ge = th["gamma_extra_GHz"]
    if ge > 0:                               # extra optical dephasing
        for j in range(3, 9):
            J = np.zeros((N, N), complex); J[j, j] = 1
            Ls.append(np.sqrt(2 * ge) * J)
    return TWOPI * H, Ls, Vp, dp, dict(N=N, p_idx=p_idx, c_idx=c_idx)


def chi_pair_p(T, Bx, Bz, d2, th):
    H, Ls, Vp, dp, meta = build_full_p(T, Bx, Bz, d2, th)
    N, p_idx, c_idx = meta["N"], meta["p_idx"], meta["c_idx"]
    L = liouvillian(H, Ls)
    rho0, _, _ = steady_state(L)
    Hp = TWOPI * 0.5 * OP * (Vp + Vp.conj().T)
    I = np.eye(N)
    V = -1j * (np.kron(I, Hp) - np.kron(Hp.T, I))
    det = np.zeros(N * N, complex)
    for e, a in enumerate(dp):
        det[p_idx * N + (3 + e)] = np.conj(a)
    S = [c_idx * N + p_idx, p_idx * N + c_idx]
    X = [k for k in range(N * N) if k not in S]
    Lc = L.copy(); Lc[np.ix_(S, X)] = 0; Lc[np.ix_(X, S)] = 0
    xf, _ = first_order(L, V, rho0)
    xc, _ = first_order(Lc, V, rho0)
    return complex(-2 * (det.conj() @ xf) / OP), complex(-2 * (det.conj() @ xc) / OP)


def peak_contrast(T, Bx, th, half_MHz=5.0, n=61):
    """Signed peak sector contrast, two-stage fixed-schedule window."""
    for stage in range(2):
        d2s = np.linspace(-half_MHz * 1e-3, half_MHz * 1e-3, n)
        Af = np.empty(n); Ac = np.empty(n)
        for i, d2 in enumerate(d2s):
            f, c = chi_pair_p(T, Bx, p2.rp.BZ0, float(d2), th)
            Af[i] = f.imag; Ac[i] = c.imag
        if not np.all(np.isfinite(Ac)) or np.any(np.abs(Ac) < 1e-300):
            return np.nan
        C = (Ac - Af) / Ac
        i0 = int(np.argmax(np.abs(C)))
        if stage == 0:
            halfmax = abs(C[i0]) / 2
            idx = np.where(np.abs(C) >= halfmax)[0]
            fwhm = (d2s[idx[-1]] - d2s[idx[0]]) * 1e3 if len(idx) > 1 else np.nan
            if np.isfinite(fwhm) and fwhm > 0:
                half_MHz = float(np.clip(4 * fwhm, 0.02, 60.0))
            else:
                half_MHz = min(half_MHz * 6, 60.0)
            continue
        return float(C[i0])
    return np.nan


TGRID = np.array([50., 60., 70., 80., 85., 90., 95., 100., 105., 110., 120., 135.])


def sample_thresholds_full(th, Bx):
    Cs = np.array([peak_contrast(float(T), Bx, th) for T in TGRID])
    out = {}
    for name, target in THRESHOLDS.items():
        f = Cs - target
        root = np.nan
        for a, b, fa, fb in zip(TGRID[:-1], TGRID[1:], f[:-1], f[1:]):
            if not (np.isfinite(fa) and np.isfinite(fb)):
                continue
            if fa * fb < 0:                      # linear interpolation in T
                root = float(a + (b - a) * fa / (fa - fb))
                break
        out[name] = root
    return out, Cs


def _worker(a):
    idx, Bx, seed = a
    rng = np.random.default_rng(seed)
    th = p2.draw(rng, Bx)
    try:
        r, Cs = sample_thresholds_full(th, Bx)
    except Exception as exc:                                # pragma: no cover
        return dict(sample=idx, Bx_T=Bx, error=f"{type(exc).__name__}: {exc}",
                    **{k: np.nan for k in THRESHOLDS})
    return dict(sample=idx, Bx_T=float(Bx),
                **{k: float(v) for k, v in r.items()},
                scale_orb=th["scale_orb"], d_GHz=th["d_GHz"],
                gg_GHz=th["gg_GHz"], Oc_GHz=th["Oc_GHz"],
                gamma_extra_GHz=th["gamma_extra_GHz"],
                phonon_model=th["phonon_model"])


def _q(a):
    a = np.asarray(a, float); ok = np.isfinite(a)
    if not ok.sum():
        return dict(median=None, q16=None, q84=None, n=0)
    return dict(median=float(np.percentile(a[ok], 50)),
                q16=float(np.percentile(a[ok], 16)),
                q84=float(np.percentile(a[ok], 84)), n=int(ok.sum()))


def main(quick=False, jobs=4, samples=None):
    Bx = rp.BX0
    n = samples if samples else (12 if quick else 80)
    print(f"P6: {n} full-Liouvillian Monte-Carlo samples at B_perp={Bx:.5f} T, "
          f"{len(TGRID)} temperatures each, jobs={jobs}")
    tasks = [(i, Bx, SEED + 1000 + i) for i in range(n)]
    t0 = time.time()
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(jobs) as pool:
            rows = []
            for i, r in enumerate(pool.imap(_worker, tasks), 1):
                rows.append(r)
                if i % 5 == 0 or i == n:
                    print(f"  {i}/{n}  ({time.time()-t0:.0f}s)", flush=True)
    else:
        rows = [_worker(t) for t in tasks]
    elapsed = time.time() - t0

    tabdir = RES / "tables"
    with (tabdir / "p6_full_uncertainty.csv").open("w", newline="") as fh:
        keys = ["sample", "Bx_T", "T_1pct", "T_sign", "scale_orb", "d_GHz",
                "gg_GHz", "Oc_GHz", "gamma_extra_GHz", "phonon_model"]
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    full = {k: _q([r.get(k, np.nan) for r in rows]) for k in THRESHOLDS}

    # reduced-model band at the same field, from P2
    red = {}
    try:
        p2sum = json.loads((tabdir / "p2_summary.json").read_text())
        key = f"{Bx:.5f}"
        red["T_1pct"] = p2sum["T_1pct_band_vs_field"].get(key)
        red["T_sign"] = p2sum["T_sign_band_vs_field"].get(key)
    except Exception:
        red = {}

    def overlaps(a, b):
        if not a or not b or a["median"] is None or b.get("median") is None:
            return None
        return bool(a["q16"] <= b["q84"] and b["q16"] <= a["q84"])

    comparison = {}
    for k in THRESHOLDS:
        comparison[k] = dict(
            full=full[k], reduced=red.get(k),
            bands_overlap=overlaps(full[k], red.get(k)),
            median_shift_K=(None if not red.get(k) or full[k]["median"] is None
                            else float(full[k]["median"] - red[k]["median"])))

    gates = dict(
        both_thresholds_resolved=bool(all(full[k]["n"] > 0.5 * n
                                          for k in THRESHOLDS)),
        sign_reversal_band_overlaps_reduced=bool(
            comparison["T_sign"]["bands_overlap"]),
        one_pct_band_overlaps_reduced=bool(
            comparison["T_1pct"]["bands_overlap"]),
    )
    summary = dict(
        what="P6 full-Liouvillian uncertainty check on the P2 reduced-model bands",
        why="P4 finds the reduced kernel invalid above ~90 K, which is where "
            "the P2 sign-reversal band sits",
        model="parameter-exposed 9-level Liouvillian (build_full_p)",
        priors="identical to P2 (p2_threshold_bands.draw)",
        B_perp_T=float(Bx), n_samples=n, temperature_grid=[float(t) for t in TGRID],
        comparison=comparison, gates=gates,
        all_gates_pass=bool(all(v for v in gates.values() if v is not None)),
        seconds=elapsed, quick=bool(quick),
    )
    with (tabdir / "p6_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    for k in THRESHOLDS:
        c = comparison[k]
        f, r = c["full"], c["reduced"]
        fs = (f"{f['median']:.1f} [{f['q16']:.1f}, {f['q84']:.1f}] (n={f['n']})"
              if f["median"] is not None else "unresolved")
        rs = (f"{r['median']:.1f} [{r['q16']:.1f}, {r['q84']:.1f}]"
              if r and r.get("median") is not None else "n/a")
        print(f"\n{k}:  full {fs}")
        print(f"{' '*len(k)}   reduced {rs}   overlap={c['bands_overlap']}  "
              f"shift={c['median_shift_K']}")
    print(f"\ngates: {gates}")
    print(f"elapsed {elapsed:.0f}s")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--samples", type=int, default=None)
    a = ap.parse_args()
    main(quick=a.quick, jobs=a.jobs, samples=a.samples)
