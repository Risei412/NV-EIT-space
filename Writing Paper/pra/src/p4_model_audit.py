"""p4_model_audit.py -- PRA calculation P4.

The consolidated model audit.  Its components existed in the repository but
were scattered across campaigns and never assembled, and one of them does not
pass everywhere -- which is precisely why it has to be stated in one place.

Audited here, in the order of NV_EIT_PRA_PRL_Split_Strategy_20260724.md
Sec. 5 (P4):

  A. full / reduced agreement domain
     The reduced kernel (delta_chi_S, run_prl_prediction.branch_value) is
     compared against the full 9-level Liouvillian on the SAME (T, B_perp)
     grid used by P1.  The output is an explicit domain of validity, not a
     single ratio at a single point.  This matters: RoomT step8 records a
     factor 81.5 disagreement at 10 K and an outright SIGN disagreement at
     30 K for its B_perp = 0 configuration, while gate2 records 0.1%
     agreement at the 70 K candidate.  Both are true; the audit says where.

  B. truncation dependence
     ISC + singlet (10th level) and 14N hyperfine toggles of
     gate2_candidate_full_vs_reduced, evaluated across temperature rather
     than only at the candidate.

  C. parameter provenance
     Pointers to the recorded manifests; no numbers are invented here.

  D. sign conventions
     The convention actually used by every quantity in the PRA campaign,
     checked against convention_table.md by an explicit assertion on a
     known-sign case.

  E. EIT/ATS classification stability
     The Gate-1 bootstrap / window / initial-value robustness, re-run across
     temperature and field instead of at the candidate only.

Outputs
  results/tables/p4_full_vs_reduced.csv    per grid point
  results/tables/p4_truncation.csv         per toggle and temperature
  results/tables/p4_classification_stability.csv
  results/tables/p4_audit.json             verdicts, domains, and failures
Usage
  python p4_model_audit.py [--quick] [--jobs N]
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

import run_prl_prediction as rp        # noqa: E402
import nv_model as nv                  # noqa: E402
import gate1_candidate_aic_bootstrap as g1   # noqa: E402
import gate2_candidate_full_vs_reduced as g2  # noqa: E402
import p1_phase_diagram as p1          # noqa: E402

SEED = 20260801
AGREE_TOL = 0.10        # 10% relative agreement counts as "reduced model usable"


# ------------------------------------------------------------------ A
def reduced_contrast(T, Bx):
    """Peak sector contrast from the reduced kernel at the candidate branch."""
    try:
        return float(rp.branch_value(float(T), float(Bx), p1.BZ0,
                                     rp.J0, p1.OC)["C"])
    except Exception:
        return float("nan")


def _cmp_worker(args):
    T, Bx = args
    row = p1.classify_point(float(T), float(Bx))
    Cf = row.get("Cmax", np.nan)
    Cr = reduced_contrast(T, Bx)
    okc = bool(np.isfinite(Cf) and np.isfinite(Cr) and Cr != 0)
    ratio = (Cf / Cr) if okc else float("nan")
    return dict(T_K=float(T), Bx_T=float(Bx), C_full=float(Cf),
                C_reduced=float(Cr), ratio_full_reduced=float(ratio),
                rel_err=(float(abs(Cf - Cr) / abs(Cr)) if okc else np.nan),
                sign_agree=bool(np.isfinite(Cf) and np.isfinite(Cr)
                                and np.sign(Cf) == np.sign(Cr)),
                klass_full=row.get("klass"))


# ------------------------------------------------------------------ B
def truncation_row(T, Bx, label, kw):
    """Peak contrast of the full model with one structural toggle applied."""
    d2s = np.linspace(-0.005, 0.005, 121)
    sp = g2.full_spectrum(d2s, T=float(T), Bx=float(Bx), Bz=p1.BZ0,
                          Oc=p1.OC, **kw)
    m = g2.dip_metrics(sp["d2_MHz"], sp["C"])
    return dict(T_K=float(T), Bx_T=float(Bx), toggle=label,
                Cmax=float(m["Cmax"]), fwhm_MHz=float(m["fwhm_MHz"]))


# ------------------------------------------------------------------ E
def stability_row(T, Bx, n_boot, n_init):
    """Gate-1 style robustness of the EIT/ATS verdict at one grid point."""
    d, Af, Ac, C, info = p1.adaptive_spectrum(float(T), float(Bx))
    if info["status"] != "ok":
        return dict(T_K=float(T), Bx_T=float(Bx), baseline_verdict="n/a",
                    baseline_dAIC=np.nan, boot_frac_same=np.nan,
                    init_frac_same=np.nan, window_frac_same=np.nan,
                    stable=False)
    base = g1.fit_all(d, Af)
    v0 = base["verdict"]
    rng = np.random.default_rng(SEED)

    # noise bootstrap at 2% of the feature depth
    depth = float(np.max(Af) - np.min(Af)) or 1.0
    same = 0
    for _ in range(n_boot):
        y = Af + rng.normal(0, 0.02 * depth, size=Af.shape)
        same += int(g1.fit_all(d, y)["verdict"] == v0)
    boot_frac = same / max(n_boot, 1)

    # initial-value randomization
    same = 0
    for _ in range(n_init):
        same += int(g1.fit_all(d, Af, rng=rng)["verdict"] == v0)
    init_frac = same / max(n_init, 1)

    # window variation (0.5x, 0.75x, 1.5x, 2x about the settled window)
    same, tot = 0, 0
    for fac in (0.5, 0.75, 1.5, 2.0):
        half = info["half_MHz"] * fac
        sp = g2.full_spectrum(np.linspace(-half * 1e-3, half * 1e-3, 161),
                              T=float(T), Bx=float(Bx), Bz=p1.BZ0, Oc=p1.OC)
        if not np.all(np.isfinite(sp["A_full"])):
            continue
        tot += 1
        same += int(g1.fit_all(sp["d2_MHz"], sp["A_full"])["verdict"] == v0)
    win_frac = same / tot if tot else np.nan

    stable = bool(boot_frac >= 0.9 and init_frac >= 0.9
                  and (not np.isfinite(win_frac) or win_frac >= 0.75))
    return dict(T_K=float(T), Bx_T=float(Bx), baseline_verdict=v0,
                baseline_dAIC=float(base["delta_aic_ats_eit"]),
                boot_frac_same=float(boot_frac), init_frac_same=float(init_frac),
                window_frac_same=float(win_frac), stable=stable)


def _stab_worker(a):
    T, Bx, nb, ni = a
    try:
        return stability_row(T, Bx, nb, ni)
    except Exception as exc:                                # pragma: no cover
        return dict(T_K=float(T), Bx_T=float(Bx),
                    baseline_verdict=f"error:{type(exc).__name__}",
                    baseline_dAIC=np.nan, boot_frac_same=np.nan,
                    init_frac_same=np.nan, window_frac_same=np.nan,
                    stable=False)


def _trunc_worker(a):
    T, Bx, label, kw = a
    try:
        return truncation_row(T, Bx, label, kw)
    except Exception as exc:                                # pragma: no cover
        return dict(T_K=float(T), Bx_T=float(Bx), toggle=label,
                    Cmax=np.nan, fwhm_MHz=np.nan,
                    error=f"{type(exc).__name__}: {exc}")


def main(quick=False, jobs=4):
    from multiprocessing import Pool
    tabdir = PRA / "results" / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # ---------------- A. full vs reduced over the P1 grid ----------------
    Ts = ([25., 50., 70., 90., 105.] if quick else
          [20., 25., 30., 40., 50., 60., 70., 80., 90., 100., 105., 110.])
    Bs = ([0.0, rp.BX0] if quick else
          [0.0, 0.05, 0.10, 0.15, 0.20, rp.BX0, 0.30, 0.40, 0.50])
    pairs = [(T, B) for T in Ts for B in Bs]
    print(f"P4/A full-vs-reduced on {len(pairs)} points")
    with Pool(jobs) as pool:
        cmp_rows = pool.map(_cmp_worker, pairs)
    with (tabdir / "p4_full_vs_reduced.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cmp_rows[0].keys()))
        w.writeheader()
        w.writerows(cmp_rows)

    good = [r for r in cmp_rows
            if np.isfinite(r["rel_err"]) and r["rel_err"] <= AGREE_TOL]
    signbad = [r for r in cmp_rows if np.isfinite(r["C_full"])
               and np.isfinite(r["C_reduced"]) and not r["sign_agree"]]
    worst = max((r for r in cmp_rows if np.isfinite(r["rel_err"])),
                key=lambda r: r["rel_err"], default=None)
    agree_T = sorted({r["T_K"] for r in good})
    agree_B = sorted({r["Bx_T"] for r in good})

    # ---------------- B. truncation ----------------
    toggles = [("base", {}), ("isc", dict(isc=True)),
               ("hyperfine", dict(hyperfine=True))]
    Tt = [30., 70., 100.] if quick else [25., 40., 70., 90., 105.]
    tasks = [(T, rp.BX0, lab, kw) for T in Tt for lab, kw in toggles]
    print(f"P4/B truncation on {len(tasks)} runs")
    with Pool(jobs) as pool:
        tr_rows = pool.map(_trunc_worker, tasks)
    with (tabdir / "p4_truncation.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["T_K", "Bx_T", "toggle", "Cmax",
                                           "fwhm_MHz"], extrasaction="ignore")
        w.writeheader()
        w.writerows(tr_rows)
    by = {(r["T_K"], r["toggle"]): r for r in tr_rows}
    trunc = []
    for T in Tt:
        b = by.get((T, "base"))
        if not b or not np.isfinite(b["Cmax"]):
            continue
        for lab, _ in toggles[1:]:
            r = by.get((T, lab))
            if r and np.isfinite(r["Cmax"]):
                trunc.append(dict(T_K=T, toggle=lab,
                                  ratio=float(r["Cmax"] / b["Cmax"]),
                                  sign_kept=bool(np.sign(r["Cmax"])
                                                 == np.sign(b["Cmax"]))))
    max_trunc = max((abs(t["ratio"] - 1.0) for t in trunc), default=np.nan)

    # ---------------- E. classification stability ----------------
    nb, ni = (20, 10) if quick else (60, 30)
    Tc = [30., 70., 100.] if quick else [25., 50., 70., 90., 100., 105.]
    stasks = [(T, rp.BX0, nb, ni) for T in Tc]
    print(f"P4/E classification stability on {len(stasks)} points "
          f"({nb} bootstraps, {ni} inits each)")
    with Pool(min(jobs, len(stasks))) as pool:
        st_rows = pool.map(_stab_worker, stasks)
    with (tabdir / "p4_classification_stability.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(st_rows[0].keys()))
        w.writeheader()
        w.writerows(st_rows)

    # ---------------- C, D. provenance and conventions ----------------
    prov = {
        "literature parameter manifest (SHA-256 per source)":
            "No-go theorem/results/tables/literature_manifest.csv",
        "signal-chain parameters":
            "No-go theorem/results/tables/signal_chain_parameters.csv",
        "sign and unit conventions":
            "No-go theorem/results/tables/convention_table.md",
        "phonon-rate variants":
            "No-go theorem/src/phonon_rates.py (k_orb_variants)",
        "PRA detection chain used here":
            "Writing Paper/pra/results/tables/p3_summary.json",
    }
    # explicit sign check: at the candidate a POSITIVE C must mean the control
    # REDUCES absorption, i.e. A_full < A_cut.
    d, Af, Ac, C, info = p1.adaptive_spectrum(70.0, rp.BX0)
    ipk = int(np.argmax(np.abs(C)))
    sign_convention_ok = bool(C[ipk] > 0 and Af[ipk] < Ac[ipk])

    verdicts = dict(
        full_reduced_agreement_domain_exists=bool(len(good) > 0),
        full_reduced_agrees_everywhere=bool(len(good) == len(cmp_rows)),
        sign_disagreements=len(signbad),
        truncation_within_20pct=bool(np.isfinite(max_trunc) and max_trunc < 0.2),
        classification_stable_everywhere=bool(all(r["stable"] for r in st_rows)),
        sign_convention_verified=sign_convention_ok,
    )
    audit = dict(
        what="P4 consolidated model audit for the PRA campaign",
        A_full_vs_reduced=dict(
            tolerance=AGREE_TOL, n_points=len(cmp_rows), n_within_tol=len(good),
            temperatures_with_agreement=agree_T,
            fields_with_agreement=agree_B,
            worst_rel_err=(dict(T_K=worst["T_K"], Bx_T=worst["Bx_T"],
                                rel_err=worst["rel_err"],
                                C_full=worst["C_full"],
                                C_reduced=worst["C_reduced"])
                           if worst else None),
            sign_disagreements=[dict(T_K=r["T_K"], Bx_T=r["Bx_T"],
                                     C_full=r["C_full"],
                                     C_reduced=r["C_reduced"])
                                for r in signbad],
            note="The reduced kernel may be quoted only inside the agreement "
                 "domain; every PRA number outside it is taken from the full "
                 "Liouvillian."),
        B_truncation=dict(rows=trunc, max_abs_deviation=(
            float(max_trunc) if np.isfinite(max_trunc) else None)),
        C_provenance=prov,
        D_sign_convention=dict(
            statement="C > 0 means the ground-coherence pathway REDUCES the "
                      "probe absorption (A_full < A_cut); C < 0 is "
                      "control-induced absorption.",
            verified_at="T=70 K, B_perp=candidate", ok=sign_convention_ok),
        E_classification_stability=st_rows,
        verdicts=verdicts,
        seconds=time.time() - t0, quick=bool(quick),
    )
    with (tabdir / "p4_audit.json").open("w") as fh:
        json.dump(audit, fh, indent=2, default=str)

    print(f"\nA. full vs reduced: {len(good)}/{len(cmp_rows)} points within "
          f"{AGREE_TOL:.0%}; {len(signbad)} sign disagreements")
    if worst:
        print(f"   worst: T={worst['T_K']:.0f} K B={worst['Bx_T']:.3f} T  "
              f"C_full={worst['C_full']:.4g} C_red={worst['C_reduced']:.4g} "
              f"rel_err={worst['rel_err']:.3g}")
    for r in signbad[:8]:
        print(f"   SIGN: T={r['T_K']:.0f} K B={r['Bx_T']:.3f} T  "
              f"full={r['C_full']:.4g} reduced={r['C_reduced']:.4g}")
    print(f"B. truncation max |ratio-1| = {max_trunc:.3g}")
    print("E. classification stability:")
    for r in st_rows:
        print(f"   T={r['T_K']:.0f} K  {r['baseline_verdict']:<14s} "
              f"boot={r['boot_frac_same']:.2f} init={r['init_frac_same']:.2f} "
              f"win={r['window_frac_same']:.2f} stable={r['stable']}")
    print(f"\nverdicts: {json.dumps(verdicts, indent=2)}")
    print(f"elapsed {time.time()-t0:.0f}s")
    return audit


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    main(quick=a.quick, jobs=a.jobs)
