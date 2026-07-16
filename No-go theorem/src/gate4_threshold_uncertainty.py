"""gate4_threshold_uncertainty.py -- Gate 4 of SIMULATION_PLAN.md:
temperature thresholds of the candidate EIT contrast as intervals, not
single numbers.

Monte-Carlo propagation of physical-parameter uncertainty through the
candidate contrast curve C(T) (same model as run_prl_prediction.branch_value,
re-implemented with explicit parameters so each can be varied):

  prior                                   parameter
  lognormal(x/ 1.35)                      orbital-hopping rate scale (k_orb)
  lognormal(x/ 1.2)                       radiative rate scale (GRAD)
  normal(1.683, 0.34) GHz, >0.3           transverse strain d
  fixed 0 (calibrated)                    strain azimuth phi
  normal(BX0, 1%)                         transverse field Bx
  normal(BZ0, 5%)                         axial bias Bz
  loguniform(2e-5, 2e-4) GHz              ground decoherence gg
  lognormal(x/ 1.1)                       control Rabi Oc
  normal(0, 3 deg)                        probe/control polarization angles

Thresholds extracted per sample by bisection on a T grid:
  T_1% (C=1e-2), T_0.1% (1e-3), T_0.01% (1e-4), T_sign (C=0).
Outputs: results/tables/gate4_threshold_bands.csv (quantiles),
         results/tables/gate4_sensitivity.json (Spearman ranking),
         results/figures/gate4_threshold_bands.png.
Usage: python gate4_threshold_uncertainty.py [--quick]
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import spearmanr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE))
import run_prl_prediction as rp
import nv_model as nv

TWOPI = 2*np.pi
SEED = 20260716
TGRID = np.linspace(20.0, 200.0, 61)
THRESHOLDS = {'T_1pct': 1e-2, 'T_0.1pct': 1e-3, 'T_0.01pct': 1e-4, 'T_sign': 0.0}
PARAM_NAMES = ['scale_orb', 'scale_rad', 'd_GHz', 'Bx_T', 'Bz_T',
               'gg_GHz', 'Oc_GHz', 'pol_p_deg', 'pol_c_deg']

def draw(rng):
    return dict(
        scale_orb=float(np.exp(rng.normal(0, np.log(1.35)))),
        scale_rad=float(np.exp(rng.normal(0, np.log(1.2)))),
        d_GHz=float(max(0.3, rng.normal(1.683, 0.34))),
        Bx_T=float(rng.normal(rp.BX0, 0.01*rp.BX0)),
        Bz_T=float(rng.normal(rp.BZ0, 0.05*rp.BZ0)),
        gg_GHz=float(np.exp(rng.uniform(np.log(2e-5), np.log(2e-4)))),
        Oc_GHz=float(rp.OC*np.exp(rng.normal(0, np.log(1.1)))),
        pol_p_deg=float(rng.normal(0, 3)),
        pol_c_deg=float(rng.normal(0, 3)),
    )

BASELINE = dict(scale_orb=1.0, scale_rad=1.0, d_GHz=rp.D, Bx_T=rp.BX0,
                Bz_T=rp.BZ0, gg_GHz=6.3e-5, Oc_GHz=rp.OC,
                pol_p_deg=0.0, pol_c_deg=0.0)

def contrast(T, th):
    """Candidate-branch contrast with explicit parameters (mirrors
    rp.branch_value / nv_model.response, but gamma and gg are exposed)."""
    B = (th['Bx_T'], 0.0, th['Bz_T'])
    H = nv.Hes(B, th['d_GHz'], 0.0)
    _, U = nv.dressed_ground(B)
    ap, ac = np.deg2rad(th['pol_p_deg']), np.deg2rad(th['pol_c_deg'])
    ppol = np.array([np.cos(ap), np.sin(ap)], complex)
    cpol = np.array([np.cos(ac), np.sin(ac)], complex)
    dp = np.kron(ppol, U[:, 1]); dc = np.kron(cpol, U[:, 2])
    w = np.linalg.eigvalsh(H); z = float(w[rp.J0])
    gamma = 0.5*th['scale_orb']*nv.korb_GHz(T, th['d_GHz']) \
        + 0.5*th['scale_rad']*nv.GRAD
    beta = (TWOPI*th['Oc_GHz'])**2/4
    geff = 2*th['gg_GHz'] + 2e-6
    G = np.linalg.inv(gamma*np.eye(6) + 1j*TWOPI*(H - z*np.eye(6)))
    K12 = np.vdot(dp, G@dc); K21 = np.vdot(dc, G@dp)
    S1 = np.vdot(dp, G@dp); S2 = np.vdot(dc, G@dc)
    dXi = -beta*K12*K21/(geff + beta*S2)
    Acut = float(np.real(S1))
    return float(-np.real(dXi))/Acut if abs(Acut) > 1e-300 else np.nan

def sample_thresholds(th):
    Cs = np.array([contrast(float(T), th) for T in TGRID])
    out = {}
    for name, target in THRESHOLDS.items():
        f = Cs - target
        root = np.nan
        for a, b, fa, fb in zip(TGRID[:-1], TGRID[1:], f[:-1], f[1:]):
            if fa == 0:
                root = float(a); break
            if fa*fb < 0:
                root = float(brentq(lambda t: contrast(t, th) - target,
                                    float(a), float(b)))
                break                      # first (descending) crossing
        out[name] = root
    return out, Cs

def main(quick=False):
    tabdir = ROOT/'results'/'tables'; figdir = ROOT/'results'/'figures'
    tabdir.mkdir(parents=True, exist_ok=True); figdir.mkdir(parents=True, exist_ok=True)
    n_samples = 60 if quick else 500
    rng = np.random.default_rng(SEED)

    base_th, base_C = sample_thresholds(BASELINE)
    samples, curves, roots = [], [], {k: [] for k in THRESHOLDS}
    for _ in range(n_samples):
        th = draw(rng)
        r, Cs = sample_thresholds(th)
        samples.append(th); curves.append(Cs)
        for k in THRESHOLDS: roots[k].append(r[k])
    curves = np.array(curves)

    rows = []
    for k in THRESHOLDS:
        arr = np.array(roots[k]); ok = np.isfinite(arr)
        q = (np.percentile(arr[ok], [2.5, 16, 50, 84, 97.5])
             if ok.sum() else [np.nan]*5)
        rows.append(dict(threshold=k, target=THRESHOLDS[k],
                         baseline_K=base_th[k], n_valid=int(ok.sum()),
                         n_samples=n_samples,
                         q02_5=q[0], q16=q[1], median=q[2], q84=q[3], q97_5=q[4],
                         band_68_K=q[3]-q[1], band_95_K=q[4]-q[0]))
        print(f"{k:10s} baseline={base_th[k]:.1f} K  "
              f"median={q[2]:.1f}  68%=[{q[1]:.1f},{q[3]:.1f}]  "
              f"95%=[{q[0]:.1f},{q[4]:.1f}]  valid={ok.sum()}/{n_samples}")
    with (tabdir/'gate4_threshold_bands.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # Spearman sensitivity: |rho(threshold, parameter)| ranking
    sens = {}
    X = {p: np.array([s[p] for s in samples]) for p in PARAM_NAMES}
    for k in THRESHOLDS:
        arr = np.array(roots[k]); ok = np.isfinite(arr)
        if ok.sum() < 10: continue
        rhos = {p: float(spearmanr(X[p][ok], arr[ok]).statistic)
                for p in PARAM_NAMES}
        order = sorted(rhos, key=lambda p: -abs(rhos[p]))
        sens[k] = dict(spearman=rhos, dominant=order[0],
                       ranking=order)
    with (tabdir/'gate4_sensitivity.json').open('w') as fh:
        json.dump(dict(seed=SEED, n_samples=n_samples, priors=draw.__doc__,
                       baseline=base_th, sensitivity=sens,
                       gate4_pass=dict(
                           all_thresholds_are_intervals=all(
                               np.isfinite([r['q16'], r['q84']]).all()
                               for r in rows),
                           dominant_identified=len(sens) == len(THRESHOLDS)),
                       ), fh, indent=2)

    lo, med, hi = np.nanpercentile(curves, [16, 50, 84], axis=0)
    plt.figure(figsize=(7.6, 5.2))
    plt.fill_between(TGRID, np.maximum(lo, 1e-12), np.maximum(hi, 1e-12),
                     alpha=0.3, label='16-84% band')
    plt.semilogy(TGRID, np.maximum(med, 1e-12), label='median C(T)')
    plt.semilogy(TGRID, np.maximum(base_C, 1e-12), 'k--', label='baseline')
    for k, tgt in THRESHOLDS.items():
        if tgt > 0: plt.axhline(tgt, ls=':', lw=0.7)
    for r in rows:
        if np.isfinite(r['median']):
            plt.axvspan(r['q16'], r['q84'], color='gray', alpha=0.12)
    plt.xlabel('Temperature (K)'); plt.ylabel('candidate contrast C')
    plt.title('Gate 4: threshold uncertainty bands '
              f'(n={n_samples}, seed={SEED})')
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(figdir/'gate4_threshold_bands.png', dpi=220); plt.close()
    print('dominant parameters:',
          {k: sens[k]['dominant'] for k in sens})
    return rows, sens

if __name__ == '__main__':
    main(quick='--quick' in sys.argv[1:])
