"""gate1_candidate_aic_bootstrap.py -- Gate 1 of SIMULATION_PLAN.md.

EIT/ATS/Lorentzian/Fano model comparison at the fixed PRL candidate point
(T=70 K, Bx=BX0, Bz=BZ0, branch J0, Oc=0.1 GHz), with bootstrap noise
resampling, frequency-window variation, fit-initial-value randomization,
and a control-power sweep mapping the EIT->ATS crossover.

Spectrum source: the SAME machinery as run_prl_prediction.py (model_at +
nv_model conventions), but with the full z-dependence of the resolvent
kept while scanning the two-photon detuning: z = z0 + d2 and the
ground-coherence denominator extended as den -> geff + i*2*pi*d2 + beta*S2(z).
A_on(d2) = Re S1(z) + Re dXi(z)  (Acut - dA, convention_table.md).

Models (fit over d2 in MHz; all include an explicit constant background b,
which deviates from bare Anisimov-Dowling-Sanders PRL 107,163604 because at
the candidate the one-photon line, gamma ~ 72 GHz, is flat across the
MHz-wide window and would otherwise force gp to a degenerate large value):
  M1 EIT :  b + Cp^2/(gp^2+d^2) - Cm^2/(gm^2+d^2)                 (k=5)
  M2 ATS :  b + C^2 [1/(g^2+(d-d0)^2) + 1/(g^2+(d+d0)^2)]         (k=4)
  M3 Lor :  b + a/(g^2+d^2)                                       (k=3)
  M4 Fano:  b + a (q g + d)^2 / ((1+q^2)(g^2+d^2))                (k=4)
AIC = N ln(RSS/N) + 2k (repo convention, eit_ats_classifier.py); AICc and
BIC = N ln(RSS/N) + k ln N alongside; fixed robust gate on
Delta_AIC = IC_ATS - IC_EIT: >= +6 robust EIT, <= -6 robust ATS.

Outputs (results/ relative to package root):
  tables/gate1_aic_bootstrap.json
  tables/gate1_window_sweep.csv
  tables/gate1_oc_sweep.csv
  figures/gate1_model_selection_stability.png
Usage: python gate1_candidate_aic_bootstrap.py [--quick]
"""
from __future__ import annotations
import csv, json, sys, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE))
import run_prl_prediction as rp
import nv_model as nv

TWOPI = 2*np.pi
GG = 6.3e-5
GEFF = 2*GG + 2e-6          # identical to nv_model.response
ROBUST = 6.0                # |Delta_AIC| gate (SiV_SnV_phonon_AIC_parameters.md 10.2)
SEED = 20260716

# ---------------------------------------------------------------- spectrum
def probe_scan(T=70.0, Bx=rp.BX0, Bz=rp.BZ0, j=rp.J0, Oc=rp.OC,
               half_GHz=0.02, n=1601):
    """Control-on absorption vs two-photon detuning d2 (probe scanned,
    control fixed), full z-dependence retained."""
    H, dp, dc, w, V, pp, cc = rp.model_at(Bx, Bz)
    z0 = float(w[j]); gamma = nv.gamma_oc_GHz(T, rp.D)
    beta = (TWOPI*Oc)**2/4
    d2s = np.linspace(-half_GHz, half_GHz, n)
    A = np.empty(n); Acut = np.empty(n)
    I6 = np.eye(6)
    for i, d2 in enumerate(d2s):
        z = z0 + d2
        G = np.linalg.inv(gamma*I6 + 1j*TWOPI*(H - z*I6))
        K12 = np.vdot(dp, G@dc); K21 = np.vdot(dc, G@dp)
        S1 = np.vdot(dp, G@dp); S2 = np.vdot(dc, G@dc)
        den = GEFF + 1j*TWOPI*d2 + beta*S2
        dXi = -beta*K12*K21/den
        Acut[i] = float(np.real(S1))
        A[i] = Acut[i] + float(np.real(dXi))       # Acut - dA
    return d2s*1e3, A, Acut                        # detuning in MHz

# ---------------------------------------------------------------- models
def m_eit(d, b, Cp, gp, Cm, gm):
    return b + Cp**2/(gp**2 + d**2) - Cm**2/(gm**2 + d**2)

def m_ats(d, b, C, g, d0):
    return b + C**2*(1.0/(g**2 + (d - d0)**2) + 1.0/(g**2 + (d + d0)**2))

def m_lor(d, b, a, g):
    return b + a/(g**2 + d**2)

def m_fano(d, b, a, q, g):
    return b + a*(q*g + d)**2/((1.0 + q**2)*(g**2 + d**2))

MODELS = dict(EIT=(m_eit, 5), ATS=(m_ats, 4), Lorentzian=(m_lor, 3), Fano=(m_fano, 4))

def _feature(d, A):
    """Background, extremum depth and FWHM of the narrow feature."""
    edge = int(max(4, 0.05*len(d)))
    b = float(np.median(np.concatenate([A[:edge], A[-edge:]])))
    dev = A - b
    i0 = int(np.argmax(np.abs(dev)))
    depth = float(dev[i0])                        # signed (dip < 0)
    half = b + depth/2
    inside = np.abs(dev) >= abs(depth)/2
    idx = np.where(inside)[0]
    fwhm = float(d[idx[-1]] - d[idx[0]]) if len(idx) > 1 else float(d[1]-d[0])
    return dict(b=b, center=float(d[i0]), depth=depth, fwhm=max(fwhm, float(d[1]-d[0])))

def _p0(name, f):
    b, dep, g = f['b'], f['depth'], max(f['fwhm']/2, 1e-3)
    amp = abs(dep)*g**2
    if name == 'EIT':
        # broad positive minus narrow negative reproduces a dip on background
        return [b, np.sqrt(amp*4), 20*g, np.sqrt(max(amp, 1e-300)), g]
    if name == 'ATS':
        return [b, np.sqrt(max(amp, 1e-300)), g, 2*g]
    if name == 'Lorentzian':
        return [b, np.sign(dep)*amp if dep else -amp, g]
    if name == 'Fano':
        return [b, np.sign(dep)*abs(dep) if dep else -abs(dep), 0.1, g]
    raise KeyError(name)

_BOUNDS = dict(
    EIT=([-np.inf, 0, 1e-4, 0, 1e-4], [np.inf]*5),
    ATS=([-np.inf, 0, 1e-4, 0], [np.inf]*4),
    Lorentzian=([-np.inf, -np.inf, 1e-4], [np.inf]*3),
    Fano=([-np.inf, -np.inf, -np.inf, 1e-4], [np.inf]*4),
)

def _ics(rss, n, k):
    rss = max(rss, 1e-300)
    aic = n*np.log(rss/n) + 2*k
    aicc = aic + (2*k*(k+1)/(n-k-1) if n-k-1 > 0 else np.inf)
    bic = n*np.log(rss/n) + k*np.log(n)
    return aic, aicc, bic

def fit_all(d, A, rng=None, maxfev=20000):
    """Fit all four models; return per-model ICs plus the fixed-gate verdict
    on Delta_AIC(ATS-EIT). Optional rng jitters the initial values."""
    scale = float(np.max(np.abs(A))) or 1.0
    y = A/scale
    f = _feature(d, y)
    out = {}
    n = len(d)
    for name, (fn, k) in MODELS.items():
        p0 = np.array(_p0(name, f), float)
        if rng is not None:
            jit = np.exp(rng.uniform(np.log(1/3), np.log(3), size=p0.size))
            p0 = p0*jit
        lo, hi = _BOUNDS[name]
        p0 = np.clip(p0, [l if np.isfinite(l) else -1e30 for l in lo],
                         [h if np.isfinite(h) else 1e30 for h in hi])
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', OptimizeWarning)
                p, _ = curve_fit(fn, d, y, p0=p0, bounds=(lo, hi), maxfev=maxfev)
            rss = float(np.sum((y - fn(d, *p))**2))
        except Exception:
            p, rss = p0, np.inf
        aic, aicc, bic = _ics(rss, n, k)
        out[name] = dict(rss=rss, k=k, aic=aic, aicc=aicc, bic=bic,
                         params=[float(v) for v in p])
    use_aicc = (n/max(m[1] for m in MODELS.values()) < 40)
    key = 'aicc' if use_aicc else 'aic'
    vals = np.array([out[m][key] for m in MODELS])
    w = np.exp(-(vals - vals.min())/2); w /= w.sum()
    for i, m in enumerate(MODELS):
        out[m]['akaike_weight'] = float(w[i])
    dAIC = out['ATS'][key] - out['EIT'][key]
    verdict = ('robust EIT' if dAIC >= ROBUST else
               'robust ATS' if dAIC <= -ROBUST else 'inconclusive')
    return dict(models=out, delta_aic_ats_eit=float(dAIC), verdict=verdict,
                used_aicc=bool(use_aicc), scale=scale,
                best=str(min(MODELS, key=lambda m: out[m][key])))

# ---------------------------------------------------------------- stages
def stage_baseline(d, A):
    r = fit_all(d, A)
    f = _feature(d, A)
    r['feature'] = f
    return r

def stage_noise_bootstrap(d, A, n_boot, sigma_rels, seed=SEED):
    f = _feature(d, A)
    depth = abs(f['depth'])
    rows = []
    for s in sigma_rels:
        rng = np.random.default_rng(seed + int(1e6*s))
        sigma = s*depth
        counts = dict(); daics = []
        for _ in range(n_boot):
            r = fit_all(d, A + rng.normal(0, sigma, size=A.size), maxfev=5000)
            counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
            daics.append(r['delta_aic_ats_eit'])
        daics = np.array(daics)
        rows.append(dict(sigma_rel=float(s), sigma_abs=float(sigma), n_boot=n_boot,
                         frac_robust_eit=counts.get('robust EIT', 0)/n_boot,
                         frac_robust_ats=counts.get('robust ATS', 0)/n_boot,
                         frac_inconclusive=counts.get('inconclusive', 0)/n_boot,
                         delta_aic_median=float(np.median(daics)),
                         delta_aic_p05=float(np.percentile(daics, 5)),
                         delta_aic_p95=float(np.percentile(daics, 95)),
                         all_daic=[float(x) for x in daics]))
    return rows

def stage_window_sweep(base_half, factors, n_boot, sigma_rel, seed=SEED, **scan_kw):
    rows = []
    for fac in factors:
        d, A, _ = probe_scan(half_GHz=base_half*fac, **scan_kw)
        r0 = fit_all(d, A)
        f = _feature(d, A); depth = abs(f['depth'])
        rng = np.random.default_rng(seed + int(1000*fac))
        n_eit = 0
        for _ in range(n_boot):
            rb = fit_all(d, A + rng.normal(0, sigma_rel*depth, size=A.size), maxfev=5000)
            n_eit += (rb['verdict'] == 'robust EIT')
        rows.append(dict(window_factor=float(fac), half_MHz=float(base_half*fac*1e3),
                         delta_aic=r0['delta_aic_ats_eit'], verdict=r0['verdict'],
                         best_model=r0['best'],
                         frac_robust_eit_boot=n_eit/n_boot))
    return rows

def stage_init_robustness(d, A, n_draws, seed=SEED):
    r0 = fit_all(d, A)
    rng = np.random.default_rng(seed)
    daics, verdicts = [], []
    for _ in range(n_draws):
        r = fit_all(d, A, rng=rng)
        daics.append(r['delta_aic_ats_eit']); verdicts.append(r['verdict'])
    daics = np.array(daics)
    return dict(n_draws=n_draws, baseline_delta_aic=r0['delta_aic_ats_eit'],
                delta_aic_min=float(daics.min()), delta_aic_max=float(daics.max()),
                delta_aic_spread=float(daics.max() - daics.min()),
                all_verdicts_equal=bool(len(set(verdicts)) == 1),
                verdict_counts={v: verdicts.count(v) for v in set(verdicts)})

def stage_oc_sweep(ocs, seed=SEED, n=1201):
    """Adaptive-window classification vs control power: coarse pass to
    locate the feature width, refined pass for the fit."""
    H, dp, dc, w, V, pp, cc = rp.model_at(rp.BX0, rp.BZ0)
    z0 = float(w[rp.J0]); gamma = nv.gamma_oc_GHz(70.0, rp.D)
    I6 = np.eye(6)
    G0 = np.linalg.inv(gamma*I6 + 1j*TWOPI*(H - z0*I6))
    S2_0 = np.vdot(dc, G0@dc)
    rows = []
    for Oc in ocs:
        beta = (TWOPI*Oc)**2/4
        guess = max(0.02, 12*(GEFF + beta*abs(np.real(S2_0)))/TWOPI)
        d, A, _ = probe_scan(Oc=float(Oc), half_GHz=guess, n=801)
        f = _feature(d, A)
        half = max(0.02, 6*f['fwhm']*1e-3)
        d, A, _ = probe_scan(Oc=float(Oc), half_GHz=half, n=n)
        r = fit_all(d, A)
        f = _feature(d, A)
        rows.append(dict(Oc_GHz=float(Oc), half_MHz=float(half*1e3),
                         fwhm_MHz=f['fwhm'], depth=f['depth'],
                         delta_aic=r['delta_aic_ats_eit'], verdict=r['verdict'],
                         best_model=r['best'],
                         w_eit=r['models']['EIT']['akaike_weight'],
                         w_ats=r['models']['ATS']['akaike_weight'],
                         w_lor=r['models']['Lorentzian']['akaike_weight'],
                         w_fano=r['models']['Fano']['akaike_weight']))
    return rows

# ---------------------------------------------------------------- io/figure
def write_csv(path, rows):
    if not rows: return
    with path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

def figure(figpath, d, A, base, boot_rows, win_rows, oc_rows):
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    a = ax[0, 0]
    sc = base['scale']
    a.plot(d, A, 'k.', ms=2, label='model spectrum')
    pe = base['models']['EIT']['params']; pa = base['models']['ATS']['params']
    a.plot(d, sc*m_eit(d, *pe), '-', label='EIT fit')
    a.plot(d, sc*m_ats(d, *pa), '--', label='ATS fit')
    a.set_xlabel('two-photon detuning (MHz)'); a.set_ylabel('absorption (arb.)')
    a.set_title(f"candidate 70 K: {base['verdict']} "
                f"(ΔAIC={base['delta_aic_ats_eit']:.1f})")
    a.legend(fontsize=8)
    a = ax[0, 1]
    for row in boot_rows:
        a.hist(row['all_daic'], bins=30, alpha=0.5,
               label=f"σ={row['sigma_rel']:.2f}·depth")
    a.axvline(ROBUST, ls=':', c='k'); a.axvline(-ROBUST, ls=':', c='k')
    a.set_xlabel('ΔAIC (ATS−EIT)'); a.set_ylabel('bootstrap count')
    a.set_title('noise bootstrap'); a.legend(fontsize=8)
    a = ax[1, 0]
    xf = [r['window_factor'] for r in win_rows]
    a.plot(xf, [r['delta_aic'] for r in win_rows], 'o-')
    a.axhline(ROBUST, ls=':', c='k'); a.axhline(-ROBUST, ls=':', c='k')
    a.set_xlabel('window factor (× ±20 MHz)'); a.set_ylabel('ΔAIC (ATS−EIT)')
    a.set_title('frequency-window robustness')
    a = ax[1, 1]
    oc = [r['Oc_GHz'] for r in oc_rows]
    a.semilogx(oc, [r['delta_aic'] for r in oc_rows], 'o-')
    a.axhline(ROBUST, ls=':', c='k'); a.axhline(-ROBUST, ls=':', c='k')
    a.axvline(rp.OC, ls='--', c='gray', label='candidate Ωc')
    a.set_xlabel('Ω$_c$ (GHz)'); a.set_ylabel('ΔAIC (ATS−EIT)')
    a.set_title('EIT→ATS crossover'); a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figpath, dpi=220); plt.close(fig)

def main(quick=False):
    tabdir = ROOT/'results'/'tables'; figdir = ROOT/'results'/'figures'
    tabdir.mkdir(parents=True, exist_ok=True); figdir.mkdir(parents=True, exist_ok=True)
    n_boot = 40 if quick else 200
    n_init = 10 if quick else 50
    n_winboot = 10 if quick else 50
    ocs = np.logspace(np.log10(0.02), np.log10(20.0), 9 if quick else 25)

    d, A, Acut = probe_scan()
    base = stage_baseline(d, A)
    print(f"baseline: verdict={base['verdict']} dAIC={base['delta_aic_ats_eit']:.2f} "
          f"best={base['best']} weights=" +
          ' '.join(f"{m}:{base['models'][m]['akaike_weight']:.3f}" for m in MODELS))

    boot = stage_noise_bootstrap(d, A, n_boot, (0.02, 0.05, 0.1, 0.2))
    for r in boot:
        print(f"bootstrap σ={r['sigma_rel']:.2f}: EIT {r['frac_robust_eit']:.2%} "
              f"ATS {r['frac_robust_ats']:.2%} inconcl {r['frac_inconclusive']:.2%}")

    win = stage_window_sweep(0.02, (0.5, 0.75, 1.0, 1.5, 2.0), n_winboot, 0.05)
    init = stage_init_robustness(d, A, n_init)
    print(f"init robustness: spread={init['delta_aic_spread']:.3f} "
          f"all_equal={init['all_verdicts_equal']}")
    oc_rows = stage_oc_sweep(ocs)

    # EIT->ATS crossover must be physically ordered: verdict rank
    # (EIT=1, inconclusive=0, ATS=-1) non-increasing with control power,
    # and the candidate Oc must sit inside the EIT region.
    rank = {'robust EIT': 1, 'inconclusive': 0, 'robust ATS': -1}
    seq = [rank[r['verdict']] for r in oc_rows]
    cand_rows = [r for r in oc_rows if r['Oc_GHz'] <= rp.OC]
    passed = dict(
        noise_95pct_eit=all(r['frac_robust_eit'] >= 0.95 for r in boot
                            if r['sigma_rel'] <= 0.05),
        window_no_flip=len({r['verdict'] for r in win}) == 1,
        init_stable=init['all_verdicts_equal'],
        crossover_continuous=bool(np.all(np.diff(seq) <= 0)),
        candidate_in_eit_region=bool(cand_rows and
                                     all(r['verdict'] == 'robust EIT'
                                         for r in cand_rows)),
    )
    summary = dict(
        candidate=dict(T_K=70.0, Bx_T=rp.BX0, Bz_T=rp.BZ0, branch=rp.J0,
                       Oc_GHz=rp.OC, seed=SEED, quick=quick),
        baseline=base, noise_bootstrap=[{k: v for k, v in r.items()
                                         if k != 'all_daic'} for r in boot],
        window_sweep=win, init_robustness=init,
        gate1_pass=passed, gate1_all_passed=all(passed.values()),
    )
    with (tabdir/'gate1_aic_bootstrap.json').open('w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2, default=str)
    write_csv(tabdir/'gate1_window_sweep.csv', win)
    write_csv(tabdir/'gate1_oc_sweep.csv', oc_rows)
    figure(figdir/'gate1_model_selection_stability.png', d, A, base, boot, win, oc_rows)
    print('gate1_pass:', passed, '->', all(passed.values()))
    return summary

if __name__ == '__main__':
    main(quick='--quick' in sys.argv[1:])
