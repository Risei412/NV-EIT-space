"""gate5_ensemble_average.py -- Gate 5 of SIMULATION_PLAN.md: does the
70 K candidate window survive ensemble averaging?

Fixed lasers (the candidate probe/control frequencies) are applied to an
inhomogeneous ensemble; each member's spectrum is evaluated with the
reduced-kernel machinery generalized to off-candidate members:
  z_probe(d2) = z0 + d2 + [eps_p(theta) - eps_p(cand)] - delta_opt
  delta_2(theta) = d2 - ([eps_c-eps_p](theta) - [eps_c-eps_p](cand))
(eps_a = dressed ground energies; delta_opt = member's optical-frequency
offset, which cancels in the two-photon detuning). The ensemble
absorption is the equal-weight average and
  C_ens(d2) = (<A_off> - <A_on>)/<A_off>.

Distributions: the four tetrahedral NV orientations under the fixed lab
field (Bz0 along a1, Bx0 along a perpendicular axis), static-field spread,
optical-transition (strain/spectral-diffusion) offsets, transverse-strain
spread, and control-intensity spread. Members of foreign orientations are
Raman-detuned by ~GHz, so they contribute (exactly, via the same formula)
background absorption that dilutes the contrast -- the main washout channel
alongside the ground-splitting spread from field inhomogeneity.

Scenarios:
  single         orientation-selected single defect (reference)
  low_density    4 orientations, sigma_B=0.1 G, sigma_opt=5 GHz, sigma_d=0.1
  high_density   4 orientations, sigma_B=1 G,   sigma_opt=30 GHz, sigma_d=0.3
  post_selected  high-density widths, but orientation 1 only and
                 |delta_opt| < 5 GHz (spectral selection)

Outputs: results/tables/gate5_ensemble_contrast.csv,
         results/figures/gate5_ensemble_spectra.png,
         results/tables/gate5_summary.json
Usage: python gate5_ensemble_average.py [--quick]
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np

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
GEFF = 2*GG + 2e-6
SEED = 20260716
GAUSS_T = 1e-4                     # 1 gauss in tesla

AXES = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
                float)/np.sqrt(3)
XLAB = np.array([1, -1, 0], float)/np.sqrt(2)   # perpendicular to AXES[0]
B_LAB = rp.BZ0*AXES[0] + rp.BX0*XLAB

def frame_field(k, dB_lab):
    """(Bx, Bz) of the lab field (+ vector noise, tesla) in NV frame k."""
    B = B_LAB + dB_lab
    Bz = float(B@AXES[k])
    Bx = float(np.linalg.norm(B - Bz*AXES[k]))
    return Bx, Bz

def member_context(Bx, Bz, d):
    B = (Bx, 0.0, Bz)
    H = nv.Hes(B, d, 0.0)
    eg, U = nv.dressed_ground(B)
    dp = np.kron(rp.Y, U[:, 1]); dc = np.kron(rp.Y, U[:, 2])
    return H, dp, dc, float(eg[1]), float(eg[2])

# candidate reference (member the lasers are tuned to)
_H0, _dp0, _dc0, EPS_P0, EPS_C0 = member_context(rp.BX0, rp.BZ0, rp.D)
Z0 = float(np.linalg.eigvalsh(_H0)[rp.J0])

def member_spectra(d2s, Bx, Bz, d, Oc, delta_opt, T=70.0):
    """(A_on, A_off) arrays for one ensemble member under the fixed lasers."""
    H, dp, dc, eps_p, eps_c = member_context(Bx, Bz, d)
    gamma = nv.gamma_oc_GHz(T, d)
    beta = (TWOPI*Oc)**2/4
    dg = (eps_c - eps_p) - (EPS_C0 - EPS_P0)
    zoff = (eps_p - EPS_P0) - delta_opt
    I6 = np.eye(6)
    A_on = np.empty(len(d2s)); A_off = np.empty(len(d2s))
    for i, d2 in enumerate(d2s):
        z = Z0 + d2 + zoff
        G = np.linalg.inv(gamma*I6 + 1j*TWOPI*(H - z*I6))
        K12 = np.vdot(dp, G@dc); K21 = np.vdot(dc, G@dp)
        S1 = np.vdot(dp, G@dp); S2 = np.vdot(dc, G@dc)
        den = GEFF + 1j*TWOPI*(d2 - dg) + beta*S2
        dXi = -beta*K12*K21/den
        A_off[i] = float(np.real(S1))
        A_on[i] = A_off[i] + float(np.real(dXi))
    return A_on, A_off

SCENARIOS = dict(
    single=dict(orientations=[0], n_draws=1, sigma_B_G=0.0,
                sigma_opt_GHz=0.0, sigma_d_GHz=0.0, sigma_Oc_rel=0.0,
                opt_cut_GHz=None),
    low_density=dict(orientations=[0, 1, 2, 3], n_draws=60, sigma_B_G=0.1,
                     sigma_opt_GHz=5.0, sigma_d_GHz=0.1, sigma_Oc_rel=0.05,
                     opt_cut_GHz=None),
    high_density=dict(orientations=[0, 1, 2, 3], n_draws=60, sigma_B_G=1.0,
                      sigma_opt_GHz=30.0, sigma_d_GHz=0.3, sigma_Oc_rel=0.10,
                      opt_cut_GHz=None),
    post_selected=dict(orientations=[0], n_draws=60, sigma_B_G=1.0,
                       sigma_opt_GHz=30.0, sigma_d_GHz=0.3, sigma_Oc_rel=0.10,
                       opt_cut_GHz=5.0),
    # orientation + spectral selection AND field shimmed to 0.1 G: the
    # ground-splitting spread from a 1 G inhomogeneity smears the ~0.3 MHz
    # dip by several MHz, so field homogeneity is itself a required
    # experimental selection condition (main Gate 5 finding).
    post_selected_shimmed=dict(orientations=[0], n_draws=60, sigma_B_G=0.1,
                               sigma_opt_GHz=30.0, sigma_d_GHz=0.3,
                               sigma_Oc_rel=0.10, opt_cut_GHz=5.0),
)

def run_scenario(d2s, cfg, rng):
    A_on = np.zeros(len(d2s)); A_off = np.zeros(len(d2s)); n = 0
    for k in cfg['orientations']:
        for _ in range(cfg['n_draws']):
            dB = rng.normal(0, cfg['sigma_B_G']*GAUSS_T, size=3) \
                if cfg['sigma_B_G'] else np.zeros(3)
            Bx, Bz = frame_field(k, dB)
            d = max(0.3, rp.D + (rng.normal(0, cfg['sigma_d_GHz'])
                                 if cfg['sigma_d_GHz'] else 0.0))
            Oc = rp.OC*np.exp(rng.normal(0, cfg['sigma_Oc_rel'])) \
                if cfg['sigma_Oc_rel'] else rp.OC
            dopt = rng.normal(0, cfg['sigma_opt_GHz']) \
                if cfg['sigma_opt_GHz'] else 0.0
            if cfg['opt_cut_GHz'] is not None:
                while abs(dopt) > cfg['opt_cut_GHz']:
                    dopt = rng.normal(0, cfg['sigma_opt_GHz'])
            on, off = member_spectra(d2s, Bx, Bz, float(d), float(Oc),
                                     float(dopt))
            A_on += on; A_off += off; n += 1
    return A_on/n, A_off/n

def metrics(d_MHz, C):
    i0 = int(np.argmax(np.abs(C)))
    cmax = float(C[i0]); half = abs(cmax)/2
    idx = np.where(np.abs(C) >= half)[0]
    fwhm = float(d_MHz[idx[-1]] - d_MHz[idx[0]]) if len(idx) > 1 else np.nan
    return dict(Cmax=cmax, center_MHz=float(d_MHz[i0]), fwhm_MHz=fwhm)

def main(quick=False):
    tabdir = ROOT/'results'/'tables'; figdir = ROOT/'results'/'figures'
    tabdir.mkdir(parents=True, exist_ok=True); figdir.mkdir(parents=True, exist_ok=True)
    d2s = np.linspace(-0.005, 0.005, 81 if quick else 161)
    rng = np.random.default_rng(SEED)
    rows, spectra = [], {}
    scen = dict(SCENARIOS)
    if quick:
        for k in scen: scen[k] = dict(scen[k], n_draws=min(scen[k]['n_draws'], 10))
    for name, cfg in scen.items():
        A_on, A_off = run_scenario(d2s, cfg, rng)
        C = (A_off - A_on)/A_off
        m = metrics(d2s*1e3, C)
        rows.append(dict(scenario=name, **{k: (v if not isinstance(v, list)
                                               else len(v))
                                           for k, v in cfg.items()
                                           if k != 'orientations'},
                         n_orientations=len(cfg['orientations']),
                         Cmax=m['Cmax'], fwhm_MHz=m['fwhm_MHz'],
                         center_MHz=m['center_MHz']))
        spectra[name] = C
        print(f"{name:14s} Cmax={m['Cmax']:.4g} fwhm={m['fwhm_MHz']:.3f} MHz")
    ref = next(r for r in rows if r['scenario'] == 'single')['Cmax']
    for r in rows:
        r['washout_factor'] = r['Cmax']/ref

    order = {r['scenario']: r['Cmax'] for r in rows}
    passed = dict(
        post_selection_survives=bool(order['post_selected_shimmed'] > 1e-3),
        washout_quantified=all(np.isfinite(r['washout_factor']) for r in rows),
        regimes_distinguished=bool(order['single'] > order['low_density']
                                   > order['high_density']),
    )
    with (tabdir/'gate5_ensemble_contrast.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with (tabdir/'gate5_summary.json').open('w') as fh:
        json.dump(dict(seed=SEED, quick=quick, rows=rows, gate5_pass=passed,
                       gate5_all_passed=all(passed.values())), fh, indent=2)

    plt.figure(figsize=(7.6, 5.0))
    for name, C in spectra.items():
        plt.plot(d2s*1e3, C, label=name)
    plt.xlabel('two-photon detuning (MHz)'); plt.ylabel('ensemble contrast C')
    plt.title('Gate 5: ensemble averaging of the 70 K candidate window')
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(figdir/'gate5_ensemble_spectra.png', dpi=220); plt.close()
    print('gate5_pass:', passed, '->', all(passed.values()))
    return rows, passed

if __name__ == '__main__':
    main(quick='--quick' in sys.argv[1:])
