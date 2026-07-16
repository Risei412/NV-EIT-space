"""gate3_snr_map.py -- Gate 3 of SIMULATION_PLAN.md: convert the candidate
EIT contrast into laboratory observables and map the detectable region.

Chain (signal_chain.py): delta_chi -> delta_alpha -> Delta OD -> Delta T/T
-> Delta N_ph -> SNR, applied at the fixed PRL candidate (70 K, BX0, BZ0,
Oc=0.1 GHz). The candidate contrast is taken live from run_prl_prediction
machinery (same source as Gates 1-2), gamma_h from nv_model.gamma_oc_GHz.

Outputs:
  results/tables/signal_chain_parameters.csv   (value/unit/source/uncertainty)
  results/tables/gate3_required_conditions.csv (per target contrast)
  results/figures/gate3_snr_map.png            (SNR over density x time)
  results/tables/gate3_summary.json            (pass flags)
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE))
import run_prl_prediction as rp
import nv_model as nv
import signal_chain as sc

TARGET_SNR = 5.0
PPM = 1.76e17            # NV cm^-3 per ppm in diamond

# parameter table: name -> (value, unit, source, uncertainty note)
def parameter_table(gamma_h):
    return dict(
        lambda_nm=(637.0, 'nm', 'NV- ZPL', 'exact'),
        n_refr=(2.41, '-', 'diamond @637 nm', '<1%'),
        debye_waller=(0.035, '-', 'Santori 2010 / Doherty review (3-5%)', '+-0.01'),
        gamma_rad_GHz=(nv.GRAD, 'GHz', 'repo GRAD (Happacher convention)', '+-20%'),
        gamma_h_GHz=(gamma_h, 'GHz', 'nv_model.gamma_oc_GHz(70 K, d=1.683)', 'model'),
        gamma_inh_GHz=(30.0, 'GHz', 'ensemble ZPL strain broadening (typ. 10-100)', 'x/ 3'),
        f_orient=(0.25, '-', 'one of four NV orientations', 'exact'),
        f_spin=(1/3, '-', 'unpolarized ground spin (no MW init)', 'up to ~0.9 with init'),
        eta=(0.1, '-', 'collection+detector efficiency', 'x/ 3'),
        power_W=(1e-6, 'W', 'weak probe, below sector saturation', 'design'),
        sigma_tech=(1e-6, '-', 'relative technical noise per sqrt(sample)', 'x/ 10'),
        L_cm=(0.05, 'cm', '0.5 mm sample', 'design'),
        n_nv_max_cm3=(PPM, 'cm^-3', '1 ppm (coherence-preserving upper end)', 'x/ 3'),
        tau_max_s=(3600.0, 's', 'one hour practical integration', 'design'),
    )

def main():
    tabdir = ROOT/'results'/'tables'; figdir = ROOT/'results'/'figures'
    tabdir.mkdir(parents=True, exist_ok=True); figdir.mkdir(parents=True, exist_ok=True)

    cand = rp.branch_value(70.0, rp.BX0, rp.BZ0, rp.J0, rp.OC)
    gamma_h = nv.gamma_oc_GHz(70.0, rp.D)
    C_cand = float(cand['C'])
    P = parameter_table(gamma_h)
    v = {k: p[0] for k, p in P.items()}

    with (tabdir/'signal_chain_parameters.csv').open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['parameter', 'value', 'unit', 'source', 'uncertainty'])
        for k, (val, unit, src, unc) in P.items():
            w.writerow([k, f'{val:.6g}', unit, src, unc])

    sigma = sc.sigma_zpl_cm2(v['lambda_nm'], v['n_refr'], v['debye_waller'],
                             v['gamma_rad_GHz'], v['gamma_h_GHz'])
    fspec = sc.spectral_fraction(v['gamma_h_GHz'], v['gamma_inh_GHz'])

    def od_sector(n_nv):
        return sc.od(sc.alpha_cm(sigma, n_nv, v['f_orient'], v['f_spin'], fspec),
                     v['L_cm'])

    # SNR map over density x integration time at the candidate contrast
    n_grid = np.logspace(14, np.log10(3*PPM), 121)
    tau_grid = np.logspace(-3, 5, 121)
    S = np.zeros((len(tau_grid), len(n_grid)))
    for j, n_nv in enumerate(n_grid):
        ods = od_sector(n_nv)
        dod = sc.delta_od(ods, C_cand)
        for i, tau in enumerate(tau_grid):
            S[i, j] = sc.snr(dod, ods, v['power_W'], v['lambda_nm'],
                             tau, v['eta'], v['sigma_tech'])

    # required conditions per target contrast
    rows = []
    for C_t in (1e-2, 1e-3, 1e-4):
        for n_nv in (0.01*PPM, 0.1*PPM, PPM):
            ods = od_sector(n_nv)
            dod = sc.delta_od(ods, C_t)
            tau = sc.required_tau_s(TARGET_SNR, dod, ods, v['power_W'],
                                    v['lambda_nm'], v['eta'], v['sigma_tech'])
            rows.append(dict(contrast=C_t, n_nv_cm3=n_nv, n_nv_ppm=n_nv/PPM,
                             od_sector=ods, delta_od=dod,
                             delta_T_over_T=sc.delta_T_over_T(dod),
                             tau_required_s=tau,
                             feasible=bool(tau <= v['tau_max_s'])))
    with (tabdir/'gate3_required_conditions.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # candidate-point summary + pass flags
    ods_c = od_sector(v['n_nv_max_cm3']); dod_c = sc.delta_od(ods_c, C_cand)
    tau_c = sc.required_tau_s(TARGET_SNR, dod_c, ods_c, v['power_W'],
                              v['lambda_nm'], v['eta'], v['sigma_tech'])
    feasible_region = S >= TARGET_SNR
    ok_n = n_grid <= v['n_nv_max_cm3']; ok_t = tau_grid <= v['tau_max_s']
    passed = dict(
        conditions_listed_all_targets=bool(len(rows) == 9),
        realistic_detection_region=bool(
            feasible_region[np.ix_(ok_t, ok_n)].any()),
        candidate_feasible=bool(tau_c <= v['tau_max_s']),
    )
    summary = dict(candidate_contrast=C_cand, gamma_h_GHz=gamma_h,
                   sigma_zpl_cm2=sigma, spectral_fraction=fspec,
                   od_sector_at_1ppm=ods_c, delta_od_at_1ppm=dod_c,
                   tau_required_s_at_1ppm=tau_c,
                   target_snr=TARGET_SNR, gate3_pass=passed,
                   gate3_all_passed=all(passed.values()))
    with (tabdir/'gate3_summary.json').open('w') as fh:
        json.dump(summary, fh, indent=2)

    plt.figure(figsize=(7.6, 5.4))
    pcm = plt.pcolormesh(n_grid/PPM, tau_grid, np.maximum(S, 1e-6),
                         norm=LogNorm(vmin=1e-2, vmax=max(S.max(), 10)),
                         shading='auto')
    plt.colorbar(pcm, label='SNR')
    cs = plt.contour(n_grid/PPM, tau_grid, S, levels=[TARGET_SNR],
                     colors='w', linewidths=2)
    plt.clabel(cs, fmt=lambda x: f'SNR={x:g}')
    plt.axvline(1.0, ls='--', c='r', label='1 ppm')
    plt.axhline(v['tau_max_s'], ls=':', c='r', label='1 hour')
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('NV density (ppm)'); plt.ylabel('integration time (s)')
    plt.title(f'Gate 3: shot-limited SNR of the 70 K candidate '
              f'(C={C_cand:.3%}, L=0.5 mm)')
    plt.legend(loc='lower left', fontsize=8); plt.tight_layout()
    plt.savefig(figdir/'gate3_snr_map.png', dpi=220); plt.close()

    print(f"sigma_ZPL={sigma:.3e} cm^2, OD_sector(1ppm)={ods_c:.3e}, "
          f"dOD={dod_c:.3e}, tau(SNR=5)={tau_c:.3g} s")
    print('gate3_pass:', passed, '->', all(passed.values()))
    return summary

if __name__ == '__main__':
    main()
