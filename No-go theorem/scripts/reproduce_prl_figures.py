"""reproduce_prl_figures.py -- single entry point for the PRL numerical
campaign (SIMULATION_PLAN.md Sec. 6): runs the candidate prediction and
verification gates 1-5 in order and asserts every expected output exists.

Usage:
  python scripts/reproduce_prl_figures.py          # full run
  python scripts/reproduce_prl_figures.py --quick  # CI-sized run
"""
from __future__ import annotations
import subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT/'src'
RES = ROOT/'results'

STEPS = [
    ('run_prl_prediction.py', []),          # candidate, phase maps, figs 1-6
    ('gate1_candidate_aic_bootstrap.py', ['--quick']),
    ('gate2_candidate_full_vs_reduced.py', ['--quick']),
    ('gate3_snr_map.py', []),
    ('gate4_threshold_uncertainty.py', ['--quick']),
    ('gate5_ensemble_average.py', ['--quick']),
]

EXPECTED = [
    RES/'tables'/'gate1_aic_bootstrap.json',
    RES/'tables'/'gate1_window_sweep.csv',
    RES/'tables'/'gate1_oc_sweep.csv',
    RES/'figures'/'gate1_model_selection_stability.png',
    RES/'tables'/'gate2_candidate_comparison.csv',
    RES/'figures'/'gate2_full_vs_reduced_spectrum.png',
    RES/'tables'/'signal_chain_parameters.csv',
    RES/'tables'/'gate3_required_conditions.csv',
    RES/'figures'/'gate3_snr_map.png',
    RES/'tables'/'gate4_threshold_bands.csv',
    RES/'figures'/'gate4_threshold_bands.png',
    RES/'tables'/'gate5_ensemble_contrast.csv',
    RES/'figures'/'gate5_ensemble_spectra.png',
]

def main(quick=False):
    for script, quick_args in STEPS:
        args = [sys.executable, str(SRC/script)] + (quick_args if quick else [])
        t0 = time.time()
        print(f'== {script} {"(quick)" if quick and quick_args else ""}',
              flush=True)
        subprocess.run(args, cwd=str(SRC), check=True)
        print(f'   done in {time.time()-t0:.1f} s', flush=True)
    missing = [str(p) for p in EXPECTED if not p.exists()]
    if missing:
        raise SystemExit('missing outputs:\n  ' + '\n  '.join(missing))
    print(f'all {len(EXPECTED)} expected outputs present under {RES}')

if __name__ == '__main__':
    main(quick='--quick' in sys.argv[1:])
