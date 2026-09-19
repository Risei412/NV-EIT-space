"""Bounded full-model search using the repository's four-model EIT gate.

Finite-grid maximum only; never a certified global temperature ceiling.
Preserves the existing model and classifiers to make the baseline reproducible.
"""
from pathlib import Path
import sys
import json
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT/'results/audits/eit_temperature_limit'


def worker(task):
    T, Bx, Bz, Oc = task
    sys.path.insert(0,str(ROOT/'calculations/numerics/manuscript_figures/pra/src'))
    import p0_2_two_layer_classification as p0
    import p1_phase_diagram as p1
    original = p1.adaptive_spectrum
    p1.BZ0, p1.OC = Bz, Oc
    p1.adaptive_spectrum = lambda t,b: original(t,b,Oc=Oc)
    try:
        row = p0.evaluate_point((T,Bx))
        row['observable_EIT_1pct'] = (row.get('joint_class') == 'spectroscopic_EIT'
            and row.get('Cmax',0) >= .01)
        return row
    finally:
        p1.adaptive_spectrum = original


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    refine = '--refine' in sys.argv
    prefix = 'operational_refined' if refine else 'operational'
    tasks = ([(t,b,z,.1) for t in (31.,32.,33.,34.,35.)
              for b in (.175,.2,.225) for z in (.005,.05)] if refine else
             [(t,b,z,o) for t in (30.,40.,50.,60.,70.)
              for b in (.15,.2,.35,.5) for z in (.005,.05) for o in (.03,.1)])
    rows = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        pending = {pool.submit(worker,t):t for t in tasks}
        for f in as_completed(pending):
            row = f.result(); rows.append(row)
            (OUT/(prefix+'_grid.json')).write_text(json.dumps(rows,indent=2),encoding='utf-8')
            print(json.dumps({k:row.get(k) for k in ('T_K','Bx_T','Bz_T','Oc_GHz','Cmax',
                'spectral_class','best_runner_margin','observable_EIT_1pct','seconds')}),flush=True)
    valid = [r for r in rows if r['observable_EIT_1pct']]
    summary = dict(n_points=len(rows), highest_found_K=max((r['T_K'] for r in valid),default=None),
        certified_global_maximum=False, EIT_points=valid,
        limitations=['fixed branch 3, Y/Y, strain 1.683 GHz, no ISC/hyperfine',
                     'legacy four-model fit; no noise or sampling convergence audit',
                     'coarse grid, no exclusion proof between samples'])
    (OUT/(prefix+'_summary.json')).write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    main()
