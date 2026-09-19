"""Check the highest grid EIT label against fit initialization and sampling."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'calculations/numerics/manuscript_figures/pra/src'))
import p1_phase_diagram as p1
import p0_2_two_layer_classification as p0
import gate1_candidate_aic_bootstrap as g1


def main():
    rows=[]
    for n in (161,321):
        p1.N_FINE=n
        d,A,Ac,C,info=p1.adaptive_spectrum(31.,.2)
        for seed in (None,1,2,3,4):
            fit=g1.fit_all(d,A,rng=None if seed is None else np.random.default_rng(seed))
            row=dict(n=n,seed=seed,half_MHz=info['half_MHz'],**p0.classify_spectrum(fit))
            rows.append(row)
            print(json.dumps(row),flush=True)
    path=ROOT/'results/audits/eit_temperature_limit/boundary_stability.json'
    path.write_text(json.dumps(rows,indent=2),encoding='utf-8')


if __name__=='__main__':main()
