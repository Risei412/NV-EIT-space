"""Effect of ISC/singlet shelving and 14N hyperfine on the full-model 1% temperature.

Evaluates the column-vectorized full model at the 70 K candidate and at the
certified reduced-model witnesses (refined_search.json). For each toggle the
two-photon detuning is re-optimized, since hyperfine shifts the mI dips.
Witness geometries are NOT re-optimized with ISC/hyperfine: the reported
temperatures are lower bounds on the full-model optimum at each geometry.
"""
import json
import numpy as np
from scipy.optimize import brentq,minimize_scalar
import nv_model as nv
import run_prl_prediction as rp
import full_model_interference_check as fm
from search_interference_temperature import OUT

CONFIGS={'base':dict(isc=False,hyperfine=False),'isc':dict(isc=True,hyperfine=False),
         'hyperfine':dict(isc=False,hyperfine=True),
         'isc+hyperfine':dict(isc=True,hyperfine=True)}
HALF=.003       # minimum two-photon detuning window, GHz


def half_window(Oc):
    # light-shifted dip centre grows ~Oc^2; window must contain it
    return max(HALF,.005*Oc**2)


def witness_args(r):
    q=r['parameters'];d=r['strain_d_GHz']
    B=(q['Bx'],0.,q['Bz'])
    w=np.linalg.eigvalsh(nv.Hes(B,d,q['strain_phi']))
    return dict(Bx=q['Bx'],Bz=q['Bz'],Oc=q['Oc'],d=d,phi=q['strain_phi'],
        ppol=np.array([np.cos(q['probe_theta']),np.exp(1j*q['probe_phase'])*np.sin(q['probe_theta'])]),
        cpol=np.array([np.cos(q['control_theta']),np.exp(1j*q['control_phase'])*np.sin(q['control_theta'])]),
        ctrl='+1' if r['channel'] else '-1',z=float(w[0]+q['z_fraction']*(w[-1]-w[0])))


def candidate_args():
    w=np.linalg.eigvalsh(nv.Hes((rp.BX0,0,rp.BZ0),rp.D,rp.PHI))
    return dict(Bx=rp.BX0,Bz=rp.BZ0,Oc=rp.OC,d=rp.D,phi=rp.PHI,ppol=rp.Y,cpol=rp.Y,
                ctrl='+1',z=float(w[rp.J0]))


def best_delta(T,a,kw):
    f=lambda x:fm.point(T,**a,delta=x,**kw)['C_sector']
    h=half_window(a['Oc']);grid=np.linspace(-h,h,61);vals=[f(x) for x in grid]
    i=int(np.argmax(vals));lo=grid[max(i-1,0)];hi=grid[min(i+1,60)]
    op=minimize_scalar(lambda x:-f(x),bounds=(lo,hi),method='bounded',options={'xatol':1e-7})
    x=float(op.x) if -op.fun>=vals[i] else float(grid[i])
    return x,f(x),bool(i in (0,60))


def threshold(a,kw,Tlo=40.,Thi=110.):
    g=lambda T:best_delta(T,a,kw)[1]-.01
    grid=np.linspace(Tlo,Thi,15);v=[g(t) for t in grid]
    br=[(s,e) for s,e,fs,fe in zip(grid[:-1],grid[1:],v[:-1],v[1:]) if fs>=0>fe]
    if not br:return None
    return float(brentq(g,*br[-1],xtol=1e-3))


def main(only=None):
    """only: case names to recompute; other rows of the existing JSON are kept."""
    rows=json.loads((OUT/'refined_search.json').read_text())
    cases=[('candidate_70K',candidate_args())]
    for key in [(1.683,.1),(1.683,1.),(1.683,10.),(.5,.1),(10.,.1),(100.,.3)]:
        r=next(r for r in rows if (r['strain_d_GHz'],r['control_cap_GHz'])==key)
        cases.append((f'witness_d{key[0]}_Oc{key[1]}',witness_args(r)))
    path=OUT/'full_model_isc_hyperfine.json'
    out=[r for r in json.loads(path.read_text()) if r['case'] not in only] if only else []
    for name,a in cases:
        if only and name not in only:continue
        for cfg,kw in CONFIGS.items():
            T=threshold(a,kw)
            x,C,edge=best_delta(T,a,kw) if T else (None,None,None)
            pt=fm.point(T,**a,delta=x,**kw) if T else {}
            row=dict(case=name,config=cfg,T_1pct_K=T,delta_opt_MHz=None if x is None else 1e3*x,
                     delta_at_window_edge=edge,C_at_70K=best_delta(70.,a,kw)[1],
                     singlet_population=pt.get('singlet_population'),
                     max_residual=max(pt.get('residual_full',0),pt.get('residual_cut',0)),
                     d_GHz=a['d'],Oc_GHz=a['Oc'],delta_window_MHz=1e3*half_window(a['Oc']))
            out.append(row);print(json.dumps(row),flush=True)
            path.write_text(json.dumps(out,indent=2),encoding='utf-8')


if __name__=='__main__':
    import sys
    main(sys.argv[1:] or None)
