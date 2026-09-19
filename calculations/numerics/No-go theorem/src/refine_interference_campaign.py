"""Root-based warm-start refinement; avoids infeasible SLSQP endpoints."""
import json
import numpy as np
from scipy.optimize import brentq,minimize
import eit_interference_certificate as ec
import nv_model as nv
from search_interference_temperature import OUT,geometry


def refine(d,cap,start):
    channel=start['channel']
    x0=np.array(list(start['parameters'].values()))
    x0[8]=cap
    bounds=[(0,.5),(.005,.05),(0,np.pi/2),(0,np.pi/2),(-np.pi,np.pi),
            (-np.pi,np.pi),(0,2*np.pi),(-4,5),(1e-5,cap)]
    upper=brentq(lambda T:ec.contrast_upper_bound(T,d,cap)-.01,20,300)
    def temperature(x):
        H,dp,dc,z=geometry(np.r_[x,channel],d)
        def contrast(T):
            return ec.certificate(H,dp,dc,z,nv.gamma_oc_GHz(T,d),x[8])['C_at_slow_center']-.01
        # Explicit bracket with a positive point; zero is not presumed monotone.
        grid=np.linspace(40,upper,9)
        vals=[contrast(t) for t in grid]
        brackets=[(a,b) for a,b,fa,fb in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]) if fa>=0 and fb<0]
        if not brackets:return 40.
        a,b=brackets[-1]
        return brentq(contrast,a,b,xtol=1e-7)
    op=minimize(lambda x:-temperature(x),x0,method='L-BFGS-B',bounds=bounds,
        options={'maxiter':100,'ftol':1e-10,'gtol':1e-5,'eps':1e-5,'maxls':30})
    T=temperature(op.x)-.001
    H,dp,dc,z=geometry(np.r_[op.x,channel],d)
    cert=ec.certificate(H,dp,dc,z,nv.gamma_oc_GHz(T,d),op.x[8])
    return dict(strain_d_GHz=d,control_cap_GHz=cap,witness_T_K=float(T),
        parameters=dict(zip(start['parameters'].keys(),map(float,op.x))),channel=channel,
        witness_certificate=cert,local_success=bool(op.success),message=str(op.message),
        certified_uniform_exclusion_above_K=upper,global_optimum_proven=False)


def main():
    original=json.loads((OUT/'extended_search.json').read_text())
    valid=[r for r in original if r['witness_certificate']['passed']]
    rows=[]
    for d in (.5,1.683,5.,10.,30.,100.):
        for cap in (.1,.3,1.,3.,10.):
            # Include a successful geometry at the same strain, if available.
            choices=sorted(valid+rows,key=lambda r:abs(r['strain_d_GHz']-d)+abs(np.log(r['control_cap_GHz']/cap)))
            start=next(r for r in choices if r['witness_certificate']['passed'])
            r=refine(d,cap,start);rows.append(r)
            (OUT/'refined_search.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
            print(json.dumps(dict(d=d,cap=cap,T=r['witness_T_K'],passed=r['witness_certificate']['passed'],
                success=r['local_success'])),flush=True)


if __name__=='__main__':main()
