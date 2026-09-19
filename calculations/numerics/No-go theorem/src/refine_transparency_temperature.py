"""Local constrained refinement of a transparency witness (not EIT proof)."""
import json
import numpy as np
from scipy.optimize import minimize
from audit_eit_temperature_limit import OUT
import run_prl_prediction as rp
import nv_model as nv


def values(y):
    T,Bx,Bz,ap,ac,pp,pc,Oc=y
    p=np.array([np.cos(ap),np.exp(1j*pp)*np.sin(ap)])
    c=np.array([np.cos(ac),np.exp(1j*pc)*np.sin(ac)])
    H,dp,dc,w,V,wp,wc=rp.model_at(Bx,Bz,p,c,'+1')
    C=nv.response(H,dp,dc,w[3],nv.gamma_oc_GHz(T,rp.D),Oc)['C']
    return np.array([C-.01,wp[3]-.1,wc[3]-.1])


def main():
    rows=json.loads((OUT/'audit.json').read_text())['witnesses']
    results=[]
    for r in rows:
        if r['T_K']!=70: continue
        p=r['parameters']
        x=[71,p['Bx_T'],p['Bz_T'],p['probe_theta'],p['control_theta'],
           p['probe_phase'],p['control_phase'],p['Oc_GHz']]
        opt=minimize(lambda y:-y[0],x,method='SLSQP',bounds=[(20,120),(0,.5),(.005,.05),
                     (0,np.pi/2),(0,np.pi/2),(-np.pi,np.pi),(-np.pi,np.pi),(0,.1)],
                     constraints=[{'type':'ineq','fun':values}],
                     options={'ftol':1e-10,'maxiter':300})
        out=dict(T_K=float(opt.x[0]),parameters=list(map(float,opt.x)),
                 constraint_slack=values(opt.x).tolist(),success=bool(opt.success),
                 message=str(opt.message),is_EIT_certified=False,is_global_maximum=False)
        results.append(out)
    (OUT/'transparency_refinement.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results),flush=True)


if __name__=='__main__': main()
