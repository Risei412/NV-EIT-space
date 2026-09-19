"""Independent numerical checks for the temperature audit; writes only new audit output."""
import json
import numpy as np
from scipy.optimize import brentq
from audit_eit_temperature_limit import OUT, vectorization_audit
import nv_model as nv
import run_prl_prediction as rp
import gate2_candidate_full_vs_reduced as g2
import liouvillian_core as lc


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    R=2*np.pi*(np.ptp(np.linalg.eigvalsh(nv.Hes((0,0,0))))+2*nv.GE*np.hypot(.5,.05))
    beta=(2*np.pi*.1)**2/4
    geff=2*6.3e-5+2e-6
    def power_bound(T):
        gamma=nv.gamma_oc_GHz(T,rp.D)
        return beta*R**2/(gamma**3*(geff+beta*gamma/(gamma**2+R**2)))
    rng=np.random.default_rng(123)
    for _ in range(200):
        T=rng.uniform(20,300)
        p=rng.normal(size=2)+1j*rng.normal(size=2); p/=np.linalg.norm(p)
        c=rng.normal(size=2)+1j*rng.normal(size=2); c/=np.linalg.norm(c)
        H,dp,dc,w,V,pp,cc=rp.model_at(rng.uniform(0,.5),rng.uniform(.005,.05),p,c)
        z=rng.uniform(w.min(),w.max()); gamma=nv.gamma_oc_GHz(T,rp.D)
        r=nv.response(H,dp,dc,z,gamma,rng.uniform(0,.1))
        assert abs(np.vdot(dp,dc))<1e-12
        assert r['C']<=power_bound(T)+1e-12
        assert nv.response(H,dp,dc,z,gamma,0)['C']==0
    args=(70.,rp.BX0,.005,0.)
    old=g2.chi_pair(*args)
    original_vec,original_unvec=lc.vec,lc.unvec
    try:
        lc.vec=lambda rho:rho.reshape(-1,order='F')
        lc.unvec=lambda v,n:v.reshape((n,n),order='F')
        corrected=g2.chi_pair(*args)
    finally:
        lc.vec,lc.unvec=original_vec,original_unvec
    delta=max(abs(old[0]-corrected[0]),abs(old[1]-corrected[1]))
    assert delta<1e-12
    fixed=brentq(lambda t:rp.branch_value(t,rp.BX0,rp.BZ0,rp.J0)['C']-.01,60,80)
    assert abs(fixed-71.31455440980206)<1e-8
    result=dict(random_bound_and_zero_control_checks=200,
                vectorization=vectorization_audit(),
                candidate_column_change_abs=delta,
                fixed_candidate_T_1pct_K=fixed,
                power_limited_contrast_bound_crossing_K=brentq(lambda t:power_bound(t)-.01,20,300))
    (OUT/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    # Add the tighter bound to the existing audit without repeating the DE search.
    path=OUT/'audit.json'
    if path.exists():
        audit=json.loads(path.read_text())
        audit['power_limited_contrast_bound_crossing_K']=result['power_limited_contrast_bound_crossing_K']
        audit['upper_temperature_monotonicity']='fixed-energy Bose integral increases with T; see audit note; scope <=300 K'
        path.write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(result),flush=True)


if __name__=='__main__': main()
