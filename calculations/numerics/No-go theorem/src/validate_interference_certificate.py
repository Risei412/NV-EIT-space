"""Validation of pole certificate, Schur response and uniform mixing bound."""
import json
import numpy as np
import nv_model as nv
import run_prl_prediction as rp
import eit_interference_certificate as ec
from search_interference_temperature import OUT


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    H=np.diag(np.arange(6)*10.).astype(complex)
    dp=np.eye(6,dtype=complex)[:,0]; dc=dp.copy()
    eit=ec.certificate(H,dp,dc,0,1.,.02)
    ats=ec.certificate(H,dp,dc,0,1.,1.)
    assert eit['passed'] and not ats['passed']
    rng=np.random.default_rng(20260920)
    max_error=0.;max_bound_ratio=0.
    for _ in range(500):
        T=rng.uniform(20,300); d=rng.uniform(.5,100); Oc=10**rng.uniform(-4,1)
        B=(rng.uniform(0,.5),0,rng.uniform(.005,.05)); phi=rng.uniform(0,2*np.pi)
        _,U=nv.dressed_ground(B); H=nv.Hes(B,d,phi)
        p=rng.normal(size=2)+1j*rng.normal(size=2);p/=np.linalg.norm(p)
        c=rng.normal(size=2)+1j*rng.normal(size=2);c/=np.linalg.norm(c)
        dp=np.kron(p,U[:,1]);dc=np.kron(c,U[:,0 if rng.random()<.5 else 2])
        gamma=nv.gamma_oc_GHz(T,d);z=rng.uniform(-10,10)*gamma/(2*np.pi)
        r=nv.response(H,dp,dc,z,gamma,Oc)
        f=np.r_[dp,0j]; M=ec.matrix(H,dc,z,gamma,Oc)
        A=float(np.vdot(f,np.linalg.solve(M,f)).real)
        err=abs(A-r['Acut']*(1-r['C']))
        max_error=max(max_error,err)
        assert err<1e-10
        bound=ec.contrast_upper_bound(T,d,Oc)
        max_bound_ratio=max(max_bound_ratio,r['C']/bound)
        assert r['C']<=bound+1e-11
        assert np.linalg.eigvalsh((M+M.conj().T)/2).min()>0
    H,dp,dc,w,*_=rp.model_at(rp.BX0,rp.BZ0)
    candidate=ec.certificate(H,dp,dc,w[3],nv.gamma_oc_GHz(70,rp.D),.1)
    assert candidate['passed']
    result=dict(ideal_lambda_EIT=eit,ideal_lambda_ATS=ats,
        random_cases=500,max_schur_error=max_error,max_C_over_bound=max_bound_ratio,
        mixing_bound_rate=ec.mixing_bound(),candidate_70K=candidate)
    (OUT/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
