"""Search certified interference witnesses, with a separate global exclusion bound.

Search maximum is NOT a proven global optimum. Analytical upper bound is
uniform in B, polarizations, strain azimuth and both optical detunings,
within the stated reduced-model assumptions.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import differential_evolution,minimize,brentq
import nv_model as nv
import eit_interference_certificate as ec

OUT=Path(__file__).resolve().parents[4]/'results/audits/eit_interference_temperature'


def geometry(x,d):
    Bx,Bz,ap,ac,pp,pc,phi,zfraction,Oc,channel=x
    B=(Bx,0,Bz)
    H=nv.Hes(B,d,phi); _,U=nv.dressed_ground(B)
    p=np.array([np.cos(ap),np.exp(1j*pp)*np.sin(ap)])
    c=np.array([np.cos(ac),np.exp(1j*pc)*np.sin(ac)])
    dp=np.kron(p,U[:,1]); dc=np.kron(c,U[:,2 if round(channel) else 0])
    w=np.linalg.eigvalsh(H)
    z=float(w[0]+zfraction*(w[-1]-w[0]))
    return H,dp,dc,z


def optimize(d,cap,seed,zrange=(0,1),maxiter=160):
    bounds=[(0,.5),(.005,.05),(0,np.pi/2),(0,np.pi/2),(-np.pi,np.pi),
            (-np.pi,np.pi),(0,2*np.pi),zrange,(0,cap),(0,1)]
    Tstart=70.
    gamma=nv.gamma_oc_GHz(Tstart,d)
    def loss(x):
        H,dp,dc,z=geometry(x,d)
        return -nv.response(H,dp,dc,z,gamma,x[8])['C']
    de=differential_evolution(loss,bounds,seed=seed,popsize=10,maxiter=maxiter,
                              polish=False,integrality=[False]*9+[True],tol=1e-7)
    channel=round(de.x[-1])
    def cert_at(y):
        T=y[0]; x=np.r_[y[1:],channel]
        H,dp,dc,z=geometry(x,d)
        return ec.certificate(H,dp,dc,z,nv.gamma_oc_GHz(T,d),x[8])
    def constraints(y):
        c=cert_at(y)
        return np.array([c['C_at_slow_center']-.01,c['decay_disk_margin'],
                         -c['slow_residue_real']*1e6-1e-6])
    y0=np.r_[70.,de.x[:-1]]
    op=minimize(lambda y:-y[0],y0,method='SLSQP',bounds=[(20,180)]+bounds[:-1],
                constraints=[{'type':'ineq','fun':constraints}],
                options={'ftol':1e-8,'maxiter':200})
    c=cert_at(op.x)
    # Back away from the equality boundary to provide a verified feasible point.
    witness_y=op.x.copy(); witness_y[0]-=.001
    wc=cert_at(witness_y)
    upper=brentq(lambda T:ec.contrast_upper_bound(T,d,cap)-.01,20,300)
    return dict(strain_d_GHz=d,control_cap_GHz=cap,seed=seed,
        searched_z_fraction=list(zrange),
        found_boundary_K=float(op.x[0]),witness_T_K=float(witness_y[0]),
        local_success=bool(op.success),local_message=str(op.message),
        DE_success=bool(de.success),DE_nfev=de.nfev,
        parameters=dict(zip(['Bx','Bz','probe_theta','control_theta','probe_phase',
            'control_phase','strain_phi','z_fraction','Oc'],map(float,op.x[1:]))),
        channel=int(channel),boundary_certificate=c,witness_certificate=wc,
        certified_uniform_exclusion_above_K=float(upper),
        global_optimum_proven=False)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for d in (.5,1.683,5.,10.):
        for cap in (.1,.3,1.):
            for seed in (20260920,20260921):
                r=optimize(d,cap,seed); rows.append(r)
                (OUT/'search.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
                print(json.dumps(r),flush=True)


if __name__=='__main__':main()
