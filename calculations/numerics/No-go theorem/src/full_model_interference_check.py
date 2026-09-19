"""Direct column-vectorized full-model check for arbitrary one-photon z.

Uses explicit matrix readout; leaves legacy shared utilities untouched.
Includes control-driven steady populations, ISC/singlet shelving and secular
14N hyperfine (absorptions averaged over mI=-1,0,+1 before forming C), with
the gate2 rates. Both are ON by default: hyperfine lowers C by up to ~40% at
Oc<=0.3 GHz (full_model_isc_hyperfine.json); ISC changes C by <0.4%.
Pass isc=False, hyperfine=False for the bare 9-level model.
"""
import numpy as np
import gate2_candidate_full_vs_reduced as g2
import liouvillian_core as lc


def trace_solve(L,rhs,trace):
    N=int(round(np.sqrt(len(rhs))))
    Q=L.copy();b=rhs.copy()
    Q[0,:]=0;Q[0,np.arange(N)*(N+1)]=1;b[0]=trace
    x=np.linalg.solve(Q,b)
    return x,float(np.linalg.norm(L@x-rhs))


def point(T,Bx,Bz,Oc,d,phi,ppol,cpol,ctrl,z,delta=0,isc=True,hyperfine=True):
    if hyperfine:
        rows=[_point(T,Bx,Bz,Oc,d,phi,ppol,cpol,ctrl,z,delta,isc,mI=m) for m in (-1,0,1)]
        af=float(np.mean([r['Afull'] for r in rows]));ac=float(np.mean([r['Acut'] for r in rows]))
        return dict(C_sector=(ac-af)/ac,Afull=af,Acut=ac,
            residual_steady=max(r['residual_steady'] for r in rows),
            residual_full=max(r['residual_full'] for r in rows),
            residual_cut=max(r['residual_cut'] for r in rows),
            ground_probe_population=float(np.mean([r['ground_probe_population'] for r in rows])),
            singlet_population=float(np.mean([r['singlet_population'] for r in rows])),
            rho_min_eigenvalue=min(r['rho_min_eigenvalue'] for r in rows))
    return _point(T,Bx,Bz,Oc,d,phi,ppol,cpol,ctrl,z,delta,isc)


def _point(T,Bx,Bz,Oc,d,phi,ppol,cpol,ctrl,z,delta=0,isc=False,mI=None):
    H,Ls,Vp,dp,meta=g2.build_full(T,Bx,Bz,delta,Oc=Oc,d=d,phi=phi,
        ppol=ppol,cpol=cpol,ctrl=ctrl,j0=0,isc=isc,mI=mI)
    H[3:9,3:9]-=2*np.pi*(z-meta['z0'])*np.eye(6)
    L=lc.liouvillian(H,Ls);N=meta['N'];p=meta['p_idx'];c=meta['c_idx']
    rho0,r0=trace_solve(L,np.zeros(N*N,complex),1)
    Hp=np.pi*g2.OP*(Vp+Vp.conj().T)
    V=-1j*(np.kron(np.eye(N),Hp)-np.kron(Hp.T,np.eye(N)))
    rhs=-V@rho0
    S=[c+N*p,p+N*c];X=[k for k in range(N*N) if k not in S]
    cut=L.copy();cut[np.ix_(S,X)]=0;cut[np.ix_(X,S)]=0
    xf,rf=trace_solve(L,rhs,0);xc,rc=trace_solve(cut,rhs,0)
    af=float((-2*np.vdot(dp,xf.reshape((N,N),order='F')[3:9,p])/g2.OP).imag)
    ac=float((-2*np.vdot(dp,xc.reshape((N,N),order='F')[3:9,p])/g2.OP).imag)
    return dict(C_sector=(ac-af)/ac,Afull=af,Acut=ac,residual_steady=r0,
        residual_full=rf,residual_cut=rc,ground_probe_population=float(rho0[p+N*p].real),
        singlet_population=float(rho0[9+N*9].real) if N==10 else 0.,
        rho_min_eigenvalue=float(np.linalg.eigvalsh(rho0.reshape((N,N),order='F')).min()))
