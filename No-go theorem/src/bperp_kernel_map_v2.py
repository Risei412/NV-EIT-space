import os,sys,json
from functools import lru_cache
import numpy as np
from scipy.optimize import linear_sum_assignment
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
from happacher_rate import k_orb,ALPHA_PH

TWOPI=2*np.pi
s=1/np.sqrt(2)
Sz=np.diag([-1.,0.,1.]).astype(complex)
Sx=s*np.array([[0,1,0],[1,0,1],[0,1,0]],complex)
Sy=s*np.array([[0,1j,0],[-1j,0,1j],[0,-1j,0]],complex)
I3=np.eye(3,dtype=complex)
sz_o=np.array([[-1,0],[0,1]],complex)
sx_o=np.array([[0,1],[1,0]],complex)
sy_o=np.array([[0,1j],[-1j,0]],complex)
I2=np.eye(2,dtype=complex)
OY=np.array([1,0],complex); OX=np.array([0,1],complex)
DGS=2.877; GRAD=0.0668; GE=28.02495164
P=dict(Dpar=1.42,Lpar=5.33,Dperp=1.55/2,Lperp=0.20/np.sqrt(2))

@lru_cache(maxsize=None)
def korb(T,d): return k_orb(T,d,ALPHA_PH)*1e-9

def Hgs(BxT,BzT):
    return DGS*(Sz@Sz-(2/3)*I3)+GE*(BxT*Sx+BzT*Sz)

@lru_cache(maxsize=None)
def Hes(BxT,BzT,d=1.683,phi=0.0):
    dx=d*np.cos(phi); dy=d*np.sin(phi)
    H=np.kron(I2,P['Dpar']*(Sz@Sz-(2/3)*I3))
    H += -P['Lpar']*np.kron(sy_o,Sz)
    H += P['Dperp']*(np.kron(sz_o,Sy@Sy-Sx@Sx)-np.kron(sx_o,Sx@Sy+Sy@Sx))
    H += P['Lperp']*(np.kron(sz_o,Sx@Sz+Sz@Sx)-np.kron(sx_o,Sy@Sz+Sz@Sy))
    H += dx*np.kron(sz_o,I3)-dy*np.kron(sx_o,I3)
    H += GE*(BxT*np.kron(I2,Sx)+BzT*np.kron(I2,Sz))
    return H

def track(Bxs,Bz):
    vals=[]; vecs=[]; Uprev=None
    for k,Bx in enumerate(Bxs):
        w,U=np.linalg.eigh(Hgs(float(Bx),float(Bz)))
        if k==0:
            # assign eigenvectors to bare (-1,0,+1)
            ov=np.abs(U.conj().T@np.eye(3))**2
            r,c=linear_sum_assignment(-ov); idx=np.empty(3,dtype=int)
            for rr,cc in zip(r,c): idx[cc]=rr
            U=U[:,idx]; w=w[idx]
        else:
            ov=np.abs(Uprev.conj().T@U)**2
            r,c=linear_sum_assignment(-ov); idx=np.empty(3,dtype=int)
            for rr,cc in zip(r,c): idx[rr]=cc
            U=U[:,idx]; w=w[idx]
        for j in range(3):
            if Uprev is None:
                q=np.argmax(np.abs(U[:,j])); z=U[q,j]
            else: z=np.vdot(Uprev[:,j],U[:,j])
            if abs(z): U[:,j]*=np.exp(-1j*np.angle(z))
        Uprev=U; vals.append(w); vecs.append(U)
    return np.array(vals),np.array(vecs)

def context(T,Bx,Bz,U,d=1.683,Oc=1.0,gg=6.3e-5):
    H=Hes(float(Bx),float(Bz),float(d),0.0)
    dp=np.kron(OX,U[:,1])  # adiabatic ms=0 branch
    dc=np.kron(OY,U[:,0])  # adiabatic ms=-1 branch
    w,V=np.linalg.eigh(H); ov=np.abs(V.conj().T@dp)**2; j=int(np.argmax(ov))
    rate=korb(float(T),float(d)); gamma=.5*(rate+GRAD)
    beta=(TWOPI*Oc)**2/4; geff=2*gg+2e-6
    return H,dp,dc,w,float(w[j]),float(ov[j]),rate,gamma,beta,geff

def eval_z(ctx,z):
    H,dp,dc,w,Eprobe,fid,rate,gamma,beta,geff=ctx
    G=np.linalg.inv(gamma*np.eye(6)+1j*TWOPI*(H-z*np.eye(6)))
    K12=np.vdot(dp,G@dc); K21=np.vdot(dc,G@dp)
    S1=np.vdot(dp,G@dp); S2=np.vdot(dc,G@dc)
    dXi=-beta*K12*K21/(geff+beta*S2)
    chi0=1j*S1; dchi=1j*dXi
    Acut=float(np.imag(chi0)); dA=float(-np.imag(dchi))
    C=dA/Acut if abs(Acut)>1e-300 else np.nan
    eta=abs(K12*K21)/max(abs(S1*S2),1e-300)
    return dict(C=C,Acut=Acut,dA=dA,eta=eta,K12=K12,K21=K21,S1=S1,S2=S2,
                z=float(z),Eprobe=Eprobe,fid=fid,rate=rate,gamma=gamma)

def eval_res(T,Bx,Bz,U,**kw):
    ctx=context(T,Bx,Bz,U,**kw)
    return eval_z(ctx,ctx[4])

def scan_line(T,Bx,Bz,U,**kw):
    ctx=context(T,Bx,Bz,U,**kw); w=ctx[3]
    # dense around manifold; room-T broad but same grid used for audit
    zz=np.linspace(float(w.min()-4),float(w.max()+4),601)
    rr=[eval_z(ctx,z) for z in zz]
    iA=int(np.argmax([abs(r['dA']) for r in rr]))
    iC=int(np.argmax([abs(r['C']) if abs(r['Acut'])>.01*max(abs(q['Acut']) for q in rr) else -1 for r in rr]))
    Amax=max(abs(r['Acut']) for r in rr)
    out=rr[iA].copy(); out['Amax']=Amax; out['dA_norm']=abs(out['dA'])/max(Amax,1e-300)
    out['C_absorption_constrained_max']=abs(rr[iC]['C']); out['z_Cmax']=rr[iC]['z']
    return out

def main(outdir):
    os.makedirs(outdir,exist_ok=True)
    Bz=0.02
    Bxs=np.linspace(0,1.0,101); _,Us=track(Bxs,Bz)
    Ts=np.array([5,10,20,30,50,77,100,150,200,250,300.])
    C=np.zeros((len(Ts),len(Bxs))); ETA=np.zeros_like(C); DA=np.zeros_like(C); FID=np.zeros_like(C)
    for i,T in enumerate(Ts):
        for j,Bx in enumerate(Bxs):
            r=eval_res(T,Bx,Bz,Us[j]); C[i,j]=abs(r['C']); ETA[i,j]=r['eta']; DA[i,j]=abs(r['dA']); FID[i,j]=r['fid']
    # room-temperature line scan at selected fields
    selected=[]
    for Bx in [0,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0]:
        j=int(np.argmin(abs(Bxs-Bx))); rr=scan_line(300,float(Bxs[j]),Bz,Us[j])
        rr2=eval_res(300,float(Bxs[j]),Bz,Us[j])
        selected.append(dict(Bx_T=float(Bxs[j]),C_res=abs(rr2['C']),eta_res=rr2['eta'],
                             dA_res=abs(rr2['dA']),fid=rr2['fid'],
                             dA_norm_max=rr['dA_norm'],C_at_max_dA=abs(rr['C']),z_max_dA=rr['z'],
                             C_absorption_constrained_max=rr['C_absorption_constrained_max']))
    # low-field scaling with axial bias, use excess over Bx=0
    Bsmall=np.linspace(0,0.03,31); _,Usm=track(Bsmall,Bz)
    sm=[eval_res(300,float(B),Bz,U) for B,U in zip(Bsmall,Usm)]
    def slope(metric):
        y=np.array([metric(r) for r in sm]); dy=np.abs(y-y[0]); x=Bsmall
        m=(x>0)&(dy>max(np.max(dy)*1e-10,1e-30))
        return float(np.polyfit(np.log10(x[m]),np.log10(dy[m]),1)[0]) if m.sum()>4 else np.nan
    slopes={'delta_eta_power':slope(lambda r:r['eta']),
            'delta_C_power':slope(lambda r:abs(r['C'])),
            'delta_dA_power':slope(lambda r:abs(r['dA']))}
    # robustness vs axial bias at 300 K
    bias_rows=[]
    for bz in [0.005,0.01,0.02,0.05]:
        bxv=np.linspace(0,0.3,61); _,UU=track(bxv,bz)
        vals=[eval_res(300,float(bx),bz,u) for bx,u in zip(bxv,UU)]
        k=int(np.argmax([abs(r['C']) for r in vals]))
        bias_rows.append(dict(Bz_T=bz,Bx_at_max=float(bxv[k]),C0=abs(vals[0]['C']),Cmax=abs(vals[k]['C']),
                              enhancement=abs(vals[k]['C'])/max(abs(vals[0]['C']),1e-300),
                              eta0=vals[0]['eta'],etamax=max(r['eta'] for r in vals)))
    summary={'model':'corrected 3E, dressed ground states, orbital-selective Model A1, Bz channel-selection bias',
             'Bz_T':Bz,'Omega_c_GHz':1.0,'gamma_g_GHz':6.3e-5,'delta_perp_GHz':1.683,
             'korb_300K_GHz':korb(300.,1.683),'gamma_optcoh_300K_GHz':.5*(korb(300.,1.683)+GRAD),
             'low_field_scaling':slopes,'selected_300K':selected,'bias_robustness':bias_rows}
    np.savez(os.path.join(outdir,'bperp_kernel_maps_v2.npz'),Bxs=Bxs,Ts=Ts,C=C,ETA=ETA,DA=DA,FID=FID,Bz=Bz,
             Bsmall=Bsmall,Csmall=np.array([abs(r['C']) for r in sm]),
             Etasmall=np.array([r['eta'] for r in sm]),DAsmall=np.array([abs(r['dA']) for r in sm]))
    with open(os.path.join(outdir,'bperp_kernel_summary_v2.json'),'w') as f: json.dump(summary,f,indent=2)
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main(sys.argv[1] if len(sys.argv)>1 else 'bperp_results_v2')
