import os,sys,json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
import bperp_kernel_map_v2 as km
from weak_probe_response import steady_state,first_order,vec
from nv_system import liouvillian
N=9; DIM=81; TWOPI=2*np.pi

def leg_matrix(d,ground_index):
    V=np.zeros((N,N),complex); V[3:9,ground_index]=d; return V

def build(T,Bx,Bz,U,Oc=1.0,gg=6.3e-5,d=1.683,delta2=0.0):
    ctx=km.context(T,Bx,Bz,U,d=d,Oc=Oc,gg=gg)
    H6,dp,dc,w,Eprobe,fid,rate,gamma,beta,geff=ctx
    H=np.zeros((N,N),complex); H[3:9,3:9]=H6-Eprobe*np.eye(6)
    # ground dressed basis ordered adiabatically (-1,0,+1); probe=1, control=0
    Eg=np.linalg.eigvalsh(km.Hgs(Bx,Bz))
    # spectator only; control shifted by two-photon frame
    H[0,0]=-delta2
    H[2,2]=float(Eg.max()-Eg.min()) # spectator far; exact value immaterial for selected sector
    Vc=leg_matrix(dc,0); Vp=leg_matrix(dp,1)
    H += 0.5*Oc*(Vc+Vc.conj().T)
    Ls=[]
    # orbital hopping
    for m in range(3):
        up=np.zeros((N,N),complex); dn=np.zeros((N,N),complex)
        up[6+m,3+m]=1; dn[3+m,6+m]=1
        Ls += [np.sqrt(rate)*up,np.sqrt(rate)*dn]
    # radiative: spin-conserving in bare basis, expressed in dressed ground basis
    for orb0 in (3,6):
        for m in range(3):
            J=np.zeros((N,N),complex)
            for a in range(3): J[a,orb0+m]=np.conj(U[m,a])
            Ls.append(np.sqrt(km.GRAD)*J)
    # dressed-basis symmetric T1 and dephasing, same phenomenological convention
    T1=1e-6
    for a in range(3):
        for b in range(3):
            if a!=b:
                J=np.zeros((N,N),complex); J[b,a]=1; Ls.append(np.sqrt(T1)*J)
    for a in range(3):
        J=np.zeros((N,N),complex); J[a,a]=1; Ls.append(np.sqrt(2*gg)*J)
    return TWOPI*H,Ls,Vp,dp,ctx

def probe_super(Vp,Op=1e-5):
    Hp=TWOPI*0.5*Op*(Vp+Vp.conj().T); I=np.eye(N)
    return -1j*(np.kron(I,Hp)-np.kron(Hp.T,I))

def detector(dp):
    d=np.zeros(DIM,complex)
    for e,a in enumerate(dp): d[(3+e)+N*1]=np.conj(a)
    return d

def cut(L):
    S=[1+N*0,0+N*1]; X=[k for k in range(DIM) if k not in S]
    M=L.copy(); M[np.ix_(S,X)]=0; M[np.ix_(X,S)]=0; return M

def run_one(T,Bx,Bz,U):
    H,Ls,Vp,dp,ctx=build(T,Bx,Bz,U)
    L=liouvillian(H,Ls); rho0,r0,gap=steady_state(L)
    V=probe_super(Vp); b=-(V@vec(rho0)); d=detector(dp)
    xf,rf=first_order(L,V,rho0); Mc=cut(L); xc,rc=first_order(Mc,V,rho0)
    Op=1e-5
    chif=-2*(d.conj()@xf)/Op; chic=-2*(d.conj()@xc)/Op
    C=(np.imag(chic)-np.imag(chif))/np.imag(chic)
    red=km.eval_z(ctx,ctx[4])
    return dict(T=T,Bx=Bx,C_full=float(np.real(C)),C_red=float(np.real(red['C'])),
                ratio=float(np.real(C/red['C'])) if red['C']!=0 else np.nan,
                res0=float(r0),res_full=float(rf),res_cut=float(rc),gap=float(gap),fid=float(ctx[5]))

def main(out):
    Bz=.02; Bxs=np.array([0,.01,.02,.05,.1,.2,.5,1.0]); _,Us=km.track(Bxs,Bz)
    rows=[]
    for T in [100.,150.,300.]:
        for B,U in zip(Bxs,Us):
            r=run_one(T,float(B),Bz,U); rows.append(r); print(r)
    with open(out,'w') as f: json.dump(rows,f,indent=2)
if __name__=='__main__': main(sys.argv[1] if len(sys.argv)>1 else 'bperp_full_validation.json')
