"""nv_model.py -- NV 3E model, Model-A1 dipole legs, kernel & susceptibility.
UNITS: all Hamiltonians in GHz (ordinary frequency). Angular conversion 2*pi
applied ONLY inside the resolvent: G = [gamma*I + i*2pi*(H - z I)]^-1 with
gamma in GHz(ordinary, s^-1 *1e-9) -- matching the archived v2 convention.
Magnetic field input in Tesla. Rates from phonon_rates in Hz -> *1e-9 GHz.
gamma_oc = k_orb/2 (k_orb = one-way X<->Y jump rate, Happacher Eq.23)."""
import numpy as np
from functools import lru_cache
from scipy.optimize import linear_sum_assignment
import phonon_rates as pr

TWOPI=2*np.pi
s2=1/np.sqrt(2)
Sz=np.diag([-1.,0.,1.]).astype(complex)
Sx=s2*np.array([[0,1,0],[1,0,1],[0,1,0]],complex)
Sy=s2*np.array([[0,1j,0],[-1j,0,1j],[0,-1j,0]],complex)
I3=np.eye(3,dtype=complex)
sz_o=np.array([[-1,0],[0,1]],complex); sx_o=np.array([[0,1],[1,0]],complex)
sy_o=np.array([[0,1j],[-1j,0]],complex); I2=np.eye(2,dtype=complex)
OY=np.array([1,0],complex); OX=np.array([0,1],complex)
DGS=2.877; GRAD=0.0668; GE=28.02495164  # GHz, GHz, GHz/T
P=dict(Dpar=1.42,Lpar=5.33,Dperp=1.55/2,Lperp=0.20/np.sqrt(2))

def korb_GHz(T,d): return pr.k_orb(float(T),float(d))*1e-9
def gamma_oc_GHz(T,d,scale=1.0): return scale*0.5*korb_GHz(T,d)+0.5*GRAD

def Hgs(Bvec):
    Bx,By,Bz=Bvec
    return DGS*(Sz@Sz-(2/3)*I3)+GE*(Bx*Sx+By*Sy+Bz*Sz)

def Hes(Bvec,d=1.683,phi=0.0):
    Bx,By,Bz=Bvec
    dx=d*np.cos(phi); dy=d*np.sin(phi)
    H=np.kron(I2,P['Dpar']*(Sz@Sz-(2/3)*I3))
    H+=-P['Lpar']*np.kron(sy_o,Sz)
    H+=P['Dperp']*(np.kron(sz_o,Sy@Sy-Sx@Sx)-np.kron(sx_o,Sx@Sy+Sy@Sx))
    H+=P['Lperp']*(np.kron(sz_o,Sx@Sz+Sz@Sx)-np.kron(sx_o,Sy@Sz+Sz@Sy))
    H+=dx*np.kron(sz_o,I3)-dy*np.kron(sx_o,I3)
    H+=GE*(Bx*np.kron(I2,Sx)+By*np.kron(I2,Sy)+Bz*np.kron(I2,Sz))
    return H

def dressed_ground(Bvec):
    """eigvecs assigned to bare (ms=-1,0,+1) by max overlap; deterministic phase."""
    w,U=np.linalg.eigh(Hgs(Bvec))
    ov=np.abs(U.conj().T@np.eye(3))**2
    r,c=linear_sum_assignment(-ov); idx=np.empty(3,dtype=int)
    for rr,cc in zip(r,c): idx[cc]=rr
    U=U[:,idx]; w=w[idx]
    for j in range(3):
        q=int(np.argmax(np.abs(U[:,j]))); z=U[q,j]
        if abs(z): U[:,j]*=np.exp(-1j*np.angle(z))
    return w,U

def legs(U,ctrl='-1',probe_orb=OX,ctrl_orb=OY):
    gp=U[:,1]; gc=U[:,0] if ctrl=='-1' else U[:,2]
    return np.kron(probe_orb,gp),np.kron(ctrl_orb,gc)

def probe_line(H,dp):
    w,V=np.linalg.eigh(H); j=int(np.argmax(np.abs(V.conj().T@dp)**2))
    return float(w[j]),float((np.abs(V.conj().T@dp)**2)[j]),w

def response(H,dp,dc,z,gamma,Oc=1.0,gg=6.3e-5):
    beta=(TWOPI*Oc)**2/4; geff=2*gg+2e-6
    G=np.linalg.inv(gamma*np.eye(6)+1j*TWOPI*(H-z*np.eye(6)))
    K12=np.vdot(dp,G@dc);K21=np.vdot(dc,G@dp)
    S1=np.vdot(dp,G@dp);S2=np.vdot(dc,G@dc)
    den=geff+beta*S2
    dXi=-beta*K12*K21/den
    Acut=float(np.imag(1j*S1)); dA=float(-np.imag(1j*dXi))
    C=dA/Acut if abs(Acut)>1e-300 else np.nan
    return dict(C=C,Acut=Acut,K12=K12,K21=K21,S1=S1,S2=S2,den=den,z=z)

def window(H,dp,protocol,z_near=2.0,fid_min=0.01,z_glob=3000.,n=601):
    Ep,fid,w=probe_line(H,dp)
    if protocol=='A': return np.linspace(Ep-z_near,Ep+z_near,n),Ep
    if protocol=='B': return np.linspace(float(w.min()-4),float(w.max()+4),n),Ep
    if protocol=='C': return np.linspace(-z_glob,z_glob,4001),Ep
    raise ValueError

def scanC(H,dp,dc,gamma,protocol='B',Oc=1.0,gg=6.3e-5,**kw):
    zz,Ep=window(H,dp,protocol,**kw)
    rr=[response(H,dp,dc,z,gamma,Oc,gg) for z in zz]
    Amax=max(abs(r['Acut']) for r in rr)
    keep=[r for r in rr if abs(r['Acut'])>0.01*Amax] if protocol!='C' else rr
    Cs=np.array([r['C'] for r in keep]); zs=np.array([r['z'] for r in keep])
    iM=int(np.argmax(Cs)); im=int(np.argmin(Cs))
    edge=bool(iM in (0,len(Cs)-1))
    return dict(Cmax=float(Cs[iM]),Cmin=float(Cs[im]),maxabs=float(np.max(np.abs(Cs))),
                z_at_Cmax=float(zs[iM]),edge_max=edge,Ep=Ep,Cres=response(H,dp,dc,Ep,gamma,Oc,gg)['C'])
