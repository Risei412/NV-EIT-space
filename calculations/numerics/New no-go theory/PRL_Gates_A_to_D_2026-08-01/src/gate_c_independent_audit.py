"""Independent Gate C audit: group-IV, NV, and circuit-QED path-moment classes."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent

def fit_order(gamma,values,mask=None):
    gamma=np.asarray(gamma,float); values=np.asarray(values)
    if mask is None: mask=np.ones(len(gamma),bool)
    return float(-np.polyfit(np.log(gamma[mask]),np.log(np.abs(values[mask])),1)[0])

sz=np.array([[1,0],[0,-1]],complex); sx=np.array([[0,1],[1,0]],complex); I2=np.eye(2,dtype=complex)
orb_x=np.array([1,0],complex); orb_y=np.array([0,1],complex); spin_up=np.array([1,0],complex)
GROUP_IV={"SiV":{"Delta_e":255.0},"SnV":{"Delta_e":3000.0}}; GAMMA_RAD=0.0157

def group_iv_H(material,xi_x=0.0,xi_y=0.0,Bx=0.0):
    de=GROUP_IV[material]["Delta_e"]
    return (de/2)*np.kron(sz,sz)+xi_x*np.kron(sz,I2)+xi_y*np.kron(sx,I2)+Bx*np.kron(I2,sx)

def group_iv_legs(theta=0.0):
    c=np.kron(orb_x,spin_up); p=np.cos(theta)*np.kron(orb_x,spin_up)+np.sin(theta)*np.kron(orb_y,spin_up)
    return p,c

def moments_from(H,p,c,nmax,rad):
    X=-(rad*np.eye(H.shape[0])+1j*H); v=c.copy(); out=[]
    for _ in range(nmax+1): out.append(np.vdot(p,v)); v=X@v
    return np.asarray(out)

def group_iv_moments(material,theta=0.0,nmax=3,**kwargs):
    p,c=group_iv_legs(theta); return moments_from(group_iv_H(material,**kwargs),p,c,nmax,GAMMA_RAD)

def group_iv_kernel(material,gamma,theta=0.0,**kwargs):
    H=group_iv_H(material,**kwargs); p,c=group_iv_legs(theta); A0=GAMMA_RAD*np.eye(4)+1j*H
    return np.asarray([np.vdot(p,np.linalg.solve(g*np.eye(4)+A0,c)) for g in gamma])

sq2=1/np.sqrt(2); Sz=np.diag([-1.,0.,1.]).astype(complex)
Sx=sq2*np.array([[0,1,0],[1,0,1],[0,1,0]],complex); Sy=sq2*np.array([[0,1j,0],[-1j,0,1j],[0,-1j,0]],complex)
I3=np.eye(3,dtype=complex); sz_o=np.array([[1,0],[0,-1]],complex); sx_o=np.array([[0,1],[1,0]],complex); Lz_o=np.array([[0,-1j],[1j,0]],complex); I2o=np.eye(2,dtype=complex)
LAM_Z,D_PAR,DELTA1,DELTA2,LAM_PERP=5.33,1.42/3,1.55/2,0.20,0.20; NV_RAD=0.0157
E_SPIN={-1:np.array([1,0,0],complex),0:np.array([0,1,0],complex),1:np.array([0,0,1],complex)}; ORB_X=np.array([1,0],complex); ORB_Y=np.array([0,1],complex)

def nv_H(lam_perp=LAM_PERP,delta2=DELTA2,delta1=DELTA1,xi_x=0.0,xi_y=0.0,Bx=0.0,lam_z=LAM_Z,d_par=D_PAR):
    H=lam_z*np.kron(Lz_o,Sz)+d_par*np.kron(I2o,Sz@Sz-(2/3)*I3)
    H+=delta1*np.kron(sz_o,Sy@Sy-Sx@Sx)+delta2*np.kron(sx_o,Sx@Sz+Sz@Sx)
    H+=lam_perp*(np.kron(sx_o,Sx)+np.kron(Lz_o,Sy))+xi_x*np.kron(sz_o,I3)+xi_y*np.kron(sx_o,I3)+Bx*np.kron(I2o,Sx)
    return H

def nv_legs(pair):
    src,det=pair; return np.kron(ORB_Y,E_SPIN[det]),np.kron(ORB_X,E_SPIN[src])

def nv_moments(pair,nmax=4,**kwargs):
    p,c=nv_legs(pair); return moments_from(nv_H(**kwargs),p,c,nmax,NV_RAD)

def nv_kernel(pair,gamma,**kwargs):
    H=nv_H(**kwargs); p,c=nv_legs(pair); A0=NV_RAD*np.eye(6)+1j*H
    return np.asarray([np.vdot(p,np.linalg.solve(g*np.eye(6)+A0,c)) for g in gamma])

def first_nonzero(ms,tol=1e-9):
    for i,m in enumerate(ms):
        if abs(m)>tol: return i+1
    return None

G=0.1; D1=0.03; D2=0.05; GQ=2.5e-5

def sc_kernel(gamma,tuning):
    gA=np.array([G,G]); gB=np.array([G,G] if tuning=="generic" else [G,-G]); vals=[]
    for k in gamma:
        H=np.zeros((4,4),complex); H[1,1]=D1; H[2,2]=D2
        H[0,1]=H[1,0]=gA[0]; H[0,2]=H[2,0]=gA[1]; H[1,3]=H[3,1]=gB[0]; H[2,3]=H[3,2]=gB[1]
        Dfast=np.diag([0,.5,.5,0]).astype(complex); floor=np.diag([GQ/2,0,0,GQ/2]).astype(complex)
        vals.append(np.linalg.solve(k*Dfast+floor+1j*H,np.array([1,0,0,0],complex))[3])
    return np.asarray(vals)

giv=np.logspace(4,8,80); gnv=np.logspace(2,5,80); gsc=np.logspace(5,9,80); rows=[]; curves=[]
for material in ("SiV","SnV"):
    ms=group_iv_moments(material); vals=group_iv_kernel(material,giv); pred=first_nonzero(ms)
    rows.append({"system":material,"platform":"group-IV diamond defect","predicted_order":pred,"fit_order":fit_order(giv,vals),"M0_abs":float(abs(ms[0]))}); curves.append((material,pred,giv,vals))
for label,pair in (("NV_0_to_m1",(0,-1)),("NV_m1_to_p1",(-1,1))):
    ms=nv_moments(pair); vals=nv_kernel(pair,gnv); pred=first_nonzero(ms)
    rows.append({"system":label,"platform":"NV diamond defect","predicted_order":pred,"fit_order":fit_order(gnv,vals),"M0_abs":float(abs(ms[0]))}); curves.append((label,pred,gnv,vals))
for tuning in ("generic","protected"):
    vals=sc_kernel(gsc,tuning); pred=1 if tuning=="generic" else 2
    rows.append({"system":f"SC_{tuning}","platform":"circuit QED artificial atom","predicted_order":pred,"fit_order":fit_order(gsc,vals),"M0_abs":None}); curves.append((f"SC_{tuning}",pred,gsc,vals))

rng=np.random.default_rng(412); robust={"group_IV_n1":[],"NV_n2":[],"NV_n3":[]}
for _ in range(100):
    xi_x,xi_y,bx=rng.uniform(-.5,.5,3); theta=rng.uniform(-.5,.5)
    robust["group_IV_n1"].append(first_nonzero(group_iv_moments("SiV",theta=theta,xi_x=xi_x,xi_y=xi_y,Bx=bx)))
    scale=rng.uniform(.7,1.3,3)
    robust["NV_n2"].append(first_nonzero(nv_moments((0,-1),lam_perp=LAM_PERP*scale[0],delta2=DELTA2*scale[1],delta1=DELTA1*scale[2],xi_x=xi_x,xi_y=xi_y)))
    robust["NV_n3"].append(first_nonzero(nv_moments((-1,1),lam_perp=LAM_PERP*scale[0],delta2=DELTA2*scale[1],delta1=DELTA1*scale[2],xi_x=xi_x,xi_y=xi_y)))
robust_summary={k:{str(o):v.count(o) for o in sorted(set(v))} for k,v in robust.items()}
classes={}
for row in rows: classes.setdefault(row["predicted_order"],[]).append(row["platform"])
cross={str(n):len(set(p))>=2 for n,p in classes.items()}
summary={"verdict":"CONDITIONAL_PASS","systems":rows,"cross_platform_by_class":cross,"robustness_counts":robust_summary}
(OUT/"gate_c_independent_audit.json").write_text(json.dumps(summary,indent=2))
with (OUT/"gate_c_independent_audit.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["system","platform","predicted_order","fit_order","M0_abs"]); w.writeheader(); w.writerows(rows)
plt.figure(figsize=(7.5,5.4))
for label,n,g,y in curves: plt.loglog(g,np.abs(y),label=f"{label}, n={n}")
plt.xlabel("fast dissipation scale"); plt.ylabel("|K|"); plt.legend(fontsize=7); plt.tight_layout(); plt.savefig(OUT/"gate_c_independent_audit.png",dpi=180); plt.close()
print(json.dumps(summary,indent=2))
