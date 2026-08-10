"""P0-6: unresolved four-orientation ODMR rank-4 gate.

Tests whether [100] and [110] equal-|B| fields remain distinguishable when the
four NV orientations cannot be individually assigned.  Uses the exact spin-1
ground Hamiltonian, 14N hyperfine triplets, Lorentzian broadening and nuisance
fits allowing field magnitude, linewidth, contrast, baseline, global frequency
shift, +/- transition imbalance and up to +/-20% orientation-weight imbalance.

The analytic low-field statement for an equal-weight four-axis ensemble with an
even single-line profile L is

 A(f,B)=L + (gamma^2 B^2/6)L''
        + gamma^4[(B^2)^2/72 - S4/108]L'''' + O(B^6),
 S4=Bx^4+By^4+Bz^4.

Thus direction dependence is absent through rank 2 and first appears at rank 4.
For white/background-dominated noise, orientation Fisher information therefore
scales as B^8 in the unresolved low-field limit.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks

D=2.877; GAMMA=28.02495164; A14=0.00216
FWHM=0.030; CONTRAST=0.03; NPTS=601
S2=1/np.sqrt(2)
Sx=S2*np.array([[0,1,0],[1,0,1],[0,1,0]],complex)
Sy=S2*np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]],complex)
Sz=np.diag([1,0,-1]).astype(complex); Sz2=Sz@Sz
AXES=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],float)/np.sqrt(3)
U100=np.array([1.,0,0]); U110=np.array([1.,1.,0])/np.sqrt(2)
FGRID=np.linspace(D-0.15,D+0.15,NPTS)


def frame(n):
    r=np.array([1.,0,0]) if abs(n[0])<0.9 else np.array([0.,1.,0])
    x=r-(r@n)*n; x/=np.linalg.norm(x); y=np.cross(n,x); y/=np.linalg.norm(y)
    return np.vstack([x,y,n])
FRAMES=[frame(n) for n in AXES]

# deterministic coherent-orientation strain stress test, rms scale ~ few MHz
EVALS=np.array([0.001586603208264,0.002243008247593,0.004533035733277,0.006710221224016])
PHIS=np.array([1.582884144999216,1.809953919601361,2.328707045371354,2.028753175842423])


def transitions(Bglob,Evals=EVALS,phis=PHIS):
    out=[]
    for a,F in enumerate(FRAMES):
        B=F@np.asarray(Bglob,float); E=float(Evals[a]); ph=float(phis[a])
        Ex=E*np.cos(2*ph); Ey=E*np.sin(2*ph)
        H=D*Sz2+GAMMA*(B[0]*Sx+B[1]*Sy+B[2]*Sz) \
          +Ex*(Sx@Sx-Sy@Sy)+Ey*(Sx@Sy+Sy@Sx)
        w,V=np.linalg.eigh(H); j0=int(np.argmax(np.abs(V[1,:])**2))
        out.append(sorted(abs(float(w[j]-w[j0])) for j in range(3) if j!=j0))
    return np.asarray(out)


def lorentz(f,f0,w):
    g=w/2; return g*g/((f-f0)**2+g*g)


def spectrum(Bglob,width=FWHM,ow=None,asym=0.,shift=0.,amp=CONTRAST,base=1.):
    fs=transitions(Bglob); ow=np.ones(4)/4 if ow is None else np.asarray(ow,float)
    tw=np.array([1-asym,1+asym]); y=np.zeros_like(FGRID)
    for a in range(4):
        for s in range(2):
            for m in (-1,0,1):
                y += ow[a]*tw[s]/3*lorentz(FGRID-shift,fs[a,s]+m*A14,width)
    y/=np.max(y)
    return base-amp*y


def fit_wrong(Btrue,rel_or=0.2,trans_rel=0.2):
    target=spectrum(Btrue*U100)
    low=.25*(1-rel_or); high=.25*(1+rel_or)
    def unpack(p):
        w=np.r_[p[5:8],1-np.sum(p[5:8])]; return w,float(p[8])
    def residual(p):
        w,a=unpack(p)
        model=spectrum(p[0]*U110,p[1],w,a,p[3],p[2],p[4])
        pen=max(low-w[3],0.)+max(w[3]-high,0.)
        return np.r_[model-target,100*pen]
    p0=np.r_[[Btrue,FWHM,CONTRAST,0.,1.],[.25,.25,.25],0.]
    lo=np.r_[[1e-5,.005,.001,-.02,.98],[low]*3,-trans_rel]
    hi=np.r_[[.01,.12,.1,.02,1.02],[high]*3,trans_rel]
    r=least_squares(residual,p0,bounds=(lo,hi),max_nfev=2000)
    w,a=unpack(r.x); rr=spectrum(r.x[0]*U110,r.x[1],w,a,r.x[3],r.x[2],r.x[4])-target
    return r.x,w,rr


def raw_difference(B):
    return float(np.sqrt(np.mean((spectrum(B*U100)-spectrum(B*U110))**2)))


def main():
    rows=[]
    for B in (0.0002,0.0003,0.0005,0.0007,0.001):
        p,w,rr=fit_wrong(B)
        y100=1-spectrum(B*U100); y110=1-spectrum(B*U110)
        n100=len(find_peaks(y100,prominence=0.001)[0]); n110=len(find_peaks(y110,prominence=0.001)[0])
        rms=float(np.sqrt(np.mean(rr**2)))
        rows.append(dict(B_mT=B*1e3,peaks_100=n100,peaks_110=n110,
            wrong_direction_rms=rms,max_residual=float(np.max(abs(rr))),
            fitted_orientation_weights=w.tolist(),fitted_B_mT=float(p[0]*1e3),
            fitted_fwhm_MHz=float(p[1]*1e3),fitted_transition_asym=float(p[8]),
            per_point_noise_for_5sigma=float(rms*np.sqrt(NPTS)/5),
            dprime_if_noise_1e3=float(rms/1e-3*np.sqrt(NPTS)),
            dprime_if_noise_1e4=float(rms/1e-4*np.sqrt(NPTS))))
    Bs=np.array([5e-5,7.5e-5,1e-4,1.5e-4,2e-4,2.5e-4,3e-4])
    dr=np.array([raw_difference(b) for b in Bs])
    amp_exp=float(np.polyfit(np.log(Bs),np.log(dr),1)[0])
    dist2_exp=float(np.polyfit(np.log(Bs),np.log(dr**2),1)[0])
    out=dict(run='RCI-20260810-P0-6',linewidth_MHz=30,contrast=CONTRAST,
        hyperfine_MHz=A14*1e3,strain_GHz=EVALS.tolist(),rows=rows,
        low_field_scaling=dict(spectral_difference_exponent=amp_exp,
            squared_distance_exponent=dist2_exp,
            analytic='directional spectrum term = -gamma^4 S4 L''''/108 + O(B^6); Fisher information ~ B^8'),
        verdict=dict(witness='PASS WITH NARROWING',prl_novelty='HOLD',
            caveat='requires near-balanced/calibrated orientation weights; not full-vector reconstruction'))
    path=Path(__file__).resolve().parents[1]/'results'/'tables'/'p0_6_unresolved_rank4_odmr_gate.json'
    path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
