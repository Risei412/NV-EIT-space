from __future__ import annotations
import csv, json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import nv_model as nv

Y=np.array([1.0,0.0],complex)
X=np.array([0.0,1.0],complex)
D=1.683
PHI=0.0
CTRL='+1'
OC=0.1
BZ0=0.005
BX0=0.23226114680059184
J0=3
MINW=0.1


def model_at(Bx:float,Bz:float,ppol=Y,cpol=Y,ctrl=CTRL,d=D,phi=PHI):
    B=(float(Bx),0.0,float(Bz))
    H=nv.Hes(B,d,phi)
    _,U=nv.dressed_ground(B)
    gp=U[:,1]
    gc=U[:,0] if ctrl=='-1' else U[:,2]
    dp=np.kron(ppol,gp); dc=np.kron(cpol,gc)
    w,V=np.linalg.eigh(H)
    pp=np.abs(V.conj().T@dp)**2
    cc=np.abs(V.conj().T@dc)**2
    return H,dp,dc,w,V,pp,cc


def branch_value(T:float,Bx:float,Bz:float,j:int,Oc:float=OC,ppol=Y,cpol=Y,ctrl=CTRL,d=D,phi=PHI):
    H,dp,dc,w,V,pp,cc=model_at(Bx,Bz,ppol,cpol,ctrl,d,phi)
    r=nv.response(H,dp,dc,float(w[j]),nv.gamma_oc_GHz(T,d),Oc)
    return dict(T=float(T),Bx=float(Bx),Bz=float(Bz),Oc=float(Oc),j=int(j),z=float(w[j]),
                pp=float(pp[j]),cc=float(cc[j]),joint=float(np.sqrt(pp[j]*cc[j])),
                C=float(r['C']),Acut=float(r['Acut']),dA=float(r['C']*r['Acut']))


def envelope(T:float,Bx:float,Bz:float,Oc:float=OC,ppol=Y,cpol=Y,ctrl=CTRL,d=D,phi=PHI,minw=MINW):
    H,dp,dc,w,V,pp,cc=model_at(Bx,Bz,ppol,cpol,ctrl,d,phi)
    best=None
    for j,z in enumerate(w):
        if pp[j] < minw or cc[j] < minw:
            continue
        r=nv.response(H,dp,dc,float(z),nv.gamma_oc_GHz(T,d),Oc)
        q=dict(T=float(T),Bx=float(Bx),Bz=float(Bz),Oc=float(Oc),j=int(j),z=float(z),
               pp=float(pp[j]),cc=float(cc[j]),joint=float(np.sqrt(pp[j]*cc[j])),
               C=float(r['C']),Acut=float(r['Acut']),dA=float(r['C']*r['Acut']))
        if best is None or q['C'] > best['C']:
            best=q
    return best


def write_csv(path:Path,rows:list[dict]):
    if not rows: return
    keys=list(rows[0].keys())
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def thresholds(Ts,Cs,targets):
    out={}
    for target in targets:
        roots=[]
        fvals=Cs-target
        for a,b,fa,fb in zip(Ts[:-1],Ts[1:],fvals[:-1],fvals[1:]):
            if fa==0: roots.append(float(a))
            elif fa*fb<0:
                roots.append(float(brentq(lambda t: branch_value(t,BX0,BZ0,J0,OC)['C']-target,float(a),float(b))))
        out[str(target)]=roots
    return out


def main():
    outdir=ROOT/'outputs'; figdir=ROOT/'figures'; tabdir=ROOT/'tables'
    outdir.mkdir(exist_ok=True);figdir.mkdir(exist_ok=True);tabdir.mkdir(exist_ok=True)

    # Fixed experimental candidate and exact thresholds.
    candidate=branch_value(70.0,BX0,BZ0,J0,OC)
    Ts=np.linspace(20,140,241)
    trows=[branch_value(float(T),BX0,BZ0,J0,OC) for T in Ts]
    Cs=np.array([r['C'] for r in trows]); dAs=np.array([r['dA'] for r in trows])
    th=thresholds(Ts,Cs,[0.01,0.001,0.0001,0.0])
    for target in [1e-4,1e-5,1e-6]:
        roots=[]; vals=dAs-target
        for a,b,fa,fb in zip(Ts[:-1],Ts[1:],vals[:-1],vals[1:]):
            if fa*fb<0:
                roots.append(float(brentq(lambda t: branch_value(t,BX0,BZ0,J0,OC)['dA']-target,float(a),float(b))))
        th['dA_'+str(target)]=roots
    write_csv(tabdir/'fixed_candidate_temperature_sweep.csv',trows)

    H,dp,dc,w,V,pp,cc=model_at(BX0,BZ0)
    spread=float(w.max()-w.min())
    crossings=[]
    for key,roots in th.items():
        for rt in roots:
            crossings.append(dict(quantity=key,T_K=rt,gamma_oc_GHz=float(nv.gamma_oc_GHz(rt,D)),
                                  gamma_over_excited_spread=float(nv.gamma_oc_GHz(rt,D)/spread)))
    write_csv(tabdir/'threshold_crossings.csv',crossings)

    # Branch-resolved phase map: exact eigenenergy, both optical weights >=10%.
    Tgrid=np.arange(20.0,140.0001,1.0)
    Bgrid=np.arange(0.0,0.50001,0.0025)
    Cmap=np.full((len(Tgrid),len(Bgrid)),np.nan)
    dAmap=np.full_like(Cmap,np.nan); Jmap=np.full_like(Cmap,np.nan)
    PPmap=np.full_like(Cmap,np.nan); CCmap=np.full_like(Cmap,np.nan)
    for ib,Bx in enumerate(Bgrid):
        H,dp,dc,w,V,pp,cc=model_at(float(Bx),BZ0)
        eligible=np.where((pp>=MINW)&(cc>=MINW))[0]
        for it,T in enumerate(Tgrid):
            best=None
            for j in eligible:
                r=nv.response(H,dp,dc,float(w[j]),nv.gamma_oc_GHz(float(T),D),OC)
                q=(float(r['C']),float(r['C']*r['Acut']),int(j),float(pp[j]),float(cc[j]))
                if best is None or q[0]>best[0]: best=q
            if best:
                Cmap[it,ib],dAmap[it,ib],Jmap[it,ib],PPmap[it,ib],CCmap[it,ib]=best
    np.savez_compressed(outdir/'branch_resolved_phase_map.npz',T_K=Tgrid,Bx_T=Bgrid,C=Cmap,dA=dAmap,
                        branch=Jmap,probe_weight=PPmap,control_weight=CCmap)

    # Bx-Bz map at 70 K.
    Bxs=np.arange(0.0,0.5001,0.005)
    Bzs=np.arange(0.005,0.0501,0.001)
    Fmap=np.full((len(Bzs),len(Bxs)),np.nan); Fbranch=np.full_like(Fmap,np.nan)
    for iz,Bz in enumerate(Bzs):
        for ix,Bx in enumerate(Bxs):
            q=envelope(70.0,float(Bx),float(Bz),OC)
            if q:
                Fmap[iz,ix]=q['C'];Fbranch[iz,ix]=q['j']
    np.savez_compressed(outdir/'field_robustness_map_70K.npz',Bx_T=Bxs,Bz_T=Bzs,C=Fmap,branch=Fbranch)

    # Control sweep.
    Ocs=np.logspace(-3,1,161)
    orows=[branch_value(70.0,BX0,BZ0,J0,float(o)) for o in Ocs]
    write_csv(tabdir/'control_Rabi_sweep_70K.csv',orows)

    # Linear polarization tolerance around Y on the fixed physical branch.
    deg=np.linspace(-15,15,121)
    Pmap=np.full((len(deg),len(deg)),np.nan); PWp=np.full_like(Pmap,np.nan);PWc=np.full_like(Pmap,np.nan)
    for i,ap in enumerate(np.deg2rad(deg)):
        ppol=np.array([np.cos(ap),np.sin(ap)],complex)
        for k,ac in enumerate(np.deg2rad(deg)):
            cpol=np.array([np.cos(ac),np.sin(ac)],complex)
            q=branch_value(70.0,BX0,BZ0,J0,OC,ppol,cpol)
            Pmap[i,k]=q['C'];PWp[i,k]=q['pp'];PWc[i,k]=q['cc']
    np.savez_compressed(outdir/'polarization_tolerance_70K.npz',probe_angle_deg=deg,control_angle_deg=deg,C=Pmap,
                        probe_weight=PWp,control_weight=PWc)

    # Two-photon spectrum at fixed one-photon resonance.  The original response()
    # is evaluated at delta_2ph=0; here the ground-coherence denominator is
    # extended as den -> den + i*2*pi*delta_2ph.
    Hs,dps,dcs,ws,Vs,pps,ccs=model_at(BX0,BZ0)
    z0=float(ws[J0]); gamma=nv.gamma_oc_GHz(70,D)
    G=np.linalg.inv(gamma*np.eye(6)+1j*2*np.pi*(Hs-z0*np.eye(6)))
    K12=np.vdot(dps,G@dcs);K21=np.vdot(dcs,G@dps)
    S1=np.vdot(dps,G@dps);S2=np.vdot(dcs,G@dcs)
    beta=(2*np.pi*OC)**2/4; geff=2*6.3e-5+2e-6
    Acut=float(np.real(S1))
    d2s=np.linspace(-0.02,0.02,1601)  # GHz = +/-20 MHz
    spec=[]
    for d2 in d2s:
        den=geff+1j*2*np.pi*d2+beta*S2
        dXi=-beta*K12*K21/den
        dA=float(-np.real(dXi)); C=float(dA/Acut)
        spec.append(dict(two_photon_detuning_GHz=float(d2),two_photon_detuning_MHz=float(d2*1e3),
                         A_control_off=Acut,transparency_reduction=dA,A_control_on=float(Acut-dA),C=C))
    write_csv(tabdir/'two_photon_spectrum_70K.csv',spec)
    sx=np.array([r['two_photon_detuning_MHz'] for r in spec]); sa=np.array([r['A_control_on'] for r in spec])
    imin=int(np.argmin(sa)); depth=float(Acut-sa[imin]); half=float(Acut-depth/2)
    li=np.where(sa[:imin]>=half)[0]; ri=np.where(sa[imin:]>=half)[0]
    def interp(i,k): return float(sx[i]+(half-sa[i])*(sx[k]-sx[i])/(sa[k]-sa[i]))
    left=interp(int(li[-1]),int(li[-1]+1)); jr=int(imin+ri[0]); right=interp(jr-1,jr)
    two_photon_metrics=dict(center_MHz=float(sx[imin]),FWHM_MHz=float(right-left),
                            maximum_depth=float(depth),maximum_relative_depth=float(depth/Acut))

    # Strain robustness: take common-branch envelope at exact eigenenergies.
    srows=[]
    for d in np.linspace(0.5,5.0,19):
        for phideg in range(0,91,15):
            q=envelope(70.0,BX0,BZ0,OC,d=float(d),phi=np.deg2rad(phideg))
            if q:
                q.update(strain_GHz=float(d),strain_azimuth_deg=phideg)
                srows.append(q)
    write_csv(tabdir/'strain_robustness_70K.csv',srows)

    # Room-temperature hard check over basis channels/polarizations/fields at exact common branches.
    basis={'Y':Y,'X':X}; rt=[]
    for Oc in (0.01,0.1,1.0,10.0):
        for ctrl in ('-1','+1'):
            for pn,pv in basis.items():
                for cn,cv in basis.items():
                    for Bx in np.linspace(0,1,41):
                        for Bz in np.linspace(0.005,0.05,7):
                            q=envelope(300.0,float(Bx),float(Bz),Oc,pv,cv,ctrl,minw=MINW)
                            if q:
                                q.update(ctrl=ctrl,polarization=pn+cn)
                                rt.append(q)
    rt_sorted=sorted(rt,key=lambda q:q['C'])
    room_summary=dict(min_C=rt_sorted[0],max_C=max(rt,key=lambda q:q['C']),max_abs=max(rt,key=lambda q:abs(q['C'])),n=len(rt))

    # Kernel moment and high-gamma slope for candidate.
    H,dp,dc,w,V,pp,cc=model_at(BX0,BZ0)
    z=float(w[J0]); Xop=2*np.pi*(H-z*np.eye(6))
    moments=[];v=dc.copy()
    for n in range(4):
        moments.append(np.vdot(dp,v));v=Xop@v
    gammas=np.logspace(2,5,61);K=[]
    for g in gammas:
        G=np.linalg.inv(g*np.eye(6)+1j*Xop);K.append(abs(np.vdot(dp,G@dc)))
    slope=float(np.polyfit(np.log10(gammas[-20:]),np.log10(K[-20:]),1)[0])

    summary=dict(
        primary_claim='Transverse-field brightening of a shared excited branch extends observable NV spin-Lambda EIT to approximately 90 K, followed by a dissipation-driven sign reversal at 101.44 K and room-temperature practical no-go.',
        candidate=candidate,
        excited_manifold_spread_GHz=spread,
        thresholds=th,
        threshold_rows=crossings,
        room_temperature_common_branch_scan=room_summary,
        candidate_moments_angular=[{'n':n,'real':float(m.real),'imag':float(m.imag),'abs':float(abs(m))} for n,m in enumerate(moments)],
        candidate_kernel_high_gamma_slope=slope,
        two_photon_resonance_70K=two_photon_metrics,
        conditions=dict(channel='ms=0 to ms=+1',probe_polarization='Y linear',control_polarization='Y linear',
                        Omega_c_GHz=OC,Bz_T=BZ0,Bx_T=BX0,strain_GHz=D,branch_index=J0,
                        branch_rule='exact excited eigenenergy; probe and control oscillator weights each >= 0.1')
    )
    with (outdir/'summary.json').open('w',encoding='utf-8') as f:json.dump(summary,f,indent=2)

    # Figures: separate plots, no subplots.
    plt.figure(figsize=(7.2,5.2))
    pcm=plt.pcolormesh(Bgrid,Tgrid,Cmap,shading='auto',norm=SymLogNorm(linthresh=1e-5,vmin=np.nanmin(Cmap),vmax=np.nanmax(Cmap)))
    plt.colorbar(pcm,label='Local EIT contrast C')
    plt.contour(Bgrid,Tgrid,Cmap,levels=[1e-4,1e-3,1e-2],linewidths=1.0)
    plt.scatter([BX0],[70],marker='x',s=60,label='Experimental candidate')
    plt.xlabel(r'$B_\perp$ (T)');plt.ylabel('Temperature (K)');plt.title('Branch-resolved NV EIT envelope (exact eigenenergies)');plt.legend();plt.tight_layout()
    plt.savefig(figdir/'fig1_branch_resolved_phase_map.png',dpi=220);plt.savefig(figdir/'fig1_branch_resolved_phase_map.pdf');plt.close()

    plt.figure(figsize=(7.2,4.8))
    plt.semilogy(Ts,np.maximum(Cs,1e-12),label='Positive local contrast')
    for y in (1e-2,1e-3,1e-4): plt.axhline(y,linestyle='--',linewidth=.8)
    plt.axvline(th['0.0'][0],linestyle=':',label=f"sign reversal {th['0.0'][0]:.2f} K")
    plt.xlabel('Temperature (K)');plt.ylabel('Local EIT contrast C');plt.title('Dissipation-driven transparency-to-absorption transition');plt.legend();plt.tight_layout()
    plt.savefig(figdir/'fig2_temperature_thresholds.png',dpi=220);plt.savefig(figdir/'fig2_temperature_thresholds.pdf');plt.close()

    plt.figure(figsize=(7.2,4.8))
    sp=np.array([r['two_photon_detuning_MHz'] for r in spec]);a0=np.array([r['A_control_off'] for r in spec]);a1=np.array([r['A_control_on'] for r in spec])
    plt.plot(sp,a0,label='Control off');plt.plot(sp,a1,label='Control on')
    plt.xlabel('Two-photon detuning (MHz)');plt.ylabel('Local absorption (arb.)');plt.title('Predicted two-photon EIT resonance at 70 K');plt.legend();plt.tight_layout()
    plt.savefig(figdir/'fig3_candidate_spectrum_70K.png',dpi=220);plt.savefig(figdir/'fig3_candidate_spectrum_70K.pdf');plt.close()

    plt.figure(figsize=(7.2,5.2))
    pc=plt.pcolormesh(Bxs,Bzs*1e3,Fmap,shading='auto',norm=SymLogNorm(linthresh=1e-5,vmin=np.nanmin(Fmap),vmax=np.nanmax(Fmap)))
    plt.colorbar(pc,label='Local EIT contrast C');plt.contour(Bxs,Bzs*1e3,Fmap,levels=[1e-4,1e-3,1e-2],linewidths=1.0)
    plt.scatter([BX0],[BZ0*1e3],marker='x',s=60);plt.xlabel(r'$B_\perp$ (T)');plt.ylabel(r'$B_z$ (mT)');plt.title('Magnetic-field robustness at 70 K');plt.tight_layout()
    plt.savefig(figdir/'fig4_field_robustness_70K.png',dpi=220);plt.savefig(figdir/'fig4_field_robustness_70K.pdf');plt.close()

    plt.figure(figsize=(7.2,5.2))
    pc=plt.pcolormesh(deg,deg,Pmap,shading='auto',norm=SymLogNorm(linthresh=1e-5,vmin=np.nanmin(Pmap),vmax=np.nanmax(Pmap)))
    plt.colorbar(pc,label='Local EIT contrast C');plt.contour(deg,deg,Pmap,levels=[1e-4,1e-3,1e-2],linewidths=1.0)
    plt.xlabel('Control linear-polarization error (deg)');plt.ylabel('Probe linear-polarization error (deg)');plt.title('Polarization tolerance at 70 K');plt.tight_layout()
    plt.savefig(figdir/'fig5_polarization_tolerance_70K.png',dpi=220);plt.savefig(figdir/'fig5_polarization_tolerance_70K.pdf');plt.close()

    plt.figure(figsize=(7.2,4.8))
    plt.loglog(Ocs,np.maximum([r['C'] for r in orows],1e-12))
    plt.axhline(1e-4,linestyle='--',linewidth=.8);plt.axhline(1e-3,linestyle='--',linewidth=.8)
    plt.xlabel(r'Control Rabi frequency $\Omega_c/2\pi$ (GHz)');plt.ylabel('Local EIT contrast C');plt.title('Control-power requirement at 70 K');plt.tight_layout()
    plt.savefig(figdir/'fig6_control_sweep_70K.png',dpi=220);plt.savefig(figdir/'fig6_control_sweep_70K.pdf');plt.close()

    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
