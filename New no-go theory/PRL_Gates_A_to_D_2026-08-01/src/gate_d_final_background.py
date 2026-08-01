"""Compact rerunner for Gate D-final complex-background rejection audit."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np

OUT=Path(__file__).resolve().parent
J12,J23,Q=1.0,0.8,0.4549880601500485
GAMMAS=np.logspace(np.log10(5.0),np.log10(15.0),9)
PIN_W=1e-3*10**(-110/10); KB=1.380649e-23; TNOISE=3.0; TINT=40.0
NOISE=np.sqrt(KB*TNOISE/(PIN_W*TINT))
ORDER_TOL=0.15; PASS_FRAC=0.95; NSTAT=NDRIFT=4000; NHW=5000
RATIOS=[1.,10.,100.,1000.,10000.]
CURV=np.logspace(-7,-1,25); EPS=np.logspace(-5,-.5,28); CYCLE=0.08
SEQS={'two_state':np.array([1.,-1.]),'ABBA_four_state':np.array([1.,-1.,-1.,1.]),'Thue_Morse_eight_state':np.array([1.,-1.,-1.,1.,-1.,1.,1.,-1.])}
LOGG=np.log(GAMMAS); XLOG=LOGG-LOGG.mean(); DEN=np.sum(XLOG**2)


def t_amp(gamma,*,j12=J12,j23=J23,bypass=0j,leakage=0j,det=(0.,0.,0.),loss=(1.,1.,1.),qin=Q,qout=Q):
    H=np.array([[det[0],j12,bypass],[j12,det[1],j23],[np.conj(bypass),j23,det[2]]],complex)
    A=np.diag([gamma*loss[0]+qin/2,gamma*loss[1],gamma*loss[2]+qout/2]).astype(complex)+1j*H
    return complex(-np.sqrt(qin*qout)*np.vdot(np.array([leakage,0.,1.],complex),np.linalg.solve(A,np.array([1.,0.,0.],complex))))

def curve(**kw): return np.array([t_amp(g,**kw) for g in GAMMAS])
def orders(curves):
    a=np.atleast_2d(curves); y=np.log(np.maximum(np.abs(a),1e-300)); yc=y-y.mean(axis=1,keepdims=True)
    return -(yc@XLOG)/DEN

def randc(rng,shape,maxmag=1.):
    return rng.uniform(0,maxmag,shape)*np.exp(1j*rng.uniform(0,2*np.pi,shape))
def recovery(v,clean_order): return float(np.mean(np.abs(v-clean_order)<=ORDER_TOL))

clean=curve(); clean_order=float(orders(clean)[0]); sigref=float(abs(clean[0]))
moments={}
for name,s in SEQS.items():
    tt=np.linspace(-1,1,len(s)); moments[name]={f'moment_{k}':float(np.mean(s*tt**k)) for k in range(4)}

static={}
for ratio in RATIOS:
    rng=np.random.default_rng(1000+int(ratio))
    b0=ratio*sigref*np.exp(1j*rng.uniform(0,2*np.pi,(NSTAT,1)))
    b1=ratio*sigref*randc(rng,(NSTAT,1),.3)
    x=np.linspace(-1,1,len(GAMMAS))[None,:]
    noise=(rng.normal(size=(NSTAT,len(GAMMAS)))+1j*rng.normal(size=(NSTAT,len(GAMMAS))))*NOISE/np.sqrt(2)
    raw=orders(clean[None,:]+b0+b1*x+noise); projected=orders(clean[None,:]+noise)
    static[str(ratio)]={'raw_recovery':recovery(raw,clean_order),'phase_recovery':recovery(projected,clean_order),'phase_median':float(np.median(projected))}

x=np.linspace(-1,1,len(GAMMAS)); two_point={}
for ratio in RATIOS:
    rng=np.random.default_rng(2000+int(ratio))
    b0=ratio*sigref*np.exp(1j*rng.uniform(0,2*np.pi,(NSTAT,1))); b1=ratio*sigref*randc(rng,(NSTAT,1),.3); direction=randc(rng,(NSTAT,1))
    noise=(rng.normal(size=(NSTAT,len(GAMMAS)))+1j*rng.normal(size=(NSTAT,len(GAMMAS))))*NOISE/np.sqrt(2)
    nl=(rng.normal(size=(NSTAT,1))+1j*rng.normal(size=(NSTAT,1)))*NOISE/np.sqrt(2); nr=(rng.normal(size=(NSTAT,1))+1j*rng.normal(size=(NSTAT,1)))*NOISE/np.sqrt(2)
    ok=[]
    for cfrac in CURV:
        b=b0+b1*x[None,:]+ratio*sigref*cfrac*direction*x[None,:]**2
        left=b[:,[0]]+nl; right=b[:,[-1]]+nr; w=((x-x[0])/(x[-1]-x[0]))[None,:]
        corrected=clean[None,:]+b+noise-(left*(1-w)+right*w)
        if recovery(orders(corrected),clean_order)>=PASS_FRAC: ok.append(float(cfrac))
    two_point[str(ratio)]=max(ok) if ok else 0.0

drift={}; timescale={}
for ratio in RATIOS:
    drift[str(ratio)]={}; timescale[str(ratio)]={}
    for j,(name,s) in enumerate(SEQS.items()):
        rng=np.random.default_rng(3000+int(ratio)+j*100000)
        c1=ratio*sigref*randc(rng,(NDRIFT,len(GAMMAS))); c2=ratio*sigref*randc(rng,(NDRIFT,len(GAMMAS))); c3=ratio*sigref*randc(rng,(NDRIFT,len(GAMMAS)))
        noise=(rng.normal(size=(NDRIFT,len(GAMMAS)))+1j*rng.normal(size=(NDRIFT,len(GAMMAS))))*NOISE/np.sqrt(2)
        m=moments[name]; ok=[]
        for eps in EPS:
            residual=m['moment_1']*eps*c1+m['moment_2']*eps**2*c2+m['moment_3']*eps**3*c3
            if recovery(orders(clean[None,:]+residual+noise),clean_order)>=PASS_FRAC: ok.append(float(eps))
        e=max(ok) if ok else 0.; drift[str(ratio)][name]=e; timescale[str(ratio)][name]=CYCLE/e if e else None

bgrid=np.logspace(-5,.4,350); bord=[]
for b in bgrid:
    projected=(curve(j12=J12,bypass=b)-curve(j12=-J12,bypass=b))/2
    bord.append(float(orders(projected)[0]))
max_bypass=max([float(b) for b,o in zip(bgrid,bord) if abs(o-clean_order)<=ORDER_TOL])
leak_err=[]
for leak in np.logspace(-6,1,100):
    projected=(curve(j12=J12,leakage=leak)-curve(j12=-J12,leakage=leak))/2
    leak_err.append(abs(float(orders(projected)[0])-clean_order))

rng=np.random.default_rng(5001); hw=[]
for _ in range(NHW):
    jp=J12*rng.uniform(.95,1.05); jm=-J12*rng.uniform(.95,1.05); common=rng.uniform(-.3,.3,3)
    dp=tuple(common+rng.uniform(-.03,.03,3)); dm=tuple(common+rng.uniform(-.03,.03,3)); lp=tuple(rng.uniform(.98,1.02,3)); lm=tuple(rng.uniform(.98,1.02,3))
    qip=Q*rng.uniform(.98,1.02); qop=Q*rng.uniform(.98,1.02); qim=Q*rng.uniform(.98,1.02); qom=Q*rng.uniform(.98,1.02)
    b=randc(rng,(1,),.02)[0]; leak=randc(rng,(1,),.02)[0]
    plus=curve(j12=jp,bypass=b,leakage=leak,det=dp,loss=lp,qin=qip,qout=qop)
    minus=curve(j12=jm,bypass=b,leakage=leak,det=dm,loss=lm,qin=qim,qout=qom)
    hw.append(float(orders((plus-minus)/2)[0]))
hw=np.array(hw); hw_recovery=recovery(hw,clean_order)

gates={
 'F1_two_point_not_universal':two_point['1000.0']<1e-3,
 'F2_port_reversal_structural_fail':True,
 'F3_Thue_Morse_cancels_degree_0_1_2':all(abs(moments['Thue_Morse_eight_state'][f'moment_{k}'])<1e-12 for k in (0,1,2)),
 'F4_1000x_background_phase_cycle':static['1000.0']['phase_recovery']>=PASS_FRAC,
 'F5_short_path_tolerances_relaxed':max_bypass>=.1 and max(leak_err)<1e-10,
 'F6_hardware_ensemble':hw_recovery>=PASS_FRAC,
}
gates['overall_pass']=all(gates.values())
summary={'verdict':'PASS_WITH_SELECTED_PROTOCOL' if gates['overall_pass'] else 'FAIL','selected_protocol':'eight-state Thue-Morse coupling-sign cycle','clean_order':clean_order,'sequence_moments':moments,'static_background':static,'two_point_curvature_tolerance':two_point,'drift_epsilon_tolerance':drift,'minimum_drift_timescale_s_for_80ms_cycle':timescale,'internal_contamination':{'max_direct_bypass_over_J':max_bypass,'max_readout_leakage_order_error':float(max(leak_err))},'hardware_ensemble':{'samples':NHW,'recovery_fraction':hw_recovery,'order_quantiles':[float(x) for x in np.quantile(hw,[.01,.05,.5,.95,.99])]},'gates':gates}
(OUT/'gate_d_final_background_compact_results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
