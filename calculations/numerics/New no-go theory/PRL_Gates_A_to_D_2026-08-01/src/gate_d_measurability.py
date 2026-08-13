"""Compact rerunner for NV-EIT PRL Gate D measurability audit."""
from __future__ import annotations
from pathlib import Path
import json, math
import numpy as np

OUT = Path(__file__).resolve().parent
J12, J23 = 1.0, 0.8
GAMMAS = np.logspace(np.log10(5.0), np.log10(15.0), 9)
J12_MHZ = 3.0
PIN_DBM = -110.0
PIN_W = 1e-3 * 10 ** (PIN_DBM / 10.0)
TNOISE, TARGET_SNR, TINT, KB = 3.0, 10.0, 40.0, 1.380649e-23
CRITERIA = dict(global_min=2.90, global_max=3.10, local_min=2.80, local_max=3.20)


def transmission(gamma, qin, qout, *, j12=J12, j23=J23,
                 det=(0.0, 0.0, 0.0), loss=(1.0, 1.0, 1.0),
                 bypass=0j, leakage=0j):
    H = np.array([[det[0], j12, bypass],
                  [j12, det[1], j23],
                  [np.conj(bypass), j23, det[2]]], complex)
    A = np.diag([gamma*loss[0]+qin/2, gamma*loss[1], gamma*loss[2]+qout/2]).astype(complex)
    A += 1j*H
    c = np.array([1.0, 0.0, 0.0], complex)
    p = np.array([leakage, 0.0, 1.0], complex)
    return complex(-np.sqrt(qin*qout) * np.vdot(p, np.linalg.solve(A, c)))


def metrics(qin, qout, **kwargs):
    vals = np.array([transmission(g, qin, qout, **kwargs) for g in GAMMAS])
    x, y = np.log(GAMMAS), np.log(np.abs(vals))
    global_order = float(-np.polyfit(x, y, 1)[0])
    local = -np.diff(y)/np.diff(x)
    return vals, global_order, local


def monte_carlo_slopes(values, seed, n=10_000):
    rng = np.random.default_rng(seed)
    pout = PIN_W*np.abs(values)**2
    rel = np.sqrt(KB*TNOISE/np.maximum(pout*TINT, 1e-300))
    out = np.empty(n)
    for i in range(n):
        z = (rng.normal(size=len(values))+1j*rng.normal(size=len(values)))/np.sqrt(2)
        measured = values + z*rel*np.abs(values)
        out[i] = -np.polyfit(np.log(GAMMAS), np.log(np.abs(measured)), 1)[0]
    return out

candidates=[]
for q in np.logspace(-3,1,2001):
    vals, nu, local = metrics(q,q)
    if CRITERIA['global_min'] <= nu <= CRITERIA['global_max'] and local.min() >= CRITERIA['local_min'] and local.max() <= CRITERIA['local_max']:
        candidates.append((float(abs(vals[-1])**2), float(q), nu, float(local.min()), float(local.max())))
if not candidates:
    raise RuntimeError('No admissible port coupling.')
_, Q, nu3, local_min, local_max = max(candidates)
clean, _, _ = metrics(Q,Q)
n2, nu2, _ = metrics(Q,Q,bypass=0.2)

power = np.abs(clean)**2
pout = PIN_W*power
required_t = TARGET_SNR**2*KB*TNOISE/np.maximum(pout,1e-300)
snr = np.sqrt(pout*TINT/(KB*TNOISE))
slopes3 = monte_carlo_slopes(clean, 412)
slopes2 = monte_carlo_slopes(n2, 413)

rng=np.random.default_rng(414); ensemble=[]
for _ in range(5000):
    qin=Q*rng.uniform(.9,1.1); qout=Q*rng.uniform(.9,1.1)
    det=tuple(rng.uniform(-.5,.5,3)); loss=tuple(rng.uniform(.9,1.1,3))
    _,nu,local=metrics(qin,qout,j12=rng.uniform(.9,1.1),j23=J23*rng.uniform(.9,1.1),det=det,loss=loss)
    ensemble.append((nu, bool(2.8<=nu<=3.2 and local.min()>=2.5 and local.max()<=3.5)))
pass_fraction=float(np.mean([p for _,p in ensemble]))

phases=np.linspace(0,2*np.pi,73,endpoint=False)
def acceptable(kind,mag):
    for phase in phases:
        v=mag*np.exp(1j*phase)
        kw={'bypass':v} if kind=='bypass' else {'leakage':v}
        vals,nu,local=metrics(Q,Q,**kw)
        if not (2.8<=nu<=3.2 and local.min()>=2.5 and local.max()<=3.5 and np.all(np.diff(np.abs(vals))<0)):
            return False
    return True

def largest(kind,grid):
    ok=[float(v) for v in grid if acceptable(kind,float(v))]
    return max(ok) if ok else 0.0

max_bypass=largest('bypass',np.logspace(-6,-1,101))
max_leak=largest('leakage',np.logspace(-7,-1,121))
summary={
 'verdict':'CONDITIONAL_PASS',
 'design':{'Gamma_over_J':[5.0,15.0],'port_energy_decay_over_J':Q,'J12_over_2pi_MHz':J12_MHZ,'J23_over_2pi_MHz':J23*J12_MHZ},
 'clean_n3':{'global_order':nu3,'local_min':local_min,'local_max':local_max,'worst_power_transmission':float(power[-1]),'worst_output_dBm':float(10*np.log10(pout[-1]/1e-3)),'required_s_for_SNR10':float(required_t[-1]),'SNR_at_40s':float(snr[-1]),'slope_95pct':[float(np.quantile(slopes3,.025)),float(np.quantile(slopes3,.975))]},
 'comparison_n2':{'global_order':nu2,'slope_95pct':[float(np.quantile(slopes2,.025)),float(np.quantile(slopes2,.975))]},
 'fabrication':{'samples':5000,'pass_fraction':pass_fraction,'order_quantiles':[float(x) for x in np.quantile([x for x,_ in ensemble],[.01,.05,.5,.95,.99])]},
 'tolerances':{'direct_bypass_over_J':max_bypass,'direct_bypass_over_2pi_kHz':max_bypass*J12_MHZ*1000,'readout_leakage':max_leak,'readout_leakage_dB':20*math.log10(max_leak)},
}
(OUT/'gate_d_measurability_compact_results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
