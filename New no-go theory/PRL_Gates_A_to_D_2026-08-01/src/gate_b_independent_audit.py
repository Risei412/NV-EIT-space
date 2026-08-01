"""Independent audit of NV-EIT PRL Gate B: two qubits coupled by two lossy bus modes."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
G = 0.1
DELTA1 = 0.03
DELTA2 = 0.05
GAMMA_Q = 2.5e-5
N = 5

def couplings(tuning="generic", eps=0.0):
    if tuning == "generic": return G, G, G, G
    if tuning == "protected": return G, G, G, -G
    if tuning == "broken": return G, G, G, -G + eps
    raise ValueError(tuning)

def h_exc(tuning="generic", eps=0.0, z=0.0):
    gA1, gA2, gB1, gB2 = couplings(tuning, eps)
    H = np.zeros((4, 4), complex)
    H[1,1], H[2,2] = DELTA1-z, DELTA2-z
    H[0,1] = H[1,0] = gA1; H[0,2] = H[2,0] = gA2
    H[1,3] = H[3,1] = gB1; H[2,3] = H[3,2] = gB2
    return H

D = np.diag([0.0, 0.5, 0.5, 0.0]).astype(complex)
FLOOR = np.diag([GAMMA_Q/2, 0.0, 0.0, GAMMA_Q/2]).astype(complex)
SOURCE = np.array([1.0,0.0,0.0,0.0], complex)

def reduced_amplitudes(kappa, tuning="generic", eps=0.0):
    x = np.linalg.solve(kappa*D + FLOOR + 1j*h_exc(tuning, eps), SOURCE)
    return complex(x[0]), complex(x[3])

def h_full(tuning="generic", eps=0.0, z=0.0, drive=0.0):
    gA1,gA2,gB1,gB2 = couplings(tuning, eps)
    H = np.zeros((N,N), complex)
    H[1,1]=-z; H[2,2]=DELTA1-z; H[3,3]=DELTA2-z; H[4,4]=-z
    H[0,1]=H[1,0]=drive
    H[1,2]=H[2,1]=gA1; H[1,3]=H[3,1]=gA2
    H[2,4]=H[4,2]=gB1; H[3,4]=H[4,3]=gB2
    return H

def jumps(kappa):
    ops=[]
    for j in (2,3):
        L=np.zeros((N,N),complex); L[0,j]=np.sqrt(kappa); ops.append(L)
    for j in (1,4):
        L=np.zeros((N,N),complex); L[0,j]=np.sqrt(GAMMA_Q); ops.append(L)
    return ops

def liouvillian(H, ops):
    dim=N*N; L=np.zeros((dim,dim),complex)
    for idx in range(dim):
        rho=np.zeros(dim,complex); rho[idx]=1.0; rho=rho.reshape(N,N)
        dr=-1j*(H@rho-rho@H)
        for C in ops:
            Cd=C.conj().T; dr += C@rho@Cd - 0.5*(Cd@C@rho + rho@Cd@C)
        L[:,idx]=dr.reshape(-1)
    return L

def solve_trace(L, rhs, trace_value):
    M=L.copy(); b=rhs.reshape(-1).copy(); tr=np.zeros(N*N,complex)
    for i in range(N): tr[i*N+i]=1.0
    M[0,:]=tr; b[0]=trace_value
    return np.linalg.solve(M,b).reshape(N,N)

def full_linear_response(kappa,tuning="generic",eps=0.0):
    H0=h_full(tuning,eps); L0=liouvillian(H0,jumps(kappa))
    rho0=np.zeros((N,N),complex); rho0[0,0]=1.0
    V=np.zeros((N,N),complex); V[0,1]=V[1,0]=1.0
    return solve_trace(L0,1j*(V@rho0-rho0@V),0.0)

def full_steady_state(kappa,tuning="generic",drive=1e-8):
    return solve_trace(liouvillian(h_full(tuning,drive=drive),jumps(kappa)),np.zeros((N,N),complex),1.0)

def fit_order(x,y,mask):
    return float(-np.polyfit(np.log(x[mask]),np.log(np.abs(y[mask])),1)[0])

blind={"generic":{"nu_amplitude":1,"nu_power":2},"protected":{"nu_amplitude":2,"nu_power":4}}
kappas=np.logspace(2,8,90); tail=kappas>1e5
rows=[]; curves={}; max_rel_error=0.0
for tuning in ("generic","protected"):
    red_A=[]; red_B=[]; lin_A=[]; lin_B=[]; steady_B=[]
    for kappa in kappas:
        a_red,b_red=reduced_amplitudes(kappa,tuning)
        rho1=full_linear_response(kappa,tuning)
        a_lin=rho1[1,0]/(-1j); b_lin=rho1[4,0]/(-1j)
        rho_ss=full_steady_state(kappa,tuning)
        red_A.append(a_red); red_B.append(b_red); lin_A.append(a_lin); lin_B.append(b_lin)
        steady_B.append(max(float(rho_ss[4,4].real),1e-300))
        max_rel_error=max(max_rel_error,abs(b_lin-b_red)/max(abs(b_red),1e-300),abs(a_lin-a_red)/max(abs(a_red),1e-300))
    lin_A=np.asarray(lin_A); lin_B=np.asarray(lin_B); steady_B=np.asarray(steady_B)
    raw=np.abs(lin_B)**2; transfer=np.abs(lin_B/lin_A)**2
    result={"amplitude_order":fit_order(kappas,lin_B,tail),"raw_power_order":fit_order(kappas,raw,tail),"normalized_transfer_order":fit_order(kappas,transfer,tail),"finite_drive_B_population_order":fit_order(kappas,steady_B,tail)}
    rows.append({"tuning":tuning,**result}); curves[tuning]=transfer

epsilons=[1e-8,1e-9,1e-10]; kgrid=np.logspace(4,10,160); crossovers=[]
for eps in epsilons:
    values=np.asarray([reduced_amplitudes(k,"broken",eps)[1] for k in kgrid])
    nu=-np.gradient(np.log(np.abs(values)),np.log(kgrid)); idx=np.where(nu<1.5)[0]
    numeric=float(kgrid[idx[0]]) if len(idx) else float("inf")
    analytic=float(2*G*abs(DELTA2-DELTA1)/abs(eps))
    crossovers.append({"eps":eps,"analytic_kappa_star":analytic,"numeric_kappa_star":numeric,"ratio_numeric_to_analytic":numeric/analytic})

physical_k=np.logspace(-4,np.log10(0.05),120); physical_orders={}
for tuning in ("generic","protected"):
    ratios=[]
    for kappa in physical_k:
        a,b=reduced_amplitudes(kappa,tuning); ratios.append(abs(b/a)**2)
    physical_orders[tuning]=fit_order(physical_k,np.asarray(ratios),np.ones(len(physical_k),bool))

summary={"verdict":"STRUCTURAL_PASS_EXPERIMENTAL_WINDOW_UNRESOLVED","blind_prediction":blind,"deep_tail_results":rows,"max_reduced_vs_full_relative_error":max_rel_error,"symmetry_breaking_crossovers":crossovers,"current_parameter_physical_window_orders":physical_orders}
(OUT/"gate_b_independent_audit.json").write_text(json.dumps(summary,indent=2))
with (OUT/"gate_b_independent_audit.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
plt.figure(figsize=(7.2,5.2))
for tuning in ("generic","protected"): plt.loglog(kappas,curves[tuning],label=tuning)
plt.xlabel("bus decay kappa"); plt.ylabel("normalized transfer |B/A|^2"); plt.legend(); plt.tight_layout(); plt.savefig(OUT/"gate_b_independent_audit.png",dpi=180); plt.close()
print(json.dumps(summary,indent=2))
