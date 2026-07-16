"""gate2_candidate_full_vs_reduced.py -- Gate 2 of SIMULATION_PLAN.md.

Spectrum-level comparison of the reduced kernel (delta_chi_S, gate1
probe_scan) against the full Lindblad master equation at the fixed PRL
candidate point (T=70 K, Bx=BX0, Bz=BZ0, branch J0, Oc=0.1 GHz, Y/Y
polarization, control on ms=+1).

Full model: generalization of the validated 9-level pipeline of
bperp_full_validation.py (same Lindblad structure: orbital hopping at
k_orb, spin-conserving radiative decay GRAD expressed in the dressed
ground basis, symmetric ground T1=1e-6 GHz, ground dephasing 2*gg) to
  * arbitrary control channel / polarizations (candidate: Y/Y, ctrl=+1),
  * a genuine probe scan: excited block shifted by -(z0+d2) AND control
    ground by -d2 (bperp_full_validation scanned only the Raman axis),
  * optional ISC + singlet shelving (10th level; Gamma_ISC(ms=+-1)=0.098,
    Gamma_ISC(0)=0.011 GHz per Goldman/Robledo-range values, singlet
    lifetime 371 ns with branching 0.6/0.2/0.2 into ms=0/+-1),
  * optional 14N hyperfine (secular A_par: ground -2.16 MHz, excited
    -40 MHz); mI is conserved by every jump operator used here, so the
    hyperfine ensemble is the exact average of three shifted models.

Observable: chi from the probe-ground->excited coherence of the
first-order response (weak_probe_response conventions); the EIT part is
isolated by cutting the ground-coherence sector (rho_pc, rho_cp) out of
the Liouvillian: C(d2) = (Im chi_cut - Im chi_full)/Im chi_cut.

Outputs: results/tables/gate2_candidate_comparison.csv,
         results/figures/gate2_full_vs_reduced_spectrum.png,
         results/tables/gate2_candidate_comparison.json
Usage: python gate2_candidate_full_vs_reduced.py [--quick]
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE))
import run_prl_prediction as rp
import nv_model as nv
import gate1_candidate_aic_bootstrap as g1
from liouvillian_core import liouvillian, steady_state, first_order, vec

TWOPI = 2*np.pi
GG = 6.3e-5
T1_GROUND = 1e-6
OP = 1e-5                       # perturbative probe Rabi (GHz)
# ISC / singlet (order-of-magnitude literature values; see module docstring)
K_ISC = {0: 0.011, 1: 0.098, -1: 0.098}     # GHz, bare-ms resolved, per orbital
K_SINGLET = 1e-9/371e-9                      # GHz (371 ns lifetime)
P_SINGLET = {0: 0.6, 1: 0.2, -1: 0.2}        # branching into bare ms
A_GS = -2.16e-3                              # GHz, 14N A_par ground
A_ES = -40e-3                                # GHz, 14N A_par excited
MS_OF = [-1, 0, 1]                           # bare spin label of index m

def dressed_from(Hg):
    """nv_model.dressed_ground for an explicit ground Hamiltonian."""
    w, U = np.linalg.eigh(Hg)
    ov = np.abs(U.conj().T@np.eye(3))**2
    r, c = linear_sum_assignment(-ov)
    idx = np.empty(3, dtype=int)
    for rr, cc in zip(r, c): idx[cc] = rr
    U = U[:, idx]; w = w[idx]
    for j in range(3):
        q = int(np.argmax(np.abs(U[:, j]))); z = U[q, j]
        if abs(z): U[:, j] *= np.exp(-1j*np.angle(z))
    return w, U

def build_full(T, Bx, Bz, d2, Oc=rp.OC, ctrl='+1', ppol=rp.Y, cpol=rp.Y,
               isc=False, mI=None, d=rp.D, phi=rp.PHI, j0=rp.J0):
    """Frame: probe laser fixed at the mI=0 candidate resonance z0+d2,
    control laser fixed at z0; hyperfine shifts enter the detunings."""
    N = 10 if isc else 9
    Bvec = (float(Bx), 0.0, float(Bz))
    Sz = nv.Sz
    Hg0 = nv.Hgs(Bvec)
    He0 = nv.Hes(Bvec, d, phi)
    eg0, U0 = dressed_from(Hg0)
    z0 = float(np.linalg.eigvalsh(He0)[j0])
    if mI is None:
        Hg, He, eg, U = Hg0, He0, eg0, U0
    else:
        Hg = Hg0 + mI*A_GS*Sz
        He = He0 + mI*A_ES*np.kron(nv.I2, Sz)
        eg, U = dressed_from(Hg)
    p_idx, c_idx = 1, (2 if ctrl == '+1' else 0)
    s_idx = 3 - p_idx - c_idx
    dp = np.kron(np.asarray(ppol, complex), U[:, p_idx])
    dc = np.kron(np.asarray(cpol, complex), U[:, c_idx])
    H = np.zeros((N, N), complex)
    # probe-scan rotating frame (probe ground at 0)
    H[3:9, 3:9] = He - (z0 + d2)*np.eye(6) + (eg[p_idx] - eg0[p_idx])*np.eye(6)
    H[c_idx, c_idx] = -d2 + (eg[c_idx]-eg0[c_idx]) - (eg[p_idx]-eg0[p_idx])
    H[s_idx, s_idx] = float(eg[s_idx] - eg[p_idx])   # spectator, uncoupled
    Vc = np.zeros((N, N), complex); Vc[3:9, c_idx] = dc
    Vp = np.zeros((N, N), complex); Vp[3:9, p_idx] = dp
    H += 0.5*Oc*(Vc + Vc.conj().T)
    rate = nv.korb_GHz(T, d)
    Ls = []
    for m in range(3):                       # orbital hopping X<->Y
        up = np.zeros((N, N), complex); dn = np.zeros((N, N), complex)
        up[6+m, 3+m] = 1; dn[3+m, 6+m] = 1
        Ls += [np.sqrt(rate)*up, np.sqrt(rate)*dn]
    for orb0 in (3, 6):                      # radiative, spin-conserving
        for m in range(3):
            J = np.zeros((N, N), complex)
            for a in range(3): J[a, orb0+m] = np.conj(U[m, a])
            Ls.append(np.sqrt(nv.GRAD)*J)
    for a in range(3):                       # ground T1 + dephasing
        for b in range(3):
            if a != b:
                J = np.zeros((N, N), complex); J[b, a] = 1
                Ls.append(np.sqrt(T1_GROUND)*J)
    for a in range(3):
        J = np.zeros((N, N), complex); J[a, a] = 1
        Ls.append(np.sqrt(2*GG)*J)
    if isc:
        for orb0 in (3, 6):                  # 3E -> singlet, bare-ms resolved
            for m in range(3):
                J = np.zeros((N, N), complex); J[9, orb0+m] = 1
                Ls.append(np.sqrt(K_ISC[MS_OF[m]])*J)
        for m in range(3):                   # singlet -> dressed grounds
            J = np.zeros((N, N), complex)
            for a in range(3): J[a, 9] = np.conj(U[m, a])
            Ls.append(np.sqrt(K_SINGLET*P_SINGLET[MS_OF[m]])*J)
    return TWOPI*H, Ls, Vp, dp, dict(N=N, p_idx=p_idx, c_idx=c_idx, z0=z0)

def chi_pair(T, Bx, Bz, d2, **kw):
    """(chi_full, chi_cut, diagnostics) at one two-photon detuning."""
    H, Ls, Vp, dp, meta = build_full(T, Bx, Bz, d2, **kw)
    N, p_idx, c_idx = meta['N'], meta['p_idx'], meta['c_idx']
    L = liouvillian(H, Ls)
    rho0, res0, gap = steady_state(L)
    Hp = TWOPI*0.5*OP*(Vp + Vp.conj().T)
    I = np.eye(N)
    V = -1j*(np.kron(I, Hp) - np.kron(Hp.T, I))
    det = np.zeros(N*N, complex)
    for e, a in enumerate(dp): det[p_idx*N + (3+e)] = np.conj(a)
    S = [c_idx*N + p_idx, p_idx*N + c_idx]
    X = [k for k in range(N*N) if k not in S]
    Lc = L.copy(); Lc[np.ix_(S, X)] = 0; Lc[np.ix_(X, S)] = 0
    xf, rf = first_order(L, V, rho0)
    xc, rc = first_order(Lc, V, rho0)
    chif = -2*(det.conj()@xf)/OP
    chic = -2*(det.conj()@xc)/OP
    rho0m = rho0.reshape(N, N)
    diag = dict(res0=res0, res_full=rf, res_cut=rc, gap=gap,
                exc_pop=float(np.real(np.trace(rho0m[3:9, 3:9]))),
                singlet_pop=float(np.real(rho0m[9, 9])) if N == 10 else 0.0,
                rho_cp_over_Op=float(abs(xf[c_idx*N + p_idx])/OP))
    return complex(chif), complex(chic), diag

def full_spectrum(d2s_GHz, T=70.0, Bx=rp.BX0, Bz=rp.BZ0, isc=False,
                  hyperfine=False, **kw):
    """C(d2) and absorption spectra; hyperfine = exact average over mI."""
    mIs = (-1, 0, 1) if hyperfine else (None,)
    chif = np.zeros(len(d2s_GHz), complex)
    chic = np.zeros(len(d2s_GHz), complex)
    diags = []
    for mI in mIs:
        for i, d2 in enumerate(d2s_GHz):
            f, c, g = chi_pair(T, Bx, Bz, float(d2), isc=isc, mI=mI, **kw)
            chif[i] += f/len(mIs); chic[i] += c/len(mIs)
            if i == len(d2s_GHz)//2: diags.append(g)
    Af = np.imag(chif); Ac = np.imag(chic)
    C = (Ac - Af)/Ac
    return dict(d2_MHz=np.asarray(d2s_GHz)*1e3, A_full=Af, A_cut=Ac, C=C,
                diag=diags[-1])

def dip_metrics(d_MHz, C):
    """Peak contrast, center and FWHM of the C(d2) transparency feature."""
    i0 = int(np.argmax(np.abs(C)))
    cmax = float(C[i0]); half = abs(cmax)/2
    idx = np.where(np.abs(C) >= half)[0]
    fwhm = float(d_MHz[idx[-1]] - d_MHz[idx[0]]) if len(idx) > 1 else np.nan
    return dict(Cmax=cmax, center_MHz=float(d_MHz[i0]), fwhm_MHz=fwhm)

def main(quick=False):
    tabdir = ROOT/'results'/'tables'; figdir = ROOT/'results'/'figures'
    tabdir.mkdir(parents=True, exist_ok=True); figdir.mkdir(parents=True, exist_ok=True)
    n = 61 if quick else 241
    half = 0.005                                  # +/-5 MHz around the dip
    d2s = np.linspace(-half, half, n)

    # reduced model on the identical grid (gate1 conventions)
    dr, Ar, Acutr = g1.probe_scan(half_GHz=half, n=n)
    Cr = (Acutr - Ar)/Acutr
    red = dip_metrics(dr, Cr)

    configs = [('base', dict()),
               ('isc', dict(isc=True)),
               ('hyperfine', dict(hyperfine=True)),
               ('isc+hyperfine', dict(isc=True, hyperfine=True))]
    if quick: configs = configs[:2]
    rows, spectra = [], {}
    for name, kw in configs:
        sp = full_spectrum(d2s, **kw)
        m = dip_metrics(sp['d2_MHz'], sp['C'])
        rows.append(dict(config=name,
                         C_full=m['Cmax'], C_red=red['Cmax'],
                         ratio_full_red=m['Cmax']/red['Cmax'],
                         fwhm_full_MHz=m['fwhm_MHz'], fwhm_red_MHz=red['fwhm_MHz'],
                         fwhm_ratio=m['fwhm_MHz']/red['fwhm_MHz'],
                         center_full_MHz=m['center_MHz'], center_red_MHz=red['center_MHz'],
                         sign_preserved=bool(np.sign(m['Cmax']) == np.sign(red['Cmax'])),
                         exc_pop=sp['diag']['exc_pop'],
                         singlet_pop=sp['diag']['singlet_pop'],
                         rho_cp_over_Op=sp['diag']['rho_cp_over_Op'],
                         res_steady=sp['diag']['res0'], res_full=sp['diag']['res_full']))
        spectra[name] = sp
        print(f"{name:14s} C_full={m['Cmax']:.4g} C_red={red['Cmax']:.4g} "
              f"ratio={rows[-1]['ratio_full_red']:.3f} "
              f"fwhm={m['fwhm_MHz']:.3f}/{red['fwhm_MHz']:.3f} MHz "
              f"sign_ok={rows[-1]['sign_preserved']}")

    passed = dict(
        contrast_same_order=all(0.1 <= r['ratio_full_red'] <= 10 for r in rows),
        linewidth_same_order=all(0.1 <= r['fwhm_ratio'] <= 10 for r in rows),
        no_sign_flip=all(r['sign_preserved'] for r in rows),
        survives_all_toggles=all(abs(r['C_full']) > 0.1*abs(red['Cmax'])
                                 for r in rows),
    )
    with (tabdir/'gate2_candidate_comparison.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with (tabdir/'gate2_candidate_comparison.json').open('w') as fh:
        json.dump(dict(candidate=dict(T_K=70.0, Bx_T=rp.BX0, Bz_T=rp.BZ0,
                                      Oc_GHz=rp.OC, branch=rp.J0),
                       reduced=red, rows=rows, gate2_pass=passed,
                       gate2_all_passed=all(passed.values()), quick=quick),
                  fh, indent=2)

    plt.figure(figsize=(7.6, 5.0))
    plt.plot(dr, Cr, 'k-', lw=2, label='reduced kernel δχ_S')
    for name, sp in spectra.items():
        plt.plot(sp['d2_MHz'], sp['C'], '--', label=f'full Lindblad ({name})')
    plt.xlabel('two-photon detuning (MHz)')
    plt.ylabel('EIT contrast C(δ₂)')
    plt.title('Gate 2: full Liouvillian vs reduced kernel, 70 K candidate')
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(figdir/'gate2_full_vs_reduced_spectrum.png', dpi=220)
    plt.close()
    print('gate2_pass:', passed, '->', all(passed.values()))
    return rows, passed

if __name__ == '__main__':
    main(quick='--quick' in sys.argv[1:])
