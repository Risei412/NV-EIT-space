"""Reproducible reduced-response bound and implementation audit.

This does NOT certify a maximum EIT observation temperature. It separates
transparency witnesses, a conservative necessary-condition bound, and
implementation failures. Existing certificates are never overwritten.
Run this file with python -B and a quoted path.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import differential_evolution, brentq
import nv_model as nv
import run_prl_prediction as rp
import liouvillian_core as lc

OUT = Path(__file__).resolve().parents[4] / 'results' / 'audits' / 'eit_temperature_limit'


def vectorization_audit():
    H = np.array([[0, 1-2j], [1+2j, .7]], complex)
    rho = np.array([[.7, .1+.2j], [.1-.2j, .3]], complex)
    J = np.array([[0, .3+.1j], [0, 0]], complex)
    Q = J.conj().T @ J
    direct = -1j*(H@rho-rho@H) + J@rho@J.conj().T - (Q@rho+rho@Q)/2
    L = lc.liouvillian(H, [J])
    row = float(np.linalg.norm(L@lc.vec(rho)-lc.vec(direct)))
    col = float(np.linalg.norm(L@rho.ravel(order='F')-direct.ravel(order='F')))
    assert row > 1e-3 and col < 1e-12
    return dict(row_residual=row, column_residual=col,
                conclusion='liouvillian uses column stacking; vec/unvec use row stacking')


def witness(T, seed):
    """Search all six eigenenergy branches and both control ground states.

    Elliptic Jones vectors, fixed strain, exact two-photon resonance.
    A found maximum is a numerical witness, NOT a certified global optimum.
    """
    gamma = nv.gamma_oc_GHz(T, rp.D)
    def evaluate(x):
        Bx, Bz, ap, ac, pp, pc, Oc, jj, channel = x
        j = int(round(jj)); ctrl = '+1' if round(channel) else '-1'
        p = np.array([np.cos(ap), np.exp(1j*pp)*np.sin(ap)])
        c = np.array([np.cos(ac), np.exp(1j*pc)*np.sin(ac)])
        H, dp, dc, w, V, wp, wc = rp.model_at(Bx, Bz, p, c, ctrl)
        r = nv.response(H, dp, dc, w[j], gamma, Oc)
        return r['C'], min(wp[j], wc[j]), j, ctrl
    def loss(x):
        C, weight, _, _ = evaluate(x)
        return -C if weight >= .1 else 10 + 100*(.1-weight)
    bounds = [(0,.5),(.005,.05),(0,np.pi/2),(0,np.pi/2),
              (-np.pi,np.pi),(-np.pi,np.pi),(0,.1),(0,5),(0,1)]
    opt = differential_evolution(loss, bounds, seed=seed, popsize=10,
            maxiter=90, polish=False, integrality=[False]*7+[True,True], tol=1e-6)
    C, weight, branch, ctrl = evaluate(opt.x)
    keys = ['Bx_T','Bz_T','probe_theta','control_theta','probe_phase',
            'control_phase','Oc_GHz','branch','channel_index']
    return dict(T_K=T, C=float(C), minimum_branch_weight=float(weight),
                branch=branch, ctrl=ctrl, seed=seed, nfev=opt.nfev,
                optimizer_converged=bool(opt.success),
                parameters=dict(zip(keys, map(float,opt.x))),
                is_eit_certified=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    audit = vectorization_audit()
    # Weyl bound: span(H0 + Zeeman) <= span(H0) + 2*g_e*|B|.
    span0 = float(np.ptp(np.linalg.eigvalsh(nv.Hes((0,0,0), rp.D, 0))))
    span_bound = span0 + 2*nv.GE*np.hypot(.5,.05)
    R = 2*np.pi*span_bound
    def bound(T):
        q = (R/nv.gamma_oc_GHz(T,rp.D))**2
        return q*(1+q)
    upper_crossing = brentq(lambda T:bound(T)-.01,20,300)
    beta_max = (2*np.pi*.1)**2/4
    geff = 2*6.3e-5+2e-6
    def power_bound(T):
        gamma = nv.gamma_oc_GHz(T,rp.D)
        return beta_max*R**2/(gamma**3*(geff+beta_max*gamma/(gamma**2+R**2)))
    power_crossing = brentq(lambda T:power_bound(T)-.01,20,300)
    fixed = brentq(lambda T:rp.branch_value(T,rp.BX0,rp.BZ0,rp.J0)['C']-.01,60,80)
    result = dict(vectorization=audit, fixed_candidate_T_1pct_K=fixed,
        contrast_floor=.01, strain_GHz=rp.D, strain_azimuth=0,
        spectral_span_bound_GHz=span_bound, R_rate_units=R,
        conservative_contrast_bound_crossing_K=upper_crossing,
        power_limited_contrast_bound_crossing_K=power_crossing,
        bound_formula='C <= q*(1+q); q=(2*pi*span_bound/gamma(T))**2',
        bound_scope='reduced response; orthogonal normalized optical legs; z inside excited eigenvalue range; fixed strain; positive ground damping; given field bounds',
        upper_temperature_monotonicity='fixed-energy Bose integral increases with T; see audit note; model claims restricted to <=300 K',
        highest_EIT_temperature_K=None, witnesses=[])
    path = OUT/'audit.json'
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='witnesses'}),flush=True)
    for T in (70.,80.,90.,100.,120.):
        for seed in (20260920,20260921):
            w = witness(T,seed)
            result['witnesses'].append(w)
            path.write_text(json.dumps(result,indent=2),encoding='utf-8')
            print(json.dumps(w),flush=True)
    print(str(path),flush=True)


if __name__ == '__main__':
    main()
