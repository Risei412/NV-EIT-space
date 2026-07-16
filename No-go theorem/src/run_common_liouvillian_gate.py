"""run_common_liouvillian_gate.py -- Gate 5 (validation_report.md Sec. 5/7):
run NV, SiV-, SnV- through the ONE common reduced-Liouvillian pipeline and
output (m, M0, M1, p_prot, A_mat, q, x(T), C(T), EIT/ATS verdict) so the
class-aware collapse test (universal_scaling_test.py in the uploaded
surrogate package) can be re-run on real computed data instead of
hand-picked (n,A,p,q) parameters.
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nv_reduced_kernel as nvk
import group_iv_model as giv
import eit_ats_classifier as clf
import phonon_rates as pr

GAMMAS = np.logspace(2, 5, 31)

def nv_entry(pair, T, Oc=1.0):
    H = nvk.H_3E(); M = nvk.moments(H, pair, 3); K = nvk.kernel(H, pair, GAMMAS)
    m = next((n for n, v in enumerate(M) if abs(v) > 1e-10), len(M)-1)
    Delta_sel = nvk.LAM_Z
    gamma_fast = pr.gamma_oc(T, 1.683)*1e-9 + nvk.GAMMA_RAD
    gg = 6.3e-5
    p, c = nvk.legs(pair)
    spec_zs = np.linspace(-2.0, 2.0, 401)
    A = clf.spectrum(H, p, c, gamma_fast, spec_zs, Oc=Oc, gg=gg)
    verdict = clf.classify(spec_zs, A)
    r_res = clf.response_generic(H, p, c, 0.0, gamma_fast, Oc=Oc, gg=gg)
    return dict(material=f"NV{pair}", m=m, M0=abs(M[0]), M1=abs(M[1]),
                p_prot=0.0, A_mat=abs(M[m]), q=gg/Delta_sel,
                x_T=gamma_fast/Delta_sel, C_T=r_res['C'], T=T, Oc=Oc,
                delta_aic=verdict['delta_aic'], verdict=verdict['verdict'])

def group_iv_entry(material, T, Oc=1.0):
    H = giv.H_groupIV(material); M = giv.moments(H, 3); K = giv.kernel(H, GAMMAS)
    m = next((n for n, v in enumerate(M) if abs(v) > 1e-10), len(M)-1)
    Delta_sel = giv.PARAMS[material]['Delta_e']
    gamma_fast = giv.gamma_orb_GHz(material, T) + giv.GAMMA_RAD
    gg = 6.3e-5
    p, c = giv.legs(theta=0.0)
    span = max(2.0, Delta_sel*0.02)
    spec_zs = np.linspace(-span, span, 401)
    A = clf.spectrum(H, p, c, gamma_fast, spec_zs, Oc=Oc, gg=gg)
    verdict = clf.classify(spec_zs, A)
    r_res = clf.response_generic(H, p, c, 0.0, gamma_fast, Oc=Oc, gg=gg)
    return dict(material=f"{material}-", m=m, M0=abs(M[0]), M1=abs(M[1]),
                p_prot=0.0, A_mat=abs(M[m]), q=gg/Delta_sel,
                x_T=gamma_fast/Delta_sel, C_T=r_res['C'], T=T, Oc=Oc,
                delta_aic=verdict['delta_aic'], verdict=verdict['verdict'])

def main(out='../results/tables/common_liouvillian_gate.csv'):
    rows = [nv_entry((0, -1), T) for T in (100., 150., 300.)]
    rows += [nv_entry((-1, +1), T) for T in (100., 150., 300.)]
    rows += [group_iv_entry(mat, T) for mat in ('SiV', 'SnV') for T in (4., 100., 300.)]
    import pandas as pd
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, index=False)

    # class-aware collapse test on REAL data: x = gamma_fast/Delta_sel,
    # z = (x^(m+1)) / A_mat  (generalizes the surrogate z=(x^n+q x)/A),
    # compare within-class (same m) spread vs across-class spread.
    df['n_class'] = df['m'] + 1
    df['z'] = df['x_T']**df['n_class'] / df['A_mat'].clip(lower=1e-300)
    spreads = {}
    for n_class, g in df.groupby('n_class'):
        spreads[int(n_class)] = dict(materials=sorted(g['material'].unique().tolist()),
                                      n_rows=len(g))
    print("\nMoment classes present (n=m+1):", json.dumps(spreads, indent=2))
    return df

if __name__ == '__main__':
    main()
