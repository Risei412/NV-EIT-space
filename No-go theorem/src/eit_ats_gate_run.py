"""eit_ats_gate_run.py -- Gate 3 driver: (a) sanity-check the AIC classifier
on toy pure-EIT / pure-ATS lineshapes, (b) apply it to NV and group-IV
(SiV-, SnV-) absorption spectra generated from the SAME reduced-kernel
Hamiltonians used in moment_order_common_pipeline.py, connecting delta_chi_S
all the way to a verdict (per SiV_SnV_phonon_AIC_parameters.md Sec. 8-11).
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nv_reduced_kernel as nvk
import group_iv_model as giv
import eit_ats_classifier as clf

def sanity_check():
    d = np.linspace(-20, 20, 400)
    A_eit_toy = clf.A_EIT(d, 5.0, 1.0, 4.7, 0.9) + 0.002*np.random.default_rng(0).standard_normal(d.size)
    r_eit = clf.classify(d, A_eit_toy)
    A_ats_toy = clf.A_ATS(d, 3.0, 1.0, 8.0) + 0.002*np.random.default_rng(1).standard_normal(d.size)
    r_ats = clf.classify(d, A_ats_toy)
    print("sanity: toy EIT ->", r_eit['verdict'], f"(dAIC={r_eit['delta_aic']:.2f})")
    print("sanity: toy ATS ->", r_ats['verdict'], f"(dAIC={r_ats['delta_aic']:.2f})")
    assert r_eit['verdict'] == 'robust EIT', "classifier failed toy EIT sanity check"
    assert r_ats['verdict'] == 'robust ATS', "classifier failed toy ATS sanity check"
    return dict(toy_eit=r_eit, toy_ats=r_ats)

def nv_row(T, pair=(0, -1), Oc=1.0):
    import phonon_rates as pr
    H = nvk.H_3E()
    p, c = nvk.legs(pair)
    gamma = pr.gamma_oc(T, 1.683)*1e-9 + nvk.GAMMA_RAD
    zs = np.linspace(-2.0, 2.0, 401)
    A = clf.spectrum(H, p, c, gamma, zs, Oc=Oc, gg=6.3e-5)
    r = clf.classify(zs, A)
    r.update(system=f"NV {pair}", T=T, Oc=Oc, gamma_GHz=gamma)
    return r

def group_iv_row(material, T, Oc=1.0):
    H = giv.H_groupIV(material)
    p, c = giv.legs(theta=0.0)
    gamma = giv.gamma_orb_GHz(material, T) + giv.GAMMA_RAD
    span = max(2.0, giv.PARAMS[material]['Delta_e']*0.02)
    zs = np.linspace(-span, span, 401)
    A = clf.spectrum(H, p, c, gamma, zs, Oc=Oc, gg=6.3e-5)
    r = clf.classify(zs, A)
    r.update(system=f"{material}- orbital-Lambda", T=T, Oc=Oc, gamma_GHz=gamma,
              k_B_T_over_h_Delta_e=giv.thermal_regime(material, T))
    return r

def main(out='../results/tables/eit_ats_gate.json'):
    sanity = sanity_check()
    rows = []
    for T in (100.0, 150.0, 300.0):
        for Oc in (0.05, 1.0, 20.0):
            rows.append(nv_row(T, Oc=Oc))
    for mat in ('SiV', 'SnV'):
        for T in (4.0, 100.0, 300.0):
            for Oc in (0.05, 1.0, 20.0):
                rows.append(group_iv_row(mat, T, Oc=Oc))
    print(f"\n{'system':22s} {'T':>6s} {'Oc':>6s} {'C':>12s} {'dAIC':>9s}  verdict")
    for r in rows:
        print(f"{r['system']:22s} {r['T']:6.1f} {r['Oc']:6.2f} "
              f"{'-':>12s} {r['delta_aic']:9.2f}  {r['verdict']}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(dict(sanity_check_passed=True, rows=rows), f, indent=2, default=str)
    return sanity, rows

if __name__ == '__main__':
    main()
