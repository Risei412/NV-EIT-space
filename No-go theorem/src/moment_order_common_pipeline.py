"""moment_order_common_pipeline.py -- Gate 1: does the Raman-transfer moment
order (first nonzero moment m, K ~ Gamma^-(m+1)) actually differ between NV
(spin-Lambda) and group-IV (SiV-/SnV-, same-spin orbital-Lambda) in the SAME
computational pipeline (identical A_Gamma=Gamma*D+A0, D=I, moments/kernel
machinery; see nv_reduced_kernel.py / group_iv_model.py)?
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nv_reduced_kernel as nvk
import group_iv_model as giv

GAMMAS = np.logspace(2, 5, 31)

def fit_slope(Gammas, K):
    return float(np.polyfit(np.log10(Gammas[-15:]), np.log10(np.abs(K[-15:])), 1)[0])

def main(out_json='../results/tables/moment_order_common_pipeline.json',
         out_fig='../results/figures/moment_order_common_pipeline.png'):
    rows = []
    H_nv = nvk.H_3E()
    for label, pair in [("NV (0,-1)", (0, -1)), ("NV (0,+1)", (0, +1)),
                        ("NV (-1,+1)", (-1, +1))]:
        M = nvk.moments(H_nv, pair, 3)
        K = nvk.kernel(H_nv, pair, GAMMAS)
        rows.append(dict(system=label, M0=complex(M[0]), M1=complex(M[1]),
                          slope=fit_slope(GAMMAS, K), K_abs=np.abs(K).tolist()))
    for mat in ('SiV', 'SnV'):
        H = giv.H_groupIV(mat)
        M = giv.moments(H, 3)
        K = giv.kernel(H, GAMMAS)
        rows.append(dict(system=f"{mat}- (orbital-Lambda)", M0=complex(M[0]),
                          M1=complex(M[1]), slope=fit_slope(GAMMAS, K),
                          K_abs=np.abs(K).tolist()))

    print(f"{'system':18s} {'|M0|':>10s} {'slope':>8s} {'predicted':>10s}")
    for r in rows:
        m = abs(r['M0'])
        pred = -1.0 if m > 1e-8 else (-2.0 if abs(complex(r['M1'])) > 1e-8 else -3.0)
        print(f"{r['system']:18s} {m:10.4g} {r['slope']:8.4f} {pred:10.1f}")
        r['predicted_slope'] = pred

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.figure(figsize=(7.2, 5.2))
    for r in rows:
        plt.loglog(GAMMAS, r['K_abs'], 'o-', ms=3,
                    label=f"{r['system']}: slope={r['slope']:.2f}")
    plt.xlabel(r'$\Gamma$ (GHz)'); plt.ylabel(r'$|K(\Gamma)|$')
    plt.title('NV vs group-IV moment order, common reduced-kernel pipeline')
    plt.legend(fontsize=8); plt.tight_layout()
    os.makedirs(os.path.dirname(out_fig), exist_ok=True)
    plt.savefig(out_fig, dpi=200)

    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    with open(out_json, 'w') as f:
        json.dump([{k: (str(v) if isinstance(v, complex) else v) for k, v in r.items()}
                    for r in rows], f, indent=2)
    return rows

if __name__ == '__main__':
    main()
