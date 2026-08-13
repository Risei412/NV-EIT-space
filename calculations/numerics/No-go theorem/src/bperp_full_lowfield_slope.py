"""bperp_full_lowfield_slope.py -- Gate: does the B_perp^2 opening law survive
in the genuine full-Liouvillian (Lindblad) NV model, not just the reduced
6x6 non-Hermitian resolvent kernel (bperp_kernel_map_v2.py)?

Runs the same low-field Bx sweep through bperp_full_validation.run_one
(full 9-level Lindbladian, weak_probe_response steady-state + linear
response) and fits the power law of |C_full(Bx)| near Bx=0, then compares
directly against the already-fitted reduced-kernel exponent
(bperp_kernel_map_v2.py: delta_eta_power/delta_C_power/delta_dA_power ~ 2.03).
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bperp_kernel_map_v2 as km
import bperp_full_validation as bfv

def slope(x, y):
    m = (x > 0) & (np.abs(y) > 0)
    return float(np.polyfit(np.log10(x[m]), np.log10(np.abs(y[m])), 1)[0])

def main(out='../results/tables/bperp_full_vs_reduced_slopes.json'):
    Bz = 0.02; T = 300.0
    Bsmall = np.linspace(0, 0.03, 16)
    _, Us = km.track(Bsmall, Bz)
    rows = []
    for Bx, U in zip(Bsmall, Us):
        r = bfv.run_one(T, float(Bx), Bz, U)
        rows.append(r)
    C_full = np.array([r['C_full'] for r in rows])
    C_red = np.array([r['C_red'] for r in rows])
    dC_full = np.abs(C_full - C_full[0])
    dC_red = np.abs(C_red - C_red[0])
    slope_full = slope(Bsmall, dC_full)
    slope_red = slope(Bsmall, dC_red)
    reduced_reference = km.main.__module__  # marker only
    result = dict(T=T, Bz_T=Bz, Bx_grid=Bsmall.tolist(),
                  C_full=C_full.tolist(), C_red=C_red.tolist(),
                  slope_full_liouvillian=slope_full,
                  slope_reduced_kernel_this_run=slope_red,
                  predicted_power_law=2.0,
                  note='reduced-kernel independent cross-check via bperp_kernel_map_v2.main() '
                       'gives delta_C_power~2.026 (see bperp_kernel_summary_v2.json)')
    print(json.dumps({k: v for k, v in result.items() if k not in ('Bx_grid', 'C_full', 'C_red')}, indent=2))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f: json.dump(result, f, indent=2)

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '../results/tables/bperp_full_vs_reduced_slopes.json')
