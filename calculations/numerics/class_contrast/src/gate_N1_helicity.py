"""
Gate N1 - helicity asymmetry: is the (0,-1)/(0,+1) chirality null robust?

Background. verify_nv_3E_graph_distance_PRL.py test T6 shows that lam_perp
alone opens the ms=0<->+1 leg only, leaving ms=0<->-1 exactly zero, and calls
this a falsifiable polarization/helicity null test. This gate asks whether the
null survives (a) a finite nonsecular spin-spin term Delta2, (b) strain, (c)
the transverse field required to sit inside the transparency island, and (d) a
realistic ensemble strain distribution.

Observable:  A = (|K(0,+1)| - |K(0,-1)|) / (|K(0,+1)| + |K(0,-1)|)
"""
import csv
from pathlib import Path
import numpy as np
from nv3e_loader import load

m = load()
H_3E, kernel, moments = m.H_3E, m.kernel, m.moments
LAM_PERP, DELTA2 = m.LAM_PERP, m.DELTA2

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "No-go theorem" / "src"))
import phonon_rates as pr

GE = 28.0              # GHz/T, electron gyromagnetic ratio
D_STRAIN = 1.683       # GHz, manuscript transverse strain
OUT = Path(__file__).resolve().parents[1] / "out"


def gamma_oc_GHz(T, d=D_STRAIN):
    return pr.gamma_oc(T, d) / 1e9          # phonon_rates returns Hz


def asym_moment(lam_perp, delta2, **kw):
    """Asymmetry of the leading moment M1 (asymptotic Gamma -> inf)."""
    a = abs(moments(H_3E(lam_perp=lam_perp, delta2=delta2, **kw), (0, +1), 1)[1])
    b = abs(moments(H_3E(lam_perp=lam_perp, delta2=delta2, **kw), (0, -1), 1)[1])
    return (a - b) / (a + b)


def asym_kernel(T, B_T, xi_x=0.0, xi_y=0.0, lam_perp=LAM_PERP, delta2=DELTA2):
    """Asymmetry of the finite-Gamma kernel at temperature T and field B_T."""
    H = H_3E(lam_perp=lam_perp, delta2=delta2, xi_x=xi_x, xi_y=xi_y, Bx=GE * B_T)
    G = gamma_oc_GHz(T)
    a = abs(kernel(H, (0, +1), [G])[0])
    b = abs(kernel(H, (0, -1), [G])[0])
    return (a - b) / (a + b)


def main():
    OUT.mkdir(exist_ok=True)
    print("=" * 78)
    print("N1-1  closed form of the asymptotic asymmetry")
    print("=" * 78)
    worst = 0.0
    for lp in (0.05, 0.20, 0.50):
        for d2 in (0.0, 0.10, 0.20, 0.60):
            worst = max(worst, abs(asym_moment(lp, d2) - lp / (lp + d2)))
    print(f"  A_inf = lam_perp/(lam_perp+Delta2)   max deviation = {worst:.2e}")
    print(f"  literature lam_perp={LAM_PERP}, Delta2={DELTA2}"
          f"  ->  A_inf = {LAM_PERP/(LAM_PERP+DELTA2):.4f}")

    print()
    print("=" * 78)
    print("N1-2  at zero field the moment-level null is protected")
    print("=" * 78)
    print(f"  {'perturbation':<28} {'A':>10}")
    for label, kw in [("none", {}),
                      ("strain xi_x = 5 GHz", dict(xi_x=5.0)),
                      ("strain xi_y = 5 GHz", dict(xi_y=5.0)),
                      ("strain xi = 20 GHz", dict(xi_x=20.0, xi_y=14.0)),
                      ("transverse field 3 GHz", dict(Bx=3.0))]:
        print(f"  {label:<28} {asym_moment(LAM_PERP, DELTA2, **kw):10.6f}")
    print("  M1 = -i p^dag H c samples only the orbital-off-diagonal block of H;")
    print("  xi_x and Bx are orbital-diagonal, xi_y is orbital-off-diagonal but")
    print("  spin-scalar, so all three drop out for the Delta m_s != 0 pairs.")

    print()
    print("=" * 78)
    print("N1-3  inside the island the asymmetry is field- and T-dependent")
    print("=" * 78)
    Ts = [60.0, 70.0, 80.0, 90.0, 101.0, 110.0]
    Bs = [0.0, 0.15, 0.232, 0.35, 0.50]
    print(f"  {'B_perp(T)':>10} " + " ".join(f"{T:>9.0f}K" for T in Ts))
    rows = []
    for B in Bs:
        vals = [asym_kernel(T, B) for T in Ts]
        rows.append((B, vals))
        print(f"  {B:10.3f} " + " ".join(f"{v:10.4f}" for v in vals))
    with (OUT / "N1_asymmetry_vs_T_B.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["B_perp_T"] + [f"A_at_{T:g}K" for T in Ts])
        for B, vals in rows:
            w.writerow([B] + [f"{v:.6f}" for v in vals])

    print()
    print("=" * 78)
    print("N1-4  ensemble strain spread at the operating field B_perp = 0.232 T")
    print("=" * 78)
    rng = np.random.default_rng(20260820)
    print(f"  {'T(K)':>6} {'sigma_strain(GHz)':>18} {'<A>':>10} {'std':>10} {'rel':>8}")
    ens = []
    for T in (70.0, 80.0, 90.0, 101.0):
        for sig in (1.683, 5.0):
            xs, ys = rng.normal(0, sig, 200), rng.normal(0, sig, 200)
            As = np.array([asym_kernel(T, 0.232, xi_x=x, xi_y=y)
                           for x, y in zip(xs, ys)])
            ens.append((T, sig, As.mean(), As.std(), As.std() / As.mean()))
            print(f"  {T:6.1f} {sig:18.3f} {As.mean():10.6f} {As.std():10.2e}"
                  f" {As.std()/As.mean():8.1e}")
    with (OUT / "N1_ensemble_spread.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["T_K", "sigma_strain_GHz", "A_mean", "A_std", "A_rel_std"])
        for r in ens:
            w.writerow([f"{v:.6g}" for v in r])

    print()
    print("  VERDICT: the null is exact only at zero transverse field, which is")
    print("  outside the transparency island (B_perp >~ 0.15 T is required).")
    print("  Inside the island A is still fully predicted by (T, B_perp) with no")
    print("  free parameters, but it acquires a sample-dependent strain spread")
    print("  that falls from ~11% at 70 K to ~2% at 101 K for sigma = 1.683 GHz.")
    print("  Gate N1 therefore does NOT pass as a null test; it passes only in")
    print("  the weaker form of a predicted curve requiring strain characterization.")


if __name__ == "__main__":
    main()
