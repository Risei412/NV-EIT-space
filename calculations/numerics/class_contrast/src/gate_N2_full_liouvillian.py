"""
Gate N2, full-Liouvillian recomputation.

The reduced-kernel pass (gate_N2_channel_ratio.py) predicted the
calibration-free relation R(T)*Gamma(T) = M2/M1 = 58.4176 GHz, holding to
better than 1% for T >= 70 K, and reported PASS. NON_CLAIMS N7.4 requires that
to be repeated against the full nine-level Liouvillian before it is quoted.
This is that repeat. It overturns the result.

Channels, on dressed ground indices (MS_OF = [-1, 0, +1]):
    class 2 : probe ms=0,  control ms=+1
    class 2 : probe ms=0,  control ms=-1
    class 3 : probe ms=-1, control ms=+1
Observable: dA = Im chi_cut - Im chi_full at the peak of the two-photon
feature, i.e. the unnormalised sector-resolved correction, which is what
inherits the kernel's class exponent (the normalised contrast C carries an
extra normalisation and a different exponent).
"""
import csv
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from full_pair import chi_pair_idx, IDX_OF_MS, test_matches_canonical
import run_prl_prediction as rp
import nv_model as nv
import phonon_rates as pr
import gate2_candidate_full_vs_reduced as g2

OUT = Path(__file__).resolve().parents[1] / "out"
D_STRAIN = 1.683
D_GS = 2.877          # GHz, ground zero-field splitting
GE = 28.0             # GHz/T
D2_SCAN = np.linspace(-0.02, 0.02, 41)      # GHz, two-photon detuning
REDUCED_PREDICTION = 58.4176                # GHz, M2/M1

CHANNELS = [("class2_0_p1", "class 2  probe 0, ctrl +1", IDX_OF_MS[0], IDX_OF_MS[+1]),
            ("class2_0_m1", "class 2  probe 0, ctrl -1", IDX_OF_MS[0], IDX_OF_MS[-1]),
            ("class3_m1_p1", "class 3  probe -1, ctrl +1", IDX_OF_MS[-1], IDX_OF_MS[+1])]


def gamma_oc_GHz(T):
    return pr.gamma_oc(T, D_STRAIN) / 1e9


def peak_dA(T, p_idx, c_idx, B):
    """Signed dA at the extremum of |dA| over the two-photon scan."""
    best_dA, best_Ac = 0.0, np.nan
    for d2 in D2_SCAN:
        f, c, _ = chi_pair_idx(T, B, rp.BZ0, float(d2), p_idx, c_idx)
        dA = c.imag - f.imag
        if abs(dA) > abs(best_dA):
            best_dA, best_Ac = dA, c.imag
    return best_dA, best_Ac


def main():
    OUT.mkdir(exist_ok=True)
    print("=" * 80)
    print("N2F-0  regression: generalised builder vs canonical build_full")
    print("=" * 80)
    test_matches_canonical()

    print()
    print("=" * 80)
    print("N2F-1  why the class label may not survive: ground-manifold mixing")
    print("=" * 80)
    print(f"  bare spin sectors are good labels only while ge*B_perp << D_gs"
          f" = {D_GS} GHz,")
    print(f"  i.e. B_perp << {D_GS/GE:.3f} T. The transparency island needs"
          f" B_perp >~ 0.15 T.")
    print(f"  {'B_perp(T)':>10} {'ge*B(GHz)':>10}   |<ms|dressed>|^2  per dressed state")
    mix_rows = []
    for B in (0.0, 0.05, 0.10, 0.15, 0.232, 0.50):
        eg, U = g2.dressed_from(nv.Hgs((B, 0.0, rp.BZ0)))
        w = np.abs(U) ** 2
        s = " | ".join(" ".join(f"{w[m, a]:.2f}" for m in range(3)) for a in range(3))
        print(f"  {B:10.3f} {GE*B:10.2f}   {s}")
        mix_rows.append([B, GE * B] + [f"{w[m, a]:.4f}" for a in range(3) for m in range(3)])
    with (OUT / "N2full_ground_mixing.csv").open("w", newline="") as f:
        w_ = csv.writer(f)
        w_.writerow(["B_perp_T", "ge_B_GHz"] +
                    [f"P_ms{ms}_in_dressed{a}" for a in range(3) for ms in (-1, 0, 1)])
        w_.writerows(mix_rows)

    print()
    print("=" * 80)
    print("N2F-2  full-Liouvillian channel sweep at the operating field"
          f" B_perp = {rp.BX0:.3f} T")
    print("=" * 80)
    Ts = [50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0]
    print(f"  {'T(K)':>6} {'Gamma(GHz)':>11} " +
          " ".join(f"{lbl.split()[1]+lbl.split()[-1]:>13}" for _, lbl, _, _ in CHANNELS))
    rows = []
    for T in Ts:
        G = gamma_oc_GHz(T)
        rec = dict(T_K=T, Gamma_GHz=G)
        for key, _, p, c in CHANNELS:
            dA, Ac = peak_dA(T, p, c, rp.BX0)
            rec[f"dA_{key}"] = dA
            rec[f"C_{key}"] = dA / Ac
        rows.append(rec)
        print(f"  {T:6.1f} {G:11.2f} " +
              " ".join(f"{rec['dA_'+k]:13.4e}" for k, _, _, _ in CHANNELS))

    print()
    print("=" * 80)
    print("N2F-3  the calibration-free prediction, tested")
    print("=" * 80)
    print(f"  reduced-kernel prediction:  |R|*Gamma = {REDUCED_PREDICTION:.4f} GHz, constant")
    print(f"  {'T(K)':>6} {'Gamma':>9} {'|R|':>11} {'|R|*Gamma':>12} {'sign c2':>8} {'sign c3':>8}")
    for r in rows:
        a, b = r["dA_class2_0_p1"], r["dA_class3_m1_p1"]
        R = abs(b / a)
        r["R_abs"] = R
        r["R_times_Gamma"] = R * r["Gamma_GHz"]
        print(f"  {r['T_K']:6.1f} {r['Gamma_GHz']:9.2f} {R:11.3e} {R*r['Gamma_GHz']:12.3e}"
              f" {np.sign(a):8.0f} {np.sign(b):8.0f}")
    vals = [r["R_times_Gamma"] for r in rows]
    print(f"  spread of |R|*Gamma over 50-120 K: {min(vals):.2e} .. {max(vals):.2e}"
          f"  ({max(vals)/min(vals):.1e} x)")

    with (OUT / "N2full_channel_sweep.csv").open("w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w_.writeheader()
        w_.writerows(rows)

    print()
    print("=" * 80)
    print("N2F-4  is any transverse field better? (weaker field = less mixing)")
    print("=" * 80)
    field_rows = []
    for B in (0.01, 0.05, rp.BX0):
        print(f"  --- B_perp = {B:.3f} T (ge*B = {GE*B:.2f} GHz) ---")
        print(f"  {'T(K)':>6} {'Gamma':>9} {'dA class2':>12} {'dA class3':>12}"
              f" {'|R|':>11} {'|R|*Gamma':>12}")
        for T in (70.0, 80.0, 90.0, 100.0):
            G = gamma_oc_GHz(T)
            a, _ = peak_dA(T, IDX_OF_MS[0], IDX_OF_MS[+1], B)
            b, _ = peak_dA(T, IDX_OF_MS[-1], IDX_OF_MS[+1], B)
            R = abs(b / a)
            field_rows.append(dict(B_perp_T=B, T_K=T, Gamma_GHz=G, dA_class2=a,
                                   dA_class3=b, R_abs=R, R_times_Gamma=R * G))
            print(f"  {T:6.1f} {G:9.2f} {a:12.3e} {b:12.3e} {R:11.3e} {R*G:12.3e}")
    with (OUT / "N2full_field_sweep.csv").open("w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=list(field_rows[0].keys()))
        w_.writeheader()
        w_.writerows(field_rows)

    print()
    print("=" * 80)
    print("N2F-5  helicity asymmetry: full vs reduced (Gate N1 cross-check)")
    print("=" * 80)
    reduced_A = {70.0: 0.2238, 80.0: 0.3042, 90.0: 0.3640, 101.0: 0.4082}
    print(f"  {'T(K)':>6} {'|dA(0,+1)|':>12} {'|dA(0,-1)|':>12} {'A_full':>9} {'A_reduced':>10}")
    asym_rows = []
    for T in (70.0, 80.0, 90.0, 101.0):
        a = abs(peak_dA(T, IDX_OF_MS[0], IDX_OF_MS[+1], rp.BX0)[0])
        b = abs(peak_dA(T, IDX_OF_MS[0], IDX_OF_MS[-1], rp.BX0)[0])
        A = (a - b) / (a + b)
        asym_rows.append(dict(T_K=T, dA_plus=a, dA_minus=b, A_full=A,
                              A_reduced=reduced_A[T]))
        print(f"  {T:6.1f} {a:12.4e} {b:12.4e} {A:9.4f} {reduced_A[T]:10.4f}")
    with (OUT / "N2full_helicity_asymmetry.csv").open("w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=list(asym_rows[0].keys()))
        w_.writeheader()
        w_.writerows(asym_rows)

    print()
    print("  VERDICT: the reduced-kernel PASS does not survive. |R|*Gamma is not")
    print("  constant at any field tried; it moves by orders of magnitude, and the")
    print("  class-3 channel carries the opposite sign to class 2 throughout.")
    print("  Diagnosis in N2F-1: the transverse field that opens the Raman channel")
    print("  (B_perp >~ 0.15 T) is already past the field that dissolves the bare")
    print("  spin sectors (ge*B_perp = D_gs at B_perp = 0.103 T), and the class")
    print("  label is defined on those sectors. The two windows do not overlap.")


if __name__ == "__main__":
    main()
