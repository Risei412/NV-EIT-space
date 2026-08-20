"""
Gate N4b - which orbital basis are the Lambda's two ground states in?

N4 derived M0 = A^2 (e_c x e_p) for ground states |g,x>, |g,y> and concluded
theta is just the relative polarization of probe and control. That assumed the
two ground states are the LINEAR orbital states. They need not be: the ground
spin-orbit interaction is lambda*L_z*S_z, L_z is the A_2 component of E (x) E,
and its eigenstates are the CIRCULAR combinations.

With d(e) = e_x D_X + e_y D_Y, D_X = sigma_z, D_Y = sigma_x (the E components
in the linear basis),

    d(e_c)^dag d(e_p) = (e_c . e_p) I + C (i sigma_y),
    C = conj(e_cx) e_py - conj(e_cy) e_px

so for orthogonal ground states the identity piece drops and

    M0 = C * <g_c| i sigma_y |g_p> .

The second factor is the ground manifold's A_2 (angular-momentum) matrix
element. It vanishes when the ground states are L_z eigenstates.

CONVENTIONS. sigma_y is the Hermitian A_2 generator; the physical spin-orbit
term is (lambda/2) sigma_y at fixed spin. i*sigma_y is anti-Hermitian and
appears only inside the operator product above, never as a Hamiltonian.
"""
import csv
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
_NOGO = Path(__file__).resolve().parents[2] / "No-go theorem" / "src"
sys.path.insert(0, str(_NOGO))
import group_iv_model as giv

OUT = Path(__file__).resolve().parents[1] / "out"
sz = np.array([[1, 0], [0, -1]], complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)     # Hermitian A_2 generator = L_z
D_X, D_Y = sz, sx


def M0(e_p, e_c, g_p, g_c):
    c = (e_p[0] * D_X + e_p[1] * D_Y) @ g_p
    p = (e_c[0] * D_X + e_c[1] * D_Y) @ g_c
    return complex(np.vdot(p, c))


def main():
    OUT.mkdir(exist_ok=True)
    X, Y = np.array([1, 0], complex), np.array([0, 1], complex)
    PLUS = np.array([1, 1j], complex) / np.sqrt(2)
    MINUS = np.array([1, -1j], complex) / np.sqrt(2)

    print("=" * 78)
    print("N4b-1  the factorisation, verified with the correct sign")
    print("=" * 78)
    rng = np.random.default_rng(11)
    worst = 0.0
    for _ in range(2000):
        ep = rng.normal(size=2) + 1j * rng.normal(size=2); ep /= np.linalg.norm(ep)
        ec = rng.normal(size=2) + 1j * rng.normal(size=2); ec /= np.linalg.norm(ec)
        C = np.conj(ec[0]) * ep[1] - np.conj(ec[1]) * ep[0]
        for gp, gc in ((X, Y), (PLUS, MINUS)):
            rhs = C * np.vdot(gc, (1j * sy) @ gp)
            worst = max(worst, abs(M0(ep, ec, gp, gc) - rhs))
    print(f"  M0 = C * <g_c| i sigma_y |g_p>   max|diff| = {worst:.2e}")

    print()
    print("=" * 78)
    print("N4b-2  the ground-manifold factor decides the class")
    print("=" * 78)
    print(f"  {'ground basis':>26} {'|<g_c|sigma_y|g_p>|':>21} {'|M0|':>8} {'class':>7}")
    for lbl, gp, gc in (("linear   |x>, |y>", X, Y),
                        ("circular |+>, |->", PLUS, MINUS)):
        ov = abs(np.vdot(gc, sy @ gp))
        m = abs(M0(X, Y, gp, gc))
        print(f"  {lbl:>26} {ov:21.4f} {m:8.4f} {1 if m > 1e-9 else '>=2':>7}")
    print("  Spin-orbit eigenstates are the circular ones, and there M0 vanishes.")

    print()
    print("=" * 78)
    print("N4b-3  class 1 requires strain to dominate the ground spin-orbit")
    print("=" * 78)
    print("  H_g (fixed spin) = (lambda/2) sigma_y + xi sigma_x")
    print(f"  {'xi/lambda':>11} {'|M0| (crossed pol)':>20} {'class':>10}")
    rows = []
    for r in (0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0):
        Hg = 0.5 * sy + r * sx
        w, U = np.linalg.eigh(Hg)
        m = abs(M0(X, Y, U[:, 0], U[:, 1]))
        rows.append(dict(xi_over_lambda=r, absM0=m))
        print(f"  {r:11.3g} {m:20.6f} {('1' if m > 0.5 else '>=2 (suppressed)'):>10}")

    def xi_for(target, lam):
        lo, hi = 1e-6, 1e9
        for _ in range(300):
            mid = np.sqrt(lo * hi)
            w, U = np.linalg.eigh(0.5 * lam * sy + mid * sx)
            if abs(M0(X, Y, U[:, 0], U[:, 1])) < target:
                lo = mid
            else:
                hi = mid
        return np.sqrt(lo * hi)

    print()
    print("  Ground spin-orbit (SiV_SnV_phonon_AIC_parameters.md):"
          " SiV ~48 GHz, SnV ~850 GHz")
    print(f"  {'material':>8} {'lambda_g(GHz)':>14} {'xi for |M0|=0.5':>17}"
          f" {'xi for |M0|=0.9':>17}")
    need_rows = []
    for mat, lam_g in (("SiV", 48.0), ("SnV", 850.0)):
        a, b = xi_for(0.5, lam_g), xi_for(0.9, lam_g)
        need_rows.append(dict(material=mat, lambda_g_GHz=lam_g,
                              xi_M0_half_GHz=a, xi_M0_p9_GHz=b))
        print(f"  {mat:>8} {lam_g:14.1f} {a:17.1f} {b:17.1f}")

    print()
    print("=" * 78)
    print("N4b-4  the repository's group-IV model cannot answer this")
    print("=" * 78)
    print("  group_iv_model.H_groupIV writes")
    print("    H_SO     = (Delta_e/2) kron(sigma_z_orb, sigma_z_spin)")
    print("    H_strain = xi_x kron(sigma_z_orb, I) + xi_y kron(sigma_x_orb, I)")
    print("  Spin-orbit is A_2 and should be sigma_y; strain is E and is correctly")
    print("  (sigma_z, sigma_x). Writing spin-orbit as sigma_z puts it in the same")
    print("  representation slot as one strain component, so no orbital basis makes")
    print("  both terms simultaneously correct. group_iv_model.legs then sends both")
    print("  legs to a common excited ket -- the linear-basis answer -- while H_SO is")
    print("  written as though the basis were already diagonal in spin-orbit.")
    print(f"  model |M0| at theta = 0: {abs(giv.moments(giv.H_groupIV('SiV'), 0)[0]):.3f}")
    print("  This is the geometry NON_CLAIMS N4 flags as schematic. The class-1")
    print("  assignment is a consequence of that schematic choice, not of D3d.")

    with (OUT / "N4b_class_vs_strain.csv").open("w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w_.writeheader(); w_.writerows(rows)
    with (OUT / "N4b_strain_required.csv").open("w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=list(need_rows[0].keys()))
        w_.writeheader(); w_.writerows(need_rows)


if __name__ == "__main__":
    main()
