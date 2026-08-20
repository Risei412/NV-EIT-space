"""
Gate N4 - what is theta?

group_iv_model.legs documents theta as "the (unmeasured) ratio of cross-branch
to same-branch dipole matrix elements", sets it to 0 as the representative
case, and Gate N3 showed the whole class-1 argument collapses if the real value
is near pi/2. This gate derives it.

Method: Wigner-Eckart in D3d. Both orbital doublets and the in-plane dipole
transform as E, and E (x) E = A1 + A2 + E contains E exactly once, so a single
reduced matrix element fixes all four orbital transitions. In the real orbital
basis {x, y} the E component of E (x) E is the symmetric traceless pair

    {(xx - yy)/sqrt2, (xy + yx)/sqrt2}

so the two dipole components are, up to one common reduced element A,

    d_x = A * sigma_z(orbital),      d_y = A * sigma_x(orbital)

The same E structure appears in group_iv_model.H_groupIV as
xi_x*kron(sz_orb, I) + xi_y*kron(sx_orb, I), which is checked below.
"""
import csv
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
_NOGO = Path(__file__).resolve().parents[2] / "No-go theorem" / "src"
sys.path.insert(0, str(_NOGO))
import group_iv_model as giv
from nv3e_loader import load

nv3e = load()
OUT = Path(__file__).resolve().parents[1] / "out"

sz = np.array([[1, 0], [0, -1]], complex)      # orbital, real {x,y} basis
sx = np.array([[0, 1], [1, 0]], complex)
I2 = np.eye(2, dtype=complex)

# E-type dipole components, one common reduced element A (set to 1)
D_X, D_Y = sz, sx


def excited_dest(ground_orb, pol):
    """Excited orbital state reached from `ground_orb` by polarization `pol`."""
    ex, ey = pol
    return (ex * D_X + ey * D_Y) @ ground_orb


def M0_of(pol_probe, pol_ctrl, g_probe, g_ctrl):
    """M0 = <control leg | probe leg> in the excited orbital space."""
    c = excited_dest(g_probe, pol_probe)         # source
    p = excited_dest(g_ctrl, pol_ctrl)           # readout
    return complex(np.vdot(p, c))


def main():
    OUT.mkdir(exist_ok=True)
    X = np.array([1, 0], complex)
    Y = np.array([0, 1], complex)
    sig_p = np.array([1, 1j], complex) / np.sqrt(2)
    sig_m = np.array([1, -1j], complex) / np.sqrt(2)

    print("=" * 78)
    print("N4-1  the E (x) E dipole structure, checked against the model's strain term")
    print("=" * 78)
    Hs = giv.H_groupIV("SiV", xi_x=1.0, xi_y=0.0) - giv.H_groupIV("SiV")
    Hc = giv.H_groupIV("SiV", xi_x=0.0, xi_y=1.0) - giv.H_groupIV("SiV")
    print(f"  d(H)/d(xi_x) == kron(sigma_z_orb, I) : "
          f"{np.max(np.abs(Hs - np.kron(sz, I2))):.1e}")
    print(f"  d(H)/d(xi_y) == kron(sigma_x_orb, I) : "
          f"{np.max(np.abs(Hc - np.kron(sx, I2))):.1e}")
    print("  strain and dipole carry the same E structure, as they must.")

    print()
    print("=" * 78)
    print("N4-2  M0 for the orbital Lambda: probe on |g,x>, control on |g,y>")
    print("=" * 78)
    print(f"  {'probe pol':>12} {'control pol':>12} {'M0':>22} {'|M0|':>7} {'theta(deg)':>11}")
    rows = []
    for lp, ep in (("x", X), ("y", Y), ("sigma+", sig_p), ("sigma-", sig_m)):
        for lc, ec in (("x", X), ("y", Y), ("sigma+", sig_p), ("sigma-", sig_m)):
            m0 = M0_of(ep, ec, X, Y)
            th = np.degrees(np.arccos(min(1.0, abs(m0))))
            rows.append(dict(probe_pol=lp, ctrl_pol=lc, M0_re=m0.real,
                             M0_im=m0.imag, absM0=abs(m0), theta_deg=th))
            print(f"  {lp:>12} {lc:>12} {str(np.round(m0, 6)):>22}"
                  f" {abs(m0):7.4f} {th:11.2f}")

    print()
    print("=" * 78)
    print("N4-3  the closed form: M0 is the polarization cross product")
    print("=" * 78)
    print("  probe on |g,x>:   c = A ( e_px |e,x> + e_py |e,y> )")
    print("  control on |g,y>: p = A ( e_cy |e,x> - e_cx |e,y> )")
    print("  so  M0 = p^dag c = A^2 ( e_cy* e_px - e_cx* e_py )  =  A^2 (e_c x e_p)_z")
    worst = 0.0
    rng = np.random.default_rng(7)
    for _ in range(2000):
        ep = rng.normal(size=2) + 1j * rng.normal(size=2); ep /= np.linalg.norm(ep)
        ec = rng.normal(size=2) + 1j * rng.normal(size=2); ec /= np.linalg.norm(ec)
        lhs = M0_of(ep, ec, X, Y)
        rhs = np.conj(ec[1]) * ep[0] - np.conj(ec[0]) * ep[1]
        worst = max(worst, abs(lhs - rhs))
    print(f"  verified on 2000 random polarization pairs, max|diff| = {worst:.2e}")
    print("  |M0| = 1 when the two polarizations are orthogonal, 0 when parallel.")
    print("  theta is therefore NOT a material constant. It is the relative")
    print("  polarization of probe and control, and theta = 0 is the ordinary")
    print("  crossed-polarization configuration.")

    print()
    print("=" * 78)
    print("N4-4  the same construction for NV, for contrast")
    print("=" * 78)
    print("  NV's optical dipole is orbital and spin-scalar, and the two ground")
    print("  states of its Lambda differ in SPIN, not in orbital branch. The legs")
    print("  therefore carry a spin overlap factor <s_ctrl|s_probe> that no choice")
    print("  of polarization can change:")
    for pair in ((0, -1), (0, +1), (-1, +1)):
        m0 = abs(nv3e.moments(nv3e.H_3E(), pair, 0)[0])
        print(f"    NV pair {str(pair):>9}:  |M0| = {m0:.1e}  (spin overlap = 0)")
    print("  For group-IV the two ground states differ in ORBITAL branch and share")
    print("  the spin, so the overlap factor is <s|s> = 1 and only the orbital")
    print("  cross product survives.")

    print()
    print("=" * 78)
    print("N4-5  robustness of theta = 0 against a misset polarization angle")
    print("=" * 78)
    print(f"  {'misalignment (deg)':>20} {'|M0|':>9} {'theta(deg)':>11} {'class':>7}")
    for da in (0, 1, 5, 10, 30, 60, 89, 90):
        a = np.deg2rad(da)
        ec = np.array([np.sin(a), np.cos(a)], complex)   # 90 deg from probe at da=0
        m0 = M0_of(X, ec, X, Y)
        th = np.degrees(np.arccos(min(1.0, abs(m0))))
        print(f"  {da:20d} {abs(m0):9.4f} {th:11.2f} {1 if abs(m0) > 1e-9 else '>=2':>7}")
    print("  A 10 degree polarization error costs 1.5% of M0; class 1 survives any")
    print("  misalignment short of exactly parallel.")

    with (OUT / "N4_theta_polarization.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
