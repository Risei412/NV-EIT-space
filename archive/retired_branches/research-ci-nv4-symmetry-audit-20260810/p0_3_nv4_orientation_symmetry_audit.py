"""P0-3: four-orientation NV symmetry audit.

Purpose
-------
Test what the four <111> NV orientations actually forbid/generate once the
low-temperature 3E orbital frame, transverse strain, optical polarization and
spin-channel choice are retained.  The observable is deliberately branch-order
independent: sum the raw signed sector response Delta A = A_cut-A_full over all
six excited branches at two-photon resonance, then average over an orientation
frame orbit.

Two frame sets are compared:
  V4  : four marked frames obtained by the pi rotations diag(+--), diag(-+-),
        diag(--+).  These reproduce the four NV axes but keep one transverse
        frame mark per axis.
  T12 : V4 completed by the three C3-related transverse frames per axis.  This
        is the 12-element proper-tetrahedral group orbit.

The full calculation uses the existing validated nine-level Lindblad structure
from gate2_candidate_full_vs_reduced.py, generalized only to arbitrary local
(Bx,By,Bz).  No ISC or hyperfine terms are added in this audit.

Outputs
-------
Writing Paper/pra/results/tables/p0_3_nv4_orientation_symmetry_audit.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PRA = HERE.parent
REPO = PRA.parents[1]
NOGO_SRC = REPO / "No-go theorem" / "src"
sys.path.insert(0, str(NOGO_SRC))

import nv_model as nv  # noqa: E402
from liouvillian_core import liouvillian, steady_state, first_order  # noqa: E402

TWOPI = 2*np.pi
GG = 6.3e-5
T1_GROUND = 1e-6
OP = 1e-5
T_K = 70.0
OC = 0.1
D_STRAIN = 1.683
PHI = 0.0
POL_Y = np.array([1.0, 0.0], complex)
H_STEP = 1e-4  # tesla


def build_full_vec(Bvec, d2=0.0, Oc=OC, ctrl="+1", ppol=POL_Y,
                   cpol=POL_Y, d=D_STRAIN, phi=PHI, j0=3):
    """Nine-level full model with arbitrary local B=(Bx,By,Bz)."""
    N = 9
    Bvec = tuple(float(x) for x in Bvec)
    Hg = nv.Hgs(Bvec)
    He = nv.Hes(Bvec, d, phi)
    eg, U = nv.dressed_ground(Bvec)
    z0 = float(np.linalg.eigvalsh(He)[j0])

    p_idx = 1
    c_idx = 2 if ctrl == "+1" else 0
    s_idx = 3 - p_idx - c_idx
    dp = np.kron(np.asarray(ppol, complex), U[:, p_idx])
    dc = np.kron(np.asarray(cpol, complex), U[:, c_idx])

    H = np.zeros((N, N), complex)
    H[3:9, 3:9] = He - (z0+d2)*np.eye(6)
    H[c_idx, c_idx] = -d2
    H[s_idx, s_idx] = float(eg[s_idx] - eg[p_idx])

    Vc = np.zeros((N, N), complex)
    Vp = np.zeros((N, N), complex)
    Vc[3:9, c_idx] = dc
    Vp[3:9, p_idx] = dp
    H += 0.5*Oc*(Vc + Vc.conj().T)

    rate = nv.korb_GHz(T_K, d)
    Ls = []
    for m in range(3):
        up = np.zeros((N, N), complex)
        dn = np.zeros((N, N), complex)
        up[6+m, 3+m] = 1
        dn[3+m, 6+m] = 1
        Ls += [np.sqrt(rate)*up, np.sqrt(rate)*dn]

    for orb0 in (3, 6):
        for m in range(3):
            J = np.zeros((N, N), complex)
            for a in range(3):
                J[a, orb0+m] = np.conj(U[m, a])
            Ls.append(np.sqrt(nv.GRAD)*J)

    for a in range(3):
        for b in range(3):
            if a != b:
                J = np.zeros((N, N), complex)
                J[b, a] = 1
                Ls.append(np.sqrt(T1_GROUND)*J)
    for a in range(3):
        J = np.zeros((N, N), complex)
        J[a, a] = 1
        Ls.append(np.sqrt(2*GG)*J)

    return TWOPI*H, Ls, Vp, dp, dict(N=N, p_idx=p_idx, c_idx=c_idx)


def chi_pair_vec(Bvec, j0, ctrl="+1"):
    H, Ls, Vp, dp, meta = build_full_vec(Bvec, j0=j0, ctrl=ctrl)
    N, p_idx, c_idx = meta["N"], meta["p_idx"], meta["c_idx"]
    L = liouvillian(H, Ls)
    rho0, _, _ = steady_state(L)

    Hp = TWOPI*0.5*OP*(Vp + Vp.conj().T)
    I = np.eye(N)
    V = -1j*(np.kron(I, Hp) - np.kron(Hp.T, I))

    det = np.zeros(N*N, complex)
    for e, a in enumerate(dp):
        det[p_idx*N + (3+e)] = np.conj(a)

    S = [c_idx*N+p_idx, p_idx*N+c_idx]
    X = [k for k in range(N*N) if k not in S]
    Lc = L.copy()
    Lc[np.ix_(S, X)] = 0
    Lc[np.ix_(X, S)] = 0

    xf, _ = first_order(L, V, rho0)
    xc, _ = first_order(Lc, V, rho0)
    chi_full = -2*(det.conj()@xf)/OP
    chi_cut = -2*(det.conj()@xc)/OP
    return complex(chi_full), complex(chi_cut)


def local_branch_sum(Bvec, ctrl="+1"):
    """Raw signed response summed over all excited branches, order-independent."""
    total = 0.0
    for j0 in range(6):
        chi_full, chi_cut = chi_pair_vec(Bvec, j0=j0, ctrl=ctrl)
        total += float(np.imag(chi_cut) - np.imag(chi_full))
    return total


def frame_sets():
    n0 = np.array([1., 1., 1.])/np.sqrt(3)
    x0 = np.array([1., -1., 0.])/np.sqrt(2)
    y0 = np.cross(n0, x0)
    y0 /= np.linalg.norm(y0)
    F0 = np.vstack([x0, y0, n0])

    Qs = [np.diag([1.,1.,1.]), np.diag([1.,-1.,-1.]),
          np.diag([-1.,1.,-1.]), np.diag([-1.,-1.,1.])]
    v4 = [F0@Q for Q in Qs]

    t12 = []
    for F in v4:
        for theta in (0., 2*np.pi/3, 4*np.pi/3):
            c, s = np.cos(theta), np.sin(theta)
            Rloc = np.array([[c, s, 0.], [-s, c, 0.], [0.,0.,1.]])
            t12.append(Rloc@F)
    return v4, t12


def orbit_response(Bglob, frames, ctrl="+1"):
    Bglob = np.asarray(Bglob, float)
    return float(np.mean([local_branch_sum(F@Bglob, ctrl=ctrl) for F in frames]))


def hessian(func, h=H_STEP):
    z = np.zeros(3)
    f0 = func(z)
    H = np.zeros((3,3))
    for i in range(3):
        e = np.zeros(3); e[i] = h
        H[i,i] = (func(e)+func(-e)-2*f0)/h**2
    for i in range(3):
        for j in range(i+1,3):
            ei = np.zeros(3); ej = np.zeros(3)
            ei[i] = h; ej[j] = h
            H[i,j] = H[j,i] = (
                func(ei+ej)-func(ei-ej)-func(-ei+ej)+func(-ei-ej)
            )/(4*h*h)
    return f0, H


def dxyz(func, h=H_STEP):
    val = 0.0
    for sx in (-1,1):
        for sy in (-1,1):
            for sz in (-1,1):
                val += sx*sy*sz*func(h*np.array([sx,sy,sz], float))
    return float(val/(8*h**3))


def tensor_metrics(H):
    eig = np.linalg.eigvalsh(H)
    iso = np.trace(H)/3
    traceless = H - iso*np.eye(3)
    return dict(
        hessian=H.tolist(),
        eigenvalues=eig.tolist(),
        relative_eigen_spread=float((eig.max()-eig.min())/abs(eig.mean())),
        relative_traceless_norm=float(np.linalg.norm(traceless)/np.linalg.norm(H)),
    )


def main():
    v4, t12 = frame_sets()
    f0_v4, H_v4 = hessian(lambda B: orbit_response(B, v4, "+1"))
    f0_t12, H_t12 = hessian(lambda B: orbit_response(B, t12, "+1"))

    cubic_plus = dxyz(lambda B: orbit_response(B, t12, "+1"))
    cubic_minus = dxyz(lambda B: orbit_response(B, t12, "-1"))
    cubic_pm = 0.5*(cubic_plus+cubic_minus)

    out = dict(
        base_commit="e5273393fe6fcfb7962c8ac37d014a0075b04126",
        model="9-level full Lindblad, base configuration, branch-summed raw sector response",
        T_K=T_K, Omega_c_GHz=OC, transverse_strain_GHz=D_STRAIN,
        finite_difference_step_T=H_STEP,
        V4_four_marked_frames=dict(f0=f0_v4, **tensor_metrics(H_v4)),
        T12_tetrahedral_completion=dict(f0=f0_t12, **tensor_metrics(H_t12)),
        cubic_mixed_derivative_d3R_dBxdBydBz=dict(
            ctrl_plus=cubic_plus,
            ctrl_minus=cubic_minus,
            symmetric_pm=cubic_pm,
            plus_to_symmetric_suppression=float(abs(cubic_plus)/max(abs(cubic_pm),1e-300)),
        ),
        analytic_selection=dict(
            V4_degree_le_3="linear forbidden; diagonal quadratic anisotropy allowed; xyz cubic allowed",
            T12_degree_le_3="linear forbidden; traceless quadratic forbidden; isotropic B^2 allowed; xyz cubic allowed",
        ),
    )

    path = PRA/"results"/"tables"/"p0_3_nv4_orientation_symmetry_audit.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
