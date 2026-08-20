"""
Generalised probe/control channel selection for the full nine-level NV model.

`gate2_candidate_full_vs_reduced.build_full` hard-codes the probe on dressed
ground index 1 (ms = 0) and lets only the control move, so it can build the
class-2 Lambda (0,+1) and (0,-1) but not the class-3 Lambda (-1,+1). The body
below is that function with the single change that both indices are arguments.

Everything else -- Hamiltonian, jump operators, frame, ISC and hyperfine
options -- is reproduced verbatim. `test_matches_canonical()` asserts that this
version with (p_idx, c_idx) = (1, 2) is byte-identical to the canonical one
with ctrl='+1', so the two cannot drift apart silently.
"""
from pathlib import Path
import sys
import numpy as np

_NOGO = Path(__file__).resolve().parents[2] / "No-go theorem" / "src"
sys.path.insert(0, str(_NOGO))

import matplotlib
matplotlib.use("Agg")
import run_prl_prediction as rp
import nv_model as nv
import gate2_candidate_full_vs_reduced as g2
from liouvillian_core import liouvillian, steady_state, first_order

TWOPI = 2 * np.pi
MS_OF = g2.MS_OF                     # [-1, 0, 1] -> dressed ground index order
IDX_OF_MS = {ms: i for i, ms in enumerate(MS_OF)}


def build_full_pair(T, Bx, Bz, d2, p_idx, c_idx, Oc=rp.OC, ppol=rp.Y,
                    cpol=rp.Y, isc=False, mI=None, d=rp.D, phi=rp.PHI,
                    j0=rp.J0):
    """As gate2_candidate_full_vs_reduced.build_full, with p_idx exposed."""
    N = 10 if isc else 9
    Bvec = (float(Bx), 0.0, float(Bz))
    Sz = nv.Sz
    Hg0 = nv.Hgs(Bvec)
    He0 = nv.Hes(Bvec, d, phi)
    eg0, U0 = g2.dressed_from(Hg0)
    z0 = float(np.linalg.eigvalsh(He0)[j0])
    if mI is None:
        Hg, He, eg, U = Hg0, He0, eg0, U0
    else:
        Hg = Hg0 + mI * g2.A_GS * Sz
        He = He0 + mI * g2.A_ES * np.kron(nv.I2, Sz)
        eg, U = g2.dressed_from(Hg)
    s_idx = 3 - p_idx - c_idx
    dp = np.kron(np.asarray(ppol, complex), U[:, p_idx])
    dc = np.kron(np.asarray(cpol, complex), U[:, c_idx])
    H = np.zeros((N, N), complex)
    H[3:9, 3:9] = He - (z0 + d2) * np.eye(6) + (eg[p_idx] - eg0[p_idx]) * np.eye(6)
    H[c_idx, c_idx] = -d2 + (eg[c_idx] - eg0[c_idx]) - (eg[p_idx] - eg0[p_idx])
    H[s_idx, s_idx] = float(eg[s_idx] - eg[p_idx])
    Vc = np.zeros((N, N), complex); Vc[3:9, c_idx] = dc
    Vp = np.zeros((N, N), complex); Vp[3:9, p_idx] = dp
    H += 0.5 * Oc * (Vc + Vc.conj().T)
    rate = nv.korb_GHz(T, d)
    Ls = []
    for m in range(3):
        up = np.zeros((N, N), complex); dn = np.zeros((N, N), complex)
        up[6 + m, 3 + m] = 1; dn[3 + m, 6 + m] = 1
        Ls += [np.sqrt(rate) * up, np.sqrt(rate) * dn]
    for orb0 in (3, 6):
        for m in range(3):
            J = np.zeros((N, N), complex)
            for a in range(3):
                J[a, orb0 + m] = np.conj(U[m, a])
            Ls.append(np.sqrt(nv.GRAD) * J)
    for a in range(3):
        for b in range(3):
            if a != b:
                J = np.zeros((N, N), complex); J[b, a] = 1
                Ls.append(np.sqrt(g2.T1_GROUND) * J)
    for a in range(3):
        J = np.zeros((N, N), complex); J[a, a] = 1
        Ls.append(np.sqrt(2 * g2.GG) * J)
    if isc:
        for orb0 in (3, 6):
            for m in range(3):
                J = np.zeros((N, N), complex); J[9, orb0 + m] = 1
                Ls.append(np.sqrt(g2.K_ISC[MS_OF[m]]) * J)
        for m in range(3):
            J = np.zeros((N, N), complex)
            for a in range(3):
                J[a, 9] = np.conj(U[m, a])
            Ls.append(np.sqrt(g2.K_SINGLET * g2.P_SINGLET[MS_OF[m]]) * J)
    return TWOPI * H, Ls, Vp, dp, dict(N=N, p_idx=p_idx, c_idx=c_idx, z0=z0)


def chi_pair_idx(T, Bx, Bz, d2, p_idx, c_idx, **kw):
    """(chi_full, chi_cut, diagnostics); cut removes the ground-coherence sector."""
    H, Ls, Vp, dp, meta = build_full_pair(T, Bx, Bz, d2, p_idx, c_idx, **kw)
    N, p_idx, c_idx = meta["N"], meta["p_idx"], meta["c_idx"]
    L = liouvillian(H, Ls)
    rho0, res0, gap = steady_state(L)
    Hp = TWOPI * 0.5 * g2.OP * (Vp + Vp.conj().T)
    I = np.eye(N)
    V = -1j * (np.kron(I, Hp) - np.kron(Hp.T, I))
    det = np.zeros(N * N, complex)
    for e, a in enumerate(dp):
        det[p_idx * N + (3 + e)] = np.conj(a)
    S = [c_idx * N + p_idx, p_idx * N + c_idx]
    X = [k for k in range(N * N) if k not in S]
    Lc = L.copy(); Lc[np.ix_(S, X)] = 0; Lc[np.ix_(X, S)] = 0
    xf, rf = first_order(L, V, rho0)
    xc, rc = first_order(Lc, V, rho0)
    chif = -2 * (det.conj() @ xf) / g2.OP
    chic = -2 * (det.conj() @ xc) / g2.OP
    return complex(chif), complex(chic), dict(res0=res0, res_full=rf,
                                              res_cut=rc, gap=gap)


def test_matches_canonical(verbose=True):
    """(p_idx,c_idx)=(1,2) must reproduce the canonical build_full(ctrl='+1')."""
    args = (70.0, rp.BX0, rp.BZ0, 0.017)
    Ha, Lsa, Vpa, dpa, ma = g2.build_full(*args, ctrl="+1")
    Hb, Lsb, Vpb, dpb, mb = build_full_pair(*args, p_idx=1, c_idx=2)
    checks = {
        "H": np.max(np.abs(Ha - Hb)),
        "Vp": np.max(np.abs(Vpa - Vpb)),
        "dp": np.max(np.abs(dpa - dpb)),
        "n_jumps": abs(len(Lsa) - len(Lsb)),
        "Ls": max((np.max(np.abs(x - y)) for x, y in zip(Lsa, Lsb)), default=0.0),
    }
    fa = g2.chi_pair(*args, ctrl="+1")[0]
    fb = chi_pair_idx(*args, p_idx=1, c_idx=2)[0]
    checks["chi"] = abs(fa - fb)
    ok = all(v == 0 or v < 1e-13 for v in checks.values()) and ma == mb
    if verbose:
        for k, v in checks.items():
            print(f"    {k:9s} max|diff| = {v:.3e}")
        print(f"  regression vs canonical build_full: {'PASS' if ok else 'FAIL'}")
    assert ok, checks
    return ok


if __name__ == "__main__":
    print("  regression test of full_pair.build_full_pair")
    test_matches_canonical()
