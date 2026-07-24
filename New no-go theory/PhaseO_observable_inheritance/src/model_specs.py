"""model_specs.py -- SectorSpec builders for the Gate A observable-inheritance
verification. Every model reuses an existing, already-validated pencil or
Lindbladian in the repository; nothing physical is re-parameterized here.

Bilinear (EIT-type susceptibility) SectorSpecs:
  * three synthetic 3-level pencils spanning the generic / symmetry-protected
    / doubly-protected classes (nu_obs = 2 / 4 / 6) purely by leg selection
    rules -- the clean controlled demonstration of the inheritance law;
  * the physical NV EIT pencil (No-go theorem/src/nv_model.py), the canonical
    nu_obs = 4 witness (M0 = dp^dag dc = 0 by orbital-branch orthogonality).

Quadratic (QFI-type) and frozen-source-difference (Phase P/N) observables are
handled directly in run_gate_a.py against their own modules.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
_NOGO_SRC = os.path.join(_HERE, "..", "..", "..", "No-go theorem", "src")
for _p in (_HERE, _NOGO_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import gate_a_observable as gao  # noqa: E402


# ----------------------------------------------------------------------
# Synthetic bilinear family: generic / protected / doubly protected
# ----------------------------------------------------------------------
# A fixed 3-level Hermitian Hamiltonian (exact rationals). The three cases
# differ ONLY in the control leg dc and (for the doubly-protected case) in
# whether the direct 0-1 coupling is present -- i.e. purely by selection
# rule, not by changing the dissipation scale.
_DIAG = [sp.Integer(0), sp.Rational(13, 10), sp.Rational(17, 10)]


def _H_sym(h01):
    """3x3 Hermitian H with tunable |0>-|1> coupling h01 (rational)."""
    H = sp.zeros(3, 3)
    for i in range(3):
        H[i, i] = _DIAG[i]
    H[0, 1] = sp.I * h01
    H[1, 0] = -sp.I * h01
    H[0, 2] = sp.I * sp.Rational(1, 3)
    H[2, 0] = -sp.I * sp.Rational(1, 3)
    H[1, 2] = sp.I * sp.Rational(1, 5)
    H[2, 1] = -sp.I * sp.Rational(1, 5)
    return H


def _spec_from_symbolic(H_sym, dp_sym, dc_sym, g_eff, beta, label):
    """Build a SectorSpec with matched numeric and symbolic pencils
    A(Gamma,z) = Gamma*I + i*(H - z*I), D = I."""
    H_num = np.array(H_sym.evalf(), dtype=complex)
    n = H_num.shape[0]
    D = np.eye(n, dtype=complex)
    dp = np.array(dp_sym.evalf(), dtype=complex).flatten()
    dc = np.array(dc_sym.evalf(), dtype=complex).flatten()

    def B_of_z(z, _H=H_num, _n=n):
        return 1j * (_H - z * np.eye(_n))

    D_sym = sp.eye(n)

    def B_sym_of_z(z, _H=H_sym, _n=n):
        return sp.I * (_H - z * sp.eye(_n))

    return gao.SectorSpec(
        D=D, B_of_z=B_of_z, dp=dp, dc=dc, g_eff=g_eff, beta=beta,
        readout_mode="bilinear", label=label,
        D_sym=D_sym, B_sym_of_z=B_sym_of_z, dp_sym=dp_sym, dc_sym=dc_sym,
    )


def synthetic_specs(g_eff=1e-3, beta=1.0):
    """Return [generic (nu=2), protected (nu=4), doubly-protected (nu=6)]."""
    e0 = sp.Matrix([1, 0, 0])
    e1 = sp.Matrix([0, 1, 0])

    # generic: dp^dag dc = 3/5 != 0  -> n12 = 1 -> nu_obs = 2
    dc_generic = sp.Rational(3, 5) * e0 + sp.Rational(4, 5) * e1
    generic = _spec_from_symbolic(
        _H_sym(sp.Rational(1, 2)), e0, dc_generic, g_eff, beta,
        "synthetic-generic (nu_obs=2)")

    # protected: dp^dag dc = 0, H[0,1] != 0 -> M1 != 0 -> n12 = 2 -> nu_obs = 4
    protected = _spec_from_symbolic(
        _H_sym(sp.Rational(1, 2)), e0, e1, g_eff, beta,
        "synthetic-protected (nu_obs=4)")

    # doubly protected: dp^dag dc = 0 AND H[0,1] = 0, bridged by |0>-|2>-|1>
    # -> M0 = M1 = 0, M2 != 0 -> n12 = 3 -> nu_obs = 6
    doubly = _spec_from_symbolic(
        _H_sym(sp.Integer(0)), e0, e1, g_eff, beta,
        "synthetic-doubly-protected (nu_obs=6)")

    return [generic, protected, doubly]


# ----------------------------------------------------------------------
# Physical NV EIT bilinear pencil (No-go theorem/src/nv_model.py)
# ----------------------------------------------------------------------
def nv_spec(T_boundary_K=300.0):
    """NV(-) ZPL spin-Lambda EIT reduced pencil, matching RoomT step3's
    numeric_cross_check convention: A = Gamma*I + i*2pi*(H - z*I), Gamma the
    native excited dissipation scale. Probe on orbital X / ms=0, control on
    orbital Y / ms=-1 (dp^dag dc = 0 -> n12 = 2 -> nu_obs = 4)."""
    import nv_model as nm

    Bvec = (0.0, 0.0, 0.02)
    D_STRAIN = 1.683
    REFERENCE_OC = 1.0
    REFERENCE_GG = 6.3e-5

    _w, U = nm.dressed_ground(Bvec)
    H = nm.Hes(Bvec, d=D_STRAIN)
    dp, dc = nm.legs(U)
    dp = np.asarray(dp, dtype=complex).flatten()
    dc = np.asarray(dc, dtype=complex).flatten()
    n = H.shape[0]

    beta = (nm.TWOPI * REFERENCE_OC) ** 2 / 4
    g_eff = 2 * REFERENCE_GG + 2e-6

    def B_of_z(z, _H=H, _n=n):
        return 1j * nm.TWOPI * (_H - z * np.eye(_n))

    spec = gao.SectorSpec(
        D=np.eye(n, dtype=complex), B_of_z=B_of_z, dp=dp, dc=dc,
        g_eff=g_eff, beta=beta, readout_mode="bilinear",
        label="NV EIT (physical, nu_obs=4)",
    )
    spec.meta["gamma_phys_300K"] = float(nm.gamma_oc_GHz(T_boundary_K, D_STRAIN))
    spec.meta["gamma_phys_10K"] = float(nm.gamma_oc_GHz(10.0, D_STRAIN))
    return spec
