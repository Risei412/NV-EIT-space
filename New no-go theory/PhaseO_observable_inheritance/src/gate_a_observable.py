"""gate_a_observable.py -- Gate A (Priority 1): observable-order inheritance.

Model-independent library that PREDICTS the suppression order nu_obs of a
physical observable from the orders of its constituent path-kernels, and
verifies the prediction three independent ways (moment predictor, exact
rational-degree certificate, full log-log slope), separating the true
asymptotic index from a pre-asymptotic effective exponent.

Central law (bilinear / EIT-type sector response). With the fast (optical)
resolvent G(Gamma,z) = [Gamma*D + A0(z)]^-1, probe leg dp, control leg dc,
cross-kernels

    K12 = dp^dag G dc,   K21 = dc^dag G dp,   S2 = dc^dag G dc,

and the adiabatically-eliminated sector observable

    R_obs(Gamma,z) = -beta * K12 * K21 / (g_eff + beta*S2),

the suppression order inherits as

    nu_obs = n12 + n21 - nu_den,
    nu_den = 0        (generic: g_eff != 0, slow denominator S_g -> g_eff),
    nu_den = nu_S2    (protected floor: g_eff = 0, S_g -> beta*S2).

n12, n21 are the cross-kernel orders (n1j = index of first nonzero path
moment + 1). Selection-rule cancellation (dp^dag dc = 0, ...) promotes n12,
n21 above the generic value 1; this is what turns nu_obs from the naive
single-kernel nu_K into the observable order.

Reconciles the three known repository data points:
  * NV EIT (RoomT step3):  n12=n21=2  (M0 = dp^dag dc = 0 selection rule)
                           => nu_obs = 4  (not the kernel order nu_K = 2).
  * QFI (Phase M):         quadratic readout |x_S|^2, n12 = n21 = nu
                           => nu_obs = 2*nu.
  * Phase N q-fan / z*:    intervention scaling and resonance zeros modify
                           nu_den / promote the order (handled by the exact
                           Phase N cores, cross-checked here).

Everything Gamma-asymptotic is decided by the EXACT rational-degree
certificate (SymPy), never by a float slope alone; the log-log fit is an
independent numeric corroboration and the source of the pre-asymptotic
crossover diagnostics. This mirrors the repository convention (RoomT step3,
New no-go theory/src/core.py).

This module reuses New no-go theory/src/core.py and does not modify it.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
_PHASE_SRC = os.path.join(_HERE, "..", "..", "src")
if _PHASE_SRC not in sys.path:
    sys.path.insert(0, _PHASE_SRC)
import core  # noqa: E402  (New no-go theory/src/core.py)

TOL = 1e-9


# ----------------------------------------------------------------------
# Sector specification
# ----------------------------------------------------------------------
@dataclass
class SectorSpec:
    """A reduced optical pencil A(Gamma,z) = Gamma*D + B_of_z(z) with a probe
    leg `dp`, a control leg `dc`, and the slow-sector floor (g_eff) and
    cross coupling (beta) of the adiabatically-eliminated observable.

    readout_mode:
      'bilinear'  -> R_obs = -beta*K12*K21/(g_eff+beta*S2)   (EIT susceptibility)
      'quadratic' -> R_obs = |K12|^2                          (QFI-type readout)

    The *_sym fields (SymPy) are optional; when present the exact
    rational-degree certificate can be evaluated. D must be invertible on
    the optical sector (the fast block), which is the case for every reduced
    pencil in this repository (radiative D = diag(rates) > 0).
    """

    D: np.ndarray
    B_of_z: Callable[[float], np.ndarray]
    dp: np.ndarray
    dc: np.ndarray
    g_eff: float = 0.0
    beta: float = 1.0
    readout_mode: str = "bilinear"
    label: str = ""

    D_sym: Optional[sp.Matrix] = None
    B_sym_of_z: Optional[Callable[[float], sp.Matrix]] = None
    dp_sym: Optional[sp.Matrix] = None
    dc_sym: Optional[sp.Matrix] = None
    meta: dict = field(default_factory=dict)


# ----------------------------------------------------------------------
# (i) Predictor -- cross-kernel orders and nu_obs from path moments
# ----------------------------------------------------------------------
def cross_kernel_orders(spec: SectorSpec, z: float = 0.0, kmax: int = 8,
                        tol: float = TOL) -> dict:
    """First-nonzero-moment orders of K12, K21, S2 (numeric moment ladder,
    Theorem II). Also flags an EXACT leading-moment cancellation (M0 == 0
    within tol) as a protected selection rule.

    n1j = (index of first nonzero moment) + 1  ==  nu of that kernel.
    """
    n12, mu12 = core.nu_from_moments(spec.D, spec.B_of_z, spec.dc, spec.dp, z, kmax, tol)
    n21, mu21 = core.nu_from_moments(spec.D, spec.B_of_z, spec.dp, spec.dc, z, kmax, tol)
    nS2, muS2 = core.nu_from_moments(spec.D, spec.B_of_z, spec.dc, spec.dc, z, kmax, tol)

    m0_overlap = complex(np.conj(spec.dp) @ spec.dc)
    protected = abs(m0_overlap) < tol  # dp^dag dc = 0 => leading moment of K12 vanishes

    return dict(
        n12=n12, n21=n21, nu_S2=nS2,
        mu12=mu12, mu21=mu21, muS2=muS2,
        m0_overlap=m0_overlap,
        protected_leading_cancellation=bool(protected),
    )


def predict_nu_obs(spec: SectorSpec, z: float = 0.0, kmax: int = 8,
                   tol: float = TOL) -> dict:
    """Predicted observable order from the inheritance law."""
    ck = cross_kernel_orders(spec, z=z, kmax=kmax, tol=tol)
    n12, n21, nu_S2 = ck["n12"], ck["n21"], ck["nu_S2"]

    if spec.readout_mode == "quadratic":
        # |K12|^2 ~ Gamma^{-2 n12}
        nu_den = 0
        nu_obs = 2 * n12 if n12 is not None else None
    else:
        geff_nonzero = abs(spec.g_eff) > 1e-15
        nu_den = 0 if geff_nonzero else nu_S2
        if n12 is None or n21 is None:
            nu_obs = None  # Class I candidate (all checked moments vanish)
        else:
            nu_obs = n12 + n21 - (nu_den or 0)

    if n12 is not None and n12 <= 1 and (n21 is None or n21 <= 1):
        mechanism = "generic"
    else:
        mechanism = "symmetry-protected"

    return dict(
        n12=n12, n21=n21, nu_S2=nu_S2, nu_den=nu_den,
        nu_obs_pred=nu_obs, mechanism=mechanism,
        protected=bool(ck["protected_leading_cancellation"]),
        m0_overlap=ck["m0_overlap"],
    )


# ----------------------------------------------------------------------
# Numeric observable and (ii) log-log verifier
# ----------------------------------------------------------------------
def observable_value(spec: SectorSpec, Gamma: float, z: float = 0.0) -> complex:
    A = Gamma * spec.D + spec.B_of_z(z)
    G = np.linalg.inv(A)
    K12 = np.conj(spec.dp) @ G @ spec.dc
    if spec.readout_mode == "quadratic":
        return complex(abs(K12) ** 2)
    K21 = np.conj(spec.dc) @ G @ spec.dp
    S2 = np.conj(spec.dc) @ G @ spec.dc
    den = spec.g_eff + spec.beta * S2
    return complex(-spec.beta * K12 * K21 / den)


def verify_nu_obs_loglog(spec: SectorSpec, gammas: np.ndarray, z: float = 0.0,
                         tail_frac: float = 0.4) -> dict:
    """Full numeric observable across a Gamma sweep, with a global slope and a
    deep-tail slope (the tail is the true asymptotic estimate, robust to the
    pre-asymptotic crossover)."""
    gammas = np.asarray(gammas, dtype=float)
    vals = np.array([observable_value(spec, g, z) for g in gammas])
    fit = core.fit_nu_loglog(gammas, vals)
    ntail = max(3, int(len(gammas) * tail_frac))
    fit_tail = core.fit_nu_loglog(gammas[-ntail:], vals[-ntail:])
    return dict(
        nu_global=fit["nu_global"],
        nu_tail=fit_tail["nu_global"],
        nu_eff=fit["nu_eff"],
        gamma_mid=fit["gamma_mid"],
        vals=vals,
        gammas=gammas,
    )


# ----------------------------------------------------------------------
# (iii) Exact rational-degree certificate
# ----------------------------------------------------------------------
def certify_nu_obs_exact(spec: SectorSpec, z_val: float = 0.0) -> dict:
    """Exact Gamma-degree certificate of nu_obs from the symbolic pencil.

    Builds Q = det A, N12 = dp^dag adj(A) dc, N21, N22 as exact SymPy
    polynomials in Gamma (A = Gamma*D_sym + B_sym(z)), then:

      bilinear, g_eff != 0:  nu_obs = 2*deg Q - deg N12 - deg N21
                                    = nu_K12 + nu_K21          (generic law)
      bilinear, g_eff == 0:  nu_obs = deg Q + deg N22 - deg N12 - deg N21
                                    = nu_K12 + nu_K21 - nu_S2  (protected floor)
      quadratic:             nu_obs = 2*(deg Q - deg N12)

    Requires spec.D_sym, spec.B_sym_of_z, spec.dp_sym, spec.dc_sym.
    """
    if spec.D_sym is None or spec.B_sym_of_z is None:
        raise ValueError("exact certificate requires symbolic pencil (D_sym, B_sym_of_z, dp_sym, dc_sym)")
    Gamma = sp.symbols("Gamma")
    A = Gamma * spec.D_sym + spec.B_sym_of_z(z_val)
    Q = sp.Poly(sp.expand(A.det()), Gamma)
    adjA = A.adjugate()
    dp, dc = spec.dp_sym, spec.dc_sym
    N12 = sp.Poly(sp.expand((dp.conjugate().T * adjA * dc)[0, 0]), Gamma)
    N21 = sp.Poly(sp.expand((dc.conjugate().T * adjA * dp)[0, 0]), Gamma)
    N22 = sp.Poly(sp.expand((dc.conjugate().T * adjA * dc)[0, 0]), Gamma)

    degQ = Q.degree()
    nu_K12 = degQ - (N12.degree() if N12.as_expr() != 0 else -sp.oo)
    nu_K21 = degQ - (N21.degree() if N21.as_expr() != 0 else -sp.oo)
    nu_S2 = degQ - (N22.degree() if N22.as_expr() != 0 else -sp.oo)

    if spec.readout_mode == "quadratic":
        nu = 2 * nu_K12
    else:
        deg_num = N12.degree() + N21.degree()
        if abs(spec.g_eff) > 1e-15:
            deg_den = degQ + max(degQ, N22.degree())
        else:
            deg_den = degQ + N22.degree()
        nu = deg_den - deg_num

    return dict(
        nu_cert=int(nu), nu_K12=int(nu_K12), nu_K21=int(nu_K21), nu_S2=int(nu_S2),
        deg_Q=int(degQ), deg_N12=int(N12.degree()), deg_N21=int(N21.degree()),
        deg_N22=int(N22.degree()),
    )


# ----------------------------------------------------------------------
# (iv) Pre-asymptotic vs asymptotic separator
# ----------------------------------------------------------------------
def separate_regimes(spec: SectorSpec, gammas: np.ndarray, z: float = 0.0,
                     gamma_phys: Optional[float] = None,
                     generic_order: Optional[float] = None,
                     asymptotic_order: Optional[float] = None) -> dict:
    """Local effective index nu_eff(Gamma) and the crossover scale where it
    passes the midpoint between the pre-asymptotic and true asymptotic
    orders. If gamma_phys is given, report nu_eff there and whether it is in
    the pre-asymptotic regime (gamma_phys < gamma_cross)."""
    v = verify_nu_obs_loglog(spec, gammas, z=z)
    nu_eff, gmid = v["nu_eff"], v["gamma_mid"]
    nu_asymptotic = v["nu_tail"]

    gamma_cross = None
    if generic_order is not None and asymptotic_order is not None and len(nu_eff):
        mid = 0.5 * (generic_order + asymptotic_order)
        idx = int(np.argmin(np.abs(nu_eff - mid)))
        gamma_cross = float(gmid[idx])

    nu_eff_at_phys = None
    is_preasymptotic = None
    if gamma_phys is not None and len(gmid):
        j = int(np.argmin(np.abs(gmid - gamma_phys)))
        nu_eff_at_phys = float(nu_eff[j])
        if gamma_cross is not None:
            is_preasymptotic = bool(gamma_phys < gamma_cross)

    return dict(
        nu_asymptotic=float(nu_asymptotic),
        nu_eff_small=float(nu_eff[0]) if len(nu_eff) else None,
        nu_eff_large=float(nu_eff[-1]) if len(nu_eff) else None,
        gamma_cross=gamma_cross,
        nu_eff_at_phys=nu_eff_at_phys,
        is_preasymptotic=is_preasymptotic,
    )


# ----------------------------------------------------------------------
# Top-level per-model consistency check (predict vs certificate vs fit)
# ----------------------------------------------------------------------
def check_model(spec: SectorSpec, gammas: np.ndarray, z: float = 0.0,
                kmax: int = 8, slope_tol: float = 0.06) -> dict:
    """Run all three routes on one SectorSpec and report agreement.

    'agree' requires: predictor == exact certificate (integers), and the
    deep-tail log-log slope within slope_tol of that integer.
    """
    pred = predict_nu_obs(spec, z=z, kmax=kmax)
    cert = None
    if spec.D_sym is not None:
        cert = certify_nu_obs_exact(spec, z_val=z)
    fit = verify_nu_obs_loglog(spec, gammas, z=z)

    nu_pred = pred["nu_obs_pred"]
    nu_cert = cert["nu_cert"] if cert else None
    nu_tail = fit["nu_tail"]

    pred_cert_agree = (nu_cert is None) or (nu_pred == nu_cert)
    ref = nu_cert if nu_cert is not None else nu_pred
    fit_agree = (ref is not None) and (abs(nu_tail - ref) < slope_tol)

    return dict(
        label=spec.label,
        readout_mode=spec.readout_mode,
        mechanism=pred["mechanism"],
        n12=pred["n12"], n21=pred["n21"], nu_den=pred["nu_den"],
        nu_obs_pred=nu_pred, nu_obs_cert=nu_cert, nu_obs_tail_fit=nu_tail,
        pred_cert_agree=bool(pred_cert_agree),
        fit_agree=bool(fit_agree),
        agree=bool(pred_cert_agree and fit_agree),
    )
