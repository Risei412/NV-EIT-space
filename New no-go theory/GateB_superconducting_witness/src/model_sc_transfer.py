"""model_sc_transfer.py -- superconducting dissipative-state-transfer witness
for PRL Gate B (non-EIT, non-diamond physical witness + blind prediction).

Physical system (circuit QED). An input transmon qubit A and an output
transmon qubit B are coupled through TWO lossy bus modes R1, R2 (a lossy
resonator with two accessible modes, or two lossy resonators). The bus decay
rate kappa is the fast dissipation scale Gamma. The observable is the A -> B
state-TRANSFER amplitude / efficiency -- a conversion/transport quantity, NOT
an EIT dark-state transparency.

Coherent graph (single-excitation manifold):
    A --gA1-- R1 --gB1-- B          two interfering transfer paths
    A --gA2-- R2 --gB2-- B
Fast dissipation acts ONLY on the bus (jump sqrt(kappa)|g><R1|, |g><R2|), so
the reduced dissipation operator D is SINGULAR (rank 2 on the bus, 0 on the
qubits). Weak fixed qubit relaxation gamma acts on A, B (the Gamma-independent
floor). Adiabatic elimination of the lossy bus gives an effective A->B
coupling (self-energy)

    Sigma(kappa) = sum_i gAi gBi / (kappa + i Delta_i + gamma_i),

so the transfer amplitude K(kappa) = p_B^dag [kappa D + A0]^-1 c_A scales as:
  * generic  (sum_i gAi gBi != 0):     K ~ kappa^-1   (nu_K = 1, efficiency 2)
  * protected(sum_i gAi gBi  = 0,       K ~ kappa^-2   (nu_K = 2, efficiency 4)
              Delta_1 != Delta_2):
  * broken   (small imbalance eps):     nu_K = 1 asymptotically, with a
                                        high->low crossover (Gamma_cross ~ 1/eps).
(Fully degenerate Delta_1 = Delta_2 with sum gAi gBi = 0 gives an EXACT
all-order cancellation -> Class I, nu = infinity; recorded as an optional
schematic data point, not the main gate.)

Literature-grounded parameters (GHz, ordinary frequency; provenance below):
  qubit-resonator coupling  g/2pi ~ 0.1 GHz   (50-300 MHz)
  bus decay (Gamma knob)    kappa/2pi          tunable, 0.1-50 MHz physical
  qubit relaxation/dephasing gamma ~ 2.5e-5 GHz (T1 ~ 40 us => 1/T1 ~ 25 kHz)
  detunings                 Delta_i            few tens of MHz
Sources: Krantz et al., Appl. Phys. Rev. 6, 021318 (2019); Blais et al.,
Rev. Mod. Phys. 93, 025005 (2021); lossy-bus adiabatic-elimination transfer
Burkhart et al., PRX Quantum 2, 030321 (2021), arXiv:2005.12334, 2403.02203.

This module mirrors the code structure of New no-go theory/src/model_physical.py
(reduced pencil + full vectorized Liouvillian) and model_metro_lindblad.py
(convention-free basis-loop liouvillian). It uses numpy/sympy only (no qutip).

Conventions:
  * Reduced pencil basis (single-excitation amplitudes): (A, R1, R2, B).
  * The coherence-sector decay rate is HALF the population jump rate
    (standard: a jump sqrt(r)|g><j| damps the coherence rho_{j,g} at r/2), so
    the reduced pencil uses kappa/2 on the bus and gamma/2 on the qubits; this
    makes the reduced amplitude matrix EQUAL the full Liouvillian's
    coherence-sector block (verified to ~1e-12 in run_gate_b.py).
  * Full GKSL basis: |g>=0, |A>=1, |R1>=2, |R2>=3, |B>=4 (single-excitation
    hardcore truncation; the linear response closes on one excitation).
"""
from __future__ import annotations

import numpy as np
import sympy as sp

# ----------------------------------------------------------------------
# parameters (GHz)
# ----------------------------------------------------------------------
G = 0.1          # base qubit-bus coupling, ~100 MHz
DELTA1 = 0.03    # R1 detuning, ~30 MHz
DELTA2 = 0.05    # R2 detuning, ~50 MHz (Delta1 != Delta2 keeps protected finite)
GAMMA_Q = 2.5e-5  # qubit relaxation floor, T1 ~ 40 us
DELTA_A = 0.0
DELTA_B = 0.0


def couplings(tuning="generic", eps=0.0):
    """Return (gA1, gA2, gB1, gB2) for the requested selection-rule class.

    generic   : sum gAi gBi = 2 g^2 != 0            -> nu_K = 1
    protected : gB2 = -gB1 so sum gAi gBi = 0        -> nu_K = 2
    broken    : protected + imbalance eps on gB2     -> nu_K = 1, crossover
    """
    if tuning == "generic":
        return (G, G, G, G)
    if tuning == "protected":
        return (G, G, G, -G)
    if tuning == "broken":
        return (G, G, G, -G + eps)
    raise ValueError(f"unknown tuning {tuning!r}")


# ----------------------------------------------------------------------
# reduced single-excitation amplitude pencil  A(kappa,z) = kappa*D + A0(z)
# basis (A, R1, R2, B); coherence-sector rates = half the jump rates
# ----------------------------------------------------------------------
def _H_exc(gA1, gA2, gB1, gB2):
    """4x4 Hermitian single-excitation Hamiltonian, basis (A,R1,R2,B)."""
    H = np.zeros((4, 4), dtype=complex)
    H[0, 0] = DELTA_A
    H[1, 1] = DELTA1
    H[2, 2] = DELTA2
    H[3, 3] = DELTA_B
    H[0, 1] = gA1; H[1, 0] = gA1
    H[0, 2] = gA2; H[2, 0] = gA2
    H[1, 3] = gB1; H[3, 1] = gB1
    H[2, 3] = gB2; H[3, 2] = gB2
    return H


def D_reduced():
    """Singular fast-dissipation operator (coherence rate = kappa/2 on bus)."""
    return np.diag([0.0, 0.5, 0.5, 0.0]).astype(complex)


def floor_reduced():
    """Gamma-independent qubit decay floor (coherence rate = gamma/2)."""
    return np.diag([GAMMA_Q / 2, 0.0, 0.0, GAMMA_Q / 2]).astype(complex)


def c_source():
    """Drive leg: weak probe creates the input-qubit A amplitude."""
    return np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)


def p_readout():
    """Readout leg: output-qubit B transfer amplitude."""
    return np.array([0.0, 0.0, 0.0, 1.0], dtype=complex)


def A0_reduced(z, tuning="generic", eps=0.0):
    gA1, gA2, gB1, gB2 = couplings(tuning, eps)
    H = _H_exc(gA1, gA2, gB1, gB2)
    return floor_reduced() + 1j * (H - z * np.eye(4))


def A_reduced(kappa, z=0.0, tuning="generic", eps=0.0):
    return kappa * D_reduced() + A0_reduced(z, tuning, eps)


def transfer_kernel(kappa, z=0.0, tuning="generic", eps=0.0):
    """K(kappa) = p_B^dag [kappa D + A0]^-1 c_A (transfer amplitude)."""
    A = A_reduced(kappa, z, tuning, eps)
    return complex(np.conj(p_readout()) @ np.linalg.solve(A, c_source()))


def B_of_z_factory(tuning="generic", eps=0.0):
    """Return a callable z -> A0(z) for core.certificate/moment helpers, and
    the matched D, so K = p^dag [kappa D + B_of_z(z)]^-1 c."""
    def B_of_z(z):
        return A0_reduced(z, tuning, eps)
    return B_of_z


# ----------------------------------------------------------------------
# exact symbolic pencil (for core.certificate_deg_nu; valid for singular D)
# ----------------------------------------------------------------------
def _rat(x):
    return sp.nsimplify(x, rational=True)


def symbolic_pencil(tuning="generic", eps=0.0):
    """Return (D_sym, B_sym_of_z, c_sym, p_sym) with exact rational entries."""
    gA1, gA2, gB1, gB2 = (float(v) for v in couplings(tuning, eps))
    gA1, gA2, gB1, gB2 = map(_rat, (gA1, gA2, gB1, gB2))
    d1, d2, gq = _rat(DELTA1), _rat(DELTA2), _rat(GAMMA_Q)

    H = sp.zeros(4, 4)
    H[1, 1] = d1; H[2, 2] = d2
    H[0, 1] = gA1; H[1, 0] = gA1
    H[0, 2] = gA2; H[2, 0] = gA2
    H[1, 3] = gB1; H[3, 1] = gB1
    H[2, 3] = gB2; H[3, 2] = gB2

    D_sym = sp.diag(0, sp.Rational(1, 2), sp.Rational(1, 2), 0)
    floor = sp.diag(gq / 2, 0, 0, gq / 2)

    def B_sym_of_z(z):
        zz = _rat(z) if z != 0 else sp.Integer(0)
        return floor + sp.I * (H - zz * sp.eye(4))

    c_sym = sp.Matrix([1, 0, 0, 0])
    p_sym = sp.Matrix([0, 0, 0, 1])
    return D_sym, B_sym_of_z, c_sym, p_sym


# ----------------------------------------------------------------------
# full GKSL (5-level single-excitation Lindblad)
# basis |g>=0, |A>=1, |R1>=2, |R2>=3, |B>=4
# ----------------------------------------------------------------------
N_FULL = 5


def vec(rho):
    return rho.reshape(-1)


def unvec(v):
    return v.reshape(N_FULL, N_FULL)


def H_full(kappa_unused=None, z=0.0, tuning="generic", eps=0.0, drive=0.0):
    gA1, gA2, gB1, gB2 = couplings(tuning, eps)
    H = np.zeros((N_FULL, N_FULL), dtype=complex)
    # excited-state (drive-frame) detunings, shifted by the probe detuning z
    H[1, 1] = DELTA_A - z
    H[2, 2] = DELTA1 - z
    H[3, 3] = DELTA2 - z
    H[4, 4] = DELTA_B - z
    # weak coherent drive on the input qubit A: drive*(|g><A| + h.c.)
    H[0, 1] = drive; H[1, 0] = drive
    # transfer couplings
    H[1, 2] = gA1; H[2, 1] = gA1
    H[1, 3] = gA2; H[3, 1] = gA2
    H[2, 4] = gB1; H[4, 2] = gB1
    H[3, 4] = gB2; H[4, 3] = gB2
    return H


def jump_operators(kappa):
    """Fast bus decay sqrt(kappa)|g><Rj| plus weak fixed qubit relaxation
    sqrt(gamma)|g><A|, |g><B| (all return to ground)."""
    ops = []
    for j in (2, 3):  # R1, R2
        L = np.zeros((N_FULL, N_FULL), dtype=complex)
        L[0, j] = np.sqrt(kappa)
        ops.append(L)
    for j in (1, 4):  # A, B qubit relaxation
        L = np.zeros((N_FULL, N_FULL), dtype=complex)
        L[0, j] = np.sqrt(GAMMA_Q)
        ops.append(L)
    return ops


def liouvillian_block(H, jump_ops):
    """Convention-free basis-loop superoperator (copied from
    model_metro_lindblad.py): d/dt vec(rho) = L vec(rho), row-major vec."""
    dim = N_FULL * N_FULL
    L = np.zeros((dim, dim), dtype=complex)
    for idx in range(dim):
        e = np.zeros(dim, dtype=complex)
        e[idx] = 1.0
        rho = unvec(e)
        drho = -1j * (H @ rho - rho @ H)
        for C in jump_ops:
            Cd = C.conj().T
            drho += C @ rho @ Cd - 0.5 * (Cd @ C @ rho + rho @ Cd @ C)
        L[:, idx] = vec(drho)
    return L


def implicit_linear_response(kappa, z=0.0, tuning="generic", eps=0.0):
    """O(drive) linear response without finite-eps subtraction (mirrors
    model_physical.implicit_linear_response). rho0 = |g><g|; probe V = |g><A|
    + h.c.; solve L0 rho1 = i[V, rho0] with a trace(rho1)=0 constraint.

    Returns rho1 (5x5, traceless). The transfer amplitude to B is rho1[4,0]
    (the |B><g| coherence); the input amplitude is rho1[1,0]."""
    H0 = H_full(z=z, tuning=tuning, eps=eps, drive=0.0)
    Ls = jump_operators(kappa)
    L0 = liouvillian_block(H0, Ls)

    rho0 = np.zeros((N_FULL, N_FULL), dtype=complex)
    rho0[0, 0] = 1.0
    Vp = np.zeros((N_FULL, N_FULL), dtype=complex)
    Vp[0, 1] = 1.0; Vp[1, 0] = 1.0
    rhs = 1j * (Vp @ rho0 - rho0 @ Vp)

    M = L0.copy()
    b = rhs.reshape(-1).copy()
    trace_row = np.zeros(N_FULL * N_FULL, dtype=complex)
    for i in range(N_FULL):
        trace_row[i * N_FULL + i] = 1.0
    M[0, :] = trace_row
    b[0] = 0.0
    x = np.linalg.solve(M, b)
    return x.reshape(N_FULL, N_FULL)


def full_transfer_amplitude(kappa, z=0.0, tuning="generic", eps=0.0):
    """Full-GKSL A->B transfer amplitude = |B><g| coherence at O(drive)."""
    rho1 = implicit_linear_response(kappa, z, tuning, eps)
    return complex(rho1[4, 0])


def verify_rho0_steady(kappa, tuning="generic", eps=0.0):
    H0 = H_full(z=0.0, tuning=tuning, eps=eps, drive=0.0)
    Ls = jump_operators(kappa)
    L0 = liouvillian_block(H0, Ls)
    rho0 = np.zeros((N_FULL, N_FULL), dtype=complex); rho0[0, 0] = 1.0
    return float(np.max(np.abs(L0 @ rho0.reshape(-1))))
