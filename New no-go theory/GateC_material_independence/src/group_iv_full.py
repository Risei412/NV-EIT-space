"""group_iv_full.py -- FULL physical GKSL (Lindblad) master equation for the
group-IV orbital-Lambda color centers SiV-/SnV-, for PRL Gate C Priority 4.

Certifies the Gamma^-1 suppression class (M0 = p^dag c != 0, "same-spin
orbital-Lambda") in a genuine density-matrix master equation, validated to
machine precision against the existing reduced kernel
(No-go theorem/src/group_iv_model.py).

Design (Gate-B-isomorphic). Instead of a steady-state EIT spectrum, we use the
fixed-rho0 weak-probe linear response of GateB_superconducting_witness, whose
coherence sector equals the reduced amplitude pencil exactly. The reduced
kernel K(Gamma) = p^dag [Gamma D + A0]^-1 c (group_iv_model) reappears as the
optical-coherence sector rho_{e,g} of the full Lindbladian.

Hilbert space (Layout A, N=5):
    |g>=0, and the 4 excited states 1..4 = group_iv_model.H_groupIV in its
    kron(orbital{e_x,e_y}, spin{up,dn}) ordering:
        1 = e_x,up   2 = e_x,dn   3 = e_y,up   4 = e_y,dn

Fast dissipation (the swept Gamma):
    mode='dephasing' : L_j = sqrt(2 Gamma) |e_j><e_j|  for j=1..4.
        Pure dephasing damps every coherence rho_{e_j,g} at (2 Gamma)/2 = Gamma,
        uniformly and WITHOUT transfer, so the coherence sub-block is exactly
        Gamma*I + A0 == group_iv_model's A_Gamma with D = I. (uniform D, unlike
        the singular bus D of Gate B.)
    mode='hopping'   : L = sqrt(Gamma) |e_y,s><e_x,s| + reverse (spin-conserving
        phonon orbital hopping, the physical process). D_hop is full-rank, so
        the leading order is still Gamma^-1 (coefficient shifts to
        p^dag D_hop^-1 c) -- confirming the class is a LEG property (M0 != 0),
        not an artifact of D = I.

Radiative decay: L = sqrt(2 GAMMA_RAD) |g><e_j|, so the coherence damping equals
group_iv_model.GAMMA_RAD = 0.0157 (an HWHM). The GKSL refill L rho L^dag lands
on |g><g| population, outside the coherence sector, leaving the kernel match
exact.

Physical parameter anchors (GHz, ordinary frequency; used only as guide lines
on the Gamma sweep -- the Gamma^-1 exponent and the M0 coefficient are
independent of the absolute phonon rate, which is swept over decades):
    SiV-: excited SO splitting Delta_e = 255 GHz; ground SO splitting ~48 GHz;
          ground orbital relaxation T1 ~ 39 ns @ 5 K  -> ~0.026 GHz;
          excited radiative lifetime ~1.7 ns          -> ~0.59 GHz.
    SnV-: Delta_e = 3000 GHz; ground SO splitting ~850 GHz.
Sources: Pingault et al., Nat. Commun. 8, 15579 (2017); Jahnke et al.,
New J. Phys. 17, 043011 (2015) [arXiv:1411.2871]; Trusheim et al., PRX 11,
041041 (2021); Meesala et al., PRB 97, 205444 (2018). (The parameter file
group_iv_model.py cites, SiV_SnV_phonon_AIC_parameters.md, is absent from the
tree; the relevant numbers are the Delta_e in group_iv_model.PARAMS plus the
anchors above.)

Reuses group_iv_model (H_groupIV, legs, moments, kernel, GAMMA_RAD) and
liouvillian_core conventions; numpy only (no qutip).
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_NOGO_SRC = os.path.join(_HERE, "..", "..", "..", "No-go theorem", "src")
if _NOGO_SRC not in sys.path:
    sys.path.insert(0, _NOGO_SRC)

import group_iv_model as giv  # noqa: E402

N_FULL = 5  # |g> + 4 excited

# physical guide-line anchors (GHz); see module docstring
ANCHORS = {
    "SiV": dict(Delta_e=255.0, ground_orbital_relax=0.026, excited_radiative=0.59),
    "SnV": dict(Delta_e=3000.0, ground_orbital_relax=None, excited_radiative=None),
}

# excited full-index map for hopping (orbital x spin)
_EX = {"e_x": {0: 1, 1: 2}, "e_y": {0: 3, 1: 4}}  # spin 0=up,1=dn


def vec(rho):
    return rho.reshape(-1)


def unvec(v):
    return v.reshape(N_FULL, N_FULL)


def build_H_full(material, xi_x=0.0, xi_y=0.0, Bx=0.0, z=0.0):
    """5x5 Hamiltonian: |g> plus the group_iv_model excited manifold (no 2pi,
    matching group_iv_model.A0_of = 1j*(H - z*I))."""
    Hg = giv.H_groupIV(material, xi_x=xi_x, xi_y=xi_y, Bx=Bx)
    H = np.zeros((N_FULL, N_FULL), dtype=complex)
    H[1:5, 1:5] = Hg - z * np.eye(4)
    return H


def dipole_legs_full(theta=0.0, spin=None):
    """Embed group_iv_model.legs(theta) into the N=5 excited block.
    Returns (p_full, c_full, Vp_probe): p/c length-5 readout/source on the
    excited coherence, Vp_probe = |c-leg><g| + h.c."""
    spin = giv.SPIN_UP if spin is None else spin
    pg, cg = giv.legs(theta, spin)  # length-4 excited vectors
    p_full = np.zeros(N_FULL, dtype=complex); p_full[1:5] = pg
    c_full = np.zeros(N_FULL, dtype=complex); c_full[1:5] = cg
    # probe couples the source excited ket (c, = e_x,up at theta=0) to |g>
    Vp = np.zeros((N_FULL, N_FULL), dtype=complex)
    Vp[1:5, 0] = cg
    Vp[0, 1:5] = np.conj(cg)
    return p_full, c_full, Vp


def jump_operators(Gamma, material=None, mode="dephasing", grad_pop=None):
    """Fast excited dissipation (dephasing or hopping) + radiative decay."""
    if grad_pop is None:
        grad_pop = 2 * giv.GAMMA_RAD
    ops = []
    if mode == "dephasing":
        for j in range(1, 5):
            L = np.zeros((N_FULL, N_FULL), dtype=complex)
            L[j, j] = np.sqrt(2 * Gamma)
            ops.append(L)
    elif mode == "hopping":
        for s in (0, 1):  # spin up, dn (spin-conserving hop)
            ex, ey = _EX["e_x"][s], _EX["e_y"][s]
            Lxy = np.zeros((N_FULL, N_FULL), dtype=complex); Lxy[ey, ex] = np.sqrt(Gamma)
            Lyx = np.zeros((N_FULL, N_FULL), dtype=complex); Lyx[ex, ey] = np.sqrt(Gamma)
            ops.extend([Lxy, Lyx])
    else:
        raise ValueError(f"unknown mode {mode!r}")
    for j in range(1, 5):  # radiative excited -> ground
        L = np.zeros((N_FULL, N_FULL), dtype=complex)
        L[0, j] = np.sqrt(grad_pop)
        ops.append(L)
    return ops


def liouvillian_block(H, jump_ops):
    """Convention-free basis-loop superoperator (as in Gate B /
    model_metro_lindblad): d/dt vec(rho) = L vec(rho), row-major vec."""
    dim = N_FULL * N_FULL
    L = np.zeros((dim, dim), dtype=complex)
    for idx in range(dim):
        e = np.zeros(dim, dtype=complex); e[idx] = 1.0
        rho = unvec(e)
        drho = -1j * (H @ rho - rho @ H)
        for C in jump_ops:
            Cd = C.conj().T
            drho += C @ rho @ Cd - 0.5 * (Cd @ C @ rho + rho @ Cd @ C)
        L[:, idx] = vec(drho)
    return L


def full_response(Gamma, material, theta=0.0, z=0.0, Bx=0.0,
                  xi_x=0.0, xi_y=0.0, mode="dephasing"):
    """Full-GKSL orbital-Lambda transfer response R(Gamma) = p^dag rho1[e,g],
    via the fixed-rho0 (=|g><g|) weak-probe linear response. Scales as
    Gamma^-1 because M0 = p^dag c != 0. (Equals the reduced kernel up to the
    trivial -i drive factor; see full_response returns and reduced_kernel_response.)"""
    H = build_H_full(material, xi_x=xi_x, xi_y=xi_y, Bx=Bx, z=z)
    p_full, c_full, Vp = dipole_legs_full(theta)
    L0 = liouvillian_block(H, jump_operators(Gamma, material, mode=mode))

    rho0 = np.zeros((N_FULL, N_FULL), dtype=complex); rho0[0, 0] = 1.0
    rhs = 1j * (Vp @ rho0 - rho0 @ Vp)

    M = L0.copy(); b = rhs.reshape(-1).copy()
    tr = np.zeros(N_FULL * N_FULL, dtype=complex)
    for i in range(N_FULL):
        tr[i * N_FULL + i] = 1.0
    M[0, :] = tr; b[0] = 0.0
    x = np.linalg.solve(M, b).reshape(N_FULL, N_FULL)
    return complex(np.conj(p_full[1:5]) @ x[1:5, 0])


def reduced_kernel_response(Gamma, material, theta=0.0, z=0.0, Bx=0.0,
                            xi_x=0.0, xi_y=0.0):
    """Reduced amplitude kernel K(Gamma) = p^dag [Gamma D + A0]^-1 c
    (group_iv_model.kernel), D = I."""
    H = giv.H_groupIV(material, xi_x=xi_x, xi_y=xi_y, Bx=Bx)
    return complex(giv.kernel(H, np.array([Gamma]), theta=theta, z=z)[0])


def M0(material, theta=0.0):
    """Leading path moment M0 = p^dag c (the Gamma^-1 coefficient); != 0 for
    the same-spin orbital-Lambda (theta small)."""
    H = giv.H_groupIV(material)
    return complex(giv.moments(H, 0, theta=theta)[0])


def verify_rho0_steady(Gamma, material, mode="dephasing"):
    H = build_H_full(material)
    L0 = liouvillian_block(H, jump_operators(Gamma, material, mode=mode))
    rho0 = np.zeros((N_FULL, N_FULL), dtype=complex); rho0[0, 0] = 1.0
    return float(np.max(np.abs(L0 @ rho0.reshape(-1))))
