"""Step 2 of the room-temperature SMRT no-go plan (plan Sec. 5, Step 2):
physical audit of the EIT operational sector cut on the FULL NV Liouvillian
(not the reduced 6-level kernel of Step 1), at the same low-temperature
positive-control geometry validated there (T=10 K, Bz=0.02 T, Bx=0, strain
d=1.683 GHz, probe on orbital branch X / ms=0, control on orbital branch Y /
ms=-1 -- No-go theorem/src/nv_model.py's convention).

The cut sector S = {rho_pc, rho_cp} (the two ground-Raman-coherence Liouville
coordinates) is the same one already used, on an ad hoc "zero the S<->X
Liouvillian blocks" basis, by gate2_candidate_full_vs_reduced.py, and shown
in tests/test_operational_cut_equivalence.py to be numerically identical to
the algebraic D_S = diag(0, I_S) construction of Theorem 0A / Corollary
"Algebraic-operational equivalence for the block cut" (New no-go
theory/Sector/theorem2B_operational_realization.tex), implemented in
New no-go theory/Sector/src/operational_cut.py. This script reuses that
same operational_cut.py machinery and runs the three checks the plan's
Step 2 requires (mirroring the already-validated Gates U1/U2/U5 of
Sector/src/run_gates_step3.py, now applied to the full many-level NV
Liouvillian instead of a 2-3 level toy model):

  Gate 2.1 (~U1): finite-kappa intervention chi_op(kappa) = p^dag (A +
    kappa D_S)^-1 c converges to the ideal Riesz-projector cut chi_cut^(S)
    as O(kappa^-1), and chi_cut^(S) agrees with gate2's own ad hoc cut.
  Gate 2.2 (~U2): kernel/implementation universality -- two admissible cut
    generators (same retained projector P_S, different nonzero eigenvalues)
    give the SAME ideal cut but DIFFERENT finite-kappa responses.
  Gate 2.3 (~U5): non-invasiveness (condition C4) -- D_S annihilates the
    UNDRIVEN (Oc=0, no probe) reference steady state, and the steady state
    of L0 + kappa*D_S stays pinned to that reference for every kappa (the
    cut does not perturb the physical state it is applied to).
  Gate 2.4: the cut does not remove the one-photon absorption background --
    chi_cut^(S) with the control field off (Oc=0) still reproduces the bare
    probe absorption, i.e. cutting S has no effect when there is no
    control-induced coherence to cut in the first place.

A(z) regularization: the full vectorized Liouvillian L is exactly singular
(one trace-preserving steady state), so a generic finite-kappa/ideal-cut
matrix inverse needs an invertible response operator. Following the plan's
own A_Gamma(z) = zI - L^(1) convention (Sec. 2.1), we use A(z) = i*z*I - L
at a small representative z (regularization only, not a physical spectral
scan -- Steps elsewhere scan the physical two-photon detuning d2) and
cross-validate against gate2's own DC (z=0, lstsq + trace-gauge) solver in
Gate 2.4 to confirm the regularization does not distort the physics.

Outputs: RoomT/results/gates_summary_step2.json,
         RoomT/results/figures/fig_step2_operational_cut_audit.png
"""
from __future__ import annotations
import os
import sys
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOMT_ROOT = os.path.dirname(HERE)
NOGO_SRC = os.path.join(ROOMT_ROOT, "..", "..", "No-go theorem", "src")
SECTOR_SRC = os.path.join(ROOMT_ROOT, "..", "Sector", "src")
sys.path.insert(0, NOGO_SRC)
sys.path.insert(0, SECTOR_SRC)

import gate2_candidate_full_vs_reduced as g2  # noqa: E402
from liouvillian_core import liouvillian, steady_state  # noqa: E402
import operational_cut as oc  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Same geometry as Step 1 (step1_low_temperature_validation.py): Bz=0.02 T
# splits ms=-1,0,+1 unambiguously, Bx=0 (no B_perp-induced mixing -- this is
# the ORDINARY branch-resolved spin-Lambda, not the B_perp PRL candidate of
# SIMULATION_PLAN.md), strain d=1.683 GHz, control on ms=-1.
BZ = 0.02
D_STRAIN = 1.683
CTRL_BRANCH = "-1"
PROBE_ORB = np.array([1.0, 0.0], complex)  # matches nv_model.ORB_X
CTRL_ORB = np.array([0.0, 1.0], complex)   # matches nv_model.ORB_Y
REFERENCE_T = 10.0     # K, matches Step 1's REFERENCE_T
REFERENCE_OC = 1.0     # GHz, matches Step 1's REFERENCE_OC
Z_REG = 1e-3           # GHz, regularizing shift (>> L's spectral gap ~3e-6)


def build(T, Oc, d2=0.0):
    H, Ls, Vp, dp, meta = g2.build_full(
        T, 0.0, BZ, d2, Oc=Oc, ctrl=CTRL_BRANCH,
        ppol=PROBE_ORB, cpol=CTRL_ORB, isc=False, d=D_STRAIN, phi=0.0, j0=3,
    )
    L = liouvillian(H, Ls)
    return H, Ls, Vp, dp, meta, L


def probe_source_readout(Vp, dp, rho0_vec, N, p_idx, OP=1e-5):
    """Matches gate2_candidate_full_vs_reduced.chi_pair()'s own convention:
    Hp = weak-probe perturbation Hamiltonian, V its commutator superoperator,
    c = V @ rho0 (probe-induced source on the driven steady state), p = det
    (readout on the probe ground -> excited coherence rho_{p_idx, e})."""
    TWOPI = g2.TWOPI
    Hp = TWOPI * 0.5 * OP * (Vp + Vp.conj().T)
    I = np.eye(N)
    V = -1j * (np.kron(I, Hp) - np.kron(Hp.T, I))
    c = V @ rho0_vec
    det = np.zeros(N * N, complex)
    for e, a in enumerate(dp):
        det[p_idx * N + (3 + e)] = np.conj(a)
    return c, det, OP


def sector_D_S(dim, S, values):
    D_S = np.zeros((dim, dim), complex)
    for s, v in zip(S, values):
        D_S[s, s] = v
    return D_S


def run():
    H, Ls, Vp, dp, meta, L_driven = build(REFERENCE_T, REFERENCE_OC)
    N, p_idx, c_idx = meta["N"], meta["p_idx"], meta["c_idx"]
    dim = N * N
    S = [c_idx * N + p_idx, p_idx * N + c_idx]
    X = [k for k in range(dim) if k not in S]

    rho0_driven_vec, res0, gap = steady_state(L_driven, n=N)
    c, p_readout, OP = probe_source_readout(Vp, dp, rho0_driven_vec, N, p_idx)

    A = 1j * Z_REG * np.eye(dim) - L_driven
    cond_A = np.linalg.cond(A)

    # --- Gate 2.4a: A(z) regularization cross-check against gate2's own
    # DC (z=0, lstsq + trace gauge) solver, at the SAME (Oc, d2) point.
    from liouvillian_core import first_order
    Hp = g2.TWOPI * 0.5 * OP * (Vp + Vp.conj().T)
    I = np.eye(N)
    Vsuper = -1j * (np.kron(I, Hp) - np.kron(Hp.T, I))
    x_lstsq, res_lstsq = first_order(L_driven, Vsuper, rho0_driven_vec)
    chi_full_lstsq = -2 * (p_readout.conj() @ x_lstsq) / OP
    chi_full_regularized = -2 * oc.core.transfer(A, p_readout, c) / OP
    regularization_rel_err = abs(chi_full_regularized - chi_full_lstsq) / abs(chi_full_lstsq)

    # --- Gate 2.1 (~U1): O(kappa^-1) convergence + agreement with gate2's
    # own ad hoc block-zeroed cut (test_operational_cut_equivalence.py's
    # construction, now on THIS system rather than the B_perp candidate).
    D_S = sector_D_S(dim, S, [1.0, 1.0])
    v = oc.validate_cut_generator(D_S)
    kappas = np.logspace(2, 6, 20)
    res_u1 = oc.compare_cut_equivalence(A, D_S, c, p_readout, kappas)
    chi_ideal = res_u1["chi_ideal"]

    Lc_gate2 = L_driven.copy()
    Lc_gate2[np.ix_(S, X)] = 0
    Lc_gate2[np.ix_(X, S)] = 0
    A_cut_gate2 = 1j * Z_REG * np.eye(dim) - Lc_gate2
    chi_cut_gate2 = np.conj(p_readout) @ np.linalg.solve(A_cut_gate2, c)
    gate2_agreement_rel_err = abs(chi_cut_gate2 - chi_ideal) / abs(chi_ideal)

    gate_2_1 = dict(
        semisimple=bool(v["semisimple"]),
        selective_damping_ok=bool(v["selective_damping_ok"]),
        fit_slope=res_u1["fit_slope"],
        gate2_agreement_rel_err=float(gate2_agreement_rel_err),
        pass_=bool(
            v["semisimple"] and v["selective_damping_ok"]
            and abs(res_u1["fit_slope"] + 1.0) < 0.05
            and gate2_agreement_rel_err < 1e-6
        ),
    )

    # --- Gate 2.2 (~U2): kernel/implementation universality.
    D_S2 = sector_D_S(dim, S, [3.0 + 0.7j, 3.0 - 0.7j])
    v2 = oc.validate_cut_generator(D_S2)
    same_P_S = np.allclose(v["P_S"], v2["P_S"], atol=1e-10)
    ideal_1 = oc.ideal_cut_response(A, D_S, c, p_readout)
    ideal_2 = oc.ideal_cut_response(A, D_S2, c, p_readout)
    ideal_diff = abs(ideal_1 - ideal_2)
    kappas_u2 = np.logspace(1, 4, 10)
    op1 = np.array([oc.operational_cut_response(A, D_S, c, p_readout, k) for k in kappas_u2])
    op2 = np.array([oc.operational_cut_response(A, D_S2, c, p_readout, k) for k in kappas_u2])
    finite_diff = np.abs(op1 - op2)
    gate_2_2 = dict(
        same_P_S=bool(same_P_S),
        ideal_limit_difference=float(ideal_diff),
        finite_kappa_difference_at_kappa10=float(finite_diff[0]),
        pass_=bool(same_P_S and ideal_diff < 1e-9 and finite_diff[0] > 1e-6 * abs(ideal_1)),
    )

    # --- Gate 2.3 (~U5): non-invasiveness on the UNDRIVEN reference state.
    _, _, _, _, meta0, L0 = build(REFERENCE_T, 0.0)
    rho0_undriven_vec, res0u, gap0 = steady_state(L0, n=N)
    C_S_rho0 = D_S @ rho0_undriven_vec
    noninvasive_norm = float(np.linalg.norm(C_S_rho0))

    A0 = 1j * Z_REG * np.eye(dim) - L0
    max_steady_state_dev = 0.0
    for k in [0.0, 1.0, 10.0, 100.0, 1e4, 1e6]:
        L0_cut = L0 + k * D_S
        rho_k_vec, res_k, gap_k = steady_state(L0_cut, n=N)
        rho_k_vec = rho_k_vec / np.trace(rho_k_vec.reshape(N, N))
        max_steady_state_dev = max(
            max_steady_state_dev,
            float(np.linalg.norm(rho_k_vec - rho0_undriven_vec)),
        )
    gate_2_3 = dict(
        noninvasive_C_S_rho0_norm=noninvasive_norm,
        max_steady_state_deviation_across_kappa=max_steady_state_dev,
        pass_=bool(noninvasive_norm < 1e-8 and max_steady_state_dev < 1e-6),
    )

    # --- Gate 2.4b: one-photon absorption background survives the cut
    # when there is no control-induced coherence to remove (Oc=0): the
    # cut should be a near-identity operation on the bare (Oc=0) response.
    H0, Ls0, Vp0, dp0, meta_b, L_bare = build(REFERENCE_T, 0.0)
    rho0_bare_vec, _, _ = steady_state(L_bare, n=N)
    c_bare, p_bare, _ = probe_source_readout(Vp0, dp0, rho0_bare_vec, N, p_idx)
    A_bare = 1j * Z_REG * np.eye(dim) - L_bare
    chi_full_bare = oc.core.transfer(A_bare, p_bare, c_bare)
    chi_cut_bare = oc.ideal_cut_response(A_bare, D_S, c_bare, p_bare)
    bare_background_rel_err = abs(chi_cut_bare - chi_full_bare) / abs(chi_full_bare)
    gate_2_4 = dict(
        regularization_rel_err_vs_gate2_lstsq=float(regularization_rel_err),
        bare_background_rel_err=float(bare_background_rel_err),
        cond_A=float(cond_A),
        pass_=bool(regularization_rel_err < 1e-3 and bare_background_rel_err < 1e-6),
    )

    gates = dict(
        gate_2_1_kappa_inverse_convergence=gate_2_1,
        gate_2_2_kernel_universality=gate_2_2,
        gate_2_3_noninvasiveness=gate_2_3,
        gate_2_4_regularization_and_background=gate_2_4,
    )
    overall_pass = all(g["pass_"] for g in gates.values())

    summary = dict(
        model="No-go theorem/src/gate2_candidate_full_vs_reduced.py build_full "
              "(full NV Liouvillian, N=9), Step 1 low-T geometry",
        T_K=REFERENCE_T, Bz_T=BZ, d_strain_GHz=D_STRAIN, Oc_GHz=REFERENCE_OC,
        z_reg_GHz=Z_REG, spectral_gap_of_L_GHz=float(gap),
        sector_S_liouville_indices=S,
        gates=gates,
        overall_pass=bool(overall_pass),
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step2.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure: kappa convergence (Gate 2.1/2.2) ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    diff1 = np.abs(res_u1["chi_op"] - chi_ideal)
    ax.loglog(kappas, diff1, "o-", label="|chi_op(kappa) - chi_ideal|")
    ref = diff1[0] * (kappas[0] / kappas)
    ax.loglog(kappas, ref, "k--", label="O(kappa^-1) reference")
    ax.set_xlabel("kappa"); ax.set_ylabel("deviation from ideal cut")
    ax.set_title("Gate 2.1: finite-kappa convergence"); ax.legend(fontsize=8)

    ax = axes[1]
    ax.loglog(kappas_u2, finite_diff, "s-")
    ax.set_xlabel("kappa"); ax.set_ylabel("|chi_op^(1) - chi_op^(2)|")
    ax.set_title("Gate 2.2: two cut generators, same P_S")
    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step2_operational_cut_audit.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps({k: v["pass_"] for k, v in gates.items()}, indent=2))
    print("overall_pass:", overall_pass)
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
