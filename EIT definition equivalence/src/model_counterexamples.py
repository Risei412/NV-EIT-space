"""
Gate E2: minimal counterexamples to the four non-implications of
Theorem III (strategy document, Part IV):

  C1. optical dark vector  =/=>  stationary pure dark state
  C2. stationary dark state  =/=>  observable EIT
  C3. CPT  =/=>  perfect response-zero EIT
  C4. observable coherent transparency  =/=>  dark state

C1 and C4 are not rebuilt here: they are already established, with
witnesses, by other gates in this project --
  - C1 is exactly the gamma_g>0 branch of gate E1
    (model_lambda_chain.py): ker(Omega) stays nontrivial for every
    gamma_g (witness_i always holds), while the stationary-dark-state
    residual (witness_ii) becomes nonzero as soon as gamma_g>0.
  - C4 is the 2g+2e dark-state-free model reproduced in gate E0
    (2g2e_package/): rank([d_p,d_c])=2 (no optical dark vector) yet a
    79.03% absorption suppression is observed.
This script builds C2 and C3 fresh.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

gamma_e = 1.0


def basis_op(dim, i, j):
    op = np.zeros((dim, dim), dtype=complex)
    op[i, j] = 1.0
    return op


def liouvillian(dim, H, jumps):
    I = np.eye(dim, dtype=complex)
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for Lk in jumps:
        LkdLk = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk) - 0.5 * (np.kron(I, LkdLk) + np.kron(LkdLk.T, I))
    return L


def steady_state(dim, L):
    w, V = np.linalg.eig(L)
    idx = np.argmin(np.abs(w))
    rho = V[:, idx].reshape(dim, dim, order="F")
    rho = (rho + rho.conj().T) / 2
    rho = rho / np.trace(rho).real
    return rho, w[idx]


def solve_linear_response(dim, L0, H_p_op, rho0, delta):
    I = np.eye(dim * dim, dtype=complex)
    comm = H_p_op @ rho0 - rho0 @ H_p_op
    rhs = 1j * comm.reshape(-1, order="F")
    Aop = L0 + 1j * delta * I
    vec_drho, *_ = np.linalg.lstsq(Aop, rhs, rcond=1e-12)
    return vec_drho.reshape(dim, dim, order="F")


# ---------------------------------------------------------------------
# C2: stationary dark state exists, but a genuinely decoupled probe
# transition sees no EIT feature at all (detector/channel insensitive).
#
# Basis: 0=g1, 1=g2, 2=e2 (the actual Lambda system: probe2 on g1-e2,
# control on g2-e2, forms the dark state), 3=e1 (a SEPARATE excited
# state probed by an independent "witness" field Omega_w on g1-e1, with
# NO coupling at all to g2 or the control field -- oblivious to the
# dark-state physics happening in the g1-g2-e2 sector).
# ---------------------------------------------------------------------

def run_C2(Omega_c=0.8, Omega_p2=0.8, Omega_w=0.3, gamma_g=0.0):
    dim = 4  # g1, g2, e2, e1
    H0 = np.zeros((dim, dim), dtype=complex)
    H0[2, 0] += Omega_p2 / 2; H0[0, 2] += np.conj(Omega_p2) / 2   # probe2: g1-e2
    H0[2, 1] += Omega_c / 2;  H0[1, 2] += np.conj(Omega_c) / 2    # control: g2-e2
    H0[3, 0] += Omega_w / 2;  H0[0, 3] += np.conj(Omega_w) / 2    # witness probe: g1-e1

    jumps = [
        np.sqrt(gamma_e / 2) * basis_op(dim, 0, 2),
        np.sqrt(gamma_e / 2) * basis_op(dim, 1, 2),
        np.sqrt(gamma_e) * basis_op(dim, 0, 3),   # e1 decays only to g1
    ]
    if gamma_g > 0:
        jumps.append(np.sqrt(gamma_g) * (basis_op(dim, 0, 0) - basis_op(dim, 1, 1)) / np.sqrt(2))

    # Dark-state / CPT diagnostics must be evaluated on the Lambda
    # subsystem ALONE (witness field off): the witness field, though it
    # never touches g2/control, still shares the g1 leg and so its
    # back-action would spoil the Lambda dark state if left on --
    # exactly the perturbation this counterexample is NOT about (we want
    # an undisturbed dark state whose EXISTENCE the witness detector
    # simply fails to reveal, not a witness field that destroys it).
    H0_bg = H0.copy()
    H0_bg[3, 0] = 0; H0_bg[0, 3] = 0  # remove witness field for background state
    L0_bg = liouvillian(dim, H0_bg, jumps)
    rho_bg, _ = steady_state(dim, L0_bg)

    rho_e2e2 = rho_bg[2, 2].real
    D = np.array([Omega_c, -Omega_p2, 0.0, 0.0], dtype=complex)
    D = D / np.linalg.norm(D)
    rho_D_full = np.outer(D, D.conj())
    dark_residual = np.linalg.norm(L0_bg @ rho_D_full.reshape(-1, order="F"))

    # witness probe response on the UNRELATED e1-g1 transition, linearized
    # on top of the (undisturbed) Lambda-subsystem dark state: should stay
    # a plain passive Lorentzian, unaffected by the dark-state formation.
    H_w_unit = basis_op(dim, 3, 0) * (Omega_w / 2) + basis_op(dim, 0, 3) * (np.conj(Omega_w) / 2)
    drho = solve_linear_response(dim, L0_bg, H_w_unit, rho_bg, 0.0)
    chi_witness = (1j / (Omega_w / 2)) * drho[3, 0]
    chi_bare = 1.0 / gamma_e  # plain passive Lorentzian at line center, Gamma_e1=gamma_e

    return {
        "rho_e2e2_ss": float(rho_e2e2),
        "dark_state_residual": float(dark_residual),
        "chi_witness": [float(chi_witness.real), float(chi_witness.imag)],
        "chi_bare_expected": float(chi_bare),
        "witness_matches_bare_Lorentzian": bool(abs(chi_witness.real - chi_bare) < 1e-6),
    }


# ---------------------------------------------------------------------
# C3: CPT holds for the Lambda-system excited state, but the total probe
# response stays finite because of an EXTRA, off-resonant background
# excited state e' coupled directly to g1 (bypassing the Raman
# interference entirely).
#
# Basis: 0=g1, 1=g2, 2=e (Lambda excited state), 3=e' (background state,
# large detuning, coupled to g1 only, no coupling to g2/control at all).
# The PROBE couples to BOTH e (resonantly) and e' (off-resonantly): the
# readout is the total g1-polarization, Tr[(|e><g1|+|e'><g1|) rho].
# ---------------------------------------------------------------------

def run_C3(Omega_c=0.8, Omega_p=0.8, Delta_bg=6.0, Omega_p_bg=0.5, gamma_g=0.0):
    dim = 4  # g1, g2, e, e'
    H0 = np.zeros((dim, dim), dtype=complex)
    H0[2, 0] += Omega_p / 2; H0[0, 2] += np.conj(Omega_p) / 2       # probe: g1-e
    H0[2, 1] += Omega_c / 2; H0[1, 2] += np.conj(Omega_c) / 2       # control: g2-e
    H0[3, 3] = Delta_bg                                             # e' detuning (background channel)
    H0[3, 0] += Omega_p_bg / 2; H0[0, 3] += np.conj(Omega_p_bg) / 2  # probe also weakly excites e'

    jumps = [
        np.sqrt(gamma_e / 2) * basis_op(dim, 0, 2),
        np.sqrt(gamma_e / 2) * basis_op(dim, 1, 2),
        np.sqrt(gamma_e) * basis_op(dim, 0, 3),  # e' decays back to g1 only
    ]
    if gamma_g > 0:
        jumps.append(np.sqrt(gamma_g) * (basis_op(dim, 0, 0) - basis_op(dim, 1, 1)) / np.sqrt(2))

    # background (probe off) for the linear-response calculation
    H0_bg = H0.copy()
    H0_bg[2, 0] = 0; H0_bg[0, 2] = 0
    H0_bg[3, 0] = 0; H0_bg[0, 3] = 0
    L0_bg = liouvillian(dim, H0_bg, jumps)
    rho_bg, _ = steady_state(dim, L0_bg)

    # CPT check: population in the Lambda-system excited state e (index 2)
    # should be ~0 at gamma_g=0 when the probe/control are on (true
    # nonlinear steady state including both probe and control coherently).
    H0_full = H0.copy()
    L0_full = liouvillian(dim, H0_full, jumps)
    rho_full, _ = steady_state(dim, L0_full)
    rho_ee = rho_full[2, 2].real

    # total probe response: sum of the two coherences (e-g1 and e'-g1),
    # each linearized in its own coupling amplitude.
    H_p_unit = (basis_op(dim, 2, 0) * (Omega_p / 2) + basis_op(dim, 0, 2) * (np.conj(Omega_p) / 2)
                + basis_op(dim, 3, 0) * (Omega_p_bg / 2) + basis_op(dim, 0, 3) * (np.conj(Omega_p_bg) / 2))
    drho = solve_linear_response(dim, L0_bg, H_p_unit, rho_bg, 0.0)
    chi_e = (1j / (Omega_p / 2)) * drho[2, 0] if Omega_p != 0 else 0.0
    chi_ebg = (1j / (Omega_p_bg / 2)) * drho[3, 0] if Omega_p_bg != 0 else 0.0
    # total normalized response (weighting by relative coupling strengths,
    # since both channels contribute additively to the same g1 coherence
    # when normalized by a COMMON reference amplitude, taken here as Omega_p)
    chi_total = drho[2, 0] * (1j / (Omega_p / 2)) + drho[3, 0] * (1j / (Omega_p / 2))

    return {
        "rho_ee_lambda_subsystem": float(rho_ee),
        "chi_e_lambda_channel": [float(chi_e.real), float(chi_e.imag)] if not isinstance(chi_e, float) else [0.0, 0.0],
        "chi_ebg_background_channel": [float(chi_ebg.real), float(chi_ebg.imag)] if not isinstance(chi_ebg, float) else [0.0, 0.0],
        "chi_total": [float(chi_total.real), float(chi_total.imag)],
        "CPT_holds_but_response_nonzero": bool(rho_ee < 1e-6 and abs(chi_total) > 1e-3),
    }


def main():
    print("=== C2: stationary dark state exists, witness channel sees no EIT ===")
    c2 = run_C2(gamma_g=0.0)
    print(json.dumps(c2, indent=2))
    assert c2["dark_state_residual"] < 1e-9 and abs(c2["rho_e2e2_ss"]) < 1e-9
    assert c2["witness_matches_bare_Lorentzian"]

    print("\n=== C3: CPT holds increasingly well as Omega_p_bg -> 0, "
          "while the response stays proportionally finite ===")
    c3_scaling = []
    for Op_bg in [0.5, 0.2, 0.1, 0.05, 0.02]:
        r = run_C3(Omega_p_bg=Op_bg, gamma_g=0.0)
        chi_mag = abs(complex(*r["chi_total"]))
        c3_scaling.append({
            "Omega_p_bg": Op_bg,
            "rho_ee_lambda": r["rho_ee_lambda_subsystem"],
            "rho_ee_over_Opbg_sq": r["rho_ee_lambda_subsystem"] / Op_bg ** 2,
            "chi_total_mag": chi_mag,
            "chi_total_over_Opbg": chi_mag / Op_bg,
        })
        print(f"  Omega_p_bg={Op_bg:.3f}: rho_ee(lambda)={r['rho_ee_lambda_subsystem']:.3e} "
              f"(~ Omega_p_bg^2, ratio {r['rho_ee_lambda_subsystem']/Op_bg**2:.4f}), "
              f"|chi_total|={chi_mag:.3e} (~ Omega_p_bg, ratio {chi_mag/Op_bg:.4f})")
    c3 = run_C3(Omega_p_bg=0.5, gamma_g=0.0)
    c3["scaling_with_Omega_p_bg"] = c3_scaling
    print(f"\n  As Omega_p_bg->0: rho_ee ~ Omega_p_bg^2 -> 0 (CPT increasingly exact) "
          f"while |chi_total| ~ Omega_p_bg -> 0 only linearly, i.e. the response/CPT-"
          f"violation RATIO diverges -- CPT does not imply a vanishing response.")

    report = {"C2": c2, "C3": c3,
              "C1_reference": "see results/gate_E1_lambda_chain.json (gamma_g>0 rows): "
                               "witness_i (ker Omega) always nontrivial, "
                               "witness_ii (dark-state residual) nonzero for gamma_g>0",
              "C4_reference": "see 2g2e_package/ (gate E0): rank([d_p,d_c])=2 "
                               "(no optical dark vector), 79.03% absorption suppression observed"}
    out_path = RESULTS_DIR / "gate_E2_counterexamples.json"
    out_path.write_text(json.dumps(report, indent=2))
    print("\nReport written to", out_path)


if __name__ == "__main__":
    main()
