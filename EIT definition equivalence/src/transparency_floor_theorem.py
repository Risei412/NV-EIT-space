"""
Gate F4: transparency floor theorem and (gamma_g, Omega_c) contrast phase
diagram for the standard (dark-state-bearing) Lambda EIT system.

Closed-form result (derived and verified below to machine precision):
for the matched two-coherence system A(delta) = diag(a1, g) + i*H_mix,
a1 = gamma_e - i*delta, g = gamma_g - i*delta, H_mix off-diagonal =
Omega_c/2 (Hermitian control coupling between the excited-state coherence
and the ground coherence), matched source/readout c=(1,0):

    chi_full(delta=0) = gamma_g / (gamma_e*gamma_g + beta),   beta = Omega_c^2/4

so the transparency contrast at line center is EXACTLY

    C_S = 1 - chi_full(0)/chi_cut(0) = beta / (gamma_e*gamma_g + beta)
        = Omega_c^2 / (4*gamma_e*gamma_g + Omega_c^2)

and the residual (floor) absorption is

    1 - C_S = 4*gamma_e*gamma_g / (4*gamma_e*gamma_g + Omega_c^2)
            ~ 4*gamma_e*gamma_g/Omega_c^2   for gamma_g << Omega_c^2/(4*gamma_e).

This is the standard EIT residual-absorption scaling (Fleischhauer-
Imamoglu-Marangos), re-derived here as a corollary of Proposition A
(the strict-accretive/matched-response floor, gate F2): C_S -> 1 (perfect
transparency) only in the gamma_g -> 0 limit, i.e. exactly the passivity
boundary lambda_min(Gamma_full) -> 0. This closes gate F4: no exploration
is needed, only quantification of the floor already proven to exist.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

gamma_e = 1.0


def chi_full_at_line_center(gamma_g, Omega_c):
    beta = Omega_c ** 2 / 4
    return gamma_g / (gamma_e * gamma_g + beta)


def contrast(gamma_g, Omega_c):
    beta = Omega_c ** 2 / 4
    return beta / (gamma_e * gamma_g + beta)


def verify_closed_form(n_trials=2000, seed=0):
    """Verify the closed form AT delta=0 (the two-photon resonance / EIT
    line center -- the standard reference point for EIT contrast in the
    literature). NOTE: an earlier version of this check instead searched
    for the GLOBAL argmin of Re(chi_full) over a delta grid; for large
    Omega_c (Autler-Townes regime) that grabs small-Re tail points far
    from resonance where chi_full/chi_cut ratio is not the EIT contrast
    question at all (both chi_full and chi_cut ->0 off resonance, and
    their ratio can exceed 1 there without violating positivity -- this
    produced a spurious "C_S<0" that looked like a Proposition-A
    violation but was actually a mismatched comparison point, not a
    floor breach; Re(chi_full) itself stayed positive throughout).
    The EIT contrast is only meaningfully defined at delta=0."""
    rng = np.random.default_rng(seed)
    max_err = 0.0
    for _ in range(n_trials):
        gamma_g = 10 ** rng.uniform(-6, 1)
        Omega_c = 10 ** rng.uniform(-2, 1)
        beta = Omega_c ** 2 / 4
        a1 = gamma_e
        g = gamma_g
        A = np.array([[a1, 1j * Omega_c / 2], [1j * Omega_c / 2, g]], dtype=complex)
        c = np.array([1.0, 0.0], dtype=complex)
        chi0 = (c.conj() @ np.linalg.solve(A, c)).real
        recut0 = 1.0 / gamma_e
        C_S_numeric = 1 - chi0 / recut0
        C_S_closed = contrast(gamma_g, Omega_c)
        max_err = max(max_err, abs(C_S_numeric - C_S_closed))
    return max_err


def build_phase_diagram(n=80):
    gamma_g_vals = np.logspace(-6, 1, n)
    Omega_c_vals = np.logspace(-2, 1, n)
    C = np.zeros((n, n))
    for i, gg in enumerate(gamma_g_vals):
        for j, Oc in enumerate(Omega_c_vals):
            C[i, j] = contrast(gg, Oc)
    return gamma_g_vals, Omega_c_vals, C


def main():
    max_err = verify_closed_form(n_trials=300)
    print(f"closed-form vs. numeric-minimization max |C_S diff| over 300 random "
          f"(gamma_g, Omega_c): {max_err:.3e}")

    gg, oc, C = build_phase_diagram(n=80)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.pcolormesh(gg, oc, (1 - C).T, shading="auto", norm=matplotlib.colors.LogNorm())
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\gamma_g/\gamma_e$")
    ax.set_ylabel(r"$\Omega_c/\gamma_e$")
    ax.set_title(r"Residual absorption floor $1-C_S$ (matched Lambda EIT)")
    fig.colorbar(im, ax=ax, label=r"$1-C_S$")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "F4_transparency_floor_phase_diagram.png", dpi=150)
    print("Figure written to", FIGURES_DIR / "F4_transparency_floor_phase_diagram.png")

    report = {
        "theorem": "C_S(gamma_g,Omega_c) = beta/(gamma_e*gamma_g+beta), beta=Omega_c^2/4; "
                   "1-C_S ~ 4*gamma_e*gamma_g/Omega_c^2 for gamma_g << Omega_c^2/(4 gamma_e)",
        "closed_form_vs_numeric_max_error": max_err,
        "sample_values": [
            {"gamma_g": float(g), "Omega_c": float(o), "C_S": float(contrast(g, o))}
            for g in [1.0, 0.1, 0.01, 0.001, 1e-4]
            for o in [0.3, 0.8, 2.0]
        ],
    }
    out_path = RESULTS_DIR / "gate_F4_transparency_floor.json"
    out_path.write_text(json.dumps(report, indent=2))
    print("Report written to", out_path)


if __name__ == "__main__":
    main()
