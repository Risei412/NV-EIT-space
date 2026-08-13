"""Gate F5-C: upgrade the matched-response floor from a random scan to a
proved identity plus an adversarial (directed) search.

Motivation
----------
`results/F5B_findings.md` certifies "no counterexample to the matched-response
floor" from 150 *random* parameter draws (`min Re chi_matched = 0.0029`).
Random sampling is weak evidence for a no-go: a positivity claim is only as
strong as the worst case, and random draws do not look for the worst case.
This gate replaces that evidence on both sides of the model boundary.

Part 1 --- where the floor is a theorem (coherence-block model).
For the reduced coherence block used by F0/F1/F4,

    A(delta) = Gamma + i (H - delta I),   Gamma = Gamma^dag >= 0,  H = H^dag,

the matched response is chi = e^dag A^{-1} e with source = readout = e.
Writing x = A^{-1} e,

    chi = (A x)^dag x       (because e = A x, and e^dag = (A x)^dag)
        = x^dag A^dag x
        = x^dag Gamma x - i x^dag (H - delta I) x,

so, since the second term is real-multiplied by -i,

    Re chi = x^dag Gamma x >= 0,

with strict inequality whenever Gamma > 0, and equality possible only if
x lies in ker Gamma.  This is a one-line proof of the F4/F5B "floor":
matched readout of a passive Markovian coherence block can never be
gain-like, and an exact transparency zero requires the response vector to
enter the kernel of the dissipator.  It also predicts exactly where the
no-go can be evaded: unmatched readout (source != readout), which is the
regime where F5-B found genuine exact zeros in the Raman channel.

Part 1 verifies this identity numerically, including the singular-Gamma
boundary where the bound is saturated.

Part 2 --- where it is not a theorem (full GKSL with populations).
The full density-matrix response is not of the form e^dag A^{-1} e, so the
identity does not transfer, and F5-B's evidence there is the random scan.
Part 2 runs a directed minimization of min_delta Re chi_matched over the
full GKSL parameter space (Nelder-Mead from many seeds, rates bounded below
by the same floor F5-B used) and reports the adversarial infimum next to the
random-scan baseline.  A directed search that still cannot go below the
random baseline's order of magnitude is much stronger evidence than the
random baseline alone; a directed search that reaches zero would falsify the
floor outside the coherence-block model, which is exactly the useful outcome.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

import full_gksl_2g3e as fg

RESULTS_DIR = Path(__file__).resolve().parents[4] / "results" / "certificates"

RATE_FLOOR = 1e-3          # same lower bound on gamma_g / gamma_relax as F5-B
DELTA_SCAN = np.linspace(-6.0, 6.0, 241)


# ----------------------------------------------------------------------
# Part 1: the accretivity identity on the coherence block
# ----------------------------------------------------------------------
def part1_identity(n_trials: int = 4000, seed: int = 20260731) -> dict:
    rng = np.random.default_rng(seed)
    max_identity_dev = 0.0
    min_re_chi = np.inf
    min_re_chi_singular = np.inf
    for _ in range(n_trials):
        n = int(rng.integers(2, 6))
        G = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        Gamma = G @ G.conj().T                       # >= 0 by construction
        Hh = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        H = 0.5 * (Hh + Hh.conj().T)
        delta = float(rng.normal(scale=3.0))
        e = np.zeros(n, dtype=complex)
        e[0] = 1.0

        for singular in (False, True):
            Gam = Gamma.copy()
            if singular:
                # project out one eigendirection to make Gamma singular
                w, V = np.linalg.eigh(Gam)
                k = int(np.argmin(w))
                Gam = Gam - w[k] * np.outer(V[:, k], V[:, k].conj())
            else:
                Gam = Gam + 1e-6 * np.eye(n)
            A = Gam + 1j * (H - delta * np.eye(n))
            if np.linalg.cond(A) > 1e12:
                continue
            x = np.linalg.solve(A, e)
            chi = np.vdot(e, x)                      # e^dag A^{-1} e
            quad = float(np.real(np.vdot(x, Gam @ x)))
            scale = max(abs(chi), 1e-12)
            max_identity_dev = max(max_identity_dev,
                                   abs(chi.real - quad) / scale)
            if singular:
                min_re_chi_singular = min(min_re_chi_singular, chi.real)
            else:
                min_re_chi = min(min_re_chi, chi.real)
    return {
        "n_trials": n_trials,
        "max_relative_deviation_from_identity": max_identity_dev,
        "min_Re_chi_strictly_positive_Gamma": min_re_chi,
        "min_Re_chi_singular_Gamma": min_re_chi_singular,
        "verdict": ("PASS" if (max_identity_dev < 1e-8 and min_re_chi >= 0.0
                               and min_re_chi_singular >= -1e-12) else "FAIL"),
        "statement": ("Re chi_matched = x^dag Gamma x with x = A^{-1} e; "
                      "nonnegative for every passive coherence block, and "
                      "zero only if x in ker Gamma."),
    }


# ----------------------------------------------------------------------
# Part 2: adversarial search in the full GKSL model
# ----------------------------------------------------------------------
N_E = 3


def unpack(v: np.ndarray):
    """Map an unconstrained vector to full-GKSL parameters."""
    i = 0

    def take(k):
        nonlocal i
        out = v[i:i + k]
        i += k
        return out

    Delta = take(N_E) * 3.0
    Omega_c = float(np.exp(take(1)[0]))                      # > 0
    dp = take(2 * N_E).reshape(N_E, 2) @ np.array([1.0, 1j])
    dc = take(2 * N_E).reshape(N_E, 2) @ np.array([1.0, 1j])
    Jm = take(3)
    Jmix = np.zeros((N_E, N_E))
    Jmix[0, 1] = Jmix[1, 0] = Jm[0]
    Jmix[1, 2] = Jmix[2, 1] = Jm[1]
    Jmix[0, 2] = Jmix[2, 0] = Jm[2]
    gamma_rad = (RATE_FLOOR + np.exp(take(2 * N_E))).reshape(N_E, 2).tolist()
    gamma_g = float(RATE_FLOOR + np.exp(take(1)[0]))
    gamma_dephase_e = (RATE_FLOOR + np.exp(take(N_E))).tolist()
    gamma_relax = float(RATE_FLOOR + np.exp(take(1)[0]))
    assert i == v.size, (i, v.size)
    # the probe leg must be coupled, else chi_matched is trivially zero
    if abs(dp[0]) < 1e-3:
        dp[0] = 1e-3
    return (Delta, Omega_c, dp, dc, Jmix, gamma_rad, gamma_g,
            gamma_dephase_e, gamma_relax)


NDIM = N_E + 1 + 2 * N_E + 2 * N_E + 3 + 2 * N_E + 1 + N_E + 1


def objective(v: np.ndarray) -> float:
    """min over the detuning scan of Re chi_matched (to be minimized)."""
    try:
        (Delta, Omega_c, dp, dc, Jmix, gamma_rad, gamma_g,
         gamma_dephase_e, gamma_relax) = unpack(v)
        res, _, eig0 = fg.run_point(Delta, Omega_c, dp, dc, Jmix, gamma_rad,
                                    gamma_g, gamma_dephase_e, DELTA_SCAN,
                                    gamma_relax=gamma_relax)
    except Exception:
        return 1e3
    if abs(eig0) > 1e-6:
        return 1e3                     # steady state not resolved; reject
    re = np.array([r["chi_matched"][0] for r in res])
    if not np.all(np.isfinite(re)):
        return 1e3
    # normalize by scale so the optimizer cannot "win" by shrinking chi
    scale = max(float(np.max(np.abs(re))), 1e-12)
    return float(re.min() / scale)


def part2_adversarial(n_starts: int = 24, maxiter: int = 600,
                      seed: int = 20260731) -> dict:
    rng = np.random.default_rng(seed)
    best = {"value": np.inf, "x": None}
    trace = []
    for k in range(n_starts):
        x0 = rng.normal(scale=0.8, size=NDIM)
        r = minimize(objective, x0, method="Nelder-Mead",
                     options={"maxiter": maxiter, "xatol": 1e-4,
                              "fatol": 1e-8, "adaptive": True})
        trace.append({"start": k, "f0": float(objective(x0)),
                      "f_final": float(r.fun), "nit": int(r.nit)})
        if r.fun < best["value"]:
            best = {"value": float(r.fun), "x": r.x.tolist()}
    return {
        "n_starts": n_starts,
        "normalized_objective_definition":
            "min_delta Re chi_matched / max_delta |Re chi_matched|",
        "adversarial_infimum": best["value"],
        "random_scan_baseline_F5B_min_Re_chi": 0.0029,
        "gain_found": bool(best["value"] < -1e-8),
        "verdict": "PASS" if best["value"] >= -1e-8 else "FAIL",
        "per_start": trace,
        "best_x": best["x"],
    }


def main() -> None:
    report = {"gate": "F5-C",
              "part1_coherence_block_identity": part1_identity(),
              "part2_full_gksl_adversarial": part2_adversarial()}
    p1, p2 = report["part1_coherence_block_identity"], report["part2_full_gksl_adversarial"]
    print("Part 1 (identity Re chi = x^dag Gamma x):")
    print(f"  max relative deviation : {p1['max_relative_deviation_from_identity']:.3e}")
    print(f"  min Re chi (Gamma > 0) : {p1['min_Re_chi_strictly_positive_Gamma']:.3e}")
    print(f"  min Re chi (Gamma sing): {p1['min_Re_chi_singular_Gamma']:.3e}")
    print(f"  verdict                : {p1['verdict']}")
    print("Part 2 (adversarial search, full GKSL):")
    print(f"  adversarial infimum of normalized min Re chi : "
          f"{p2['adversarial_infimum']:.6e}")
    print(f"  gain-like counterexample found              : {p2['gain_found']}")
    print(f"  verdict                                     : {p2['verdict']}")
    out = RESULTS_DIR / "gate_F5C_matched_floor_accretive.json"
    out.write_text(json.dumps(report, indent=2))
    print("Report written to", out)


if __name__ == "__main__":
    main()
