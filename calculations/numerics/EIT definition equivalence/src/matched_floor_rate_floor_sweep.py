"""Gate F5-C follow-up: is the adversarial infimum physics, or the rate floor?

Why this exists
---------------
`results/certificates/gate_F5C_matched_floor_accretive.json` reports

    adversarial_infimum                 = 0.0031270095844016937
    random_scan_baseline_F5B_min_Re_chi = 0.0029

side by side, and Part 2's docstring reads the near-coincidence as "a directed
search that still cannot go below the random baseline's order of magnitude is
much stronger evidence than the random baseline alone".

Two things are wrong with reading it that way.

1. **They are not the same quantity.** The certificate's own
   `normalized_objective_definition` says the first is
   `min_delta Re chi_matched / max_delta |Re chi_matched|` — dimensionless, and
   deliberately so: `objective()` divides by `scale` precisely so the optimizer
   cannot win by shrinking chi. F5-B's 0.0029 is an unnormalized
   `min Re chi_matched`. Their agreeing to one significant figure is a
   coincidence of units, not a corroboration.

2. **The search was floored.** `unpack()` builds every rate as
   `RATE_FLOOR + exp(.)`, so no rate can go below `RATE_FLOOR = 1e-3` — the same
   order of magnitude as the reported infimum. A minimizer that stops at
   3.1e-3 while its parameters are pinned against a 1e-3 floor has not
   necessarily found a physical bound; it may have found its own constraint.

What this script does
---------------------
Reruns Part 2 at several values of `RATE_FLOOR`, holding `n_starts`,
`maxiter`, the RNG seed and `DELTA_SCAN` fixed, and records for each floor:

* the normalized infimum (the quantity Part 2 actually minimizes),
* the **unnormalized** `min_delta Re chi_matched` at the same point, so the two
  are never again confused,
* the response scale `max_delta |Re chi_matched|` that relates them,
* the smallest rate at the optimum, and how close it sits to the floor.

The last of these is the decisive one, and it is sharper than a log-log slope:
if the optimizer's rates are pinned at the floor, the constraint is binding and
the infimum is an artifact of it. If the optimum sits well above the floor and
the infimum stops moving, there is a scale-free bound to explain — one that
neither passivity nor the coherence-block identity `Re chi = x^dag Gamma x`
predicts, since both permit the normalized ratio to reach zero.

Reading the output
------------------
* infimum tracks RATE_FLOOR (slope ~1) and rates pinned  -> Part 2's headline
  conclusion is the search hitting its constraint. The comparison against
  F5-B's 0.0029 should be removed from the certificate.
* infimum saturates and rates unpinned                   -> a relative floor
  stronger than passivity, and a real question.

Neither outcome is assumed here. This script reports; it does not conclude.
"""
from __future__ import annotations

import json
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

import matched_floor_accretive_gate as gate
import full_gksl_2g3e as fg

RESULTS_DIR = Path(__file__).resolve().parents[4] / "results" / "certificates"
OUT_PATH = RESULTS_DIR / "gate_F5C_rate_floor_sweep.json"

#: Spanning three decades below and one above the committed 1e-3. 1e-3 is kept
#: so the sweep contains the published setting and the trend is internal to one
#: consistent (n_starts, maxiter, seed).
RATE_FLOORS = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]

N_STARTS = 12
MAXITER = 600
SEED = 20260731


def _rates_at(v: np.ndarray) -> dict:
    """Every dissipative rate at a parameter point, and its distance to the floor."""
    (_, _, _, _, _, gamma_rad, gamma_g, gamma_dephase_e, gamma_relax) = gate.unpack(v)
    rates = [r for leg in gamma_rad for r in leg]
    rates += [gamma_g, gamma_relax, *gamma_dephase_e]
    smallest = float(min(rates))
    return {
        "smallest_rate": smallest,
        "smallest_rate_over_floor": smallest / gate.RATE_FLOOR,
        "n_rates_within_1pct_of_floor": int(
            sum(1 for r in rates if r < gate.RATE_FLOOR * 1.01)
        ),
        "n_rates": len(rates),
    }


def _unnormalized_at(v: np.ndarray) -> dict:
    """The absolute min Re chi_matched and the scale the objective divides by.

    Recomputed rather than inferred: conflating these two is the defect this
    script exists to settle, so it must not repeat it.
    """
    (Delta, Omega_c, dp, dc, Jmix, gamma_rad, gamma_g,
     gamma_dephase_e, gamma_relax) = gate.unpack(v)
    res, _, eig0 = fg.run_point(
        Delta, Omega_c, dp, dc, Jmix, gamma_rad, gamma_g,
        gamma_dephase_e, gate.DELTA_SCAN, gamma_relax=gamma_relax,
    )
    re = np.array([r["chi_matched"][0] for r in res])
    scale = float(np.max(np.abs(re)))
    return {
        "min_Re_chi_matched_unnormalized": float(re.min()),
        "response_scale_max_abs_Re_chi": scale,
        "steady_state_eigenvalue": float(abs(eig0)),
    }


def _one_start(args) -> dict:
    """One Nelder-Mead descent at one rate floor. Runs in its own process."""
    rate_floor, k = args
    gate.RATE_FLOOR = rate_floor          # module-global; unpack() reads it
    rng = np.random.default_rng(SEED)
    for _ in range(k):                    # same start sequence as part2_adversarial
        rng.normal(scale=0.8, size=gate.NDIM)
    x0 = rng.normal(scale=0.8, size=gate.NDIM)
    t0 = time.time()
    r = minimize(
        gate.objective, x0, method="Nelder-Mead",
        options={"maxiter": MAXITER, "xatol": 1e-4, "fatol": 1e-8, "adaptive": True},
    )
    return {
        "rate_floor": rate_floor,
        "start": k,
        "f_final": float(r.fun),
        "nit": int(r.nit),
        "x": r.x.tolist(),
        "seconds": round(time.time() - t0, 2),
    }


def main() -> None:
    tasks = [(f, k) for f in RATE_FLOORS for k in range(N_STARTS)]
    started = time.time()
    print(f"{len(tasks)} descents ({len(RATE_FLOORS)} floors x {N_STARTS} starts)", flush=True)

    with ProcessPoolExecutor() as pool:
        results = list(pool.map(_one_start, tasks))

    per_floor = []
    for floor in RATE_FLOORS:
        rows = [r for r in results if r["rate_floor"] == floor]
        best = min(rows, key=lambda r: r["f_final"])
        gate.RATE_FLOOR = floor
        x = np.array(best["x"])
        entry = {
            "rate_floor": floor,
            "n_starts": len(rows),
            "adversarial_infimum_normalized": best["f_final"],
            "best_start": best["start"],
            "f_final_all": sorted(r["f_final"] for r in rows),
            "seconds_total": round(sum(r["seconds"] for r in rows), 1),
            **_unnormalized_at(x),
            **_rates_at(x),
        }
        per_floor.append(entry)
        print(
            f"  floor {floor:.0e}: normalized {entry['adversarial_infimum_normalized']:.6e}"
            f"  unnormalized {entry['min_Re_chi_matched_unnormalized']:.6e}"
            f"  smallest rate/floor {entry['smallest_rate_over_floor']:.2f}",
            flush=True,
        )

    floors = np.array([e["rate_floor"] for e in per_floor])
    infima = np.array([e["adversarial_infimum_normalized"] for e in per_floor])
    positive = infima > 0
    slope = None
    if positive.sum() >= 2:
        slope = float(
            np.polyfit(np.log10(floors[positive]), np.log10(infima[positive]), 1)[0]
        )

    report = {
        "gate": "F5-C rate-floor sweep",
        "question": (
            "Does the Part 2 adversarial infimum track RATE_FLOOR (a constraint "
            "artifact) or saturate (a scale-free bound stronger than passivity)?"
        ),
        "raised_by": "Research CI run_dceae5b2570e47a3, candidate matched-floor-scale-invariance",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config": {
            "rate_floors": RATE_FLOORS,
            "n_starts": N_STARTS,
            "maxiter": MAXITER,
            "seed": SEED,
            "delta_scan": [
                float(gate.DELTA_SCAN[0]), float(gate.DELTA_SCAN[-1]), int(gate.DELTA_SCAN.size)
            ],
        },
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "committed_reference": {
            "note": (
                "From gate_F5C_matched_floor_accretive.json. n_starts there was 24, "
                "so it is not directly comparable to a sweep row; the trend across "
                "floors below is internal to one fixed (n_starts, maxiter, seed)."
            ),
            "adversarial_infimum_normalized": 0.0031270095844016937,
            "rate_floor": 1e-3,
            "n_starts": 24,
            "random_scan_baseline_F5B_min_Re_chi_UNNORMALIZED": 0.0029,
            "why_the_two_must_not_be_compared": (
                "The first is min_delta Re chi / max_delta |Re chi| (dimensionless); "
                "the second is an unnormalized min Re chi. Agreement to one "
                "significant figure is a coincidence of units."
            ),
        },
        "per_floor": per_floor,
        "log_log_slope_infimum_vs_rate_floor": slope,
        "reading": (
            "slope near 1 with rates pinned at the floor => the infimum is the "
            "search terminating on its own constraint; slope near 0 with rates "
            "well above the floor => a floor-independent relative bound."
        ),
        "wall_clock_seconds": round(time.time() - started, 1),
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nlog-log slope (infimum vs RATE_FLOOR): {slope}")
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
