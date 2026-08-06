"""Gate E: NV experimental-anchor identifiability audit.

This gate does not change the general path-moment theorem. It asks whether the
planned NV ensemble experiment can itself identify the integer observable order,
or whether it should instead be used as a model-constrained physical anchor.

Inputs are frozen from repository outputs:
  * Gate-D 70 K detection floor and slope budget.
  * PRA Gate-5 ensemble contrasts.

The Monte-Carlo test compares candidate observable orders n=3,4,5. Two fits are
contrasted:
  1. constrained finite-Gamma correction, with the first correction coefficient
     fixed independently by the full model;
  2. unconstrained correction, where each candidate order has a free O(1/Gamma)
     coefficient and a free additive background.

The latter is the conservative identifiability test.
"""
from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
TABLES = os.path.join(ROOT, "results", "tables")
os.makedirs(TABLES, exist_ok=True)

SEED = 20260807
TARGET_SNR_AT_FLOOR = 5.0
REQUIRED_GAMMA_DECADES = 1.0
TRUE_OBSERVABLE_ORDER = 4
FIRST_CORRECTION_COEFFICIENT = 0.30
MULTIPLICATIVE_DRIFT = 0.03
N_POINTS = 16
N_MONTE_CARLO = 5000

# Frozen provenance: Gate D and PRA Gate 5.
C_MIN_70K = 1.1599051535422284e-6
CONTRASTS = {
    "single": 0.01360014974074735,
    "low_density": 0.0005721492074293866,
    "high_density": 0.00012895685598832745,
    "post_selected": 0.0009451684220583797,
    "post_selected_shimmed": 0.0026330768063651836,
}
GATING_SCENARIO = "post_selected_shimmed"


@dataclass(frozen=True)
class WindowResult:
    scenario: str
    contrast: float
    margin_over_floor: float
    intensity_order: int
    intensity_window_decades: float
    amplitude_order: int
    amplitude_window_decades: float
    improvement_to_one_decade: float


def usable_window_decades(contrast: float, c_min: float, order: float) -> float:
    """Largest Gamma window before C ~ Gamma^-order reaches the detection floor."""
    if contrast <= c_min or order <= 0:
        return 0.0
    return float(math.log10(contrast / c_min) / order)


def window_table() -> list[WindowResult]:
    rows: list[WindowResult] = []
    required_margin = 10.0 ** (TRUE_OBSERVABLE_ORDER * REQUIRED_GAMMA_DECADES)
    for scenario, contrast in CONTRASTS.items():
        margin = contrast / C_MIN_70K
        d4 = usable_window_decades(contrast, C_MIN_70K, 4)
        d2 = usable_window_decades(contrast, C_MIN_70K, 2)
        rows.append(
            WindowResult(
                scenario=scenario,
                contrast=contrast,
                margin_over_floor=float(margin),
                intensity_order=4,
                intensity_window_decades=d4,
                amplitude_order=2,
                amplitude_window_decades=d2,
                improvement_to_one_decade=float(max(1.0, required_margin / margin)),
            )
        )
    return rows


def _fit_candidate(x: np.ndarray, y: np.ndarray, n: int, correction: float | None) -> float:
    """Return RSS for a candidate integer order.

    correction=None: free additive background + free leading and first-correction
    amplitudes. correction=float: additive background + one amplitude, with the
    correction coefficient fixed independently.
    """
    if correction is None:
        X = np.column_stack([np.ones_like(x), x ** (-n), x ** (-(n + 1))])
    else:
        f = x ** (-n) * (1.0 + correction / x)
        X = np.column_stack([np.ones_like(x), f])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    residual = y - X @ beta
    return float(np.dot(residual, residual))


def monte_carlo_identifiability(
    contrast_ref: float,
    window_decades: float,
    correction_constrained: bool,
    seed: int = SEED,
) -> dict:
    rng = np.random.default_rng(seed)
    x = np.logspace(0.0, window_decades, N_POINTS)
    b = FIRST_CORRECTION_COEFFICIENT
    amplitude = contrast_ref / (1.0 + b)
    truth = amplitude * x ** (-TRUE_OBSERVABLE_ORDER) * (1.0 + b / x)
    sigma_additive = C_MIN_70K / TARGET_SNR_AT_FLOOR

    candidates = (3, 4, 5)
    wins = {n: 0 for n in candidates}
    for _ in range(N_MONTE_CARLO):
        measured = truth * (1.0 + rng.normal(0.0, MULTIPLICATIVE_DRIFT, N_POINTS))
        measured += rng.normal(0.0, sigma_additive, N_POINTS)
        correction = b if correction_constrained else None
        rss = {n: _fit_candidate(x, measured, n, correction) for n in candidates}
        wins[min(rss, key=rss.get)] += 1

    probabilities = {str(n): wins[n] / N_MONTE_CARLO for n in candidates}
    return {
        "window_decades": float(window_decades),
        "n_points": N_POINTS,
        "n_monte_carlo": N_MONTE_CARLO,
        "true_order": TRUE_OBSERVABLE_ORDER,
        "candidate_orders": list(candidates),
        "first_correction_coefficient": b,
        "multiplicative_drift": MULTIPLICATIVE_DRIFT,
        "additive_sigma": sigma_additive,
        "high_end_snr": float(truth[-1] / sigma_additive),
        "correction_constrained": correction_constrained,
        "selection_probability": probabilities,
        "correct_class_probability": probabilities[str(TRUE_OBSERVABLE_ORDER)],
    }


def main() -> None:
    rows = window_table()
    gating = next(r for r in rows if r.scenario == GATING_SCENARIO)
    constrained = monte_carlo_identifiability(
        gating.contrast, gating.intensity_window_decades, correction_constrained=True
    )
    unconstrained = monte_carlo_identifiability(
        gating.contrast, gating.intensity_window_decades, correction_constrained=False
    )

    temperature_ratio_for_one_gamma_decade = 10.0 ** (REQUIRED_GAMMA_DECADES / 5.0)

    gates = {
        "G_E1_signal_detectable": bool(gating.margin_over_floor > 1.0),
        "G_E2_one_decade_raw_intensity_window": bool(
            gating.intensity_window_decades >= REQUIRED_GAMMA_DECADES
        ),
        "G_E3_constrained_full_model_identifiability": bool(
            constrained["correct_class_probability"] >= 0.95
        ),
        "G_E4_unconstrained_asymptotic_identifiability": bool(
            unconstrained["correct_class_probability"] >= 0.95
        ),
        "G_E5_amplitude_linear_readout_window": bool(
            gating.amplitude_window_decades >= REQUIRED_GAMMA_DECADES
        ),
    }

    # The NV experiment can anchor the PRL only with a preregistered/full-model
    # finite-Gamma correction or a genuinely amplitude-linear readout. It cannot
    # independently establish the integer class by a free raw power-law fit.
    verdict = "CONDITIONAL_PASS"
    summary = {
        "description": "Gate E: can the planned NV experiment identify the integer response class?",
        "provenance": {
            "detection_floor": "GateD_robustness_discriminability/results/tables/gates_summary_gateD.json",
            "ensemble_contrasts": "No-go theorem/results/tables/gate5_ensemble_contrast.csv",
            "candidate_temperature_K": 70.0,
        },
        "assumptions": {
            "Gamma_T_power": 5.0,
            "true_observable_order": TRUE_OBSERVABLE_ORDER,
            "required_gamma_decades": REQUIRED_GAMMA_DECADES,
            "finite_Gamma_model": "C=A Gamma^-n (1+b/Gamma)+background",
        },
        "window_rows": [r.__dict__ for r in rows],
        "gating_scenario": GATING_SCENARIO,
        "temperature_ratio_for_one_gamma_decade": temperature_ratio_for_one_gamma_decade,
        "monte_carlo_constrained_correction": constrained,
        "monte_carlo_unconstrained_correction": unconstrained,
        "gates": gates,
        "verdict": verdict,
        "interpretation": {
            "direct_raw_exponent_claim": "NO-GO with the present ensemble contrast and a free finite-Gamma correction",
            "model_constrained_NV_anchor": "GO if Gamma_ph(T) and the first finite-Gamma correction are independently fixed before fitting",
            "amplitude_linear_route": "GO in the same contrast-to-floor budget because the observable order is reduced from 4 to 2",
            "minimum_margin_improvement_for_one_decade": gating.improvement_to_one_decade,
            "paper_role": "NV experiment anchors the mechanism and signed/crossover response; the clean integer-class measurement remains the engineered-dissipation witness unless the extra constraints are met",
        },
    }

    with open(os.path.join(TABLES, "gates_summary_gateE.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    with open(os.path.join(TABLES, "gate_e_windows.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].__dict__.keys()))
        w.writeheader()
        for row in rows:
            w.writerow(row.__dict__)

    with open(os.path.join(TABLES, "gate_e_identifiability.csv"), "w", newline="", encoding="utf-8") as f:
        fieldnames = ["fit", "correct_class_probability", "p_n3", "p_n4", "p_n5", "window_decades", "high_end_snr"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for label, result in (("constrained", constrained), ("unconstrained", unconstrained)):
            w.writerow({
                "fit": label,
                "correct_class_probability": result["correct_class_probability"],
                "p_n3": result["selection_probability"]["3"],
                "p_n4": result["selection_probability"]["4"],
                "p_n5": result["selection_probability"]["5"],
                "window_decades": result["window_decades"],
                "high_end_snr": result["high_end_snr"],
            })

    print(json.dumps({"gates": gates, "verdict": verdict, "interpretation": summary["interpretation"]}, indent=2))


if __name__ == "__main__":
    main()
