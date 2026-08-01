"""experimental_budget.py -- the measurement model behind Gate D's P8.

P8 asks whether the integer classes are experimentally DISCRIMINABLE. What has
to be measurable is the exponent, not the absolute response: two adjacent
classes differ by Delta nu = 1, so the question is whether a log-log fit over
the reachable Gamma window resolves the slope to better than 0.5 (a 2-sigma
separation of adjacent classes). Absolute contrast enters only as the
precondition that each point on that sweep clears its detection threshold.

Splitting that out of run_gate_d.py made three previously hidden assumptions
explicit:

  * the representative contrast was hard-coded to 1e-2, which is essentially
    the SINGLE-DEFECT value from the PRA campaign's ensemble study (0.0136).
    Every ensemble-averaged scenario in that same table is one to two orders
    of magnitude smaller, and the worst is level with the detection floor. The
    contrast is now read from that table, per scenario.
  * the slope precision came from fitting a NOISELESS analytic kernel, which
    unsurprisingly said one decade suffices. Noise is now injected and the
    fitted slope's spread is what decides.
  * the superconducting dynamic range was the literal 9.0, attributed to a
    sweep. Reach is now taken from the engineerable kappa range in
    model_sc_transfer, not from how far a numerical sweep happened to run.
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_NOGO_RESULTS = os.path.join(_HERE, "..", "..", "..", "No-go theorem", "results", "tables")
GATE5_CSV = os.path.join(_NOGO_RESULTS, "gate5_ensemble_contrast.csv")
GATE5_JSON = os.path.join(_NOGO_RESULTS, "gate5_summary.json")

# Adjacent classes differ by 1, so resolving one needs the fitted slope to be
# tighter than half that gap; |bias| + 2*std < 0.5 is a 2-sigma separation.
SLOPE_RESOLUTION = 0.5


def load_gate5_contrast(path=GATE5_CSV, summary_path=GATE5_JSON):
    """Ensemble-averaged EIT contrasts from the PRA campaign's Gate 5.

    Refuses a --quick table: the gating number must not come from a 10-draw
    ensemble average.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Generate it first with:\n"
            f"    python \"No-go theorem/src/gate5_ensemble_average.py\"\n"
            f"(full run, not --quick: Gate D gates on these numbers)."
        )
    if os.path.exists(summary_path):
        with open(summary_path) as f:
            meta = json.load(f)
        if meta.get("quick"):
            raise RuntimeError(
                f"{summary_path} was produced by a --quick run; its contrasts are "
                f"too coarse to gate on. Re-run gate5_ensemble_average.py in full."
            )
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            rows[r["scenario"]] = dict(
                Cmax=float(r["Cmax"]),
                washout_factor=float(r["washout_factor"]),
                fwhm_MHz=float(r["fwhm_MHz"]),
                n_draws=int(r["n_draws"]),
                n_orientations=int(r["n_orientations"]),
            )
    return rows


def slope_precision(gammas, values, sigma_rel, n_trials=400, seed=0):
    """Spread of the fitted log-log slope under multiplicative measurement noise.

    Returns the Monte-Carlo bias and standard deviation, plus the closed-form
    estimate std ~ sigma_rel / (sqrt(N) * std(ln Gamma)), which is the standard
    error of an ordinary-least-squares slope when the log-residuals are
    sigma_rel. Reporting both follows the repository convention of pairing an
    exact and a numerical estimator so a mistake in either is visible.
    """
    gammas = np.asarray(gammas, float)
    values = np.abs(np.asarray(values))
    rng = np.random.default_rng(seed)
    x = np.log(gammas)
    truth = -np.polyfit(x, np.log(values), 1)[0]

    noise = rng.normal(size=(n_trials, len(gammas))) * sigma_rel
    measured = values[None, :] * (1.0 + noise)
    measured = np.abs(measured) + 1e-300
    slopes = np.array([-np.polyfit(x, np.log(m), 1)[0] for m in measured])

    closed_form = sigma_rel / (np.sqrt(len(gammas)) * np.std(x))
    return dict(
        truth=float(truth),
        bias=float(np.mean(slopes) - truth),
        std=float(np.std(slopes)),
        closed_form_std=float(closed_form),
        resolves_adjacent_class=bool(
            abs(np.mean(slopes) - truth) + 2 * np.std(slopes) < SLOPE_RESOLUTION),
    )


def decades_needed(kernel_fn, gamma_lo, sigma_rel, widths=(1, 2, 3, 4, 5, 6),
                   points_per_decade=8, seed=0):
    """Smallest sweep width (in decades) whose noisy slope resolves Delta nu = 1.

    kernel_fn maps an array of Gamma to the response amplitude.
    """
    per = []
    needed = None
    for w in widths:
        gs = np.logspace(np.log10(gamma_lo), np.log10(gamma_lo) + w,
                         int(points_per_decade * w))
        prec = slope_precision(gs, kernel_fn(gs), sigma_rel, seed=seed)
        row = dict(window_decades=w, **prec)
        per.append(row)
        if needed is None and prec["resolves_adjacent_class"]:
            needed = w
    return dict(sigma_rel=sigma_rel, per_window=per, decades_needed=needed)


def sigma_rel_from_contrast(contrast, c_min, target_snr=5.0):
    """Relative precision of one contrast measurement.

    c_min is the contrast that reaches target_snr, so a measurement of size
    `contrast` has SNR = target_snr * contrast / c_min and relative precision
    1/SNR. Below c_min the point is not measurable at all.
    """
    if contrast <= 0 or not np.isfinite(c_min) or c_min <= 0:
        return np.inf
    snr = target_snr * contrast / c_min
    return float("inf") if snr <= 0 else float(1.0 / snr)


def usable_gamma_window(nu, contrast_ref, gamma_ref, c_min, target_snr=5.0,
                        gamma_max=None):
    """Gamma range over which the response still clears the detection floor.

    The response falls as (Gamma/Gamma_ref)^-nu from the anchor point, so the
    window ends where contrast reaches c_min. This is what makes a steep class
    hard to measure optically: at nu = 4 a single decade of Gamma costs four
    decades of signal.
    """
    if contrast_ref <= 0 or not np.isfinite(c_min) or c_min <= 0 or nu <= 0:
        return dict(gamma_lo=gamma_ref, gamma_hi=gamma_ref, decades=0.0)
    ratio = contrast_ref / c_min
    if ratio <= 1:
        return dict(gamma_lo=gamma_ref, gamma_hi=gamma_ref, decades=0.0)
    dec = float(np.log10(ratio) / nu)
    hi = gamma_ref * 10 ** dec
    if gamma_max is not None and hi > gamma_max:
        hi = float(gamma_max)
        dec = float(np.log10(hi / gamma_ref))
    return dict(gamma_lo=float(gamma_ref), gamma_hi=float(hi), decades=max(dec, 0.0))
