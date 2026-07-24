"""Phase Z runner: physical detuning signature z*Z (gates Z0-Z6).

Certifies that the Phase N path-order results survive when the common
resolvent shift -izI is replaced by a physical diagonal detuning
signature -izZ derived from swept laser frequencies.

Usage:
    python run_phase_z.py [--smoke]
"""

import argparse
import json
import math
import os
import time
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import phase_n_exact_core as base_core
import phase_n_frequency_core as freq_core
import phase_z_detuning_core as zc
from phase_n_exact_core import (
    D_RATES,
    GaussianRational,
    J45_STAR,
    S34,
    ZERO,
    as_gq,
    edges,
    linear_fit,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "results"))
FIGDIR = os.path.join(OUT, "figures")

Q_MAX = Fraction(6)


# ---------------------------------------------------------------------------
# serialization helpers
# ---------------------------------------------------------------------------

def clean(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, GaussianRational):
        return value.text()
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, float):
        return value
    return value


def poly3_equal(left, right):
    return left == right


# ---------------------------------------------------------------------------
# expected generic fan (Priority 1-3 certified): 4-q, 2+q, 4
# ---------------------------------------------------------------------------

def expected_generic_fan():
    return [
        {"q_lo": Fraction(0), "q_hi": Fraction(1),
         "slope": Fraction(-1), "intercept": Fraction(4)},
        {"q_lo": Fraction(1), "q_hi": Fraction(2),
         "slope": Fraction(1), "intercept": Fraction(2)},
        {"q_lo": Fraction(2), "q_hi": Q_MAX,
         "slope": Fraction(0), "intercept": Fraction(4)},
    ]


# ---------------------------------------------------------------------------
# Gate Z0: exact regression against Priority 3 and Phase N polynomials
# ---------------------------------------------------------------------------

def gate_z0():
    n3, q3 = zc.detuned_master_polynomials(zc.Z_IDENTITY)
    n3_ref, q3_ref = freq_core.frequency_master_polynomials()
    n0, q0 = base_core.master_polynomials()
    moments = zc.ideal_moment_polynomials(zc.Z_IDENTITY)
    moments_ref = freq_core.ideal_moment_polynomials()
    checks = {
        "numerator3_matches_priority3": poly3_equal(n3, n3_ref),
        "denominator3_matches_priority3": poly3_equal(q3, q3_ref),
        "numerator_z0_matches_phaseN": zc.substitute_z(n3, 0) == n0,
        "denominator_z0_matches_phaseN": zc.substitute_z(q3, 0) == q0,
        "ideal_moments_match_priority3": moments == moments_ref,
    }
    return {
        **checks,
        "overall_pass": all(checks.values()),
    }


# ---------------------------------------------------------------------------
# Gate Z1: exact path-order fan per detuning signature
# ---------------------------------------------------------------------------

GENERIC_Z0S = [Fraction(1, 3), Fraction(2, 7), Fraction(5, 11)]


def fan_for(zeta, z0):
    n3, q3 = zc.detuned_master_polynomials(zeta)
    return zc.exact_fan(
        zc.substitute_z(n3, z0), zc.substitute_z(q3, z0), q_max=Q_MAX
    ), n3, q3


def gate_z1(smoke=False):
    rows = {}
    all_pass = True
    decimal_checks = []
    generic = expected_generic_fan()
    for name, zeta in zc.Z_SIGNATURES.items():
        fan_a, n3, q3 = fan_for(zeta, GENERIC_Z0S[0])
        fan_b, _, _ = fan_for(zeta, GENERIC_Z0S[1])
        z0_independent = zc.fans_equal(fan_a, fan_b)
        if not z0_independent:
            fan_c, _, _ = fan_for(zeta, GENERIC_Z0S[2])
            if zc.fans_equal(fan_a, fan_c):
                z0_independent = True
            elif zc.fans_equal(fan_b, fan_c):
                fan_a = fan_b
                z0_independent = True
        preserved = zc.fans_equal(fan_a, generic)
        moments = zc.ideal_moment_polynomials(zeta, kmax=8)
        k0 = zc.first_nonzero_moment_index(moments)
        plateau = zc.fan_order_at(fan_a, Q_MAX)
        plateau_matches = (k0 is not None and plateau == Fraction(k0))
        rows[name] = {
            "zeta": [str(Fraction(v)) for v in zeta],
            "fan": zc.fan_text(fan_a),
            "z0_independent": bool(z0_independent),
            "generic_fan_preserved": bool(preserved),
            "ideal_moment_first_nonzero_index": k0,
            "plateau_order_at_qmax": str(plateau),
            "plateau_matches_ideal_moment_index": bool(plateau_matches),
        }
        all_pass = all_pass and z0_independent and plateau_matches
        qs = [Fraction(1, 2)] if smoke else [
            Fraction(1, 2), Fraction(3, 2), Fraction(5, 2)
        ]
        for q in qs:
            exact_order = zc.fan_order_at(fan_a, q)
            _, measured = zc.point_logabs_and_order(
                n3, q3, 10.0, float(q), 1.0, float(GENERIC_Z0S[0]),
                dps=60,
            )
            error = abs(measured - float(exact_order))
            decimal_checks.append({
                "signature": name,
                "q": str(q),
                "exact_order": str(exact_order),
                "measured_order": measured,
                "abs_error": error,
            })
            all_pass = all_pass and error < 1e-3
    identity_ok = rows["probe_common"]["generic_fan_preserved"]
    all_pass = all_pass and identity_ok
    return {
        "signatures": rows,
        "decimal_cross_checks": decimal_checks,
        "identity_signature_reproduces_generic_fan": bool(identity_ok),
        "overall_pass": bool(all_pass),
    }


# ---------------------------------------------------------------------------
# Gate Z2: frequency-promotion roots per signature (exact)
# ---------------------------------------------------------------------------

WINDOW = (Fraction(-3), Fraction(3))


def analyze_promotion(zeta, kmax=9):
    moments = zc.ideal_moment_polynomials(zeta, kmax=kmax)
    k0 = zc.first_nonzero_moment_index(moments)
    result = {"first_nonzero_moment_index": k0}
    if k0 is None or k0 >= kmax:
        result["verdict"] = "undetermined"
        return result
    constraint = zc.real_root_constraint(moments[k0 - 1])
    if constraint is None or len(constraint) <= 1:
        result["verdict"] = "no_real_promotion_root"
        result["constraint_degree"] = 0 if constraint is None else max(
            len(constraint) - 1, 0
        )
        return result
    result["constraint_degree"] = len(constraint) - 1
    result["constraint_coefficients"] = [str(c) for c in constraint]
    total_roots = zc.count_real_roots(constraint, Fraction(-10 ** 6),
                                      Fraction(10 ** 6))
    in_window = zc.count_real_roots(constraint, WINDOW[0], WINDOW[1])
    result["real_roots_total"] = total_roots
    result["real_roots_in_window"] = in_window
    rational_root = zc.exact_rational_root_if_linear(constraint)
    roots = []
    if rational_root is not None:
        roots.append({
            "exact": str(rational_root),
            "float": float(rational_root),
            "in_window": bool(WINDOW[0] <= rational_root <= WINDOW[1]),
        })
    else:
        for lo, hi in zc.isolate_real_roots(constraint, Fraction(-10 ** 6),
                                            Fraction(10 ** 6)):
            mid = (lo + hi) / 2
            roots.append({
                "isolating_interval": [str(lo), str(hi)],
                "float": float(mid),
                "in_window": bool(WINDOW[0] <= mid <= WINDOW[1]),
            })
    result["roots"] = roots
    if total_roots > 0:
        result["promotion_exactly_one_order"] = zc.promotion_is_exactly_one_order(
            constraint, moments[k0]
        )
        result["verdict"] = "promotion_root_found"
    else:
        result["verdict"] = "no_real_promotion_root"
    return result


def gate_z2():
    rows = {}
    for name, zeta in zc.Z_SIGNATURES.items():
        rows[name] = analyze_promotion(zeta)
    identity = rows["probe_common"]
    identity_ok = (
        identity["first_nonzero_moment_index"] == 4
        and identity.get("roots")
        and identity["roots"][0].get("exact") == "543/280"
        and identity.get("promotion_exactly_one_order") is True
    )
    verdict_defined = all(
        row.get("verdict") in (
            "promotion_root_found", "no_real_promotion_root"
        )
        for row in rows.values()
    )
    representative = None
    for name in ("coupling23", "coupling24", "coupling35", "two_tone_23_24"):
        row = rows[name]
        if row.get("verdict") != "promotion_root_found":
            continue
        for root in row.get("roots", []):
            if root.get("in_window") and "exact" in root:
                representative = {
                    "signature": name,
                    "z_star_exact": root["exact"],
                    "z_star_float": root["float"],
                }
                break
        if representative:
            break
    if representative is None:
        representative = {
            "signature": "probe_common",
            "z_star_exact": "543/280",
            "z_star_float": float(Fraction(543, 280)),
            "note": "no physical signature with exact in-window root; "
                    "falling back to the probe_common signature",
        }
    return {
        "signatures": rows,
        "identity_regression_z_star_543_280": bool(identity_ok),
        "all_verdicts_defined": bool(verdict_defined),
        "representative": representative,
        "overall_pass": bool(identity_ok and verdict_defined),
    }


# ---------------------------------------------------------------------------
# Gate Z3: frequency unfolding at the representative promotion root
# ---------------------------------------------------------------------------

def find_crossover(n3, q3, q, z, target, log_min, log_max,
                   dps=80, grid=501):
    xs = [log_min + (log_max - log_min) * idx / (grid - 1)
          for idx in range(grid)]
    values = [
        zc.point_logabs_and_order(n3, q3, x, q, 1.0, z, dps=dps)[1] - target
        for x in xs
    ]
    brackets = []
    for left, right, fl, fr in zip(xs[:-1], xs[1:], values[:-1], values[1:]):
        if fl == 0 or fl * fr < 0:
            brackets.append((left, right))
    if not brackets:
        raise RuntimeError(f"no crossover bracket q={q} z={z}")
    left, right = brackets[-1]
    for _ in range(60):
        mid = 0.5 * (left + right)
        fl = zc.point_logabs_and_order(n3, q3, left, q, 1.0, z, dps=dps)[1] - target
        fm = zc.point_logabs_and_order(n3, q3, mid, q, 1.0, z, dps=dps)[1] - target
        if fl * fm <= 0:
            right = mid
        else:
            left = mid
    return 0.5 * (left + right)


def gate_z3(representative, smoke=False):
    zeta = zc.Z_SIGNATURES[representative["signature"]]
    z_star = Fraction(representative["z_star_exact"])
    fan_gen, n3, q3 = fan_for(zeta, GENERIC_Z0S[0])
    fan_star = zc.exact_fan(
        zc.substitute_z(n3, z_star), zc.substitute_z(q3, z_star),
        q_max=Q_MAX,
    )
    q_candidates = [Fraction(k, 4) for k in range(1, int(Q_MAX * 4))]
    promoted_qs = [
        q for q in q_candidates
        if zc.fan_order_at(fan_star, q) > zc.fan_order_at(fan_gen, q)
    ]
    if not promoted_qs:
        return {
            "fan_generic": zc.fan_text(fan_gen),
            "fan_at_z_star": zc.fan_text(fan_star),
            "overall_pass": False,
            "note": "no q with promoted order at z_star",
        }
    # Sample the unfolding inside the promoted window: from just above the
    # promotion onset up to half a unit past the saturation breakpoint
    # (deep-plateau q adds no new scaling law).
    deltas = {
        q: zc.fan_order_at(fan_star, q) - zc.fan_order_at(fan_gen, q)
        for q in promoted_qs
    }
    max_delta = max(deltas.values())
    q_sat = min(q for q in promoted_qs if deltas[q] == max_delta)
    q_start = min(promoted_qs) - Fraction(1, 4)
    ramp = [q for q in promoted_qs if q < q_sat]
    if smoke:
        qs = ([ramp[0]] if ramp else []) + [q_sat + Fraction(1, 2)]
        epsilons = [1e-5, 1e-6]
        grid = 161
        dps = 50
    else:
        qs = []
        if ramp:
            qs.append(ramp[0])
            if len(ramp) > 2:
                qs.append(ramp[len(ramp) // 2])
            qs.append(ramp[-1])
        qs.extend([q_sat, q_sat + Fraction(1, 2)])
        qs = sorted(set(q for q in qs if q <= Q_MAX))
        epsilons = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
        grid = 301
        dps = 80
    _ = q_start
    rows = []
    slope_pass = True
    collapse_spreads = []
    precision_checks = []
    for q in qs:
        nu_gen = zc.fan_order_at(fan_gen, q)
        nu_star = zc.fan_order_at(fan_star, q)
        delta_nu = nu_star - nu_gen
        predicted_slope = float(-1 / delta_nu)
        target = 0.5 * float(nu_gen + nu_star)
        log_max = min(
            8.0 + 1.2 * (-math.log10(min(epsilons))) / float(delta_nu), 38.0
        )
        crossovers = []
        for eps in epsilons:
            z_val = float(z_star) + eps
            crossovers.append(find_crossover(
                n3, q3, float(q), z_val, target, 0.4, log_max,
                dps=dps, grid=grid,
            ))
        n_fit = min(3, len(epsilons))
        slope, intercept, rms = linear_fit(
            [math.log10(eps) for eps in epsilons[-n_fit:]],
            crossovers[-n_fit:],
        )
        error = abs(slope - predicted_slope)
        tolerance = max(0.05, 0.02 * abs(predicted_slope))
        slope_pass = slope_pass and error < tolerance
        offsets = [-1.25 + 0.25 * k for k in range(12)]
        curves = []
        for eps, log_cross in zip(epsilons[-2:], crossovers[-2:]):
            z_val = float(z_star) + eps
            curves.append([
                zc.point_logabs_and_order(
                    n3, q3, log_cross + off, float(q), 1.0, z_val, dps=dps
                )[1]
                for off in offsets
            ])
        spread = max(
            max(curve[idx] for curve in curves)
            - min(curve[idx] for curve in curves)
            for idx in range(len(offsets))
        )
        collapse_spreads.append(spread)
        high = find_crossover(
            n3, q3, float(q), float(z_star) + epsilons[-1], target,
            0.4, log_max, dps=dps + 40, grid=grid,
        )
        precision_checks.append({
            "q": str(q),
            "log10_crossover_dps": crossovers[-1],
            "log10_crossover_dps_plus_40": high,
            "abs_difference": abs(crossovers[-1] - high),
        })
        rows.append({
            "q": str(q),
            "nu_generic": str(nu_gen),
            "nu_promoted": str(nu_star),
            "predicted_slope": predicted_slope,
            "fitted_slope": slope,
            "abs_error": error,
            "fit_rms": rms,
            "log10_crossovers": crossovers,
            "collapse_spread": spread,
        })
    collapse_pass = max(collapse_spreads) < 0.06
    precision_pass = all(
        row["abs_difference"] < 1e-6 for row in precision_checks
    )
    return {
        "representative": representative,
        "fan_generic": zc.fan_text(fan_gen),
        "fan_at_z_star": zc.fan_text(fan_star),
        "rows": rows,
        "precision_checks": precision_checks,
        "slope_pass": bool(slope_pass),
        "collapse_pass": bool(collapse_pass),
        "precision_pass": bool(precision_pass),
        "overall_pass": bool(slope_pass and collapse_pass and precision_pass),
    }


# ---------------------------------------------------------------------------
# Gate Z4: broadening / finite spectral resolution
# ---------------------------------------------------------------------------
# Fast float evaluator: same max-rescaled summation as the Decimal core, in
# float64 (adequate for the 1e-2 order tolerances of the window gates).

def poly_float_terms(poly):
    keys = list(poly)
    gs = np.array([g for g, _, _ in keys], dtype=float)
    ks = np.array([k for _, k, _ in keys], dtype=float)
    pzs = np.array([pz for _, _, pz in keys], dtype=float)
    coefs = np.array([
        complex(float(poly[key].re), float(poly[key].im)) for key in keys
    ])
    return gs, ks, pzs, coefs


def logabs_float(terms, log10_gamma, q, z):
    gs, ks, pzs, coefs = terms
    ln_gamma = log10_gamma * math.log(10.0)
    log_scale = (gs + q * ks) * ln_gamma
    if z == 0.0:
        mask = pzs == 0
        log_scale = log_scale[mask]
        coefs = coefs[mask]
        signs = np.ones(len(coefs))
    else:
        log_scale = log_scale + pzs * math.log(abs(z))
        signs = np.where((z < 0) & (pzs % 2 == 1), -1.0, 1.0)
    peak = log_scale.max()
    value = (coefs * signs * np.exp(log_scale - peak)).sum()
    return math.log(abs(value)) + peak


def window_log_m2(nterms, qterms, log10_gamma, q, z_values, log_weights):
    log_terms = [
        2.0 * (logabs_float(nterms, log10_gamma, q, z_val)
               - logabs_float(qterms, log10_gamma, q, z_val)) + log_w
        for z_val, log_w in zip(z_values, log_weights)
    ]
    return 0.5 * (zc.logsumexp(log_terms) - zc.logsumexp(log_weights))


def window_grid(z_star, sigma, points, shape):
    zs = np.linspace(z_star - 10 * sigma, z_star + 10 * sigma, points)
    eps = (zs - z_star) / sigma
    if shape == "gaussian":
        logw = -0.5 * eps ** 2
    else:
        logw = -np.log1p(eps ** 2)
    return list(map(float, zs)), list(map(float, logw))


def window_order(nterms, qterms, q, z_values, log_weights, log_gammas):
    values = [
        window_log_m2(nterms, qterms, lg, q, z_values, log_weights)
        for lg in log_gammas
    ]
    slope, _, _ = linear_fit(
        [lg * math.log(10.0) for lg in log_gammas], values
    )
    return -slope


def gate_z4(representative, z3_result, smoke=False):
    zeta = zc.Z_SIGNATURES[representative["signature"]]
    z_star = float(Fraction(representative["z_star_exact"]))
    fan_gen, n3, q3 = fan_for(zeta, GENERIC_Z0S[0])
    fan_star = zc.exact_fan(
        zc.substitute_z(n3, Fraction(representative["z_star_exact"])),
        zc.substitute_z(q3, Fraction(representative["z_star_exact"])),
        q_max=Q_MAX,
    )
    nterms = poly_float_terms(n3)
    qterms = poly_float_terms(q3)
    sigmas = [1e-2] if smoke else [1e-3, 1e-2, 1e-1]
    qs = [Fraction(1, 2), Fraction(5, 2)] if smoke else [
        Fraction(1, 2), Fraction(3, 2), Fraction(5, 2), Fraction(7, 2)
    ]
    points = 61 if smoke else 241
    log_gammas = [20.0 + k for k in range(4)] if smoke else [
        20.0 + k for k in range(7)
    ]
    rows = []
    order_pass = True
    for sigma in sigmas:
        zs, logw = window_grid(z_star, sigma, points, "gaussian")
        for q in qs:
            expected = float(zc.fan_order_at(fan_gen, q))
            measured = window_order(nterms, qterms, float(q), zs, logw,
                                    log_gammas)
            error = abs(measured - expected)
            order_pass = order_pass and error < 1e-2
            rows.append({
                "sigma": sigma,
                "q": str(q),
                "shape": "gaussian",
                "expected_generic_order": expected,
                "measured_window_order": measured,
                "abs_error": error,
            })
    zs, logw = window_grid(z_star, 1e-2, points, "lorentzian")
    lorentz_rows = []
    for q in ([qs[-1]] if smoke else qs):
        expected = float(zc.fan_order_at(fan_gen, q))
        measured = window_order(nterms, qterms, float(q), zs, logw,
                                log_gammas)
        error = abs(measured - expected)
        order_pass = order_pass and error < 1e-2
        lorentz_rows.append({
            "sigma": 1e-2,
            "q": str(q),
            "shape": "lorentzian",
            "expected_generic_order": expected,
            "measured_window_order": measured,
            "abs_error": error,
        })
    refine = []
    for pts in ([61, 121] if smoke else [121, 241, 481]):
        zs_r, logw_r = window_grid(z_star, 1e-2, pts, "gaussian")
        refine.append(window_order(
            nterms, qterms, 2.5, zs_r, logw_r, log_gammas
        ))
    refine_change = max(abs(a - b) for a, b in zip(refine[:-1], refine[1:]))
    refine_pass = refine_change < 1e-6

    # crossover Gamma_x(sigma) at a q with unit promotion gap
    q_unit = None
    for q in [Fraction(k, 4) for k in range(1, int(Q_MAX * 4))]:
        if zc.fan_order_at(fan_star, q) - zc.fan_order_at(fan_gen, q) == 1:
            q_unit = q
    crossover_rows = []
    slope_sigma = None
    if q_unit is not None and not smoke:
        nu_gen = float(zc.fan_order_at(fan_gen, q_unit))
        target = nu_gen + 0.5
        for sigma in sigmas:
            zs_c, logw_c = window_grid(z_star, sigma, 121, "gaussian")

            def local_order(lg):
                dlg = 0.125
                lo = window_log_m2(nterms, qterms, lg - dlg,
                                   float(q_unit), zs_c, logw_c)
                hi = window_log_m2(nterms, qterms, lg + dlg,
                                   float(q_unit), zs_c, logw_c)
                return -(hi - lo) / (2 * dlg * math.log(10.0))

            grid_lg = [0.5 + 0.25 * k for k in range(95)]
            vals = [local_order(lg) - target for lg in grid_lg]
            bracket = None
            for left, right, fl, fr in zip(grid_lg[:-1], grid_lg[1:],
                                           vals[:-1], vals[1:]):
                if fl == 0 or fl * fr < 0:
                    bracket = (left, right)
            if bracket is None:
                continue
            left, right = bracket
            for _ in range(30):
                mid = 0.5 * (left + right)
                if (local_order(left) - target) * (local_order(mid) - target) <= 0:
                    right = mid
                else:
                    left = mid
            crossover_rows.append({
                "sigma": sigma,
                "q": str(q_unit),
                "log10_Gamma_x": 0.5 * (left + right),
            })
        if len(crossover_rows) >= 2:
            slope_sigma, _, _ = linear_fit(
                [math.log10(row["sigma"]) for row in crossover_rows],
                [row["log10_Gamma_x"] for row in crossover_rows],
            )
    slope_pass = True if smoke or slope_sigma is None else (
        abs(slope_sigma + 1.0) < 0.2
    )
    return {
        "gaussian_rows": rows,
        "lorentzian_rows": lorentz_rows,
        "grid_refinement_orders": refine,
        "grid_refinement_max_change": refine_change,
        "sigma_crossover_rows": crossover_rows,
        "sigma_crossover_slope": slope_sigma,
        "expected_sigma_slope": -1.0,
        "order_pass": bool(order_pass),
        "refine_pass": bool(refine_pass),
        "sigma_slope_pass": bool(slope_pass),
        "overall_pass": bool(order_pass and refine_pass and slope_pass),
    }


# ---------------------------------------------------------------------------
# Gate Z5: full 25x25 GKSL Liouvillian consistency (exact + float)
# ---------------------------------------------------------------------------

N_LEVELS = 5


def gq_matrix(n):
    return [[ZERO for _ in range(n)] for _ in range(n)]


def gq_matmul(a, b):
    n = len(a)
    out = gq_matrix(n)
    for i in range(n):
        for k in range(n):
            if a[i][k].is_zero():
                continue
            for j in range(n):
                if b[k][j].is_zero():
                    continue
                out[i][j] = out[i][j] + a[i][k] * b[k][j]
    return out


def gq_dagger(a):
    n = len(a)
    return [[a[j][i].conjugate() for j in range(n)] for i in range(n)]


def exact_hamiltonian(zeta, z):
    h = gq_matrix(N_LEVELS)
    for row, col, coupling, _ in edges(J45_STAR):
        h[row + 1][col + 1] = coupling
        h[col + 1][row + 1] = coupling.conjugate()
    for reduced_idx in range(4):
        shift = -Fraction(z) * Fraction(zeta[reduced_idx])
        if shift != 0:
            h[reduced_idx + 1][reduced_idx + 1] = as_gq(shift)
    return h


def exact_dissipators(gamma, kappa):
    """List of (rate, jump matrix E) with E = |0><p| unit matrices."""
    out = []
    for physical, rate in zip(range(1, N_LEVELS), D_RATES):
        e = gq_matrix(N_LEVELS)
        e[0][physical] = as_gq(1)
        out.append((Fraction(gamma) * rate, e))
    for reduced_idx in S34:
        e = gq_matrix(N_LEVELS)
        e[0][reduced_idx + 1] = as_gq(1)
        out.append((Fraction(kappa), e))
    return out


def exact_liouvillian(h, dissipators):
    n2 = N_LEVELS * N_LEVELS
    columns = []
    minus_i = GaussianRational(Fraction(0), Fraction(-1))
    for index in range(n2):
        rho = gq_matrix(N_LEVELS)
        rho[index // N_LEVELS][index % N_LEVELS] = as_gq(1)
        comm = gq_matmul(h, rho)
        rhs = gq_matmul(rho, h)
        drho = [[minus_i * (comm[i][j] - rhs[i][j])
                 for j in range(N_LEVELS)] for i in range(N_LEVELS)]
        for rate, e in dissipators:
            ed = gq_dagger(e)
            gain = gq_matmul(gq_matmul(e, rho), ed)
            ede = gq_matmul(ed, e)
            loss1 = gq_matmul(ede, rho)
            loss2 = gq_matmul(rho, ede)
            half = as_gq(Fraction(1, 2))
            for i in range(N_LEVELS):
                for j in range(N_LEVELS):
                    drho[i][j] = drho[i][j] + as_gq(rate) * (
                        gain[i][j] - half * (loss1[i][j] + loss2[i][j])
                    )
        columns.append([drho[i][j] for i in range(N_LEVELS)
                        for j in range(N_LEVELS)])
    return [[columns[col][row] for col in range(n2)] for row in range(n2)]


def gq_solve(matrix, vector):
    n = len(matrix)
    a = [row[:] + [vector[idx]] for idx, row in enumerate(matrix)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if not a[row][col].is_zero():
                pivot = row
                break
        if pivot is None:
            raise RuntimeError("singular exact system")
        a[col], a[pivot] = a[pivot], a[col]
        inv = zc.gq_div(as_gq(1), a[col][col])
        a[col] = [entry * inv for entry in a[col]]
        for row in range(n):
            if row == col or a[row][col].is_zero():
                continue
            factor = a[row][col]
            a[row] = [
                entry - factor * ref
                for entry, ref in zip(a[row], a[col])
            ]
    return [a[row][n] for row in range(n)]


def exact_reduced_response(zeta, gamma, kappa, z):
    a = gq_matrix(4)
    for idx, rate in enumerate(D_RATES):
        diag = as_gq(Fraction(gamma) * rate / 2)
        if idx in S34:
            diag = diag + as_gq(Fraction(kappa) / 2)
        shift = Fraction(z) * Fraction(zeta[idx])
        if shift != 0:
            diag = diag + GaussianRational(Fraction(0), -shift)
        a[idx][idx] = diag
    for row, col, coupling, _ in edges(J45_STAR):
        a[row][col] = a[row][col] + zc.I * coupling
        a[col][row] = a[col][row] + zc.I * coupling.conjugate()
    source = [GaussianRational(Fraction(0), Fraction(-1)),
              ZERO, ZERO, ZERO]
    return gq_solve(a, source)[3]


def exact_full_response(zeta, gamma, kappa, z):
    h = exact_hamiltonian(zeta, z)
    generator = exact_liouvillian(h, exact_dissipators(gamma, kappa))
    n2 = N_LEVELS * N_LEVELS
    # rhs = i (V rho0 - rho0 V) with V = |0><1| + |1><0|, rho0 = |0><0|:
    # (V rho0 - rho0 V) = |1><0| - |0><1|.
    rhs = [ZERO] * n2
    rhs[1 * N_LEVELS + 0] = GaussianRational(Fraction(0), Fraction(1))
    rhs[0 * N_LEVELS + 1] = GaussianRational(Fraction(0), Fraction(-1))
    for col in range(n2):
        generator[0][col] = as_gq(
            1 if col % (N_LEVELS + 1) == 0 else 0
        )
    rhs[0] = ZERO
    solution = gq_solve(generator, rhs)
    return solution[4 * N_LEVELS + 0]


def numeric_full_response(zeta, gamma, kappa, z):
    h = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    for row, col, coupling, _ in edges(J45_STAR):
        value = complex(float(coupling.re), float(coupling.im))
        h[row + 1, col + 1] = value
        h[col + 1, row + 1] = value.conjugate()
    for reduced_idx in range(4):
        h[reduced_idx + 1, reduced_idx + 1] = -z * float(zeta[reduced_idx])
    jumps = []
    for physical, rate in zip(range(1, N_LEVELS), D_RATES):
        op = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
        op[0, physical] = math.sqrt(gamma * float(rate))
        jumps.append(op)
    for reduced_idx in S34:
        op = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
        op[0, reduced_idx + 1] = math.sqrt(kappa)
        jumps.append(op)
    n2 = N_LEVELS * N_LEVELS
    generator = np.zeros((n2, n2), dtype=complex)
    for index in range(n2):
        rho = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
        rho.flat[index] = 1.0
        drho = -1j * (h @ rho - rho @ h)
        for op in jumps:
            opd = op.conj().T
            drho += op @ rho @ opd - 0.5 * (
                opd @ op @ rho + rho @ opd @ op
            )
        generator[:, index] = drho.reshape(-1)
    hermiticity = float(np.linalg.norm(h - h.conj().T))
    trace_row = np.zeros(n2, dtype=complex)
    trace_row[::N_LEVELS + 1] = 1.0
    trace_residual = float(np.linalg.norm(trace_row @ generator))
    rho0 = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    rho0[0, 0] = 1.0
    probe = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    probe[0, 1] = probe[1, 0] = 1.0
    rhs = (1j * (probe @ rho0 - rho0 @ probe)).reshape(-1)
    matrix = generator.copy()
    matrix[0, :] = trace_row
    rhs[0] = 0.0
    rho1 = np.linalg.solve(matrix, rhs).reshape(N_LEVELS, N_LEVELS)
    return rho1[4, 0], trace_residual, hermiticity


def numeric_reduced_response(zeta, gamma, kappa, z):
    a = np.zeros((4, 4), dtype=complex)
    for row, col, coupling, _ in edges(J45_STAR):
        value = complex(float(coupling.re), float(coupling.im))
        a[row, col] = 1j * value
        a[col, row] = 1j * value.conjugate()
    a += np.diag([gamma * float(rate) / 2 for rate in D_RATES])
    for idx in S34:
        a[idx, idx] += kappa / 2
    a -= 1j * z * np.diag([float(v) for v in zeta])
    source = np.array([-1j, 0j, 0j, 0j])
    return np.linalg.solve(a, source)[3]


def gate_z5(representative, smoke=False):
    exact_rows = []
    exact_all_zero = True
    exact_points = [(Fraction(3), Fraction(27), Fraction(1, 2))]
    if not smoke:
        exact_points.append((Fraction(10), Fraction(1000), Fraction(-2, 3)))
    signatures = (
        {"probe_common": zc.Z_SIGNATURES["probe_common"],
         representative["signature"]:
             zc.Z_SIGNATURES[representative["signature"]]}
        if smoke else zc.Z_SIGNATURES
    )
    for name, zeta in signatures.items():
        for gamma, kappa, z in exact_points:
            reduced = exact_reduced_response(zeta, gamma, kappa, z)
            full = exact_full_response(zeta, gamma, kappa, z)
            difference = full - reduced
            exact_all_zero = exact_all_zero and difference.is_zero()
            exact_rows.append({
                "signature": name,
                "Gamma": str(gamma),
                "kappa": str(kappa),
                "z": str(z),
                "reduced": reduced.text(),
                "full_liouvillian": full.text(),
                "difference": difference.text(),
                "exactly_zero": bool(difference.is_zero()),
            })
    float_rows = []
    max_error = 0.0
    max_trace = 0.0
    max_herm = 0.0
    z_rep = float(Fraction(representative["z_star_exact"]))
    z_values = [0.0, -1.0, z_rep, 2.5]
    cases = [(3.0, 3.0 ** 1.5)] if smoke else [
        (3.0, 3.0 ** 1.5), (10.0, 10.0 ** 2.5)
    ]
    for name, zeta in signatures.items():
        for gamma, kappa in cases:
            for z in z_values:
                full, trace_res, herm = numeric_full_response(
                    zeta, gamma, kappa, z
                )
                reduced = numeric_reduced_response(zeta, gamma, kappa, z)
                error = abs(full - reduced)
                max_error = max(max_error, error)
                max_trace = max(max_trace, trace_res)
                max_herm = max(max_herm, herm)
                float_rows.append({
                    "signature": name,
                    "Gamma": gamma,
                    "kappa": kappa,
                    "z": z,
                    "reduced_vs_full_error": error,
                })
    return {
        "exact_rows": exact_rows,
        "exact_all_differences_zero": bool(exact_all_zero),
        "float_rows": float_rows,
        "max_float_reduced_vs_full_error": max_error,
        "max_trace_preservation_residual": max_trace,
        "max_hamiltonian_hermiticity_residual": max_herm,
        "overall_pass": bool(
            exact_all_zero
            and max_error < 1e-11
            and max_trace < 1e-12
            and max_herm < 1e-12
        ),
    }


# ---------------------------------------------------------------------------
# Gate Z6: dimensional estimate (report-level numbers)
# ---------------------------------------------------------------------------

def gate_z6(representative):
    gamma0_mhz = 10.0  # unit rate gamma_0 = 2*pi x 10 MHz
    z_star = float(Fraction(representative["z_star_exact"]))
    zeta = zc.Z_SIGNATURES[representative["signature"]]
    _, n3, q3 = fan_for(zeta, GENERIC_Z0S[0])
    magnitudes = []
    for log10_gamma in [1.0, 1.5, 2.0]:
        logabs, order = zc.point_logabs_and_order(
            n3, q3, log10_gamma, 2.5, 1.0, z_star, dps=60
        )
        magnitudes.append({
            "Gamma_over_gamma0": 10.0 ** log10_gamma,
            "Gamma_physical_MHz_2pi": gamma0_mhz * 10.0 ** log10_gamma,
            "abs_R_S": math.exp(logabs),
            "local_effective_order": order,
        })
    return {
        "unit_rate_gamma0_MHz_2pi": gamma0_mhz,
        "representative_signature": representative["signature"],
        "z_star_model_units": z_star,
        "z_star_physical_MHz_2pi": z_star * gamma0_mhz,
        "resolution_requirement_note": (
            "Z4 shows any Gaussian broadening sigma >= 1e-3 model units "
            f"(2*pi x {1e-3 * gamma0_mhz * 1e3:.0f} kHz) restores the "
            "generic fan; resolving the pointwise promotion needs "
            "spectral selectivity below that scale."
        ),
        "R_S_magnitude_at_q_2p5": magnitudes,
        "overall_pass": True,
    }


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def make_fig_z1(z1):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    qs = np.linspace(0, float(Q_MAX), 481)
    for ax, (name, row) in zip(axes.flat, z1["signatures"].items()):
        orders = []
        for q in qs:
            for cell in row["fan"]:
                if Fraction(cell["q_lo"]) <= Fraction(str(round(q, 6))) \
                        <= Fraction(cell["q_hi"]):
                    orders.append(float(Fraction(cell["slope"])) * q
                                  + float(Fraction(cell["intercept"])))
                    break
        color = "tab:blue" if row["generic_fan_preserved"] else "tab:red"
        ax.plot(qs[:len(orders)], orders, color=color, lw=2)
        status = "generic fan" if row["generic_fan_preserved"] else "modified"
        ax.set_title(f"{name}\nZ=({','.join(row['zeta'])}) [{status}]",
                     fontsize=9)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("q  (kappa = Gamma^q)")
    for ax in axes[:, 0]:
        ax.set_ylabel("exact path order nu(q)")
    fig.suptitle("Phase Z: exact path-order fan per detuning signature Z")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figZ1_fan_per_signature.png"), dpi=180)
    plt.close(fig)


def make_fig_z2(representative):
    zeta_target = zc.Z_SIGNATURES[representative["signature"]]
    ts = np.linspace(0.0, 1.0, 41)
    branches = []
    for t in ts:
        tt = Fraction(str(round(float(t), 6)))
        zeta = tuple(
            (1 - tt) * Fraction(a) + tt * Fraction(b)
            for a, b in zip(zc.Z_IDENTITY, zeta_target)
        )
        moments = zc.ideal_moment_polynomials(zeta, kmax=8)
        k0 = zc.first_nonzero_moment_index(moments)
        roots_here = []
        if k0 is not None and k0 < 8:
            constraint = zc.real_root_constraint(moments[k0 - 1])
            if constraint and len(constraint) > 1:
                coeffs = [float(c) for c in reversed(constraint)]
                for root in np.roots(coeffs):
                    if abs(root.imag) < 1e-9 and -4 <= root.real <= 4:
                        roots_here.append(float(root.real))
        branches.append((float(t), roots_here, k0))
    fig, ax = plt.subplots(figsize=(7, 5))
    for t, roots_here, k0 in branches:
        for root in roots_here:
            ax.plot(t, root, "o", color="tab:blue", ms=3)
    ax.axhline(float(Fraction(543, 280)), color="gray", ls="--", lw=1,
               label="z* = 543/280 (Z = I)")
    ax.set_xlabel("t  in  Z(t) = (1-t) I + t Z_phys")
    ax.set_ylabel("real promotion roots z* in [-4, 4]")
    ax.set_title(
        "Phase Z: promotion-root locus toward "
        f"Z_phys = {representative['signature']}"
    )
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figZ2_root_locus.png"), dpi=180)
    plt.close(fig)


def make_fig_z3(z3):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for row in z3["rows"]:
        eps = [1e-2, 1e-3, 1e-4, 1e-5][:len(row["log10_crossovers"])]
        axes[0].plot(
            [math.log10(e) for e in eps], row["log10_crossovers"],
            "o-", label=(f"q={row['q']}  fit {row['fitted_slope']:.3f} "
                         f"(pred {row['predicted_slope']:.3f})"),
        )
    axes[0].set_xlabel("log10 |epsilon|")
    axes[0].set_ylabel("log10 Gamma_x")
    axes[0].set_title("frequency-unfolding crossover scale")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    axes[1].bar(
        [row["q"] for row in z3["rows"]],
        [row["collapse_spread"] for row in z3["rows"]],
        color="tab:green",
    )
    axes[1].axhline(0.06, color="red", ls="--", label="gate threshold")
    axes[1].set_xlabel("q")
    axes[1].set_ylabel("max collapse spread of nu_eff curves")
    axes[1].set_title("crossover collapse quality")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.suptitle(
        f"Phase Z gate Z3: unfolding at z* for {z3['representative']['signature']}"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figZ3_unfolding.png"), dpi=180)
    plt.close(fig)


def make_fig_z4(z4):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sigmas = sorted({row["sigma"] for row in z4["gaussian_rows"]})
    for sigma in sigmas:
        qs = [float(Fraction(row["q"])) for row in z4["gaussian_rows"]
              if row["sigma"] == sigma]
        measured = [row["measured_window_order"]
                    for row in z4["gaussian_rows"] if row["sigma"] == sigma]
        axes[0].plot(qs, measured, "o-", label=f"sigma={sigma:g}")
    expected = [(row["q"], row["expected_generic_order"])
                for row in z4["gaussian_rows"]
                if row["sigma"] == sigmas[0]]
    axes[0].plot([float(Fraction(q)) for q, _ in expected],
                 [v for _, v in expected], "k--", label="generic fan")
    axes[0].set_xlabel("q")
    axes[0].set_ylabel("window (L2) order")
    axes[0].set_title("broadened-window order vs generic fan")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    if z4["sigma_crossover_rows"]:
        axes[1].plot(
            [math.log10(row["sigma"]) for row in z4["sigma_crossover_rows"]],
            [row["log10_Gamma_x"] for row in z4["sigma_crossover_rows"]],
            "o-",
            label=(f"fit slope {z4['sigma_crossover_slope']:.3f} "
                   "(pred -1)"),
        )
        axes[1].legend()
    axes[1].set_xlabel("log10 sigma")
    axes[1].set_ylabel("log10 Gamma_x")
    axes[1].set_title("false-order window crossover vs broadening")
    axes[1].grid(alpha=0.3)
    fig.suptitle("Phase Z gate Z4: finite spectral resolution")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figZ4_broadening.png"), dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def run(smoke=False):
    start = time.time()
    os.makedirs(FIGDIR, exist_ok=True)
    gates = {}
    gates["Z0_exact_regression"] = gate_z0()
    gates["Z1_fan_per_signature"] = gate_z1(smoke=smoke)
    gates["Z2_promotion_roots"] = gate_z2()
    representative = gates["Z2_promotion_roots"]["representative"]
    gates["Z3_frequency_unfolding"] = gate_z3(representative, smoke=smoke)
    gates["Z4_broadening"] = gate_z4(
        representative, gates["Z3_frequency_unfolding"], smoke=smoke
    )
    gates["Z5_full_liouvillian"] = gate_z5(representative, smoke=smoke)
    gates["Z6_dimensional_estimate"] = gate_z6(representative)
    gates["phase_z_overall_pass"] = all(
        gate["overall_pass"]
        for gate in gates.values()
        if isinstance(gate, dict) and "overall_pass" in gate
    )
    gates["runtime_seconds"] = time.time() - start
    if not smoke:
        make_fig_z1(gates["Z1_fan_per_signature"])
        make_fig_z2(representative)
        make_fig_z3(gates["Z3_frequency_unfolding"])
        make_fig_z4(gates["Z4_broadening"])
    output = os.path.join(
        OUT,
        "gates_summary_phaseZ_smoke.json" if smoke
        else "gates_summary_phaseZ.json",
    )
    with open(output, "w") as handle:
        json.dump(clean(gates), handle, indent=2)
    print(json.dumps({
        "mode": "smoke" if smoke else "production",
        "overall_pass": gates["phase_z_overall_pass"],
        "runtime_seconds": round(gates["runtime_seconds"], 2),
        "representative": clean(representative),
        "output": output,
    }, indent=2))
    return gates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    run(smoke=args.smoke)
