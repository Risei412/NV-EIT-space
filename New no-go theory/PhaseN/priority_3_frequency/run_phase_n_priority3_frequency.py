"""Phase N Priority 3: finite-frequency and frequency-unfolding certification."""

from fractions import Fraction
import argparse
import json
import math
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-smrt-phase-n-p3")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from phase_n_exact_core import (
    D_RATES,
    J45_STAR,
    S34,
    effective_order,
    find_midpoint_crossover,
    linear_fit,
    master_polynomials,
    path_order,
    upper_breakpoints,
    GaussianRational,
)
from phase_n_frequency_core import (
    exact_linear_real_root,
    first_nonzero_moment_at,
    frequency_master_polynomials,
    ideal_moment_polynomials,
    logsumexp,
    point_logabs_and_order,
    point_scaled_complex,
    substitute_z,
)


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "phase_n_frequency_results")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)


def expected_generic(q):
    if q <= 1:
        return 4 - q
    if q <= 2:
        return 2 + q
    return 4


def expected_special(q):
    if q <= 1:
        return 4 - q
    if q <= 3:
        return 2 + q
    return 5


def fraction_rows(poly):
    return [
        {
            "power": power,
            "coefficient": coefficient.text(),
        }
        for power, coefficient in sorted(poly.items())
    ]


def gate_f0_f2_exact(smoke=False):
    numerator3, denominator3 = frequency_master_polynomials()
    old_numerator, old_denominator = master_polynomials()
    zero_numerator = substitute_z(numerator3, Fraction(0))
    zero_denominator = substitute_z(denominator3, Fraction(0))
    regression_pass = (
        zero_numerator == old_numerator
        and zero_denominator == old_denominator
    )

    moments = ideal_moment_polynomials(kmax=7)
    z_star = exact_linear_real_root(moments[3])
    root_pass = (
        z_star == Fraction(543, 280)
        and first_nonzero_moment_at(moments, Fraction(0)) == 4
        and first_nonzero_moment_at(moments, z_star) == 5
    )

    q_values = (
        [Fraction(0), Fraction(1), Fraction(2), Fraction(3), Fraction(7, 2)]
        if smoke else
        [Fraction(index, 8) for index in range(29)]
    )
    fan_rows = []
    fan_pass = True
    breakpoint_rows = []
    for label, z, expected, breakpoints_expected in [
        ("generic_z0", Fraction(0), expected_generic, [Fraction(1), Fraction(2)]),
        ("special_zstar", z_star, expected_special, [Fraction(1), Fraction(3)]),
    ]:
        numerator = substitute_z(numerator3, z)
        denominator = substitute_z(denominator3, z)
        breakpoints = sorted(set(
            upper_breakpoints(numerator) + upper_breakpoints(denominator)
        ))
        breakpoint_rows.append({
            "frequency": label,
            "z": str(z),
            "breakpoints": [str(value) for value in breakpoints],
            "expected": [str(value) for value in breakpoints_expected],
        })
        fan_pass &= breakpoints == breakpoints_expected
        for q in q_values:
            order, coefficient_n, coefficient_d = path_order(
                q, numerator, denominator
            )
            target = expected(q)
            fan_pass &= order == target
            fan_rows.append({
                "frequency": label,
                "z": str(z),
                "q": str(q),
                "order": str(order),
                "expected": str(target),
                "leading_numerator": coefficient_n.text(),
                "leading_denominator": coefficient_d.text(),
            })

    return {
        "F0_zero_frequency_exact_regression": {
            "numerator_exact_match": regression_pass,
            "denominator_exact_match": regression_pass,
            "overall_pass": bool(regression_pass),
        },
        "F1_frequency_selective_promotion": {
            "ideal_moment_polynomials": [
                fraction_rows(poly) for poly in moments[:6]
            ],
            "z_star": str(z_star),
            "z_star_float": float(z_star),
            "order_at_z0": first_nonzero_moment_at(moments, Fraction(0)),
            "order_at_z_star": first_nonzero_moment_at(moments, z_star),
            "overall_pass": bool(root_pass),
        },
        "F2_exact_generic_and_special_fans": {
            "breakpoint_rows": breakpoint_rows,
            "fan_rows": fan_rows,
            "overall_pass": bool(fan_pass),
        },
        "z_star": z_star,
        "numerator3": numerator3,
        "denominator3": denominator3,
    }


def frequency_crossover_gate(numerator3, denominator3, z_star, smoke=False):
    q_values = [2.5, 3.0, 3.5] if smoke else [2.25, 2.5, 2.75, 3.0, 3.5]
    epsilon_powers = [2, 3, 4] if smoke else [2, 3, 4, 5, 6]
    rows = []
    collapse_rows = []
    slope_pass = True
    precision_pass = True
    collapse_pass = True

    for q in q_values:
        delta_nu = q - 2 if q < 3 else 1
        xs = []
        ys = []
        crossovers = []
        for power in epsilon_powers:
            epsilon = Fraction(1, 10 ** power)
            z = z_star + epsilon
            numerator = substitute_z(numerator3, z)
            denominator = substitute_z(denominator3, z)
            predicted_log = power / delta_nu
            cross80 = find_midpoint_crossover(
                numerator,
                denominator,
                q,
                1.0,
                expected_special(q),
                expected_generic(q),
                predicted_log - 8,
                predicted_log + 8,
                dps=80,
            )
            cross120 = find_midpoint_crossover(
                numerator,
                denominator,
                q,
                1.0,
                expected_special(q),
                expected_generic(q),
                predicted_log - 8,
                predicted_log + 8,
                dps=120,
            )
            precision_error = abs(cross80 - cross120)
            precision_pass &= precision_error < 1e-8
            xs.append(-float(power))
            ys.append(cross120)
            crossovers.append((epsilon, numerator, denominator, cross120))
            rows.append({
                "q": q,
                "epsilon": f"1e-{power}",
                "z": f"{z.numerator}/{z.denominator}",
                "log10_Gamma_cross": cross120,
                "precision_abs_difference_log10": precision_error,
            })

        fit_count = 2 if smoke else 3
        slope, intercept, rms = linear_fit(xs[-fit_count:], ys[-fit_count:])
        expected_slope = -1 / delta_nu
        slope_error = abs(slope - expected_slope)
        tolerance = 0.10 if q == 2.25 else 0.05
        slope_pass &= slope_error < tolerance
        rows.append({
            "q": q,
            "fit_slope_dlogGamma_dlogepsilon": slope,
            "expected_slope": expected_slope,
            "slope_abs_error": slope_error,
            "fit_intercept": intercept,
            "fit_rms_log10": rms,
        })

        offsets = np.linspace(-1.0, 1.0, 31)
        curves = []
        selected = crossovers[-fit_count:]
        for _, numerator, denominator, crossover in selected:
            curve = [
                effective_order(
                    numerator,
                    denominator,
                    crossover + float(offset),
                    q,
                    1.0,
                    dps=100,
                )
                for offset in offsets
            ]
            curves.append(curve)
        spread = max(
            max(curve[index] for curve in curves)
            - min(curve[index] for curve in curves)
            for index in range(len(offsets))
        )
        collapse_pass &= spread < 0.05
        collapse_rows.append({"q": q, "max_order_spread": spread})

        if not smoke:
            for curve, (epsilon, _, _, _) in zip(curves, selected):
                plt.plot(
                    offsets,
                    curve,
                    linewidth=1.2,
                    label=rf"$\epsilon={float(epsilon):.0e}$",
                )
            plt.axhline(
                expected_special(q), color="0.5", linestyle="--", linewidth=0.8
            )
            plt.axhline(
                expected_generic(q), color="0.5", linestyle=":", linewidth=0.8
            )
            plt.xlabel(r"$\log_{10}(\Gamma/\Gamma_{\times,z})$")
            plt.ylabel(r"$\nu_{\rm eff}$")
            plt.title(f"Frequency crossover collapse, q={q:g}")
            plt.legend(fontsize=6)
            plt.tight_layout()
            plt.savefig(
                os.path.join(FIG, f"frequency_collapse_q_{q:g}.png"), dpi=180
            )
            plt.close()

    return {
        "crossover_rows": rows,
        "collapse_rows": collapse_rows,
        "slope_scaling_pass": bool(slope_pass),
        "precision_pass": bool(precision_pass),
        "collapse_pass": bool(collapse_pass),
        "overall_pass": bool(slope_pass and precision_pass and collapse_pass),
    }


def trapezoid_weights(size, width):
    weights = np.full(size, width / (size - 1), dtype=float)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    return weights


def window_metrics(
    numerator3,
    denominator3,
    log10_gamma,
    q,
    alpha,
    z_min,
    z_max,
    grid_size,
    dps=60,
):
    zs = np.linspace(z_min, z_max, grid_size)
    weights = trapezoid_weights(grid_size, z_max - z_min)
    logabs = []
    orders = []
    for z in zs:
        value, order = point_logabs_and_order(
            numerator3,
            denominator3,
            log10_gamma,
            q,
            alpha,
            float(z),
            dps=dps,
        )
        logabs.append(value)
        orders.append(order)

    max_index = int(np.argmax(logabs))
    result = {
        "M_inf": {
            "log_norm": logabs[max_index],
            "effective_order": orders[max_index],
            "argmax_z": float(zs[max_index]),
        }
    }
    for label, power in [("M_1", 1), ("M_2", 2)]:
        log_terms = [
            power * value + math.log(weight)
            for value, weight in zip(logabs, weights)
        ]
        log_integral = logsumexp(log_terms)
        normalized = [
            math.exp(value - log_integral) for value in log_terms
        ]
        effective_order = sum(
            weight * order for weight, order in zip(normalized, orders)
        )
        result[label] = {
            "log_norm": log_integral / power,
            "effective_order": effective_order,
        }
    return result


def fit_order(log10_gammas, log_norms):
    ln_gammas = [value * math.log(10) for value in log10_gammas]
    slope, intercept, rms = linear_fit(ln_gammas, log_norms)
    return -slope, intercept, rms


def finite_window_gate(numerator3, denominator3, smoke=False):
    windows = {
        "near": (-1.0, 1.0),
        "wide": (-3.0, 3.0),
    }
    q_values = [1.5, 2.5] if smoke else [0.5, 1.5, 2.5, 3.5]
    log_gammas = [8.0, 10.0, 12.0] if smoke else [10.0, 12.0, 14.0, 16.0]
    grid_size = 81 if smoke else 241
    rows = []
    all_pass = True
    for window, (z_min, z_max) in windows.items():
        for q in q_values:
            by_norm = {"M_inf": [], "M_1": [], "M_2": []}
            for log_gamma in log_gammas:
                metrics = window_metrics(
                    numerator3,
                    denominator3,
                    log_gamma,
                    q,
                    1.0,
                    z_min,
                    z_max,
                    grid_size,
                    dps=70,
                )
                for label in by_norm:
                    by_norm[label].append(metrics[label]["log_norm"])
            target = expected_generic(q)
            for label, log_norms in by_norm.items():
                order, _, rms = fit_order(log_gammas[-3:], log_norms[-3:])
                error = abs(order - target)
                all_pass &= error < 0.025
                rows.append({
                    "window": window,
                    "z_range": [z_min, z_max],
                    "q": q,
                    "norm": label,
                    "fitted_order": order,
                    "expected_generic_order": target,
                    "absolute_error": error,
                    "fit_rms_natural_log": rms,
                    "grid_size": grid_size,
                })
    return {
        "rows": rows,
        "interpretation": (
            "An isolated frequency-selective promotion does not promote "
            "finite-band L1, L2, or supremum norms."
        ),
        "overall_pass": bool(all_pass),
    }


def find_window_crossover(
    numerator3,
    denominator3,
    q,
    tuned_order,
    generic_order,
    predicted_log,
    z_min=-1.0,
    z_max=1.0,
    grid_size=51,
):
    target = 0.5 * (tuned_order + generic_order)

    def function(log_gamma):
        return (
            window_metrics(
                numerator3,
                denominator3,
                log_gamma,
                q,
                1.0,
                z_min,
                z_max,
                grid_size,
                dps=60,
            )["M_2"]["effective_order"]
            - target
        )

    grid = np.linspace(predicted_log - 7, predicted_log + 7, 71)
    values = [function(float(value)) for value in grid]
    brackets = []
    for left, right, fleft, fright in zip(
        grid[:-1], grid[1:], values[:-1], values[1:]
    ):
        if fleft == 0 or fleft * fright < 0:
            brackets.append((float(left), float(right)))
    if not brackets:
        raise RuntimeError("No finite-window crossover bracket")
    left, right = brackets[-1]
    for _ in range(45):
        mid = 0.5 * (left + right)
        if function(left) * function(mid) <= 0:
            right = mid
        else:
            left = mid
    return 0.5 * (left + right)


def finite_band_interference_unfolding(smoke=False):
    q = 1.5
    # The largest perturbation (1e-2) is intentionally retained only in
    # production diagnostics.  It is pre-asymptotic and must not define the
    # smoke-test scaling slope.
    powers = [3, 4] if smoke else [2, 3, 4, 5]
    rows = []
    crossovers = []
    for power in powers:
        delta = Fraction(1, 10 ** power)
        numerator3, denominator3 = frequency_master_polynomials(
            J45_STAR + GaussianRational(delta, 0)
        )
        crossover = find_window_crossover(
            numerator3,
            denominator3,
            q,
            tuned_order=3.5,
            generic_order=3.0,
            predicted_log=2 * power,
            grid_size=31 if smoke else 51,
        )
        rows.append({
            "q": q,
            "delta": f"1e-{power}",
            "log10_Gamma_cross_M2_near_window": crossover,
        })
        crossovers.append((delta, numerator3, denominator3, crossover))

    fit_count = 2 if smoke else 3
    xs = [-float(power) for power in powers[-fit_count:]]
    ys = [row["log10_Gamma_cross_M2_near_window"] for row in rows[-fit_count:]]
    slope, intercept, rms = linear_fit(xs, ys)
    slope_pass = abs(slope + 2.0) < 0.08

    offsets = np.linspace(-0.8, 0.8, 17)
    curves = []
    for _, numerator3, denominator3, crossover in crossovers[-fit_count:]:
        curves.append([
            window_metrics(
                numerator3,
                denominator3,
                crossover + float(offset),
                q,
                1.0,
                -1.0,
                1.0,
                51,
                dps=60,
            )["M_2"]["effective_order"]
            for offset in offsets
        ])
    spread = max(
        max(curve[index] for curve in curves)
        - min(curve[index] for curve in curves)
        for index in range(len(offsets))
    )
    collapse_pass = spread < 0.06
    return {
        "rows": rows,
        "fit_slope_dlogGamma_dlogdelta": slope,
        "expected_slope": -2.0,
        "fit_intercept": intercept,
        "fit_rms_log10": rms,
        "collapse_max_order_spread": spread,
        "slope_pass": bool(slope_pass),
        "collapse_pass": bool(collapse_pass),
        "overall_pass": bool(slope_pass and collapse_pass),
    }


def physical_frequency_gate(z_star, smoke=False):
    n_levels = 5

    def hamiltonian():
        matrix = np.zeros((n_levels, n_levels), dtype=complex)
        from phase_n_exact_core import edges
        for row, col, coupling, _ in edges(J45_STAR):
            value = complex(float(coupling.re), float(coupling.im))
            matrix[row + 1, col + 1] = value
            matrix[col + 1, row + 1] = value.conjugate()
        return matrix

    def jumps(gamma, kappa):
        operators = []
        for physical, rate in zip(range(1, 5), D_RATES):
            operator = np.zeros((n_levels, n_levels), dtype=complex)
            operator[0, physical] = np.sqrt(gamma * float(rate))
            operators.append(operator)
        for reduced_idx in S34:
            operator = np.zeros((n_levels, n_levels), dtype=complex)
            operator[0, reduced_idx + 1] = np.sqrt(kappa)
            operators.append(operator)
        return operators

    def liouvillian(h, operators):
        output = np.zeros((n_levels ** 2, n_levels ** 2), dtype=complex)
        for index in range(n_levels ** 2):
            rho = np.zeros((n_levels, n_levels), dtype=complex)
            rho.flat[index] = 1
            drho = -1j * (h @ rho - rho @ h)
            for operator in operators:
                product = operator.conj().T @ operator
                drho += (
                    operator @ rho @ operator.conj().T
                    - 0.5 * (product @ rho + rho @ product)
                )
            output[:, index] = drho.reshape(-1)
        return output

    def reduced_response(gamma, kappa, z):
        matrix = np.zeros((4, 4), dtype=complex)
        from phase_n_exact_core import edges
        for row, col, coupling, _ in edges(J45_STAR):
            value = complex(float(coupling.re), float(coupling.im))
            matrix[row, col] = 1j * value
            matrix[col, row] = 1j * value.conjugate()
        matrix += np.diag([gamma * float(rate) / 2 for rate in D_RATES])
        for index in S34:
            matrix[index, index] += kappa / 2
        matrix -= 1j * z * np.eye(4)
        source = np.array([-1j, 0j, 0j, 0j])
        return np.linalg.solve(matrix, source)[3]

    h = hamiltonian()
    rows = []
    maximum_error = 0.0
    maximum_trace = 0.0
    maximum_hermiticity = 0.0
    frequencies = [0.0, -1.0, float(z_star), 2.5]
    cases = [(3.0, 3.0 ** 1.5)] if smoke else [
        (3.0, 3.0 ** 1.5),
        (10.0, 10.0 ** 2.5),
    ]
    for gamma, kappa in cases:
        generator = liouvillian(h, jumps(gamma, kappa))
        trace = np.zeros(n_levels ** 2, dtype=complex)
        trace[::n_levels + 1] = 1
        trace_residual = float(np.linalg.norm(trace @ generator))
        hermiticity = float(np.linalg.norm(h - h.conj().T))
        rho0 = np.zeros((n_levels, n_levels), dtype=complex)
        rho0[0, 0] = 1
        probe = np.zeros((n_levels, n_levels), dtype=complex)
        probe[0, 1] = probe[1, 0] = 1
        rhs = 1j * (probe @ rho0 - rho0 @ probe)
        for z in frequencies:
            matrix = generator + 1j * z * np.eye(n_levels ** 2)
            vector = rhs.reshape(-1)
            if z == 0:
                matrix = matrix.copy()
                vector = vector.copy()
                matrix[0] = trace
                vector[0] = 0
            rho1 = np.linalg.solve(matrix, vector).reshape(n_levels, n_levels)
            reduced = reduced_response(gamma, kappa, z)
            error = abs(rho1[4, 0] - reduced)
            maximum_error = max(maximum_error, error)
            maximum_trace = max(maximum_trace, trace_residual)
            maximum_hermiticity = max(maximum_hermiticity, hermiticity)
            rows.append({
                "Gamma": gamma,
                "kappa": kappa,
                "z": z,
                "reduced_vs_full_error": error,
                "trace_preservation_residual": trace_residual,
                "hamiltonian_hermiticity_residual": hermiticity,
            })
    return {
        "rows": rows,
        "max_reduced_vs_full_error": maximum_error,
        "max_trace_preservation_residual": maximum_trace,
        "max_hamiltonian_hermiticity_residual": maximum_hermiticity,
        "overall_pass": bool(
            maximum_error < 1e-11
            and maximum_trace < 1e-12
            and maximum_hermiticity < 1e-12
        ),
    }


def no_real_pole_gate(z_star, smoke=False):
    # The Hermitian part of every reduced pencil is diagonal and bounded below
    # by Gamma*min(D_rates)/2 = Gamma/2 for Gamma>0. This excludes real-axis
    # zeros independently of z and kappa>=0.
    rows = []
    minimum_singular = float("inf")
    maximum_condition = 0.0
    gammas = [1.0, 10.0] if smoke else [1.0, 10.0, 100.0]
    qs = [1.5, 2.5]
    zs = np.linspace(-3, 3, 61 if smoke else 241)
    from phase_n_exact_core import edges
    b = np.zeros((4, 4), dtype=complex)
    for row, col, coupling, _ in edges(J45_STAR):
        value = complex(float(coupling.re), float(coupling.im))
        b[row, col] = 1j * value
        b[col, row] = 1j * value.conjugate()
    d = np.diag([float(rate) / 2 for rate in D_RATES])
    ds = np.zeros((4, 4))
    for index in S34:
        ds[index, index] = 0.5
    for gamma in gammas:
        for q in qs:
            kappa = gamma ** q
            local_min = float("inf")
            local_max_condition = 0.0
            for z in zs:
                for matrix in [
                    gamma * d + b - 1j * z * np.eye(4),
                    gamma * d + b + kappa * ds - 1j * z * np.eye(4),
                ]:
                    singular = np.linalg.svd(matrix, compute_uv=False)
                    local_min = min(local_min, float(singular[-1]))
                    local_max_condition = max(
                        local_max_condition, float(singular[0] / singular[-1])
                    )
            minimum_singular = min(minimum_singular, local_min)
            maximum_condition = max(maximum_condition, local_max_condition)
            rows.append({
                "Gamma": gamma,
                "q": q,
                "kappa": kappa,
                "minimum_singular_value": local_min,
                "maximum_condition_number": local_max_condition,
            })
    return {
        "analytic_accretivity_bound": "sigma_min(A)>=Gamma/2 for real z and kappa>=0",
        "frequency_window": [-3.0, 3.0],
        "rows": rows,
        "minimum_sampled_singular_value": minimum_singular,
        "maximum_sampled_condition_number": maximum_condition,
        "z_star_inside_window": -3 < float(z_star) < 3,
        "overall_pass": bool(minimum_singular > 0.49),
    }


def grid_refinement_gate(numerator3, denominator3, smoke=False):
    log_gammas = [10.0, 12.0, 14.0]
    sizes = [61, 121] if smoke else [121, 241, 481]
    rows = []
    orders = []
    for size in sizes:
        logs = []
        for log_gamma in log_gammas:
            metric = window_metrics(
                numerator3,
                denominator3,
                log_gamma,
                2.5,
                1.0,
                -3.0,
                3.0,
                size,
                dps=70,
            )
            logs.append(metric["M_2"]["log_norm"])
        order, _, rms = fit_order(log_gammas, logs)
        orders.append(order)
        rows.append({"grid_size": size, "M2_fitted_order": order, "fit_rms": rms})
    maximum_difference = max(orders) - min(orders)
    return {
        "rows": rows,
        "maximum_order_difference": maximum_difference,
        "overall_pass": bool(maximum_difference < 0.005),
    }


def make_figures(numerator3, denominator3, z_star):
    # Exact pointwise map, explicitly retaining the isolated z=z_star row.
    q_values = [Fraction(index, 40) for index in range(141)]
    z_values = [Fraction(-3) + Fraction(index, 20) for index in range(121)]
    closest = min(range(len(z_values)), key=lambda index: abs(z_values[index] - z_star))
    z_values[closest] = z_star
    z_values = sorted(set(z_values))
    values = np.empty((len(z_values), len(q_values)))
    for row, z in enumerate(z_values):
        numerator = substitute_z(numerator3, z)
        denominator = substitute_z(denominator3, z)
        for col, q in enumerate(q_values):
            values[row, col] = float(path_order(q, numerator, denominator)[0])
    plt.figure(figsize=(7.0, 4.6))
    image = plt.imshow(
        values,
        origin="lower",
        aspect="auto",
        extent=[0, 3.5, -3, 3],
        cmap="viridis",
        vmin=3,
        vmax=5,
        interpolation="nearest",
    )
    plt.colorbar(image, label=r"pointwise order $\nu(q,z)$")
    plt.axhline(float(z_star), color="white", linewidth=0.8, linestyle="--")
    plt.xlabel(r"path exponent $q$")
    plt.ylabel(r"detuning $z$")
    plt.title("Frequency-resolved path-order map")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "frequency_pointwise_order_map.png"), dpi=180)
    plt.close()

    # Effective-order heat map around the special frequency.
    q = 2.5
    log_gammas = np.linspace(2, 15, 100)
    offsets = np.linspace(-0.12, 0.12, 91)
    heat = np.empty((len(offsets), len(log_gammas)))
    for row, offset in enumerate(offsets):
        z = float(z_star) + float(offset)
        for col, log_gamma in enumerate(log_gammas):
            _, heat[row, col] = point_logabs_and_order(
                numerator3,
                denominator3,
                float(log_gamma),
                q,
                1.0,
                z,
                dps=60,
            )
    plt.figure(figsize=(7.0, 4.6))
    image = plt.imshow(
        heat,
        origin="lower",
        aspect="auto",
        extent=[log_gammas[0], log_gammas[-1], offsets[0], offsets[-1]],
        cmap="magma",
        vmin=4.0,
        vmax=4.5,
    )
    plt.colorbar(image, label=r"$\nu_{\rm eff}$")
    plt.xlabel(r"$\log_{10}\Gamma$")
    plt.ylabel(r"$z-z_\star$")
    plt.title(r"Frequency-selective false-high-order window ($q=5/2$)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "frequency_effective_order_heatmap.png"), dpi=180)
    plt.close()

    # Rescaled complex lineshape.
    zs = np.linspace(float(z_star) - 0.15, float(z_star) + 0.15, 301)
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.5), sharex=True)
    for log_gamma in [3.0, 5.0, 7.0]:
        curve = np.array([
            point_scaled_complex(
                numerator3,
                denominator3,
                log_gamma,
                2.5,
                1.0,
                float(z),
                rescale_order=4.0,
                dps=70,
            )
            for z in zs
        ])
        label = rf"$\Gamma=10^{{{int(log_gamma)}}}$"
        axes[0].plot(zs, curve.real, label=label)
        axes[1].plot(zs, curve.imag, label=label)
        axes[2].plot(zs, np.abs(curve), label=label)
    for axis, title in zip(axes, ["Real part", "Imaginary part", "Magnitude"]):
        axis.axvline(float(z_star), color="0.5", linestyle="--", linewidth=0.8)
        axis.set_title(title)
        axis.set_xlabel(r"$z$")
    axes[0].set_ylabel(r"$\Gamma^4 R_{S_{34}}$")
    axes[0].legend(fontsize=7)
    fig.suptitle(r"Rescaled frequency response near $z_\star$ ($q=5/2$)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "frequency_rescaled_lineshape.png"), dpi=180)
    plt.close(fig)


def run(smoke=False):
    exact = gate_f0_f2_exact(smoke=smoke)
    z_star = exact.pop("z_star")
    numerator3 = exact.pop("numerator3")
    denominator3 = exact.pop("denominator3")
    gates = dict(exact)
    gates["F3_F4_frequency_unfolding"] = frequency_crossover_gate(
        numerator3, denominator3, z_star, smoke=smoke
    )
    gates["F5_finite_band_interference_unfolding"] = (
        finite_band_interference_unfolding(smoke=smoke)
    )
    gates["F6_finite_window_norms"] = finite_window_gate(
        numerator3, denominator3, smoke=smoke
    )
    gates["F7_no_real_frequency_poles"] = no_real_pole_gate(
        z_star, smoke=smoke
    )
    gates["F8_full_liouvillian_frequency_consistency"] = (
        physical_frequency_gate(z_star, smoke=smoke)
    )
    gates["F9_frequency_grid_refinement"] = grid_refinement_gate(
        numerator3, denominator3, smoke=smoke
    )
    gates["priority3_overall_pass"] = all(
        gate["overall_pass"]
        for gate in gates.values()
        if isinstance(gate, dict) and "overall_pass" in gate
    )
    if not smoke:
        make_figures(numerator3, denominator3, z_star)
    output = os.path.join(
        OUT, "priority3_frequency_smoke.json"
        if smoke else "priority3_frequency_production.json"
    )
    with open(output, "w") as handle:
        json.dump(gates, handle, indent=2)
    print(json.dumps({
        "mode": "smoke" if smoke else "production",
        "overall_pass": gates["priority3_overall_pass"],
        "output": output,
    }))
    return gates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    run(smoke=args.smoke)
