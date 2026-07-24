"""Phase N Priority 1: exact unfolding and finite-window crossover."""

from fractions import Fraction
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-smrt-phase-n-p1")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from phase_n_exact_core import (
    D_RATES,
    GaussianRational,
    J45_STAR,
    S34,
    edges,
    effective_order,
    exact_moments,
    find_midpoint_crossover,
    linear_fit,
    master_polynomials,
    path_order,
    suppression_index,
)


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "phase_n_priority_results")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)


def expected_tuned(q):
    if q <= 1:
        return 4 - q
    if q <= 2:
        return 2 + q
    return 4


def expected_generic(q):
    return 4 - q if q <= 1 else 3


def moment_gate():
    perturbations = {
        "tuned": GaussianRational(),
        "real_1e-3": GaussianRational(Fraction(1, 1000), 0),
        "imag_1e-3": GaussianRational(0, Fraction(1, 1000)),
    }
    rows = {}
    for name, delta in perturbations.items():
        moments = exact_moments(S34, J45_STAR + delta)
        rows[name] = {
            "delta": delta.text(),
            "moments": [value.text() for value in moments[:5]],
            "order": suppression_index(moments),
        }
    passed = (
        rows["tuned"]["order"] == 4
        and rows["real_1e-3"]["order"] == 3
        and rows["imag_1e-3"]["order"] == 3
    )
    return {"rows": rows, "overall_pass": passed}


def fan_gate():
    q_values = [Fraction(idx, 4) for idx in range(13)]
    perturbations = {
        "tuned": GaussianRational(),
        "real_1e-3": GaussianRational(Fraction(1, 1000), 0),
        "imag_1e-3": GaussianRational(0, Fraction(1, 1000)),
    }
    rows = {}
    passed = True
    for label, delta in perturbations.items():
        numerator, denominator = master_polynomials(J45_STAR + delta)
        data = []
        for q in q_values:
            order, coefficient_n, coefficient_d = path_order(
                q, numerator, denominator
            )
            expected = expected_tuned(q) if label == "tuned" else expected_generic(q)
            passed &= order == expected
            data.append({
                "q": str(q),
                "order": str(order),
                "expected": str(expected),
                "leading_numerator": coefficient_n.text(),
                "leading_denominator": coefficient_d.text(),
            })
        rows[label] = data
    return {"rows": rows, "overall_pass": bool(passed)}


def physical_liouvillian_gate():
    n_levels = 5

    def hamiltonian(j45, eliminated=()):
        removed = {idx + 1 for idx in eliminated}
        matrix = np.zeros((n_levels, n_levels), dtype=complex)
        for row, col, coupling, _ in edges(j45):
            physical_row, physical_col = row + 1, col + 1
            if physical_row in removed or physical_col in removed:
                continue
            value = complex(float(coupling.re), float(coupling.im))
            matrix[physical_row, physical_col] = value
            matrix[physical_col, physical_row] = value.conjugate()
        return matrix

    def jumps(gamma, eliminated=(), kappa=None):
        operators = []
        for physical, rate in zip(range(1, 5), D_RATES):
            operator = np.zeros((n_levels, n_levels), dtype=complex)
            operator[0, physical] = np.sqrt(gamma * float(rate))
            operators.append(operator)
        if kappa is not None:
            for reduced_idx in eliminated:
                operator = np.zeros((n_levels, n_levels), dtype=complex)
                operator[0, reduced_idx + 1] = np.sqrt(kappa)
                operators.append(operator)
        return operators

    def liouvillian(h, operators):
        output = np.zeros((n_levels ** 2, n_levels ** 2), dtype=complex)
        for idx in range(n_levels ** 2):
            rho = np.zeros((n_levels, n_levels), dtype=complex)
            rho.flat[idx] = 1
            drho = -1j * (h @ rho - rho @ h)
            for operator in operators:
                product = operator.conj().T @ operator
                drho += (
                    operator @ rho @ operator.conj().T
                    - 0.5 * (product @ rho + rho @ product)
                )
            output[:, idx] = drho.reshape(-1)
        return output

    def reduced_response(j45, gamma, kappa=None):
        b = np.zeros((4, 4), dtype=complex)
        for row, col, coupling, _ in edges(j45):
            value = complex(float(coupling.re), float(coupling.im))
            b[row, col] = 1j * value
            b[col, row] = 1j * value.conjugate()
        diagonal = np.diag([float(rate) / 2 for rate in D_RATES])
        cut = np.zeros((4, 4))
        for idx in S34:
            cut[idx, idx] = 0.5
        matrix = gamma * diagonal + b
        if kappa is not None:
            matrix = matrix + kappa * cut
        source = np.array([-1j, 0j, 0j, 0j])
        return np.linalg.solve(matrix, source)[3]

    def full_response(j45, gamma, kappa=None):
        h = hamiltonian(j45)
        generator = liouvillian(h, jumps(gamma, S34, kappa))
        rho0 = np.zeros((n_levels, n_levels), dtype=complex)
        rho0[0, 0] = 1
        probe = np.zeros((n_levels, n_levels), dtype=complex)
        probe[0, 1] = probe[1, 0] = 1
        rhs = 1j * (probe @ rho0 - rho0 @ probe)
        matrix = generator.copy()
        vector = rhs.reshape(-1)
        trace = np.zeros(n_levels ** 2, dtype=complex)
        trace[::n_levels + 1] = 1
        matrix[0] = trace
        vector[0] = 0
        rho1 = np.linalg.solve(matrix, vector).reshape(n_levels, n_levels)
        trace_residual = float(np.linalg.norm(trace @ generator))
        return rho1[4, 0], trace_residual, float(np.linalg.norm(h - h.conj().T))

    rows = []
    maximum_error = 0.0
    maximum_trace = 0.0
    maximum_hermiticity = 0.0
    perturbations = [
        GaussianRational(Fraction(1, 1000), 0),
        GaussianRational(0, Fraction(1, 1000)),
    ]
    for delta in perturbations:
        j45 = J45_STAR + delta
        for gamma, kappa in [(10.0, 50.0), (100.0, 500.0)]:
            reduced = reduced_response(j45, gamma, kappa)
            full, trace_residual, hermiticity = full_response(j45, gamma, kappa)
            error = abs(reduced - full)
            maximum_error = max(maximum_error, error)
            maximum_trace = max(maximum_trace, trace_residual)
            maximum_hermiticity = max(maximum_hermiticity, hermiticity)
            rows.append({
                "delta": delta.text(),
                "Gamma": gamma,
                "kappa": kappa,
                "reduced_vs_full_error": error,
                "trace_preservation_residual": trace_residual,
                "hamiltonian_hermiticity_residual": hermiticity,
            })
    passed = maximum_error < 1e-12 and maximum_trace < 1e-12 and maximum_hermiticity < 1e-12
    return {
        "rows": rows,
        "max_reduced_vs_full_error": maximum_error,
        "max_trace_preservation_residual": maximum_trace,
        "max_hamiltonian_hermiticity_residual": maximum_hermiticity,
        "overall_pass": bool(passed),
    }


def crossover_gate(smoke=False):
    q_values = [1.5, 2.0, 2.5] if smoke else [1.25, 1.5, 1.75, 2.0, 2.5, 3.0]
    delta_powers = [2, 3, 4] if smoke else [2, 3, 4, 5, 6]
    rows = []
    collapse_rows = []
    slopes_pass = True
    precision_pass = True
    collapse_pass = True
    for q in q_values:
        delta_nu = q - 1 if q < 2 else 1
        xs, ys = [], []
        q_crossovers = []
        for power in delta_powers:
            delta = Fraction(1, 10 ** power)
            numerator, denominator = master_polynomials(
                J45_STAR + GaussianRational(delta, 0)
            )
            predicted_log = power / delta_nu
            log_max = max(12.0, predicted_log + 8.0)
            cross80 = find_midpoint_crossover(
                numerator, denominator, q, 1.0,
                expected_tuned(q), expected_generic(q),
                predicted_log - 8.0, log_max, dps=80,
            )
            cross120 = find_midpoint_crossover(
                numerator, denominator, q, 1.0,
                expected_tuned(q), expected_generic(q),
                predicted_log - 8.0, log_max, dps=120,
            )
            precision_error = abs(cross80 - cross120)
            precision_pass &= precision_error < 1e-8
            xs.append(-float(power))
            ys.append(cross120)
            q_crossovers.append((delta, numerator, denominator, cross120))
            rows.append({
                "q": q,
                "delta": f"1e-{power}",
                "log10_Gamma_cross": cross120,
                "precision_abs_difference_log10": precision_error,
            })
        fit_count = 2 if smoke else 3
        slope, intercept, rms = linear_fit(xs[-fit_count:], ys[-fit_count:])
        expected_slope = -1 / delta_nu
        slope_error = abs(slope - expected_slope)
        slopes_pass &= slope_error < (0.08 if q == 1.25 else 0.04)
        rows.append({
            "q": q,
            "fit_slope_dlogGamma_dlogdelta": slope,
            "expected_slope": expected_slope,
            "slope_abs_error": slope_error,
            "fit_rms_log10": rms,
        })

        scaled_offsets = np.linspace(-1.0, 1.0, 31)
        curves = []
        collapse_crossovers = q_crossovers[-fit_count:]
        for _, numerator, denominator, crossover in collapse_crossovers:
            curve = [
                effective_order(
                    numerator, denominator, crossover + float(offset),
                    q, 1.0, dps=100,
                )
                for offset in scaled_offsets
            ]
            curves.append(curve)
        spread = max(
            max(curve[idx] for curve in curves) - min(curve[idx] for curve in curves)
            for idx in range(len(scaled_offsets))
        )
        collapse_pass &= spread < 0.05
        collapse_rows.append({"q": q, "max_order_spread": spread})

        if not smoke:
            for curve, (delta, _, _, _) in zip(curves, collapse_crossovers):
                plt.plot(
                    scaled_offsets, curve, linewidth=1.3,
                    label=f"delta={float(delta):.0e}",
                )
            plt.axhline(expected_tuned(q), color="0.5", linestyle="--", linewidth=0.8)
            plt.axhline(expected_generic(q), color="0.5", linestyle=":", linewidth=0.8)
            plt.xlabel(r"$\log_{10}(\Gamma/\Gamma_\times)$")
            plt.ylabel(r"$\nu_{\rm eff}$")
            plt.title(f"Crossover collapse, q={q:g}")
            plt.legend(fontsize=6, ncol=2)
            plt.tight_layout()
            plt.savefig(os.path.join(FIG, f"p1_collapse_q_{q:g}.png"), dpi=180)
            plt.close()
    return {
        "crossover_rows": rows,
        "collapse_rows": collapse_rows,
        "slope_scaling_pass": bool(slopes_pass),
        "precision_pass": bool(precision_pass),
        "collapse_pass": bool(collapse_pass),
        "overall_pass": bool(slopes_pass and precision_pass and collapse_pass),
    }


def make_heatmap():
    q = 1.5
    delta_powers = np.linspace(2, 7, 70)
    log_gammas = np.linspace(2, 18, 100)
    values = np.empty((len(delta_powers), len(log_gammas)))
    for row, power in enumerate(delta_powers):
        delta = Fraction(1, 10 ** 7) if power == 7 else None
        if delta is None:
            delta_float = 10 ** (-float(power))
            delta = Fraction(delta_float).limit_denominator(10 ** 12)
        numerator, denominator = master_polynomials(
            J45_STAR + GaussianRational(delta, 0)
        )
        values[row] = [
            effective_order(numerator, denominator, point, q, 1.0, dps=80)
            for point in log_gammas
        ]
    plt.figure(figsize=(6.5, 4.5))
    image = plt.imshow(
        values, origin="lower", aspect="auto",
        extent=[log_gammas[0], log_gammas[-1], delta_powers[0], delta_powers[-1]],
        vmin=3.0, vmax=3.5, cmap="viridis",
    )
    plt.colorbar(image, label=r"$\nu_{\rm eff}$")
    plt.xlabel(r"$\log_{10}\Gamma$")
    plt.ylabel(r"$-\log_{10}|\delta|$")
    plt.title(r"Finite-window false high-order region ($q=3/2$)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "p1_effective_order_heatmap_q_1p5.png"), dpi=180)
    plt.close()


def run(smoke=False):
    gates = {
        "U1_U2_exact_moments": moment_gate(),
        "U3_exact_unfolded_fan": fan_gate(),
        "U4_U5_U6_crossover": crossover_gate(smoke=smoke),
        "U7_physical_liouvillian": physical_liouvillian_gate(),
    }
    gates["priority1_overall_pass"] = all(
        gate["overall_pass"] for gate in gates.values()
        if isinstance(gate, dict) and "overall_pass" in gate
    )
    if not smoke:
        make_heatmap()
    output = os.path.join(
        OUT, "priority1_smoke.json" if smoke else "priority1_production.json"
    )
    with open(output, "w") as handle:
        json.dump(gates, handle, indent=2)
    print(json.dumps({
        "mode": "smoke" if smoke else "production",
        "overall_pass": gates["priority1_overall_pass"],
        "output": output,
    }))
    return gates


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    run(smoke=args.smoke)
