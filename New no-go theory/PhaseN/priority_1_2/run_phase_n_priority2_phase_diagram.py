"""Phase N Priority 2: exact (q, alpha) phase diagram."""

from fractions import Fraction
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-smrt-phase-n-p2")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from phase_n_exact_core import (
    J45_STAR,
    effective_order,
    leading_face,
    master_polynomials,
    path_order,
    upper_breakpoints,
)


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "phase_n_priority_results")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)


def expected_order(q):
    if q <= 1:
        return 4 - q
    if q <= 2:
        return 2 + q
    return Fraction(4)


def face_text(face):
    return {str(power): coefficient.text() for power, coefficient in sorted(face.items())}


def common_positive_roots(face):
    """Exact positive-common-root count using a rational Sturm sequence."""

    def trim(poly):
        poly = list(poly)
        while poly and poly[-1] == 0:
            poly.pop()
        return poly

    def divmod_poly(left, right):
        left = trim(left)
        right = trim(right)
        quotient = [Fraction(0)] * max(1, len(left) - len(right) + 1)
        while left and len(left) >= len(right):
            degree = len(left) - len(right)
            factor = left[-1] / right[-1]
            quotient[degree] += factor
            for idx, value in enumerate(right):
                left[idx + degree] -= factor * value
            left = trim(left)
        return trim(quotient), left

    def monic(poly):
        poly = trim(poly)
        if not poly:
            return []
        lead = poly[-1]
        return [value / lead for value in poly]

    def gcd_poly(left, right):
        left, right = trim(left), trim(right)
        while right:
            _, remainder = divmod_poly(left, right)
            left, right = right, remainder
        return monic(left)

    def derivative(poly):
        return trim([Fraction(idx) * poly[idx] for idx in range(1, len(poly))])

    def sign(value):
        return 1 if value > 0 else -1 if value < 0 else 0

    def variations(signs):
        signs = [value for value in signs if value]
        return sum(left != right for left, right in zip(signs[:-1], signs[1:]))

    def positive_root_count(poly):
        poly = trim(poly)
        if len(poly) <= 1:
            return 0
        sequence = [poly, derivative(poly)]
        while sequence[-1]:
            _, remainder = divmod_poly(sequence[-2], sequence[-1])
            if not remainder:
                break
            sequence.append([-value for value in remainder])
        at_zero_plus = []
        at_infinity = []
        for item in sequence:
            first_nonzero = next((value for value in item if value), Fraction(0))
            at_zero_plus.append(sign(first_nonzero))
            at_infinity.append(sign(item[-1]))
        return variations(at_zero_plus) - variations(at_infinity)

    max_power = max(face) if face else 0
    real = [Fraction(0)] * (max_power + 1)
    imag = [Fraction(0)] * (max_power + 1)
    for power, coefficient in face.items():
        real[power] = coefficient.re
        imag[power] = coefficient.im
    real = trim(real)
    imag = trim(imag)
    if not real:
        common = monic(imag)
    elif not imag:
        common = monic(real)
    else:
        common = gcd_poly(real, imag)
    count = positive_root_count(common)
    return [] if count == 0 else [{"exact_positive_root_count": count}]


def exact_phase_gate(smoke=False):
    numerator, denominator = master_polynomials(J45_STAR)
    breakpoints_n = upper_breakpoints(numerator)
    breakpoints_d = upper_breakpoints(denominator)
    breakpoints = sorted(set(breakpoints_n + breakpoints_d))
    q_tests = (
        [Fraction(0), Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3)]
        if smoke else
        [Fraction(idx, 8) for idx in range(25)]
    )
    alphas = (
        [Fraction(1, 10), Fraction(1), Fraction(10)]
        if smoke else
        [Fraction(1, 100), Fraction(1, 10), Fraction(1), Fraction(10), Fraction(100)]
    )
    rows = []
    cells_pass = True
    for q in q_tests:
        for alpha in alphas:
            order, coefficient_n, coefficient_d = path_order(
                q, numerator, denominator, alpha=alpha
            )
            expected = expected_order(q)
            cells_pass &= order == expected
            rows.append({
                "q": str(q),
                "alpha": str(alpha),
                "order": str(order),
                "expected": str(expected),
                "leading_numerator": coefficient_n.text(),
                "leading_denominator": coefficient_d.text(),
            })
    face_rows = []
    no_positive_roots = True
    for q in sorted(set(breakpoints)):
        for label, poly in [("numerator", numerator), ("denominator", denominator)]:
            face = leading_face(poly, q)
            roots = common_positive_roots(face)
            no_positive_roots &= not roots
            face_rows.append({
                "q": str(q),
                "polynomial": label,
                "face_coefficients_by_alpha_power": face_text(face),
                "positive_common_roots": roots,
            })
    exact_breakpoints_pass = breakpoints == [Fraction(1), Fraction(2)]
    return {
        "upper_breakpoints_numerator": [str(value) for value in breakpoints_n],
        "upper_breakpoints_denominator": [str(value) for value in breakpoints_d],
        "response_breakpoints": [str(value) for value in breakpoints],
        "face_rows": face_rows,
        "cell_rows": rows,
        "exact_breakpoints_pass": exact_breakpoints_pass,
        "no_positive_alpha_face_roots": bool(no_positive_roots),
        "all_cells_match_fan": bool(cells_pass),
        "overall_pass": bool(exact_breakpoints_pass and no_positive_roots and cells_pass),
    }


def unit_invariance_gate():
    numerator, denominator = master_polynomials(J45_STAR)
    rows = []
    passed = True
    # Omega_ref' = c Omega_ref implies alpha' = alpha*c^(q-1).
    for q in [Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2), Fraction(5, 2)]:
        for alpha in [Fraction(1, 10), Fraction(1), Fraction(10)]:
            baseline = path_order(q, numerator, denominator, alpha)[0]
            for scale in [Fraction(1, 10), Fraction(10)]:
                transformed = alpha * (scale ** (q - 1))
                order = path_order(q, numerator, denominator, transformed)[0]
                passed &= order == baseline
                rows.append({
                    "q": str(q),
                    "alpha": str(alpha),
                    "unit_scale": str(scale),
                    "transformed_alpha": str(transformed),
                    "order": str(order),
                })
    return {"rows": rows, "overall_pass": bool(passed)}


def finite_window_map(smoke=False):
    numerator, denominator = master_polynomials(J45_STAR)
    qs = np.linspace(0, 3, 31 if smoke else 121)
    log_alphas = np.linspace(-2, 2, 21 if smoke else 101)
    log_gamma = 6.0
    values = np.empty((len(log_alphas), len(qs)))
    for row, log_alpha in enumerate(log_alphas):
        alpha = 10 ** float(log_alpha)
        values[row] = [
            effective_order(numerator, denominator, log_gamma, float(q), alpha, dps=80)
            for q in qs
        ]
    finite = bool(np.isfinite(values).all())
    if not smoke:
        plt.figure(figsize=(7.0, 4.5))
        image = plt.imshow(
            values, origin="lower", aspect="auto",
            extent=[qs[0], qs[-1], log_alphas[0], log_alphas[-1]],
            cmap="coolwarm", vmin=2.8, vmax=4.2,
        )
        plt.colorbar(image, label=r"$\nu_{\rm eff}(\Gamma=10^6)$")
        plt.axvline(1, color="black", linestyle="--", linewidth=0.8)
        plt.axvline(2, color="black", linestyle="--", linewidth=0.8)
        plt.xlabel(r"path exponent $q$")
        plt.ylabel(r"$\log_{10}\alpha$")
        plt.title("Finite-window path-order map")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG, "p2_finite_window_phase_map.png"), dpi=180)
        plt.close()
    return {
        "Gamma": 1e6,
        "q_range": [float(qs[0]), float(qs[-1])],
        "log10_alpha_range": [float(log_alphas[0]), float(log_alphas[-1])],
        "minimum_effective_order": float(np.min(values)),
        "maximum_effective_order": float(np.max(values)),
        "all_finite": finite,
        "overall_pass": finite,
    }


def run(smoke=False):
    gates = {
        "P2_exact_phase_diagram": exact_phase_gate(smoke=smoke),
        "P2_unit_invariance": unit_invariance_gate(),
        "P2_finite_window_map": finite_window_map(smoke=smoke),
    }
    gates["priority2_overall_pass"] = all(
        gate["overall_pass"] for gate in gates.values()
        if isinstance(gate, dict) and "overall_pass" in gate
    )
    output = os.path.join(
        OUT, "priority2_smoke.json" if smoke else "priority2_production.json"
    )
    with open(output, "w") as handle:
        json.dump(gates, handle, indent=2)
    print(json.dumps({
        "mode": "smoke" if smoke else "production",
        "overall_pass": gates["priority2_overall_pass"],
        "output": output,
    }))
    return gates


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    run(smoke=args.smoke)
