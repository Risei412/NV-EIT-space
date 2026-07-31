"""Exact and stable numerical core for Phase N frequency certification.

The frequency-dependent reduced pencil is

    A(Gamma, kappa, z) = Gamma D + B - i z I + kappa D_S.

Sparse exact polynomials use keys ``(power_Gamma, power_kappa, power_z)``.
"""

from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import permutations
import math

from phase_n_exact_core import (
    D_INV,
    D_RATES,
    GaussianRational,
    I,
    J45_STAR,
    S34,
    ZERO,
    edges,
)


def gq_div(left, right):
    """Exact division of Gaussian rationals."""
    left = left if isinstance(left, GaussianRational) else GaussianRational(Fraction(left))
    right = right if isinstance(right, GaussianRational) else GaussianRational(Fraction(right))
    norm = right.re * right.re + right.im * right.im
    if norm == 0:
        raise ZeroDivisionError("Gaussian rational division by zero")
    return GaussianRational(
        (left.re * right.re + left.im * right.im) / norm,
        (left.im * right.re - left.re * right.im) / norm,
    )


def poly_add(left, right):
    out = dict(left)
    for key, value in right.items():
        out[key] = out.get(key, ZERO) + value
        if out[key].is_zero():
            del out[key]
    return out


def poly_neg(poly):
    return {key: -value for key, value in poly.items()}


def poly_mul(left, right):
    out = {}
    for (ag, ak, az), value in left.items():
        for (bg, bk, bz), other in right.items():
            key = (ag + bg, ak + bk, az + bz)
            out[key] = out.get(key, ZERO) + value * other
    return {key: value for key, value in out.items() if not value.is_zero()}


def permutation_sign(perm):
    inversions = sum(
        perm[i] > perm[j]
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
    )
    return -1 if inversions % 2 else 1


def poly_det(matrix):
    out = {}
    for perm in permutations(range(len(matrix))):
        term = {(0, 0, 0): GaussianRational(Fraction(permutation_sign(perm)))}
        for row in range(len(matrix)):
            term = poly_mul(term, matrix[row][perm[row]])
        out = poly_add(out, term)
    return out


def exact_frequency_pencil(operational=False, j45=J45_STAR):
    matrix = [[{} for _ in range(4)] for _ in range(4)]
    for idx, rate in enumerate(D_RATES):
        matrix[idx][idx] = {
            (1, 0, 0): GaussianRational(rate / 2),
            (0, 0, 1): -I,
        }
        if operational and idx in S34:
            matrix[idx][idx] = poly_add(
                matrix[idx][idx],
                {(0, 1, 0): GaussianRational(Fraction(1, 2))},
            )
    for row, col, coupling, _ in edges(j45):
        matrix[row][col] = poly_add(
            matrix[row][col], {(0, 0, 0): I * coupling}
        )
        matrix[col][row] = poly_add(
            matrix[col][row], {(0, 0, 0): I * coupling.conjugate()}
        )
    return matrix


def response_numerator(matrix):
    minor = [
        [matrix[row][col] for col in range(4) if col != 3]
        for row in range(4) if row != 0
    ]
    cofactor = poly_neg(poly_det(minor))
    return {
        key: value * GaussianRational(0, -1)
        for key, value in cofactor.items()
    }


def frequency_master_polynomials(j45=J45_STAR):
    full = exact_frequency_pencil(False, j45)
    operational = exact_frequency_pencil(True, j45)
    q_full = poly_det(full)
    q_operational = poly_det(operational)
    n_full = response_numerator(full)
    n_operational = response_numerator(operational)
    numerator = poly_add(
        poly_mul(n_full, q_operational),
        poly_neg(poly_mul(n_operational, q_full)),
    )
    return numerator, poly_mul(q_full, q_operational)


def substitute_z(poly, z):
    """Substitute an exact rational z and return a two-variable polynomial."""
    z = Fraction(z)
    out = {}
    for (g, k, power_z), coefficient in poly.items():
        key = (g, k)
        out[key] = out.get(key, ZERO) + coefficient * (z ** power_z)
    return {key: value for key, value in out.items() if not value.is_zero()}


def zpoly_add(left, right):
    out = dict(left)
    for power, value in right.items():
        out[power] = out.get(power, ZERO) + value
        if out[power].is_zero():
            del out[power]
    return out


def zpoly_mul(left, right):
    out = {}
    for power, value in left.items():
        for other_power, other in right.items():
            key = power + other_power
            out[key] = out.get(key, ZERO) + value * other
    return {key: value for key, value in out.items() if not value.is_zero()}


def zmatvec(matrix, vector):
    output = []
    for row in range(4):
        value = {}
        for col in range(4):
            value = zpoly_add(value, zpoly_mul(matrix[row][col], vector[col]))
        output.append(value)
    return output


def zdinv(vector):
    return [
        {power: D_INV[idx] * value for power, value in vector[idx].items()}
        for idx in range(4)
    ]


def ideal_moment_polynomials(kmax=7, j45=J45_STAR):
    """Return exact z-polynomials for full-minus-ideal-cut moments."""

    def bmatrix(eliminated):
        eliminated = set(eliminated)
        matrix = [[{} for _ in range(4)] for _ in range(4)]
        for idx in range(4):
            matrix[idx][idx] = {1: -I}
        for row, col, coupling, _ in edges(j45):
            if row in eliminated or col in eliminated:
                continue
            matrix[row][col] = {0: I * coupling}
            matrix[col][row] = {0: I * coupling.conjugate()}
        return matrix

    source = [{0: GaussianRational(0, -1)}, {}, {}, {}]
    full = bmatrix(())
    cut = bmatrix(S34)
    vf = zdinv(source)
    vc = zdinv(source)
    moments = []
    for _ in range(kmax):
        moments.append(zpoly_add(vf[3], {p: -v for p, v in vc[3].items()}))
        vf = zdinv(zmatvec(full, vf))
        vc = zdinv(zmatvec(cut, vc))
    return moments


def evaluate_zpoly(poly, z):
    z = Fraction(z)
    value = ZERO
    for power, coefficient in poly.items():
        value += coefficient * (z ** power)
    return value


def first_nonzero_moment_at(moments, z):
    for idx, poly in enumerate(moments):
        if not evaluate_zpoly(poly, z).is_zero():
            return idx + 1
    return None


def exact_linear_real_root(poly):
    """Return the exact real root of a nonconstant linear GQ polynomial."""
    if set(poly).issubset({0, 1}) and 1 in poly and not poly[1].is_zero():
        root = gq_div(-poly.get(0, ZERO), poly[1])
        if root.im == 0:
            return root.re
    return None


def _decimal_scaled_value(poly, log10_gamma, q, alpha, z, dps):
    """Scaled complex value and Gamma logarithmic derivative."""
    with localcontext() as context:
        context.prec = dps
        qd = Decimal(str(q))
        logg = Decimal(str(log10_gamma)) * Decimal(10).ln()
        loga = Decimal(str(alpha)).ln()
        zd = Decimal(str(z))
        active = []
        max_log_scale = None
        for (g, k, power_z), coefficient in poly.items():
            if power_z and zd == 0:
                continue
            weight = Decimal(g) + qd * Decimal(k)
            log_scale = weight * logg + Decimal(k) * loga
            z_sign = Decimal(1)
            if power_z:
                if zd < 0 and power_z % 2:
                    z_sign = Decimal(-1)
                log_scale += Decimal(power_z) * abs(zd).ln()
            if max_log_scale is None or log_scale > max_log_scale:
                max_log_scale = log_scale
            active.append((weight, log_scale, z_sign, coefficient))
        if not active:
            return (Decimal(0), Decimal(0)), (Decimal(0), Decimal(0)), Decimal(0)
        value = [Decimal(0), Decimal(0)]
        derivative = [Decimal(0), Decimal(0)]
        for weight, log_scale, z_sign, coefficient in active:
            scale = (log_scale - max_log_scale).exp() * z_sign
            real = (
                Decimal(coefficient.re.numerator)
                / Decimal(coefficient.re.denominator)
                * scale
            )
            imag = (
                Decimal(coefficient.im.numerator)
                / Decimal(coefficient.im.denominator)
                * scale
            )
            value[0] += real
            value[1] += imag
            derivative[0] += weight * real
            derivative[1] += weight * imag
        return value, derivative, max_log_scale


def point_logabs_and_order(
    numerator,
    denominator,
    log10_gamma,
    q,
    alpha,
    z,
    dps=80,
):
    """Return natural log |R| and -d log|R|/d log Gamma."""
    n, dn, nscale = _decimal_scaled_value(
        numerator, log10_gamma, q, alpha, z, dps
    )
    d, dd, dscale = _decimal_scaled_value(
        denominator, log10_gamma, q, alpha, z, dps
    )

    def norm2(value):
        return value[0] * value[0] + value[1] * value[1]

    def real_ratio(left, right):
        return (
            left[0] * right[0] + left[1] * right[1]
        ) / norm2(right)

    with localcontext() as context:
        context.prec = dps
        logabs = (
            (norm2(n).ln() - norm2(d).ln()) / Decimal(2)
            + nscale - dscale
        )
        order = real_ratio(dd, d) - real_ratio(dn, n)
        return float(logabs), float(order)


def point_scaled_complex(
    numerator,
    denominator,
    log10_gamma,
    q,
    alpha,
    z,
    rescale_order=0.0,
    dps=80,
):
    """Return Gamma**rescale_order * R as a Python complex number."""
    n, _, nscale = _decimal_scaled_value(
        numerator, log10_gamma, q, alpha, z, dps
    )
    d, _, dscale = _decimal_scaled_value(
        denominator, log10_gamma, q, alpha, z, dps
    )
    with localcontext() as context:
        context.prec = dps
        denom = d[0] * d[0] + d[1] * d[1]
        ratio_re = (n[0] * d[0] + n[1] * d[1]) / denom
        ratio_im = (n[1] * d[0] - n[0] * d[1]) / denom
        exponent = (
            nscale - dscale
            + Decimal(str(rescale_order))
            * Decimal(str(log10_gamma))
            * Decimal(10).ln()
        )
        scale = exponent.exp()
        return complex(float(ratio_re * scale), float(ratio_im * scale))


def logsumexp(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))

