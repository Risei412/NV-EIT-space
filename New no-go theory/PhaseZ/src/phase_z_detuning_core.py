"""Exact core for Phase Z: physical detuning signature z*Z.

Generalizes the Priority 3 common-shift pencil

    A(Gamma, kappa, z) = Gamma D + B + kappa D_S - i z I

to a diagonal physical detuning signature Z = diag(zeta_1..zeta_4),

    A(Gamma, kappa, z; Z) = Gamma D + B + kappa D_S - i z Z,

acting on the optical coherences (rho21, rho31, rho41, rho51) of the
five-level Phase N diamond.  The physical origin of Z is a static level
shift  delta_H = -z * diag(0, zeta_1, .., zeta_4)  on states |2>..|5> in
the rotating frame, i.e. a co-sweep of laser detunings (Gate Z5 certifies
this correspondence against the full 25x25 Liouvillian, exactly).

Sparse exact polynomials use keys ``(power_Gamma, power_kappa, power_z)``
with GaussianRational coefficients, as in ``phase_n_frequency_core``.
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
    as_gq,
    edges,
)

# --- detuning signatures -------------------------------------------------
# Reduced coordinate k <-> coherence rho_{(k+2),1} (physical labels |1..5>).
# Rotating-frame spanning tree: eps1=0, eps2=wp, eps3=wp+w23, eps4=wp+w24,
# eps5=wp+w23+w35, with the loop constraint w45 = w23+w35-w24 (the J45 laser
# is co-swept to keep the frame static).  Sweeping one generator frequency
# by z shifts each coherence rho_{k1} by z times the entries below.
Z_SIGNATURES = {
    "probe_common": (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    "coupling23": (Fraction(0), Fraction(1), Fraction(0), Fraction(1)),
    "coupling24": (Fraction(0), Fraction(0), Fraction(1), Fraction(0)),
    "coupling35": (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
    "two_tone_23_24": (Fraction(0), Fraction(1), Fraction(1), Fraction(1)),
    "generic_control": (
        Fraction(1), Fraction(7, 5), Fraction(3, 10), Fraction(2)
    ),
}

Z_IDENTITY = Z_SIGNATURES["probe_common"]


def gq_div(left, right):
    """Exact division of Gaussian rationals."""
    left = as_gq(left)
    right = as_gq(right)
    norm = right.re * right.re + right.im * right.im
    if norm == 0:
        raise ZeroDivisionError("Gaussian rational division by zero")
    return GaussianRational(
        (left.re * right.re + left.im * right.im) / norm,
        (left.im * right.re - left.re * right.im) / norm,
    )


# --- three-variable sparse polynomial algebra ----------------------------

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


# --- Z-weighted pencil and master polynomials ----------------------------

def exact_detuned_pencil(zeta, operational=False, j45=J45_STAR):
    """4x4 pencil Gamma*D + B + kappa*D_S - i z Z as a polynomial matrix."""
    matrix = [[{} for _ in range(4)] for _ in range(4)]
    for idx, rate in enumerate(D_RATES):
        entry = {(1, 0, 0): GaussianRational(rate / 2)}
        weight = Fraction(zeta[idx])
        if weight != 0:
            entry[(0, 0, 1)] = GaussianRational(Fraction(0), -weight)
        matrix[idx][idx] = entry
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


def detuned_master_polynomials(zeta, j45=J45_STAR):
    """Exact numerator and denominator of R_S(Gamma, kappa, z; Z)."""
    full = exact_detuned_pencil(zeta, False, j45)
    operational = exact_detuned_pencil(zeta, True, j45)
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
    """Substitute an exact rational z; return a (Gamma, kappa) polynomial."""
    z = Fraction(z)
    out = {}
    for (g, k, power_z), coefficient in poly.items():
        key = (g, k)
        out[key] = out.get(key, ZERO) + coefficient * (z ** power_z)
    return {key: value for key, value in out.items() if not value.is_zero()}


# --- ideal-cut moment ladder (z-polynomials) ------------------------------

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


def ideal_moment_polynomials(zeta, kmax=7, j45=J45_STAR):
    """Exact z-polynomials for the full-minus-ideal-cut moment ladder."""

    def bmatrix(eliminated):
        eliminated = set(eliminated)
        matrix = [[{} for _ in range(4)] for _ in range(4)]
        for idx in range(4):
            weight = Fraction(zeta[idx])
            if weight != 0:
                matrix[idx][idx] = {
                    1: GaussianRational(Fraction(0), -weight)
                }
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


def first_nonzero_moment_index(moments):
    """1-based index of the first moment that is nonzero as a polynomial."""
    for idx, poly in enumerate(moments):
        if poly:
            return idx + 1
    return None


# --- exact real-root machinery over Fraction[z] ---------------------------
# A real polynomial is a list of Fraction coefficients [c0, c1, ...].

def zpoly_real_imag(poly):
    degree = max(poly) if poly else 0
    real = [Fraction(0)] * (degree + 1)
    imag = [Fraction(0)] * (degree + 1)
    for power, coefficient in poly.items():
        real[power] = coefficient.re
        imag[power] = coefficient.im
    return trim(real), trim(imag)


def trim(coeffs):
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def rpoly_eval(coeffs, x):
    value = Fraction(0)
    for coefficient in reversed(coeffs):
        value = value * x + coefficient
    return value


def rpoly_divmod(num, den):
    num = list(num)
    quotient = [Fraction(0)] * max(len(num) - len(den) + 1, 0)
    while len(num) >= len(den) and num:
        factor = num[-1] / den[-1]
        shift = len(num) - len(den)
        quotient[shift] = factor
        for idx, coefficient in enumerate(den):
            num[shift + idx] -= factor * coefficient
        num = trim(num)
    return trim(quotient), num


def rpoly_gcd(left, right):
    left, right = trim(list(left)), trim(list(right))
    while right:
        _, remainder = rpoly_divmod(left, right)
        left, right = right, remainder
    if left:
        lead = left[-1]
        left = [coefficient / lead for coefficient in left]
    return left


def rpoly_derivative(coeffs):
    return trim([
        coefficient * idx for idx, coefficient in enumerate(coeffs)
    ][1:])


def squarefree_part(coeffs):
    if len(coeffs) <= 1:
        return list(coeffs)
    g = rpoly_gcd(coeffs, rpoly_derivative(coeffs))
    if len(g) <= 1:
        return list(coeffs)
    quotient, _ = rpoly_divmod(coeffs, g)
    return quotient


def sturm_chain(coeffs):
    chain = [trim(list(coeffs)), rpoly_derivative(coeffs)]
    while chain[-1]:
        _, remainder = rpoly_divmod(chain[-2], chain[-1])
        if not remainder:
            break
        chain.append([-c for c in remainder])
    return [entry for entry in chain if entry]


def sign_changes(chain, x):
    signs = []
    for entry in chain:
        value = rpoly_eval(entry, x)
        if value != 0:
            signs.append(1 if value > 0 else -1)
    return sum(
        1 for left, right in zip(signs, signs[1:]) if left != right
    )


def count_real_roots(coeffs, lo, hi):
    """Distinct real roots of coeffs in (lo, hi], exact Sturm count."""
    sf = squarefree_part(coeffs)
    if len(sf) <= 1:
        return 0
    chain = sturm_chain(sf)
    return sign_changes(chain, Fraction(lo)) - sign_changes(chain, Fraction(hi))


def isolate_real_roots(coeffs, lo, hi, width=Fraction(1, 10 ** 12)):
    """Disjoint rational isolating intervals for distinct real roots."""
    sf = squarefree_part(coeffs)
    if len(sf) <= 1:
        return []
    chain = sturm_chain(sf)

    def count(a, b):
        return sign_changes(chain, a) - sign_changes(chain, b)

    stack = [(Fraction(lo), Fraction(hi))]
    intervals = []
    while stack:
        a, b = stack.pop()
        n = count(a, b)
        if n == 0:
            continue
        if n == 1 and b - a <= width:
            intervals.append((a, b))
            continue
        mid = (a + b) / 2
        if rpoly_eval(sf, mid) == 0:
            intervals.append((mid, mid))
            shrink = min(width, (b - a) / 4)
            stack.append((a, mid - shrink))
            stack.append((mid + shrink, b))
            continue
        stack.append((a, mid))
        stack.append((mid, b))
    return sorted(intervals)


def exact_rational_root_if_linear(coeffs):
    if len(coeffs) == 2:
        return -coeffs[0] / coeffs[1]
    return None


def real_root_constraint(poly):
    """Real polynomial whose real roots are exactly the real roots of the
    GaussianRational z-polynomial ``poly`` (gcd of real and imaginary parts).
    Returns None when poly is identically zero."""
    real, imag = zpoly_real_imag(poly)
    if not real and not imag:
        return None
    if not real:
        return imag
    if not imag:
        return real
    return rpoly_gcd(real, imag)


def promotion_is_exactly_one_order(constraint, next_moment):
    """True when no root of ``constraint`` is also a root of the next
    moment polynomial: gcd(constraint, Re next, Im next) is constant."""
    real, imag = zpoly_real_imag(next_moment)
    g = list(constraint)
    if real:
        g = rpoly_gcd(g, real)
    if imag:
        g = rpoly_gcd(g, imag)
    if not real and not imag:
        return False
    return len(g) <= 1


# --- weighted Newton data on (Gamma, kappa) after z substitution ----------

def weighted_degree2(poly, q, alpha=Fraction(1)):
    buckets = {}
    for (power_gamma, power_kappa), coefficient in poly.items():
        weight = Fraction(power_gamma) + Fraction(q) * Fraction(power_kappa)
        scaled = coefficient * as_gq(Fraction(alpha) ** power_kappa)
        buckets[weight] = buckets.get(weight, ZERO) + scaled
    buckets = {w: v for w, v in buckets.items() if not v.is_zero()}
    degree = max(buckets)
    return degree, buckets[degree]


def upper_breakpoints2(poly):
    support = list(poly)
    candidates = set()
    for (ag, ak), (bg, bk) in permutations(support, 2):
        if ak == bk:
            continue
        q = Fraction(bg - ag, ak - bk)
        if q < 0:
            continue
        weights = [Fraction(g) + q * Fraction(k) for g, k in support]
        target = Fraction(ag) + q * Fraction(ak)
        if target == max(weights) and weights.count(target) >= 2:
            candidates.add(q)
    return sorted(candidates)


def exact_fan(numerator2, denominator2, q_max=Fraction(4)):
    """Exact piecewise-affine path order nu(q) = d_q(Q) - d_q(N) on [0, q_max].

    Returns a list of cells {"q_lo", "q_hi", "slope", "intercept"} with exact
    Fractions, merged where adjacent cells continue the same affine law."""
    points = sorted({
        Fraction(0), Fraction(q_max),
        *[q for q in upper_breakpoints2(numerator2) if q <= q_max],
        *[q for q in upper_breakpoints2(denominator2) if q <= q_max],
    })
    cells = []
    for lo, hi in zip(points[:-1], points[1:]):
        q1 = lo + (hi - lo) / 3
        q2 = lo + 2 * (hi - lo) / 3
        order1 = (weighted_degree2(denominator2, q1)[0]
                  - weighted_degree2(numerator2, q1)[0])
        order2 = (weighted_degree2(denominator2, q2)[0]
                  - weighted_degree2(numerator2, q2)[0])
        slope = (order2 - order1) / (q2 - q1)
        intercept = order1 - slope * q1
        cells.append({"q_lo": lo, "q_hi": hi,
                      "slope": slope, "intercept": intercept})
    merged = [cells[0]]
    for cell in cells[1:]:
        last = merged[-1]
        if (cell["slope"] == last["slope"]
                and cell["intercept"] == last["intercept"]):
            last["q_hi"] = cell["q_hi"]
        else:
            merged.append(cell)
    return merged


def fan_order_at(fan, q):
    q = Fraction(q)
    for cell in fan:
        if cell["q_lo"] <= q <= cell["q_hi"]:
            return cell["slope"] * q + cell["intercept"]
    raise ValueError(f"q={q} outside fan range")


def fan_text(fan):
    return [
        {
            "q_lo": str(cell["q_lo"]),
            "q_hi": str(cell["q_hi"]),
            "slope": str(cell["slope"]),
            "intercept": str(cell["intercept"]),
            "law": f"{cell['intercept']} + ({cell['slope']})*q",
        }
        for cell in fan
    ]


def fans_equal(left, right):
    return fan_text(left) == fan_text(right)


# --- high-precision decimal evaluation along paths ------------------------

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
    numerator, denominator, log10_gamma, q, alpha, z, dps=80,
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
        return (left[0] * right[0] + left[1] * right[1]) / norm2(right)

    with localcontext() as context:
        context.prec = dps
        logabs = (
            (norm2(n).ln() - norm2(d).ln()) / Decimal(2)
            + nscale - dscale
        )
        order = real_ratio(dd, d) - real_ratio(dn, n)
        return float(logabs), float(order)


def logsumexp(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(v - maximum) for v in values))
