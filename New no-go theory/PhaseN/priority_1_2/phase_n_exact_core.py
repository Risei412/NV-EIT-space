"""Exact polynomial core for Phase N intervention-path calculations."""

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import permutations
import math


@dataclass(frozen=True)
class GaussianRational:
    re: Fraction = Fraction(0)
    im: Fraction = Fraction(0)

    def __add__(self, other):
        other = as_gq(other)
        return GaussianRational(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return GaussianRational(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-as_gq(other))

    def __mul__(self, other):
        other = as_gq(other)
        return GaussianRational(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    __rmul__ = __mul__

    def conjugate(self):
        return GaussianRational(self.re, -self.im)

    def is_zero(self):
        return self.re == 0 and self.im == 0

    def text(self):
        if self.im == 0:
            return str(self.re)
        if self.re == 0:
            return f"{self.im}*I"
        sign = "+" if self.im >= 0 else "-"
        return f"{self.re}{sign}{abs(self.im)}*I"

def as_gq(value):
    if isinstance(value, GaussianRational):
        return value
    return GaussianRational(Fraction(value), Fraction(0))


ZERO = GaussianRational()
I = GaussianRational(0, 1)
J45_STAR = GaussianRational(Fraction(-51, 65), 0)
D_RATES = [Fraction(1), Fraction(13, 10), Fraction(17, 10), Fraction(21, 10)]
D_INV = [Fraction(2) / rate for rate in D_RATES]
S34 = (1, 2)


def edges(j45=J45_STAR):
    return [
        (0, 1, GaussianRational(Fraction(1)), "J23"),
        (0, 2, GaussianRational(Fraction(1)), "J24"),
        (0, 3, GaussianRational(Fraction(1, 10)), "J25"),
        (1, 2, GaussianRational(Fraction(1, 20)), "J34"),
        (1, 3, GaussianRational(Fraction(3, 5)), "J35"),
        (2, 3, j45, "J45"),
    ]


def exact_B(eliminated=(), j45=J45_STAR):
    eliminated = set(eliminated)
    out = [[ZERO for _ in range(4)] for _ in range(4)]
    for row, col, coupling, _ in edges(j45):
        if row in eliminated or col in eliminated:
            continue
        out[row][col] = I * coupling
        out[col][row] = I * coupling.conjugate()
    return out


def matvec(matrix, vector):
    return [
        sum((matrix[row][col] * vector[col] for col in range(4)), ZERO)
        for row in range(4)
    ]


def dinv(vector):
    return [D_INV[idx] * vector[idx] for idx in range(4)]


def exact_moments(eliminated=S34, j45=J45_STAR, kmax=7):
    full = exact_B((), j45)
    cut = exact_B(eliminated, j45)
    source = [GaussianRational(0, -1), ZERO, ZERO, ZERO]
    vf = dinv(source)
    vc = dinv(source)
    moments = []
    for _ in range(kmax):
        moments.append(vf[3] - vc[3])
        vf = dinv(matvec(full, vf))
        vc = dinv(matvec(cut, vc))
    return moments


def suppression_index(moments):
    for idx, value in enumerate(moments):
        if not value.is_zero():
            return idx + 1
    return None


# A polynomial is {(power_Gamma, power_kappa): GaussianRational}.
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
    for (ag, ak), value in left.items():
        for (bg, bk), other in right.items():
            key = (ag + bg, ak + bk)
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
        term = {(0, 0): GaussianRational(Fraction(permutation_sign(perm)))}
        for row in range(len(matrix)):
            term = poly_mul(term, matrix[row][perm[row]])
        out = poly_add(out, term)
    return out


def exact_pencil(operational=False, j45=J45_STAR):
    matrix = [[{} for _ in range(4)] for _ in range(4)]
    for idx, rate in enumerate(D_RATES):
        matrix[idx][idx] = {(1, 0): GaussianRational(rate / 2)}
        if operational and idx in S34:
            matrix[idx][idx] = poly_add(
                matrix[idx][idx],
                {(0, 1): GaussianRational(Fraction(1, 2))},
            )
    for row, col, coupling, _ in edges(j45):
        matrix[row][col] = poly_add(
            matrix[row][col], {(0, 0): I * coupling}
        )
        matrix[col][row] = poly_add(
            matrix[col][row], {(0, 0): I * coupling.conjugate()}
        )
    return matrix


def response_numerator(matrix):
    minor = [
        [matrix[row][col] for col in range(4) if col != 3]
        for row in range(4) if row != 0
    ]
    cofactor = poly_neg(poly_det(minor))
    return {key: value * GaussianRational(0, -1) for key, value in cofactor.items()}


def master_polynomials(j45=J45_STAR):
    full = exact_pencil(False, j45)
    operational = exact_pencil(True, j45)
    q_full = poly_det(full)
    q_operational = poly_det(operational)
    n_full = response_numerator(full)
    n_operational = response_numerator(operational)
    numerator = poly_add(
        poly_mul(n_full, q_operational),
        poly_neg(poly_mul(n_operational, q_full)),
    )
    return numerator, poly_mul(q_full, q_operational)


def weighted_buckets(poly, q, alpha=Fraction(1)):
    buckets = {}
    for (power_gamma, power_kappa), coefficient in poly.items():
        weight = Fraction(power_gamma) + q * Fraction(power_kappa)
        scaled = coefficient * (alpha ** power_kappa)
        buckets[weight] = buckets.get(weight, ZERO) + scaled
    return {weight: value for weight, value in buckets.items() if not value.is_zero()}


def weighted_degree(poly, q, alpha=Fraction(1)):
    buckets = weighted_buckets(poly, q, alpha)
    degree = max(buckets)
    return degree, buckets[degree]


def path_order(q, numerator, denominator, alpha=Fraction(1)):
    degree_n, coefficient_n = weighted_degree(numerator, q, alpha)
    degree_d, coefficient_d = weighted_degree(denominator, q, alpha)
    return degree_d - degree_n, coefficient_n, coefficient_d


def upper_breakpoints(poly):
    """Return nonnegative q where the leading Newton face can change."""
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


def leading_face(poly, q):
    max_weight = max(Fraction(g) + q * Fraction(k) for g, k in poly)
    face = {}
    for (g, k), coefficient in poly.items():
        if Fraction(g) + q * Fraction(k) == max_weight:
            face[k] = face.get(k, ZERO) + coefficient
    return {power: coefficient for power, coefficient in face.items()
            if not coefficient.is_zero()}


def effective_order(poly_n, poly_d, log10_gamma, q, alpha=1.0, dps=80):
    """Analytic logarithmic slope of |N/Q| along kappa=alpha*Gamma**q."""
    with localcontext() as context:
        context.prec = dps
        log_gamma = Decimal(str(log10_gamma)) * Decimal(10).ln()
        log_alpha = Decimal(str(alpha)).ln()
        q_decimal = Decimal(str(q))

        def value_and_derivative(poly):
            value = [Decimal(0), Decimal(0)]
            derivative = [Decimal(0), Decimal(0)]
            for (g, k), coefficient in poly.items():
                weight = Decimal(g) + q_decimal * Decimal(k)
                scale = (Decimal(k) * log_alpha + weight * log_gamma).exp()
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
            return value, derivative

        def real_ratio(left, right):
            denominator = right[0] * right[0] + right[1] * right[1]
            return (left[0] * right[0] + left[1] * right[1]) / denominator

        n, dn = value_and_derivative(poly_n)
        d, dd = value_and_derivative(poly_d)
        return float(real_ratio(dd, d) - real_ratio(dn, n))


def find_midpoint_crossover(
    poly_n,
    poly_d,
    q,
    alpha,
    tuned_order,
    generic_order,
    log_min,
    log_max,
    dps=80,
):
    target = 0.5 * (tuned_order + generic_order)
    grid = [
        log_min + (log_max - log_min) * idx / 500
        for idx in range(501)
    ]
    values = [
        effective_order(poly_n, poly_d, point, q, alpha, dps=dps) - target
        for point in grid
    ]
    brackets = []
    for left, right, fleft, fright in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if fleft == 0 or fleft * fright < 0:
            brackets.append((left, right))
    if not brackets:
        raise RuntimeError(f"no crossover bracket for q={q}, alpha={alpha}")
    left, right = brackets[-1]
    for _ in range(70):
        mid = 0.5 * (left + right)
        fleft = effective_order(poly_n, poly_d, left, q, alpha, dps=dps) - target
        fmid = effective_order(poly_n, poly_d, mid, q, alpha, dps=dps) - target
        if fleft * fmid <= 0:
            right = mid
        else:
            left = mid
    return 0.5 * (left + right)


def fraction_text(value):
    return str(value) if isinstance(value, Fraction) else f"{value:g}"


def linear_fit(xs, ys):
    count = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    slope = (count * sxy - sx * sy) / (count * sxx - sx * sx)
    intercept = (sy - slope * sx) / count
    residuals = [y - (slope * x + intercept) for x, y in zip(xs, ys)]
    rms = math.sqrt(sum(value * value for value in residuals) / count)
    return slope, intercept, rms
