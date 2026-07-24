"""Phase N: full-response non-identifiability of a sector order.

This calculation uses one fixed physical five-level GKSL model.  No part of
the uncut generator, probe, readout, or native dissipation is changed between
the cases.  The only additional datum is the operational sector intervention:

  S3  : add a strong decay from intermediate |3> to the ground;
  S4  : add a strong decay from intermediate |4> to the ground;
  S34 : add both strong decays.

The ideal kappa -> infinity cuts eliminate the corresponding intermediate
coherences.  Exact rational moments give nu[S3] = nu[S4] = 3 and nu[S34] = 4,
while the common uncut response has nu[chi_full] = 2.  Hence the uncut
response alone -- even the complete uncut dynamics -- cannot return a unique
sector order unless the intervention S is supplied as part of the question.

The theorem does NOT claim non-identifiability after S and a tomographically
complete generator are supplied.
"""

from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-smrt-phase-n")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
FIGDIR = os.path.join(RESULTS, "figures")
os.makedirs(FIGDIR, exist_ok=True)


# Reduced weak-probe basis: (rho21, rho31, rho41, rho51).
D_RATES = [Fraction(1), Fraction(13, 10), Fraction(17, 10), Fraction(21, 10)]
D_INV = [Fraction(2, 1) / d for d in D_RATES]

# At phi=pi the physical J45 exp(i phi) coupling is -51/65.
EDGES = [
    (0, 1, Fraction(1), "J23"),
    (0, 2, Fraction(1), "J24"),
    (0, 3, Fraction(1, 10), "J25"),
    (1, 2, Fraction(1, 20), "J34"),
    (1, 3, Fraction(3, 5), "J35"),
    (2, 3, Fraction(-51, 65), "J45_phase_pi"),
]

SECTORS = {
    "S3_single_intermediate": (1,),
    "S4_single_intermediate": (2,),
    "S34_two_intermediate": (1, 2),
}


@dataclass(frozen=True)
class GaussianRational:
    """A minimal exact complex number with rational real/imaginary parts."""

    re: Fraction = Fraction(0)
    im: Fraction = Fraction(0)

    def __add__(self, other):
        other = as_gq(other)
        return GaussianRational(self.re + other.re, self.im + other.im)

    def __radd__(self, other):
        return self + other

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

    def __rmul__(self, other):
        return self * other

    def is_zero(self):
        return self.re == 0 and self.im == 0

    def to_complex(self):
        return complex(float(self.re), float(self.im))

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


def exact_B(eliminated=()):
    """Exact B matrix for the ideal retained-block cut."""
    eliminated = set(eliminated)
    out = [[ZERO for _ in range(4)] for _ in range(4)]
    for i, j, coupling, _ in EDGES:
        if i in eliminated or j in eliminated:
            continue
        value = GaussianRational(Fraction(0), coupling)  # B = i H
        out[i][j] = value
        out[j][i] = value
    return out


def exact_matvec(matrix, vector):
    return [sum((matrix[i][j] * vector[j] for j in range(4)), ZERO)
            for i in range(4)]


def exact_Dinv(vector):
    return [D_INV[i] * vector[i] for i in range(4)]


def exact_moments(eliminated=(), kmax=6):
    """Moments of R_S = chi_full - chi_cut on the doubled realization."""
    bf = exact_B(())
    bc = exact_B(eliminated)
    source = [GaussianRational(0, -1), ZERO, ZERO, ZERO]
    vf = exact_Dinv(source)
    vc = exact_Dinv(source)
    moments = []
    for _ in range(kmax):
        moments.append(vf[3] - vc[3])
        vf = exact_Dinv(exact_matvec(bf, vf))
        vc = exact_Dinv(exact_matvec(bc, vc))
    return moments


def exact_full_moments(kmax=4):
    b = exact_B(())
    vector = exact_Dinv([GaussianRational(0, -1), ZERO, ZERO, ZERO])
    moments = []
    for _ in range(kmax):
        moments.append(vector[3])
        vector = exact_Dinv(exact_matvec(b, vector))
    return moments


def suppression_index(moments):
    for k, value in enumerate(moments):
        if not value.is_zero():
            return k + 1
    return None


# Exact bivariate-polynomial machinery for the path-order theorem. A
# polynomial is a dict {(power_Gamma, power_kappa): GaussianRational}.
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
    for (ig, ik), value in left.items():
        for (jg, jk), other in right.items():
            key = (ig + jg, ik + jk)
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
    size = len(matrix)
    out = {}
    for perm in permutations(range(size)):
        term = {(0, 0): GaussianRational(Fraction(permutation_sign(perm)), 0)}
        for row in range(size):
            term = poly_mul(term, matrix[row][perm[row]])
        out = poly_add(out, term)
    return out


def exact_pencil(operational=False):
    """A(Gamma,kappa)=Gamma D+B+kappa D_S for the joint S34 cut."""
    matrix = [[{} for _ in range(4)] for _ in range(4)]
    for idx, rate in enumerate(D_RATES):
        matrix[idx][idx] = {(1, 0): GaussianRational(rate / 2, 0)}
        if operational and idx in SECTORS["S34_two_intermediate"]:
            matrix[idx][idx] = poly_add(
                matrix[idx][idx],
                {(0, 1): GaussianRational(Fraction(1, 2), 0)},
            )
    for i, j, coupling, _ in EDGES:
        entry = {(0, 0): GaussianRational(0, coupling)}
        matrix[i][j] = poly_add(matrix[i][j], entry)
        matrix[j][i] = poly_add(matrix[j][i], entry)
    return matrix


def response_numerator(matrix):
    """Numerator of e3^T A^-1 (-i e0), before division by det(A)."""
    # (A^-1)_{3,0}=cofactor_{0,3}/det(A). Remove row 0, column 3;
    # cofactor sign is (-1)^3=-1, then multiply by the source -i.
    minor = [
        [matrix[row][col] for col in range(4) if col != 3]
        for row in range(4) if row != 0
    ]
    cofactor = poly_neg(poly_det(minor))
    return {
        key: value * GaussianRational(0, -1)
        for key, value in cofactor.items()
    }


def exact_master_rational_polynomials():
    """Return exact N,Q for R=chi_full-chi_operational=N/Q."""
    full = exact_pencil(operational=False)
    operational = exact_pencil(operational=True)
    q_full = poly_det(full)
    q_op = poly_det(operational)
    n_full = response_numerator(full)
    n_op = response_numerator(operational)
    numerator = poly_add(poly_mul(n_full, q_op), poly_neg(poly_mul(n_op, q_full)))
    denominator = poly_mul(q_full, q_op)
    return numerator, denominator


def weighted_degree(poly, q):
    """Exact degree after kappa=Gamma^q, including face cancellation."""
    buckets = {}
    for (power_gamma, power_kappa), coefficient in poly.items():
        weight = Fraction(power_gamma) + q * Fraction(power_kappa)
        buckets[weight] = buckets.get(weight, ZERO) + coefficient
    for weight in sorted(buckets, reverse=True):
        if not buckets[weight].is_zero():
            return weight, buckets[weight]
    raise ValueError("zero polynomial has no weighted degree")


def exact_path_order(q, numerator=None, denominator=None):
    if numerator is None or denominator is None:
        numerator, denominator = exact_master_rational_polynomials()
    degree_n, coefficient_n = weighted_degree(numerator, q)
    degree_q, coefficient_q = weighted_degree(denominator, q)
    return degree_q - degree_n, coefficient_n, coefficient_q


def evaluate_poly(poly, Gamma, kappa):
    return sum(
        coefficient.to_complex() * Gamma ** power_gamma * kappa ** power_kappa
        for (power_gamma, power_kappa), coefficient in poly.items()
    )


def numeric_D():
    return np.diag([float(d) / 2.0 for d in D_RATES])


def numeric_B(eliminated=()):
    eliminated = set(eliminated)
    out = np.zeros((4, 4), dtype=complex)
    for i, j, coupling, _ in EDGES:
        if i in eliminated or j in eliminated:
            continue
        value = 1j * float(coupling)
        out[i, j] = value
        out[j, i] = value
    return out


def source_vector():
    return np.array([-1j, 0.0j, 0.0j, 0.0j])


def readout_vector():
    return np.array([0.0, 0.0, 0.0, 1.0])


def response(Gamma, eliminated=(), kappa=None):
    """Reduced response: ideal cut if kappa=None, finite operational cut otherwise."""
    d = numeric_D()
    if kappa is None:
        a = Gamma * d + numeric_B(eliminated)
    else:
        ds = np.zeros((4, 4))
        for idx in eliminated:
            ds[idx, idx] = 0.5
        a = Gamma * d + numeric_B(()) + kappa * ds
    return readout_vector() @ np.linalg.solve(a, source_vector())


def fit_nu(gammas, values):
    slope = np.polyfit(np.log(gammas), np.log(np.abs(values)), 1)[0]
    return float(-slope)


# Full 5-level GKSL validation. Physical states are |1>...|5>, represented by
# Python indices 0...4. Reduced coherence indices 1 and 2 correspond to the
# physical intermediate states |3> and |4>, hence +1 below.
N_LEVELS = 5


def hamiltonian(eliminated=()):
    eliminated_physical = {idx + 1 for idx in eliminated}
    out = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    for i, j, coupling, _ in EDGES:
        pi, pj = i + 1, j + 1
        if pi in eliminated_physical or pj in eliminated_physical:
            continue
        out[pi, pj] = float(coupling)
        out[pj, pi] = float(coupling)
    return out


def full_hamiltonian():
    return hamiltonian(())


def jump_operators(Gamma, eliminated=(), kappa=None):
    jumps = []
    for physical, rate in zip(range(1, 5), D_RATES):
        op = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
        op[0, physical] = np.sqrt(Gamma * float(rate))
        jumps.append(op)
    if kappa is not None:
        for reduced_idx in eliminated:
            op = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
            op[0, reduced_idx + 1] = np.sqrt(kappa)
            jumps.append(op)
    return jumps


def liouvillian(h, jumps):
    n2 = N_LEVELS * N_LEVELS
    out = np.zeros((n2, n2), dtype=complex)
    for idx in range(n2):
        rho = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
        rho.flat[idx] = 1.0
        drho = -1j * (h @ rho - rho @ h)
        for op in jumps:
            opd = op.conj().T
            drho += op @ rho @ opd - 0.5 * (opd @ op @ rho + rho @ opd @ op)
        out[:, idx] = drho.reshape(-1)
    return out


def trace_preservation_residual(generator):
    trace_row = np.zeros(N_LEVELS * N_LEVELS, dtype=complex)
    trace_row[::N_LEVELS + 1] = 1.0
    return float(np.linalg.norm(trace_row @ generator))


def liouvillian_linear_response(Gamma, eliminated=(), kappa=None):
    # Ideal algebraic cut uses a retained-block Hamiltonian. Finite kappa uses
    # the unchanged full Hamiltonian plus valid extra GKSL jumps.
    h = hamiltonian(eliminated) if kappa is None else full_hamiltonian()
    jumps = jump_operators(Gamma, eliminated=eliminated, kappa=kappa)
    l0 = liouvillian(h, jumps)

    rho0 = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    rho0[0, 0] = 1.0
    probe = np.zeros((N_LEVELS, N_LEVELS), dtype=complex)
    probe[0, 1] = 1.0
    probe[1, 0] = 1.0
    rhs = 1j * (probe @ rho0 - rho0 @ probe)

    matrix = l0.copy()
    vector = rhs.reshape(-1)
    trace_row = np.zeros(N_LEVELS * N_LEVELS, dtype=complex)
    trace_row[::N_LEVELS + 1] = 1.0
    matrix[0, :] = trace_row
    vector[0] = 0.0
    rho1 = np.linalg.solve(matrix, vector).reshape(N_LEVELS, N_LEVELS)
    return rho1[4, 0], l0


def gate_N1_exact_nonidentifiability():
    full_moments = exact_full_moments()
    table = {}
    for name, eliminated in SECTORS.items():
        moments = exact_moments(eliminated)
        table[name] = {
            "eliminated_reduced_indices": list(eliminated),
            "moments": [value.text() for value in moments],
            "nu_R_S": suppression_index(moments),
        }
    orders = {row["nu_R_S"] for row in table.values()}
    gate = {
        "common_full_model": True,
        "common_source": "-i e_rho21",
        "common_readout": "e_rho51",
        "full_moments": [value.text() for value in full_moments],
        "nu_chi_full": suppression_index(full_moments),
        "sector_table": table,
        "distinct_sector_orders": sorted(orders),
        "theorem": (
            "No function of the uncut response alone can return a unique "
            "intervention-resolved order when S is not supplied."
        ),
        "scope_caveat": (
            "The result does not apply after the intervention S and a "
            "tomographically complete generator are supplied."
        ),
    }
    gate["overall_pass"] = bool(
        gate["nu_chi_full"] == 2
        and table["S3_single_intermediate"]["nu_R_S"] == 3
        and table["S4_single_intermediate"]["nu_R_S"] == 3
        and table["S34_two_intermediate"]["nu_R_S"] == 4
    )
    return gate


def gate_N2_direct_fits():
    gammas = np.logspace(2, 3, 31)
    full = np.array([response(g) for g in gammas])
    table = {}
    curves = {}
    for name, eliminated in SECTORS.items():
        cut = np.array([response(g, eliminated=eliminated) for g in gammas])
        master = full - cut
        table[name] = {
            "nu_R_S_fit": fit_nu(gammas, master),
            "nu_R_S_exact": suppression_index(exact_moments(eliminated)),
        }
        curves[name] = master
    gate = {
        "gammas": gammas.tolist(),
        "nu_chi_full_fit": fit_nu(gammas, full),
        "nu_chi_full_exact": suppression_index(exact_full_moments()),
        "max_full_response_spread_between_sector_labels": 0.0,
        "table": table,
    }
    gate["overall_pass"] = bool(
        abs(gate["nu_chi_full_fit"] - 2.0) < 0.01
        and all(abs(row["nu_R_S_fit"] - row["nu_R_S_exact"]) < 0.01
                for row in table.values())
    )
    return gate, gammas, full, curves


def gate_N3_operational_admissibility():
    rows = []
    max_trace = 0.0
    max_hermiticity = 0.0
    max_reduced_full = 0.0
    for name, eliminated in SECTORS.items():
        for mode, kappa in [("ideal", None), ("finite_kappa", 100.0)]:
            h = hamiltonian(eliminated) if kappa is None else full_hamiltonian()
            jumps = jump_operators(17.0, eliminated=eliminated, kappa=kappa)
            generator = liouvillian(h, jumps)
            hermiticity = float(np.linalg.norm(h - h.conj().T))
            trace_residual = trace_preservation_residual(generator)
            max_hermiticity = max(max_hermiticity, hermiticity)
            max_trace = max(max_trace, trace_residual)
            rows.append({
                "sector": name,
                "mode": mode,
                "hamiltonian_hermiticity_residual": hermiticity,
                "trace_preservation_residual": trace_residual,
                "number_of_jumps": len(jumps),
            })
        for gamma in [10.0, 100.0, 1000.0]:
            reduced = response(gamma, eliminated=eliminated)
            full_liouville, _ = liouvillian_linear_response(
                gamma, eliminated=eliminated, kappa=None
            )
            max_reduced_full = max(max_reduced_full, abs(reduced - full_liouville))
        for gamma, kappa in [(10.0, 50.0), (100.0, 500.0)]:
            reduced = response(gamma, eliminated=eliminated, kappa=kappa)
            full_liouville, _ = liouvillian_linear_response(
                gamma, eliminated=eliminated, kappa=kappa
            )
            max_reduced_full = max(max_reduced_full, abs(reduced - full_liouville))
    gate = {
        "all_native_and_cut_rates_nonnegative": True,
        "noninvasive_reference_state": "All added jumps annihilate |1><1|.",
        "fixed_probe_and_readout": True,
        "semisimple_cut_zero_mode": True,
        "rows": rows,
        "max_hamiltonian_hermiticity_residual": max_hermiticity,
        "max_trace_preservation_residual": max_trace,
        "max_reduced_vs_full_liouvillian_error": max_reduced_full,
    }
    gate["overall_pass"] = bool(
        max_hermiticity < 1e-12
        and max_trace < 1e-12
        and max_reduced_full < 1e-12
    )
    return gate


def gate_N4_operational_limit_and_path():
    fixed_gamma = 10.0
    kappas = np.logspace(2, 7, 31)
    convergence = {}
    for name, eliminated in SECTORS.items():
        ideal = response(fixed_gamma, eliminated=eliminated)
        errors = np.array([
            abs(response(fixed_gamma, eliminated=eliminated, kappa=kappa) - ideal)
            for kappa in kappas
        ])
        slope = float(np.polyfit(np.log(kappas[-15:]), np.log(errors[-15:]), 1)[0])
        convergence[name] = {
            "slope": slope,
            "error_at_kappa_max": float(errors[-1]),
        }

    # The promoted S34 order is recovered only if the intervention scale grows
    # sufficiently rapidly relative to native Gamma. At kappa~Gamma a residual
    # finite-cut term demotes the observed order from 4 to 3.
    gammas = np.logspace(1, 2.5, 31)
    path_table = {}
    for name in ["S3_single_intermediate", "S34_two_intermediate"]:
        eliminated = SECTORS[name]
        path_table[name] = {}
        for power in [1, 2, 3]:
            values = np.array([
                response(g) - response(g, eliminated=eliminated, kappa=g ** power)
                for g in gammas
            ])
            path_table[name][f"kappa_equals_Gamma_power_{power}"] = fit_nu(
                gammas[-15:], values[-15:]
            )

    gate = {
        "fixed_Gamma_for_kappa_limit": fixed_gamma,
        "kappas": kappas.tolist(),
        "convergence": convergence,
        "path_table": path_table,
        "interpretation": (
            "For S34, kappa proportional to Gamma leaves a Gamma^-3 "
            "finite-cut term; kappa growing as Gamma^2 or faster recovers "
            "the promoted exponent 4."
        ),
    }
    gate["overall_pass"] = bool(
        all(abs(row["slope"] + 1.0) < 0.01 for row in convergence.values())
        and abs(path_table["S3_single_intermediate"]["kappa_equals_Gamma_power_1"] - 3.0) < 0.03
        and abs(path_table["S34_two_intermediate"]["kappa_equals_Gamma_power_1"] - 3.0) < 0.03
        and abs(path_table["S34_two_intermediate"]["kappa_equals_Gamma_power_2"] - 4.0) < 0.03
        and abs(path_table["S34_two_intermediate"]["kappa_equals_Gamma_power_3"] - 4.0) < 0.03
    )
    return gate


def gate_N5_exact_path_order_fan():
    numerator, denominator = exact_master_rational_polynomials()
    q_values = [Fraction(k, 4) for k in range(0, 13)]
    exact_rows = []
    direct_rows = []
    # Evaluate the exact rational response deep in the asymptotic regime.
    # Matrix subtraction would lose the promoted Gamma^-4 difference here;
    # evaluating the already-certified N/Q avoids that cancellation error.
    gammas = np.logspace(8, 12, 61)
    for q in q_values:
        nu, leading_n, leading_q = exact_path_order(q, numerator, denominator)
        exact_rows.append({
            "q": str(q),
            "q_float": float(q),
            "nu_exact": str(nu),
            "nu_exact_float": float(nu),
            "leading_numerator_coefficient": leading_n.text(),
            "leading_denominator_coefficient": leading_q.text(),
        })
        values = np.array([
            evaluate_poly(numerator, g, g ** float(q))
            / evaluate_poly(denominator, g, g ** float(q))
            for g in gammas
        ])
        nu_fit = fit_nu(gammas[-25:], values[-25:])
        direct_rows.append({
            "q": float(q),
            "nu_fit": nu_fit,
            "nu_exact": float(nu),
            "abs_error": abs(nu_fit - float(nu)),
        })

    numerator_terms = [
        {
            "power_Gamma": powers[0],
            "power_kappa": powers[1],
            "coefficient": coefficient.text(),
        }
        for powers, coefficient in sorted(numerator.items())
    ]
    gate = {
        "theorem": (
            "For finite-dimensional R=N/Q and kappa=Gamma^q, "
            "nu(q)=deg_q(Q)-deg_q(N), with deg_q the exact weighted "
            "Newton degree after face cancellation."
        ),
        "consequence": "nu(q) is finitely decidable and piecewise linear in q.",
        "phase_n_fan": (
            "nu(q)=4-q for 0<=q<=1; nu(q)=2+q for 1<=q<=2; "
            "nu(q)=4 for q>=2."
        ),
        "breakpoints": [1.0, 2.0],
        "minimum_order": 3.0,
        "minimum_at_q": 1.0,
        "numerator_terms": numerator_terms,
        "denominator_term_count": len(denominator),
        "exact_rows": exact_rows,
        "exact_rational_fit_rows": direct_rows,
        "max_exact_rational_fit_error": max(row["abs_error"] for row in direct_rows),
    }
    expected = []
    for q in q_values:
        if q <= 1:
            expected.append(Fraction(4) - q)
        elif q <= 2:
            expected.append(Fraction(2) + q)
        else:
            expected.append(Fraction(4))
    obtained = [Fraction(row["nu_exact"]) for row in exact_rows]
    gate["overall_pass"] = bool(
        obtained == expected and gate["max_exact_rational_fit_error"] < 0.002
    )
    return gate


def make_figure(gammas, full, curves, gate_n2):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    # There is only one full curve: the intervention is not part of any uncut
    # measurement. State the exact overlap rather than hiding it under three
    # visually indistinguishable colored lines.
    axes[0].loglog(gammas, np.abs(full), color="black", linewidth=2.4,
                   label="common full response")
    axes[0].text(0.05, 0.08, "S3, S4, S34 overlap exactly\nmax spread = 0",
                 transform=axes[0].transAxes, fontsize=9)
    axes[0].set_xlabel(r"$\Gamma$")
    axes[0].set_ylabel(r"$|\chi_{\mathrm{full}}|$")
    axes[0].set_title(r"Identical uncut response ($\nu_{\rm full}=2$)")
    axes[0].legend(fontsize=7)

    nu_single = gate_n2["table"]["S3_single_intermediate"]["nu_R_S_fit"]
    axes[1].loglog(
        gammas, np.abs(curves["S3_single_intermediate"]), linewidth=2.4,
        label=f"S3, S4: nu={nu_single:.3f} (opposite phase)",
    )
    nu_joint = gate_n2["table"]["S34_two_intermediate"]["nu_R_S_fit"]
    axes[1].loglog(
        gammas, np.abs(curves["S34_two_intermediate"]), linewidth=2.4,
        label=f"S34 joint: nu={nu_joint:.3f}",
    )
    axes[1].set_xlabel(r"$\Gamma$")
    axes[1].set_ylabel(r"$|R_S|$")
    axes[1].set_title("Different operational-sector orders")
    axes[1].legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figN1_full_response_nonidentifiability.png"), dpi=180)
    plt.close(fig)


def make_path_fan_figure(gate_n5):
    exact = gate_n5["exact_rows"]
    direct = gate_n5["exact_rational_fit_rows"]
    q_exact = [row["q_float"] for row in exact]
    nu_exact = [row["nu_exact_float"] for row in exact]
    q_fit = [row["q"] for row in direct]
    nu_fit = [row["nu_fit"] for row in direct]

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(q_exact, nu_exact, color="black", linewidth=2.4,
            label="exact Newton-degree law")
    ax.scatter(q_fit, nu_fit, color="tab:red", s=28, zorder=3,
               label="exact-rational response fits")
    ax.axvline(1.0, color="0.6", linestyle="--", linewidth=1)
    ax.axvline(2.0, color="0.6", linestyle="--", linewidth=1)
    ax.set_xlabel(r"intervention path exponent $q$ in $\kappa=\Gamma^q$")
    ax.set_ylabel(r"observed sector order $\nu(q)$")
    ax.set_title("Intervention-scaling fan for the promoted S34 sector")
    ax.set_xlim(0, 3)
    ax.set_ylim(2.85, 4.12)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figN2_intervention_scaling_fan.png"), dpi=180)
    plt.close(fig)


def main():
    gates = {}
    gates["N1_exact_nonidentifiability"] = gate_N1_exact_nonidentifiability()
    gates["N2_direct_fits"], gammas, full, curves = gate_N2_direct_fits()
    gates["N3_operational_admissibility"] = gate_N3_operational_admissibility()
    gates["N4_operational_limit_and_path"] = gate_N4_operational_limit_and_path()
    gates["N5_exact_path_order_fan"] = gate_N5_exact_path_order_fan()
    gates["phase_n_overall_pass"] = bool(all(
        value["overall_pass"] for value in gates.values()
        if isinstance(value, dict) and "overall_pass" in value
    ))

    make_figure(gammas, full, curves, gates["N2_direct_fits"])
    make_path_fan_figure(gates["N5_exact_path_order_fan"])
    output = os.path.join(RESULTS, "gates_summary_phaseN.json")
    with open(output, "w") as handle:
        json.dump(gates, handle, indent=2)

    for name, gate in gates.items():
        if isinstance(gate, dict) and "overall_pass" in gate:
            print(f"{name}: {gate['overall_pass']}")
    print(f"phase_n_overall_pass: {gates['phase_n_overall_pass']}")
    print(output)


if __name__ == "__main__":
    main()
