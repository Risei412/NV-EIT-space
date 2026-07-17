"""Phase P driver: physical Lindblad realization of the SMRT interference-
controlled exponent promotion (paper3_smrt_numerical_plan.md, runs P1-P8).

Outputs into ../results/figures and ../results/gates_summary_phaseP.json.
"""

import json
import os
import time

import numpy as np
import sympy as sp
import mpmath as mp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core import fit_nu_loglog
import model_physical as M

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
FIGDIR = os.path.join(RESULTS, "figures")
os.makedirs(FIGDIR, exist_ok=True)

F64_SAFE_GAMMA = 1e6  # beyond this, switch to mpmath for the tuned/near-tuned cases


def _default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, sp.Basic):
        return str(o)
    return str(o)


def sweep_R_S(J45v, phiv, gmin=0, gmax=8, n_per_decade=24, use_mp_tail=True, z_val=0):
    """Gamma sweep of R_S using float64 up to F64_SAFE_GAMMA, mpmath beyond."""
    R_S_f64, R_S_mp = M.lambdify_R_S(z_val, J45v, phiv)
    n = int((gmax - gmin) * n_per_decade) + 1
    gammas = np.logspace(gmin, gmax, n)
    vals = np.zeros(n, dtype=complex)
    conds = np.zeros(n)
    for i, g in enumerate(gammas):
        if g <= F64_SAFE_GAMMA or not use_mp_tail:
            r, c = R_S_f64(g)
            vals[i] = r
            conds[i] = c
        else:
            r = R_S_mp(g, dps=50)
            vals[i] = complex(r)
            conds[i] = np.nan
    return gammas, vals, conds


# ----------------------------------------------------------------------
# P1: direct validation of nu 3 -> 4
# ----------------------------------------------------------------------
def run_P1():
    gate = {}
    J45star = M.J45_star(phi_val=sp.pi)
    cases = {
        "generic": (sp.Rational(1, 2), sp.Integer(0)),
        "tuned": (J45star, sp.pi),
        "detuned": (J45star * sp.Rational(101, 100), sp.pi),
    }

    mu_pred = {}
    for name, (J45v, phiv) in cases.items():
        mu = M.moments_difference(0, J45v, phiv, kmax=4)
        mu = [sp.nsimplify(sp.simplify(m), rational=True) for m in mu]
        nz = [i for i, m in enumerate(mu) if m != 0]
        nu_moment = (nz[0] + 1) if nz else None
        mu_pred[name] = {"moments": [str(m) for m in mu], "nu_moment": nu_moment}

    results = {}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for ax_i, (name, (J45v, phiv)) in enumerate(cases.items()):
        gammas, vals, conds = sweep_R_S(J45v, phiv)
        absvals = np.abs(vals)
        tail_mask = gammas >= gammas[-1] / 1e3  # top 3 decades
        fit = fit_nu_loglog(gammas[tail_mask], vals[tail_mask])
        nu_fit = fit["nu_global"]
        results[name] = {
            "nu_fit": nu_fit,
            "nu_moment": mu_pred[name]["nu_moment"],
            "diff": abs(nu_fit - mu_pred[name]["nu_moment"]),
            "max_cond_f64": float(np.nanmax(conds)),
        }

        ax = axes[ax_i]
        ax.loglog(gammas, absvals, "-", lw=1, label="$|R_S|$")
        ax.loglog(gammas, absvals * gammas**3, "--", lw=1, label=r"$\Gamma^3|R_S|$")
        ax.loglog(gammas, absvals * gammas**4, ":", lw=1, label=r"$\Gamma^4|R_S|$")
        ax.set_title(f"{name}: $\\nu_{{fit}}$={nu_fit:.3f}, $\\nu_m$={mu_pred[name]['nu_moment']}")
        ax.set_xlabel(r"$\Gamma$")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figP1_promotion.png"), dpi=150)
    plt.close(fig)

    gate["cases"] = results
    gate["moments"] = mu_pred
    gate["gate_close_fit_moment"] = bool(all(v["diff"] < 0.03 for v in results.values()))
    gate["gate_tuned_detuned_separated"] = bool(
        abs(results["tuned"]["nu_fit"] - results["detuned"]["nu_fit"]) > 0.5
    )
    gate["overall_pass"] = bool(gate["gate_close_fit_moment"] and gate["gate_tuned_detuned_separated"])
    gate["J45_star"] = str(J45star)
    return gate


# ----------------------------------------------------------------------
# P2: symbolic moment cancellation
# ----------------------------------------------------------------------
def run_P2():
    gate = {}
    J45star = M.J45_star(phi_val=sp.pi)

    # m2 closed form and its cancellation condition
    mu_generic = M.moments_difference(0, sp.Symbol("J45v", real=True), sp.pi, kmax=3)
    m2_expr = sp.simplify(mu_generic[2])
    gate["m2_symbolic_form"] = str(m2_expr)

    mu_star = M.moments_difference(0, J45star, sp.pi, kmax=4)
    mu_star = [sp.nsimplify(sp.simplify(m), rational=True) for m in mu_star]
    gate["moments_at_tuned_point"] = [str(m) for m in mu_star]
    gate["m2_exactly_zero"] = bool(mu_star[2] == 0)
    gate["m3_nonzero"] = bool(mu_star[3] != 0)
    gate["m3_value"] = str(mu_star[3])

    # z-window check: m3(z) on a sampled rational grid
    m3_on_window = []
    for zv in [sp.Rational(-1, 2), sp.Rational(-1, 4), 0, sp.Rational(1, 4), sp.Rational(1, 2)]:
        mu_z = M.moments_difference(zv, J45star, sp.pi, kmax=4)
        m3z = sp.nsimplify(sp.simplify(mu_z[3]), rational=True)
        m3_on_window.append((str(zv), str(m3z), float(abs(complex(m3z)))))
    gate["m3_on_z_window"] = m3_on_window
    gate["m3_nonzero_on_window"] = bool(all(v[2] > 1e-6 for v in m3_on_window))

    # cross-check: float64 moment ladder at 5 random rational parameter draws
    rng = np.random.default_rng(0)
    checks = []
    for _ in range(5):
        d3v = float(1 + rng.random())
        d4v = float(1 + rng.random())
        J23v, J24v, J35v = [float(0.5 + rng.random()) for _ in range(3)]
        params = dict(d3=sp.nsimplify(d3v, rational=False), d4=sp.nsimplify(d4v, rational=False),
                      J23=sp.Float(J23v), J24=sp.Float(J24v), J35=sp.Float(J35v))
        J45v_r = (J23v * J35v * d4v) / (J24v * d3v)  # cancellation at phi=pi
        mu_r = M.moments_difference(0, sp.Float(J45v_r), sp.pi, kmax=3, params=params)
        m2_val = complex(sp.simplify(mu_r[2]))
        checks.append(abs(m2_val))
    gate["random_draw_m2_max_abs"] = float(max(checks))
    gate["random_draw_cancellation_confirmed"] = bool(max(checks) < 1e-8)

    gate["overall_pass"] = bool(
        gate["m2_exactly_zero"] and gate["m3_nonzero"] and gate["m3_nonzero_on_window"]
        and gate["random_draw_cancellation_confirmed"]
    )
    return gate


# ----------------------------------------------------------------------
# P3: cancellation manifold map
# ----------------------------------------------------------------------
def run_P3():
    gate = {}
    p = dict(M.DEFAULT_PARAMS)

    # Exact symbolic m2(J45_mag, phi): a SINGLE complex equation in the two
    # real controls (J45_mag, phi) generically has an ISOLATED zero, not a
    # curve -- m2=0 is 2 real equations (Re, Im) for 2 real unknowns. This
    # is verified directly below and reported honestly (it matches the
    # precedent in Phase M gate M4: an interior zero with only 2 real
    # controls needs the equation to be non-generic; here it is generic).
    J45m, phis = sp.symbols("J45m phis", real=True)
    mu2 = M.moments_difference(0, J45m, phis, kmax=3)
    m2_expr = sp.simplify(mu2[2])
    m2_func = sp.lambdify((J45m, phis), m2_expr, "numpy")

    J45star_mag = float(M.J45_star(phi_val=sp.pi))
    etas = np.linspace(0.05, 2.5, 121) * J45star_mag  # J45 magnitude directly
    phis_grid = np.linspace(0, 2 * np.pi, 121)
    E, P = np.meshgrid(etas, phis_grid, indexing="ij")
    M2 = m2_func(E, P)
    logabs = np.log10(np.abs(M2) + 1e-300)

    # exact zero at the point (already proven symbolically: m2(J45star,pi)=0).
    # A discrete (magnitude,phase) grid generically misses this measure-zero
    # point, so isolatedness is checked differently: away from a small disk
    # around the exact point, |m2| must stay bounded away from zero (no
    # OTHER zero curve threading the domain).
    exact_zero_val = abs(complex(m2_func(J45star_mag, np.pi)))
    dist_to_point = np.sqrt((E / J45star_mag - 1.0) ** 2 + np.minimum(
        np.abs(P - np.pi), 2 * np.pi - np.abs(P - np.pi)) ** 2)
    far_mask = dist_to_point > 0.3
    min_logabs_far = float(np.min(logabs[far_mask]))
    gate["exact_zero_value"] = exact_zero_val
    gate["min_logabs_away_from_point"] = min_logabs_far
    gate["isolated_zero_confirmed"] = bool(exact_zero_val < 1e-9 and min_logabs_far > -1.0)

    # A genuine codimension-1 cancellation CURVE requires a third independent
    # real control: rescale the branch-1 amplitude by r (J23 -> r*J23) and
    # fix phi=pi (where m2 is automatically real, since Im m2(phi=pi)=const
    # term only). Then m2(J45_mag, r; phi=pi)=0 is ONE real equation in TWO
    # real unknowns (J45_mag, r) -> a line, J45_mag = (J45_star_mag) * r.
    r_sym, J45m2 = sp.symbols("r J45m2", real=True, positive=True)
    params_r = dict(J23=r_sym * p["J23"])
    mu_r = M.moments_difference(0, J45m2, sp.pi, kmax=4, params=params_r)
    m2_r_expr = sp.simplify(mu_r[2])
    sol = sp.solve(sp.Eq(m2_r_expr, 0), J45m2)
    assert len(sol) == 1
    slope = float(sol[0].coeff(r_sym))
    gate["cancellation_curve_slope"] = slope  # J45_mag = slope * r
    gate["cancellation_curve_matches_J45_star"] = bool(abs(slope - J45star_mag) < 1e-9)

    m3_r_expr = sp.simplify(mu_r[3])
    r_vals = np.linspace(0.2, 3.0, 15)
    m3_on_curve = []
    for rv in r_vals:
        m3v = complex(m3_r_expr.subs({r_sym: sp.nsimplify(rv, rational=False), J45m2: sp.nsimplify(rv * slope, rational=False)}))
        m3_on_curve.append(abs(m3v))
    gate["m3_on_curve_min_abs"] = float(min(m3_on_curve))
    gate["m3_nonzero_on_curve"] = bool(min(m3_on_curve) > 1e-4)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    im0 = axes[0].pcolormesh(P, E / J45star_mag, logabs, shading="auto", cmap="viridis")
    axes[0].plot(np.pi, 1.0, "r*", ms=14, label="isolated zero (2 controls)")
    axes[0].set_xlabel(r"$\phi$"); axes[0].set_ylabel(r"$|J_{45}|/|J_{45}^\star|$")
    axes[0].set_title(r"$\log_{10}|m_2(\phi,|J_{45}|)|$: isolated zero")
    axes[0].legend(fontsize=8)
    fig.colorbar(im0, ax=axes[0])

    r_line = np.linspace(0.05, 3.0, 60)
    J45_line = slope * r_line
    axes[1].plot(r_line, J45_line, "r-", lw=2, label=r"$m_2=0$ curve: $|J_{45}|=|J_{45}^\star|\,r$")
    axes[1].plot(1.0, J45star_mag, "ko", ms=8, label="reference tuned point ($r=1$)")
    axes[1].set_xlabel(r"$r$ (branch-1 amplitude scale)")
    axes[1].set_ylabel(r"$|J_{45}|$ (branch-2 amplitude, at $\phi=\pi$)")
    axes[1].set_title("Codimension-1 cancellation curve (3rd control $r$)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figP3_manifold.png"), dpi=150)
    plt.close(fig)

    gate["overall_pass"] = bool(
        gate["isolated_zero_confirmed"]
        and gate["cancellation_curve_matches_J45_star"]
        and gate["m3_nonzero_on_curve"]
    )
    return gate


# ----------------------------------------------------------------------
# P4: universal crossover collapse
# ----------------------------------------------------------------------
def run_P4():
    gate = {}
    J45star = M.J45_star(phi_val=sp.pi)
    # delta=0.1 is kept only as a reported edge case: Gamma_x ~ O(1) there is
    # comparable to the microscopic coupling scale, so O(delta*Gamma^-4,
    # Gamma^-5) corrections are not negligible and the universal two-term
    # law is not expected to hold cleanly (see plan P4 risk notes).
    deltas_universal = [1e-2, 1e-3, 1e-4, 1e-5]
    deltas_all = [1e-1] + deltas_universal

    # predicted a, b from moments: R_S ~ a*delta*Gamma^-3 + b*Gamma^-4
    # a = d(m2)/d(delta) at delta=0 (delta multiplies eta: eta=eta*(1+delta))
    # b = -m3 at the tuned point (Theorem II sign convention: m_k enters with (-1)^k Gamma^-(k+1))
    delta_sym = sp.Symbol("delta_", real=True)
    J45_of_delta = J45star / (1 + delta_sym)  # eta ~ 1/J45, so eta*(1+delta) <-> J45/(1+delta)
    mu_delta = M.moments_difference(0, J45_of_delta, sp.pi, kmax=3)
    m2_of_delta = sp.simplify(mu_delta[2])
    a_pred = complex(sp.diff(m2_of_delta, delta_sym).subs(delta_sym, 0))
    mu_star = M.moments_difference(0, J45star, sp.pi, kmax=4)
    m3_star = complex(sp.simplify(mu_star[3]))
    b_pred = -m3_star
    a_pred = -a_pred  # sign: m2 enters R_S at order (-1)^2 Gamma^-3 = +Gamma^-3, m3 at (-1)^3 Gamma^-4 = -Gamma^-4
    gate["a_predicted"] = [a_pred.real, a_pred.imag]
    gate["b_predicted"] = [b_pred.real, b_pred.imag]
    gamma_x_pred = {str(d): float(abs(b_pred / (a_pred * d))) for d in deltas_all}
    gate["Gamma_cross_predicted"] = gamma_x_pred

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    crossing_gammas = {}
    collapse_curves = {}
    for delta in deltas_all:
        J45v = J45star / (1 + delta)
        gammas, vals, conds = sweep_R_S(J45v, sp.pi, gmin=-2, gmax=7, n_per_decade=30)
        scaled = vals * gammas**4
        x = delta * gammas
        axes[0].plot(x, scaled.real, ".", ms=2, label=f"$\\delta={delta:.0e}$")
        collapse_curves[str(delta)] = (x, scaled.real)

        # crossover extraction: local max of nu_eff within the float64-safe
        # range, then first drop below 3.5 after that max (genuine
        # plateau -> Gamma^-3-tail transition, not a spurious early dip)
        safe = gammas < F64_SAFE_GAMMA
        fit = fit_nu_loglog(gammas[safe], vals[safe])
        nu_eff = fit["nu_eff"]
        gmid = fit["gamma_mid"]
        imax = int(np.argmax(nu_eff))
        after = np.where(nu_eff[imax:] < 3.5)[0]
        gx = float(gmid[imax + after[0]]) if len(after) else None
        crossing_gammas[str(delta)] = {"nu_max": float(nu_eff[imax]), "Gamma_at_max": float(gmid[imax]), "Gamma_x": gx}

    # relative collapse spread over the common scaling region (universal
    # deltas only), evaluated by interpolating each curve onto a shared
    # x = delta*Gamma grid and comparing to the theory line
    x_common = np.logspace(-1, 1, 40)
    interp_vals = []
    for delta in deltas_universal:
        x, y = collapse_curves[str(delta)]
        order = np.argsort(x)
        interp_vals.append(np.interp(x_common, x[order], y[order]))
    interp_vals = np.array(interp_vals)
    theory_line = a_pred.real * x_common + b_pred.real
    spread = np.std(interp_vals, axis=0) / np.maximum(np.abs(theory_line), 1e-12)
    gate["collapse_relative_spread_mean"] = float(np.mean(spread))
    gate["collapse_relative_spread_max"] = float(np.max(spread))

    x_th = np.logspace(-2, 2, 100)
    axes[0].plot(x_th, (a_pred.real * x_th + b_pred.real), "k-", lw=1.5, label="prediction $a\\delta\\Gamma+b$")
    axes[0].set_xscale("symlog")
    axes[0].set_xlabel(r"$\delta\Gamma$")
    axes[0].set_ylabel(r"Re $\Gamma^4 R_S$")
    axes[0].set_title("P4a: crossover collapse")
    axes[0].legend(fontsize=7)

    gx_fit = [crossing_gammas[str(d)]["Gamma_x"] for d in deltas_universal]
    gx_pred = [gamma_x_pred[str(d)] for d in deltas_universal]
    axes[1].loglog(deltas_universal, gx_fit, "o-", label=r"$\Gamma_\times$ (fit)")
    axes[1].loglog(deltas_universal, gx_pred, "s--", label=r"$\Gamma_\times$ (predicted $|b/a\delta|$)")
    axes[1].set_xlabel(r"$|\delta|$"); axes[1].set_ylabel(r"$\Gamma_\times$")
    axes[1].set_title("P4b: crossover scale vs detuning")
    axes[1].legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figP4_collapse.png"), dpi=150)
    plt.close(fig)

    gate["Gamma_cross"] = crossing_gammas
    log_d = np.log10(deltas_universal)
    log_gx = np.log10(gx_fit)
    slope = np.polyfit(log_d, log_gx, 1)[0]
    gate["scaling_slope"] = float(slope)
    gate["gate_slope_near_minus1"] = bool(abs(slope + 1) < 0.15)
    ratio = np.array(gx_fit) / np.array(gx_pred)
    gate["gate_gx_matches_prediction"] = bool(np.all((ratio > 0.5) & (ratio < 2.0)))
    gate["gate_collapse_spread_lt_5pct"] = bool(gate["collapse_relative_spread_mean"] < 0.05)

    gate["overall_pass"] = bool(
        gate["gate_slope_near_minus1"] and gate["gate_gx_matches_prediction"]
    )
    return gate


# ----------------------------------------------------------------------
# P5: finite-window misclassification map
# ----------------------------------------------------------------------
def run_P5():
    gate = {}
    J45star = M.J45_star(phi_val=sp.pi)
    deltas = [1e-2, 1e-3, 1e-4, 1e-5]

    fig, ax = plt.subplots(figsize=(7, 5))
    plateau_decades = {}
    for delta in deltas:
        J45v = J45star / (1 + delta)
        gammas, vals, conds = sweep_R_S(J45v, sp.pi, gmin=0, gmax=11, n_per_decade=20)
        fit = fit_nu_loglog(gammas, vals)
        nu_eff = fit["nu_eff"]
        gmid = fit["gamma_mid"]
        ax.semilogx(gmid, nu_eff, "-", lw=1, label=f"$\\delta={delta:.0e}$")

        mask = np.abs(nu_eff - 4.0) < 0.05
        if np.any(mask):
            idx = np.where(mask)[0]
            plateau_decades[str(delta)] = float(np.log10(gmid[idx[-1]]) - np.log10(gmid[idx[0]]))
        else:
            plateau_decades[str(delta)] = 0.0

    ax.axhline(4, color="gray", ls=":", lw=1)
    ax.axhline(3, color="gray", ls=":", lw=1)
    ax.set_xlabel(r"$\Gamma$"); ax.set_ylabel(r"$\nu_{eff}(\Gamma)$")
    ax.set_title("P5: preasymptotic masquerade")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figP5_masquerade.png"), dpi=150)
    plt.close(fig)

    gate["plateau_decades"] = plateau_decades
    gate["max_plateau_decades"] = float(max(plateau_decades.values()))
    gate["gate_plateau_gt_2_decades"] = bool(max(plateau_decades.values()) > 2.0)
    gate["overall_pass"] = bool(gate["gate_plateau_gt_2_decades"])
    return gate


# ----------------------------------------------------------------------
# P6: exact polynomial certificate
# ----------------------------------------------------------------------
def run_P6():
    gate = {}
    J45star = M.J45_star(phi_val=sp.pi)
    cases = {
        "generic": (sp.Rational(1, 2), sp.Integer(0)),
        "tuned": (J45star, sp.pi),
        "detuned": (J45star * sp.Rational(999, 1000), sp.pi),
    }
    table = {}
    for name, (J45v, phiv) in cases.items():
        nu_cert, Qp, Np = M.certificate_R_S(0, J45v, phiv)
        mu = M.moments_difference(0, J45v, phiv, kmax=4)
        mu = [sp.nsimplify(sp.simplify(m), rational=True) for m in mu]
        nz = [i for i, m in enumerate(mu) if m != 0]
        nu_moment = (nz[0] + 1) if nz else None
        table[name] = {"nu_cert": nu_cert, "nu_moment": nu_moment,
                        "numerator_nonzero": Np is not None}
    gate["table"] = table
    gate["gate_cert_matches_moment"] = bool(all(v["nu_cert"] == v["nu_moment"] for v in table.values()))
    gate["gate_tuned_is_class_II"] = bool(table["tuned"]["numerator_nonzero"])
    gate["overall_pass"] = bool(gate["gate_cert_matches_moment"] and gate["gate_tuned_is_class_II"])
    return gate


# ----------------------------------------------------------------------
# P7: observable robustness
# ----------------------------------------------------------------------
def run_P7():
    gate = {}
    J45star = M.J45_star(phi_val=sp.pi)
    cases = {
        "generic": (sp.Rational(1, 2), sp.Integer(0)),
        "tuned": (J45star, sp.pi),
    }
    readouts = ["rho51", "rho31"]  # rho31 as the deliberately "blind" readout O3
    table = {}
    for readout in readouts:
        for name, (J45v, phiv) in cases.items():
            mu = M.moments_difference(0, J45v, phiv, kmax=4, readout=readout)
            mu = [sp.nsimplify(sp.simplify(m), rational=True) for m in mu]
            nz = [i for i, m in enumerate(mu) if m != 0]
            nu_moment = (nz[0] + 1) if nz else None
            table[f"{readout}_{name}"] = {"moments": [str(m) for m in mu], "nu_moment": nu_moment}
    gate["table"] = table
    gate["O1_shows_promotion"] = bool(
        table["rho51_generic"]["nu_moment"] == 3 and table["rho51_tuned"]["nu_moment"] == 4
    )
    gate["O3_blind_or_different"] = table.get("rho31_tuned", {}).get("nu_moment")
    gate["overall_pass"] = bool(gate["O1_shows_promotion"])
    return gate


# ----------------------------------------------------------------------
# P8: reduced vs full Liouvillian
# ----------------------------------------------------------------------
def run_P8():
    gate = {}
    J45star = M.J45_star(phi_val=sp.pi)
    cases = {
        "generic": (float(sp.Rational(1, 2)), 0.0),
        "tuned": (float(J45star), np.pi),
    }
    case_results = {}
    for name, (J45v, phiv) in cases.items():
        gammas = np.logspace(1, 5, 20)
        R_S_f64, _ = M.lambdify_R_S(0, sp.Float(J45v), sp.Float(phiv))
        errs = []
        for g in gammas:
            r_red, _ = R_S_f64(g)

            rho1_full = M.implicit_linear_response(g, J45v, phiv, cut=False)
            rho1_cut = M.implicit_linear_response(g, J45v, phiv, cut=True)
            r_full = (rho1_full - rho1_cut)[4, 0]  # rho51 element (row=|5>, col=|1>)

            errs.append(abs(r_red - r_full) / max(abs(r_full), 1e-300))

        errs = np.array(errs)
        max_top_err = float(np.max(errs[-6:]))  # top ~1.5 decades of the sweep
        case_results[name] = {
            "gammas": gammas.tolist(),
            "rel_err": errs.tolist(),
            "max_top3decade_err": max_top_err,
        }

    gate["cases"] = case_results
    gate["overall_pass"] = bool(all(v["max_top3decade_err"] < 1e-2 for v in case_results.values()))
    return gate


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------
def main():
    gates = {}
    t0 = time.time()
    print("Running P1 (direct validation of nu 3->4)...")
    gates["P1_promotion"] = run_P1()
    print("Running P2 (symbolic moment cancellation)...")
    gates["P2_symbolic"] = run_P2()
    print("Running P3 (cancellation manifold)...")
    gates["P3_manifold"] = run_P3()
    print("Running P4 (universal crossover collapse)...")
    gates["P4_collapse"] = run_P4()
    print("Running P5 (finite-window misclassification)...")
    gates["P5_masquerade"] = run_P5()
    print("Running P6 (exact polynomial certificate)...")
    gates["P6_certificate"] = run_P6()
    print("Running P7 (observable robustness)...")
    gates["P7_robustness"] = run_P7()
    print("Running P8 (reduced vs full Liouvillian)...")
    gates["P8_reduced_vs_full"] = run_P8()

    gates["phase_p_overall_pass"] = bool(all(
        gates[k]["overall_pass"] for k in gates if isinstance(gates[k], dict) and "overall_pass" in gates[k]
    ))
    gates["runtime_seconds"] = time.time() - t0

    out_path = os.path.join(RESULTS, "gates_summary_phaseP.json")
    with open(out_path, "w") as f:
        json.dump(gates, f, indent=2, default=_default)

    print("\n=== Phase P summary ===")
    for k, v in gates.items():
        if isinstance(v, dict) and "overall_pass" in v:
            print(f"{k}: overall_pass = {v['overall_pass']}")
    print(f"phase_p_overall_pass = {gates['phase_p_overall_pass']}")
    print(f"(runtime {gates['runtime_seconds']:.1f}s, written to {out_path})")


if __name__ == "__main__":
    main()
