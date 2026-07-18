"""Step 4 of the room-temperature SMRT no-go plan (plan Sec. 5, Step 4):
does 300 K sit deep in the merged-manifold asymptotic regime, using
conservative (least-suppressive) Gamma_XY(T) models?

Two DISTINCT scales matter here, and conflating them was the open flag
left by Step 3 (RoomT/README.md):

  x1(T) = Gamma(T) / d_strain -- the Hamiltonian-level "branches merged"
      scale (plan's x(T) = Gamma_XY/Delta_pair). Gate 4 as literally
      stated wants x1(300K) >> 1.

  Gamma_cross ~ beta/geff -- the SEPARATE denominator-dominance scale
      found in Step 3, below which R_EIT's exact leading behavior is the
      PRE-crossover Gamma^-3 asymptote (not the Gamma^-4 one), because the
      control-power-broadening term beta*S2 in the denominator still
      dominates the Gamma-independent ground decoherence geff.

Both exact leading coefficients come from Step 3's symbolic certificate
(reused directly, not rederived):

  R_asym_3(Gamma) = -2*pi^2*(Dperp^2+Lperp^2) / Gamma^3   (Gamma << Gamma_cross)
  R_asym_4(Gamma) = -2*pi^2*beta*(Dperp^2+Lperp^2)/geff / Gamma^4   (Gamma >> Gamma_cross)

so this script checks x1(300K) >> 1 (manifold IS Hamiltonian-merged) AND
reports which of the two asymptotes actually matches the exact numeric
R_EIT(300K) -- resolving, rather than assuming past, the Step 3 flag.

Gamma_XY(T) models: reuses phonon_rates.k_orb_variants (4 models: full
Happacher, literature-uncertainty conservative lower bound, saturation,
naive T^5 extrapolation) and phonon_rates.k_orb_most_conservative (picks
the smallest dissipation in [4,300] K, plan Sec. 4.1's "no-goに最も不利な
...モデル"), guarded by tests/test_phonon_rate_variants.py.

Outputs: RoomT/results/gates_summary_step4.json,
         RoomT/results/figures/fig_step4_temperature_scaling.png
"""
from __future__ import annotations
import os
import sys
import json

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOMT_ROOT = os.path.dirname(HERE)
NOGO_SRC = os.path.join(ROOMT_ROOT, "..", "..", "No-go theorem", "src")
sys.path.insert(0, NOGO_SRC)
sys.path.insert(0, HERE)

import nv_model as nm  # noqa: E402
import phonon_rates as pr  # noqa: E402
from step3_merged_manifold_moments import build_symbolic_certificate  # noqa: E402

RESULTS_DIR = os.path.join(ROOMT_ROOT, "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

D_STRAIN = 1.683
BVEC = (0.0, 0.0, 0.02)
REFERENCE_OC = 1.0
REFERENCE_GG = 6.3e-5

CAMPAIGN_T = [4.0, 10.0, 30.0, 50.0, 77.0, 150.0, 300.0]  # plan Sec. 4.1 focus points
MODELS = ["full_happacher", "conservative_lower_bound", "saturation", "naive_T5_extrapolation"]
MAIN_MODEL = "conservative_lower_bound"  # picked by k_orb_most_conservative in [4,300] K


def gamma_oc_variant(T, d, label):
    """Same combination nv_model.gamma_oc_GHz uses (0.5*k_orb + 0.5*GRAD),
    with k_orb replaced by the requested Gamma_XY(T) model variant."""
    k_orb_variant_Hz = pr.k_orb_variants(T, d)[label]
    return 0.5 * k_orb_variant_Hz * 1e-9 + 0.5 * nm.GRAD


def asymptotic_coefficients(Ep):
    """Exact leading coefficients of both R_EIT asymptotes (pre- and
    post-denominator-crossover), reusing Step 3's certificate -- rebuilt
    at z=Ep (the actual probe resonance, matching nv_model.probe_line())
    rather than Step 3's default z=0, so the numeric comparison below is
    physically meaningful (on-resonance), not at an arbitrary detuning."""
    cert = build_symbolic_certificate(z_val=Ep)
    Qp, N12, N21, N22 = cert["Qp"], cert["N12"], cert["N21"], cert["N22"]
    Dperp_val = sp.Rational(155, 200)            # 1.55/2
    Lperp_val = sp.sqrt(2) * sp.Rational(1, 10)   # 0.20/sqrt(2)
    subs = {cert["Dperp"]: Dperp_val, cert["Lperp"]: Lperp_val}

    coeff3 = complex(sp.N((-N12.LC() * N21.LC() / (Qp.LC() * N22.LC())).subs(subs), 30))
    beta_s, geff_s = sp.symbols("beta geff", positive=True)
    coeff4_expr = -beta_s * N12.LC() * N21.LC() / (geff_s * Qp.LC() ** 2)
    coeff4_fn = sp.lambdify((beta_s, geff_s), coeff4_expr.subs(subs), "numpy")
    return coeff3, coeff4_fn


def run():
    Bvec = BVEC
    w, U = nm.dressed_ground(Bvec)
    H = nm.Hes(Bvec, d=D_STRAIN)
    dp, dc = nm.legs(U)
    Ep, fid, _ = nm.probe_line(H, dp)  # on-resonance point, matches Step 1's convention

    coeff3, coeff4_fn = asymptotic_coefficients(Ep)
    beta = (nm.TWOPI * REFERENCE_OC) ** 2 / 4
    geff = 2 * REFERENCE_GG + 2e-6
    coeff4 = complex(coeff4_fn(beta, geff))

    rows = []
    for T in CAMPAIGN_T:
        row = dict(T=T)
        for label in MODELS:
            g = gamma_oc_variant(T, D_STRAIN, label)
            r = nm.response(H, dp, dc, Ep, g, Oc=REFERENCE_OC, gg=REFERENCE_GG)
            # nv_model.response() doesn't return dXi directly; recompute
            # from its own K12,K21,den fields (same formula, Step 3's sign).
            dXi = -beta * r["K12"] * r["K21"] / r["den"]
            row[f"Gamma_{label}"] = g
            row[f"C_EIT_{label}"] = r["C"]
            row[f"R_EIT_{label}"] = complex(dXi)
        rows.append(row)

    # main-model (conservative_lower_bound) asymptotic comparison at 300 K
    g300 = gamma_oc_variant(300.0, D_STRAIN, MAIN_MODEL)
    x1_300 = g300 / D_STRAIN
    R_num_300 = rows[-1][f"R_EIT_{MAIN_MODEL}"]
    R_asym3_300 = coeff3 / g300 ** 3
    R_asym4_300 = coeff4 / g300 ** 4
    err3 = abs(R_num_300 - R_asym3_300) / abs(R_num_300)
    err4 = abs(R_num_300 - R_asym4_300) / abs(R_num_300)
    gamma_cross = abs(beta / geff)

    # C_EIT(T) monotonicity over the campaign grid, main model
    C_main = [abs(row[f"C_EIT_{MAIN_MODEL}"]) for row in rows]
    monotone_decreasing = all(C_main[i + 1] <= C_main[i] * 1.05 for i in range(len(C_main) - 1))

    # Oc scan at 300 K: does lowering Oc push 300K into the true Gamma^-4
    # asymptote (smaller Gamma_cross ~ Oc^2/geff), and where is |C_EIT| max?
    Oc_grid = np.logspace(-3, 1, 25)  # 1 MHz .. 10 GHz, plan Sec. 4.2 range (+stress test)
    Oc_scan = []
    for Oc in Oc_grid:
        r = nm.response(H, dp, dc, Ep, g300, Oc=Oc, gg=REFERENCE_GG)
        beta_o = (nm.TWOPI * Oc) ** 2 / 4
        gcross_o = abs(beta_o / geff)
        Oc_scan.append(dict(Oc=float(Oc), C_EIT=float(abs(r["C"])), gamma_cross=float(gcross_o)))
    best = max(Oc_scan, key=lambda d: d["C_EIT"])

    gates = dict(
        x1_300K_much_greater_than_1=bool(x1_300 > 10),
        C_EIT_monotone_decreasing_main_model=bool(monotone_decreasing),
        crossover_diagnosis_self_consistent=bool(min(err3, err4) < 0.3),
        conservative_model_still_tiny_at_300K=bool(C_main[-1] < 1e-3),
    )
    gates["overall_pass"] = bool(all(gates.values()))

    which_regime = "pre-crossover (Gamma^-3)" if err3 < err4 else "post-crossover (Gamma^-4)"

    summary = dict(
        model="No-go theorem/src/nv_model.py, Step 1-3 geometry; Step 3 exact asymptotic coefficients",
        Oc_GHz=REFERENCE_OC, gg_GHz=REFERENCE_GG, d_strain_GHz=D_STRAIN,
        gamma_cross_GHz_at_reference_Oc_gg=gamma_cross,
        campaign_T_K=CAMPAIGN_T,
        rows=[{k: (v if not isinstance(v, complex) else [v.real, v.imag]) for k, v in row.items()}
              for row in rows],
        at_300K=dict(
            main_model=MAIN_MODEL,
            Gamma_GHz=g300,
            x1_Gamma_over_d_strain=x1_300,
            R_num=[R_num_300.real, R_num_300.imag],
            R_asym_Gamma3=[R_asym3_300.real, R_asym3_300.imag],
            R_asym_Gamma4=[R_asym4_300.real, R_asym4_300.imag],
            rel_err_vs_Gamma3_asymptote=err3,
            rel_err_vs_Gamma4_asymptote=err4,
            regime_300K_actually_in=which_regime,
        ),
        Oc_scan_at_300K=Oc_scan,
        best_Oc_at_300K=best,
        gates=gates,
    )

    out_json = os.path.join(RESULTS_DIR, "gates_summary_step4.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # --- figure ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    for label in MODELS:
        Ts = CAMPAIGN_T
        Cs = [abs(row[f"C_EIT_{label}"]) for row in rows]
        ax.semilogy(Ts, Cs, "o-", label=label)
    ax.set_xlabel("T (K)"); ax.set_ylabel("|C_EIT|")
    ax.set_title("C_EIT(T), Oc=1 GHz, gg=6.3e-5 GHz"); ax.legend(fontsize=7)

    ax = axes[1]
    Oc_arr = np.array([d["Oc"] for d in Oc_scan])
    C_arr = np.array([d["C_EIT"] for d in Oc_scan])
    gc_arr = np.array([d["gamma_cross"] for d in Oc_scan])
    ax.loglog(Oc_arr, C_arr, "o-", label="|C_EIT|(300K)")
    ax2 = ax.twinx()
    ax2.loglog(Oc_arr, gc_arr, "r--", label="Gamma_cross(Oc)")
    ax2.axhline(g300, color="gray", linestyle=":", label="Gamma(300K)")
    ax.set_xlabel("Omega_c (GHz)"); ax.set_ylabel("|C_EIT|(300K)")
    ax2.set_ylabel("Gamma_cross (GHz)")
    ax.set_title("300K: control-power scan")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=6, loc="best")
    fig.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_step4_temperature_scaling.png")
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)

    print(json.dumps(gates, indent=2))
    print(f"300K regime: {which_regime} (err3={err3:.3f}, err4={err4:.3f})")
    print(f"x1(300K) = {x1_300:.2f}")
    print(f"best Oc at 300K: {best}")
    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    return summary


if __name__ == "__main__":
    run()
