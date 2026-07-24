"""run_gate_c.py -- Gate C driver: material independence (PRL Priorities 4 & 5).

Two new results plus a material-independence overlay:
  P4  group-IV SiV-/SnV- Gamma^-1 class certified in a FULL physical GKSL
      (group_iv_full.py), matched to the reduced kernel to machine precision.
  P5  three physical suppression classes n = 1, 2, 3 collapse under Gamma^n:
        n=1  group-IV orbital-Lambda (M0 != 0)
        n=2  NV ms=0<->-1        (graph distance d=1)
        n=3  NV ms=-1<->+1       (graph distance d=2)
  Material independence: the SAME integer class appears in diamond AND in the
  non-diamond superconducting witness of Gate B --
        class 1: group-IV (diamond)  &  SC generic (non-diamond)
        class 2: NV (diamond)        &  SC protected (non-diamond)
  The plateau VALUE (the first path moment) is material-specific; the integer
  EXPONENT is shared -- that is the universality claim.

Gates (strategy doc Sec.13 P4/P5, Sec.14 Gate C):
  G-C1  group-IV full-GKSL slope -1; reduced==full to <1e-7 (dephasing);
        M0 != 0; Gamma*R -> M0; hopping mode also slope -1.
  G-C2  three physical classes: first-nonzero-moment (graph-distance)
        predicted exponent matches the log-log slope for n=1,2,3.
  G-C3  Gamma^n |K| collapses to a plateau (rel spread < 5%), exponent stable
        across fit sub-windows; finite-Gamma correction reported.
  G-C4  material independence: class 1 slope=1 in group-IV & SC-generic;
        class 2 slope=2 in NV & SC-protected.
  G-C5  full/reduced agree (G-C1) and pre-asymptotic vs asymptotic separated.

Usage:  python run_gate_c.py [--quick] [--smoke]
Outputs: results/tables/gates_summary_gateC.json, gate_c_collapse.csv,
         results/figures/fig_gateC_three_class_collapse.png (+ .pdf)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE_SRC = os.path.join(HERE, "..", "..", "src")
NOGO_SRC = os.path.join(HERE, "..", "..", "..", "No-go theorem", "src")
GATEB_SRC = os.path.join(HERE, "..", "..", "GateB_superconducting_witness", "src")
for _p in (HERE, PHASE_SRC, NOGO_SRC, GATEB_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import core                              # noqa: E402
import group_iv_full as gf               # noqa: E402
import group_iv_model as giv             # noqa: E402
import nv_reduced_kernel as nvk          # noqa: E402
import model_sc_transfer as sc           # noqa: E402

RESULTS = os.path.join(HERE, "..", "results")
TABLES = os.path.join(RESULTS, "tables")
FIGS = os.path.join(RESULTS, "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

_NV_H = nvk.H_3E()


# ---- kernel callables (all return |K(Gamma)| amplitude) -----------------
def k_groupiv_full(g, material="SiV", mode="dephasing"):
    return gf.full_response(g, material, mode=mode)


def k_groupiv_reduced(g, material="SiV"):
    return gf.reduced_kernel_response(g, material)


def k_nv(g, pair):
    return nvk.kernel(_NV_H, pair, np.array([g]))[0]


def k_sc(g, tuning):
    return sc.transfer_kernel(g, tuning=tuning)


# ---- helpers ------------------------------------------------------------
def slope(gammas, vals):
    return core.fit_nu_loglog(np.asarray(gammas, float), np.asarray(vals))["nu_global"]


def compensated_spread(gammas, vals, n):
    """rel spread of Gamma^n |K| over the window (plateau flatness)."""
    y = np.asarray(gammas, float) ** n * np.abs(vals)
    return float(np.std(y) / np.mean(y)), float(np.mean(y))


def window_stability(gammas, vals):
    """slope on lower vs upper half of the window (fit-window stability)."""
    g = np.asarray(gammas, float); v = np.asarray(vals)
    h = len(g) // 2
    return slope(g[:h], v[:h]), slope(g[h:], v[h:])


# ---- P4: group-IV full GKSL --------------------------------------------
def run_p4(quick=False):
    rows = {}
    npts = 18 if quick else 30
    for material, (lo, hi) in (("SiV", (4, 8)), ("SnV", (5, 10))):
        gammas = np.logspace(lo, hi, npts)
        Kf = np.array([k_groupiv_full(g, material, "dephasing") for g in gammas])
        Kh = np.array([k_groupiv_full(g, material, "hopping") for g in gammas])
        Kr = np.array([k_groupiv_reduced(g, material) for g in gammas])
        m0 = gf.M0(material)
        # reduced==full up to the -i drive factor
        red_full_err = float(max(abs(Kf[i] / (-1j) - Kr[i]) / abs(Kr[i]) for i in range(len(gammas))))
        gamma_R = complex(gammas[-1] * Kf[-1] / (-1j))  # -> M0
        rows[material] = dict(
            slope_dephasing=slope(gammas, Kf),
            slope_hopping=slope(gammas, Kh),
            reduced_vs_full_max_rel_err=red_full_err,
            abs_M0=float(abs(m0)),
            gamma_times_R_over_M0=float(abs(gamma_R) / abs(m0)),
            rho0_residual=float(gf.verify_rho0_steady(gammas[len(gammas) // 2], material)),
        )
    return rows


# ---- P5 + material independence: three classes -------------------------
def run_p5(quick=False):
    npts = 20 if quick else 32
    # (class n, systems) ; each system: (label, material_type, kernel fn, window)
    classes = [
        (1, [
            ("group-IV SiV", "diamond", lambda g: k_groupiv_full(g, "SiV"), (4, 8)),
            ("SC generic", "non-diamond", lambda g: k_sc(g, "generic"), (5, 9)),
        ]),
        (2, [
            ("NV 0<->-1", "diamond", lambda g: k_nv(g, (0, -1)), (2, 5)),
            ("SC protected", "non-diamond", lambda g: k_sc(g, "protected"), (5, 9)),
        ]),
        (3, [
            ("NV -1<->+1", "diamond", lambda g: k_nv(g, (-1, 1)), (2, 5)),
        ]),
    ]
    out = []
    curves = []
    for n, systems in classes:
        for label, mtype, kfn, (lo, hi) in systems:
            gammas = np.logspace(lo, hi, npts)
            vals = np.array([kfn(g) for g in gammas])
            # predicted exponent = first nonzero moment index + 1 (graph distance)
            sl = slope(gammas, vals)
            sp, plateau = compensated_spread(gammas, vals, n)
            slo, shi = window_stability(gammas, vals)
            out.append(dict(
                n_class=n, system=label, material_type=mtype,
                slope=sl, predicted=n, slope_match=bool(abs(sl - n) < 0.05),
                compensated_spread=sp, plateau_value=plateau,
                collapse_ok=bool(sp < 0.05),
                slope_lo=slo, slope_hi=shi,
                window_stable=bool(abs(slo - shi) < 0.05),
            ))
            curves.append((n, label, mtype, gammas, np.abs(vals)))
    return out, curves


# ---- predicted exponents from moments (graph-distance predictor) -------
def predicted_exponents():
    giv_m0 = abs(gf.M0("SiV"))
    nv_pred = {}
    for pair in [(0, -1), (-1, 1)]:
        M = nvk.moments(_NV_H, pair, 3)
        nz = [i for i, m in enumerate(M) if abs(m) > 1e-8]
        nv_pred[str(pair)] = (nz[0] + 1) if nz else None
    return dict(group_IV=(1 if giv_m0 > 1e-8 else None),
                NV_0_m1=nv_pred["(0, -1)"], NV_m1_p1=nv_pred["(-1, 1)"])


# ---- figure -------------------------------------------------------------
def make_figure(curves, path_png):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.7))
    colors = {1: "C0", 2: "C1", 3: "C2"}
    styles = {"diamond": "o-", "non-diamond": "s--"}
    for n, label, mtype, g, y in curves:
        ax1.loglog(g, y, styles[mtype], ms=3, lw=1.0, color=colors[n],
                   label=f"{label} (n={n})")
        ax2.loglog(g, g ** n * y, styles[mtype], ms=3, lw=1.0, color=colors[n],
                   label=f"{label}")
    ax1.set_xlabel(r"$\Gamma$"); ax1.set_ylabel(r"$|K(\Gamma)|$")
    ax1.set_title("Gate C: three suppression classes")
    ax1.legend(fontsize=7)
    ax2.set_xlabel(r"$\Gamma$"); ax2.set_ylabel(r"$\Gamma^{n}\,|K(\Gamma)|$ (compensated)")
    ax2.set_title("scaling collapse to plateaux (diamond o / non-diamond s)")
    ax2.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path_png, dpi=200)
    fig.savefig(path_png.replace(".png", ".pdf"))
    plt.close(fig)


# ---- main ---------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Gate C: material independence")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    quick = args.quick or args.smoke

    t0 = time.time()
    p4 = run_p4(quick=quick)
    pred = predicted_exponents()
    p5, curves = run_p5(quick=quick)

    # -- gates --------------------------------------------------------
    g_c1 = bool(all(
        abs(r["slope_dephasing"] - 1) < 1e-2 and abs(r["slope_hopping"] - 1) < 1e-2
        and r["reduced_vs_full_max_rel_err"] < 1e-7 and r["abs_M0"] > 1e-8
        and abs(r["gamma_times_R_over_M0"] - 1) < 1e-3 and r["rho0_residual"] < 1e-9
        for r in p4.values()))

    g_c2 = bool(pred["group_IV"] == 1 and pred["NV_0_m1"] == 2 and pred["NV_m1_p1"] == 3
                and all(r["slope_match"] for r in p5))

    g_c3 = bool(all(r["collapse_ok"] and r["window_stable"] for r in p5))

    diamond1 = [r for r in p5 if r["n_class"] == 1 and r["material_type"] == "diamond"]
    nond1 = [r for r in p5 if r["n_class"] == 1 and r["material_type"] == "non-diamond"]
    diamond2 = [r for r in p5 if r["n_class"] == 2 and r["material_type"] == "diamond"]
    nond2 = [r for r in p5 if r["n_class"] == 2 and r["material_type"] == "non-diamond"]
    g_c4 = bool(
        diamond1 and nond1 and diamond2 and nond2
        and all(abs(r["slope"] - 1) < 0.05 for r in diamond1 + nond1)
        and all(abs(r["slope"] - 2) < 0.05 for r in diamond2 + nond2)
    )

    g_c5 = bool(all(r["reduced_vs_full_max_rel_err"] < 1e-7 for r in p4.values())
                and all(r["window_stable"] for r in p5))

    gates = dict(
        G_C1_groupIV_full_gksl=g_c1,
        G_C2_three_physical_classes=g_c2,
        G_C3_scaling_collapse=g_c3,
        G_C4_material_independence=g_c4,
        G_C5_full_reduced_and_finite_gamma=g_c5,
    )
    gates["overall_pass"] = bool(all(gates.values()))

    summary = dict(
        description="Gate C: material independence -- group-IV full GKSL (P4) and "
                    "three-class scaling collapse n=1,2,3 (P5), diamond & non-diamond.",
        p4_group_iv_full_gksl=p4,
        predicted_exponents=pred,
        p5_three_class=p5,
        quick=quick,
        gates=gates,
        runtime_s=round(time.time() - t0, 2),
    )

    out_json = os.path.join(TABLES, "gates_summary_gateC.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    out_csv = os.path.join(TABLES, "gate_c_collapse.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["n_class", "system", "material_type", "slope",
                                          "predicted", "compensated_spread", "plateau_value",
                                          "window_stable"])
        w.writeheader()
        for r in p5:
            w.writerow({k: r[k] for k in w.fieldnames})

    fig_png = os.path.join(FIGS, "fig_gateC_three_class_collapse.png")
    make_figure(curves, fig_png)

    print(json.dumps(gates, indent=2))
    print("\nP4 group-IV full GKSL:")
    for mat, r in p4.items():
        print(f"  {mat}: slope_deph={r['slope_dephasing']:.3f} slope_hop={r['slope_hopping']:.3f} "
              f"reduced==full={r['reduced_vs_full_max_rel_err']:.1e} |M0|={r['abs_M0']:.3f}")
    print(f"predicted exponents: {pred}")
    print("P5 three classes (diamond & non-diamond):")
    for r in p5:
        print(f"  n={r['n_class']} {r['system']:16s} [{r['material_type']:11s}] "
              f"slope={r['slope']:.3f} collapse_spread={r['compensated_spread']:.1e} "
              f"stable={r['window_stable']}")
    print(f"\nwrote {out_json}\nwrote {out_csv}\nwrote {fig_png}")
    print(f"runtime {summary['runtime_s']} s")
    return summary


if __name__ == "__main__":
    main()
