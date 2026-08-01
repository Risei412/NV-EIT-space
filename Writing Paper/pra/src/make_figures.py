"""make_figures.py -- the six PRA main-text figures.

Presentation only: every panel reads a recorded result of the P1-P6 campaigns
or calls their functions unchanged.  No new physics is introduced here.

  Fig. 1  NV level scheme and the temperature-dependent orbital mixing path
  Fig. 2  Representative spectra: low-T transparency to high-T absorption
  Fig. 3  Full-Liouvillian T-B_perp phase diagram (P1)
  Fig. 4  Signed contrast vs temperature with the uncertainty band (P1/P2/P6)
  Fig. 5  Transverse-field opening, C ~ B_perp^2 (P5)
  Fig. 6  Observables and the detectable region (P3)

Style is inherited from the PRL figure package so the two papers look like
one body of work.

Usage: python make_figures.py [--quick]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRA = HERE.parent
REPO = PRA.parents[1]
sys.path.insert(0, str(REPO / "No-go theorem" / "src"))
sys.path.insert(0, str(REPO / "Writing Paper" / "prl_figures" / "src"))
sys.path.insert(0, str(HERE))

import matplotlib                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

import prl_style as st                 # noqa: E402

TAB = PRA / "results" / "tables"
FIG = PRA / "results" / "figures"
C = st.COLORS


def _save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=300)
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


def _load_json(p):
    return json.loads((TAB / p).read_text())


def _load_csv(p):
    with (TAB / p).open() as fh:
        return list(csv.DictReader(fh))


def _f(x, default=np.nan):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------- Fig. 1
def fig1_level_scheme():
    fig, axes = plt.subplots(1, 2, figsize=(st.COL_2, 2.7))

    ax = axes[0]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    # ground triplet
    for y, lab in [(1.6, r"$m_s=-1$"), (1.0, r"$m_s=0$"), (2.2, r"$m_s=+1$")]:
        ax.plot([1.2, 3.6], [y, y], color=C["black"], lw=1.6)
        ax.text(3.75, y, lab, va="center", fontsize=7.5)
    ax.text(0.9, 1.6, r"$^3A_2$", ha="right", va="center", fontsize=8.5)
    # excited orbital branches
    for y, lab, col in [(7.2, r"$E_x$", C["blue"]), (8.2, r"$E_y$", C["vermillion"])]:
        ax.plot([1.2, 3.6], [y, y], color=col, lw=1.6)
        ax.text(3.75, y, lab, va="center", fontsize=7.5, color=col)
    ax.text(0.9, 7.7, r"$^3E$", ha="right", va="center", fontsize=8.5)
    # probe / control legs
    ax.annotate("", xy=(1.9, 7.15), xytext=(1.9, 1.05),
                arrowprops=dict(arrowstyle="<->", color=C["green"], lw=1.3))
    ax.text(1.55, 4.1, "probe", rotation=90, va="center", ha="center",
            fontsize=7.5, color=C["green"])
    ax.annotate("", xy=(2.9, 7.15), xytext=(2.9, 2.25),
                arrowprops=dict(arrowstyle="<->", color=C["purple"], lw=1.3))
    ax.text(3.25, 4.7, "control", rotation=90, va="center", ha="center",
            fontsize=7.5, color=C["purple"])
    # orbital hopping between branches
    ax.add_patch(FancyArrowPatch((2.4, 7.25), (2.4, 8.15),
                                 arrowstyle="<->", mutation_scale=9,
                                 color=C["orange"], lw=1.4))
    ax.text(2.65, 7.7, r"$\Gamma_{XY}(T)$", fontsize=7.5, color=C["orange"])
    # transverse field mixing in the ground manifold
    ax.add_patch(FancyArrowPatch((1.35, 1.05), (1.35, 2.2),
                                 arrowstyle="<->", mutation_scale=8,
                                 color=C["sky"], lw=1.2, linestyle="--"))
    ax.text(0.35, 1.65, r"$B_\perp$", fontsize=7.5, color=C["sky"])
    st.panel_label(ax, "a", x=0.02, y=0.94)
    ax.set_title("spin-$\\Lambda$ channel and the paths that close it",
                 fontsize=8.5)

    # optical-coherence damping vs temperature
    ax = axes[1]
    import nv_model as nv
    import run_prl_prediction as rp
    Ts = np.linspace(4, 300, 300)
    g = np.array([nv.gamma_oc_GHz(float(t), rp.D) for t in Ts])
    ax.semilogy(Ts, g, color=C["orange"], lw=1.5)
    ax.axhline(rp.D, color=C["gray"], ls=":", lw=1.0)
    ax.text(298, rp.D * 1.35, r"$\delta_{\rm strain}$", ha="right",
            fontsize=7, color=C["gray"])
    for T0, lab in [(25, "island\nopens"), (95, "island\ncloses"),
                    (103, None)]:
        if lab:
            ax.axvline(T0, color=C["gray"], ls="--", lw=0.7)
            ax.text(T0 + 4, 3e-4, lab, fontsize=6.5, color=C["gray"])
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel(r"optical-coherence damping $\gamma_{\rm oc}$ (GHz)")
    ax.set_xlim(0, 300)
    st.panel_label(ax, "b")
    ax.set_title("phonon-driven damping of the Raman path", fontsize=8.5)

    fig.tight_layout()
    _save(fig, "fig1_level_scheme")


# --------------------------------------------------------------- Fig. 2
def fig2_spectra(quick=False):
    import p1_phase_diagram as p1
    import run_prl_prediction as rp
    fig, axes = plt.subplots(1, 3, figsize=(st.COL_2, 2.4), sharey=False)
    shown = [(30.0, "transparency"), (70.0, "transparency"),
             (105.0, "control-induced absorption")]
    for ax, (T, tag) in zip(axes, shown):
        d, Af, Ac, Cc, info = p1.adaptive_spectrum(T, rp.BX0)
        ax.plot(d, Ac, color=C["gray"], lw=1.1, ls="--",
                label=r"pathway cut, $A_{\rm cut}$")
        ax.plot(d, Af, color=C["blue"], lw=1.4, label=r"full, $A_{\rm full}$")
        ax.set_xlabel(r"two-photon detuning $\delta_2$ (MHz)")
        ipk = int(np.argmax(np.abs(Cc)))
        ax.set_title(f"$T$ = {T:.0f} K\n$C$ = {Cc[ipk]:+.3g}", fontsize=8)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    axes[0].set_ylabel("absorption (arb.)")
    axes[0].legend(loc="best", fontsize=6.5)
    for i, ax in enumerate(axes):
        st.panel_label(ax, "abc"[i])
    fig.tight_layout()
    _save(fig, "fig2_spectra")


# --------------------------------------------------------------- Fig. 3
def fig3_phase_diagram():
    z = np.load(TAB / "p1_phase_diagram.npz", allow_pickle=True)
    T, B, K, Cg = z["T_K"], z["Bx_T"], z["klass"], z["C"]
    code = json.loads(str(z["class_code"]))
    inv = {v: k for k, v in code.items()}

    fig, axes = plt.subplots(1, 2, figsize=(st.COL_2, 3.0))

    # (a) classification
    ax = axes[0]
    order = ["transparency", "ATS", "absorption", "unresolved"]
    cols = [C["blue"], C["orange"], C["vermillion"], "#DDDDDD"]
    remap = np.full(K.shape, 3)
    for i, name in enumerate(order):
        remap[K == code[name]] = i
    cmap = ListedColormap(cols)
    norm = BoundaryNorm(np.arange(-0.5, 4.5, 1), cmap.N)
    dB = (B[1] - B[0]) if len(B) > 1 else 0.05
    dT = (T[1] - T[0]) if len(T) > 1 else 5.0
    ax.pcolormesh(np.append(B, B[-1] + dB) - dB / 2,
                  np.append(T, T[-1] + dT) - dT / 2,
                  remap, cmap=cmap, norm=norm, shading="flat")
    ax.set_xlabel(r"transverse field $B_\perp$ (T)")
    ax.set_ylabel("temperature (K)")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="none") for c in cols]
    labels = ["genuine transparency", "Autler-Townes",
              "control-induced absorption", "unresolved"]
    ax.legend(handles, labels, fontsize=6, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=2)
    st.panel_label(ax, "a")

    # (b) signed contrast, log magnitude, sign by colour
    ax = axes[1]
    mag = np.abs(Cg)
    mag[~np.isfinite(mag)] = np.nan
    pos = np.where(Cg > 0, mag, np.nan)
    im = ax.pcolormesh(np.append(B, B[-1] + dB) - dB / 2,
                       np.append(T, T[-1] + dT) - dT / 2,
                       pos, cmap="viridis",
                       norm=LogNorm(vmin=1e-6, vmax=1.0), shading="flat")
    neg = np.isfinite(Cg) & (Cg < 0)
    Bg, Tg = np.meshgrid(B, T)
    ax.plot(Bg[neg], Tg[neg], ".", color=C["vermillion"], ms=2.4,
            label=r"$C<0$")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label(r"$|C|$  (where $C>0$)", fontsize=7.5)
    ax.set_xlabel(r"transverse field $B_\perp$ (T)")
    ax.set_ylabel("temperature (K)")
    ax.legend(fontsize=6.5, loc="lower left")
    st.panel_label(ax, "b")

    fig.tight_layout()
    _save(fig, "fig3_phase_diagram")


# --------------------------------------------------------------- Fig. 4
def fig4_contrast_vs_T():
    z = np.load(TAB / "p1_phase_diagram.npz", allow_pickle=True)
    T, B, Cg = z["T_K"], z["Bx_T"], z["C"]
    import run_prl_prediction as rp
    j = int(np.argmin(np.abs(B - rp.BX0)))
    Cc = Cg[:, j]

    p2 = _load_json("p2_summary.json")
    key = f"{rp.BX0:.5f}"
    b1 = p2["T_1pct_band_vs_field"].get(key, {})
    bs = p2["T_sign_band_vs_field"].get(key, {})
    try:
        p6 = _load_json("p6_summary.json")
        f1 = p6["comparison"]["T_1pct"]["full"]
        fs = p6["comparison"]["T_sign"]["full"]
    except Exception:
        f1 = fs = None

    fig, ax = plt.subplots(figsize=(st.COL_1 * 1.5, 2.9))
    m = np.isfinite(Cc)

    def runs(mask):
        """Contiguous index runs of a boolean mask.

        The sign-reversed points are not contiguous -- there is a separate
        low-temperature lobe as well as the high-temperature reversal -- so
        joining them with one line would draw a spurious segment straight
        through the transparency region.
        """
        out, cur = [], []
        for i, v in enumerate(mask):
            if v:
                cur.append(i)
            elif cur:
                out.append(cur); cur = []
        if cur:
            out.append(cur)
        return out

    for k, idx in enumerate(runs(m & (Cc > 0))):
        ax.semilogy(T[idx], np.abs(Cc[idx]), "o-", color=C["blue"], ms=3.4,
                    lw=1.3, label=r"$C>0$ (transparency)" if k == 0 else None)
    for k, idx in enumerate(runs(m & (Cc < 0))):
        ax.semilogy(T[idx], np.abs(Cc[idx]), "s--", color=C["vermillion"],
                    ms=3.4, lw=1.3,
                    label=r"$|C|$, $C<0$ (absorption)" if k == 0 else None)

    if b1.get("q16") is not None:
        ax.axvspan(b1["q16"], b1["q84"], color=C["green"], alpha=0.16, lw=0)
        ax.text(b1["median"], 3.2e-7, r"$T_{1\%}$", fontsize=7, ha="center",
                color=C["green"])
    if bs.get("q16") is not None:
        ax.axvspan(bs["q16"], bs["q84"], color=C["vermillion"], alpha=0.14, lw=0)
        ax.text(bs["median"], 3.2e-7, r"$T_{\rm sign}$", fontsize=7,
                ha="center", color=C["vermillion"])
    for band, col, dy in ((f1, C["green"], 0.55), (fs, C["vermillion"], 0.30)):
        if band and band.get("median") is not None:
            ax.errorbar([band["median"]], [10 ** (-6 + dy)],
                        xerr=[[band["median"] - band["q16"]],
                              [band["q84"] - band["median"]]],
                        fmt="D", color=col, ms=3.2, lw=1.1, capsize=2)
    ax.text(0.02, 0.04, "diamonds: full-Liouvillian bands (P6)\n"
                        "shaded: reduced-model bands (P2)",
            transform=ax.transAxes, fontsize=6, va="bottom")

    ax.axhline(1.534e-7, color=C["gray"], ls=":", lw=0.9)
    ax.text(T[-1], 1.9e-7, "detection floor", fontsize=6.3, ha="right",
            color=C["gray"])
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel(r"sector contrast $|C|$")
    ax.set_ylim(1e-7, 12)
    ax.legend(fontsize=6.8, loc="upper right", ncol=1)
    ax.set_title(rf"$B_\perp$ = {rp.BX0:.3f} T, full Liouvillian", fontsize=8.5)
    fig.tight_layout()
    _save(fig, "fig4_contrast_vs_T")


# --------------------------------------------------------------- Fig. 5
def fig5_bperp_scaling():
    rows = _load_csv("p5_bperp_scaling.csv")
    summ = _load_json("p5_summary.json")
    B_pert = summ["physical_window"]["B_pert_T"]
    B_mix = summ["physical_window"]["B_full_mixing_T"]
    fits = summ["temperature_fits"]
    oc0 = 0.1
    Ts = sorted({_f(r["T_K"]) for r in rows if _f(r["Oc_GHz"]) == oc0})

    fig, axes = plt.subplots(1, 2, figsize=(st.COL_2, 2.9))

    # (a) the data and the fitted curves
    ax = axes[0]
    for i, T in enumerate(Ts):
        f = fits[f"{T:.0f}"]["contrast"]
        sub = sorted([r for r in rows if _f(r["T_K"]) == T
                      and _f(r["Oc_GHz"]) == oc0],
                     key=lambda r: _f(r["Bx_T"]))
        Bv = np.array([_f(r["Bx_T"]) for r in sub])
        Cv = np.array([_f(r["C"]) for r in sub])
        col = C[st.CAT_ORDER[i % len(st.CAT_ORDER)]]
        m = np.isfinite(Cv) & (Cv > 0) & (Bv > 0)
        ax.loglog(Bv[m], Cv[m], st.MARKERS[i % len(st.MARKERS)], color=col,
                  ms=3.6, ls="none", label=f"{T:.0f} K")
        rng = f.get("fit_range_T")
        if rng and np.isfinite(f.get("n", np.nan)):
            bb = np.logspace(np.log10(rng[0]), np.log10(rng[1]), 60)
            yy = f["y_res"] + f["a"] * bb ** f["n"]
            g = yy > 0
            ax.loglog(bb[g], yy[g], "-", color=col, lw=1.0, alpha=0.9)
    ax.axvline(B_pert, color=C["gray"], ls="--", lw=0.8)
    ax.axvline(B_mix, color=C["gray"], ls=":", lw=0.8)
    ax.text(B_pert * 0.92, 3e-6, r"$B_{\rm pert}$", fontsize=6.3,
            color=C["gray"], ha="right", rotation=90)
    ax.text(B_mix * 1.08, 3e-6, r"$\gamma_e B = D_{gs}$", fontsize=6.3,
            color=C["gray"], ha="left", rotation=90)
    ax.set_xlabel(r"transverse field $B_\perp$ (T)")
    ax.set_ylabel(r"sector contrast $C$")
    ax.legend(fontsize=6.5, loc="lower right")
    st.panel_label(ax, "a")

    # (b) exponent against fitting cutoff -- plateau vs drift
    ax = axes[1]
    for i, T in enumerate(Ts):
        f = fits[f"{T:.0f}"]["contrast"]
        w = f.get("n_by_window", {})
        if not w:
            continue
        cuts = np.array([v["b_max"] for v in w.values()])
        ns = np.array([v["n"] for v in w.values()])
        er = np.array([v["n_err"] for v in w.values()])
        o = np.argsort(cuts)
        col = C[st.CAT_ORDER[i % len(st.CAT_ORDER)]]
        rho = f.get("drift_rho")
        tag = ("drifts" if f.get("monotone_drift") else "stable")
        ax.errorbar(cuts[o], ns[o], yerr=er[o], fmt=st.MARKERS[i % 6] + "-",
                    color=col, ms=3.4, lw=1.1, capsize=2,
                    label=f"{T:.0f} K ({tag})")
    ax.axhline(2.0, color=C["black"], ls=":", lw=1.0)
    ax.text(ax.get_xlim()[1], 2.04, r"$n=2$", fontsize=6.5, ha="right",
            va="bottom")
    ax.set_xlabel(r"upper fitting cutoff $B_{\rm max}$ (T)")
    ax.set_ylabel(r"fitted exponent $n$")
    ax.legend(fontsize=6.3, loc="upper right")
    st.panel_label(ax, "b")

    fig.tight_layout()
    _save(fig, "fig5_bperp_scaling")


# --------------------------------------------------------------- Fig. 6
def fig6_observables():
    rows = _load_csv("p3_observables_vs_T.csv")
    rows = [r for r in rows if np.isfinite(_f(r["C"]))]
    T = np.array([_f(r["T_K"]) for r in rows])
    dTT = np.array([_f(r["dT_over_T_m"]) for r in rows])
    dPL = np.array([_f(r["dPL_over_PL_m"]) for r in rows])
    tt = np.array([_f(r["tau_trans_m_s"]) for r in rows])
    tp = np.array([_f(r["tau_pl_m_s"]) for r in rows])
    ceil = _load_json("p3_summary.json")["detection_chain"]["tau_ceiling_s"]["value"]

    fig, axes = plt.subplots(1, 2, figsize=(st.COL_2, 2.8))

    ax = axes[0]
    for y, col, mk, lab in ((dTT, C["blue"], "o", r"$\Delta T/T$"),
                            (dPL, C["purple"], "s", r"$\Delta {\rm PL}/{\rm PL}$")):
        p = y > 0; n = y < 0
        ax.semilogy(T[p], np.abs(y[p]), mk, color=col, ms=3.4, ls="-", lw=1.1,
                    label=lab + " > 0")
        ax.semilogy(T[n], np.abs(y[n]), mk, color=col, ms=3.4, ls="--", lw=1.1,
                    mfc="none", label=lab + " < 0")
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("fractional signal (OD-matched sample)")
    ax.legend(fontsize=6, loc="lower left")
    st.panel_label(ax, "a")

    ax = axes[1]
    ax.semilogy(T, np.where(np.isfinite(tt), tt, np.nan), "o-",
                color=C["blue"], ms=3.4, lw=1.1, label="transmission")
    ax.semilogy(T, np.where(np.isfinite(tp), tp, np.nan), "s--",
                color=C["purple"], ms=3.4, lw=1.1, label="fluorescence")
    ax.axhline(ceil, color=C["vermillion"], ls=":", lw=1.0)
    ax.text(T.max(), ceil * 1.5, "24 h", ha="right", fontsize=6.5,
            color=C["vermillion"])
    und = ~np.isfinite(tt)
    if und.any():
        ax.axvspan(T[und].min(), T.max(), color=C["gray"], alpha=0.13, lw=0)
        ax.text(T[und].min() + 3, 1e-6, "undetectable", fontsize=6.5,
                color=C["gray"])
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel(r"integration time for SNR = 5 (s)")
    ax.legend(fontsize=6.5, loc="upper left")
    st.panel_label(ax, "b")

    fig.tight_layout()
    _save(fig, "fig6_observables")


def main(quick=False):
    st.apply_style()
    print("building PRA figures ...")
    fig1_level_scheme()
    fig2_spectra(quick)
    fig3_phase_diagram()
    fig4_contrast_vs_T()
    try:
        fig5_bperp_scaling()
    except FileNotFoundError:
        print("  skipping fig5 (p5 results not present yet)")
    fig6_observables()
    print(f"figures in {FIG}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    main(quick=a.quick)
