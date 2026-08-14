"""Slide-sized versions of the manuscript figures.

The manuscript figures are drawn for a two-column journal page: 7 pt tick
labels in a 7-inch-wide canvas.  Dropped onto a slide they scale up by at most
1.3x, so their smallest text lands near 9 pt on the projector -- half the floor
this deck works to.

Rescaling the PNG cannot fix that, because enlarging the canvas enlarges the
data area with it.  The only way to get 18-20 pt text is to redraw at the size
the figure will actually occupy on the slide, with the font sizes set to what
they should be there.  Each function below therefore takes its width and height
in inches straight from the slide placement in build.py, so one matplotlib
point is one slide point and nothing is scaled afterwards.

That budget is unforgiving: an 18 pt tick label is a quarter of an inch tall,
so a panel that carried eight ticks in the paper carries three here.  Panels
whose detail cannot survive that are not shrunk -- they are left in the
manuscript, and the slide shows the one panel that carries the message.

Usage: python slide_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PRA = REPO / "calculations" / "numerics" / "manuscript_figures" / "pra" / "src"
RES = REPO / "results"
OUT = HERE / "figures"
sys.path.insert(0, str(PRA))
sys.path.insert(0, str(REPO / "calculations" / "numerics" / "No-go theorem" / "src"))
sys.path.insert(0, str(REPO / "calculations" / "numerics" / "manuscript_figures"
                      / "prl" / "src"))

import matplotlib                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt         # noqa: E402
import numpy as np                      # noqa: E402
from matplotlib.colors import (BoundaryNorm, LinearSegmentedColormap,  # noqa: E402
                               ListedColormap, SymLogNorm)
from matplotlib.patches import FancyArrowPatch  # noqa: E402

# The deck's palette, and no more: navy carries structure, blue means the
# pathway produces transparency, red means it produces absorption, amber is
# the one Autler-Townes corner, grey is everything unresolved or auxiliary.
NAVY = "#1C3077"
BLUE = "#2E5AAC"
RED = "#C00000"
AMBER = "#B88000"
GREY = "#767676"
PALE = "#DDDDDD"

BASE = 18          # ticks, legends, in-axes annotation
LABEL = 20         # axis labels
DPI = 200


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": BASE,
        "axes.labelsize": LABEL,
        "axes.titlesize": LABEL,
        "xtick.labelsize": BASE,
        "ytick.labelsize": BASE,
        "legend.fontsize": BASE,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "lines.linewidth": 2.4,
        "lines.markersize": 7,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    p = OUT / f"{name}.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  {p.name}")


def _json(name):
    with (RES / "tables" / name).open() as fh:
        return json.load(fh)


def _csv(name):
    import csv
    with (RES / "tables" / name).open() as fh:
        return list(csv.DictReader(fh))


def _f(x, default=np.nan):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ------------------------------------------------------------------ Lambda
def lambda_scheme(w=5.00, h=3.50):
    """The three-level channel and what closes it."""
    fig, ax = plt.subplots(figsize=(w, h))
    # the labels sit outside the level bars, which leaves the strip under the
    # ground pair free for the coherence arrow and its caption
    ax.set_xlim(-1.4, 11.0)
    ax.set_ylim(-0.6, 10)
    ax.axis("off")

    ax.plot([1.0, 4.2], [1.2, 1.2], color="black", lw=3)
    ax.plot([5.8, 9.0], [2.6, 2.6], color="black", lw=3)
    ax.plot([3.4, 6.6], [8.4, 8.4], color=NAVY, lw=3)
    ax.text(0.75, 1.2, r"$|1\rangle$", ha="right", va="center", fontsize=LABEL)
    ax.text(9.25, 2.6, r"$|2\rangle$", ha="left", va="center", fontsize=LABEL)
    ax.text(5.0, 9.05, r"$|e\rangle$", ha="center", fontsize=LABEL)

    ax.annotate("", xy=(4.35, 8.25), xytext=(2.75, 1.35),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.6,
                                mutation_scale=22))
    ax.annotate("", xy=(5.75, 8.25), xytext=(7.25, 2.75),
                arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.6,
                                mutation_scale=22))
    ax.text(2.55, 5.2, "probe", color=BLUE, fontsize=BASE, ha="right",
            rotation=77, va="center")
    ax.text(7.55, 5.6, "control", color=RED, fontsize=BASE, ha="left",
            rotation=-75, va="center")

    ax.annotate("", xy=(5.7, 2.45), xytext=(4.3, 1.35),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=1.8, ls="--",
                                mutation_scale=16))
    ax.text(5.0, 0.55, "ground coherence", color=GREY, fontsize=BASE,
            ha="center", va="center")
    save(fig, "sf_lambda")


# ------------------------------------------------------------------ island
def island(w=9.55, h=3.62):
    """Classification map and signed contrast over the (T, B) plane."""
    z = np.load(RES / "raw" / "p1_phase_diagram.npz", allow_pickle=True)
    T, B, K, Cg = z["T_K"], z["Bx_T"], z["klass"], z["C"]
    code = json.loads(str(z["class_code"]))
    dB = (B[1] - B[0]) if len(B) > 1 else 0.05
    dT = (T[1] - T[0]) if len(T) > 1 else 5.0
    be = np.append(B, B[-1] + dB) - dB / 2
    te = np.append(T, T[-1] + dT) - dT / 2

    fig, axes = plt.subplots(1, 2, figsize=(w, h))

    ax = axes[0]
    order = ["transparency", "ATS", "absorption", "unresolved"]
    cols = [BLUE, AMBER, RED, PALE]
    remap = np.full(K.shape, 3)
    for i, nm in enumerate(order):
        remap[K == code[nm]] = i
    ax.pcolormesh(be, te, remap, cmap=ListedColormap(cols),
                  norm=BoundaryNorm(np.arange(-0.5, 4.5, 1), 4),
                  shading="flat")
    ax.set_xlabel(r"$B_\perp$ (T)")
    ax.set_ylabel("temperature (K)")
    ax.set_xticks([0.0, 0.25, 0.5])
    ax.set_yticks([20, 60, 100])
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="none")
               for c in (BLUE, RED, AMBER)]
    ax.legend(handles, ["transparency", "absorption", "Autler-Townes"],
              fontsize=BASE - 3, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), ncol=3, frameon=False,
              handlelength=1.0, columnspacing=1.0, handletextpad=0.4)

    ax = axes[1]
    cmap = LinearSegmentedColormap.from_list("ta", [RED, "#F7F7F7", BLUE])
    im = ax.pcolormesh(be, te, np.where(np.isfinite(Cg), Cg, 0.0), cmap=cmap,
                       norm=SymLogNorm(linthresh=1e-6, vmin=-1.0, vmax=1.0),
                       shading="flat")
    cb = fig.colorbar(im, ax=ax, pad=0.03, ticks=[-1, -1e-3, 0, 1e-3, 1])
    cb.ax.set_yticklabels([r"$-1$", r"$-10^{-3}$", "0", r"$10^{-3}$", "1"])
    cb.ax.tick_params(labelsize=BASE - 2)
    cb.set_label("contrast $C$", fontsize=BASE)
    ax.set_xlabel(r"$B_\perp$ (T)")
    ax.set_xticks([0.0, 0.25, 0.5])
    ax.set_yticks([20, 60, 100])
    ax.set_ylabel("temperature (K)")

    fig.tight_layout()
    save(fig, "sf_island")


# --------------------------------------------------------- contrast vs T
def contrast_vs_T(w=6.30, h=3.62):
    """Sector contrast down the candidate field, with the two boundaries."""
    z = np.load(RES / "raw" / "p1_phase_diagram.npz", allow_pickle=True)
    T, B, Cg = z["T_K"], z["Bx_T"], z["C"]
    import run_prl_prediction as rp
    Cc = Cg[:, int(np.argmin(np.abs(B - rp.BX0)))]

    p2 = _json("p2_summary.json")
    key = f"{rp.BX0:.5f}"
    b1 = p2["T_1pct_band_vs_field"].get(key, {})
    bs = p2["T_sign_band_vs_field"].get(key, {})

    fig, ax = plt.subplots(figsize=(w, h))
    m = np.isfinite(Cc)

    def runs(mask):
        out, cur = [], []
        for i, v in enumerate(mask):
            if v:
                cur.append(i)
            elif cur:
                out.append(cur)
                cur = []
        if cur:
            out.append(cur)
        return out

    for k, idx in enumerate(runs(m & (Cc > 0))):
        ax.semilogy(T[idx], np.abs(Cc[idx]), "o-", color=BLUE,
                    label="transparency" if k == 0 else None)
    for k, idx in enumerate(runs(m & (Cc < 0))):
        ax.semilogy(T[idx], np.abs(Cc[idx]), "s--", color=RED,
                    label="absorption" if k == 0 else None)

    for band, col in ((b1, GREY), (bs, RED)):
        if band.get("q16") is not None:
            ax.axvspan(band["q16"], band["q84"], color=col, alpha=0.15, lw=0)
    ax.axhline(1.534e-7, color=GREY, ls=":")
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel(r"contrast $|C|$")
    ax.set_ylim(1e-7, 30)
    ax.set_yticks([1e-6, 1e-4, 1e-2, 1e0])
    ax.set_xticks([20, 60, 100])
    ax.legend(fontsize=BASE, loc="upper right", frameon=False,
              handlelength=1.4, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, "sf_contrast_T")


# ------------------------------------------------------------- lineshape
def lineshape(w=5.35, h=3.62, T=105.0):
    """The reversed lineshape: the slide's whole point is the sign."""
    import p1_phase_diagram as p1
    import run_prl_prediction as rp
    cache = OUT / f"spectrum_{T:.0f}K.npz"
    if cache.exists():
        d = np.load(cache)
        det, Af, Ac = d["det"], d["Af"], d["Ac"]
    else:
        det, Af, Ac, _, _ = p1.adaptive_spectrum(T, rp.BX0)
        OUT.mkdir(exist_ok=True)
        np.savez(cache, det=det, Af=Af, Ac=Ac)

    fig, ax = plt.subplots(figsize=(w, h))
    base = float(np.nanmedian(Ac))
    ax.plot(det, (Ac - base) * 1e7, ls="--", color=GREY, label="pathway cut")
    ax.plot(det, (Af - base) * 1e7, color=RED, label="full")
    ax.axhline(0, color="black", lw=0.9)
    ax.set_xlabel(r"$\delta_2$ (MHz)")
    ax.set_ylabel(r"absorption  ($10^{-7}$, offset)")
    ax.set_xticks([-1, 0, 1])
    ax.set_ylim(top=float(np.nanmax((Af - base) * 1e7)) * 1.75)
    ax.legend(fontsize=BASE, loc="upper left", frameon=False,
              handlelength=1.4, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, "sf_lineshape")


# ----------------------------------------------------------- observables
def observables(w=5.35, h=3.62):
    """Integration time to SNR = 5, and where it runs out."""
    rows = [r for r in _csv("p3_observables_vs_T.csv")
            if np.isfinite(_f(r["C"]))]
    T = np.array([_f(r["T_K"]) for r in rows])
    tt = np.array([_f(r["tau_trans_m_s"]) for r in rows])
    tp = np.array([_f(r["tau_pl_m_s"]) for r in rows])
    ceil = _json("p3_summary.json")["detection_chain"]["tau_ceiling_s"]["value"]

    fig, ax = plt.subplots(figsize=(w, h))
    ax.semilogy(T, np.where(np.isfinite(tt), tt, np.nan), "o-", color=BLUE,
                label="transmission")
    ax.semilogy(T, np.where(np.isfinite(tp), tp, np.nan), "s--", color=NAVY,
                label="fluorescence")
    ax.axhline(ceil, color=RED, ls=":")
    ax.text(T.max(), ceil * 3, "24 h", ha="right", fontsize=BASE, color=RED)
    und = ~np.isfinite(tt)
    if und.any():
        ax.axvspan(T[und].min(), T.max(), color=GREY, alpha=0.15, lw=0)
        ax.text(T[und].min() + 6, 1e-7, "undetectable", fontsize=BASE,
                color=GREY)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("time to SNR = 5 (s)")
    ax.set_xticks([0, 100, 200, 300])
    ax.set_yticks([1e-8, 1e-4, 1e0, 1e4])
    ax.legend(fontsize=BASE, loc="upper left", frameon=False,
              handlelength=1.4, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, "sf_observables")


# --------------------------------------------------------------- B-perp
def bperp(w=9.55, h=3.62):
    """The field dependence, and the exponent's dependence on the window."""
    rows = _csv("p5_bperp_scaling.csv")
    summ = _json("p5_summary.json")
    B_pert = summ["physical_window"]["B_pert_T"]
    fits = summ["temperature_fits"]
    oc0 = 0.1
    Ts = sorted({_f(r["T_K"]) for r in rows if _f(r["Oc_GHz"]) == oc0})
    cols = [NAVY, BLUE, GREY]
    marks = ["o", "s", "^"]

    fig, axes = plt.subplots(1, 2, figsize=(w, h))

    ax = axes[0]
    for i, T in enumerate(Ts):
        f = fits[f"{T:.0f}"]["contrast"]
        sub = sorted([r for r in rows if _f(r["T_K"]) == T
                      and _f(r["Oc_GHz"]) == oc0],
                     key=lambda r: _f(r["Bx_T"]))
        Bv = np.array([_f(r["Bx_T"]) for r in sub])
        Cv = np.array([_f(r["C"]) for r in sub])
        m = np.isfinite(Cv) & (Cv > 0) & (Bv > 0)
        ax.loglog(Bv[m], Cv[m], marks[i % 3], color=cols[i % 3], ls="none",
                  label=f"{T:.0f} K")
        rng = f.get("fit_range_T")
        if rng and np.isfinite(f.get("n", np.nan)):
            bb = np.logspace(np.log10(rng[0]), np.log10(rng[1]), 60)
            yy = f["y_res"] + f["a"] * bb ** f["n"]
            ax.loglog(bb[yy > 0], yy[yy > 0], "-", color=cols[i % 3], lw=1.8)
    ax.axvline(B_pert, color=GREY, ls="--", lw=1.4)
    ax.set_xlabel(r"$B_\perp$ (T)")
    ax.set_ylabel(r"contrast $C$")
    ax.set_yticks([1e-6, 1e-3, 1e0])
    ax.set_ylim(top=8.0)
    ax.legend(fontsize=BASE, loc="upper left", frameon=False,
              handlelength=1.0, borderaxespad=0.2)

    ax = axes[1]
    for i, T in enumerate(Ts):
        f = fits[f"{T:.0f}"]["contrast"]
        w_ = f.get("n_by_window", {})
        if not w_:
            continue
        cuts = np.array([v["b_max"] for v in w_.values()])
        ns = np.array([v["n"] for v in w_.values()])
        er = np.array([v["n_err"] for v in w_.values()])
        o = np.argsort(cuts)
        tag = "drifts" if f.get("monotone_drift") else "stable"
        ax.errorbar(cuts[o], ns[o], yerr=er[o], fmt=marks[i % 3] + "-",
                    color=cols[i % 3], capsize=3, lw=1.8,
                    label=f"{T:.0f} K, {tag}")
    ax.axhline(2.0, color=RED, ls=":")
    ax.set_xlabel(r"fitting cutoff $B_{\max}$ (T)")
    ax.set_ylabel(r"exponent $n$")
    ax.set_xticks([0.025, 0.035, 0.045])
    ax.set_yticks([2.0, 2.5, 3.0])
    ax.set_ylim(top=3.9)          # headroom so the legend clears the curves
    ax.legend(fontsize=BASE - 3, loc="upper right", frameon=False,
              handlelength=1.2, borderaxespad=0.2)

    fig.tight_layout()
    save(fig, "sf_bperp")


def main():
    style()
    print("building slide figures ...")
    lambda_scheme()
    island()
    contrast_vs_T()
    lineshape()
    observables()
    bperp()
    print(f"figures in {OUT}")


if __name__ == "__main__":
    main()
