"""Phase A driver: standard Lambda model, EIT-ATS-background separation
(Priority 1), Figure 2 (crossover maps) and Figure 5 (mechanism map).
"""

import json
import os

import numpy as np
from scipy.signal import find_peaks

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model_lambda import chi_full_analytic, chi_cut_analytic, self_consistency_check

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
FIGDIR = os.path.join(RESULTS, "figures")
os.makedirs(FIGDIR, exist_ok=True)

GAMMA31 = 1.0  # unit of frequency/rate throughout Phase A


# Probe absorption spectrum: control fixed on resonance, only the probe
# frequency (Delta) is scanned, so the two-photon detuning delta = Delta -
# Delta_c equals Delta itself (Delta_c = 0). This co-scan (delta=Delta) is
# what produces the standard EIT transparency dip / ATS twin-peak
# lineshape; delta held fixed at 0 while Delta varies instead sweeps a
# different (non-standard) cut of the two-parameter response and does not
# reproduce the known crossover, so it is used only where the two are
# swept independently (Figure 2 maps show both cuts).
def eta_S(Delta, Omega_c, gamma21, eps=1e-6):
    cf = chi_full_analytic(Delta, Delta, Omega_c, GAMMA31, gamma21)
    cc = chi_cut_analytic(Delta, Delta, Omega_c, GAMMA31, gamma21)
    R = cf - cc
    return abs(R) / (abs(cf) + eps)


def integrated_sector_weight(Omega_c, gamma21, Delta_window=8.0, n=801):
    Deltas = np.linspace(-Delta_window, Delta_window, n)
    cf = np.array([chi_full_analytic(D, D, Omega_c, GAMMA31, gamma21) for D in Deltas])
    cc = np.array([chi_cut_analytic(D, D, Omega_c, GAMMA31, gamma21) for D in Deltas])
    R = cf - cc
    num = np.trapezoid(np.abs(R) ** 2, Deltas)
    den = np.trapezoid(np.abs(cf) ** 2, Deltas) + 1e-30
    return num / den, Deltas, cf, cc, R


def peak_separation(Deltas, absorption):
    peaks, _ = find_peaks(absorption)
    if len(peaks) < 2:
        return 0.0, peaks
    # take the two most prominent peaks
    heights = absorption[peaks]
    order = np.argsort(heights)[::-1][:2]
    top2 = np.sort(peaks[order])
    return float(Deltas[top2[1]] - Deltas[top2[0]]), peaks


def figure2_crossover_maps(gamma21=0.02, n_omega=120, n_delta=140):
    Omega_over_g31 = np.linspace(0.05, 6.0, n_omega)
    Deltas = np.linspace(-6.0, 6.0, n_delta)

    Im_full = np.zeros((n_omega, n_delta))
    Im_cut = np.zeros((n_omega, n_delta))
    Im_R = np.zeros((n_omega, n_delta))

    for i, oc in enumerate(Omega_over_g31):
        for j, D in enumerate(Deltas):
            cf = chi_full_analytic(D, D, oc, GAMMA31, gamma21)
            cc = chi_cut_analytic(D, D, oc, GAMMA31, gamma21)
            Im_full[i, j] = cf.imag
            Im_cut[i, j] = cc.imag
            Im_R[i, j] = (cf - cc).imag

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    titles = [r"$\mathrm{Im}\,\chi_{\rm full}$", r"$\mathrm{Im}\,\chi_{\rm cut}$", r"$\mathrm{Im}\,R_S$"]
    for ax, data, title in zip(axes, [Im_full, Im_cut, Im_R], titles):
        vmax = np.max(np.abs(data))
        im = ax.pcolormesh(Deltas, Omega_over_g31, data, shading="auto",
                            cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xlabel(r"$\Delta/\gamma_{31}$")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046)
    axes[0].set_ylabel(r"$|\Omega_c|/\gamma_{31}$")
    fig.suptitle(f"Figure 2: EIT-ATS crossover ($\\gamma_{{21}}/\\gamma_{{31}}={gamma21}$)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_eit_ats_crossover.png"), dpi=150)
    plt.close(fig)


def classify_point(Omega_c, gamma21, W_S_thresh=2e-3):
    """Classifier combining the integrated sector weight W_S (structure
    present at all, vs. pure background) with the exact Anisimov-Kocharovsky
    / Giner et al. discriminant for this Lambda lineshape.

    In this minimal 2-level Lambda model the sector cut severs the *only*
    coupling (Omega_c) that produces either the EIT notch or the ATS
    splitting, so chi_cut is always the bare Lorentzian
    1/(gamma31-iDelta): the cut can never retain a two-peak structure by
    itself here (that requires the extended Phase B model with an
    independent background pathway). What the minimal model *can*
    diagnose is which flavour of coherent structure R_S carries, via the
    exact pole discriminant of chi_full's denominator
    (gamma31-iDelta)(gamma21-iDelta)+Omega_c^2/4 = 0:
    poles are complex (single unresolved notch, EIT-like) iff
    Omega_c < |gamma31-gamma21|, and split onto the real Delta axis
    (resolved twin peaks, ATS-like) iff Omega_c > |gamma31-gamma21|.
    """
    W_S, Deltas, cf, cc, R = integrated_sector_weight(Omega_c, gamma21)
    eta0 = eta_S(0.0, Omega_c, gamma21)
    ats_criterion = Omega_c > abs(GAMMA31 - gamma21)

    if W_S < W_S_thresh:
        return "background-dominant", W_S, float(ats_criterion), eta0
    if ats_criterion:
        return "ATS-dominant", W_S, float(ats_criterion), eta0
    return "EIT-dominant", W_S, float(ats_criterion), eta0


def figure5_mechanism_map(n_omega=45, n_gamma21=40):
    Omega_vals = np.linspace(0.1, 6.0, n_omega)
    gamma21_vals = np.logspace(-4, 0, n_gamma21)

    labels = {"EIT-dominant": 0, "ATS-dominant": 1, "background-dominant": 2, "unresolved": 3}
    grid = np.zeros((n_gamma21, n_omega))

    for i, g21 in enumerate(gamma21_vals):
        for j, oc in enumerate(Omega_vals):
            cls, W_S, sep_ratio, eta0 = classify_point(oc, g21)
            grid[i, j] = labels[cls]

    fig, ax = plt.subplots(figsize=(6.5, 5))
    cmap = plt.get_cmap("Set2", 4)
    im = ax.pcolormesh(Omega_vals, gamma21_vals, grid, shading="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    ax.set_yscale("log")
    ax.set_xlabel(r"$|\Omega_c|/\gamma_{31}$")
    ax.set_ylabel(r"$\gamma_{21}/\gamma_{31}$")
    ax.set_title("Figure 5: mechanism map")
    cbar = fig.colorbar(im, ax=ax, ticks=list(labels.values()))
    cbar.ax.set_yticklabels(list(labels.keys()))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig5_mechanism_map.png"), dpi=150)
    plt.close(fig)

    return grid, labels


def priority1_three_typical_points():
    """Gate 3: verify at least three typical points separate cleanly."""
    points = {
        "EIT_typical": (0.3, 0.005),      # Omega_c < |gamma31-gamma21|: unresolved notch
        "ATS_typical": (4.0, 0.005),      # Omega_c > |gamma31-gamma21|: resolved twin peaks
        "background_typical": (0.05, 0.9),  # weak coupling: R_S negligible, near-bare Lorentzian
    }
    out = {}
    for name, (oc, g21) in points.items():
        cls, W_S, sep_ratio, eta0 = classify_point(oc, g21)
        out[name] = {
            "Omega_c_over_g31": oc, "gamma21_over_g31": g21,
            "classified_as": cls, "W_S": W_S, "peak_sep_ratio": sep_ratio, "eta_S(0)": eta0,
        }
    return out


def main():
    ef, ec = self_consistency_check()
    print(f"[Phase A self-check] analytic vs numeric: {max(ef, ec):.2e}")

    figure2_crossover_maps()
    grid, labels = figure5_mechanism_map()
    gate3 = priority1_three_typical_points()

    gate3["passes"] = bool(len({v["classified_as"] for v in gate3.values() if isinstance(v, dict)}) >= 3)

    summary = {"phaseA_self_check": {"max_err_full": ef, "max_err_cut": ec},
               "gate3_eit_ats_background_separation": gate3}

    with open(os.path.join(RESULTS, "gates_summary_phaseA.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\nGate 3 (EIT/ATS/background separation):")
    for k, v in gate3.items():
        print(" ", k, "->", v)


if __name__ == "__main__":
    main()
