#!/usr/bin/env python3
"""Phenomenological hBN EIT no-go/go scan.

This is not a microscopic certification of K12*K21. It maps experimentally
reported linewidth/coherence scales onto an ideal-Lambda visibility budget.
All rates are handled in ordinary frequency units (MHz); the common 2*pi
factor cancels in the dimensionless ratios.
"""
from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
FIGURES = ROOT / "figures"
ANALYSIS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)


def ideal_lambda_contrast(omega_mhz: np.ndarray | float,
                          gamma_e_mhz: np.ndarray | float,
                          gamma_g_mhz: float,
                          eta: float = 1.0) -> np.ndarray:
    """Conditional normalized absorption reduction.

    C = eta * Omega^2 / (Omega^2 + 4 gamma_e gamma_g).
    eta is a phenomenological kernel/branching quality, not measured for hBN.
    """
    omega = np.asarray(omega_mhz, dtype=float)
    gamma_e = np.asarray(gamma_e_mhz, dtype=float)
    return eta * omega**2 / (omega**2 + 4.0 * gamma_e * gamma_g_mhz)


def required_control(gamma_e_mhz: float,
                     gamma_g_mhz: float,
                     target_contrast: float,
                     eta: float = 1.0) -> float:
    """Omega required by the conditional visibility proxy, in MHz."""
    if not (0.0 < target_contrast < eta <= 1.0):
        return math.inf
    return 2.0 * math.sqrt(
        gamma_e_mhz * gamma_g_mhz * target_contrast / (eta - target_contrast)
    )


def resolvedness_proxy(delta_mhz: np.ndarray | float,
                       gamma_e_mhz: np.ndarray | float) -> np.ndarray:
    """Heuristic branch-resolution factor eta_res=r^2/(1+r^2), r=Delta/gamma_e.

    This is a visualization aid only; it is not the theorem's exact kernel.
    """
    delta = np.asarray(delta_mhz, dtype=float)
    gamma_e = np.asarray(gamma_e_mhz, dtype=float)
    r = delta / gamma_e
    return r**2 / (1.0 + r**2)


# Ground-state spin-coherence proxies from Stern et al. 2024.
t2star_ns = 106.0
gamma_g_t2star_mhz = 1000.0 / (math.pi * t2star_ns)  # 1/(pi T2*) in MHz
gamma_g_conservative_mhz = 10.0  # unsaturated ODMR linewidth scale

# Experiment-derived optical linewidth scenarios.
# FWHM values are used conservatively as gamma_e proxies. If the theory convention
# requires HWHM, required controls decrease by sqrt(2).
scenarios = [
    {
        "configuration": "Mechanically decoupled emitter (Dietrich 2020)",
        "temperature_K": 300,
        "optical_width_MHz": 65.0,
        "ground_width_MHz": gamma_g_t2star_mhz,
        "same_defect_inputs": False,
        "evidence_status": "best-case cross-defect synthesis",
        "classification": "Class 0; conditional go candidate",
        "note": "Room-temperature near-FT optical line, but no verified spin Lambda in that emitter.",
    },
    {
        "configuration": "Mechanically decoupled emitter + conservative spin width",
        "temperature_K": 300,
        "optical_width_MHz": 65.0,
        "ground_width_MHz": gamma_g_conservative_mhz,
        "same_defect_inputs": False,
        "evidence_status": "best-case cross-defect synthesis",
        "classification": "Class 0; conditional go candidate",
        "note": "Uses 10 MHz spin-dephasing proxy; Lambda/dipoles remain unknown.",
    },
    {
        "configuration": "Akbari T^3 stress-test extrapolation",
        "temperature_K": 300,
        "optical_width_MHz": np.nan,
        "ground_width_MHz": gamma_g_t2star_mhz,
        "same_defect_inputs": False,
        "evidence_status": "extrapolation beyond measured 120 K",
        "classification": "Class 0 / likely practical no at fixed power",
        "note": "Temperature model is a stress test, not a measured 300 K linewidth.",
    },
    {
        "configuration": "Whitefield representative room-temperature PL band",
        "temperature_K": 300,
        "optical_width_MHz": 1_356_500.0,
        "ground_width_MHz": gamma_g_t2star_mhz,
        "same_defect_inputs": True,
        "evidence_status": "PL bandwidth proxy, not homogeneous PLE linewidth",
        "classification": "Not certifiable; practical no under this proxy",
        "note": "Optical spin readout exists, but PL FWHM cannot be substituted for homogeneous linewidth in an exact test.",
    },
]

# Akbari: fit anchored to 89 MHz at 6.5 K and ~1.15 GHz at 120 K.
T0 = 6.5
gamma0_ghz = 0.089
A_akbari_ghz_per_K3 = (1.15 - gamma0_ghz) / (120.0**3 - T0**3)

def akbari_width_mhz(T: np.ndarray | float) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    return 1000.0 * (gamma0_ghz + A_akbari_ghz_per_K3 * (T**3 - T0**3))

scenarios[2]["optical_width_MHz"] = float(akbari_width_mhz(300.0))

# Numerical visibility results.
rows: list[dict[str, object]] = []
for sc in scenarios:
    for eta in (1.0, 0.5, 0.2, 0.1):
        for target in (0.10, 0.50, 0.90):
            omega = required_control(
                float(sc["optical_width_MHz"]),
                float(sc["ground_width_MHz"]),
                target,
                eta,
            )
            rows.append({
                **sc,
                "eta_assumed": eta,
                "target_contrast": target,
                "required_control_MHz": omega,
                "required_control_GHz": omega / 1000.0 if math.isfinite(omega) else math.inf,
            })
results = pd.DataFrame(rows)
results.to_csv(ANALYSIS / "hbn_numerical_results.csv", index=False)

# A compact room-temperature comparison table at eta=1.
rt_rows = []
for sc in scenarios:
    rt_rows.append({
        "configuration": sc["configuration"],
        "T_K": sc["temperature_K"],
        "gamma_e_proxy_MHz": sc["optical_width_MHz"],
        "gamma_g_proxy_MHz": sc["ground_width_MHz"],
        "Omega_for_50pct_MHz_eta1": required_control(float(sc["optical_width_MHz"]), float(sc["ground_width_MHz"]), 0.50, 1.0),
        "Omega_for_90pct_MHz_eta1": required_control(float(sc["optical_width_MHz"]), float(sc["ground_width_MHz"]), 0.90, 1.0),
        "classification": sc["classification"],
        "limitation": sc["note"],
    })
pd.DataFrame(rt_rows).to_csv(ANALYSIS / "room_temperature_scenarios.csv", index=False)

# Parameter table.
parameter_rows = [
    {
        "platform_or_emitter": "Mechanically decoupled emitter",
        "reference": "Dietrich et al., Phys. Rev. B 101, 081401(R) (2020)",
        "temperature_range_K": "3-300",
        "optical_linewidth": "64.7, 75.5, 65.5, 62.8 MHz at 3, 100, 200, 300 K",
        "spin_parameter": "not established for the same emitter",
        "usable_nogo_input": "gamma_e(T) only",
        "assessment": "optical gate passes at 300 K; Lambda/kernel unknown",
    },
    {
        "platform_or_emitter": "Mechanically decoupled emitter",
        "reference": "Hoese et al., Sci. Adv. 6, eaba6038 (2020)",
        "temperature_range_K": "cryogenic to 300",
        "optical_linewidth": "61 +/- 10 MHz homogeneous; low-frequency phonon gap",
        "spin_parameter": "not established",
        "usable_nogo_input": "small optical damping / phonon-gap mechanism",
        "assessment": "does not by itself prove ker(D) or EIT capability",
    },
    {
        "platform_or_emitter": "Mechanically decoupled emitter",
        "reference": "Koch et al., Commun. Mater. 5, 240 (2024)",
        "temperature_range_K": "5-150",
        "optical_linewidth": "FTL 109(9) MHz; full homogeneous 324(16) MHz; gap 2.42(11) THz",
        "spin_parameter": "not established",
        "usable_nogo_input": "gamma_e and coherent-drive temperature boundary",
        "assessment": "FTL-like to about 50 K; coherent Rabi limit between 20-30 K for this emitter",
    },
    {
        "platform_or_emitter": "Electric-field-stabilized emitter",
        "reference": "Akbari et al., Nano Lett. 22, 7798 (2022)",
        "temperature_range_K": "6.5-120",
        "optical_linewidth": "89 MHz at 6.5 K; approximately T^3 broadening",
        "spin_parameter": "not established",
        "usable_nogo_input": "gamma_e(T) empirical fit",
        "assessment": "no measured 300 K input; extrapolation only",
    },
    {
        "platform_or_emitter": "Typical visible emitter",
        "reference": "White et al., Optica 8, 1153 (2021)",
        "temperature_range_K": "4-40",
        "optical_linewidth": "about 1.10 GHz at 5 K; 36+0.32 T^3 MHz above 20 K",
        "spin_parameter": "not established",
        "usable_nogo_input": "low-temperature gamma_e(T)",
        "assessment": "strong temperature broadening; fit must not be extrapolated to 300 K as data",
    },
    {
        "platform_or_emitter": "Carbon-related single spin",
        "reference": "Stern et al., Nat. Mater. 23, 1379 (2024)",
        "temperature_range_K": "ambient",
        "optical_linewidth": "resonant homogeneous width not reported",
        "spin_parameter": "S=1, D/h=1.959 GHz, T2*=106(12) ns, DD coherence 1.08 us",
        "usable_nogo_input": "gamma_g only",
        "assessment": "spin gate passes; optical Lambda and gamma_e unknown",
    },
    {
        "platform_or_emitter": "Narrowband optically addressable spin complexes",
        "reference": "Whitefield et al., arXiv:2501.15341 / Nat. Mater. (2026)",
        "temperature_range_K": "300",
        "optical_linewidth": "representative PL FWHM 2.76 nm at 781 nm (about 1.36 THz)",
        "spin_parameter": "S=1 and S=1/2; >25% optical spin-readout yield",
        "usable_nogo_input": "spin existence and broad PL proxy only",
        "assessment": "homogeneous PLE linewidth and coherent optical Lambda not established",
    },
    {
        "platform_or_emitter": "B center",
        "reference": "Horder et al., ACS Photonics 12, 1284 (2025); Gérard et al., Nat. Commun. (2026)",
        "temperature_range_K": "cryogenic",
        "optical_linewidth": "PLE 0.88-1.15 GHz in favorable samples; coherent two-level drive",
        "spin_parameter": "long-lived lower-state pair not established",
        "usable_nogo_input": "optical coherence only",
        "assessment": "two-level optical go, but not an EIT Lambda configuration",
    },
    {
        "platform_or_emitter": "V_B^- ensemble",
        "reference": "Mathur et al., Nat. Commun. 13, 3233 (2022)",
        "temperature_range_K": "300",
        "optical_linewidth": "broad PL/strong ISC; resonant optical homogeneous width unavailable",
        "spin_parameter": "ground ZFS about 3.5 GHz; excited ZFS 2.1 GHz",
        "usable_nogo_input": "spin Hamiltonian only",
        "assessment": "no resolved shared optical Lambda; presently practical no/undetermined",
    },
]
pd.DataFrame(parameter_rows).to_csv(ANALYSIS / "hbn_parameter_table.csv", index=False)

# Figure 1: linewidth models vs temperature.
T = np.linspace(5.0, 300.0, 500)
plt.figure(figsize=(8, 5.2))
plt.plot(T, np.full_like(T, 65.0), label="Dietrich best case: ~65 MHz")
plt.plot(T, akbari_width_mhz(T), label="Akbari T^3 stress-test")
white = np.where((T >= 20.0) & (T <= 40.0), 36.0 + 0.32 * T**3, np.nan)
plt.plot(T, white, label="White fit (valid 20-40 K)")
plt.axvline(120.0, linestyle="--", linewidth=1.0, label="Akbari measured limit")
plt.yscale("log")
plt.xlabel("Temperature (K)")
plt.ylabel("Optical linewidth proxy (MHz)")
plt.title("Reported and stress-test hBN linewidth scenarios")
plt.grid(True, which="both", alpha=0.25)
plt.legend()
plt.tight_layout()
plt.savefig(FIGURES / "linewidth_models_vs_temperature.png", dpi=220)
plt.close()

# Figure 2: required control vs optical linewidth.
gamma_e_grid = np.logspace(np.log10(50.0), np.log10(1.0e7), 500)
plt.figure(figsize=(8, 5.2))
for gg in (gamma_g_t2star_mhz, gamma_g_conservative_mhz):
    for target in (0.50, 0.90):
        omega = 2.0 * np.sqrt(gamma_e_grid * gg * target / (1.0 - target))
        plt.plot(gamma_e_grid, omega, label=f"gamma_g={gg:.1f} MHz, C={target:.0%}")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Optical linewidth proxy gamma_e (MHz)")
plt.ylabel("Required control Omega_c (MHz), eta=1")
plt.title("Conditional ideal-Lambda control budget")
plt.grid(True, which="both", alpha=0.25)
plt.legend()
plt.tight_layout()
plt.savefig(FIGURES / "required_control_vs_optical_linewidth.png", dpi=220)
plt.close()

# Figure 3: conditional contrast vs temperature at fixed control.
omega_fixed = 100.0
plt.figure(figsize=(8, 5.2))
plt.plot(T, ideal_lambda_contrast(omega_fixed, np.full_like(T, 65.0), gamma_g_t2star_mhz), label="Dietrich best case")
plt.plot(T, ideal_lambda_contrast(omega_fixed, akbari_width_mhz(T), gamma_g_t2star_mhz), label="Akbari T^3 stress-test")
white_contrast = ideal_lambda_contrast(omega_fixed, 36.0 + 0.32*T**3, gamma_g_t2star_mhz)
white_contrast[(T < 20.0) | (T > 40.0)] = np.nan
plt.plot(T, white_contrast, label="White fit (valid 20-40 K)")
plt.axvline(120.0, linestyle="--", linewidth=1.0, label="Akbari measured limit")
plt.ylim(0.0, 1.02)
plt.xlabel("Temperature (K)")
plt.ylabel("Conditional contrast proxy")
plt.title("Ideal-Lambda visibility at Omega_c = 100 MHz, eta = 1")
plt.grid(True, alpha=0.25)
plt.legend()
plt.tight_layout()
plt.savefig(FIGURES / "conditional_contrast_vs_temperature.png", dpi=220)
plt.close()

# Figure 4: heuristic unresolved-manifold penalty.
gamma_grid = np.logspace(1.0, 7.0, 500)
plt.figure(figsize=(8, 5.2))
for delta_ghz in (0.85, 2.0, 3.5):
    eta_res = resolvedness_proxy(delta_ghz * 1000.0, gamma_grid)
    plt.plot(gamma_grid, eta_res, label=f"Delta_sel={delta_ghz:g} GHz")
plt.xscale("log")
plt.ylim(0.0, 1.02)
plt.xlabel("Optical linewidth proxy gamma_e (MHz)")
plt.ylabel("Heuristic resolution factor eta_res")
plt.title("Unresolved-manifold penalty (heuristic, not theorem)")
plt.grid(True, which="both", alpha=0.25)
plt.legend()
plt.tight_layout()
plt.savefig(FIGURES / "resolvedness_kernel_proxy.png", dpi=220)
plt.close()

print(f"gamma_g(T2*) = {gamma_g_t2star_mhz:.6f} MHz")
print(f"Akbari 300 K stress-test width = {float(akbari_width_mhz(300.0)):.3f} MHz")
for sc in scenarios:
    ge = float(sc["optical_width_MHz"])
    gg = float(sc["ground_width_MHz"])
    print(sc["configuration"])
    print("  Omega50 =", required_control(ge, gg, 0.5, 1.0), "MHz")
    print("  Omega90 =", required_control(ge, gg, 0.9, 1.0), "MHz")
