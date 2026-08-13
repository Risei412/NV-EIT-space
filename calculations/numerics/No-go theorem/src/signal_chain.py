"""signal_chain.py -- Gate 3 of SIMULATION_PLAN.md: pure conversion
functions from the local susceptibility contrast to a measurable signal,

  delta_chi -> delta_alpha -> Delta OD -> Delta T/T -> Delta N_ph -> SNR.

All functions are side-effect free and unit-annotated so they can be
unit-tested (tests/test_core.py). Wavelengths in nm, lengths in cm,
rates/linewidths in GHz (ordinary frequency, repo convention), powers in W,
times in s. "sector OD" is the optical depth of the probed (post-selected)
ZPL sub-ensemble; the EIT contrast C from the kernel/Liouvillian models is
a relative change of exactly that absorption, so Delta OD = OD_sector * C.
"""
from __future__ import annotations
import numpy as np

H_PLANCK = 6.62607015e-34          # J s
C_LIGHT = 2.99792458e10            # cm / s

def photon_energy_J(lambda_nm: float) -> float:
    return H_PLANCK*C_LIGHT/(lambda_nm*1e-7)

def sigma_zpl_cm2(lambda_nm: float, n_refr: float, debye_waller: float,
                  gamma_rad_GHz: float, gamma_h_GHz: float) -> float:
    """Peak ZPL absorption cross-section: two-level resonant cross-section
    3 lambda_medium^2 / 2 pi, scaled by the Debye-Waller factor and by the
    ratio of radiative to total homogeneous linewidth."""
    lam_cm = lambda_nm*1e-7/n_refr
    return 3*lam_cm**2/(2*np.pi)*debye_waller*(gamma_rad_GHz/gamma_h_GHz)

def spectral_fraction(gamma_h_GHz: float, gamma_inh_GHz: float) -> float:
    """Fraction of the inhomogeneous ensemble resonant with the probe."""
    return min(1.0, gamma_h_GHz/max(gamma_inh_GHz, 1e-300))

def alpha_cm(sigma_cm2: float, n_nv_cm3: float, f_orient: float = 0.25,
             f_spin: float = 1/3, f_spectral: float = 1.0) -> float:
    """Absorption coefficient of the probed sector (cm^-1)."""
    return sigma_cm2*n_nv_cm3*f_orient*f_spin*f_spectral

def od(alpha_percm: float, L_cm: float) -> float:
    return alpha_percm*L_cm

def delta_od(od_sector: float, contrast: float) -> float:
    """EIT-induced change of optical depth (contrast = relative change of
    the sector absorption, dA/Acut of the models)."""
    return od_sector*contrast

def transmission(od_total: float) -> float:
    return float(np.exp(-od_total))

def delta_T_over_T(d_od: float) -> float:
    """Relative transmission change; exact for constant background."""
    return float(np.expm1(d_od))

def detected_photons(power_W: float, lambda_nm: float, tau_s: float,
                     eta: float, od_total: float) -> float:
    return power_W/photon_energy_J(lambda_nm)*tau_s*eta*transmission(od_total)

def snr(d_od: float, od_total: float, power_W: float, lambda_nm: float,
        tau_s: float, eta: float, sigma_tech: float = 0.0) -> float:
    """Shot-noise + relative technical-noise SNR of the transparency signal."""
    N = detected_photons(power_W, lambda_nm, tau_s, eta, od_total)
    if N <= 0: return 0.0
    rel = abs(delta_T_over_T(d_od))
    return rel/np.sqrt(1.0/N + sigma_tech**2)

def required_tau_s(target_snr: float, d_od: float, od_total: float,
                   power_W: float, lambda_nm: float, eta: float,
                   sigma_tech: float = 0.0) -> float:
    """Integration time for target SNR; inf if the technical-noise floor
    already exceeds the signal."""
    rel = abs(delta_T_over_T(d_od))
    denom = rel**2 - (target_snr*sigma_tech)**2
    if denom <= 0: return float('inf')
    rate = power_W/photon_energy_J(lambda_nm)*eta*transmission(od_total)
    return target_snr**2/(denom*rate)

def min_detectable_contrast(target_snr: float, od_sector: float, od_total: float,
                            power_W: float, lambda_nm: float, tau_s: float,
                            eta: float, sigma_tech: float = 0.0) -> float:
    """Inverts snr() at fixed (realistic) integration time to give the
    smallest sector-EIT contrast C = dOD/OD_sector detectable at
    target_snr, in the weak-signal (d_od << 1) linear regime where
    delta_T_over_T(d_od) ~= d_od. This is the epsilon_th of the room-
    temperature no-go plan (Sec. 0, "sup C_EIT(300K) < epsilon_th"): it
    must be fixed from the detection chain BEFORE the 300 K global
    optimization (Step 5) is run, not chosen after the fact to make the
    optimization result look conclusive.

    Returns inf if even a d_od = 1 (order-unity transmission change)
    signal cannot reach target_snr at this tau_s -- i.e. no achievable
    contrast is detectable under these conditions."""
    N = detected_photons(power_W, lambda_nm, tau_s, eta, od_total)
    if N <= 0:
        return float('inf')
    d_od_min = target_snr*np.sqrt(1.0/N + sigma_tech**2)
    if d_od_min >= 1.0:
        return float('inf')
    return d_od_min/od_sector
