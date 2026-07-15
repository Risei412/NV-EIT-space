"""
happacher_rate.py -- Gate 1: phonon-induced orbital hopping rate.

SOURCE (verified against arXiv:2302.00011 / PRL 131, 086904 (2023), SI):
  Eq. (21):  I(w_c,T,delta_perp) = int_0^{x_c} dx  x(x+x_D)[x^2+(x+x_D)^2]
                                   / [(e^x - 1)(1 - e^{-(x+x_D)})]
             x_c = hbar*w_c/(k_B T),  x = hbar*w/(k_B T),
             x_D = Delta_xy/(k_B T),  Delta_xy = E_x - E_y = 2*delta_perp
  Eq. (22):  Gamma_XY = alpha_ph * I(w_c = 60 meV, T, delta_perp) * T^5
  Eq. (23):  k_{Y->X} = k_{X->Y} = Gamma_XY        <-- UNIDIRECTIONAL jump rate
  alpha_ph = 1.70 +/- 0.08  Hz K^-5   (SI, "Phonon-coupling parameter")
  Cross-check (main text, NV#1): Gamma_XY = (1+eps)*gamma_NV*(1.10+/-0.05e-6 K^-5)*T^5,
  gamma_NV = (12.5 ns)^-1, eps << 1 for T ~ 10..100 K.

RATE CONVENTION (fixed for this project):
  k_orb     := Gamma_XY  (one-directional X<->Y jump rate, Eq. 23)
  Gamma_pop := 2 k_orb   (population-imbalance relaxation rate)
  gamma_oc  := k_orb / 2 (phonon damping of each optical coherence)
  Jump operators in the full Liouvillian:
      L_{X->Y} = sqrt(k_orb) |Y><X| (x) I_spin  and reverse.
  gamma_oc is NEVER added separately in the full Liouvillian (no double count);
  it is used only as the fast scale of the reduced kernel (D = I).

UNITS:
  delta_perp input in GHz (ordinary frequency).  Energy conversion uses h:
      Delta_xy [J] = h * (2 * delta_perp[Hz]).
  x_D = h*2*delta_perp / (k_B T).   Output rates in Hz (ordinary rate, s^-1).
  Conversion table (see also unit tests):
      1 GHz  = 1e9 Hz ; 1 meV/h = 241.7990 GHz ; 1 meV/hbar = 1519.267 rad/ns
      radiative: lifetime tau -> population decay k_r = 1/tau
                 -> optical HWHM = k_r/2 (amplitude damping).
"""
import numpy as np
from scipy.integrate import quad

H_PLANCK = 6.62607015e-34   # J s
HBAR = 1.054571817e-34      # J s
KB = 1.380649e-23           # J / K
MEV = 1.602176634e-22       # J
ALPHA_PH = 1.70             # Hz / K^5
ALPHA_PH_ERR = 0.08
OMEGA_C_MEV = 60.0          # hbar*omega_c = 60 meV (fixed, SI)


def _integrand(x, xD):
    # numerically stable using expm1: (e^x-1) = expm1(x); (1-e^{-(x+xD)}) = -expm1(-(x+xD))
    num = x * (x + xD) * (x * x + (x + xD) ** 2)
    den = np.expm1(x) * (-np.expm1(-(x + xD)))
    return num / den


def I_integral(T, delta_perp_GHz, omega_c_meV=OMEGA_C_MEV):
    """Dimensionless Eq. (21). delta_perp in GHz (ordinary frequency)."""
    if T <= 0:
        return 0.0
    x_c = omega_c_meV * MEV / (KB * T)
    xD = H_PLANCK * 2.0 * delta_perp_GHz * 1e9 / (KB * T)
    # integrand -> x*(x)*(2x^2)/(x * x) ~ 2x^2 as x->0 when xD=0 : finite; quad handles it,
    # but shift lower limit to tiny epsilon for xD=0 safety.
    val, err = quad(_integrand, 1e-12, x_c, args=(xD,), limit=400)
    return val


def k_orb(T, delta_perp_GHz, alpha_ph=ALPHA_PH):
    """Unidirectional X<->Y jump rate, Hz (Eq. 22-23)."""
    return alpha_ph * I_integral(T, delta_perp_GHz) * T ** 5


def gamma_pop(T, d, a=ALPHA_PH):
    return 2.0 * k_orb(T, d, a)


def gamma_oc(T, d, a=ALPHA_PH):
    return 0.5 * k_orb(T, d, a)


def k_orb_naive_T5(T, delta_perp_GHz, T_ref=30.0, alpha_ph=ALPHA_PH):
    """Naive extrapolation: freeze I at its low-T (T_ref) value, pure T^5."""
    I_ref = I_integral(T_ref, delta_perp_GHz)
    return alpha_ph * I_ref * T ** 5


def maintext_fit_NV1(T):
    """Main-text NV#1 cross-check: gamma_NV * 1.10e-6 K^-5 * T^5 (eps neglected), Hz."""
    gamma_NV = 1.0 / 12.5e-9
    return gamma_NV * 1.10e-6 * T ** 5
