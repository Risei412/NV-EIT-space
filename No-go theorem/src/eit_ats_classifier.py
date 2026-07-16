"""eit_ats_classifier.py -- Gate 3: connect delta_chi_S to an observed
absorption spectrum and objectively discern EIT from Autler-Townes
splitting (ATS), following Anisimov, Dowling & Sanders, PRL 107, 163604
(2011), with the AIC/AICc thresholds fixed in SiV_SnV_phonon_AIC_parameters.md.

Pipeline: generic (dimension-agnostic) excited-manifold response, identical
convention to nv_model.response() (Acut = Re(S1), dA = -Re(dXi),
C = dA/Acut, per convention_table.md), applied to whichever reduced-kernel
Hamiltonian (NV: nv_reduced_kernel.H_3E; group-IV: group_iv_model.H_groupIV)
is passed in -- this is delta_chi_S connected all the way to a probe-
frequency absorption spectrum A(z) = Acut(z) - dA(z).

Model comparison (Anisimov et al. 2011):
  A_EIT(d) = C+^2/(g+^2+d^2) - C-^2/(g-^2+d^2)
  A_ATS(d) = C^2 [ 1/(g^2+(d-d0)^2) + 1/(g^2+(d+d0)^2) ]
AIC_i = -2 ln L_i + 2 K_i  (Gaussian likelihood, unknown variance:
  AIC = N ln(RSS/N) + 2K); Delta_AIC = AIC_ATS - AIC_EIT.
Robust decision gate (fixed, not a universal physical constant):
  Delta_AIC >= +6 -> robust EIT ; <= -6 -> robust ATS ; else inconclusive.
"""
import numpy as np
from scipy.optimize import curve_fit

TWOPI = 2*np.pi

def response_generic(H, dp, dc, z, gamma, Oc=1.0, gg=6.3e-5):
    """Dimension-agnostic copy of nv_model.response's convention (Acut, dA, C),
    generalized to any excited-manifold dimension N=H.shape[0]."""
    N = H.shape[0]
    beta = (TWOPI*Oc)**2/4; geff = 2*gg + 2e-6
    G = np.linalg.inv(gamma*np.eye(N) + 1j*TWOPI*(H - z*np.eye(N)))
    K12 = np.vdot(dp, G@dc); K21 = np.vdot(dc, G@dp)
    S1 = np.vdot(dp, G@dp); S2 = np.vdot(dc, G@dc)
    den = geff + beta*S2
    dXi = -beta*K12*K21/den
    Acut = float(np.real(S1)); dA = float(-np.real(dXi))
    C = dA/Acut if abs(Acut) > 1e-300 else np.nan
    return dict(C=C, Acut=Acut, dA=dA, K12=K12, K21=K21, S1=S1, S2=S2, den=den, z=z)

def spectrum(H, dp, dc, gamma, zs, Oc=1.0, gg=6.3e-5):
    Acut = np.empty_like(zs); dA = np.empty_like(zs)
    for i, z in enumerate(zs):
        r = response_generic(H, dp, dc, z, gamma, Oc, gg)
        Acut[i] = r['Acut']; dA[i] = r['dA']
    return Acut - dA  # full (control-on) absorption line shape

# ---- Anisimov et al. (2011) model-comparison lineshapes ----
def A_EIT(d, Cp, gp, Cm, gm):
    return Cp**2/(gp**2 + d**2) - Cm**2/(gm**2 + d**2)

def A_ATS(d, C, g, d0):
    return C**2*(1.0/(g**2 + (d - d0)**2) + 1.0/(g**2 + (d + d0)**2))

def _aic(rss, n, k):
    rss = max(rss, 1e-300)
    return n*np.log(rss/n) + 2*k

def _aicc(aic, n, k):
    if n - k - 1 <= 0: return np.inf
    return aic + 2*k*(k+1)/(n - k - 1)

def classify(delta, A, robust_threshold=6.0, use_aicc=None):
    """Fit A_EIT and A_ATS to (delta, A); return dict with AIC, Delta_AIC,
    Akaike weight and a fixed-threshold verdict (Sec. 10.2 of
    SiV_SnV_phonon_AIC_parameters.md: |Delta_AIC|>=6 => robust)."""
    n = len(delta)
    A0 = float(np.max(np.abs(A))) or 1.0
    g0 = float((delta.max() - delta.min())/10.0) or 1.0
    p_eit0 = (A0, g0, A0*0.5, g0*0.3)
    p_ats0 = (A0**0.5, g0, g0)
    try:
        p_eit, _ = curve_fit(A_EIT, delta, A, p0=p_eit0, maxfev=20000)
        rss_eit = float(np.sum((A - A_EIT(delta, *p_eit))**2))
    except Exception:
        p_eit, rss_eit = p_eit0, np.inf
    try:
        p_ats, _ = curve_fit(A_ATS, delta, A, p0=p_ats0, maxfev=20000)
        rss_ats = float(np.sum((A - A_ATS(delta, *p_ats))**2))
    except Exception:
        p_ats, rss_ats = p_ats0, np.inf
    k_eit, k_ats = 4, 3
    if use_aicc is None: use_aicc = (n/max(k_eit, k_ats) < 40)
    aic_eit = _aic(rss_eit, n, k_eit); aic_ats = _aic(rss_ats, n, k_ats)
    if use_aicc:
        aic_eit = _aicc(aic_eit, n, k_eit); aic_ats = _aicc(aic_ats, n, k_ats)
    dAIC = aic_ats - aic_eit
    from scipy.special import expit
    w_eit = float(expit(dAIC/2.0))
    if dAIC >= robust_threshold: verdict = 'robust EIT'
    elif dAIC <= -robust_threshold: verdict = 'robust ATS'
    else: verdict = 'inconclusive'
    return dict(rss_eit=rss_eit, rss_ats=rss_ats, params_eit=list(p_eit),
                params_ats=list(p_ats), aic_eit=aic_eit, aic_ats=aic_ats,
                delta_aic=float(dAIC), weight_eit=float(w_eit),
                used_aicc=bool(use_aicc), verdict=verdict)
