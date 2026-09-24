"""P0-7: symmetry-breaking rank switch for unresolved four-orientation NV ODMR."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

D = 2.877
GAMMA = 28.02495164
FWHM = 0.030
A14 = 0.00216
CONTRAST = 0.03
AXES = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) / np.sqrt(3)
U100 = np.array([1.0, 0, 0])
U110 = np.array([1.0, 1.0, 0]) / np.sqrt(2)
F = np.linspace(D - 0.15, D + 0.15, 601)
REPO = Path(__file__).resolve().parents[3]


def leakage_tensor(eps):
    return sum(e * np.outer(n, n) for e, n in zip(np.asarray(eps, float), AXES))


def lorentz(f, f0, w):
    g = w / 2
    return g * g / ((f - f0) ** 2 + g * g)


def spectrum(B, u, w):
    y = np.zeros_like(F)
    for wa, n in zip(w, AXES):
        d = GAMMA * B * np.dot(n, u)
        for s in (-1, 1):
            for m in (-1, 0, 1):
                y += wa / 6 * lorentz(F, D + s * d + m * A14, FWHM)
    y /= y.max()
    return 1 - CONTRAST * y


def rms_pair(B, w):
    return float(np.sqrt(np.mean((spectrum(B, U100, w) - spectrum(B, U110, w)) ** 2)))


def main():
    eps = np.array([0.05, -0.05, -0.05, 0.05])
    w = 0.25 + eps
    E = leakage_tensor(eps)
    epsnorm = float(np.linalg.norm(eps))
    Efnorm = float(np.linalg.norm(E))
    Exy = float(E[0, 1])
    Bc = (FWHM / GAMMA) * np.sqrt(27 * abs(Exy) / np.sqrt(105))
    Bc_bound = (FWHM / GAMMA) * np.sqrt(18 * epsnorm / np.sqrt(105))
    small = np.array([5e-6, 7.5e-6, 1e-5, 1.5e-5, 2e-5, 3e-5])
    req = np.array([rms_pair(B, np.ones(4) / 4) for B in small])
    rim = np.array([rms_pair(B, w) for B in small])
    rows = []
    for BmT in (0.1, 0.2, 0.3, 0.5):
        B = BmT * 1e-3
        a = rms_pair(B, np.ones(4) / 4)
        b = rms_pair(B, w)
        rows.append(dict(B_mT=BmT, equal_weight_RMS=a, controlled_imbalance_RMS=b, ratio_imbalance_to_equal=b / a))
    ratios = {}
    for BmT in (0.2, 0.3, 0.5, 0.7, 1.0):
        B = BmT * 1e-3
        ratios[str(BmT)] = float((np.sqrt(105) / 27) * (GAMMA * B / FWHM) ** 2 / abs(Exy))
    out = dict(
        run="RCI-20260810-P0-7",
        theorem=dict(
            spectrum="A=sum_a w_a [L(f-delta_a)+L(f+delta_a)]/2, delta_a=gamma n_a.B",
            weights="w_a=1/4+epsilon_a, sum epsilon_a=0",
            quadratic_tensor="Q2=I/3+E, E=sum_a epsilon_a n_a n_a^T",
            exact_norm_identity="||E||_F^2=(8/9)||epsilon||_2^2",
            necessary_sufficient="E=0 iff all epsilon_a=0",
            rank_switch="equal weights -> first directional rank 4; generic known imbalance -> rank 2 at O(epsilon B^2)",
            fisher_scaling="equal weights F_theta=O(B^8); calibrated generic imbalance F_theta=O(epsilon^2 B^4)",
        ),
        matched_pair_100_110=dict(
            quadratic_leakage="Delta A2=-(gamma^2/2) E_xy B^2 L''",
            intrinsic_quartic="Delta A4=-(gamma^4/216) B^4 L''''",
            lorentzian_norm_ratio="||L''''||_2/||L''||_2=sqrt(105)/(Gamma/2)^2",
            quartic_visibility="|E_xy| < (sqrt(105)/27)(gamma B/Gamma)^2",
            worst_case_sufficient="||epsilon||_2 < (sqrt(105)/18)(gamma B/Gamma)^2",
        ),
        stress_pattern=dict(
            epsilon=eps.tolist(),
            weights=w.tolist(),
            epsilon_l2=epsnorm,
            E=E.tolist(),
            E_fro=Efnorm,
            norm_identity_ratio=(Efnorm**2) / (epsnorm**2),
            E_xy=Exy,
            crossover_B_mT=float(Bc * 1e3),
            crossover_from_eps_bound_mT=float(Bc_bound * 1e3),
            quartic_to_rank2_L2_ratio=ratios,
        ),
        numerical_validation=dict(
            linewidth_MHz=FWHM * 1e3,
            hyperfine_MHz=A14 * 1e3,
            equal_weight_spectral_difference_exponent=float(np.polyfit(np.log(small), np.log(req), 1)[0]),
            imbalanced_spectral_difference_exponent=float(np.polyfit(np.log(small), np.log(rim), 1)[0]),
            benchmarks=rows,
        ),
    )
    path = REPO / "results" / "tables" / "p0_7_symmetry_breaking_rank_switch.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
