# P0-7 — Symmetry-breaking rank switch audit

Date: 2026-08-10  
Parent: P0-6 unresolved rank-4 ODMR gate  
Original branch: `research-ci/nv4-symmetry-audit-20260810`  
Status: **witness PASS / generalization PASS / PRL novelty FAIL**

## 1. Exact leakage law

For symmetric `ms=0 <-> +/-1` transitions with a common even line profile `L`, write

\[
A(f,\mathbf B)=\sum_a w_a\frac{L(f-\gamma\mathbf n_a\cdot\mathbf B)+L(f+\gamma\mathbf n_a\cdot\mathbf B)}{2},
\qquad
w_a=\frac14+\epsilon_a,
\qquad
\sum_a\epsilon_a=0.
\]

The quadratic response tensor is

\[
Q^{(2)}=\frac13 I+E,\qquad E=\sum_a\epsilon_a\mathbf n_a\mathbf n_a^T.
\]

For the four tetrahedral NV axes,

\[
\boxed{\|E\|_F^2=\frac89\sum_a\epsilon_a^2}.
\]

Hence equal orientation weights are necessary and sufficient to eliminate anisotropic rank-2 response over all directions. Generic imbalance produces an `O(epsilon B^2)` directional term.

The information hierarchy becomes

- equal weights: directional amplitude `O(B^4)`, Fisher information `O(B^8)`;
- calibrated generic imbalance: amplitude `O(epsilon B^2)`, Fisher information `O(epsilon^2 B^4)`.

## 2. Crossover bound

For `[100]` versus `[110]`,

\[
\Delta A_2=-\frac{\gamma^2}{2}E_{xy}B^2L'',
\qquad
\Delta A_4=-\frac{\gamma^4}{216}B^4L''''.
\]

For a Lorentzian of FWHM `Gamma`,

\[
\frac{\|L''''\|_2}{\|L''\|_2}=\frac{\sqrt{105}}{(\Gamma/2)^2},
\]

so the quartic sector dominates when

\[
|E_{xy}|<\frac{\sqrt{105}}{27}\left(\frac{\gamma B}{\Gamma}\right)^2.
\]

For weights `(0.30,0.20,0.20,0.30)` and `Gamma=30 MHz`, the predicted crossover is about `0.449 mT`. Numerical exponents are `3.9996` for equal weights and `1.9976` for the imbalanced case.

## 3. Novelty verdict

The exact leakage law and crossover are retained as internal design rules, but the PRL candidate fails: weighted spherical-design mathematics is established, and deliberate orientation weighting using optical or microwave polarization is already established in NV vector magnetometry.

Do not rescue this as a new magnetometer claim. Retain it for systematics and experiment design: equal weighting protects rank-4-only access; unknown imbalance creates rank-2 leakage; calibrated imbalance can be used as a rank-2 resource.
