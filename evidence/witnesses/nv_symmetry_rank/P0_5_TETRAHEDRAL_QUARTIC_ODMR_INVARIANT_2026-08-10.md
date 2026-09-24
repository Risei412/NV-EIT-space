# P0-5 — Tetrahedral quartic ODMR invariant

Date: 2026-08-10  
Parent: P0-4 tensor-rank generation  
Original branch: `research-ci/nv4-symmetry-audit-20260810`  
Status: **clean internal witness / PRL HOLD**

## 1. Observable

For one NV orientation `a`, use the ground spin-1 Hamiltonian

\[
H_a=D S_z^2+\gamma\mathbf B_a\cdot\mathbf S.
\]

Let `lambda_0,a` denote the nondegenerate eigenvalue continuously connected to `ms=0`. The sum of the two ODMR transition frequencies is

\[
\Omega_{\Sigma,a}=\operatorname{Tr}H_a-3\lambda_{0,a}=2D-3\lambda_{0,a}.
\]

Define

\[
\bar\Omega_\Sigma=\frac14\sum_{a=1}^4\Omega_{\Sigma,a},
\qquad
f_{\rm cent}=\frac12\bar\Omega_\Sigma.
\]

This quantity is invariant under permutation of orientations and exchange of the two transitions inside each orientation.

## 2. Small-field result

With `q=gamma(n_a.B)` and `p^2=gamma^2(B^2-(n_a.B)^2)`, the `ms=0` branch is

\[
\lambda_0=-\frac{p^2}{D}+\frac{p^4-p^2q^2}{D^3}+O(B^6).
\]

Tetrahedral averaging gives

\[
\boxed{\bar\Omega_\Sigma=2D+\frac{2\gamma^2}{D}B^2-\frac{2\gamma^4}{D^3}(B^2)^2+\frac{4\gamma^4}{3D^3}(B_x^4+B_y^4+B_z^4)+O(B^6)}
\]

and

\[
\boxed{f_{\rm cent}=D+\frac{\gamma^2}{D}B^2-\frac{\gamma^4}{D^3}(B^2)^2+\frac{2\gamma^4}{3D^3}(B_x^4+B_y^4+B_z^4)+O(B^6).}
\]

Thus rank 1 is forbidden, rank 2 is isotropic, rank 3 is absent in the field-even spectral sum, and rank 4 is the first direction-dependent tensor.

## 3. Matched-pair witness

For equal-magnitude fields `[100]` and `[110]`, `B^2` is identical and `Bx By Bz=0` for both, but

\[
S_4[100]=B^4,\qquad S_4[110]=\frac12B^4.
\]

Hence

\[
\bar\Omega_\Sigma[100]-\bar\Omega_\Sigma[110]=\frac{2\gamma^4}{3D^3}B^4+O(B^6).
\]

Direct diagonalization agrees with the analytic coefficient `2 gamma^4/(3D^3)=1.7269058e4 GHz/T^4` in the small-field regime.

## 4. Claim status

The mathematical witness is verified. The bare quartic formula remains at natural-corollary risk and is not promoted as a standalone PRL claim. Its value is as the clean, basis-invariant foundation for the unresolved-spectrum gate P0-6.
