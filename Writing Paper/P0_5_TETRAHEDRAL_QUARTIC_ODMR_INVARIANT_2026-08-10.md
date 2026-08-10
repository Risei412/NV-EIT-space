# P0-5 — Tetrahedral quartic ODMR invariant

Date: 2026-08-10  
Parent: P0-4 tensor-rank generation  
Branch: `research-ci/nv4-symmetry-audit-20260810`  
Status: **clean internal witness / PRL HOLD**

## 1. Why P0-5 is needed

P0-4 exposed a rank hierarchy, but a one-sided `ms=+1` or `ms=-1` control channel at exactly zero magnetic field retains a labeling subtlety because the `ms=±1` ground doublet is degenerate.  That makes it unsuitable as the cleanest intrinsic tensor-rank witness.

P0-5 removes that ambiguity completely by using only a basis-invariant spectral quantity: the sum of both physical ODMR transition frequencies for each NV orientation, then the average over the four tetrahedral orientations.

## 2. Observable

For one NV orientation `a`, with local NV axis `n_a`, use the ground spin-1 Hamiltonian

\[
H_a=D S_z^2+\gamma\mathbf B_a\cdot\mathbf S.
\]

Let `lambda_0,a` denote the nondegenerate eigenvalue continuously connected to `ms=0` in the small-field regime.  If the other two eigenvalues are `lambda_1,a` and `lambda_2,a`, the sum of the two ODMR transition frequencies is

\[
\Omega_{\Sigma,a}=(\lambda_1-\lambda_0)+(\lambda_2-\lambda_0)
=\operatorname{Tr}H_a-3\lambda_{0,a}
=2D-3\lambda_{0,a}.
\]

No labeling of the degenerate `ms=±1` doublet is required.

Define

\[
\bar\Omega_\Sigma=\frac14\sum_{a=1}^4\Omega_{\Sigma,a}.
\]

The mean frequency of all eight orientation-times-spin ODMR lines is simply

\[
f_{\rm cent}=\frac12\bar\Omega_\Sigma.
\]

This centroid is permutation invariant under exchange of the four orientations and under exchange of the two transitions within an orientation.

## 3. Local small-field expansion

Write

\[
q=\gamma(\mathbf n_a\cdot\mathbf B),
\qquad
p^2=\gamma^2\left(B^2-(\mathbf n_a\cdot\mathbf B)^2\right).
\]

The exact characteristic polynomial is

\[
\lambda^3-2D\lambda^2+(D^2-p^2-q^2)\lambda+Dp^2=0.
\]

The root connected to `ms=0` is

\[
\lambda_0=-\frac{p^2}{D}
+\frac{p^4-p^2q^2}{D^3}+O(B^6).
\]

## 4. Tetrahedral averaging

For the four axes

\[
\mathbf n_a\in\frac1{\sqrt3}
\{(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)\},
\]

\[
\langle(\mathbf n_a\cdot\mathbf B)^2\rangle=\frac{B^2}{3},
\]

and

\[
\langle(\mathbf n_a\cdot\mathbf B)^4\rangle
=\frac13(B^2)^2-rac29(B_x^4+B_y^4+B_z^4).
\]

Substitution gives the closed small-field formula

\[
\boxed{
\bar\Omega_\Sigma
=2D+\frac{2\gamma^2}{D}B^2
-\frac{2\gamma^4}{D^3}(B^2)^2
+\frac{4\gamma^4}{3D^3}(B_x^4+B_y^4+B_z^4)
+O(B^6)
}
\]

and therefore

\[
\boxed{
f_{\rm cent}
=D+\frac{\gamma^2}{D}B^2
-\frac{\gamma^4}{D^3}(B^2)^2
+\frac{2\gamma^4}{3D^3}(B_x^4+B_y^4+B_z^4)
+O(B^6).}
\]

## 5. Tensor-rank selection

For this observable:

- rank 1 is forbidden;
- rank 2 survives only as the isotropic scalar `B^2`;
- rank 3 is absent because the spectral sum is field-even;
- rank 4 is the first direction-dependent magnetic tensor.

Thus the four-orientation ensemble converts direction information into a quartic invariant even after orientation labels are discarded at the level of the observable.

## 6. Matched-pair witness

Choose two fields with equal magnitude:

\[
\mathbf B_{100}=B(1,0,0),
\qquad
\mathbf B_{110}=\frac{B}{\sqrt2}(1,1,0).
\]

They have the same lower aggregate invariant `B^2`, and both have tetrahedral cubic invariant `Bx By Bz=0`.  But

\[
S_4=B_x^4+B_y^4+B_z^4
\]

is

\[
S_4[100]=B^4,
\qquad
S_4[110]=\frac12B^4.
\]

Hence

\[
\boxed{
\bar\Omega_\Sigma[100]-\bar\Omega_\Sigma[110]
=\frac{2\gamma^4}{3D^3}B^4+O(B^6)
}
\]

and

\[
\boxed{
f_{\rm cent}[100]-f_{\rm cent}[110]
=\frac{\gamma^4}{3D^3}B^4+O(B^6).}
\]

With `D=2.877 GHz` and `gamma=28.02495164 GHz/T`, direct diagonalization gives for the transition-sum difference:

| |B| | numerical Delta Omega | Delta Omega / B^4 |
|---:|---:|---:|
| 0.1 mT | 0.00173 Hz | 1.7266e4 GHz/T^4 |
| 0.5 mT | 1.079 Hz | 1.7267e4 GHz/T^4 |
| 1 mT | 17.263 Hz | 1.7263e4 GHz/T^4 |
| 2 mT | 275.886 Hz | 1.7243e4 GHz/T^4 |
| 5 mT | 10.691 kHz | 1.7106e4 GHz/T^4 |

The analytic coefficient is

\[
\frac{2\gamma^4}{3D^3}=1.7269058\times10^4\;{\rm GHz/T^4}.
\]

The all-eight-line centroid difference is one half of these frequency differences: about 8.63 Hz at 1 mT, 138 Hz at 2 mT, and 5.35 kHz at 5 mT.

## 7. Novelty audit

The mathematical pieces alone are not a PRL claim:

- third-rank tetrahedral/octupolar order is standard;
- four NV orientations are standard for vector magnetometry;
- nonlinear orientation dependence from the spin-1 Hamiltonian is standard;
- recent NV vector-magnetometry work explicitly uses both `ms=0↔±1` transition frequencies to recover magnetic-field magnitude and angle after assigning each orientation.

The targeted literature search did **not** find a direct primary-paper match for using the permutation-invariant all-eight-line centroid as a quartic directional observable.  That absence is not sufficient for novelty because the formula is a natural consequence of the standard Hamiltonian plus tetrahedral averaging.

Therefore:

- mathematical witness: **VERIFIED**;
- novelty of bare quartic formula: **weak / natural-corollary risk**;
- PRL status: **HOLD**, not PASS.

## 8. Decisive next gate: unresolved-spectrum capability

The route with real leverage is operational, not algebraic.

Construct a realistic four-orientation ODMR spectrum with linewidth, contrast, strain, hyperfine and orientation-population imbalance, then ask whether an **orientation-unassigned spectral statistic** can recover the quartic invariant when conventional peak assignment fails or requires an external bias field.

The clean capability target is:

> Two equal-|B| fields that are indistinguishable to the lower aggregate moments of an unresolved spectrum are separated by a permutation-invariant rank-4 ODMR statistic, without assigning the four NV orientations.

The first matched pair is `[100]` versus `[110]`.  If the quartic separation survives realistic linewidth and nuisance parameters and can be extracted without peak assignment, the central claim becomes a sensing capability rather than a restatement of tetrahedral invariant theory.

That is the appropriate next Gate.
