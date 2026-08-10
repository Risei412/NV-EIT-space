# P0-7 — Symmetry-breaking rank switch audit

Date: 2026-08-10  
Parent: P0-6 unresolved rank-4 ODMR gate  
Branch: `research-ci/nv4-symmetry-audit-20260810`  
Status: **witness PASS / generalization PASS / PRL novelty FAIL**

## 1. Question

If the four tetrahedral NV orientations are not exactly equally weighted, how does the leading directional tensor rank of an unresolved, spin-summed ODMR spectrum change, and can the rank-4 sector still dominate?

## 2. Exact theorem

For symmetric `m_s=0 <-> +/-1` transitions with a common even single-line profile `L`, write

\[
A(f,\mathbf B)=\sum_a w_a\frac{L(f-\gamma\mathbf n_a\!\cdot\!\mathbf B)+L(f+\gamma\mathbf n_a\!\cdot\!\mathbf B)}{2},
\]

with

\[
w_a=\frac14+\epsilon_a,\qquad \sum_a\epsilon_a=0.
\]

The quadratic response tensor is

\[
Q^{(2)}=\frac13 I+E,\qquad E=\sum_a\epsilon_a\,\mathbf n_a\mathbf n_a^T.
\]

For the four tetrahedral NV axes,

\[
\boxed{\|E\|_F^2=\frac89\sum_a\epsilon_a^2}.
\]

Therefore

\[
E=0\iff \epsilon_1=\epsilon_2=\epsilon_3=\epsilon_4=0.
\]

Equal orientation weights are necessary and sufficient to eliminate anisotropic rank-2 response over all directions.  Generic weight imbalance generates a rank-2 directional term at `O(epsilon B^2)`.

This changes the information hierarchy:

- equal weights: directional amplitude `O(B^4)`, hence background-dominated directional Fisher information `O(B^8)`;
- calibrated generic imbalance: directional amplitude `O(epsilon B^2)`, hence Fisher information `O(epsilon^2 B^4)`.

Thus controlled symmetry breaking can lower the accessible tensor rank from 4 to 2.

## 3. Matched-pair bound

For equal-|B| fields `[100]` and `[110]`, the weight-induced quadratic difference is

\[
\Delta A_2=-\frac{\gamma^2}{2}E_{xy}B^2L'',
\]

while the intrinsic equal-weight quartic difference is

\[
\Delta A_4=-\frac{\gamma^4}{216}B^4L''''.
\]

For a Lorentzian of FWHM `Gamma`,

\[
\frac{\|L''''\|_2}{\|L''\|_2}=\frac{\sqrt{105}}{(\Gamma/2)^2}.
\]

Therefore the quartic sector exceeds the rank-2 leakage in L2 norm when

\[
\boxed{|E_{xy}|<\frac{\sqrt{105}}{27}\left(\frac{\gamma B}{\Gamma}\right)^2}.
\]

Using `|E_xy| <= (2/3)||epsilon||_2` gives the sufficient condition

\[
\boxed{\|\epsilon\|_2<\frac{\sqrt{105}}{18}\left(\frac{\gamma B}{\Gamma}\right)^2}.
\]

For the stress pattern `(w1,w2,w3,w4)=(0.30,0.20,0.20,0.30)`, `||epsilon||_2=0.1`, `E_xy=1/15`.  With `Gamma=30 MHz`, the predicted crossover is

\[
B_c\simeq0.449\,\mathrm{mT}.
\]

Below this scale the calibrated rank-2 term dominates; above it the intrinsic quartic pair difference dominates in the asymptotic L2 comparison.

## 4. Numerical validation

A 30 MHz Lorentzian ODMR model including a common `14N` hyperfine triplet (`2.16 MHz`) gives, in the small-field regime,

- equal weights: `[100]-[110]` spectral RMS exponent `3.9996`,
- controlled imbalance: exponent `1.9976`.

Thus the expected rank switch `B^4 -> B^2` is numerically verified.

For the same stress pattern, the controlled imbalance increases the `[100]-[110]` RMS difference relative to equal weights by about `16.2x` at `0.1 mT`, `3.52x` at `0.2 mT`, and `1.21x` at `0.3 mT`.  Around and above the predicted crossover, the intrinsic quartic term becomes competitive and the simple low-field enhancement interpretation ceases to be monotone.

## 5. Novelty kill

### Known mathematics

The tetrahedral equal-weight four-axis set is an instance of spherical-design / equal-weight quadrature structure.  Weighted and approximate spherical designs are also established.  Therefore the statement that equal weights enforce low-order isotropy and weight perturbations create lower-order moment errors is not independent new mathematics.

The exact coefficient `8/9` and the Lorentzian crossover formula are useful NV-specialized diagnostics, but under Research CI novelty rules they are natural specializations, not a standalone PRL claim.

### Known NV physics

Controlled optical polarization and microwave polarization are already used to change the relative contribution of crystallographic NV orientations.  Polarization-assisted vector magnetometry explicitly exploits orientation-dependent optical excitation, and prior work further suppresses selected NV orientations with microwave polarization to reconstruct small magnetic fields without a bias field.

Low-field overlapping-orientation CW-ODMR and its field-orientation dependence have also been studied directly.

Hence the broad physical statement

> deliberately break the four-orientation symmetry to restore lower-rank directional information at low field

is already occupied in mechanism and purpose, even if the existing papers do not phrase it as a tensor-rank switch or Fisher-scaling law.

## 6. Gate verdict

- Exact leakage theorem: **PASS**.
- Crossover bound: **PASS**.
- Numerical `B^4 -> B^2` rank-switch witness: **PASS**.
- New low-field sensing capability from controlled orientation weighting: **FAIL, absorbed by established polarization-selective NV vector magnetometry**.
- `B^8 -> epsilon^2 B^4` Fisher-scaling interpretation: **useful explanatory residue but insufficient independent PRL novelty**.
- Final P0-7 PRL status: **FAIL**.

Do not rescue this candidate by renaming it as a new magnetometer.

## 7. Retained value

P0-7 should be retained as an internal design law:

1. equal orientation weighting is a symmetry resource that protects rank-4-only directional access;
2. unknown imbalance is a nuisance that generates rank-2 leakage;
3. known/calibrated imbalance is a resource that restores rank-2 low-field information;
4. the crossover between those regimes is quantitatively predictable.

This is useful for experimental design and for interpreting any future genuinely new tensor-rank protocol, but it is not promoted as an active PRL question.
