# P0-6 — Unresolved rank-4 ODMR gate

Date: 2026-08-10  
Parent: P0-5 tetrahedral quartic ODMR invariant  
Branch: `research-ci/nv4-symmetry-audit-20260810`  
Status: **PASS WITH NARROWING as a witness; PRL novelty HOLD**

## 1. Question

Can a four-orientation NV ensemble distinguish equal-magnitude magnetic-field directions even when its individual crystallographic ODMR families cannot be assigned or resolved, and if so what tensor rank carries the first directional information?

Matched pair: `[100]` versus `[110]` at the same `|B|`.

## 2. Exact unresolved-spectrum selection rule

For equal occupation of the four tetrahedral axes and symmetric `m_s=0 <-> +/-1` transition weights, let `L(f)` be any common even single-transition line profile.  In the low-field linear-Zeeman regime the orientation-unassigned normalized spectrum has expansion

\[
A(f,\mathbf B)
=L+\frac{\gamma^2B^2}{6}L''
+\gamma^4\left[\frac{(B^2)^2}{72}-\frac{S_4}{108}\right]L''''+O(B^6),
\]

with

\[
S_4=B_x^4+B_y^4+B_z^4.
\]

Therefore:

- rank 1 is forbidden,
- rank 2 is isotropic,
- rank 3 is forbidden by field-even spin summation,
- rank 4 is the first direction-dependent term.

For background-dominated independent spectral noise, `dA/dtheta = O(B^4)`, so classical directional Fisher information scales as `O(B^8)`.

The numerical ideal-spectrum audit gives power `3.959` for the `[100]-[110]` RMS spectral difference and `7.918` for its squared distance.

## 3. Realistic unresolved witness

Spectrum model:

- exact spin-1 ground-state diagonalization,
- four tetrahedral NV axes,
- `14N` hyperfine triplet, `A=2.16 MHz`,
- Lorentzian FWHM `30 MHz`,
- total ODMR contrast `3%`,
- 601 frequency samples across 300 MHz,
- deterministic orientation-dependent transverse-strain stress values `1.59, 2.24, 4.53, 6.71 MHz`.

The wrong `[110]` hypothesis is allowed to optimize nuisance parameters: field magnitude, linewidth, overall contrast, baseline, common frequency shift, common `+/-` transition imbalance up to `20%`, and each orientation population within `0.20..0.30` (plus/minus `20%` around equal `0.25`).

At `B=0.5 mT`, both `[100]` and `[110]` spectra have a single visible ODMR dip under the stated 30 MHz linewidth.  Even after the nuisance fit, the wrong-direction residual is

- RMS `1.378e-4` of normalized fluorescence,
- maximum residual `6.39e-4`,
- best wrong-fit orientation weights approximately `(0.284,0.200,0.216,0.300)`.

For 601 statistically independent points, the Euclidean template separation reaches `d'=5` when per-point normalized-fluorescence noise is approximately `6.76e-4` or smaller.  This is a detectability benchmark, not a full experimental sensitivity claim.

At `0.3 mT` the corresponding noise requirement is about `1.02e-4`; at `0.7 mT` it relaxes to about `2.45e-3`.

## 4. Robustness boundary

The witness survives the imposed hyperfine, several-MHz strain, moderate line-amplitude asymmetry and plus/minus 20% orientation-population nuisance.

However, if the four orientation weights are allowed to vary essentially freely, a wrong direction can mimic the target spectrum to near numerical precision by adopting extreme orientation populations.  Therefore near-equal or independently calibrated orientation populations are a load-bearing resource.

This kills any claim of a generic calibration-free vector magnetometer.

## 5. External novelty audit

### Established and therefore killed

1. Four-axis NV vector magnetometry from resolved resonance projections is established.
2. Bias-free and polarization-assisted schemes that label the four NV families are established.
3. Full ODMR spectra and machine-learning inference are already used to infer magnetic quantities; a 2026 conference contribution publicly claims full-vector reconstruction directly from raw CW-ODMR spectra with strain, temperature and contrast fluctuations.
4. Low-field/overlapping-transition ODMR and ensemble line-shape modeling are established.
5. Complete transition overlap has been intentionally engineered for sensitivity enhancement.

Thus the broad claim `raw unresolved ODMR contains vector information` is not novel enough.

### Residue not directly matched in the targeted search

The specific structural statement

> in an equal four-orientation single-crystal NV ensemble with orientation-unassigned, spin-summed symmetric ODMR readout, directional information is symmetry-forbidden through quadratic order and first appears as a tetrahedral quartic spectral distortion, enforcing B^4 signal and B^8 Fisher-information scaling at low field

was not directly matched in the targeted primary-literature search.

This is more valuable as a **selection/no-go law for low-field orientation information** than as a new magnetometer architecture.

## 6. Gate verdict

- Minimal matched-pair witness: **PASS**.
- Realistic unresolved-spectrum witness: **PASS WITH NARROWING**.
- Broad new-sensor claim: **FAIL / absorbed by existing raw-spectrum and vector-inference work**.
- Structural rank-4 / B^8 law: **HOLD for PRL novelty**.

The next decisive gate should test whether the rank-4/B^8 law remains exact or becomes a controlled inequality under realistic but calibrated orientation-weight, microwave-polarization and strain distributions.  A theorem or bound showing when lower-rank leakage is absent, and how leakage scales when symmetry is weakly broken, would create an independent claim.  Without that, the result remains an elegant explanation of an existing sensing possibility rather than a PRL-level new capability.
