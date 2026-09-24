# P0-6 — Unresolved rank-4 ODMR gate

Date: 2026-08-10  
Parent: P0-5 tetrahedral quartic ODMR invariant  
Original branch: `research-ci/nv4-symmetry-audit-20260810`  
Status: **PASS WITH NARROWING as a witness; PRL novelty HOLD**

## 1. Exact unresolved-spectrum selection rule

For equal occupation of the four tetrahedral axes and symmetric `ms=0 <-> +/-1` transition weights, let `L(f)` be any common even single-transition line profile. In the low-field linear-Zeeman regime,

\[
A(f,\mathbf B)=L+\frac{\gamma^2B^2}{6}L''+\gamma^4\left[\frac{(B^2)^2}{72}-\frac{S_4}{108}\right]L''''+O(B^6),
\]

where `S4=Bx^4+By^4+Bz^4`.

Therefore rank 1 is forbidden, rank 2 is isotropic, rank 3 is forbidden by field-even spin summation, and rank 4 is the first direction-dependent term. For background-dominated independent spectral noise, `dA/dtheta=O(B^4)`, so classical directional Fisher information scales as `O(B^8)`.

The numerical ideal-spectrum audit gives power `3.959` for the `[100]-[110]` RMS spectral difference and `7.918` for its squared distance.

## 2. Realistic unresolved witness

The witness uses exact spin-1 diagonalization, four tetrahedral axes, `14N` hyperfine triplet (`2.16 MHz`), Lorentzian FWHM `30 MHz`, total ODMR contrast `3%`, 601 samples, and orientation-dependent transverse-strain stress values `1.59, 2.24, 4.53, 6.71 MHz`.

The wrong `[110]` hypothesis may optimize field magnitude, linewidth, contrast, baseline, common frequency shift, `+/-` transition imbalance up to 20%, and orientation populations within `0.20..0.30`.

At `B=0.5 mT`, both spectra have a single visible ODMR dip. After the nuisance fit, the wrong-direction residual is RMS `1.378e-4`, maximum `6.39e-4`, with fitted weights approximately `(0.284,0.200,0.216,0.300)`.

## 3. Robustness boundary

The witness survives hyperfine, several-MHz strain, moderate line-amplitude asymmetry, and ±20% orientation-population nuisance. If orientation weights are allowed to vary essentially freely, however, a wrong direction can mimic the target to near numerical precision. Near-equal or independently calibrated orientation populations are therefore load-bearing.

## 4. Gate verdict

- Minimal matched-pair witness: **PASS**.
- Realistic unresolved-spectrum witness: **PASS WITH NARROWING**.
- Broad new-sensor claim: **FAIL / absorbed by existing raw-spectrum and vector-inference work**.
- Structural rank-4 / B^8 law: **HOLD for PRL novelty**.

The retained claim is a structural low-field selection law, not a generic calibration-free magnetometer architecture.
