# ASSUMPTIONS.md — regime and scope

Every claim in `../CLAIMS.md` is conditional on these. Quoting a result
outside this regime is a misquotation.

## Standing regime for the classification

1. **Finite-dimensional.** All spaces are finite-dimensional complex
   inner-product spaces.
2. **Markovian.** GKSL dynamics, time-independent generator.
3. **Weak probe.** Linear response only.
4. **Stationary rotating frame.**
5. **Passive optical response, no gain.**
6. **A well-defined source and readout channel.**
7. **A physically specified long-lived coherence sector `S`.** The sector cut
   depends on this choice (T1 Lemma L1.2) — it is not canonical.
8. **Regular observation frequencies**, away from poles.

The trichotomy in `ν` is stated in this regime and **nowhere else**.

## Explicitly outside scope

- Strong-probe nonlinearities
- Non-Markovian memory
- Continuum or propagation instabilities
- Gain media
- Many-body collective effects
- Time-dependent feedback or Floquet driving

## Additional assumptions per result

### Matched-readout floor (`CLAIMS.md` C3)

- Passive coherence block `A(δ) = Γ + i(H − δI)` with `Γ = Γ† ≥ 0`.
- `Γ > 0` gives strict positivity; `Γ ≥ 0` gives the non-strict floor.
- "Matched" means source = readout. The identity says nothing about unmatched
  readout, where exact zeros do occur.

### Group-IV (SiV/SnV) material independence (Gate C)

**Load-bearing and known to be schematic.** The phonon normalization and
dipole geometry are schematic: single rate scale, orthogonal orbital basis
legs, `|M0|` normalization. Constants are recorded in
`../calculations/numerics/New no-go theory/GateC_material_independence/SiV_SnV_phonon_AIC_parameters.md`
but were **not refit**. Integer exponents are unaffected; group-IV curves must
not be read as quantitative predictions of absolute response
(`../NON_CLAIMS.md` N4).

### Gate E 3-mode conditional PASS (`CLAIMS.md` C7)

The PASS carries design assumptions that are **not demonstrated**:

- 15–45 MHz tunable engineered loss
- coupling-preserving independent loss sweep
- −190 dBm-class amplification and calibration
- sign-switching transients
- a measured feedthrough drift spectrum

### NV campaigns

- Candidate Λ geometry with small `B_z` (0.02 T) so the (m_s = −1, m_s = +1)
  branch assignment is unambiguous; strain `d = 1.683` fixed across the gate
  scripts.
- Phonon rate variants are treated as a model class, not a single model; the
  threshold bands quoted are quantiles over that class.
