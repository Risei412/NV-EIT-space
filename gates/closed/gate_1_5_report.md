# Gates 1–5 execution report (SIMULATION_PLAN.md)

Run date: 2026-07-16. All five PRL verification gates of `SIMULATION_PLAN.md`
were implemented and executed at the fixed candidate point
(T = 70 K, Bx = 0.23226 T, Bz = 0.005 T, branch J0 = 3, Ωc = 0.1 GHz,
Y/Y polarization, control on ms = +1). Code: `src/gate{1..5}_*.py`,
`src/signal_chain.py`; entry point `scripts/reproduce_prl_figures.py`;
seeds fixed (20260716). **All gates pass.** Test suite: 19/19
(`python tests/test_core.py`).

| Gate | Verdict | Key numbers |
|---|---|---|
| 1 EIT/ATS statistics | **pass** | ΔAIC(ATS−EIT) = +1902; 100 % robust-EIT over 200 bootstraps at every noise level up to σ = 0.2·depth; no window/init flip; EIT→ATS crossover continuous (EIT ≤ 0.63 GHz, inconclusive ≈ 1.5 GHz, ATS ≥ 3.6 GHz) |
| 2 full vs reduced | **pass** | contrast ratio full/reduced = 1.001 (base), 1.004 (+ISC), 0.866 (+hyperfine avg); FWHM 0.33 → 0.46 MHz under hyperfine; no sign flip in any toggle |
| 3 signal chain | **pass** | OD_sector(1 ppm, 0.5 mm) ≈ 7.9, ΔOD ≈ 0.11; τ(SNR = 5) sub-second for all contrast targets 1 %/0.1 %/0.01 % at 0.01–1 ppm; margin ≥ 10⁴ against cross-section uncertainty |
| 4 threshold bands | **pass** | T_1% = 71.5 K [67.2, 75.8]; T_0.1% = 81.1 K [76.1, 86.3]; T_0.01% = 90.5 K [84.7, 96.4]; T_sign = 101.9 K [95.1, 108.9] (68 %; n = 500, all valid); dominant parameter everywhere: orbital-hopping rate scale |
| 5 ensemble | **pass** | single 1.36 % → unselected low-density 0.057 %, high-density 0.013 % (washed out); orientation + spectral selection alone 0.095 %; **+ field shimmed to 0.1 G: 0.26 %, FWHM 1.1 MHz** (detectable in ms per Gate 3) |

## Gate 1 — model comparison at the candidate (robust EIT)

`gate1_candidate_aic_bootstrap.py`. The probe-scan spectrum (full
z-dependence, den → geff + i2πδ₂ + βS₂) was fit with four models (EIT,
ATS, Lorentzian, Fano — all with explicit constant background, since the
one-photon line, γ ≈ 72 GHz, is flat across the MHz window). The fixed
Anisimov–Dowling–Sanders gate |ΔAIC| ≥ 6 gives **robust EIT** with
ΔAIC ≈ 1.9 × 10³, stable under: 200 noise bootstraps at σ/depth ∈
{0.02, 0.05, 0.1, 0.2} (100 % robust EIT each), window factors 0.5–2×,
50 randomized initial values (verdict unchanged, ΔAIC spread 90 ≪ 1902).
The control sweep shows a physically ordered, monotonic EIT → inconclusive
→ ATS transition, with the candidate deep in the EIT region.

Caveat recorded for the manuscript: among the *four* models the noiseless
spectrum slightly prefers Fano/Lorentzian over the ADS-EIT form (Akaike
weights 0.47/0.40/0.13) — the interference dip carries the small asymmetry
expected from the complex kernel phase arg(K₁₂K₂₁). This does not affect
the EIT-vs-ATS discrimination (the PRL claim), but "Fano-compatible
interference window" is the accurate fine-grained description.

## Gate 2 — full Lindblad confirms the reduced kernel

`gate2_candidate_full_vs_reduced.py`. The validated 9-level pipeline of
`bperp_full_validation.py` was generalized to the candidate configuration
and a genuine probe scan, plus toggles: ISC/singlet shelving (10th level;
Γ_ISC(±1) = 0.098, Γ_ISC(0) = 0.011, singlet 371 ns), and ¹⁴N hyperfine
(secular; mI conserved by every jump operator, so the ensemble is the exact
average of three shifted models). Results: the full CPTP model reproduces
the reduced δχ_S spectrum to 0.1 % at the candidate; ISC shifts the
contrast by +0.4 %; hyperfine averaging retains 86 % of the contrast with
moderate broadening (0.33 → 0.46 MHz). No feature loss or sign flip in any
configuration.

## Gate 3 — the window is measurable

`signal_chain.py` + `gate3_snr_map.py`. Chain δχ → δα → ΔOD → ΔT/T →
ΔN_ph → SNR with a fully sourced parameter table
(`signal_chain_parameters.csv`). At 1 ppm NV density, 0.5 mm sample, 1 µW
probe, η = 0.1: sector OD ≈ 7.9 (computed σ_ZPL ≈ 1.1 × 10⁻¹⁴ cm² at the
70 K homogeneous width), candidate ΔOD ≈ 0.11, shot-limited τ(SNR=5) ≈
2 × 10⁻⁵ s. All nine (contrast × density) rows of
`gate3_required_conditions.csv` are feasible within one hour — even a
100× smaller cross-section leaves seconds-scale integration times.

## Gate 4 — thresholds are intervals, quoted to matching precision

`gate4_threshold_uncertainty.py`. 500 Monte-Carlo samples over priors on
orbital-hopping scale (×/1.35), radiative scale (×/1.2), strain
(1.68 ± 0.34 GHz), Bx (1 %), Bz (5 %), gg (2×10⁻⁵–2×10⁻⁴), Ωc (×/1.1),
polarization angles (3°). All 500 samples yield all four thresholds.
Spearman ranking identifies the orbital-hopping (phonon) rate as the
dominant uncertainty for every threshold — consistent with the
dissipation-driven interpretation. Manuscript should quote e.g.
T_sign ≈ 102 K with 68 % band [95, 109] K (no more than 2 significant
figures).

## Gate 5 — ensemble washout quantified; selection conditions identified

`gate5_ensemble_average.py`. Fixed lasers applied to inhomogeneous
ensembles (4 orientations under the lab field, field spread, optical
offsets, strain spread, control-intensity spread; foreign orientations
enter exactly as Raman-detuned background absorbers). Findings:

- Unselected ensembles wash the window out by ×24 (low density) to ×105
  (high density, where the residue also distorts).
- Orientation + spectral post-selection alone is NOT sufficient: a 1 G
  field inhomogeneity smears the ~0.4 MHz dip to ~7 MHz (ground-splitting
  spread ≫ window width). **Field homogeneity ≲ 0.1 G over the probed
  volume is a required experimental selection condition.**
- With orientation + spectral selection + 0.1 G shimming, 19 % of the
  single-defect contrast survives (0.26 %, FWHM 1.1 MHz) — detectable in
  milliseconds by Gate 3.

## Decision points (SIMULATION_PLAN.md §8)

- **Decision A (after Gates 1–2): keep the EIT claim.** The verdict is
  statistically robust and the full Liouvillian confirms the reduced
  theory at the candidate point.
- **Decision B (after Gates 3, 5): keep the experimental-feasibility claim
  in the main text**, conditioned on the Gate 5 selection requirements
  (single orientation class, spectral selection, ≲0.1 G homogeneity),
  which should be stated explicitly in the experimental-prediction section.

## Reproduction

```bash
cd "No-go theorem"
python tests/test_core.py                      # 19 assertions
python scripts/reproduce_prl_figures.py --quick   # CI-sized full chain
python scripts/reproduce_prl_figures.py           # full statistics
```
