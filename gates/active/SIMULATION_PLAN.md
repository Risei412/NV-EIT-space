# Simulation Plan — NV EIT No-Go (PRL Verification Gates)

This document is the working plan for the numerical campaign that decides
whether the PRL claim can be submitted. It maps the five verification gates
of the submission roadmap onto the existing code base under
`No-go theorem/src/`, marks what is already done, and specifies the scripts,
inputs, outputs, and pass criteria for the remaining work.

**Central claim to be verified numerically:**

> In NV centers, a transverse magnetic field opens a symmetry-suppressed
> Raman pathway quadratically (K ∝ B⊥²), producing an EIT-type response at
> finite temperature, while thermal dissipation systematically suppresses it
> and defines a practical no-go boundary.

**Fixed candidate point (all gates use the same point unless stated):**

| Parameter | Value | Where fixed |
|---|---|---|
| Temperature | 70 K | roadmap §16 |
| B⊥ (Bx) | 0.23226 T | `run_prl_prediction.py` (`BX0`) |
| Bz | 0.005 T | `run_prl_prediction.py` (`BZ0`) |
| Ω_c | 0.1 GHz | `run_prl_prediction.py` (`OC`) |
| Control branch | ms = +1 | `run_prl_prediction.py` (`CTRL`) |

---

> **Status update (2026-07-16): Gates 1–5 implemented, executed, and ALL
> PASSED** — see `No-go theorem/results/gate_1_5_report.md` for numbers and
> caveats. Code: `src/gate{1..5}_*.py`, `src/signal_chain.py`; entry point
> `scripts/reproduce_prl_figures.py`; tests 19/19 green. Decision A: EIT
> claim kept (ΔAIC = +1902, full Liouvillian matches reduced kernel to
> 0.1 %). Decision B: feasibility claim kept, conditioned on the Gate 5
> selection requirements (single orientation, spectral selection, ≤0.1 G
> field homogeneity). Remaining before submission: compose final PRL
> Figures 1–4 from the gate outputs (§6) and the submission freeze (§8).

## 0. Current status (as of this plan)

Already completed and archived in `No-go theorem/results/`
(see `results/common_liouvillian_gate_report.md`):

- **Moment-order pipeline** — NV spin-Λ legs have |M0| = 0 (order −2 / −3),
  group-IV orbital-Λ legs |M0| = 1 (order −1); fitted slopes match predictions
  (`moment_order_common_pipeline.py`, `run_common_liouvillian_gate.py`).
- **B⊥² opening in the full Liouvillian** — reduced-kernel slope 2.16 vs
  full 9-level Lindblad slope 2.18 on the same Bx grid; C_full/C_red = 0.96–1.00
  at 100/150/300 K (`bperp_full_lowfield_slope.py`, `bperp_full_validation.py`).
  → **Gate 2 largely satisfied**; remaining item is the spectrum-level
  comparison at the candidate point itself (§2 below).
- **EIT/ATS classifier** — Anisimov–Dowling–Sanders A_EIT/A_ATS comparison
  with AIC/AICc and the |ΔAIC| ≥ 6 decision gate exists and passes sanity
  tests (`eit_ats_classifier.py`, `eit_ats_gate_run.py`). It has **not** yet
  been run with bootstrap/noise/window robustness at the fixed candidate point.
- **Candidate spectrum machinery** — `run_prl_prediction.py` produces the
  70 K candidate spectrum and Figures (candidate spectrum, field robustness,
  polarization tolerance, control sweep).

**Not yet implemented** (verified by grep over `src/`): photon-level
signal conversion (OD / ΔT/T / SNR), ensemble averaging over orientation /
strain / inhomogeneity, bootstrap resampling for the AIC verdict, and
parameter-uncertainty bands on the temperature thresholds.

---

## 1. Gate 1 — EIT/ATS statistical discrimination at the candidate point

**Goal:** show that the 70 K candidate window is identified as EIT (not ATS,
not a single Lorentzian, not Fano) robustly, not just at one fit.

**Existing code:** `eit_ats_classifier.py`, `eit_ats_gate_run.py`.

**New script:** `src/gate1_candidate_aic_bootstrap.py`

1. Generate the candidate spectrum from `run_prl_prediction.py` machinery
   (same Hamiltonian and Lindbladian as Gate 2 — one source of truth).
2. Fit four models: EIT, ATS, single Lorentzian + background, Fano.
3. Compute AIC, AICc, BIC, Akaike weights for each.
4. Robustness sweeps, each re-running the full model comparison:
   - additive noise realizations (≥ 200 bootstrap resamples, fixed seed),
   - frequency-window variations (±25 %, ±50 % of the fit window),
   - randomized fit initial values (≥ 50 draws),
   - control-power sweep Ω_c ∈ [0.02, 20] GHz to map the EIT→ATS crossover.

**Outputs:** `results/tables/gate1_aic_bootstrap.json`,
`results/tables/gate1_window_sweep.csv`,
`results/figures/gate1_model_selection_stability.png`.

**Pass criteria** (roadmap Gate 1):
- EIT beats ATS with |ΔAIC| ≥ 6 in ≥ 95 % of bootstrap resamples;
- verdict does not flip under any frequency-window choice;
- EIT→ATS transition vs Ω_c is continuous and physically ordered;
- verdict independent of initial values / fit constraints.

**Fallback if failed:** re-label the claim as
"transverse-field-induced narrow interference response" (roadmap §13,
Phase 1 fallback) and re-run this plan with the revised observable.

---

## 2. Gate 2 — Reduced kernel vs full Liouvillian at the candidate point

**Goal:** the opening exponent is already confirmed (status §0); what remains
is the **spectrum-level** comparison at the fixed candidate point.

**Existing code:** `liouvillian_core.py`, `nv_system.py`,
`weak_probe_response.py`, `nv_reduced_kernel.py`, `bperp_full_validation.py`.

**New script:** `src/gate2_candidate_full_vs_reduced.py`

Compare, at identical parameters, over the full two-photon detuning axis:
- full-model absorption spectrum vs reduced-model δχ_S;
- EIT contrast, linewidth, center frequency (fitted the same way for both);
- ground-state coherence ρ_12 and excited-state population from the full
  steady state;
- with/without: singlet shelving, ISC, hyperfine, both orbital branches
  (toggle flags, one row per configuration).

**Outputs:** `results/tables/gate2_candidate_comparison.csv`,
`results/figures/gate2_full_vs_reduced_spectrum.png`.

**Pass criteria:**
- opening exponent identical (already met, re-confirm at candidate point);
- contrast and linewidth agree within the same order of magnitude;
- no feature seen in the reduced model vanishes or flips sign in the full
  model, for every toggle configuration.

---

## 3. Gate 3 — From δχ to measurable signal (SNR) — *new*

**Goal:** convert the theoretical contrast chain
δχ → δα → ΔOD → ΔT/T → ΔN_ph → SNR and show a realistic detection window.

**New scripts:**
- `src/signal_chain.py` — pure functions for each conversion step, unit-tested.
- `src/gate3_snr_map.py` — applies the chain to the candidate spectrum.

**Physical inputs (all tabulated with literature source + uncertainty in
`results/tables/signal_chain_parameters.csv`):** NV density, sample
thickness, optical depth, ZPL branching ratio (Debye–Waller ≈ 3–4 %), phonon
sideband background, inhomogeneous broadening, probe/control intensities,
detection efficiency, shot noise, technical-noise floor, integration time.

**Outputs:**
- `results/tables/gate3_required_conditions.csv` — for each target contrast
  (1 %, 0.1 %, 0.01 %): required NV density × thickness × integration time;
- `results/figures/gate3_snr_map.png` — SNR over (density, integration time)
  at the candidate point, with the detectability frontier marked.

**Pass criteria:**
- required conditions listed for all three contrast levels;
- at least one region of (realistic density, ≤ mm thickness, ≤ hours
  integration) achieves SNR ≥ 5;
- the theoretical optimum used in Figures 3–4 lies inside the feasible region.

---

## 4. Gate 4 — Temperature-threshold uncertainty bands — *new analysis*

**Goal:** report T_1%, T_0.1%, T_0.01%, T_sign, T_no-go as intervals, not
single numbers.

**Existing code:** temperature sweeps exist
(`fixed_candidate_temperature_sweep.csv`, `threshold_crossings.csv`,
`phonon_rates.py`, `happacher_rate.py`).

**New script:** `src/gate4_threshold_uncertainty.py`

1. Define priors (interval or distribution, with literature citation) for:
   optical decoherence rate, phonon relaxation rate, strain, B⊥, Bz, dipole
   matrix elements, ISC branching ratio, ground-coherence time T2,
   control Rabi frequency, inhomogeneous width.
2. Latin-hypercube or Monte-Carlo sampling (≥ 500 samples, fixed seed);
   re-solve the temperature sweep per sample; extract each threshold via the
   existing `brentq` crossing logic.
3. Report median and 16–84 % / 2.5–97.5 % quantiles per threshold.

**Outputs:** `results/tables/gate4_threshold_bands.csv`,
`results/figures/gate4_threshold_bands.png` (thresholds vs T with error bands).

**Pass criteria:**
- every threshold quoted as an interval;
- dominant parameter per threshold identified (rank by Spearman correlation);
- significant digits in the paper trimmed to match the band widths.

---

## 5. Gate 5 — Ensemble averaging — *new*

**Goal:** verify the candidate window survives (or quantify how much it is
washed out by) averaging over realistic distributions.

**New script:** `src/gate5_ensemble_average.py`

Average the candidate spectrum over:
- the four NV orientations (fixed lab-frame B),
- local strain distribution (Gaussian, width from literature),
- static-field inhomogeneity,
- optical-transition frequency distribution (inhomogeneous broadening),
- control-intensity and MW-intensity spatial profiles,
- temperature gradient (small perturbation).

Run three scenarios: single defect (no averaging — reference), low-density
ensemble (orientation + mild strain), high-density ensemble (all
distributions at full width).

**Outputs:** `results/tables/gate5_ensemble_contrast.csv`,
`results/figures/gate5_ensemble_spectra.png` (three scenarios overlaid).

**Pass criteria:**
- with experimentally implementable post-selection (single orientation class,
  spectral hole selection) the feature survives;
- the unselected washout factor is quantified;
- the three density regimes are explicitly distinguished in the outputs.

---

## 6. Figure production (PRL Figures 1–4)

Single entry point: `scripts/reproduce_prl_figures.py` (new; thin wrapper
that calls the scripts below in order and asserts the output files exist).

| Figure | Content | Source script(s) | Status |
|---|---|---|---|
| Fig. 1 | Liouvillian block structure, sector S, Schur complement, forbidden path at B⊥=0 vs escape channel at B⊥≠0 | schematic + `moment_order_common_pipeline.py` inset | partial |
| Fig. 2 | \|K\| vs B⊥ (reduced) and contrast vs B⊥ (full), log–log fits, low/high-field regimes | `bperp_kernel_map_v2.py`, `bperp_full_lowfield_slope.py` | data done, figure to compose |
| Fig. 3 | (T, B⊥) phase map, contrast color, 1 %/0.1 %/0.01 % thresholds, sign-reversal line, EIT/ATS boundary, visibility boundary (from Gate 3) | `run_prl_prediction.py` + Gate 3/4 outputs | needs Gates 3–4 |
| Fig. 4 | Candidate spectrum on/off control, EIT vs ATS fits, AICc/BIC panel, transmission + SNR panel | Gate 1 + Gate 3 outputs | needs Gates 1, 3 |

---

## 7. Automated tests (extend `tests/test_core.py`)

Add, in this order of priority:

1. trace / Hermiticity preservation of the full Liouvillian steady state;
2. positivity audit (min eigenvalue of ρ_ss ≥ −1e−10);
3. zero-control limit → bare absorption; zero-field limit → K ≡ 0;
   cut-sector limit → δχ_S ≡ 0;
4. finite-difference check of the weak-probe linear response;
5. frequency-grid convergence (halving the grid changes contrast < 1 %);
6. `signal_chain.py` unit conversions (round-trip and known reference values);
7. fixed seeds everywhere; classifier verdict reproducibility across runs.

CI note: every `gate*_*.py` script must be runnable headless
(`matplotlib.use('Agg')`) and finish < 30 min on a laptop-class CPU; longer
sweeps (Gate 4 Monte Carlo) accept an `--n-samples` flag with a small
CI default.

---

## 8. Execution order and decision points

```
Week 1   Gate 1 (AIC + bootstrap at candidate point)      ── decision A
Week 1-2 Gate 2 spectrum-level comparison                 ── decision A
Week 2-3 Gate 3 signal chain + SNR map                    ── decision B
Week 3-4 Gate 4 threshold uncertainty bands
Week 4   Gate 5 ensemble averaging                        ── decision B
Week 5   Figures 1-4 frozen, reproduce_prl_figures.py, tests green
```

- **Decision A (after Gates 1–2):** if the EIT verdict is unstable or ATS
  wins, or the window vanishes in the full model → switch the claim to the
  interference-response wording and retarget PRA/PRB; Gates 3–5 still run,
  with the revised observable.
- **Decision B (after Gates 3, 5):** if the converted signal is unrealistically
  small or fully washed out without post-selection → keep the theory claim,
  demote the experimental-feasibility claim to the Supplement.

**Submission freeze (when all gates pass):** tag the commit, create a GitHub
release, mint a Zenodo DOI, record the commit hash in the manuscript, and
verify `scripts/reproduce_prl_figures.py` regenerates every paper figure from
a clean clone.
