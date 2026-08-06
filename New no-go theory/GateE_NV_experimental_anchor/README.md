# Gate E — NV experimental anchor and exponent identifiability

This gate is the next step after the general-theory Gates A–D. It asks a more
specific question required by the revised PRL plan:

> Can the planned NV ensemble experiment itself identify the integer response
> class, or should it be used as a model-constrained experimental anchor while
> the clean exponent measurement is supplied by the engineered-dissipation
> witness?

## Frozen inputs

No new contrast value is invented here. The calculation reuses:

- `No-go theorem/results/tables/gate5_ensemble_contrast.csv`
- `GateD_robustness_discriminability/results/tables/gates_summary_gateD.json`

At the 70 K candidate point, the post-selected and field-shimmed ensemble has

- signed-feature magnitude `C_ref = 2.6331e-3`,
- minimum detectable contrast `C_min = 1.1599e-6`,
- contrast-to-floor margin `C_ref/C_min = 2.2701e3`.

The observable inherited from the two-sided sector response is taken to have
order `nu_obs = 4`, consistent with Gate A and Gate D.

## Calculation

For a response `C(Gamma) ~ Gamma^-nu`, the usable dissipation window before the
signal reaches the detection floor is

```text
D_Gamma = log10(C_ref/C_min) / nu .
```

For the gating ensemble scenario:

- intensity/contrast readout, `nu=4`: `D_Gamma = 0.839` decade,
- amplitude-linear readout, `nu=2`: `D_Gamma = 1.678` decades.

A one-decade intensity window therefore needs a further factor

```text
4.405
```

in contrast-to-floor margin. Because `Gamma(T) ~ T^5`, one decade in Gamma
requires only a temperature ratio `T_max/T_min = 10^(1/5) = 1.585`; the main
bottleneck is not temperature reach but observable identifiability.

## Finite-Gamma identifiability test

The synthetic data use

```text
C(Gamma) = A Gamma^-n (1 + b/Gamma) + background,
```

with `n=4`, `b=0.30`, 16 log-spaced points across the available 0.839-decade
window, 3% multiplicative drift, and additive noise fixed by the Gate-D SNR=5
floor. Candidate classes are `n=3,4,5`.

Two protocols are compared over 5000 deterministic Monte-Carlo repetitions:

| fit protocol | probability of selecting n=4 |
|---|---:|
| first finite-Gamma correction fixed independently | 0.9996 |
| correction and background freely refitted for every class | 0.5998 |

The second protocol is the conservative one. It shows that a visually clean
log-log trace is not enough: over the current window, adjacent integer classes
can absorb each other through finite-Gamma corrections.

## Gate verdict

**CONDITIONAL PASS**

- `G-E1` signal detectability: **PASS**.
- `G-E2` one-decade raw intensity window: **FAIL** (`0.839 < 1`).
- `G-E3` model-constrained class identification: **PASS**.
- `G-E4` unconstrained raw asymptotic identification: **FAIL**.
- `G-E5` amplitude-linear readout window: **PASS**.

## Consequence for the PRL story

The present NV ensemble experiment should not be advertised as a standalone
measurement of the integer exponent from an unconstrained power-law fit.
Instead it should serve as the experimental anchor for the mechanism, signed
response, and intervention/crossover behavior. The clean integer-class
measurement remains the superconducting or other engineered-dissipation
witness unless one of the following is achieved:

1. `Gamma_ph(T)` and the leading finite-Gamma correction are independently
   calibrated and frozen before fitting the NV data;
2. an amplitude-linear readout is implemented;
3. the contrast-to-floor margin is improved by at least `4.405x`, together with
   sufficient control of finite-Gamma corrections.

The safest PRL claim is therefore:

> The general integer-class law is established by theorem and blind physical
> witnesses, while the NV experiment reveals its selection-rule-controlled
> loss of coherent response in a solid-state optical system.

## Contents

- `src/run_gate_e.py` — self-contained window and Monte-Carlo audit.
- `tests/test_gate_e.py` — deterministic regression checks.
- `results/tables/gates_summary_gateE.json` — complete frozen result.
- `results/tables/gate_e_windows.csv` — readout windows and required margin.
- `results/tables/gate_e_identifiability.csv` — constrained/free model selection.

Run:

```bash
python src/run_gate_e.py
python tests/test_gate_e.py
```
