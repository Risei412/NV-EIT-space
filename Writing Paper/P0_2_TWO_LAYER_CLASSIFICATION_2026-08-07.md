# P0-2 — sector response and spectroscopic classification must be separated

Date: 2026-08-07

## Question

The old P1 phase diagram assigned one joint label from three ingredients:

1. the sign of the pathway-cut response,
2. the EIT-versus-ATS information-criterion difference,
3. the best fit among EIT, ATS, Fano and Lorentzian.

The implementation used item 2 for the final label and merely recorded item 3.
Consequently, a spectrum could be called `transparency` while its best model was
Fano. P0-2 asks whether the planned NV anchor is spectroscopically identified as
EIT once all four models are treated symmetrically.

## Frozen two-layer definition

### Mechanistic layer

The pathway-sector object is

\[
\Delta A=A_{\rm cut}-A_{\rm full}.
\]

- \(\Delta A>0\): `sector_transparency`, the selected coherence pathway lowers absorption;
- \(\Delta A<0\): `sector_absorption`, the pathway raises absorption;
- below the detection/containment gate: `sector_unresolved`.

This layer is a counterfactual mechanism statement. It is not, by itself, a
spectroscopic EIT diagnosis.

### Spectroscopic layer

The measured control-on line \(A_{\rm full}(\delta_2)\) is fitted to

- EIT,
- ATS,
- Fano,
- single Lorentzian.

A robust label requires the best model to beat the runner-up by at least six
AIC/AICc units. EIT and ATS must additionally have the appropriate sign of the
EIT-versus-ATS difference. In particular,

> rejecting ATS does not establish EIT.

## Decisive 70 K candidate result

At the planned candidate

\[
T=70~\mathrm{K},\qquad B_\perp=0.232261~\mathrm{T},
\]

the full-Liouvillian phase table gives

\[
C_{\max}=1.3871\times10^{-2},\qquad
\Delta A=3.7789\times10^{-4}>0.
\]

The point is therefore a clear positive sector response.

The archived four-model comparison gives

| model | AIC |
|---|---:|
| Fano | -23021.9700 |
| Lorentzian | -23021.6432 |
| EIT | -23019.3940 |
| ATS | -21117.2349 |

Hence:

- ATS is very strongly rejected relative to EIT;
- Fano is the numerical best model;
- Fano beats Lorentzian by only 0.327 AIC;
- EIT trails Fano by 2.576 AIC;
- no model beats its runner-up by the required six units.

The strict candidate label is therefore

> **sector transparency + spectrally ambiguous**, not spectroscopic EIT.

The frequency-window sweep reinforces the verdict. The best model is Fano for
half-windows 10, 15, 20 and 40 MHz and Lorentzian for 30 MHz. EIT is never the
best model in that sweep.

## Consequence for the old phase diagram

Many rows previously labelled `transparency` have `best_model=Fano`. Those rows
must be renamed `Fano-shaped sector transparency` or, when the Fano advantage
is below the robust threshold, `sector transparency / spectrally ambiguous`.

The old headline

> the genuine EIT region is a closed island

is not supported by the current four-model analysis.

The supported replacement is

> the coherence-mediated signed response occupies a bounded region of the
> temperature-field plane and undergoes a high-temperature sign reversal.

Whether a smaller subset qualifies as robust spectroscopic EIT is the output
of the full P0-2 grid rerun, not an assumption.

## Transverse-field statement

The archived grid contains small positive sector responses at \(B_\perp=0\),
including 55 and 60 K. Therefore the transverse field is not an absolute binary
switch that creates the pathway from an exact zero in this full model.

The defensible statement is narrower:

> transverse field enhances and stabilizes the experimentally useful positive
> sector response, while residual zero-field pathways remain in the complete
> model.

An exact-zero switch claim would require an additional operational cut or a
symmetry-restricted observable that removes those residual paths.

## Implementation

PR #25 adds `p0_2_two_layer_classification.py`, which recomputes the full
240-point grid and exports separate sector, spectral and joint maps. Its gates
require:

- both classification layers to be present;
- the four-model runner-up margin to be enforced;
- no `spectroscopic_EIT` row with a non-EIT best model;
- an explicit \(B_\perp=0\) audit;
- a resolved 70 K candidate label.

The branch is intentionally not merged until the full-grid outputs are
available and inspected.

## PRL decision

### Allowed now

- NV demonstrates a selection-rule-controlled signed coherence-sector response.
- The response is strongly distinct from ATS at the candidate.
- The response loses strength and reverses sign with increasing phonon mixing.
- NV can anchor the pathway-moment mechanism without claiming a textbook EIT lineshape.

### Not allowed now

- the 70 K candidate is statistically identified as EIT against all alternatives;
- the full positive-response region is a genuine EIT island;
- \(B_\perp\) opens an exactly forbidden response from zero.

## Next computation

After the P0-2 grid is frozen, P0-3 should generate end-to-end pseudo-data and
ask which of the two experimentally relevant claims is identifiable:

1. the normalized dissipation class, or
2. the intervention-controlled signed-response crossover.
