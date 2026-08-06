# P0-2 — two-layer NV phase classification

Date: 2026-08-07  
Status: **PASS, full 240-point grid completed**

## Why P0-2 was required

The old P1 classifier combined three logically different statements:

1. the sign of the pathway-cut response;
2. the EIT-versus-ATS information-criterion difference;
3. the best fit among EIT, ATS, Fano and Lorentzian.

Its final label used item 2 even when item 3 selected Fano. Thus `robust EIT`
could mean only “EIT beats ATS”, not “EIT beats every relevant alternative”.
P0-2 separates the mechanism and spectroscopy layers.

## Frozen definitions

The mechanistic observable is

\[
\Delta A=A_{\rm cut}-A_{\rm full}.
\]

- \(\Delta A>0\): `sector_transparency`;
- \(\Delta A<0\): `sector_absorption`;
- failed containment or detection gate: `sector_unresolved`.

The measured line \(A_{\rm full}(\delta_2)\) is independently fitted to EIT,
ATS, Fano and a single Lorentzian. A robust spectral label requires the best
model to beat the runner-up by at least six AIC/AICc units. EIT and ATS also
require the corresponding sign of the ATS-minus-EIT difference.

> Rejecting ATS is not sufficient to establish EIT.

## Full-grid result

The calculation uses 20 temperatures from 20 to 115 K and 12 transverse-field
values from 0 to 0.5 T, for 240 full-Liouvillian adaptive spectra.

### Mechanistic layer

| class | points |
|---|---:|
| positive sector response | 169 |
| control-induced absorption | 68 |
| unresolved | 3 |

### Spectroscopic layer

| best robust class | points |
|---|---:|
| Fano | 220 |
| EIT | 9 |
| ATS | 8 |
| ambiguous | 3 |

### Joint layer

| joint class | points |
|---|---:|
| Fano-shaped sector transparency | 154 |
| spectroscopic EIT | 9 |
| ATS with positive sector response | 3 |
| positive sector response, spectrally ambiguous | 3 |
| control-induced absorption | 68 |
| unresolved | 3 |

The old classifier labelled 148 points as `transparency`. Of these, 139 had a
non-EIT best model. The mixed-label rate was therefore 93.9%.

## Surviving strict EIT region

A genuine four-model EIT region survives, but it is small and confined to low
temperature and strong transverse field:

- 25 K at \(B_\perp=0.20\) through 0.50 T;
- 30 K at \(B_\perp=0.20\) T.

There are nine strict EIT grid points in total. No strict EIT point exists at
\(B_\perp=0\), and no strict EIT survives at the planned 70 K operating point.

Thus the old phrase “the genuine EIT region is a closed island” was too broad.
The corrected statement is:

> a small low-temperature EIT subset lies inside a much larger bounded region
> of positive coherence-sector response.

## The 70 K candidate

At

\[
T=70~\mathrm{K},\qquad B_\perp=0.232261~\mathrm{T},
\]

the full calculation gives

\[
C_{\max}=1.3871069\times10^{-2},\qquad
\Delta A=3.7788842\times10^{-4}>0,
\]

with a 0.375 MHz feature. The adaptive 1.5 MHz half-window selects Fano over
Lorentzian by 89.45 AICc units and rejects ATS relative to EIT by 327.12 units.
The joint label is therefore

> **Fano-shaped sector transparency**, not spectroscopic EIT.

A broad-window audit remains important: at half-windows 10, 15, 20 and 40 MHz
the best model is Fano, whereas at 30 MHz it is Lorentzian. The 20 MHz Fano-to-
Lorentzian margin is only 0.327 AIC. Hence the manuscript may say
“Fano-shaped narrow feature”, but should not claim that a unique microscopic
Fano mechanism has been identified independently of the background window.

## Zero-field audit

At \(B_\perp=0\), positive sector response appears at 50, 55, 60, 65 and 70 K.
It is not strict EIT at any of these points. Consequently, transverse field is
not an exact binary switch from a mathematically zero response in the complete
model.

The defensible statement is:

> transverse field enhances and stabilizes the experimentally useful positive
> sector response, while residual zero-field pathways remain.

An exact-zero claim would require an additional operational cut or a
symmetry-restricted observable that removes those residual pathways.

## PRL consequence

### Allowed

- NV exhibits a selection-rule-controlled signed coherence-sector response.
- The useful positive response occupies a bounded temperature-field region.
- The response weakens and eventually reverses sign under stronger phonon mixing.
- The 70 K candidate is strongly distinct from ATS.
- NV can serve as the experimental mechanism anchor without being called EIT.
- A small strict-EIT subset exists near 25–30 K at sufficiently large transverse field.

### Not allowed

- the entire positive-response region is a genuine EIT island;
- the 70 K candidate is spectroscopically identified as EIT;
- a large positive ATS-minus-EIT score alone proves EIT;
- \(B_\perp\) creates the response from an exact zero in the full model;
- the broad-window data uniquely identify a microscopic Fano mechanism.

## Reproducibility

The quick grid reproduced the archived P1 values, including the 70 K candidate,
to numerical precision. The full run used four workers with one BLAS thread per
worker and completed in 274.2 s. The numeric full-grid CSV generated in the run
had SHA-256

`29acae87046b02c9efbdc2dad3389aaa3a1cad3272fd1ef9340007110f974f3a`.

Frozen repository outputs are:

- `p0_2_two_layer_summary.json`;
- `p0_2_joint_class_grid.csv`;
- `p0_2_boundaries.csv`;
- `p0_2_candidate_preflight.json`.

## Gate verdict

All P0-2 gates pass:

- mechanism and spectroscopy are stored separately;
- all four spectral models enter the runner-up test;
- no EIT label has a non-EIT best model;
- \(B_\perp=0\) is explicitly audited;
- the 70 K candidate is resolved and relabelled;
- the old mixed-label population is quantified.

P0-2 is complete. The next computational gate is P0-3: an end-to-end
pseudo-experiment comparing integer-class recovery against signed-crossover
recovery under realistic nuisance parameters.
