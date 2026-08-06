# Gate E — NV experimental anchor and exponent identifiability

This gate asks whether the planned NV ensemble experiment can identify the
integer response class directly, or whether it should remain a
model-constrained experimental anchor.

## P0-1 observable freeze

The previous Gate E used the ensemble **normalized contrast** as its numerical
input but assigned it the asymptotic order of the **raw signed response**. P0-1
separates the two quantities throughout the code, tables, floors, and fits:

\[
\Delta A=A_{\rm cut}-A_{\rm full},
\qquad
C=\frac{\Delta A}{A_{\rm cut}}.
\]

Observable inheritance gives

\[
\Delta A\sim\Gamma^{-4},
\qquad
A_{\rm cut}\sim\Gamma^{-1},
\qquad
C\sim\Gamma^{-3}.
\]

Thus the theorem-facing observable is the raw signed absorption difference
\(\Delta A\) with \(\nu=4\), while the directly normalized experimental
observable is \(C\) with \(\nu=3\). The order of \(C\) is derived in
`run_gate_e.py` as numerator order minus denominator order rather than entered
as an independent hard-coded class.

## Frozen inputs

The calculation reuses:

- `No-go theorem/results/tables/gate5_ensemble_contrast.csv`
- `GateD_robustness_discriminability/results/tables/gates_summary_gateD.json`

Gate 5 now exports, at the same spectral feature, `Cmax`, `Aoff_at_Cmax`, and
`dA_at_Cmax`. For the 70 K post-selected and field-shimmed ensemble:

| quantity | value |
|---|---:|
| normalized contrast \(C_{\rm ref}\) | `2.6330768e-3` |
| cut absorption \(A_{\rm cut}\) | `9.0307923e-3` |
| raw signed difference \(\Delta A_{\rm ref}\) | `2.3778770e-5` |
| contrast floor \(C_{\min}\) | `1.1599052e-6` |
| raw floor \(\Delta A_{\min}=C_{\min}A_{\rm cut}\) | `1.0474863e-8` |
| signal-to-floor margin | `2.2700794e3` |

The raw and normalized quantities have the same signal-to-floor margin because
the floor is expressed on the same detection chain. They have different usable
\(\Gamma\) windows because their asymptotic orders differ.

## Dissipation windows

For an observable \(O\sim\Gamma^{-\nu}\),

\[
D_\Gamma=\frac{\log_{10}(|O_{\rm ref}|/O_{\min})}{\nu}.
\]

For the gating ensemble scenario:

| observable | order | usable window | one-decade gate |
|---|---:|---:|---:|
| raw \(\Delta A\) | 4 | `0.8390` decade | **FAIL** |
| normalized \(C\) | 3 | `1.1187` decades | **PASS** |

The old conclusion that the normalized contrast had only `0.839` decade was an
observable-order mismatch. Correcting it removes the raw window-size blocker
for \(C\), but does not by itself establish free exponent identifiability.

## Finite-\(\Gamma\) identifiability

The synthetic stress test uses

\[
O(\Gamma)=A\Gamma^{-n}(1+b/\Gamma)+B,
\]

with `b=0.30`, 16 logarithmically spaced points, 3% multiplicative drift, and
additive noise fixed by the Gate-D SNR=5 floor. The candidate classes are the
true integer and its two neighbours.

| observable | fit protocol | correct-class probability |
|---|---|---:|
| \(\Delta A\), true \(n=4\) | correction fixed independently | `0.9996` |
| \(\Delta A\), true \(n=4\) | correction and background free | `0.5998` |
| \(C\), true \(n=3\) | correction fixed independently | `0.9998` |
| \(C\), true \(n=3\) | correction and background free | `0.6538` |

The conservative free fit still fails. Even with the corrected `1.119`-decade
contrast window, adjacent classes absorb one another through finite-\(\Gamma\)
corrections and background freedom.

## Gate verdict

**CONDITIONAL PASS**

- both observables are above their matched detection floors: **PASS**;
- normalized contrast has a one-decade window: **PASS**;
- raw signed difference has a one-decade window: **FAIL**;
- model-constrained identification succeeds for both: **PASS**;
- unconstrained identification succeeds for both: **FAIL**;
- observable orders are derived consistently: **PASS**.

## Consequence for the PRL story

The paper should use the observables in distinct roles:

- **theorem statement and exact class certificate:** raw signed response
  \(\Delta A\), \(\nu=4\);
- **direct NV experimental plot:** normalized contrast \(C\), \(\nu=3\);
- **bridge between them:** the measured or modeled normalization
  \(A_{\rm cut}(\Gamma)\sim\Gamma^{-1}\).

The NV experiment should still not be advertised as a standalone discovery of
an integer exponent from an unconstrained power-law fit. It becomes a direct
one-decade **signal-window** candidate, but remains a model-constrained class
measurement until P0-3 generates the finite-\(\Gamma\) correction and nuisance
structure from the full NV Liouvillian rather than from the synthetic
`b=0.30` stress model.

The safe PRL sentence is:

> The raw sector response obeys the theorem-level fourth-order suppression,
> while the experimentally normalized NV contrast inherits a third-order law;
> NV data test the mechanism and the normalization-linked crossover, with the
> finite-response correction fixed independently.

## Contents

- `src/run_gate_e.py` — observable definitions, matched floors, windows, and
  Monte-Carlo identifiability.
- `tests/test_gate_e.py` — regression checks that pin \(\nu_{\Delta A}=4\),
  \(\nu_C=3\), and the two different window verdicts.
- `results/tables/gates_summary_gateE.json` — complete frozen result.
- `results/tables/gate_e_windows.csv` — both observables for all ensemble
  scenarios.
- `results/tables/gate_e_identifiability.csv` — constrained/free selection for
  both integer classes.

Run:

```bash
python src/run_gate_e.py
python tests/test_gate_e.py
```
