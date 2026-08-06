# P0-1 — NV observable freeze and Gate E rerun

Date: 2026-08-07

## Decision

The PRL must not assign a single dissipation order to both the raw signed
absorption response and the normalized experimental contrast.

Freeze the two observables as

\[
\Delta A=A_{\rm cut}-A_{\rm full},
\qquad
C=\frac{\Delta A}{A_{\rm cut}}.
\]

The observable-inheritance result gives

\[
\Delta A\sim\Gamma^{-4},
\qquad
A_{\rm cut}\sim\Gamma^{-1},
\qquad
C\sim\Gamma^{-3}.
\]

Therefore:

- theorem-facing observable: raw signed response \(\Delta A\), \(\nu=4\);
- direct NV experimental observable: normalized contrast \(C\), \(\nu=3\);
- the manuscript must show the normalization bridge explicitly.

## Recomputed 70 K ensemble anchor

For the post-selected and field-shimmed ensemble:

| quantity | value |
|---|---:|
| \(C_{\rm ref}\) | \(2.6330768\times10^{-3}\) |
| \(A_{\rm cut}\) | \(9.0307923\times10^{-3}\) |
| \(\Delta A_{\rm ref}\) | \(2.3778770\times10^{-5}\) |
| \(C_{\min}\) | \(1.1599052\times10^{-6}\) |
| \(\Delta A_{\min}=C_{\min}A_{\rm cut}\) | \(1.0474863\times10^{-8}\) |
| signal-to-floor margin | \(2.2700794\times10^3\) |

The Gate-5 ensemble calculation now exports `Cmax`, `Aoff_at_Cmax`, and
`dA_at_Cmax` at the same spectral feature.

## Window verdict

\[
D_\Gamma=\frac{\log_{10}(|O_{\rm ref}|/O_{\min})}{\nu}.
\]

| observable | order | usable window | verdict |
|---|---:|---:|---:|
| raw \(\Delta A\) | 4 | 0.8390 decade | FAIL |
| normalized \(C\) | 3 | 1.1187 decades | PASS |

The previous 0.839-decade result for normalized contrast was caused by using
\(\nu=4\) on a \(\nu=3\) quantity.

## Identifiability verdict

A deterministic 5000-run Monte-Carlo stress test was rerun separately for the
two observables.

| observable | constrained correction | free correction/background |
|---|---:|---:|
| \(\Delta A\), true \(n=4\) | 0.9996 | 0.5998 |
| \(C\), true \(n=3\) | 0.9998 | 0.6538 |

Thus the normalized contrast now supplies enough raw dynamic range for a
one-decade scan, but a free fit still cannot separate adjacent integer classes
reliably. The present Gate E remains a **conditional pass**.

## PRL consequence

Allowed:

> The raw sector response carries the exact fourth-order class, while the
> normalized NV contrast inherits a third-order law through the first-order
> cut absorption. With the leading finite-response correction calibrated
> independently, the NV experiment can test the corresponding class.

Not yet allowed:

> A free power-law fit to the NV ensemble data independently discovers the
> integer class.

## Files changed

- `No-go theorem/src/gate5_ensemble_average.py`
- `No-go theorem/results/tables/gate5_ensemble_contrast.csv`
- `No-go theorem/results/tables/gate5_summary.json`
- `New no-go theory/GateE_NV_experimental_anchor/src/run_gate_e.py`
- `New no-go theory/GateE_NV_experimental_anchor/tests/test_gate_e.py`
- `New no-go theory/GateE_NV_experimental_anchor/results/tables/*`
- `New no-go theory/GateE_NV_experimental_anchor/README.md`

## Next gate

P0-2 remains the next logical blocker for the manuscript: separate
sector-induced signed transparency from spectroscopic EIT/Fano/ATS model
classification. P0-3 then replaces the synthetic finite-\(\Gamma\) correction
with end-to-end full-Liouvillian pseudo-data.
