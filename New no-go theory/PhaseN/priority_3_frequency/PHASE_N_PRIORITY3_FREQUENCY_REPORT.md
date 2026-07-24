# Phase N Priority 3 Frequency Certification

Date: 2026-07-24  
Model: fixed five-level GKSL Phase N witness  
Sector: \(S_{34}\)  
Production runtime: 435.11 s  
Status: **F0–F9 PASS**

## Executive conclusion

The frequency-dependent operational response

\[
R_{S_{34}}(\Gamma,\kappa,z)
=
\frac{N(\Gamma,\kappa,z)}{Q(\Gamma,\kappa,z)},
\qquad
\kappa=\alpha\Gamma^q,
\]

was constructed as an exact three-variable rational function for

\[
A(\Gamma,\kappa,z)
=\Gamma D+B_0-izI+\kappa D_{S_{34}}.
\]

At \(z=0\), the numerator and denominator exactly reproduce the previously
certified Phase N polynomials. For a generic fixed frequency, the original
path-order fan remains

\[
\nu_{\mathrm{gen}}(q)=
\begin{cases}
4-q,&0\le q\le1,\\
2+q,&1\le q\le2,\\
4,&q\ge2.
\end{cases}
\]

An isolated exact frequency,

\[
\boxed{z_\star=\frac{543}{280}\simeq1.939285714,}
\]

cancels the fourth-order ideal-cut coefficient. The pointwise ideal-cut order
is promoted from \(4\) to \(5\), and the path-order fan becomes

\[
\nu_\star(q)=
\begin{cases}
4-q,&0\le q\le1,\\
2+q,&1\le q\le3,\\
5,&q\ge3.
\end{cases}
\]

Thus the Newton breakpoints change exactly from

\[
\{1,2\}\longrightarrow\{1,3\}.
\]

This promotion is pointwise in frequency. It does not promote the asymptotic
order of fixed nonzero-width \(L^1\), \(L^2\), or supremum frequency-window
norms, because generic neighboring frequencies have the lower order.

## Exact frequency promotion

The first nonzero ideal-cut moment at generic frequency is

\[
\mathfrak m_4(z)
=-\frac{173760}{140777}
+\frac{12800}{20111}z.
\]

Its unique real root is \(z_\star=543/280\). At that root,
\(\mathfrak m_5(z_\star)\ne0\), proving exact promotion by one order rather
than an unresolved higher cancellation.

The special frequency lies inside the certified wide window
\([-3,3]\). Reversing the Fourier-transform sign convention reverses the sign
of \(z_\star\) without changing any order statement.

## Frequency unfolding

Writing

\[
\varepsilon=z-z_\star,
\]

the certified crossover law is

\[
\Gamma_{\times,z}\propto
\begin{cases}
|\varepsilon|^{-1/(q-2)},&2<q<3,\\
|\varepsilon|^{-1},&q\ge3.
\end{cases}
\]

Production fits were:

| \(q\) | fitted \(d\log\Gamma_{\times,z}/d\log|\varepsilon|\) | predicted | absolute error |
|---:|---:|---:|---:|
| 2.25 | -4.0000000000 | -4 | \(3.98\times10^{-12}\) |
| 2.50 | -1.9998793861 | -2 | \(1.21\times10^{-4}\) |
| 2.75 | -1.3140844265 | -1.3333333333 | \(1.92\times10^{-2}\) |
| 3.00 | -0.9999866550 | -1 | \(1.33\times10^{-5}\) |
| 3.50 | -0.9985196809 | -1 | \(1.48\times10^{-3}\) |

The maximum collapse spread ranged from \(3.85\times10^{-12}\) to
\(2.07\times10^{-2}\), below the production threshold in every tested cell.
The 80-digit and 120-digit crossover locations agreed within the stored
tolerance.

## Finite-band interference unfolding

The interference-coupling perturbation from Priority 1 was repeated after
integrating over the nonzero-width frequency window \(K=[-1,1]\). For
\(q=3/2\), the \(L^2\)-window crossover obeyed

\[
\Gamma_{\times,\delta}\propto|\delta|^{-2}.
\]

The fitted slope over the asymptotic perturbations
\(\delta=10^{-3},10^{-4},10^{-5}\) was

\[
-2.0151487,
\]

and the maximum rescaled-curve spread was \(8.15\times10^{-3}\). Therefore
the Priority 1 interference-unfolding law is not a single-frequency artifact.

## Finite frequency windows

The norms

\[
M_\infty=\sup_{z\in K}|R(z)|,\qquad
M_1=\int_K|R(z)|\,dz,\qquad
M_2=\left(\int_K|R(z)|^2\,dz\right)^{1/2}
\]

were evaluated for \(K=[-1,1]\) and \(K=[-3,3]\), with the latter containing
\(z_\star\). All three norms converged to the generic fan for
\(q=0.5,1.5,2.5,3.5\).

The largest absolute fitted-order error was

\[
3.23\times10^{-7}.
\]

Grid refinement from 121 to 241 to 481 frequency points changed the fitted
wide-window \(L^2\) order by at most

\[
2.84\times10^{-14}.
\]

Consequently, the frequency-selective order-five result must be described as
an isolated pointwise promotion with a parametrically long narrow-band
crossover, not as an order-five finite-band phase.

## Physical and numerical certification

For real \(z\) and \(\kappa\ge0\), the Hermitian part of the reduced pencil is
positive:

\[
\operatorname{Re}A
=\Gamma D+\kappa D_{S_{34}}
\succeq\frac{\Gamma}{2}I.
\]

Hence no real-frequency pole can occur for \(\Gamma>0\). The smallest sampled
singular value over the production grid was \(0.617671\).

The maximum frequency-domain reduced-response versus full-Liouvillian error
was

\[
6.99\times10^{-18}.
\]

Trace-preservation and Hamiltonian-Hermiticity residuals were zero at machine
precision.

## Gate summary

| Gate | Result |
|---|---|
| F0 \(z=0\) exact Phase N regression | PASS |
| F1 exact \(z_\star\) and order \(4\to5\) | PASS |
| F2 generic and special fans | PASS |
| F3 frequency-unfolding slope | PASS |
| F4 crossover collapse | PASS |
| F5 nonzero-band Priority 1 collapse | PASS |
| F6 \(L^1,L^2,L^\infty\) window orders | PASS |
| F7 no real-axis pole | PASS |
| F8 full-Liouvillian frequency consistency | PASS |
| F9 precision and grid refinement | PASS |

## Consequence for the central claim

The result supports a stronger operational statement:

> The sector suppression order is not a scalar invariant of the unperturbed
> GKSL generator or of the sector label. It is a valuation selected by the
> experimental path through intervention strength and probe frequency.
> Interference manifolds create exact promoted orders, while their
> neighborhoods exhibit universal crossover laws.

The model-specific statement must remain narrower:

> In the certified five-level witness, the tuned interference manifold carries
> the exact V-shaped intervention fan. A second, frequency-selective
> codimension-one cancellation promotes the pointwise order from four to five,
> but finite-band norms retain the generic fan.

## Next numerical priority

The common resolvent shift \(-izI\) is mathematically controlled, but a
physical EIT implementation generally places laser detunings with different
weights on different coherences. The next calculation should therefore
replace \(zI\) by a physically derived diagonal detuning signature \(zZ\) and
classify:

1. which physical detuning directions preserve the generic V-shaped fan;
2. whether a real, accessible frequency-promotion root survives;
3. how inhomogeneous broadening and finite spectral resolution modify the
   crossover;
4. the dimensional parameter range and signal magnitude in one concrete
   double-\(\Lambda\), cavity-EIT, ion, or solid-state realization.

Noise and exponent-estimation calculations should follow this physical
embedding. A second abstract GKSL architecture is lower priority until the
detuning signature and observability of the first architecture are fixed.

