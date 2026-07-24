# Phase Z (Phase N Priority 4): Physical Detuning Signature zZ

Date: 2026-07-24
Model: fixed five-level GKSL Phase N witness
Sector: \(S_{34}\)
Production runtime: 266.9 s
Status: **Z0–Z6 PASS** (`results/gates_summary_phaseZ.json`)

## Executive conclusion

The Priority 3 frequency certification used the common resolvent shift
\(-izI\). A physical EIT implementation instead places laser detunings with
different weights on different coherences. Phase Z replaces the shift by a
diagonal physical detuning signature,

\[
A(\Gamma,\kappa,z;Z)=\Gamma D+B+\kappa D_{S_{34}}-izZ,
\qquad Z=\operatorname{diag}(\zeta_1,\dots,\zeta_4),
\]

and certifies, exactly where possible:

1. **The generic V-shaped intervention fan is direction-independent.** For
   every physical detuning signature tested (and for a generic rational
   control signature), the exact path-order fan is unchanged:
   \(\nu(q)=4-q\ (0\le q\le1),\ 2+q\ (1\le q\le2),\ 4\ (q\ge2)\).
   The path-order theorem's central claim therefore does not depend on the
   idealized common-shift frequency variable.
2. **The frequency-promotion root is direction-dependent, and physically
   steerable.** The first nonzero ideal-cut moment is
   \(\mathfrak m_4(z;Z)=-\tfrac{173760}{140777}+c_1(Z)\,z\) with a
   direction-dependent slope \(c_1(Z)\). Sweeping the branch-3 coupling
   laser moves the promotion root from \(z_\star=543/280\simeq1.94\) down to
   \(z_\star=543/1190\simeq0.456\) — well inside the low-detuning window —
   while sweeping the closing laser \(\omega_{35}\) gives \(c_1=0\): **no**
   real promotion root exists in that direction.
3. **The reduced Z-pencil is the exact weak-probe block of the physical
   Liouvillian.** With detunings implemented as rotating-frame level shifts
   \(\delta H=-z\sum_k\zeta_k|k{+}1\rangle\langle k{+}1|\), the reduced
   response equals the full 25×25 GKSL linear response **exactly** — the
   difference is zero in Gaussian-rational arithmetic, not merely small in
   floating point.

## Physical detuning signatures

Reduced coordinate \(k\) is the optical coherence \(\rho_{(k+2),1}\) of the
diamond \(|1\rangle\to|2\rangle\to\{|3\rangle,|4\rangle\}\to|5\rangle\).
Rotating-frame spanning tree: \(\varepsilon_1=0\),
\(\varepsilon_2=\omega_p\), \(\varepsilon_3=\omega_p+\omega_{23}\),
\(\varepsilon_4=\omega_p+\omega_{24}\),
\(\varepsilon_5=\omega_p+\omega_{23}+\omega_{35}\), with the loop laser
co-swept to satisfy \(\omega_{45}=\omega_{23}+\omega_{35}-\omega_{24}\)
(keeping the frame static; sweeping \(\omega_{45}\) *alone* is not
representable as a diagonal \(Z\) and is out of scope). Sweeping one
generator frequency by \(z\) shifts coherence \(\rho_{k1}\) by \(z\) times:

| swept laser | \(Z=(\zeta_{\rho21},\zeta_{\rho31},\zeta_{\rho41},\zeta_{\rho51})\) |
|---|---|
| probe \(\omega_p\) (Priority 3 convention) | \((1,1,1,1)=I\) |
| branch-3 coupling \(\omega_{23}\) | \((0,1,0,1)\) |
| branch-4 coupling \(\omega_{24}\) | \((0,0,1,0)\) |
| closing laser \(\omega_{35}\) | \((0,0,0,1)\) |
| two-tone \(\omega_{23}+\omega_{24}\) | \((0,1,1,1)\) |
| generic rational control | \((1,7/5,3/10,2)\) |

Gate Z5 certifies this table against the full Liouvillian rather than
assuming it.

## Z0: exact regression

At \(Z=I\), the three-variable numerator and denominator polynomials, and
the ideal-cut moment ladder, reproduce the Priority 3 objects **bit-exactly**
(dictionary equality of Gaussian-rational coefficients). At \(z=0\) they
reproduce the certified Phase N two-variable polynomials bit-exactly.

## Z1: the fan is exactly direction-independent

For each signature, the exact \((\Gamma,\kappa)\)-Newton upper hulls at two
independent generic rational frequencies (\(z_0=1/3,\ 2/7\)) give identical
piecewise-affine fans, and every signature returns

\[
\nu(q)=
\begin{cases}
4-q,&0\le q\le1,\\
2+q,&1\le q\le2,\\
4,&q\ge2,
\end{cases}
\]

with breakpoints exactly \(\{1,2\}\). The \(q\to\infty\) plateau equals the
first nonzero ideal-moment index (4) for every signature. High-precision
decimal cross-checks at \(\Gamma=10^{10}\), \(q\in\{1/2,3/2,5/2\}\) agree
with the exact orders to \(\le1.5\times10^{-5}\).

**Consequence.** The interference cancellation responsible for the promoted
ideal order (\(m_1=m_2=m_3\equiv0\) as \(z\)-polynomials) survives every
physical detuning direction: it is a property of the coupling graph, not of
the common-shift idealization.

## Z2: the promotion root moves with the detuning direction

Exact first nonzero moments (all signatures have \(k_0=4\)):

| signature | \(\mathfrak m_4(z;Z)\) | promotion root \(z_\star\) | in \([-3,3]\) |
|---|---|---|---|
| probe \(I\) | \(-\frac{173760}{140777}+\frac{12800}{20111}z\) | \(543/280\simeq1.939\) | yes |
| \(\omega_{23}\) sweep | \(-\frac{173760}{140777}+\frac{3200}{1183}z\) | \(543/1190\simeq0.456\) | yes |
| \(\omega_{24}\) sweep | \(-\frac{173760}{140777}-\frac{3200}{1547}z\) | \(-543/910\simeq-0.597\) | yes |
| \(\omega_{35}\) sweep | \(-\frac{173760}{140777}\) (constant) | **none** | — |
| two-tone \(\omega_{23}{+}\omega_{24}\) | \(-\frac{173760}{140777}+\frac{12800}{20111}z\) | \(543/280\) | yes |
| generic control | \(-\frac{173760}{140777}+\frac{63680}{20111}z\) | \(543/1393\simeq0.390\) | yes |

Structural observations, all exact:

- the constant term is the certified Phase N value \(m_4(z{=}0)\) and is
  **signature-independent**;
- the slope \(c_1(Z)\) is linear in \(\zeta\) and independent of
  \(\zeta_{\rho21}\) (probe and two-tone rows coincide);
- \(c_1=0\) for the \(\omega_{35}\) sweep: that detuning direction cannot
  reach the pointwise order-5 promotion at any real frequency;
- at every existing root, \(\mathfrak m_5(z_\star)\ne0\) is certified by an
  exact polynomial-gcd argument (`promotion_exactly_one_order`), so the
  promotion is by exactly one order, with no deeper cancellation.

The physically most convenient direction is the \(\omega_{23}\) sweep:
\(z_\star=543/1190\simeq0.456\) sits deep inside the certified window and
close to line center. It is used as the representative for Z3, Z4, Z6.

## Z3: frequency unfolding at the physical root

At \(z_\star=543/1190\) the exact fan becomes

\[
\nu_\star(q)=
\begin{cases}
4-q,&0\le q\le1,\\
2+q,&1\le q\le3,\\
5,&q\ge3,
\end{cases}
\]

i.e. Newton breakpoints \(\{1,2\}\to\{1,3\}\) — the same promotion pattern
as Priority 3, now along a physical detuning direction. Writing
\(\varepsilon=z-z_\star\), the certified crossover law
\(\Gamma_\times\propto|\varepsilon|^{-1/(\nu_\star(q)-\nu_{\rm gen}(q))}\)
was fitted over the asymptotic perturbations
\(\varepsilon\in\{10^{-4},10^{-5},10^{-6}\}\):

| \(q\) | fitted slope | predicted | abs. error | collapse spread |
|---:|---:|---:|---:|---:|
| 2.25 | -4.000000 | -4 | 0.000000 | 0.000000 |
| 2.50 | -2.000169 | -2 | 0.000169 | 0.000037 |
| 2.75 | -1.345200 | -4/3 | 0.011870 | 0.005133 |
| 3.00 | -0.999986 | -1 | 0.000014 | 0.000021 |
| 3.50 | -1.018056 | -1 | 0.018056 | 0.027611 |

All within the 0.05 slope gate; rescaled effective-order curves for the two
smallest \(|\varepsilon|\) collapse with spread \(\le2.8\times10^{-2}\)
(<0.06 gate). The 80-digit and 120-digit crossover locations agree to
machine bisection resolution (recorded difference 0.0) at all five \(q\).

## Z4: broadening / finite spectral resolution

\(L^2\) window norms weighted by normalized Gaussian profiles of width
\(\sigma\in\{10^{-3},10^{-2},10^{-1}\}\) centered **at** \(z_\star\) (and a
Lorentzian control at \(\sigma=10^{-2}\)), fitted over
\(\Gamma\in[10^{20},10^{26}]\):

- every window order matches the generic fan value at
  \(q\in\{0.5,1.5,2.5,3.5\}\) to \(\le1.7\times10^{-11}\) — **any** finite
  broadening restores the generic fan asymptotically, confirming that the
  order-5 result is a pointwise promotion, not a finite-band phase;
- grid refinement \(121\to241\to481\) points changes the fitted order by
  0.0 at the reported precision;
- the false-order window shrinks with resolution as
  \(\Gamma_\times(\sigma)\propto\sigma^{-1}\): fitted slope \(-1.045\)
  (predicted \(-1\); the \(\sigma=10^{-1}\) cell has no crossover bracket
  because the promoted window is already absorbed by the O(1) transient
  scale there, and it is excluded from the fit).

## Z5: exact full-Liouvillian certification

Detunings are implemented physically as level shifts
\(\delta H=-z\operatorname{diag}(0,\zeta_1,\dots,\zeta_4)\) in the 5-level
GKSL model (Hermitian Hamiltonian, nonnegative jump rates, trace
preserving), with the finite intervention \(\kappa\) as additional jumps to
the ground state, and the weak-probe steady-state response solved on the
full 25×25 Liouvillian.

- **Exact check:** at rational points
  \((\Gamma,\kappa,z)=(3,27,1/2)\) and \((10,1000,-2/3)\), for **all six**
  signatures, the full-Liouvillian response minus the reduced Z-pencil
  response is **exactly zero** in Gaussian-rational arithmetic
  (12 rows, all `difference = 0`).
- Float grid (\(z\in\{0,\pm1,z_\star,2.5\}\), two \((\Gamma,\kappa)\)
  cases, all signatures): maximum error \(1.07\times10^{-17}\); trace and
  Hermiticity residuals 0.0.

This upgrades the P8/F8-style consistency gate from "agrees to \(10^{-12}\)"
to an exact algebraic identity, and simultaneously certifies the
swept-laser table above.

## Z6: first dimensional estimate (NV-like scale)

With the model unit rate \(\gamma_0=2\pi\times10\) MHz (optical-linewidth
class):

- representative promotion root: \(z_\star=0.456\gamma_0
  \simeq2\pi\times4.6\) MHz of \(\omega_{23}\) detuning;
- resolving the pointwise promotion requires spectral selectivity below the
  \(\sigma\sim10^{-3}\gamma_0\simeq2\pi\times10\) kHz class (Z4);
- sector-response magnitude on the path \(q=2.5\) at \(z_\star\):
  \(|R_S|\simeq1.7\times10^{-5}\) at \(\Gamma=10\gamma_0\),
  \(1.1\times10^{-7}\) at \(\Gamma=31.6\gamma_0\),
  \(6.7\times10^{-10}\) at \(\Gamma=100\gamma_0\) (probe-normalized model
  units) — the measurable window sits at moderate \(\Gamma/\gamma_0\),
  consistent with the crossover-based (not deep-asymptotic) detection
  strategy of Priority 1.

## Consequence for the PRL claim

The referee-facing gap "your frequency variable is not a physical detuning"
is closed at the level of this witness:

> The path-order fan \(\nu_{S_{34}}(q)=4-q,\ 2+q,\ 4\) is exactly invariant
> under every physical detuning signature of the five-level witness, and
> the reduced pencil with signature \(Z\) is the exact weak-probe block of
> the corresponding GKSL Liouvillian with swept-laser level shifts. The
> frequency-selective order-five promotion is real, physical, and
> steerable: its root location is an exact rational function of the
> detuning direction, reaching \(z_\star\simeq0.46\gamma_0\) for a
> branch-coupling sweep, and it is absent for the closing-laser sweep.

The narrower statements remain as in Priority 3: the promotion is pointwise
in frequency (finite-band norms keep the generic fan), and the V-shaped fan
itself is a tuned-interference property with certified crossover laws away
from the manifold.

## Remaining physics tasks

1. Inhomogeneous ensemble averaging over *static parameter* disorder
   (rates, couplings), beyond the spectral broadening in Z4.
2. Noise-floor and exponent-estimation analysis for the operational
   full-minus-cut signal at the \(q=2.5\) magnitudes quoted in Z6.
3. Mapping onto one concrete double-\(\Lambda\)/NV level scheme with
   selection rules and realistic \(\kappa\) implementations (auxiliary
   dephasing beam) — prerequisite for quoting laboratory numbers beyond
   the \(\gamma_0\) scaling of Z6.

## Reproduction

```bash
pip install numpy matplotlib
python src/run_phase_z.py --smoke   # ~25 s, coarse gates
python src/run_phase_z.py           # ~270 s production, figures + JSON
```

Outputs: `results/gates_summary_phaseZ.json`,
`results/figures/figZ1_fan_per_signature.png` (fan per signature),
`figZ2_root_locus.png` (root locus along \(Z(t)=(1-t)I+tZ_{23}\)),
`figZ3_unfolding.png` (crossover slopes and collapse),
`figZ4_broadening.png` (window orders and \(\Gamma_\times(\sigma)\)).
