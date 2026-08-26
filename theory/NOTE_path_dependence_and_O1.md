# Note — the declared resource path as an unused knob in O1

Written: 2026-08-26.
Status: **note, not a claim.** Nothing here changes `CLAIMS.md`.
External result: `Sector-Master-Resolved-Theory`,
`paper_eit_smrt/gates/SM_R1/` (verdict PASS) and
`paper_eit_smrt/03_GENERAL_LAW_AND_REACH.md` section 2.3.

## Why this note exists

`NEXT_GATES.md` records open question **O1**:

> Is there any platform where `ν` is both large enough to matter and small
> enough to measure? `ν = 4` costs four decades of contrast per decade of `Γ`,
> which is what closes the NV window. The two requirements pull opposite ways;
> whether the tension is fundamental or an artifact of these platforms is open.

O1 is posed as a **platform-selection** question. A result obtained in the SMRT
repository shows that it is not only that: for a fixed platform, sector and
observable, the **declared resource path** `κ = κ₀Γ^q` already moves `ν` by an
integer amount. This note records that, and records equally clearly why it does
**not** reopen `NON_CLAIMS.md` N2.

## The external result, stated exactly

For the two-scale operational response `R_S^op = N_S/Q_S`, write `b_N(q)` and
`b_Q(q)` for the `κ`-degree of the monomial selected by the exact weighted
Newton degree at `q`. Then on each open stratum

```
ν_S'(q) = b_Q(q) − b_N(q),
```

and for a diagonal nonnegative cut generator on a nonsingular diagonal native
damping, `b_Q` jumps exactly once, at `q = 1`, from `0` to `m = deg_κ Q_S`. The
jump cannot be removed by tuning, because the coefficients on that Newton face
are strictly positive. Consequently `ν_S` is non-increasing on `(0,1)` and
non-decreasing on `(1,∞)`: **the fan is minimised as `q → 1`.**

Verified there against the four already-certified fans (Phase N, N6, N8 base,
N8 equal-weight control), reproduced exactly, plus 732 models with zero
violations.

## What it means for O1

In every fan computed to date the minimum sits **one integer order below the
plateau**:

| witness | plateau (`= ν_S^ideal`) | value as `q → 1` |
|---|---:|---:|
| Phase N (five-level) | 4 | 3 |
| N6 (six-level blind) | 5 | 4 |
| N8 base | 5 | 4 |

The NV quantities that O1 is about — `ν_obs = 4` for the signed absorption
difference — are **ideal-cut** quantities. In the two-scale language they are
plateau values (`q = ∞` corner; see the SMRT-side coordinate audit SM-J1, which
established that no NV certificate sweeps the `q` axis at all).

So O1's slope budget has been computed only at the corner of the fan where the
exponent is largest. Whether NV has a fan at all, and whether its interior
branches reach a smaller order, is exactly what the SMRT-side gate **SM-J2**
(NV physical model on the fan) is defined to decide. It is **unresolved**.

## What it does not mean

- **N2 stands unchanged.** "Optical NV can measure its own suppression
  exponent" remains a non-claim. The result above is about asymptotic *orders*;
  it says nothing about the detection floor, the usable window (0.51–1.02
  decade against 1 decade required), or SNR. A smaller exponent would relax the
  slope budget, but the budget is only one of the three intersecting conditions
  in Gate E (`I_engineerable ∩ I_asymptotic ∩ I_detectable`).
- **N1 stands unchanged.** Nothing here is a measurable design rule across
  platforms.
- **`q = 1` is not asserted to be optimal in general.** The unimodality result
  carries a nonsingular-damping hypothesis, and the SMRT-side gate records two
  boundaries: a singular native damping can move the minimiser off `q = 1`
  (6 of 228 cases), and an isolated exceptional value can sit at `q = 1` itself
  where the Newton face cancels exactly, so the fan need not be continuous
  there.
- **No NV fan exists yet.** No certificate in this repository sweeps `q`. The
  two `q`-resolved blocks that appear in `gates_summary_gateD.json` and
  `gates_summary_gateO.json` both evaluate the abstract Phase-N model imported
  via `phase_n_exact_core`, not an NV model.

## A candidate gate, not a permitted one

This suggests one gate, recorded here so it is not lost. It is **not** in
`gates/active/` and is therefore **not permitted work** under `NEXT_GATES.md`.

> **Candidate — path-resolved slope budget.** Re-run the Gate E window analysis
> with the exponent taken from the interior of the fan rather than from the
> ideal-cut corner, for whichever platform SM-J2 first supplies a fan for.
> Pass condition would have to be stated before running, and would have to
> include the engineerability of the path itself: holding `κ ∝ Γ` across a
> decade of `Γ` is a design assumption of exactly the kind `theory/ASSUMPTIONS.md`
> already lists as undemonstrated.

Its prerequisite is SM-J2, which has not been run.

## A finite test this repository can use now

The same result gives a decidable criterion for a question SM-J2 has to answer
anyway. The SMRT integrated theorem defines the ideal-cut order as an
ordered-limit quantity and **does not assume** that the strong-intervention
limit commutes with `Γ → ∞`. SM-J2's pass condition (b) asks whether
`lim_{q→∞} ν(q)` equals the Gate A ideal-cut value `4`. That is now decidable
without any numerical limit:

```
the fan has a plateau  ⟺  deg_κ N_S = deg_κ Q_S,
```

so comparing two `κ`-degrees of the NV operational response settles it.
