# Gates N1 / N2 — in-material class contrast

_Opened 2026-08-20. **Status after the full-Liouvillian repeat: N2 FAIL,
N1 undetermined.** The reduced-kernel pass reported N2 PASS; that verdict is
withdrawn below._

## Why these gates exist

Gate E rejected the general form of a measurable design rule across
engineered-dissipation platforms (`NON_CLAIMS.md` N1), and N2 records that
optical NV cannot resolve its own suppression exponent. Both failures share a
shape: they ask an experiment to measure an exponent over a decade of dynamic
range that the detection chain does not have.

Gate C appeared to offer a way around that shape. The integer classes it
certifies are not one-per-material — a single NV centre carries two of them
(`results/tables/gate_c_collapse.csv`): `ms = 0 <-> -1` at graph distance 1
(`Gamma^-2`, class 2) and `ms = -1 <-> +1` at graph distance 2 (`Gamma^-3`,
class 3). Same crystal, same bath, same `Gamma(T)`. If the class label is real,
the two channels differ by exactly one power of `Gamma`, and that difference
shows up in a *ratio*, which needs no absolute calibration and no decade of
dynamic range.

Two observables were tested: the helicity asymmetry of T6 (**N1**) and the
class-2 / class-3 channel ratio (**N2**).

Code: `calculations/numerics/class_contrast/src/`. Outputs: `../out/`.
The reduced-kernel 3E model is exec'd verbatim from
`No-go theorem/src/verify_nv_3E_graph_distance_PRL.py` by `nv3e_loader.py`.
The full nine-level model is `full_pair.build_full_pair`, which is
`gate2_candidate_full_vs_reduced.build_full` with the probe index exposed;
`full_pair.test_matches_canonical()` asserts bit-identical agreement with the
canonical builder at `(p_idx, c_idx) = (1, 2)` and passes.

---

## The result that decides both gates

**The transverse field that opens the Raman channel is already past the field
that dissolves the sectors the class label is defined on.**

The bare spin sectors are good labels only while `ge*B_perp << D_gs`, i.e.
`B_perp << D_gs/ge = 2.877/28 = 0.103 T`. The transparency island requires
`B_perp >~ 0.15 T` (manuscript Sec. 3); below `0.10 T` both thresholds collapse
to 27-33 K and there is no usable window at all.

Ground-manifold composition (`out/N2full_ground_mixing.csv`), `|<ms|dressed>|^2`
for the dressed state that is `ms = -1` at zero field:

| `B_perp` (T) | `ge*B` (GHz) | `ms=-1` | `ms=0` | `ms=+1` |
|---|---|---|---|---|
| 0.000 | 0.00 | 1.00 | 0.00 | 0.00 |
| 0.100 | 2.80 | 0.55 | 0.00 | 0.45 |
| 0.150 | 4.20 | 0.52 | 0.00 | 0.48 |
| 0.232 | 6.50 | 0.51 | 0.00 | 0.49 |
| 0.500 | 14.00 | 0.50 | 0.00 | 0.50 |

At the operating field the state is a 51/49 mixture of `ms = -1` and `ms = +1`.
There is no dressed pair that realises the `(-1,+1)` Lambda, so there is no
channel for the class-3 label to attach to.

`I_class` (`B_perp <~ 0.10 T`) and `I_island` (`B_perp >~ 0.15 T`) do not
intersect. This is the Gate E failure shape again — non-overlapping windows —
but here it follows from two fixed NV constants, `D_gs` and `ge`, rather than
from a detection budget.

---

## N2 — class ratio. Reduced kernel PASS, **full Liouvillian FAIL**.

### What the reduced kernel said

`gate_N2_channel_ratio.py`, `out/N2_channel_ratio.csv`:
`R(T)*Gamma(T) = M2/M1 = 58.4176 GHz`, constant to better than 1% for
`T >= 70 K`, with `R` swinging by a factor 5.3 across 70-101 K and both
channels far inside the integration-time ceiling. Reported PASS.

### What the full Liouvillian says

`gate_N2_full_liouvillian.py`, `out/N2full_channel_sweep.csv`. Observable is
`dA = Im chi_cut - Im chi_full` at the peak of the two-photon feature — the
unnormalised sector correction, which is what inherits the kernel exponent.

| T (K) | `Gamma` (GHz) | `dA` class 2 | `dA` class 3 | `\|R\|` | `\|R\|*Gamma` |
|---|---|---|---|---|---|
| 50 | 13.80 | +1.932e-02 | -7.39e-06 | 3.83e-04 | 5.28e-03 |
| 70 | 72.23 | +3.710e-04 | -3.37e-06 | 9.09e-03 | 6.57e-01 |
| 90 | 233.08 | +1.372e-06 | -1.22e-07 | 8.91e-02 | 2.08e+01 |
| 100 | 369.89 | +2.44e-08 | -2.74e-08 | 1.12e+00 | 4.14e+02 |
| 120 | 781.64 | -1.02e-08 | -2.13e-09 | 2.10e-01 | 1.64e+02 |

`|R|*Gamma` moves by a factor **7.8e4** across 50-120 K. It is not constant,
and `|R|` *increases* with `Gamma` where the reduced kernel required it to
fall. The class-3 channel carries the **opposite sign** to class 2 at every
temperature — control-induced absorption throughout, never transparency.

Weakening the field does not help (`out/N2full_field_sweep.csv`): at
`B_perp = 0.01 T` and `0.05 T` the spread is worse, not better, because those
fields are below the island and the class-2 channel itself collapses.

Two independent validations say the full model is the one to believe: it
reproduces the manuscript's `C = 1.4e-2` at 70 K for the `(0,+1)` channel
(computed 1.3618e-2), and it reproduces the sign reversal between 100 K and
110 K against the manuscript's `T_sign = 103 K`.

**Verdict: FAIL.** The calibration-free prediction does not survive. The
reduced-kernel PASS is withdrawn.

---

## N1 — helicity asymmetry. **Undetermined; reduced-kernel numbers withdrawn.**

The reduced kernel gave a closed form for the asymptotic asymmetry, verified to
1.1e-16 over a `(lam_perp, Delta2)` grid:

```
A_inf = lam_perp / (lam_perp + Delta2)   ->  exactly 1/2 at the literature constants
```

and showed that at zero field the moment-level null is structurally protected
against strain to 20 GHz and transverse field to 3 GHz, because `p^dag c = 0`
makes `M1 = -i p^dag H c` sample only the orbital-off-diagonal block of `H`,
from which `xi_x` and `Bx` (orbital-diagonal) and `xi_y` (orbital-off-diagonal
but spin-scalar) all drop out for `Delta m_s != 0` pairs. That algebra stands.

Everything the reduced kernel said about the *operating point* does not
(`out/N2full_helicity_asymmetry.csv`):

| T (K) | `A` full | `A` reduced |
|---|---|---|
| 70 | 0.9860 | 0.2238 |
| 80 | 0.9883 | 0.3042 |
| 90 | 0.8742 | 0.3640 |
| 101 | **-0.1542** | 0.4082 |

The two disagree by a factor of four and then by sign. The reduced-kernel
conclusion recorded in the first pass — "N1 fails as a null test, survives as a
predicted curve with 11% strain spread at 70 K" — rested on those numbers and
is withdrawn with them.

The full model is closer to a null than the reduced kernel suggested (0.986 at
70 K), but the interpretation is compromised by the mixing result above: at
`B_perp = 0.232 T` the two channels are dressed-state channels, not bare-spin
ones, so the asymmetry is real and measurable but it is not the T6 chirality
null. Whether a strain-robust null survives in the full model at the operating
field has not been computed.

**Verdict: undetermined.** Not a pass, not a clean fail.

---

## Corrections to the first pass

1. **The contrast anchor was attached to the wrong channel.** The first pass
   assumed the manuscript's 1.4e-2 at 70 K was the class-3 `(-1,+1)` channel,
   on the strength of the "experimentally standard pair" phrase in
   `verify_nv_3E_graph_distance_PRL.py`. It is not: the manuscript candidate is
   probe on `ms = 0`, control on `ms = +1`, which is class 2. Confirmed by
   direct computation (`C = 1.3618e-2` at 70 K). This was caveat 2 of the first
   pass and it turned out to matter.
2. **Reduced-kernel amplitudes were caveat 1 and they were the fatal one.**
   `NON_CLAIMS.md` N7.4 records 13 sign disagreements out of 108 points
   concentrated above 90 K. The disagreement here is larger than that record
   suggests: it is not confined to the tails, and it reverses a trend.

---

## What this means for the in-material class-contrast idea

The idea is blocked in NV, and blocked structurally rather than budgetarily.
The channel-opening mechanism and the class-labelling mechanism are the same
degree of freedom pulling in opposite directions: `B_perp` supplies the spin
overlap the orbital dipole cannot, and in doing so destroys the spin sectors
that make "class 2" and "class 3" different names for different things.

That obstruction is itself a result, and a sharper one than a failed gate: it
says the NV suppression-class structure is not experimentally separable in the
only regime where NV EIT exists at all. It belongs in the PRA's discussion of
the room-temperature no-go, not in a PRL about universality.

Nothing here is a claim. Neither gate is closed and neither appears in
`CLAIMS.md`.

## Next, if this is pursued

- The obstruction is specific to sectors labelled by *ground spin*. Group-IV
  orbital-Lambda channels are same-spin (`M0 != 0`, class 1) and need no
  transverse field, so they do not hit it — but they carry no NV experiment,
  which is what the original proposal wanted.
- A class contrast built on *orbital* rather than spin sectors would avoid the
  `ge*B_perp` vs `D_gs` collision. Nothing in the repository computes one.
- Remaining unfixed caveats from the first pass: the detection chain used
  `OD_total = OD_sector = 1`, and the strain distribution widths were scaled
  from the manuscript's `delta = 1.683 GHz` rather than measured.
