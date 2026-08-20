# Gates N1 / N2 — in-material class contrast

_Opened 2026-08-20. Status: N1 FAIL (as posed) / N2 PASS (reduced kernel)._

## Why these gates exist

Gate E rejected the general form of a measurable design rule across
engineered-dissipation platforms (`NON_CLAIMS.md` N1), and N2 records that
optical NV cannot resolve its own suppression exponent. Both failures share a
shape: they ask an experiment to measure an exponent over a decade of dynamic
range that the detection chain does not have.

Gate C already contains a way around that shape. The integer classes it
certifies are not one-per-material — **a single NV centre carries two of them**
(`results/tables/gate_c_collapse.csv`):

| class | channel | mechanism | scaling |
|---|---|---|---|
| 2 | `ms = 0 <-> -1` | graph distance `d = 1` | `Gamma^-2` |
| 3 | `ms = -1 <-> +1` | graph distance `d = 2` | `Gamma^-3` |

Same crystal, same phonon bath, same `Gamma(T)`, same detection chain. If the
class label is real, the two channels must differ by exactly one power of
`Gamma`, and that difference is visible in a *ratio*, which needs no absolute
calibration and no decade of dynamic range.

These two gates test the two observables that would carry such a claim.

- **N1** — the helicity null of T6 in `verify_nv_3E_graph_distance_PRL.py`:
  `lam_perp` alone opens `(0,+1)` and leaves `(0,-1)` exactly zero.
- **N2** — the class-2 / class-3 channel ratio.

Code: `calculations/numerics/class_contrast/src/`. Outputs: `../out/`.
The NV 3E model is loaded verbatim from
`calculations/numerics/No-go theorem/src/verify_nv_3E_graph_distance_PRL.py`
via `nv3e_loader.py`; nothing is duplicated.

---

## N1 — helicity asymmetry. **Does not pass as posed.**

Observable: `A = (|K(0,+1)| - |K(0,-1)|) / (|K(0,+1)| + |K(0,-1)|)`.

### What holds

**A closed form for the asymptotic asymmetry**, verified to 1.1e-16 over a
grid in `(lam_perp, Delta2)`:

```
A_inf = lam_perp / (lam_perp + Delta2)
```

At the literature constants `lam_perp = Delta2 = 0.20 GHz` this gives
`A_inf = 1/2` exactly. Both constants are known, so this is parameter-free.

**At zero field the moment-level null is structurally protected.** `A` is
unchanged by strain up to 20 GHz on either component and by a transverse field
up to 3 GHz. The mechanism is not accidental: since `p` and `c` sit on
different orbital branches, `p^dag c = 0` and

```
M1 = -i p^dag H c
```

samples only the orbital-off-diagonal block of `H`. Of the perturbations,
`xi_x` (`sz_o (x) I3`) and `Bx` (`I2 (x) Sx`) are orbital-diagonal and drop
out; `xi_y` (`sx_o (x) I3`) is orbital-off-diagonal but spin-scalar, so it
drops out for the `Delta m_s != 0` pairs.

### What breaks it

**The null is a zero-field statement, and zero field is outside the island.**
Genuine transparency requires `B_perp >~ 0.15 T` (manuscript Sec. 3). At the
operating field the finite-`Gamma` asymmetry departs strongly from `1/2`
(`out/N1_asymmetry_vs_T_B.csv`):

| `B_perp` (T) | 60 K | 70 K | 80 K | 90 K | 101 K | 110 K |
|---|---|---|---|---|---|---|
| 0.000 | 0.4991 | 0.4998 | 0.4999 | 0.5000 | 0.5000 | 0.5000 |
| 0.150 | 0.1835 | 0.2798 | 0.3536 | 0.4029 | 0.4366 | 0.4539 |
| 0.232 | 0.1314 | 0.2238 | 0.3042 | 0.3640 | 0.4082 | 0.4321 |
| 0.500 | 0.0553 | 0.1295 | 0.2064 | 0.2757 | 0.3362 | 0.3731 |

**And once the field is on, strain matters.** Ensemble spread at
`B_perp = 0.232 T` (`out/N1_ensemble_spread.csv`):

| T | `sigma_strain = 1.683 GHz` | `sigma_strain = 5 GHz` |
|---|---|---|
| 70 K | 0.205 +- 0.022 (11%) | 0.158 +- 0.062 (39%) |
| 80 K | 0.287 +- 0.023 (8.0%) | 0.229 +- 0.068 (29%) |
| 90 K | 0.354 +- 0.013 (3.5%) | 0.301 +- 0.061 (20%) |
| 101 K | 0.401 +- 0.009 (2.3%) | 0.368 +- 0.043 (12%) |

### Verdict

**FAIL as a null test.** What survives is weaker but not worthless: `A(T,
B_perp)` remains a free-parameter-free prediction, and the strain spread
shrinks toward the warm end of the island — 2.3% at 101 K for
`sigma = 1.683 GHz`, which is also where the sign reversal sits. Testing it
requires characterising the sample's strain distribution, which is an added
experimental burden the original null-test framing did not carry.

---

## N2 — class-2 / class-3 channel ratio. **PASS (reduced kernel).**

Because the two channels differ by exactly one graph-distance step, their
ratio carries exactly one power of `Gamma`:

```
R(T) * Gamma(T) = M2 / M1 = 58.4176 GHz     (parameter-free)
```

with `M1 = 0.141421 GHz^2` (class 2) and `M2 = 8.261500 GHz^3` (class 3).

`out/N2_channel_ratio.csv`:

| T (K) | `Gamma` (GHz) | `R` | `R*Gamma` | deviation |
|---|---|---|---|---|
| 40 | 4.53 | 5.436 | 24.63 | -57.8% |
| 60 | 34.07 | 1.668 | 56.83 | -2.7% |
| 70 | 72.23 | 0.804 | 58.05 | -0.63% |
| 80 | 136.04 | 0.429 | 58.31 | -0.18% |
| 90 | 233.08 | 0.251 | 58.38 | -0.07% |
| 101 | 385.97 | 0.151 | 58.40 | -0.03% |
| 120 | 781.64 | 0.075 | 58.41 | -0.01% |

Two things matter here. `R*Gamma` is constant to better than **1% for
T >= 70 K**, and `R` itself swings by a factor of **5.3 across 70-101 K** — a
large, calibration-free signature. Below 60 K the ratio is pre-asymptotic and
the prediction does not apply.

**Both channels clear the floor simultaneously.** Using the repo's own signal
chain (`signal_chain.py`, parameters from
`results/tables/signal_chain_parameters.csv`), OD-matched sample, SNR = 5,
technical floor 1e-6, contrast anchored to the manuscript's 1.4e-2 at 70 K for
the `(-1,+1)` channel:

| T (K) | `tau` class 2 | `tau` class 3 |
|---|---|---|
| 70 | 0.69 us | 1.1 us |
| 90 | 75 us | 1.2 ms |
| 101 | 0.56 ms | 25 ms |
| 120 | 9.5 ms | 2.1 s |

Everything is orders of magnitude inside the one-hour ceiling, across the whole
island and beyond it.

### Verdict

**PASS**, subject to the caveats below. This is the observable that carries an
in-material class-contrast claim: it needs no absolute calibration, no decade
of dynamic range, and it is a ratio, so optical depth, defect density,
collection efficiency and probe power cancel.

---

## Caveats that must be cleared before any of this is quoted

1. **Reduced kernel only.** Both gates use the v6.2 reduced amplitude kernel
   (`D = I`). `NON_CLAIMS.md` and the manuscript audit record that the reduced
   kernel disagrees in sign with the full Liouvillian at 13 of 108 points,
   concentrated below 40 K and above 90 K. The `R*Gamma` result must be
   repeated with the full nine-level Liouvillian, in particular across
   90-101 K where the prediction is quoted as tightest.
2. **One-point contrast anchor.** The class-3 contrast scale is anchored to the
   single manuscript value 1.4e-2 at 70 K, on the assumption that the
   manuscript's "experimentally standard pair" is `(-1,+1)`. That assumption is
   taken from the docstring of `verify_nv_3E_graph_distance_PRL.py` and has not
   been re-derived here.
3. **Detection chain simplified.** `OD_total = OD_sector = 1` was used; the
   manuscript's chain carries more structure.
4. **Strain distribution assumed.** The `sigma` values are scaled from the
   manuscript's `delta = 1.683 GHz`; a real ensemble distribution has to come
   from the sample.
5. **Nothing here is a claim yet.** Neither gate is closed and neither appears
   in `CLAIMS.md`.

## Next

- Repeat N2 with the full Liouvillian across 70-110 K (clears caveat 1, the
  only one that can overturn the result).
- Decide whether N1 is worth carrying in its weakened form, or dropped in
  favour of N2 plus the already-computed sign reversal.
