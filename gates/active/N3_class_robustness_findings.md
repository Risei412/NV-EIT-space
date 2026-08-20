# Gate N3 — two kinds of suppression class

_Opened 2026-08-20, after gates N1/N2 failed. Status: the structural result is
computed and holds; the experimental consequence rests on one unmeasured
parameter._

Code: `calculations/numerics/class_contrast/src/gate_N3_window_law.py`.
Outputs: `../out/N3_window_law.csv`, `N3_requirement_table.csv`.

## What is not new here

`GateE_NV_experimental_anchor/src/run_gate_e.py` defines

```python
usable_window_decades(reference, floor, order) = log10(|ref|/|floor|) / order
```

so `window_decades * order = D` is an **identity of Gate E's own definition**,
not something discovered in its output. The check in N3-1 confirms it to 1e-15
across all five scenarios, which it must.

The only thing worth restating is the consequence, which the repository draws
nowhere as a general statement: **the requirement `window >= 1 decade` is the
requirement `D >= nu`.** The detection headroom a system needs grows linearly
with its class index. Gate E reports this per scenario (1.02, 0.84, 0.51
decades) rather than as a scaling in `nu`.

Applying that to NV's own measured headroom (`out/N3_requirement_table.csv`):

| scenario | `D` | `nu=1` | `nu=2` | `nu=3` | `nu=4` |
|---|---|---|---|---|---|
| single | 4.069 | 4.07 ok | 2.03 ok | 1.36 ok | 1.02 ok |
| post_selected_shimmed | 3.356 | 3.36 ok | 1.68 ok | 1.12 ok | 0.84 **NO** |
| post_selected | 2.911 | 2.91 ok | 1.46 ok | 0.97 **NO** | 0.73 **NO** |
| low_density | 2.693 | 2.69 ok | 1.35 ok | 0.90 **NO** | 0.67 **NO** |
| high_density | 2.046 | 2.05 ok | 1.02 ok | 0.68 **NO** | 0.51 **NO** |

A class-1 system clears every scenario on the chain that fails for the orders
NV actually has.

## What is new: the class index is set in one of two ways, and they behave differently

**NV.** Probe and control sit on orthogonal orbital branches, so `p^dag c = 0`
identically and `M0 = 0`. The class is therefore fixed by the first nonzero
moment, i.e. by a path *through* `H`, with the order given by the sector-graph
distance. Any perturbation that changes `H` can change the class — and
`gate_N2_full_liouvillian.py` shows that the transverse field required to open
that path is exactly such a perturbation. At `B_perp = 0.232 T` the dressed
ground state that is `ms = -1` at zero field is a 51/49 mixture of `ms = -1`
and `ms = +1`, so the sectors the class is defined on no longer exist.

**Group-IV.** The orbital-Lambda is same-spin and probe and control share a
bright excited sublevel, so `M0 = p^dag c = cos(theta)`, where `theta` is the
leg-mixing angle. **No term of `H` appears in it.** Verified directly:

| perturbation | NV `M0` | SiV `M0` | SnV `M0` |
|---|---|---|---|
| none | 0.000e+00 | 1.000 | 1.000 |
| strain `xi_x = 100 GHz` | 0.000e+00 | 1.000 | 1.000 |
| strain `xi_y = 100 GHz` | 0.000e+00 | 1.000 | 1.000 |
| strain `xi_y = 5000 GHz` | 0.000e+00 | 1.000 | 1.000 |
| transverse field `Bx = 50 GHz` | 0.000e+00 | 1.000 | 1.000 |

Strain twenty times the SiV excited splitting, and a field an order of
magnitude past anything an experiment would apply, leave `M0` at exactly 1.
The only thing that removes group-IV from class 1 is `theta -> pi/2`, a dipole
geometry with no same-branch overlap:

| `theta` (deg) | `M0` | class |
|---|---|---|
| 0 | 1.000 | 1 |
| 60 | 0.500 | 1 |
| 85 | 8.72e-02 | 1 |
| 89.9 | 1.75e-03 | 1 |
| 90 | 6.1e-17 | >= 2 |

### The dichotomy

> A class fixed by the **dipole geometry** (`M0 != 0`) is invariant under every
> Hamiltonian perturbation and needs the least detection headroom, `D >= 1`.
> A class fixed by a **path through the Hamiltonian** (`M0 = 0`, order = graph
> distance) is destroyable by the same perturbation that opens the path, and its
> headroom requirement grows as `D >= nu`.
>
> The classes that selection rules make interesting are the fragile ones.

The window verification (N3-2) confirms `window = D/nu` on the actual kernels
and shows the asymptotic entry point drops out of it: `Gamma_a` is 15.8 GHz for
NV class 2, 23.4 GHz for NV class 3, 264 GHz for SiV and 3.1e3 GHz for SnV —
roughly `Delta_e` in the group-IV cases — but the window depends only on `D`
and `nu`.

## What this does not establish

1. **`theta` is unmeasured.** `group_iv_model.legs` documents it as "the
   (unmeasured) ratio of cross-branch to same-branch dipole matrix elements",
   and `theta = 0` is the representative choice, not a measured one. If the real
   D3d dipole tensor puts `theta` near `pi/2`, group-IV loses class 1 and the
   argument collapses. **This is the single measurement that decides the
   claim.**
2. **Group-IV's own `D` is not computed.** The table above uses NV's measured
   headroom to illustrate the `D >= nu` requirement. Group-IV's actual headroom
   depends on its dipole strength, optical depth and phonon rate, none of which
   are refit (`NON_CLAIMS.md` N4).
3. **`Gamma -> T` needs the refit.** The asymptotic entry points are stated in
   `Gamma`, which is N4-safe. Converting them to temperatures requires the
   schematic phonon normalisation that N4 forbids quoting.
4. **No experimental leg.** Nothing here is measured; there is no group-IV
   experiment in this repository.
5. **Not yet checked against the Gate E audit.** `NON_CLAIMS.md` N1 forbids
   resurrecting the rejected general design rule by rewording. The claim here is
   close to its negation — an obstruction to measurability rather than a rule
   for measuring — and the evidence that killed N1 supports it. That reading is
   the author's and has not been checked against
   `gates/closed/GENERAL_PRL_GATE_AUDIT.md` line by line.

## Next

Measure or compute `theta` from the D3d dipole tensor. Everything else waits on
it.
