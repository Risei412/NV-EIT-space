# Gate N4 — what theta is, and what it does to Gate N3

_Opened 2026-08-20. Status: theta is resolved. **It falsifies the central claim
of Gate N3**, which is corrected below._

Code: `calculations/numerics/class_contrast/src/gate_N4_theta.py`,
`gate_N4b_basis.py`. Outputs: `../out/N4_theta_polarization.csv`,
`N4b_class_vs_strain.csv`, `N4b_strain_required.csv`.

## The question

`group_iv_model.legs` documents `theta` as "the (unmeasured) ratio of
cross-branch to same-branch dipole matrix elements" and sets it to 0 as the
representative case. Gate N3 rested its whole argument on `theta != pi/2`.

## The answer: theta is not one number, it factorises

Both orbital doublets and the in-plane dipole transform as `E` in D3d, and
`E (x) E = A1 + A2 + E` contains `E` once, so one reduced matrix element fixes
all four orbital transitions. In the real orbital basis the two dipole
components are `d_x = A sigma_z`, `d_y = A sigma_x` — the same `E` structure the
model already uses for strain, checked directly. Then

```
d(e_c)^dag d(e_p) = (e_c . e_p) I + C (i sigma_y),   C = conj(e_cx) e_py - conj(e_cy) e_px
```

and for orthogonal ground states the identity piece drops, leaving

```
M0 = C * <g_ctrl| i sigma_y |g_probe>
```

verified on 2000 random polarization pairs in both candidate bases to 2.3e-16.
`M0` is a product of two independent factors:

1. **A polarization factor `C`**, the cross product of the probe and control
   polarizations. It is maximal for crossed polarizations and zero for parallel
   ones — a free experimental knob, not a material property. A 10 degree
   misalignment costs 1.5% (`out/N4_theta_polarization.csv`).
2. **A ground-manifold factor**, the `A_2` (angular-momentum) matrix element
   between the two ground states. **This is the one that decides the class.**

## The ground factor kills class 1 unless strain dominates spin-orbit

`i sigma_y` is the angular-momentum operator `L_z`, and the ground spin-orbit
interaction is `lambda L_z S_z`. So the spin-orbit eigenstates are exactly the
states that make the second factor vanish:

| ground basis | `\|<g_c\|sigma_y\|g_p>\|` | `\|M0\|` | class |
|---|---|---|---|
| linear `\|x>, \|y>` | 1.0000 | 1.0000 | 1 |
| circular `\|+>, \|->` (spin-orbit eigenstates) | 0.0000 | 0.0000 | >= 2 |

Strain is `E`-type and off-diagonal in the `L_z` basis, so it competes. With
`H_g = (lambda/2) sigma_y + xi sigma_x` at fixed spin
(`out/N4b_class_vs_strain.csv`):

| `xi/lambda` | 0 | 0.03 | 0.1 | 0.3 | 1 | 3 | 10 | 100 |
|---|---|---|---|---|---|---|---|---|
| `\|M0\|` | 0.000 | 0.060 | 0.196 | 0.514 | 0.894 | 0.986 | 0.999 | 1.000 |

Using the recorded ground spin-orbit splittings
(`SiV_SnV_phonon_AIC_parameters.md`):

| material | `lambda_g` (GHz) | `xi` for `\|M0\|=0.5` | `xi` for `\|M0\|=0.9` |
|---|---|---|---|
| SiV | 48 | **13.9 GHz** | **49.6 GHz** |
| SnV | 850 | 245 GHz | 878 GHz |

**Class 1 in group-IV is a strained-sample phenomenon.** Unstrained SiV and
SnV, where spin-orbit sets the ground eigenstates, are not class 1.

## This falsifies the central claim of Gate N3

Gate N3 asserted:

> group-IV: `M0 = p^dag c = cos(theta)`, a property of the dipole geometry
> alone. No term in `H` appears in it.

**That is wrong.** The legs `p` and `c` are the excited states reached *from the
two ground states*, and the ground states are eigenstates of `H_g`. The
dependence on `H` is there; `group_iv_model.legs` hides it by hard-coding the
destinations, which is why the numerical check in N3-4 showed `M0` unmoved by
strain of 5000 GHz. That check varied `H` in the model's excited manifold while
the ground states, which is where the dependence lives, were never represented.

The N3 dichotomy "class fixed by dipole geometry vs class fixed by a path
through the Hamiltonian" therefore does not stand. Both are Hamiltonian
dependent.

### What survives, in weaker form

The two systems still differ, but in the *coupling of the knobs*, not in
whether `H` matters:

| | what sets the class | what opens the channel | independent? |
|---|---|---|---|
| NV | `ge*B_perp` vs `D_gs` | `B_perp` | **no — same knob** |
| group-IV | `xi` vs `lambda_g` | nothing (channel is open) | **yes** |

NV has no free parameter: the field that supplies the spin overlap is the field
that mixes the sectors, so `I_class` and `I_island` cannot both be satisfied
(gate N2, full Liouvillian). Group-IV has one: strain is set by growth and
fabrication, independently of anything the channel needs. That is a real
difference and it is testable — **strained SiV should be class 1, unstrained
SiV and SnV should not** — but it is a statement about experimental
controllability, not about Hamiltonian invariance.

## The repository's group-IV model cannot settle this on its own

`group_iv_model.H_groupIV` writes spin-orbit as
`(Delta_e/2) kron(sigma_z_orb, sigma_z_spin)` and strain as
`xi_x kron(sigma_z_orb, I) + xi_y kron(sigma_x_orb, I)`. Spin-orbit is `A_2`
and belongs on `sigma_y`; strain is `E` and is correctly `(sigma_z, sigma_x)`.
Writing spin-orbit on `sigma_z` puts it in the same representation slot as one
strain component, so **no single orbital basis makes both terms correct at
once**. `legs` then sends both legs to a common excited ket, which is the
linear-basis answer, while `H_SO` is written as though the basis were already
diagonal in spin-orbit.

This is the dipole geometry `NON_CLAIMS.md` N4 flags as schematic. The
consequence is sharper than N4 states: it is not only the absolute response
levels that are unreliable — **the class-1 assignment itself follows from the
schematic choice rather than from D3d symmetry.** Gate C's group-IV row is one
of the two legs of its material-independence claim.

## What would settle it

1. Rebuild the group-IV ground manifold with spin-orbit on `sigma_y` and strain
   on `(sigma_z, sigma_x)`, and recompute the Gate C class assignment as a
   function of `xi/lambda_g`. This is a small change to `group_iv_model` and it
   decides whether Gate C's class-1 row survives.
2. Compare against measured strain in the SiV samples of interest. The
   threshold `xi ~ 14 GHz` for half-maximal `M0` is well inside the range that
   strain-tuning experiments reach; `xi ~ 245 GHz` for SnV is not obviously so.
3. Only then is the N3 comparison worth redoing.

## Caveats

- This is a symmetry-level derivation with one reduced matrix element. Real SiV
  carries Jahn-Teller quenching and higher-order terms that can reduce the
  off-diagonal element below its Wigner-Eckart value.
- The spin-orbit values 48 GHz and 850 GHz are recorded in the parameter file as
  "context only (not a model input)"; they have not been used in a computation
  in this repository before now.
- Nothing here is measured, and nothing here is a claim in `CLAIMS.md`.
