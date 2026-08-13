# Dark-State-Free Coherent Transparency in a Minimal 2g+2e Markovian Model

## 1. Main result

A minimal four-level model with two lower states and two excited states was constructed with

- no optical dark state (`rank Ω = 2`);
- no stationary pure Lindblad dark state for nonzero probe and control fields;
- a strong sector-mediated absorption minimum;
- no Autler–Townes pole splitting near the transparency window.

For the representative parameter set, the absorption is reduced by **79.03%** relative to the sector-cut response. However, the exact transfer zero remains below the real-frequency axis. A determinant identity proves that, in a regular ideal 2g+2e model with zero ground-coherence decay, a perfect real-frequency transparency zero is impossible when `rank Ω = 2`.

Thus the minimal model supports **dark-state-free EIT-like coherent transparency**, but not perfect dark-state-free EIT.

## 2. Hilbert space, Hamiltonian and jumps

Use the basis

`|g1>, |g2>, |e1>, |e2>`.

The probe and control dipole vectors in the excited manifold are

- `d_p = (1,0)^T`,
- `d_c = (cos θ, sin θ)^T`.

The coupling matrix is

`Ω = [Ω_p d_p, Ω_c d_c]`.

For `θ = π/4`,

- `det[d_p,d_c] = 0.707107`,
- singular values are `1.306563` and `0.541196`,
- therefore `rank Ω = 2` for every nonzero `Ω_p Ω_c`.

The Markovian jump set may be chosen as

`L_(jα) = sqrt(γ_(jα)) |g_α><e_j|`

for all `j=1,2` and `α=1,2`, supplemented by ground dephasing and the Γ-scaled excited-state dephasing jumps used below.

### No stationary pure dark state

Every radiative jump is nilpotent and has only eigenvalue zero. A common eigenvector of all radiative jumps must therefore lie entirely in the ground subspace. For a stationary pure Lindblad state, the effective-Hamiltonian condition then requires

`P_e H_eff |D> = Ω |D>/2 = 0`.

Because `rank Ω = 2`, its kernel is trivial. Hence no nonzero pure stationary state exists for nonzero probe and control fields.

At the strict reference point `Ω_p=0`, `|g1>` is a trivial control-dark state. The dark-state-free statement concerns every finite probe amplitude and its weak-probe limit. This order-of-limits caveat must be stated in a paper.

## 3. Exact weak-probe response

Define

`a1 = γ_e + Γ_1 - iδ`,

`a2 = γ_e + Γ_2 + iΔ_e - iδ`,

`g = γ_g + Γ_g - iδ`,

`β = |Ω_c|^2/4`,

`c = cos θ`, `s = sin θ`.

The exact reduced response is

`χ_cut^(S) = 1/a1`,

`χ_full = (g a2 + β s^2) / [g a1 a2 + β(c^2 a2 + s^2 a1)]`,

`R_S = χ_full - χ_cut^(S)`

`    = -β c^2 a2 / {a1 [g a1 a2 + β(c^2 a2 + s^2 a1)]}`.

The sector cut removes the feedback through the ground-coherence variable while preserving the optical block, source and detector.

## 4. A new minimal no-zero theorem

Let `D=[d_p,d_c]` and `G=A^-1`. With

`S_p=d_p†Gd_p`, `S_c=d_c†Gd_c`,
`K_pc=d_p†Gd_c`, `K_cp=d_c†Gd_p`,

one has the exact 2×2 identity

`S_p S_c - K_pc K_cp = det(D† G D)`
`                         = |det D|^2 det G`.

If `rank D=2` and `A` is invertible, the right-hand side is nonzero. At ideal two-photon resonance `g=0`,

`χ_full = (S_p S_c-K_pc K_cp)/S_c`.

Therefore a regular perfect transfer zero requires `det D=0`, which is precisely the optical-dark-state condition. This obstruction is special to the minimal two-dimensional excited manifold and does not automatically extend to three or more excited states.

## 5. Numerical coherent-transparency example

Parameters, in units of `γ_e`, are

- `γ_e = 1.0`,
- `γ_g = 0.02`,
- `Δ_e = 8.0`,
- `Ω_c = 0.8`,
- `θ = π/4`.

The local minimum occurs at

- `δ_min/γ_e = -0.010000`,
- `Re χ_full = 0.209709460`,
- `Re χ_cut = 0.999900010`,
- relative absorption suppression = `79.027%`.

The nearest complex transfer zero is

`z_0/γ_e = -0.009840765 -0.021202535 i`.

It is not on the real axis, but its small imaginary distance produces the deep absorption minimum.

## 6. Pole/residue analysis and ATS exclusion

The two near-resonant poles are

1. `z_1/γ_e = 0.001123215 -0.910127975 i`,
2. `z_2/γ_e = -0.010975891 -0.111077512 i`.

Their real-frequency separation is

`|Re z_1-Re z_2| = 0.012099105 γ_e`,

whereas their halfwidths are `0.910128 γ_e` and `0.111078 γ_e`. The poles are therefore not spectrally split.

The corresponding residues are

1. `r_1 = 0.003123414 +1.112428284 i`,
2. `r_2 = -0.003123045 -0.112429773 i`.

At the absorption minimum their real contributions are

- broad/near pole contribution: `1.222052296`,
- narrow/near pole contribution: `-1.012342903`.

The contributions have opposite signs and cancel. This is a Fano/EIT interference mechanism, not a pair of independently absorbing Autler–Townes peaks.

## 7. Γ-dependent Class I–III map

The master response is classified, not the total susceptibility.

### Class I: exact structural zero

Choose `θ=π/2`, so `d_p=(1,0)` and `d_c=(0,1)`. The coupling matrix still has rank two, but `c=0`, hence

`R_S ≡ 0`

for every frequency and Γ. Thus `ν=∞`.

### Class II: algebraic dissipative suppression

Damp both optical coherences:

`Γ_1=Γ_2=Γ`, `Γ_g=0`.

For fixed nonzero `g`,

`R_S = -β c^2/(g Γ^2) + O(Γ^-3)`,

so `ν=2`. The numerical log-log slope is `-1.999099`.

If all three response coherences are scaled, including the ground coherence,

`R_S=O(Γ^-3)`,

with numerical slope `-2.999817`.

### Class III: protected survival

Damp only the orthogonal excited channel:

`Γ_2=Γ`, `Γ_1=Γ_g=0`.

Then

`R_S -> -β c^2 / [a1(g a1+β c^2)]`,

which is nonzero. The numerical slope is `0.000072`, i.e. order Γ^0. The source, readout and common Raman route live in the kernel of the scaled dissipator.

## 8. Mapping to the NV center

A local NV mapping is

- `|g1>,|g2>`: two ground-spin sublevels;
- `|e1>,|e2>`: two selected excited orbital-spin branches;
- `d_p,d_c`: polarization- and mixing-dependent dipole vectors;
- `c=cos θ`: the common excited-manifold component opened by spin-orbit, spin-spin, hyperfine, strain or transverse-field mixing;
- `Γ`: optical-coherence damping generated by excited-state orbital hopping.

In the symmetry-preserving zero-field limit, the two spin legs can be orthogonal (`c=0`), giving the simplified Class-I result. A transverse field or another mixing term produces `c∝B_perp`, and the full-optical Class-II asymptotic law gives

`R_S ∝ B_perp^2 Γ_oc^-2`

for fixed ground-coherence damping.

For symmetric NV orbital hopping, the optical-coherence rate is `Γ_oc=Γ_XY/4`. Both selected optical coherences are damped, so the standard NV mapping is Class II rather than Class III. A Class-III NV realization would require an engineered dissipative kernel, such as a genuinely protected optical-coherence combination.

The actual NV excited manifold is six dimensional. This 2e mapping is therefore a local mechanism and scaling model, not a final material-level proof. The next validation step is the full six-level excited-manifold Liouvillian with the same sector-cut and pole/residue audit.

## 9. Final verdict

The requested program succeeds with an important qualification:

1. `rank Ω=2` and absence of a stationary pure dark state are proved.
2. A strong dark-state-free coherent transparency dip is obtained.
3. Pole positions and opposite-sign residues exclude ATS for the representative regime.
4. The same minimal model realizes Class I, II and III under distinct coupling/dissipation paths.
5. A new obstruction theorem shows that perfect dark-state-free EIT is impossible in the regular ideal 2g+2e model.

Consequently, an exact dark-state-free EIT theorem should next be sought in a `2g+3e` model, a closed-loop model, or a full NV excited manifold. The present 2g+2e result supplies both the benchmark and the no-go boundary.
