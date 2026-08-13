# DEFINITIONS.md — symbols and objects

Canonical symbol table. Where a symbol is defined formally in a source
document, that document is named; this file is the index, not the authority.

## Spaces and operators

| Symbol | Meaning |
|---|---|
| `H_g`, `H_e` | Lower (ground) and excited manifolds; `dim H_g = N_g`. |
| `P`, `P_e`, `P_g` | Orthogonal projectors, `P² = P = P†`. |
| `H` | Hamiltonian, `H = H†`. |
| `L_μ` | GKSL jump operators. |
| `𝓛` | Liouvillian: `ρ̇ = 𝓛(ρ) = −i[H,ρ] + Σ_μ (L_μ ρ L_μ† − ½{L_μ†L_μ, ρ})`. |
| `V₊ = P_e H_drive P_g = ½Ω` | Rotating-frame drive coupling; `Ω: H_g → H_e` has columns `Ω_a d_a`. |
| `d_a` | Dipole legs, `(d_a)_j = ⟨e_j| ε_a·d̂ |g_a⟩`. |
| `Γ` | Fast-dissipative rate scale; also the dissipative block `Γ = Γ† ≥ 0` of the coherence generator. |
| `A(δ) = Γ + i(H − δI)` | Passive coherence block at probe detuning `δ`. |

Source: `proofs/eit_nogo_proofs.tex` §"Notation and standing assumptions".

## Response objects

| Symbol | Meaning |
|---|---|
| `χ(δ)` | Full weak-probe linear response (susceptibility) at detuning `δ`. |
| `S` | The physically specified long-lived coherence sector. |
| `χ^(S)_cut` | Response with sector `S` cut out (the *sector cut* construction). |
| `Δχ_S = χ − χ^(S)_cut` | Sector-mediated response — the object of the no-go statement. |
| `R_S` | Sector-resolved response used in the Phase A/B campaigns. |
| `η_S(0)`, `W_S` | Sector transparency and width diagnostics reported by the gates. |
| `χ_matched` | Response with source = readout (matched readout). |
| `x` | Coherence amplitude vector solving `A(δ)x = source`. |

Source: `proofs/T1_sector_cut_axiomatization.tex` (Definition: sector cut);
`proofs/prl_eit_equivalence_conditions.md` §3, §4.

## Classification objects

| Symbol | Meaning |
|---|---|
| `ν` | Suppression index: `‖Δχ_S‖ ~ Γ^(−ν)`. |
| `ν = ∞` | Class I — exact structural no-go. |
| `ν ∈ (0,∞)` | Class II — asymptotic no-go; integer-valued. |
| `ν = 0` | Class III — protected go. |
| `μ_k` | Moments of the sector-mediated transfer; the first nonvanishing one fixes `ν`. |
| `F_Q,S` | Sector-mediated quantum Fisher information; scales as `Γ^(−2ν)` (Phase M). |
| `x_S` | Tangent-vector difference `∂_θ ρ_full − ∂_θ ρ_cut`; `‖x_S‖ ~ Γ^(−ν)`. |

## NV-specific quantities

| Symbol | Meaning |
|---|---|
| `B⊥`, `B_z` | Transverse and axial magnetic field. |
| `K ∝ B⊥²` | Quadratically opened symmetry-suppressed Raman pathway. |
| `d` | Strain parameter (`D_STRAIN = 1.683` in the candidate configuration). |
| `Ω_c` | Control Rabi frequency. |
| `γ_oc`, `γ_pop` | Orbital-coherence and population phonon rates; `γ_oc = Γ_XY/4` (rate-map corollary). |
| `k_orb(T, d)` | Orbital hopping rate at temperature `T`. |
| `C` | Contrast — the reported EIT observable. |

Source: `calculations/numerics/No-go theorem/src/nv_model.py`,
`phonon_rates.py`.
