# New No-Go Theory: Three-Tier Response Classification

This directory contains the new no-go theory that upgrades the binary EIT
go/no-go criterion (see `../No-go theorem/`) to a **three-tier classification**
of sector-mediated coherent response in finite-dimensional Markovian
weak-probe systems, following the strategy in
`EIT_new_nogo_PRX_taxonomy_and_quantum_network_strategy.md`.

## The three central theorems

| Theorem | Class | Index | Content |
|---|---|---|---|
| **I** — Sector-Resolved Transfer-Function No-Go | I: Exact Structural No-Go | ν = ∞ | δχ_S(z) ≡ 0 ⇔ ℓ†Mᵏr = 0 (k < n) ⇔ 𝒦(M,r) ⊆ 𝒩_obs(M,ℓ); finite Krylov certificate |
| **II** — Dissipative Asymptotic Hierarchy | II: Asymptotic No-Go | 0 < ν < ∞ | Uniform-in-z expansion F_Γ = Σ(−1)ⁿΓ^−(n+1)μ_n(z); suppression index ν_diss = m+1 ≤ dim, basis independent |
| **III** — Protected Channel under Singular Dissipation | III: Protected Go | ν = 0 | Non-Hermitian D, semisimple kernel, Riesz projection P; F_Γ = F₀ + Γ⁻¹F₁ + O(Γ⁻²), F₀ = p†P B_P⁻¹ Pc; EIT corollary δχ_{S,0} |

Plus: basis/realization independence (Kalman), the exclusive trichotomy
ν ∈ {∞} ∪ (0,∞) ∪ {0}, the counterexample showing Pc ≠ 0 and p†P ≠ 0 do
**not** suffice for protection, and a proof that the 3-step decision
algorithm terminates and is correct.

## Contents

- `Theorem and proofs/three_theorems_proofs.tex` — self-contained proofs
  document (cites the companion `../No-go theorem/Theorem and proofs/eit_nogo_proofs.tex`
  for reused building blocks: Theorems 2A/2B, 3, 4, 5, 6 there).
- `src/three_theorem_verification.py` — numerical verification of all three
  theorems and the counterexample lemma.

## Reproduce

```bash
pip install numpy
python "src/three_theorem_verification.py"
```

Expected output: four PASS lines (Theorem I machine-precision zero;
Theorem II fitted slope −ν_diss; Theorem III O(1) plateau at F₀ with 1/Γ
correction; counterexample F₀ = 0 with nonzero endpoint projections).
