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

## Sector-resolved response calculation (`new_nogo_numerical_priorities.md`)

A second, independent calculation implements the priorities in
`new_nogo_numerical_priorities.md`: it computes χ_full, the frozen-source
sector-cut χ_cut, and the difference R_S = χ_full − χ_cut together (never
χ_full alone), extracts the suppression index ν three independent ways
(direct log–log fit, moment method, protected-coefficient method), searches
for a *hidden class transition* — a parameter point where ν[χ_full] stays 0
while ν[R_S] jumps 0→1 — and checks the associated Γ(λ−λ_c) scaling
collapse.

- `src/core.py` — shared infrastructure (transfer function, Krylov
  certificate, moment method, Riesz-projection protected coefficient,
  log–log ν fit).
- `src/model_lambda.py` — Phase A: standard 2×2 Λ model (analytic and
  numeric cross-checked).
- `src/model_protected.py` — Phase B: Class I/II/III unit-test models and
  the singular-D full/cut transition model.
- `src/run_phase_a.py` — Figure 2 (EIT–ATS crossover maps), Figure 5
  (mechanism map), Gate 3.
- `src/run_phase_b.py` — Figure 1 (unit test), Figure 3 (hidden class
  transition), Figure 4 (scaling collapse), Gates 1/2/4.
- `src/report.py` — assembles `results/summary.md` from the two gate JSON
  files.

```bash
pip install numpy scipy matplotlib
python "src/run_phase_a.py"
python "src/run_phase_b.py"
python "src/report.py"
```

Results (all four automated gates currently PASS): see
[`results/summary.md`](results/summary.md) and `results/figures/`.

Key finding: in the minimal 2×2 Λ model, the sector cut severs the only
coupling responsible for both the EIT notch and the ATS splitting, so
χ_cut is always the bare Lorentzian — R_S alone carries the coherent
structure, and the Anisimov–Kocharovsky/Giner discriminant
Ω_c ≷ |γ31−γ21| separates the EIT-dominant and ATS-dominant regions (Figure
5). In the extended singular-D model (Phase B), a hidden class transition
is realized explicitly: χ_full remains Class III (ν=0, O(1) plateau) on
both sides of a tuned control parameter λ_c, while R_S transitions from
ν≈0 to ν≈1 exactly at the root of r₀(λ) = δχ_{S,0}(λ), with Γ·R_S vs.
Γ(λ−λ_c) collapsing onto a single curve across a decade of Γ (Figure 4).
