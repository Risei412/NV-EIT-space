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

## Metrological (QFI) extension: does ν also classify sensing information?

A further calculation tests the candidate prediction that the tangent-vector
difference x_S = ∂_θρ_full − ∂_θρ_cut obeys ‖x_S‖ ∼ Γ^−ν, and that the
sector-mediated quantum Fisher information F_{Q,S} = x_S^† G_ρ x_S ∼ Γ^−2ν
— i.e. that the same three-tier classification governs not just response
magnitude but how much parameter-estimation information a sector carries
under strong dissipation.

- `src/model_metro_linear.py` — abstract vector-valued generalization of
  the theorems, plus a rank-1 (Sherman–Morrison) sector cut giving an
  **exact** vector hidden-class transition in closed form (no root search).
- `src/model_metro_lindblad.py` — a genuine 3-level Λ Lindblad master
  equation: vectorized steady state, implicit-differentiation tangent
  vector, and SLD quantum Fisher information.
- `src/run_phase_m.py` — Gates M1–M5, Figures M1–M5.

```bash
python "src/run_phase_m.py"
python "src/report.py"
```

Results (Gates M1, M2, M3, M5 PASS; M4 is an informative negative result
— see below): appended to
[`results/summary.md`](results/summary.md).

Key findings:
- **The ν → 2ν QFI translation holds** in a genuine 3-level Lindblad model:
  at a generic point, ‖x_S‖ ∼ Γ^{−2.00} and F_{Q,S} ∼ Γ^{−3.99} — the
  ratio ν_F/(2ν_x) = 1.0007 (Gate M3, Figure M3).
- **The vector hidden class transition is exact** in the abstract
  (non-Hermitian, rank-1 cut) arena: ν[x_S] jumps 0 → 1 at a
  closed-form λ_c, with Γ·‖x_S‖ vs. Γ(λ−λ_c) collapsing across two decades
  of Γ (Gate M2, Figure M2).
- **A genuine physical (Hermitian) coupling cut does *not* exhibit an
  interior transition** with only a 2-parameter (amplitude, phase)
  control — a 2D search finds no dip beyond the trivial edge (Gate M4,
  Figure M4). This is explained, not contradicted, by the theory: the
  leading-order response lives in a 3-real-dimensional traceless-Hermitian
  slice of the protected block, so codimension counting predicts an
  interior zero needs ≥3 independent real controls. The theory correctly
  predicts *when* a hidden transition is physically achievable.
- **The full-rank / non-singular-SLD-metric assumption is load-bearing**:
  forcing the steady state toward a near-pure (rank-deficient) regime
  (ε → 10⁻⁶) breaks the clean ν_F = 2ν_x relation (ratio drifts to 1.49),
  confirming Gate M5 as a genuine negative control rather than a
  tautology.
