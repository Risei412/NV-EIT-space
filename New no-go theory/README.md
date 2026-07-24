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
- `Theorem and proofs/no_go_theory_undergrad_guide.tex` — pedagogical guide
  to the whole theory at undergraduate level (only linear algebra and
  geometric series assumed): intuition, analogies, exactly solvable 2×2
  examples for each class, the trichotomy, the decision algorithm, scope,
  and a jargon dictionary.
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

## Phase P: physical Lindblad realization of interference-controlled exponent promotion (Paper III)

A physical CPTP realization of the abstract Phase-D order-promotion result
(`ν: 3→4`), implementing the plan in `paper3_smrt_numerical_plan.md`.

- `src/model_physical.py` — a 5-level "diamond" Lindblad model (two coherent
  paths into a readout state, sector cut on both closing couplings), with
  reduced 4×4 weak-probe machinery, exact symbolic moments/certificate, and
  the full vectorized 24×24 Liouvillian with an implicit linear-response
  solver.
- `src/run_phase_p.py` — runs P1–P8, writes
  `results/gates_summary_phaseP.json` and `results/figures/figP*.png`.

```bash
python "src/run_phase_p.py"
python "src/report.py"
```

Results (all eight gates PASS): appended to
[`results/summary.md`](results/summary.md).

Key findings: the cancellation condition `J45* = J23 J35 d4/(J24 d3)` at
`φ=π` gives exact `m2=0, m3≠0` (ν: 3→4), matching the direct large-Γ fit to
`9×10⁻⁷`; the universal crossover collapse (`Γ⁴R_S` vs `δΓ`) holds within a
`3%` spread and `Γ_×∝|δ|⁻¹` to slope `−1.06`; a false `ν≈4` plateau spans
`2.15` decades at `δ=10⁻⁵`; and the 4×4 reduced system agrees with the
full 24×24 Liouvillian to `~10⁻¹²`–`10⁻⁶`. One correction to the original
plan: with only 2 real controls `(|J45|,φ)`, the cancellation `m2=0` is an
*isolated* point, not a codimension-1 curve (matching the precedent in gate
M4 above) — a genuine cancellation curve needs a 3rd real control, shown
explicitly by also rescaling the first branch's amplitude.

## Phase Z: physical detuning signature zZ (Phase N Priority 4)

Closes the gap flagged by the Phase N Priority 3 report: a physical EIT
implementation places laser detunings with different weights on different
coherences, so the common resolvent shift `-izI` must be replaced by a
diagonal signature `-izZ`, `Z = diag(zeta)`, derived from which laser is
swept (rotating-frame spanning tree of the 5-level diamond; the loop laser
is co-swept). Gates Z0-Z6 certify, exactly where possible:

- **Z0** bit-exact regression to the Priority 3 polynomials at `Z = I` and
  to the Phase N polynomials at `z = 0`;
- **Z1** the generic V-shaped path-order fan `nu(q) = 4-q, 2+q, 4` is
  preserved *exactly* for every physical detuning signature tested;
- **Z2** the frequency-promotion root `z*` moves with the detuning
  direction (exact rational values per signature; sweeping the closing
  laser `omega35` has *no* real promotion root — its `m4(z)` is constant);
- **Z3** frequency-unfolding crossover laws and collapse at the physical
  root `z* = 543/1190` of the `omega23` sweep;
- **Z4** finite broadening/resolution restores the generic fan
  (Gaussian/Lorentzian window norms), with `Gamma_x ∝ sigma^-1`;
- **Z5** reduced Z-pencil vs. full 25x25 GKSL Liouvillian with explicit
  level-shift detunings: difference **exactly zero** in Gaussian-rational
  arithmetic (plus float grid `<1e-11`);
- **Z6** first dimensional (NV-like) estimate of `z*` and the required
  spectral selectivity.

- `PhaseZ/src/phase_z_detuning_core.py` — Z-weighted exact 3-variable
  pencil, moment ladder, Sturm real-root machinery, weighted Newton fans.
- `PhaseZ/src/run_phase_z.py` — gates Z0-Z6, figures, JSON summary.
- `PhaseZ/PHASE_Z_DETUNING_REPORT.md` — production report.

```bash
pip install numpy matplotlib
python "PhaseZ/src/run_phase_z.py" --smoke   # ~25 s
python "PhaseZ/src/run_phase_z.py"           # production
```
