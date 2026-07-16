# PRX Theorem Plan — Upgrading the No-Go Theory to Three Central Theorems

Companion to `SIMULATION_PLAN.md` (PRL campaign, Gates 1–5 all passed).
This plan specifies how to reorganize and extend the existing theorem
package in `No-go theorem/Theorem and proofs/` into the three PRX central
theorems of the strategy roadmap (§7), what is already proven, what must be
newly proven, and how each claim is audited numerically with the existing
code base.

**PRX acceptance criteria this plan must satisfy (roadmap §14):**

1. The three central theorems are mathematically rigorous.
2. The formulation is basis-independent.
3. At least one realistic physical example of a protected channel exists.
4. Several materials are explained by one classification (A–D).
5. The result is a conceptual classification, not just a computational method.
6. The difference from dark-state conditions, control theory, and
   decoherence-free-subspace theory is explicit.

---

## 0. Inventory: what the current package already proves

`eit_nogo_proofs.tex` (929 lines, elementary self-contained proofs) and
`EIT_no_go_go_theory_v6_2_English.tex` (§1–23) contain:

| Existing result | Statement (abridged) | Feeds PRX Theorem |
|---|---|---|
| Thm 1A/1B | optical dark-subspace rank; stationary pure Lindblad dark state | I (contrast: what no-go is NOT) |
| Thm 2A/2B | reduced susceptibility = Schur complement; sector-cut block formula δΞ_S = −βK₁₂K₂₁/(γ_g+βS₂) | I (core object) |
| Thm 3 + sector-graph remark | K_Γ = M_m Γ^{−(m+1)} + O(Γ^{−m−2}) for scalar damping ΓI; path-length selection rule M_n = 0 for n < d | II (scalar-damping special case) |
| Thm 4 | exact transfer zero ⇔ first N moments vanish ⇔ Krylov orthogonality (Cayley–Hamilton), at fixed z | I (finite certificate, fixed z) |
| Thm 5 | symmetry-protected transfer zero via commuting projectors | I & II (coefficient-vanishing mechanism) |
| Thm 6 | singular damping: K_Γ = p†P(PA₀P)⁻¹Pc + O(Γ⁻¹), O(1) channel iff leading scalar ≠ 0 | III (leading order) |
| Thm 7 | adjugate identity, analytic domains, singular points | I (analyticity input) |
| v6_2 §9–10 | full-Liouvillian weak-probe susceptibility and observable no-go | I (connects δχ_S to observables) |
| v6_2 §11–13, 18 | symmetry/group analysis, taxonomy, decision procedure, "final theorem package" | classification chapter |

Numerical anchors already archived (reuse as theorem audits):
moment classes n = m+1 ∈ {1, 2, 3} across NV/SiV/SnV in one pipeline
(`moment_order_common_pipeline.py`, slopes −1.000/−2.000/−3.000),
B⊥² opening in reduced and full Lindblad models, and the Gate 1–5 outputs.

---

## 1. Theorem I — Sector-Resolved Transfer-Function No-Go Theorem

**Target statement.** For a finite-dimensional Markovian weak-probe system,
with sector S and the sector-mediated transfer function
T_S(z) = K₁₂(z)K₂₁(z) (block form: the full Schur-complement transfer map),
the following are equivalent on the analytic domain:
(i) δχ_S(z) ≡ 0; (ii) T_S(z) ≡ 0; (iii) a finite Krylov certificate holds;
(iv) the (input, output) pair fails joint reachability/observability of the
sector in the sense of a minimal realization.

### Lemmas to prove (new work)

- **I.1 Rationality and identity theorem.** K₁₂(z), K₂₁(z), δχ_S(z) are
  rational in z with degree bounded by the dimension (input: Thm 7 adjugate
  identity). Hence "≡ 0 on the domain" ⇔ "= 0 at any 2N+1 distinct regular
  points" — this turns exact no-go into a finitely checkable property and
  cleanly separates **identical zeros (no-go)** from **isolated spectral
  zeros (accidental dark points)**. *Proof route:* adjugate/cofactor
  expansion + fundamental theorem of algebra. Difficulty: low.
- **I.2 z-resolved finite-moment certificate.** Extend Thm 4 (fixed z,
  Γ-resolvent) to the z-resolvent: δχ_S(z) ≡ 0 iff the first N moments of
  the z-expansion vanish at one regular point (derivatives up to order N−1),
  N = sector dimension. *Proof route:* Cayley–Hamilton exactly as Thm 4
  applied to X(z₀) plus I.1. Difficulty: low–medium.
- **I.3 Control-theory correspondence.** Recognize K(z) = C(zI−A)⁻¹B as a
  transfer function; show the Krylov condition of I.2 is exactly Kalman
  non-reachability/non-observability, and state the no-go as "the sector is
  removed in the minimal realization of (A, B, C)". Prove equivalence with
  the Kalman decomposition; identify controllability/observability Gramian
  rank conditions valid for stable (Hurwitz) Liouvillian blocks.
  *Proof route:* standard linear-systems theory, restated and proved
  self-containedly in the repo's elementary style. Difficulty: medium
  (mostly careful bookkeeping; the payoff is criterion 6 — explicit
  difference from dark-state conditions: no-go = non-minimality, not
  kernel existence).
- **I.4 Basis independence.** All conditions in I.1–I.3 are invariant under
  (a) unitary changes of system basis, (b) similarity transformations of the
  Liouvillian, (c) rescaling of input/output vectors; the sector S enters
  only through its projector. *Proof route:* transformation bookkeeping;
  restate theorem purely in terms of (P_S, L, |in⟩, ⟨out|). Difficulty: low.
- **I.5 Dark-state demarcation.** A corollary contrasting Thm 1A/1B:
  examples where a dark state exists but T_S ≢ 0 (EIT possible) and where
  no dark state exists but T_S ≡ 0 (structural no-go regardless).
  NV at B⊥ = 0 is the physical instance of the latter half. Difficulty: low
  (constructive examples, one 3-level and NV itself).

### Numerical audit (new script `src/prx_theorem1_certificate.py`)

- Verify I.1/I.2 on NV and group-IV Hamiltonians: evaluate the N-point /
  N-derivative certificate and compare against a dense z-scan of |T_S(z)|
  (agreement = certificate never misclassifies; include the B⊥ = 0 NV case
  as the "certified identical zero" instance and B⊥ ≠ 0 as certified
  nonzero).
- Verify I.3: rank of the Kalman controllability/observability matrices vs
  the certificate verdict on the same models + 20 random sector models
  (seeded) as a falsification sweep.

---

## 2. Theorem II — Dissipative Asymptotic Hierarchy Theorem

**Target statement.** With A_Γ(z) = ΓD + A₀(z), D full rank, the sector
response obeys δχ_S = c_m(z) Γ^{−(m+1)} + O(Γ^{−m−2}) where m is the first
index with nonvanishing generalized moment
M_n = p†(D⁻¹A₀)ⁿD⁻¹c, and the hierarchy classifies materials by decay
order (NV spin-Λ: m = 1, 2 → Γ⁻², Γ⁻³; group-IV orbital-Λ: m = 0 → Γ⁻¹).

### Lemmas to prove (new work)

- **II.1 General full-rank damping.** Thm 3 currently covers ΓI only.
  Redo the Neumann expansion for ΓD + A₀ with D ≻ 0 (or merely invertible
  with numerical range in the right half plane): moments become words in
  D⁻¹A₀ acting on D⁻¹c; convergence for Γ > ‖D⁻¹A₀‖; identical hierarchy
  statement. *Proof route:* factor A_Γ = ΓD(I + Γ⁻¹D⁻¹A₀). Difficulty: low.
- **II.2 Non-normal and Jordan corrections.** When D⁻¹A₀ is
  non-diagonalizable, the remainder constants involve the Jordan condition
  number; prove the hierarchy exponent is unchanged and only the constant
  degrades (bound via ‖(D⁻¹A₀)ⁿ‖ ≤ C n^{k−1} ρⁿ with k = largest Jordan
  block). State explicitly when a Jordan block *at the moment index* can
  produce log-free polynomial enhancement of the coefficient. Difficulty:
  medium.
- **II.3 Joint limits.** Uniformity of the expansion in z on compact subsets
  of the regular domain, and the two-parameter regime z ~ Γ (detuning
  scaling with dissipation): show the hierarchy holds pointwise and
  uniformly, and characterize the crossover scale. *Proof route:* resolvent
  bounds from I.1 rationality. Difficulty: medium.
- **II.4 Symmetry-forced coefficient vanishing.** If a symmetry (Thm 5
  projectors commuting with D and A_diag) kills M_m predicted by the
  path-length rule, the order increases to the next symmetry-allowed moment;
  give the group-theoretic selection rule (connect to v6_2 §11).
  Difficulty: medium.
- **II.5 Multi-sector dominance.** For δχ = Σ_S δχ_S with different orders
  m_S: the observable decay order is min(m_S + 1) *unless* leading
  coefficients cancel; prove genericity of non-cancellation and state the
  measure-zero cancellation condition. Difficulty: low–medium.

### Numerical audit (extend existing + new `src/prx_theorem2_hierarchy.py`)

- Reuse `moment_order_common_pipeline.py` (NV −2/−3, SiV/SnV −1 already
  confirmed) as the physical-model verification; add an atomic Λ reference
  system (analytic standard EIT, m = 0).
- New checks: anisotropic D (random positive-definite, seeded) slope vs
  predicted m; a constructed Jordan-block model confirming II.2; a
  two-sector model confirming II.5 dominance and an engineered cancellation.
- All physical models, no random-matrix-only claims in the paper (roadmap
  requirement); random models are used only as internal falsification.

---

## 3. Theorem III — Singular-Dissipation Protected-Channel Theorem

**Target statement.** If ker D ≠ {0} with projector P, the response retains
an O(1) part iff (a) input projects onto ker D, (b) output projects onto
ker D, and (c) the in-kernel effective operator PA₀P (after Schur
elimination of the damped complement) connects them; otherwise the response
is O(Γ⁻¹) or higher, with a complete order classification.

### Lemmas to prove (new work)

- **III.1 Order classification.** Thm 6 gives the O(1) criterion and the
  O(Γ⁻¹) remainder. Extend to the full hierarchy: when the O(1) scalar
  vanishes, iterate the kernel decomposition on the Schur-corrected
  effective operator to classify O(Γ⁻¹), O(Γ⁻²), … channels (nested
  protection). *Proof route:* recursive application of the Thm 6 block
  elimination. Difficulty: medium.
- **III.2 Leakage bound.** Quantify how Hamiltonian coupling PA₀Q "leaks"
  the protected subspace: the O(Γ⁻¹) correction constant is
  ‖PA₀Q D_QQ^{−1/2}‖²-controlled; a no-leak condition ([A₀, P] = 0 on the
  relevant vectors) makes the protected value exact to all orders.
  Difficulty: low–medium.
- **III.3 PA₀P singular case.** Thm 6 assumes PA₀P invertible; treat the
  degenerate case (resonances inside the protected subspace) via the
  pseudo-inverse on the reachable subspace — needed for clock-transition
  examples. Difficulty: medium.
- **III.4 Dictionary lemmas.** Prove the formal correspondence claimed in
  the roadmap: decoherence-free subspace = ker D for collective dephasing;
  dark polariton = protected channel of the cavity-Λ Liouvillian;
  subradiant manifold = ker of the collective radiative dissipator.
  Each as a short lemma identifying (D, P, A₀) for that literature model.
  Difficulty: medium (mostly careful model translation).

### Physical example (mandatory for PRX — pick primary + backup)

- **Primary (recommended): nuclear-spin-assisted protected Raman channel in
  NV.** The Gate 2 machinery already implements ¹⁴N hyperfine levels with
  mI conserved by *every* jump operator — i.e. the nuclear coherence sector
  is exactly ker D of the implemented dissipator. Build the Λ scheme on
  nuclear-spin-split ground states and show the O(1)/O(Γ⁻¹) plateau of the
  window as orbital dissipation Γ(T) is increased (temperature sweep).
  Advantage: continuous with the PRL paper, all parameters sourced, code
  exists (`gate2_candidate_full_vs_reduced.py` extension).
- **Backup: two-emitter subradiant Λ.** Two NV-like defects with a common
  radiative bath; the antisymmetric collective state spans ker D. Small
  dimension (≤ 12 levels), directly buildable on `liouvillian_core.py`.
  Use if the nuclear-channel plateau turns out to be spoiled by the
  nuclear-spin-flipping terms once added realistically (decision point
  below).
- Rare-earth clock transitions remain a citation-level discussion, not the
  computed example (parameter provenance too weak for a first-principles
  plateau plot).

### Numerical audit (new script `src/prx_theorem3_protected.py`)

- Plateau test: |δχ_S| vs Γ over 4+ decades for (i) the protected model
  (expect O(1) or the classified O(Γ⁻ᵏ)), (ii) the same model with the
  protection explicitly broken (expect Thm II hierarchy), (iii) leakage
  scan validating the III.2 bound.

---

## 4. Cross-material classification chapter (Classes A–D)

Assemble the classification (roadmap §8) as a *computed* table, one row per
(material, channel), using the theorem-level quantities:

| Input | Tool (exists) | Class test |
|---|---|---|
| NV spin-Λ, B⊥ = 0 | `nv_reduced_kernel.py` + I-certificate | A (exact structural no-go) |
| NV spin-Λ, B⊥ ≠ 0, T sweep | Gate 4 outputs | B (asymptotic no-go, order 2/3) |
| SiV/SnV orbital-Λ | `group_iv_model.py`, AIC gate | C (conditional go, order 1) |
| hBN candidates | `hbn_nogo_scan.py` | C conditional (parameter-gap statement, roadmap §10) |
| nuclear-assisted NV / subradiant pair | Theorem III example | D (protected go) |
| atomic Λ (reference) | new analytic module | go reference, order 1 with small γ_g |

New script `src/prx_classification_table.py` emits
`results/tables/prx_class_table.csv` — this table *is* PRX acceptance
criterion 4 and the paper's central figure input.

---

## 5. Manuscript and repository layout

- New TeX: `Theorem and proofs/prx/theorem_I.tex`, `theorem_II.tex`,
  `theorem_III.tex`, `classification.tex` — each theorem file contains
  statement, lemmas, proofs, and a "numerical audit" subsection referencing
  the scripts; keep the elementary, self-contained proof style of
  `eit_nogo_proofs.tex` and import (do not re-prove) Thms 2A/2B/3–7 by
  reference.
- Scope statement fixed once, used verbatim everywhere:
  *"finite-dimensional, Markovian, weak-probe open quantum systems"* —
  with the roadmap §9 exclusion list as a single table in the discussion.
- PRX figure plan: Fig 1 concept (transfer function vs dark state, Kalman
  picture); Fig 2 hierarchy (slopes −1/−2/−3, three materials, one
  pipeline — data exists); Fig 3 protected-channel plateau (Theorem III
  example); Fig 4 class map A–D across materials.
- Tests: one `tests/test_prx_theorems.py` with a numerical check per lemma
  tag (I.1–III.3), same style as the Gate tests, fixed seed.

---

## 6. Execution order, decision points, milestones

```
M1  Theorem II first (lowest risk, anchors exist):
    II.1 -> II.5 proofs + prx_theorem2_hierarchy.py            (~1 week)
M2  Theorem I: I.1, I.2, I.4 proofs + certificate script;
    then I.3 control-theory correspondence; I.5 examples       (~1-2 weeks)
M3  Theorem III: III.1-III.3 proofs; DICTIONARY III.4;
    primary physical example computed                          (~2 weeks)
      Decision point D-III: if the nuclear-spin channel shows no usable
      plateau under realistic nuclear relaxation, switch to the
      two-emitter subradiant example (backup path, same scripts).
M4  Classification table + class map figure                    (~3 days)
M5  TeX assembly, audit subsections, tests green, plan update  (~1 week)
```

- Theorem numbering discipline: PRX Theorems I–III are *new statements*;
  existing Thms 1A–7 keep their numbers and are cited as lemma inputs, so
  the PRL supplement and the PRX manuscript never conflict.
- Every lemma lands only with its numerical audit passing; a lemma whose
  audit fails is demoted to a conjecture remark, not silently dropped.
- PRL referee feedback (submission in progress) is triaged into this plan
  at M4 (roadmap §15: reviews feed the PRX theory).

## 7. Acceptance checklist (maps to roadmap §14)

- [ ] I, II, III each: statement + complete proofs + passing numerical audit
- [ ] Basis independence proven (I.4) and stated for II/III via projector form
- [ ] One computed realistic protected channel (III example) with plateau figure
- [ ] `prx_class_table.csv` covers NV, SiV, SnV, hBN(conditional), atomic Λ,
      protected example under one criterion set
- [ ] Demarcation results included: dark-state ≠ no-go (I.5), Kalman
      correspondence (I.3), DFS dictionary (III.4)
- [ ] Scope limited to finite-dimensional Markovian weak-probe systems in
      every general claim
