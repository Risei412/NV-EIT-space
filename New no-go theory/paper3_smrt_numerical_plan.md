# Paper III numerical plan: interference-controlled exponent promotion in a physical Lindblad model (Phase P)

Status: **EXECUTED.** All eight gates (P1–P8) PASS —
see `results/summary.md` (Phase P section) and
`results/gates_summary_phaseP.json`. Code: `src/model_physical.py`,
`src/run_phase_p.py`. One correction to this plan surfaced during
execution: P3's assumption of a codimension-1 curve directly in the
2-real-parameter (η, φ) plane was wrong (2 real controls generically give
an *isolated* zero, not a curve — see the P3 section of `results/summary.md`
for the corrected 3-control construction). Everything else below matches
what was executed.

Source strategy: `SMRT_Paper_III_numerical_strategy.md` (user-provided, not committed
here), which defines Calculations A–H, the acceptance gates, and the figure plan for

> **Interference-Controlled Dissipative Response Exponents in Open Quantum Systems**

This document turns that strategy into an executable numerical campaign
(Phase P, runs P0–P8) on top of the existing infrastructure in this directory:

- `src/core.py` — moments (`moments_invertible_D`), ν estimators
  (`fit_nu_loglog`, `nu_from_moments`), exact polynomial certificate
  (`certificate_deg_nu`), window norm (`max_abs_on_grid`). All reused as-is.
- `src/model_cancel.py`, `src/model_masquerade.py`, `src/run_phase_d.py` —
  the abstract D2 (order promotion, ν 3→4, v23\* = −0.7846…) and D4
  (masquerade, 2.43-decade false ν≈4 plateau, m₃ = 895/2688) results already
  PASS (`results/gates_summary_phaseD.json`). They are the Stage-1 baseline;
  Phase P embeds the same mechanism in a CPTP Lindblad model.

---

## P0. Design decision: the physical model (revised from the strategy §4.2)

### P0.1 Why the strategy's 4-level scheme is modified

The strategy suggests states |1..4⟩ with paths 1→2→4 and 1→3→4 and the cut
S = {J₁₃, J₃₄} (remove the second branch). Writing the weak-probe reduced
system on the optical coherences and expanding the **master response of the
difference** R_S = χ_full − χ_cut in moments m_k = μ_k^full − μ_k^cut shows a
structural problem:

- the k-hop moment of the *difference* contains only paths that use at least
  one cut edge;
- if the cut removes one entire branch, the leading moment of R_S is the
  amplitude of **that branch alone**, m_lead = 𝒜₂ — a single product of
  couplings that cannot be cancelled by tuning (η, φ) against 𝒜₁.

The two-path interference m_lead = 𝒜₁ + 𝒜₂ (strategy §2.1) lives in R_S only
if **both** interfering routes pass through the cut sector. The minimal fix:
cut the two couplings that close both routes into the readout state.

### P0.2 Adopted model: 5-level diamond ("double-Λ ladder")

States: ground |1⟩; probe-excited |2⟩; two intermediates |3⟩, |4⟩;
readout state |5⟩.

Rotating-frame Hamiltonian (ħ = 1):

    H = Σ_{j=2..5} Δ_j |j⟩⟨j|
        + ε (|2⟩⟨1| + h.c.)                        (weak probe, amplitude ε)
        + [ J23 |3⟩⟨2| + J24 |4⟩⟨2|                 (branch splitters, kept)
          + J25 |5⟩⟨2|                              (weak direct route, kept)
          + J35 |5⟩⟨3| + J45 e^{iφ} |5⟩⟨4|          (sector S: route closers)
          + J34 |4⟩⟨3|                              (weak cross-coupling)
          + h.c. ]

CPTP dissipation, all rates ≥ 0 (Lindblad jump operators):

    L_j = sqrt(Γ d_j) |1⟩⟨j| ,  j = 2..5            (radiative decay to ground)
    optional: L_j^φ = sqrt(Γ κ_j) |j⟩⟨j|            (pure dephasing, if needed)

Sector cut (strategy §4.3, revised): **S = {J35, J45 e^{iφ}}**. The cut model
sets J35 = J45 = 0 while preserving D, c, p — physically, switching off the
control field(s) that connect both coherent branches to |5⟩. R_S then measures
exactly "the response added by coherently closing the two routes into the
readout state", and its leading moment is the coherent **sum of both branch
amplitudes** — the tunable interference the paper needs. The weak J25 keeps
χ_cut ≠ 0 so the intervention is a genuine full/cut pair (χ_cut ≡ 0 variant,
J25 = 0, is kept as a degenerate check).

### P0.3 Weak-probe reduced linear system (the classified object)

To O(ε), with ρ⁽⁰⁾ = |1⟩⟨1| (all decay returns to |1⟩), the response closes on
the optical coherences x = (ρ21, ρ31, ρ41, ρ51). Steady state at probe
frequency z:

    A(z) x = c ,   A(z) = Γ D + B(z)
    D  = diag(d2, d3, d4, d5)/2            (invertible ⇔ all d_j > 0)
    B(z) = i diag(Δ_j − z) + i J           (J = coupling adjacency above)
    c  = −i ε e_{ρ21} ,   p = e_{ρ51}      (primary readout below)

This 4×4 system is graph-isomorphic to the validated abstract D2 model
(`model_cancel.build_promotion_model`): node ρ21 ↔ site 0, {ρ31, ρ41} ↔
{1, 2}, ρ51 ↔ site 3; J34 ↔ v12. All Phase-D machinery therefore transfers
verbatim.

Readouts (Calculation G): O1 = Im ρ51 (absorption on the 5-branch; linear in
ε, exactly the SMRT scalar p†A⁻¹c), O1′ = Re ρ51 (dispersion), O2 = ρ55
(emitted intensity Γ d5 ρ55; O(ε²), computed from the full Liouvillian, with
the Phase-M ν→2ν translation as prediction), O3 = a deliberately blind
readout p = e_{ρ31} (expected ν = 2, explained by observability).

### P0.4 Moment structure and the analytic cancellation condition

With D, B as above, m_k = p†(D⁻¹B)^k D⁻¹ c on the full/cut difference:

- m₀ = m₁ = 0 structurally (c on ρ21, p on ρ51, direct J25 term cancels in
  the full−cut difference);
- m₂ ∝ J23·J35 e^{iφ}/d3 + J24·J45/d4 — **both** branch amplitudes, because
  both 2-hop routes end on a cut edge; z-independent, so the cancellation is
  uniform in z (as Theorem II requires);
- cancellation manifold:  η e^{iφ} = −1,  η ≔ (J23 J35 d4)/(J24 J45 d3);
  i.e. φ⋆ = π at η⋆ = 1 — closed form, no root search needed;
- m₃ at the cancellation point: dressed 2-hop paths ∝ (Δ_a − z)/d_a asymmetry
  (a = 3, 4) plus the 3-hop cross-path via J34. Choosing Δ3 = Δ4 = 0 and
  J34 ≠ 0 makes m₃ z-independent and manifestly nonzero at z = 0; the
  z-window check inf_K |m₃(z)| > 0 is an automated gate.

Hence generically ν = 3, at the tuned point ν = 4 — the strategy's ν: 3→4 in
a CPTP model with all four requirements of §4.1 met.

### P0.5 Baseline parameter set

Exact rationals throughout (needed for the symbolic gates):

    d = (1, 13/10, 17/10, 21/10, 9/10)/... → d2..d5 = (1, 13/10, 17/10, 21/10)
    Δ2..Δ5 = (0, 0, 0, 0) primary;  robustness variant (1/5, 0, 0, −1/10)
    J23 = 1, J24 = 1, J35 = 3/5, J45 = tuned (≈ −(J23 J35 d4)/(J24 d3) at φ=0,
          or |J45| fixed and φ tuned to φ⋆), J25 = 1/10, J34 = 1/20
    probe window K: z ∈ [−1/2, 1/2], 21-point grid; primary evaluation z = 0
    ε (full-Liouvillian runs only): 10⁻⁶ with an ε/2 linearity check

Before any production run, an automated pre-check verifies: no eigenvalue of
A(z) approaches 0 on K for Γ ∈ [Γ_min, Γ_max] (pole exclusion, strategy §9
"promoted behavior is not caused by a nearby pole").

---

## P1–P8: run specifications (= strategy Calculations A–H)

### Precision policy (applies to every run)

The reduced system is 4×4, so **exact rational arithmetic is affordable
everywhere**: with rational parameters and Γ = 10^k, R_S(Γ) is an exact
`sympy.Rational`. Policy:

1. dense Γ grids in float64 while the a-posteriori bound
   `eps · cond(A) · ‖x‖ / |R_S|` < 10⁻³ (cond monitored at every point);
2. beyond that (large Γ near cancellation, where |R_S| ~ Γ⁻⁴ underflows the
   float64 *solve accuracy*, not the exponent range), switch to `mpmath`
   (mp.dps = 50) on the dense grid;
3. exact `sympy` rational solves on a sparse grid (every half-decade) as the
   ground truth that certifies both float64 and mpmath branches;
4. all moments and the polynomial certificate: exact symbolic only — a float
   |m₂| < 10⁻¹² is recorded as a cross-check but never used as evidence
   (strategy Calculation B gate).

The full 5-level Liouvillian (25×25 vectorized; 24×24 on the traceless
subspace) is solved in float64 with an mpmath (dps = 50) branch for
Γ > 10⁴ at the tuned point; the probe is applied as an **implicit linear
response** (solve L₀ρ₁ = −i[V_probe, ρ₀] on the zero-trace subspace, ρ₀ from
the null space of L₀) so no finite-ε subtraction ever limits precision;
finite-ε steady states at moderate Γ validate the implicit solve.

### P1 = Calculation A: direct validation of ν 3→4

For the three parameter sets — generic (φ = 0), tuned (φ = φ⋆, η = η⋆),
detuned (η = η⋆(1+δ), δ = 10⁻²) —

1. build full and cut A(z); compute R_S(Γ, z) on Γ ∈ 10^[0, 8] (dense, 33/decade
   ≥ 8 decades; exact-rational spine per the precision policy);
2. window norm ‖R_S‖_K via `max_abs_on_grid`;
3. fit ν two independent ways: global log–log LSQ over the top 3 decades and
   the local effective index ν_eff(Γ) (`fit_nu_loglog`), plus the moment
   prediction (`nu_from_moments`, exact);
4. compensated plots Γ³|R_S| and Γ⁴|R_S| (plateau diagnosis).

Gates (strategy §5A): |ν_fit − ν_moment| < 0.03 per case; tuned vs detuned
separated by more than the fit spread; generic plateau in Γ³|R_S|, tuned
plateau in Γ⁴|R_S|, detuned Γ⁴-plateau → Γ³-plateau crossover visible.

Output: `figP1_promotion.png` (paper Figure 2 panels a–c),
gate block `P1_promotion` in `results/gates_summary_phaseP.json`.

### P2 = Calculation B: symbolic moment cancellation

Sympy, exact rationals, on the reduced 4×4 difference system:

1. derive m₂ = α₁ J24 J45 e^{iφ} + α₂ J23 J35 with α_j explicit functions of
   (d_j); verify α-structure against §P0.4 by expanding p†(D⁻¹B)²D⁻¹c;
2. solve m₂ = 0 in closed form for φ⋆(η) and η⋆; confirm the solved point
   reproduces the P1 tuned parameters exactly;
3. prove m₃(φ⋆, η⋆) ≠ 0 symbolically at z = 0 and show inf_K |m₃(z)| > 0 on
   the sampled window (rational z samples);
4. float64 moment ladder cross-check (`moments_invertible_D`) at 10 random
   rational parameter draws: symbolic and float values agree to 10⁻¹⁰.

Gate: exact m₂ = 0 and exact m₃ ≠ 0 (never a small float). Output: symbolic
audit appended to `results/summary.md`, gate block `P2_symbolic`.

### P3 = Calculation C: cancellation manifold map

Scan (η, φ) ∈ [η⋆/4, 4η⋆] × [0, 2π] (log × linear grid, 81×81, symbolic m₂
evaluated as float — cheap since m₂ is a closed-form expression):

1. map log₁₀|m₂(η, φ)|; overlay the analytic curve η e^{iφ} = −1 … i.e.
   φ = π, η = η⋆ point and, with the robustness detuning variant Δ ≠ 0, the
   full codimension-1 curve φ⋆(η);
2. verify m₃ ≠ 0 along the sampled curve;
3. ν_eff(η, φ) from a 2-decade Γ fit at 9 points on/off the manifold
   (on-manifold → ≈ 4, off → ≈ 3).

Gate: the zero set is a codimension-1 curve (not an isolated accident);
on-manifold ν_eff > 3.7, off-manifold |ν_eff − 3| < 0.1.
Output: `figP3_manifold.png` (paper Figure 3, panels a–b).

### P4 = Calculation D: universal crossover collapse

Detuning δ ∈ {±10⁻¹, ±10⁻², ±10⁻³, ±10⁻⁴, ±10⁻⁵} applied to η (and, as a
variant, to φ):

1. symbolic a, b: R_S = a δ Γ⁻³ + b Γ⁻⁴ + O(δΓ⁻⁴, Γ⁻⁵) with a = ∂m₂/∂δ|₀ and
   b = −m₃ (signs per Theorem II expansion) — predicted, not fitted;
2. plot Γ⁴ R_S vs δΓ for all δ: collapse onto a(δΓ) + b;
3. extract Γ×(δ) two ways: ν_eff crossing 3.5, and equality |aδ|Γ⁻³ = |b|Γ⁻⁴
   (predicts Γ× = |b/(aδ)|);
4. log–log fit of Γ× vs |δ|.

Gates (strategy §5D): |d log Γ×/d log|δ| + 1| < 0.1; relative collapse spread
< 5% over the common scaling region; fitted Γ× within a factor 2 of |b/(aδ)|.
Output: `figP4_collapse.png` (paper Figure 3, panels c–d).

### P5 = Calculation E: finite-window misclassification map

1. ν_eff(Γ) curves for each δ of P4 (exact-rational spine so the far tail,
   Γ up to 10¹², is trustworthy);
2. sliding-window fits: window width 2 decades, centers sweeping
   Γ ∈ 10^[0.5, 11]; heat map ν_fit(δ, window center);
3. false-plateau width: decades with |ν_eff − 4| < 0.05 at each δ; compare
   with the D4 abstract reference (2.43 decades at δ = 10⁻⁵);
4. side-by-side: finite-window ν_fit vs the exact certificate ν_cert (P6).

Gate: at least one δ with a false ν≈4 plateau > 2 decades while ν_cert = 3.
Output: `figP5_masquerade.png` (paper Figure 4), gate block `P5_masquerade`.

### P6 = Calculation F: exact polynomial certificate

`certificate_deg_nu` on the reduced difference system with rational
parameters (numerator of R_S = N_full Q_cut − N_cut Q_full over Q_full Q_cut,
assembled symbolically in Γ at z = 0 and at two more rational z ∈ K):

| Case | ν_cert | ν_moment | tail fit | 2-decade window fit |
|---|---:|---:|---:|---:|
| Generic | 3 | 3 | ≈3 | ≈3 |
| Tuned | 4 | 4 | ≈4 | ≈4 |
| Detuned 10⁻³ | 3 | 3 | ≈3 | can read ≈4 |

Gate: table reproduced exactly in the first two columns; numerator not
identically zero at the tuned point (Class II, not Class I — strategy §9).

### P7 = Calculation G: observable robustness

For O1, O1′, O3 (linear readouts): repeat P1's three cases by changing p
only. For O2 = ρ55 (second order): full-Liouvillian implicit response at
O(ε²) on a sparse Γ grid.

Expected/gated: O1 and O1′ both show ν 3→4 (coefficients differ); O3 shows
ν = 2 with the moment-level explanation (p†D⁻¹B D⁻¹c ≠ 0 — an observability
statement, reported, not hidden); O2 consistent with the ν→2ν translation
(slope 6↔8 at generic/tuned; Phase M provides the theory reference).

### P8 = Calculation H: reduced vs full Liouvillian

1. same Γ grid, tuned + generic cases, z = 0 and one z ≠ 0;
2. R_S^reduced (4×4) vs R_S^full (24×24 implicit linear response, ε-free);
3. relative deviation plotted across the sweep; a finite-ε (10⁻⁶, 5·10⁻⁷)
   steady-state check at Γ ∈ {10¹, 10³}.

Gate: |R_S^red − R_S^full|/|R_S^full| < 10⁻² over the top 3 decades (< 10⁻³
expected since the reduced system is exact at O(ε) for this model; any excess
must be explained, e.g. dephasing variants coupling population sectors).

---

## Deliverables and file layout

New code (mirroring the Phase-D layout):

    src/model_physical.py   5-level model: reduced (D, B(z), c, p) builder for
                            full/cut, Liouvillian builder (vectorized, CPTP
                            check: Kossakowski matrix ⪰ 0, trace preservation),
                            implicit linear-response solver, mpmath/rational
                            solve branches
    src/moments_p.py        symbolic m_k, cancellation manifold, a & b, certificate
    src/run_phase_p.py      runs P1–P8, writes results/gates_summary_phaseP.json
                            and figures figP1..figP5 (+ supplementary figP6..P8)
    src/report.py           extended to append the Phase-P section to summary.md

Automated gates: one JSON block per run, `overall_pass` aggregating the gate
conditions quoted above — same convention as Phases A/B/M/D, so `report.py`
needs only a new section loop.

Figure → paper mapping: figP1 → Fig. 2; figP3+figP4 → Fig. 3; figP5 → Fig. 4;
figP6–P8 + P2 audit → Supplement. Fig. 1 (schematic) is a drawing task, not a
computation.

## Execution order and effort

| Stage (strategy §11) | Runs | Status / effort |
|---|---|---|
| 1. abstract mechanism | D2, D4 | **done** (Phase D, all gates PASS) |
| 2. physical model | P0 pre-checks, model_physical.py | ~1 session |
| 3. main result | P1, P2, P6 | ~1 session |
| 4. universal scaling | P3, P4 | ~1 session |
| 5. finite-window failure | P5 | short |
| 6. robustness + figures | P7, P8, report | ~1 session |

All matrices are ≤ 25×25; total runtime is minutes. The only nontrivial cost
is symbolic (4×4 adjugate/determinant in Γ — trivial; the 24×24 Liouvillian is
**never** treated symbolically; the certificate is defined on the reduced
system, which is the classified object, and P8 ties it to the full dynamics).

## Risks and fallbacks

1. **m₃ accidentally small at the tuned point** (weak J34): raise J34 or turn
   on the Δ3 ≠ Δ4 variant; gate P2.3 catches this before production runs.
2. **z-window zero of m₃(z)**: shrink K or use the Δ = 0, J34-only variant
   (m₃ z-independent); gated in P2.
3. **Float64 breakdown near cancellation at large Γ**: by design — the
   rational/mpmath spine is the evidence; float64 is only a fast first pass.
   (This breakdown is itself reportable: it is the numerical face of the
   masquerade warning.)
4. **Pole intrusion into K** at some (η, φ) scan corner: pre-check spectral
   distance of A(z) to singularity; mask such points in P3 maps.
5. **ρ55 readout not showing clean 2ν doubling** at small Γ (population
   sector mixing): restrict the fit window to the top decades; Phase M Gate
   M5 documents when the translation legitimately fails.
6. **Cut interpretation challenged** (S = two couplings rather than one
   branch): the χ_cut ≡ 0 single-knob variant (J25 = 0, S closes the only
   access to |5⟩) is retained as a supplementary check, and the strategy's
   branch-cut variant is reported with its (non-tunable) m₂ = 𝒜₂ as the
   explicit counter-illustration of why the sector must gate both routes.

## Acceptance checklist (strategy §9, automated where possible)

- [ ] m₂ = 0 exact (symbolic), m₃ ≠ 0 exact, at the tuned point (P2)
- [ ] numerator ≢ 0 at the tuned point — Class II, not Class I (P6)
- [ ] detuned certificate ν = 3; tail fit agrees with certificate (P5/P6)
- [ ] D invertible on the response-relevant block; CPTP verified; rates ≥ 0 (P0)
- [ ] weak-probe linearity verified (ε vs ε/2) (P8)
- [ ] no nearby pole in K across the sweep (P0 pre-check)
- [ ] condition numbers logged across every Γ sweep; window-stable fits (P1–P5)
- [ ] two independent ν estimators per claim; symbolic ↔ float moments agree (P1/P2)
- [ ] reduced ↔ full Liouvillian agree < 1% (P8)
