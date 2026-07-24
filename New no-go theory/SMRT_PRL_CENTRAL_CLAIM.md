# SMRT PRL central claim — target statement and closure plan

**Status: DRAFT / planning document.** This is not a finished theorem. It
fixes what the PRL paper must prove, separates that from what is already
proven, and specifies the blind-prediction protocol that will be run on a
second, independent architecture before the claim is considered closed.

## 0. Journal split

| Paper | Journal | Character |
|---|---|---|
| No-go theorem (Theorem I–III, operational sector cut, taxonomy) | PRA | Exhaustive structural classification; full rigor, full scope |
| **SMRT** (this document) | **PRL** | One law that closes EIT, ATS, and CPT as sector-response regimes of a single classification, demonstrated by blind prediction on an independent second architecture |

PRL requires a single, broad-impact claim inside a 4-page Letter (see prior
discussion in this session on PRL submission criteria). The PRA paper carries
the proofs and the full taxonomy; the PRL paper carries **one closed law**
plus the cross-phenomenon unification plus the blind-prediction result.

## 1. What is already proven (building blocks, not the PRL claim itself)

These results exist in this repository and are reused, not re-derived:

1. **Three-tier trichotomy** (`Theorem and proofs/three_theorems_proofs.tex`,
   `src/three_theorem_verification.py`): every sector-mediated response falls
   into exactly one of Class I (ν = ∞, exact structural no-go), Class II
   (0 < ν < ∞, asymptotic no-go), Class III (ν = 0, protected channel). This
   is basis-independent (Kalman) and the trichotomy is exclusive and
   exhaustive.

2. **Path-resolved order theorem** (`SMRT_scaling_assumption_revision.md`
   §6): for a joint native/intervention scaling path κ = κ₀Γ^q, the
   suppression order is a weighted-Newton-degree difference,
   ν_S(q; κ₀) = d_{q,κ₀}(Q) − d_{q,κ₀}(N), hence finitely decidable and
   piecewise linear in q. This is the mathematical spine the PRL claim sits
   on: **suppression order is not an invariant of the system alone but a
   valuation selected by (system, probe, readout, sector intervention, joint
   scaling path).**

3. **Exact witness of path-dependence** (Phase N, `PhaseN/`): a single fixed
   five-level GKSL diamond exhibits ν_{S34}(q) = 4−q (0≤q≤1), 2+q (1≤q≤2), 4
   (q≥2) — proving the order is genuinely path-selected, not a system
   constant, in a concrete physical model with exact rational arithmetic.

4. **EIT/ATS discriminant already connected to the sector-cut structure**
   (`New no-go theory/README.md` §"Key finding", Phase A/B `src/run_phase_a.py`):
   in the minimal 2×2 Λ model, the sector cut severs the coupling responsible
   for both the EIT notch and the ATS splitting; R_S alone carries the
   coherent structure, and the Anisimov–Kocharovsky/Giner discriminant
   Ω_c ≷ |γ31−γ21| separates the EIT-dominant and ATS-dominant regions. The
   AIC/AICc statistical classifier for EIT vs. ATS vs. Lorentzian vs. Fano
   already exists (`No-go theorem/src/eit_ats_classifier.py`).

5. **Physical hidden-transition realization** (Phase H): a fixed scalar
   full-cut response can stay Class III while its sector-attributable part
   transitions Class III → Class II, invisible to single-readout
   classification.

None of these, individually, is the PRL claim. They are the parts that must
be assembled into one law.

## 2. The unification target: EIT, ATS, CPT as regimes of one sector-response law

EIT, ATS, and CPT are conventionally treated as three distinct phenomena
(spectral transparency window, resonance splitting, and steady-state dark-state
population trapping, respectively), each with its own discriminant literature
(Anisimov–Kocharovsky/Giner for EIT vs. ATS; separate treatments for CPT as a
DC/zero-detuning limit). The PRL claim is that **all three are the same
sector-mediated response function R_S(Γ, κ; z), read out at different points
of the (q, z) plane fixed by the path-resolved order theorem**:

- **EIT** ↔ Class III / protected regime of χ_full at finite probe detuning
  z ≠ 0, where the sector-cut interference produces the notch (existing Phase
  A/B result, item 4 above).
- **ATS** ↔ the same response once the discriminant Ω_c ≷ |γ31 − γ21| crosses,
  i.e., the Class III protected structure gives way to resolved-splitting
  behavior in the same R_S(z) at large control coupling — a different region
  of the same discriminant, not a different theory.
- **CPT** ↔ the z → 0 (steady-state) limit of the same sector response: the
  dark-state population trapping is the DC value of the same protected
  channel whose finite-z behavior is the EIT notch. This connects the
  spectral classification (Theorem I–III, path-order ν_S(q)) to a
  steady-state observable via the existing z-dependence machinery
  (`PhaseN/priority_3_frequency/`).

**What this requires that does not yet exist:** an explicit derivation
connecting the z → 0 limit of R_S to the CPT dark-state population, and a
single statement of the discriminant (in terms of Ω_c, |γ31−γ21|, and the
path variable q) that reproduces the EIT/ATS boundary as one section and the
CPT limit as another section of the same law. This is the theorem to be
closed for the PRL submission — not a new numerical architecture, but a
closure argument over existing machinery (items 1–5).

## 3. Target theorem statement (to be proven, not yet final)

> **SMRT unification theorem (target).** For a finite-dimensional Markovian
> weak-probe system with a designated coherent sector S, the sector-mediated
> response R_S(Γ, κ₀Γ^q; z) determines, from a single discriminant built from
> (Ω_c, |γ31 − γ21|, q, z):
> (a) an EIT notch at fixed z ≠ 0 in the Class III / small-q regime;
> (b) an Autler–Townes splitting of the same resonance once the discriminant
>     crosses Ω_c ≷ |γ31 − γ21|;
> (c) coherent population trapping as the z → 0 limit of the same response,
>     with trapped-population fraction determined by the Class III plateau
>     coefficient F₀ of Theorem III.
> All three are sections of one path-resolved, finitely-decidable
> classification; none requires a separate theory.

This is a *target* — the derivation connecting (c) to (a)/(b) through the
existing z-dependence and path-order machinery is the open item.

## 4. Blind-prediction protocol on a second architecture

To make the PRL claim falsifiable rather than a post-hoc fit to the Phase N
witness, the unification theorem above must be **pre-registered as a
prediction on a second, independent architecture** before that architecture
is numerically run.

### 4.1 Requirements for the second architecture

- Must be a genuinely different level scheme from the Phase N five-level
  diamond (different topology, not a relabeling) — e.g., a different Λ/V/ladder
  combination, or a different sector-cut structure.
- Must support all three observables: a finite-z spectral response (for
  EIT/ATS), a z → 0 steady state (for CPT), and a sector cut admissible under
  the existing GKSL intervention definition.
- Must be specified completely (Hamiltonian, dissipators, source, readout,
  sector S, scaling path family) **before** any numerical run against it.

### 4.2 Pre-registration procedure

1. Freeze the architecture specification (Hamiltonian, jump operators,
   source/readout, candidate sector S) in a committed file, timestamped by
   git commit, **before** running any code against it.
2. From the target theorem (§3) and the existing path-order machinery, derive
   and commit numerical predictions in the same commit or an immediately
   following one, still before execution:
   - the predicted ν_S(q) fan (breakpoints and slopes),
   - the predicted EIT/ATS discriminant crossing point,
   - the predicted CPT trapped-population value at z = 0.
3. Only after the prediction commit lands, run the numerical pipeline against
   the second architecture and record the actual gate outputs.
4. Report both the predicted and observed values with the two commit hashes
   (prediction-before, result-after) as the audit trail. Git commit order is
   the blindness guarantee — no prediction may be edited after the result
   commit exists.

### 4.3 Acceptance gate for the PRL claim

The PRL unification claim is considered closed only if:
- the predicted ν_S(q) fan breakpoints match the observed ones exactly (as in
  Phase N, Priority 2's exact breakpoint matching), and
- the predicted EIT/ATS discriminant crossing agrees with the observed AIC/AICc
  classification boundary within the existing statistical tolerance
  (`No-go theorem/src/eit_ats_classifier.py` gate criteria), and
- the predicted CPT z → 0 value matches the observed steady-state trapped
  population within the numerical precision already demonstrated in Phase N
  (rational/high-precision arithmetic where available, else the existing
  smoke/production tolerance).

A partial match (e.g., correct fan but wrong discriminant crossing) is a
negative result to be reported honestly, not adjusted after the fact — the
same standard already applied to Gate M4/P3's corrected predictions elsewhere
in this repository.

## 5. Remaining work before the PRL draft can be written

1. Derive the z → 0 CPT connection to the Class III plateau coefficient F₀
   (currently missing — see §2).
2. Write the single discriminant statement combining Ω_c ≷ |γ31−γ21| with the
   path variable q and frequency z (currently the EIT/ATS discriminant and
   the path-order theorem exist as separate results).
3. Specify and freeze the second architecture (§4.1–4.2).
4. Run the pre-registration/blind-prediction procedure and record the
   outcome, positive or negative.
5. Only then compress the result into the PRL 4-page format with the
   mandatory 100-word justification (broad impact: one closed law replacing
   three separately-discriminated phenomena, demonstrated by a genuine blind
   prediction).

## 6. Related files

- `SMRT_scaling_assumption_revision.md` — path-resolved order theorem (item 2
  above).
- `PhaseN/` — exact witness of path-dependence (item 3), frequency-dependence
  machinery needed for the CPT z → 0 limit (item 2 requirement).
- `PhaseH/` — hidden-transition mechanism relevant to why a single readout can
  miss the sector-level classification.
- `Theorem and proofs/three_theorems_proofs.tex` — Class I/II/III trichotomy
  (item 1).
- `README.md` §"Key finding" — existing EIT/ATS sector-cut connection (item
  4), Anisimov–Kocharovsky/Giner discriminant.
- `../No-go theorem/src/eit_ats_classifier.py` — AIC/AICc statistical
  classifier reused for the acceptance gate in §4.3.
