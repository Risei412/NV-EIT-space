# SMRT PRL central claim — target statement and closure plan

**Status: DRAFT / planning document.** This is not a finished theorem. It
fixes what the PRL paper must prove, separates that from what is already
proven, and specifies the blind-prediction protocol that will be run on a
second, independent architecture before the claim is considered closed.

## 0. Journal split

| Paper | Journal | Character |
|---|---|---|
| No-go theorem (Theorem I–III, operational sector cut, taxonomy) | PRA | Exhaustive structural classification; full rigor, full scope |
| **SMRT** (this document) | **PRL** | One label-free law — sector suppression order is a path-selected valuation, not a system invariant — demonstrated by blind prediction on an independent second architecture, with EIT/ATS/CPT recovered as named instances |

PRL requires a single, broad-impact claim inside a 4-page Letter (see prior
discussion in this session on PRL submission criteria). The PRA paper carries
the proofs and the full taxonomy; the PRL paper carries **one closed,
label-free law**, its blind-prediction demonstration, and the recovery of
EIT/ATS/CPT as corollaries that make the claim legible to the AMO/quantum-optics
audience.

**Load-bearing vs. demonstration layer.** The theorem itself is stated and
proved entirely in terms of (sector, source, readout, intervention path) — it
does not name EIT, ATS, or CPT, and its truth does not depend on reproducing
any spectroscopic discriminant from the literature. EIT, ATS, and CPT are a
*corollary layer*: named instances recovered by choosing a specific readout,
frequency, and path from the label-free classification. This separation
matters for two reasons. First, it keeps the PRL claim from being hostage to
matching a specific external discriminant (Anisimov–Kocharovsky/Giner) — if
that identification later needs correction, the theorem is untouched.
Second, it matches the deepest existing statement in this repository
(`SMRT_scaling_assumption_revision.md` §10, §13: "intervention-path
non-intrinsicness" — stated with no reference to EIT/ATS/CPT at all).

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

## 2. The theorem body: label-free path-selected valuation

The PRL claim, stated without reference to any spectroscopic phenomenon:

> Sector suppression order is not an invariant of a system or of a sector
> label. For a specified (sector S, source c, readout p, intervention path
> κ = κ₀Γ^q), the order ν_S(q; κ₀) is a well-defined, finitely-decidable
> valuation — a weighted-Newton-degree difference — and, in general, changes
> value across finitely many path breakpoints. Two different specified paths
> through the same system and the same sector can certify different orders;
> this is not an inconsistency, because the classification was never a
> property of the system alone.

This is exactly `SMRT_scaling_assumption_revision.md` §10's
"intervention-path non-intrinsicness" statement and §6's path-order theorem,
restated as the PRL headline. It requires nothing beyond (system, probe,
readout, sector intervention, joint scaling path) as input — no discriminant
literature, no spectroscopic identification. Phase N (item 3, §1) is its
exact witness: one fixed five-level GKSL diamond, one fixed sector S₃₄, and
the *same* system certifies ν = 3 along one path and ν = 4 along another.

The theorem is falsifiable on its own terms: given a second architecture,
predict the ν_S(q) fan (breakpoints, slopes, exact values at each breakpoint)
before running it, and check for exact agreement (§5).

## 3. Corollary layer: EIT, ATS, CPT as named instances

EIT, ATS, and CPT are not separately hypothesized — they are recovered as
specific (readout, frequency, path) choices inside the label-free
classification of §2, using machinery that already exists in this
repository:

- **EIT** ↔ the Class III / protected regime of the full response at finite
  probe detuning z ≠ 0, where sector-cut interference produces the notch
  (existing Phase A/B result, item 4 of §1: the sector cut severs the
  coupling responsible for both the EIT notch and the ATS splitting).
- **ATS** ↔ the same response once the existing Anisimov–Kocharovsky/Giner
  discriminant Ω_c ≷ |γ31 − γ21| crosses — a different readout region of the
  same R_S(z), not a different theory (item 4).
- **CPT** ↔ the z → 0 (steady-state) limit of the same sector response —
  dark-state population trapping as the DC value of the protected channel
  whose finite-z behavior is the EIT notch.

**Open item for this layer (not for the theorem itself):** the explicit
derivation connecting the z → 0 limit of R_S to the CPT dark-state
population is not yet written down, and the existing Ω_c ≷ |γ31 − γ21|
discriminant has not yet been re-expressed in terms of the path variable q.
Closing this is corollary work — useful for making the PRL claim legible to
an AMO audience, and worth doing, but its status is independent of whether
§2's theorem is true. A gap here would narrow the demonstrated scope of the
corollary, not falsify the theorem.

## 4. What the theorem does *not* claim

Per `SMRT_scaling_assumption_revision.md` §10: given the sector S, the
intervention path, and a tomographically complete generator, the order is
computable — not unidentifiable. The theorem's content is that the order
depends on these specified inputs, not that it is unrecoverable once they are
given. This boundary must be kept in the PRL text; overstating it as "no
observation can ever determine the order" is false and would be corrected in
review.

## 5. Blind-prediction protocol on a second architecture

To make the PRL claim falsifiable rather than a post-hoc fit to the Phase N
witness, the §2 theorem must be **pre-registered as a prediction on a
second, independent architecture** before that architecture is numerically
run. The primary prediction target is the label-free ν_S(q) fan; the
EIT/ATS/CPT identification (§3) is checked as a secondary, corollary-layer
prediction and does not gate the theorem's acceptance.

### 5.1 Requirements for the second architecture

- Must be a genuinely different level scheme from the Phase N five-level
  diamond (different topology, not a relabeling) — e.g., a different Λ/V/ladder
  combination, or a different sector-cut structure.
- Must be specified completely (Hamiltonian, dissipators, source, readout,
  sector S, scaling path family) **before** any numerical run against it.
- If the corollary layer (§3) is also to be tested, the architecture must
  additionally support a finite-z spectral response and a z → 0 steady state
  — but this is only required for the secondary check, not for the theorem
  itself.

### 5.2 Pre-registration procedure

1. Freeze the architecture specification (Hamiltonian, jump operators,
   source/readout, candidate sector S) in a committed file, timestamped by
   git commit, **before** running any code against it.
2. From the §2 theorem and the existing path-order machinery, derive and
   commit numerical predictions in the same commit or an immediately
   following one, still before execution:
   - **primary:** the predicted ν_S(q) fan (breakpoints and slopes) for the
     specified sector and path family;
   - **secondary, if applicable:** the predicted EIT/ATS discriminant
     crossing point and the predicted CPT trapped-population value at z = 0.
3. Only after the prediction commit lands, run the numerical pipeline against
   the second architecture and record the actual gate outputs.
4. Report both the predicted and observed values with the two commit hashes
   (prediction-before, result-after) as the audit trail. Git commit order is
   the blindness guarantee — no prediction may be edited after the result
   commit exists.

### 5.3 Acceptance gate

**Theorem (§2), load-bearing:** closed only if the predicted ν_S(q) fan
breakpoints and values match the observed ones exactly (as in Phase N,
Priority 2's exact breakpoint matching).

**Corollary layer (§3), secondary:** reported separately. If tested, closed
only if the predicted EIT/ATS discriminant crossing agrees with the observed
AIC/AICc classification boundary within the existing statistical tolerance
(`No-go theorem/src/eit_ats_classifier.py` gate criteria) and the predicted
CPT z → 0 value matches the observed steady-state trapped population within
the numerical precision already demonstrated in Phase N. A miss here narrows
the demonstrated scope of the corollary (per §3) without affecting whether
the theorem itself is accepted.

Any partial match is a negative result to be reported honestly, not adjusted
after the fact — the same standard already applied to Gate M4/P3's corrected
predictions elsewhere in this repository.

## 6. Remaining work before the PRL draft can be written

1. Specify and freeze the second architecture (§5.1–5.2).
2. Run the pre-registration/blind-prediction procedure for the §2 theorem
   and record the outcome, positive or negative.
3. Derive the z → 0 CPT connection to the Class III plateau coefficient F₀,
   and write the single discriminant statement combining Ω_c ≷ |γ31−γ21|
   with the path variable q and frequency z (§3's open item) — needed only
   for the corollary layer, can proceed in parallel with or after item 2.
4. If the corollary layer is ready in time, run its secondary blind
   prediction alongside item 2.
5. Compress the result into the PRL 4-page format with the mandatory 100-word
   justification (broad impact: one label-free law governing sector
   suppression order, demonstrated by a genuine blind prediction, with EIT,
   ATS, and CPT recovered as named instances).

## 7. Related files

- `SMRT_scaling_assumption_revision.md` §10, §13 — the label-free
  "intervention-path non-intrinsicness" statement the §2 theorem restates
  for PRL.
- `PhaseN/` — exact witness of path-dependence (§1 item 3), frequency-dependence
  machinery needed for the CPT z → 0 limit (§3's open item).
- `PhaseH/` — hidden-transition mechanism relevant to why a single readout can
  miss the sector-level classification.
- `Theorem and proofs/three_theorems_proofs.tex` — Class I/II/III trichotomy
  (§1 item 1).
- `README.md` §"Key finding" — existing EIT/ATS sector-cut connection (§1
  item 4), Anisimov–Kocharovsky/Giner discriminant.
- `../No-go theorem/src/eit_ats_classifier.py` — AIC/AICc statistical
  classifier reused for the corollary-layer acceptance gate in §5.3.
