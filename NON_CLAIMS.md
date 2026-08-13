# Non-claims

Statements this repository must **not** make. Each entry names why, and where
the refutation or the gap is recorded. A claim that appears here may not be
resurrected by rewording; it can only be replaced by new evidence that closes
the specific gap named.

## N1 — "A measurable design rule across engineered-dissipation platforms"

**Rejected in its general form by Gate E (2026-08-02).** Two of three
platforms fail outright. The surviving statement is the narrow one in
`CLAIMS.md` C7, and it is model-level and conditional.

Record: `results/certificates/gates_summary_gateE.json`;
`gates/closed/GENERAL_PRL_GATE_AUDIT.md`.

## N2 — "Optical NV can measure its own suppression exponent"

**It cannot.** The signal reaches the detection floor before the asymptotic
regime: usable window 1.02 decade (single), 0.84 (post-selected shimmed),
0.51 (high density), against 1 decade required. This is the single most
important finding of the 2026-08-02 audit (G-D7 slope budget).

Record: `evidence/failures/` and `results/tables/gate_e_windows.csv`.

## N3 — "Gate B gives a measurement procedure for a superconducting device"

**No.** At realistic `κ` the effective exponent is ≈ 0 at 0.1–50 MHz for both
generic and protected cases; integer asymptotics are readable only 8+ orders
above physical bath damping. Gate B's claim is *structural universality*, not
a measurement procedure.

Record: `results/certificates/gates_summary_gateB.json`.

## N4 — Group-IV (SiV/SnV) curves as quantitative predictions

**Not quantitative.** The group-IV phonon normalization and dipole geometry
are schematic — single rate scale, orthogonal orbital basis legs, `|M0|`
normalization — and were **not** refit. Integer exponents are unaffected;
absolute response levels must not be quoted.

Record:
`calculations/numerics/New no-go theory/GateC_material_independence/SiV_SnV_phonon_AIC_parameters.md`.

## N5 — Gate E's 3-mode PASS as an experimentally demonstrated result

**It is conditional on undemonstrated design assumptions:** 15–45 MHz tunable
engineered loss, coupling-preserving independent loss sweep, −190 dBm-class
amplification and calibration, sign-switching transients, and a measured
feedthrough drift spectrum.

Record: `theory/ASSUMPTIONS.md`.

## N6 — T1 Theorem 2.1 as a result of this repository

**It is not proved here.** It is applied verbatim from a package that is not
present in this repository. Theorem 2.2's counterexample numbers are generated
by no code here either.

Record: `theory/LIMITATIONS.md`.

## N7 — T1 Theorem 2.3(ii) as a genericity statement

**The genericity is not defined.** The result is argued for a single
Lorentzian pole but asserted for "generic `K(z)`" with no definition of the
dense set. Until that is fixed, the defensible form is "any passive response
with a single Lorentzian pole".

Record: `theory/LIMITATIONS.md`; open question O3 in `NEXT_GATES.md`.

## N8 — Manuscript Fig. 4 as a `quick=False` figure

**It is a `quick=True` figure** while the README instructs `quick=False`,
which prints different values in the figure. Do not quote its numbers as the
full-run numbers until it is regenerated.

Record: `manuscript/submission_notes.md`.

## N9 — "150 random draws show no counterexample" as evidence for positivity

**Withdrawn as evidence.** F5B's random scan was the weakest possible support
for a positivity claim, and it came within two implementation bugs (kron
argument order, linear-response sign) of reporting the wrong conclusion. It
was replaced by the one-line identity proof in F5C, which is what `CLAIMS.md`
C3 rests on.

Record: `gates/closed/F5B_findings.md`, `gates/closed/F5C_findings.md`;
`archive/superseded_claims/`.

## N10 — Absence from the Zotero library as evidence of novelty

Not evidence. The novelty boundary this repository asserts is the deliberately
narrow one in `literature/NOVELTY_AUDIT.md`; several theorems are known
results in other clothing and are labeled as such.
