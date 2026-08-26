# Next gates

The work that is currently permitted, and nothing else. A gate is *binding*
while it sits in `gates/active/`; its plan document is the authority on what
the gate computes, its inputs, outputs and pass criteria. When a gate closes,
its document moves to `gates/closed/` and the outcome is recorded in
`gates/GATE_LEDGER.md`.

## Binding gate plans (in `gates/active/`)

| Plan | Scope |
|---|---|
| `SIMULATION_PLAN.md` | The PRL verification-gate campaign: what each gate computes, its inputs, outputs and pass criteria. |
| `EXECUTION_PLAN.md` | Execution plan for the EIT-definition-equivalence campaign (Phase E/F). |
| `NEXT_NUMERICAL_PLAN.md` | The next numerical steps for the equivalence campaign. |
| `new_nogo_numerical_priorities.md` | Priority ordering for the three-tier classification campaign (Gates 1–5). |
| `phase_d_beyond_old_theory_plan.md` | Phase D: what goes beyond the binary no-go theory. |

## Open questions that gate the next claim

These are the questions whose answers would change `CLAIMS.md`. They are open;
none of them may be assumed answered.

- **O1 — Is there any platform where `ν` is both large enough to matter and
  small enough to measure?** `ν = 4` costs four decades of contrast per decade
  of `Γ`, which is what closes the NV window. The two requirements pull
  opposite ways; whether the tension is fundamental or an artifact of these
  platforms is open.
  O1 is posed as a platform-selection question, but the declared resource path
  `κ = κ₀Γ^q` is a second, unexploited axis: every `ν` this repository reports
  is an ideal-cut value, which in two-scale language is the corner of the fan
  where the exponent is largest. See `theory/NOTE_path_dependence_and_O1.md`.
  This does **not** reopen `NON_CLAIMS.md` N2, and the gate it suggests is a
  candidate only — it is not in `gates/active/` and is not permitted work.
- **O2 — Does the sector-cut classification have an observable signature that
  is not the exponent?** If NV can only ever serve as a witness certifying an
  exact class, what certifies it?
- **O3 — What is the right genericity statement for T1 Theorem 2.3(ii)?**
  Either the claim weakens to "any passive response with a single Lorentzian
  pole", or genericity gets defined (the dense set must be named).
- **O4 — Can the unmatched escape route be characterized?** The identity says
  exact zeros need `x ∈ ker Γ`, reachable only when source ≠ readout, and F5B
  did find exact zeros in Raman-basis coherence readout. What is the general
  condition?

## Housekeeping items that are permitted without a new gate

- Regenerate manuscript Fig. 4 with `quick=False` and re-check the printed
  values (closes `NON_CLAIMS.md` N8).
- Write out the Theorem 1B rearrangement via
  `H_eff := H − (i/2) Σ_μ L_μ†L_μ` (closes the unwritten step recorded in
  `theory/LIMITATIONS.md`).
- Import, or explicitly cite as external, the T1 Theorem 2.1 proof and the
  Theorem 2.2 counterexample generator (closes `NON_CLAIMS.md` N6).
