# Changelog

## 2026-08-13 — canonical layout restructure

Reorganized the repository from campaign-shaped directories into the canonical
theory-repository layout. No physics changed; no numerical result changed.

**Moved**

| From | To |
|---|---|
| `No-go theorem/`, `New no-go theory/`, `EIT definition equivalence/` (code) | `calculations/numerics/<same name>/` |
| `Writing Paper/pra/src`, `Writing Paper/prl_figures/src` | `calculations/numerics/manuscript_figures/{pra,prl}/src` |
| `master-equation/` | `calculations/analytic/master-equation/` |
| every campaign's `tests/` | `calculations/tests/{nogo,gate_a…gate_d,pra,prl_figures}/` |
| every campaign's `results/` | flat `results/{tables,figures,raw,certificates}/` |
| `No-go theorem/Theorem and proofs/`, `EIT definition equivalence/tex`, `prl_eit_equivalence_conditions.md` | `theory/proofs/` |
| `SIMULATION_PLAN.md`, `EXECUTION_PLAN.md`, `NEXT_NUMERICAL_PLAN.md`, `new_nogo_numerical_priorities.md`, `phase_d_beyond_old_theory_plan.md` | `gates/active/` |
| gate findings and reports (`E1_E3`, `E4`, `F4`, `F5A–C`, `F_preview_audit`, `gate_1_5_report`, `common_liouvillian_gate_report`, `hbn_nogo_EIT_report`, `GENERAL_PRL_GATE_AUDIT`, `P0_1`, `P0_2`) | `gates/closed/` |
| `2g2e_package/` (intact) | `evidence/witnesses/` |
| `AUDIT_theorem_and_code_2026-08-01.md` | `evidence/benchmarks/` |
| `PRL_CLAIM_GAP_AUDIT.md` | `evidence/failures/` |
| `EIT_general_theory_literature_2026-07-13/`, `hBN papers/` | `literature/` |
| `Writing Paper/pra/manuscript/nv_eit_pra.{tex,pdf}` | `manuscript/main.{tex,pdf}` |
| `Writing Paper/{README,split strategy,drafts}` | `manuscript/{submission_notes.md, …, drafts/}` |

**Added** — `summary.json`, `REPOSITORY_STATE.md`, `NEXT_GATES.md`,
`CLAIMS.md`, `NON_CLAIMS.md`; `theory/{THEORY,DEFINITIONS,THEOREMS,PROOFS,ASSUMPTIONS,LIMITATIONS}.md`;
`calculations/README.md`; `results/summary.json`; `gates/GATE_LEDGER.md`;
`literature/{NOVELTY_AUDIT,PRIOR_ART_BOUNDARY}.md`;
`provenance/{MANIFEST.json,PROVENANCE.json,SHA256SUMS.txt,CHANGELOG.md}`;
directory READMEs under `evidence/`, `archive/`, `calculations/`.

**Code changes** — output and import path constants only:

- Every runner now writes into the repository-level `results/` tree, routed by
  artifact kind. Gate JSONs go to `results/certificates/`.
- Test modules resolve their sources from `calculations/tests/<x>/` to
  `calculations/numerics/…`.
- Three scripts that wrote to directories that did not exist in the repository
  (`outputs/`, `analysis/`, and an absolute `/home/claude/` path) now write
  into the canonical `results/` tree.

**Verification** — `pytest calculations/tests`: 70 passed. The PRL figure
regression tests regenerate `results/figures/fig{1..4}_*.{png,pdf}`
identically to the committed artifacts.

## Earlier

Dated audit history is in `PROVENANCE.json` (`audits`) and, in narrative form,
in `../archive/historical_notes/README.md`.

## Note on path references in dated documents

Navigational documents (`README.md` files, `seed.md`) had their path references
rewritten to the new layout. **Dated records were not rewritten** — audits,
gate plans and findings documents (`evidence/benchmarks/AUDIT_*`,
`evidence/failures/PRL_CLAIM_GAP_AUDIT.md`, `gates/active/*`,
`gates/closed/*`) still cite the paths that existed when they were written.
Editing them would falsify the record; the table above is the old → new map.
