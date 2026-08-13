# PROOFS.md — where the arguments are

The proofs are kept in their original typeset form rather than paraphrased here. This file maps statement to argument and records which source is authoritative.

## Sources and authority

| File | Contents |
|---|---|
| `proofs/eit_nogo_proofs.tex` | Compatibility entry point. It loads the canonical v7 source below. |
| `proofs/eit_nogo_proofs_v7.tex` | **Authoritative no-go proof source.** Laurent-first asymptotics, sector-graph lower-bound theorem, corrected stationary-pure-state proof, Krylov exact-zero criterion, observable-order inheritance, and scoped singular-damping results. |
| `../archive/historical_notes/theory/eit_nogo_proofs_v6_2_archive.tex` | Verbatim pre-v7 `eit_nogo_proofs.tex`, retained for provenance and comparison. |
| `proofs/EIT_no_go_go_theory_v6_2_English.tex` (+ `.pdf`) | Historical full binary NV no-go/go theory document, v6.2. |
| `proofs/eit_nogo_lecture.tex` | Lecture-form exposition predating the v7 proof revision. |
| `proofs/T1_sector_cut_axiomatization.tex` | Independent sector-cut axiomatization package. |
| `proofs/prl_eit_equivalence_conditions.md` | Response-theoretic EIT-equivalence package. |

## v7 revision status

The v7 source changes the proof package in the following publication-critical ways:

- the large-dissipation order theorem is derived from the **Laurent expansion at infinity** of the finite-dimensional rational resolvent; the Neumann series is retained as a quantitative convergence and remainder lemma;
- the former sector-graph remark is promoted to the **Sector-graph lower-bound theorem**: graph distance gives a lower bound, while equality additionally requires a nonzero summed shortest-path amplitude;
- the stationary pure Lindblad-state proof is rewritten without the invalid shortcut `ρ_D L = λ* ρ_D`;
- the exact transfer-zero theorem is formulated using the reachable Krylov dimension `r` and no longer needs the Neumann norm condition for the exact-zero implication;
- observable-order inheritance and the normalization caveat are explicit;
- the singular-damping theorem is explicitly restricted to Hermitian positive-semidefinite `D`; a general nonnormal Liouvillian extension is not claimed;
- non-semisimple zero modes are handled through the Laurent principal part.

The August 1 audit in `../evidence/benchmarks/AUDIT_theorem_and_code_2026-08-01.md` was performed against the pre-v7 source. It remains provenance for that version. Where v6.2 and v7 differ, v7 controls. PR #23 records successful local TeX compilation, static theorem-contract checks, and a three-sector numerical check performed during the v7 revision.

## Remaining independent gaps

The T1 package is separate from the v7 cleanup. T1 Theorem 2.1 is not proved in this repository, Theorem 2.2's counterexample numbers are not generated here, and the genericity language in Theorem 2.3(ii) remains scoped in `LIMITATIONS.md`.

## Numerical corroboration

| Statement | Check |
|---|---|
| Laurent first-nonzero-moment rule | Existing moment-classification gates plus `calculations/tests/nogo/test_theorem_tex_revision.py`. |
| Sector-graph lower bound | Three-sector v7 witness plus `calculations/numerics/No-go theorem/src/verify_nv_3E_graph_distance_PRL.py`. |
| Sector cut / Schur response | Operational-cut tests and hidden-class-transition certificates. |
| Exact transfer zero / protected channels | Phase A tests and certificates. |
| Observable-order inheritance | Phase M and observable-inheritance campaign results. |
| Matched-readout floor | `results/certificates/gate_F5C_matched_floor_accretive.json` and `gate_F5B_full_gksl.json`. |
| NV rate map `γ_oc = Γ_XY/4` | `calculations/tests/nogo/test_phonon_rate_variants.py`. |

Compile the canonical package from `theory/proofs/` with `pdflatex eit_nogo_proofs.tex`.
