# Gate ledger

Every gate this repository has run, and how it ended. A gate is a
pre-registered computation with a stated pass criterion: the criterion is
fixed before the run, and the verdict is whatever the run returns.

Certificates are in `../results/certificates/`. Plans for gates still binding
are in `active/`; documents for gates that have closed are in `closed/`.

## Binary NV no-go campaign (Gates 1–5)

| Gate | Question | Verdict | Certificate |
|---|---|---|---|
| 1 | Candidate model selection stability (AIC bootstrap) | PASS | `gate1_aic_bootstrap.json` |
| 2 | Full vs. reduced Liouvillian agreement at the candidate | PASS | `gate2_candidate_comparison.json` |
| 3 | SNR map — is the predicted contrast detectable? | PASS | `gate3_summary.json` |
| 4 | Threshold uncertainty bands | PASS | `gate4_sensitivity.json` |
| 5 | Ensemble averaging survival | PASS | `gate5_summary.json` |
| 1.5 | Interim consolidation | closed | `closed/gate_1_5_report.md` |
| 21–33 | Detuning protocols, convention and contour audits | closed | `gate21_detuning_protocols.json`, `gates_22_23_26.json`, `gates_24_28_31_33.json` |
| — | Common-Liouvillian gate | closed | `closed/common_liouvillian_gate_report.md` |
| — | hBN no-go scan | closed | `closed/hbn_nogo_EIT_report.md` |

## Three-tier classification campaign

| Gate | Question | Verdict | Certificate |
|---|---|---|---|
| Phase A | Class I/II/III unit tests; EIT/ATS/background separation | PASS | `gates_summary_phaseA.json` |
| Phase B | Classifier consistency; hidden class transition; scaling collapse | PASS | `gates_summary_phaseB.json` |
| Phase D | Beyond the binary theory: oblique protection, cancellation promotion, masquerade | closed | `gates_summary_phaseD.json` |
| Phase M | `ν → 2ν` metrological (QFI) translation, with rank-deficiency negative control | PASS | `gates_summary_phaseM.json` |
| Phase P | Promotion / manifold / collapse / masquerade | closed | `gates_summary_phaseP.json` |
| Phase O | Observable inheritance | PASS | `gates_summary_gateO.json` |
| Gate B | Superconducting transfer witness | structural PASS; **experimental window unresolved** | `gates_summary_gateB.json` |
| Gate C | Material independence (group-IV three-class collapse) | PASS, on schematic parameters | `gates_summary_gateC.json` |
| Gate D | Robustness and discriminability | PASS | `gates_summary_gateD.json` |
| Gate E | Is there a non-empty `I_engineerable ∩ I_asymptotic ∩ I_detectable`? | **optical NV FAIL, SC transfer FAIL, 3-mode conditional PASS** | `gates_summary_gateE.json` |
| A–D consolidated | Independent re-audit of Gates A–D | closed | `gates_a_to_d_consolidated_results.json` |
| RoomT 1–9 | Room-temperature campaign: low-T positive control, cut audit, moments, temperature scaling, global adversarial optimization, dip discrimination, correction bounds, reduced-vs-full, signal conversion | closed | `gates_summary_step1.json` … `step9.json` |

## EIT-definition-equivalence campaign (Phases E, F)

| Gate | Question | Verdict | Certificate |
|---|---|---|---|
| E1 | Λ-chain equivalence | closed | `gate_E1_lambda_chain.json`, `closed/E1_E3_findings.md` |
| E2 | Counterexample search | closed | `gate_E2_counterexamples.json` |
| E3 | 2g2e boundary | closed | `gate_E3_2g2e_boundary.json` |
| E4 | 2g3e search (+ fast scan) | closed | `gate_E4_2g3e_search.json`, `gate_E4_2g3e_fast_scan.json`, `gate_E4_FINAL.json`, `closed/E4_findings.md` |
| F0/F1/F3 | Convention audit | closed | `gate_F0_F1_F3_convention_audit.json` |
| F2 | Boundary scaling | closed | `gate_F2_boundary_scaling.json` |
| F4 | Transparency floor | closed | `gate_F4_transparency_floor.json`, `closed/F4_findings.md` |
| F5A | Unmatched readout | closed | `gate_F5A_unmatched_readout.json`, `closed/F5A_findings.md` |
| F5B | Full GKSL 2g3e | closed; **artifact was stale until 2026-08-10** | `gate_F5B_full_gksl.json`, `closed/F5B_findings.md` |
| F5C | Matched-floor accretive identity | PASS — replaced F5B's random scan with a proof | `gate_F5C_matched_floor_accretive.json`, `closed/F5C_findings.md` |

## Manuscript-side gates

| Gate | Question | Verdict | Certificate |
|---|---|---|---|
| P0-1 | NV observable freeze | closed | `closed/P0_1_NV_OBSERVABLE_FREEZE_2026-08-07.md` |
| P0-2 | Two-layer classification | PASS | `p0_2_two_layer_summary.json` (in `results/tables/`), `closed/P0_2_TWO_LAYER_CLASSIFICATION_2026-08-07.md` |
| P1–P6 | PRA phase diagram, threshold bands, observables, model audit, `B⊥` scaling, full uncertainty | closed | `p1_summary.json` … `p6_summary.json` (in `results/tables/`) |
| PRL gate audit | Do the assembled PRL claims survive their own gates? | closed | `closed/GENERAL_PRL_GATE_AUDIT.md` |

## Standing rules

- A gate's verdict is what its certificate says, not what its prose says. Where
  they disagree, the certificate wins and the prose is a defect to be fixed
  (this is exactly what happened with F5B).
- Random scans are not accepted as evidence for a positivity claim.
- A FAIL is an asset. Failures are kept in `../evidence/failures/` and cited in
  `../NON_CLAIMS.md`; they are not deleted when a later gate succeeds.
