# evidence/counterexamples/

Systems that break a candidate statement.

| Counterexample | Where |
|---|---|
| E2 counterexample search — systems defeating the naive equivalence | `../../results/certificates/gate_E2_counterexamples.json`; generator: `../../calculations/numerics/EIT definition equivalence/src/model_counterexamples.py` |
| F5A / F5B unmatched-readout exact zeros — matched readout cannot reach `x ∈ ker Γ`, unmatched Raman-basis coherence readout does | `../../results/certificates/gate_F5A_unmatched_readout.json`, `gate_F5B_full_gksl.json`; findings in `../../gates/closed/F5A_findings.md`, `F5B_findings.md` |
| Hidden class transition — a class change invisible in the total spectrum, so `χ_full` is a counterexample to itself as a diagnostic | `../../results/certificates/gates_summary_phaseB.json` (Gate 2) |
| Phase D / P masquerade — a system that imitates another class's signature | `../../results/certificates/gates_summary_phaseD.json`, `gates_summary_phaseP.json`; figures `figD4_masquerade.png`, `figP5_masquerade.png` |
| Phase M rank-deficiency negative control — as `ε → 0` the clean `ν_F = 2ν_x` relation degrades, showing the full-rank assumption is load-bearing | `../../results/certificates/gates_summary_phaseM.json` (Gate M5) |

The artifacts live in `results/` because they are gate outputs; this file is
the index that says *which* of them refute *what*.
