# evidence/failures/

Failures are primary assets here. Each one narrows `../../CLAIMS.md`, and each
is cited by an entry in `../../NON_CLAIMS.md`. Nothing in this directory gets
deleted because a later result succeeded.

| Failure | Record |
|---|---|
| **Optical NV cannot resolve its own exponent.** Usable window 1.02 / 0.84 / 0.51 decade (single / post-selected shimmed / high density) against 1 decade required. G-D7 slope budget, 2026-08-02. | `../../results/certificates/gates_summary_gateE.json`, `../../results/tables/gate_e_windows.csv` |
| **Gate B's experimental window is unresolved.** Effective exponent ≈ 0 at 0.1–50 MHz; integer asymptotics 8+ orders above physical bath damping. | `../../results/certificates/gates_summary_gateB.json` |
| **The general "measurable design rule" claim fails Gate E.** Two of three platforms FAIL outright. | `../../gates/closed/GENERAL_PRL_GATE_AUDIT.md` |
| **Claim/evidence gaps found in the PRL package.** | `PRL_CLAIM_GAP_AUDIT.md` (this directory) |
| **`gate_F5B_full_gksl.json` was stale and self-contradictory** until 2026-08-10 — the committed verdict came from a pre-bugfix run. | `../../gates/closed/F5B_findings.md` §0; `../../theory/LIMITATIONS.md` |
| **F5B's 150-draw random scan was inadequate evidence** for a positivity claim and was superseded by F5C's identity proof. | `../../archive/superseded_claims/` |
