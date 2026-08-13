# evidence/benchmarks/

Independent checks of the repository against itself.

| Benchmark | Contents |
|---|---|
| `AUDIT_theorem_and_code_2026-08-01.md` | Line-level audit of `theory/proofs/eit_nogo_proofs.tex` (Theorems 1A–7B) **and** of the code: no mathematical error found; the novelty boundary drawn in §I-4; predicted re-execution values in §A-1 that later caught the stale F5B artifact. |
| Cross-library reproduction | 47 tests passed on a library generation different from the author's (numpy 2.4 / scipy 1.17 / pandas 3.0). Gate A–D, RoomT step1–9, Phase A/B/D/M/P and No-go gate1–5 JSON/CSV matched exactly apart from `runtime_s`; PNGs matched pixel-for-pixel. |
| Post-restructure reproduction (2026-08-13) | `pytest calculations/tests` — 70 tests pass against the reorganized layout; PRL figure regression tests reproduce the committed figures from `results/figures/`. |
