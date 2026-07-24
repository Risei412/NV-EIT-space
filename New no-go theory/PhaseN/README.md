# Phase N: SMRT Foundation — Sector-Mediated Response Theory

**Status:** All gates (Phase N N1–N5, Priority 1 U1–U7, Priority 2, Priority 3 F0–F9, Phase H H0–H5) **PASS**.

This directory contains the core SMRT numerical campaign verifying that the sector-mediated suppression order is determined by the joint scaling path \(κ=κ(Γ)\) through intervention strength and dissipation, not a fixed invariant of the system.

## Contents

### phase_n_base/
- `run_phase_n.py` — five-level GKSL witness model runner
- `gates_summary_phaseN.json` — gates N1–N5 results (all PASS)

Key result: full model has ν=2; three sectors have ν_S3=3, ν_S4=3, ν_S34=4 — demonstrating non-identifiability without sector intervention.

### priority_1_2/
- `phase_n_exact_core.py` — exact rational-moment machinery
- `run_phase_n_priority1_unfolding.py` — interference-point unfolding
- `run_phase_n_priority2_phase_diagram.py` — (q,α) phase diagram
- `PHASE_N_PRIORITY1_PRIORITY2_REPORT.md` — production report
- Results JSON (smoke + production) and 8 figures

Key result: V-shaped fan ν(q)=4-q, 2+q, 4 is an exact property of the tuned interference manifold, with universal crossover laws ∝ |δ|^{-1/(q-1)}.

### priority_3_frequency/
- `phase_n_frequency_core.py` — frequency-dependent response
- `run_phase_n_priority3_frequency.py` — production runner
- `PHASE_N_PRIORITY3_FREQUENCY_REPORT.md` — production report
- Results JSON (smoke + production) and 8 figures

Key result: isolated frequency z★=543/280 promotes pointwise order 4→5. Finite-band norms retain generic fan.

## Reproduce

```bash
pip install numpy scipy sympy matplotlib mpmath
python phase_n_base/run_phase_n.py
python priority_1_2/run_phase_n_priority1_unfolding.py
python priority_1_2/run_phase_n_priority2_phase_diagram.py
python priority_3_frequency/run_phase_n_priority3_frequency.py
```

## Publication status

These results are primary evidence for the SMRT central claim:

> The sector suppression order is not an invariant of the unperturbed generator or sector label. It is a valuation selected by the experimental path through intervention strength and probe frequency. Interference manifolds create exact promoted orders, while their neighborhoods exhibit universal crossover laws.
