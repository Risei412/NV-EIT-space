# Phase H: Physical GKSL Hidden Class Transition

**Status:** All gates (H0–H5) **PASS**.

This directory contains the physical realization of the central SMRT phenomenon in a minimal 4-level GKSL model: the sector-mediated response transitions from Class III (protected, ν=0) to Class II (ν=1) while the full observable response remains in Class III.

## Contents

- `run_phase_h.py` — four-level Lindblad model with tuned hidden point
- `physical_hidden_transition_gksl.md` — complete report of gates H0–H5

## Key result

At a non-trivial physical coupling J₁₂★ = 0.115 + 0.205i:
- The full response χ_full remains Class III (ν=0, finite O(1) plateau)
- The sector response R_S = χ_full - χ_cut transitions to Class II (ν=1, ∝ Γ⁻¹)
- The transition is invisible to single-readout classification
- Requires at least two independent observables or a multi-sector decomposition to detect

This is the physical foundation of SMRT: the hidden mechanism of order suppression can exist even when the measured coherence appears protected.

## Reproduce

```bash
pip install numpy scipy matplotlib
python run_phase_h.py
```

Expected outputs:
- `results/gates_summary_phaseH.json`
- `results/figures/figH1_physical_hidden_transition.png`

**Note:** The core model file (`model_hidden_physical.py`) was 0 bytes in the original bundle. The report and runner are included, but full reproduction requires the missing Lindblad coefficient definitions (can be reconstructed from the .md report section 2–3).
