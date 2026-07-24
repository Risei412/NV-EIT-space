# PRL Main-Text Figures

Publication figures assembled from the PRL Gate A–D results (see
`Writing Paper/NV_EIT_PRA_PRL_Split_Strategy_20260724.md`, §17: "PRL本文では
4図程度へ圧縮"). **No new physics** — every panel reuses the existing gate
campaigns' functions unchanged; this directory is presentation-only.

## Figures

- **Fig. 1 — general structure + three integer classes.** (a) schematic
  `c → D → A0 → p` with the selection-rule chain; (b) `|K(Γ)|` for the three
  physical classes n=1 (group-IV), n=2, n=3 (NV), slopes 1.00/2.00/3.00.
- **Fig. 2 — observable-order inheritance + crossover.** (a) `ν_obs =
  n₁₂+n₂₁−ν_den` on the synthetic generic/protected/doubly-protected family
  (2/4/6); (b) NV pre-asymptotic vs asymptotic crossover, `Γ_cross` and
  `Γ(300 K)` marked.
- **Fig. 3 — material independence + blind prediction.** (a) `Γⁿ|K|`
  normalized to its own plateau, overlaying diamond (group-IV, NV) and
  non-diamond (superconducting) systems at the same integer n; (b)
  superconducting transfer efficiency, blind-predicted `κ⁻²`/`κ⁻⁴`.
- **Fig. 4 — robustness + experimental discriminability.** (a) exact (NV,
  unbreakable) vs approximate (superconducting, tuned) class; `ν_eff(Γ)`
  family and crossover `Γ*(ε) ∝ 1/ε` (inset); (b) `Γ(T)` platform reach (NV
  ~8.5 decades, group-IV narrow) against the required Γ range.

## Contents

- `src/prl_style.py` — shared style: Okabe–Ito colorblind-safe palette
  (validated with the dataviz skill's `validate_palette.js`, all six checks
  PASS), PRL column widths, panel labels, PDF+PNG export.
- `src/fig1_classes.py` … `src/fig4_robustness.py` — one figure per file.
- `src/make_figures.py` — driver, `--quick` for a faster sweep.
- `tests/test_figures.py` — regression tests (outputs exist; key slopes/values
  match the certified Gate A–D numbers).

Data sources (imported unchanged): `No-go theorem/src/nv_reduced_kernel.py`,
`New no-go theory/GateC_material_independence/src/group_iv_full.py`,
`New no-go theory/PhaseO_observable_inheritance/src/{gate_a_observable,
model_specs}.py`, `New no-go theory/GateB_superconducting_witness/src/
model_sc_transfer.py`, `New no-go theory/GateD_robustness_discriminability/
src/run_gate_d.py`, `New no-go theory/src/core.py`.

## Run

```
python src/make_figures.py            # full
python src/make_figures.py --quick    # fast
python tests/test_figures.py          # regression tests
```

Outputs: `figures/fig{1..4}_*.pdf` (vector, for submission) and
`figures/fig{1..4}_*.png` (dpi=300, preview).
