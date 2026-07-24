# Gate D — Physical Significance (robustness/crossover + experimental discriminability)

PRL Priorities 6 & 8 (see `Writing Paper/NV_EIT_PRA_PRL_Split_Strategy_20260724.md`,
§13–16). This is the **final gate**; passing it completes the PRL central claim
(§14: "Gate D合格 = PRLの中心主張が完成する").

## P6 — Robustness / crossover

**Exact vs approximate class.** The integer response classes come in two kinds:

- **Exact (symmetry-protected).** NV ms=−1↔+1 (n=3): the leading moments
  `M0 = M1 = 0` are *structural* zeros (probe/control on orthogonal orbital
  branches, graph distance 2). Strain (`xi_x, xi_y`), transverse field `Bx`,
  and detuning `z` do **not** lift them — the order stays 3 with no crossover
  (`Γ* = ∞`). The high-order class is robust.
- **Approximate (tuned).** The superconducting protected class (ν=2) is a tuned
  cancellation `g_A1 g_B1 = −g_A2 g_B2`. A small coupling imbalance ε breaks it
  to ν=1 with a **crossover `Γ*(ε) ∝ 1/ε`** (numerically −1.00 in log Γ* vs
  log ε). The high-order class is observable only for `Γ < Γ*(ε)`, a window
  that opens to all Γ as ε→0.

The effective-exponent map `ν_eff(Γ, ε)` and the crossover scale `Γ*(ε)` are the
robustness fan. The intervention-scaling `ν(q) = 4−q, 2+q, 4` fan (Phase N,
`path_order`) corroborates that the crossover is a general Newton-polygon
phenomenon.

## P8 — Experimental discriminability

- **Required Γ dynamic range.** ~1–2 decades of clean asymptotic Γ pin a
  power-law exponent to ±0.1, i.e. resolve an adjacent class `Δν = 1`.
- **Platform reach** (relative to that requirement):
  - **Superconducting**: bus decay κ is directly engineerable over ≥9 decades →
    the exponent (class) is measured directly. Best platform.
  - **NV**: the phonon rate `Γ(T) = k_orb ∝ T⁵` spans **~8.5 decades over
    4–300 K** → the exponent is accessible by temperature tuning (in the
    asymptotic regime Γ ≫ internal splitting).
  - **group-IV**: the single-phonon Bose-law `Γ(T)` is narrow (SiV ~1.6, SnV
    ~0.6 decade) → temperature tuning alone is insufficient.
- **Optical read-out** (diamond): the ZPL OD→SNR chain (`signal_chain.py`) at a
  feasible point (1 ppm, 0.5 mm, 1 h, η=0.1, SNR=5) gives a minimum detectable
  contrast `C_min ≈ 1×10⁻⁴`, well below an order-10⁻² response — so the log-log
  slope (the class) is measurable.
- **Design rule.** Dissipation-tolerant response is engineered by inspecting the
  path moments / selection rules and interface geometry, not merely by seeking a
  small linewidth (§15.5) — a measurable rule, not just asymptotic mathematics.

## Contents

- `src/run_gate_d.py` — P6 (robustness/crossover) + P8 (discriminability) driver.
- `tests/test_gate_d.py` — regression tests.

Reuses, unchanged: Gate B `model_sc_transfer.py` (approximate-class crossover),
`No-go theorem/src/{nv_reduced_kernel, phonon_rates, group_iv_model, nv_model,
signal_chain}`, `New no-go theory/src/core.py`, and
`PhaseN/priority_1_2/phase_n_exact_core.py` (ν(q) fan). No new physics model —
this is a robustness/experimental-conversion analysis of the existing kernels.
numpy/scipy/sympy only (no qutip).

## Run

```
python src/run_gate_d.py            # full
python src/run_gate_d.py --quick    # fast
python tests/test_gate_d.py         # regression tests
```

Outputs: `results/tables/gates_summary_gateD.json`, `gate_d_robustness.csv`,
`gate_d_discriminability.csv`, `results/figures/fig_gateD_robustness.png` and
`fig_gateD_gammaT.png` (+`.pdf`).

## Gate criteria (all currently PASS)

- **G-D1** exact class (NV, unbreakable) vs approximate class (SC, breaks with ε)
  distinguished.
- **G-D2** crossover scale law `Γ*(ε) ∝ 1/ε` (fitted power −1.0).
- **G-D3** effective-exponent map + observable window; ν(q) fan `4−q, 2+q, 4`.
- **G-D4** required Γ dynamic range to resolve `Δν = 1` computed.
- **G-D5** Γ(T)/platform reach: NV (~8.5 dec) and SC (≥9 dec) exceed the
  requirement; group-IV is narrow.
- **G-D6** optical SNR: finite `C_min`, response detectable → slope measurable.

**Gate D verdict: PASS.**

## PRL program complete

With Gate D, the four gates are all green:

| gate | result | campaign |
|---|---|---|
| A | observable-order inheritance `ν_obs = n12+n21−ν_den` | `PhaseO_observable_inheritance/` |
| B | non-EIT non-diamond superconducting witness + blind prediction | `GateB_superconducting_witness/` |
| C | material independence (group-IV full-GKSL + 3-class collapse) | `GateC_material_independence/` |
| D | robustness (exact vs approximate) + experimental discriminability | `GateD_robustness_discriminability/` |

The strategy's minimum PRL conditions (§21) — observable-level order law,
non-EIT full physical witness, non-diamond blind prediction, plus
robustness/crossover and experimental discriminability — are all satisfied. What
remains is the writing (≈4 main figures + Supplement) and the separate PRA line.
