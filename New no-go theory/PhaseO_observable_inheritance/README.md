# Phase O — Observable-Order Inheritance (PRL Gate A)

Priority 1 of the PRL program (see
`Writing Paper/NV_EIT_PRA_PRL_Split_Strategy_20260724.md`, §13–14). This is
the make-or-break gate: if the kernel suppression order `ν_K` does not
inherit to the physical observable order `ν_obs` as a general law, the PRL
claim collapses into PRA. Here it is derived, certified, and verified.

## The law

For a bilinear (EIT-type) sector response with fast optical resolvent
`G(Γ,z) = [ΓD + A₀(z)]⁻¹`, probe leg `dp`, control leg `dc`,

```
R_obs(Γ,z) = −β · K12 · K21 / (g_eff + β·S2),
K12 = dp† G dc,   K21 = dc† G dp,   S2 = dc† G dc,
```

the suppression order inherits as

```
ν_obs = n12 + n21 − ν_den,
ν_den = 0        (generic:  g_eff ≠ 0, slow denominator → g_eff)
ν_den = ν_S2     (protected floor: g_eff = 0, slow denominator → β·S2)
```

where `n1j` = index of the first nonzero path moment of that cross-kernel
`+ 1`. **Selection-rule cancellation** (`dp†dc = 0`, …) promotes `n12`, `n21`
above the generic value `1`, which is exactly what lifts `ν_obs` above the
naive single-kernel order `ν_K`.

Three known repository data points are reconciled:

| observable | mechanism | result |
|---|---|---|
| NV EIT susceptibility (RoomT step3) | `M0 = dp†dc = 0` selection rule | `n12=n21=2 ⇒ ν_obs=4` |
| Quantum Fisher information (Phase M) | quadratic readout `|x_S|²` | `ν_obs = 2ν` |
| Phase N intervention scaling / resonance | `κ=κ₀Γ^q` fan; isolated `z★` | `ν(q)=4−q,2+q,4`; `z★` promotes 4→5 |

## Contents

- `src/gate_a_observable.py` — model-independent library: `SectorSpec`,
  `predict_nu_obs` (moment predictor), `certify_nu_obs_exact` (exact SymPy
  Γ-degree certificate), `verify_nu_obs_loglog` (full numeric slope),
  `separate_regimes` (pre-asymptotic vs asymptotic + `Γ_cross`), `check_model`.
  Reuses `New no-go theory/src/core.py` unchanged.
- `src/model_specs.py` — `SectorSpec` builders: a synthetic three-class family
  (`ν_obs = 2/4/6` by pure leg selection rules) and the physical NV EIT pencil
  (`No-go theorem/src/nv_model.py`).
- `src/run_gate_a.py` — driver; runs every observable class and writes the
  gate summary.
- `tests/test_gate_a.py` — fast regression tests pinning the known orders.

## Run

```
python src/run_gate_a.py            # full
python src/run_gate_a.py --quick    # fast CI-sized sweep
python tests/test_gate_a.py         # regression tests
```

Outputs: `results/tables/gates_summary_gateO.json`,
`results/tables/gate_a_models.csv`,
`results/figures/fig_gateO_inheritance.png` (+`.pdf`).

## Gate criteria (all currently PASS)

- **G-O1** inheritance `ν_obs = n12 + n21 − ν_den` confirmed for ≥3 distinct
  observables (bilinear susceptibility, quadratic QFI, frozen-source Raman/
  difference).
- **G-O2** generic vs symmetry-protected distinguished by the **exact**
  rational-degree certificate (not a float slope): `ν_obs = 2/4/6` for the
  synthetic family; NV `M0 = 0` exactly.
- **G-O3** predicted order matches the deep-tail log-log slope.
- **G-O4** pre-asymptotic effective exponent separated from the true
  asymptotic index; on physical NV `Γ_cross ≈ 7.5×10⁴ GHz > Γ(300 K)`, so the
  300 K response is pre-asymptotic (`ν_eff ≈ 3`) while the true index is 4.
- **G-O5** predictor / exact certificate / log-log fit agree across all
  models, and the Phase N q-fan and `z★` resonance corroborate the modifiers.

**Gate A verdict: PASS** → observable-order inheritance closes; the PRL may
proceed to Gate B (non-EIT physical witness — superconducting circuit).

## Conventions

Follows the repository conventions: numpy/scipy/sympy (no qutip), GHz /
`2π` / sign conventions per `No-go theorem/results/tables/convention_table.md`,
`Γ`-asymptotics decided by exact arithmetic and only corroborated by float
log-log fits.
