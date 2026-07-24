# Gate B — Non-EIT / Non-Diamond Physical Witness (Superconducting)

PRL Priorities 2 & 7 (see `Writing Paper/NV_EIT_PRA_PRL_Split_Strategy_20260724.md`,
§13–14). Gate A (`../PhaseO_observable_inheritance/`) closed the observable-order
inheritance law; Gate B shows the **same path-moment law predicts a different
integer order in a physical system that is neither EIT nor diamond** — a
superconducting circuit — with the prediction made **blind** (before the full
model is solved).

## Witness

Circuit-QED **dissipative state transfer through a lossy bus**. An input
transmon qubit A and output qubit B are coupled through two lossy bus modes
R1, R2 (bus decay `κ` = the fast dissipation `Γ`). The observable is the
**A→B transfer amplitude / efficiency** — a conversion/transport quantity, not
an EIT dark-state transparency.

```
A --gA1-- R1(κ) --gB1-- B        two interfering transfer paths
A --gA2-- R2(κ) --gB2-- B
```

Adiabatic elimination of the lossy bus gives the effective coupling
`Σ(κ) = Σ_i gAi gBi/(κ + iΔ_i)`, so the transfer amplitude
`K(κ) = p_B†[κD + A0]⁻¹c_A` (with **singular** D acting only on the bus) has a
**selection-rule-controlled integer order**:

| tuning | selection rule | amplitude | efficiency `|K|²` |
|---|---|---|---|
| generic | `Σ gAi gBi ≠ 0` | `κ⁻¹` (ν_K=1) | `κ⁻²` |
| protected | `gA1gB1 = −gA2gB2`, `Δ1≠Δ2` | `κ⁻²` (ν_K=2) | `κ⁻⁴` |
| broken (ε) | small imbalance | ν_K=1 asymptotically | crossover 2→1 |

Fully degenerate (`Δ1=Δ2`) gives an exact all-order cancellation → Class I
(ν=∞); recorded as a schematic aside, not the main gate.

Literature-grounded parameters (GHz): `g/2π ≈ 0.1`, `α/2π ≈ −0.3`,
`κ` tunable, `γ ≈ 2.5×10⁻⁵` (T1 ≈ 40 µs). Sources: Krantz et al., *Appl. Phys.
Rev.* **6**, 021318 (2019); Blais et al., *Rev. Mod. Phys.* **93**, 025005
(2021); lossy-bus transfer — Burkhart et al., *PRX Quantum* **2**, 030321
(2021), arXiv:2005.12334, arXiv:2403.02203.

## Contents

- `src/model_sc_transfer.py` — witness model: reduced single-excitation pencil
  `A(κ,z)=κD+A0` (singular D), full 5-level GKSL (convention-free basis-loop
  Liouvillian), `implicit_linear_response`, `tuning` switch
  generic/protected/broken. numpy/sympy only (no qutip).
- `src/run_gate_b.py` — driver with the explicit **blind protocol** and gates.
- `tests/test_gate_b.py` — regression tests pinning ν_K = 1/2 and reduced==full.

Reuses `New no-go theory/src/core.py` (`certificate_deg_nu` — valid for
singular D, `fit_nu_loglog`) and `PhaseO_observable_inheritance/src/gate_a_observable.py`
(`SectorSpec`/`certify_nu_obs_exact`/`verify_nu_obs_loglog` in quadratic
readout mode = efficiency `|K|²`), unchanged. Note: the moment-ladder path
(`nu_from_moments`/`predict_nu_obs`) is **not** used — it assumes invertible D,
whereas the bus dissipation D is singular; the blind order comes from the exact
symbolic κ-degree certificate instead.

## Run

```
python src/run_gate_b.py            # full
python src/run_gate_b.py --quick    # fast
python tests/test_gate_b.py         # regression tests
```

Outputs: `results/tables/gates_summary_gateB.json` (with the `blind_prediction`
recorded first), `results/tables/gate_b_scaling.csv`,
`results/figures/fig_gateB_transfer_scaling.png` (+`.pdf`).

## Blind protocol (P7)

1. Build the reduced pencil from `(H, jumps, input, readout)` only.
2. `certificate_deg_nu` → ν_K (generic 1, protected 2); efficiency `2ν_K` (2, 4).
3. Record in `blind_prediction` **before** any full-model fit; class definition
   is not changed afterward.
4. Full GKSL κ-sweep; log-log deep-tail slope confirms `2ν_K`.

## Gate criteria (all currently PASS)

- **G-B1** non-EIT witness well-posed: reduced pencil == full-Liouvillian
  coherence sector (max rel err ≈ 3×10⁻¹⁰; `|g><g|` steady residual 0).
- **G-B2** blind prediction fixed before the full sweep: generic ν_K=1,
  protected ν_K=2 (efficiency 2, 4).
- **G-B3** full-GKSL deep-tail slope matches the blind prediction: generic
  efficiency 1.99, protected 3.99.
- **G-B4** selection rule changes the order, distinguished by the exact
  rational-degree certificate (efficiency 2 vs 4).
- **G-B5** small symmetry breaking ε → high→low crossover; `κ*(ε) ∝ 1/ε`
  (grows as ε→0); protected (ε=0) stays ν=2.

**Gate B verdict: PASS** → a non-EIT, non-diamond physical witness confirms the
path-moment law predicts (blind) a selection-rule-controlled integer order. PRL
may proceed to Gate C (further non-diamond witnesses, group-IV full-GKSL,
three-class scaling collapse). Since this witness is itself non-diamond, it also
advances part of Gate C.
