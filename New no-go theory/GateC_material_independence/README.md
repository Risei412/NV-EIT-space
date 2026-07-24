# Gate C — Material Independence (group-IV full-GKSL + three-class collapse)

PRL Priorities 4 & 5 (see `Writing Paper/NV_EIT_PRA_PRL_Split_Strategy_20260724.md`,
§13–14). Gate A closed observable-order inheritance; Gate B gave a non-EIT,
non-diamond (superconducting) witness with a blind prediction. Gate C
establishes **material independence**: the integer suppression classes appear
across diamond and non-diamond systems, and the group-IV Γ⁻¹ class is certified
in a **full physical GKSL** rather than only a reduced kernel.

## Results

**P4 — group-IV SiV⁻/SnV⁻ full GKSL.** `group_iv_full.py` builds a genuine
5-level density-matrix Lindblad master equation for the orbital-Λ (|g⟩ + the
4-state excited manifold of `group_iv_model.H_groupIV`). The swept fast
dissipation Γ is realized as excited-state pure dephasing
`L_j = √(2Γ)|e_j⟩⟨e_j|` (uniform D = I, exactly the reduced kernel's `A_Γ =
ΓI + A0`) or physical phonon orbital hopping `√Γ|e_y⟩⟨e_x|` (full-rank D). The
fixed-ρ₀ weak-probe response (Gate-B pattern) reproduces the reduced amplitude
kernel `K(Γ) = p†[ΓD+A0]⁻¹c` **to machine precision** (~2×10⁻¹⁶), with slope −1
because `M0 = p†c ≠ 0` (same-spin orbital-Λ).

**P5 — three physical suppression classes collapse under Γⁿ.**

| class | system | mechanism | Γ-scaling |
|---|---|---|---|
| n=1 | group-IV orbital-Λ | `M0 ≠ 0` | Γ⁻¹ |
| n=2 | NV ms=0↔−1 | graph distance d=1 | Γ⁻² |
| n=3 | NV ms=−1↔+1 | graph distance d=2 | Γ⁻³ |

The first-nonzero-moment (graph-distance) predictor fixes the integer exponent
before any fit; `Γⁿ|K(Γ)|` collapses to a plateau (the first path moment).

**Material independence.** The same integer class appears in **diamond and
non-diamond**:
- class 1: group-IV (diamond) **and** the Gate-B superconducting generic
  witness (non-diamond),
- class 2: NV (diamond) **and** the superconducting protected witness
  (non-diamond).

The plateau *value* (the path moment) is material-specific; the integer
*exponent* is shared — that is the universality claim. (P3, a diamond-external
witness, was already delivered in Gate B; this overlay reuses it.)

Literature anchors (guide lines only; the exponent is swept-Γ independent of the
absolute phonon rate): SiV⁻ Δ_e=255 GHz, ground SO ~48 GHz, ground orbital T1≈39
ns@5K (~0.026 GHz), excited lifetime ~1.7 ns; SnV⁻ Δ_e=3000 GHz, ground SO ~850
GHz. Sources: Pingault et al., *Nat. Commun.* **8**, 15579 (2017); Jahnke et al.,
*New J. Phys.* **17**, 043011 (2015)/arXiv:1411.2871; Trusheim et al., *PRX*
**11**, 041041 (2021); Meesala et al., *PRB* **97**, 205444 (2018).

## Contents

- `src/group_iv_full.py` — group-IV full GKSL (P4): `build_H_full`,
  `dipole_legs_full`, `jump_operators` (dephasing/hopping + radiative),
  `full_response` (fixed-ρ₀ linear response), `reduced_kernel_response`, `M0`.
- `src/run_gate_c.py` — driver: P4 certification, P5 three-class collapse,
  diamond/non-diamond overlay, gates.
- `tests/test_gate_c.py` — regression tests.

Reuses `No-go theorem/src/{group_iv_model, nv_reduced_kernel,
verify_nv_3E_graph_distance_PRL, moment_order_common_pipeline, liouvillian_core,
phonon_rates}`, `New no-go theory/src/core.py`, and Gate B's
`model_sc_transfer.py`, unchanged. numpy/sympy only (no qutip). The parameter
file `group_iv_model.py` cites (`SiV_SnV_phonon_AIC_parameters.md`) is absent
from the tree; the anchors above live as constants in `group_iv_full.py`.

## Run

```
python src/run_gate_c.py            # full
python src/run_gate_c.py --quick    # fast
python tests/test_gate_c.py         # regression tests
```

Outputs: `results/tables/gates_summary_gateC.json`, `results/tables/gate_c_collapse.csv`,
`results/figures/fig_gateC_three_class_collapse.png` (+`.pdf`).

## Gate criteria (all currently PASS)

- **G-C1** group-IV full-GKSL slope −1 (dephasing *and* hopping); reduced==full
  to ~3×10⁻¹⁶; `M0≠0`; `Γ·R→M0`; ρ₀ steady.
- **G-C2** three physical classes: graph-distance predicted exponents (1,2,3)
  match the log-log slopes.
- **G-C3** `Γⁿ|K|` collapses to a plateau (rel spread <5%), exponent stable
  across fit sub-windows.
- **G-C4** material independence: class 1 slope=1 in group-IV & SC-generic;
  class 2 slope=2 in NV & SC-protected.
- **G-C5** full/reduced agree; pre-asymptotic vs asymptotic separated.

**Gate C verdict: PASS** → the integer classes are material-independent (diamond
and non-diamond) and group-IV Γ⁻¹ is certified in a full physical GKSL. PRL may
proceed to Gate D (robustness/crossover map, experimental discriminability).
