# NV-Side Campaigns of the Three-Tier Response Classification

The three-tier classification upgrades the binary EIT go/no-go criterion (see
`../No-go theorem/`) to an exclusive trichotomy of the suppression index
ν ∈ {∞} ∪ (0,∞) ∪ {0} - exact structural no-go, asymptotic no-go, protected go
- for sector-mediated coherent response in finite-dimensional Markovian
weak-probe systems.

> **The theory itself now lives in a separate repository.**
> Theorems I-III and their proofs, the operational cut, and the campaigns that
> are independent of NV physics were moved to
> **[Sector-Master-Resolved-Theory](https://github.com/Risei412/Sector-Master-Resolved-Theory)**
> (SMRT), which is where the theory is developed.
>
> What remains here is the material that is bound to NV-center physics or that
> feeds this repository's manuscript figures.

## What moved to the SMRT repository

| Moved | Now at |
|---|---|
| Theorem documents and proofs | `theory/` |
| `three_theorem_verification.py` | `src/` |
| Operational cut, gates U1-U8 | `Sector/` |
| Exact path fan | `PhaseN/` |
| Detuning campaign | `PhaseZ/` |
| Physical hidden-class transition | `PhaseH/` |
| SMRT plan and revision notes | `docs/` |

`PhaseH` was incomplete here - `model_hidden_physical.py` and
`results/gates_summary_phaseH.json` were never committed - and is restored and
runnable in the SMRT repository.

## What stays here

| Directory | Contents | Why it stays |
|---|---|---|
| `PhaseO_observable_inheritance/` | Gate A: observable-order inheritance `ν_obs = n₁₂+n₂₁−ν_den`, and the NV pre-asymptotic/asymptotic crossover | Imports `nv_model` from `../calculations/numerics/No-go theorem/src/`; produces PRL Fig. 2 |
| `GateB_superconducting_witness/` | Non-diamond witness: superconducting transfer efficiency, blind-predicted κ⁻²/κ⁻⁴ | Produces PRL Fig. 3(b) |
| `GateC_material_independence/` | Material independence: group-IV and NV at the same integer n | Imports `group_iv_model`; produces PRL Fig. 3(a) |
| `GateD_robustness_discriminability/` | Exact vs approximate class, `ν_eff(Γ)`, crossover `Γ*(ε) ∝ 1/ε`, platform reach | Imports `signal_chain`; produces PRL Fig. 4 |
| `GateE_NV_experimental_anchor/` | NV ensemble contrast-to-floor window and finite-Γ exponent-identifiability audit | Determines whether the NV experiment is a direct exponent witness or a model-constrained PRL anchor |
| `RoomT/` | Steps 1-9 of the NV room-temperature no-go: low-temperature validation, merged-manifold moments, temperature scaling, adversarial optimization, dip discrimination, correction mechanisms, reduced-vs-full Liouvillian, signal conversion | Imports `liouvillian_core`, `nv_model`, `phonon_rates`, `gate2_candidate_full_vs_reduced` from `../calculations/numerics/No-go theorem/src/` |
| `src/` | Shared models and the Phase A/B/D/M/P sector-response campaigns | Imported by all of the above and by `../calculations/numerics/manuscript_figures/prl/` |
| `results/` | Archived gate summaries and figures for those campaigns | Outputs of `src/run_phase_*.py` |

## Modules duplicated across the two repositories

These exist in both repositories by design - each side needs them and neither
should depend on the other's checkout:

- `src/core.py` - transfer function, Krylov certificate, moment method,
  Riesz-projection protected coefficient, log-log ν fit.
- `src/model_metro_lindblad.py`
- `src/operational_cut.py` - used here by `RoomT/src/step2_operational_cut_audit.py`
  and by `../calculations/tests/nogo/test_operational_cut_equivalence.py`, which
  checks that gate 2's block-zeroed cut agrees with the ideal operational cut.
- `src/phase_n_exact_core.py`, `src/phase_n_frequency_core.py` - used here by
  `PhaseO_observable_inheritance/`.

## Reproduce

```bash
pip install -r requirements.txt

python PhaseO_observable_inheritance/src/run_gate_a.py
python GateB_superconducting_witness/src/run_gate_b.py
python GateC_material_independence/src/run_gate_c.py
python GateD_robustness_discriminability/src/run_gate_d.py
python GateE_NV_experimental_anchor/src/run_gate_e.py
python GateE_NV_experimental_anchor/tests/test_gate_e.py
python src/run_phase_a.py    # also run_phase_{b,d,m,p}.py
```

Tests live in each campaign's `tests/` directory. Gate results are archived as
JSON/CSV under each campaign's `results/`, and `results/summary.md` collects the
Phase A/B/D/M/P outcomes.

## Planning documents

- `new_nogo_numerical_priorities.md` - the priorities implemented by the
  sector-resolved response calculation in `src/`: computing χ_full, the
  frozen-source sector cut χ_cut, and the difference R_S = χ_full − χ_cut
  together (never χ_full alone), extracting ν three independent ways, and
  searching for a hidden class transition where ν[χ_full] stays 0 while
  ν[R_S] jumps 0→1.
- `phase_d_beyond_old_theory_plan.md` - the Phase D plan.
