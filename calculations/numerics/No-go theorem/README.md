# No-Go Theorem: EIT in NV Centers

The original, binary go/no-go theory for this project. It asks whether an
EIT-type coherent response can survive in the NV center in diamond once
thermal dissipation is included, and answers it as a **practical no-go
boundary** in temperature.

**Central claim.** A transverse magnetic field opens a symmetry-suppressed
Raman pathway quadratically (K ∝ B⊥²), producing an EIT-type response at
finite temperature, while thermal dissipation systematically suppresses it and
defines the boundary above which the response is no longer observable.

The three-tier generalization of this binary criterion lives in
`../New no-go theory/`; the redefinition of EIT itself lives in
`../EIT definition equivalence/`.

## Verification status

All five PRL verification gates were implemented, executed, and **passed** at
the fixed candidate point (T = 70 K, B⊥ = 0.23226 T, Bz = 0.005 T, Ωc = 0.1 GHz,
control on ms = +1). Numbers, caveats, and per-gate discussion:
`results/gate_1_5_report.md`. The gate definitions and pass criteria are in
`../gates/active/SIMULATION_PLAN.md`.

| Gate | Question | Script |
|---|---|---|
| 1 | Is the candidate spectrum robustly EIT rather than ATS? | `src/gate1_candidate_aic_bootstrap.py` |
| 2 | Does the reduced kernel agree with the full Liouvillian? | `src/gate2_candidate_full_vs_reduced.py` |
| 3 | Is the signal detectable through a realistic signal chain? | `src/gate3_snr_map.py`, `src/signal_chain.py` |
| 4 | Where are the temperature thresholds, with uncertainty? | `src/gate4_threshold_uncertainty.py` |
| 5 | Does the contrast survive ensemble averaging? | `src/gate5_ensemble_average.py` |

## Layout

- `Theorem and proofs/` — the proof documents. `eit_nogo_proofs.tex` is the
  technical package (Theorems 2A/2B, 3–6, reused by the newer theories);
  `EIT_no_go_go_theory_v6_2_English.tex` (with built PDF) is the full theory
  write-up; `eit_nogo_lecture.tex` is the lecture-format presentation.
- `src/` — the models and gate implementations: NV Hamiltonian and reduced
  kernel (`nv_model.py`, `nv_system.py`, `nv_reduced_kernel.py`), phonon rates
  (`phonon_rates.py`, `happacher_rate.py`), Liouvillian core
  (`liouvillian_core.py`), weak-probe response (`weak_probe_response.py`),
  EIT/ATS classification (`eit_ats_classifier.py`), transverse-field scans
  (`bperp_*.py`), the group-IV comparison (`group_iv_model.py`), and the hBN
  scan (`hbn_nogo_scan.py`).
- `scripts/reproduce_prl_figures.py` — single entry point that reruns the gate
  campaign and regenerates the figures.
- `tests/` — `pytest` suite covering the core kernel, the operational-cut
  equivalence, and the phonon-rate variants.
- `results/` — archived outputs: `tables/` (CSV/JSON gate summaries),
  `figures/`, `metadatas/` (NPZ parameter maps and run metadata), plus the
  written reports `gate_1_5_report.md`, `common_liouvillian_gate_report.md`,
  `hbn_nogo_EIT_report.md`, and `literature_notes.md`.

## Reproduce

```bash
pip install -r requirements.txt
python scripts/reproduce_prl_figures.py
python -m pytest tests/
```

Random seeds are fixed (20260716), so the archived numbers in `results/` are
reproducible run to run.

## Literature inputs

`results/literature_notes.md` groups the hBN and NV parameter literature by
the role it plays in the argument, and
`results/tables/literature_manifest.csv` records the bibliographic details and
the SHA-256 of every source file a number was read from. See
`../hBN papers/README.md` for the hBN bibliography.
