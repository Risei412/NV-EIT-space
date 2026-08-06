# No-Go Theorem: EIT in NV Centers

The original, binary go/no-go theory for this project. It asks whether an
EIT-type coherent response can survive in the NV center in diamond once
thermal dissipation is included, and answers it as a **practical no-go
boundary** in temperature.

**Central claim.** A transverse magnetic field opens a symmetry-suppressed
Raman pathway perturbatively, while thermal dissipation suppresses it and
defines the boundary above which the response is no longer observable. The
field exponent is temperature- and fit-window dependent in the full model;
only the perturbative regime should be quoted as a power law.

The three-tier generalization of this binary criterion lives in
`../New no-go theory/`; the redefinition of EIT itself lives in
`../EIT definition equivalence/`.

## Theorem package revision (2026-08-06)

`Theorem and proofs/eit_nogo_proofs.tex` is the compatibility entry point for
the canonical revised source
`Theorem and proofs/eit_nogo_proofs_v7.tex`.

The revision makes four publication-critical changes:

1. The large-dissipation order theorem is proved by the **Laurent expansion at
   infinity of the finite-dimensional rational resolvent**. It therefore does
   not require the conservative norm hypothesis `||X||/Gamma < 1`.
2. The **Neumann expansion is retained** as a quantitative lemma because it
   gives a directly checkable convergence region and an explicit truncation
   error, including a uniform-on-compact parameter version.
3. The former sector-graph remark is promoted to a theorem. Graph distance
   gives a rigorous lower bound on the first nonzero moment; equality also
   requires that the summed shortest-path amplitude does not cancel.
4. Previously abbreviated points are closed or scoped explicitly: the pure
   stationary-state algebra, the Neumann-free Krylov transfer-zero proof,
   observable-order inheritance, the Hermitian restriction of singular
   damping, and the higher-pole test for non-semisimple zero modes.

The previous v6.2 technical source remains recoverable from Git history. The
full v6.2 English write-up and its PDF are retained as historical theory
artifacts; theorem statements used in a new manuscript should be taken from
v7 where the two differ.

## Verification status

All five original verification gates were implemented and executed at the
fixed candidate point (T = 70 K, B⊥ = 0.23226 T, Bz = 0.005 T, Ωc = 0.1 GHz,
control on ms = +1). Numbers, caveats, and per-gate discussion:
`results/gate_1_5_report.md`. The gate definitions and pass criteria are in
`../SIMULATION_PLAN.md`.

| Gate | Question | Script |
|---|---|---|
| 1 | Is the candidate spectrum robustly EIT rather than ATS? | `src/gate1_candidate_aic_bootstrap.py` |
| 2 | Does the reduced kernel agree with the full Liouvillian? | `src/gate2_candidate_full_vs_reduced.py` |
| 3 | Is the signal detectable through a realistic signal chain? | `src/gate3_snr_map.py`, `src/signal_chain.py` |
| 4 | Where are the temperature thresholds, with uncertainty? | `src/gate4_threshold_uncertainty.py` |
| 5 | Does the contrast survive ensemble averaging? | `src/gate5_ensemble_average.py` |

## Layout

- `Theorem and proofs/` — the proof documents.
  `eit_nogo_proofs.tex` routes to the canonical v7 technical package;
  `eit_nogo_proofs_v7.tex` contains the Laurent, sector-graph, Krylov,
  observable-inheritance, and singular-point proofs;
  `EIT_no_go_go_theory_v6_2_English.tex` (with built PDF) is the historical
  full theory write-up; `eit_nogo_lecture.tex` is the lecture-format
  presentation.
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

Compile the revised theorem package from its directory:

```bash
cd "Theorem and proofs"
pdflatex eit_nogo_proofs.tex
```

Random seeds are fixed (20260716), so the archived numerical outputs are
reproducible run to run apart from explicitly recorded runtime metadata.

## Literature inputs

`results/literature_notes.md` groups the hBN and NV parameter literature by
the role it plays in the argument, and
`results/tables/literature_manifest.csv` records the bibliographic details and
the SHA-256 of every source file a number was read from. See
`../hBN papers/README.md` for the hBN bibliography.
