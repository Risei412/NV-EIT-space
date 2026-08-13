# PRA campaign — NV EIT: where it lives and where it inverts

The Physical Review A paper of the two-paper split described in
`../NV_EIT_PRA_PRL_Split_Strategy_20260724.md`: the microscopic mechanism by
which NV$^-$ spin-Λ EIT collapses on warming, the temperature window it
occupies, the transverse-field pathway, and the room-temperature no-go —
quantified.

The strategy document §5 lists four calculations that had to exist before the
manuscript could be closed (P1–P4). P1 did not exist, P2 existed only at one
field, P3 existed as two mutually incommensurate halves, and P4 had never been
assembled. This directory contains all four, plus two that the audit forced:

| Script | Delivers |
|---|---|
| `src/p1_phase_diagram.py` | Full-Liouvillian *T*–*B*⊥ phase diagram, every point classified by sign + lineshape + EIT/ATS model comparison |
| `src/p2_threshold_bands.py` | Boundary-temperature uncertainty **as a function of *B*⊥**, with the optical linewidth given its own prior and the four phonon-rate models as a systematic |
| `src/p3_observables.py` | One detection chain across the whole temperature range, transmission **and** fluorescence, on an optical-depth-matched common axis |
| `src/p4_model_audit.py` | The consolidated audit: reduced-vs-full agreement domain, truncation, provenance, sign conventions, classification stability |
| `src/p5_bperp_scaling.py` | The *B*⊥² opening exponent, with the zero-field residual pathway subtracted |
| `src/p6_full_uncertainty.py` | Full-Liouvillian re-derivation of the P2 bands, required because P4 found the reduced kernel invalid where the sign reversal sits |
| `src/make_figures.py` | The six main-text figures |

Nothing here re-derives physics that already existed. P1 composes
`gate2_candidate_full_vs_reduced.build_full` (the validated nine-level
Lindblad pipeline) with `gate1_candidate_aic_bootstrap.fit_all` (the four-model
AIC comparison) and adds an adaptive two-photon window; P2 extends
`gate4_threshold_uncertainty`; P3 reuses `signal_chain.py` and reconciles
`gate3_snr_map.py` with `RoomT/step9_signal_conversion.py`. Only two functions
are parameter-exposed copies rather than direct reuse — `p2.contrast` and
`p6.build_full_p` — following the precedent set by
`gate4_threshold_uncertainty.contrast`, whose docstring records the same
choice.

## Results

Every run writes a `*_summary.json` carrying its own pass gates alongside the
numbers, and the per-point CSV it was computed from, in `results/tables/`.

The headline findings:

- **The transparency region is a closed island**, not a half-space below a
  ceiling. It needs *B*⊥ ≳ 0.15 T to exist at all, closes below ≈22 K through
  an Autler–Townes crossover, and closes above 90–95 K through phonon-driven
  collapse of the leading Raman path.
- **Contrast falls 4.5 decades** from 0.99 at 30 K to 2.1×10⁻⁵ at 100 K
  (*B*⊥ = 0.232 T), then **changes sign** at 100 [94, 108] K.
- **The transverse field is a switch, not a dial**: it opens the pathway as
  *B*⊥², but moving from 0.15 T to 0.50 T shifts the boundary by under 2 K.
- **The sign reversal is measurable** — 1.6 s of integration at 105 K — while
  the room-temperature residual (1.1×10⁻⁹) is nine orders below the floor.
- **The reduced kernel reproduces boundary temperatures but not magnitudes**,
  and fails qualitatively (including in sign) below 40 K and above 90 K.

## Reproducing

```bash
pip install -r "../../No-go theorem/requirements.txt"
cd src
# Pin BLAS to one thread per worker: the workers are already parallel, and
# oversubscribing turns a 3-minute run into a 30-minute one.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python p1_phase_diagram.py --jobs 4      # ~5 min   (240 grid points)
python p2_threshold_bands.py --jobs 4    # ~1 min   (10 fields x 500 samples)
python p3_observables.py --jobs 4        # ~30 s
python p4_model_audit.py --jobs 4        # ~4 min
python p5_bperp_scaling.py --jobs 4      # ~1 min
python p6_full_uncertainty.py --jobs 4   # ~8 min   (80 full-Liouvillian samples)
python make_figures.py
```

Add `--quick` to any of them for a reduced grid. `p1` validates itself against
the archived 70 K candidate (`gate2_candidate_comparison.json`, C = 0.013836)
and reports the relative difference; it currently reproduces it to 0.25%.

## Manuscript

`manuscript/nv_eit_pra.tex` (REVTeX 4.2, `aps,pra,reprint`), with
`references.bib` copied from the curated literature package. Build:

```bash
cd manuscript
pdflatex nv_eit_pra && bibtex nv_eit_pra && pdflatex nv_eit_pra && pdflatex nv_eit_pra
```

Figures are read from `../results/figures/`, so run `make_figures.py` first.

## Scope kept out

Per strategy §20 and §22, the general pathway-order classification — observable
inheritance, path-moment orders, non-EIT and non-diamond witnesses — stays in
the companion PRL and is referred to, not developed, here. Group-IV centres
appear as a design comparison only.
