# Common full-Liouvillian gate: NV vs group-IV (SiV-/SnV-)

Reproduces the "next hard gate" of `VALIDATION_REPORT.md` from the uploaded
surrogate package: run NV, SiV-, SnV- through one common computational
pipeline and report moment order, B_perp opening, and delta_chi_S -> EIT/ATS
verdict. Code: `src/{nv_reduced_kernel,group_iv_model,liouvillian_core,
eit_ats_classifier,moment_order_common_pipeline,run_common_liouvillian_gate,
bperp_full_lowfield_slope}.py`. Reproduce with
`python -m unittest tests/test_core.py` and the scripts listed above.

## 1. Moment order differs between NV and group-IV (same pipeline)

`moment_order_common_pipeline.py`, identical `A_Gamma=Gamma*D+A0` reduced-
kernel machinery for both:

| system | \|M0\| | fitted slope | predicted |
|---|---:|---:|---:|
| NV (0,-1) | 0 | -2.000 | -2 |
| NV (0,+1) | 0 | -2.000 | -2 |
| NV (-1,+1) | 0 | -3.000 | -3 |
| SiV- (orbital-Lambda) | 1 | -1.000 | -1 |
| SnV- (orbital-Lambda) | 1 | -0.985 | -1 |

NV's spin-Lambda legs are orbital-orthogonal (M0=0 identically); group-IV's
same-spin orbital-Lambda legs share a common excited ket at zeroth order
(M0=1). Three distinct moment classes appear (n=m+1 = 1, 2, 3) -- confirmed
in `run_common_liouvillian_gate.py` output (`common_liouvillian_gate.csv`).

## 2. B_perp^2 opening survives in the full Liouvillian

- Reduced kernel (`bperp_kernel_map_v2.py`, run for the first time --
  previously coded but never executed/archived): `delta_C_power = 2.026`.
- Full 9-level Lindblad model (`bperp_full_lowfield_slope.py`, using the
  completed `liouvillian_core.py`/`nv_system.py`/`weak_probe_response.py`):
  `slope_full_liouvillian = 2.178` vs `slope_reduced_kernel_this_run = 2.162`
  (same run, same Bx grid).
- `bperp_full_validation.py` (T=100/150/300 K, Bx up to 1 T) gives
  `C_full/C_red` ratios of 0.96-1.00 throughout, i.e. the full CPTP
  Lindbladian and the phenomenological resolvent agree to a few percent.

Conclusion: the B_perp^2 opening is not an artifact of the reduced
resolvent -- it survives (same order of magnitude exponent, ~2) in the
genuine Lindblad master-equation model.

## 3. delta_chi_S -> absorption -> EIT/ATS verdict

`eit_ats_classifier.py` implements the Anisimov-Dowling-Sanders (PRL 107,
163604, 2011) A_EIT/A_ATS model comparison with AIC/AICc and the fixed
robust-decision gate `|Delta_AIC|>=6` (`SiV_SnV_phonon_AIC_parameters.md`
Sec. 10.2). Sanity check on toy lineshapes: correctly recovers "robust EIT"
and "robust ATS" (`test_eit_ats_classifier_sanity`).

Applied to NV and group-IV absorption spectra generated from the SAME
Hamiltonians as Sec. 1 (`eit_ats_gate_run.py`, `common_liouvillian_gate.csv`):
at the default ground-coherence protection (gg=6.3e-5 GHz) used throughout
this repo, no configuration reaches a "robust EIT" verdict at Oc=1 GHz;
weak control is AIC-inconclusive (no resolvable delta_chi_S feature above
background), strong control (Oc=20 GHz) is "robust ATS" for several
configurations. This is consistent with -- and now makes quantitative --
the existing qualitative verdict of `hbn_nogo_EIT_report.md` ("no
configuration is currently a confirmed room-temperature go").

## Caveats

- Group-IV Hamiltonian parameters are fixed to the values in
  `SiV_SnV_phonon_AIC_parameters.md`; the orbital-Lambda leg-mixing angle
  theta=0 (maximal same-ket overlap) is a modeling choice, not a measured
  dipole ratio.
- `p_prot` is reported as 0 throughout: none of these concrete channels
  realize the abstract protected-subspace theorem (only demonstrated on
  block-random matrices in the uploaded surrogate package) -- this gate does
  not claim a protected group-IV or NV channel exists.
- phonon-rate normalization constants for SiV/SnV (`NORM_GHZ` in
  `group_iv_model.py`) are schematic (order-of-magnitude), since per-sample
  fit constants are not given in the source material.
