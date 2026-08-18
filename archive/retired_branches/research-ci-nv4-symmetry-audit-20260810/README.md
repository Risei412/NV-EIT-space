# Retired branch: research-ci/nv4-symmetry-audit-20260810

This directory preserves the P0-3 and P0-4 stage of the 2026-08-10 NV four-orientation symmetry audit.

The branch started from `e5273393fe6fcfb7962c8ac37d014a0075b04126`. P0-3 identified the distinction between a four-axis orbit and a fully completed tetrahedral frame orbit, but the broad PRL novelty claim was killed. P0-4 then corrected the zero-field `ms=±1` labeling subtlety and separated protocol-resolved rank-3 structure from basis-invariant spin-summed rank-4 structure.

These files are retained as research history, not as current claim authority. The cleaned basis-invariant continuation was promoted into the canonical repository layout as:

- `calculations/numerics/nv_symmetry/p0_5_tetrahedral_quartic_odmr_invariant.py`
- `calculations/numerics/nv_symmetry/p0_6_unresolved_rank4_odmr_gate.py`
- `calculations/numerics/nv_symmetry/p0_7_symmetry_breaking_rank_switch.py`
- corresponding evidence under `evidence/witnesses/nv_symmetry_rank/` and `evidence/failures/`

The archived scripts retain their original paths/import assumptions and are preserved verbatim for provenance; they are not part of the current runnable canonical pipeline.
