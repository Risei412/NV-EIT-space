# LIMITATIONS.md

This file lists the current boundaries of the canonical v7 theory package. The detailed pre-v7 limitation record remains recoverable from Git history before the v7 migration.

## Current theory boundaries

- The main asymptotic results assume a finite-dimensional weak-probe response with fixed scaling `A_Gamma = Gamma D + A_0`.
- A second parameter scale that grows with `Gamma` is outside that fixed-scaling theorem.
- The sector-graph distance is a lower bound on the response order. Equality requires the leading path sum to remain nonzero.
- The reported exponent belongs to a named observable. A later normalization can change that exponent.
- The singular-damping result in v7 is proved for Hermitian positive-semidefinite `D`; a general nonnormal extension is not included.
- The T1 sector-cut package still has independent documentation gaps: its Theorem 2.1 is not proved in this repository, its Theorem 2.2 numerical example has no local generating code, and the genericity language in Theorem 2.3(ii) remains undefined.
- The sector cut depends on the physical choice of sector and is therefore not unique across different declared cuts.

## Experimental boundaries retained from earlier audits

- Optical NV does not provide a clean experimental window for resolving its own asymptotic exponent under the audited measurement budget.
- The Gate B device witness supports structural behavior but not a practical exponent-measurement protocol in its tested physical window.
- The broad measurable-design-rule claim remains a non-claim; see `../NON_CLAIMS.md`.

## v7 closure note

The pre-v7 source contained an omitted algebraic step in the stationary pure-state argument. The canonical `proofs/eit_nogo_proofs_v7.tex` replaces that passage with an explicit calculation, so that item is not an active v7 proof gap.

## Artifact notes

The regenerated F5B certificate supersedes its stale predecessor. Manuscript Fig. 4 still needs regeneration with the documented non-quick setting.
