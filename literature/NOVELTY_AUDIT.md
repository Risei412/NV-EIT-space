# Novelty audit

What is new here, what is a known result in other clothing, and where the
boundary was drawn. The boundary is deliberately **narrow**: several theorems
in this repository are restatements of established results, and are labeled as
such rather than defended.

Primary source: `../evidence/benchmarks/AUDIT_theorem_and_code_2026-08-01.md`
§I-4. Per-paper notes, the claim-to-citation map and the search method are in
`EIT_general_theory_2026-07-13/`.

## Known results in other clothing

| Statement | Established as |
|---|---|
| Theorem 1A (optical dark-subspace rank) | Rank–nullity; Morris–Shore transformation (Morris & Shore 1983, PRA **27**, 906; extended by Rangelov, Vitanov & Shore 2006, PRA **74**, 053402). |
| Theorem 1B (stationary pure Lindblad dark state) | Kraus-form argument; a special case of Albert & Jiang, PRA **89**, 022118 (2014), dark-state structure. |
| Theorem 2A (reduced susceptibility as a Schur complement) | Feshbach–Fano projection / Schur complement. |
| Theorem 2B block formula (as an algebraic identity) | Schur complement. Its *use* — that the sector cut, not `χ_full`, is the object of the no-go — is the new part. |
| Theorem 5 (symmetry-protected transfer zero) | Selection rule from a weak symmetry. |
| `rank[p c] = 1` / nullspace / SVD criteria for an exact dark state | Not new. Standard multilevel EIT theory. |

## Claimed as new

1. **The exclusive trichotomy in `ν`** — `ν ∈ {∞} ∪ (0,∞) ∪ {0}` as an
   exhaustive and mutually exclusive classification for finite-dimensional
   Markovian weak-probe systems, with the finite exponent integer-valued and
   fixed by the first nonvanishing moment.
2. **Identifying the sector cut `χ − χ^(S)_cut` as the correct object of a
   no-go statement.** `χ_full = 0` is not a no-go; a class transition can be
   invisible in the total spectrum and visible only in the cut.
3. **The material-independent classification that follows** — the same three
   tiers in NV, group-IV color centers, and an engineered-loss chain.
4. **The response-theoretic redefinition of EIT** and the conditions
   equivalent to it (absorption cancellation, complex response-zero,
   Schur-complement criterion), together with the dark-state-free 2g+2e
   witness that separates the definition from the dark-state diagnosis.
5. **The matched-readout floor as an identity** — `Re χ_matched = x†Γx ≥ 0`,
   with exact zeros requiring `x ∈ ker Γ`, unreachable under matched readout.

The earlier general-theory framing also identified the defensible core as the
combination `K(ω,T) = p† G_e(ω,T) c` with a controlled merged-manifold
expansion, the `d_p† Π_E d_c` point-group selection rule, and the separation
of exact dark state / finite EIT / zero kernel / observable contrast
(`EIT_general_theory_2026-07-13/00_README.md`).

## What is explicitly not claimed as novelty

- Absence of a result from the author's Zotero library is not evidence that it
  is new (`../NON_CLAIMS.md` N10).
- "Multiple excited levels destroy EIT" as a phenomenology is well covered in
  the atomic-gas literature (e.g. Mishina *et al.* 2011, PRA **83**, 053809)
  and is not claimed here.
- EIT/ATS discrimination as a method is established (Anisimov *et al.*); this
  repository uses it, and does not claim it.
