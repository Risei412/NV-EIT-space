# Prior-art boundary

Where this work stops and the literature begins, stated per topic so a
reviewer can check the line rather than infer it. The per-paper notes,
including the "competition / novelty risk" field for each entry, are in
`EIT_general_theory_2026-07-13/paper_notes/`; the claim-to-citation map is
`EIT_general_theory_2026-07-13/05_claim_to_citation_map.md`.

## Boundary by topic

| Topic | Prior art (theirs) | Boundary (ours) |
|---|---|---|
| Standard 3-level EIT, linear susceptibility, dark-state polaritons | Fleischhauer, Imamoglu & Marangos, RMP **77**, 633 (2005) | Used as the notation and convention baseline. We extend which assumptions, not which formulas. |
| Bright/dark decomposition of a degenerate coupling matrix | Morris & Shore (1983); Rangelov, Vitanov & Shore (2006) | Enters as a *lemma*, not as a theorem of ours. Neither handles dissipation, frequency-dependent Green functions, or thermally merged manifolds. |
| Dark-state classification in dissipative multilevel systems | "Classification of dark states in multilevel dissipative systems"; "Simple determination of dark states"; "General dark-state theory for arbitrary multilevel quantum systems" (C05–C07) | These classify *states*. We classify *response*: the sector-mediated susceptibility and its suppression exponent. A system with no dark state can still be transparent (see the 2g2e witness). |
| Multiple excited levels + inhomogeneous broadening degrade EIT | Mishina *et al.*, PRA **83**, 053809 (2011); C09 | Their result is phenomenological degradation. Ours is a structural statement about *which* class the degradation lands in, and an integer exponent. |
| EIT vs. ATS discrimination | Anisimov *et al.* (C10) | Used as a tool inside the gates; not claimed. |
| Temporal buildup of EIT/EIA | C08 | Outside our weak-probe stationary scope; cited as an adjacent regime, not competed with. |
| Stationary dark states of Lindblad dynamics | Albert & Jiang, PRA **89**, 022118 (2014) | Theorem 1B is a special case of theirs and is labeled as such. |
| NV / group-IV microscopic structure, phonons, strain, temperature dependence | S01–S14 (reinforcement set) | Inputs to our models. We do not claim NV spectroscopy results; we claim what the response classification does with them. |
| Applications: CPT and EIT in SiV/SnV/PbV, SiC, rare-earth, quantum dots, hBN | A01–A17 (application set) | Targets of the classification, not results of it. |

## Rules for maintaining this boundary

- Every claim in `../CLAIMS.md` must be traceable to a row above, or must add
  a new row.
- A statement that a competitor's result "does not apply" needs the specific
  assumption that breaks, named in `../theory/ASSUMPTIONS.md`.
- New references go in `references.bib`, with a note in
  `EIT_general_theory_2026-07-13/paper_notes/` recording the novelty risk.
