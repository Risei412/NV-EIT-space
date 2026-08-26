# Submission notes

Manuscript work: how the results in this repository are split across papers,
the figures assembled for submission, and the drafts themselves.

## Contents

- `NV_EIT_PRA_PRL_Split_Strategy_20260724.md` — the strategy document. It
  argues that the current results are clearest as two papers with different
  questions:
  - **PRA** — the microscopic mechanism by which low-temperature EIT in NV
    centers collapses as temperature rises: the temperature window where it
    holds, the transverse-field rescue pathway, and the room-temperature
    no-go, quantified.
  - **PRL** — EIT demoted to one instance of a general response law: for
    linear response mediated by a fast-dissipative subspace, selection-rule
    path cancellation selects an *integer* dissipative suppression order,
    across materials and across phenomena.

  The document also fixes what each paper must show before submission, and
  §17 specifies the compression of the PRL main text to four figures.

- `../calculations/numerics/manuscript_figures/prl/` — those four figures, assembled from the Gate A–D campaigns.
  Presentation only: every panel reuses the gate functions unchanged, with no
  new physics. See `../calculations/numerics/manuscript_figures/prl/README.md` for panel-by-panel description,
  data sources, and how to regenerate.

- `drafts/` — manuscript sources.
  `The theory of EIT in Nitrogen Vacancy.tex` is the working LaTeX draft.

## Regenerating the figures

```bash
python ../calculations/numerics/manuscript_figures/prl/src/make_figures.py          # full
python ../calculations/numerics/manuscript_figures/prl/src/make_figures.py --quick  # fast sweep
pytest ../calculations/tests/prl_figures                                            # regression tests
```

Outputs land in `../results/figures/` as both PDF (vector, for submission)
and PNG (300 dpi, preview). The regression tests check the figures against the
certified Gate A–D numbers, so a figure cannot silently drift from the result
it depicts.


## Layout after the 2026-08-13 restructure

- `main.tex` / `main.pdf` — the PRA main text (formerly `nv_eit_pra.tex`).
  `\graphicspath{{../results/figures/}}`, so the figures it includes are the
  generated ones under `results/figures/` — build from this directory.
- `references.bib` — the manuscript bibliography. The literature package's
  master bibliography is separate, at `../literature/references.bib`.
- `drafts/` — working LaTeX drafts.
- `NV_EIT_PRA_PRL_Split_Strategy_20260724.md` — the PRA/PRL split strategy.
- `figures/` — see that directory's README: manuscript figures are not
  duplicated here, they are read from `../results/figures/`.

## Closed defect

**Fig. 4 was a `quick=True` figure** while the instructions above say
`quick=False`, which prints different values in the figure. Its numbers were
not to be quoted as full-run numbers until it was regenerated
(`../NON_CLAIMS.md` N8).

**Closed 2026-08-26.** `results/figures/fig4_robustness.{pdf,png}` were
regenerated with `make_figures.py` and no `--quick`. The inset label moved from
`-1.01` to `-1.00`, and both printed quantities now equal the Gate D
certificate exactly:

| quantity | figure (`quick=False`) | `gates_summary_gateD.json` |
|---|---|---|
| `crossover_power` | `-0.9971223021582747` | `-0.9971223021582747` |
| `nu_protected_eps0` | `1.9989602590799627` | `1.9989602590799627` |

Figures 1-3 were regenerated in the same run and are byte-identical to the
committed versions apart from the PDF `/CreationDate`, so they were left
untouched.

**Root cause.** Regenerating the figure alone would not have closed this. The
test suite itself rebuilt Fig. 4: `test_fig4_outputs_and_crossover` called
`fig4_robustness.build(quick=True)`, and `build()` writes into
`results/figures/`, so every `pytest` run silently overwrote the shipped
figure with the 90-sample version. That call is now `quick=False`, matching
`make_figures.py`, with a comment saying why it must not be changed back.

**Observed while closing this, not fixed here.** `build()` writes a fresh
`/CreationDate` into every PDF, so `results/figures/*.pdf` change on each run
even when the rendered content is identical (the PNGs are byte-stable). The
root README's claim that the figure regression tests reproduce the committed
figures byte-identically therefore holds for the PNGs but not the PDFs.
Setting `SOURCE_DATE_EPOCH` before building would make the PDFs deterministic.

The regression that would have caught this is now in place:
`calculations/tests/prl_figures/test_figures.py` runs
`test_fig4a_crossover_power_label` at `quick=False` (the setting the committed
figure is built with) and pins the printed label exactly, and
`test_fig4a_matches_gate_d_certificate` asserts the figure's two numbers equal
the certificate's.
