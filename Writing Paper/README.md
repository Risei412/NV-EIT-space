# Writing Paper

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

- `prl_figures/` — those four figures, assembled from the Gate A–D campaigns.
  Presentation only: every panel reuses the gate functions unchanged, with no
  new physics. See `prl_figures/README.md` for panel-by-panel description,
  data sources, and how to regenerate.

- `drafts/` — manuscript sources.
  `The theory of EIT in Nitrogen Vacancy.tex` is the working LaTeX draft.

## Regenerating the figures

```bash
python prl_figures/src/make_figures.py            # full
python prl_figures/src/make_figures.py --quick    # fast sweep
python prl_figures/tests/test_figures.py          # regression tests
```

Outputs land in `prl_figures/figures/` as both PDF (vector, for submission)
and PNG (300 dpi, preview). The regression tests check the figures against the
certified Gate A–D numbers, so a figure cannot silently drift from the result
it depicts.
