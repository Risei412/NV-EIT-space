# Lab seminar deck — NV EIT phase diagram and sign reversal

26 slides, English, 15–20 minutes, built on the Institute of Science Tokyo 4:3
template. IMRAD order, with the anticipated-but-secondary material moved to a
supplementary block at the end.

## Design contract

Every content slide carries two sentences that must, on their own, tell the
whole story:

* a **head message** at the top — interpretive, subject + verb, 20 pt bold navy;
* a **bottom message** in the tinted band at the base — cold and factual, 20 pt.

Body text never drops below 16 pt. Navy and white throughout; red is reserved
for a single emphasised value per slide.

Head and bottom messages are written to fit **two rendered lines**. The boxes
are sized for exactly that (`HEAD_H`, `BOT_H` in `build.py`), so a longer
sentence overflows rather than shrinking — `qa_geometry.py` reports it.

## Files

| file | what it does |
| --- | --- |
| `build.py` | writes `nv_eit_seminar.pptx` from `template.pptx` |
| `qa_geometry.py` | geometric QA: off-canvas shapes, overlaps, text that cannot fit its box, fonts under 16 pt |
| `preview.py` | approximate PIL raster of each slide (LibreOffice does not run in this environment) |
| `make_lambda_fig.py` | draws `lambda_scheme.png`, the schematic on the first content slide |

Figures come from `results/figures/`, produced by
`calculations/numerics/manuscript_figures/pra/src/make_figures.py`.

## Rebuild

```
python build.py
python qa_geometry.py nv_eit_seminar.pptx      # expect 2 known template hits
python preview.py nv_eit_seminar.pptx prev 13 15   # eyeball selected slides
```

`qa_geometry.py` reports two issues that are the template's own geometry, not
ours: the title-slide footer at y = 7.10 and the layout's title placeholder on
the outline slide.
