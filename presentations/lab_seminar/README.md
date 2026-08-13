# Lab seminar deck — NV EIT phase diagram and sign reversal

26 slides, English, 15–20 minutes, built on the Institute of Science Tokyo 4:3
template (layout `コンテンツ1-3`). IMRAD order, with the anticipated-but-secondary
material in an appendix block at the end.

## The five bands

Every content slide places the same five things in the same place, so the eye
does not have to search:

| band | y (in) | what |
| --- | --- | --- |
| head message | 0.11 – 0.87 | interpretive, subject + verb, 20 pt bold navy |
| title line | 0.94 – 1.44 | `Result 2: …` — section and number, 22 pt bold |
| body | 1.62 – 5.28 | figure and/or cards |
| bottom message | 5.38 – 6.18 | cold and factual, 20 pt, in a tinted band |
| citations | 6.42 – 7.07 | hairline, then 10 pt italic sources |

Head and bottom messages are written to fit **two rendered lines**, the title
to fit **one**. The boxes are sized for exactly that, so a longer sentence
overflows rather than shrinking — `qa_geometry.py` reports it.

## Conventions

**Position.** Every content slide's title names its section and its number, so
the audience can place itself without a progress bar.

**Colour, spent sparingly.** Structure is navy `1C3077`. A card's tint is its
meaning and there are only four: `EEF4FF` neutral, `F6F8FC` secondary,
`FFF6E8` a caveat, `FFF1F1` a failure. Red `C00000` is the emphasis of last
resort — it appears eight times in the whole deck, always on the one number
that is the point of the slide.

**Figure text at 18–20 pt.** The manuscript figures cannot be reused: they are
drawn 7 pt in a 7-inch canvas, so on a slide their text lands near 9 pt.
`slide_figures.py` redraws them at the exact inch size they occupy on the
slide, with the fonts set to what they should be there, so nothing is scaled
afterwards. That budget only fits about three ticks an axis, which is why some
slides carry one panel where the paper carries two.

**Citations on the same screen as the claim**, under a hairline at the foot,
numbered `[1]`–`[15]` consistently across the deck. Slides whose numbers are
ours say so rather than citing nothing.

## Files

| file | what it does |
| --- | --- |
| `build.py` | writes `nv_eit_seminar.pptx` from `template.pptx` |
| `slide_figures.py` | writes `figures/sf_*.png` at slide size with 18–20 pt text |
| `qa_geometry.py` | off-canvas shapes, overlaps, text that cannot fit its box, body fonts under 16 pt |
| `preview.py` | approximate PIL raster of each slide (LibreOffice does not run in this environment) |

`slide_figures.py` reads `results/` in the repository root — the same tables and
arrays the manuscript figures use, so a slide number can never drift from the
paper.

## Rebuild

```
python slide_figures.py                        # only when the numbers change
python build.py
python qa_geometry.py nv_eit_seminar.pptx      # expect 0 issues
python preview.py nv_eit_seminar.pptx prev 13 15   # eyeball selected slides
```

`preview.py` draws the layout's placeholder prompt text (「見出し 游ゴシック」 and
friends); PowerPoint does not. Ignore it in the previews.
