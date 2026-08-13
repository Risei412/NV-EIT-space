"""Build the NV-EIT lab seminar deck on the supplied 4:3 corporate template.

Design rules the deck follows (from the requester):
  * IMRAD order
  * every content slide carries a HEAD message (one interpretive sentence,
    subject + verb) and a BOTTOM message (one cold factual sentence); reading
    only those two lines on every slide must give the whole story
  * body text never below 16 pt
  * navy + white, red only for emphasis
"""
from __future__ import annotations

import copy
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

HERE = Path(__file__).resolve().parent
FIG = HERE.parents[1] / "results" / "figures"
OUT = HERE / "nv_eit_seminar.pptx"

NAVY = RGBColor(0x1C, 0x30, 0x77)      # template lt2 -- the centre colour
BLUE = RGBColor(0x7F, 0x96, 0xC2)      # template dk2
TINT = RGBColor(0xEE, 0xF2, 0xFA)      # very light navy, for the bottom band
RED = RGBColor(0xC0, 0x00, 0x00)       # emphasis only
INK = RGBColor(0x1A, 0x1A, 0x1A)
GREY = RGBColor(0x55, 0x55, 0x55)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

FONT = "Arial"

L_TITLE, L_SECTION, L_CONTENT, L_BLANK = 0, 1, 13, 21

SW, SH = 10.0, 7.5
HEAD_Y, HEAD_H = 0.32, 0.80
BODY_TOP, BODY_BOT = 1.32, 6.16
BOT_Y, BOT_H = 6.28, 0.80      # bottom edge 7.08 -- clears the slide number
MARGIN = 0.45


# --------------------------------------------------------------- helpers
def strip_placeholders(slide, keep_idx=()):
    for ph in list(slide.placeholders):
        if ph.placeholder_format.idx not in keep_idx:
            ph._element.getparent().remove(ph._element)


def textbox(slide, x, y, w, h, *, anchor=MSO_ANCHOR.TOP, wrap=True):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = wrap
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0)
    tf.margin_top = tf.margin_bottom = Inches(0.02)
    return tb, tf


def para(tf, first=False):
    return tf.paragraphs[0] if first else tf.add_paragraph()


def run(p, text, *, size=16, bold=False, color=INK, italic=False, font=FONT):
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.italic = italic
    r.font.color.rgb = color
    r.font.name = font
    return r


def rich(p, parts, *, size=16):
    """parts: list of (text, {bold, color, italic})"""
    for text, opt in parts:
        run(p, text, size=opt.get("size", size), bold=opt.get("bold", False),
            color=opt.get("color", INK), italic=opt.get("italic", False))


HEAD_W = 7.00      # stops clear of the template logo at x = 7.59
HEAD_PT = 20       # head and bottom messages carry the deck, so they lead
BOT_PT = 20


def head(slide, text, parts=None):
    """Interpretive one-sentence message across the top."""
    tb, tf = textbox(slide, MARGIN, HEAD_Y, HEAD_W, HEAD_H,
                     anchor=MSO_ANCHOR.MIDDLE)
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    if parts:
        rich(p, parts, size=HEAD_PT)
    else:
        run(p, text, size=HEAD_PT, bold=True, color=NAVY)
    return tb


def bottom(slide, text, parts=None):
    """Cold factual one-sentence message along the base."""
    card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(MARGIN),
                                  Inches(BOT_Y), Inches(SW - 2 * MARGIN),
                                  Inches(BOT_H))
    card.fill.solid()
    card.fill.fore_color.rgb = TINT
    card.line.fill.background()
    card.shadow.inherit = False
    try:
        card.adjustments[0] = 0.12
    except Exception:
        pass
    tf = card.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = tf.margin_right = Inches(0.16)
    tf.margin_top = tf.margin_bottom = Inches(0.04)
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    if parts:
        rich(p, parts, size=BOT_PT)
    else:
        run(p, text, size=BOT_PT, color=INK)
    return card


def picture(slide, name, *, top=None, max_h=None, max_w=None, cx=SW / 2):
    path = FIG / name
    im = Image.open(path)
    ar = im.width / im.height
    mw = max_w if max_w else SW - 2 * MARGIN
    mh = max_h if max_h else (BODY_BOT - (top if top else BODY_TOP))
    w = mw
    h = w / ar
    if h > mh:
        h = mh
        w = h * ar
    x = cx - w / 2
    y = top if top else (BODY_TOP + ((BODY_BOT - BODY_TOP) - h) / 2)
    return slide.shapes.add_picture(str(path), Inches(x), Inches(y),
                                    Inches(w), Inches(h))


def card(slide, x, y, w, h, *, fill=TINT, line=None):
    sh = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x),
                                Inches(y), Inches(w), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = fill
    if line is None:
        sh.line.fill.background()
    else:
        sh.line.color.rgb = line
        sh.line.width = Pt(1)
    sh.shadow.inherit = False
    try:
        sh.adjustments[0] = 0.06
    except Exception:
        pass
    sh.text_frame.word_wrap = True
    sh.text_frame.margin_left = sh.text_frame.margin_right = Inches(0.14)
    sh.text_frame.margin_top = sh.text_frame.margin_bottom = Inches(0.10)
    return sh


def bullets(slide, x, y, w, h, items, *, size=16, space=8):
    tb, tf = textbox(slide, x, y, w, h)
    for i, it in enumerate(items):
        p = para(tf, first=(i == 0))
        p.space_after = Pt(space)
        if isinstance(it, str):
            run(p, "•  " + it, size=size)
        else:
            run(p, "•  ", size=size)
            rich(p, it, size=size)
    return tb


# --------------------------------------------------------------- deck
prs = Presentation(str(HERE / "template.pptx"))
sld_lst = prs.slides._sldIdLst
for sid in list(sld_lst):
    rId = sid.get(
        "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
    prs.part.drop_rel(rId)
    sld_lst.remove(sid)

master = prs.slide_masters[0]


def add(layout_idx, keep=()):
    s = prs.slides.add_slide(master.slide_layouts[layout_idx])
    strip_placeholders(s, keep_idx=keep)
    return s


def section(title, n):
    s = add(L_SECTION, keep=(10,))
    tb, tf = textbox(s, 0.73, 3.35, 8.66, 1.1, anchor=MSO_ANCHOR.MIDDLE)
    p = tf.paragraphs[0]
    run(p, f"{n}.  ", size=28, bold=True, color=BLUE)
    run(p, title, size=28, bold=True, color=NAVY)
    return s


# ---------------------------------------------------------------- 1 title
s = add(L_TITLE)
tb, tf = textbox(s, 0.71, 1.55, 8.81, 1.6, anchor=MSO_ANCHOR.MIDDLE)
p = tf.paragraphs[0]
run(p, "Where NV electromagnetically induced transparency lives,",
    size=28, bold=True, color=NAVY)
p2 = para(tf)
run(p2, "and where it inverts", size=28, bold=True, color=NAVY)
tb, tf = textbox(s, 0.75, 3.60, 6.03, 0.4)
run(tf.paragraphs[0], "Risei Abe", size=18, bold=True, color=INK)
tb, tf = textbox(s, 0.75, 4.05, 6.5, 1.0)
p = tf.paragraphs[0]
run(p, "Department of Electrical and Electronic Engineering", size=16,
    color=GREY)
p = para(tf)
run(p, "Institute of Science Tokyo", size=16, color=GREY)
tb, tf = textbox(s, 5.69, 7.10, 4.15, 0.32)
p = tf.paragraphs[0]
p.alignment = PP_ALIGN.RIGHT
run(p, "Group seminar / progress report", size=16, color=GREY)

# ---------------------------------------------------------------- 2 outline
s = add(L_CONTENT, keep=(0, 11))
s.shapes.title.text_frame.paragraphs[0].runs
tt = s.shapes.title.text_frame
tt.paragraphs[0].runs
for r in list(tt.paragraphs[0].runs):
    r._r.getparent().remove(r._r)
run(tt.paragraphs[0], "Outline", size=22, bold=True, color=NAVY)
strip_placeholders(s, keep_idx=(0, 11))
items = [
    ("Introduction", "What EIT is, why NV, and what we ask"),
    ("Method", "The no-go criterion and how each operating point is judged"),
    ("Results", "Sanity check, the island, room temperature, the field"),
    ("Discussion", "What to measure, and why room temperature fails"),
    ("Summary", "Operating window fixed; experiment this year"),
]
y = 1.48
for i, (t, d) in enumerate(items, 1):
    c = card(s, MARGIN, y, SW - 2 * MARGIN, 0.88)
    tf = c.text_frame
    p = tf.paragraphs[0]
    run(p, f"{i}.  ", size=18, bold=True, color=NAVY)
    run(p, t, size=18, bold=True, color=NAVY)
    p2 = para(tf)
    run(p2, d, size=16, color=INK)
    y += 0.97

# ================================================== INTRODUCTION
section("Introduction", "I")

# ---- I-1 what EIT is
s = add(L_CONTENT)
head(s, 'A control laser makes two excitation paths interfere, turning an absorber transparent.')
s.shapes.add_picture(str(HERE / "lambda_scheme.png"), Inches(0.55),
                     Inches(1.85), Inches(5.05), Inches(3.43))
bullets(s, 5.80, 1.85, 3.75, 4.3, [
    [("Two ground states, one excited state: the ", {}),
     ("Λ", {"bold": True}), (" system.", {})],
    [("The two absorption paths cancel at two-photon resonance.", {})],
    [("Population is trapped in a ", {}),
     ("dark state", {"bold": True, "color": NAVY}),
     (" that the light cannot excite.", {})],
    [("Requires a long-lived ground coherence: it is destroyed by "
      "dephasing.", {})],
])
bottom(s, "In the ideal Λ system probe absorption "
          "vanishes at two-photon resonance when the ground coherence "
          "is long-lived.")

# ---- I-2 why NV
s = add(L_CONTENT)
head(s, 'We turned to NV because the group-IV centres that motivated us have no settled model.')
c1 = card(s, MARGIN, 1.50, 4.35, 3.5, fill=RGBColor(0xF4, 0xF4, 0xF4))
tf = c1.text_frame
run(tf.paragraphs[0], "SiV / SnV — where we started", size=18, bold=True,
    color=GREY)
for t in ["Inversion symmetric: optical line robust against electric noise",
          "Orbital Λ channel, no spin-overlap obstruction",
          "Room-temperature EIT was the hope",
          "But: excited-state Hamiltonian and phonon rates not settled"]:
    p = para(tf); p.space_before = Pt(7)
    run(p, "•  " + t, size=16, color=INK)
c2 = card(s, MARGIN + 4.65, 1.50, 4.35, 3.5, fill=TINT)
tf = c2.text_frame
run(tf.paragraphs[0], "NV — where we can compute", size=18, bold=True,
    color=NAVY)
for t in ["³E Hamiltonian measured and published",
          "Phonon-induced orbital hopping rate Γₓᵧ(T) measured",
          "Cryogenic ensemble EIT already demonstrated",
          "So the theory can be tested against something known"]:
    p = para(tf); p.space_before = Pt(7)
    run(p, "•  " + t, size=16, color=INK)
tb, tf = textbox(s, MARGIN, 5.25, SW - 2 * MARGIN, 0.95)
p = tf.paragraphs[0]
rich(p, [("Strategy: ", {"bold": True, "color": NAVY}),
         ("build the criterion on NV, where every input is known, then carry "
          "the same criterion back to group-IV.", {})], size=16)
bottom(s, "NV wins as a test bed not on material grounds but because its "
          "excited-state model and phonon rates are "
          "independently measured.")

# ---- I-3 objective
s = add(L_CONTENT)
head(s, "We ask whether NV EIT survives to room temperature, and if not, "
        "exactly where it stops.")
q = card(s, MARGIN, 1.55, SW - 2 * MARGIN, 1.35, fill=NAVY)
tf = q.text_frame
tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p = tf.paragraphs[0]
p.alignment = PP_ALIGN.CENTER
run(p, "Is room-temperature NV EIT possible?", size=22, bold=True, color=WHITE)
p = para(tf)
p.alignment = PP_ALIGN.CENTER
run(p, "If not — up to what temperature can it be observed at all?",
    size=18, color=RGBColor(0xD8, 0xE2, 0xF3))
bullets(s, MARGIN, 3.20, SW - 2 * MARGIN, 3.0, [
    [("Known: ", {"bold": True, "color": NAVY}),
     ("NV ensemble EIT is demonstrated at cryogenic temperature.", {})],
    [("Missing: ", {"bold": True, "color": NAVY}),
     ("any microscopic account of why it disappears on warming.", {})],
    [("The usual answer — “dephasing grows” — predicts "
      "only a monotone fade. It predicts no lower bound, no field threshold, "
      "and no change of sign.", {})],
    [("All three occur.", {"bold": True, "color": RED})],
])
bottom(s, "Cryogenic NV ensemble EIT is established; no "
          "published work fixes its upper temperature limit or its "
          "collapse mechanism.")

# ================================================== METHOD
section("Method", "II")

# ---- M-1 no-go / go
s = add(L_CONTENT)
head(s, 'A no-go says the targeted pathway contributes nothing — not that the total response vanishes.')
w1 = card(s, MARGIN, 1.45, 4.35, 1.5, fill=RGBColor(0xFD, 0xEE, 0xEE))
tf = w1.text_frame
run(tf.paragraphs[0], "✗   The tempting definition", size=17, bold=True,
    color=RED)
p = para(tf); p.space_before = Pt(6)
run(p, "“EIT is forbidden when χ_full = 0”", size=16, color=INK)
w2 = card(s, MARGIN + 4.65, 1.45, 4.35, 1.5, fill=TINT)
tf = w2.text_frame
run(tf.paragraphs[0], "✓   The correct object", size=17, bold=True,
    color=NAVY)
p = para(tf); p.space_before = Pt(6)
run(p, "δχₛ ≡ χ_full − χ_cut  ≡  0",
    size=16, color=INK)
tb, tf = textbox(s, MARGIN, 3.12, SW - 2 * MARGIN, 0.62)
rich(tf.paragraphs[0], [
    ("χ_full = 0 is what ", {}),
    ("perfect EIT", {"bold": True, "color": RED}),
    (" looks like.", {})], size=17)
c = card(s, MARGIN, 3.80, SW - 2 * MARGIN, 2.25, fill=RGBColor(0xF7, 0xF9, 0xFC))
tf = c.text_frame
run(tf.paragraphs[0], "Ideal Λ system, worked through", size=17,
    bold=True, color=NAVY)
for t, v, col in [
        ("cut  (pathway removed)", "Ξ_cut = 1", INK),
        ("full (pathway intact)", "Ξ_full = 0    ← response vanishes", INK),
        ("sector correction", "δΞₛ = −1    ← maximally non-zero", RED)]:
    p = para(tf); p.space_before = Pt(7)
    run(p, f"{t:<26s}", size=16, color=GREY, font="Courier New")
    run(p, v, size=16, bold=(col == RED), color=col, font="Courier New")
bottom(s, "The pathway cut is an algebraic counterfactual (B = C = 0 with the "
          "steady state held fixed), not a control-off measurement.")

# ---- M-2 material independence
s = add(L_CONTENT)
head(s, 'Built from the response kernel, the criterion holds for any Markovian emitter.')
steps = [("Material", "H, jump operators,\ndipoles"),
         ("Kernel", "Kᵢⱼ = dᵢ† G(z) dⱼ"),
         ("Test", "δχₛ = 0\n⇔  K₁₂K₂₁ = 0"),
         ("Verdict", "no-go / suppressed /\ngo")]
x = MARGIN
wid = (SW - 2 * MARGIN - 3 * 0.32) / 4
for i, (t, d) in enumerate(steps):
    fill = NAVY if i == 3 else TINT
    fg = WHITE if i == 3 else NAVY
    c = card(s, x, 1.55, wid, 1.75, fill=fill)
    tf = c.text_frame
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    run(p, t, size=18, bold=True, color=fg)
    for line in d.split("\n"):
        p = para(tf); p.alignment = PP_ALIGN.CENTER
        run(p, line, size=16, color=(WHITE if i == 3 else INK))
    if i < 3:
        cn = s.shapes.add_connector(1, Inches(x + wid + 0.05), Inches(2.42),
                                    Inches(x + wid + 0.27), Inches(2.42))
        cn.line.color.rgb = BLUE
        cn.line.width = Pt(2.25)
    x += wid + 0.32
bullets(s, MARGIN, 3.60, SW - 2 * MARGIN, 2.6, [
    [("Nothing in the test is specific to NV, to diamond, or even to "
      "optics.", {})],
    [("Already applied to group-IV centres (SiV, SnV) and to a "
      "non-optical witness.", {})],
    [("What changes between materials is only which pathway is closed and "
      "by what.", {})],
])
bottom(s, "For the scalar block, δχₛ vanishes identically if "
          "and only if the product of the two Raman vertices "
          "K₁₂K₂₁ vanishes.")

# ---- M-3 how we classify
s = add(L_CONTENT)
head(s, 'We judge every point by three criteria at once, because a dip alone identifies nothing.')
crit = [("1. Sign", "Does the pathway reduce or\nincrease absorption?",
         "C > 0  transparency\nC < 0  induced absorption"),
        ("2. Containment", "Is the feature actually inside\nthe scan window?",
         "adaptive two-photon window\n(widths span 3 decades)"),
        ("3. Lineshape", "EIT interference, or\nAutler–Townes doublet?",
         "AIC model comparison\n|ΔAIC| ≥ 6")]
x = MARGIN
wid = (SW - 2 * MARGIN - 2 * 0.30) / 3
for t, q_, a in crit:
    c = card(s, x, 1.50, wid, 2.55, fill=TINT)
    tf = c.text_frame
    run(tf.paragraphs[0], t, size=18, bold=True, color=NAVY)
    for line in q_.split("\n"):
        p = para(tf); p.space_before = Pt(4)
        run(p, line, size=16, color=INK)
    for line in a.split("\n"):
        p = para(tf); p.space_before = Pt(4)
        run(p, line, size=16, color=GREY, italic=True)
    x += wid + 0.30
bullets(s, MARGIN, 4.35, SW - 2 * MARGIN, 1.9, [
    [("Physics from the ", {}),
     ("full 9-level Lindblad generator", {"bold": True, "color": NAVY}),
     (" — no adiabatic elimination.", {})],
    [("240 grid points over T = 20–115 K and B⊥ = 0–0.5 T; each "
      "point gets one of four labels.", {})],
])
bottom(s, "A strong control gives an Autler–Townes "
          "doublet even with zero ground coherence, so a dip alone cannot "
          "establish EIT.")

# ================================================== RESULTS
section("Results", "III")

# ---- R-1 sanity check
s = add(L_CONTENT)
head(s, "The general formula reproduces the textbook three-level result as a "
        "one-line special case.")
c = card(s, MARGIN, 1.45, SW - 2 * MARGIN, 1.9, fill=RGBColor(0xF7, 0xF9, 0xFC))
tf = c.text_frame
run(tf.paragraphs[0], "General (Schur complement, any number of excited "
                      "branches)", size=16, bold=True, color=NAVY)
p = para(tf); p.space_before = Pt(6)
run(p, "Ξ  =  S₁  −  β K₁₂K₂₁ / "
       "(γ_g + β S₂)", size=18, color=INK)
p = para(tf); p.space_before = Pt(10)
rich(p, [("One excited state, scalar dipoles ", {"color": GREY}),
         ("→", {"color": BLUE, "bold": True}),
         ("   Ξ = γ_g / (Aγ_g + β)   = the standard "
          "Λ-system result", {})], size=16)
bullets(s, MARGIN, 3.60, SW - 2 * MARGIN, 2.6, [
    [("Recovers the published three-level susceptibility exactly — no "
      "adiabatic approximation is used anywhere.", {})],
    [("Holds for multi-branch excited manifolds: K₁₂ becomes a "
      "coherent sum over branches.", {})],
    [("Full Liouvillian reproduces the archived 70 K benchmark to ", {}),
     ("0.25 %", {"bold": True, "color": RED}), (".", {})],
])
bottom(s, "Textbook Λ susceptibility recovered as a special case; "
          "70 K contrast C = 0.01387 against the archived 0.01384.")

# ---- R-2 phase diagram
s = add(L_CONTENT)
head(s, 'NV EIT does not simply fade — it occupies a closed island in the temperature–field plane.')
picture(s, "fig3_phase_diagram.png", top=1.35, max_h=3.55)
bullets(s, MARGIN, 5.02, SW - 2 * MARGIN, 1.20, [
    [("Needs ", {}), ("B⊥ ≳ 0.15 T", {"bold": True, "color": RED}),
     (" to open the Raman pathway at all — at zero field the two "
      "Λ legs reach orthogonal spin states.", {})],
    [("Closed ", {}), ("below ≈ 22 K", {"bold": True, "color": RED}),
     (" by an Autler–Townes crossover, and ", {}),
     ("above 90–95 K", {"bold": True, "color": RED}),
     (" by phonon-driven collapse.", {})],
], size=16, space=4)
bottom(s, "148 of 240 grid points are classified genuine transparency, 68 "
          "control-induced absorption, 3 Autler–Townes, 21 unresolved.")

# ---- R-3a room temperature
s = add(L_CONTENT)
head(s, "Room-temperature NV EIT is out of reach by nine orders "
        "of magnitude.")
picture(s, "fig4_contrast_vs_T.png", top=1.42, max_h=4.30, max_w=6.35, cx=3.55)
rows = [("30 K", "0.99", "instant"), ("70 K", "1.4×10⁻²", "1 μs"),
        ("90 K", "1.5×10⁻⁴", "10 ms"),
        ("300 K", "−1.1×10⁻⁹", "never")]
c = card(s, 6.95, 1.55, 2.60, 4.05, fill=TINT)
tf = c.text_frame
run(tf.paragraphs[0], "contrast → time", size=17, bold=True, color=NAVY)
for t, cc, tau in rows:
    p = para(tf); p.space_before = Pt(9)
    col = RED if t == "300 K" else INK
    run(p, f"{t}   ", size=16, bold=True, color=col)
    run(p, f"{cc}", size=16, color=col)
    p2 = para(tf)
    run(p2, f"          {tau}", size=16, italic=True, color=GREY)
bottom(s, "At 300 K the residual contrast is 1.1×10⁻⁹ and of "
          "the wrong sign; no integration time reaches SNR = 5.")

# ---- R-3b sign reversal
s = add(L_CONTENT)
head(s, "Above 103 K the control stops inducing transparency and "
        "starts inducing absorption.")
picture(s, "fig2_spectra.png", top=1.38, max_h=2.35)
c = card(s, MARGIN, 3.95, 4.35, 2.25, fill=TINT)
tf = c.text_frame
run(tf.paragraphs[0], "Why this matters", size=17, bold=True, color=NAVY)
for t in ["No Autler–Townes mechanism gives a sign change",
          "So it directly tests the interference account",
          "Detectable: 1.6 s at 105 K"]:
    p = para(tf); p.space_before = Pt(6)
    run(p, "•  " + t, size=16, color=INK)
c = card(s, MARGIN + 4.65, 3.95, 4.35, 2.25, fill=RGBColor(0xF7, 0xF9, 0xFC))
tf = c.text_frame
run(tf.paragraphs[0], "Where it happens", size=17, bold=True, color=NAVY)
p = para(tf); p.space_before = Pt(8)
run(p, "T_sign = ", size=16, color=INK)
run(p, "103 [97, 108] K", size=18, bold=True, color=RED)
p = para(tf); p.space_before = Pt(6)
run(p, "full Liouvillian, 68 % interval,", size=16, color=GREY)
p = para(tf)
run(p, "80 Monte-Carlo samples", size=16, color=GREY)
bottom(s, "Peak contrast turns negative between 100 and 105 K above 0.15 T; "
          "at 105 K the transmission signal is −1.2×10⁻⁵.")

# ---- R-4 transverse field
s = add(L_CONTENT)
head(s, "The transverse field is a switch that opens the pathway, not a dial "
        "that buys temperature.")
tbl = [("B⊥ (T)", "T₁% (K)", "verdict"),
       ("0 – 0.10", "27 – 33", "no usable window"),
       ("0.15", "70.3", "window opens"),
       ("0.2323", "71.1", "—"),
       ("0.50", "69.7", "no further gain")]
x0, y0 = MARGIN, 1.50
cw = [1.75, 1.65, 2.85]
for r, row in enumerate(tbl):
    xx = x0
    for i, cell in enumerate(row):
        isf = (r == 0)
        c = card(s, xx, y0 + r * 0.72, cw[i], 0.66,
                 fill=(NAVY if isf else (TINT if r % 2 else WHITE)))
        if not isf and r % 2 == 0:
            c.line.color.rgb = RGBColor(0xDD, 0xDD, 0xDD)
            c.line.width = Pt(0.75)
        tf = c.text_frame
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        col = WHITE if isf else (RED if (r in (2, 4) and i == 2) else INK)
        run(p, cell, size=16, bold=isf or (i == 1 and not isf), color=col)
        xx += cw[i] + 0.10
bullets(s, 7.62, 1.50, 1.95, 4.6, [
    [("0 → 0.15 T", {"bold": True, "color": NAVY})],
    [("boundary jumps by 40 K", {})],
    [("0.15 → 0.50 T", {"bold": True, "color": NAVY})],
    [("boundary moves < 2 K", {"bold": True, "color": RED})],
], size=16, space=10)
bottom(s, "The field is required to have a channel at all; once open, the "
          "upper temperature is set by dissipation, not by the field.")

# ================================================== DISCUSSION
section("Discussion", "IV")

# ---- D-1 what to measure
s = add(L_CONTENT)
head(s, 'The calculation says where to look and which single measurement is decisive.')
picture(s, "fig6_observables.png", top=1.38, max_h=2.65)
c1 = card(s, MARGIN, 4.15, 4.35, 2.05, fill=TINT)
tf = c1.text_frame
run(tf.paragraphs[0], "Where to look", size=17, bold=True, color=NAVY)
for t in ["25 – 90 K, B⊥ ≳ 0.15 T",
          "Transmission or fluorescence both work",
          "PL dips where transmission rises"]:
    p = para(tf); p.space_before = Pt(6)
    run(p, "•  " + t, size=16, color=INK)
c2 = card(s, MARGIN + 4.65, 4.15, 4.35, 2.05, fill=RGBColor(0xFD, 0xEE, 0xEE))
tf = c2.text_frame
run(tf.paragraphs[0], "The decisive measurement", size=17, bold=True, color=RED)
for t in ["Follow the fringe through zero at 103 K",
          "A sign flip cannot come from ATS",
          "1.6 s at 105 K, 5.6 s at 110 K"]:
    p = para(tf); p.space_before = Pt(6)
    run(p, "•  " + t, size=16, color=INK)
bottom(s, "Transmission is detectable to 110 K, fluorescence to 105 K; "
          "past 120 K technical noise blocks detection outright.")

# ---- D-2 why RT fails
s = add(L_CONTENT)
head(s, 'Room temperature fails because a thermalising orbital mode destroys the Raman path.')
c = card(s, MARGIN, 1.42, SW - 2 * MARGIN, 1.90, fill=RGBColor(0xF7, 0xF9, 0xFC))
tf = c.text_frame
run(tf.paragraphs[0], "The mechanism", size=17, bold=True, color=NAVY)
p = para(tf); p.space_before = Pt(6)
rich(p, [("A jump inside the excited manifold conserves population but acts "
          "on every optical coherence as ", {}),
         ("pure loss", {"bold": True, "color": RED}),
         (".  Phonon hopping Eₓ ↔ E_y therefore decoheres the "
          "pathway while removing no population.", {})], size=16)
p = para(tf); p.space_before = Pt(4)
run(p, "Rate entering the coherence equation:  Γₓᵧ / 4.",
    size=16, color=GREY)
c1 = card(s, MARGIN, 3.52, 4.35, 2.68, fill=TINT)
tf = c1.text_frame
run(tf.paragraphs[0], "What NV would need", size=17, bold=True, color=NAVY)
for t in ["Freeze out the orbital dynamics, or",
          "Push the branch splitting far above kT",
          "Neither is available in bulk NV at 300 K"]:
    p = para(tf); p.space_before = Pt(7)
    run(p, "•  " + t, size=16, color=INK)
c2 = card(s, MARGIN + 4.65, 3.52, 4.35, 2.68, fill=RGBColor(0xF4, 0xF4, 0xF4))
tf = c2.text_frame
run(tf.paragraphs[0], "Back to group-IV", size=17, bold=True, color=GREY)
for t in ["Orbital Λ: no spin-overlap obstruction",
          "But orbital coherence damped by single phonons",
          "Strain and spin–orbit splitting are the levers"]:
    p = para(tf); p.space_before = Pt(7)
    run(p, "•  " + t, size=16, color=INK)
bottom(s, "The obstruction is a thermalising orbital degree of freedom in the "
          "excited manifold, not the spin dephasing usually blamed.")

# ================================================== SUMMARY
s = add(L_CONTENT)
head(s, 'We have fixed the operating window and found a falsifiable prediction to test this year.')
findings = [
    ("Closed island", "Transparency exists only for 25–90 K and "
                      "B⊥ ≳ 0.15 T — bounded on both "
                      "temperature sides."),
    ("Room temperature: no", "Residual contrast 1.1×10⁻⁹, wrong "
                             "sign, undetectable at any integration time."),
    ("Sign reversal at 103 K", "Control-induced absorption above "
                               "103 [97, 108] K — measurable in seconds."),
]
y = 1.50
for i, (t, d) in enumerate(findings, 1):
    c = card(s, MARGIN, y, SW - 2 * MARGIN, 1.08,
             fill=(TINT if i != 3 else RGBColor(0xFD, 0xEE, 0xEE)))
    tf = c.text_frame
    p = tf.paragraphs[0]
    run(p, f"{i}.  ", size=18, bold=True, color=(NAVY if i != 3 else RED))
    run(p, t, size=18, bold=True, color=(NAVY if i != 3 else RED))
    p2 = para(tf)
    run(p2, d, size=16, color=INK)
    y += 1.20
c = card(s, MARGIN, y + 0.02, SW - 2 * MARGIN, 0.95, fill=NAVY)
tf = c.text_frame
tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p = tf.paragraphs[0]
run(p, "Next:  ", size=18, bold=True, color=RGBColor(0xD8, 0xE2, 0xF3))
run(p, "measure within the year and test the predicted window and the sign "
       "reversal against experiment.", size=18, bold=True, color=WHITE)
bottom(s, "PRA manuscript complete; six figures and all numerical gates "
          "recorded and reproducible from the repository.")

# ================================================== SUPPLEMENTARY
section("Supplementary", "S")

# ---- S-1 theorems
s = add(L_CONTENT)
head(s, "The criterion rests on a small set of theorems, each provable with "
        "first-course linear algebra.")
th = [("1A", "Optical dark subspace has dimension N_g − rank Ω; "
              "basis independent."),
      ("1B", "A pure dark state is stationary iff it is an eigenvector of "
             "every jump operator and of an effective Hamiltonian."),
      ("2A", "The weak-probe response is an exact Schur complement — no "
             "adiabatic elimination."),
      ("2B", "δχₛ = 0 ⇔ K₁₂K₂₁ = 0.  "
             "This, not χ_full = 0, is the no-go."),
      ("3",  "Under fast damping the first non-zero path moment fixes the "
             "suppression power of 1/Γ."),
      ("4",  "Finitely many moments certify an all-orders exact zero "
             "(Cayley–Hamilton)."),
      ("5",  "Sector symmetry forces an identically vanishing transfer."),
      ("6",  "Singular damping leaves an O(1) protected channel in ker D.")]
y = 1.36
for k, d in th:
    tb, tf = textbox(s, MARGIN, y, SW - 2 * MARGIN, 0.58)
    p = tf.paragraphs[0]
    run(p, f"Thm {k:<3s} ", size=16, bold=True, color=NAVY)
    run(p, d, size=16, color=INK)
    y += 0.605
bottom(s, "Full statements and proofs are in the group's theory note; the "
          "manuscript uses 2A, 2B and the rate map only.")

# ---- S-2 model audit
s = add(L_CONTENT)
head(s, 'The reduced kernel reproduces boundary temperatures but not contrast magnitudes.')
bullets(s, MARGIN, 1.45, SW - 2 * MARGIN, 2.5, [
    [("Compared against the full Liouvillian on 108 grid points: ", {}),
     ("47 agree within 10 %", {"bold": True}),
     (", and ", {}), ("13 disagree in sign", {"bold": True, "color": RED}),
     (".", {})],
    [("Failures are structured, not random: below 40 K the reduced kernel "
      "reports transparency where the full model reports induced "
      "absorption; above 90 K no field agrees within 10 %.", {})],
    [("Threshold temperatures still agree — the full-model bands "
      "overlap the reduced ones with medians shifted by ≈ 3 K.", {})],
])
c = card(s, MARGIN, 4.10, SW - 2 * MARGIN, 2.1, fill=TINT)
tf = c.text_frame
run(tf.paragraphs[0], "How the manuscript uses this", size=17, bold=True,
    color=NAVY)
for t in ["Every contrast magnitude quoted comes from the full Liouvillian",
          "The reduced kernel is used only for threshold locations, and those "
          "were re-derived with the full model as a check",
          "Truncation (singlet, hyperfine) shifts contrast by ≤ 30 %, "
          "never the sign"]:
    p = para(tf); p.space_before = Pt(7)
    run(p, "•  " + t, size=16, color=INK)
bottom(s, "Reduced-vs-full agreement is 47/108 points within 10 % with 13 "
          "sign disagreements, concentrated below 40 K and above 90 K.")

# ---- S-3 B-perp exponent
s = add(L_CONTENT)
head(s, 'The quadratic law is confirmed at the warm end of the island, not at the cold end.')
picture(s, "fig5_bperp_scaling.png", top=1.38, max_h=3.30)
bullets(s, MARGIN, 4.78, SW - 2 * MARGIN, 1.48, [
    [("Expected C = C_res + a B⊥²: one symmetry-breaking insertion "
      "per Λ leg.", {})],
    [("85 K: ", {"bold": True, "color": NAVY}),
     ("n = 2.11 ± 0.08", {"bold": True}),
     (", stable against the fitting cutoff.  55 / 70 K: the exponent drifts "
      "monotonically and is not determined.", {})],
    [("Independent of control power (×10) and identical for the ratio "
      "and the absolute signal — so the drift is the pathway, not the "
      "probe.", {})],
], size=16, space=3)
bottom(s, "Between the residual floor and mixing saturation at "
          "γₑB = D_gs = 0.103 T there is under a decade of field.")

# ---- S-4 EIT vs ATS
s = add(L_CONTENT)
head(s, "We separate EIT from Autler–Townes by model comparison rather "
        "than by eye.")
c = card(s, MARGIN, 1.45, 4.35, 2.4, fill=TINT)
tf = c.text_frame
run(tf.paragraphs[0], "EIT model (interference)", size=17, bold=True,
    color=NAVY)
p = para(tf); p.space_before = Pt(8)
run(p, "A = C₊²/(γ₊²+δ²) − "
       "C₋²/(γ₋²+δ²)", size=16, color=INK)
p = para(tf); p.space_before = Pt(8)
run(p, "a narrow window inside a broad line", size=16, italic=True, color=GREY)
c = card(s, MARGIN + 4.65, 1.45, 4.35, 2.4, fill=RGBColor(0xF4, 0xF4, 0xF4))
tf = c.text_frame
run(tf.paragraphs[0], "ATS model (doublet)", size=17, bold=True, color=GREY)
p = para(tf); p.space_before = Pt(8)
run(p, "A = C²[1/(γ²+(δ−δ₀)²) + "
       "1/(γ²+(δ+δ₀)²)]", size=16, color=INK)
p = para(tf); p.space_before = Pt(8)
run(p, "two dressed-state peaks, no interference", size=16, italic=True,
    color=GREY)
bullets(s, MARGIN, 4.05, SW - 2 * MARGIN, 2.2, [
    [("Both fitted to the same computed spectrum, together with Lorentzian "
      "and Fano alternatives.", {})],
    [("Verdict from ΔAIC = IC_ATS − IC_EIT with a fixed robust gate "
      "at |ΔAIC| ≥ 6.", {})],
    [("Verdict unchanged under 60 noise bootstraps, 30 randomised fit "
      "initialisations, and window changes of 0.5× to 2×.", {})],
])
bottom(s, "Following Anisimov, Dowling and Sanders, PRL 107, 163604 (2011); "
          "the classification is stable at every temperature tested.")

# ---- S-5 detection chain
s = add(L_CONTENT)
head(s, 'One detection chain at every temperature makes the comparison measure physics, not geometry.')
bullets(s, MARGIN, 1.45, SW - 2 * MARGIN, 2.0, [
    [("The sector cross-section grows as the optical line narrows: a fixed "
      "sample is OD ≈ 16 at 10 K and OD ≈ 0.04 at 300 K.", {})],
    [("Primary axis is therefore the ", {}),
     ("OD-matched sample", {"bold": True, "color": NAVY}),
     (" (sector OD = 1); the fixed 1 ppm / 0.5 mm sample is reported "
      "alongside.", {})],
])
par = [("λ", "637 nm"), ("Debye–Waller", "0.030"),
       ("γ_inh", "30 GHz"), ("probe power", "1 μW"),
       ("η detect", "0.1"), ("σ_tech", "10⁻⁶"),
       ("target SNR", "5"), ("ceiling", "24 h")]
x, y = MARGIN, 3.50
for i, (k, v) in enumerate(par):
    c = card(s, x, y, 2.22, 0.86, fill=TINT)
    tf = c.text_frame
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]
    run(p, k + "   ", size=16, color=GREY)
    run(p, v, size=16, bold=True, color=NAVY)
    x += 2.32
    if (i + 1) % 4 == 0:
        x = MARGIN
        y += 0.98
bottom(s, "Where the archived conversions disagreed (70 K and 300 K), the "
          "conservative value was taken and the disagreement recorded.")

prs.save(str(OUT))
print(f"wrote {OUT}  ({len(prs.slides.__iter__.__self__._sldIdLst)} slides)")
