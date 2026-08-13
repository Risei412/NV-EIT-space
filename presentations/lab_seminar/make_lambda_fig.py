"""Generic three-level Lambda schematic for the introduction slide.

Drawn with matplotlib rather than with PowerPoint shapes: this environment
cannot render a .pptx, so hand-placed shapes and their labels cannot be
checked visually, whereas a figure produced here can be looked at directly.
Styling follows the deck (template navy, light blue, red emphasis).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAVY = "#1C3077"
BLUE = "#7F96C2"
RED = "#C00000"
GREY = "#555555"

fig, ax = plt.subplots(figsize=(5.0, 3.4))
ax.set_xlim(0, 10)
ax.set_ylim(0.6, 9.6)
ax.axis("off")

# levels
ax.plot([1.1, 3.1], [2.2, 2.2], color=NAVY, lw=3, solid_capstyle="butt")
ax.plot([6.9, 8.9], [3.0, 3.0], color=NAVY, lw=3, solid_capstyle="butt")
ax.plot([4.0, 6.0], [8.4, 8.4], color=NAVY, lw=3, solid_capstyle="butt")

ax.text(0.95, 2.2, r"$|g_1\rangle$", ha="right", va="center",
        fontsize=15, color=NAVY)
ax.text(9.05, 3.0, r"$|g_2\rangle$", ha="left", va="center",
        fontsize=15, color=NAVY)
ax.text(5.0, 8.85, r"$|e\rangle$", ha="center", va="bottom",
        fontsize=15, color=NAVY)

# drive legs
ax.annotate("", xy=(4.45, 8.3), xytext=(2.35, 2.35),
            arrowprops=dict(arrowstyle="<->", color=NAVY, lw=2.2))
ax.annotate("", xy=(5.6, 8.3), xytext=(7.7, 3.15),
            arrowprops=dict(arrowstyle="<->", color=BLUE, lw=2.2))
ax.text(2.55, 5.7, r"probe  $\Omega_p$", fontsize=14, color=NAVY,
        rotation=68, ha="center", va="center")
ax.text(7.55, 5.9, r"control  $\Omega_c$", fontsize=14, color=BLUE,
        rotation=-68, ha="center", va="center")

# the ground coherence that carries the interference
ax.annotate("", xy=(6.85, 2.85), xytext=(3.15, 2.35),
            arrowprops=dict(arrowstyle="<->", color=RED, lw=1.8,
                            linestyle=(0, (4, 3))))
ax.text(5.0, 1.62, "ground coherence\n(must be long-lived)", fontsize=13,
        color=RED, ha="center", va="top")

fig.tight_layout(pad=0.2)
out = Path(__file__).resolve().parent / "lambda_scheme.png"
fig.savefig(out, dpi=300, facecolor="white")
print("wrote", out)
