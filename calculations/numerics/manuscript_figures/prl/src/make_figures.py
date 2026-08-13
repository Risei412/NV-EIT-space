"""make_figures.py -- generate all four PRL main-text figures.

Usage:
    python make_figures.py [--quick]

Outputs (figures/): fig1_classes, fig2_inheritance, fig3_material_independence,
fig4_robustness, each as .pdf (vector, for submission) and .png (dpi=300, preview).

No new physics: every figure module reuses the Gate A-D campaign code
unchanged (see each fig*.py docstring for its data sources).
"""
from __future__ import annotations

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import fig1_classes
import fig2_inheritance
import fig3_material_independence
import fig4_robustness


def main():
    ap = argparse.ArgumentParser(description="Generate PRL main-text figures")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    outputs = []
    for mod, name in [(fig1_classes, "Fig. 1 (classes)"),
                      (fig2_inheritance, "Fig. 2 (inheritance)"),
                      (fig3_material_independence, "Fig. 3 (material independence)"),
                      (fig4_robustness, "Fig. 4 (robustness)")]:
        t1 = time.time()
        pdf, png = mod.build(quick=args.quick)
        print(f"{name}: {pdf}  ({time.time()-t1:.2f}s)")
        outputs.append((pdf, png))

    print(f"\n{len(outputs)} figures written in {time.time()-t0:.2f}s")
    return outputs


if __name__ == "__main__":
    main()
