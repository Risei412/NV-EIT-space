"""Regression checks for P0-2 two-layer classification."""
from __future__ import annotations

import csv
import importlib.util
import json
import pathlib
import sys

ROOT = (pathlib.Path(__file__).resolve().parents[2]
        / "numerics" / "manuscript_figures" / "pra")
SCRIPT = ROOT / "src" / "p0_2_two_layer_classification.py"
spec = importlib.util.spec_from_file_location("p0_2", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
assert spec.loader is not None
spec.loader.exec_module(mod)


def fake_fit(best: str, runner: str, margin: float, d_ats_eit: float):
    order = [best, runner] + [m for m in ("EIT", "ATS", "Fano", "Lorentzian")
                              if m not in (best, runner)]
    values = {order[0]: 0.0, order[1]: margin, order[2]: margin + 10.0,
              order[3]: margin + 20.0}
    models = {
        name: {"aic": values[name], "aicc": values[name],
               "akaike_weight": 1.0 if name == best else 0.0}
        for name in values
    }
    return {
        "models": models,
        "used_aicc": False,
        "delta_aic_ats_eit": d_ats_eit,
        "verdict": ("robust EIT" if d_ats_eit >= mod.ROBUST else
                    "robust ATS" if d_ats_eit <= -mod.ROBUST else "inconclusive"),
    }


# EIT must be best over all four models and robust against the runner-up.
r = mod.classify_spectrum(fake_fit("EIT", "Fano", 7.0, 10.0))
assert r["spectral_class"] == "EIT"
r = mod.classify_spectrum(fake_fit("Fano", "EIT", 7.0, 10.0))
assert r["spectral_class"] == "Fano"
r = mod.classify_spectrum(fake_fit("EIT", "Fano", 2.0, 10.0))
assert r["spectral_class"] == "spectral_ambiguous"

# Sector sign and spectroscopy are deliberately independent.
assert mod.joint_class("sector_transparency", "Fano") == \
    "Fano_shaped_sector_transparency"
assert mod.joint_class("sector_transparency", "EIT") == "spectroscopic_EIT"
assert mod.joint_class("sector_absorption", "EIT") == "control_induced_absorption"

# Frozen full-grid result must remain internally consistent.
tab = pathlib.Path(__file__).resolve().parents[2].parent / "results" / "tables"
summary = json.loads((tab / "p0_2_two_layer_summary.json").read_text(encoding="utf-8"))
assert summary["verdict"] == "PASS"
assert all(summary["gates"].values())
assert summary["grid"]["n_points"] == 240
assert sum(summary["sector_counts"].values()) == 240
assert sum(summary["spectral_counts"].values()) == 240
assert sum(summary["joint_counts"].values()) == 240
assert summary["joint_counts"]["spectroscopic_EIT"] == 9
assert summary["legacy_classifier_audit"]["old_transparency_with_nonEIT_best_model"] == 139
assert summary["candidate_70K"]["joint_class"] == "Fano_shaped_sector_transparency"
assert summary["candidate_70K"]["best_model"] == "Fano"
assert not summary["Bperp_zero_audit"]["has_robust_spectroscopic_EIT"]

# Compact grid encodes exactly 20 x 12 classified points.
grid_path = tab / "p0_2_joint_class_grid.csv"
with grid_path.open(encoding="utf-8") as fh:
    comment = fh.readline().strip()
    rows = list(csv.DictReader(fh))
assert comment.startswith("# E=spectroscopic_EIT")
assert len(rows) == 20
assert all(len(row) == 13 for row in rows)
codes = [value for row in rows for key, value in row.items() if key != "T_K"]
assert len(codes) == 240
assert codes.count("E") == 9
assert codes.count("F") == 154
assert codes.count("A") == 3
assert codes.count("?") == 3
assert codes.count("N") == 68
assert codes.count("U") == 3

# A freshly generated full CSV, when present, must obey the same no-mixing rule.
csv_path = tab / "p0_2_two_layer_phase_diagram.csv"
if csv_path.exists():
    with csv_path.open(encoding="utf-8") as fh:
        full_rows = list(csv.DictReader(fh))
    assert len(full_rows) == 240
    assert not any(r["joint_class"] == "spectroscopic_EIT"
                   and r["best_model"] != "EIT" for r in full_rows)
    assert not any(r["joint_class"] == "spectroscopic_EIT"
                   and float(r["best_runner_margin"]) < mod.ROBUST for r in full_rows)

print("P0-2 two-layer classification regression checks: PASS")
