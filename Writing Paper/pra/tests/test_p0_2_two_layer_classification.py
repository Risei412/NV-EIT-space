"""Regression checks for P0-2 two-layer classification."""
from __future__ import annotations

import csv
import importlib.util
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
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

summary_path = ROOT / "results" / "tables" / "p0_2_two_layer_summary.json"
csv_path = ROOT / "results" / "tables" / "p0_2_two_layer_phase_diagram.csv"
if summary_path.exists() and csv_path.exists():
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["verdict"] == "PASS"
    assert summary["new_mixed_EIT_nonEIT_best_count"] == 0
    assert all(summary["gates"].values())
    with csv_path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    assert not any(r["joint_class"] == "spectroscopic_EIT"
                   and r["best_model"] != "EIT" for r in rows)
    assert not any(r["joint_class"] == "spectroscopic_EIT"
                   and float(r["best_runner_margin"]) < mod.ROBUST for r in rows)

print("P0-2 two-layer classification regression checks: PASS")
