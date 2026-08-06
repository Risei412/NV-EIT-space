"""Regression checks for Gate E."""
import importlib.util
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = HERE / "src" / "run_gate_e.py"
spec = importlib.util.spec_from_file_location("gate_e", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
assert spec.loader is not None
spec.loader.exec_module(mod)

rows = mod.window_table()
gating = next(r for r in rows if r.scenario == mod.GATING_SCENARIO)
assert gating.intensity_window_decades < 1.0
assert gating.amplitude_window_decades > 1.0
assert 4.3 < gating.improvement_to_one_decade < 4.5

constrained = mod.monte_carlo_identifiability(
    gating.contrast, gating.intensity_window_decades, True
)
unconstrained = mod.monte_carlo_identifiability(
    gating.contrast, gating.intensity_window_decades, False
)
assert constrained["correct_class_probability"] > 0.95
assert unconstrained["correct_class_probability"] < 0.95
print("Gate E regression checks: PASS")
