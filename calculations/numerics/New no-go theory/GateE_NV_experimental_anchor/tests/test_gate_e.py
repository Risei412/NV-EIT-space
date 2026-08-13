"""Regression checks for Gate E P0-1 observable freeze."""
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

points = mod.load_ensemble_points()
rows = mod.window_table(points)
gating = {
    row.observable: row
    for row in rows
    if row.scenario == mod.GATING_SCENARIO
}
raw = gating["signed_absorption_difference"]
contrast = gating["normalized_contrast"]

assert raw.order == 4
assert contrast.order == 3
assert abs(raw.margin_over_floor - contrast.margin_over_floor) < 1e-9
assert raw.window_decades < 1.0
assert contrast.window_decades > 1.0
assert 4.3 < raw.improvement_to_one_decade < 4.5
assert contrast.improvement_to_one_decade == 1.0

for row in (raw, contrast):
    constrained = mod.monte_carlo_identifiability(
        row.reference_value,
        row.detection_floor,
        row.order,
        row.window_decades,
        True,
    )
    unconstrained = mod.monte_carlo_identifiability(
        row.reference_value,
        row.detection_floor,
        row.order,
        row.window_decades,
        False,
    )
    assert constrained["correct_class_probability"] > 0.95
    assert unconstrained["correct_class_probability"] < 0.95

print("Gate E P0-1 regression checks: PASS")
