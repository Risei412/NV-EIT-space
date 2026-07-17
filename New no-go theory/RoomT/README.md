# Room-temperature SMRT no-go campaign (NV ZPL spin-Lambda)

Source strategy: `SMRT_NV_room_temperature_EIT_no_go_numerical_plan.md`
(user-provided, not committed here), which defines Steps 1-9, the four
EIT criteria (E1-E4), the SMRT class certificate, and the figure/table
plan for the room-temperature (300 K) no-go claim. This directory turns
that plan into an executable numerical campaign on top of the existing
infrastructure in `No-go theorem/` and `New no-go theory/`.

## Status

**Step 1 (positive control): EXECUTED, all 6 gates PASS.**
See `results/gates_summary_step1.json` and
`results/figures/fig_step1_low_T_positive_control.png`.
Code: `src/step1_low_temperature_validation.py`.

Reuses two already-validated models, no new physics:

- `No-go theorem/src/nv_model.py` -- the physical NV branch-resolved
  spin-Lambda (ground |ms=0>, |ms=-1> connected through the thermally-
  split 3E orbital branches X/Y), with `gamma_oc_GHz(T,d)` the phonon-
  driven orbital-hopping rate (Happacher et al., `phonon_rates.py`).
  Used for: nonzero contrast at low T, the sector-cut null test (Oc->0
  kills the feature identically), the qualitative low-T temperature
  trend, and the ground-decoherence trend.
- `New no-go theory/src/model_lambda.py` -- the clean three-level Lambda
  (Phase A positive control, plan Sec. 3.2 item 1), used only for the
  "EIT linewidth broadens with control power" check, because
  `nv_model.response()` ties probe and control to a single shared
  detuning (no independent two-photon/Raman-detuning axis to scan) --
  see the code comment in `step1_low_temperature_validation.py` for why
  that check cannot be done directly on `nv_model.py`.

Gates certified (`gates_summary_step1.json`):

| Gate | Requirement (plan Sec. 5, Step 1) | Result |
|---|---|---|
| `nonzero_at_low_T` | R_EIT != 0 at T=4 K | PASS (\|C\|=0.89) |
| `cut_kills_feature` | sector cut (Oc=0) removes the window exactly | PASS (\|C_cut\|<1e-12 at every T) |
| `monotone_decreasing_in_T` | qualitative low-T suppression trend | PASS (4->30 K) |
| `width_increasing_with_control_power` | EIT FWHM broadens with Omega_c | PASS (~Omega_c^2 scaling, 0.53->842 GHz over Omega_c in [0.2,8] GHz) |
| `no_window_edge_clipping` | broadening is genuine, not a fixed-window artifact | PASS (window scaled adaptively per Omega_c) |
| `monotone_decreasing_in_ground_decoherence` | contrast falls as gamma_ground grows | PASS |

Two real bugs were caught and fixed during this run (see git history):
an absorption-quadrature sign error (used Im(chi) instead of Re(chi) for
the Lorentzian absorption line in `model_lambda`'s 1/(gamma-i*Delta)
convention, which silently zeroed the linewidth-vs-power check via a
divide-by-zero), and a fixed-width detuning window that clipped the
FWHM measurement at large Omega_c and made the broadening look like it
saturated. Both are guarded against regressing via the
`no_window_edge_clipping` gate and by scanning `Omega_c` over almost two
decades.

Gate 1 (plan: "Gate 1を通過しないcodeで室温no-goを論じてはならない")
is therefore satisfied; Step 2 (operational cut audit on the full NV
Liouvillian, matching `New no-go theory/Sector/src/operational_cut.py`'s
GKSL-admissible D_S construction) is next.

## Layout

    src/step1_low_temperature_validation.py   Step 1 (this run)
    results/gates_summary_step1.json
    results/figures/fig_step1_low_T_positive_control.png

Future steps (2-9) will follow the same convention: one `run_stepN_*.py`
per step, one `gates_summary_stepN.json`, figures under `results/figures/`.
