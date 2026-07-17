# Room-temperature SMRT no-go campaign (NV ZPL spin-Lambda)

Source strategy: `SMRT_NV_room_temperature_EIT_no_go_numerical_plan.md`
(user-provided, not committed here), which defines Steps 1-9, the four
EIT criteria (E1-E4), the SMRT class certificate, and the figure/table
plan for the room-temperature (300 K) no-go claim. This directory turns
that plan into an executable numerical campaign on top of the existing
infrastructure in `No-go theorem/` and `New no-go theory/`.

## Status

**Step 2 (operational cut audit, full NV Liouvillian): EXECUTED, all 4 gates PASS.**
See `results/gates_summary_step2.json` and
`results/figures/fig_step2_operational_cut_audit.png`.
Code: `src/step2_operational_cut_audit.py`.

Builds the FULL (N=9, dim=81) NV Liouvillian at Step 1's exact geometry
(T=10 K, Bz=0.02 T, Bx=0, control on ms=-1, probe/control on orbital
branches X/Y), reusing `No-go theorem/src/gate2_candidate_full_vs_reduced.py`'s
`build_full()` for the Lindbladian construction and
`New no-go theory/Sector/src/operational_cut.py` for the operational-cut
machinery already validated on toy models in
`New no-go theory/Sector/src/run_gates_step3.py` (Gates U1/U2/U5). The
sector cut S = {rho_pc, rho_cp} is the same one `tests/
test_operational_cut_equivalence.py` proved is numerically identical to
gate2's own ad hoc "zero the S<->X Liouvillian blocks" construction.

Since the full vectorized Liouvillian is exactly singular (one
trace-preserving steady state), the response operator is regularized as
A(z) = i*z*I - L at a small representative z (0.001 GHz, chosen well
above L's spectral gap of ~3e-6 GHz) -- a technical device to get an
invertible matrix for the operational-cut formulas, cross-validated in
Gate 2.4 against gate2's own DC (lstsq + trace-gauge) solver (agreement to
5e-5 relative).

Gates certified (`gates_summary_step2.json`):

| Gate | Requirement (plan Sec. 5, Step 2) | Result |
|---|---|---|
| 2.1 (~U1) | O(kappa^-1) convergence to the ideal Riesz-projector cut; agreement with gate2's own cut | PASS (fit slope -0.99998; exact agreement with gate2's cut) |
| 2.2 (~U2) | implementation universality: two admissible cut generators with the same retained projector P_S converge to the identical ideal cut but differ at finite kappa | PASS (ideal-limit difference = 0 exactly; finite-kappa difference nonzero and itself O(kappa^-1)) |
| 2.3 (~U5) | non-invasiveness (C4): the cut generator annihilates the undriven reference steady state, and L0+kappa*D_S's steady state stays pinned to it for every kappa | PASS (both norms at or below machine precision) |
| 2.4 | the z-regularization doesn't distort the physics; the cut doesn't remove the one-photon absorption background when there's no control-induced coherence to cut (Oc=0) | PASS (5e-5 regularization error; background error at machine precision) |

Step 2 (plan: two independent GKSL-admissible cut implementations must
agree in the kappa -> infinity limit, matching the already-validated
Sector methodology) is therefore satisfied on the full many-level NV
model, not just the toy Lambda/3-level models it was first proven on.
Step 3 (merged-manifold moment analysis: symbolic M0=0 and the SMRT class
of the NV spin-Lambda sector) is next.

---

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

    src/step1_low_temperature_validation.py   Step 1
    src/step2_operational_cut_audit.py         Step 2
    results/gates_summary_step1.json
    results/gates_summary_step2.json
    results/figures/fig_step1_low_T_positive_control.png
    results/figures/fig_step2_operational_cut_audit.png

Future steps (3-9) will follow the same convention: one `stepN_*.py` per
step, one `gates_summary_stepN.json`, figures under `results/figures/`.
