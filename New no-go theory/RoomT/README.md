# Room-temperature SMRT no-go campaign (NV ZPL spin-Lambda)

Source strategy: `SMRT_NV_room_temperature_EIT_no_go_numerical_plan.md`
(user-provided, not committed here), which defines Steps 1-9, the four
EIT criteria (E1-E4), the SMRT class certificate, and the figure/table
plan for the room-temperature (300 K) no-go claim. This directory turns
that plan into an executable numerical campaign on top of the existing
infrastructure in `No-go theorem/` and `New no-go theory/`.

## Status

**Step 8 (reduced-kernel vs full-Liouvillian agreement): EXECUTED, all 4 gates PASS.**
See `results/gates_summary_step8.json` and
`results/figures/fig_step8_model_agreement.png`.
Code: `src/step8_reduced_vs_full_liouvillian.py`.

Compares Model R (`nv_model.py`'s reduced 6-level kernel, used for Steps
1/3/4/5/7's asymptotic and optimization work) against Model L (the full
N=9 Lindblad master equation, `gate2_candidate_full_vs_reduced.py`, used
in Steps 2/6) at the SAME reference geometry, T=10 K through 300 K, both
evaluated at a single unambiguous point (Model R at z=Ep, Model L at
d2=0) rather than via a peak search -- an exploratory attempt at a
window-max peak search hit the same ambiguity Steps 5/6 already found,
worse in the 30-70 K crossover region where multiple spectral features of
comparable size coexist (a +/-50 MHz window was needed to find a
"non-edge" extremum at T=30 K, and it did not correspond to the same
physical feature as the near-resonance dip). Evaluating at z=Ep/d2=0
sidesteps this entirely.

**Key finding**: Model R systematically OVERESTIMATES |C| relative to
Model L, by a factor shrinking monotonically with T: ~80x at 10 K, ~28x
at 30 K (also sign-disagreeing there), ~8x at 50 K, ~3x at 70 K, then
1.00-1.04x for T >= 100 K. This is a genuine, self-diagnosing
perturbative-breakdown effect, not a bug: Model R's adiabatic-elimination
formula (`dXi = -beta*K12*K21/den`) assumes the Raman correction is a
SMALL perturbation to the bare susceptibility. At T=10 K, Model R's own
`|C|~0.87` is manifestly not small -- signaling its own breakdown -- while
at T>=100 K, `|C|` has collapsed to ~1e-6 or below, confirming the
perturbative assumption holds. The correlation between log|C_modelR| and
log-deviation-from-Model-L across the whole T range is 0.90, confirming
Model R's own magnitude is a reliable self-diagnostic for when to trust
it quantitatively.

**Consequence for the rest of this campaign** (why this doesn't undercut
anything already concluded): Step 1's low-T contrast VALUES (0.89 at 4 K)
were only ever used QUALITATIVELY there (EIT works, cut kills it, correct
trends) -- unaffected. Steps 3/4's exact symbolic Gamma-scaling analysis
(nu_K=2, nu_R=4) is derived in and applies to the LARGE-Gamma regime,
exactly where this step shows Model R and Model L agree to a few percent
-- so the exponent result is now validated on the model that actually
matters (Model L), not merely self-consistent within Model R. Steps
5/6/7's 300 K conclusions were already close to model-independent (both
agree to ~4% there); two spot-checks with Model L directly (the Bx=0
reference point and Step 5's rejected Bx=0.2T candidate) both confirm
negative/negligible C at 300 K, matching Model R's verdict.

Gates certified (`gates_summary_step8.json`): `high_T_agreement_within_
10pct` (max relative error 4.0% for T>=100K), `low_T_discrepancy_
explained_by_perturbative_breakdown` (correlation 0.90), `model_L_
confirms_negative_or_negligible_at_300K_spotchecks`, `sign_agreement_
from_100K_up`.

Step 8 (plan: "Model Rが指数・機構を正しく再現し、Model Lが有限温度での上界を
確認すること") is satisfied, with the added value of an explicit,
quantitative, self-diagnosing account of where and why the two models
diverge. Step 9 (local response to detectable-signal conversion, tying
together `signal_chain.py`'s OD/SNR machinery already built during the
preflight work) is the last step in the plan.

---

**Step 7 (correction-mechanism upper bounds): EXECUTED, all 4 gates PASS.**
See `results/gates_summary_step7.json` and
`results/figures/fig_step7_correction_bounds.png`.
Code: `src/step7_correction_mechanisms.py`.

Four corrections evaluated at T=300 K on the reference geometry, each
reported as `eta_max_deviation` (how much the correction changes C from
baseline) and `C_EIT_max` (the largest POSITIVE/transparency value found,
per Step 6's Criterion-E1 lesson -- max(C,0), not abs(C)):

- **Hyperfine (14N secular, exact literature A_GS/A_ES)**: mI=-1,+1
  branches show |C| roughly 200x SMALLER than the mI=None baseline
  (-6.5e-12, -8.3e-12 vs -1.82e-9); the mI-ensemble average is
  -6.1e-10 -- hyperfine mixing SUPPRESSES the already-negative feature
  further, it does not restore transparency (`eta_max_deviation`~1.8e-9,
  `C_EIT_max`=0).
- **ISC / singlet shelving (exact literature branching)**: negligible
  effect (`eta_max_deviation`~6.9e-11).
- **Branch-linewidth asymmetry** (free parameter, up to a generous 50%
  radiative-rate difference between orbital branches X and Y -- not
  already in the model, needed a modified jump-operator construction,
  `build_full_asym`): C stays completely FLAT at -1.82e-9 across the
  entire 0-50% range (`eta_max_deviation`~1.7e-15) -- at 300 K, Gamma
  dominates so completely that this O(1) fractional rate asymmetry is
  irrelevant.
- **Polarization impurity** (free parameter, probe/control orbital
  polarization mixed up to 45 degrees off pure orthogonal X/Y): equally
  negligible (`eta_max_deviation`~4.5e-15), confirming Step 5's structural
  argument (M0=0 is protected by ground-state orthogonality regardless of
  polarization, so this can only perturb the subleading response, and
  here it barely does even that).
- **Combined worst case** (all four simultaneously, each at its most
  adversarial individual setting): C_EIT_max = 0 -- no combination
  restores positive (transparency-signed) contrast.

Transverse/axial field, strain magnitude/azimuth, and "off-axis
orientation" are NOT re-done here -- Step 5's global optimization already
searched Bx, Bz, d_strain, phi over the full realistic domain (a
superset of what a dedicated small-perturbation sweep would show) and
found no correctly-signed transparency anywhere.

A quick robustness check during development: since Step 6 found d2=0 is
a special (non-representative) point for one particular small-Oc
expansion, the most adversarial correction setting (branch asymmetry
eta=0.5) was re-scanned across a +/-20 MHz window to rule out a hidden
positive peak nearby that a d2=0-only evaluation might miss. The window
scan confirms the d2=0 evaluation is representative: the only positive
values found are at the window edge and utterly negligible (~1e-14, five
orders of magnitude below the detection threshold).

Gates certified (`gates_summary_step7.json`): `each_correction_below_
threshold`, `combined_worst_case_below_threshold`, `combined_worst_
case_tiny_absolute`, `combined_not_worse_than_sum_of_individual`.

Step 7 (plan: show small symmetry-breaking corrections, individually and
combined, do not restore room-temperature EIT) is satisfied. Step 8
(reduced-kernel vs full-Liouvillian agreement across the full parameter
range) and Step 9 (local-response to detectable-signal conversion) remain.

---

**Step 6 (separating apparent dips from genuine EIT, Criteria E1-E4): EXECUTED, all 8 gates PASS.**
See `results/gates_summary_step6.json` and
`results/figures/fig_step6_dip_discrimination.png`.
Code: `src/step6_dip_discrimination.py`.

Reuses `No-go theorem/src/gate2_candidate_full_vs_reduced.py`'s
`full_spectrum()` (the FULL N=9 Lindblad master equation with a genuine
independent two-photon detuning d2 -- unlike the reduced kernel's single
shared z, per Step 4/5's finding) to compute full vs EIT-sector-cut
absorption spectra at both T=10 K (Step 1's positive control) and T=300 K,
on the identical reference geometry used throughout. Four criteria (plan
Sec. 6) are checked at BOTH temperatures, not asserted only at 300 K:

- **T=10 K (positive control): genuine EIT, all four criteria satisfied.**
  E1 transparency (C(d2)>0, peak +1.07e-2, full spectrum visibly below cut
  spectrum); E2 sector causality (R_EIT=0 exactly at Oc=0); E4 continuous
  quadratic emergence from Oc=0 (fit slope ~2.0, verified over five
  decades of Oc AT THE ACTUAL PEAK location, not an arbitrary d2=0 -- see
  bug note below); ground coherence proxy 0.0127 (order-1%, substantial).
- **T=300 K: NOT EIT -- fails E1 outright.** C(d2) is negative
  EVERYWHERE in the scanned +/-20 MHz window (peak -1.82e-9 exactly at
  d2=0), meaning the control field INCREASES absorption there -- the
  opposite of transparency. This confirms, now across the full two-photon
  detuning window and the FULL (not reduced) Liouvillian, exactly what
  Step 5 found at the single resonance point. The full and cut spectra
  are visually indistinguishable at any practical plotting resolution
  (differ only at the 1e-9 level on a ~1.6e-4 background) -- there is no
  appreciable dip of ANY kind (real or fake) to worry about at 300 K, not
  merely a fake one that fails discrimination. Ground coherence proxy is
  ~1e-5, about 970x smaller than at 10 K.

**Bug caught during development**: the first version of the E4
(continuous-emergence) check evaluated C at the arbitrary fixed point
d2=0 for small Oc. At T=10 K this gave a spurious SIGN FLIP and an
apparent cubic (not quadratic) power law -- d2=0 turned out to be a
non-representative, near-degenerate point for the small-Oc expansion of
this particular model. Evaluating instead at the ACTUAL feature peak
(d2~0.2-0.25 MHz, found from the reference-Oc spectrum) gives a clean,
single-signed Oc^2 scaling across five decades (1e-4 to 1 GHz) before
saturating at large Oc -- the expected, unambiguous EIT power-broadening
signature. This does not affect T=300 K's conclusion (which fails E1
regardless, independent of the E4 check), but was necessary to get a
trustworthy positive-control confirmation at 10 K.

Gates certified (`gates_summary_step6.json`): `low_T_passes_E1_transparency`,
`low_T_passes_E2_sector_causality`, `low_T_E4_continuous_emergence`,
`room_T_fails_E1_transparency` (yes, this is the PASS direction -- the
gate is that room-T genuinely fails E1, confirming no-go), `room_T_passes_
E2_sector_causality`, `room_T_dip_not_edge_artifact`,
`ground_coherence_far_smaller_at_room_T`.

Step 6 (plan: "dip exists AND R_EIT~=0 => not EIT") is satisfied with a
stronger result than the plan's minimum bar: at 300 K there is no
appreciable dip of any kind, real or fake, and what negligible residual
exists fails Criterion E1 by sign, not merely by magnitude. Step 7
(correction mechanisms: hyperfine, strain, ensemble, linewidth asymmetry)
is next.

---

**Step 5 (300 K global adversarial optimization): EXECUTED, all 6 gates PASS.**
See `results/gates_summary_step5.json` and
`results/figures/fig_step5_global_optimization.png`.
Code: `src/step5_global_adversarial_optimization.py`.

Latin-hypercube pre-sampling (800 points) + local refinement + an
independent differential-evolution search over Oc, gg, Bx, Bz, d_strain,
phi, at T=300 K, using the `conservative_lower_bound` Gamma_XY(T,d) model
(Step 4's pick). Two separate domains are searched: MAIN (plan Sec. 4's
primary declared realistic range -- what Gate 5 is checked against) and
STRESS (plan Sec. 4's explicit non-physical-extreme upper bounds,
reported separately per plan Sec 1.2/4.2/4.3/16 -- a positive finding
there does not count against the main claim).

**Two real bugs were caught during this run, both material to the
conclusion, not just cosmetic:**

1. **Wrong sign convention in the objective (the important one).** The
   first version scored `abs(C)`. The adversarial search promptly found
   Bx=0.2 T -- well inside the MAIN range -- giving `|C|~1.85e-6`,
   apparently exceeding the illustrative detection threshold
   (`epsilon_th~1.53e-7`) by ~12x. Investigating (mechanism check: does
   the signal vanish as Oc->0? yes, correctly) initially suggested this
   was a genuine, large EIT-type signal -- and it is the SAME B_perp-
   induced branch-mixing mechanism already studied at T=70 K elsewhere in
   this repository (`No-go theorem/`'s SIMULATION_PLAN.md / PRL
   candidate, BX0~0.23 T). Checking Criterion E1 (transparency: Im
   chi_full < Im chi_cut) directly -- computing Acut, dA, Afull
   explicitly -- showed C<0 at this point: Afull > Acut, i.e. the control
   field INCREASES absorption there. That is the opposite of EIT. Scoring
   `abs(C)` let the adversarial search "win" by finding absorption
   increase instead of transparency, which is not a valid escape route.
   Fixed to `max(C,0)`, matching the plan's own objective
   `C_EIT(theta) = max[...]_+` (Sec. 5) exactly. After the fix, the MAIN
   and STRESS range global optima are both consistent with 0 (no
   detectable transparency anywhere sampled in either domain).
2. **A companion structural finding surfaced by chasing bug #1**: at the
   Step 1/3/4 reference configuration itself (Bx=0), C(T) does not merely
   shrink in magnitude with T -- it changes SIGN, from positive
   (transparency, matching Step 1's positive control) at T <~ 50 K to
   negative (absorption-increase) from T >~ 77 K onward (crossover
   bracketed at [50, 77] K; see `sign_crossover_scan` and the bottom-right
   panel of the figure). By 300 K, even the "conventional" reference
   configuration is not weakly-EIT -- it is not EIT-signed at all. Step
   1/3/4's `abs(C)`-based monotonicity checks did not catch this (a
   magnitude check can't detect a sign flip); it does not change any of
   their PASS verdicts (the magnitude claims all still hold) but sharpens
   the physical picture.

Gates certified (`gates_summary_step5.json`): `M0_robust_under_Bx` (ground-
state orthogonality protects M0=0 for ANY Bx, Bz -- verified directly, not
just at Bx=0 -- so opening the transverse field cannot reopen the leading-
order channel, only perturb the already-suppressed subleading one),
`main_range_max_below_detection_threshold`, `main_range_max_tiny_absolute`,
`stress_range_max_still_far_below_unity`, `lhs_and_DE_consistent`,
`reference_config_sign_crossover_found`.

Step 5 (plan: search the full realistic domain simultaneously, not just
convenient points) is satisfied: no combination found anywhere in the
declared MAIN domain gives detectable, correctly-signed (transparency)
contrast at 300 K, and the same holds even in the STRESS-test domain.
Step 6 (separating apparent dips from genuine EIT via Criteria E1-E4) is
largely already anticipated by this step's Criterion-E1 check; the
remaining Step 6 work is characterizing the OTHER failure modes (ATS,
saturation, pumping) for whatever residual features exist.

---

**Step 4 (temperature scaling, crossover resolution): EXECUTED, all 4 gates PASS.**
See `results/gates_summary_step4.json` and
`results/figures/fig_step4_temperature_scaling.png`.
Code: `src/step4_temperature_scaling.py`.

Reuses Step 3's exact symbolic certificate to get closed-form leading
coefficients of BOTH R_EIT asymptotes (pre- and post-denominator-crossover,
see Step 3 section below), then checks numerically, using the actual
`nv_model.response()` (not just the asymptotic formulas) evaluated on
resonance (z=Ep, matching Step 1's convention) across the campaign
temperature grid and four `phonon_rates.py` Gamma_XY(T) models:

- **x1(300K) = Gamma(300K)/d_strain ~ 7559 >> 1**: the excited manifold IS
  deep in the Hamiltonian-merged regime at 300 K (the plan's own x(T)
  criterion), confirmed for the most-conservative (`conservative_lower_bound`)
  dissipation model.
- **300 K sits in the PRE-crossover Gamma^-3 regime, not yet Gamma^-4**:
  the Step 3 flag is now resolved quantitatively -- R_EIT(300K) matches
  the exact Gamma^-3 asymptote to 16.5% relative error (vs 606% for the
  Gamma^-4 one), consistent with Gamma(300K)~1.27e4 GHz sitting below
  Gamma_cross~7.7e4 GHz at the Step 1/3 reference point (Oc=1 GHz,
  gg=6.3e-5 GHz).
- **|C_EIT|(T) collapses monotonically and enormously**: 0.89 (4 K) ->
  0.026 (30 K) -> 1.3e-7 (150 K) -> ~2e-9 (300 K), for every one of the
  four Gamma_XY(T) models tried (full Happacher, literature-uncertainty
  conservative lower bound, saturation, naive T^5 extrapolation) -- they
  agree closely up to ~100 K and differ only by factors of a few by 300 K,
  never changing the qualitative conclusion.
- **Control-power scan at 300 K (correcting a Step-3 speculation)**:
  |C_EIT|(300K) increases MONOTONICALLY with Omega_c across the plan's
  entire declared range (1 MHz .. 10 GHz) and plateaus at only ~2.4e-9 at
  the top of that range -- there is no interior optimum at small Omega_c.
  Step 3's remark that smaller (more realistic) control powers would
  "make the asymptotic regime easier to reach, strengthening the no-go"
  was correct about WHICH regime (Gamma_cross falls with Omega_c^2, so
  small Omega_c does sit deeper in the Gamma^-4 tail) but misleading about
  practical consequence: raw signal magnitude is smallER at smaller
  Omega_c regardless of regime, so more control power never helps recover
  the signal in the declared realistic range, and even the best case found
  is 8+ orders of magnitude below anything detectable.

Gates certified (`gates_summary_step4.json`): `x1_300K_much_greater_than_1`,
`C_EIT_monotone_decreasing_main_model`, `crossover_diagnosis_self_consistent`
(the regime R_EIT is actually in matches its own asymptote to <30% error),
`conservative_model_still_tiny_at_300K` (|C_EIT| < 1e-3, in fact ~2e-9).

Step 4 (plan: verify 300 K is deep in the applicable asymptotic regime
using conservative phonon models) is satisfied, with the caveat resolved
honestly rather than assumed. Step 5 (300 K global adversarial
optimization over the full realistic parameter domain) is next.

---

**Step 3 (merged-manifold moment analysis): EXECUTED, all 8 gates PASS.**
See `results/gates_summary_step3.json` and
`results/figures/fig_step3_moment_scaling.png`.
Code: `src/step3_merged_manifold_moments.py`.

At Bx=By=0 (Step 1/2's geometry), `nv_model.Hgs` is exactly diagonal, so
the probe/control legs dp, dc are EXACT closed-form basis vectors (no
eigendecomposition needed) -- this makes the whole moment analysis exact
symbolic algebra rather than a numerical approximation. Using sympy's
adjugate/determinant method (`New no-go theory/src/core.py`'s Theorem 8.1
certificate style) on the 6x6 excited-manifold response operator:

- **M0 = dp^T dc = 0 exactly** -- trivial orthogonality of the probe (orbital
  X) and control (orbital Y) branches, a group-theoretic selection rule
  holding for ANY Hamiltonian parameter values, not a numerical coincidence.
- **nu_K = 2**: the first nonzero kernel moment is at order Gamma^-2, from
  the exact closed form **H[1,3] = -i\*Lperp/sqrt(2)** (Dpar, Lpar, Dperp all
  drop out of this specific matrix element; Lperp alone generates the
  branch-mixing at this order, proven nonzero generically by keeping Lperp
  free in the certificate).
- **nu_R = 4**: R_EIT (the actual EIT observable, not just the kernel)
  inherits nu_K12 + nu_K21 = 2+2 = 4 asymptotically, NOT nu_K = 2 -- exactly
  the nu_K != nu_R distinction the plan's Step 3 warns about, since R_EIT
  enters as the PRODUCT K12*K21 divided by a denominator whose
  Gamma-independent part (ground decoherence geff) eventually dominates its
  Gamma-dependent correction (S2 ~ Gamma^-1, nu_S2=1).
- SMRT class: **Class II** (finite suppression order, not Class I exact
  zero, not Class III protected) -- the correct outcome for the no-go
  claim (Class III would reject it).
- Exact leading coefficients (both reused directly by Step 4):
  `R_EIT ~ -2*pi^2*Lperp^2 / Gamma^3` (Gamma << Gamma_cross) and
  `R_EIT ~ -2*pi^2*Lperp^2*beta/geff / Gamma^4` (Gamma >> Gamma_cross).

**Bug caught and fixed during development**: the first version of the
symbolic Hamiltonian transcribed the Dperp term as
`kron(sz_o,SySy-SxSx) - kron(sx_o,SxSz+SzSx)` where the real
`nv_model.Hes` has `kron(sz_o,SySy-SxSx) - kron(sx_o,SxSy+SySx)` (Sx@Sy,
not Sx@Sz, in the second piece). This gave a spurious Dperp contribution
to H[1,3] and a leading coefficient of `(Dperp^2+Lperp^2)` instead of the
correct `Lperp^2` alone -- caught by cross-checking the exact symbolic K12
against `nv_model`'s numeric K12 at Gamma(300K) (Step 4), which disagreed
by a factor of ~5 before the fix and agree to 5e-6 (float64 sqrt(2)
rounding only) after it. All eight Step 3 gates were re-verified and still
PASS with the corrected Hamiltonian (expected: the bug was in a
branch-orthogonal, Gamma-independent diagonal-adjacent term that could
only change the coefficient, not the M0=0 result or the nu_K=2/nu_R=4
degrees). The corrected coefficient is smaller (Lperp^2=0.02 vs the
erroneous Dperp^2+Lperp^2=0.621), so the fix makes the no-go conclusion
slightly STRONGER, not weaker.

A second, related fix: the leading coefficients were originally verified
at the certificate's default z=0 rather than the actual probe resonance
z=Ep, giving numerically correct but not directly Step-1-comparable
contrast values (C_EIT(4K) came out as -0.03 instead of Step 1's 0.89).
`build_symbolic_certificate` now takes a `z_val` parameter; Step 4 rebuilds
it at z=Ep. The leading coefficients turned out to be exactly
z-independent regardless (a consequence of M0=0: the z-term in the moment
expansion is multiplied by M0 and drops out), so only the raw C_EIT(T)
trend needed the z=Ep fix, not the asymptotic degree/coefficient analysis.

Gates certified (`gates_summary_step3.json`):

| Gate | Requirement (plan Sec. 5/7, Step 3) | Result |
|---|---|---|
| `M0_exact_zero` | M0=0 via exact identity, not a small float | PASS (dp^T dc = 0 by orbital orthogonality) |
| `nu_K12_symbolic_equals_2` | exact polynomial-degree certificate for the kernel | PASS (deg_Q=6, deg_N12=4) |
| `nu_R_symbolic_equals_4` | exact degree certificate for the actual R_EIT observable | PASS (deg_num=8, deg_den=12) |
| `nu_K_neq_nu_R` | kernel and observable exponents reported separately, not conflated | PASS (2 != 4) |
| `numeric_fit_matches_symbolic_K12` / `_R` | independent numeric log-log fit (deep asymptotic tail) agrees with the exact degree | PASS |
| `precision_stable` | exponent unchanged between float64 and 50-digit mpmath | PASS |
| `not_class_III` | reject the no-go if a protected (Class III) response were found | PASS (Class II confirmed) |

Step 3 (plan: "M0=0をsymbolicに証明", "最初の非零momentを特定",
"nu_Kとnu_Rを分離") is therefore satisfied with an exact, non-numerical
certificate. See the Step 4 section above for the crossover resolution.

---

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
    src/step3_merged_manifold_moments.py       Step 3
    src/step4_temperature_scaling.py           Step 4
    src/step5_global_adversarial_optimization.py  Step 5
    src/step6_dip_discrimination.py             Step 6
    src/step7_correction_mechanisms.py          Step 7
    src/step8_reduced_vs_full_liouvillian.py    Step 8
    results/gates_summary_step1.json
    results/gates_summary_step2.json
    results/gates_summary_step3.json
    results/gates_summary_step4.json
    results/gates_summary_step5.json
    results/gates_summary_step6.json
    results/gates_summary_step7.json
    results/gates_summary_step8.json
    results/figures/fig_step1_low_T_positive_control.png
    results/figures/fig_step2_operational_cut_audit.png
    results/figures/fig_step3_moment_scaling.png
    results/figures/fig_step4_temperature_scaling.png
    results/figures/fig_step5_global_optimization.png
    results/figures/fig_step6_dip_discrimination.png
    results/figures/fig_step7_correction_bounds.png
    results/figures/fig_step8_model_agreement.png

Step 9 (local response to detectable-signal conversion) will follow the
same convention: one `stepN_*.py`, one `gates_summary_stepN.json`, a
figure under `results/figures/`.
