# Phase D plan: phenomena invisible to the old no-go framework

Status: **PLAN ONLY — not yet executed.** Awaiting approval before any code is run.

This document identifies theoretical results / phenomena that exist in the *new*
master-response theorem package (`master_response_no_go_theorem_package.pdf`,
`Theorem and proofs/three_theorems_proofs.tex`) but are structurally invisible to
the *old* no-go framework (`eit_nogo_lecture.pdf`: per-kernel criterion
`K12*K21 = 0`, invertible-D moment hierarchy, Hermitian orthogonal-projector
protected channels, class attached to the configuration alone). It then defines
a numerical campaign (Phase D, runs D1–D4) to demonstrate each one.

Already covered by Phases A/B/M (do **not** repeat): hidden class transition
(nu[chi_full]=0 while nu[R_S]: 0->1), scaling collapse, QFI nu->2nu translation,
EIT/ATS/background mechanism map.

---

## 0. What exactly is "new" in the new theory

| Ingredient | Old framework (lecture) | New package (master response) |
|---|---|---|
| Object classified | kernels `K12`, `K21` separately; scalar criterion `delta_chi_S = 0 <=> K12 K21 = 0` (valid only in the 2-lower-state scalar block) | the difference `R_S = chi_full - chi_cut` itself, via master moments `m_k = p†[(D^-1 B)^k D^-1]c` on the *augmented* (doubled) space — full-vs-cut cancellation is inside the classified object (Thm 5.1, Eq. 12) |
| Singular dissipation | only `D = D† >= 0`, orthogonal projector onto `ker D` (lecture Eq. 23) | any semisimple zero of a generally **non-Hermitian** `D`, **Riesz (oblique) projection** (Assumption 6.1, Thm 6.2) |
| What the class is attached to | the configuration tuple `C` | the triple (response, `D`, scaling path `Gamma -> Gamma D + B`); reassigning terms between `D` and `B` legitimately changes the class (Remark 3.2) |
| Validity range of the classification | Neumann series, `Gamma > ||X||` | finite polynomial certificate `nu = deg_Gamma Q - deg_Gamma N` from adjugate/determinant — exact, valid **beyond** the convergence radius (Thm 8.1 Step 5, Remark 8.2) |

Each row is one blind spot of the old theory, and each yields one demonstrable
phenomenon (D1–D4 below).

---

## D1. Oblique protection and dissipation-enhanced sector response

**Phenomenon.** When the fast dissipator `D` is non-normal with a semisimple
zero eigenvalue whose kernel is *not orthogonal* to its complement, the correct
leading protected coefficient is `R_{S,0} = p†P[B_P,full^-1 - B_P,cut^-1]Pc`
with the **Riesz** projector `P`. Two consequences:

1. *Misclassification by the old formula.* The old orthogonal projector
   `P_orth = Proj(ker D)` predicts a wrong `r0` — quantitatively wrong by a
   factor growing with the non-normality angle, and, at tuned parameters,
   *qualitatively* wrong: `r0_orth = 0` (old theory says Class II, suppressed)
   while the true Riesz `r0 != 0` (Class III, protected) — or the reverse
   ("false protection").
2. *Non-monotonic `|R_S(Gamma)|`.* The oblique projector has condition number
   `kappa(P) > 1`; the subleading terms carry `kappa(P)`-enhanced coefficients,
   so `|R_S|` can *rise* with increasing dissipation over an intermediate
   window before settling on the protected plateau — "dissipation-enhanced
   sector response". The old asymptotic taxonomy (a single suppression order
   per kernel) has no slot for a non-monotonic dissipation dependence; the new
   package predicts both the bump scale (via `kappa(P)`, cf. lecture Remark 6.3
   which flags but does not resolve pseudospectral amplification) and the exact
   plateau.

**Why the old theory cannot see it.** Lecture Eq. (23) is derived only for
`D = D† >= 0`; a Lindblad-space dissipator is generically non-Hermitian, so the
orthogonal-projector formula is simply the wrong formula there, and the old
class table has no observable that distinguishes "wrong r0" from "right r0".

**Minimal model (3-dim optical block + sector, N=4 augmented 8).**

```
D(a) = [[0, a, 0],
        [0, 1, 0],
        [0, 0, 2]]          # semisimple zero; ker D = span{e1 - ... oblique for a != 0}
```

`a` is the non-normality knob (`a = 0` recovers the Hermitian old case —
built-in control experiment). `B(z)` a generic analytic coupling block; frozen-
source cut on the sector coupling exactly as in Phase B.

**Compute / plot.**
- exact Riesz `P(a)` (closed form for the 3x3 above; cross-check
  `core.riesz_projection_diag`), `kappa(P)(a)`;
- `r0_Riesz(a)` vs `r0_orth(a)` vs directly measured plateau
  `lim |R_S(Gamma)|`;
- misclassification region in the `(a, lambda)` plane where
  `sign/zero structure` of `r0_orth` and `r0_Riesz` disagree;
- `|R_S(Gamma)|` curves showing the intermediate bump; bump height vs
  `kappa(P)`.

**Gate D1 (pass criteria).**
- at `a = 0`: Riesz = orthogonal (control, machine precision);
- for `a != 0`: measured plateau matches `r0_Riesz` to <1e-6 relative, and
  deviates from `r0_orth` by the predicted factor;
- exhibit at least one parameter point with `r0_orth = 0` but plateau
  `= r0_Riesz != 0` (old theory predicts suppression, response survives), and
  one point with the reverse;
- exhibit non-monotonic `|R_S(Gamma)|` with bump height growing monotonically
  with `kappa(P)`.

**Figures.** `figD1a_oblique_r0.png`, `figD1b_dissipation_enhanced.png`,
`figD1c_misclassification_map.png`.

---

## D2. Cancellation-promoted suppression and the symmetry-free exact zero

**Phenomenon.** The master moments classify the *difference*; full and cut
contributions can cancel order by order. Therefore:

1. *Order promotion:* `nu[R_S]` can strictly exceed the old sector-graph bound
   `d+1` (minimal coupling-insertion path length), because `m_k(full)` and
   `m_k(cut)` are individually nonzero but identical for `k < m` with `m > d`.
   The old theory, classifying kernels separately and taking the minimal path,
   predicts `nu = d+1` and is wrong.
2. *Extreme case — symmetry-free exact Class I:* in a multi-excited-branch
   (Mishina-type) model, couplings can be tuned so that *every* master moment
   vanishes identically on a parameter domain — an exact all-order no-go —
   while `K12 K21 != 0` pointwise and no reducing symmetry exists. The old
   scalar criterion (18) declares "go"; the old symmetry audit finds nothing;
   only the Krylov certificate on the *augmented* space (Thm 4.2 + Cor. 4.3)
   detects it.

**Why the old theory cannot see it.** Its exact-zero tests are (a) the scalar
`K12 K21 = 0` criterion — false here — and (b) reducing symmetry — absent here.
Its asymptotic order is the per-kernel first moment — too small here.

**Minimal model.** Two lower states, `Ne = 2..3` excited branches, diagonal
`A = diag(gamma_j - i*Delta_j)` so all kernel sums are closed-form
(`K12 = sum_j d1j* d2j / a_j`); sector cut = frozen-source cut of the ground-
coherence coupling. Tune complex dipoles `{d_ij}` such that on the augmented
space the first `m` master moments cancel (linear conditions on the `d`'s —
solvable exactly). For the Class-I endpoint, impose all `N_resp` moment
conditions and verify with **exact rational/symbolic arithmetic** (sympy),
never floats (Assumption 1.3 / Remark 1.4 discipline; floating-point smallness
is inadmissible as a Class-I proof).

**Compute / plot.**
- `m_k(full)`, `m_k(cut)`, `m_k(master)` side by side vs `k` — the cancellation
  ladder;
- `nu_eff(Gamma)` fit vs sector-graph bound `d+1` vs master-moment prediction;
- for the Class-I point: symbolic certificate output (all `N_resp` moments
  identically zero as rational functions) + numeric `|R_S|` at machine-noise
  floor, alongside `|K12 K21| = O(1)`.

**Gate D2.**
- at least one point with `nu[R_S] >= d+2` while per-kernel analysis gives
  `d+1`, confirmed by moment method AND direct fit (Gate-1-style consistency);
- one symbolically certified Class-I point with `K12 K21 != 0` on the domain
  and no reducing projector (audit table included);
- control: detuning a single dipole phase by `epsilon` demotes the point back
  to the generic order, with `|R_S|` reopening `∝ epsilon` (shows the zero is
  cancellation, not hidden symmetry).

**Figures.** `figD2a_cancellation_ladder.png`, `figD2b_promotion_fit.png`,
`figD2c_classI_no_symmetry.png`.

---

## D3. Scaling-path dependence: the class belongs to the knob, not the material

**Phenomenon.** Remark 3.2: the classification is attached to the triple
(response, `D`, path `Gamma -> Gamma D + B`). The *same* physical generator at
the *same* physical parameter point carries **different classes along different
physically meaningful scaling paths**. Concretely: scaling only the phonon
(temperature-driven) channel — `D_phonon` singular, radiative decay left in
`B` — can give Class III (protected plateau vs temperature), while scaling all
dissipation uniformly — `D_all` invertible — gives Class II with some
`nu = m >= 1`. Same model, same `chi_full(z)` at every matched point; different
experiment (temperature scan vs pressure/broadening scan) ⇒ different no-go
verdict and different collapse law.

**Why the old theory cannot see it.** The old taxonomy assigns one class to a
configuration `C`; "Gamma" is a single unlabeled fast scale read off the
optical-coherence equation. There is no slot for two inequivalent
decompositions of the same generator, so the old framework cannot even state —
let alone predict — that the temperature scan and the broadening scan of the
same sample disagree in class.

**Minimal model.** 4–5 dim optical block with two dissipation channels:
`A(z, Gamma) = Gamma * D_path + B_path(z)` with
- Path T (temperature-like): `D = D_phonon` (rank-deficient, semisimple zero on
  a phonon-dark subspace), `B` contains radiative + detunings;
- Path U (uniform): `D = D_phonon + D_rad` (full rank), `B` = detunings only.
Choose rates so both paths pass through the identical physical point at
`Gamma = Gamma_0` (calibration constraint: `Gamma_0 D_T + B_T = Gamma_0 D_U + B_U`).

**Compute / plot.**
- `nu_eff(Gamma)` along both paths from the shared physical point; classifier
  (moments / protected coefficient / fit) verdicts per path;
- verify `chi_full`, `chi_cut`, `R_S` numerically identical at the calibration
  point (sanity: it is one physical system);
- 2-channel mixing sweep `D(theta) = cos(theta) D_phonon + sin(theta) D_rad`
  scaled path: class as a function of path direction `theta` — a "class fan"
  around one physical point, with the critical `theta_c` where `D(theta)` loses
  rank.

**Gate D3.**
- the two paths return different classes (target: III vs II) from the same
  calibrated point, each internally consistent across the three classifier
  methods;
- the class flip coincides with the analytically computed rank change /
  protected-block criterion of `D(theta)`;
- physical translation stated in the summary: predicted observable difference
  between a T-scan and a broadening-scan collapse exponent in an NV-like
  effective model.

**Figures.** `figD3a_two_paths.png`, `figD3b_class_fan.png`.

---

## D4. Preasymptotic masquerade and the finite polynomial certificate

**Phenomenon.** The Neumann expansion converges only for `Gamma > X_K`; with
two widely separated internal scales the effective exponent `nu_eff(Gamma)`
can sit on a **false integer plateau** (e.g. `nu_eff ≈ 1` over 2+ decades)
before crossing over to the true asymptotic order (e.g. `nu = 3`). Any
finite-window fit — i.e. any experiment, and the old theory's only asymptotic
instrument — misassigns the class. The new package gives (i) the exact class
at zero asymptotic cost via `nu = deg_Gamma det - deg_Gamma adj` (polynomials
of bounded degree, exact arithmetic), and (ii) a *predicted crossover scale*
from the explicit remainder bound (Thm 5.1, Eq. 14): the false plateau ends at
`Gamma ~ X_K` computable a priori.

**Why the old theory cannot see it.** It has the moment hierarchy but no
degree certificate and no uniform remainder bound, so it cannot distinguish a
true plateau from a preasymptotic one, and cannot say *where* the crossover
must occur.

**Minimal model.** Block model with a weakly coupled far sector:
`B = B_near + eta * B_far` with `||B_far|| >> ||B_near||`, `eta` small, tuned so
the first master moment is dominated by a near-path term that cancels only at
higher order. `X_K = ||D^-1 B||` made large (~1e3) so the false plateau spans
the "experimentally accessible" window.

**Compute / plot.**
- `nu_eff(Gamma)` over 6+ decades with the false plateau and the true tail;
- exact certificate: construct `N_S(Gamma, z)`, `Q_S(Gamma, z)` symbolically
  (sympy, augmented space, dimension <= 10) and report
  `deg Q - deg N`;
- overlay the Eq.-(14) remainder bound and the predicted crossover
  `Gamma_* = X_K`; show fitted-window verdicts ("what a finite experiment would
  conclude") vs certificate verdict, as a function of window position.

**Gate D4.**
- a >= 2-decade false plateau at an integer `!=` true `nu`, followed by
  convergence to true `nu` (fit within 1%);
- certificate returns the true `nu` exactly and cheaply;
- measured crossover location within a factor ~3 of the Eq.-(14) prediction;
- a "misassignment table": for each fit window, the class a window-limited fit
  would report — demonstrating that only the certificate is window-independent.

**Figures.** `figD4a_masquerade.png`, `figD4b_certificate_vs_fit.png`.

---

## Implementation plan

```
New no-go theory/src/
  model_oblique.py       # D1 model + closed-form Riesz projector
  model_cancel.py        # D2 multi-branch cancellation model (+ sympy exact part)
  model_paths.py         # D3 two-path / theta-fan construction
  model_masquerade.py    # D4 two-scale model + symbolic N,Q degree certificate
  run_phase_d.py         # runs D1–D4, writes gates_summary_phaseD.json,
                         # appends to results/summary.md, saves figures
```

- Reuse `core.py` (`moments_invertible_D`, `nu_from_moments`,
  `riesz_projection_diag`, `fit_nu_loglog`, `max_abs_on_grid`); add
  `riesz_projection_exact` (contour-free, via eigenprojector algebra, with an
  orthogonal-projector counterpart `proj_orth_kernel` for the D1 comparison)
  and `certificate_deg_nu` (sympy `Matrix.adjugate`/`det` degrees in `Gamma`).
- All symbolic Class-I / certificate checks in exact arithmetic (sympy
  Rational); floats used only for fits and figures.
- Style, gate reporting, and JSON summaries follow Phase A/B/M conventions;
  everything headless (`matplotlib.use('Agg')`), fixed seeds, target runtime
  < 10 min total on laptop CPU.

**Execution order & effort:** D1 → D3 → D2 → D4 (D1 and D3 are the two most
publication-relevant claims: "old formula misclassifies non-Hermitian
protection" and "the class belongs to the scaling path"). Estimated ~1 session
of implementation + runs.

**Overall success statement (Gate D-final).** Each of D1–D4 exhibits a
quantity where the old framework's prediction is either undefined or wrong,
and the new package's prediction is verified by at least two independent
methods (analytic/symbolic + direct numerics).
