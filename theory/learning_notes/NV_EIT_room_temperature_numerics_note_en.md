# Why EIT in NV Centers Disappears at Room Temperature — An Answer from Numerics

> Unless stated otherwise, all numbers come from an unreduced nine-level full Liouvillian (GKSL) calculation.
> Fixed parameters: $B\_z = 5\,\mathrm{mT}$, $\Omega\_c = 0.1\,\mathrm{GHz}$, transverse strain $\delta = 1.683\,\mathrm{GHz}$, $D\_\mathrm{gs} = 2.877\,\mathrm{GHz}$

## 0. Conclusions First

The dominant reason that EIT in the NV$^-$ center cannot be observed at room temperature is that phonon-induced orbital hopping destroys the optical coherence and thereby destroys the Raman interference pathway itself. The numerics give five conclusions.

- The region where EIT occurs is a closed island in the plane of temperature $T$ and transverse magnetic field $B\_\perp$. It closes not only on the high-temperature side but also on the low-temperature side.
- At $B\_\perp = 0.232\,\mathrm{T}$, the temperature at which the contrast falls below 1% is $T\_{1\mathrm{pct}}=74\,[70,78]\,\mathrm{K}$.
- The response reverses sign at $T\_\mathrm{sign}=103\,[97,108]\,\mathrm{K}$, where the control field turns to increasing the probe absorption.
- The transverse field is a switch that opens the EIT pathway, not a dial that moves the temperature ceiling.
- The residual contrast at room temperature is $1.1\times10^{-9}$, nine orders of magnitude below the assumed detection limit. Longer integration does not reach it.

## 1. The Four Questions This Note Answers

1. Why is EIT hard to observe in NV centers at room temperature?
2. What is the theoretical ceiling temperature for observing EIT?
3. Why can alkali atoms show EIT at room temperature while NV centers cannot?
4. What is the transverse magnetic field actually doing?

Coherent population trapping in single NV spins [S02] and EIT in NV ensembles [S01] have both been demonstrated at cryogenic temperatures, so the question is not whether the pathway exists at all, but what closes it as the temperature rises.

Below, the common microscopic mechanism is presented first, and the four questions are answered in order afterwards.

---

## 2. Microscopic Mechanism: Orbital Hopping Destroys Coherence

The excited state $^3E$ of NV$^-$ is a six-dimensional orbital doublet carrying spin-orbit, spin-spin, transverse strain, and Zeeman terms [S07][S03][S08]. As the temperature rises, the phonon-induced orbital hopping $E\_x\leftrightarrow E\_y$ becomes stronger. This process conserves the population within the excited state, while acting as pure loss on the optical coherence.

Writing the optical coherence as $\sigma=P\_e\rho P\_g$, the jump operators inside the excited state $L\_\mu=P\_eL\_\mu P\_e$ satisfy

```math
P_{e}\,\mathcal{D}[L_{\mu}](\rho)\,P_{g}
=-\tfrac12 L_{\mu}^{\dagger}L_{\mu}\,\sigma
```

The recycling term $L\rho L^\dagger$ does not return to the optical-coherence block. Orbital hopping therefore damps the Raman pathway without removing population.

Denoting the symmetric two-way hopping rate by $k$, the population-imbalance relaxation rate commonly quoted in the literature [S04][S05] is $\Gamma\_{xy}=2k$, while the optical-coherence decay rate is $k/2=\Gamma\_{xy}/4$. When feeding values into a calculation, it is important not to confuse this factor of $1/4$.

### 2.1 Basic Equations of the No-Go Theory

The no-go criterion does not look at the magnitude of the total response alone. The decision quantity is the difference between the response that retains the long-lived coherence sector $S$ of interest and a counterfactual in which only that pathway is cut,

```math
\delta\chi_S
:=\chi_{\mathrm{full}}-\chi_{\mathrm{cut}}^{(S)}
```

The condition $\chi\_\mathrm{full}=0$ is also the outcome of perfect EIT, and by itself does not imply a no-go. By contrast, $\delta\chi\_S=0$ states that the observed response is unchanged whether the sector $S$ is retained or cut.

Splitting the weak-probe linear response into the fast optical variables $x$ and the long-lived variables $s$,

```math
\begin{pmatrix}
A & B\\
C & G_g
\end{pmatrix}
\begin{pmatrix}
x\\s
\end{pmatrix}
=
\begin{pmatrix}
b_p\\0
\end{pmatrix}
```

Eliminating the fast variables, the Schur complement of the long-lived sector is

```math
S_g=G_g-CA^{-1}B
```

and the response difference through the sector is given by

```math
\delta\chi_S
=\chi_0 d^{\dagger}A^{-1}B\,S_g^{-1}CA^{-1}b_p
```

For a closed two-level ground-state block this reduces to the form of Theorem 2A,

```math
\Xi
=S_1-\frac{\beta K_{12}K_{21}}
{\gamma_g+\beta S_2},
\qquad
\beta=\frac{|\Omega_c|^2}{4}
```

where

```math
S_i=\sum_j\frac{|d_{ij}|^2}{a_j},
\qquad
K_{12}=\sum_j\frac{d_{1j}^{*}d_{2j}}{a_j}
```

and $K\_{12}$ and $K\_{21}$ are coherent sums of the Raman amplitudes over all excited branches. If

```math
\delta\chi_S\equiv0
\quad\Longleftrightarrow\quad
K_{12}K_{21}\equiv0
```

holds over the whole region of interest, the system is an exact no-go.

Three statements must be kept apart, and conflating them is the most common error in reading this theory.

| Type | Defining condition | What it asserts |
|---|---|---|
| Exact no-go | $\delta\chi\_S\equiv0$ on the regular domain; for a closed two-level ground block, $K\_{12}K\_{21}\equiv0$ | the transfer is algebraically zero. A selection rule, not a small number |
| Asymptotic no-go | $\delta\chi\_S=O(\Gamma^{-\nu})$ with $0<\nu<\infty$ as $\Gamma\to\infty$ | the pathway exists and is nonzero at finite $\Gamma$, but is suppressed as an integer power of the dissipation |
| Practical no-go | $K\_{12}K\_{21}\neq0$, yet the engineerable, asymptotic, and detectable parameter windows have empty intersection | a property of a measurement budget, not an algebraic identity |

Only the first is a transfer zero. Note also that $\chi\_\mathrm{full}=0$ is none of the three. The third type appears as "experimental no-go" in the proof package and as "practical no-go" in the claim files; the two names denote the same category.

The NV center is placed as follows. At $B\_\perp=0$ the spin-$\Lambda$ channel is an asymptotic no-go, not an exact one: in the physical six-level $^3E$ manifold, transverse spin-orbit coupling, nonsecular spin-spin terms, and hyperfine coupling generically break the spin sector, so the exact statement holds only in a simplified sector-preserving model. With phonon-induced orbital hopping supplying full-rank fast damping $\Gamma(T)$ on the optical-coherence block, the first moment vanishes and

```math
M_0=0,\qquad K=O(\Gamma^{-2}),\qquad
\delta\Xi_S=O(\Gamma^{-3})\ \text{or}\ O(\Gamma^{-4})
```

the branch depending on $|\beta S\_2/\gamma\_g|$. A transverse field does not change the type; it opens the vertex quadratically, giving $K\_{12}K\_{21}=O(\lambda^2\Gamma^{-4})$ and a tunable escape rather than a new class of no-go.

The room-temperature conclusion of this note is therefore a statement of the third kind. It is a numerical practical no-go, not a claim that the response is exactly zero.

Simply thinking that "dephasing washes out the line" explains neither the low-temperature boundary, nor the transverse-field threshold, nor the sign reversal of the response. Figure 1 shows the selection rules of the pathway, and Figure 2 shows how the transparency collapses with temperature and eventually reverses into increased absorption.

![Figure 1: NV$^-$ spin-$\Lambda$ channel and a scale comparison of the optical-coherence decay](https://raw.githubusercontent.com/Risei412/NV-EIT-space/main/results/figures/fig1_level_scheme.png)

![Figure 2: Absorption spectra at $B_\perp = 0.232\,\mathrm{T}$ (30 K / 70 K / 105 K)](https://raw.githubusercontent.com/Risei412/NV-EIT-space/main/results/figures/fig2_spectra.png)

---

## 3. Question 1: Why Is Observation Hard at Room Temperature?

The reason is that rising temperature destroys the optical coherence inside each individual NV center, so that the interference pathway itself disappears. This is not merely a question of ensemble averaging, and it is not recovered by homogenizing the sample or by applying echo sequences.

At $B\_\perp=0.232\,\mathrm{T}$ the contrast falls from $0.99$ at $30\,\mathrm{K}$ to $1.4\times10^{-2}$ at $70\,\mathrm{K}$, $1.5\times10^{-4}$ at $90\,\mathrm{K}$, and $2.1\times10^{-5}$ at $100\,\mathrm{K}$. Beyond that the sign reverses, and only $1.1\times10^{-9}$ remains at room temperature.

Converting the contrast into an experimental signal requires the same detection chain and the same optical-density (OD) reference at every temperature. Because the optical linewidth changes with temperature, the estimated OD for the same sample can vary from about $0.04$ near room temperature to about $16$ at $10\,\mathrm{K}$. The present calculation compares at a matched $\mathrm{OD}\_\mathrm{sector}=1$.

| $T$ | Transmission signal | Integration time for SNR 5 |
|---|---:|---:|
| $70\,\mathrm{K}$ | $1.4\times10^{-2}$ | $\mathrm{\mu s}$ |
| $90\,\mathrm{K}$ | $1.5\times10^{-4}$ | $10\,\mathrm{ms}$ |
| $105\,\mathrm{K}$ | $-1.2\times10^{-5}$ | $1.6\,\mathrm{s}$ |
| $110\,\mathrm{K}$ | — | $5.6\,\mathrm{s}$ |
| above $120\,\mathrm{K}$ | — | unreachable |
| $300\,\mathrm{K}$ | $1.1\times10^{-9}$ | nine orders below the detection limit |

The sign reversal near $100\,\mathrm{K}$ is therefore a measurable target, whereas the residual room-temperature signal cannot be rescued by integration.

---

## 4. Question 2: Where Is the Theoretical Ceiling Temperature?

The region in which EIT appears is not a simple "colder is better" region, but a closed island in the $(T,B\_\perp)$ plane.

| Boundary | Location | Main mechanism |
|---|---:|---|
| Field side | $B\_\perp\gtrsim0.15\,\mathrm{T}$ | without a transverse field the main Raman pathway is closed by symmetry |
| Low-temperature side | below about $22\,\mathrm{K}$ | the linewidth narrows and, for a fixed control field, the system crosses over to Autler–Townes splitting (ATS); the two are separated by the model-comparison test of [C10] |
| High-temperature side | $90$–$95\,\mathrm{K}$ | phonon-driven orbital hopping destroys the Raman pathway |

The Monte Carlo 68% intervals at $B\_\perp=0.232\,\mathrm{T}$ are as follows.

| Quantity | Value |
|---|---:|
| $T\_{1\mathrm{pct}}$ (1% contrast threshold) | $74\,[70,78]\,\mathrm{K}$ |
| $T\_\mathrm{sign}$ | $103\,[97,108]\,\mathrm{K}$ |

The hopping rate $k\_\mathrm{orb}(T)$ is taken from the measured photophysics of single NV centers [S06]. The boundaries obtained from four different phonon-rate models differ by $0.8$–$2.7\,\mathrm{K}$, which lies inside the statistical band of about $9\,\mathrm{K}$.

![Figure 3: Classification of the $(T, B_\perp)$ plane (nine-level Liouvillian)](https://raw.githubusercontent.com/Risei412/NV-EIT-space/main/results/figures/fig3_phase_diagram.png)

Figure 4 was generated with `quick=True`, so the full-run values in the table above are cited rather than the numbers in the figure.

![Figure 4: Temperature dependence of the signed contrast (generated with quick=True; numbers not citable)](https://raw.githubusercontent.com/Risei412/NV-EIT-space/main/results/figures/fig4_contrast_vs_T.png)

---

## 5. Question 3: Why Are Alkali Atoms Visible at Room Temperature?

The difference is not the simple contrast that "alkali atoms have no inhomogeneous broadening while NV centers have phonons." Inhomogeneous broadening degrades EIT as well. What matters is whether the pathway is open to begin with, and how temperature acts on that pathway.

|  | Alkali vapor | NV$^-$ |
|---|---|---|
| $\Lambda$ pathway | both legs terminate on the same pair of ground levels and are open in principle | at $B\_\perp=0$ the two legs terminate on orthogonal spin states and the main Raman vertex vanishes |
| Main degradation with temperature | ensemble averaging such as Doppler broadening; part of it can be cancelled by the geometry | homogeneous decoherence inside each NV center; the Raman pathway itself is damped |
| Extra condition to open the pathway | in principle none | $B\_\perp\gtrsim0.15\,\mathrm{T}$ required |

This organization is based on the standard three-level EIT literature [C01], on Mishina et al. [C04] treating multiple excited levels and Doppler broadening, and on [C09] treating inhomogeneous broadening in solids. Since Theorem 2A of the present theory yields the susceptibility corresponding to the $\Lambda$ block of Mishina Eq. (2.11), the two are connected at the level of the susceptibility. This is not, however, a numerical comparison of alkali vapor and NV centers under matched conditions.

---

## 6. Question 4: What Is the Transverse Field Doing?

$B\_\perp$ mixes the ground spin sublevels and creates a spin overlap that the optical dipole alone cannot supply. This opens the Raman pathway that is closed at $B\_\perp=0$. On the weak-field side,

```math
C\simeq C_{\mathrm{res}}+aB_{\perp}^{\,n},\qquad n\approx2
```

is expected, and $n=2.11\pm0.08$ was obtained at $85\,\mathrm{K}$. At $55$–$70\,\mathrm{K}$, by contrast, the exponent moves with the fit window and is not determined as a single power law. Moreover, the dynamic range needed to measure the exponent experimentally cannot be secured, so the exponent itself is not taken as a verification target.

What matters in practice is that once $B\_\perp$ crosses the threshold it opens the EIT window, while barely moving the temperature ceiling.

| $B\_\perp$ | $T\_{1\mathrm{pct}}$ | $T\_\mathrm{sign}$ | State |
|---:|---:|---:|---|
| $\le0.10\,\mathrm{T}$ | $27$–$33\,\mathrm{K}$ | $27$–$33\,\mathrm{K}$ | no EIT window |
| $0.15\,\mathrm{T}$ | $70.3\,\mathrm{K}$ | $98.5\,\mathrm{K}$ | EIT window present |
| $0.50\,\mathrm{T}$ | $69.7\,\mathrm{K}$ | $103.8\,\mathrm{K}$ | EIT window present |

Increasing the field from $0.15$ to $0.50\,\mathrm{T}$ changes $T\_{1\mathrm{pct}}$ by only $0.6\,\mathrm{K}$. The transverse field is a switch for entering the inside of the island, not a dial for recovering high-temperature EIT.

![Figure 5: $B_\perp$ dependence of the sector contrast and the exponent versus the upper fit cutoff](https://raw.githubusercontent.com/Risei412/NV-EIT-space/main/results/figures/fig5_bperp_scaling.png)

---

## 7. What to Measure in an Experimental Test

The decisive target is the sign reversal near $103\,\mathrm{K}$. A dip at a single temperature can also arise from ATS, but a sign reversal under a temperature sweep is a discriminator specific to the present mechanism.

Before measuring, calibrate the sample-specific transverse strain $\delta$, the control Rabi frequency $\Omega\_c$, the local temperature, and the OD, and recompute the prediction with those values. The following order of measurement is recommended.

| Stage | Condition | What to look at |
|---|---|---|
| E1 | $30\,\mathrm{K}$, fixed field | confirm that a large EIT signal appears |
| E2 | $70\,\mathrm{K}$, $B\_\perp=0$–$0.3\,\mathrm{T}$ | confirm that no window exists below $0.10\,\mathrm{T}$ and that it opens above $0.15\,\mathrm{T}$ |
| E3 | $30/50/70/90\,\mathrm{K}$, matched OD | follow the four-and-a-half decade collapse of the contrast |
| E4 | fine sweep over $95$–$115\,\mathrm{K}$ | see whether the fringe crosses zero and reverses near $103\,\mathrm{K}$ |

In E4, record not only the sign but also the spectral shape and the $\Omega\_c$ dependence, and distinguish the interference window from the dressed-state doublet. E3 and E4 use the same readout scheme and the same detection chain.

![Figure 6: Fractional transmission and fluorescence signals for an OD-matched sample, and the required integration times](https://raw.githubusercontent.com/Risei412/NV-EIT-space/main/results/figures/fig6_observables.png)

---

## 8. Scope, Open Items, and Reproduction

The claims of this note presuppose a finite-dimensional Markovian GKSL generator, weak-probe linear response, a stationary rotating frame, passive response, a unique steady state, and a specified long-lived coherence sector. Strong probes, non-Markovian memory, many-body and propagation effects, gain media, and time-dependent feedback are not treated.

The reduced kernel agrees with the full Liouvillian to within 10% at only 47 of 108 points, and differs even in sign at 13 points. For this reason the full Liouvillian values are used for the contrasts and the main boundaries in the text. Adding the singlet or the $^{14}$N hyperfine manifold can change the peak contrast by up to 30%, so the absolute values carry a structural systematic error of this order.

The open items are five: (1) a full-only scan over $95$–$115\,\mathrm{K}$, (2) the shift of the sign-reversal temperature due to the singlet and the hyperfine structure, (3) recomputation with the $\delta$ and $\Omega\_c$ of a real sample, (4) regeneration of Figure 4 with `quick=False`, and (5) whether ground-state transverse strain can open a different $\Lambda$ pathway in practice. On strain in particular, the excited-state transverse strain $\delta$ held fixed in this note changes the dissipation rates, but does not open the main Raman vertex the way $B\_\perp$ does.

| Kind | Location |
|---|---|
| Figures | `results/figures/` |
| Numerical tables | `results/tables/` |
| Gate certificates | `results/certificates/` |
| Boundaries of the claims | `CLAIMS.md`, `NON_CLAIMS.md` |
| Assumptions and constraints | `theory/ASSUMPTIONS.md`, `theory/LIMITATIONS.md` |
| Theory | `theory/proofs/eit_nogo_lecture.tex` |
| Manuscript | `manuscript/main.tex` |

The current calculation suite passes 70 tests under `pytest calculations/tests`.

---

## References

Keys follow the repository's literature index (`literature/EIT_general_theory_2026-07-13/`, `literature/references.bib`).

### EIT theory and comparison systems

- **[C01]** M. Fleischhauer, A. Imamoglu, J. P. Marangos, "Electromagnetically induced transparency: Optics in coherent media," *Rev. Mod. Phys.* **77**, 633 (2005). [doi:10.1103/RevModPhys.77.633](https://doi.org/10.1103/RevModPhys.77.633) — reference treatment of three-level EIT, the linear susceptibility, and dark-state polaritons.
- **[C04]** O. S. Mishina *et al.*, "Electromagnetically induced transparency in an inhomogeneously broadened Lambda transition with multiple excited levels," *Phys. Rev. A* **83**, 053809 (2011). [doi:10.1103/PhysRevA.83.053809](https://doi.org/10.1103/PhysRevA.83.053809) — multiple excited levels plus Doppler broadening in alkali vapor; Eq. (2.11) is the susceptibility connected to Theorem 2A in §5.
- **[C09]** H. Q. Fan *et al.*, "Electromagnetically induced transparency in inhomogeneously broadened rare-earth-doped solids," *Phys. Rev. A* **99**, 053821 (2019). [doi:10.1103/PhysRevA.99.053821](https://doi.org/10.1103/PhysRevA.99.053821) — the same degradation problem in a solid host.
- **[C10]** P. M. Anisimov, J. P. Dowling, B. C. Sanders, "Objectively discerning Autler-Townes splitting from electromagnetically induced transparency," *Phys. Rev. Lett.* **107**, 163604 (2011). [doi:10.1103/PhysRevLett.107.163604](https://doi.org/10.1103/PhysRevLett.107.163604) — the model-comparison criterion used to label the EIT/ATS boundary in Figure 3.

### NV center structure and phonon physics

- **[S03]** A. Batalov *et al.*, "Low temperature studies of the excited-state structure of negatively charged nitrogen-vacancy color centers in diamond," *Phys. Rev. Lett.* **102**, 195506 (2009). [doi:10.1103/PhysRevLett.102.195506](https://doi.org/10.1103/PhysRevLett.102.195506)
- **[S04]** K.-M. C. Fu *et al.*, "Observation of the dynamic Jahn-Teller effect in the excited states of nitrogen-vacancy centers in diamond," *Phys. Rev. Lett.* **103**, 256404 (2009). [doi:10.1103/PhysRevLett.103.256404](https://doi.org/10.1103/PhysRevLett.103.256404) — orbital averaging by phonons; source of the $\Gamma\_{xy}$ convention discussed in §2.
- **[S05]** M. L. Goldman *et al.*, "Phonon-induced population dynamics and intersystem crossing in nitrogen-vacancy centers," *Phys. Rev. Lett.* **114**, 145502 (2015). [doi:10.1103/PhysRevLett.114.145502](https://doi.org/10.1103/PhysRevLett.114.145502)
- **[S06]** J. Happacher *et al.*, "Temperature-dependent photophysics of single NV centers in diamond," *Phys. Rev. Lett.* **131**, 086904 (2023). [doi:10.1103/PhysRevLett.131.086904](https://doi.org/10.1103/PhysRevLett.131.086904) — the measured rates behind $k\_\mathrm{orb}(T)$.
- **[S07]** M. W. Doherty *et al.*, "The nitrogen-vacancy colour centre in diamond," *Phys. Rep.* **528**, 1 (2013). [doi:10.1016/j.physrep.2013.02.001](https://doi.org/10.1016/j.physrep.2013.02.001)
- **[S08]** G. Thiering, A. Gali, "*Ab initio* calculation of spin-orbit coupling for an NV center in diamond exhibiting dynamic Jahn-Teller effect," *Phys. Rev. B* **96**, 081115 (2017). [doi:10.1103/PhysRevB.96.081115](https://doi.org/10.1103/PhysRevB.96.081115)

### EIT and CPT demonstrated in diamond

- **[S01]** V. M. Acosta *et al.*, "Electromagnetically induced transparency in a diamond spin ensemble enables all-optical electromagnetic field sensing," *Phys. Rev. Lett.* **110**, 213605 (2013). [doi:10.1103/PhysRevLett.110.213605](https://doi.org/10.1103/PhysRevLett.110.213605)
- **[S02]** C. Santori *et al.*, "Coherent population trapping of single spins in diamond under optical excitation," *Phys. Rev. Lett.* **97**, 247401 (2006). [doi:10.1103/PhysRevLett.97.247401](https://doi.org/10.1103/PhysRevLett.97.247401)
