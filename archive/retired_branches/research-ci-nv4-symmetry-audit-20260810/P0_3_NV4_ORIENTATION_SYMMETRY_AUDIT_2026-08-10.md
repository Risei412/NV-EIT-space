# P0-3 — NV four-orientation symmetry audit

Date: 2026-08-10  
Run: RCI-20260810-PRL14  
Official baseline: RCI-20260810-PRL12  
Interactive precursor: PRL13 pilot/no-write  
Base repository state: `e5273393fe6fcfb7962c8ac37d014a0075b04126`  
Branch: `research-ci/nv4-symmetry-audit-20260810`  
Literature mode: web only  
Status: **PRL NOVELTY KILL / informative internal result retained**

## 1. Internal anomaly

The existing Gate-5 ensemble code already contains the four tetrahedral NV axes

\[
(1,1,1),\;(1,-1,-1),\;(-1,1,-1),\;(-1,-1,1)
\]

but used them only to quantify ensemble washout of the 70 K response.  It did not ask what the orientation set forbids or permits as a response-selection object.

The internal anomaly is that the real NV optical model is not an axis-only scalar response.  The \(^3E\) orbital doublet, transverse strain and optical polarization mark a local transverse frame.  Therefore equal occupation of the four crystallographic axes does not automatically imply a full tetrahedral average of the measured response.

## 2. Audit observable

To avoid branch-index switching at the zero-field excited-state degeneracy, the audit uses a branch-order-independent raw observable

\[
R_G(\mathbf B)=\frac1{|G|}\sum_{F\in G}\sum_{j=0}^{5}
\left[\operatorname{Im}\chi^{F,j}_{\rm cut}(\mathbf B)
-\operatorname{Im}\chi^{F,j}_{\rm full}(\mathbf B)\right].
\]

The calculation retains the existing nine-level full Lindblad structure, phonon orbital hopping, radiative decay, ground relaxation/dephasing, transverse strain and Y/Y optical protocol at 70 K.  Only the local magnetic field is generalized from `(Bx,0,Bz)` to arbitrary `(Bx,By,Bz)`.

## 3. Two orientation-frame orbits

### V4: four marked frames

Take one transverse frame for each of the four NV axes, related by the four proper \(\pi\)-rotations

\[
I,\quad {\rm diag}(1,-1,-1),\quad {\rm diag}(-1,1,-1),\quad {\rm diag}(-1,-1,1).
\]

This orbit represents four crystallographic axes while retaining a definite transverse orbital/polarization mark.

The exact invariant polynomials through degree three permit

\[
1,\quad B_x^2,\;B_y^2,\;B_z^2,\quad B_xB_yB_z.
\]

Hence V4 forbids the vector term and off-diagonal quadratic terms, but it does **not** forbid diagonal traceless rank-2 anisotropy.

### T12: proper-tetrahedral frame completion

Complete each NV axis by its three \(C_3\)-related transverse frames.  The resulting 12-frame orbit is the proper tetrahedral rotation group.

Its invariant polynomials through degree three are

\[
1,\quad B_x^2+B_y^2+B_z^2,\quad B_xB_yB_z.
\]

Therefore the complete tetrahedral frame orbit forbids rank-1 response and traceless rank-2 anisotropy, while allowing an isotropic quadratic term and a tetrahedral cubic invariant.

This is an exact group-averaging statement; the numerical Lindblad calculation is a physical witness that the distinction matters for the present NV optical protocol.

## 4. Full-Lindblad witness

Finite-difference step: \(10^{-4}\,\mathrm T\).

### Four marked frames V4

The Hessian of the branch-summed raw response at \(\mathbf B=0\) is

\[
H_{V4}\approx
\begin{pmatrix}
1.05026135&0&0\\
0&0.10844234&0\\
0&0&2.04972040
\end{pmatrix}.
\]

Eigenvalues:

\[
(0.10844234,\;1.05026135,\;2.04972040).
\]

Relative eigenvalue spread: `1.81517`.  
Relative traceless-Hessian norm: `0.59544`.

**Result:** four equal NV axes alone leave a large rank-2 anisotropic response in this marked optical model.

### Twelve-frame tetrahedral completion T12

\[
H_{T12}\approx
\begin{pmatrix}
1.06947370&0&0\\
0&1.06947070&0\\
0&0&1.06948411
\end{pmatrix}.
\]

Eigenvalues:

\[
(1.06947070,\;1.06947370,\;1.06948411).
\]

Relative eigenvalue spread: `1.25e-5`.  
Relative traceless-Hessian norm: `5.37e-6`.

**Result:** completing the local \(C_3\) frame orbit removes the rank-2 anisotropy to numerical precision, exactly as required by tetrahedral invariant theory.

## 5. Generated cubic response and protocol test

The mixed cubic derivative

\[
\partial_{B_x}\partial_{B_y}\partial_{B_z}R
\]

is nonzero for the one-sided \(m_s=0\leftrightarrow +1\) protocol:

- V4, `ctrl=+1`: `-16.0280`
- T12, `ctrl=+1`: `-16.0342`
- T12, `ctrl=-1`: `+15.9559`

A symmetric average of the \(+1\) and \(-1\) control channels gives `-0.03918`, a suppression factor of about `409`, but **not an exact zero** in the present full Hamiltonian/readout model.

Therefore the stronger statement that spin-channel symmetrization exactly switches the leading tetrahedral cubic term off is **HOLD / not established**.

## 6. Candidate question

> Are the four crystallographic NV axes by themselves sufficient to enforce tetrahedral response selection in an optical NV experiment, or does a marked transverse orbital/polarization frame leave symmetry-allowed lower-rank anisotropy that disappears only after the local \(C_3\) frame orbit is completed?

Conceptual contrast:

> The same four NV axes can exhibit or forbid rank-2 anisotropy depending on whether the measurement protocol resolves a transverse frame that the axis-only description discards.

## 7. External novelty audit

Closest absorbers:

1. **Four-axis NV vector magnetometry.** Simultaneous use of the four tetrahedral NV axes as a vector coordinate system is established.  Thus “four NV orientations contain directional information” is not novel.
2. **NV \(C_{3v}\) group theory and optical selection rules.** The local symmetry, orbital structure and polarization-dependent transition rules are established.  Thus introducing a transverse optical frame is not a new primitive.
3. **Polarization-assisted orientation discrimination.** Optical polarization is already used to distinguish NV orientations and exploit transition-dipole anisotropy.  Thus “readout can mark the transverse/orientation sector” is known physics.
4. **Four-orientation tensor averaging in nonlinear NV response.** Existing NV work explicitly averages anisotropic tensors/nonlinear response over the four possible NV axes to obtain macroscopic symmetry constraints.  Thus the broad claim “four orientations select allowed tensor components” is already occupied.
5. **Finite-group / Reynolds averaging.** The V4-versus-tetrahedral invariant calculation is standard group averaging.  The rank-selection theorem is therefore not a new mathematical theorem by itself.

Primary papers checked:

- J. M. Schloss et al., *Phys. Rev. Applied* **10**, 034044 (2018), simultaneous vector magnetometry using all four NV orientations.
- F. Münzhuber et al., *Phys. Rev. Applied* **14**, 014055 (2020), polarization-assisted vector magnetometry with NV ensembles.
- J. R. Maze et al., *New J. Phys.* **13**, 025025 (2011), group-theoretic description of NV states and selection rules.
- Recent single-NV polarization-anisotropy work assigning crystallographic orientations from polarization-dependent optical response.
- M. Ichikawa et al., *Nature Communications* **15**, 7174 (2024), nonlinear/tensor response constrained by averaging over the four possible NV axes.

## 8. Novelty residue

The internal result that survives the broad collision is **not** a new tetrahedral selection theorem.  It is the physically useful distinction

> `axis orbit` \(\neq\) `frame orbit` when the NV optical protocol carries a transverse orbital/strain/polarization mark.

The present calculation shows that replacing the full marked-frame problem by an axis-only tetrahedral average can create a false rank-2 selection rule.

However, after external audit this residue remains a natural physical specialization of known NV optical anisotropy plus standard group averaging.  No independent theorem, new observable class, or previously excluded physical regime has yet been established.

Novelty confidence for **PRL-level novelty**: **High confidence KILL**.

## 9. Gate ledger

- Gate 0 Internal Seed Mining: **PASS**. Existing four-orientation Gate-5 implementation contained an unused symmetry question.
- Gate 1 External Collision: **PASS WITH ONE NARROWING**. Broad “four orientations generate tensor selection” is killed; narrowed axis-orbit versus marked-frame-orbit distinction retained for audit.
- Gate 2 Question Synthesis: **PASS**.
- Gate 3 Preliminary Novelty Kill: **FAIL: known mathematics / natural embedding**.
- Gate 4 Minimal Witness: executed for information only; **witness verified** for V4 rank-2 leakage and T12 rank-2 cancellation. This does not rescue the Gate-3-failed PRL candidate.
- Alternate seed, exact spin-channel rank switch: **HOLD / witness not exact**, so it is not promoted.
- Gates 5–8: **not opened for the killed PRL candidate**.

## 10. Final verdict

**FAIL: PRL novelty.**

Do not create a new active PRL Research Question from the broad four-orientation selection claim and do not rescue it by renaming it.  Retain this branch as a negative/design-rule result: four crystallographic axes do not automatically realize full tetrahedral frame averaging when the optical protocol marks the transverse frame.

A genuinely new future seed would require a fresh internal anomaly beyond group averaging, for example a non-factorizable NV-internal selection rule that survives after the full frame symmetry has already been accounted for.  The present \(m_s=\pm1\) cubic suppression is not yet such a result because the residual is nonzero.
