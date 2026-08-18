# P0-4 — Tensor-rank generation in the four-orientation NV ensemble

Date: 2026-08-10  
Parent: P0-3 NV four-orientation symmetry audit  
Base repository state: `e5273393fe6fcfb7962c8ac37d014a0075b04126`  
Branch: `research-ci/nv4-symmetry-audit-20260810`  
Status: **promising internal seed / PRL not yet passed**

## 1. Exact tetrahedral moment hierarchy

For the four directed tetrahedral axes

\[
\mathbf n_a\in\frac1{\sqrt3}\{(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)\},
\]

\[
\sum_a n_{ai}=0,
\qquad
\sum_a n_{ai}n_{aj}=\frac43\delta_{ij},
\]

but

\[
T_{ijk}=\sum_a n_{ai}n_{aj}n_{ak}\neq0.
\]

The only nonzero components of `T_ijk` are permutations of xyz, each equal to

\[
\frac{4}{3\sqrt3}.
\]

Hence

\[
\sum_a(\mathbf n_a\!\cdot\!\mathbf B)^3
=\frac8{\sqrt3}B_xB_yB_z.
\]

The four-orientation set is therefore dipole-free and quadrupole-isotropic but octupolar.

The Molien coefficients of the proper tetrahedral group for scalar homogeneous polynomials of ranks 0...10 are

`1, 0, 1, 1, 2, 1, 4, 2, 5, 4, 7`.

Thus:

- rank 1: no invariant,
- rank 2: one invariant, `B^2`, hence no rank-2 anisotropy,
- rank 3: one invariant, `Bx By Bz`,
- rank 4: two invariants, so even anisotropy becomes possible.

## 2. Correction to P0-3

At `B=0`, the ground `ms=+1` and `ms=-1` states are degenerate.  A one-sided `+1` or `-1` channel label obtained by overlap assignment is therefore not a basis-invariant object under arbitrary transverse perturbations.

P0-3 correctly removed excited-branch ordering by summing all six excited branches, but it did not remove this ground-doublet labeling subtlety.

Therefore the one-sided cubic derivative is retained as a **protocol-resolved witness**, not as an intrinsic property of an unresolved four-orientation ensemble.

## 3. Basis-invariant spin-summed observable

Define

\[
R_{\pm}(\mathbf B)=\frac12[R_+(\mathbf B)+R_-(\mathbf B)].
\]

Using the same T12 proper-tetrahedral frame orbit and nine-level full Lindblad model as P0-3, the Hessian at zero field (`h=10^-4 T`) is

\[
H_{\pm}\approx
\begin{pmatrix}
1.00707882&0&0\\
0&1.00707544&0\\
0&0&1.00708560
\end{pmatrix},
\]

with eigenvalues

\[
(1.00707544,1.00707882,1.00708560).
\]

Hence the complete rank-2 response is isotropic to numerical precision.

Direct field-reversal tests give relative deviations between `R_pm(B)` and `R_pm(-B)` of approximately `1.7e-5`, `4.0e-5`, and `1.7e-4` for the tested 1–2 mT vectors.  The spin-summed response is therefore numerically field-even in the small-field regime.

## 4. Rank-4 generation

For equal field magnitude r, compare `[100]` and `[111]` directions.  Their isotropic rank-2 contributions are identical.  In the spin-summed nearly field-even sector, their first symmetry-allowed difference is rank 4.

Using

\[
c_4(r)=\frac{R_{\pm}(r\hat x)-R_{\pm}(r\hat n_{111})}{(2/3)r^4},
\]

we obtain

| r (T) | Delta R | c4 |
|---|---:|---:|
| 5e-4 | -1.074e-11 | -257.74 |
| 1e-3 | -1.567e-10 | -235.01 |
| 2e-3 | -2.353e-9 | -220.62 |
| 3e-3 | -1.084e-8 | -200.82 |

A small-field extrapolation using the first three points gives approximately

\[
c_4(0)\simeq -2.53\times10^2.
\]

Thus the basis-invariant spin-summed ensemble has a verified **rank-4 leading anisotropic response** after its isotropic rank-2 sector.

## 5. Protocol-resolved rank-3 sector

For one-sided spin control, the T12 cubic mixed derivative remains large:

\[
\partial_x\partial_y\partial_z R_+\approx -16.03,
\qquad
\partial_x\partial_y\partial_z R_-\approx +15.96.
\]

The robust feature is the opposite-sign O(16) cubic scale.  The tiny symmetric residual is finite-difference-step dependent and is not promoted to a microscopic symmetry-breaking claim.

The reduced-kernel model gives an essentially exact antisymmetry,

\[
d_{xyz}^{(+)}=+1.0009303647,
\qquad
d_{xyz}^{(-)}=-1.0009303659,
\]

with symmetric residual `~6e-10`.

## 6. Interpretation: protocol-selected tensor-rank accessibility

The same spatial four-orientation structure supports different leading anisotropic tensor ranks depending on what the protocol resolves:

- orientation/spin-resolved odd sector: tetrahedral rank 3 (`Bx By Bz`) can be exposed;
- spin-summed, field-even sector: odd rank is removed, rank 2 is isotropic, so rank 4 is the first anisotropic sector.

This is best viewed as **tensor-rank accessibility selected jointly by spatial symmetry and protocol parity**, not as a new tetrahedral invariant theorem.

A parity caveat is essential: magnetic field is an axial vector.  A cubic contraction with the polar tetrahedral octupole has pseudoscalar character under improper operations.  Therefore a physical rank-3 scalar signal requires the full measurement/control protocol to supply the corresponding handedness or to reduce the effective symmetry from an achiral tetrahedral group to its proper-rotation subgroup.  The present one-sided spin-channel protocol is a candidate realization, but that identification is not yet proved as an exact symmetry theorem.

## 7. Novelty audit

### Killed components

1. Tetrahedral order described by a symmetric rank-3/octupolar tensor is established physics.
2. Nonlinear responses as probes of octupolar order are established in liquid-crystal, multipolar and nonlinear-optical contexts.
3. Four NV orientations as a vector-sensing basis are established.
4. Optical polarization can resolve NV orientation and transverse optical structure.
5. Averaging tensor response over four NV orientations is established in NV-related nonlinear-response work.

### Surviving residue

No direct primary-literature match was found in the targeted search for the specific NV statement:

> the same four-orientation NV ensemble admits a protocol-resolved tetrahedral rank-3 magnetic sector while its basis-invariant spin-summed small-field response suppresses rank-1 and anisotropic rank-2 and exposes rank-4 as the first anisotropic sector.

This remains **HOLD**, not PASS, because it may still be a natural application of standard invariant projection plus ordinary spin-channel selection.

## 8. Decisive next gate

Do not spend the next gate proving group theory already known.  The decisive question is whether rank selection creates a capability unavailable to ordinary four-axis vector reconstruction.

The strongest candidate is a matched-pair sensing test:

> Can two magnetic-field configurations be identical in all rank-0, rank-1 and rank-2 ensemble invariants accessible to a conventional unresolved readout, yet be distinguished with finite contrast by the rank-3 protocol channel or rank-4 spin-summed channel?

Required witness:

1. construct a matched pair with identical lower-rank invariants,
2. show conventional total/second-order readout is degenerate,
3. show a rank-3 or rank-4 channel separates the pair in the full Lindblad model,
4. show separation survives strain, hyperfine, four-orientation population imbalance and realistic polarization geometry,
5. audit against nonlinear NV magnetometry and octupolar sensing literature.

If this witness exists, the claim changes from “tetrahedral group theory appears in NV” to “NV symmetry provides an operational high-rank discriminator unavailable to lower-rank sensing.”  That is the route with plausible PRL leverage.
