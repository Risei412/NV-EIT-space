# Group-IV (SiV⁻ / SnV⁻) parameter provenance

`group_iv_model.py` referred to this file for the SiV⁻/SnV⁻ constants, but it did
not exist in the tree, so the numbers had no traceable source outside a docstring
and both `group_iv_full.py` and the Gate C README carried a note saying so. This
file records where each constant comes from. **No value is changed here** — this
is provenance only, written so a reader can check the inputs rather than take
them on trust.

## Why these values are guide lines, not the result

Gate C's claim is an *exponent*: the group-IV orbital-Λ kernel falls as Γ⁻¹
because `M0 ≠ 0`. That exponent is independent of the absolute phonon rate, which
the gate sweeps over four decades, and independent of the overall prefactor. The
constants below set where on the Γ axis the physical operating point sits and
what the plateau value is; they do not set the class. A reader who disagrees with
a factor of two in Δ_e should still get slope 1.

## Constants

| Symbol | Value | Where used | Source |
|---|---|---|---|
| SiV⁻ Δ_e (excited spin-orbit splitting) | 255 GHz | `group_iv_full.PARAMS["SiV"]`, `group_iv_model.PARAMS` | Pingault *et al.*, Nat. Commun. **8**, 15579 (2017) |
| SiV⁻ ground spin-orbit splitting | ≈48 GHz | context only (not a model input) | Pingault *et al.* (2017); Jahnke *et al.*, New J. Phys. **17**, 043011 (2015), arXiv:1411.2871 |
| SiV⁻ ground orbital relaxation | T₁ ≈ 39 ns @ 5 K → 0.026 GHz | `PARAMS["SiV"]["ground_orbital_relax"]` | Jahnke *et al.* (2015) |
| SiV⁻ excited radiative lifetime | ≈1.7 ns → 0.59 GHz | `PARAMS["SiV"]["excited_radiative"]` | Jahnke *et al.* (2015) |
| SnV⁻ Δ_e | 3000 GHz | `group_iv_full.PARAMS["SnV"]` | Trusheim *et al.*, PRX **11**, 041041 (2021) |
| SnV⁻ ground spin-orbit splitting | ≈850 GHz | context only | Trusheim *et al.* (2021) |
| Radiative damping `GAMMA_RAD` | 0.0157 GHz (HWHM) | `group_iv_model.GAMMA_RAD` | Meesala *et al.*, PRB **97**, 205444 (2018); consistent with the ~1.7 ns lifetime above |

SnV⁻ has no entry for the ground orbital relaxation or the excited radiative
rate (`None` in `PARAMS`): the gate does not use them, and rather than carry an
invented number the field is left empty.

## What is *not* pinned down

The docstring of `group_iv_full.py` and §12.3 of the split-strategy document both
flag that the phonon normalization and the dipole geometry of the group-IV model
remain schematic. That is still true and this file does not change it. What is
schematic:

- the phonon spectral density is represented by a single rate scale swept over
  decades, not by a measured density of states;
- the dipole legs are taken as orthogonal orbital basis vectors rather than from
  a computed transition-dipole geometry;
- the Gate C figures normalize `|M0|` to 1, so the plateau height is a
  normalized quantity and no absolute prefactor is claimed.

None of these affect the integer exponent, which is what Gate C certifies, but
they do mean the group-IV curves should not be read as quantitative predictions
of an absolute response.
