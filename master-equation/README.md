# Master Equation

Jupyter notebooks for the density-matrix formulation underlying the EIT
analysis. These are the hand-derivation and symbolic-algebra stage that the
numerical campaigns in `../No-go theorem/` and `../New no-go theory/` build
on — exploratory rather than part of the verification pipeline.

## Notebooks

- **`EIT in 6 Levels.ipynb`** — how the EIT response changes when the system
  has excited levels and transition pathways beyond the standard three-level
  Λ scheme. Sets up the parameter definitions, works through the density-matrix
  equations of motion by hand, and checks that the six-level result reduces to
  the familiar three-level form in the Ω_p,2 → 0 limit.

- **`Mathmatical formula.ipynb`** — the derivation chain from the model to the
  observable, kept in notebook form so the LaTeX for each step can be lifted
  directly into a manuscript without re-typesetting: the four-level model with
  phonon transitions, the density operator and dipole-moment operator in
  matrix form, the electric polarization P, and the absorption (EIT signal) as
  Im[χ] via the Kramers–Kronig relation.

## Notes

Notebook prose and annotations are in Japanese. `.ipynb_checkpoints/` is
gitignored; committed notebooks retain their outputs so the derivations can be
read without executing them.
