# calculations/

All computation. Code lives here; its outputs live in `../results/`.

```
analytic/      symbolic / analytic derivations (Jupyter, sympy)
numerics/      the numerical campaigns
simulations/   dynamics / Monte Carlo campaigns
tests/         theorem and regression tests
```

## Layout note (read before moving anything)

Campaign code keeps its original internal structure — `<campaign>/src/`, and
in the sector-response tree one directory per gate. What changed in the
2026-08-13 restructure is that **outputs no longer live next to the code**:
every script writes into the repository-level `results/` tree

| output | destination |
|---|---|
| `.csv`, non-gate `.json` | `results/tables/` |
| `.png`, `.pdf` | `results/figures/` |
| `.npz` | `results/raw/` |
| `gate*.json`, `gates_summary_*.json` | `results/certificates/` |

so a result is found by *what it is*, not by which campaign produced it. All
filenames across campaigns are distinct, which is what makes the flat layout
safe; keep it that way when adding outputs.

## The campaigns

| Path | Contents |
|---|---|
| `numerics/No-go theorem/` | The binary NV EIT no-go theory: NV/phonon models, verification gates 1–5, the PRL prediction runner, and `scripts/reproduce_prl_figures.py`. |
| `numerics/New no-go theory/` | NV-side campaigns of the three-tier (Class I–III) classification: Phases A, B, D, M, P (`src/`), Gates B–E, Phase O observable inheritance, and the room-temperature campaign (`RoomT/`). The material-independent theory itself lives in [Sector-Master-Resolved-Theory](https://github.com/Risei412/Sector-Master-Resolved-Theory). |
| `numerics/EIT definition equivalence/` | Phase E/F: the redefinition of EIT and the conditions equivalent to it — Λ-chain, counterexamples, 2g2e boundary, 2g3e search, transparency floor, full GKSL, matched-floor accretive gate. |
| `numerics/manuscript_figures/pra/` | PRA figure and analysis pipeline (P0-2, P1–P6). |
| `numerics/manuscript_figures/prl/` | The four PRL main-text figures, assembled from the Gate A–D campaigns. Presentation only — every panel reuses the gate functions unchanged. |
| `analytic/master-equation/` | Jupyter notebooks for the master-equation formulation (6-level EIT, symbolic derivations). |
| `simulations/` | Reserved for dynamics / Monte Carlo campaigns. Currently empty: the room-temperature campaign is a parameter-sweep campaign and stayed under `numerics/New no-go theory/RoomT/`. |

## Running

Each campaign pins its own dependencies:

```bash
pip install -r "calculations/numerics/No-go theorem/requirements.txt"
python "calculations/numerics/No-go theorem/scripts/reproduce_prl_figures.py"

pip install -r "calculations/numerics/New no-go theory/requirements.txt"
python "calculations/numerics/New no-go theory/PhaseO_observable_inheritance/src/run_gate_a.py"

python calculations/numerics/manuscript_figures/prl/src/make_figures.py
```

## Tests

```bash
pytest calculations/tests            # 70 tests
```

`tests/` is organized by campaign (`nogo/`, `gate_a/` … `gate_d/`, `pra/`,
`prl_figures/`). The figure tests check each figure against the certified gate
numbers, so a figure cannot silently drift from the result it depicts.
