# Dark-State-Free EIT 2g+2e Numerical Package

## Purpose

This package contains the analytical summary, numerical code, raw data, and figures for the minimal \(2g+2e\) Markovian model showing:

- full-rank optical coupling with no optical dark state;
- absence of a stationary pure Lindblad dark state;
- strong dark-state-free coherent transparency;
- pole/residue exclusion of Autler–Townes splitting;
- a determinant no-go theorem for a perfect real-axis zero;
- Class I–III scaling under different dissipation paths;
- a local mapping to the NV center.

## Directory structure

- `docs/`: Japanese summary and English report
- `code/`: reproduction script and requirements
- `data/`: spectrum and dissipation-scaling CSV files
- `figures/`: absorption, pole-zero, and Class I–III plots
- `metadata/`: parameters and SHA-256 checksums

## Reproduction

```bash
python -m pip install -r code/requirements.txt
python code/dark_state_free_2g2e_analysis.py
```

The central result is a strong dark-state-free coherent-transparency dip together with a proof that a perfect regular real-axis zero is forbidden for full-rank \(2g+2e\) optical coupling.
