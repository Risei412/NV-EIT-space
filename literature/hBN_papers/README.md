# hBN literature — bibliographic index

This directory previously held publisher PDFs of the hexagonal-boron-nitride
(hBN) papers used as parameter input for the no-go analysis. The PDFs are **no
longer tracked in git**: this repository is public, and the articles are
copyrighted by their publishers.

What is kept instead is everything the analysis actually depends on — the full
bibliographic record, the DOI/arXiv identifier, and the SHA-256 of the exact
file each number was read from, so any result can be traced back to a specific
version of a specific paper.

## How to restore the PDFs locally

Download each article from its DOI below into this directory using the file
name in the first column. `.gitignore` excludes `hBN papers/*.pdf`, so local
copies will never be committed. Verify with:

```bash
sha256sum "hBN papers"/*.pdf
```

and compare against the `sha256` column of
`../No-go theorem/results/tables/literature_manifest.csv`, which remains the
machine-readable source of truth for this table.

## Related documents

- `../No-go theorem/results/literature_notes.md` — how each paper feeds the
  no-go argument, grouped by role (temperature/linewidth, room-temperature spin
  coherence, high-coherence emitters without an established Λ, defect
  identification).
- `../No-go theorem/results/hbn_nogo_EIT_report.md` — the analysis these
  parameters were used for.
- `../EIT_general_theory_literature_2026-07-13/references.bib` — BibTeX entries
  for the wider literature survey.

## Index

| File name | Reference | Year | Venue | DOI / ID | Role in the analysis |
|---|---|---|---|---|---|
| `Akbari_2022_lifetime_limited_tunable_hBN.pdf` | Akbari et al., *Lifetime-Limited and Tunable Quantum Light Emission in h-BN via Electric Field Modulation* | 2022 | Nano Letters 22, 7798–7803 | [10.1021/acs.nanolett.2c02163](https://doi.org/10.1021/acs.nanolett.2c02163) | Electric-field suppression of spectral diffusion; 6.5–120 K linewidth |
| `Bhat_2026_hybrid_jump_diffusion_hBN.pdf` | Bhat et al., *A Hybrid Jump-Diffusion Model for Coherent Optical Control of Quantum Emitters in hBN* | 2026 | arXiv preprint | arXiv:2601.20587 | Stochastic spectral diffusion/dephasing model calibrated to 5–30 K data |
| `Cholsuk_2024_hBN_defects_database.pdf` | Cholsuk et al., *The hBN Defects Database: A Theoretical Compilation of Color Centers in Hexagonal Boron Nitride* | 2024 | J. Phys. Chem. C | [10.1021/acs.jpcc.4c03404](https://doi.org/10.1021/acs.jpcc.4c03404) | DFT database of structures, lifetimes, dipoles, spin multiplicities and quantum-memory candidates |
| `Dietrich_2018_Fourier_limited_lines_hBN.pdf` | Dietrich et al., *Observation of Fourier Transform Limited Lines in Hexagonal Boron Nitride* | 2018 | Physical Review B 98, 081414(R) | [10.1103/PhysRevB.98.081414](https://doi.org/10.1103/PhysRevB.98.081414) | Cryogenic near-lifetime-limited optical transitions |
| `Dietrich_2020_room_temperature_FT_limited_hBN.pdf` | Dietrich et al., *Solid-State Single Photon Source with Fourier Transform Limited Lines at Room Temperature* | 2020 | Physical Review B 101, 081401(R) | [10.1103/PhysRevB.101.081401](https://doi.org/10.1103/PhysRevB.101.081401) | 3–300 K sub-100 MHz PLE linewidth in a mechanically decoupled emitter |
| `Gerard_2026_B_center_resonance_fluorescence.pdf` | Gérard et al., *Resonance Fluorescence and Indistinguishable Photons from a Coherently Driven B Centre in hBN* | 2026 | Nature Communications | [10.1038/s41467-026-68555-5](https://doi.org/10.1038/s41467-026-68555-5) | Coherent two-level optical drive, Mollow triplet and photon indistinguishability |
| `Hoese_2020_mechanical_decoupling_hBN.pdf` | Hoese et al., *Mechanical Decoupling of Quantum Emitters in Hexagonal Boron Nitride from Low-Energy Phonon Modes* | 2020 | Science Advances 6, eaba6038 | [10.1126/sciadv.aba6038](https://doi.org/10.1126/sciadv.aba6038) | Low-frequency phonon gap and about 61 MHz homogeneous linewidth |
| `Horder_2025_B_center_decoherence.pdf` | Horder et al., *Optical Coherence of B Center Quantum Emitters in Hexagonal Boron Nitride* | 2025 | ACS Photonics 12, 1284–1290 | [10.1021/acsphotonics.4c02088](https://doi.org/10.1021/acsphotonics.4c02088) | B-center PLE linewidths and processing-induced charge-noise decoherence |
| `Jungwirth_2016_temperature_dependence_hBN.pdf` | Jungwirth et al., *Temperature Dependence of Wavelength Selectable Zero-Phonon Emission from Single Defects in Hexagonal Boron Nitride* | 2016 | Nano Letters 16 | [10.1021/acs.nanolett.6b01987](https://doi.org/10.1021/acs.nanolett.6b01987) | Temperature shifts and broadening of visible hBN emitters |
| `Koch_2024_limits_coherent_optical_control_hBN.pdf` | Koch et al., *Probing the Limits for Coherent Optical Control of a Mechanically Decoupled Defect Center in Hexagonal Boron Nitride* | 2024 | Communications Materials 5, 240 | [10.1038/s43246-024-00686-y](https://doi.org/10.1038/s43246-024-00686-y) | Separates lifetime limit, homogeneous dephasing and spectral diffusion versus temperature |
| `Mathur_2022_excited_state_spin_VBminus_hBN.pdf` | Mathur et al., *Excited-State Spin-Resonance Spectroscopy of V_B^- Defect Centers in Hexagonal Boron Nitride* | 2022 | Nature Communications 13 | [10.1038/s41467-022-30772-z](https://doi.org/10.1038/s41467-022-30772-z) | Ground/excited spin splittings and room-temperature spin Hamiltonian |
| `Stern_2022_room_temperature_single_spin_ODMR_hBN.pdf` | Stern et al., *Room-Temperature Optically Detected Magnetic Resonance of Single Defects in Hexagonal Boron Nitride* | 2022 | Nature Communications 13 | [10.1038/s41467-022-28169-z](https://doi.org/10.1038/s41467-022-28169-z) | Single-defect room-temperature ODMR and optical/spin dynamics |
| `Stern_2024_quantum_coherent_spin_hBN.pdf` | Stern et al., *A Quantum Coherent Spin in Hexagonal Boron Nitride at Ambient Conditions* | 2024 | Nature Materials 23 | [10.1038/s41563-024-01887-z](https://doi.org/10.1038/s41563-024-01887-z) | Single carbon-related S=1 spin, T2*, echo and dynamical decoupling |
| `White_2021_phonon_dephasing_spectral_diffusion_hBN.pdf` | White et al., *Phonon Dephasing and Spectral Diffusion of Quantum Emitters in Hexagonal Boron Nitride* | 2021 | Optica 8, 1153–1158 | [10.1364/OPTICA.431262](https://doi.org/10.1364/OPTICA.431262) | 4–40 K homogeneous broadening and spectral diffusion |
| `Whitefield_2025_narrowband_spin_emitters_hBN.pdf` | Whitefield et al., *Generation of Narrowband Quantum Emitters in hBN with Optically Addressable Spins* | 2025 | arXiv preprint; Nature Materials 2026 | arXiv:2501.15341 | Room-temperature optical spin readout; S=1 and S=1/2 spin complexes |

`Dietrich_2018_Fourier_limited_lines_hBN.pdf` is listed in the manifest but was
not among the files stored in this directory; it is retained in the table
because the analysis cites it.
