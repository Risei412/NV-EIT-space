# Convention table (v2.1 / audit package; authoritative)

| Item | Convention |
|---|---|
| Hamiltonian units | GHz, ordinary frequency |
| Resolvent | G = [γ·I + i·2π·(H − z·I)]⁻¹ ; the single 2π lives here |
| Zeeman | γ_e = 28.02495164 GHz/T ; H += γ_e(B_x S_x + B_y S_y + B_z S_z), B in Tesla |
| Excited fine structure | D_par=1.42, λ_par=5.33, D_perp=1.55/2, λ_perp=0.20/√2 GHz (operator forms in nv_model.py) |
| k_orb | one-way X↔Y jump rate = Γ_XY (Happacher Eq. 23), Hz |
| Γ_pop | population-imbalance rate = 2·k_orb |
| γ_oc | phonon damping of each optical coherence = k_orb/2 = Γ_pop/4 (consistent with theory v6.2 §5.2 corollary γ_oc = Γ_XY^{imb}/4 where Γ_XY^{imb}=Γ_pop) |
| γ used in G | γ = γ_oc + GRAD/2, GRAD = 0.0668 GHz radiative |
| Optical legs | Model A1: d = |orb⟩ ⊗ |dressed ground eigvec⟩ (field-dressed, NOT bare spin) |
| K12, K21 | K12 = ⟨d_p|G|d_c⟩, K21 = ⟨d_c|G|d_p⟩ (computed independently; G non-symmetric under †) |
| β | (2π·Ω_c)²/4 with Ω_c in GHz ordinary |
| γ_g,eff | 2γ_g + 2×10⁻⁶ GHz |
| Contrast C | C = −Im(δΞ·i)/Im(S₁·i) evaluated as dA/A_cut; C > 0 = transparency (EIT-like), C < 0 = absorption-adding |
| Warning | v6 numerical package not supplied to this audit; any v6-derived numbers must be converted to the above before reuse |
