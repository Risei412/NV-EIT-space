# 材料横断の検証ワークフロー

## Stage 0: 対象channelの固定

- lower states `|g1>, |g2>`
- operative excited manifold `E`
- probe/control polarizationと伝搬方向
- electric dipole / magnetic dipole / Raman / cavity-mediatedの別
- resolved、partially merged、fully mergedの温度領域

曖昧なまま `eta_geom` を計算しない。

## Stage 1: Symmetry pre-screen

1. point groupとdouble groupを確定
2. lower/excited irrepsを割り当て
3. dipole irrepとpolarization tensorを構築
4. spin overlapを判定
5. `Gamma(g1)^* x Gamma(d†Pi_Ed) x Gamma(g2)` がtotally symmetric irrepを含むか判定

出力:
- exact symmetry zero
- allowed but reduced matrix element未確定
- polarization-selective go

## Stage 2: Microscopic kernel

`H_e`, dipole matrices, strain, magnetic field, spin-orbit, spin-spin, Jahn-Tellerを含め、

`K(omega,T)=p†[omega I-H_e-Sigma_e]^{-1}c`

を計算する。複素位相を保持する。

出力:
- `K(omega,T)`
- `eta_geom`
- singular values of coupling map
- selected-subspace `eta_S`
- non-scalar correction `p†delta M c`

## Stage 3: Liouvillian benchmark

full Lindblad/Redfield modelで、次を含める。

- spontaneous emission
- pure dephasing
- phonon transitions
- intersystem crossing
- population pumping/shelving
- ground-state relaxation
- branch-dependent linewidth

弱probe kernelと一致する領域を確認する。

## Stage 4: Observable definition

最低限、次を別々に出す。

- baseline absorption `A_0`
- control-on absorption `A_c`
- EIT contrast `C=(A_0-A_c)/A_0`
- linewidth
- delay/dispersion
- dark-state fidelity
- `eta_geom`
- ATS/EIT model score

## Stage 5: 材料分類

1. **Symmetry no-go**: leading overlapが厳密zero
2. **Asymptotic no-go**: merged leading term zero、有限補正あり
3. **Partial EIT**: overlap非零だがcoupling vectors非平行
4. **Exact dark-state**: coupling vectorsがcollinear
5. **Tunable**: field/strain/cavity/polarizationでclassが変わる
6. **Practical no-go**: geometryはgoだがcoherence/power budgetを満たさない

## 優先適用順

1. SiC divacancy — 既報EITとLindblad比較があり、外部検証に最適
2. SiV — 微視的HamiltonianとCPTデータが豊富
3. Rare-earth — empirical goとして理論の一般性を確認
4. PbV/SnV/GeV — group-IV内の予測比較
5. Quantum dot/GaAs — 異なる固体platformへの外挿
6. hBN — level identification不足のため探索的
