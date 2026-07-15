# hBN物性パラメータ文献の用途別整理

## A. 温度依存・均一幅・スペクトル拡散

- **Dietrich et al. (2020)** — 機械的脱結合エミッタで3–300 Kのsub-100 MHz PLE幅。室温conditional-go計算の最良光学入力。
- **Hoese et al. (2020)** — 約61 MHzの均一幅と低周波フォノンgap。full-rank dampingを弱める機構候補。ただし厳密な `ker(D)` の証明ではない。
- **Koch et al. (2024)** — lifetime limit、pure dephasing、spectral diffusionを分離。別エミッタでは約50 K以上で広がり、coherent driveは20–30 Kでoverdampedへ移行。
- **Akbari et al. (2022)** — 電場でspectral diffusionを抑制し、6.5–120 Kでほぼ立方則の温度広がり。300 K値は未測定。
- **White et al. (2021)** — 典型的エミッタは5 Kでも約1 GHz。4–40 Kのphonon dephasingとspectral diffusionを定量化。
- **Jungwirth et al. (2016)** — 二種類のZPLの温度shift・linewidth・intensityを測定し、面内phononとのpiezoelectric couplingで記述。
- **Bhat et al. (2026)** — Koch型データをjump-diffusion模型で再現する理論。約25.9 Kのoverdamped crossoverはモデル予測。

## B. 室温spin coherence・spin Hamiltonian

- **Stern et al. (2024)** — 炭素関連単一欠陥、`S=1`、`D/h=1.959 GHz`、`T2*=106(12) ns`、DDで約1.08 µs。ground-coherence budgetの主要入力。
- **Stern et al. (2022)** — 室温単一欠陥ODMR。spin-dependent optical dynamicsはあるが、resonant optical Λと均一幅は未確定。
- **Mathur et al. (2022)** — `V_B^-`の基底ZFS約3.5 GHz、励起ZFS約2.1 GHz。光学EITに必要な狭いresonant transitionとdipole tensorは未取得。
- **Whitefield et al. (2025/2026)** — 室温でoptical spin readoutを示すnarrowband emitter complexes。PL幅はhomogeneous PLE幅ではない。

## C. 高コヒーレンス光学エミッタだがΛ未確立

- **Dietrich et al. (2018)** — 低温で約50 MHzのFourier-limit線。
- **Horder et al. (2025)** — B-centerの材料処理依存coherence。pristine/e-beam生成で最小幅、annealing/C dopingでcharge-noise decoherence。
- **Gérard et al. (2026)** — B-centerのresonance fluorescence、Mollow triplet、indistinguishable photon。現状はcoherent two-level systemでありEIT Λではない。

## D. 欠陥同定・理論入力候補

- **Cholsuk et al. (2024)** — 数百のsinglet/triplet configurationについてZPL、lifetime、quantum efficiency、transition dipole、spin multiplicityを収録。実験欠陥の同定とΛ候補screeningに有用。ただしDFT候補を実験同定なしにno-go入力へ直結させてはならない。
