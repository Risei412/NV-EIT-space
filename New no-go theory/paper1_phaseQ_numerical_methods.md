# Paper I / Phase Q 数値計算の実装方針
## Scaling-Path-Dependent Class II/III Transition — 物理 CPTP Lindblad モデル

対象: 二エミッタ bright/dark モデル(3準位 Hilbert 空間、9×9 Liouvillian)で
Path T (ν=0, Class III) / Path U (ν=2, Class II) の経路依存クラス転移を実証する計算
(計画書 `paper1_scaling_path_numerical_plan.md` 相当、戦略文書の Paper I)。

位置づけ:

- 抽象 4×4 行列モデルでは Phase D3 (`src/model_paths.py`, `run_phase_d.py`) で
  ν_T = −0.0022, ν_U = 1.961、class fan まで **確認済み**
  (`results/gates_summary_phaseD.json`)。
- 今回の新規性は「完全に CPTP な物理 Lindblad モデル」で同じ転移を示すこと。
  理論側は theorem package (Theorems I–III, 有限判定定理) で閉じており、
  本計算は Theorem II (invertible D̂ → ν = m+1) と
  Theorem III (singular semisimple D̂ → Riesz 射影, ν = 0) の物理実証にあたる。

---

## 1. 最重要の数値的問題: full − cut 減算の桁落ち

計画書は R_S = ⟨D|(ρ_ss^full − ρ_ss^cut)|0⟩ を Γ/Γ₀ ∈ [1, 10⁸] で要求するが、
**full と cut の定常状態を別々に解いて引く素朴な実装は Γ ≳ 10⁵–10⁶ で破綻する**。

オーダー見積り (ε = 10⁻³, J = 0.2):

| Γ | \|R_U\| ≈ εJ/Γ² | 状況 |
|---|---|---|
| 10³ | 2×10⁻¹⁰ | 可 |
| 10⁵ | 2×10⁻¹⁴ | 危険域 (κ(𝕃) 増大でさらに悪化) |
| 10⁶ | 2×10⁻¹⁶ | 機械精度の床 (ρ₀₀ ≈ 1 の状態を引くため誤差床 ~ε_mach·‖ρ‖ ≈ 10⁻¹⁶) |
| 10⁸ | 2×10⁻²⁰ | 素朴な減算では 4 桁ノイズに沈む |

漸近 fitting window [10⁴,10⁶], [10⁵,10⁷], [10⁶,10⁸] はすべてこの危険域にあり、
対策なしでは Gate Q4 (ν_U = 2) が **偽の失敗**を起こす。

### 対策: 差分を直接解く (difference solve)

δρ = ρ_full − ρ_cut は厳密に

```
𝕃_full |δρ⟩⟩ = −Δ𝕃 |ρ_cut⟩⟩,   Δ𝕃 = 𝕃_full − 𝕃_cut = −i[J|D⟩⟨D|, ·],
Tr δρ = 0
```

を満たす。cut の定常状態を解いた後、この線形方程式で δρ を直接求めれば
O(1) 同士の減算が発生せず、相対精度が保たれる。R_S = ⟨D|δρ|0⟩ を直接読む。

2×2 プロトタイプでの検証結果 (このドキュメント作成時に実測):
差分ソルバは Γ = 10⁸ でも Γ²|R_U|/ε = 0.200000 (= J) を machine precision で維持。
素朴な減算は Γ ~ 10⁶ で |R_U| ~ 2×10⁻¹⁶ となり破綻する。

---

## 2. ソルバの 3 層構成

### Level 1 (主計算): O(ε) 線形応答・コヒーレンス小ブロック

probe ε の 0 次で定常状態は ρ⁰ = |0⟩⟨0| に厳密に落ちる。O(ε) では
コヒーレンス (ρ_B0, ρ_D0) だけが立ち、9×9 の vectorized Liouvillian は
2×2 複素ブロックに厳密に閉じる:

```
M(J) v = −i e_D,   v = (ρ_B0, ρ_D0)/ε,
M(J) = [ γ_B + iΔ_B      iδ/2            ]
       [ iδ/2            γ_D + i(Δ_D+J)  ]
```

(δ は §3 のエミッタ非対称項。計画書の最小モデルは δ = 0。)

- ε がスケールアウトするので weak-probe 線形性は構成上厳密。
  R_S/ε を直接得る。
- 差分ソルバは `M(J) d = −iJ (0, v_cut,D)` の 2×2 解。ここまで小さいと
  条件数の問題は消滅し、Γ = 10⁸ どころか任意の Γ で安全。
- sweep (Γ 300–500 点 × θ 20 点 × robustness) が全て秒オーダーで回る。

### Level 2 (検証): full 9×9 非線形 Lindblad 定常解

計画書 §7.1 の通り column-stacking で 𝕃 を組み、trace 行置換
(冗長行を vec(I)† に差し替え、RHS = e_trace) で ρ_ss を解く。

- probe ε を H に入れた**非線形 (全次数) 解**であり、Level 1 との比較が
  Gate Q5 (reduced vs full) の実体になる。一致は O(ε²) ~ 10⁻⁶ 相対まで
  期待でき、計画書の目標 10⁻³ を大きく上回る。
- full−cut は必ず §1 の差分ソルバ (𝕃_full δρ = −Δ𝕃 ρ_cut, Tr δρ = 0 を
  同じ trace 行置換で課す) を使う。素朴な減算は Γ ≤ 10⁴ の
  クロスチェック専用とする。
- 行スケーリング (equilibration) を入れ、residual ‖𝕃|ρ_ss⟩⟩‖₂ を毎点記録
  (Gate Q2)。
- 計算 C (ε ∈ {10⁻²,…,10⁻⁵} の線形性チェック) はこの Level で実施し、
  「Level 1 が厳密な ε→0 極限である」ことの数値確認と読み替える。

### Level 3 (確定): 厳密証明書 — 有限窓 fitting に依存しない ν

theorem package の有限判定定理 (Theorem 8.1) をそのまま実装する:

1. **多項式証明書**: R_S(Γ) = N(Γ)/Q(Γ) は Γ の有理関数。sympy で
   2×2 (または 9×9) ブロックの adj/det を Γ シンボリックに計算し、
   `ν = deg_Γ Q − deg_Γ N` を**厳密に**得る。
   Path U: deg Q = 4, deg N = 2 → ν = 2。Path T: deg Q = 2, deg N = 2 → ν = 0。
   これで preasymptotic masquerade (Paper III の主張 C) の懸念を原理的に排除。
2. **moment 法** (Theorem II, `core.py` 実装済み): Path U では D̂ が可逆なので
   m₀ ≡ 0, m₁ ≠ 0 を数値確認 → ν = m+1 = 2 が整数として出る。
3. **Riesz 射影** (Theorem III, `core.py` 実装済み): Path T では scaled dissipator
   が dark コヒーレンス方向に semisimple zero を持つ。R_{S,0} を式 (20) の
   protected-block 逆行列で評価し、Level 1/2 の plateau と照合 (相対誤差 < 10⁻⁶
   を目標; Phase D では ~10⁻⁷ を達成済み)。

fitting (局所有効指数 + 複数 window 回帰) は Level 1/2 の結果に対する
「見え方」の検証であり、クラス判定の根拠は Level 3 に置く。
論文の主張構造 (有限窓では class が偽装され得る → 有限判定が必要) とも整合する。

---

## 3. モデル上の注意: B–D 非結合の自明性と δ 項の追加

計画書の H は probe が |D⟩ のみを駆動し、B と D を結ぶ項がない。このままだと
bright セクターは応答から**完全に切り離され**、R_S は γ_B に厳密に非依存
— ν_T = 0 は成立するが「Path T は何もしていない」という自明性批判を受ける。

対策: エミッタ周波数非対称 δ (物理的には 2 エミッタの detuning 差 ±δ/2) を加える:

```
H → H + (δ/2)(|B⟩⟨D| + |D⟩⟨B|)
```

- δ ≠ 0 なら dark コヒーレンスは強く減衰する bright モードと結合しつつ、
  γ_B → ∞ で Schur 補正 (δ/2)²/(γ_B + iΔ_B) → 0 により**非自明に**保護される
  (Zeno 型 decoupling; Theorem III の Riesz/Schur 機構そのもの)。
- 既定値 δ = 0.1 Γ₀。δ = 0 は解析式との厳密照合用リミットとして残す。
- δ は H の項なので両 path 共通であり、同一物理点条件 (Gate Q1) は不変。
- 2×2 プロトタイプ実測: δ = 0.1 で |R_T|/ε → 0.1713412 の非自明 plateau
  (γ_B 依存の過渡を経て収束)、Path U は Γ²|R_U|/ε → J を維持。

robustness (§12) に δ sweep ([0, 0.5]Γ₀) を追加し、δ = 0 (自明保護) と
δ > 0 (Schur 保護) で ν_T = 0 が共通であることを示す。

---

## 4. 実装構成 — 既存インフラへの追加

新規ディレクトリではなく、既存の Phase A/B/D/M/P と同じ規約で
`New no-go theory/` に追加する (計画書 §14 の standalone 構成は採らない):

```
New no-go theory/
├── src/
│   ├── model_bright_dark.py   # H(δ, J), L_B, L_D, L_φ; 9×9 𝕃 と 2×2 線形応答ブロック
│   │                          #   γ_B(Γ,θ), γ_D(Γ,θ) の path 定義 (T: θ=0, U: θ=π/2)
│   ├── run_phase_q.py         # Gates Q1–Q7 実行、figQ*.png、gates_summary_phaseQ.json
│   └── (core.py, report.py)   # 再利用: moment 法 / Riesz / 局所 fit / gate 集計
├── results/
│   ├── gates_summary_phaseQ.json
│   ├── phaseQ_baseline.csv          # Γ, |R_T|, |R_U|, ν_eff (Level 1/2 両方)
│   ├── phaseQ_class_fan.npz         # (Γ, θ) → ν_eff
│   └── figures/figQ1..figQ6_*.png   # 計画書 Figure 1–6 対応
└── tests/ 相当は run_phase_q.py 内の gate 関数として実装 (既存 phase の流儀)
```

パラメータ既定値 (無次元, 計画書 §7.2 + δ):
Γ₀ = 1, J = 0.2, Δ_D = 0.3, Δ_B = 0, δ = 0.1, ε = 10⁻³ (Level 2 のみ),
Γ/Γ₀ ∈ [10⁰, 10⁸] 対数 400 点。

---

## 5. グリッドと fitting の具体値

- **ν_eff**: 対数グリッド上の中心差分 (`core.py` の local fit)。端点除外。
- **漸近 window**: [10⁴,10⁶], [10⁵,10⁷], [10⁶,10⁸] の 3 window 線形回帰、
  max−min < 0.03 (Gate Q3/Q4)。Level 1 では全域で有効。
- **class fan (θ)**: γ_D(Γ,θ) = Γ₀ + sin²θ (Γ−Γ₀)。θ = 0 と、
  sin²θ を [10⁻⁶, 1] で対数 20 点。
  crossover Γ_× は ν_eff(Γ) = 1 の交差で定義。
  Γ_× ~ Γ₀/sin²θ ≤ 10⁶ となるのは sin²θ ≥ 10⁻⁶ の範囲で、
  fitting 可能な Γ_× ∈ [10¹, 10⁶] に収まるよう **slope fit には
  sin²θ ∈ [10⁻⁵, 10⁻¹] の点のみ**使う (それ以外は fan 図の表示用)。
  合格: |s_fit + 1| < 0.05 (Gate 相当は Phase D の fan 判定を流用)。
- **collapse**: x = Γ sin²θ/Γ₀ に対する ν_eff の master curve。
  Phase D3/P4 の collapse 判定 (相対 spread) を流用。

---

## 6. Gates (計画書 Q1–Q5 + 追加 2 件)

| Gate | 内容 | 実装 |
|---|---|---|
| Q1 | 同一物理点: ‖𝕃_T(Γ₀)−𝕃_U(Γ₀)‖_F / ‖𝕃_T‖_F < 10⁻¹², ρ_ss/R の一致 < 10⁻¹⁰ | 構成上厳密のはずだが必ず実測 |
| Q2 | CPTP 監査: Tr, Hermiticity, λ_min > −10⁻¹⁰, residual < 10⁻¹⁰ | Level 2 の全 Γ 点で記録 |
| Q3 | Path T: \|ν_T\| < 0.03, plateau 非零 | Level 1/2 + Riesz R_{S,0} 照合 |
| Q4 | Path U: \|ν_U−2\| < 0.03, std/mean(Γ²\|R_U\|) < 0.02 (最終 decade) | Level 1/2 + moment 法 (m₀≡0, m₁≠0) |
| Q5 | reduced (Level 1) vs full (Level 2): 相対差 < 10⁻³ (目標 O(ε²)~10⁻⁶), \|Δν\| < 0.02 | 同一 Γ グリッドで直接比較 |
| **Q6 (追加)** | 精度監査: 差分ソルバ vs 素朴減算 (Γ ≤ 10⁴ で一致確認)、mpmath 50 桁での spot check (Γ ∈ {10⁴,10⁶,10⁸} 各 path) | 桁落ち対策の妥当性を明示的に閉じる |
| **Q7 (追加)** | 厳密証明書: sympy で deg Q − deg N (T: 0, U: 2)、Path U の moment 整数判定 | 有限窓 fitting 非依存のクラス確定 |

---

## 7. 実行順序

1. **Q1 段階**: `model_bright_dark.py` (δ = 0) → Level 1 で Path T/U sweep →
   ν_T = 0, ν_U = 2 を確認。解析式 R_S ∝ −iεJ/[(γ_D+i(Δ_D+J))(γ_D+iΔ_D)]
   (δ=0 の閉形式) と machine precision で照合。
2. **Q2 段階**: Level 2 (9×9) 実装 → Gates Q1, Q2, Q5, Q6。ε sweep。
3. δ = 0.1 に切替えて 1–2 を再実行 (非自明保護版が論文の主結果)。
4. **Q3 段階**: class fan + crossover collapse (Level 1 で全グリッド、
   Level 2 で代表 θ 3 点を検証) → Figure 4, 5。
5. **Q4 段階**: robustness (Δ_D ∈ [−1,1], J ∈ [10⁻³,1], δ ∈ [0,0.5],
   dephasing γ_φ 固定/スケール) + Gate Q7 → Figure 6。
   - dephasing は L_φ = √γ_φ(|D⟩⟨D|−|0⟩⟨0|) で 0–D コヒーレンス減衰 2γ_φ。
     γ_φ ~ Γ で ν_T: 0 → 2 への移行を確認 (「どのチャネルを scale するかが
     class を決める」の直接実証)。
   - Δ_D sweep では Δ_D = 0 および Δ_D + J = 0 の近傍点を明示的に含め、
     γ_D0 = 1 のため特異でないことを確認しつつ prefactor の谷を記録。
6. 図の組版 (Figure 1–6) と `results/summary.md` への追記。

期待値の再掲 (Level 1 実測済みプロトタイプ、δ = 0.1):
|R_T|/ε → 1.713×10⁻¹ (plateau), Γ²|R_U|/ε → 0.2000 (= J), Γ_× ~ Γ₀/sin²θ。

全計算は 9×9 以下の線形代数のみで、robustness 込みでも数分以内に完了する。
