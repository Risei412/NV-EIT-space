# Phase F 先行監査結果(F0/F1/F2/F4 のスポットチェック)

`NEXT_NUMERICAL_PLAN.md` の受領直後に、本格実装の前へ分岐点を確定させるため
4 つのスポットチェックを実行した。**結論: 監査指摘は全面的に正しく、
命題 A・B はいずれも成立。さらに命題 A は基底コヒーレンス(γ_g > 0)を含む
全コヒーレンス空間へそのまま拡張され、F4 の結末も理論的に確定した。**

---

## 1. F0/F3 先行: 旧 238 解はすべて非受動(命題 B 成立)

旧 closed-loop モデル(実数非対角 `A_23 = J_23`)の全 238 解について
`m_Γ = λ_min((A+A†)/2)` を数値計算した結果:

```
min |J23| over all 238 solutions : 1.1214   (> γ_e = 1)
max m_Γ  over all 238 solutions  : -0.1214  (< 0)
count with m_Γ ≥ 0               : 0
```

**全解が受動性を破って初めて存在する零点**であり、旧解は
「passivity-breaking boundary を越えると exact zero が生成される」対照例
(計画書 Outcome IV)として再分類する。`E4_findings.md` の「構成成功」は
この意味に格下げして訂正する。

## 2. F2 先行: 命題 A(strict-accretive no-go)の数値証明書

恒等式 `(G+G†)/2 = G†ΓG` と `Re M = D†G†ΓGD > 0` を、ランダムな
Hermitian 混成・非可換 Γ・non-normal A を含む 20,000 試行
(N_e = 2〜6)で検証:

```
恒等式残差の最大値            : 1.7e-14  (機械精度)
λ_min(Re M) の全試行最小値    : 3.5e-05  (> 0、反例なし)
```

対角性は本質ではなく、**strict passivity(Γ > 0)+ matched response が障害**
である、という計画書の読み替えが正しい。

## 3. F1 先行: 物理的 Hermitian 混成(A_23 = iJ23)では零点が消える

旧探索と同一の機構で `iJ` 実装の零点を 600 始動で探索した結果、
solver は 46 件の「成功」を返したが、**全件が逃走解**
(|J23| または |δ| → ∞ で det M がアンダーフロー、σ_min(M) ~ 1e-76)であり、
**有限パラメータの真の零点は 0 件**。計画書 §6.3 の警告
(「root solver の成功表示を信用しない」)がそのまま現実になった。
命題 A の帰結と完全に整合。

## 4. F4 先行: 命題 A は γ_g > 0 の全コヒーレンス空間へ拡張される(重要)

弱プローブ線形化の変数を
(光学コヒーレンス σ_{e_j g1}, j=1..N_e; 基底コヒーレンス σ_{g2 g1})
に取ると、係数行列は

```
A_full = Γ_full + i K_full,
Γ_full = diag(γ_1,...,γ_Ne, γ_g),
K_full = Hermitian(離調、Hermitian 励起混成 J、制御結合 (Ω_c/2)d_c)
```

の形になり、プローブの source と readout は同一ベクトル
`c = p = (d_p, 0)`(**matched**)。検証:

- この matched 全空間形式は 2g+2e パッケージの閉形式 χ_full を
  機械精度(2.8e-16)で再現する。すなわちこれは正しい物理表現である。
- したがって命題 A がそのまま適用され、γ_g > 0 なら Γ_full > 0 で

  ```
  Re χ_full = c† G† Γ_full G c > 0   (厳密に正)
  ```

- 30,000 のランダム matched モデル(N_e=1..4、γ_g = 1e-6..1、
  ランダム Hermitian K、実軸離調)で `min Re χ_full = 1.7e-05 > 0`、反例なし。

### 帰結(F4 は理論決着)

1. **matched passive Markovian 応答では、有限 γ_g の exact zero どころか
   absorption-only zero(Re χ = 0)すら不可能。** 透明化はすべて有限コントラスト
   dip(C_S < 1)である。F4 の分類作業は「exact / absorption-only / dip」の
   探索ではなく、**コントラスト上限の定量化**に置き換わる。
2. **教科書の完全 EIT(γ_g = 0)は、まさに passivity boundary
   λ_min(Γ_full) = 0 上の零点**である。dark state による完全 EIT とは
   「基底コヒーレンスという lossless mode を持つ系の境界現象」と
   統一的に再解釈できる。
3. 新しい定理候補(**透明化フロア定理**):passivity margin から

   ```
   Re χ_full(δ) ≥ λ_min(Γ_full) · |G(δ)c|²  >  0
   ```

   の形の下限が従い、1 − C_S の定量的下限(透明化の到達不可能領域)を
   γ_g で表せる。これは PRL の主張を「no-go + 定量境界」へ強化する。

## 5. Phase F 優先順位への影響

- F0–F4 は「探索」から「証明書化+図表化」へ性格が変わった(結末確定)。
- 真の研究フロンティアは **F5-A(unmatched readout / 直接経路)と
  F5-B(完全 GKSL)** である。決定的な観点:完全 GKSL の線形応答は
  `c = L_probe ρ_0`(制御でドレスされた定常状態への作用)、
  `p = 観測双極子` であり、**一般に p ≠ c(unmatched)**。
  命題 A は適用されず、受動系のまま零点が生じる余地はここにしかない。
- したがって F5-B の中心的問いは「制御場による ρ_0 のドレスが応答を
  どの程度 unmatch し、それが gain なしで numerator zero を生めるか」。

再現: 本ノートの数値は使い捨てスクリプトによる先行チェックであり、
正式な証明書は `audit_closedloop_convention.py` / `model_passive_general.py`
(Gate F0–F2)として実装し直す。
