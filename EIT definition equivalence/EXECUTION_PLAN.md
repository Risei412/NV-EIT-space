# 実行プラン — EIT 定義変更プログラム(同値条件の定式化)

**Phase E: Equivalence-conditions campaign**

このディレクトリは、EIT を

> **「観測応答の中に制御場が生成するコヒーレントな相殺現象」**

として再定義し、その定義と**同値になる条件群**を定理として定式化・検証する
プログラムの作業場所である。戦略文書は
`prl_eit_equivalence_conditions.md`(以下「戦略MD」)、
数学的基盤は Class I–III master-response no-go 定理パッケージ
(`../New no-go theory/Theorem and proofs/three_theorems_proofs.tex` および
アップロード版 PDF)、数値的出発点は `2g2e_package/`
(dark-state-free coherent transparency、79.03% 吸収抑制、2g+2e no-go)である。

---

## 0. 何を証明すれば「定義の変更」が成立するか

新定義の骨子(戦略MD Part I–II):

```
完全 EIT ⟺ χ_full(ω0) = 0 (正則点での transfer zero)
        ⟺ R_S(ω0) = −χ_cut^(S)(ω0), R_S(ω0) ≠ 0 (セクター応答による厳密相殺)
```

ここで `R_S = χ_full − χ_cut^(S)` は master sector response(no-go 理論と同一の量)。

定義変更が PRL レベルで成立するための論理構造は 4 本柱:

| 柱 | 内容 | 状態 |
|---|---|---|
| (a) 一般同値定理 | 応答ゼロ ⟺ セクター相殺(Theorem 1, 2) | 骨子あり。厳密化が必要(T1, T2) |
| (b) 教科書との整合 | 理想 Λ 極限で dark state ⟺ CPT ⟺ 応答ゼロ(Corollary 1) | 主張のみ。証明が必要(T3) |
| (c) 一般には非同値 | 4 つの非含意それぞれに反例(Theorem III) | 主張のみ。反例構成が必要(T4) |
| (d) 構成的存在 | dark state なしで **厳密な** 実軸応答ゼロを実現するモデル | **未達(最大リスク)**。2g+2e では不可能が証明済 → 2g+3e で構成(T6, E4) |

(d) が成功して初めて「dark state は EIT の定義ではなく一機構にすぎない」
という主張が完結する(戦略MD §18「This constructive result would complete
the PRL-level argument」)。(d) の成否がこのプログラム全体のクリティカルパス。

---

## 1. 入力資産

| 資産 | 場所 | 使い方 |
|---|---|---|
| 戦略MD(定理列・主張・証拠要件) | `prl_eit_equivalence_conditions.md` | 定理番号・主張文の原案 |
| Class I–III closed theorem package | `../New no-go theory/Theorem and proofs/` | R_S の定義・不変性(Thm 3.1)・トリコトミー(Thm 7.2)・有限判定(Thm 8.1)をそのまま引用 |
| 2g+2e 数値パッケージ | `2g2e_package/` | (c) の反例 C4(近似版)、(d) の no-go 境界、E3 の再現元 |
| 数値基盤 | `../New no-go theory/src/core.py` | 伝達関数・Krylov 証明書・moment 法・Riesz 射影・ν フィット(再利用) |
| EIT/ATS 分離指標 | Phase A Gate 3(`run_phase_a.py`)の W_S, η_S | E5 の判定器 |

---

## 2. 理論タスク(T1–T7)

### T1. セクターカットの公理化と well-definedness 【基礎・最初にやる】

Theorem 1, 2 の証明自体は分解の線形性から自明なので、この理論の非自明性は
**「χ_cut^(S) が物理的に一意に定義できること」** に集中する。査読で最初に
突かれるのはここ。

定式化: 線形応答系 `A(z) x = c`, `χ = p† x` に対しセクター射影 Π_S を導入し

```
B_cut = B − Q B Π_S − Π_S B Q   (Q = I − Π_S)
```

すなわち「S への出入りの結合のみ除去、S 内部および Q 内部のブロック
(直接光学ブロック・対角散逸・源 c・検出 p・零次定常状態)は保持」を
**セクターカットの公理**として明文化する。

証明すべき補題:
- **L1.1** カットの基底不変性・実現不変性(no-go パッケージ Thm 3.1 (i)(ii) を継承)。
- **L1.2** カット表現不変性の限界の明示: Thm 3.1 (iii) は「同じ cut 伝達関数を
  与える構成は同値」までしか言わない。物理的に異なるカット(異なる S の選択)は
  異なる R_S を与える。→ 定義は常に「対 (χ, S)」に付随することを明記。
- **L1.3** 吸収汎関数 A の規約非依存性: 受動性(散逸性、no-gain)の仮定から
  A[χ_cut^(S)] ≥ 0 を導出し、規約(Im χ / Re χ)の選択が結果を変えないこと。

### T2. Theorem 1(吸収相殺同値)と Theorem 2(複素 transfer-zero 同値)の完成

戦略MD §6, §7 の証明骨子は完成している。追加すべき厳密化:

- **正則性条件**: ω0 が Q_full の零点でない(pole でない)ことの検証可能な形。
- **control-induced 条件**: `N_full(ω0; Ω_c = 0) ≠ 0` かつ
  `N_full(ω0; Ω_c) = 0` — 零点が制御場で移動してきたことの定式化
  (pre-existing zero の除外、戦略MD §7 定義の条件 4)。
- **部分 EIT の定量化**: コントラスト
  `C_S(ω0) = 1 − A[χ_full]/A[χ_cut^(S)] ∈ (0, 1]`
  を導入し、完全 EIT = C_S = 1、2g+2e 数値例 = C_S ≈ 0.79 と統一的に書く。

### T3. Corollary 1(教科書 Λ 同値連鎖)の完全証明

理想 Λ 系で (i) Hamiltonian dark state ⟺ (ii) 定常純 Lindblad dark state
⟺ (iii) CPT ⟺ (iv) 基底コヒーレンス ≠ 0 ⟺ (v) セクター相殺 ⟺ (vi) 応答ゼロ。

証明の落とし穴(必ず処理する):
- **極限の順序**: 弱プローブ極限 Ω_p → 0 と t → ∞ の順序。γ_g = 0 では
  定常状態が初期条件依存になり得るため、「Ω_p ≠ 0 で optical pumping により
  一意化 → 弱プローブ極限」の順で取る(2g+2e レポート §2 の
  order-of-limits caveat と同じ構造)。
- (vi) の ω0 は二光子共鳴点に固定し、γ_g = 0 でのみ厳密ゼロになることを明示。
- 各含意を独立に証明せず、(i)→(ii)→(iii)→(iv)→(v)→(vi)→(i) の巡回で証明。

### T4. Theorem III(一般非同値)— 4 つの反例の最小モデル構成

各非含意に対し、最小次元の具体的 GKSL モデル + 数値検証(E2)を与える:

| 反例 | 非含意 | 最小モデル案 |
|---|---|---|
| C1 | optical dark vector ⇏ 定常純 dark state | Λ + 基底ディフェージング γ_g > 0(ker Ω ≠ 0 だが定常状態は混合) |
| C2 | 定常 dark state ⇏ 観測可能 EIT | 4 準位: dark subspace は存在するがプローブ源との overlap がゼロ(偏光直交)/検出器不感 |
| C3 | CPT ⇏ 完全応答ゼロ | Λ + 第 2 励起状態への off-resonant 結合(ρ_ee ≈ 0 でも別チャネルの吸収で χ_full(ω0) ≠ 0) |
| C4 | 観測可能なコヒーレント透明化 ⇏ dark state | `2g2e_package` の 79% 抑制例(近似版)+ T6 の厳密版 |

### T5. Theorem IV(2g+2e no-go)の境界の明確化

行列式恒等式 `S_p S_c − K_pc K_cp = |det D|² det G` による no-go
(`2g2e_package/docs/` §4)は **g = 0(理想二光子共鳴・γ_g = 0)** で証明
されている。詰めるべきこと:

- 有限 γ_g では Schur 条件(戦略MD Theorem 3 の kernel identity)
  `S_1(γ_g + β S_2) − β K_12 K_21 = 0` は full-rank でも満たし得るか?
  満たせる場合、それは「散逸に助けられた零点」であり、位相コヒーレント相殺
  という定義とどう整合するかを解釈として整理(おそらく γ_g > 0 の零点は
  Condition 6〔位相コヒーレンス〕テストで区別できる)。
- no-go の結論を「γ_g = 0 の理想境界における障害」と正確にスコープする。

### T6. 構成的存在定理 — 2g+3e での厳密 dark-state-free EIT 【最重要・最大リスク】

**数学的機構(このプランの核心)**: 励起多様体を 3 次元にすると D ∈ C^{3×2},
G(δ) ∈ C^{3×3} となり、Cauchy–Binet 恒等式により

```
S_p S_c − K_pc K_cp = det(D† G D)
                    = Σ_{|I|=|J|=2} det((D†)_I) · det(G_{IJ}) · det(D_J)
```

これは G の 2×2 小行列式の**複数項の和**になり、2e の場合の
`|det D|² det G`(単項、rank D = 2 なら零になれない)と違って、
**rank D = 2 のまま項間の複素位相相殺で和がゼロになれる**。
つまり dark state(det D = 0)なしの厳密実軸零点の障害が消える。

構成手順:
1. G を対角(3 つの励起共鳴 a_1, a_2, a_3)+ 必要なら閉ループ位相で
   パラメトライズし、零点条件を Re/Im の実 2 方程式とみなす。
2. 実パラメータ 2 つ(例: 第 3 双極子角 θ₂ と離調 Δ_e2)で陰関数定理により
   解の存在を示す(codim 2 → ジェネリックに可解)。まず sympy で厳密解を探す。
3. 解の点で 7 条件(§3 Gate E4 のチェックリスト)を全て検証。
4. 代替ルート(2g+3e が詰まった場合): 閉ループ結合(loop 位相 Φ を
   相殺パラメータに使う)/相関散逸(correlated jumps)/NV 6 準位励起多様体。

**Fallback(重要)**: もし 2g+3e でも厳密零点が不可能なら、それは
Cauchy–Binet 和の符号構造による**拡張 no-go 定理**であり、それ自体が
主結果になり得る(「dark-state-free 厳密 EIT の最小次元はさらに上」)。
どちらに転んでも論文は成立するように、E4 は「存在の構成 or 障害の証明」
の二値ゲートとして設計する。

### T7. Class I–III との統合定理(EIT verdict)

戦略MD Part VII の

```
class(R_S) + (χ_cut に対する R_S の相対位相・振幅) = 観測 EIT の判定
```

を定理化する。具体的には:
- Class III でも位相が合わなければ EIT にならない例、Class II でも有限 Γ で
  EIT になる例を構成(E5)。
- ATS 除外(Condition 4)の定量化: 極分離 |Re z₁ − Re z₂| と半値幅
  min(|Im z₁|, |Im z₂|) の比、および留数の符号反転
  (`2g2e_package` §6 の判定と Phase A Gate 3 の W_S, η_S を統一)。

---

## 3. 数値検証キャンペーン(Phase E ゲート)

実装は `../New no-go theory/src/core.py` の基盤(伝達関数、Krylov 証明書、
moment 法、Riesz 射影、ν フィット)を再利用。新規コードは本ディレクトリの
`src/` に置く:

```
EIT definition equivalence/
├── EXECUTION_PLAN.md            (本書)
├── prl_eit_equivalence_conditions.md
├── 2g2e_package/                (取り込み済み数値パッケージ)
├── src/
│   ├── model_lambda_chain.py    (E1: Λ 同値連鎖)
│   ├── model_counterexamples.py (E2: C1–C4)
│   ├── model_2g2e_boundary.py   (E3: no-go 境界)
│   ├── model_2g3e_zero.py       (E4: 構成的零点探索; sympy + numpy)
│   ├── model_class_verdict.py   (E5: クラス×位相 → verdict)
│   └── run_phase_e.py           (全ゲート実行 + gates_summary_phaseE.json)
├── results/                     (figures, JSON サマリ)
└── tex/equivalence_theorems.tex (定理パッケージ; three_theorems_proofs.tex と対)
```

| ゲート | 内容 | PASS 基準 |
|---|---|---|
| **E0** | 資産取り込み(2g2e_package)と `dark_state_free_2g2e_analysis.py` の再現実行 | 公表値(δ_min = −0.01, 抑制 79.027%, z₀ = −0.00984 − 0.02120i)を再現 |
| **E1** | Λ 同値連鎖 (i)–(vi) の機械検証: γ_g = 0 で全鎖一致、γ_g > 0 で全鎖が同時に破れ始めること | 各同値量が machine precision で一致 / γ_g スキャンで破れの連続発生 |
| **E2** | 反例 C1–C4 の witness 計算(ker Ω, 定常状態純度, ρ_ee, χ_full(ω0), R_S) | 4 つの非含意すべてで witness が定義基準を満たす |
| **E3** | 2g+2e no-go 境界: g = 0 で行列式恒等式と \|Im z₀\| > 0 の整合、γ_g スキャンで実軸零点の可否 | 恒等式残差 < 1e−12; γ_g > 0 零点の存在有無を確定しレポート |
| **E4** | **2g+3e 厳密零点の構成(or 拡張 no-go の証明)** | 下記 7 条件チェックリスト、または不可能性の記号的証明 |
| **E5** | クラス×位相 verdict テーブル(Class III で EIT なし / Class II で有限 Γ EIT) | verdict が定理 T7 の予言と全点一致 |
| **E6** (stretch) | NV 6 準位励起多様体への写像、R_S ∝ B⊥² Γ_oc⁻² スケーリング(既存 `SIMULATION_PLAN.md` Gate 資産再利用) | Class II スケーリング指数一致 |

**E4 の 7 条件チェックリスト**(戦略MD §15 の操作的条件に対応):

1. `|χ_full(ω0)| < 1e−12` を実軸上の ω0 で(sympy による厳密ゼロなら最良)
2. `σ_min(Ω) > 0`(full rank、optical dark vector なし)
3. 定常純 Lindblad dark state の不存在(jump 共通固有ベクトル判定)
4. `R_S(ω0) ≠ 0` かつ `R_S(ω0) = −χ_cut^(S)(ω0)`
5. セクターカットで透明化が消える(`A[χ_cut](ω0) > 0`)
6. 極・留数分析で ATS 除外(極分離 < 半値幅、留数符号反転)
7. 受動性(観測域で gain なし)+ control-induced(Ω_c → 0 で零点消失)

---

## 4. 執筆物

1. `tex/equivalence_theorems.tex` — 定理パッケージ本体。構成は戦略MD Part VIII:
   - Theorem I(一般応答 EIT 同値)/ Theorem II(複素 transfer-zero 同値)
   - Corollary(教科書 Λ 同値)/ Theorem III(一般非同値 + 反例 4 種)
   - Theorem IV(2g+2e 障害)/ **構成的定理(2g+3e 厳密零点)or 拡張 no-go**
   - Class I–III 統合定理(verdict)
2. PRL 本文への反映は `../writing paper/The theory of EIT in Nitrogen Vacancy.tex`
   系列とは独立の新規原稿として開始(定義変更論文は NV 特化論文と分離)。
3. Definition of Done は戦略MD §23 の 10 項目チェックリストをそのまま採用。

---

## 5. 実行順序・依存関係・リスク

```
T1 ──→ T2 ──→ T3 ──→ E1
 │              │
 ├──→ T4 ──→ E2 (C1–C3), E0/E3 (C4 近似版)
 ├──→ T5 ──→ E3
 └──→ T6 ══→ E4  ←― クリティカルパス(存在 or 拡張 no-go の二値ゲート)
        │
T7 ────┴──→ E5 ──→ tex 執筆 ──→ (stretch) E6
```

推奨実行順: **E0 → T1 → T6/E4(リスク最大なので早期着手)→ T2–T5/E1–E3(並行)→ T7/E5 → 執筆**。

主要リスクと対処:

| リスク | 対処 |
|---|---|
| E4 で 2g+3e 厳密零点が存在しない | Cauchy–Binet 構造から拡張 no-go を証明し主結果を反転(§2 T6 Fallback)。閉ループ・相関散逸ルートも並行検討 |
| セクターカットの一意性への査読反論 | T1 の公理化 + L1.2 で「定義は (χ, S) 対に付随」と先回りして明示 |
| γ_g > 0 零点の解釈(散逸支援 vs コヒーレント) | Condition 6(位相コヒーレンステスト)を判定器として定式化(T5) |
| EIT/ATS 分離の恣意性批判 | Phase A Gate 3 の定量指標(W_S, η_S)と極・留数判定を統合し閾値を明記(T7) |

---

## 6. 再現・実行手順(実装後)

```bash
cd "EIT definition equivalence"
pip install numpy scipy sympy matplotlib
python src/run_phase_e.py            # 全ゲート実行
# 出力: results/gates_summary_phaseE.json, results/figures/figE*.png
```

各ゲートは `../New no-go theory/results/gates_summary_phase*.json` と同形式の
PASS/FAIL サマリを出す。全ゲート PASS(E4 は「構成」または「拡張 no-go」の
いずれかで PASS)が執筆開始の条件。
