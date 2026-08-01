# 論文化前監査レポート — 定理証明および計算コード

対象コミット: `20de02a` / 監査日: 2026-08-01
監査環境: Linux 6.18.5, Python 3.11.15, numpy 2.4.6, scipy 1.17.1,
matplotlib 3.11.1, sympy 1.14.0, mpmath 1.3.0, pandas 3.0.5

監査の依頼事項は二点。
1. 各定理の証明に他の方法がないか
2. 計算コードが PC 環境の異なる他者の手元で再現するか

---

# 総括

**証明**: `No-go theorem/Theorem and proofs/eit_nogo_proofs.tex` の定理 1A–7B、
補題、系をすべて検証した。**数学的な誤りは見つからなかった。** ただし
(a) 定理 1B の証明中に一箇所、書き下されていない代数変形があり、
(b) 定理 3・4・6・7A は、現在採用されている証明よりも**仮定の弱い / 主張の強い
別証明**が存在する。特に定理 3 と定理 4 は、別証明を採ると仮定 `q<1` と
`N` が不要になり、結果が真に強くなる。
(c) 定理 1A/1B/2A/2B/5 はいずれも既知結果の別表現であり、査読では
**新規性の切り分け**を明示的に書かないと priority を突かれる。

**再現性**: 全体としては**非常に良好**。RNG は例外なく `default_rng(固定 seed)`、
matplotlib は `Agg` 固定、出力先は `__file__` 相対、`differential_evolution` は
`workers=1` 明示。著者環境とは世代の違うライブラリ (numpy 2.4 / scipy 1.17 /
pandas 3.0) で全 47 テストが通り、Gate A–D・RoomT step1–9・Phase A/B/D/M/P・
No-go theorem gate1–5 の JSON/CSV は `runtime_s` を除き**完全一致**、PNG は
**ピクセル完全一致**した。

ただし公表前に必ず潰すべき欠陥が **4 件** ある。

| # | 深刻度 | 内容 |
|---|---|---|
| A-1 | **High** | `gate_F5B_full_gksl.json` が古い(バグ修正前の)結果のまま。gate 判定 `any_gain_in_matched` が再実行で `true → false` に反転し、同梱の `F5B_findings.md` の記述とも矛盾する |
| A-2 | **High** | 原稿 Fig. 4 の committed PNG/PDF が `quick=True` 版。README が指示する `python make_figures.py`(= `quick=False`)は図中に印字される値まで異なる図を出す |
| A-3 | **Medium** | `2g2e_package` の「再現パッケージ」の生成スクリプトが、同梱の CSV も PNG も生成しない。加えて唯一 `Agg` 未指定で `plt.show()` を呼ぶ |
| A-4 | **Medium** | T1 定理 2.2 の反例の数値が、リポジトリのどのコードからも生成されない |

以下、詳細。

---

# 第 I 部 定理の証明の監査

対象: `No-go theorem/Theorem and proofs/eit_nogo_proofs.tex`(929 行、定理 1A–7B)
および `EIT definition equivalence/tex/T1_sector_cut_axiomatization.tex`。

## I-1. 検証結果(正しさ)

各証明を独立に追跡した。結論として、**すべての定理の主張は正しく、証明も
(下記 1 点を除いて)完結している。**

補助的に確認した具体的整合性:

- 定理 2A の注(教科書 Λ 系への帰着): `Ξ = γ_g/(Aγ_g+β)` の代数を再計算し一致。
- 定理 2B の注(`χ_full = 0` は no-go ではない): `A=1, d₁=d₂=1, β=1, γ_g=0` で
  `Ξ_full = 0`, `δΞ_S = −1` を再計算し一致。
- 定理 4 の `B(Γ) = adj(ΓI−X)` の構成: telescoping と Cayley–Hamilton による
  検証手順は正しい。
- 系(hopping rate map)の `γ_oc = Γ_XY/4`: 対称二方向 hopping で
  `Γ_XY = 2k`、各光学コヒーレンスの減衰が `k/2` になることを再計算し一致。

## I-2. 唯一の証明上のギャップ — 定理 1B

証明 (L. 240–242) の

> A short rearrangement shows this equals `−i(H̃ρ_D − ρ_D H̃†) + c ρ_D`

が書き下されていない。ここは実際に非自明で、`L_μ|D⟩ = λ_μ|D⟩` からは
`ρ_D L_μ† = λ_μ* ρ_D` は従うが `ρ_D L_μ = λ_μ* ρ_D` は従わない
(`L_μ†|D⟩` は `|D⟩` の固有ベクトルとは仮定されていない)ため、
「短い変形」で片づけると査読で必ず止まる。

**推奨する書き換え**(結論は変わらない、経路だけ差し替える):
`H_eff := H − (i/2)Σ_μ L_μ†L_μ` を導入して

```
L(ρ_D) = −i(H_eff ρ_D − ρ_D H_eff†) + Σ_μ L_μ ρ_D L_μ†
```

と非エルミート形に書く。(i) より `Σ_μ L_μ ρ_D L_μ† = (Σ_μ|λ_μ|²) ρ_D`。
`L(ρ_D)=0` の両辺を右から `|D⟩` に作用させると `ρ_D|D⟩ = |D⟩` により

```
H_eff|D⟩ = ( ⟨D|H_eff†|D⟩ − i Σ_μ|λ_μ|² ) |D⟩
```

が直ちに出る。左辺の反エルミート部を分離すれば (ii) が得られる。
`ρ_D L_μ` を経由しないので、上記の落とし穴を回避できる。

**先行研究**: この定理は Kraus 型の標準結果であり、
Albert & Jiang, *Phys. Rev. A* **89**, 022118 (2014)
("Symmetries and conserved quantities in Lindblad master equations") の
dark-state 構造定理の特別な場合である。初等証明を自前で置くのは
appendix として妥当だが、本文では出典を明示すべき。

## I-3. 別証明の可能性(依頼事項 1 への回答)

「他に方法がないか」に対して、定理ごとに実在する代替経路を挙げる。
★ は**採用を推奨する**(仮定が弱くなる/主張が強くなる/文献接続が良くなる)もの。

### 定理 1A — 光学 dark subspace の rank
現行: rank–nullity。

- **SVD 版**: `Ω = UΣV†` の零特異値の重複度が `dim D_opt`。数え上げは同じだが、
  最小非零特異値 `σ_min` が「どれだけ dark に近いか」の**定量指標**を同時に与える。
  実験の不完全性(偏光誤差、strain)に対する頑健性を論じる節があるので、
  こちらの方が後段と接続が良い。
- ★ **Morris–Shore 変換**: 縮退二準位系を独立 Λ 系 + dark 状態群にブロック
  対角化する標準手法。リポジトリの文献ノート `C03`(Morris–Shore の多準位拡張)
  および `C02` が既にこれを扱っている。定理 1A は本質的に MS 分解の dark ブロック
  の次元計算であり、**そう書けば競合文献 C02/C03/C05/C07 との関係が自明になる**。
  現状の書き方だと「独立に再導出した」ように読め、新規性の主張が過剰に見える。

### 定理 1B — 定常純粋 Lindblad dark 状態
現行: purity `Tr ρ²` の時間微分 → 分散和 → 各分散が消える。

- ★ 上記 I-2 の `H_eff` 経路(初等性を保ったまま穴を塞ぐ)。
- **散逸子の正値性による直接証明**: `⟨D|L(ρ_D)|D⟩ = 0` を Spohn/Frigerio の
  H-定理型の議論で処理する。purity ではなく相対エントロピーの単調性を使う版。
- **半群の不動点構造から**: CPTP 半群の不動点集合は commutant 構造を持つ
  (Baumgartner–Narnhofer)。純粋不動点はその極点であり、(i)(ii) が構造定理から
  一行で出る。ただし「学部レベル自己完結」の方針とは相容れない。

### 定理 2A / 2B — Schur 補元としての被約感受率
現行: ブロック代入(ブロック Gauss 消去)。

- **Feshbach–Fano 射影法**: 物理側の標準名称。同じ式だが、`A^{-1}` が
  有効ハミルトニアン/自己エネルギーとして解釈でき、後段 §14.3 の
  「周波数依存自己エネルギー」の議論と用語が揃う。**改名だけで文献接続が改善する。**
- ★ **Sherman–Morrison(rank-1 更新)**: 定理 2A の結合ブロックは
  `B = (iΩ_c/2)d₂`、`C = (iΩ_c*/2)d₂†` で **rank 1**。したがって Schur 補元を
  経由せず Sherman–Morrison 公式一行で
  `δΞ_S = −β K₁₂K₂₁/(γ_g + βS₂)` が出る。
  「`δΞ_S = 0 ⟺ K₁₂K₂₁ = 0`」という判定条件が
  **rank-1 更新の分子がゼロになる条件**として構造的に見えるようになり、
  なぜ積 `K₁₂K₂₁` なのかが自明になる。定理 2B の一般ブロック版は現行のままでよい。
- **Cramer / 行列式比**: `χ = N/Δ` の形。定理 7A で結局この形に戻るので、
  2A の段階でこの形も併記しておくと 7A が独立の定理でなく系になる。

### 定理 3 — モーメント展開と最初の非零モーメント ★重要
現行: Neumann 級数、仮定 `q = ‖X‖/Γ < 1`。

- ★ **有理関数の Laurent 展開による証明**:
  `K_Γ = p†(ΓI−X)^{-1}ν = p† adj(ΓI−X) ν / det(ΓI−X)` は `Γ` の**有理関数**で、
  分子の次数 ≤ `N−1`、分母の次数 `N`。したがって `Γ → ∞` での Laurent 展開
  `K_Γ = Σ_{n≥0} M_n Γ^{−n−1}` は**ノルム条件なしに**成立する
  (収束半径は `1/ρ(X)`、スペクトル半径であって作用素ノルムではない)。
  よって `M_0 = … = M_{m−1} = 0, M_m ≠ 0` ⇒ `K_Γ = M_m Γ^{−m−1} + O(Γ^{−m−2})`
  は **`q<1` を仮定せずに**従う。
  現行の Neumann 証明は `Γ > ‖X‖` を要求するが、非正規な `X`(Liouvillian は
  一般に非正規)では `‖X‖ ≫ ρ(X)` になり得るため、この仮定は**実際に効く制約**。
  §5.9(F) で "nonnormal fast" を別枠扱いしているのは、まさにこの制約の帰結であり、
  Laurent 版に差し替えれば**その別枠自体が不要になる可能性がある**。
  Neumann 版は「明示的な打ち切り誤差評価」を与える点に価値があるので、
  補題として残し、**主定理の証明は Laurent 版に差し替える**ことを推奨する。
- **Schrieffer–Wolff / 断熱消去**: 物理側の同値な導出。次数勘定が
  「速い部分空間の消去の次数」として読めるので §5.2 の微視的マッピングと接続が良い。
- **Cauchy 積分による評価**: `M_n = (1/2πi)∮ Γ^n K_Γ dΓ` から、
  作用素ノルムでなくレゾルベントの一様評価で誤差限界が出る。非正規性に強い。

### 定理 4 — 厳密 transfer zero(Krylov / Cayley–Hamilton)★重要
現行: (i)⇔(ii) 自明、(i)⇒(iii) Cayley–Hamilton、(iii)⇒(iv) adjugate、
(iv)⇒(i) Neumann 級数 + べき級数係数の一意性。

- ★ **最小多項式による鋭化**: Cayley–Hamilton の代わりに
  `ν` に関する `X` の**最小多項式**(次数 `d ≤ N`)を使うと、
  検証すべきモーメントは `M_0, …, M_{d−1}` の `d` 個で足りる。
  現行の主張(`N` 個)は正しいが**鋭くない**。数値実装上も
  Krylov 空間の実効次元で止められる。
- ★ **Kalman 分解 / 可制御-可観測性**:
  (i)–(iv) は制御理論では「実現 `(X, ν, p)` の伝達関数が恒等的に零
  ⟺ 可到達部分空間が不可観測」という **Kalman 分解定理そのもの**である。
  これを採ると:
  - 証明が一行になる、
  - `d = deg(最小多項式)` の鋭い形が自動的に出る、
  - **`T1_sector_cut_axiomatization.tex` の定理 2.1(ii)「realization invariance
    ― 同じスカラー伝達関数の二つの実現は最小 Kalman 実現を共有する」が、
    外部から引用した仮定ではなく本パッケージ内で閉じた帰結になる。**
    現状は T1 定理 2.1 が別パッケージ (master-response no-go package) の
    定理 3.1 を「verbatim 適用」しているだけで、本リポジトリ内に証明がない。
    ここを閉じられるのは論文構成上大きい。
- ★ **(iv)⇒(i) の Neumann 依存の除去**: 定理 3 と同じく有理関数論法で置き換えれば、
  `K_Γ` が開集合上で消える有理関数 ⇒ 恒等的に零 ⇒ 無限遠での Laurent 係数が
  すべて零、で終わる。これにより**定理 4 が定理 3 に依存しなくなる**
  (929 行目の依存関係表の "4 ← 3" が消える)。定理の独立性が上がる。

### 定理 5 — 対称性に保護された transfer zero
現行: 可換射影子。

- ★ **Schur の補題 / Wigner–Eckart**: 対称群 `G` の既約表現が異なる
  `p` と `c` の間の行列要素は消える、という標準論法。現行の射影子版は
  **可換(アーベル)対称性に限定**されている。§11 の点群選択則は
  非可換群(`C₃ᵥ` 等)を扱うのだから、非可換版を主定理にしておくべき。
  現行版はその可換特殊ケースとして残る。
- **弱対称性 (weak symmetry) の言語**: 注の「微視的十分条件」
  (`UH_cU† = H_c`, `UL_μU† = Σ V_{μν}L_ν`)は Lindbladian の
  **weak symmetry** の定義そのもの(Buča–Prosen, *New J. Phys.* **14**, 073007 (2012);
  Albert–Jiang 2014)。strong symmetry(各 `L_μ` が個別に可換)との区別は
  定常状態の一意性に効くので、系(Full-Liouvillian version)が
  「`L_c` が一意な定常状態を持つ」を仮定している以上、**明示的に区別すべき**。

### 定理 6 — 特異減衰と保護されたチャネル
現行: `ker D ⊕ (ker D)^⊥` のブロック分解 + Schur + Neumann。

- **Drazin 逆の極限公式**: `lim_{Γ→∞} (ΓD + A_0)^{-1} = P(PA_0P)^{-1}P`
  は特異摂動の標準結果。証明が短くなる。
- **Kato の解析的摂動論**: `Γ^{-1}` を摂動パラメータとする正則摂動として
  レゾルベントを展開。高次項も系統的に取れる。
- ⚠ **スコープの注意(欠落している仮定)**: 現行証明は
  `D = D† ≥ 0` を使って `ker D ⊥ Ran D`(直交分解)を導いている。
  しかし**物理的な Liouvillian の減衰部は一般に非正規**であり、
  その場合 `ker D` と `Ran D` は直交しない。一般の場合は
  **斜交(スペクトル)射影子**が必要で、`P = Proj(ker D)` を
  直交射影子と書いている現行の主張は**そのままでは非エルミート `D` に拡張できない**。
  定理文に「`D` エルミート半正定値」を仮定として書いてあるので誤りではないが、
  §5.9(F)・§8.2 が非正規ケースを扱う以上、**適用範囲の注記を本文に上げるべき**。

### 定理 7A — adjugate 恒等式と解析接続
現行: 解析関数の一致の定理(ブラックボックスとして引用)。

- ★ **多項式版で十分**: 本プロジェクトで実際に使うモデル
  (NV 3E kernel、2g+2e/2g+3e、SC transfer、group-IV)はすべて
  `A(θ)` の成分が**パラメータの多項式**である。多項式の場合
  「非零多項式の零点集合は内点を持たない」は初等的に示せるので、
  **解析関数の一致の定理というブラックボックスを外せる**。
  一般解析族版は remark として残す。929 行目の「二つのブラックボックス」の
  一方が消える。
- **Zariski 位相の言葉**: 代数幾何的には「`N ≡ 0` は既約多様体上の閉条件」。
  複数パラメータ族を扱う §15(外部摂動による no-go の解除)と接続が良い。

### 定理 7B — 半単純特異点と極の相殺
現行: 不変部分空間ブロック形 + Neumann。

- **Drazin 逆の Laurent 展開**: `(A+sI)^{-1} = P_0/s + A^D + O(s)` は
  Drazin 逆の標準的な性質として引用できる。
- ★ **非半単純(Jordan)ケースの明示**: remark (3) で
  「Jordan ブロックでは `s^{-2}` 以上の極が出るので各 Laurent 係数を個別に調べよ」
  と述べているが、**一般の Laurent 展開を定理として書いていない**。
  index `k` の固有値 0 に対して
  `(A+sI)^{-1} = Σ_{j=1}^{k} (−1)^{j−1} N^{j−1} P_0 / s^{j} + A^D + O(s)`
  (`N` はニルポテント部)を定理 7C として書けば、
  remark が「注意書き」から「証明済みの手順」に格上げされる。
  数値的な `0/0` の扱いは本理論の実運用上の要なので、ここは書き切る価値がある。

### 系 — hopping rate map `γ_oc = Γ_XY/4`
現行: rank-1 の直接計算。

- **一般公式から**: 任意の jump `L` について、コヒーレンス `σ_ab` の減衰率は
  `(1/2)(⟨a|L†L|a⟩ + ⟨b|L†L|b⟩)` から recycling 項を引いたもの。
  これを一度書けば、対称・非対称・detailed balance の全ケース
  (現行では注で「同様に計算すれば」と済ませている非対称ケースを含む)が
  **一つの式から出る**。§5.2 の微視的マッピングが 3 ケース別々に書かれているのは
  この一般公式がないため。

## I-4. 新規性の切り分け(査読対策)

上の整理から明らかなように、**定理 1A・1B・2A・2B・5 はいずれも既知結果の
別表現**である(順に rank–nullity/Morris–Shore、Kraus/Albert–Jiang、
Feshbach–Fano/Schur、Schur 補元、weak symmetry による選択則)。
これは欠陥ではないが、現在の `eit_nogo_proofs.tex` は出典を一つも引かずに
すべてを自前で証明しているため、**「既知結果を再発見している」と読まれるリスクが高い**。

本当に新しいのは

- 抑制指数 `ν` による**排他的三分律** (`ν ∈ {∞} ∪ (0,∞) ∪ {0}`)、
- sector cut `χ − χ^(S)_cut` を **no-go の正しい対象として同定**したこと
  (定理 2B の注「`χ_full = 0` は no-go ではない」がこの主張の核)、
- そこから導かれる**材料非依存の分類**、

の 3 点である。証明文書の冒頭に「定理 1A–2B・5 は標準結果の本枠組みでの
再定式化であり、新規性は定理 3・4・6 と分類定理にある」と一段落置くだけで、
priority 論争をほぼ回避できる。

## I-5. T1(sector cut 公理化)固有の指摘

- **定理 2.1 (invariance) が本リポジトリ内で証明されていない。**
  別パッケージ ("master-response no-go closed package") の定理 3.1 を
  "applied verbatim" とするのみ。その文書はこのリポジトリに存在しない。
  論文添付時には自己完結しないので、I-3 の定理 4 の Kalman 版を採用して
  内部で閉じるか、当該文書を同梱するかの**どちらかが必須**。
- **定理 2.2 (L1.2, 非一意性) の反例が検証不能** → 後述 A-4。
- 定理 2.3 (L1.3, 吸収汎関数の符号) の (ii) は「generic な `K(z)` では」と
  述べつつ、証明は単一 Lorentzian の具体例のみ。`generic` の意味
  (どの集合で稠密か)が未定義。**主張を「単一 Lorentzian 極を持つ任意の
  受動応答について」に弱めるか、genericity を定義する**必要がある。
  結論(`θ=0` が正準)は正しいので、主張の書き方だけの問題。

---

# 第 II 部 計算コードの再現性監査

## II-1. 実施した検証

著者環境とはライブラリ世代の異なる clean 環境
(numpy 2.4.6 / scipy 1.17.1 / pandas 3.0.5)で、以下をすべて実行した。

| 対象 | 結果 |
|---|---|
| 全 pytest(リポジトリ root から / 各キャンペーン dir から) | **47 passed**(両方) |
| `New no-go theory` Phase A/B/D/M/P | 全 exit 0、JSON 一致(`runtime_seconds` のみ差) |
| Gate A / B / C / D | 全 exit 0、JSON 一致(`runtime_s` のみ差) |
| `RoomT` step1–step9 | 全 exit 0、JSON 一致(step5 の `runtime_s` のみ差) |
| `No-go theorem/scripts/reproduce_prl_figures.py`(full run) | exit 0、13 expected outputs 全て存在、JSON/CSV **完全一致** |
| `EIT definition equivalence/src/*.py`(9 本) | 全 exit 0、**1 本が不一致**(→ A-1) |
| `2g2e_package` の SHA256SUMS | 11/11 OK |
| 生成 PNG のピクセル比較(委員会全体) | Fig.4 を除き**すべてピクセル完全一致** |

**この一致度は特筆に値する。** 疑似乱数は例外なく
`np.random.default_rng(固定 seed)`(legacy `np.random.seed` の使用は皆無)、
`differential_evolution` は `workers=1` を明示、matplotlib は 30 箇所で
`use("Agg")`、出力パスはすべて `os.path.dirname(os.path.abspath(__file__))` 相対。
`scipy.stats.qmc.LatinHypercube` + `differential_evolution` を使う
RoomT step5 まで scipy 1.17 でビット単位一致した。

## II-2. 公表前に必ず直すべき欠陥

### A-1 【High】`gate_F5B_full_gksl.json` が修正前の古い結果

`EIT definition equivalence/src/full_gksl_2g3e.py` を再実行すると、
committed JSON と以下が食い違う。

| キー | committed | 再実行 |
|---|---|---|
| `min_Re_chi_matched` | **−0.17246** | **+0.012866** |
| `min_Re_chi_matched_at_delta` | 0.0 | −6.0 |
| `any_gain_in_matched` | **true** | **false** |

**gate 判定そのものが反転している。** 同一環境での再実行は完全に決定論的
(md5 一致)であり、環境差ではない。原因は明白で、
JSON とソースを同時に追加したコミット `f7c472e` のメッセージ自身が

> Before the fix, the simplest possible sanity check ... showed Re(chi) going
> negative — which would have been reported as a breakdown of Proposition A
> had it not been caught against the analytic Lorentzian.
> After the fix: the matched-response floor (Re(chi) >= 0) is confirmed ...

と書いており、**JSON は「fix 前」の値のまま commit されている**。
同梱の `results/F5B_findings.md` も
「min Re(χ_matched) = 0.0029 (>0、反例なし)」「floor は完全 GKSL でも成立」
と述べており、JSON とだけ矛盾する。

→ **JSON が誤り。再生成して差し替えること。** 論文に「受動性の floor は
完全 GKSL でも破れない」と書く根拠がこの gate である以上、
アーカイブが逆の値を記録したまま公開されるのは致命的。

なお本ファイルの数値経路自体には別の潜在的脆弱性がある
(`np.linalg.eig` の `argmin|w|` による定常状態選択、および特異行列への
`np.linalg.lstsq(rcond=1e-12)`)。今回の条件では
第 2 固有値が `−0.0816`、特異値ギャップが `9.8e−18` 対 `7.9e−02` と
十分離れているため安全だが、**パラメータを変えた再実行では LAPACK 実装差
(OpenBLAS / MKL / Accelerate)で選択が変わり得る**。
`steady_state()` に「零固有値が一意であること」の assert
(例: `|w[2nd]| / |w[1st]| > 1e6`)を入れておくことを推奨する。

### A-2 【High】原稿 Fig. 4 が `quick` モードの図

`Writing Paper/prl_figures/figures/fig4_robustness.{png,pdf}` の committed 版は
`build(quick=True)` の出力である。根拠:

- `tests/test_figures.py::test_fig4_outputs_and_crossover` だけが
  `build(quick=True)` を呼んでおり、pytest 実行後は committed PNG と
  ピクセル一致する。
- README が案内する `python make_figures.py`(既定 `quick=False`)を実行すると、
  パネル (a) 領域の **12,472 px が変化する**(他の 3 図は完全一致)。

差は見た目だけではない。`quick` は `run_gate_d.sc_approximate_class` の
`Γ` グリッドを 140 点 → 90 点に落とすため:

| 量 | quick=True(= committed 図) | quick=False(= README の手順) |
|---|---|---|
| `crossover_power`(**図中に印字される**) | **−1.0112** | **−0.9971** |
| `gamma_star` (最小 ε) | 4.6616e5 | 4.4779e5 |
| `gamma_star` (最大 ε) | 5.0379e8 | 4.2609e8(**−18%**) |

理論の主張は `Γ* ∝ ε^{−1}` なので、**full 実行の −0.9971 の方が主張に近い**。
committed 図は精度の低い方を載せている。

→ **`make_figures.py`(full)で再生成して差し替える**こと。併せて
`test_fig4` が `quick=True` で committed 成果物を上書きしてしまう構造も直す
(テストは一時ディレクトリに出力するか、`quick` 版を別名にする)。

### A-3 【Medium】`2g2e_package` が「再現パッケージ」として成立していない

`EIT definition equivalence/2g2e_package/` は SHA256SUMS 付きの配布パッケージだが、

- 同梱の `code/dark_state_free_2g2e_analysis.py` は全 53 行で、
  **`data/*.csv`(2 本)も `figures/*.png`(3 枚)も一切生成しない**。
  2 つの数値を print して `plt.show()` するだけ。
  `from pathlib import Path` は import されているが未使用で、
  ファイル出力コードが失われた痕跡に見える。
- したがって SHA256SUMS が保証しているのは**同一性であって再現性ではない**。
  第三者はチェックサムを検証できるが、CSV/PNG を作り直せない。
- さらにこのスクリプトは**リポジトリ内で唯一 `matplotlib.use("Agg")` を
  持たず `plt.show()` を呼ぶ**。GUI バックエンドが入った headless 環境では
  警告・ブロック・例外のいずれかになり、まさに「PC 環境の異なる他者」で
  壊れる典型形。

→ CSV/PNG を生成する完全版スクリプトを同梱し、`Agg` を指定し、
生成後に SHA256SUMS を張り直すこと。論文の Data Availability に
このパッケージを挙げるなら必須。

### A-4 【Medium】T1 定理 2.2 の反例数値が検証不能

`T1_sector_cut_axiomatization.tex` の定理 2.2 は
"Proof by explicit construction" として

```
χ_full(δ₀) = 0.8129 + 0.4262i
R_{S1}(δ₀) = −0.0498 + 0.0623i
R_{S2}(δ₀) = −0.0071 − 0.0008i
```

を提示するが、**これらの数値はリポジトリ内のコード・結果 JSON・
findings のいずれにも存在しない**(全文検索で 0 ヒット)。
定理本文は「gate F4/F5 で検証した 2g+3e matched coherence-space model」と
述べるが、その model から上記パラメータでこの 3 値を再現するスクリプトがない。

この定理は「理論は `χ` 単独でなく組 `(χ, S)` についての主張である」という
**理論の意味づけを決める中心的な主張**であり、その唯一の証明が
再現できない数値では査読を通らない。

→ 反例を生成する 20 行程度のスクリプトを
`EIT definition equivalence/src/` に追加し、結果を
`gate_T1_sector_nonuniqueness.json` として archive すること。

## II-3. 中〜低リスクの指摘

### B-1 依存パッケージにバージョン指定が一切ない
3 つの `requirements.txt` はすべてパッケージ名のみ(`numpy`, `scipy`, …)。
Python の下限も未指定。今回たまたま numpy 2.4 / scipy 1.17 で完全一致したが、
これは保証ではない。特に

- `scipy.optimize.differential_evolution` の内部戦略、
- `scipy.stats.qmc.LatinHypercube` のサンプル生成、
- `np.linalg.lstsq` の既定 driver、

は過去に scipy/numpy のマイナー更新で挙動が変わった実績がある。

→ 検証済み版で `>=` ではなく `==` に固定した
`requirements.lock.txt` を追加し、README に
「論文の数値は本 lock で再現される」と明記する。
`python_requires >= 3.11` も明示する(`from __future__ import annotations` を
使用しているため 3.7+ で動くが、検証済みは 3.11)。

### B-2 `runtime_s` / `runtime_seconds` が成果物 JSON に入っている
`gates_summary_phaseP.json`、Gate B/C/D/O、RoomT step5 の差分は**すべてこれだけ**。
JSON を diff で照合する運用(README がまさにそれを推奨)のノイズになる。

→ 実行時間は別ファイル(`runtime.json` 等)に出すか、
照合スクリプトで除外キーとして扱う。**現状、これが唯一の非決定要素であり、
除けば全キャンペーンがビット単位で決定論的になる。**

### B-3 PDF がバイト単位で再現しない(PNG は再現する)
matplotlib PDF backend は `/CreationDate` を埋め込むため、
内容が同一でも PDF のバイト列は毎回変わる。
`SOURCE_DATE_EPOCH` を設定すれば固定できることを実測で確認した
(設定値が `/CreationDate` に反映される)。

→ `prl_style.save()` 内で `SOURCE_DATE_EPOCH` 未設定時に固定値を
`os.environ` に入れるか、README に
`SOURCE_DATE_EPOCH=<定数> python make_figures.py` と記載する。
PNG は既にピクセル完全一致するので、PDF だけの問題。

### B-4 `reproduce_prl_figures.py` が追跡外ディレクトリを生成し、その成果物が archive されていない
full run すると `No-go theorem/` 直下に `figures/`, `outputs/`, `tables/` の
3 ディレクトリが作られる(`run_prl_prediction.py` L84 が `ROOT/'outputs'` 等に書く)。
中身は fig1–fig6(PRL 候補予測図)、`*.npz` 4 本、`*.csv` 4 本。

問題点:
- `results/` 配下ではないため、リポジトリの「成果物は `results/` に archive する」
  という規約から外れている。
- `.gitignore` にも入っていないので、実行するたび `git status` が汚れる。
- **`reproduce_prl_figures.py` の `EXPECTED` リストにも含まれていない**ため、
  これらが生成されたかどうかは検証されない。
- 何より、これら fig1–fig6 は**リポジトリに commit されていない**。
  他の全キャンペーンが図を archive しているのに、ここだけ照合対象がない。

→ 出力先を `results/figures`, `results/tables`, `results/outputs` に統一し、
`EXPECTED` に追加し、生成物を commit する。

### B-5 テストのカバレッジが campaign 間で不均一
`tests/` を持つのは GateB / GateC / GateD / PhaseO / No-go theorem / prl_figures の 6 つ。
一方、

- `New no-go theory/RoomT/`(step1–9、9 スクリプト)
- `EIT definition equivalence/src/`(9 スクリプト)— **A-1 が起きた場所**
- `New no-go theory/src/` の Phase A/B/D/M/P

には **1 つもテストがない**。README は
「Tests live in the `tests/` directory of each campaign」と書いており、
現状と一致しない。

→ 最低限、各 gate の PASS/FAIL 判定値を committed JSON と照合する
regression test を追加する。A-1 は、この 1 本があれば commit 時点で検出できていた。

### B-6 モジュール名の衝突(潜在)
`make_figures.py` が `Writing Paper/prl_figures/src/` と `No-go theorem/src/` の
両方に存在する。`fig1_classes.py` は `sys.path.insert(0, ...)` で
`No-go theorem/src` を**自分のディレクトリより前に**挿入するため、
以後 `import make_figures` は No-go theorem 側を掴む。
現在は `make_figures.py` が常に `__main__` として実行されるため顕在化していないが、
将来 import 経由で使うと壊れる。

→ どちらかを改名する(例: `make_prl_figures.py`)。

### B-7 `font.family: serif` の指定
`prl_style.py` が `serif` を指定しているが具体的なフォント名がない。
フォント構成の異なるマシンでは matplotlib が別フォントにフォールバックし、
PNG がピクセル一致しなくなる(今回の環境ではたまたま一致した)。

→ 具体名(例: `DejaVu Serif`、matplotlib 同梱で全 OS 共通)を指定するか、
`mathtext.fontset` と併せて明示する。

### B-8 未 commit の gate 結果
`model_2g3e_closedloop.py` は `gate_E4_2g3e_closedloop.json` を生成するが、
このファイルは commit されていない(他の E/F gate JSON はすべて commit 済み)。

→ archive するか、生成しないようにするかを決める。

## II-4. 良好であり、そのまま維持すべき点

再現性の観点で明確に良くできている点を、後退させないために列挙する。

1. **RNG が完全に統一されている** — 全 22 箇所が `np.random.default_rng(seed)`。
   グローバル状態を汚す `np.random.seed` / `random` モジュールの使用はゼロ。
   seed は引数化され既定値が固定。
2. **`matplotlib.use("Agg")` が 30 ファイルで import 直後に呼ばれている**
   (例外は A-3 の 1 本のみ)。headless CI でそのまま動く。
3. **出力パスがすべて `__file__` 相対**。CWD 依存は B-4 の 1 箇所のみ。
   リポジトリ root からでも各 campaign dir からでも pytest が通ることを実測確認済み。
4. **`differential_evolution(workers=1)` の明示** — 並列化による
   非決定性を意図的に排除している。
5. **gate 判定が JSON に PASS/FAIL 付きで archive されている** —
   まさにこの設計のおかげで A-1 を検出できた。
6. **`2g2e_package` の SHA256SUMS** — 発想は正しい(A-3 で中身を補えば完成する)。
7. **ネットワーク・外部データ・環境変数への依存がゼロ**。全計算が自己完結。

## II-5. 推奨する修正の優先順位

| 優先 | 項目 | 見積 |
|---|---|---|
| 1 | A-1: `gate_F5B_full_gksl.json` 再生成 | 5 分 |
| 2 | A-2: Fig. 4 を full モードで再生成 | 5 分 |
| 3 | A-4: T1 反例スクリプト追加 + JSON archive | 1 時間 |
| 4 | A-3: `2g2e_package` の生成スクリプト完全化 + `Agg` + SHA 再計算 | 1–2 時間 |
| 5 | B-1: `requirements.lock.txt` + Python 下限明記 | 15 分 |
| 6 | B-5: RoomT / EIT-definition の regression test 追加 | 2–3 時間 |
| 7 | B-2 / B-3: `runtime_s` 分離、`SOURCE_DATE_EPOCH` 固定 | 30 分 |
| 8 | B-4 / B-6 / B-7 / B-8 | 各 15–30 分 |

証明側は、論文提出前に最低限:

- 定理 1B の当該一段落を `H_eff` 経路に書き換える(I-2)
- 定理 3・4 を Laurent / Kalman 版に差し替え、`q<1` と `N` の仮定を外す(I-3)
- 定理 5 を Schur の補題版に一般化する(§11 の点群選択則との整合)(I-3)
- 定理 6 に `D` 非エルミートの場合の適用範囲注記を入れる(I-3)
- 新規性の切り分け段落を冒頭に置く(I-4)
- T1 定理 2.1 を自己完結させる(I-5)

を推奨する。1 番目はギャップの解消、2–4 番目は結果の強化と
適用範囲の明確化、5–6 番目は査読対策である。
