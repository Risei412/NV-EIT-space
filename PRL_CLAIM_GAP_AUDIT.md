# PRL主張に必要な計算の網羅性監査 (2026-08-01)

`SIMULATION_PLAN.md` が定める Gate 1–5 の**仕様項目**と、`No-go theorem/src/` の
**実装および実際の出力ファイル**を1対1で突き合わせた結果。

`SIMULATION_PLAN.md` §28 の status update は「Gates 1–5 implemented, executed,
and ALL PASSED」と記録しているが、これは各スクリプトが自ら出力する
`*_all_passed` フラグに基づく。本監査は**そのフラグが検査していない仕様項目**を
洗い出したものである。各 `gate*_all_passed = true` は、下記の欠落項目を
検査対象に含めていないという意味で、仕様の充足を意味しない。

## 総括

| Gate | 仕様項目 | 完全実装 | 部分的 | 未実装 |
|---|---|---|---|---|
| 1 EIT/ATS識別 | 6 | 4 | 2 | 0 |
| 2 縮約 vs フル | 6 | 1 | 4 | 1 |
| 3 SNR変換鎖 | 5 | 0 | 5 | 0 |
| 4 閾値の不確かさ帯 | 5 | 3 | 2 | 0 |
| 5 アンサンブル平均 | 8 | 3 | 4 | 1 |
| 6 論文図 Fig.1–4 | 4 | 0 | 3 | 1 |
| 7 テスト | 7 | 7 | 0 | 0 |
| 7 CI要件 | 3 | 2 | 0 | 1 |
| 8 投稿凍結 | 5 | 0 | 1 | 4 |

投稿前に**必ず**埋めるべき欠落（主張の成立に直結）は §「Blocker」の6件。

§7 のテストは7項目すべてが実装済みで、`python -m pytest tests/ -q` は
**23 passed (1.89 s)**。仕様の中で最も充実している節である。

---

## Blocker（主張そのものが未検証）

### B1. Gate 1 — 4モデル比較の合否判定が存在しない
`gate1_candidate_aic_bootstrap.py:334-343` の `passed` 辞書は
ΔAIC = AIC(ATS) − AIC(EIT) の1組しか見ていない (`:165`)。
一方、実際の出力 `results/tables/gate1_aic_bootstrap.json` の Akaike weight は

    Fano 0.47 > Lorentzian 0.40 > EIT 0.13 > ATS 0.00

で、**EIT は最良モデルではない**。仕様 §1 の目標文
「identified as EIT (not ATS, not a single Lorentzian, not Fano)」は数値上
満たされていない。`results/gate_1_5_report.md:33-38` に caveat として
記載はあるが、合否基準としては未検査。
→ ATS を退けるだけでは「EIT である」と主張できず、より単純な単一ローレンツ／
Fano で説明できないことを示さないと「単なる非対称ディップ」に格下げされる。

### B2. Gate 2 — opening exponent の候補点での再確認が存在しない
仕様 §2 の合否基準第1項「opening exponent identical (re-confirm at candidate
point)」に対し、`gate2_candidate_full_vs_reduced.py` に exponent/slope の
計算は一切ない（grep ヒット0）。指数検証は `bperp_full_lowfield_slope.py` に
あるが、そこは **T=300 K, Bz=0.02 T, Bx∈[0,0.03] T, Ω_c=1.0 GHz**
(`:22-23`, `bperp_full_validation.py:13`) で走っており、候補点
（70 K, Bz=0.005 T, Bx=0.232 T, Ω_c=0.1 GHz）ではない。
→ 指数2の一致は「縮約が全 Liouvillian の忠実な縮約である」ことの唯一の定量的
証拠。フォノン率が桁違いで Bx が摂動域外の候補点へは外挿できない。

### B3. Gate 2 — orbital branch トグルが未実装
仕様 §2 は singlet shelving / ISC / hyperfine / **両方の orbital branch** の
4トグルを要求し、合否基準は「for every toggle configuration」。実装の
`configs` (`:196-199`) は `base / isc / hyperfine / isc+hyperfine` の4行だが
これは2トグルの直積で、`j0` 引数 (`:76`) は常に `J0=3` の単一ブランチ固定。
加えて singlet shelving と ISC は単一フラグ `isc` (`:76`) に抱き合わせで
分離不能。
→ 主張は特定の 3E 分岐に紐づく干渉なので、他方の分岐で符号が反転すれば
偏光/歪み平均で信号が相殺し得る。分岐依存性なしに「観測可能」とは言えない。

### B4. Gate 3 — δχ → δα の換算が実装されていない
仕様 §3 の鎖 δχ → δα → ΔOD → ΔT/T → ΔN_ph → SNR のうち、**初段が欠落**。
`signal_chain.py` に δχ を引数に取る関数はなく、`sigma_zpl_cm2`(`:22-28`) +
`alpha_cm`(`:34-37`) は2準位断面積から**絶対**吸収係数を作るだけ。
相対コントラスト C を `delta_od`(`:42-45`) が直接 OD に掛ける形で、
「C は δα/α と同一物」という未検証の仮定で迂回されている。
さらにフォノンサイドバンド背景が `signal_chain_parameters.csv` にも
`od_total`(`:86-87, 108-109`) にも入っていない（ZPL は全体の ~3.5 % なので
背景吸収とショットノイズを最も楽観的に見積もっている）。
→ 「理論量が実際の吸収変化になる」という第一歩が未証明で、σ_ZPL の1桁誤りが
そのまま SNR に1桁乗る。

### B5. Gate 4 — T_no-go の区間が存在しない
`gate4_threshold_uncertainty.py:47` の `THRESHOLDS` は
T_1pct / T_0.1pct / T_0.01pct / T_sign の4つのみで、**T_no-go は定義も抽出も
されていない**（リポジトリ全体で `T_no-go`/`T_nogo` の実装ヒット0）。
出力 `gate4_threshold_bands.csv` も4行。
→ no-go 主張そのものの温度境界に区間が付かないと、論文の中心主張が単一数値の
まま残り、Gate 4 の目的（区間化）を主張の核心部分で満たさない。

### B6. §6 — 論文 Figure 1–4 が1枚も合成されていない
`scripts/reproduce_prl_figures.py` は仕様どおり「薄いラッパ + assert」だが、
STEPS (`:17-24`) が呼ぶのは `run_prl_prediction.py` と `gate1`–`gate5` の
6本のみ、EXPECTED (`:26-40`) の13ファイルもすべて `gate*_*.png/csv/json` で、
**`fig1`–`fig4` という論文図は1つも assert されていない**。§6 の表が指定する
`bperp_kernel_map_v2.py` / `bperp_full_lowfield_slope.py` /
`moment_order_common_pipeline.py` は STEPS に入っておらず、前2者は
そもそも `savefig` を持たない（json/npz のみ）。
→ Fig.1（Liouvillian ブロック構造・sector S・Schur 補・B⊥=0 禁制 vs
B⊥≠0 escape channel）は **完全に未実装**：`src/` を "Schur" /"schematic" /
"block structure" で grep してヒット0。論文の中核メカニズムを図で
提示する手段が存在しない。

**混同注意**: `Writing Paper/prl_figures/` の `fig1_classes.py` …
`fig4_robustness.py` は**別論文の図**である。`NV_EIT_PRA_PRL_Split_Strategy_
20260724.md` によればこちらは PRA/PRL 分割戦略の「一般理論PRL」
（Gate A–D: 整数クラス／observable 継承／物質独立性／robustness）の
Fig.1–4 で、SIMULATION_PLAN §6 が要求する NV-EIT（Gate 1–5）の Fig.1–4 とは
中身が全く別。「Fig.1–4 は既にある」と誤認しやすいが、§6 の図は依然0枚。

---

## Gate 1 — EIT/ATS 統計的識別

| 仕様項目 | 判定 |
|---|---|
| 候補スペクトル生成（Gate 2と同一機構） | 実装済 (`:56-75`) |
| 4モデルのフィット | 実装済 (`:78-90, 133-170`) |
| AIC/AICc/BIC/Akaike weights | 実装済 (`:126-131, 159-164`) |
| ロバストネス4種 | 部分的（下記） |
| 出力3ファイル | 実装済（3点すべて実在） |
| 合否基準のハードチェック | 部分的（下記） |

- **窓変動 +25 %（factor 1.25）が未実装**。`:322` の factors は
  `(0.5, 0.75, 1.0, 1.5, 2.0)` で、仕様「±25 %」の +側が欠け、代わりに仕様外の
  +100 % が入っている。→ 窓の非対称な取り方で判定が変わらないことを示すのが
  窓ロバストネスの目的なので、片側欠けは網羅性の穴。
- **σ/depth > 0.05 が 95 % 判定から除外**されている (`:336` の
  `if r['sigma_rel'] <= 0.05`)。実データは4水準すべて 100 % なので結論は
  変わらないが、実験に近い高ノイズ域が合否式に入っていない。
  → 実測 SNR での判定安定性が「実験で識別可能」という主張の根拠。
- 窓判定 (`:337`) は verdict の集合サイズのみを見て `best_model` の反転を
  許している（実際 CSV では Fano/Lorentzian/Fano と揺れる）。
- Ω_c スイープ全行で `best_model=Fano, w_fano=1.0`。Ω_c ≳ 0.84 GHz で適応窓が
  3.5–24 GHz に膨張し、1光子線幅 γ≈72 GHz と同スケールに入るため、
  docstring `:14-17` が正当化している平坦背景近似が崩れている。
- 「bootstrap」は残差リサンプリングでなくパラメトリックなガウス加法ノイズ
  (`:188`)。仕様文言とは整合するが統計用語としては誤り。
  bootstrap/window 内の `maxfev=5000` (`:188, :211`) はベースライン
  `:152` の 20000 より緩く、収束失敗が `rss=inf` に落ちて EIT 有利側に
  偏る潜在リスクがある。
- モデル関数とAIC式が `eit_ats_classifier.py:48-60`（k=4/3, 背景なし）と
  gate1 `:78-90`（k=5/4, 背景あり）で**二重実装**され、gate1 は classifier を
  import していない。片方だけ修正すると論文中の数値が不整合になる。

## Gate 2 — 縮約カーネル vs フル Liouvillian

Blocker B2 / B3 に加えて:

- **生スペクトルが保存されていない**。比較・出力されるのは派生量
  `C(δ₂) = (A_cut − A_full)/A_cut` のみで、`A_full` と縮約 δχ_S の配列は
  CSV/JSON に入らず、図 (`:239-248`) も C しか描かない。→ Gate 3 は δχ_S の
  絶対値から検出可能性を主張するので、比を取ると消える絶対振幅の一致が
  示されないと Gate 3 の入力が正当化されない。
- **ρ_12 が定常状態由来でない**。`:156` の `rho_cp_over_Op` は一次応答ベクトル
  `xf` の要素であり、仕様の「from the full steady state」を満たさない。
  絶対値のみで位相が失われ、値はグリッド中央の1点のみ (`:170`)。
  → EIT の機構は基底コヒーレンス ρ_12 の生成そのもので、定常 ρ_12 の
  大きさと符号の一致が「同じ物理」の証明の核心。
  （励起状態 population `:154` は定常状態由来で正しい。）
- **中心周波数に合否判定がない**。CSV には `center_full_MHz`/`center_red_MHz`
  列があるが `passed` 辞書 (`:222-228`) は contrast/linewidth/sign/survival の
  4項目のみ。→ 中心ずれは二光子共鳴条件のずれ＝縮約が捨てた項によるシフトで、
  実験の走査中心を決めるため予言として不可欠。
- `dip_metrics` (`:176-182`) は argmax + 半値超え端点差の**非パラメトリック
  推定**でフィットではない。FWHM は全構成でグリッド刻み 0.08333 MHz の整数倍に
  張り付き、中心差は 1 グリッド刻みそのもの（分解能限界に埋もれている）。
- hyperfine は実装されている (`:58-59, 89-90, 166-169`, mI=−1,0,+1 の3モデル
  平均) が、核スピンは Hilbert 空間に入らず A_perp・核ゼーマン・四重極は無視、
  mI 混合は表現不能。診断値 `diags[-1]` (`:174`) は mI=+1 のみを採用。
- 符号反転チェックは自動化されている (`:211, 225-227`) が**ピーク1点の符号のみ**
  で、スペクトル全域の点ごとの符号一致や副構造の保存は未検査。

## Gate 3 — δχ から測定可能信号へ

Blocker B4 に加えて:

- **不確かさが SNR に伝播されていない**。`signal_chain_parameters.csv` の
  source 列・uncertainty 列は全15行埋まっている（空欄/TODO なし、例
  `debye_waller,0.035,-,Santori 2010 / Doherty review (3-5%),+-0.01`）が、
  値は `'x/ 3'`, `'design'`, `'+-20%'` 等の**人間可読文字列**で、
  `gate3_snr_map.py:62` は `P[k][0]`（値）しか読まない。誤差バー付き SNR も
  feasibility 領域の不確かさ帯も計算されない（Gate 4 は Monte-Carlo を
  やっているが Gate 3 にはない）。→ 「SNR≥5 達成」が誤差1桁の中で意味を
  持つか示せない。
- **制御強度 Ω_c がパラメータ表にない**（`run_prl_prediction` 側に埋め込み）。
  出典も不確かさも記録されていない。
- **候補スペクトルではなくスカラー1点への適用**。`:58` で
  `rp.branch_value(...)` からピーク値 `C_cand = 0.0136` のみを取る。
  周波数軸上の ΔOD(δ)/ΔT/T(δ) プロファイルがない。→ 線幅とプローブ帯域の
  ミスマッチ（有限分解能での平均化損失）が評価されず信号を過大評価しうる。
- **厚さが変数になっていない**。`gate3_required_conditions.csv` は
  3コントラスト × 密度3点 = 9行をカバーするが `L_cm=0.05` 固定で、
  仕様の「必要 NV密度 × **厚さ** × 積分時間」の3軸のうち厚さ軸が欠落。
  合否基準の「≤ mm 厚」もスキャンではなく前提になっている。
- 数値の妥当性に懸念: 全9行 `feasible=True` で τ が 1.3e-4〜2.7e-6 s。
  τ が密度に対し非単調（0.01ppm 1.35e-4 s / 0.1ppm 2.7e-6 s / 1ppm 3.15e-5 s）
  で、OD_sector=7.9 による透過率 e^−7.9 の減衰と競合した結果。最適 OD(≈1) を
  選ぶ最適化がなく、μs 未満の τ に対する検出帯域・プローブ飽和・最小実用 τ の
  下限も課されていないため「1時間以内で余裕」が非物理的に強く出ている。
- 合否基準の実装が弱い (`:112-117`):
  - 基準1 `conditions_listed_all_targets = bool(len(rows) == 9)` (`:113`) は
    ループが常に9行生成するので**必ず true になるトートロジー**。τ=inf でも通る。
  - 基準3「Figures 3–4 の理論最適点が feasible 領域内」は**未実装**。代替の
    `candidate_feasible = (tau_c <= tau_max)` (`:116`) は 1 ppm 候補点の
    feasibility にすぎず、図のパラメータと照合するコードは存在しない
    （`optimum`/`optimal` の grep ヒット0）。
  - 基準2 のみ実質判定されている (`:114-115`)。
- 単体テストは `tests/test_core.py:146-155` のみで、`sigma_zpl_cm2`,
  `alpha_cm`, `od`, `delta_od`, `detected_photons`, `min_detectable_contrast`
  が**未テスト**（断面積・OD 段が数値的に無検証）。
- `signal_chain.min_detectable_contrast` (`:77-98`) は 300 K no-go の ε_th を
  定義する関数だが `gate3_snr_map.py` から一度も呼ばれず、その固定値は
  どの出力ファイルにも記録されていない（docstring 自身が「300 K 大域最適化の
  前に検出チェーンから固定せねばならない」と述べている）。
- 図タイトルは "shot-limited" (`:140`) だが実際は技術ノイズも含む（表記不整合）。

## Gate 4 — 温度閾値の不確かさ帯

Blocker B5 に加えて:

- **事前分布が 10 項目中 7 項目のみ**。`PARAM_NAMES` (`:48-49`) / `draw()`
  (`:51-62`) は9パラメータ（うち2つは仕様外の偏光角）。未実装は
  - **双極子行列要素** — `nv_model.P` の Dpar/Lpar/Dperp/Lperp は固定、
    双極子脚は単位ベクトル `U[:,1]`/`U[:,2]` (`:77`)
  - **ISC 分岐比** — モデルに ISC 自体がない
  - **不均一幅** — Gate 5 で別途扱われるが Gate 4 のバンドには入らない

  光学デコヒーレンス率も `scale_rad`（放射率スケール ×/1.2, `:54, :79-80`）で
  代理されている。→ これらは光学応答の振幅と線幅を直接支配するので、除外は
  バンド幅の過小評価につながる。
- **文献引用がどのパラメータにも無い**（docstring `:9-18` は分布形のみ）。
  唯一 `phonon_rates.py:4-13` に Happacher arXiv:2302.00011 の α_ph=1.70±0.08
  があるが、Gate 4 が使う ×/1.35 の σ はこの ±0.08（≒5 %）と整合せず出典不明。
- **`gate4_sensitivity.json` の `"priors": null` はバグ**。`:150` が
  `draw.__doc__` を読むが `draw` に docstring がないため、事前分布が
  成果物に一切記録されていない。
- NaN／打ち切りの扱いが未文書化。`:95` で `root = np.nan` 初期化、符号変化が
  なければ NaN、`:123` の `ok = np.isfinite(arr)` で除外し `n_valid` を記録。
  ただし「閾値が 20–200 K に存在しない」と「数値失敗」が区別されず、
  `TGRID` 上限 200 K 超の交差は無条件に NaN（警告なし）、非単調な C(T) での
  複数交差は最初のみ採用 (`:102`) で検証なし。現状 500/500 有効なので
  結果に影響はないが、T_no-go を追加すれば必ず問題化する。
  → 「一部サンプルで閾値が存在しない」は物理的に意味のある結果（no-go が
  成立しない領域）であり、暗黙に落とすと分位が生存バイアスを受ける。
- Spearman ランク付けは実装済 (`:138-148`, 4閾値すべて `dominant=scale_orb`,
  ρ≈−0.99) だが、対象がサンプリング済み9パラメータに閉じているため、
  未サンプルの双極子要素／ISC 分岐比／不均一幅が「支配的でない」ことは
  示されていない。
- **論文の有効数字がバンド幅に合わせて丸められていない**。
  `results/gate_1_5_report.md:71-73` に「T_sign ≈ 102 K, 68 % [95,109]、
  2桁以内」という指示はあるが、`Writing Paper/drafts/*.tex` には閾値温度の
  記述が一切なく（grep ヒット0）、丸めた本文が存在しない。一方
  `run_prl_prediction.py:222` の `primary_claim` は依然
  **"sign reversal at 101.44 K"** で、±7 K のバンドに対し過剰精度。
- 実装済みで問題ないもの: 固定 seed `SEED = 20260716` (`:45`)、
  n_samples=500 実行済（CSV/JSON ともに n_valid=500）、
  分位 16–84 % と 2.5–97.5 % の**両方** (`:124`)、出力2ファイル。
  ただし図 (`:159-174`) は 68 % バンドのみ描画で 95 % は未可視化。
  サンプリングは Latin-hypercube ではなく素の Monte-Carlo（仕様は
  「または」なので可）。

## Gate 5 — アンサンブル平均

- **温度勾配が未実装**。`member_spectra` の温度は `:78` の既定引数 `T=70.0`
  のみで、呼び出し (`:136-137`) は T を渡さず、`SCENARIOS` (`:98-118`) に
  温度関連キーが存在しない。`gamma = nv.gamma_oc_GHz(T, d)` (`:81`) は T に
  強く依存するのに全メンバーで同一。→ 主張自体が「70 K 窓」で、脱位相 γ は
  温度に鋭敏（Gate 4 で T_sign の 68 % 帯が ±7 K）。数 K の勾配で
  メンバー毎に γ が変われば追加的に減衰しうるので、これがないと
  「70 K で観測可能」の頑健性が未検証のまま残る。
- **MW 強度分布が未実装、かつ「空間プロファイル」になっていない**。
  実装は制御光強度のみ (`:129-130`,
  `Oc = rp.OC*exp(N(0, sigma_Oc_rel))`, σ=5 %/10 %) で、ビーム形状を持たない
  i.i.d. 対数正規スカラー乱数。ガウシアンビーム径・強度プロファイル・
  ビーム内位置の重み付けはない。MW/スピン駆動 Rabi 周波数はモデルに
  パラメータとして存在しない（probe/control が両方**光学**遷移、`:71`）。
  → 実測はビーム内の Ω_c 勾配で平均された信号で、Ω_c 依存の EIT 幅/
  コントラストは非線形なので、i.i.d. 散らしはビーム平均の washout を
  系統的に過小評価しうる（＝予測が楽観側にずれる）。
- **全 σ 値に文献根拠がない**（ファイル全体に文献引用0件）。σ_d = 0.1/0.3 GHz
  (`:127-128`)、σ_B = 0.1/1.0 G (`:124-126`)、σ_opt = 5/30 GHz (`:131-132`)、
  σ_Oc = 5/10 % がすべてマジックナンバー。→ 「0.1 G homogeneity が実験的に
  実装可能」という Decision B の条件は、その値が現実の実験値域であることの
  引用なしには成立しない。strain は `max(0.3, ...)` (`:127`) のクリップで
  分布が非対称に切断されるが注記もない。
- **平均化の数値積分の収束確認がない**。n_draws=60 固定 (`:102-117`)、
  `--quick` で10 (`:156`)。サンプル数倍化での変化率テスト、MC 標準誤差、
  複数シードの分散評価がすべて皆無で、`tests/` に gate5 のテストもない
  （`scripts/reproduce_prl_figures.py` はファイル存在を assert するのみ）。
  周波数グリッド（`:151`, 161点）の Gate 5 での収束確認もなし。
  → 60 ドローで Cmax ~1e-4 を報告しており MC 誤差が washout factor の桁と
  競合しうる。誤差付きでないと「×105 washout」「19 % survives」は定量的
  主張として支持できない。
- `metrics()` (`:141-146`) の FWHM は `|C| >= max|C|/2` を満たすインデックスの
  **最初と最後の差**なので、非連結領域や符号反転を含む場合に幅を過大評価する。
  `high_density` の fwhm=4.25 MHz / center=+2.125 MHz はこの定義由来の疑いが強い。
- spectral hole selection は棄却サンプリング (`:133-135`,
  `while abs(dopt) > opt_cut_GHz`) として実際にモデル化されている（言葉だけ
  ではない）が、ホール幅・埋め戻し・選択効率／収率はなく、選択で失われる
  原子数（S/N 代償）も出力されていない。
- 実装済みで問題ないもの: 4配向の正しいテトラヘドラル方向余弦
  (`:55-56`) と固定 lab-frame 磁場 (`:57-58`)、`frame_field()` (`:60-65`) の
  各 NV 枠への射影、静磁場不均一 (`:124-126`) と光学遷移周波数分布
  (`:131-132`, 一光子離調のみに入り二光子で相殺 — 物理的に正しい)、
  出力2ファイル（図は仕様の3本を超えて5本重ね描き）、
  3密度領域の行としての明示的区別と判定 (`:178-179`)、
  **washout factor の数値化** (`:170-172`, CSV 最終列: low 0.0421=×24,
  high 0.00948=×105, post_selected 0.0695, shimmed 0.1936)。
  軽微な近似として `member_context` (`:67-72`) は phi=0 固定・B⊥ をノルムのみで
  扱うため、各配向で B⊥ と strain 軸の相対方位角が変化しない。

## §6 — 論文 Figure 1–4

Blocker B6 に加えて、図ごとの判定:

| 図 | 判定 | 状況 |
|---|---|---|
| Fig.1 | **未実装** | schematic がコード化されておらず、`moment_order_common_pipeline.py` の inset 合成もない |
| Fig.2 | 部分的（データのみ） | slope は `bperp_kernel_map_v2.py:122-128`, `bperp_full_lowfield_slope.py:33-34` で計算・JSON化されるが両者に `savefig` がなく図は未合成 |
| Fig.3 | 部分的 | (T,B⊥) contrast カラーマップと 1e−2/1e−3/1e−4 contour は `run_prl_prediction.py:238-244` にある。**符号反転線は同図になく**別の1D図に分離 (`:249` の `axvline` 101.44 K)、**EIT/ATS 境界は未描画**（Gate 1 の Ωc 交差は ΔAIC vs Ωc の1Dのみ）、**visibility 境界（Gate 3由来）は完全に不在**（`gate3_snr_map.py:132` の SNR=5 フロンティアは (density, τ) 面にあり (T,B⊥) 面へ写像されていない） |
| Fig.4 | 部分的 | on/off スペクトルは `run_prl_prediction.py:253-257`、EIT/ATS フィット重ね描きは `gate1_candidate_aic_bootstrap.py:269-279`（2×2の左上）。**AICc/BIC を並べるパネルはなく**（`:129-131` で計算され JSON に入るだけ）、**透過+SNR パネルも別図**（しかも透過スペクトルではなく density×τ の SNR マップ）。4図とも1枚の PRL 図に組まれていない |

→ Fig.3 の visibility 境界と符号反転線がないと「理論的に開いている領域」と
「実際に見える領域」の区別＝no-go 主張の定量的な線が読者に示せない。
Fig.4 の SNR/透過パネルがないと feasibility 主張（Decision B）が本文図で
裏づけられない。

### 図生成の二重化・出力パス不整合

- **`Writing Paper/prl_figures` はクリーン clone でビルド不能**。
  `make_figures.py` が
  `../../../New no-go theory/PhaseO_observable_inheritance/src/gate_a_observable.py`
  を import し、そこで `import sympy` が失敗する（`No-go theorem/requirements.txt`
  には sympy があるが、こちらのディレクトリに requirements がなく依存が
  跨りディレクトリでリンクしている）。
- **出力パスが二重化している**。`run_prl_prediction.py:84` は
  `No-go theorem/{outputs,figures,tables}/` に書くが、commit 済みの成果物は
  `No-go theorem/results/{figures,tables,metadatas}/` にある。クリーン clone で
  再実行すると `No-go theorem/figures/fig1_branch_resolved_phase_map.png` が
  新規生成され、`results/figures/` の同名ファイルは更新されないまま二重化する。
  `reproduce_prl_figures.py` の EXPECTED は `results/` 側の gate 出力しか
  見ないため、この不整合を検出できない（削除した `results/figures/fig1..fig6`
  が復元されないまま "all 13 expected outputs present" と成功報告した）。
  → 投稿版の図がどのディレクトリの版か一意に決まらないと、原稿の図と
  再現スクリプトの出力が食い違うリスクが残る。
- `reproduce_prl_figures.py` は `Writing Paper/prl_figures` 側を一切呼ばず、
  「単一エントリポイント」要求に反する二重系統になっている。

## §7 — 自動テスト（仕様を満たしている）

`No-go theorem/tests/test_core.py`（167行、19テスト関数）。7項目すべて実装済み:

| 要求 | テスト関数 |
|---|---|
| 1 trace/Hermiticity 保存 | `test_steady_state_trace_hermiticity_positivity` (`:92-98`) |
| 2 positivity ≥ −1e−10 | 同上 (`:97`, `min(eigvalsh(rho)) >= -1e-10`) |
| 3 zero-control / zero-field / cut-sector 極限 | `test_zero_control_limit` (`:100`), `test_zero_field_limit` (`:105`), `test_cut_sector_limit` (`:117`) |
| 4 weak-probe 有限差分 | `test_finite_difference_response` (`:125-138`) — `first_order` と `steady_state(L+V)` を相対 1e−3 で比較 |
| 5 周波数グリッド収束 | `test_frequency_grid_convergence` (`:140-144`) — n=401 vs 801 で Cmax 変化 <1 % |
| 6 signal_chain 単位変換 | `test_signal_chain_units` (`:146-155`) |
| 7 固定 seed / 判定再現性 | `test_model_comparison_reproducibility` (`:157-162`) — `fit_all` 2回で ΔAIC と verdict が bit-identical |

`python -m pytest tests/ -q` → **23 passed (1.89 s)**。
軽微な点: 項目5は「グリッドを半分に」でなく「2倍に」（数学的に等価）。
警告1件 — `test_operational_cut_equivalence.py` の関数が None でなく tuple を
return（`PytestReturnNotNoneWarning`、assert 忘れの可能性）。

### CI 要件

- `matplotlib.use('Agg')`: **全 gate スクリプトにあり**（gate1–gate5,
  `run_prl_prediction.py`, および `New no-go theory/` の `run_gate_{a,b,c,d}.py`）。
- 実行時間: クリーン clone で `reproduce_prl_figures.py --quick` が **86 秒**で
  完走（laptop 級 <30 分を満たす）。フル実行は未計測。
- **`--n-samples` フラグは未実装**。`gate4_threshold_uncertainty.py:109` は
  `n_samples = 60 if quick else 500` のハードコードで、CLI は `--quick` のみ
  (`:180`)。gate1–5 いずれにも argparse がない。→ Monte Carlo 標本数を外から
  変えられないと、閾値バンド（T_sign = 102 K [95,109] 等）の収束性を
  第三者が検証できない。

## §8 — 投稿凍結

| 要求 | 判定 |
|---|---|
| コミットの tag | **未実装**（`git tag` の出力が空） |
| GitHub release | **未実装**（tag がないため） |
| Zenodo DOI | **未実装**（`.md`/`.tex` に記載なし） |
| 原稿への commit hash 記録 | **未実装**（`Writing Paper/drafts/*.tex` に記載なし） |
| クリーン clone からの全論文図再生成 | 部分的 — 実測で clone → `results/figures`,`results/tables` 削除 → `--quick` 実行 → exit 0（86 秒、13 出力すべて存在、T_sign 101.4 K など全 Gate の pass 判定も再現）。ただし再生成されるのは **Gate 診断図のみ**で、論文 Fig.1–4 は対象外なので「every paper figure」は満たしていない |

→ §6 の図合成スクリプトが存在しない限り、§8 の「clean clone で全論文図が出る」
は原理的に達成不能。

---

## 是正の優先順

1. **B1** Gate 1 の合否に「EIT の Akaike weight が Lorentzian/Fano を上回る」を
   追加し、上回らないなら SIMULATION_PLAN §1 の fallback（"transverse-field-
   induced narrow interference response" への再ラベル、roadmap §13 Phase 1）を
   発動する。これは Decision A そのもの。
2. **B2** Gate 2 に候補点での opening exponent 再確認を追加。
3. **B5** Gate 4 に T_no-go を追加（閾値の意味論と NaN 処理の明文化を含む）。
4. **B4** Gate 3 の δχ→δα を明示的な関数にし、PSB 背景を `od_total` に入れる。
5. **B3** Gate 2 の orbital branch トグルと、singlet shelving / ISC の分離。
6. Gate 3 の不確かさ伝播（Gate 4 の Monte-Carlo を再利用）＋厚さ軸のスキャン。
7. Gate 5 の温度勾配、ビーム空間プロファイル、n_draws 収束と MC 標準誤差。
8. Gate 4 の事前分布3項目追加、全事前分布の文献引用、`priors: null` バグ修正。
9. `run_prl_prediction.py:222` を含む主張文の有効数字をバンド幅に合わせる。
10. Gate 1 の窓 +25 %、σ/depth 全水準の合否算入、モデル関数の二重実装の解消。
11. **B6** NV-EIT 用の `make_prl_fig{1,2,3,4}.py` を新規作成し
    `reproduce_prl_figures.py` の STEPS/EXPECTED に追加。Fig.1 の Liouvillian
    schematic が完全に未着手で最大の穴。Fig.3 には符号反転線・EIT/ATS 境界・
    Gate 3 visibility 境界を (T,B⊥) 面へ写像して重ね描き、Fig.4 は
    AICc/BIC パネル＋透過スペクトル＋SNR を1枚に合成。
12. 出力ディレクトリを `results/` に統一（`run_prl_prediction.py:84`）し、
    `results/figures/fig1..fig6` を EXPECTED に追加。
    `Writing Paper/prl_figures` に requirements（sympy 含む）を明記し、
    どちらの PRL の図かを README 冒頭で曖昧さなく宣言。
13. `gate4_threshold_uncertainty.py` に `--n-samples`（argparse）を追加。
14. tag / GitHub release / Zenodo DOI / 原稿への hash 記録（上記完了後）。
