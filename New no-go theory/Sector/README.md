# Sector — operational cut foundations (Step 1 output)

`section0_operational_foundations.tex` は改訂案
(`Theorem and proofs/operational_cut_revision_plan.md`)の Step 1 の成果物です。
Theorem 0A–0E / Proposition 0D / Definition 0.1–0.3 を、既存の
`Theorem and proofs/three_theorems_proofs.tex`(Theorem I–III)と
`No-go theorem/Theorem and proofs/eit_nogo_proofs.tex`(Theorem 2B, Lemma 2)
を明示的に引用する形で、独立にコンパイル可能な LaTeX 文書として作成しました。

既存ファイルへの直接編集は行っていません。内容を確定させたのち、
`three_theorems_proofs.tex` の `\section{Framework...}` の直前に
本ファイルの内容(プリアンブルを除く)を統合するのが Step 1 の最終形です。

## 収録内容

- Definition 0.1: 実験仕様 $\mathfrak E$
- Definition 0.2: admissible cut generator(GKSL 性・選択性・半単純零モード・
  noninvasive 条件・prep/readout 固定の 5 条件)
- Lemma: Riesz splitting の restatement
- Definition 0.3: operational sector(同値類 $[\mathcal C_S]$、$P_S=P_S'$ で定義)
- **Theorem 0A**: 現行 Theorem III を $(\Gamma,D,B)\to(\kappa,D_S,A_\Gamma)$
  と役割読み替えして得る cut 極限公式。$\kappa_*(\Gamma)=O(\Gamma)$ である旨と、
  二重極限の非可換性リスクを remark で明記。
- **Theorem 0B**: kernel universality。ただし一致条件は $P_S=P_S'$(斜交射影の一致)。
- **Lemma 0B$'$**: $D_S$ が Hilbert–Schmidt 内積で自己随伴(例:Hermitian jump
  operator による pure dephasing)の場合に限り、$\ker D_S$ の一致だけで
  $P_S=P_S'$ が従うことの証明。
- **Theorem 0C**: 共変表現不変性。$P_S\mapsto TP_ST^{-1}$ が Riesz 射影の
  周回積分の共変性から正しい変換則であることを明示。
- **Proposition 0D** + **Lemma 0D$'$**: frozen/self-consistent 分解。
  noninvasive 条件と定常状態の一意性から $R_S^{\mathrm{reprep}}=0$ を導出。
- **Theorem 0E**: cut ambiguity に対する class robustness。
- 二重散逸スケールの標準規約($\kappa\to\infty$ を先に取る)の明記。

## 未解決・次段階に持ち越す点

- Remark 3.2(Theorem 0A 直後)で指摘した二重極限の非可換性は証明せず、
  Gate U7 の数値探索課題として残した。
- Lemma 0B$'$ は自己随伴クラスに限定した十分条件であり、一般の admissible
  cut generator に対する kernel-only 十分条件は未証明(改訂案 §2.2 の懸念は
  完全には解消していない)。

## Step 2 出力:`theorem2B_operational_realization.tex`

Theorem 2B(`eit_nogo_proofs.tex`)の block cut $B,C\to0$ が、
$D_S=\mathrm{diag}(0,I_S)$ を用いた operational cut($\kappa\to\infty$
Zeno 極限)と厳密に一致することを示す addendum です。既存ファイルへの
直接編集は行わず、挿入予定の remark/corollary を独立文書として作成しました。

- **Remark**(admissibility):$D_S=\mathrm{diag}(0,I_S)$ が Definition 0.2
  の admissibility 条件を満たすこと(Hermitian ⇒ Lemma 0B$'$ により
  $P_S$ が直交射影 $\mathrm{diag}(I,0)$ になる)。
- **Corollary(Algebraic–operational equivalence)**:Theorem 0A を
  この $D_S$ に specialize すると、正確に Theorem 2B の
  $\chi_{\mathrm{cut}}^{(S)}=\chi_0 d^\dagger A^{-1}b_p$ が retained-block
  公式として再導出される。さらに block 消去を直接計算し、
  $O(\kappa^{-1})$ 補正係数を閉形式
  $\chi_0 d^\dagger A^{-1}BCA^{-1}b_p$ として明示(Gate U1 の定量目標値)。
- **Corollary(frozen-source 十分条件)**:Lemma 0D$'$ をこの block 設定に
  適用し、$\mathcal C_S(\rho_0)=0$ + 定常状態一意性から
  $R_S^{\mathrm{reprep}}=0$ を導出。
- 末尾に `eit_nogo_proofs.tex` への挿入位置の指示を明記(本文・証明自体は
  無変更、Theorem 2B 直前への追記のみ)。

## Step 3 出力:`src/`, `results/`

既存の `New no-go theory/src/core.py` は変更せず(Riesz 射影などの関数は
import して再利用)、`Sector/src/` に operational cut 専用のモジュールを
新規追加しました。

- `operational_cut.py`:`validate_cut_generator`(admissibility の
  半単純性・Hermitian 性・selective damping を検査)、
  `check_noninvasive`(条件 C4)、`operational_cut_response`(有限 $\kappa$)、
  `ideal_cut_response`($\kappa\to\infty$ の Theorem 0A 式 (0.2) を Riesz
  射影から直接評価)、`compare_cut_equivalence`(Gate U1/U2 用の $\kappa$
  掃引・log-log フィット)。
- `model_lambda_operational.py`:Step 2 の Corollary で使った
  $D_S=\mathrm{diag}(0,1)$ による Λ 模型 block cut。閉形式の
  $O(\kappa^{-1})$ 係数 $-|\Omega_c|^2/(4A^2)$ も実装。
- `model_ground_coherence_lindblad.py`:Gate U5 用の 3 準位 Lindblad 模型
  (basis $g_1,g_2,e$)。$e\to g_1$ decay(rate $\Gamma$)+ 弱い
  $g_2\to g_1$ 緩和(定常状態の一意性確保のため追加、rate $\epsilon$)で
  $\rho_0=|g_1\rangle\langle g_1|$ を一意な定常状態にし、
  sector cut を $g_2$ 上の pure dephasing($L_S=|g_2\rangle\langle g_2|$、
  Hermitian ⇒ Lemma 0B′ 適用)として実装。
- `run_gates_step3.py`:Gate U1・U2・U5 を実行し、
  `results/gates_summary_sector.json` に保存するスクリプト。

### 実行結果(2026-07-17 時点、全 PASS)

| Gate | 内容 | 結果 |
|---|---|---|
| U1 | $\chi_{\mathrm{cut}}^{\mathrm{op}}(\kappa)\to\chi_{\mathrm{cut}}^{(S)}$、slope $-1$、閉形式係数 $-|\Omega_c|^2/(4A^2)$ との一致 | slope $=-0.99997$、係数相対誤差 $4.4\times10^{-4}$ |
| U2 | 同じ $P_S$ を持つ異なる admissible $D_S$(実レート vs 複素レート)が理想極限で一致 | 理想極限差 $=0.0$(有限 $\kappa$ では差あり、$\kappa\to\infty$ で消える) |
| U5 | $\mathcal C_S(\rho_0)=0$ かつ $\kappa$ 全域で定常状態が $\rho_0$ に固定(Lemma 0D′) | $\|\mathcal C_S(\rho_0)\|=0$、定常状態偏差 $=0.0$(全 $\kappa$) |

再現:`python3 "New no-go theory/Sector/src/run_gates_step3.py"`
(要 `numpy`)。

## Step 4 出力:`step4_amendments.md`, `src/gate_u8_naming_audit.py`

`New no-go theory/README.md` と `new_nogo_numerical_priorities.md` は
直接編集せず、`step4_amendments.md` に適用可能な差分(コピー&ペースト用の
置換テキスト)としてまとめました。

- **README.md 冒頭の一文説明**:改訂案 §9.3 の一文
  (operational sector / universal strong-intervention limit / 分類対象は
  system 単体でなく (system, probe, readout, sector intervention) の組)
  を追加する置換案。
- **`new_nogo_numerical_priorities.md`**:frozen-source cut の規約直後に
  挿入する保存量リスト($\chi_{\mathrm{cut}}^{\mathrm{alg}}$,
  $\chi_{\mathrm{cut}}^{\mathrm{op}}(\kappa)$, $\chi_{\mathrm{cut}}^{\mathrm{sc}}$,
  $R_S^{\mathrm{path}}$, $R_S^{\mathrm{reprep}}$, $R_S^{\mathrm{total}}$)と、
  Gate U1 の検証結果への参照。
- **QFI 命名変更**:README.md の Metrological extension 節の該当段落を、
  $F_{Q,S}$ → $\mathcal I_S^{(G)}$ / `I_sector_metric`(sector-mediated
  tangent distinguishability)へ置換する案。$\Delta F_Q$ とは一般に異なる
  量であることを明記し、定理化に必要な 4 仮定(faithful $\bar\rho$、
  metric の有界非退化、基準状態差の subleading 性、rank-deficiency 回避)
  への参照を追加。

### Gate U8(命名監査)の実装と結果

`src/gate_u8_naming_audit.py` を新規実装し、2 部構成で検証:

1. **合成 qubit 反例**(PASS):$\bar\rho=\mathrm{diag}(0.7,0.3)$、
   $\partial_\theta\rho_{\mathrm{full}}$ の非対角成分 $a=0.3$、
   $\partial_\theta\rho_{\mathrm{cut}}$ の非対角成分 $b=0.1i$ を用いた
   手作りの反例。正しい双線形恒等式
   $\Delta F_Q=\mathcal I_S^{(G)}+2\times\text{cross term}$
   (`qfi_cross_term` の規約を数値的に検証して確定)の下で、
   $\Delta F_Q=0.32\neq\mathcal I_S^{(G)}=0.40$ を確認し、両者が一般には
   異なる量であることを実証。
2. **既存の物理 Lindblad 模型**(`model_metro_lindblad.py`、無改変、
   情報提供用):`kappa_cut=0` 構成では cut Hamiltonian が対角のみになり
   $\partial_\theta\rho_{\mathrm{cut}}\equiv0$ となるため cross term が
   厳密にゼロ——これは README の Gate M3($\nu_F/(2\nu_x)=1.0007$)の
   結果と整合する構造的説明であり、命名区別への反証ではない。

再現:`python3 "New no-go theory/Sector/src/gate_u8_naming_audit.py"`
(結果は `results/gate_u8_naming_audit.json` にも保存)。

## 未着手(この改訂の範囲外)

- `README.md` / `new_nogo_numerical_priorities.md` への実際の適用
  (`step4_amendments.md` の差分をレビューのうえ手動反映)。
- 改訂案 §14 の「将来課題」(Γ と κ の二重極限非可換性、複数 sector 同時
  cut、non-Markovian 拡張)。

(`New no-go theory/Theorem and proofs/operational_cut_revision_plan.md` の
§4 は Step 1–4 ですべて完了。)
