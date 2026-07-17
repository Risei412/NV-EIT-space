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

## 次のステップ(Step 2 以降、未着手)

1. `eit_nogo_proofs.tex` Theorem 2B 直前への remark 追加。
2. `core.py` への `validate_cut_generator` 等の実装と Λ 模型での
   Gate U1/U2/U5 の数値検証。
3. README・数値優先度文書・QFI 命名の改訂。

(`New no-go theory/Theorem and proofs/operational_cut_revision_plan.md` の
§4 Step 2–4 に対応。)
