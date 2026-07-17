# Operational sector-cut 改訂:着手手順と証明戦略

**対象:** `SMRT_sector_cut_physical_uniqueness_revision.md` に記載された Theorem 0 群
(0A–0E)の追加と、既存文書・数値コードの改訂
**前提文書:** `three_theorems_proofs.tex`(Theorem I–III)、
`eit_nogo_proofs.tex`(Theorem 2B)、master response package(augmented realization)

---

## 1. 中核となる観察:Theorem 0A は Theorem III の役割読み替えである

現行 Theorem III は

\[
A_\Gamma = \Gamma D + B(z)
\]

に対し、(H1) 零固有値の半単純性、(H2) $D_Q$ 可逆、(H3) protected block
$B_P = PBP$ の $K$ 上一様可逆、(H4) $\sup_K\|B\|<\infty$ の下で

\[
F_\Gamma = p^\dagger P B_P^{-1} P c + \Gamma^{-1}F_1 + O(\Gamma^{-2})
\]

を証明済みである(Riesz 分裂 + Schur complement + Neumann 級数)。

改訂案の operational cut は

\[
A_{\Gamma,\kappa}(z) = A_\Gamma(z) + \kappa D_S
\]

の $\kappa\to\infty$ 極限なので、Theorem III に次の**置換**を施すだけで
Theorem 0A が得られる:

| Theorem III | Theorem 0A |
|---|---|
| 散逸スケール $\Gamma$ | cut 強度 $\kappa$ |
| 散逸行列 $D$ | cut generator $D_S$ |
| コヒーレント部 $B(z)$ | full 生成子 $A_\Gamma(z)$(固定 $\Gamma$) |
| (H1) $0\in\mathrm{spec}(D)$ 半単純 | admissibility 条件 3(semisimple zero mode) |
| (H3) $PBP$ 可逆 | $P_S A_\Gamma(z) P_S$ の $K$ 上一様可逆性 |
| (H4) $\sup_K\|B\|<\infty$ | $A_\Gamma$ は固定 $\Gamma$ で $K$ 上有界(既存 A2 より自動) |
| $F_0$ | $\chi_{\mathrm{cut},\Gamma}^{(S)} = p^\dagger P_S [P_S A_\Gamma P_S]^{-1} P_S c$ |
| $\Gamma^{-1}F_1$ | $O(\kappa^{-1})$ 収束レート(Gate U1 の理論値) |

したがって **0A の証明は新規計算ゼロ**:Theorem III を仮定込みで引用し、
役割対応表を明示した Corollary として書く。$O(\kappa^{-1})$ の係数は
既存の $F_1$ 公式そのものであり、Gate U1 の定量予測にもなる。

### 注意点(証明を書く際の非自明箇所)

1. **誤差定数の $\Gamma$ 依存性。** Theorem III の $\Gamma_*$ に相当する
   $\kappa_*$ は $\sup_K\|A_\Gamma\|$ に比例し、これは $\Gamma$ に線形で増大する。
   よって $\kappa$-展開は「固定 $\Gamma$ で $K$ 上一様」であり、
   $\Gamma\to\infty$ と合成するには実質 $\kappa/\Gamma\to\infty$ の領域が必要。
   これは二重極限問題(Gate U7)が単なる形式的注意でなく実質的である理由で、
   Theorem 0A に remark として明記する。
2. **符号規約。** $A(z) = zI - \mathcal{L}$ 系の規約では $D_S$ の非零固有値の
   実部が正になる向きに取る(改訂案 §4.1 の規約をそのまま採用)。
3. **$z$ 依存性。** $A_\Gamma(z)$ は $z$ 依存だが、Theorem III の証明は
   $B(z)$ の $z$ 依存性を既に許容している($K$ コンパクト + 解析性)ので追加作業なし。

---

## 2. Theorem 0B(kernel universality)は 0A の系

0A の極限公式は $D_S$ に **$P_S$ を通じてのみ**依存する。よって
$P_S = P'_S$ なら極限は一致する。証明は 2 行。

### 重要な精密化:「kernel だけでは足りない」

半単純な零固有値に対する Riesz 射影は
$\mathrm{Ran}\,P_S = \ker D_S$、$\ker P_S = \mathrm{Ran}\,D_S$ で決まる
**斜交射影**である。つまり同値類は $(\ker D_S, \mathrm{Ran}\,D_S)$ の
**組**で決まり、kernel が同じでも range が異なれば $P_S \neq P'_S$ となる。
改訂案の "kernel-defined uniqueness" という表現はここを曖昧にしているので、
定義 0.3 は $P_S = P'_S$(射影そのものの一致)で書く。

そのうえで、次の **補題 0B′** を追加すると物理的に使いやすくなる:

> **補題 0B′(自己随伴 cut の場合)** $D_S$ が Hilbert–Schmidt 内積に関して
> 自己随伴(例:Hermitian jump operator による pure dephasing
> $\mathcal{D}[L]$, $L=L^\dagger$)なら、$P_S$ は $\ker D_S$ への
> **直交**射影であり、このクラス内では kernel の一致だけで射影の一致が従う。

証明:$L=L^\dagger$ のとき $\mathcal{D}[L]$ は HS-自己随伴
($\langle X, \mathcal{D}[L]Y\rangle_{HS} = \langle \mathcal{D}[L]X, Y\rangle_{HS}$
を直接計算)、よって対角化可能(半単純性が自動で従う=admissibility 条件 3 の
十分条件にもなる)かつ固有空間直交。§11 の Λ 模型の dephasing cut は
このクラスに入るので、最小数値実験では kernel 一致 ⇒ 射影一致が保証される。

---

## 3. 残る定理群の証明方針

### Theorem 0C(covariant representation invariance)

現行 Lemma 1(similarity invariance)の組 $(D,B,c,p)$ に $D_S$ と $P_S$ を
追加するだけ。$P_S \mapsto T P_S T^{-1}$ は Riesz 射影の定義
(レゾルベントの周回積分)と共変なので自動:
$(\zeta I - TD_ST^{-1})^{-1} = T(\zeta I - D_S)^{-1}T^{-1}$。
証明は既存の共役消去 1 行 + この 1 行。

### Proposition 0D(frozen / self-consistent 分解)

分解式自体は telescoping 恒等式なので証明不要。実質は次の補題:

> **補題 0D′** $\mathcal{L}_0(\rho_0)=0$ かつ $\mathcal{C}_S(\rho_0)=0$ なら、
> 全ての $\kappa\ge0$ で $(\mathcal{L}_0+\kappa\mathcal{C}_S)(\rho_0)=0$。
> さらに cut generator の定常状態が一意なら、cut model の定常状態は
> $\rho_0$ 自身であり、source $c = \mathcal{V}_p\rho_0$ が不変、
> よって $R_S^{\mathrm{reprep}}=0$。

証明は線形性 + 一意性のみ。一意性の仮定(cut 後の Lindbladian が
唯一の定常状態を持つ)は明示的な仮定として置く。$\kappa\to\infty$ 極限での
一意性の保存は自明でないため、「有限 $\kappa$ ごとに一意」を仮定に含める。

### Theorem 0E(class robustness)

漸近比較の 2 行証明:
$R_S^{(1)} - R_S^{(2)} = \chi_{\mathrm{cut}}^{(2)} - \chi_{\mathrm{cut}}^{(1)} = O(\Gamma^{-q})$、
$R_S^{(2)} = \Gamma^{-\nu}r_\nu + o(\Gamma^{-\nu})$、$q>\nu$ より
$R_S^{(1)} = \Gamma^{-\nu}r_\nu + o(\Gamma^{-\nu})$。
$\sup_K$ ノルムで書けば master response package の
Definition 7.1($\nu$ の定義)とそのまま整合する。

### 現行 block cut との一致(改訂案 §4.2)

$D_S = \mathrm{diag}(0, I_S)$ に 0A を適用すると
$P_S = \mathrm{diag}(I, 0)$(直交射影)で

\[
\chi_{\mathrm{cut}}^{\mathrm{op}} = p^\dagger P_S [P_S A P_S]^{-1} P_S c
= d^\dagger A^{-1} b_p
\]

($c, p$ が非 sector block に台を持つ場合)。これは Theorem 2B の
$B,C\to0$ cut と**転送関数として**一致する。生成子としては
$B,C\to0$ cut は decoupled な sector block を残すが、$c,p$ が触れないため
transfer function は同じ——この「生成子の一致ではなく応答の一致」である点を
remark に明記する(Kalman realization invariance と整合する主張になる)。

---

## 4. 実行順序

### Step 1:LaTeX への Theorem 0 群追加(最優先・ほぼ引用のみ)

`three_theorems_proofs.tex` の Section 1 の前に
「Section 0: Operational foundations of the sector cut」を追加:

1. Definition 0.1 実験仕様 $\mathfrak{E}$
2. Definition 0.2 admissible cut generator(GKSL / selective / semisimple /
   noninvasive / fixed prep-readout の 5 条件)
3. Definition 0.3 operational sector $[\mathcal{C}_S]$($P_S=P'_S$ で定義)
4. Theorem 0A + 役割対応表(Theorem III の Corollary として)
5. Theorem 0B + 補題 0B′(自己随伴の場合の kernel 十分性)
6. Theorem 0C(Lemma 1 の拡張)
7. Proposition 0D + 補題 0D′
8. Theorem 0E
9. Remark:誤差定数の $\Gamma$ 依存性と二重極限の非自明性(Gate U7 への接続)

### Step 2:companion への remark 追加

`eit_nogo_proofs.tex` Theorem 2B の直前に、$B,C\to0$ が
$\kappa\to\infty$ Zeno 極限の algebraic realization であるという remark
(§3 末尾の transfer-function equality を含む)を追加。

### Step 3:Λ 模型での最小数値検証(理論より先でもよい)

`core.py` に追加:

- `validate_cut_generator(L0, C_S, rho0)`:GKSL 形式チェック、
  零固有値の半単純性(固有分解 + Jordan チェック)、
  $\|\mathcal{C}_S(\rho_0)\|$ の数値検証
- `operational_cut_response(A_of_z, D_S, c, p, z, kappa)`:有限 $\kappa$ 応答
- `riesz_projection_exact` は既存のものを $D_S$ に流用
- `compare_cut_equivalence`:$\kappa$ 掃引 + $O(\kappa^{-1})$ フィット
- `self_consistent_cut_response`:cut Lindbladian の定常状態再計算 +
  `sector_response_decomposition`($R^{\mathrm{path}}$ / $R^{\mathrm{reprep}}$)

検証は §11 の 3 項目に限定:
(a) $\kappa\to\infty$ で bare Lorentzian へ収束、
(b) 誤差 $O(\kappa^{-1})$(理論係数 = $F_1$ 公式と比較)、
(c) algebraic / operational 両 cut の $\nu$ 一致。

ゲート優先順:U1 → U2(dephasing rate を変えた 2 実装 + 非自己随伴な
admissible 実装 1 つ)→ U5 → U3/U4 → U6 → U7。

### Step 4:README・数値優先度文書・QFI 名称変更

- `README.md` の一文説明を改訂案 §9.3 の文へ差し替え
- `new_nogo_numerical_priorities.md` に保存量リスト
  ($\chi^{\mathrm{alg}}, \chi^{\mathrm{op}}(\kappa), \chi^{\mathrm{sc}},
  R^{\mathrm{path}}, R^{\mathrm{reprep}}$)を追加
- metrology 系:$F_{Q,S} \to \mathcal{I}_S^{(G)}$(sector-mediated tangent
  distinguishability)への改名と $\Delta F_Q$ の併記。定理化の 4 仮定
  (faithful $\bar\rho$、metric の有界非退化、基準状態差の subleading 性、
  rank-deficiency 回避)を assumptions として明記

---

## 5. リスクと未解決点(着手前に合意しておくべき事項)

1. **$P_S$ の一致 vs kernel の一致**(§2):定義を $P_S=P'_S$ に固定し、
   kernel だけで十分なのは自己随伴クラスに限る、という構成にする。
   改訂案の "kernel-defined operational uniqueness" の表現は本文で微修正。
2. **二重極限の順序**:標準規約($\kappa\to\infty$ が先)を定義に焼き込み、
   非可換性は Gate U7 + 将来課題として切り出す。誤差定数の $\Gamma$ 線形性から
   「$\kappa \gg \Gamma$ で一様」という定量的な安全領域が言えるので、
   これを remark として先に確保しておく。
3. **$P_S A_\Gamma(z) P_S$ の可逆性**:0A の実質的仮定。Λ 模型では retained
   block が bare Lorentzian で $\Im z$ 側に極が逃げるため $K$ 上可逆だが、
   一般模型では検証必須(`validate_cut_generator` に条件数チェックを入れる)。
4. **多準位模型での選択性**:§11 の警告どおり、$\mathcal{D}[L_S]$ が他の
   optical coherence にも作用しうるため、「$P_S$ が意図した保持空間に一致するか」
   を毎回数値検証する(kernel の次元と基底の明示チェック)。
