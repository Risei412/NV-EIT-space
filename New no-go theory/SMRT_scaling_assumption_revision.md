# SMRT 再統合メモ：固定スケーリング仮定から介入経路依存分類へ

**文書の役割:** 本文書は、当初の SMRT（Sector-Mediated Response
Theory）の固定スケーリング仮定をどこまで保持し、どこを拡張すれば Phase N
の強い主張を理論へ再統合できるかを整理する。既存の Theorem I--III や
operational sector cut の定義を置き換える完成原稿ではなく、再統合のための
仕様書である。

## 1. 結論

当初の SMRT は、強散逸を一つの固定された一変数経路

\[
A_\Gamma(z)=\Gamma D+B(z)
\]

として定義する。この仮定の下では、抑制次数は固定された組

\[
(D,B,c,p,S)
\]

に対する一つの整数として分類される。

Phase N で得られたより強い結果を取り込むには、この仮定を「任意の
スケーリングを許す」という形で撤廃するのではなく、次のように拡張する。

> **改訂された固定経路原則:** 各分類ではスケーリング経路を事前に固定する。
> ただし、native dissipation の強度 \(\Gamma\) だけでなく、sector
> intervention の強度 \(\kappa\) も独立な物理スケールとして認め、固定した
> joint path \(\kappa=\kappa(\Gamma)\) ごとに分類する。異なる経路の比較は
> 許すが、データを見た後に経路を選んではならない。

したがって、変更後のクラスは系だけに付随する量ではなく、少なくとも

\[
\boxed{
(\text{system},\text{probe},\text{readout},
 \text{sector intervention},\text{joint scaling path})
}
\]

に付随する。元の固定経路はこの拡張理論の特殊な一断面として保存される。

## 2. 当初の仮定 A5

当初の A5 は次の内容である。

1. 大きくする物理量は \(\Gamma\) のみである。
2. 指定された散逸チャネルは同じ係数 \(\Gamma\) でスケールする。
3. \(D\) と \(B(z)\) の分割は固定し、\(B(z)\) は \(\Gamma\) に依存しない。
4. 何を \(D\) に含めるかを変更すれば、別の物理的問いとして扱う。

この仮定は Theorem I--III の誤りではなく、そこで分類する漸近方向を一意に
定めるための仮定である。特に

\[
F_\Gamma(z)=p^\dagger[\Gamma D+B(z)]^{-1}c
\]

の \(\Gamma^{-1}\) 展開、master moments、および
\(\deg_\Gamma Q-\deg_\Gamma N\) による有限判定はそのまま有効である。

## 3. operational cut が既に導入していた第二のスケール

物理的 sector cut は

\[
A_{\Gamma,\kappa}(z)
=A_\Gamma(z)+\kappa D_S
=\Gamma D+B(z)+\kappa D_S
\]

で実現される。ここで \(D_S\) は admissible な選択的 GKSL intervention
generator である。

既存の Theorem 0A は、\(\Gamma\) を固定したまま \(\kappa\to\infty\) として

\[
\chi_{\Gamma,\kappa}(z)
=\chi_{\mathrm{cut},\Gamma}^{(S)}(z)+O(\kappa^{-1})
\]

を与える。従来の SMRT が分類する ideal cut は、二重極限の順序を

\[
\boxed{
\lim_{\Gamma\to\infty}\left[
  \lim_{\kappa\to\infty}\chi_{\Gamma,\kappa}(z)
\right]
}
\]

と固定したものである。したがって、従来理論も形式上は二尺度を含んでいたが、
\(\kappa\) は分類変数ではなく、先に除去される正則化パラメータとして扱われて
いた。

既存の double-limit remark は、\(\kappa_*(\Gamma)=O(\Gamma)\) のため二つの
極限が自由には交換できないことを既に指摘している。Phase N は、この注意が
単なる技術的 caveat ではなく、観測される抑制次数を変える物理効果であることを
示した。

## 4. 仮定の具体的な変更

### 4.1 A5 を A5a と A5b に分割する

既存 A5 は、次の二つに分けるのが最も安全である。

**A5a（native scaling の固定）**

\[
A_\Gamma(z)=\Gamma D+B(z)
\]

における \(D\) と \(B\) の分割は物理的入力として固定する。Theorem I--III
を適用するとき、\(B\) は \(\Gamma\) に依存しない。この仮定は変更しない。

**A5b（intervention scaling protocol の明示）**

有限強度の operational sector response を分類するときは、\(D_S\) と
\(\kappa(\Gamma)\) を追加の物理的入力として指定する。基本族として

\[
\kappa(\Gamma)=\kappa_0\Gamma^q,
\qquad \kappa_0>0,
\qquad q\ge0
\]

を用い、\(q\) を極限の前に固定する。より一般には
\(\kappa(\Gamma)=\kappa_0\Gamma^q[1+o(1)]\) のような regularly varying
path まで許せるが、最初の再統合では power-law path に限定するのがよい。

### 4.2 「固定」の意味は維持する

この変更は fixed-scaling principle の放棄ではない。各 \(q\) については

\[
\Gamma\longmapsto
\Gamma D+B(z)+\kappa_0\Gamma^qD_S
\]

という一つの経路が固定されている。変更点は、唯一の固定経路だけを理論の対象
とするのではなく、事前指定された固定経路の族を比較可能にすることである。

禁止すべき操作は次のとおりである。

- 観測データに最もよく合う \(q\) を事後的に選び、それを固有クラスと呼ぶこと。
- \(\Gamma\) の掃引中に、明示せず \(D\)、\(B\)、\(D_S\)、source、readout
  を変更すること。
- 異なる経路で得た次数を、同一の固定スケーリング問題の矛盾と解釈すること。

## 5. 再定義すべき応答と次数

有限介入に対する operational sector response を

\[
R_S^{\mathrm{op}}(\Gamma,\kappa;z)
:=\chi_{\mathrm{full},\Gamma}(z)
-\chi_{\Gamma,\kappa}^{(S)}(z)
\]

と定義する。固定 path \(\kappa=\kappa_0\Gamma^q\) 上の次数は

\[
\nu_S(q)
:=-\lim_{\Gamma\to\infty}
\frac{\log|R_S^{\mathrm{op}}
(\Gamma,\kappa_0\Gamma^q;z)|}{\log\Gamma}
\]

で定義する。ただし極限が存在し、評価点が漸近的な pole 上にない場合とする。
周波数窓 \(K\) 上で分類するときは、従来どおり
\(\|R\|_K=\sup_{z\in K}|R(z)|\) を用い、分母の一様可逆性と先頭係数の
一様非消失を仮定する。

ideal-cut order は別記号で

\[
\nu_S^{\mathrm{ideal}}
:=\nu\!\left[
\chi_{\mathrm{full},\Gamma}
-\lim_{\kappa\to\infty}\chi_{\Gamma,\kappa}^{(S)}
\right]
\]

と定義する。この区別により、有限介入の path order と従来の SMRT order を
同一視する誤りを避けられる。形式的には ordered limit を \(q=\infty\) と呼んで
よいが、有限 \(q\) の通常の極限と同一ではないため、証明では
\(\nu_S^{\mathrm{ideal}}\) と明記する。

## 6. 二尺度 path-order theorem

有限次元系では、固定した \(z\) における応答は二変数有理関数として

\[
R_S^{\mathrm{op}}(\Gamma,\kappa;z)
=\frac{N(\Gamma,\kappa;z)}
       {Q(\Gamma,\kappa;z)}
\]

と書ける。多項式

\[
P(\Gamma,\kappa;z)
=\sum_{a,b}p_{ab}(z)\Gamma^a\kappa^b
\]

に対し、path \(\kappa=\kappa_0\Gamma^q\) 上の exact weighted degree を

\[
d_{q,\kappa_0}(P)
=\max\left\{
w:\sum_{a+qb=w}p_{ab}(z)\kappa_0^b
\not\equiv0
\right\}
\]

と定義する。同じ重みを持つ Newton face 上の係数は、最大値を取る前に必ず
合算する。これにより breakpoint 上の厳密な干渉消去を取り込む。

分子が恒等的にゼロでなく、分母の先頭 face が評価領域で消えないなら、

\[
\boxed{
\nu_S(q;\kappa_0)
=d_{q,\kappa_0}(Q)-d_{q,\kappa_0}(N)
}
\]

である。有限個の monomial しか存在しないため、\(\nu_S(q)\) は有限判定可能で、
Newton face の交替または厳密消去が起こる点を除いて区分線形である。

この定理は Theorem I--III を否定しない。Theorem I--III は一変数断面の
係数・moment 構造を与え、path-order theorem は二変数応答のどの Newton face
が指定経路で支配的になるかを選ぶ。

## 7. Phase N による厳密な witness

同一の五準位 GKSL diamond、同一の source/readout、同一の native scaling、
同一の joint sector \(S_{34}\) に対して、Phase N は

\[
\boxed{
\nu_{S_{34}}(q)=
\begin{cases}
4-q,&0\le q\le1,\\
2+q,&1\le q\le2,\\
4,&q\ge2
\end{cases}
}
\]

を exact Gaussian-rational arithmetic で与える。

この結果の意味は次のとおりである。

- ideal cut では、二つの branch の leading moment が逆位相で消去し、次数は
  \(\nu_{S_{34}}^{\mathrm{ideal}}=4\) に昇格する。
- matched scaling \(\kappa\sim\Gamma\) では有限介入補正がその消去を覆い、
  観測次数は \(3\) になる。
- \(q\ge2\) では理想的な次数 \(4\) が回復する。
- よって次数は有限介入強度に対して単調とは限らず、同じ sector に対しても
  joint scaling path を指定しなければ一意でない。

ここで重要なのは、\(\kappa/\Gamma\to\infty\) と
「ideal sector order を保存する」が同値ではないことである。前者は operational
cut response の**絶対誤差**を小さくするための条件である。一方、ideal sector
response 自体が干渉により \(O(\Gamma^{-4})\) まで小さくなる場合、その高次構造を
保存するには有限介入誤差がさらに小さくなければならない。Phase N ではその閾値が
\(q=2\) として現れる。

## 8. 変更前後の理論構造

| 項目 | 当初の SMRT | 再統合後の SMRT |
|---|---|---|
| native generator | \(\Gamma D+B(z)\) | 変更なし |
| sector cut | \(\kappa\to\infty\) の ideal cut | ideal cut に加えて有限 \(\kappa\) を分類対象化 |
| 極限 | \(\kappa\to\infty\) を先に取る ordered limit | ordered limit と \(\kappa=\kappa_0\Gamma^q\) を区別 |
| class の帰属先 | system/probe/readout/sector | system/probe/readout/sector/joint path |
| suppression order | 一つの整数 \(\nu_S^{\mathrm{ideal}}\) | path-indexed family \(\nu_S(q;\kappa_0)\) |
| 判定法 | moments または一変数 polynomial degree | 二変数 weighted Newton degree を追加 |
| 干渉の役割 | moment 消去による次数昇格 | 昇格に加え、有限介入 face との競合を決定 |

## 9. 何を保持し、何を修正するか

### そのまま保持するもの

- 有限次元、Markovian、weak-probe、passive という基本 scope。
- native generator の固定分解 \(A_\Gamma=\Gamma D+B\)。
- Theorem I の exact-zero/Krylov certificate。
- Theorem II の moment hierarchy と cancellation promotion。
- Theorem III の protected-channel expansion と semisimple-kernel 仮定。
- sector intervention の GKSL admissibility、noninvasiveness、固定 source/readout。
- ideal cut を operational Zeno limit として定義する Theorem 0A。

### 修正または追加するもの

- A5 を A5a/A5b に分割し、joint scaling protocol を明示的入力に加える。
- ideal order と finite-intervention path order の記号を分離する。
- class の帰属先に joint scaling path を追加する。
- 二変数 rational response と weighted Newton degree の定理を追加する。
- double-limit remark を、単なる注意から path-resolved classification の導入部へ
  昇格させる。
- 「sector order は一つの系固有量である」という表現を避ける。

## 10. 強くなった主張と論理的境界

再統合後に主張できる中心命題は次である。

> **Intervention-path non-intrinsicness.** 干渉で昇格した sector 抑制次数は、
> 無介入の全応答または sector label だけに内在する一意な整数ではない。有限強度
> の物理的 GKSL intervention では、その値は事前指定された native/intervention
> joint scaling path に依存し、有限個の Newton-face transition によって変化する。

その帰結として、無介入の全応答だけから
\(\nu_S^{\mathrm{ideal}}\) または \(\nu_S(q)\) を一意に推定することはできない。
sector intervention と scaling protocol が入力に含まれていないからである。

ただし、次のより強すぎる命題は主張しない。

> sector intervention、scaling protocol、および tomographically complete な
> generator がすべて既知でも次数は同定不能である。

この命題は一般には偽である。必要な入力がすべて与えられれば、上記の有限
polynomial certificate により次数を計算できる。

## 11. RISEI との境界

この拡張は、複数の介入強度とその極限順序を扱うため RISEI と接するが、直ちに
RISEI へ統合されるわけではない。ここで扱うのは一つの指定 sector に対する
native/intervention 二尺度応答であり、ordered intervention composition、
irreducible residual、\(k_{\min}\)、functional hierarchy、observable rotation
は用いていない。

したがって当面は、名称を例えば

- **path-resolved SMRT**、または
- **two-scale operational SMRT**

とし、RISEI の多介入順序構造とは分離して記述するのが安全である。

## 12. 再統合の推奨順序

1. assumptions 節の A5 を A5a/A5b に置換する。
2. operational cut 節で ideal order と finite path order を定義し分ける。
3. Theorem 0A の double-limit remark の直後に path-order theorem を置く。
4. Theorem II の cancellation promotion を、二変数 numerator の face
   cancellation と接続する corollary を追加する。
5. Phase N の \(S_{34}\) fan を最小の物理的 witness として提示する。
6. taxonomy の表を \(\nu_S\) から \(\nu_S(q)\) へ拡張し、ideal-cut 列を別に残す。
7. abstract と conclusion では「固定経路を廃止した」と書かず、
   「classification を事前指定された joint paths の族へ拡張した」と書く。

## 13. 原稿で使える最短の定式化

> Original SMRT assigns a suppression order to a fixed native-dissipation
> path, \(A_\Gamma=\Gamma D+B\), with the sector cut defined by an ordered
> strong-intervention limit. We retain that theory as the ideal-cut slice and
> extend its fixed-path principle to jointly specified native/intervention
> paths, \(\kappa=\kappa_0\Gamma^q\). The resulting path-resolved order is a
> weighted Newton-degree difference and is therefore finitely decidable and
> piecewise linear. A physical GKSL realization exhibits a nonmonotone
> V-shaped order fan, proving that interference-promoted sector suppression is
> not an intrinsic integer of the uncut response.

## 14. 関連する既存ファイル

- `Theorem and proofs/three_theorems_proofs.tex`: 元の一変数 Theorem I--III。
- `Theorem and proofs/no_go_theory_undergrad_guide.tex`: A5 の明示的記述。
- `Sector/section0_operational_foundations.tex`: Theorem 0A と double-limit
  remark。
- `full_response_nonidentifiability.md`: Phase N の非同定命題、path-order
  theorem、厳密な V-shaped fan。
- `src/run_phase_n.py`: exact two-variable polynomial certificate と Gates
  N1--N5。

