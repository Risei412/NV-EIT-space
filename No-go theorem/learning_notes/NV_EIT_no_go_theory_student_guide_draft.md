# NV–EIT No-Go理論 学習ガイド

> **Draft v0.1**  
> 対象読者：大学学部レベルの線形代数と量子力学を一通り学習済みの学生  
> 対応する定理パッケージ：`../Theorem and proofs/eit_nogo_proofs_v7.tex`

---

## 0. このノートの目的

この理論が答えようとしている問いは、ひとことで言えば次である。

> **強い散逸の下で、ある入力からある観測量へ情報や振幅がどの程度伝わるのか。さらに、その伝達がゼロになる理由を、対称性・経路・実験検出限界に分けて理解できるか。**

NV中心のEITは重要な具体例だが、理論の数理的な骨格はNVやEITだけに依存しない。有限次元の線形応答問題を

\[
A_\Gamma=\Gamma D+A_0
\]

と書き、散逸強度 \(\Gamma\) を大きくしたときの応答

\[
K_\Gamma=p^\dagger A_\Gamma^{-1}c
\]

を調べる。

- \(c\)：外部から入るsource
- \(p\)：readout、すなわち何を測るか
- \(D\)：強くスケールする散逸部分
- \(A_0\)：コヒーレント結合、弱い散逸、detuningなど、\(\Gamma\) と共に増大しない部分

このノートでは、数式を暗記するよりも、次の流れを理解することを目標とする。

1. EITにおける「暗い状態」には複数の意味がある。
2. no-go判定には、全応答ではなく経路を切った反実仮想との差を使う。
3. resolventのLaurent展開係数が、入力から出力へ到達する経路を数える。
4. sector graphは低次係数の消失を保証する。
5. 最短経路どうしが干渉すると、graph distanceよりさらに高次へ昇格する。
6. exact no-go、asymptotic no-go、experimental no-goは別物である。

---

## 1. 前提知識と推奨する読み方

### 1.1 前提知識

このノートでは、以下を既知とする。

- 複素ベクトル空間、内積、随伴 \(A^\dagger\)
- 固有値、固有ベクトル、kernel、rank
- 射影演算子
- 行列のblock分解と逆行列
- 密度行列とLindblad方程式の基本形
- 幾何級数

Schur complement、Krylov空間、Drazin逆行列は本ノート内で必要な範囲だけ説明する。

### 1.2 推奨する読み順

初読では次の順番でよい。

1. 第2節：三種類のno-go
2. 第4節：線形応答の基本形
3. 第5節：Laurent展開とmoment
4. 第6節：sector graph
5. 第8節：NV–EITへの接続

証明を追いたい場合は、その後に第3、7、9節へ進む。

---

## 2. 最初に区別すべき三種類のno-go

この理論では、no-goという言葉を一枚岩として扱わない。

### 2.1 Exact no-go

\[
K_\Gamma=0
\]

が、regularな全パラメータ領域で恒等的に成立する場合である。

これは「有限の \(\Gamma\) では小さい」という主張ではなく、入力から出力へのtransferが代数的に完全に消えていることを意味する。

典型的な原因は次の二つである。

- sourceとreadoutが異なる保存sectorに属する
- readoutがsourceから生成されるKrylov空間全体に直交する

### 2.2 Asymptotic no-go

\[
K_\Gamma\sim \Gamma^{-r},
\qquad \Gamma\to\infty
\]

のように、強散逸極限で応答が消える場合である。

有限の \(\Gamma\) では応答は一般にゼロではない。したがって、これはexact no-goではない。

理論の中心問題は、整数 \(r\) が何によって決まるかである。

### 2.3 Experimental no-go

数理的には応答が存在しても、次の三つの窓が重ならない場合である。

\[
I_{\rm engineerable}
\cap I_{\rm asymptotic}
\cap I_{\rm detectable}
=\varnothing.
\]

- 装置で実現できる領域
- 理論の漸近則が成立する領域
- 信号が検出可能な領域

NV中心の室温EITでは、この第三の意味が特に重要になる。応答が厳密にゼロであることと、実験的に見えないことを混同してはいけない。

---

## 3. 「dark state」は一つの概念ではない

EITを学ぶとき、まず混乱しやすいのがdark stateという言葉である。

### 3.1 Optical dark vector

ground manifoldからexcited manifoldへの光学結合を

\[
\Omega:\mathcal H_g\to\mathcal H_e
\]

とする。状態 \(|D\rangle\) が

\[
\Omega|D\rangle=0
\]

を満たせば、これはoptical dark vectorである。

その次元はrank-nullityより

\[
\dim\ker\Omega
=\dim\mathcal H_g-\operatorname{rank}\Omega
\]

で与えられる。

ここでは「光がexcited manifoldへ結合しない」ということだけを述べている。

### 3.2 Stationary pure Lindblad state

密度行列

\[
\rho_D=|D\rangle\langle D|
\]

がLindblad方程式の定常状態であるためには、さらに強い条件が必要である。

各jump operatorについて

\[
L_\mu|D\rangle=\lambda_\mu|D\rangle
\]

が成立し、修正Hamiltonian

\[
\widetilde H
=H+\frac{i}{2}\sum_\mu
\left(\lambda_\mu^*L_\mu-\lambda_\mu L_\mu^\dagger\right)
\]

について

\[
\widetilde H|D\rangle=h|D\rangle
\]

でなければならない。

### 3.3 なぜ違うのか

\(\Omega|D\rangle=0\) でも、ground-state dephasingが \(|D\rangle\) の重ね合わせを壊すことがある。

したがって、

> optical dark vectorが存在する  
> \(\not\Rightarrow\)  
> 系が純粋なdark stateに定常化する

である。

この違いは、低温での理想的EITと、散逸を含む実際のNV系を区別する上で重要である。

---

## 4. no-go判定量は「全応答」ではなく「経路を切った差」

### 4.1 block線形応答

弱いprobeに対する一次応答を

\[
\begin{pmatrix}
A&B\\
C&G_g
\end{pmatrix}
\begin{pmatrix}
x\\s
\end{pmatrix}
=
\begin{pmatrix}
b_p\\0
\end{pmatrix}
\]

と書く。

- \(x\)：optical coherenceなど、速い応答変数
- \(s\)：ground coherenceなど、注目する長寿命sector
- \(B,C\)：速いsectorと長寿命sectorを結ぶ結合
- \(G_g\)：長寿命sector内部のgenerator

第一行から

\[
x=A^{-1}b_p-A^{-1}Bs
\]

である。これを第二行へ代入すると

\[
S_gs=-CA^{-1}b_p,
\qquad
S_g:=G_g-CA^{-1}B
\]

を得る。\(S_g\) がSchur complementである。

### 4.2 pathway-cut counterfactual

注目するsectorを通る経路だけを切るため、\(B=C=0\) とした反実仮想を考える。他のパラメータは固定する。

- full response：経路を残した応答
- cut response：対象経路だけを切った応答

その差を

\[
\delta\chi_S
:=\chi_{\rm full}-\chi_{\rm cut}^{(S)}
\]

と定義する。

Schur complementを使うと

\[
\delta\chi_S
=\chi_0d^\dagger A^{-1}B
S_g^{-1}CA^{-1}b_p
\]

となる。

### 4.3 なぜ全応答のゼロではいけないのか

理想的EITでは、干渉によって全吸収がゼロになることがある。

\[
\chi_{\rm full}=0
\]

は、EIT経路が存在しないことではなく、むしろEITが完全に働いた結果かもしれない。

一方、

\[
\delta\chi_S=0
\]

は、そのsectorを通る経路を残しても切っても観測量が変わらないことを意味する。

したがって、この理論でno-go判定に使うのは全応答ではなく、**pathway-resolved response contrast** \(\delta\chi_S\) である。

---

## 5. 強散逸展開：Neumann展開とLaurent展開

### 5.1 基本変形

\(D\) がinvertibleである場合、

\[
A_\Gamma=\Gamma D+A_0
=D(\Gamma I-X),
\qquad
X:=-D^{-1}A_0
\]

と書ける。

また

\[
\nu:=D^{-1}c
\]

と置けば、

\[
K_\Gamma
=p^\dagger(\Gamma I-X)^{-1}\nu
\]

となる。

ここから問題は、resolvent \((\Gamma I-X)^{-1}\) の大きな \(\Gamma\) における展開へ還元される。

---

### 5.2 Neumann展開

\[
(\Gamma I-X)^{-1}
=\frac1\Gamma
\left(I-\frac{X}{\Gamma}\right)^{-1}
\]

であり、

\[
\frac{\|X\|}{|\Gamma|}<1
\]

ならば

\[
\left(I-\frac{X}{\Gamma}\right)^{-1}
=\sum_{n=0}^\infty\frac{X^n}{\Gamma^n}
\]

と展開できる。

したがって

\[
K_\Gamma
=\sum_{n=0}^\infty
\frac{p^\dagger X^n\nu}{\Gamma^{n+1}}.
\]

Neumann展開の長所は、収束条件と打ち切り誤差を明示できる点である。

一方、非normalな行列では \(\|X\|\) が大きくなりやすく、この条件は必要以上に厳しいことがある。

---

### 5.3 Laurent展開

有限次元では

\[
K_\Gamma
=\frac{p^\dagger\operatorname{adj}(\Gamma I-X)\nu}
{\det(\Gamma I-X)}
\]

なので、\(K_\Gamma\) は \(\Gamma\) の有理関数である。

\(w=1/\Gamma\) と置くと

\[
K_{1/w}
=w\,p^\dagger(I-wX)^{-1}\nu.
\]

\(w=0\) では \(I-wX=I\) だからinvertibleであり、近傍でTaylor展開できる。

\[
(I-wX)^{-1}
=I+wX+w^2X^2+\cdots
\]

よって

\[
K_\Gamma
=\sum_{n=0}^\infty
\frac{M_n}{\Gamma^{n+1}},
\qquad
M_n:=p^\dagger X^n\nu.
\]

これが無限遠におけるLaurent展開である。

### 5.4 二つの展開の役割

- **Laurent展開**：漸近次数の存在と係数を示す主証明
- **Neumann展開**：具体的な収束領域と誤差評価を与える補助道具

したがって、Neumann展開が誤りなのではない。主定理に必要な仮定として使うと保守的すぎるため、役割を分離している。

---

## 6. projected Markov parametersとsector graph

### 6.1 momentの意味

\[
M_n=p^\dagger X^n\nu
\]

は、source \(\nu\) に \(X\) を \(n\) 回作用させた後、readout \(p\) にどれだけ到達したかを測る。

この理論ではpath momentと呼んできたが、線形システムの標準語彙へ寄せるなら **projected Markov parameter** と理解できる。

- \(M_0=p^\dagger\nu\)：直接重なり
- \(M_1=p^\dagger X\nu\)：一回の作用で到達
- \(M_2=p^\dagger X^2\nu\)：二段階経路で到達

最初に非zeroとなるものを \(M_m\) とすると、

\[
K_\Gamma
=M_m\Gamma^{-m-1}
+O(\Gamma^{-m-2}).
\]

つまり、最初の非zero momentが強散逸指数を決める。

---

### 6.2 sector分解

空間を直交sectorへ分ける。

\[
I=\sum_\alpha P_\alpha.
\]

\(D\) と \(A_{\rm diag}\) が各sectorを保存するとする。

\[
A_0=A_{\rm diag}+V.
\]

ここで

- \(A_{\rm diag}\)：sectorを変えない
- \(V\)：sector間を移動させる

と分ける。

\[
X_0=-D^{-1}A_{\rm diag},
\qquad
W=-D^{-1}V,
\qquad
X=X_0+W.
\]

sector \(\alpha\) から \(\beta\) へ

\[
P_\beta W P_\alpha\neq0
\]

なら、有向辺 \(\alpha\to\beta\) を描く。

---

### 6.3 sector-graph selection theorem

sourceがsector \(s\)、readoutがsector \(t\) にあるとする。

sector graph上の最短距離を

\[
d=\operatorname{dist}_G(s,t)
\]

とする。

このとき

\[
M_n=0
\qquad(n<d)
\]

が成立する。

理由は単純である。\(X_0\) はsectorを変えず、\(W\) 一回で高々一つの辺しか進めない。したがって、距離 \(d\) のreadoutへ到達するには最低でも \(d\) 回の \(W\) が必要である。

よって

\[
K_\Gamma=O(\Gamma^{-d-1}).
\]

これはgraph distanceが散逸次数の**下限**を与えることを意味する。

---

### 6.4 重要：graph distanceは自動的に等号を与えない

距離 \(d\) の最短経路が複数あると、それらの複素振幅が打ち消し合うことがある。

\[
M_d=0
\]

となれば、最初の非zero momentはさらに高次になる。

したがって、

> graphは「最低何回sectorを移る必要があるか」を決める。  
> 実際の係数が残るかは、経路振幅の干渉が決める。

という二段階構造になっている。

この点が、単なるgraph reachabilityと量子干渉を接続する部分である。

---

## 7. 具体例

### 7.1 三つのsectorを直列につないだ場合

basisを

\[
|s\rangle,
|a\rangle,
|t\rangle
\]

とし、

\[
X=
\begin{pmatrix}
0&0&0\\
g_1&0&0\\
0&g_2&0
\end{pmatrix},
\qquad
\nu=|s\rangle,
\qquad
p=|t\rangle
\]

とする。

sector graphは

```text
source s  ──g1──▶  a  ──g2──▶  readout t
```

であり、距離は \(d=2\) である。

実際、

\[
M_0=\langle t|s\rangle=0,
\]

\[
M_1=\langle t|X|s\rangle=0,
\]

\[
M_2=\langle t|X^2|s\rangle=g_1g_2.
\]

したがって \(g_1g_2\neq0\) なら

\[
K_\Gamma
\sim\frac{g_1g_2}{\Gamma^3}.
\]

「二段階の移動」が「\(\Gamma^{-3}\) のamplitude」に対応している。最初の \(\Gamma^{-1}\) はresolvent全体から現れ、追加の二段階がさらに \(\Gamma^{-2}\) を与える。

---

### 7.2 二本の最短経路が干渉する場合

```text
             ┌──▶ a ──▶┐
source s ────┤          ├──▶ readout t
             └──▶ b ──▶┘
```

二本の長さ2の経路振幅を

\[
A_a=v_au_a,
\qquad
A_b=v_bu_b
\]

とすると

\[
M_2=A_a+A_b.
\]

もし

\[
A_a=-A_b
\]

なら、graph distanceは \(d=2\) のままだが

\[
M_2=0
\]

となる。

この場合、\(\Gamma^{-3}\) 項が消え、次のmomentを調べなければならない。

ここで起きているのは「経路が存在しない」のではなく、「存在する複数経路が干渉して低次係数を消した」という現象である。

---

## 8. exact zeroとKrylov空間

### 8.1 reachable Krylov space

source \(\nu\) から \(X\) を繰り返し作用させて得られる空間を

\[
\mathcal K(X,\nu)
=\operatorname{span}
\{\nu,X\nu,X^2\nu,\ldots\}
\]

とする。

これはsourceから動力学によって到達可能な方向をすべて集めた空間である。

### 8.2 exact transfer-zero criterion

readout \(p\) がKrylov空間全体に直交すれば、

\[
p^\dagger X^n\nu=0
\qquad(\forall n)
\]

である。

したがってLaurent係数がすべてゼロになり、

\[
K_\Gamma=0
\]

が全regular \(\Gamma\) で成立する。

有限次元ではKrylov空間の次元を \(r\) とすると、最初の \(r\) 個だけを確認すればよい。

\[
M_0=M_1=\cdots=M_{r-1}=0
\]

なら、以後のmomentもすべてゼロになる。

### 8.3 sector graphとの関係

sourceからreadoutへの有向pathが一つも存在しないなら、すべてのwordがreadoutへ到達できない。これはKrylov直交条件の十分条件であり、exact no-goを与える。

ただしexact zeroは、graphの非到達性だけでなく、より精密な代数的相殺によっても起こり得る。

---

## 9. symmetry sectorによるexact zero

行列 \(A\) がsector projectorと可換であるとする。

\[
[A,P_\alpha]=0.
\]

このとき \(A^{-1}\) も同じsectorを保存する。

sourceとreadoutが異なる直交sectorにあれば、

\[
p^\dagger A^{-1}c=0.
\]

これは通常の選択則の線形応答版と見なせる。

- sector graphに辺がないため到達できない
- symmetryが異なるためmatrix elementが消える
- Krylov空間がreadout sectorへ入らない

という三つの言い方は、この場合には同じ構造を異なる角度から見ている。

---

## 10. kernelの次数から観測量の次数へ

PRLの主張では、個々のkernelだけでなく、最終的なobservableの指数を求める必要がある。

Schur complement応答を

\[
\delta\chi
=\chi_0L(\Gamma)S_g(\Gamma)^{-1}R(\Gamma)
\]

と書く。

各因子が

\[
L\sim a\Gamma^{-n_L},
\qquad
R\sim b\Gamma^{-n_R},
\qquad
S_g\sim c\Gamma^{-n_S}
\]

なら、leading coefficientが非zeroである限り

\[
\delta\chi
\sim
\chi_0\frac{ab}{c}
\Gamma^{-\nu_{\rm obs}},
\]

\[
\nu_{\rm obs}=n_L+n_R-n_S.
\]

これをobservable-order inheritanceと呼ぶ。

### 10.1 正規化すると指数が変わる

contrastを

\[
C=\frac{\delta\chi}{\chi_{\rm ref}}
\]

と定義し、

\[
\chi_{\rm ref}\sim d\Gamma^{-n_{\rm ref}}
\]

なら

\[
C\sim \Gamma^{-(\nu_{\rm obs}-n_{\rm ref})}.
\]

したがって、

- signed absorption difference
- normalized contrast
- population
- amplitude

では指数が同じとは限らない。

論文で「指数は3である」と書く前に、何のobservableの指数なのかを明示する必要がある。

---

## 11. singular dampingとprotected channel

これまで \(D\) はinvertibleと仮定した。しかし、散逸が全方向を抑えるとは限らない。

\[
\ker D\neq\{0\}
\]

なら、散逸から保護されたsubspaceが存在する。

\(D=D^\dagger\ge0\) とし、\(P\) を \(\ker D\) への直交射影とする。\(PA_0P\) がinvertibleなら

\[
p^\dagger(\Gamma D+A_0)^{-1}c
=
 p^\dagger P(PA_0P)^{-1}Pc
+O(\Gamma^{-1}).
\]

つまり、sourceとreadoutがprotected subspaceを通じて結ばれていれば、\(O(1)\) の応答が残り得る。

### 11.1 適用範囲の注意

この証明は

\[
D=D^\dagger\ge0
\]

を使っている。

一般のLiouvillian blockは非Hermitianかつnonnormalであり、\(\ker D\) と \(\operatorname{Ran}D\) が直交するとは限らない。

その場合はRiesz projectorやDrazin inverseを用いた別の定理が必要であり、本理論のHermitian singular-damping theoremをそのまま適用してはいけない。

---

## 12. singular pointとLaurent主部

\(A\) がzero eigenvalueを持つ点では、通常の逆行列は存在しない。

zero modeがsemisimpleなら

\[
(A+sI)^{-1}
=\frac{P_0}{s}+A^D+O(s)
\]

となる。

transferは

\[
K(s)
=\frac{p^\dagger P_0c}{s}
+p^\dagger A^Dc+O(s).
\]

したがって、poleが消える条件は

\[
p^\dagger P_0c=0
\]

である。

Jordan blockがある非semisimpleな場合には、\(s^{-2},s^{-3},\ldots\) も現れる。数値計算で \(0/0\) が出ても、それをzero responseと解釈してはいけない。Laurent主部の全係数を調べる必要がある。

---

## 13. NV–EITへどう接続するか

NV中心では、温度上昇に伴うorbital hoppingやdephasingが、光学coherenceとRaman経路を強く抑制する。

理論上は、温度依存rateを強散逸scale \(\Gamma(T)\) と見なし、probe sourceからEIT readoutまでのprojected Markov parametersを調べる。

### 13.1 横磁場の役割

横磁場 \(B_\perp\) は、対称性によって弱かったsector間経路を開く役割を持つ。

ただし、full Liouvillian計算ではfield exponentは温度とfit windowに依存する。したがって

> 「横磁場がperturbativeに経路を開く」

までは安全だが、全温度領域で一律に \(B_\perp^2\) と断言してはいけない。

### 13.2 reduced kernelとfull Liouvillian

reduced kernelは、どの経路が応答を担い、なぜmomentが消えるかを理解するために有用である。

一方、定量的なcontrast、符号反転、温度境界はfull Liouvillianを基準にすべきである。reduced modelは一部領域で大きさや符号を外すことがある。

### 13.3 室温no-goの意味

NV室温EITについて、本理論が直ちに示すのは必ずしもexact zeroではない。

より現実的な主張は、

- 応答が高い散逸次数で抑制される
- 漸近領域へ入る前後で信号が検出floorを下回る
- engineerable、asymptotic、detectableの共通窓が消える

というexperimental no-goである。

---

## 14. PRLの中心主張を学生向けに言い換える

この理論の中心を、数式を減らして言うと次のようになる。

> 強い散逸を受ける有限次元系では、入力から観測量へ至る最短のsector経路が、応答の最小抑制次数を決める。低次の経路が選択則で禁止されるほど、応答はより高い整数べきで消える。ただし、最短経路が複数ある場合は量子干渉によってその次数がさらに上がり得る。

数式では

\[
M_0=\cdots=M_{m-1}=0,
\qquad
M_m\neq0
\]

なら

\[
K_\Gamma\sim M_m\Gamma^{-m-1}.
\]

sector graphは \(m\ge d\) を保証し、\(M_d\neq0\) なら \(m=d\) となる。

---

## 15. 標準的な道具と、この理論での組み合わせ

以下の数学そのものは標準的である。

- rank-nullity
- Lindblad方程式
- Schur complement
- resolvent
- Neumann展開
- Laurent展開
- Krylov空間
- symmetry sector
- Drazin inverse

この理論での特徴は、それらを次の一本の流れに接続した点にある。

```text
symmetry / sector structure
          ↓
low-order projected Markov parameters vanish
          ↓
first nonzero coefficient selects dissipation order
          ↓
Schur complement transfers kernel orders to observable order
          ↓
finite-SNR experiment determines practical no-go
```

したがって、個々の数学用語を新概念として主張するのではなく、**既存の道具を開放量子系の強散逸応答へ統合した構成**として理解するのがよい。

---

## 16. よくある誤解

### 誤解1：\(\chi_{\rm full}=0\) ならEIT経路は存在しない

逆である。完全なEITによって全吸収がゼロになっている可能性がある。経路の有無は \(\delta\chi_S\) で判定する。

### 誤解2：graph distanceが2なら必ず \(\Gamma^{-3}\)

下限としては正しいが、最短経路間の干渉で \(M_2=0\) なら、さらに高次になる。

### 誤解3：Neumann条件を満たさないとLaurent漸近則は使えない

有限次元の有理resolventには無限遠Laurent展開が存在する。Neumann条件は便利な十分条件と誤差評価である。

### 誤解4：kernelの指数がそのまま実験contrastの指数になる

Schur complementの分母やnormalizationが指数を変える。observable全体を展開する必要がある。

### 誤解5：asymptotic no-goは有限 \(\Gamma\) でzeroという意味

有限 \(\Gamma\) では一般に非zeroである。強散逸極限で消えるという意味である。

### 誤解6：optical dark vectorがあれば定常dark stateも存在する

jump operatorの固有vector条件が別途必要である。

---

## 17. 演習問題

### 演習1：resolventの因数分解

\[
A_\Gamma=\Gamma D+A_0,
\qquad X=-D^{-1}A_0
\]

から

\[
A_\Gamma=D(\Gamma I-X)
\]

を示せ。また、\(\nu=D^{-1}c\) と置いて \(K_\Gamma\) をresolvent表示せよ。

<details>
<summary>解答の骨子</summary>

\(D(\Gamma I-X)=\Gamma D-DX\) に \(X=-D^{-1}A_0\) を代入する。

</details>

---

### 演習2：三sector chain

第7.1節の行列について、\(X^2\) と \(X^3\) を計算し、\(K_\Gamma\) を厳密に求めよ。

<details>
<summary>解答の骨子</summary>

\(X^3=0\) なのでNeumann級数は有限で止まる。

\[
(\Gamma I-X)^{-1}
=\Gamma^{-1}I+\Gamma^{-2}X+\Gamma^{-3}X^2.
\]

readout成分は \(g_1g_2\Gamma^{-3}\) のみ残る。

</details>

---

### 演習3：最短経路の相殺

二本の長さ2経路を持つgraphを行列で表し、\(M_2=A_a+A_b\) を示せ。\(A_a=-A_b\) のとき、次の非zero momentが何に依存するか考察せよ。

<details>
<summary>ヒント</summary>

sector内部の \(X_0\)、長さ3以上のwalk、または非対称な追加結合が次の係数へ寄与する。

</details>

---

### 演習4：Schur complement

block方程式から

\[
\delta\chi_S
=\chi_0d^\dagger A^{-1}BS_g^{-1}CA^{-1}b_p
\]

を導け。導出の途中で \(G_g^{-1}\) を仮定する必要がないことを確認せよ。

---

### 演習5：observable exponent

\[
L\sim\Gamma^{-2},
\qquad
R\sim\Gamma^{-3},
\qquad
S_g\sim\Gamma^{-1}
\]

のとき \(\delta\chi\) の指数を求めよ。さらに \(\chi_{\rm ref}\sim\Gamma^{-2}\) で正規化したcontrastの指数を求めよ。

<details>
<summary>解答</summary>

\[
\nu_{\rm obs}=2+3-1=4.
\]

正規化contrastは \(4-2=2\) 次である。

</details>

---

### 演習6：optical darkとstationary dark

二準位ground manifoldで

\[
|D\rangle=\frac{|g_1\rangle-|g_2\rangle}{\sqrt2}
\]

を考える。dephasing jump

\[
L=\sqrt\gamma|g_1\rangle\langle g_1|
\]

に対して、\(|D\rangle\) が \(L\) の固有vectorではないことを示せ。

---

## 18. 用語集

| 用語 | このノートでの意味 |
|---|---|
| source \(c\) | 外部probeが線形応答空間へ注入するvector |
| readout \(p\) | 最終的に測定する線形functional |
| resolvent | \((\Gamma I-X)^{-1}\) |
| projected Markov parameter | \(M_n=p^\dagger X^n\nu\) |
| invariant sector | generatorで保存される部分空間 |
| sector graph | sector間結合を有向辺で表したgraph |
| graph distance \(d\) | source sectorからreadout sectorまでの最短有向path長 |
| pathway cut | 注目sectorへの結合だけを切る反実仮想 |
| Schur complement \(S_g\) | 消去した自由度の影響を含む有効generator |
| Krylov space | \(\operatorname{span}\{\nu,X\nu,X^2\nu,\ldots\}\) |
| exact no-go | transferがregular domainで恒等的にzero |
| asymptotic no-go | \(\Gamma\to\infty\) で応答がべき的にzeroへ近づく |
| experimental no-go | 実現可能・漸近・検出可能領域の共通部分がない |

---

## 19. 次に読むべきリポジトリ内資料

1. `../Theorem and proofs/eit_nogo_proofs_v7.tex`  
   本ノートに対応する厳密な定理と証明。

2. `../Theorem and proofs/eit_nogo_proofs_v6_2_archive.tex`  
   修正前の証明パッケージ。Neumann展開中心の旧構成との比較用。

3. `../results/gate_1_5_report.md`  
   NV候補点におけるEIT/ATS、full-vs-reduced、SNR、温度閾値、ensemble平均の数値監査。

4. `../../Writing Paper/NV_EIT_PRA_PRL_Split_Strategy_20260724.md`  
   NV固有のPRAと一般則を扱うPRLの分割方針。

---

## 20. このドラフトで今後補うべき点

- NVの具体的な準位図と、各sectorがどのdensity-matrix成分に対応するか
- 三sector chainを実際の数値コードと照合する短いnotebook
- Schur complementから通常の三準位EIT感受率を得る詳細計算
- full Liouvillianとreduced kernelが一致する領域・外れる領域の図解
- experimental no-goをSNRと測定時間へ変換する例
- 先行研究の標準語彙との対応表

---

## まとめ

この理論の最小骨格は、次の四行に凝縮できる。

\[
K_\Gamma
=p^\dagger(\Gamma I-X)^{-1}\nu,
\]

\[
K_\Gamma
=\sum_{n\ge0}M_n\Gamma^{-n-1},
\qquad
M_n=p^\dagger X^n\nu,
\]

\[
M_0=\cdots=M_{m-1}=0,
\quad M_m\neq0
\Longrightarrow
K_\Gamma\sim M_m\Gamma^{-m-1},
\]

\[
\operatorname{dist}_G(s,t)=d
\Longrightarrow
m\ge d.
\]

graphは到達に必要な最小sector遷移数を与え、量子振幅はその最短経路が生き残るかを決める。そして、最終的に実験で測る指数は、個々のkernelではなくSchur complementを含むobservable全体から決まる。
