# 新 no-go 理論を用いた簡易数値計算：優先解析事項

## 0. 目的と到達点

本計算の第一目的は、新 no-go 理論の三分類

\[
\nu=\infty,\qquad 0<\nu<\infty,\qquad \nu=0
\]

を数値的に再現することではない。それだけでは、定理の例示に留まる。

優先すべき目的は、次の仮説を最小模型で検証することである。

> 全光学応答が強散逸下でも有限に残る一方、EITを担う sector-resolved response のみが高次に抑制される場合がある。この差を利用すれば、見かけ上類似するEIT、Autler–Townes splitting（ATS）、および背景透過を定量的に分離できる可能性がある。

ただし、現段階では「新 no-go 理論単独でEITとATSを完全に識別する条件式が既に証明された」とは主張しない。今回の簡易計算では、master response の散逸指数と、従来のスペクトル極・線幅・制御強度依存を接続できるかを確認する。

---

## 1. 最優先で計算する量

指定した長寿命コヒーレンス sector \(S\) に対して、

\[
R_{S,\Gamma}(z)
=
\chi_{\mathrm{full},\Gamma}(z)
-
\chi^{(S)}_{\mathrm{cut},\Gamma}(z)
\]

を直接計算する。

同時に、以下の三量を必ず保存する。

\[
\chi_{\mathrm{full},\Gamma}(z),\qquad
\chi^{(S)}_{\mathrm{cut},\Gamma}(z),\qquad
R_{S,\Gamma}(z).
\]

全応答 \(\chi_{\mathrm{full}}\) だけを計算してはならない。全応答の透明窓が残っていても、sector-mediated coherent response が残っているとは限らないためである。

### sector cut の規約

初期計算では **frozen-source cut** を用いる。

- full model と cut model で散逸行列 \(D\) を共通にする。
- source \(c\)、readout \(p\)、定常状態由来の線形化係数を固定する。
- 指定 sector と光学応答を結ぶ結合だけを切る。
- control field を単純にOFFにしたスペクトルを cut response とみなさない。

これにより、定理で分類する master response と数値実装を一致させる。

---

## 2. 最小模型

### 2.1 第一段階：標準的な \(\Lambda\) 模型

線形応答行列を

\[
M(\Delta,\delta)=
\begin{pmatrix}
\gamma_{31}-i\Delta & i\Omega_c/2\\
i\Omega_c^*/2 & \gamma_{21}-i\delta
\end{pmatrix}
\]

とし、probe optical coherence を観測する。

full response は

\[
\chi_{\mathrm{full}}
\propto
\frac{\gamma_{21}-i\delta}
{(\gamma_{31}-i\Delta)(\gamma_{21}-i\delta)+|\Omega_c|^2/4},
\]

sector-cut response は、ground coherence 経路を切った

\[
\chi_{\mathrm{cut}}
\propto
\frac{1}{\gamma_{31}-i\Delta}
\]

とする。

この模型では、既知のEIT–ATS crossoverを再現しつつ、

\[
R_S=\chi_{\mathrm{full}}-\chi_{\mathrm{cut}}
\]

がどの領域で大きいか、または小さいかを調べる。

### 2.2 第二段階：保護背景を持つ拡張模型

新理論固有の効果を出すため、少なくとも一つの散逸kernelを持つ模型へ拡張する。

\[
A_\Gamma(z)=\Gamma D+B(z),\qquad \ker D\neq\{0\}.
\]

必要条件は、

\[
\chi_{\mathrm{full},\Gamma}=a_0+\frac{a_1}{\Gamma}+\cdots,
\]

\[
\chi_{\mathrm{cut},\Gamma}=a_0+\frac{b_1}{\Gamma}+\cdots
\]

となり、共通の保護背景 \(a_0\) が差分で相殺されることである。すると、

\[
R_{S,\Gamma}
=
\frac{a_1-b_1}{\Gamma}+O(\Gamma^{-2}),
\]

となり、

\[
\nu[\chi_{\mathrm{full}}]=0,\qquad
\nu[\chi_{\mathrm{cut}}]=0,\qquad
\nu[R_S]=1
\]

が実現する。

これを **protected-background cancellation** と呼ぶ。

---

## 3. 解析課題の優先順位

## Priority 1：EIT、ATS、背景応答の数値的分離

### 目的

透明窓が現れた場合に、その起源が

1. 長寿命ground coherenceを介したEIT
2. 強制御によるATS
3. full/cutに共通する保護背景または直接透過

のどれであるかを判定する。

### 掃引パラメータ

\[
\Omega_c/\gamma_{31},\qquad
\gamma_{21}/\gamma_{31},\qquad
\Gamma/\gamma_{31},\qquad
\Delta,\delta
\]

を掃引する。

### 保存する診断量

- \(\operatorname{Im}\chi_{\mathrm{full}}\)
- \(\operatorname{Im}\chi_{\mathrm{cut}}\)
- \(\operatorname{Im}R_S\)
- transparency depth
- peak separation
- two-photon linewidth
- pole positions and residues
- \(\partial \text{linewidth}/\partial |\Omega_c|^2\)
- \(R_S\) の \(\Gamma\) 依存性

### 暫定的な識別規則

以下は最終定理ではなく、数値検証すべき候補である。

#### EIT-dominant

\[
|R_S(\Delta\simeq0)| \text{ が透明窓形成に主要な寄与を持つ},
\]

かつ、

- transparency linewidth が \(\gamma_{21}\) に敏感
- ground coherence を切ると窓が消失または大幅に減少
- \(R_S\) が無視できない
- 極の干渉またはresidue cancellationが確認できる

#### ATS-dominant

- \(\chi_{\mathrm{cut}}\) だけでも二峰構造または分裂が残る
- peak separation が概ね \(|\Omega_c|\) に比例する
- ground coherence sector を切っても中心の構造が大きく変わらない
- \(|R_S|/|\chi_{\mathrm{full}}|\) が小さい

#### Protected-background-dominant

- \(\chi_{\mathrm{full}}\) と \(\chi_{\mathrm{cut}}\) がともに \(O(1)\)
- 差分 \(R_S\) のみが \(O(\Gamma^{-\nu})\), \(\nu>0\)
- 全スペクトルでは透明または低吸収が残るが、coherent sector contribution は消失する

### 合格条件

少なくとも三つの典型点を一枚のphase map上で分離できること。

\[
\text{EIT-dominant},\qquad
\text{ATS-dominant},\qquad
\text{background-dominant}.
\]

---

## Priority 2：master-response指数 \(\nu\) の抽出

各パラメータ点で、

\[
\|R_{S,\Gamma}\|_K
=
\max_{z\in K}|R_{S,\Gamma}(z)|
\]

を計算し、

\[
\|R_{S,\Gamma}\|_K\propto \Gamma^{-\nu}
\]

をlog–log fitする。

### 数値推定

隣接点から局所有効指数

\[
\nu_{\mathrm{eff}}(\Gamma)
=
-\frac{
\log \|R_{S,\Gamma_2}\|_K-
\log \|R_{S,\Gamma_1}\|_K
}{
\log \Gamma_2-\log \Gamma_1
}
\]

を求める。

### 同時に行う厳密確認

可能な小行列模型では、数値fitだけでなく、

- master moments
- protected leading coefficient \(R_{S,0}\)
- adjugate/determinantの次数差

から理論値 \(\nu\) を求め、直接fitと照合する。

### 合格条件

- Class III：\(\nu_{\mathrm{eff}}\to0\)
- Class II：\(\nu_{\mathrm{eff}}\to m>0\)
- Class I：symbolicに \(R_S\equiv0\)

を再現する。

浮動小数点で小さいことをClass Iの証明に使ってはならない。

---

## Priority 3：隠れたclass transition

外場、歪み、偏光角、detuningなどを制御パラメータ \(\lambda\) とする。

特異散逸模型では、

\[
r_0(\lambda)
=
p^\dagger P
\left[
B_{P,\mathrm{full}}(\lambda)^{-1}
-
B_{P,\mathrm{cut}}(\lambda)^{-1}
\right]
Pc
\]

を計算する。

### 探索条件

\[
r_0(\lambda_c)=0
\]

を満たす \(\lambda_c\) を探索する。

そのとき、

\[
\lambda\neq\lambda_c:\quad
R_{S,\Gamma}=O(1),\quad \nu=0,
\]

\[
\lambda=\lambda_c:\quad
R_{S,\Gamma}=O(\Gamma^{-1}),\quad \nu=1
\]

となるか確認する。

### 最重要の新規予測

\[
\chi_{\mathrm{full}}=O(1)
\]

のままでも、

\[
\nu[R_S]:0\rightarrow1
\]

の転移が起こり得る。

これは通常の全吸収スペクトルだけでは見えず、full–cut responseを評価した場合にのみ検出される。

### 合格条件

- \(\chi_{\mathrm{full}}\) のclassは変わらない
- \(R_S\) の指数だけが変化する
- 変化点が \(r_0(\lambda_c)=0\) と一致する

---

## Priority 4：境界近傍のscaling collapse

\(r_0(\lambda)\) が単純零を持つ場合、

\[
r_0(\lambda)\simeq
a(\lambda-\lambda_c)
\]

より、

\[
R_{S,\Gamma}
\simeq
a(\lambda-\lambda_c)
+
\frac{r_1}{\Gamma}.
\]

したがって、

\[
\Gamma R_{S,\Gamma}
\simeq
a\,\Gamma(\lambda-\lambda_c)+r_1.
\]

### 作成する図

横軸：

\[
X=\Gamma(\lambda-\lambda_c)
\]

縦軸：

\[
Y=\Gamma R_{S,\Gamma}
\]

異なる \(\Gamma\) の曲線が一つのuniversal curveにcollapseするかを確認する。

### 期待される予測

- crossover width：

\[
|\lambda-\lambda_c|\sim\Gamma^{-1}
\]

- 境界上のamplitude：

\[
|R_S|\sim\Gamma^{-1}
\]

- 観測contrastが \(|R_S|^2\) に比例する場合：

\[
C_{\min}\sim\Gamma^{-2}
\]

### 合格条件

少なくとも1 decade以上の \(\Gamma\) 範囲でcollapseが改善すること。

---

## Priority 5：EIT–ATS判定量の候補構築

新理論をEIT–ATS識別に接続するため、次の無次元量を試す。

### sector fraction

\[
\eta_S(z)
=
\frac{|R_S(z)|}
{|\chi_{\mathrm{full}}(z)|+\epsilon}.
\]

透明窓近傍で \(\eta_S\) が大きければsector-mediated mechanismが支配的であり、小さければATSまたは背景応答が支配的である可能性が高い。

### integrated sector weight

\[
W_S
=
\frac{
\int_K |R_S(z)|^2\,dz
}{
\int_K |\chi_{\mathrm{full}}(z)|^2\,dz+\epsilon
}.
\]

### dissipation-resolved sector index

\[
\mathcal I_S
=
\left(
W_S,\,
\nu[R_S],\,
\frac{\partial \Delta_{\mathrm{peak}}}{\partial|\Omega_c|},\,
\frac{\partial w_{\mathrm{dip}}}{\partial\gamma_{21}}
\right).
\]

単一の閾値を最初から仮定せず、この多変量指標がEIT、ATS、背景透過を分離できるか確認する。

### 注意

\(\eta_S\) や \(W_S\) の数値閾値は普遍定数とは限らない。最初の目的は、既知のEIT領域とATS領域で明確に異なるclusterを形成するかを確認することである。

---

## 4. 最初に作成する図

### Figure 1：理論のunit test

\[
|R_{S,\Gamma}|
\quad \text{vs}\quad \Gamma
\]

をlog–log表示し、

- Class I
- Class II, \(\nu=1\)
- Class II, \(\nu=2\)
- Class III

を同一図または個別図で示す。

### Figure 2：EIT–ATS crossover

横軸を \(|\Omega_c|/\gamma_{31}\)、縦軸をdetuningとし、

- \(\operatorname{Im}\chi_{\mathrm{full}}\)
- \(\operatorname{Im}\chi_{\mathrm{cut}}\)
- \(\operatorname{Im}R_S\)

の三つのmapを並べる。

### Figure 3：hidden class transition

\[
\nu_{\mathrm{eff}}(\lambda),\qquad
|r_0(\lambda)|,\qquad
|\chi_{\mathrm{full}}|,\qquad
|R_S|
\]

を同じ \(\lambda\) 掃引で示す。

### Figure 4：scaling collapse

\[
\Gamma R_{S,\Gamma}
\quad \text{vs}\quad
\Gamma(\lambda-\lambda_c)
\]

を複数の \(\Gamma\) で重ねる。

### Figure 5：mechanism map

\[
\left(
|\Omega_c|/\gamma_{31},
\gamma_{21}/\gamma_{31}
\right)
\]

または

\[
\left(
|\Omega_c|/\gamma_{31},
\Gamma/\gamma_{31}
\right)
\]

平面上で、

- EIT-dominant
- ATS-dominant
- protected-background-dominant
- unresolved

を分類する。

---

## 5. 実装順序

### Phase A：解析式と小行列

1. 標準 \(\Lambda\) 模型を実装する。
2. full、cut、master responseを解析式と数値逆行列で照合する。
3. pole、residue、linewidth、peak separationを抽出する。
4. 既知のEIT–ATS crossoverを再現する。

### Phase B：新理論固有の模型

1. singular \(D\) を持つ3～5次元模型を構成する。
2. Riesz projection \(P\) を計算する。
3. \(R_{S,0}\) を直接評価する。
4. protected-background cancellationを実現する。
5. \(\nu=0\to1\) のclass transitionを探索する。
6. scaling collapseを確認する。

### Phase C：物理模型への移植

簡易模型で効果が確認できた後に、

1. Rbまたは標準 \(\Lambda\) 系
2. NVのzero-field spin-\(\Lambda\)
3. transverse-field escape channel
4. SiV/SnV orbital-\(\Lambda\)

へ適用する。

材料模型から始めない。最初から物性パラメータ、選択則、ensemble averagingを全部入れると、理論固有の効果とモデル依存効果を区別できなくなるためである。

---

## 6. 成功判定ゲート

### Gate 1：数値分類器の正しさ

moment法、protected coefficient法、直接fitが同じ \(\nu\) を返す。

### Gate 2：新理論固有の現象

\[
\nu[\chi_{\mathrm{full}}]=0
\]

のまま、

\[
\nu[R_S]:0\rightarrow1
\]

となるパラメータ点を一つ以上発見する。

### Gate 3：EIT–ATS接続

既知のEIT領域では \(R_S\) が透明窓に主要寄与し、ATS領域では \(\chi_{\mathrm{cut}}\) が主要構造を保持することを示す。

### Gate 4：定量的予測

少なくとも一つを確認する。

\[
|\lambda-\lambda_c|\sim\Gamma^{-1},
\]

\[
|R_S(\lambda_c)|\sim\Gamma^{-1},
\]

\[
C_{\min}\sim\Gamma^{-2},
\]

または

\[
\Gamma R_S
=
\mathcal F\!\left(\Gamma(\lambda-\lambda_c)\right).
\]

### Gate 5：論文化判断

以下をすべて満たせば、新理論を独立した主役として論文化する価値が高い。

- hidden class transitionが存在する
- total spectrumだけでは転移が判定できない
- master responseでのみ明確に検出できる
- EIT/ATS/背景応答の識別改善につながる
- 実材料または実験可能な有効模型でも再現する

---

## 7. 当面の結論

最優先課題は、材料別の大規模計算ではなく、

\[
\boxed{
\text{全応答はprotectedだが、
EIT-sector responseのみsuppressedになる模型}
}
\]

を構成し、

\[
\boxed{
\nu[\chi_{\mathrm{full}}]=0,\qquad
\nu[R_S]>0
}
\]

を数値的に示すことである。

次に、その差がEIT–ATS識別に使えるかを、標準 \(\Lambda\) 模型のpole、residue、linewidth、control-power scalingと照合する。

この二段階が成功すれば、新 no-go 理論は単なる分類の再整理ではなく、

> 見かけ上同じ透明窓の内部にあるcoherent sectorの生存・消失を識別する理論

として位置づけられる。
