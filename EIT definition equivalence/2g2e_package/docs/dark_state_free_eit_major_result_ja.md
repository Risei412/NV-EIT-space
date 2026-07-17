# Dark-State-Free Coherent Transparency in a Minimal \(2g+2e\) Markovian Model

## 最小模型におけるEIT型干渉、完全零点のno-go境界、および新no-go理論との統合

## 1. 研究上の結論

本解析では、2つの下準位と2つの励起準位からなる最小のMarkovian Lindblad模型を用いて、以下を同時に示した。

1. 結合行列を
   \[
   \operatorname{rank}\Omega=2
   \]
   に固定することで、光学dark stateを完全に排除できる。

2. 全radiative jump operatorおよび有効Hamiltonian条件を用いることで、probeとcontrolがともに非零である限り、stationary pure Lindblad dark stateも存在しないことを証明できる。

3. dark stateが存在しないにもかかわらず、sector-mediated interferenceにより、sector-cut応答と比較して約79%の吸収抑制が生じる。

4. pole/residue解析により、この吸収窓はAutler–Townes splittingではなく、反対符号のpole寄与によるEIT/Fano型干渉であると判定できる。

5. 一方で、正則な理想 \(2g+2e\) 模型では、
   \[
   \operatorname{rank}\Omega=2
   \]
   のまま完全な実軸transfer zeroを実現できないことを、\(2\times2\) 行列式恒等式から証明できる。

したがって最小模型の結論は、

\[
\boxed{
\text{dark stateなしでも強いEIT型coherent transparencyは存在し得る}
}
\]

である一方、

\[
\boxed{
\text{正則な最小 }2g+2e\text{ 模型では完全なdark-state-free EIT zeroは不可能}
}
\]

である。

この二つの結果は矛盾しない。複素transfer zeroが実軸へ近づくことで深い吸収極小は形成できるが、rank-two条件の下では零点を実軸上へ完全に置くことができない。

---

## 2. 最小 \(2g+2e\) 模型

基底を

\[
\{|g_1\rangle,|g_2\rangle,|e_1\rangle,|e_2\rangle\}
\]

とする。

励起空間におけるprobeおよびcontrolの双極子ベクトルを

\[
d_p=
\begin{pmatrix}
1\\
0
\end{pmatrix},
\qquad
d_c=
\begin{pmatrix}
\cos\theta\\
\sin\theta
\end{pmatrix}
\]

とする。

結合行列は

\[
\Omega=
\begin{pmatrix}
\Omega_p & \Omega_c\cos\theta\\
0 & \Omega_c\sin\theta
\end{pmatrix}.
\]

\(\Omega_p\Omega_c\sin\theta\neq0\)ならば、

\[
\det\Omega=\Omega_p\Omega_c\sin\theta\neq0,
\]

したがって、

\[
\operatorname{rank}\Omega=2,
\qquad
\ker\Omega=\{0\}.
\]

これは、どの下準位重ね合わせも光励起から完全には切り離されないことを意味する。代表計算では \(\theta=\pi/4\) を用いた。

---

## 3. stationary pure Lindblad dark stateの非存在証明

radiative jump operatorを

\[
L_{j\alpha}
=
\sqrt{\gamma_{j\alpha}}
|g_\alpha\rangle\langle e_j|,
\qquad
j,\alpha\in\{1,2\}
\]

とする。

各jump operatorはnilpotentであり、固有値は0のみである。stationary pure state \(|D\rangle\) が存在するためには、全jump operatorの共通固有ベクトルでなければならないため、

\[
L_{j\alpha}|D\rangle=0
\]

をすべての \(j,\alpha\) について満たす必要がある。

すべてのradiative jump operatorの共通kernelはground subspaceであるから、

\[
|D\rangle
\in
\operatorname{span}\{|g_1\rangle,|g_2\rangle\}.
\]

さらにpure stationary Lindblad stateの有効Hamiltonian条件をexcited subspaceへ射影すると、

\[
P_eH_{\mathrm{eff}}|D\rangle
=
\frac{1}{2}\Omega|D\rangle
=
0
\]

が必要である。

しかし、\(\operatorname{rank}\Omega=2\) なので \(\ker\Omega=\{0\}\) である。したがって非零の \(|D\rangle\) は存在せず、stationary pure Lindblad dark stateは存在しない。

### weak-probe極限に関する注意

厳密に \(\Omega_p=0\) と設定したreference Liouvillianでは、\(|g_1\rangle\) がcontrolに結合しない自明な状態となり得る。

したがってdark-state-freeという主張は、

1. 任意の有限な \(\Omega_p\neq0\) でrank-two模型を定義し、
2. その後にweak-probe極限を取る、

という順序に対して成立する。このorder-of-limitsは論文化の際に明示する必要がある。

---

## 4. sector-resolved response

次を定義する。

\[
a_1=\gamma_e+\Gamma_1-i\delta,
\]

\[
a_2=\gamma_e+\Gamma_2+i\Delta_e-i\delta,
\]

\[
g=\gamma_g+\Gamma_g-i\delta,
\]

\[
\beta=\frac{|\Omega_c|^2}{4},
\qquad
c=\cos\theta,
\qquad
s=\sin\theta.
\]

ground-coherence sector \(S\) を含むfull responseは、

\[
\boxed{
\chi_{\mathrm{full}}
=
\frac{ga_2+\beta s^2}
{ga_1a_2+\beta(c^2a_2+s^2a_1)}
}
\]

である。

sector \(S\) を切断したcounterfactual responseは、

\[
\boxed{
\chi_{\mathrm{cut}}^{(S)}
=
\frac{1}{a_1}
}
\]

である。

したがって、新no-go理論が分類するmaster sector responseは、

\[
R_S
=
\chi_{\mathrm{full}}
-
\chi_{\mathrm{cut}}^{(S)}
\]

であり、明示的には

\[
\boxed{
R_S
=
-\frac{\beta c^2a_2}
{a_1\left[
ga_1a_2+\beta(c^2a_2+s^2a_1)
\right]}
}
\]

となる。

重要なのは、EITの存在をtotal susceptibilityの零点だけで議論せず、

\[
\chi_{\mathrm{full}}
=
\chi_{\mathrm{cut}}^{(S)}
+
R_S
\]

という背景応答とsector-mediated correctionの干渉として解析することである。

---

## 5. 最小模型における完全零点no-go定理

励起空間のresolventを \(G=A^{-1}\) とし、

\[
D=(d_p,d_c)
\]

を2本の双極子ベクトルからなる行列とする。

次を定義する。

\[
S_p=d_p^\dagger Gd_p,
\qquad
S_c=d_c^\dagger Gd_c,
\]

\[
K_{pc}=d_p^\dagger Gd_c,
\qquad
K_{cp}=d_c^\dagger Gd_p.
\]

このとき、

\[
D^\dagger GD
=
\begin{pmatrix}
S_p & K_{pc}\\
K_{cp} & S_c
\end{pmatrix}.
\]

よって、

\[
S_pS_c-K_{pc}K_{cp}
=
\det(D^\dagger GD).
\]

行列式の乗法性から、

\[
\det(D^\dagger GD)
=
\det(D^\dagger)\det(G)\det(D)
=
|\det D|^2\det G.
\]

したがって、

\[
\boxed{
S_pS_c-K_{pc}K_{cp}
=
|\det D|^2\det G
}
\]

が成立する。

理想二光子共鳴 \(g=0\) では、

\[
\chi_{\mathrm{full}}
=
S_p-
\frac{K_{pc}K_{cp}}{S_c}
=
\frac{S_pS_c-K_{pc}K_{cp}}{S_c}.
\]

したがって、

\[
\chi_{\mathrm{full}}
=
\frac{|\det D|^2\det G}{S_c}.
\]

\(A\)が正則であり、かつ \(\operatorname{rank}D=2\) ならば \(\det D\neq0\) である。ゆえに、

\[
\boxed{
\chi_{\mathrm{full}}\neq0.
}
\]

### 定理の意味

正則な理想 \(2g+2e\) 模型では、完全な実軸transfer zeroを実現するためには

\[
\det D=0
\]

が必要である。これは \(\operatorname{rank}D<2\)、すなわち光学dark stateの存在条件と一致する。

この結果は、excited-state mixingや非Hermitian・non-normalな \(G\) を含めても成立する。証明では \(G\) を任意の正則 \(2\times2\) 行列として扱っているためである。

---

## 6. 数値結果

代表パラメータを、\(\gamma_e\) を単位として、

\[
\gamma_e=1,
\qquad
\gamma_g=0.02,
\]

\[
\Delta_e=8,
\qquad
\Omega_c=0.8,
\]

\[
\theta=\frac{\pi}{4}
\]

とした。

吸収極小は、

\[
\delta_{\min}/\gamma_e
\simeq
-0.00965
\]

で生じた。

そのとき、

\[
\operatorname{Re}\chi_{\mathrm{full}}
\simeq
0.20970,
\]

\[
\operatorname{Re}\chi_{\mathrm{cut}}^{(S)}
\simeq
0.99991.
\]

したがってsector cutに対する吸収抑制率は、

\[
1-
\frac{\operatorname{Re}\chi_{\mathrm{full}}}
{\operatorname{Re}\chi_{\mathrm{cut}}^{(S)}}
\simeq
0.7903.
\]

すなわち、

\[
\boxed{
\text{約 }79.0\%\text{ の吸収抑制}
}
\]

が得られた。

最も近い複素transfer zeroは、

\[
z_0/\gamma_e
\simeq
-0.009841
-
0.021203i
\]

に存在する。

零点は実軸上にはないが、虚部が小さいため、実周波数スペクトルに深い吸収極小を形成する。

---

## 7. pole/residue解析とATS排除

近共鳴の二つのpoleは、

\[
z_1/\gamma_e
\simeq
-0.010976
-
0.111078i,
\]

\[
z_2/\gamma_e
\simeq
0.001123
-
0.910128i.
\]

実部間隔は、

\[
|\operatorname{Re}z_1-\operatorname{Re}z_2|
\simeq
0.01210\,\gamma_e.
\]

一方、pole linewidthに対応する虚部の絶対値は、

\[
|\operatorname{Im}z_1|
\simeq
0.1111\,\gamma_e,
\]

\[
|\operatorname{Im}z_2|
\simeq
0.9101\,\gamma_e.
\]

したがって二つのpoleは、周波数方向に分離したAutler–Townes doubletを形成していない。

対応する留数は、

\[
r_1
\simeq
-0.003123
-
0.112430i,
\]

\[
r_2
\simeq
0.003123
+
1.112428i.
\]

吸収極小における各poleの実部寄与は、

\[
+1.22205,
\qquad
-1.01234.
\]

反対符号の寄与が相殺して吸収窓を作るため、この現象は独立した二つの吸収peakによるATSではなく、Fano/EIT型の干渉である。

---

## 8. \(\Gamma\) scalingと新no-go理論のClass I–III

新no-go理論では、分類は単一の \(\Gamma\) の値ではなく、指定したscaling path

\[
A_\Gamma=\Gamma D+B
\]

に沿う \(\Gamma\to\infty\) の漸近挙動として定義される。

### 8.1 Class I：exact structural zero

\[
\theta=\frac{\pi}{2}
\]

とすると、

\[
d_p=
\begin{pmatrix}
1\\0
\end{pmatrix},
\qquad
d_c=
\begin{pmatrix}
0\\1
\end{pmatrix}.
\]

この場合も結合行列はrank twoであるが、\(c=\cos\theta=0\) なので、

\[
R_S\equiv0.
\]

したがって、

\[
\boxed{
\nu=\infty,
\qquad
\mathrm{Class\ I}.
}
\]

これは光学dark stateとは異なる、source-to-readout sector transferのexact structural zeroである。

### 8.2 Class II：両optical coherenceの散逸

\[
\Gamma_1=\Gamma_2=\Gamma,
\qquad
\Gamma_g=0
\]

とする。

固定された \(g\neq0\) に対し、

\[
R_S
=
-\frac{\beta c^2}{g\Gamma^2}
+
O(\Gamma^{-3}).
\]

数値fitでは、

\[
|R_S|\propto\Gamma^{-1.9993}.
\]

したがって、

\[
\boxed{
\nu=2,
\qquad
\mathrm{Class\ II}.
}
\]

さらにground coherenceも同じscaleで散逸させると、

\[
R_S=O(\Gamma^{-3}),
\]

数値fitはほぼ \(\Gamma^{-2.9999}\) となる。

### 8.3 Class III：protected survival

\[
\Gamma_2=\Gamma,
\qquad
\Gamma_1=\Gamma_g=0
\]

とし、直交した \(e_2\) optical channelのみを強く散逸させる。

このとき、

\[
R_S
\longrightarrow
-\frac{\beta c^2}
{a_1(ga_1+\beta c^2)}
\neq0.
\]

数値fitの傾きはほぼ0である。

したがって、

\[
\boxed{
\nu=0,
\qquad
\mathrm{Class\ III}.
}
\]

source、readoutおよび共通Raman経路がscaled dissipatorのkernel内に残るため、order-one responseが保護される。

---

## 9. NV centerへの写像

局所的な対応を次のように置く。

\[
|g_1\rangle,|g_2\rangle
\leftrightarrow
\text{two NV ground-spin sublevels},
\]

\[
|e_1\rangle,|e_2\rangle
\leftrightarrow
\text{two selected excited orbital-spin branches}.
\]

\(d_p,d_c\) は、probe/controlの偏光、spin–orbit mixing、spin–spin interaction、strain、hyperfine coupling、磁場によって決まる励起空間内の双極子ベクトルである。

\(c=\cos\theta\) は、二つのoptical legが共有するexcited-manifold componentの大きさを表す。

### 9.1 symmetry-preserving limit

ゼロ磁場かつspin-preservingな簡約模型では、二つのspin legが直交し、\(c=0\) となり得る。

この場合、\(R_S=0\) であり、簡約模型ではClass Iに対応する。

### 9.2 transverse-field escape

横磁場などのsector-mixing perturbationにより、

\[
c\propto B_\perp
\]

が開くとする。

両optical coherenceがorbital hoppingで散逸するClass-II pathでは、

\[
R_S
\propto
c^2\Gamma_{\mathrm{oc}}^{-2}.
\]

したがって、

\[
\boxed{
R_S
\propto
B_\perp^2
\Gamma_{\mathrm{oc}}^{-2}.
}
\]

これは横磁場に対するquadratic escape lawである。

### 9.3 orbital hopping rate

NVの対称なorbital hoppingに対して、population imbalance rateを \(\Gamma_{XY}\) とすると、optical-coherence dampingは

\[
\Gamma_{\mathrm{oc}}
=
\frac{\Gamma_{XY}}{4}
\]

となる。

両optical coherenceが同時に散逸するため、標準的なNV mappingはClass IIである。

Class IIIをNVで実現するには、scaled dissipatorのkernel内に光学応答を運ぶ真のprotected coherence combinationが必要である。

### 9.4 六準位模型への拡張

実際のNV \(^{3}E\) excited manifoldは六準位である。

したがって今回の \(2e\) mappingは、

- 最小局所機構、
- scaling則、
- exact-zero no-go境界、

を示す模型であり、最終的なNV material-level proofではない。

次段階では六準位excited-manifold Liouvillianに対して、

1. rank audit、
2. pure dark state audit、
3. sector cut、
4. transfer-zero探索、
5. pole/residue解析、
6. \(\Gamma(T)\) scaling、

を同じ形式で実行する必要がある。

---

## 10. 理論的インパクト

この結果は、「dark stateなしでも透明窓が見える」という個別数値例だけではない。

最小次元における次の境界を与える。

### 存在側

\[
\operatorname{rank}\Omega=2
\]

かつpure Lindblad dark stateが存在しなくても、sector-mediated transferが背景吸収と干渉することで、大きな吸収抑制が可能である。

### 非存在側

正則な \(2g+2e\) excited resolventでは、完全な実軸transfer zeroは

\[
\det D=0
\]

を必要とする。

これは光学dark stateの存在条件と一致する。

### 次の理論課題

完全なdark-state-free EITを実現する最小拡張として、以下が候補となる。

1. \(2g+3e\) 模型
2. closed-loop coherent coupling
3. correlated or nontrivial Lindblad jumps
4. non-normal multichannel resolvent
5. NVの完全な六励起準位構造

とくに \(2g+3e\) 以上では、\(D\) が長方形行列となるため、

\[
\det(D^\dagger GD)
=
|\det D|^2\det G
\]

という \(2\times2\) 正方行列特有の障害はそのままでは成立しない。

したがって、完全なdark-state-free EITの最小次元が \(2g+3e\) である可能性がある。

---

## 11. 論文構成への提案

### 中心命題

> Optical dark states and input–output transparency zeros are generally distinct objects. A minimal \(2g+2e\) Markovian model supports strong dark-state-free coherent transparency, while an exact determinant identity forbids a perfect real-axis transparency zero at full optical rank.

### 推奨Figure

1. **Figure 1:** \(2g+2e\) 模型とsector cut
2. **Figure 2:** \(\chi_{\mathrm{full}}\) と \(\chi_{\mathrm{cut}}^{(S)}\) の吸収スペクトル
3. **Figure 3:** complex pole-zero map
4. **Figure 4:** 留数分解とATS排除
5. **Figure 5:** \(\Gamma\) scalingによるClass I–III
6. **Figure 6:** NVへの写像と \(B_\perp^2\Gamma_{\mathrm{oc}}^{-2}\) 予測

### 投稿上の位置づけ

この \(2g+2e\) 結果単独では、完全なdark-state-free EITの存在定理ではなく、以下の二重結果として提示するのが正確である。

- dark-state-free coherent-transparency existence
- minimal exact-zero no-go theorem

その後、\(2g+3e\) またはNV六準位模型で完全零点を実現できれば、

- 最小次元境界、
- 一般判定条件、
- 実材料への写像、

を組み合わせた強い理論論文へ発展し得る。

---

## 12. ZIP packageの構成

- `docs/`
  - 日本語研究まとめ
  - 英文解析レポート
- `code/`
  - 再現用Pythonコード
  - requirements
- `data/`
  - 吸収スペクトルCSV
  - \(\Gamma\) scaling CSV
- `figures/`
  - 吸収スペクトル図
  - pole-zero map
  - Class I–III scaling図
- `metadata/`
  - 使用パラメータJSON
  - SHA-256 checksum一覧

---

## 13. 最終結論

\[
\boxed{
\text{dark stateは内部状態のkernel条件であり、EIT型透明窓は入出力transferの干渉条件である}
}
\]

両者は同じではない。

しかし最小正則 \(2g+2e\) 模型では、完全な実軸零点に限れば両者は行列式恒等式によって再び結び付けられる。

したがって今回の結果は、

1. dark stateとcoherent transparencyの概念的分離、
2. 最小模型における強い吸収抑制の存在、
3. 完全零点に対する新しいno-go境界、
4. 新no-go理論Class I–IIIとの整合、
5. NV横磁場escapeへの具体的予測、

を同時に与える。

これは完全なdark-state-free EITを探索するための、明確な出発点と最小次元benchmarkである。
