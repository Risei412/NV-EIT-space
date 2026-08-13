# EIT 定義同値プログラム：次段階の数値計算計画

**対象:** `EIT definition equivalence/`  
**目的:** 一般化 no-go の物理的スコープを確定し、受動的 Markovian 系における dark-state-free 厳密 EIT zero の可否と、最小の正当な escape structure を数値的に決定する。  
**方針:** 最初は抽象的な \(2g+N_e e\) 最小模型で構造を確定し、最後にのみ NV 励起6準位模型へ写像する。

---

## 1. 現時点の到達点と修正すべき問題

### 1.1 確実に得られた結果

現在の計算では、対角励起応答

\[
G(\omega)=\mathrm{diag}(g_1,\ldots,g_{N_e}),
\qquad
\operatorname{Re}g_j>0
\]

と full-rank 光学結合

\[
D\in\mathbb C^{N_e\times 2},
\qquad
\operatorname{rank}D=2
\]

に対して、

\[
M(\omega)=D^\dagger G(\omega)D
\]

が実軸上で特異にならないことが示されている。したがって、励起状態数を増やすだけでは、dark state のない厳密な実軸 transfer zero は生成されない。

この結果は `results/E4_findings.md` における「任意次元・対角励起多様体の一般化 no-go」に対応する。

### 1.2 現行 closed-loop 計算の監査点

`src/model_2g3e_closedloop.py` では、励起状態間結合を

\[
A_{23}=A_{32}=J_{23}
\]

としている。一方、通常の受動的 Markovian coherence block は

\[
A(\omega)=\Gamma+i\!\left(H-\omega I\right),
\qquad
\Gamma=\Gamma^\dagger\ge0,
\qquad
H=H^\dagger
\]

と書かれるため、Hermitian なコヒーレント結合 \(J_{23}\) は原則として

\[
A_{23}=iJ_{23}
\]

に入る。

現行実装の実数非対角項は

\[
\Gamma_{\mathrm{eff}}
=
\frac{A+A^\dagger}{2}
\]

の非対角成分として働く。例として報告された \(J_{23}\simeq6.25\)、対角減衰 \(\gamma_e=1\) では、\(\Gamma_{\mathrm{eff}}\) が正半定値でなくなる。したがって、現在の厳密零点は「Hermitian mixing による passive escape」ではなく、非受動的または GKSL 条件を外れた有効散逸結合による零点である可能性が高い。

この点を隠すのではなく、次の二つに分けて再整理する。

1. **一般化 no-go の強化:** strict passive response では非対角励起混成を含めても零点が禁止されるか。
2. **旧 closed-loop 解の再分類:** passivity を破ったときに零点が出現する境界例として利用できるか。

---

## 2. 次段階で検証する中心命題

### 命題 A：strict-accretive generalized no-go

\[
A(\omega)=\Gamma+iK(\omega),
\qquad
\Gamma>0,
\qquad
K=K^\dagger
\]

とし、

\[
G=A^{-1},\qquad M=D^\dagger GD
\]

とする。このとき

\[
\frac{G+G^\dagger}{2}
=
G^\dagger\Gamma G>0
\]

であり、\(\operatorname{rank}D=2\) なら

\[
\frac{M+M^\dagger}{2}
=
D^\dagger G^\dagger\Gamma GD>0.
\]

したがって、

\[
\boxed{
\Gamma>0,\ \operatorname{rank}D=2
\Longrightarrow
\det(D^\dagger G D)\neq0
}
\]

が期待される。

重要なのは、これは以下をすべて許す点である。

- 任意の励起状態数 \(N_e\)
- 任意の Hermitian 励起状態混成
- 非可換な \(\Gamma\) と \(H\)
- 非対角・相関散逸
- non-normal な \(A\)
- 任意の実 detuning

したがって、「対角 \(G\) が障害」ではなく、より本質的には「strict positive-real / accretive response が障害」である可能性がある。

### 命題 B：零点の出現は passivity margin の消失と一致する

passivity margin を

\[
m_\Gamma=\lambda_{\min}\!\left(\frac{A+A^\dagger}{2}\right)
\]

と定義する。現行 closed-loop 模型の厳密零点は、

\[
m_\Gamma\le0
\]

となる領域でのみ出現すると予想する。

数値計算では、零点の位置、\(\sigma_{\min}(M)\)、および \(m_\Gamma\) を同時追跡し、

\[
m_\Gamma\downarrow0
\]

に伴って transfer zero がどのように実軸へ接近するかを調べる。

### 命題 C：有限基底デコヒーレンスは別の零点機構を許すか

有限 \(\gamma_g\) を含む Schur 応答を

\[
\chi_{\mathrm{full}}
=
S_p-
\frac{\beta K_{pc}K_{cp}}
{\gamma_g+\beta S_c}
\]

とし、分子を

\[
N_{\mathrm{full}}
=
\gamma_g S_p
+
\beta\left(S_pS_c-K_{pc}K_{cp}\right)
\]

とする。

\(\gamma_g=0\) では generalized no-go が直接作用するが、\(\gamma_g>0\) では二項間の複素相殺により full-rank の実軸 zero が生じる可能性を数値的に確定する。

仮に零点が存在しても、それが

- 制御場誘起か
- phase coherence に依存するか
- gain や population inversion を伴わないか
- 単なる散逸性 spectral hole でないか

を別途判定する必要がある。

---

## 3. 数値キャンペーンの全体構成

次段階を **Phase F: Passivity audit and legitimate escape campaign** とする。

---

# Gate F0：既存結果の凍結と完全再現

## 目的

現在の結果を変更前に保存し、実装修正後の結果と一対一で比較できるようにする。

## 実行内容

1. `model_2g3e_closedloop.py` を現状のまま再実行する。
2. 発見された全解について以下を保存する。
   - \(\delta_0,\Delta_2,\Delta_3,J_{23},\phi\)
   - \(|N_{\mathrm{full}}|\)
   - \(\operatorname{rank}D\)
   - \(\sigma_{\min}(D)\)
   - 極と零点
   - ATS 指標
3. 各解について
   \[
   \Gamma_{\mathrm{eff}}=(A+A^\dagger)/2
   \]
   を計算し、
   \[
   \lambda_{\min}(\Gamma_{\mathrm{eff}})
   \]
   を追加記録する。
4. コード、入力パラメータ、JSON、図の SHA-256 を保存する。

## PASS 基準

- 既報の代表解と零点残差を再現。
- 既存解の passivity margin が一覧化される。
- 以後、旧模型を `real_offdiagonal_effective_model` と明示し、物理的 Hermitian coupling 模型と混同しない。

## 出力

- `results/gate_F0_baseline_reproduction.json`
- `results/gate_F0_existing_zero_passivity_audit.csv`
- `figures/F0_zero_vs_passivity_margin.pdf`

---

# Gate F1：物理的な Hermitian mixing 実装への修正

## 目的

コヒーレント励起状態混成を、規約に依存しない形で正しく実装する。

## 標準模型

\[
A(\delta)=\Gamma+i(H-\delta I),
\]

\[
H=
\begin{pmatrix}
\Delta_1 & J_{12} & J_{13}\\
J_{12}^* & \Delta_2 & J_{23}\\
J_{13}^* & J_{23}^* & \Delta_3
\end{pmatrix},
\qquad
\Gamma=\Gamma^\dagger\ge0.
\]

単純な \(e_2-e_3\) 混成なら

\[
H_{23}=H_{32}=J_{23}\in\mathbb R
\]

とし、行列 \(A\) には \(iJ_{23}\) が入る。

## 実行内容

同一の \(\Delta_2,\Delta_3,J_{23},D\) に対し、次の三模型を比較する。

### Model R：旧実装

\[
A_{23}=J_{23}.
\]

### Model H：物理的 Hermitian mixing

\[
A_{23}=iJ_{23}.
\]

### Model D：物理的な相関散逸

\[
A_{23}=\kappa_{23},
\qquad
\Gamma=
\begin{pmatrix}
\gamma_1&0&0\\
0&\gamma_2&\kappa_{23}\\
0&\kappa_{23}^*&\gamma_3
\end{pmatrix}\ge0.
\]

Model D では Kossakowski/GKSL positivity に対応する制約を必ず課す。単純な等減衰の場合は概ね

\[
|\kappa_{23}|\le\sqrt{\gamma_2\gamma_3}
\]

を満たす範囲のみを物理領域とする。

## 計算量

- 旧解を初期値とする continuation
- ランダム多点始動
- \(\delta,J,\phi,\Delta_j\) の局所最適化
- real-axis root と complex zero の両方を追跡

## PASS 基準

以下の恒等式残差が \(10^{-11}\) 以下。

\[
\frac{G+G^\dagger}{2}-G^\dagger\Gamma G=0.
\]

物理領域では

\[
\lambda_{\min}(\Gamma)\ge0,
\qquad
\lambda_{\min}(\operatorname{Re}G)\ge0
\]

を満たす。

## 最重要出力

旧実装から物理実装へ連続変形する補間

\[
A_{23}(\alpha)
=
(1-\alpha)J_{23}+i\alpha J_{23},
\qquad
0\le\alpha\le1
\]

を導入し、零点軌跡

\[
z_0(\alpha)
\]

を複素周波数平面で描く。

期待される結果は、\(\alpha\to1\) および \(m_\Gamma>0\) への移行に伴い、零点が実軸から離れる、または消滅することである。

## 出力

- `src/audit_closedloop_convention.py`
- `results/gate_F1_convention_comparison.json`
- `figures/F1_realJ_vs_iJ_spectra.pdf`
- `figures/F1_zero_trajectory_complex_plane.pdf`
- `figures/F1_passivity_eigenvalues.pdf`

---

# Gate F2：strict-accretive generalized no-go の大規模数値検証

## 目的

一般化 no-go を、対角 \(G\) だけでなく任意の passive non-normal excited manifold に拡張できるか確認する。

## ランダム模型生成

励起状態数

\[
N_e=2,3,\ldots,8
\]

について、

\[
H=H^\dagger
\]

をランダム生成する。

散逸行列は

\[
\Gamma=LL^\dagger+\gamma_{\mathrm{floor}}I
\]

として生成し、

\[
\gamma_{\mathrm{floor}}>0
\]

により strict passivity を保証する。

光学結合 \(D\in\mathbb C^{N_e\times2}\) は複素ランダム行列から生成し、

\[
\sigma_{\min}(D)>\varepsilon_D
\]

を課す。

## 探索量

各模型で実軸 detuning を走査し、

\[
M(\delta)=D^\dagger A(\delta)^{-1}D
\]

について以下を記録する。

- \(|\det M|\)
- \(\sigma_{\min}(M)\)
- \(\lambda_{\min}(\operatorname{Re}M)\)
- \(\lambda_{\min}(\operatorname{Re}G)\)
- \(m_\Gamma\)
- \(A\) の non-normality
  \[
  \|[A,A^\dagger]\|
  \]
- condition number
- 最接近 complex zero

粗いランダム探索の後、\(\sigma_{\min}(M)\) が小さいサンプルだけを局所最適化する。

## 理論下限との比較

次の下限を数値的に検証する。

\[
\lambda_{\min}(\operatorname{Re}M)
\ge
\frac{
\lambda_{\min}(\Gamma)\,
\sigma_{\min}(D)^2
}{
\sigma_{\max}(A)^2
}.
\]

比

\[
Q=
\frac{
\lambda_{\min}(\operatorname{Re}M)
}{
\lambda_{\min}(\Gamma)\sigma_{\min}(D)^2/\sigma_{\max}(A)^2
}
\]

を計算し、数値誤差を除いて \(Q\ge1\) となることを確認する。

## 境界スキャン

以下の二つの no-go 境界へ接近する。

### Rank boundary

\[
\sigma_{\min}(D)=\epsilon_D\to0.
\]

### Passivity boundary

\[
\lambda_{\min}(\Gamma)=\epsilon_\Gamma\to0.
\]

\(\sigma_{\min}(M)\) と \(|\det M|\) の漸近スケーリングを測定する。

これにより、実軸零点へ接近するには

- optical dark vector の生成
- lossless/decoherence-free mode の生成
- passivity の破れ

のいずれかが必要かを可視化する。

## PASS 基準

- strict passive、full-rank の全サンプルで real-axis zero がゼロ件。
- 恒等式残差 \(<10^{-11}\)。
- 下限違反が数値丸めの範囲外でゼロ件。
- non-normality を大きくしても no-go が破れない。
- 小さい \(|\det M|\) が rank または passivity boundary への接近で説明される。

## 出力

- `src/model_passive_general.py`
- `results/gate_F2_passive_random_scan.json`
- `results/gate_F2_boundary_scaling.csv`
- `figures/F2_min_singular_value_vs_passivity.pdf`
- `figures/F2_boundary_scaling_rank_and_loss.pdf`
- `figures/F2_nonnormality_no_escape.pdf`

---

# Gate F3：旧 closed-loop 零点の物理的再分類

## 目的

現在得られている厳密零点を「誤り」として破棄するのではなく、passivity-breaking transition の数値例として整理する。

## 実行内容

全旧解について、

\[
m_\Gamma
=
\lambda_{\min}\!\left(\frac{A+A^\dagger}{2}\right)
\]

と以下の量を対応づける。

- \(|\chi_{\mathrm{full}}(\delta_0)|\)
- \(\operatorname{Im}z_0\)
- 最小吸収
- gain の有無
- pole stability
- \(\sigma_{\min}(M)\)

次に \(J_{23}\) を連続変化させ、相関散逸として許される領域

\[
\Gamma\ge0
\]

と禁止領域

\[
\Gamma\ngeq0
\]

の境界を跨いで零点を追跡する。

## 判定

### 期待結果 1

すべての厳密零点が

\[
m_\Gamma\le0
\]

に対応する。

この場合、旧解は

> exact dark-state-free zero is obtained only after the strict-passive condition is lost

という no-go 境界の対照例として利用する。

### 期待結果 2

\[
m_\Gamma\ge0
\]

でも零点が残る。

この場合は、\(\Gamma\) の半正定値境界、零減衰モード、特異な定常状態、または実装規約を再監査する。strict positive definite と positive semidefinite を明確に区別する。

## 出力

- `results/gate_F3_existing_solutions_reclassified.csv`
- `figures/F3_zero_emergence_at_passivity_boundary.pdf`
- `results/F3_interpretation_note.md`

---

# Gate F4：有限 \(\gamma_g\) と制御強度を含む Schur zero の探索

## 目的

\(\gamma_g=0\) における no-go が、有限の基底コヒーレンス減衰でどう変化するかを確定する。

## パラメータ

対数スキャンを用いる。

\[
\gamma_g/\gamma_e
=
10^{-6},10^{-5},\ldots,10^{1},
\]

\[
\beta/\gamma_e^2
=
10^{-4},10^{-3},\ldots,10^{4}.
\]

励起状態数は最初に \(N_e=2,3\)、必要なら \(N_e=4\) へ拡張する。

## 零点探索

分子

\[
N_{\mathrm{full}}(\delta)
=
\gamma_gS_p
+
\beta\det(D^\dagger GD)
\]

について、

\[
\operatorname{Re}N_{\mathrm{full}}=0,
\qquad
\operatorname{Im}N_{\mathrm{full}}=0
\]

を多点始動 root solver で解く。

同時に分母

\[
Q_{\mathrm{full}}
=
\gamma_g+\beta S_c
\]

がゼロでないことを確認する。

## 必ず分けて報告する量

1. **複素 transfer zero**
   \[
   \chi_{\mathrm{full}}=0.
   \]

2. **吸収のみの zero**
   \[
   \mathcal A[\chi_{\mathrm{full}}]=0,
   \qquad
   \chi_{\mathrm{full}}\neq0.
   \]

3. **有限コントラスト dip**
   \[
   0<C_S<1.
   \]

この三者を混ぜない。

## phase-coherence test

候補零点に対して、

- \(\gamma_g\) 増加
- 制御位相のランダム化
- \(K_{pc}K_{cp}\) の位相平均
- sector coupling の切断

を行い、零点または dip が破壊されるか確認する。

## control-induced test

\[
\beta\to0
\]

の零点軌跡を追跡し、

\[
N_{\mathrm{full}}(\delta_0;\beta=0)\neq0
\]

を確認する。もともと存在する background zero は除外する。

## PASS 基準

次のいずれかを確定する。

### F4-A：有限 \(\gamma_g\) でも zero なし

generalized no-go を有限基底減衰へ拡張する証明候補を作る。

### F4-B：有限 \(\gamma_g\) で zero あり

その zero が 7 条件と passivity 条件を満たすか判定する。満たさなければ「散逸補助 zero」であり response-EIT とは区別する。

## 出力

- `src/scan_finite_ground_dephasing.py`
- `results/gate_F4_ground_dephasing_scan.json`
- `figures/F4_gamma_g_beta_phase_diagram.pdf`
- `figures/F4_zero_trajectory_vs_gamma_g.pdf`
- `figures/F4_phase_randomization_test.pdf`

---

# Gate F5：正当な passive escape structure の最小探索

## 目的

strict-accretive no-go が成立する場合、どの仮定を物理的に変更すれば、gain を使わずに dark-state-free exact zero を構成できるかを調べる。

「とりあえず複雑な NV 全模型を回す」という人類にありがちな煙幕は避け、低ハードル順に進める。

---

## F5-A：非対称 input-output または直接背景経路

### 模型

一般の観測応答を

\[
\chi(\omega)
=
\chi_{\mathrm{dir}}(\omega)
+
p^\dagger G(\omega)c
\]

とする。

現在の \(D^\dagger GD\) no-go は、入力と検出が同じ双極子構造に拘束される matched input-output に強く依存する。偏光選択、干渉計測、直接散乱経路、複数の放射チャネルを含めると、

\[
p\neq c
\]

または \(\chi_{\mathrm{dir}}\neq0\) が自然に現れ、passive system でも Fano 型 transfer zero が可能になる。

### 探索条件

- \(\Gamma>0\)
- Liouvillian/pole stability
- gain なし
- optical coupling full rank
- pure dark state なし
- zero は control-induced
- sector cut で zero が消える
- ATS ではない
- 直接背景だけによる trivial zero ではない

### 判定上の注意

この構成が「EIT」と呼べるかは、位相コヒーレンスと sector response の操作的条件で判断する。単なる外部干渉計の暗ポートは除外する。

---

## F5-B：完全 GKSL \(2g+3e\) 模型

### 目的

手作業で構成した optical coherence block ではなく、正の Kossakowski 行列を持つ完全 Lindblad 方程式から線形感受率を導出する。

### Hamiltonian

\[
H=H_g+H_e+H_{\mathrm{mix}}+H_p+H_c.
\]

励起混成 \(H_{\mathrm{mix}}\) は Hermitian とする。

### jump operators

- 励起状態から基底状態への放射緩和
- 励起純粋デコヒーレンス
- 基底デコヒーレンス
- 必要なら相関放射 jump
- population inversion を生む pump は最初は導入しない

### 数値方法

1. 制御場のみで定常状態 \(\rho_0\) を求める。
2. 弱プローブに対する線形化 Liouvillian
   \[
   \mathcal L_0\,\delta\rho
   =
   -\mathcal L_p\rho_0
   \]
   を解く。
3. \(\chi_{\mathrm{full}}\)、\(\chi_{\mathrm{cut}}\)、\(R_S\) を同一の観測演算子で評価する。
4. 直接 resolvent 模型と比較する。
5. Liouvillian の complete positivity、trace preservation、steady-state uniqueness を検証する。

### 重要な確認

- reduced optical block が strict accretive か
- full Liouvillian では reduced-block no-go を回避できるか
- population dynamics や shelving が零点に不可欠か
- 零点が coherent cancellation か population hole か

---

## F5-C：補助 passive mode / pseudomode

\(2g+3e\) に passive auxiliary mode を結合し、全体は Markovian GKSL のまま、縮約応答に周波数依存 self-energy を生成する。

候補：

- cavity-like mode
- phonon pseudomode
- metastable auxiliary coherence
- singlet-like shelving coherence

目的は、active gain や非物理的散逸を用いずに、追加の伝達経路が numerator zero を生成できるか調べることである。

---

## F5-D：優先度を下げる構造

以下は F5-A〜C が失敗した後にのみ実施する。

- genuine non-Markovian memory kernel
- explicit gain
- population inversion
- PT-symmetric/non-Hermitian active model
- exceptional-point tuning

これらで零点が得られても、PRL の中心的な「passive EIT 定義」の構成例としては弱くなる。

## F5 PASS 基準

少なくとも一つの模型が、後述の拡張10条件をすべて満たす。

---

# Gate F6：拡張10条件チェックリスト

既存の7条件に、物理性を判定する3条件を追加する。

## 応答・構造条件

1. **実軸 exact zero**
   \[
   |\chi_{\mathrm{full}}(\omega_0)|<10^{-10}.
   \]

2. **optical dark vector なし**
   \[
   \sigma_{\min}(D)>10^{-6}
   \]
   または対応する完全結合行列が full rank。

3. **定常純 Lindblad dark state なし**
   jump operators の共通固有ベクトル条件と Hamiltonian 不変性を直接検査。

4. **sector cancellation**
   \[
   R_S(\omega_0)\neq0,
   \qquad
   R_S(\omega_0)=-\chi_{\mathrm{cut}}^{(S)}(\omega_0).
   \]

5. **sector-cut destruction**
   \[
   \mathcal A[\chi_{\mathrm{cut}}^{(S)}](\omega_0)>0.
   \]

6. **ATS 除外**
   極分離、半値幅、留数位相、二峰分離を用いる。

7. **control-induced**
   制御場をゼロにすると同じ周波数領域の零点が消える。

## 物理性条件

8. **passivity / no gain**
   観測帯域で負吸収または出力増幅がない。縮約模型では \(\Gamma\ge0\)、完全模型では pump-free GKSL とエネルギーフローを確認。

9. **stability and regularity**
   全 pole が安定側にあり、\(\omega_0\) は pole でない。定常状態が一意である。

10. **robustness**
    主要パラメータの \(0.1\%\)、\(1\%\)、\(5\%\) 摂動に対し、
    - zero の移動
    - 最小吸収
    - コントラスト
    - linewidth
    を評価する。機械精度だけの fine-tuned zero は構成定理の例には使えても、実験予言とは区別する。

---

# Gate F7：NV 励起6準位への写像

## 開始条件

F5 で passive かつ dark-state-free の構成が10条件を満たした場合のみ着手する。

抽象模型で機構が未確定のまま NV の全パラメータを投入すると、構造的原因が見えなくなるためである。

## 写像対象

NV 励起状態の

- spin-orbit coupling
- spin-spin coupling
- strain/electric-field mixing
- transverse magnetic-field mixing
- orbital relaxation
- radiative/ISC channels

から、F5 で得られた最小 escape structure に対応する部分を抽出する。

## 実行順

1. NV 励起6準位 Hamiltonian の固有基底を構成。
2. probe/control dipole 行列 \(D\) を完全遷移テンソルから計算。
3. GKSL dissipator を追加。
4. \(\chi_{\mathrm{full}},\chi_{\mathrm{cut}},R_S\) を計算。
5. 抽象模型の zero または near-zero を continuation で NV パラメータ領域へ追跡。
6. 温度、歪み、磁場、制御強度に対する phase diagram を作成。
7. 実軸 exact zero が失われる場合も、complex zero の虚部と最大コントラストを報告する。

## 必要な結果

- NV における passivity margin
- zero trajectory
- dark-state rank
- sector response
- ATS 指標
- 温度・磁場・歪みに対する robustness
- 実験的に必要な制御強度と linewidth

---

## 4. 実装ファイル案

```text
EIT definition equivalence/
├── NEXT_NUMERICAL_PLAN.md
├── src/
│   ├── audit_closedloop_convention.py
│   ├── model_passive_general.py
│   ├── scan_passivity_boundary.py
│   ├── scan_finite_ground_dephasing.py
│   ├── search_passive_escape.py
│   ├── full_gksl_2g3e.py
│   ├── model_auxiliary_pseudomode.py
│   ├── map_escape_to_nv6.py
│   └── run_phase_f.py
├── results/
│   ├── gate_F0_baseline_reproduction.json
│   ├── gate_F1_convention_comparison.json
│   ├── gate_F2_passive_random_scan.json
│   ├── gate_F3_existing_solutions_reclassified.csv
│   ├── gate_F4_ground_dephasing_scan.json
│   ├── gate_F5_escape_candidates.json
│   ├── gate_F6_validation_summary.json
│   └── gates_summary_phaseF.json
└── figures/
    ├── F1_realJ_vs_iJ_spectra.pdf
    ├── F1_zero_trajectory_complex_plane.pdf
    ├── F2_min_singular_value_vs_passivity.pdf
    ├── F2_boundary_scaling_rank_and_loss.pdf
    ├── F3_zero_emergence_at_passivity_boundary.pdf
    ├── F4_gamma_g_beta_phase_diagram.pdf
    ├── F5_passive_escape_spectrum.pdf
    └── F6_ten_condition_summary.pdf
```

---

## 5. 推奨 Figure 構成

### Figure 1：実装規約監査

同一パラメータで

- 実数非対角 \(J\)
- Hamiltonian \(iJ\)
- PSD 相関散逸

のスペクトルと pole-zero を比較する。

**必要な主張:** 現在の exact zero は Hermitian mixing へ置換すると実軸から離れる。

### Figure 2：strict-passive universal no-go

横軸：

\[
m_\Gamma=\lambda_{\min}(\Gamma).
\]

縦軸：

\[
\min_\delta\sigma_{\min}(D^\dagger GD).
\]

色または記号で \(N_e\)、non-normality、\(\sigma_{\min}(D)\) を表す。

**必要な主張:** 次元や non-normality では no-go を破れず、rank/passivity boundary でのみ零点へ接近する。

### Figure 3：finite-\(\gamma_g\) phase diagram

\[
(\gamma_g,\beta)
\]

平面に、

- exact zero
- absorption-only zero
- finite coherent dip
- no transparency

を分類して表示する。

### Figure 4：最小 passive escape

F5 で成功した最小模型について、

- full/cut spectra
- sector response
- pole/residue
- zero trajectory
- dark-state test
- passivity
- robustness

を一枚または二枚にまとめる。

### Figure 5：NV 写像

抽象模型の escape parameter と NV の磁場・歪み・軌道混成を対応づけ、実現可能領域を示す。

---

## 6. 計算上の注意事項

### 6.1 `det` だけで判定しない

near-singular 行列では determinant は condition number に強く依存する。必ず

\[
\sigma_{\min}(M)
\]

を併記する。

### 6.2 実軸走査だけで終えない

complex zero \(z_0\) を直接求め、

\[
\operatorname{Im}z_0
\]

を exact zero からの距離として使う。

### 6.3 root solver の成功表示を信用しない

以下を独立に検証する。

- 残差
- Jacobian rank
- condition number
- 高精度演算による再評価
- 初期値依存性
- parameter perturbation

数値計算は、成功フラグを出すことに関してだけは妙に自信満々なので、こちらが疑う必要がある。

### 6.4 passivity と pole stability を分ける

全 pole が安定でも、\(\operatorname{Re}A\) が正半定値とは限らない。逆に、縮約行列の見かけだけで完全 GKSL の物理性を判断しない。

### 6.5 exact zero と実験的 EIT を分ける

機械精度の exact zero が極端な fine tuning を必要とする場合、

- 数学的構成定理
- 実験的に観測可能な deep transparency

を別結果として報告する。

---

## 7. 優先順位

### 最優先

1. F0：既存結果の凍結
2. F1：\(J\) の実装規約修正
3. F2：strict-accretive generalized no-go の検証
4. F3：既存238解の passivity 再分類

ここまでで、副産物が本当に「任意次元 passive no-go」へ昇格するかが決まる。

### 次点

5. F4：有限 \(\gamma_g\) と制御強度
6. F5-A：非対称 input-output / direct path
7. F5-B：完全 GKSL \(2g+3e\)

### 最後

8. F5-C：補助 mode
9. F7：NV 6準位写像

---

## 8. 最終的な分岐と論文化上の意味

### Outcome I：strict-passive no-go が成立し、passive escape なし

主結果は、

\[
\boxed{
\text{strict passive Markovian matched response}
+
\operatorname{rank}D=2
\Longrightarrow
\text{no exact dark-state-free real-axis zero}
}
\]

となる。

この場合、dark-state-free EIT は原則として finite-depth coherent transparency として定義し、perfect zero は dark state、lossless boundary、direct-path interference、または非標準構造を必要とすると整理する。

### Outcome II：finite \(\gamma_g\) または full GKSL で passive escape が見つかる

その模型が10条件を満たせば、

> dark state は EIT の必要条件ではなく、passive multichannel open-system interference が厳密 EIT zero を生成できる

という構成的結果になる。

### Outcome III：direct path / asymmetric readout のみで escape

「物質内部の EIT」と「観測系を含む coherent transparency」を区別する必要がある。これは EIT の定義を系単体ではなく input-output pair に付随させる理論へ発展する。

### Outcome IV：旧解のみが exact zero

旧解は失敗例ではなく、

> passivity-breaking boundary を越えると exact zero が生成される

ことを示す対照モデルとして残す。ただし、通常の受動 EIT の構成例としては使用しない。

---

## 9. Phase F の最終判定表

| Gate | 問い | PASS 条件 |
|---|---|---|
| F0 | 旧結果を再現・保存できるか | 既報解と残差を再現 |
| F1 | 物理的 \(iJ\) 実装で零点は残るか | 規約比較と zero trajectory 完成 |
| F2 | strict-passive generalized no-go は破れるか | 反例ゼロ、恒等式・下限を確認 |
| F3 | 旧零点は非受動境界に対応するか | 全解を passivity で分類 |
| F4 | 有限 \(\gamma_g\) が零点を許すか | exact/absorption/dip を分類 |
| F5 | 正当な passive escape が存在するか | 少なくとも一例が候補化 |
| F6 | 候補は EIT の10条件を満たすか | 全条件 PASS |
| F7 | NV に写像できるか | zero/near-zero と実験条件を提示 |

---

## 10. 推奨される直近の実行単位

最初に実装すべき計算は、次の一つの統合スクリプトである。

```text
src/audit_closedloop_convention.py
```

このスクリプトは同一パラメータに対し、

1. 旧 real-\(J\) 模型
2. physical \(iJ\) 模型
3. PSD correlated-damping 模型

を計算し、

- \(\Gamma_{\mathrm{eff}}\) の固有値
- full/cut spectra
- pole-zero
- \(\sigma_{\min}(D^\dagger GD)\)
- zero trajectory
- gain/passivity 判定

を一括出力する。

この結果だけで、

- 現在の構成的存在結果が維持されるか
- 一般化 no-go がさらに強くなるか
- 次に有限 \(\gamma_g\) へ進むべきか
- passive escape 探索へ進むべきか

を決定できる。

---

## 結論

次段階の数値計算で最も重要なのは、新たなパラメータを大量に探索することではない。まず

\[
\boxed{
\text{exact zero}
\quad\text{と}\quad
\text{strict passivity}
}
\]

の両立可能性を決着させることである。

推奨する論理順序は、

\[
\boxed{
\text{実装監査}
\rightarrow
\text{generalized no-go の数値証明書}
\rightarrow
\text{有限 }\gamma_g
\rightarrow
\text{最小 passive escape}
\rightarrow
\text{完全 GKSL}
\rightarrow
\text{NV 写像}
}
\]

である。

この順序なら、現在の closed-loop 解が物理的に維持されなくても、任意次元・非対角・non-normal passive response に対するより強い no-go が主結果として残る。逆に正当な passive escape が見つかれば、それは dark-state-free EIT の構成的存在定理として、当初より明確で強い結果になる。
