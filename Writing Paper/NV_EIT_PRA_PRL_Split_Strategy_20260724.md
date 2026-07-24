# NV–EIT研究：PRA／PRL分割戦略

作成日：2026-07-24

## 0. 結論

現状の成果は、問いの異なる二報に分割するのが最も明瞭である。

- **PRA**：NV中心に、低温EITが温度上昇で崩壊する微視的機構、成立温度域、横磁場による救済経路、室温no-goを定量化する。
- **PRL**：EITを一つの実現例に格下げし、高速散逸部分空間を介する一般の線形応答について、選択則による経路相殺が整数の散逸抑制次数を選ぶ、という物質・現象横断の応答則を提示する。

二報の関係は、

\[
\text{PRA：一つの物理系を深く閉じる}
\qquad
\text{PRL：そこから抽出した一般法則を別の物理系へ広げる}
\]

である。

現状では **PRAは成立圏内**である。一方、PRLは有望だが、一般式が広いだけでなく、EIT外・diamond外の物理的witnessとobservable-levelの次数則を示す必要がある。

---

## 1. 二報を分ける理由

PRAとPRLでは、答える問いが異なる。

| 項目 | PRA | PRL |
|---|---|---|
| 中心対象 | NV\(^{-}\)のEIT | 高速散逸を介する有限次元Markovian線形応答 |
| 中心問い | NVの低温EITはなぜ失われ、どこまで回復できるか | 散逸下の応答抑制次数を、入力・選択則・読出し構造から予測できるか |
| 主役 | NVのfull physical model | 一般定理と複数の異なる物理的witness |
| 主要変数 | \(T,B_\perp,\Omega,\Delta\) | 高速率 \(\Gamma\)、経路モーメント \(M_j\)、observable |
| 結論 | NVの温度境界と室温no-go | selection-rule cancellationによる整数応答クラス |
| 適切な読者 | NV、EIT、固体量子光学 | 開放系、量子光学、物性、散逸工学 |

PRLをPRAの短縮版にしてはならない。PRAで得たNVの結果はPRLの動機・代表例にはなるが、PRLはNVを除いても独立に成立する中心命題を持つ必要がある。

---

# Part I — PRA：NVにおけるEITの成立限界

## 2. PRAが掲げる問い

> 低温で観測されるNV\(^{-}\)のspin-\(\Lambda\) EITは、なぜ温度上昇とともに失われるのか。その崩壊は外部制御でどこまで回復でき、室温では何が観測されるのか。

これは「dephasingが大きいためEITが消える」という現象論を越え、暗状態を壊す微視的経路、温度境界、回復経路、信号の符号まで明らかにする問いである。

## 3. PRAの中心主張

> NV\(^{-}\)のZPL spin-\(\Lambda\) EITは、温度上昇に伴う励起状態軌道枝の高速混合により、EITを担う先頭Raman経路が相殺されて崩壊する。横磁場は副次的経路を \(B_\perp^2\) で開くが、室温では真の透明化を回復できず、応答はcontrol-induced absorptionへ符号反転する。

一文に圧縮すると、

> NV EITの温度限界を決める微視的機構を特定し、低温での成立領域と室温no-goを定量化した。

## 4. PRAに組み込む既存結果

### 4.1 正の対照

- 低温NV EITをfull Liouvillianで再現。
- 既存の低温実験と整合するgenuine EITを確認。
- EIT、ATS、未分解構造を区別する判定系を構築。

### 4.2 崩壊機構

- 温度上昇に伴う励起状態軌道混合を導入。
- 先頭Raman経路の相殺と、高次経路への押し上げを確認。
- 「単なる線幅増大」ではなく、経路干渉構造がNVの脆弱性を決めることを示す。

### 4.3 室温no-go

- 300 Kでの大域パラメータ探索を実施。
- 補正機構、full/reduced model比較、detuning window、観測信号換算を監査。
- 以前の \(|C|\) 最大化が吸収増加を透明化と誤認し得る点を修正。
- 300 Kでは \(C<0\) が維持され、control-induced absorptionへ移ることを確認。

### 4.4 横磁場

- 横磁場による副次的経路の開口が \(B_\perp^2\) 則に従う。
- ただし、これは必ずしも室温EITの回復ではなく、抑制された応答のprefactor変化である。

### 4.5 材料比較

- SiV/SnV orbital-\(\Lambda\)はNVのspin-\(\Lambda\)と経路構造が異なる。
- SiV/SnVは構造的な厳密no-goではないが、室温ではphonon・coherence条件による実用的no-goになり得る。
- PRAでは材料比較を設計原理・展望として扱い、定量的な主役はNVに置く。

## 5. PRAに残る必須計算

### P1. Full-Liouvillian \(T\)–\(B_\perp\)相図

30–100 K程度の領域で、各点を以下に分類する。

- genuine positive transparency
- control-induced absorption
- ATS
- unresolved／inconclusive

**合格条件**：符号、lineshape、EIT/ATS判定を同時に用いて、NV EITの成立上限を決定する。

### P2. 温度境界の不確かさ

phonon rate、ground-state dephasing、Rabi frequency、strain、linewidth等の不確かさを伝播させる。

**合格条件**：「約80 Kまでgo」と固定せず、例えば「50–80 K域の構成依存境界」のように信頼区間または成立帯域として提示する。

### P3. 実験observableへの変換

理論的感受率を透過、OD、PL、contrast、必要SNRへ変換する。

**合格条件**：低温の正の信号、遷移域、300 Kの負の信号が実験でどう見えるかを同一尺度で示す。

### P4. モデル監査

- full/reduced一致域
- truncation依存性
- parameter provenance
- 符号規約
- EIT/ATS分類の安定性

を再確認する。

## 6. PRAの図構成

1. NV level schemeと温度依存の軌道混合経路
2. 低温EITから高温吸収増加への代表スペクトル
3. full-Liouvillian \(T\)–\(B_\perp\) phase diagram
4. 温度に対するsigned contrastと不確かさ帯
5. 横磁場による \(B_\perp^2\) scaling
6. 実験observableと検出可能領域

## 7. PRAのストーリー

1. NV EITは低温で既知だが、室温へ延長できない微視的理由は未整理である。
2. full physical modelにより低温EITを正しく再現する。
3. 温度上昇により軌道混合がRaman干渉経路を破壊する。
4. 横磁場は副経路を開くが、室温では透明化を救済できない。
5. 温度境界と実験信号を与え、NVでの探索可能領域とno-go領域を確定する。

---

# Part II — PRL：EITを越えた経路モーメント応答則

## 8. PRLが掲げる問い

> 高速散逸部分空間を介する応答は、散逸率や線幅だけでなく、入力・内部力学・読出し経路の構造から普遍的に分類できるか。

より物理的には、

> 同程度の高速緩和率を持つ系が、なぜ異なる整数冪でコヒーレント応答を失うのか。その次数を、詳細な数値計算に先立って選択則から予測できるか。

この問いはEIT固有ではない。EITは、入力と読出しが光学dipoleで与えられる代表的witnessである。

## 9. PRLの一般定理候補

高速部分空間を介する弱い線形応答kernelを

\[
K_{pc}(\Gamma,z)
=
p^\dagger[\Gamma D+A_0(z)]^{-1}c
\]

とする。ここで、

- \(c\)：入力・駆動ベクトル
- \(p\)：読出しベクトル
- \(D\)：高速散逸・混合演算子
- \(A_0\)：低速力学、detuning、内部結合
- \(\Gamma\)：高速散逸率

である。高速部分空間上で \(D\) が可逆で、漸近展開が許されるなら、

\[
K_{pc}(\Gamma,z)
=
\sum_{j=0}^{\infty}
\frac{M_j(z)}{\Gamma^{j+1}},
\]

\[
M_j(z)
=
(-1)^j p^\dagger
\left[D^{-1}A_0(z)\right]^jD^{-1}c .
\]

したがって、

\[
M_0=\cdots=M_{m-1}=0,\qquad M_m\neq0
\]

ならば、

\[
K_{pc}(\Gamma,z)\sim M_m(z)\Gamma^{-(m+1)}.
\]

中心命題は次のように表せる。

> 高速散逸を介する応答の抑制次数は、散逸率だけではなく、入力・内部力学・読出しを結ぶ最初の非零経路モーメントによって選択される。選択則による低次モーメントの相殺は、系を異なる整数応答クラスへ分類する。

簡潔な帰結は、

\[
\boxed{
\text{selection-rule cancellation}
\Longrightarrow
\text{first nonzero path moment}
\Longrightarrow
\text{integer dissipation order}
\Longrightarrow
\text{observable response class}
}
\]

である。

## 10. 適用範囲

この法則は「どの材料にも無条件で成立する」のではない。材料非依存だが、少なくとも以下を仮定する。

- 有限次元のMarkovian dynamics
- 線形または弱プローブ応答
- 高速部分空間と低速部分空間の分離
- 高速部分空間における適切な逆演算子またはDrazin型縮約
- 観測量が入力—内部力学—読出しのresponse kernelで表される

直接の適用範囲外または追加理論が必要な例は、

- 強い非Markovian memory
- 連続スペクトルと非孤立特異点
- 高速演算子の特異kernelが支配的な場合
- 強プローブ・非線形応答
- 相分離や熱力学極限に固有の非解析性

である。

## 11. EIT以外の適用候補

同じ構造は、以下に現れ得る。

- Raman遷移・Raman散乱
- pump–probe線形応答
- cavity-mediated spin–photon conversion
- 散逸下の状態移送・有効結合
- 多準位系の線形感受率
- 高速損失中間状態を介する輸送
- 開放系のセンサー応答

PRLでは、このうち少なくとも一つをEITとは独立したfull physical witnessとして実証する必要がある。

## 12. 現在得られているPRL候補の結果

### 12.1 経路モーメント階層

現時点の解析・共通pipelineでは、概略として、

| 系 | 最初の非零モーメント | kernelの漸近次数 |
|---|---:|---:|
| group-IV orbital-\(\Lambda\) | \(M_0\neq0\) | \(\Gamma^{-1}\) |
| NV \(m_s=0\leftrightarrow\pm1\) | \(M_0=0,\ M_1\neq0\) | \(\Gamma^{-2}\) |
| NV \(m_s=-1\leftrightarrow+1\) | \(M_0=M_1=0,\ M_2\neq0\) | \(\Gamma^{-3}\) |

が確認されている。

### 12.2 NVでの物理的認証

- \(M_0=0\)は数値的偶然ではなく、光学脚の直交性・選択則に由来する。
- reduced kernelとfull nine-level Liouvillianは高速混合域で整合する。
- kernel次数とsector observable次数が必ずしも同一でないことが判明している。
- kernel次数2に対し、sector observableでは漸近次数4、物理的有限領域ではpre-asymptotic次数3が現れ得る。

最後の点は、kernelの定理だけではPRLが閉じないことを意味する。

### 12.3 現状の限界

- group-IVモデルのphonon normalizationとdipole geometryには模式性が残る。
- SiV/SnVについてfull physical GKSLでの次数認証が不十分。
- 物理的witnessがEIT/Ramanとdiamond centerに偏っている。
- EIT外のfull physical modelによる検証が未完了。
- observable inheritanceの一般則が未完成。

したがって、現状は「一般的な候補法則と有力なEIT witness」がある段階であり、「分野横断則としてPRLに必要な認証が完了した段階」ではない。

---

## 13. PRLに必須の計算・証明

### Priority 1：Observable-order inheritance

kernelの次数が実測量へどう継承されるかを導出する。

例としてsector responseが

\[
R_{\mathrm{obs}}
\sim
K_{12}S_g^{-1}K_{21}
\]

なら、単一kernelの次数だけでなく、左右のkernel、低速Schur complement、分母の零点・相殺を考慮する必要がある。

目標は、

\[
\nu_{\mathrm{obs}}
=
F(n_{12},n_{21},\nu_{S_g},\text{cancellation data})
\]

の形の一般則または十分条件を得ることである。

**合格条件**

- 感受率、signed absorption、Raman rate等に対して次数継承を導出。
- generic caseとsymmetry-protected cancellationを区別。
- full Liouvillianで予測次数と一致。
- pre-asymptotic exponentと真の漸近指数を区別。

**失敗時の意味**

kernel分類だけではobservable分類を主張できず、PRLの中心命題は弱まる。PRAの設計原理に戻すべきである。

### Priority 2：EIT外の物理的witness

EITでない現象を一つ選び、同じ経路モーメント則から異なる整数次数を予測・検証する。

候補：

- cavity-mediated spin–photon conversion
- 高速損失中間状態を介するRaman散乱
- dissipative state transfer
- 多準位pump–probe response

**選定条件**

- EIT dark-stateの言い換えではない。
- observableが明確。
- full GKSLまたは同等の物理モデルで計算可能。
- 選択則を変えると次数が変化する。

**合格条件**

- 経路モーメントから事前に次数を予測。
- full modelのlog–log scalingが予測に一致。
- 対称性を微小に破ると高次から低次へcrossoverする。

### Priority 3：Diamond外のwitness

非diamond系を一例加える。

候補：

- SiC color center
- rare-earth ion
- quantum dot／cavity-QED emitter
- superconducting artificial atom

**合格条件**

- diamond固有の軌道構造に依存しない。
- 文献値または妥当な実験領域のパラメータを使用。
- 同じ定理で応答次数を予測できる。

EIT外witnessとdiamond外witnessを同一モデルで満たせれば最も効率がよい。

### Priority 4：Group-IV full-GKSL認証

SiV/SnVについて、

- 実測phonon rate
- dipole geometry
- polarization
- strain
- ground-state coherence

を組み込んだfull modelを構築する。

**合格条件**

- reduced kernelとfull modelが同じ漸近次数を示す。
- EIT/ATSの可視性とは独立に、signed responseのscalingを抽出できる。
- NVと異なるclassが模式値の産物でない。

### Priority 5：三クラスのscaling collapse

\(n=1,2,3\)の物理モデルについて、適切なprefactorとrateで規格化し、

\[
\Gamma^n R(\Gamma)
\]

が一定値へ収束することを示す。

**合格条件**

- 十分な漸近windowがある。
- fit windowを変えても指数が安定。
- full/reduced modelの双方で一致。
- finite-\(\Gamma\) correctionを明示。

### Priority 6：Robustnessとcrossover fan

対称性保護された高次classに小さな破れ \(\epsilon\) を加えると、

\[
R(\Gamma,\epsilon)
\sim
\frac{\epsilon a_0}{\Gamma}
+
\frac{a_1}{\Gamma^2}
+\cdots
\]

のように、有限域では高次、真の漸近域では低次へ戻り得る。

**計算内容**

- strain
- detuning
- polarization error
- magnetic-field misalignment
- dipole mixing
- ground-state dephasing

に対する指数mapを作る。

**合格条件**

- exact classとapproximate classを区別。
- crossover scale \(\Gamma_\ast(\epsilon)\)を導出・数値確認。
- 実験で高次classを観測できるparameter windowを示す。

### Priority 7：Blind prediction

定理の恣意性を抑えるため、full calculationを見る前に、

1. Hamiltonian、jump operators、input、readoutだけから \(M_j\) を計算。
2. 最初の非零モーメントから指数を予測。
3. 予測を固定。
4. full numerical sweepを実行。

というblind protocolを採用する。

**合格条件**

- 少なくとも一つの新規モデルで予測が一致。
- fit後にclass定義を変更しない。

### Priority 8：実験識別可能性

**計算内容**

- 必要な \(\Gamma\) dynamic range
- signal magnitude
- SNR
- parameter uncertainty
- competing background
- finite-temperature mapping \(\Gamma(T)\)

**合格条件**

- 二つ以上のclassを実験的に区別できる。
- 単なる漸近数学ではなく、測定可能な設計則になる。

---

## 14. PRL計算の実行順序と判定ゲート

### Gate A：論理閉包

1. observable-order inheritance
2. generic／protected caseの定理
3. pre-asymptoticとasymptoticの区別

**Aが不合格**：PRL化を止め、PRAの設計原理へ統合する。

### Gate B：EIT非依存性

1. 非EIT witnessを選定
2. blind prediction
3. full model verification

**Bが不合格**：「EITに広く成立する法則」まで主張を縮小し、PRA／PRA Letter相当を検討する。

### Gate C：物質非依存性

1. diamond外witness
2. group-IV full-GKSL
3. 三クラスscaling collapse

**Cが不合格**：diamond defectの応答分類として位置づける。PRLの分野横断性は弱い。

### Gate D：物理的重要性

1. robustness／crossover
2. measurable scaling
3. material・device design prediction

**Dが合格**：PRLの中心主張が完成する。

---

## 15. PRLの推奨ストーリー

### 15.1 導入：線幅中心の見方の限界

高速散逸は通常、線幅やdephasing rateで応答を抑えると理解される。しかし、同じ散逸率でも、入力と読出しを結ぶ低次経路が選択則で消えると、応答は異なる整数冪で失われる。

### 15.2 一般則

response resolventの経路モーメント展開を導出し、最初の非零モーメントがkernelの抑制次数を選ぶことを示す。

### 15.3 Observable theorem

kernel次数が感受率、散乱率、変換効率などへどう継承されるかを示し、generic case、symmetry-protected case、crossoverを分類する。

### 15.4 異なる物理的実現

- EIT／Raman：NVおよびgroup-IV center
- 非EIT応答：別の散乱・変換・移送モデル
- diamond外の物理系

で同じ予測則を検証する。

### 15.5 設計則

散逸耐性の探索では「小さい線幅」を探すだけでなく、

\[
M_0,\ M_1,\ M_2,\ldots
\]

をselection ruleとinterface geometryから検査し、必要な応答classを設計できる。

## 16. PRLの一文主張候補

英語：

> Selection-rule cancellations organize open-system responses into integer dissipation-order classes determined by the first nonvanishing input–dynamics–readout path moment.

日本語：

> 選択則による経路相殺は、入力・内部力学・読出しを結ぶ最初の非零モーメントを通じて、開放系応答を整数の散逸次数クラスへ分類する。

EITを含めた補助文：

> EITはこの法則の一例であり、同じ分類はRaman散乱、散逸媒介変換、pump–probe応答にも適用される。

---

## 17. PRLの図構成

1. 一般構造：input \(c\) → fast dissipative sector \(D\) → internal dynamics \(A_0\) → readout \(p\)
2. \(M_0,M_1,M_2\)相殺と \(\Gamma^{-1},\Gamma^{-2},\Gamma^{-3}\) class
3. observable inheritanceとcrossover
4. EIT witnessにおけるNV／group-IV scaling collapse
5. 非EIT・diamond外witnessのblind prediction
6. robustness mapと実験可能領域

PRL本文では4図程度へ圧縮し、詳細な材料modelと追加robustnessはSupplementへ送る。

---

## 18. 先行研究との差と新規性

数学的なresolvent／Neumann展開そのものは新しくない。したがって、新規性を展開式そのものに置いてはならない。

| 先行分野 | 既知の内容 | 本研究の新規性候補 |
|---|---|---|
| EIT/CPT標準理論 | dark state、dephasing、control強度 | 高速散逸下の整数応答classを選択則から予測 |
| Raman pathway interference | 多経路相殺 | 相殺次数を散逸応答の普遍指数へ接続 |
| adiabatic elimination／Zeno | \(1/\Gamma\)展開、有効力学 | input–dynamics–readout geometryとobservable classへ接続 |
| NV phonon研究 | orbital mixing rate、ISC | mixing rateからRaman/EIT抑制次数を導出 |
| SiV/SnV CPT | 低温dark state、phonon-limited coherence | 異なる材料を同一のpath-moment classで比較 |
| open-system response | Liouvillian resolvent、linear response | selection-rule cancellationによる有限個の整数分類とblind prediction |
| materials screening | linewidth、coherence time、dipole strength | 経路モーメントを新しいscreening descriptorとして導入 |

新規性の核は、

\[
\text{標準的な漸近展開}
\]

ではなく、

\[
\text{selection rule}
\rightarrow
\text{moment cancellation}
\rightarrow
\text{observable exponent}
\rightarrow
\text{physical universality class}
\rightarrow
\text{testable design rule}
\]

という接続全体にある。

## 19. 分野横断性

必要な認証を通過すれば、以下を横断する。

- **量子光学**：EIT、CPT、Raman、pump–probe
- **開放量子系**：Liouvillian response、強散逸極限、adiabatic elimination
- **物性・欠陥中心**：phonon mixing、dipole selection rule、strain
- **量子デバイス**：state transfer、spin–photon conversion、散逸工学
- **材料設計**：linewidth以外のscreening descriptor

ただし、EITとdiamond centerだけで検証した段階では、形式的には広くても物理的分野横断性は弱い。非EITかつdiamond外のwitnessが、PRL適性を大きく左右する。

---

## 20. PRAとPRLの重複回避

### PRAに残すもの

- NVのfull level structure
- 温度依存phonon model
- \(T\)–\(B_\perp\)相図
- 室温符号反転
- experimental signal estimate
- NVの温度境界とno-go

### PRLに残すもの

- 一般response theorem
- observable inheritance
- path-moment classification
- EIT外witness
- diamond外witness
- blind prediction
- scaling collapseとcrossover law

### 両方で使えるが役割を変えるもの

- NV：PRAでは主役、PRLでは一つのclass witness
- SiV/SnV：PRAでは比較・展望、PRLでは異なるclassの物理的認証
- 横磁場・strain：PRAでは制御変数、PRLではsymmetry-breaking perturbation

同じ図・同じ数値表を主要結果として二重使用せず、PRAは定量的物理、PRLは一般則の検証に役割を分ける。

---

## 21. 現時点の投稿判断

### PRA

**判定：成立圏内。**

full-Liouvillian \(T\)–\(B_\perp\)相図、温度境界の不確かさ、observable換算を終えれば執筆を閉じられる。

### PRL

**判定：中心概念は有望だが、未認証。**

現時点のままでは、

> 標準的なresolvent展開をdiamond EITへ適用した

と評価される危険がある。次の三条件が最低限必要である。

1. observable-levelの次数則
2. EIT外のfull physical witness
3. diamond外を含むblind prediction

これにrobustness、crossover、実験識別可能性が加われば、PRLとしての重要性と一般性が明確になる。

---

## 22. 推奨する直近の作業順

### PRAライン

1. full \(T\)–\(B_\perp\) phase diagram
2. 温度境界のuncertainty propagation
3. transmission／PL／SNR換算
4. PRA本文の執筆開始

### PRLライン

1. observable-order inheritanceの解析
2. 小規模な一般GKSL modelで定理を単体テスト
3. 非EIT witnessを一つ選定
4. blind predictorを固定
5. full modelで本番検証
6. diamond外witnessへ適用
7. 三クラスscaling collapse
8. robustness／crossover／measurability

PRAの執筆はPRL計算の完了を待つ必要はない。PRAで一般則を過度に展開せず、「より一般的なpathway-order classificationは別報で扱う」として新規性を温存する。

---

## 23. 最終的な論文像

### PRAの帰結

> NVにおける低温EITから室温control-induced absorptionへの遷移を、励起状態軌道混合と経路干渉から説明し、その成立温度域と実験的境界を確定した。

### PRLの帰結

> 高速散逸はすべてのコヒーレント応答を同じ \(1/\Gamma\) 則で消すのではない。選択則が低次経路を相殺すると、応答は \(\Gamma^{-2},\Gamma^{-3},\ldots\) の整数classへ移り、そのclassはEIT、散乱、変換など異なるobservableに予測可能な形で現れる。

この分割により、PRAは材料固有の完全性を、PRLは物理法則としての一般性を、それぞれ最大化できる。
