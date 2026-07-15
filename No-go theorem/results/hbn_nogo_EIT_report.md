# hBNへのEIT no-go/go理論の適用

作成日: 2026-07-15  
対象: **光学EIT**（強いcontrol光が、長寿命物質コヒーレンスを介して弱probe光の線形吸収を低減する現象）  
除外: 蛍光CPT、ODMR dark resonance、シェルビング、Autler–Townes分裂だけによるディップ、古典EIT analogue

## 1. 結論

### 現在の実験事実に対する判定

**室温hBNをobservable optical EITのgoとは判定できない。** hBNで室温スピン制御、室温光学スピン読出し、室温で約65 MHzの狭い光学線はそれぞれ報告されているが、同一欠陥で次の入力がそろっていない。

1. 二つの長寿命下準位と共有励起状態からなる光学Λ系
2. 複素光学双極子ベクトル \(d_1,d_2\)
3. 室温のground-coherence減衰 \(\gamma_g\)
4. 同じ欠陥の室温均一光学幅 \(\gamma_e\)
5. Raman kernel \(K_{12},K_{21}\) または完全な線形化Liouvillian
6. controlによる弱probe吸収の低減

したがって、理論分類は**Class 0（inputs insufficient）**であり、実験的には**未実証なので現在はnot-go**として扱う。ただし、厳密なClass II no-goも証明できない。

### 設計可能性に対する判定

最良の光学線幅と最良の室温スピンコヒーレンスを同一欠陥で実現でき、かつRaman kernelの品質が十分なら、室温goは数値的に可能である。したがってhBNは

> **established goではなく、room-temperature conditional-go candidate**

と位置づけるのが妥当である。

---

## 2. no-go/go理論で本来必要な判定

弱probeの縮約応答は

\[
\Xi=S_1-\frac{\beta K_{12}K_{21}}{\gamma_g+\beta S_2},
\qquad \beta=|\Omega_c|^2/4.
\]

対象セクターによる補正は

\[
\delta\Xi_S=-\frac{\beta K_{12}K_{21}}{\gamma_g+\beta S_2}.
\]

したがって正則点では、厳密なセクターno-goは

\[
K_{12}K_{21}=0
\]

で判定される。hBN文献では、特定欠陥についてこのkernelを構成できる双極子行列・励起多様体・jump operatorが未取得である。このため、今回の計算は**厳密なkernel証明ではなく、実験パラメータによる可視性budget**である。

また、フォノンギャップの観測だけからClass V protected goとは断定できない。Class Vには光学コヒーレンス減衰演算子について実際に \(\ker D\neq0\) が存在し、probe/controlの両方がそのkernelへ射影される必要がある。フォノン結合が小さいことと、厳密なdissipative kernelの存在は別である。

---

## 3. 文献から得られる主要入力

### 3.1 機械的脱結合エミッタ: 最良の室温光学候補

Dietrich et al.は3–300 Kでsub-100 MHzのPLE幅を報告し、図から約64.7 MHz（3 K）、75.5 MHz（100 K）、65.5 MHz（200 K）、62.8 MHz（300 K）が得られる。Hoese et al.は約61 MHzの均一幅と低周波フォノンスペクトル密度のギャップを報告した。

これは光学コヒーレンス側の室温gateを通過し得る。一方、そのエミッタのspin multiplicity、共有励起状態、spin-selective dipoleは未確定である。

### 3.2 Koch et al.: 同じ機構の再評価

別の機械的脱結合エミッタでは、5 KでFTL幅109(9) MHz、full homogeneous width 324(16) MHz、フォノンギャップ2.42(11) THzが得られた。固定scanのPLE幅は約50 KまでFTL付近にあるが、55 K以上で増大し、ギャップ閉鎖は約100–150 Kと推定された。また、時間領域のcoherent Rabi driveは20 Kでは確認され、30 Kではoverdampedであった。

したがって、「機械的脱結合hBNならすべて300 KでFTL」という一般化はできない。欠陥・局所配置依存が強い。

### 3.3 室温スピン候補

Stern et al.の炭素関連単一欠陥では、室温で \(S=1\)、\(D/h=1.959\) GHz、\(T_2^*=106(12)\) nsが得られている。周波数linewidthに変換すると

\[
\gamma_g^{*}\simeq\frac{1}{\pi T_2^*}=3.00\ \mathrm{MHz}.
\]

一方、同一欠陥についてresonant optical homogeneous linewidthと光学Λ双極子は報告されていない。

Whitefield et al.では室温で25%以上のエミッタに光学スピン読出しが観測された。代表例のPL幅2.76 nm at 781 nmは約1.36 THzに相当するが、これはPLE均一幅ではない。したがって厳密な \(\gamma_e\) としては使用できず、「未選別の光学bandwidthを使った悲観的proxy」としてのみ計算した。

### 3.4 B-centerと \(V_B^-\)

B-centerは共鳴蛍光、Mollow triplet、二光子干渉まで進んでおり、二準位光学コヒーレンスは高い。しかしEIT用の長寿命下準位対は確立していない。

\(V_B^-\)は室温で基底ZFS約3.5 GHz、励起ZFS約2.1 GHzが測定されているが、広い発光bandと強いISCを持ち、同じ励起準位へ結合する分解された光学Λ系およびresonant optical linewidthは未確定である。

---

## 4. 数値モデル

正確な \(K_{12}K_{21}\) がないため、理想Λ系のon-resonance吸収低減を基準にした条件付きproxyを使用した。

\[
C_{\mathrm{proxy}}
=\eta\frac{\Omega_c^2}{\Omega_c^2+4\gamma_e\gamma_g}.
\]

- \(\eta\in[0,1]\): Raman kernel、branching、偏光、共有励起状態の品質をまとめた未測定係数
- \(\eta=1\): 理想的に共通励起状態へ結合する上限
- 文献のFWHMを保守的に \(\gamma_e\) proxyとして使用
- 理論の \(\gamma_e\) がHWHM規約なら、必要 \(\Omega_c\) は表の値より \(1/\sqrt2\) 小さくなる

目標contrast \(C_*<\eta\) に必要なcontrolは

\[
\Omega_c=2\sqrt{\gamma_e\gamma_g\frac{C_*}{\eta-C_*}}.
\]

この式は可視性budgetであり、EITの存在条件そのものではない。\(K_{12}K_{21}=0\)なら、いかにcontrolを強くしてもEITは生じない。

---

## 5. 室温の計算結果

### 5.1 最良ケース: 65 MHzの室温光学線 + 3.00 MHzのspin width

異なる論文の最良値を同一欠陥で実現できると仮定する。

| 目標contrast | \(\eta=1\)で必要な \(\Omega_c\) |
|---:|---:|
| 50% | 27.9 MHz |
| 90% | 83.8 MHz |

\(\Omega_c=100\) MHzでは理想上限が

\[
C_{\rm ideal}=0.928.
\]

10%のEIT contrastを得るには \(\eta>0.108\) が必要である。したがって、kernelが理想値の約10%以上なら、光学・spin linewidth budgetだけからは室温EITは排除されない。

### 5.2 保守的spin width: \(\gamma_g=10\) MHz

| 目標contrast | \(\eta=1\)で必要な \(\Omega_c\) |
|---:|---:|
| 50% | 51.0 MHz |
| 90% | 153.0 MHz |

\(\Omega_c=100\) MHzで \(C_{\rm ideal}=0.794\)、10%観測には \(\eta>0.126\) が必要である。依然としてcontrol強度は原理的には非現実的ではない。ただしhBN欠陥の光学双極子と損傷・電離・加熱条件からレーザーpowerへ変換する入力は不足している。

### 5.3 温度広がりが残るエミッタ

Akbari et al.の6.5–120 Kの \(T^3\) 傾向を300 Kへ**stress-test extrapolation**すると、\(\gamma_e\approx16.7\) GHzとなる。この値は300 Kの測定値ではない。

- 50%: \(\Omega_c\approx0.447\) GHz
- 90%: \(\Omega_c\approx1.34\) GHz

この場合、弱いcontrolではpractical noに近づく。

Whitefield et al.の代表的PL幅1.36 THzを均一幅の悲観的proxyとして入れると、理想Λでも

- 50%: \(\Omega_c\approx4.04\) GHz
- 90%: \(\Omega_c\approx12.1\) GHz

を要する。ただしPL帯域をhomogeneous linewidthとして扱うことはできないため、この結果は上限的な警告であり、正確なEIT予測ではない。

### 5.4 励起多様体の分解に関するヒューリスティック

実際の定理のkernelではなく、未分解多様体の影響を可視化するため

\[
\eta_{\rm res}=\frac{(\Delta_{\rm sel}/\gamma_e)^2}{1+(\Delta_{\rm sel}/\gamma_e)^2}
\]

を補助的に使用した。代表値 \(\Delta_{\rm sel}=2\) GHzでは

- \(\gamma_e=65\) MHz: \(\eta_{\rm res}=0.999\)
- \(\gamma_e=16.7\) GHz: \(\eta_{\rm res}=0.014\)
- \(\gamma_e=1.36\) THz: \(\eta_{\rm res}=2.2\times10^{-6}\)

となる。ただし、ここで必要なのは実際の**励起状態の関連branch separation**であり、基底ZFSを代入してよいとは限らない。この図は定理の数値証明ではない。

---

## 6. 欠陥別のno-go/go判定

| configuration | 室温判定 | 根拠 |
|---|---|---|
| Dietrich/Hoese機械的脱結合エミッタ | Class 0 / conditional go | \(\gamma_e\sim65\) MHzは十分小さいが、spin Λとkernelが不明 |
| Stern炭素関連単一spin | Class 0 | \(\gamma_g\sim3\) MHzは良好だが、resonant \(\gamma_e\)と光学Λが不明 |
| Whitefield spin-active narrowband emitter | Class 0、現在はpractical not-go | 室温spin readoutはあるが、homogeneous PLE幅・coherent Λなし |
| \(V_B^-\) | Class VI/III候補だが未確定 | spin構造は既知、光学band・ISCが不利、kernel入力なし |
| B-center | Class 0 for EIT | coherent two-level optical emitterだが長寿命下準位対なし |

**同一欠陥として室温goが証明されたconfigurationはない。**

---

## 7. 主張の方針

### 採用すべき主張

> **Unambiguous optical EIT has not yet been demonstrated in a solid at room temperature. The no-go/go framework identifies why the presently characterized hBN configurations remain incomplete: no single defect has yet combined a verified optical Λ kernel, a narrow room-temperature homogeneous optical line, and a long-lived ground coherence. Nevertheless, the best measured hBN optical and spin parameters occupy a conditional-go region, making hBN a concrete target for the first room-temperature solid-state EIT experiment.**

日本語では次の構成がよい。

> **厳密な光学EITは室温固体で未実証である。現在同定されているhBN欠陥も、同一欠陥における光学Λ kernel、狭い室温均一幅、長寿命ground coherenceがそろっていないため、現状はgoではない。一方、別々のhBN欠陥で得られた最良パラメータを同一configurationで実現できれば、室温EITは数値的に許容される。hBNは普遍no-goの反例ではなく、no-go境界を越えるための具体的なconditional-go設計候補である。**

### 採用すべきでない主張

- 「hBNは室温goである」: 同一欠陥のΛとkernelが未実証
- 「hBNは室温no-goである」: 65 MHz級の室温光学線があり、最良値では可視性budgetを通過
- 「フォノンギャップがあるのでClass Vである」: \(\ker D\)の証明がない
- 「hBN全体に単一の \(\Gamma_e(T)\) がある」: 欠陥・局所構造依存が大きい

---

## 8. 室温goを確定する最小実験

1. 同一単一欠陥で二つのground statesと共有excited stateをPLEで同定する。
2. 偏光依存Rabi周波数から複素双極子 \(d_1,d_2\) を再構成する。
3. 5–300 Kで短時間homogeneous width、長時間spectral diffusion、\(T_2^*\)を同じ欠陥で測る。
4. 弱probe・強control条件でprobe吸収またはextinctionを測る。
5. control-offとの差ではなく、ground-coherence pathwayを切ったモデルと比較し \(\delta\chi_S\) を抽出する。
6. control-power、二光子detuning、\(\gamma_g\)依存とpole構造からATSを排除する。

単一エミッタではbulk transmissionではなく、共振器・waveguide中のprobe extinctionまたはcoherent scatteringをobservableに設定するのが現実的である。

---

## 9. ファイル構成

- `literature/`: hBNの線幅、温度依存、spin、欠陥構造に関する文献PDF
- `analysis/hbn_parameter_table.csv`: 文献パラメータ表
- `analysis/hbn_numerical_results.csv`: \(\eta\)・目標contrast別の必要control
- `analysis/room_temperature_scenarios.csv`: 室温scenario比較
- `figures/`: 温度幅、control budget、contrast、branch resolutionの図
- `code/hbn_nogo_scan.py`: 計算再現コード
- `theory/`: 適用したno-go/go理論資料

## 10. 主要参考文献

1. Dietrich et al., Phys. Rev. B 101, 081401(R) (2020), DOI: 10.1103/PhysRevB.101.081401.
2. Hoese et al., Sci. Adv. 6, eaba6038 (2020), DOI: 10.1126/sciadv.aba6038.
3. Koch et al., Commun. Mater. 5, 240 (2024), DOI: 10.1038/s43246-024-00686-y.
4. Stern et al., Nat. Mater. 23 (2024), DOI: 10.1038/s41563-024-01887-z.
5. White et al., Optica 8, 1153 (2021), DOI: 10.1364/OPTICA.431262.
6. Akbari et al., Nano Lett. 22 (2022), DOI: 10.1021/acs.nanolett.2c02163.
7. Mathur et al., Nat. Commun. 13 (2022), DOI: 10.1038/s41467-022-30772-z.
8. Whitefield et al., arXiv:2501.15341; subsequently Nature Materials (2026).
9. Gérard et al., Nat. Commun. (2026), DOI: 10.1038/s41467-026-68555-5.
