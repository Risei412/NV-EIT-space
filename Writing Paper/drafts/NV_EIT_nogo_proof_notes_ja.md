# 証明ノート — NV-EIT no-go 理論の仮定と定理・証明

**`NV_EIT_nogo_review_ja.md`（本編）の Methods 章（§2）を支える詳細ノート。学部生が独力で追えることを目標にした「講義ノート」形式。**

> 本編と同じく、リポジトリ内部の呼称（sector, SMRT class 等）は標準的な数学用語に置き換えている。対応表は本編の付録A、および本ノート末尾の付録Aを参照。
> 本ノートは3つの一次資料の要約・翻訳・平易化である：
> - `No-go theorem/Theorem and proofs/eit_nogo_lecture.tex`（学部生向け講義ノート、英語、厳密な証明つき）
> - `New no-go theory/src/core.py` とそのコメントが参照する `Sector-Master-Resolved-Theory` リポジトリの `theory/three_theorems_proofs.tex`
> - `EIT definition equivalence/prl_eit_equivalence_conditions.md` および `2g2e_package/docs/dark_state_free_eit_2g2e_report.md`
>
> **本ノートは要約であり、完全な証明の一次資料ではない。** 行間を厳密に埋めた証明が必要な場合は、必ず上記のtex/mdファイルを参照すること。特に本ノートの §3–5（暗状態・核・モーメント階層）は `eit_nogo_lecture.tex` の Section 2–5 の抄訳、§6（三分類）は SMRT の Theorem I–III の要旨、§7（EIT再定義）は equivalence 文書の要旨である。

---

## 0. このノートの読み方

理論は3段階で積み上がっている。それぞれの段階が答える問いが違う。

| 段階 | 主定理 | 答える問い |
|---|---|---|
| ① 暗状態の理論 | 定理1A・1B | 「暗状態は存在するか、それは散逸のもとで生き残るか」 |
| ② 核と厳密感受率、モーメント階層 | 核の定義、モーメント法、Krylov定理 | 「暗状態が壊れても、コヒーレント経路はどの速さ（べき）で消えるか」 |
| ③ 三分類（抑制次数） | 定理I・II・III | ①②を統一し、任意の系を3つの排他的クラスに分類する |
| ④ EITの再定義 | 透明化フロア公式、$2g{+}2e$ 存在定理 | 「暗状態がなくてもEIT的透明化と呼べる現象はあるか」 |

①②③はひとつながりの議論（§1–6）、④は独立した第2の問い（§7）である。読者は §1–6 を先に読み、必要なら §7 に進むとよい。

---

## 1. 仮定（Assumptions）— 理論が成り立つ範囲

以下の4条件（A1–A4）が、本ノート全体を通じて仮定される。本編§0の5条件と対応する（A3が本編の「弱プローブ」「利得なし」の両方を含む）。

- **(A1) 有限次元**：ヒルベルト空間、あるいは注目する応答が住む空間が有限次元である。
- **(A2) マルコフ的GKSL力学**：ゼロ次（制御場・静的な場・散逸を含む）のダイナミクスが、定常回転系における時間非依存のLindblad（GKSL）生成子 $\mathcal{L}_c$ で書ける。
- **(A3) 弱プローブ**：応答はプローブ振幅の1次まで計算する。制御場は非摂動的に（厳密に）扱う。
- **(A4) 一意な定常状態**：トレース0部分空間上で、群逆元（一般化逆）が well-defined である（縮退した定常多様体がない）。

**この範囲に入らないもの**：真に非マルコフな浴、強プローブによる飽和、伝搬が支配的な集団効果。これらは主定理の対象外であり、拡張が別途必要（`eit_nogo_lecture.tex` §1.3）。

---

## 2. 舞台設定：ヒルベルト空間の分解と光学コヒーレンスブロック

全ヒルベルト空間を

$$
\mathcal H = \mathcal H_e \oplus \mathcal H_g \oplus \mathcal H_r
$$

と分解する。$\mathcal H_e$ は注目する励起多重項、$\mathcal H_g$ は下準位（プローブ・制御が結合する側）、$\mathcal H_r$ はそれ以外（一重項、他の電荷状態、損失チャネルなど）。

プローブ応答が住むのは、密度行列の **非対角ブロック**

$$
\sigma = P_e\,\rho\,P_g \in \mathcal B(\mathcal H_g,\mathcal H_e)
$$

である（$P_e,P_g$ は $\mathcal H_e,\mathcal H_g$ への射影）。$\sigma$ を「光学コヒーレンスブロック」と呼ぶ。プローブに対する応答（吸収・分散）は、すべてこのブロックの運動から決まる。

### 補題（励起状態内のジャンプは光学コヒーレンスに対して純粋減衰として働く）

**主張**：励起多重項内に閉じたジャンプ演算子 $L_\mu = P_e L_\mu P_e$（たとえばフォノンによる軌道分枝間のホッピング）を考える。対応するLindblad散逸子

$$
\mathcal D[L_\mu](\rho) = L_\mu\rho L_\mu^\dagger - \tfrac12\{L_\mu^\dagger L_\mu,\rho\}
$$

は、光学コヒーレンスブロックに対して

$$
P_e\,\mathcal D[L_\mu](\rho)\,P_g = -\tfrac12\,L_\mu^\dagger L_\mu\,\sigma
$$

と作用する。すなわち **「補充項」（recycling term）が一切効かず、純粋な減衰としてのみ効く**。

**証明**：補充項 $L_\mu\rho L_\mu^\dagger$ を $e$–$g$ ブロックに射影すると

$$
P_e L_\mu\rho L_\mu^\dagger P_g = L_\mu(P_e\rho P_e)(P_e L_\mu^\dagger P_g) = 0
$$

となる。なぜなら $L_\mu^\dagger$ は $\mathcal H_e\to\mathcal H_e$ の写像なので $P_e L_\mu^\dagger P_g = 0$ だからである（$L_\mu$ が励起多重項に閉じているという仮定がここで効く）。反交換子の項は、$e$–$g$ ブロックに対しては左からの作用のみが残る：$P_e L_\mu^\dagger L_\mu \rho P_g = L_\mu^\dagger L_\mu\,\sigma$、一方 $P_e \rho L_\mu^\dagger L_\mu P_g = 0$（同じ理由）。$\blacksquare$

**物理的意味（なぜこれが重要か）**：たとえ「集団を移すだけ」のジャンプ（例：NVのフォノン誘起軌道混合 $|E_x\rangle\leftrightarrow|E_y\rangle$）であっても、それは光学コヒーレンスに対しては**純粋な損失**として効く。密度行列の対角ブロック（集団）では補充と流出が釣り合って定常状態が保たれるのに対し、非対角ブロック（コヒーレンス）ではそうならない。これが「フォノン力学が集団を奪わなくてもEITを殺せる」ことのミクロな起源であり、§5のモーメント展開の物理的な出発点である。

**単位の落とし穴（頻出）**：対称な双方向ホッピング $|X\rangle\leftrightarrow|Y\rangle$（各方向のレート $k$）を考えると、集団の不均衡は $\dot{(p_X-p_Y)} = -\Gamma_{XY}(p_X-p_Y)$、$\Gamma_{XY}=2k$ に従う。しかし上の補題により、各光学コヒーレンス $\sigma_{Xg}$ が受ける減衰は $k/2=\Gamma_{XY}/4$ である。「軌道緩和レート $\Gamma_{XY}$」という言葉が指すのは集団緩和レートであり、光学コヒーレンスの方程式に入るのは $\Gamma_{XY}/4$ である。この規約の取り違えは、定量予測を2〜4倍狂わせる。本リポジトリのコードは `phonon_rates.py` の冒頭コメントでこの規約を明示している（本編§2.7で引用した RATE CONVENTION の節）。

---

## 3. 暗状態：混同してはいけない2つの定理

### 定理1A（光学暗部分空間の次元）

プローブ・制御場は下準位から励起多重項への結合

$$
V_+ = P_e H_{\rm drive}P_g = \tfrac12\Omega,\qquad
\Omega = (\Omega_1 d_1,\dots,\Omega_{N_g}d_{N_g})
$$

を生む（$d_a$ は複素位相を保持した双極子ベクトル）。

**主張**：下準位のベクトル $|v\rangle\in\mathcal H_g$ が光学的に瞬時的に脱結合する（$V_+|v\rangle=0$）ことと $\Omega|v\rangle=0$ は同値。したがって光学暗部分空間は $\mathcal D_{\rm opt}=\ker\Omega$ であり

$$
\dim\mathcal D_{\rm opt} = N_g - \operatorname{rank}\Omega.
$$

**証明**：$V_+=\Omega/2$ なので同値性は自明。次元公式は階数・退化次数の定理そのもの。基底変換で $\Omega\mapsto U_e\Omega U_g^\dagger$ となるため、階数と退化次数は基底の取り方によらない。$\blacksquare$

**系**：2つの下準位（両方のRabi周波数が非零）に対して、非自明な暗ベクトルが存在するのは $\operatorname{rank}(d_1,d_2)<2$、すなわち（双極子が非零なら）$d_1$ と $d_2$ が $\mathcal H_e$ の中で**平行**であるとき、かつそのときに限る。

**物理的意味**：これがEIT干渉に対する第一の、純粋に幾何学的なフィルターである。干渉が起こるには、2つの光学経路が励起多重項の**同じ**重ね合わせに到達しなければならない。アルカリ原子で両方の経路が単一の超微細励起準位に結合する場合、$d_1,d_2$ は自明に平行（どちらもスカラー）である。NV$^-$やSiV$^-$では2つの経路が**異なる**励起軌道に到達しうるため、階数条件が非自明になる。

### 定理1B（Lindblad方程式の定常純粋暗状態）

$\dot\rho=\mathcal L(\rho)=-i[H,\rho]+\sum_\mu\big(L_\mu\rho L_\mu^\dagger-\tfrac12\{L_\mu^\dagger L_\mu,\rho\}\big)$ とする。

**主張**：規格化された $|D\rangle$ に対し、純粋状態 $\rho_D=|D\rangle\langle D|$ が定常であることと、ある $\lambda_\mu\in\mathbb C$、$h\in\mathbb R$ が存在して

$$
L_\mu|D\rangle=\lambda_\mu|D\rangle\ (\forall\mu),\qquad
\Big[H+\tfrac{i}{2}\sum_\mu(\lambda_\mu^*L_\mu-\lambda_\mu L_\mu^\dagger)\Big]|D\rangle = h|D\rangle
$$

が成り立つことは同値。

**証明の概略**：定常純粋状態は純度が一定に保たれる。直接計算により

$$
\frac{d}{dt}\operatorname{Tr}(\rho_D^2) = -2\sum_\mu\Big(\langle L_\mu^\dagger L_\mu\rangle_D - |\langle L_\mu\rangle_D|^2\Big)
$$

（分散の和、非負）。定常性はこの各項をゼロにすることを強制する。すなわち $|D\rangle$ は**すべての**ジャンプ演算子の固有ベクトルでなければならない：$L_\mu|D\rangle=\lambda_\mu|D\rangle$。これを代入し、$Q_D=I-|D\rangle\langle D|$ で射影すると、括弧内のエルミート演算子が $|D\rangle$ に平行に作用すること（実固有値 $h$）が従う。逆方向は、両条件を $\mathcal L(\rho_D)$ に代入すると非対角項がすべて消え、トレース保存とあわせて $\mathcal L(\rho_D)=0$ が得られる。$\blacksquare$

**EIT暗状態であるためにはさらに** $P_e|D\rangle=0$ かつ $\Omega|D\rangle=0$ を課す。自然放出のジャンプは自動的に $L_\mu|D\rangle=0$ を満たすが、基底状態のデフェージングや緩和、マイクロ波ノイズは一般に共通固有ベクトル条件を破る。

**重要な警告**：光学暗ベクトルが存在すること（$\ker\Omega\neq\{0\}$、定理1A）と、定常暗状態が存在すること（定理1B）は**別の問い**である。片方が成り立ってももう片方は成り立たないことがある。両方が揃って初めて、以下の応答解析と合わせて観測可能なEITが決まる。

---

## 4. 核 $K(z)$ と厳密な弱プローブ感受率

光学コヒーレンス $\sigma$ の非駆動時のダイナミクスは、補題（§2）とハミルトニアン交換子を合わせて

$$
\dot\sigma = -i(H_e\sigma-\sigma H_g) - \tfrac12 R_e\,\sigma + (\text{補正項}) + (\text{ソース項}),\qquad R_e=\sum_\mu L_\mu^\dagger L_\mu
$$

に従う。$\operatorname{vec}(AXB)=(B^{\mathsf T}\otimes A)\operatorname{vec}X$ でベクトル化すると、行列生成子 $A(z)$ が現れる。下側の添字が混ざらない（基底状態緩和がブロックに流入しない、相関ジャンプがない）限り、問題は励起多重項のレゾルベント方程式

$$
A_j(z)|x_j\rangle = |s_j\rangle,\qquad G(z)=A(z)^{-1}
$$

に帰着する。双極子行列 $C=(d_1,\dots,d_{N_g})$ を用いて

$$
K_{ab}(z) = d_a^\dagger\,G(z)\,d_b
$$

と定義する。対角成分 $S_a=K_{aa}$ は単一脚の光学自己エネルギー、非対角成分 $K_{ab}$ は励起多重項を介した**コヒーレントRaman頂点**である。これらはユニタリな基底変換に対して不変。生成子が正規行列とは限らないため、一般に $K_{21}\neq K_{12}^*$ であり、両方を計算する必要がある。

これが本編§2.2の $K_{12}=d_p^\dagger G d_c$、$K_{21}=d_c^\dagger G d_p$ に対応する。

### 比較応答の正しい定義（本編でいう $\chi_{\rm cut}$）

しばしば犯される誤りは、「EITが禁止される」ことを「全プローブ感受率がゼロになる」こと（$\chi_{\rm full}=0$）と定義してしまうことである。**これは正反対である**：完全なEITこそが $\chi_{\rm full}=0$ という姿をとる（§0の教科書公式で $\delta_2=0,\gamma_{21}\to0$ とすればわかる）。

正しい定義は、**注目部分空間を経由する寄与だけを除いた比較応答** $\chi_{\rm cut}^{(S)}$ を作ることである。除去の際、以下は固定したまま変えない：

- 直接の光学ブロック、
- プローブ源、
- 検出器、
- ゼロ次の定常状態、
- 結合していない散逸項すべて、
- 観測プロトコル。

すると、**主部分空間応答**

$$
R_S(\omega) = \chi_{\rm full}(\omega) - \chi_{\rm cut}^{(S)}(\omega)
$$

が定義され、厳密な分解 $\chi_{\rm full}=\chi_{\rm cut}^{(S)}+R_S$ が成り立つ。これは直接応答・背景応答と、コヒーレントな部分空間媒介応答とを分離する。**EITの問いは「$R_S$ が背景の吸収を破壊的に相殺するか」であって、「$\chi_{\rm full}$ がゼロかどうか」ではない。**

---

## 5. 速い散逸：モーメント階層と厳密なKrylovゼロ（三分類の土台）

### 設定

速い散逸レート $\Gamma$ を持つ励起多重項のレゾルベントを

$$
G(\Gamma,z) = [\Gamma D + A_0(z)]^{-1}
$$

と書く（$D$：散逸構造を表す行列、$A_0$：ハミルトニアンと遅い減衰）。$D$ が可逆なとき、Neumann級数展開

$$
G(\Gamma,z) = \frac{1}{\Gamma}D^{-1} - \frac{1}{\Gamma^2}D^{-1}A_0D^{-1} + \cdots
= \frac1\Gamma\sum_{k=0}^\infty\Big(-\frac{1}{\Gamma}A_0D^{-1}\Big)^k D^{-1}
$$

が成り立つ。これを用いて、プローブ脚 $p$・制御脚 $c$ に対する**モーメント**

$$
\mu_k = p^\dagger\,(-D^{-1}A_0)^k\,D^{-1}c
$$

を定義する。核は

$$
K(\Gamma,z) = p^\dagger G(\Gamma,z)c = \frac{1}{\Gamma}\sum_{k=0}^\infty \frac{\mu_k}{\Gamma^k}
$$

とモーメントで展開される。

### モーメント法（定理II の要旨）

**主張**：$\mu_0,\mu_1,\dots,\mu_{m-1}=0$ かつ $\mu_m\neq0$ なら、$K(\Gamma,z)\sim \mu_m\,\Gamma^{-(m+1)}$（$\Gamma\to\infty$）。すなわち**最初に非零になるモーメントの添字 $m$ が、抑制次数 $n_K=m+1$ を決める**。

これは展開式を見れば明らか（$\Gamma^{-1}$ の因子と $\Gamma^{-k}$ の因子を合わせると $\Gamma^{-(k+1)}$、最初に効くのは最小の $k=m$）。次数が整数になるのは、モーメントが行列の有限積のトレース（あるいは内積）として書かれる**代数的な量**であり、指数関数的減衰のような非整数のスケールを生まないからである。

### Krylov定理（定理I、厳密ゼロの判定：Class I）

**問い**：$\mu_k=0$ が**すべての** $k$ について成り立つ（応答が構造的に完全にゼロ）のはいつか。

**主張（Cayley–Hamilton的議論）**：$M=-D^{-1}A_0$ は有限次元行列（次元 $n$）であるから、Cayley–Hamilton の定理により $M$ は自身の特性多項式を満たす：$M^n = -\sum_{k=0}^{n-1}c_k M^k$。したがって、モーメント列 $\{\mu_k=p^\dagger M^k D^{-1}c\}_{k\ge0}$ は、**最初の $n$ 項** $\mu_0,\dots,\mu_{n-1}$ がすべてゼロなら、**それ以降のすべての** $\mu_k$（$k\ge n$）もゼロになる。なぜなら $\mu_n = p^\dagger M^n D^{-1}c = -\sum_{k<n}c_k\mu_k = 0$ であり、帰納的に以降も同様だからである。

**したがって、$\mu_0,\dots,\mu_{n-1}$ という有限個（行列の次元と同数）のモーメントを機械精度で計算するだけで、「全次数で厳密にゼロ」を判定できる。** これは漸近的な推測ではなく**有限の証明書（certificate）**である。本リポジトリの `New no-go theory/src/core.py` の `krylov_class1_certificate` がこれを実装しており、`max|mu_k| < 10^{-10}` を厳密ゼロの判定基準としている。

### 保護されたチャネル（$D$ が特異な場合）

$D$ が可逆でない場合（一部のモードが速い散逸を受けない）、Riesz射影 $P$（$\ker D$ への射影）を用いて系を分解する。$D$ が半単純（$D$ の像と核が互いに補空間をなす）なら、保護部分空間上の有効核は

$$
F_0 = p^\dagger P\,B_P^{-1}\,P c
$$

という**$\Gamma$ に依存しない**（$\Gamma\to\infty$ で有限値に収束する）量になる（定理III）。これが本編でいう「保護されたgo」（第3分類、$n=0$）である。

---

## 6. 三分類（抑制次数によるトリコトミー）

§5の結果をまとめると、任意の有限次元マルコフ弱プローブ系に対して、部分空間応答の抑制次数 $n$ は

$$
n \in \{\infty\}\ \cup\ (0,\infty)\ \cup\ \{0\}
$$

のいずれかに**排他的に**分類される：

| 分類 | 条件 | 判定法 |
|---|---|---|
| 厳密構造的 no-go（$n=\infty$） | $\mu_k=0$ が全 $k$ で成立 | Krylov証明書（有限個のモーメント計算） |
| 漸近的 no-go（$0<n<\infty$） | 最初の非零モーメントが $\mu_{n-1}$ | モーメント法、または対数対数フィット |
| 保護された go（$n=0$） | $D$ が特異かつ保護部分空間で $F_0\neq0$ | Riesz射影 + 定理III |

**この三分類が材料非依存である理由**：上記の議論のどこにも、系が何の物質か（ダイヤモンド、超伝導回路、アルカリ原子）は現れない。現れるのは行列 $D,A_0,p,c$ の代数的構造だけである。同じ整数 $n$ が異なる材料に現れうる（本編§3.4のGate C参照）のは、この意味においてである。ただし、モーメントの**値**（$\mu_m$ そのもの、あるいは $F_0$）は材料固有のパラメータに依存する。**次数は普遍、係数は材料依存**、というのが正確な主張である。

### 核から観測量への次数継承（本編§2.3の根拠）

物理的に測るのは核 $K_{12}$ 単体ではなく、

$$
R_{\rm obs}(\Gamma,z) = -\frac{\beta K_{12}K_{21}}{g_{\rm eff}+\beta S_2}
$$

という有理式である（$\beta=\Omega_c^2/4$）。$K_{12}\sim\Gamma^{-n_{12}}$、$K_{21}\sim\Gamma^{-n_{21}}$、分母が $g_{\rm eff}\neq0$（$\Gamma$ に依らず有限）なら $R_{\rm obs}\sim\Gamma^{-(n_{12}+n_{21})}$。分母自身が $\Gamma$ とともに消える場合（$g_{\rm eff}=0$、$S_2\sim\Gamma^{-n_{S_2}}$）は $R_{\rm obs}\sim\Gamma^{-(n_{12}+n_{21}-n_{S_2})}$。これが継承則 $n_{\rm obs}=n_{12}+n_{21}-n_{\rm den}$ であり、単なる代数（分数のべきの引き算）であって新たな物理的仮定を要しない。

---

## 7. EITの再定義とその等価条件

ここからは独立した第2の理論であり、上記①②③の枠組みの**内部**で「EITとは何か」という定義そのものを精密化する。

### 7.1 主張の形

比較応答 $\chi_{\rm cut}^{(S)}$（§4）と、実測される吸収に対応する実線形汎関数 $\mathcal A[\chi]$（規約により $\operatorname{Im}\chi$ か $\operatorname{Re}\chi$）を用いて、

$$
\boxed{\ \text{EIT} \iff \mathcal A[R_S](\omega_0) < 0\ \text{（正則な観測周波数 }\omega_0\text{ において）}\ }
$$

すなわち「部分空間を介した応答 $R_S$ が、比較応答が持つ吸収を**減らす**方向に効く」ことをEITの定義とする。これは暗状態やCPTの存在を仮定しない。

### 7.2 標準的な$\Lambda$系での閉じた形（透明化フロア公式）

暗状態を持つ標準的な $\Lambda$ 系（プローブと制御が同一の励起状態 $e$ に結合）を最小の $2\times2$ 有効モデル

$$
A(\delta) = \begin{pmatrix}a_1 & i\Omega_c/2\\ i\Omega_c/2 & g\end{pmatrix},\qquad a_1=\gamma_e-i\delta,\quad g=\gamma_g-i\delta
$$

で考える（$\gamma_e$：励起状態の減衰、$\gamma_g$：基底コヒーレンスの減衰）。二光子共鳴 $\delta=0$ での応答を直接計算すると

$$
\chi_{\rm full}(0) = \frac{\gamma_g}{\gamma_e\gamma_g+\beta},\qquad
\chi_{\rm cut}(0) = \frac{1}{\gamma_e},\qquad
\beta=\Omega_c^2/4
$$

（比較応答 $\chi_{\rm cut}$ は $2\times2$ 行列の $(1,1)$ 成分を、非対角結合＝制御場を切って求めたものに一致する）。したがって透明化コントラストは**厳密に閉じた形**で

$$
\boxed{\ C_S(\gamma_g,\Omega_c) = 1-\frac{\chi_{\rm full}(0)}{\chi_{\rm cut}(0)} = \frac{\beta}{\gamma_e\gamma_g+\beta} = \frac{\Omega_c^2}{4\gamma_e\gamma_g+\Omega_c^2}\ }
$$

**証明**：$2\times2$ 行列の逆行列公式 $A^{-1}=\frac{1}{\det A}\begin{pmatrix}g&-i\Omega_c/2\\-i\Omega_c/2&a_1\end{pmatrix}$ から $(A^{-1})_{11}=g/\det A=g/(a_1g+\beta)$。$\delta=0$ を代入して $a_1=\gamma_e,\ g=\gamma_g$ とすれば $\chi_{\rm full}(0)=(A^{-1})_{11}|_{\delta=0}=\gamma_g/(\gamma_e\gamma_g+\beta)$。制御場オフ（$\Omega_c=0\Rightarrow\beta=0$）での同じ成分は $1/\gamma_e$。あとは代入して整理するだけ。$\blacksquare$

**帰結**：$\gamma_g\to0$（理想的な基底コヒーレンス）で $C_S\to1$（完全透明化）。$\gamma_g>0$ では、どれだけ $\Omega_c$ を強くしても $C_S<1$ の**透明化フロア**が残る（$\Omega_c\to\infty$ で $C_S\to1$ に漸近するが到達しない）。これが本編§3.2で引用した閉じた形の公式である。

### 7.3 暗状態を持たない透明化の存在証明（$2g+2e$ 最小モデル）

**問い**：暗状態が原理的に存在しない系（2つの下準位、2つの励起状態、プローブと制御が異なる励起状態に結合し双極子行列 $\Omega=(\Omega_p d_p,\Omega_c d_c)$ が $\operatorname{rank}\,2$）でも、EIT的な透明化（$C_S>0$）は起こりうるか。

**構成**：基底 $|g_1\rangle,|g_2\rangle,|e_1\rangle,|e_2\rangle$、双極子 $d_p=(1,0)^{\mathsf T}$、$d_c=(\cos\theta,\sin\theta)^{\mathsf T}$（$\theta\neq0$ なら $d_p,d_c$ は非平行 → 定理1Aにより $\ker\Omega=\{0\}$、光学暗状態は存在しない）。この系の定常GKSL方程式を厳密に解き、$\chi_{\rm full}$ と $\chi_{\rm cut}$ を比較する。

**数値結果**：代表的なパラメータで、比較応答に対して吸収が **79.03%** 低減する。透明窓近傍にAutler–Townes型の極分裂も現れない（定常状態が正則にとどまる）。

**しかし完全透明化は不可能**：行列式の恒等式により、正則な理想 $2g+2e$ モデル（基底コヒーレンス減衰ゼロ、$\operatorname{rank}\,\Omega=2$）では、厳密な伝達零点（$\chi_{\rm full}=0$ となる実周波数）は複素平面の実軸上には存在しないことが証明される（`2g2e_package` の行列式解析）。すなわち：

$$
\boxed{\text{暗状態なしの「EIT的コヒーレント透明化」は存在するが、暗状態なしの「完全EIT」は存在しない}}
$$

これは §7.1 の再定義（$\mathcal A[R_S]<0$ という不等式条件）が、$\chi_{\rm full}=0$ という等式条件より真に緩い、実質的な一般化であることを示す最小の反例つき存在証明である。

---

## 8. 具体例で確認する（本編Q3への直接の根拠）

`eit_nogo_lecture.tex` §7 の3つの具体例を要約する。詳細な数値は同文書・本編Q3を参照。

### アルカリ原子（Rb）— 教科書的な go

- 定理1A：プローブ・制御は単一の超微細励起準位に結合。$d_1,d_2$ は自明に平行（両方スカラー）→ $\ker\Omega\neq\{0\}$、光学暗状態が存在。
- 定理1B：自然放出のジャンプは $L_\mu|D\rangle=0$ を自動的に満たす。バッファガス・コーティングセルにより基底状態緩和 $\gamma_g$ を極小化できるため、定常暗状態がほぼ厳密に生き残る。
- §5の言葉で言えば、抑制次数 $n$ は小さく（$M_0\neq0$、$n=1$ に相当する構造）、かつ $\Gamma$（自然放出＋ドップラー幅、室温でもMHz–GHzオーダー）が内部分裂より小さいまま保たれる。

### NV$^-$ 中心 — 教訓的な no-go

- 定理1A：プローブ・制御が異なる軌道分枝 $X,Y$ に結合する配置では、基底状態の直交性から $d_p^\dagger d_c=0$（$M_0=0$、本編§2.4）。厳密な選択則であり、任意の磁場・歪みでも破れない。
- §5–6の言葉で言えば、これは第2分類（漸近的no-go、$n_R=4$）に対応する。フォノン誘起の軌道混合 $\Gamma_{XY}(T)\propto T^5$ が速い散逸を担い、室温では励起多重項の内部分裂を4桁近く上回る（本編§3.2の $x(300\text{K})=7559$）。

### 群IV中心（SiV$^-$, SnV$^-$）— 軌道go

- プローブ・制御が同一軌道分枝（同一スピン）に結合しうる配置では $M_0\neq0$（$n=1$）。大きなスピン軌道分裂（SnV$^-$: $\Delta_e=3000$ GHz）により、$\Gamma/\delta$ を小さく保ちやすい。ただし単一フォノンBose則のため温度で $\Gamma$ を大きく動かせる範囲は狭い（本編§3.4(e)）。

**Q3への対応**：これら3系統の違いは、定理1A・1Bのレベルでは「暗状態が構造的に存在するかどうか」という**質的**な違いに見えるが、§5–6の枠組みで見ると、いずれも同じ三分類の中の点であり、実際に効いているのは（i）抑制次数 $n$（アルカリは低次、NVは高次）と（ii）$\Gamma/\delta$ の比（アルカリは室温でも $\ll1$、NVは室温で $\gg1$）の**2つの連続的パラメータ**である。本編§3.3の結論（「構造的な差ではなく桁数の差」）は、この統一的な描像から導かれる。

---

## 付録A：本ノートと一次資料の用語対応

| 一次資料（英語） | 本ノート・本編での訳語 |
|---|---|
| optical coherence block $\sigma$ | 光学コヒーレンスブロック |
| optical dark subspace $\mathcal D_{\rm opt}$ | 光学暗部分空間 |
| stationary Lindblad dark state | 定常Lindblad暗状態 |
| sector $S$ / sector-cut $\chi_{\rm cut}^{(S)}$ | 注目部分空間 $S$ / 比較応答 $\chi_{\rm cut}$ |
| sector-resolved correction $R_S$ | 部分空間応答 $R_S$ |
| moment $\mu_k$ | モーメント $\mu_k$ |
| Krylov / Cayley–Hamilton exact-zero certificate | Krylov証明書（厳密ゼロの判定） |
| suppression index $\nu$ | 抑制次数 $n$ |
| Class I / II / III | 厳密構造的no-go／漸近的no-go／保護されたgo |
| symmetry-protected transfer zero | 対称性保護された伝達零点 |
| transparency floor | 透明化フロア |

## 付録B：一次資料ファイル一覧

| 内容 | ファイル |
|---|---|
| 暗状態定理1A・1B、核とモーメント階層、Krylov定理の完全な証明 | `No-go theorem/Theorem and proofs/eit_nogo_lecture.tex`（学部生向け、本ノートの主要典拠） |
| 同内容のより簡潔な数式主体版 | `No-go theorem/Theorem and proofs/eit_nogo_proofs.tex` |
| フル版（v6.2）、応用・数値検証を含む完全版 | `No-go theorem/Theorem and proofs/EIT_no_go_go_theory_v6_2_English.tex`（4641行） |
| 三分類（定理I–III）の数値インフラとコメント中の理論要旨 | `New no-go theory/src/core.py` |
| 三分類の完全な証明（`three_theorems_proofs.tex`） | 別リポジトリ [Sector-Master-Resolved-Theory](https://github.com/Risei412/Sector-Master-Resolved-Theory) `theory/` |
| EIT再定義の枠組みと透明化フロア公式の導出 | `EIT definition equivalence/prl_eit_equivalence_conditions.md`、`results/F4_findings.md` |
| $2g+2e$ 存在定理と行列式恒等式 | `EIT definition equivalence/2g2e_package/docs/dark_state_free_eit_2g2e_report.md` |
| 具体例（アルカリ・NV・群IV）の詳細 | `No-go theorem/Theorem and proofs/eit_nogo_lecture.tex` §7、`GateC_material_independence/README.md` |

## 付録C：演習（自己確認用）

`eit_nogo_lecture.tex` §11 に学部生向け演習問題がある（レート規約の変換、階数条件の具体計算、Krylov証明書の手計算など）。本編・本ノートの内容を独力で追えたかの確認に使える。
