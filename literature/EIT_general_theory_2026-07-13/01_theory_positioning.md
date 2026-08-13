# 理論の位置づけと推奨構成

## 1. 現在の2つの理論PDFから抽出される核

第一稿は、二つのlower statesと複数excited statesをcoupling vectors `p,c`で表し、弱probe感受率がweighted overlap `p^dagger D^{-1}c`で制御されること、merged limitで `p^dagger c`へ縮退することを主張している。

改良稿はこれをbasis-invariant kernel

`K(omega,T) = p^dagger G_e(omega,T)c`

へ拡張し、次を区別している。

- exact dark state: `p || c`
- imperfect but finite EIT: `K != 0` and `p` not parallel to `c`
- no-go: `K = 0`

さらに、complete excited manifoldのclosureにより

`p^dagger c = <g1| d_p^dagger Pi_E d_c |g2>`

とし、point-group symmetryとspin-scalar electric dipoleからNV spin-Lambdaのnullを導いている。

## 2. 先行研究との境界

### 既知であり、主たる新規性にしてはいけない部分

- coupling matrixのrank/nullspaceによるdark-state存在判定
- bright/dark subspaceへのMorris-Shore reduction
- SVDによるdark-state抽出
- 任意multilevel coupling graphのdark-state分類
- Lindblad系のstationary dark-state分類
- multiple excited levelsとinhomogeneous broadeningによるEIT低下
- EITとAutler-Townes splittingの区別

### 新規性候補

1. **Thermal scalarization theorem**  
   非Hermitian excited-manifold resolventが、制御された極限でscalar operatorへ近づくことを定式化する。

2. **Symmetry collapse theorem**  
   scalarization後に干渉kernelが `d_p^dagger Pi_E d_c` のlower-manifold matrix elementへ縮退し、point-group selection ruleでzero/nonzeroが決まることを示す。

3. **Finite-EIT theorem**  
   exact dark stateの存在ではなく、observable susceptibilityのcontrol-induced termが消える条件を定理化する。

4. **Material taxonomy with falsification**  
   NV exact no-go、SiC/rare-earth empirical go、group-IV predicted/observed goを、同じkernelで再現する。

## 3. 推奨する数学的定式化

### 3.1 Hilbert-space coupling map

Lower subspace `L`、excited subspace `E`を定義し、光学結合を写像

`V: L -> E`

として扱う。二lower-stateの場合、

`V = Omega_p |p><g1| + Omega_c |c><g2|`.

Exact dark stateは `ker(V)` の非自明性で決まり、これは先行研究に帰属させる。

### 3.2 Non-Hermitian resolvent

`G_e(omega,T) = [omega I_E - H_e(T) - Sigma_e(omega,T)]^{-1}`

とする。`Sigma_e`にはradiative decay、pure dephasing、phonon-induced mixingを含める。単純な `i Gamma_e/2` はMarkov近似として明示する。

### 3.3 Controlled merged-manifold expansion

`M = omega I - H_e - Sigma_e = z I - delta M`, `Tr(delta M)=0` と分解し、

`G_e = z^{-1} I + z^{-2} delta M + O(||delta M/z||^2)`.

したがって

`K = z^{-1} p^dagger c + z^{-2} p^dagger delta M c + ...`.

これにより次を厳密に分離できる。

- strict merged no-go: `p^dagger c=0` かつ補正も対称性で禁止
- asymptotic collapse: leading termはzeroだが有限温度補正が残る
- selected-submanifold escape: `||delta M/z||`を意図的に大きくする

### 3.4 Symmetry condition

`O_pc = P_L d_p^dagger Pi_E d_c P_L`.

`<g1|O_pc|g2>` が非零となる必要十分条件を、次の積表現とreduced matrix elementで書く。

`Gamma(g1)^* tensor Gamma(O_pc) tensor Gamma(g2)` がtotally symmetric irrepを含む。

さらにelectric dipoleがspin-scalarなら `⟨s1|s2⟩ != 0` が必要。orbital degeneracyだけでは十分でない。

### 3.5 Liouville-space bridge

一般の線形応答は

`chi_pc(omega) = Tr[d_p (i omega - L_0)^{-1} V_c rho_0]`.

`p^dagger G_e c`への縮約が成立する条件を列挙する。

- weak probe
- populationがlower reference stateに固定
- jump-induced optical coherencesが無視可能
- control fieldによるrank-one self-energy reductionが有効
- propagation/optical-depth effectsを単一emitter susceptibilityから分離

これによりHilbert-space kernelとLindblad theoryを競合ではなく階層化できる。

## 4. 定理の推奨名称

避ける:
- General dark-state theory for arbitrary multilevel systems
- General existence condition of dark states

推奨:
- Symmetry criterion for EIT in thermally merged excited-state manifolds
- Thermal-manifold no-go/go theorem for coherent optical interference
- Basis-invariant response-kernel criterion for multilevel solid-state EIT

## 5. 現在の主張で修正すべき点

1. `K != 0`をEITの十分条件としない。必要条件、またはcontrol-induced interference amplitudeの非零条件とする。
2. `p^dagger c=0 => any control power cannot produce EIT`は、対象channel、電気双極子、complete manifold、単一control mode、直接ground couplingなし、という仮定を明記する。
3. `orbital-singlet => spin-Lambda no-go`は有力だが、spin-dependent dipole、magnetic dipole、spin-phonon-assisted optical pathsの上限を与える。
4. `orbital-degenerate => go`を削除し、full tensor selection ruleへ置換する。
5. `eta_geom`、EIT contrast、dark-state fidelity、kernel magnitudeを別定義にする。
6. SiV cross-polarizationの `eta_geom≈1` は、dipole phase convention、vibronic basis、remote excited statesを含む独立計算で再確認する。
7. EITとATSの判別を数値図に組み込む。

## 6. 投稿可能な中心命題

> In a thermally merged excited-state manifold, the control-induced interference kernel reduces, under a controlled scalar-resolvent limit, to a symmetry-selected lower-manifold tensor matrix element. This establishes a material-independent no-go/go criterion that is distinct from the known rank condition for exact dark states and predicts the resolved-to-merged collapse or survival of EIT across solid-state emitters.

この中心命題を守るには、NVだけでなく、少なくとも以下が必要。

- exact/no-go example: NV spin-Lambda
- empirical go example: SiC divacancyまたはrare-earth crystal
- constructive predicted go: group-IV same-spin orbital-Lambda
- correction bound: non-scalar resolvent、branch-dependent linewidth、vibronic leakage
- Liouvillian benchmark: full master equationとの一致領域
