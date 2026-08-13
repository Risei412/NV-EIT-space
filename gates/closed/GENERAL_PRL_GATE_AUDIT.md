# 一般理論PRL（経路モーメント応答則）— Gate A–E 監査 (2026-08-02)

> **この監査の対象**：`Writing Paper/NV_EIT_PRA_PRL_Split_Strategy_20260724.md`
> §13（Priority 1–8）と §14（Gate A–D）が定める**一般理論PRL**、実装は
> `New no-go theory/`。
>
> **対象外**：NV中心のEITの温度限界・室温no-go（PRAライン）。そちらは
> `SIMULATION_PLAN.md` の Gate 1–5 が仕様で、監査は `PRL_CLAIM_GAP_AUDIT.md` にある。
>
> **「Fig.1–4」という名前は2つの論文で衝突している。**
> `Writing Paper/prl_figures/` の Fig.1–4 は本監査が扱う一般理論PRLの図。
> `SIMULATION_PLAN.md` §6 が要求する NV-EIT の Fig.1–4 は依然0枚
> （`PRL_CLAIM_GAP_AUDIT.md` の Blocker B6）。

各ゲートの `overall_pass` は、**そのスクリプト自身が検査した項目**に基づく。本監査は
そのフラグが検査していなかった項目を洗い出したものである。Gate A–D は監査前から
全て `overall_pass: true` だったが、以下の3件は判定が主張を支えていなかった。

## 総括

| 項目 | 監査前 | 是正後 |
|---|---|---|
| Gate A テスト | 6 | 10 |
| Gate B テスト | 5 | 7 |
| Gate C テスト | 5 | 9 |
| Gate D テスト | 4 | 10 |
| 図テスト | 4 | 11 |
| クリーン clone で図がビルド可能 | 不可（`ModuleNotFoundError: sympy`） | 可 |

---

## Blocker（是正済み）

### B1. Gate A — 証明書が無いと合否が真空的に真になっていた
`gate_a_observable.py:302` が

```python
pred_cert_agree = (nu_cert is None) or (nu_pred == nu_cert)
```

だったため、symbolic pencil を持たない spec は証明書判定を**無条件に満たし**、
`pred_cert_agree=True` と表示された。物理NVの spec がまさにそれで、
`gates_summary_gateO.json` の NV 行は `nu_obs_cert: null` のまま PASS していた。

**是正**：RoomT step3 の厳密 pencil を `symbolic_H` / `symbolic_pencil` として共有化し
（step3 の出力は bit-identical）、NV spec に接続。`cert_present` を常時出力し、
`pred_cert_agree` は比較対象が無ければ `None`。G-O2 は NV 自身の証明書が 4 であることと
step3 ルートとの一致を要求、G-O5 は全 bilinear モデルが証明書を持つことを要求。
NV は Gate A 自身の経路で `nu_obs = 4` を証明するようになった。

**数値↔厳密の橋渡し**も追加した。NV の `dp`/`dc` は `nv_model.legs()` の固有分解由来の
float で、step3 の厳密 pencil と結びつける保証が無かった。位相を除いた脚の一致と
pencil の数値一致（8.9e-16）を検査する。なお step3 は歴史的に軌道ラベルを入れ替えた
鏡像ペア（e₁, e₃）を証明していたが、この H では両者は等価
（同じ `M0 = 0`、同じ生成項 `-√2·i·Lperp/2`、同じ次数）なので step3 の結論は有効。
Gate A は数値側と同じペア（e₄, e₀）を使う。

### B2. Gate D — 実験識別可能性の入力3つが仮定値だった
1. **代表コントラストが `1e-2` のハードコード**。これは PRA 側 Gate 5 の
   **単一欠陥値**（0.0136）に相当し、アンサンブル平均前の値。同じ表の平均後の値は
   1〜2桁小さい。→ 5シナリオを表から読み、Gate 5 自身の合格基準が指す
   `post_selected_shimmed`（2.63e-3）で判定し、残りを併記。
2. **検出限界が 50 K で計算されていた**（コントラストは 70 K）。→ 70 K に統一。
   副次的な発見として、50 K では sector OD が 19 で試料がほぼ不透明になり
   `c_min` が押し上げられていた（70 K では OD 7.9、`c_min = 1.16e-6`）。
3. **超伝導のダイナミックレンジが `9.0` のハードコード**。→ 掃引範囲と
   engineerable 範囲（κ/2π ≈ 0.1–50 MHz、2.7 decade）に分離し、判定には後者を使う。

さらに `required_gamma_range` は**ノイズ無しの解析カーネル**にフィットしており、
1 decade で `|slope−2| < 0.1` を 2e-5 で満たしていた（浮動小数点の性質であって
実験の性質ではない）。`experimental_budget.py` に測定モデルを分離し、
`|bias| + 2·std < 0.5`（隣接クラスを2σで分離）を基準にした。

### B3. Gate C — クラス3が単一材料にしか無かった
クラス1・2 は diamond と非diamond の両方で実現されていたが、クラス3は NV のみ。
中心主張が最も依存するクラスが1材料に依存していた。G-C4 はクラス1・2 しか見ていなかった
ためこれを検出できなかった。

**是正**：`chain3_witness.py`（3モード鎖、`M0 = M1 = 0`、`M2 = -J12·J23`、
振幅次数 3.000、固定読み出し population 次数 6.000、reduced==full 8.8e-16）を追加。
G-C4 は全クラスに両ホストを要求、新設 G-C6 が witness の厳密性を検査。

---

## 是正の結果として確定した主張の範囲

### Gate E — 漸近到達性・工学的到達性・検出性の結合監査（2026-08-02）

Gate D は可変幅とノイズ下の指数識別を検査したが、その可変区間が同時に
漸近域かつ検出域に入ることまでは要求していなかった。Gate E では同一の散逸区間

\[
I_E=I_{\rm engineerable}\cap I_{\rm asymptotic}\cap I_{\rm detectable}
\]

が存在し、局所指数・信号床・必要区間幅を同時に満たすことを要求した。

| プラットフォーム | Gate E | 根拠 |
|---|---|---|
| 光学NV | **FAIL** | 漸近域へ入る前に信号が検出床へ落ちる（0.51–0.84 decade） |
| Gate B SC transfer | **FAIL** | 0.1–50 MHzでは実効指数が約0、整数漸近域は8桁以上先 |
| 3モード鎖 | **MODEL-LEVEL CONDITIONAL PASS** | 15–45 MHzの有限散逸窓でclass 2/3を統計的に分離 |

3モード鎖では \(J_{12}/2\pi=3.0\,\mathrm{MHz}\)、
\(J_{23}/2\pi=2.4\,\mathrm{MHz}\)、\(\Gamma/2\pi=15\text{–}45\,\mathrm{MHz}\)
を用いた。class 3 のglobal exponentは2.9002、local exponentは
2.8163–2.9518、Monte Carlo 95%区間は[2.8286, 2.9779]。class 2の95%区間
[2.0680, 2.1090]とは非重複である。最悪点は−189.67 dBmで、3 K、入力−110 dBm、
40 s/点ならSNR=10.21（SNR 10に必要な時間38.36 s）。製造ばらつき5000標本の
合格率は100%だった。

背景除去には結合符号の8状態Thue–Morse cycle
`+ - - + - + + -`を選ぶ。定数・一次・二次時間ドリフトを相殺し、信号の1000倍の
滑らかな背景でも指数回収率100%。80 ms cycleに対し必要な背景ドリフト時間尺度は
2.53 s以上である。

このPASSはデバイス固有の実証ではない。15–45 MHzの可変engineered loss、結合を
保った独立損失掃引、−190 dBm級の増幅・校正、符号切替transient、feedthroughの
実測ドリフトスペクトルは設計仮定として残る。したがって許される主張は
「3モードengineered-loss chainにclass 2/3を識別できる有限窓が存在する」までであり、
任意のengineered-dissipation platformでの測定可能性ではない。

**光学NVでは指数が読めない。** これは新設した G-D7（slope budget）が出した結論であり、
本監査で最も重要な発見である。

| プラットフォーム | ν | 使用可能窓 | 必要 | 解像可能 |
|---|---|---|---|---|
| NV optical EIT (single) | 4 | 1.02 dec | 1 | ○（ただし平均前の値） |
| NV optical EIT (post_selected_shimmed) | 4 | 0.84 dec | 1 | **×** |
| NV optical EIT (high_density) | 4 | 0.51 dec | 1 | **×** |
| SC transfer (generic) | 1 | 2.70 dec | 1 | **×（Gate E：漸近域外）** |
| SC transfer (protected) | 2 | 2.70 dec | 1 | **×（Gate E：漸近域外）** |
| 3-mode chain (class 3) | 3 | 0.48 dec (15–45 MHz) | finite-window fit | **△（model-level）** |

ν = 4 では Γ を1桁上げるとコントラストが4桁落ちるため、検出限界までの窓が
フィットの必要幅に届かない。G-D7 は「≥1プラットフォーム」定義なので Gate D 自体は
PASS を維持するが、Gate Eによって一般的な**「測定可能な設計則」**は棄却された。
測定可能性は3モードengineered-loss chainの具体的有限窓に限定し、NV は exact class の
認証を与える witness と位置づける。
戦略文書 §14 に Gate D の fallback 節を新設してこれを明文化した（従来 Gate D だけ
不合格分岐が未定義だった）。

**Gate B の実験窓は未解決。** 現実的な κ 域では指数が generic / protected とも ~0 で、
漸近的整数次数が読めるのは物理的なバス減衰より8桁以上高い κ のみ。Gate B の主張は
構造的普遍性であり、この素子の測定手順ではない。G-B6 がこれを機械可読に記録する。

---

## 残る既知のギャップ（未着手）

- **group-IV モデルの模式性**。phonon normalization と dipole geometry は依然として
  模式的（単一レートスケール、直交軌道基底の脚、`|M0|` の正規化）。定数の出所は
  `GateC_material_independence/SiV_SnV_phonon_AIC_parameters.md` に記録したが、
  値そのものは変えていない。整数指数には影響しないが、group-IV 曲線を絶対応答の
  定量予測として読んではいけない。
- **PRL 原稿が存在しない**。`Writing Paper/drafts/` にあるのは日本語の講義ノート1本で、
  経路モーメント関連の語彙は皆無。図4枚と数値認証は揃っているが本文は0行。
  REVTeX テンプレートも Supplement の構成表も無い。
- **投稿凍結の項目が全て未達**：git tag、GitHub release、Zenodo DOI、
  原稿への commit hash 記録。

## ブランチ整理の実施結果

統合したもの（すべて本作業ブランチ `claude/prl-claim-calculation-check-vvrvx8` へ）：

| ブランチ | 内容 |
|---|---|
| `agent/prl-gates-a-d-bundle-2026-08-01`（PR #20） | Gate A–D 独立監査バンドル（9ファイル） |
| `claude/theorem-proof-code-audit-dt2ur4` | `AUDIT_theorem_and_code_2026-08-01.md` |
| `claude/notion-question-db-audit-u05gz6` | Gate F5-C の結果とランナー |

**削除できなかったもの（環境の制約）**。main の祖先であることを
`git merge-base --is-ancestor` で14本すべて確認した上で削除を試みたが、
`git push origin --delete` はいずれも **HTTP 403** で拒否された。同じ制約で
**タグの push も拒否**される（`pre-prl-cleanup-2026-08-01` はローカルにのみ存在）。
`/root/.ccr/README.md` は 403/407 を組織ポリシーによる拒否として扱い再試行しないよう
指示しているため、削除は未実施のまま残している。リモートのブランチ数は 23 のまま。

削除対象として検証済みの14本（すべて `claude/` 接頭辞）：
`eit-no-go-numerics-l2v9i6`, `theorem-proof-plan-qnouqs`,
`no-go-theorem-computation-o11kzl`, `no-go-theory-phenomena-3qqsjz`,
`no-go-theory-tex-guide-ucx404`, `numerical-computation-plan-2y83n1`,
`eit-definition-equivalence-7pcftx`, `theorem-revision-strategy-cx2yqf`,
`claude-md-documentation-jgbjor`, `nbv-eit-no-go-access-7gdyjq`,
`understand-commit-text-m4lqag`, `repository-content-review-vbly8p`,
`repository-organization-o6k79k`, `eit-no-go-theory-publication-unvqb4`。
リポジトリ管理者が GitHub 側から削除すれば済む。

**意図的に触れていないもの**：
- `claude/content-review-9afb78`（7コミット・約1900行。`Theorem and proofs/` に対する
  modify/delete 衝突を必ず起こす。SMRT 側の管轄を含むため選択的 cherry-pick が必要）
- `claude/repository-content-review-yuwiox`（対象ファイルが SMRT へ移管済み。ただし
  PR #16 マージ後に2コミット追加されているので無条件削除は不可）
- `claude/simulation-plan-markdown-pug2sy`, `claude/paper-1-numerical-methods-4lekau`
  （計画文書のみ。基点が 7/16–7/17 と古く、SMRT 移管前のパス参照を含む可能性がある）
