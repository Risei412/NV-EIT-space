# Step 4:README・数値優先度文書・QFI 命名の改訂案

Step 1–3 と同じ方針で、`New no-go theory/README.md` および
`new_nogo_numerical_priorities.md` は直接編集せず、ここに
**適用可能な差分(コピー&ペースト用の置換テキスト)** としてまとめます。
QFI 命名変更については Gate U8 の数値検証(negative control)も実行しました。

---

## 4.1 `New no-go theory/README.md` の一文説明改訂

対象:冒頭段落(3–7 行目)。

### 現行

```
This directory contains the new no-go theory that upgrades the binary EIT
go/no-go criterion (see `../No-go theorem/`) to a **three-tier classification**
of sector-mediated coherent response in finite-dimensional Markovian
weak-probe systems, following the strategy in
`EIT_new_nogo_PRX_taxonomy_and_quantum_network_strategy.md`.
```

### 置換後(改訂案 §9.3 の一文を統合)

```
This directory contains the new no-go theory that upgrades the binary EIT
go/no-go criterion (see `../No-go theorem/`) to a **three-tier classification**
of sector-mediated coherent response in finite-dimensional Markovian
weak-probe systems, following the strategy in
`EIT_new_nogo_PRX_taxonomy_and_quantum_network_strategy.md`.

SMRT classifies how an operationally specified internal sector contributes
to an observable response under strong native dissipation. A sector is
defined by an equivalence class of admissible selective interventions
sharing the same retained Riesz projector (see `Sector/`), and the ideal
sector cut is their universal strong-intervention limit. The three-tier
classification is therefore not a classification of the system alone, but
of the tuple **(system, probe, readout, sector intervention)**.
```

---

## 4.2 `new_nogo_numerical_priorities.md` への保存量追加

対象:`## 1. 最優先で計算する量` セクション末尾(現行 41–54 行目)、
frozen-source cut の規約の直後。

### 追加ブロック(挿入位置:44行目の直後、「### sector cut の規約」の前)

```markdown
### 保存すべき量(operational cut 導入後の拡張)

`Sector/section0_operational_foundations.tex`(Theorem 0A)の導入後は、
上記三量に加えて次を必須保存量とする。

\[
\chi_{\mathrm{cut}}^{\mathrm{alg}},
\quad
\chi_{\mathrm{cut}}^{\mathrm{op}}(\kappa),
\quad
\chi_{\mathrm{cut}}^{\mathrm{sc}},
\]

\[
R_S^{\mathrm{path}},
\quad
R_S^{\mathrm{reprep}},
\quad
R_S^{\mathrm{total}}.
\]

さらに、\(\kappa\to\infty\) で

\[
\chi_{\mathrm{cut}}^{\mathrm{op}}(\kappa)-\chi_{\mathrm{cut}}^{\mathrm{alg}}
=O(\kappa^{-1})
\]

となることを、Gate U1 として毎回検証する(`Sector/src/run_gates_step3.py`
に実装済み、Λ 模型で確認済み:slope \(=-0.99997\)、閉形式係数との相対誤差
\(4.4\times10^{-4}\))。
```

---

## 4.3 QFI 命名の改訂(README.md 92–99 行目)

対象:`## Metrological (QFI) extension` セクションの第一段落。

### 現行

```
A further calculation tests the candidate prediction that the tangent-vector
difference x_S = ∂_θρ_full − ∂_θρ_cut obeys ‖x_S‖ ∼ Γ^−ν, and that the
sector-mediated quantum Fisher information F_{Q,S} = x_S^† G_ρ x_S ∼ Γ^−2ν
— i.e. that the same three-tier classification governs not just response
magnitude but how much parameter-estimation information a sector carries
under strong dissipation.
```

### 置換後

```
A further calculation tests the candidate prediction that the tangent-vector
difference x_S = ∂_θρ_full − ∂_θρ_cut obeys ‖x_S‖ ∼ Γ^−ν, and that its
metric lift I_sector_metric = x_S^† G_ρ̄ x_S ∼ Γ^−2ν — i.e. that the same
three-tier classification governs not just response magnitude but how much
parameter-estimation information a sector carries under strong dissipation.

I_sector_metric (also called the **sector-mediated tangent
distinguishability**) is *not* a decomposition of the quantum Fisher
information: in general F_Q[ρ_full] − F_Q[ρ_cut] ≠ F_Q[ρ_full − ρ_cut],
because a cross term survives (Gate U8, `Sector/src/gate_u8_naming_audit.py`).
The genuine QFI contrast ΔF_Q = F_Q(ρ_full) − F_Q(ρ_cut) is tracked
separately and is not claimed to equal I_sector_metric.
```

Rationale and the four theorem-level assumptions needed to derive
\(\|x_S\|\sim\Gamma^{-\nu_x}\Rightarrow\mathcal I_S^{(G)}\sim\Gamma^{-2\nu_x}\)
rigorously (faithful reference state \(\bar\rho\); \(\Gamma\)-uniform
boundedness and non-degeneracy of the metric tensor \(G_{\bar\rho}\); the
full/cut reference-state difference not changing the leading order;
avoidance of the rank-deficient limit) are recorded in the revision note
§8 and should be added as an explicit Assumption block wherever
`I_sector_metric` is defined in the metrology code/paper text.

### Gate U8 numerical result (negative control)

Run: `python3 "New no-go theory/Sector/src/gate_u8_naming_audit.py"`
(results also saved to `results/gate_u8_naming_audit.json`).

The correct bilinear-form identity, verified numerically against `qfi`/
`qfi_cross_term` directly, is
\[
\Delta F_Q = F_Q^{\mathrm{full}}-F_Q^{\mathrm{cut}}
=\mathcal I_S^{(G)}+2\times\text{cross term},
\]
not \(\Delta F_Q=\mathcal I_S^{(G)}+\text{cross term}\) as a naive reading
of the docstring of `qfi_cross_term` might suggest.

**Part 1 — synthetic qubit counterexample** (hand-built \(2\times2\)
example, \(\bar\rho=\mathrm{diag}(0.7,0.3)\), \(\partial_\theta\rho_{\mathrm{full}}\)
off-diagonal element \(a=0.3\), \(\partial_\theta\rho_{\mathrm{cut}}\)
off-diagonal element \(b=0.1i\)), which **PASSES**:

| Quantity | Value |
|---|---|
| \(F_Q^{\mathrm{full}}\) | \(0.36\) |
| \(F_Q^{\mathrm{cut}}\) | \(0.04\) |
| \(\Delta F_Q\) | \(0.32\) |
| \(\mathcal I_S^{(G)}=g_{\bar\rho}(x_S,x_S)\) | \(0.40\) |
| cross term | \(-0.04\) |
| residual \(\Delta F_Q-\mathcal I_S^{(G)}-2\times\text{cross}\) | \(\approx0\) (machine precision) |

\(\Delta F_Q=0.32\neq\mathcal I_S^{(G)}=0.40\): the two quantities are
**numerically confirmed to differ**, closing the counterexample requested
in Sec. 8 of the revision note and justifying the rename away from
"\(F_{Q,S}\)" (which reads as a component of \(\Delta F_Q\)) to the
metric-lift name \(\mathcal I_S^{(G)}\)/`I_sector_metric`.

**Part 2 — repository's physical 3-level Lindblad model** (informative,
not a pass/fail check), at `Gamma=50, theta=0.3, lam=0.7, phi=0.4`: the
cross term is numerically zero (\(\Delta F_Q=\mathcal I_S^{(G)}=1.37\times10^{-6}\)
exactly). This is a genuine structural feature of that particular
`kappa_cut=0` construction: with the coherent coupling switched fully off,
the cut Hamiltonian is purely diagonal, so \(\partial_\theta\rho_{\mathrm{cut}}\)
vanishes identically and \(x_S=\partial_\theta\rho_{\mathrm{full}}\) exactly
— consistent with, and explaining, the near-unit ratio
\(\nu_F/(2\nu_x)=1.0007\) already reported for Gate M3 in the top-level
README. It is not evidence against the naming distinction (Part 1 already
settles that in general); it documents why the distinction happens to be
numerically invisible in that one existing model.
