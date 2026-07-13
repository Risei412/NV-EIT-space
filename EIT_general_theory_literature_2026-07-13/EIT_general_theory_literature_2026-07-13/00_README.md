# 一般化EIT理論：競合・補強・適用文献パッケージ

更新日: 2026-07-13

## 目的

添付理論を、NV-固有の議論から「熱的に広がった多励起準位系に対する一般的なEIT判定理論」へ展開するため、文献を次の3群に整理した。

1. **競合・先行理論**: dark stateのrank/SVD判定、多準位EIT、散逸系、広がり、EIT/ATS判別。
2. **補強文献**: NVおよびgroup-IV中心の対称性、微視的Hamiltonian、phonon、温度依存、strain。
3. **適用・検証文献**: SiV/PbV/SnV、SiC、希土類、半導体量子ドット、GaAs、hBN。

## 主要結論

- `rank[p c]=1`、nullspace、SVDによる**exact dark stateの存在条件は新規ではない**。
- 一般論として最も防御可能な新規核は、次の組合せである。

  `K(omega,T) = p^dagger G_e(omega,T) c`

  1. frequency/temperature-dependent coherent kernel
  2. controlled merged-manifold expansion
  3. `d_p^dagger Pi_E d_c` のpoint-group selection rule
  4. exact dark state / finite EIT / zero-kernel / observable contrast の分離
  5. no-go・partial-go・tunable・practical-goの材料分類
  6. resolved-to-merged crossoverの反証可能な予測

- 「orbital-degenerate groundならgo」は強すぎる。正しくは、**spin overlap、表現のtensor積、偏光でアクセス可能な非scalar成分、非零のreduced matrix element、coherence/power budget**をすべて満たす必要がある。
- `Gamma_opt >= Delta_ES`だけでは `G_e ≈ gI` は保証されない。branch-dependent linewidth、self-energy、非Markov phonon、incomplete manifoldを誤差評価に入れる。
- SiVの `eta_geom ≈ 1` は有力な予測だが、既報SiV-CPTは主に低温resolved spin-Lambdaであり、室温merged orbital-Lambdaを直接検証してはいない。

## 収録物

- `01_theory_positioning.md`: 競合を踏まえた理論の再定義
- `02_competing_literature.md`: 競合・先行文献一覧
- `03_reinforcement_literature.md`: 補強文献一覧
- `04_application_literature.md`: 材料別の適用・実証文献
- `05_claim_to_citation_map.md`: 主張ごとの引用先と必要な修正
- `06_validation_workflow.md`: 各材料への適用手順
- `07_priority_reading_list.md`: 読む順序
- `paper_matrix.csv`: 全文献の比較表
- `references.bib`: BibTeX
- `paper_notes/`: 論文ごとの個別ノート
- `source_theory/`: 今回与えられた2つの理論PDF
- `search_method.md`: 検索範囲と限界

## 著作権

第三者論文のPDF全文は格納していない。DOI/arXiv/publisher URLと分析ノートのみを収録した。`source_theory/`にはユーザー提供の理論PDFのみを格納した。
