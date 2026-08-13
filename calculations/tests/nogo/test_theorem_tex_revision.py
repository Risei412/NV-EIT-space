from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PROOF_DIR = ROOT / "theory" / "proofs"
ENTRY = PROOF_DIR / "eit_nogo_proofs.tex"
V7 = PROOF_DIR / "eit_nogo_proofs_v7.tex"
V6_ARCHIVE = ROOT / "archive" / "historical_notes" / "theory" / "eit_nogo_proofs_v6_2_archive.tex"


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_entrypoint_routes_to_v7_and_archive_is_retained() -> None:
    assert V7.is_file()
    assert V6_ARCHIVE.is_file()
    assert "eit_nogo_proofs_v7.tex" in _text(ENTRY)


def test_v7_contains_the_required_revision_contract() -> None:
    text = _text(V7)
    required = (
        "Laurent moment expansion",
        "Neumann expansion with explicit error",
        "Sector-graph lower bound",
        "Graph distance is a lower bound, not an automatic equality",
        "Exact transfer-zero/Krylov criterion",
        "Sector-observable inheritance",
        "Restriction of Theorem~\\ref{thm:singular-D}",
        "Non-semisimple zero mode",
        "fixed-scaling family",
    )
    for marker in required:
        assert marker in text, marker


def test_sector_graph_statement_separates_bound_from_equality() -> None:
    text = _text(V7)
    theorem_start = text.index("\\begin{theorem}[Sector-graph lower bound]")
    equality_start = text.index("\\begin{corollary}[When graph distance equals the dissipation order]")
    theorem_block = text[theorem_start:equality_start]
    assert "M_n=p^\\dagger X^n\\nu=0" in theorem_block
    assert "n<d" in theorem_block
    assert "K_\\Gamma=O(\\Gamma^{-d-1})" in theorem_block
    assert "M_d" not in theorem_block
    equality_block = text[equality_start:]
    assert "M_d=p^\\dagger(X_0+W)^d\\nu\\ne0" in equality_block
    assert "M_d\\Gamma^{-d-1}" in equality_block


def test_laurent_theorem_is_not_conditioned_on_neumann_norm() -> None:
    text = _text(V7)
    start = text.index("\\begin{theorem}[Laurent moment expansion]")
    end = text.index("\\end{theorem}", start)
    theorem_block = text[start:end]
    assert "\\|X\\|/|\\Gamma|<1" not in theorem_block
    assert "no operator-norm hypothesis" in theorem_block


def test_tex_environments_are_balanced() -> None:
    text = _text(V7)
    begins = Counter(re.findall(r"\\begin\{([^}]+)\}", text))
    ends = Counter(re.findall(r"\\end\{([^}]+)\}", text))
    assert begins == ends
    assert text.count("\\begin{document}") == 1
    assert text.count("\\end{document}") == 1
