"""
Load the NV- 3E model (Hamiltonian, legs, moments, kernel, graph_distance) from
the canonical verify_nv_3E_graph_distance_PRL.py WITHOUT running its test
battery, which executes on import and writes figures.

Single source of truth: nothing is copy-pasted; the model section of the
canonical file is exec'd verbatim, so the two cannot drift apart.
"""
from pathlib import Path
import types

_CANON = (Path(__file__).resolve().parents[2]
          / "No-go theorem" / "src" / "verify_nv_3E_graph_distance_PRL.py")


def load():
    src = _CANON.read_text().split("\n")
    # Model section = everything before the test battery header.
    end = next(i for i, l in enumerate(src) if l.startswith("# 4. Test battery"))
    mod = types.ModuleType("nv3e_model")
    mod.__dict__["__file__"] = str(_CANON)
    exec(compile("\n".join(src[:end]), str(_CANON), "exec"), mod.__dict__)
    return mod


if __name__ == "__main__":
    m = load()
    print("loaded from", _CANON)
    print("LAM_PERP =", m.LAM_PERP, " DELTA2 =", m.DELTA2, " LAM_Z =", m.LAM_Z)
