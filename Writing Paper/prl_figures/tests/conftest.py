"""Put the figure sources and the gate trees they read on sys.path.

The test module used to get `core`, `nv_reduced_kernel` and friends only
because `import fig1_classes` happened first and mutated sys.path on the way
in. Reordering the imports broke collection. Doing it here makes the
dependency explicit and independent of import order.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import _paths  # noqa: E402

_paths.add_sources("nogo", "phase", "gate_a", "gate_b", "gate_c", "gate_d", "roomt")
