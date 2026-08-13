# manuscript/figures/

Deliberately empty of binaries.

`main.tex` sets `\graphicspath{{../results/figures/}}`, so the figures the
manuscript includes are the generated artifacts in `results/figures/` — the
same files the gate regression tests check. Copying them here would create a
second copy that can silently drift from the result it depicts, which is
exactly the failure mode `NON_CLAIMS.md` N8 records.

Put a file here only if it is a hand-drawn figure with no generating script.
