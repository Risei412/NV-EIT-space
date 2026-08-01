# NV–EIT PRL Gates A–E Calculation Bundle

**Frozen:** 2026-08-01  
**Repository target:** `Risei412/NV-EIT-space`  
**Purpose:** preserve the independent calculations, kill-gate decisions, numerical tables, and measurement-design results supporting the NV–EIT PRL route.

## Decision ledger

| Stage | Decision | What was established |
|---|---|---|
| Gate A | Conditional PASS | Marked complex response and signed absorption inherit the path-moment order; normalization changes the observable order. |
| Gate B | Structural PASS | A non-EIT circuit-QED full-GKSL witness reproduces amplitude classes 1 and 2 and the expected symmetry-breaking crossover. |
| Gate C | PASS | Classes 1 and 2 recur across diamond and circuit-QED platforms; a minimal non-NV three-mode chain supplies the missing class-3 witness. |
| Gate D | Model-level PASS | A threefold finite-dissipation window retains an effective amplitude order near 3 at finite SNR. |
| Gate D-final | PASS with selected protocol | An eight-state Thue–Morse coupling-sign cycle rejects static feedthrough and smooth background drift and preserves quantitative order recovery. |
| Gate E | Model-level CONDITIONAL PASS | The joint engineerable–asymptotic–detectable window fails for optical NV and the Gate-B SC transfer, but a 15–45 MHz three-mode-chain window separates classes 2 and 3 at finite SNR. |

## Central non-NV class-3 witness

The minimal active graph is a three-mode chain,

```text
input -> mode 1 --J12-- mode 2 --J23-- mode 3 -> marked output
```

with a common ground state and engineered loss on the three active modes. For source `c=|1>` and readout `p=|3>`,

```text
M0 = 0
M1 = 0
M2 = -J12*J23 != 0
```

so the amplitude response scales as `Gamma^-3`. The full 16-dimensional GKSL calculation gives an amplitude exponent `2.9999708` and fixed-readout population exponent `5.9999417`.

## Selected measurement protocol

The final calculation selects coupling-sign phase cycling with the eight-state Thue–Morse sequence

```text
+ - - + - + + -
```

which cancels background Taylor moments of temporal degree 0, 1, and 2. In the model audit it recovers the class-3 exponent with a background up to `10^4` times the strongest marked signal, subject to a measured smooth-drift condition and switching-transient control.

## Gate E claim boundary

Gate E requires the same loss interval to be engineerable, asymptotic enough for class
identification, and detectable. Optical NV and the Gate-B superconducting transfer witness fail
this joint test. The three-mode chain gives a model-level conditional pass for
`Gamma/2pi = 15–45 MHz`: the class-3 95% slope interval `[2.8286, 2.9779]` is disjoint from
the class-2 interval `[2.0680, 2.1090]`, and the worst point reaches SNR 10.21 with 40 s
integration under the stated 3 K, −110 dBm input model.

This is not a device-level experimental pass. Tunable loss hardware, amplifier calibration near
−190 dBm, coupling-sign switching transients, and the measured feedthrough drift spectrum remain
open engineering checks.

## Directory layout

```text
reports/          consolidated human-readable Gate A–D audit
src/              standalone calculation scripts
results/tables/   compact JSON output used in the decisions
```

## Reproduction

```bash
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python src/gate_a_independent_audit.py
python src/gate_b_independent_audit.py
python src/gate_c_independent_audit.py
python src/minimal_nonNV_n3_witness.py
python src/gate_d_measurability.py
python src/gate_d_final_background.py
```

## Scope and claim boundary

These calculations support a **structural and measurement-model** universality claim for marked path-moment response orders. They do not replace a device-specific microwave-layout simulation or an experimental measurement of feedthrough drift, coupler switching transients, and amplifier-chain calibration.

The complete ZIP bundle additionally contains the generated figures and the full Monte Carlo raw tables.
