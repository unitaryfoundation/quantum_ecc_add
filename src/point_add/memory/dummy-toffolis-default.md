# DUMMY_TOFFOLIS default was inflating the scored Toffoli count

## What

`build()` ends with a `DUMMY_TOFFOLIS` block that emits `n` self-cancelling
`CCX(a,b,c); CCX(a,b,c)` pairs on a fresh ancilla. It is a **diagnostic** for
harness sensitivity tests — the pairs are logically identity (Toffoli is
self-inverse) and the ancilla returns to |0>, so peak qubits are unchanged.

The committed default was `unwrap_or(7_400)`. Every emitted CCX still counts
toward the scored *average executed Toffoli* (`sim.rs` increments
`toffoli_gates` per active shot for every CCX/CCZ regardless of whether the
controls fire). So the default silently added **2 × 7_400 = 14_800 Toffoli/shot**
of pure dead weight to the score.

## Impact

- Committed repo (default 7_400): 3_957_553 Toffoli × 2_715 qubits = 1.0745e10
- After default → 0:              3_942_753 Toffoli × 2_715 qubits = 1.0704e10

Setting the default to 0 restores the README's quoted "current main" baseline
(score 10_704_574_395) and is a strict, verified improvement over the repo as
committed. Validated end-to-end: all 9_024 shots OK, score.json confirms.

## How to apply

`unwrap_or(7_400)` → `unwrap_or(0)` in the `DUMMY_TOFFOLIS` block in `build()`.
Knob retained for diagnostics via `DUMMY_TOFFOLIS=<n>`.

## Note for the next agent: where the real wins are

This fix only ties the README baseline; beating 1.07e10 *strictly* needs a real
reduction. Measured structure (DUMMY_TOFFOLIS=0):

- **Peak qubits = 2715**, hit in `sol_addlo` (Solinas reduction inside a
  Karatsuba multiply during Kaliski inversion). Owner breakdown:
  `pair1_kaliski_forward`=916, `pair1_mul1`=770, `sol_addlo`=517, `init`=512.
  Disabling `KAL_PAIR1_MUL1_KARATSUBA` drops peak to 2712 but the binding
  constraint then becomes the Kaliski coefficient state itself
  (`bk_step4`=2712 = 916+772+512+512). So peak is gated by the ~1688 qubits of
  Kaliski forward+backward state — reducing that is the real width lever.
- **Toffoli budget** is ~60% Kaliski `step4` + `step3/step9_cswap` (forward
  `kal_*` and backward `bk_*` mirror each other). Optimize the Kaliski
  iteration to move the gate count.
- Multiply/shift knobs tested are all net-negative on the *product* score:
  `KAL_PAIR1_MUL1_KARATSUBA=0` → 3_979_605 × 2712 = 1.079e10 (worse).
  Karatsuba-on defaults are correctly tuned.
- `COMPACT_POINT_ADD=1` (Fermat-inversion, fewer qubits) is **unfinished** —
  its cleanup is a TODO and it fails validation. A working compact/Fermat path
  is the most promising direction for a large width cut if completed.
