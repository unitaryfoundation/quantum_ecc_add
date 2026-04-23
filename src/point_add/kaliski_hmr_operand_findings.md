# HMR operand finding

After the HMR-history repair, the specialized bulk-prefix step now matches the
**generic HMR count**, but it still does **not** match the HMR operand sequence.

## Diagnostic result
For iter indices `0, 1, 2, 3, 7, 15, 31`:
- generic and specialized HMR counts match exactly,
- but the common-prefix operand sequence differs in **767** positions.

The first differences are systematic:
- generic starts at `QubitId(1283), 1282, 1281, ...`
- specialized starts at `QubitId(1284), 1283, 1282, ...`

So the specialized HMR stream is effectively shifted relative to the generic one.

## Interpretation
This means the phase bug is not just “missing HMRs” anymore.
Even after restoring the count, the measurement targets themselves are still not
aligned with the generic step, so the phase history remains different.

This is now the strongest concrete residual low-level bug: the specialized bulk-
prefix step must reproduce not only the number of HMRs but the **same operand
sequence** if it is to be phase-equivalent to the generic step.
