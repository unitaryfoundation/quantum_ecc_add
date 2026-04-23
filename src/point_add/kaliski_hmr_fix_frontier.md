# Frontier after the HMR-count phase fix

After restoring the missing generic step-0 HMR history into the specialized
bulk-prefix step, the strict `main.rs` frontier changed as follows.

## Passing strict settings
- `k = 4`  → 4,392,466
- `k = 16` → 4,386,226
- `k = 72` → 4,357,106
- `k = 80` → 4,352,946
- `k = 112` → 4,336,306

## Failing strict settings
- `k = 8`   → classical mismatch
- `k = 24`  → phase
- `k = 32`  → phase
- `k = 40`  → phase
- `k = 64`  → phase
- `k = 96`  → classical mismatch
- `k = 128` → classical + phase

## Interpretation
The HMR-history repair was a real phase-bug fix:
- it repaired the strict `k = 4` failure,
- and it moved the passing frontier,
- but it did not restore the original pre-fix `k = 96` passing result.

So there are at least two distinct issues mixed together:
1. a phase-history/HMR mismatch (partially fixed),
2. a remaining classical or phase incompatibility that still blocks many larger
   prefix lengths.
