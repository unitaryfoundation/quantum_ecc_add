# Next-session plan for reaching SOTA

## Current position (committed, stable)
- **4.18M Toffoli / 2716 qubits**
- Matches Litinski 2023 frontier
- Gap to Google: 1.55× Toffoli, 2.3× qubits

## Debugging progress

### Venting adder infrastructure (DONE)
9 tested primitives in `src/point_add/venting.rs`:
- Classical-offset: carry-xor, 2-clean vented add, HRS linear-clean,
  dirty-2-clean, controlled add/sub variants (ciadd, cisub).
- Quantum-offset: 2-clean vented qadd, dirty-2-clean qadd, qsub.

### u64 shift UB bug (FIXED)
Was the root cause of the first wiring attempt's 320 phase batches.
Rust release-mode masks shift amounts modulo bit width. Fixed with
`if k >= 64 { false } else { (x >> k) & 1 != 0 }`.

### Seed-3 phase leak (IDENTIFIED via bisection)
Wiring venting halve into backward Kaliski produces 1 phase batch out
of 20480 shots at seed=3. Bisection reveals it's a **cross-call phase
interaction**, not a primitive bug:
- 0..577 calls enabled: 0 phase batches
- 0..578 calls: 0 phase batches
- 0..579 calls: 2 phase batches (!!)
- 0..1000 calls: 1 phase batch
- 0..full: 1 phase batch

Different subsets of calls produce different phase counts, meaning
Gidney's phase corrections via `cz_if(dirty, vent_keys)` aren't
fully composing across sequential calls when dirty qubits are
shared across calls.

## Concrete next steps

### 1. Debug the cross-call phase interaction (HIGH PRIORITY)
Compare my Rust port's emitted gate sequence against Gidney's Python
reference at identical inputs. Specifically: run both at n=256 with 2
sequential `cisub_dirty_2clean_classical` calls sharing dirty qubits,
check phase invariant.

If difference found: it's likely in how `vent_keys` phase corrections
sandwich the carry-xor. The Python uses explicit `broadcast_cz(dirty,
vent[1:])` with sandwich structure; my port may have subtle ordering
difference.

### 2. Wire qoffset venting into correction-3 (after 1)
`isub_dirty_2clean_qoffset` is tested and phase-clean at n=256. Wiring
into schoolbook_mul correction-3 requires:
- Handle 2n+1-wide sub via n+1-wide qoffset call + borrow ripple.
- Borrow bit captured in wide[n+1].
- Ripple decrements wide[n+2..2n+1] conditional on borrow.

Expected: pair1_mul1 peak 2716 → ~2460.

### 3. Check for other peak sites
Currently 4 phases at 2716 (pair1_mul1, pair1_mul2, bk_bulk_step6_7_8,
bk_step6_7_8). Reducing ONE doesn't reduce global. Need to attack ALL
4 simultaneously via their respective peak-triggers.

### 4. Fallback: reduce m_hist via Bennett pebble game
m_hist is 407 qubits, non-trivial fraction of base alive. Pebble game
split: store only sqrt(407) ≈ 20 checkpoints, recompute intermediate
m_hist bits on backward. Saves ~387 qubits, costs ~1M Toffoli.

Expected: peak 2716 → ~2330, Toffoli 4.18M → ~5.2M.

### 5. Ambitious: Luo 2025 register sharing
Full rewrite using location-controlled arithmetic. Reaches ~1300 qubits
but at 100-500M Toffoli. Multi-session effort.

## Files to read next session
- `src/point_add/venting.rs` — the 9 primitives
- `src/point_add/session_summary.md` — state summary
- `/tmp/gidney_venting/code/src/constadd/` — Python reference for
  venting adder
- `/tmp/zenodo_ecc/extracted/` — Google's ZKP harness (matches ours)

## Baseline commit
`8a9737b` — Primitive infrastructure + 9 tests. Baseline 4.18M/2716q
unchanged.
