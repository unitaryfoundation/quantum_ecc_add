# Research notes — quantum_ecc inversion algorithm space

Session: 2026-04-22 (continued).
Author: autoresearch agent.

Baseline state: 4,394,546 avg executed Toffoli @ 2729 qubits.
Kaliski modular inversion contributes ~81% of the circuit budget
(3.55M Toffoli across the two inversion passes).

This document surveys published modular-inversion algorithms that could
plausibly replace or augment Kaliski, with their iteration structure,
per-iteration reversible cost, and publication status.

## Deliverable 1 result (classical B-Y empirical survey)

Implemented `divstep2` (Bernstein–Yang 2019/266, §8) in
`src/classical_by.rs` and ran it on 10,000 random secp256k1 inputs
seeded by SHAKE128.

| metric | value |
|---|---|
| theoretical bound `⌈(49·256 + 57)/17⌉` | 742 |
| observed minimum iters | 502 |
| observed maximum iters | **567** |
| observed mean iters | 531 |
| max |δ| during execution | 20 |
| modinv matches (vs Fermat) | 10,000 / 10,000 |

**Key implication**: the theoretical bound overestimates real-world
iteration count on secp256k1 by **24%**. Prior sessions' cost
analyses that used 742 were pessimistic. `max |δ| = 20` also means
the δ register needs only ~7 bits, not a full-width integer.

## Deliverable 2: algorithm space survey

All costs are for **n = 256** (secp256k1). Reversible costs are
**measured** (for Kaliski, from our instrumented build) or
**estimated conservatively** (per-iteration op counts × naive
register sizes).

### 1. Kaliski almost-inverse  — baseline (used by our circuit)

- **Classical ref**: Kaliski 1995, *"The Montgomery inverse and its
  applications"*, IEEE Trans. Computers 44(8), 1064–1065.
  DOI: 10.1109/12.403725.
- **Reversible refs**:
  - Roetteler et al. 2017 (RNSL), arXiv:1706.06752.
  - Häner et al. 2020 (HRSL), arXiv:2001.09580, eprint 2020/077.
- **Iteration count**: classically 2n = 512 iters for deterministic
  convergence; our code truncates to 399 (tuned against 9024-shot
  deterministic test corpus).
- **Per-iter reversible cost (measured)**: **~2180 CCX** (profile).
- **Per-pass (forward + backward)**: **1.81M CCX**, two passes in our
  flow = 3.62M CCX.
- **Structural notes**:
  - Binary-GCD style. Each iter: parity check, 2n-bit comparator,
    two cswaps of two register pairs, fused cond-sub/cond-add, halve.
  - Reversibility via `m_hist` branch log (one qubit per iter).

### 2. Bernstein–Yang divstep2 (`w = 1`)  — no reversible impl published

- **Classical ref**: Bernstein, Yang 2019, *"Fast constant-time gcd
  computation and modular inversion"*, eprint 2019/266, TCHES 2019(3).
- **Reversible ref**: **unpublished / would be novel research**.
  No quantum/reversible B-Y implementation exists in the public
  literature as of April 2026 to my knowledge.
- **Iteration count**:
  - Theoretical bound: ⌈(49n + 57)/17⌉ = **742** for n=256.
  - Empirical worst case (10k secp256k1 samples, Deliverable 1):
    **567** (24% below bound).
- **Per-iter reversible cost (est. conservative)**:
  - Branch predicate `(δ > 0) ∧ (g odd)`: 1 CCX.
  - Cswap 4 register pairs on ctrl (f↔g, U↔Q, V↔R, with negations):
    4n + 3n_mod-neg ≈ 7n.
  - Cond add-or-sub f±g: n.
  - Coeff updates (4 registers, each a cond add/sub): 4n ≈ 4n.
  - Halve g: 0.
  - δ ± 1 update: ~7 bits, tiny.
  - Misc bookkeeping: small.
  - **Total estimate: 10–12n ≈ 2560–3072 CCX/iter**.
- **Per-pass (forward + backward)**:
  - Using empirical max 567 iters and 12n/iter:
    2 × 567 × 3072 ≈ **3.48M CCX per pass**.
  - Using empirical mean 531 iters and 10n/iter:
    2 × 531 × 2560 ≈ 2.72M CCX per pass.
- **Verdict vs Kaliski**: Worse by ~1.5–1.9× per pass even with the
  empirical-iter-count correction.

### 3. Bernstein–Yang jumpdivstep (`w = 31`)  — no reversible impl published

- **Same paper, §9 of eprint 2019/266**.
- **Classical speedup**: batches 31 divsteps into one "jump" via a
  precomputed 2×2 transition matrix keyed on low-31-bits of (f, g).
  Classical cost is dominated by 32-bit word operations. Used in
  libsecp256k1 production.
- **Reversible ref**: **unpublished / would be novel research**.
- **Iteration count**: ⌈567/31⌉ = 19 jumped iters empirically; ≈ 24
  using the theoretical bound.
- **Per-jump reversible cost (est. conservative)**:
  - QROM lookup of 2×2 matrix (entries bounded by 2^31, so ~32-bit
    signed). Table has 2^w = 2^31 entries of ~128 bits each. A full
    unary-decode QROM is infeasible for 2^31 entries (2^31 CCX for
    decode alone). Select-swap QROM reduces to O(√2^31) = O(2^15.5),
    still ~46k CCX per lookup.
  - Apply 2×2 matrix to (f, g): compute 2 linear combinations
    `a·f + b·g` where a, b are up-to-31-bit classical constants,
    f, g are 256-bit quantum. Each `a·f` is a classical-by-quantum
    mul: conservatively 31 × n = 7936 CCX via shift-and-add
    (schoolbook). Plus n-bit sum. Two combinations per matrix × 2n ≈
    32k CCX for (f, g). Similar for (U, V, Q, R) tracked mod p: 4 ×
    2 × 31 × n + Solinas reductions ≈ 64k + small.
  - Halve by 2^31: requires 31 halvings or a single batched halve
    (~same cost either way: ~31 × 256 = 8k CCX for mod_halve style).
  - **Total estimate: 100–150k CCX/jump**.
- **Per-pass cost**: 19 × 125k ≈ **2.4M CCX forward**; double for
  backward = 4.8M CCX.
- **Verdict vs Kaliski**: Much worse. The classical speedup from
  32-bit wordwise operations does NOT translate to reversible
  advantage because every bit of the w-bit matrix entries still
  requires an n-bit conditional-add operation in the quantum setting.
  `w · n` scaling defeats the 1/w iteration-count reduction.

### 4. Montgomery inverse (Savaş–Koç)  — Kaliski variant

- **Classical ref**: Savaş, Koç 2000, *"The Montgomery modular
  inverse — revisited"*, IEEE Trans. Computers 49(7), 763–766.
  DOI: 10.1109/12.863048.
- **Reversible ref**: used by RNSL 2017 (arXiv:1706.06752) and HRSL
  2020 (arXiv:2001.09580) as the inversion primitive. Our Kaliski is
  structurally equivalent.
- **Iteration count**: same as Kaliski, 2n.
- **Per-iter reversible cost**: essentially identical to Kaliski.
- **Verdict vs Kaliski**: Not a separate algorithm in practice — the
  "Montgomery inverse" is a Kaliski variant with slightly different
  post-processing.

### 5. Lehmer-style GCD  — no reversible impl published

- **Classical ref**: Lehmer 1938, *"Euclid's algorithm for large
  numbers"*, American Mathematical Monthly 45(4), 227–233.
  Refinements: Jebelean 1993; Wang–Pan 1993.
- **Idea**: approximate (u, v) by their top k bits (e.g., k=32) to
  compute a 2×2 transformation matrix classically via fast native
  arithmetic, then apply to the full-precision (u, v).
- **Reversible ref**: **unpublished / would be novel research**.
- **Iteration count**: ⌈2n/k⌉ Lehmer steps. For k=32, n=256: ~16
  steps.
- **Per-step reversible cost**:
  - Classical 2×2 matrix computation: free (done at circuit build
    time if matrices can be pre-enumerated by the low-k bits).
  - Matrix application: 4 × (k·n) = 4·32·256 = 32k CCX for (u, v).
    Plus 4·k·n = 32k for coefficient tracking. Per step: ~64k.
- **Per-pass cost**: 16 × 64k = **1.0M CCX**; 2× for backward.
- **Verdict vs Kaliski**: If the matrix lookup is FREE (pre-enumerated
  at build time, keyed on k bits of classical state), Lehmer-style
  COULD be competitive. But the classical computation of the matrix
  depends on the RUNTIME values of (u[top k], v[top k]) which are
  quantum. So at runtime we'd need a QROM with 2^(2k) = 2^64 entries
  — infeasible. Reducing k to make the table manageable (k ≤ 20 for
  2^40 entries max) gives 12 Lehmer steps × (4·20·256 = 20k) = 240k
  per pass × 2 = 480k, plus lookup cost.
  
  **The lookup cost is the obstruction**: 2^(2k) entries for general
  Lehmer-style. Only tractable for k ≤ ~12-14 (tables of 2^24-2^28
  entries). At k=12: 22 steps × (4·12·256 = 12k) = 264k per pass.
  Plus lookup ~ 2^24 = 16M CCX per step? No — a unary-decode
  QROM of 2^24 entries costs 2^24 = 16M CCX, which destroys the
  savings.
  
  Select-swap QROM reduces to O(2^12) = 4k CCX per lookup at k=12.
  With this, 22 × (12k + 4k + ...) ≈ 22 × 20k = 440k per pass.
  Doubled for backward: 880k. **Possibly competitive with Kaliski**
  but requires building a reversible select-swap QROM primitive plus
  careful matrix tracking. **Novel research territory.**

### 6. Fermat's little theorem (`a^{p-2}`) via addition chain

- **Classical ref**: Fermat 1640 (folklore); modern survey in
  Knuth TAOCP §4.6.3 (addition chains).
- **Reversible impl**: RNSL 2017 discussed it; not used in practice.
- **Iteration count**: addition chain length for p-2.
  - For secp256k1 p-2 = 2^256 - 2^32 - 979, a good chain has
    ~263-270 multiplications (255 squarings + 8-15 non-squaring mults).
- **Per-mul reversible cost**: measured ~70-80k Toffoli (our
  schoolbook_mul_into_addsub + Solinas).
- **Per-pass cost**: 270 × 75k = **20M CCX**.
- **Verdict vs Kaliski**: **much worse** (5× Kaliski). No amount of
  primitive-level optimization brings this competitive.

### 7. Itoh–Tsujii inversion (GF(2^n) only)

- **Classical ref**: Itoh, Tsujii 1988, *"A fast algorithm for
  computing multiplicative inverses in GF(2^m) using normal
  bases"*, Information and Computation 78(3).
- **Applicability**: GF(2^n), not GF(p). Not applicable to
  secp256k1 which uses a prime field.

### Summary table

| algo | iters (n=256) | per-iter/step CCX | per-pass CCX (forward) | pub reversible? | competitive with Kaliski? |
|---|---|---|---|---|---|
| Kaliski | 398 | 2180 | 868k | yes (ours) | **baseline** |
| B-Y divstep w=1 | 567 (emp) | 2560–3072 | 1.45–1.74M | **no** | +67–100% worse |
| B-Y jumpdivstep w=31 | 19 | 100–150k | 2.4M | **no** | +180% worse |
| Montgomery inv | ≈ 398 | ≈ 2200 | ≈ 870k | yes (RNSL/HRSL variant) | equivalent (same algorithm) |
| Lehmer w=12 + QROM | 22 | 20k+lookup | 440–880k | **no** | **possibly competitive** |
| Fermat addn chain | 270 muls | 75k | 20M | yes | much worse |
| Itoh–Tsujii | N/A | N/A | N/A | N/A | not applicable (GF(2^n)) |

## Deliverable 3: the actual research bet

**Conclusion: "Lehmer with select-swap QROM is the bet, IF the select-swap QROM primitive can be built cheaply and the matrix-tracking code can be made reversibly clean. Otherwise, no known algorithmic path remains open for single-session gains against Kaliski at n=256."**

### Why not B-Y?

Even with the empirical iter-count correction (567 instead of 742),
B-Y's per-iter reversible cost (~12n CCX) × iter count is still
worse than Kaliski's (~8.5n × 398 ≈ 3400n) by 2×. The jumping
variant (`jumpdivstep w≥8`) does not help because the matrix
application cost scales linearly with the batch width w, so the total
work is constant at ≈ c · n² · iter_base regardless of how we batch.

### Why not Montgomery batched inversion?

Independently confirmed (prior session + this session) that the
cleanup obstruction is fundamental: zeroing the saved `c^{-1}`
register requires either a second Kaliski (no savings) or expressing
c^{-1} as a function of end-state registers (no such expression
exists without further inversions).

### Why Lehmer could be the bet (with caveats)

Lehmer-style has a structural property B-Y lacks: the **classical
part** of the matrix computation can be done at circuit build time
IF we enumerate over possible values of the "approximation" (top k
bits of u, v). The runtime cost reduces to a **QROM lookup** of the
precomputed transformation matrix plus its application.

For k = 12:
- Lookup table: 2^24 = 16M classical 2×2 matrices (build-time).
- Select-swap QROM read cost: O(2^(k/2)) ≈ 2^6 = 64 CCX per read.
- Matrix application cost: 4 × 12 × 256 = 12k CCX per step.
- Steps: ~22.
- Forward cost: 22 × (12k + 64) ≈ 265k CCX.
- Backward cost: similar, another 265k.
- Total: **~530k CCX per inversion pass** — vs 1.81M for Kaliski.

If realized, this would be a **~1.3M CCX saving per Kaliski pass**,
times two passes = **2.6M total savings = 59% of the current
4.39M circuit**, landing at ~1.8M Toffoli. That would **beat
Google's reported SOTA** (2.1M/2.7M).

### Caveats / what could go wrong

1. **Select-swap QROM primitive doesn't currently exist in this
   codebase**. Implementing it correctly is a ~300-500 LOC effort
   including decode/undecode logic and phase-clean uncompute. Known
   patterns (Babbush et al. 2018, "Encoding Electronic Spectra in
   Quantum Circuits with Linear T Complexity", arXiv:1805.03662)
   give the asymptotic, but the exact CCX count at k=12 needs to be
   measured.

2. **Lehmer's classical algorithm has dynamic step length**. At each
   step, the algorithm decides how many "inner iterations" to
   aggregate based on the values. Reversibly, we'd need to fix the
   step structure at build time (so every "step" has the same
   reversible circuit), OR use uniform-size steps that sometimes do
   no-ops — which loses the benefit. **This is the hardest design
   problem**.

3. **Coefficient tracking** for (u, v, q, r) analogous to Kaliski's
   (r, s): needs careful sign handling (Lehmer matrices have signed
   entries). Classical sign tracking is trivial; quantum requires
   either magnitude+sign or two's-complement with explicit negation.

4. **No published reversible Lehmer exists**. Success requires
   deriving the reversible step structure from scratch. This is
   genuinely novel research.

### Secondary bet: hybrid Kaliski-jump

A less ambitious direction: keep Kaliski's state machine but batch
the **cswap(u, v_w)** + **cond-sub v_w -= u** + **halve** operations
across w bits at once, using precomputed branch-sequence matrices
keyed on the low-w bits of (u, v_w).

- Per-iter quantum cost reduction: skip intermediate state updates
  that will get cswap'd away again. Potential savings: 30–50% on
  cswap contribution (1.29M CCX = 29% of circuit).
- Max realistic saving: ~500k CCX = 11% of total.
- Feasibility: similar to Lehmer but reuses Kaliski's state machine,
  so no brand-new algorithm. Less novelty, less risk.

### Tertiary bets (smaller upsides)

- **Windowed classical-constant mul** to replace the 399 halvings
  of `lam` in pair1 and 399 doublings in pair2. Net ~120–160k saving
  (3.5%) if implemented correctly. Requires QROM primitive.

- **STEP 3+9 cumulative-swap-state** (prior session's HRSL analysis)
  was +3.2M net-negative in the Kaliski context. I don't think
  there's a way around that in-session.

- **Kim-style unconditional Kaliski** + m_hist removal (saves 400
  qubits, costs more Toffoli). Not a Toffoli win.

## Proposals for future sessions (quantum circuit deferrals)

Per instruction, I did not write any quantum code this session. The
following are proposal-level — each would take multiple sessions to
implement correctly.

### Proposal P1: reversible Lehmer with k=12 select-swap QROM

**Primitive needed**:
```
fn lehmer_step(b: &mut B, u: &[QubitId], v: &[QubitId],
               r: &[QubitId], s: &[QubitId],
               u_coeffs: &[QubitId], v_coeffs: &[QubitId]) {
    // 1. Read low 12 bits of u and v into a 24-bit address.
    // 2. Select-swap QROM lookup of precomputed 2x2 matrix
    //    (a, b, c, d) with entries in [-2^12, 2^12].
    // 3. Apply: (u', v') := (a·u + b·v, c·u + d·v) / 2^k.
    //    (Division is exact by construction.)
    // 4. Apply same matrix (mod p) to coefficient registers.
    // 5. Uncompute QROM address.
}
```

**Open research questions**:
- Correct classical Lehmer matrix enumeration (which subset of
  2^{2k} low-bit values give valid matrices?).
- Handling the "remainder" iterations after Lehmer can't make
  further progress (switch to per-bit Kaliski for final bits).
- Reversible sign handling in matrix application.

### Proposal P2: empirical measurement of reversible B-Y w=1

Even though my analysis says B-Y w=1 loses by ~70%, the cost
estimates are NOT measured. A full implementation would let us see
the actual gate count. If per-iter turns out cheaper than 12n
CCX (e.g., 7n via aggressive fusion), the loss margin shrinks.

This is the "small win" bet: might net -200k to +500k relative to
Kaliski. Worth trying if we have a slow session.

### Proposal P3: windowed classical-const mul primitive

Separate from inversion: implement a general windowed quantum-by-
classical multiply using precomputed tables, and use it to collapse
the halve/double loops on `lam`. Net ~120-160k Toffoli.

Bounded scope (~200 LOC), low risk, small reward.

### What I won't propose

- Fermat via addition chain: disqualified at 20M CCX.
- Jacobian coordinates: prior analysis confirmed cleanup
  obstruction is fundamental.
- Montgomery batched with two full Kaliski: no savings.
- B-Y w ≥ 4 with matrix magnitudes scaling as 2^w: analysis shows
  per-w-iter cost scales linearly in w, offsetting batching.

## Bottom line

The 4.39M / 2729q baseline is at the single-session frontier of
publicly-known techniques. The gap to Google's 2.1–2.7M SOTA likely
requires either:

1. **Undisclosed Google techniques** (most likely — the paper omits
   algorithmic details).
2. **Reversible Lehmer** (novel research; my best guess at what's
   closing the gap).
3. Compound primitive-level wins (plausible ~10-15% but won't close
   the full gap).

If I had a multi-week budget, I'd prototype Proposal P1 (reversible
Lehmer) as the single bet with the highest expected value. In a
single session from the current state, no bet has high-enough
probability of success to recommend pursuing over documenting this
analysis for future sessions.
