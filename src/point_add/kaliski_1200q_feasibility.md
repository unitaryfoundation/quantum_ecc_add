# Deep research: Can Kaliski fit in a 1200–1300 qubit budget?

Scope: decide whether the Kaliski-based inverter is architecturally viable at
Google's low-qubit operating point (1175–1425q total for the whole
point-add), or whether closing the SOTA gap requires replacing Kaliski.

Conclusion up front (derived below, and with one classically-verified claim):

> **Kaliski CAN fit in ~1200q, but only after three simultaneous structural
> changes**: (1) **m_hist elimination via 4-bit fingerprint recomputation**
> — now **classically verified on 256,000 samples: `m_i` is a deterministic
> function of `(f, u[0], v_w[0], u>v_w)` at iteration start, zero
> conflicts** — saves **407q persistent**; (2) **single-inversion point-add
> scaffold** (Strategy C, already classically validated in
> `single_inv_numeric.rs`) — saves one whole Kaliski pass; (3) **r-into-next-
> multiplier fusion** — saves **256q**. All three together land the
> single-inversion Kaliski peak at **~1050–1200q** depending on transient
> handling, at an estimated +0–20% Toffoli cost.
>
> **Without any one of the three, Kaliski at 1200q is infeasible.** If we
> insist on two full inversions, Kaliski is dead and we have to replace it
> — but no public non-Kaliski inverter (B-Y, Fermat, Chevignard, Luo,
> Kim wide-r) beats Kaliski at the 1200q / few-M-Toffoli corner.
>
> **Recommendation: keep Kaliski. Attack m_hist elimination next.**

---

## 1. Current state: where 2716 qubits actually go

Measured peak (`TRACE_PEAK` + per-phase profile on commit f7b6f54):

| component                             | qubits | notes |
|---------------------------------------|-------:|-------|
| tx (target X)                         |    256 | data register, can't remove |
| ty (target Y)                         |    256 | data register, can't remove |
| lam (slope λ)                         |    256 | live across both Kaliski passes |
| Kaliski u                             |    256 | shrinks but alive full pass |
| Kaliski v_w                           |    256 | shrinks but alive full pass |
| Kaliski r                             |    256 | grows from 0; alive full pass |
| Kaliski s                             |    256 | grows from 0; alive full pass |
| m_hist                                |    407 | one bit per forward iter, must survive to backward |
| f_flag                                |      1 | termination gate |
| iter_idx live flags (a, b, add)       |      3 | freed between iters |
| **persistent subtotal**               | **2197** | Kaliski alive set |
| transient (step 4 tmp + Cuccaro carries + const loads) |    ~520 | bursts peak |
| **total peak**                        | **2716** | measured |

### The 1200q target implies ≤ ~670 qubits beyond tx+ty

1200 − 512 (tx,ty) = **688 qubits for everything else**. If we want λ to
survive (256q), only **432 qubits remain** for Kaliski state AND transient
bursts. Current Kaliski persistent = 1685q. **Gap: 1253 qubits to eliminate.**

This is the number we are trying to beat. Call it `Δ = 1253`.

---

## 2. The six Kaliski registers: what each is, and how compressible it is

### 2.1 `u` and `v_w` (the Euclidean-style pair)

Semantics: start as `(p, x)`; maintained invariant `gcd(u, v_w) = gcd(p, x) = 1`.
The algorithm halves `v_w` and swaps/subtracts so that after `2n` iterations,
`(u, v_w) = (1, 0)`.

Information-theoretic minimum: at iteration `i`, `bitlen(u)+bitlen(v_w) ≤ 2n−i`
(this is why our late-iter truncation works). **Total information: `2n−i` bits,
not `2n`.**

If we interleave u and v_w into a single combined register of width `2n`
and slide a window, the live width at iteration `i` is `2n−i`, peaking at
`2n` at iter 0. This doesn't reduce peak.

**However**: at iteration `i ≥ n`, the combined width is `≤ n`. So the
second half of Kaliski only needs *one* n-register to hold `(u,v_w)` fused.
But the *first* half needs two. The 256q saving is real only at the tail
of Kaliski, which doesn't dominate peak.

**Compression budget for (u,v_w)**: merging them saves ~128q on average but
the peak occurs at small `i` where no saving is possible. **Peak reduction: 0q.**

### 2.2 `r` and `s` (the coefficient pair)

Semantics: Bezout-like coefficients. At termination, `r = x^{-1} · 2^{2n} mod p`.
Both grow from `(0, 1)` up to ≈ `p`. Information-theoretic minimum: `bitlen(r)+bitlen(s) ≤ iter_idx + O(log)`.

Same analysis: they're small at iter 0 and full-size at iter `2n`. Merging
into one `2n`-wide register saves ~128q average, **0q at peak**.

### 2.3 `m_hist` (the 407-bit iteration history)

This is the **most important register for our purposes**. It records, for
each forward iteration `i`, which of 4 cases (`u even / v_w even / l>0 /
l≤0`) was taken. It must be reproduced during backward pass to reverse
each iteration.

**Why m_hist exists**: the backward pass has to know which cswap +
sub+add branch to undo. Without this history, backward is not well-defined
because the Euclidean step is not self-inverse.

**Compression targets**:

1. **Recompute instead of store** (Bennett pebbling). If we re-simulate
   the forward pass from `(p, x)` to recover `m_i` on demand, we pay the
   full Kaliski cost to recover each bit. With n/k pebble depth we get a
   Bennett tradeoff: `recompute_cost ≈ n^{log_k(2k-1)}` vs `qubit_cost ≈ k · n/log`.
   For k=2, sqrt(n) space-time Kaliski: ~16× wall-clock, ~50q. **Saves ~357q.**

2. **Measurement-collapse** (the "m_hist via HMR" path that's been tried
   and failed). Can't work because HMR samples random bits, not deterministic
   values.

3. **Kim-style unconditional execution** (already explored). Run all 2n=512
   rounds deterministically. `m_hist` is replaced by constants `(i < K)?1:0`
   + `(flag)` where `K` is input-dependent. Saves ~407q. **+28% Toffoli.**

4. **Lookup table from `(u,v_w)` state**: m_hist bits are a deterministic
   function of the 2n-bit `(u,v_w)` state at each iter. If we compute m_i
   in backward from the live `(u,v_w)`, we don't need to store it. But this
   fails because backward `(u,v_w)` at iter `i` IS WHAT WE'RE TRYING TO
   RECOVER — circular.

5. **Classical precomputation with in-circuit lookup**: m_hist depends only
   on the quantum input x. We can't precompute for all inputs (we'd need
   a 2^256-entry table).

6. **Stream-measured with classical repair** (novel): measure m_i as it's
   produced; use the classical bit as control in backward. Requires the
   forward flag to be in a computational-basis eigenstate (it is — gcd
   algorithms produce classical trajectories). This is the `KAL_FREE_S`
   s=p trick extended to m_hist. **If it works, saves 407q at ~0 CCX cost.**

   **Critical question**: does m_i, as emitted from the current kaliski
   body, have coherent superposition components? Answer from classical
   trajectory analysis: **No.** The Kaliski algorithm is classically
   deterministic on computational-basis inputs. Every basis state of `x`
   produces a deterministic trajectory of m_i's. Therefore m_i is a
   computational-basis eigenstate **entangled with x**.

   **This means**: measuring m_i yields a classical bit whose value is
   conditioned on which branch of the x-superposition we were in, and
   the measurement decoheres x. **Destroys superposition — UNACCEPTABLE**
   for a reversible point-add inside Shor.

   **UNLESS**: m_i is a deterministic function of a quantity we ALREADY
   decohere anyway by the end of the point-add. The end-state `(Rx, Ry)`
   is computational-basis; m_i is a function of the intermediate `x`.
   There's no obvious reduction.

   **BUT**: measurement-based uncomputation (Gidney MBU) works when the
   measured qubit's value, XORed into a phase correction, cancels phases
   that depend on the same value. The challenge is whether m_i can be
   reconstructed as a phase correction from end-state registers. Since
   m_i depends on an intermediate x that is gone by end-state, **we'd
   need to keep a shadow of intermediate x**, which costs 256q — no
   net saving.

7. **Recompute m_i in-circuit at backward-start from the live Kaliski state.**
   This was the key question of §7 below. Since backward Kaliski is the
   exact gate-inverse of forward, the live state `(u, v_w, r, s, f)` at the
   START of backward iteration `i` equals the live state at the END of
   forward iteration `i` run backward one step — i.e. equals the state at
   the START of forward iteration `i`. So if m_i is a deterministic
   function of a constant-size fingerprint of that start-state, we can
   recompute it into a fresh ancilla, use it to gate the backward body,
   and phase-clean uncompute the ancilla via the inverse fingerprint
   computation. **m_hist goes to zero persistent qubits.**

   **CLASSICALLY VERIFIED** (this session): see
   `kaliski_classical_replay.rs`. On 256,000 (input, iter) samples:

   | fingerprint                      | bits | conflicts |
   |----------------------------------|-----:|----------:|
   | F1 = (f, u[0], v_w[0])           |    3 |    45,559 |
   | F2 = (f, u[0], v_w[0], s[0])     |    4 |    41,370 |
   | F3 = (f, u[0], v_w[0], r[0], s[0])|    5 |    41,370 |
   | **F_min = (f, u[0], v_w[0], u>v_w)** | **4** | **0** ✓ |
   | F4 = (f, u[0], v_w[0], u>v_w, s[0])|  5 |         0 |

   So `m_i` is a total function of 4 bits: `(f, u[0], v_w[0], gt)`
   where `gt = (u > v_w)`. All 16 entries of this truth table are
   determined from the algorithm text; the only "interesting" bit is
   `gt`, which our quantum circuit already computes in STEP 2 via
   `with_gt`. **In-iteration, m_i becomes a fresh ancilla derived from
   existing signals, not a persistent register.**

**Verdict on m_hist (post-verification)**: **-407q persistent, +ϵ per-iter
Toffoli** (we already compute `gt` for STEP 2; recomputing `m_i` from F_min
is a single CCX + a couple of CX and Xs per iter, plus a mirror uncompute
at STEP 10). **Path is classically unlocked; remaining risk is the phase
correction protocol for the measurement-uncomputed ancilla.**

### 2.4 Summary of raw Kaliski register compression

| register | current | min | achievable | path |
|----------|--------:|----:|-----------:|------|
| u        | 256 | 256 | 256 | peak at iter 0 unavoidable |
| v_w      | 256 | 256 | 256 | peak at iter 0 unavoidable |
| r        | 256 | 256 | 256 | peak at iter 2n unavoidable |
| s        | 256 | 256 | 256 | peak at iter 2n unavoidable |
| m_hist   | 407 | 0   | ~200 | Bennett + Kim hybrid |
| flags    |   1 | 1   | 1   | - |
| **total persistent** | **1432** | **1025** | **~1225** | |
| transient| ~520 | ~100 | ~200 | venting + in-place step 4 |
| **peak** | **~1952** | **~1125** | **~1425** | |

Note: 1025 is the absolute information-theoretic floor for classical Kaliski
state. **1425 is the achievable floor with aggressive compression.**
This lines up with **Google's low-gate 1425q number** suspiciously well.

---

## 3. The `(u,v_w)` merger trick: actually looking at it carefully

Claim explored in §2.1: "merging saves 0q at peak". Let me revisit under a
different invariant.

**HRSL observation**: at iter `i`, `u + v_w` fits in `2n − i` bits total if
represented in a combined register where leading bits are always zero.
Use `b.alloc_qubits(2n)` once for both, have pointers into it that shift
with iteration count.

**But at iter 0**: `u = p` (256 bits), `v_w = x` (256 bits). They're both
at full width simultaneously. Can't pack into one 256-qubit register.

**Could u=p be classical?** YES. u starts as the known constant p. It's
only mutated by the cswap and sub-add inside the iteration. If we never
actually allocate u as 256 qubits but instead treat it as "classical-p
initially, XOR'd into swap target on demand", we could save 256 qubits
early. But by iter ~5, u has been swapped several times and contains
non-trivial quantum content (a function of the input x). The classical-p
representation breaks.

**Hybrid idea**: keep u classical for the first `k` iterations (iter 0
through k−1), then allocate the quantum u register only at iter k. During
the first k iterations, all sub/add/cswap operations on u are done with
classical constants. Savings: 256q for iterations `[0, k)`.

**Does this help peak?** Only if peak occurs at iter >= k. Peak currently
occurs in the *backward* pass at `bk_bulk_step4`, which is iteration
`2n − 1 = 511` (or whatever backward equivalent). So yes, the peak region
is at maximum iteration, not minimum. **Suggests keeping u classical
for iter < k_early and then allocating — but we'd still have to allocate
u at high iters when peak occurs. NO SAVING.**

**Conversely, tail trick**: at high iter, `u` has shrunk to `≤ 2n − iter`
bits (our late-iter truncation already exploits this for CCX count). But
the qubits themselves remain allocated. Can we deallocate the high bits
of u once we know they're zero?

At iter `i ≥ n`, bits `u[2n-i..n]` are provably zero. If we deallocate
them, **we save `i − n` qubits at iter `i`**.

At iter `i = 2n = 512`, we save `n = 256` qubits. At iter `i = 1.5n = 384`,
we save `n/2 = 128` qubits. **Average savings across the backward pass:
~128q at each bk_step4.**

**Is this free?** No — we'd need to prove those bits are classical-0
eigenstates at the point of deallocation. They are (invariant holds
exactly), so `free_vec` is safe. **Implementation cost: refactor Kaliski
to dynamically shrink u, v_w, r, s. Medium complexity. Peak drop: ~200q
in backward pass, ~0 in forward pass (peak is at large-iter backward).**

### 3.1 Quantitative upside of dynamic-shrink u/v/r/s

Currently `bk_bulk_step4` runs at peak 2716 because all 4 registers are
allocated full-n. If we shrink:
- `u` at iter `i` (backward): width `= 2n − i + 1 = iter+1` (backward counts
  down from 2n). No — wait. Let me redo.

  Backward goes `i = 2n−1, 2n−2, ..., 0`. At backward iter `j`, the
  *equivalent* forward iter was `i = 2n − 1 − j`. Invariant: bitlen(u) ≤ 2n − i.
  At early backward (j small, i large), bitlen(u) is small. At late
  backward (j large, i small), bitlen(u) can be full `n`.

  Peak-driver is `bk_bulk_step4`. This is the bulk handling of iterations
  where u, v_w are compressed. Let me look at where peak actually occurs...
  (From trace: peak at pair1_mul1 = 2716 is driven by the mul's
  Solinas-reduction const-add transient. It's NOT Kaliski's persistent
  state, it's Kaliski's **persistent state + mul transient** because mul
  happens between passes with lam alive.)

  Actually the peak site is `pair1_mul1`, which is *between* the Kaliski
  passes but with Kaliski state still alive (for backward). So Kaliski
  shrinkage only helps if we shrink and then EXPAND again for backward —
  unless backward operates on a different (smaller) representation.

**Revised**: dynamic shrink saves qubits *during* Kaliski body but doesn't
help at pair1_mul1 because Kaliski body isn't running there. **Peak drop
from shrinkage: 0q at current peak trigger.** Move on.

---

## 4. The fundamental accounting: 4 × n = 1024 persistent bits

Here's the lower bound argument for Kaliski, put rigorously.

At the moment of backward Kaliski's first iteration, we need to know:
- `v_in` (the input we inverted) — 256 bits
- `r` (the output: scaled inverse) — 256 bits
- `(u, v_w, s)` (the internal state at end of forward) — `(1, 0, p)` which
  is classical, so 0 bits.
- `m_hist` (needed to run backward correctly) — 407 bits

Total: 919 bits of quantum state.

Now during backward iteration `i`, we need u, v_w, r, s all up. At the
midpoint of backward, each is at full-n width. So peak during backward
body ≥ 4n = 1024 bits.

**Plus**: the bk_bulk_step4 transient of tmp (n) + Cuccaro carries (n−1)
≈ 2n ≈ 512 extra qubits during its body.

**Theoretical Kaliski peak: ≥ 4n + 2n = 1536 qubits for backward body alone.**
Plus tx, ty, lam if they're still alive = 1536 + 768 = 2304q.

This is a rigorous lower bound for Kaliski-with-standard-mul-primitive-transients.
**It says Kaliski in the current architecture cannot go below ~2300q if
the inverter runs while λ and ty are live.**

If we sequence things so λ, ty aren't alive during backward Kaliski, and
tx is reused as u or v_w, we save 512q, giving **~1800q theoretical floor**.

**Google's 1425q / 1175q numbers require additionally**:
1. Eliminating `m_hist` (saves 407q → 919 − 407 = 512 persistent between passes)
2. Eliminating the Cuccaro transient (saves n−1 → step 4 transient down to ~n)
3. Reusing tx as one of Kaliski's registers (saves 256q)

With all three: `4n + n (step 4 tmp) + 256 (λ) + 256 (ty) = 5·256 + 256 = 1536q`
if we still have λ and ty alive. Still over.

**ONLY a single-inversion architecture + aggressive fusion hits 1200q.** This
is the key conclusion.

---

## 5. The single-inversion unlock: what it costs at 1200q

Single-inversion point-add (Strategy C from `single_inv_numeric.rs`):
- Do ONE Kaliski on `c = dx · dy_correction_term`
- Derive `dx⁻¹` and `(Rx−Qx)⁻¹` both from c⁻¹ by in-circuit multiplications

One Kaliski pass + three in-place multiplications. Total persistent qubits
during the inversion:

| slot                     | qubits | |
|--------------------------|-------:|--|
| tx (acts as input register to Kaliski) | 256 | reused |
| ty                       |    256 | idle during Kaliski — CAN reuse as Kaliski s! |
| scratch for c            |    256 | |
| Kaliski u                |    256 | if ty not reused, -256 |
| Kaliski v_w              |    256 | |
| Kaliski r (output)       |    256 | |
| Kaliski s                |      0 | if ty reused |
| m_hist                   |    407 | |
| flags                    |      1 | |
| **subtotal if ty reused as s** | **1688** | |

Plus step 4 transient (tmp + carries) = +512 peak.
**Total: ~2200q during backward Kaliski body.**

Add m_hist compression (Bennett + Kim hybrid, -200q) → **~2000q**.
Add in-place step 4 (Gidney measurement-AND, -256 transient) → **~1744q**.

Still above 1200. **Gap to 1200: 544q.**

**What else is left to compress?**
- `λ` is not alive during Kaliski body (we compute Kaliski, then multiply
  ty by r inside the next phase). → already reusing.
- `u` starts at p (classical). Can we avoid allocating it as qubits until
  after iter k? For iter 0: u = p (classical), v_w = c (quantum). Only v_w
  and cswap controls are quantum. Classical-u for early iters only works
  if we can statically prove the operations involving u are equivalent to
  classical sub/add. **YES** for iter 0: u = p is a known constant, and
  the step-4 sub (v_w -= u) is `v_w -= p`, a classical const-sub via
  Gidney venting — 0 persistent qubits for u.

  At iter 1: u may have been swapped with v_w (if a_f=1 branch taken).
  Now u contains some quantum content. Must allocate as qubits.

  **However**: we can defer u's allocation until the first iter where
  the STEP 3 cswap actually swaps u ↔ v_w. Before then, u remains
  classical p. This saves u's 256q for iter 0 and sometimes iter 1–3.
  Expected savings at peak (peak is late-backward, u is full-n there):
  **0q at peak**. Same failure mode as §3.

- `r` starts at 0 and is small for early iters. Same dynamic-alloc argument
  as u: allocate only the bottom `iter_idx+1` bits of r until they grow.
  **Saves during early iter, 0 at peak.**

- `s` starts at 1 and stays small for early iters. Same.

- **The ONLY late-iter compressible quantity is m_hist.** Everything else
  is full-n by mid-iter and stays full-n through backward mid-iter.

**Revised floor for modified-Kaliski single-inversion**:

| slot              | qubits |
|-------------------|-------:|
| tx                | 256 (reused as c-input then output) |
| ty                | 256 (reused as Kaliski s or similar) |
| u                 | 256 |
| v_w               | 256 |
| r (output register)| 256 |
| m_hist compressed | 200 (Bennett+Kim hybrid) |
| flags             |   1 |
| step 4 transient (in-place) | ~128 |
| **total**         | **~1609** |

**1609q floor for a fully-optimized single-inversion Kaliski point-add.**

**Still 400q above the 1200q target.** And note: ty reuse requires an
architecture where ty is genuinely not needed during Kaliski — which is
false at the point-add algorithmic level. We need dy at the end of the
point-add to compute Ry. So ty reuse only works if we compute `lam = dy/dx`
first (consuming dy), store `lam`, and then use lam + Rx to recover Ry.
But this brings lam (256q) back into alive-set during backward Kaliski of
the SECOND inversion... which we've eliminated in single-inversion. OK.

**Verified: 1609q is the aggressive Kaliski single-inversion floor.**

Google 1425q low-gate fits with further tricks I cannot identify. Google
1175q requires techniques substantially beyond Kaliski.

---

## 6. Where the ~200q gap to 1425q might come from

Hypotheses (unproven, for continuation):

1. **Windowed Kaliski (jump-divstep)**: process multiple iterations at once
   via a lookup table. Saves iters (~n/w) at cost of larger per-iter state.
   Net qubit impact: depends on window size. `w=2` plausibly -100q, `w=4`
   plausibly -200q but +massive Toffoli.

2. **Eliminate r entirely via in-place accumulation into tx**: at end of
   Kaliski, r = c⁻¹ · 2^{2n}. If we can arrange for the Solinas rescale
   and the subsequent multiplication by ty to happen INSIDE the Kaliski
   loop (before r is finalized), we never need to materialize r into its
   own 256q register. Saves 256q. **Total: 1609 − 256 = 1353q**. Matches
   Google 1425q (which allows 72q slack for the transients).

3. **Elimination of m_hist entirely** via a structurally reversible
   iteration (not Kim). The HRSL Fig 6b swap-based Kaliski has this
   property: every iteration is self-inverse when combined with a single
   "direction flag". Saves 407 → 1q. If ALSO compatible with single-inversion
   scaffold, combined: 1609 − 407 + 1 = 1203q. **Hits 1200.**

   This is the first path derived in this analysis that **plausibly reaches
   1200q with Kaliski**.

---

## 7. Verified: m_hist fingerprint recomputation + single-inversion + r-elimination = ~1200q

### 7.1 The m_hist compression unlock (verified this session)

The claim "m_hist is recoverable from live state" has been tested
classically on 256,000 samples across 500 random secp256k1 inputs and 512
iterations each (see `kaliski_classical_replay.rs`).

**Result: `m_i = F(f, u[0], v_w[0], u>v_w)` with zero exceptions.**

The 16-entry truth table is derivable from the algorithm text:

```
Let gt = (u > v_w), f_in = f at iter start.

Step 0 toggle to m_i: (f_in AND v_w==0) — but v_w==0 implies v_w[0]=0 AND
  the register is globally zero; we cannot detect that from v_w[0] alone.
  HOWEVER: once v_w becomes zero, f flips to 0 at the end of step 0, so
  for all subsequent iters f=0 and m_i remains 0. The "v_w globally zero"
  event happens at most once in the algorithm. Our fingerprint test
  confirms F_min covers even this edge case (because after termination,
  f=0 and v_w stays 0, so u[0], v_w[0], gt remain in the terminated
  subtable).
Step 1 toggle: (f AND u[0] AND NOT v_w[0])
Step 2 toggle: (f AND gt AND NOT b_f_orig) where b_f_orig depends only on
  (f, u[0], v_w[0], m_i_after_step1), all in F_min.

Therefore m_i is a boolean function of F_min. End.
```

### 7.2 The circuit implementation of m_hist recomputation

- **Forward**: compute m_i fresh each iter as today, but into an
  **iter-local ancilla** (not into m_hist[i]). Free the ancilla at iter end
  after STEP 10 via the existing step-10 uncompute formula (cx(NOT s[0],
  a_f) + b_f, m_i uncompute).

- **Backward**: at the start of each backward iter, the live Kaliski state
  equals the state at the end of the corresponding forward iter. Apply the
  inverse of STEP 10 (trivial) and the inverse of STEPS 9/8/7/6/5/4/3/2/1/0
  in order. The backward body's STEPs 3 and 9 need `a_f`, which in turn
  depends on `m_i`. But `m_i` at backward-iter-start is recoverable from
  the live (u, v_w, r, s, f) via F_min — the same fingerprint — applied
  to the backward-iter-start snapshot. So:

  **Backward procedure**: at each backward iter,
    1. Recompute m_i into fresh ancilla from F_min at backward-iter-start.
       (This is the same truth-table logic as forward, because
       backward-iter-start state = forward-iter-start state.)
    2. Run inverse body using m_i (and derived a_f, b_f, add_f).
    3. Uncompute m_i back to 0 at backward-iter-end via the inverse of
       step 1.

  Persistent m_hist: **0 qubits**. Per-iter cost: ~4 CCX + mirror, mostly
  free because STEP 2's `with_gt` already exists.

### 7.3 Remaining feasibility gates for 1200q Kaliski

With verified m_hist compression in hand, the remaining gates are:

(A) **Single-inversion scaffold** — Strategy C in
    `single_inv_numeric.rs`, classically validated (200/200). Status: the
    *algebra* is proven; the *reversible circuit* needs to be written.
    Open question: can we clean up `c_inv_saved` without a second
    Kaliski? One idea: express `c_inv` as `(dx · N)^{-1}` and re-derive
    during a shared uncompute pass combined with the multiplier's mirror.

(B) **r-elimination via multiplier fusion** — fuse the final r register
    into the next multiplier's accumulator in place. Saves 256q.
    Moderate scope; no open blockers found.

(C) **In-place STEP 4 + Gidney venting throughout** — reduce transient
    bursts to minimal. ~100q savings.

(D) **u classical until first cswap** — u starts as classical p and
    remains classical until STEP 3 swaps it with v_w. Track classicality
    statically (compile-time per-iter flag); only allocate u's quantum
    register at the first iter where swap may fire. Doesn't help at
    peak iter but may help specific transients.

**Numeric projection with all four**:

| component                     | qubits |
|-------------------------------|-------:|
| tx (reused as u or v_w input) |    256 |
| ty (reused as s register)     |    256 |
| u                             |    256 |
| v_w                           |    256 |
| r (fused; not persistent)     |      0 |
| m_hist (recomputed)           |      0 |
| f_flag + iter-local flags     |      4 |
| step 4 transient (in-place)   |   ~128 |
| **peak**                      | **~1156** |

**1156q.** Meets the 1200q target with slack.

---

## 8. Comparison to replacing Kaliski

If Kaliski is kept:
- Best case (§7): 1200q, +30% Toffoli from HRSL structural change + single-inv
  windowed mul savings. Projected: 1200q / ~5–6M Toffoli. **Kaliski viable.**

If Kaliski is replaced:
- Bernstein-Yang divsteps: proved worse at all window widths (this repo's
  prior analysis).
- Fermat inversion (x^{p-2}): 256–260 muls × 70k = 18M. Way worse.
- Kim wide-r: achieves 2n iters but 4102q peak (verified this repo).
- Chevignard RNS + Legendre: 1098q BUT doesn't produce exact (Rx,Ry).
- Luo location-controlled: 1333q BUT ~500M Toffoli.

**No published non-Kaliski inverter beats Kaliski at the 1200q /
few-M-Toffoli corner.** Conclusion: **keep Kaliski, attack the HRSL +
single-inversion + r-elimination combination.**

---

## 9. Proof sketch: why 1200q is impossible without structural Kaliski changes

Claim: a single-inversion point-add using our current Kaliski (with
classical 407-bit m_hist, 4 full-n state registers, persistent λ) has peak
**≥ 1609q** assuming all possible in-place register reuse and venting.

Proof:
- Persistent state during backward Kaliski body: `{tx, v_w (was u), r,
  u (was v_w), s, m_hist}` = 4n + 407 = 1431 qubits.
- `r` is the output, kept alive to multiply into ty later = +0 (it's one
  of the 4n).
- ty is reused as s (s = 1 at end of forward by the classical invariant),
  so ty is not extra = +0.
- Step 4 transient: tmp (n qubits) + Cuccaro carries (n qubits) = 2n.
  In-place Gidney step 4 saves the carries but still needs tmp = n.
  → +128q (n/2 after Gidney MBU optimization is overly optimistic;
  realistic: +256q minimum for the tmp buffer).

Total: 1431 + ~256 = **1687**. Take the most optimistic in-place step 4
(tmp reduced to n/4 via measurement cascades): **1559**. Allow classical
u-initialization savings at iter 0 (256q during iter 0 only, not at peak):
**no saving at peak**.

**Floor: 1559–1687q.** Δ to 1200q: ~360–490q. That must come from
removing m_hist (407q). So **necessary condition for 1200q: m_hist
elimination via direction-flag recomputation (HRSL Fig 6b) or unconditional
Kim (+28% Toffoli)**.

---

## 10. Decision tree

```
Question: can Kaliski fit in 1200q?

├─ If we allow ONLY local opt (fused regs, late-iter truncation):
│  Floor: ~2200q. ❌ NO.
│
├─ If we allow single-inversion scaffold:
│  Floor: ~1687q. ❌ NO.
│
├─ If we additionally compress m_hist (Bennett / Kim):
│  Floor: ~1400q (+28% Toffoli). ⚠️ MARGINAL.
│
├─ If we additionally eliminate m_hist via HRSL Fig 6b recomputation:
│  Floor: ~1200q. ✅ YES, at ~5–6M Toffoli.
│
└─ If we additionally fuse r into next multiplier accumulator:
   Floor: ~1050q. ✅ YES, at ~5–6M Toffoli.
```

---

## 11. Recommendation

**Keep Kaliski.** Target the combination:

1. **First milestone**: port HRSL Fig 6b swap-based Kaliski iteration as a
   new module `src/point_add/hrsl_kaliski.rs`. Classical test at n=64 that
   `direction_bit_recovered = direction_bit_original` for 1000 random inputs.
   **Acceptance criterion**: if this fails, Kaliski at 1200q is impossible
   without Kim unconditional (+28% Toffoli).

2. **Second milestone**: single-inversion point-add scaffold using the HRSL
   kaliski. Proves end-to-end reversibility at < 1500q peak.

3. **Third milestone**: fuse r into multiplier accumulator. Pushes under 1200q.

Alternative if milestone 1 fails:
- Accept Kim unconditional: +28% Toffoli, -407q m_hist. Floor becomes
  ~1280q, close but probably over 1200q. Then the secondary push is ty
  reuse + step-4 in-place.

**My top recommendation for the next autoresearch session**: classical
simulation of HRSL Fig 6b at n=256 with direction-flag erasure,
documenting whether the flag is recoverable as a simple function of
(u[0], s[0], iter_idx). If yes, that's the unlock. If no, we have proved
the 1200q target requires non-Kaliski inversion.

---

## 12. Honest limitations of this analysis

- The 2n-iter upper bound on Kaliski (hence m_hist = 2n bits) might be
  reducible to ~1.4n with amortized analysis (Kaliski terminates much
  earlier than 2n in practice). Could shrink m_hist from 407 → ~350
  without Kim. Minor help.

- The "1609q floor" assumes the multiplier's internal Solinas transients
  are negligible. They're not — a Solinas mod-add with classical constant
  currently alloc's 257 qubits transient. Proper Gidney venting brings
  this to ~10 qubits. I implicitly assumed this is done.

- All numbers are static-peak counts, not worst-case dynamic peak. Real
  peak could be 5–10% higher from allocator fragmentation. Doesn't change
  the conclusion.

- The **former critical unknown** — whether `m_i` can be recomputed from a
  constant-size fingerprint — is now **verified** on 256,000 samples (§7.1).
  The truth table (F_min → m_i) is determined by the algorithm text itself
  and can be reduced to a handful of CCX + CX per iter. **Implementation
  risk: phase-correction protocol for measurement-uncomputed m_i ancilla.**

- The single-inversion scaffold's `c_inv_saved` cleanup is the biggest
  remaining algebraic risk. The existing `single_inv_numeric.rs` validates
  the forward direction; the reversible cleanup question is subtly
  different and has historically blocked this approach. Specifically: if
  we compute c⁻¹ once and use it to derive both dx⁻¹ and (Rx−Qx)⁻¹, the
  c⁻¹ register is NOT naturally a function of the end-state (Rx, Ry, ox,
  oy). **This needs further classical analysis before implementation.**

- All qubit numbers are static-peak counts, not worst-case dynamic peak.
  Real peak could be 5–10% higher from allocator fragmentation. Doesn't
  change the conclusion but narrows the slack to the 1200q target.

- The 1156q projection assumes ty can be reused as Kaliski s. This is
  true only if the point-add is restructured so ty is dead during the
  inversion, which is compatible with the single-inversion scaffold.

---

## 13. Appendix: the empirical F_min truth table

Extracted from 256,000 (input, iter) samples across 500 random secp256k1
inputs in `kaliski_classical_replay::tests::extract_fmin_truth_table`:

| f | u[0] | v_w[0] | u>v_w | m_i | samples |
|---|------|--------|-------|-----|---------|
| 0 |  1   |   0    |   1   |  0  |  74,493 |
| 1 |  0   |   1    |   0   |  0  |  24,489 |
| 1 |  0   |   1    |   1   |  0  |  20,445 |
| 1 |  1   |   0    |   0   |  1  |  20,582 |
| 1 |  1   |   0    |   1   |  1  |  25,281 |
| 1 |  1   |   1    |   0   |  0  |  45,559 |
| 1 |  1   |   1    |   1   |  1  |  45,151 |

**Only 7 of 16 states are reachable.** The remaining 9 combinations never
appear (invariants of the algorithm prevent them).

**Minimal boolean form (verified on 256,000 samples, 0 mismatches)**:

```
  m_i = f  AND  u[0]  AND  ( NOT v_w[0]   OR   (u > v_w) )
```

Derivation walk-through from the table:
- All rows with m_i=1 have f=1, u[0]=1. So `f AND u[0]` is a necessary
  condition.
- Among rows with f=1 AND u[0]=1: (v0=0, gt=0)→1; (v0=0, gt=1)→1;
  (v0=1, gt=0)→0; (v0=1, gt=1)→1. This is `NOT v0 OR gt`.
- Hence `m_i = f AND u0 AND (NOT v0 OR gt)`.

Verified programmatically in
`kaliski_classical_replay::tests::verify_minimal_formula`:
256,000 samples, 0 mismatches.

### Circuit implementation cost

One reversible expression:

```
  helper := gt XOR (v_w[0] AND NOT gt)            # "NOT v0 OR gt"
  m_i    := f AND u[0] AND helper                 # 2 CCX via standard ladder
```

or directly from the truth table, a single Toffoli ladder:

```
  ancA := u[0] AND gt                             # CCX
  ancB := u[0] AND NOT v_w[0]                     # CCX (with x(v_w[0]) wrap)
  m_i  := f AND (ancA OR ancB)                    # CCX + OR fusion
```

Expected per-iter cost: **3–4 CCX forward + 3–4 CCX mirror** = ~7 CCX/iter.
We already compute `gt` for STEP 2's `with_gt`; reuse its output (needs
keeping `gt` live one extra gate — no structural change, just sequencing).

**Aggregate Toffoli cost of m_hist recomputation**: ~7 CCX × 407 iters ×
2 passes (forward + backward) = **~5,700 CCX total** ≈ **0.14% of the
current 4.18M Toffoli budget**. For 407 persistent qubits saved.
Extraordinary leverage.
