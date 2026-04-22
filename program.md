# quantum_ecc research loop

You are an autonomous research agent optimizing a reversible quantum circuit
for secp256k1 point addition. Your job is to iteratively reduce the
**average executed Toffoli count** while keeping the circuit correct and
within a qubit budget. Run continuously. Do not pause for human confirmation.

## Scope of edits

- You may ONLY modify files under `src/point_add/` (the whole folder is
  yours — add/rename/split submodules freely).
- You may NOT modify `src/main.rs`, `src/circuit.rs`, `src/sim.rs`,
  `src/weierstrass_elliptic_curve.rs`, `Cargo.toml`, `Cargo.lock`,
  `rust-toolchain`, `results.tsv` (directly) or anything else.
- You may NOT add dependencies.
- You may NOT modify the test harness or the correctness check.

## Objective

Minimize the metric `avg executed Toffoli` printed by `cargo run --release`.

### Hard constraints (run is invalid if violated)
1. `=== experiment OK ===` must print. This requires:
   - all 64 classical correctness shots pass, AND
   - `strict_apply` passes — every `R` (i.e. every `assert_zero_and_free`)
     targets a qubit whose 64-shot value is already 0, AND
   - the forward∘reverse identity check passes — after running the
     circuit and then its gate-reversed inverse, every qubit returns to
     its pre-forward snapshot.
2. `qubits` (peak live) must be ≤ **2800** (≈ current baseline).
   We need to reduce qubits over time to the below results; never exceed the current best's
   qubit count by more than 5% unless the Toffoli win is >10%.
3. `cargo build --release` must succeed with no warnings introduced by your
   edits beyond those already present on the baseline.

### Reversibility

Every ancilla must be uncomputed to |0⟩ before being freed. The standard
pattern is compute / use / uncompute. The harness enforces this two ways:

- `sim.rs` treats every `R` op (`Builder::assert_zero_and_free`) as a hard
  assertion that the target qubit is |0⟩ on every live shot. Dirty frees
  fail at the dirty op with a localized error.
- After the forward pass, the harness zeroes the output registers and
  asserts every remaining qubit is |0⟩. Lingering ancillas anywhere
  outside the four declared registers fail this check.

There are no loopholes — a Toffoli "win" from skipping uncomputation
makes the run fail, not faster.

### Tie-breakers (when Toffoli counts are within ~0.5%)
- Lower peak qubits.
- Lower total Clifford.

Code aesthetics are NOT a consideration. Long functions, hand-unrolled
loops, duplicated primitives, weird control flow, deeply nested
`emit_inverse`, hundreds of call-site-specialized helpers — all fine if
they shave Toffolis. The goal is the best-possible circuit, not the
cleanest codebase. If in doubt between a 500-line readable version and
a 2000-line ugly version at lower Toffoli count, pick the ugly one.

## Baseline (honest reversible kaliski, commit `main`)

```
avg executed Toffoli  : 101284162
avg executed Clifford : 211257273
emitted ops           : 383933667
qubits                : 3595
```

Reference targets (zenodo `zkp_ecc` Pareto frontier, for calibration —
these are aspirational, not required):

| Variant | Toffoli | Qubits |
|---|---|---|
| low-qubit | 2,700,000 | 1,175 |
| low-gate  | 2,100,000 | 1,425 |

You are ~40× above these on Toffoli and ~3× over on qubits. There is
substantial room.

## Setup

On first run only:
1. `git checkout -b autoresearch/<YYYY-MM-DD>` — work on a dated branch.
2. Read `src/point_add/mod.rs` and the module doc at its top
   (steps 1–12 of the point-add algorithm).
3. Skim `src/circuit.rs` for the `Op` IR and `src/sim.rs` for how gates
   are counted (in particular `sim.rs:102` — `executed_shots` semantics).
4. Verify the baseline runs: `cargo run --release -- --note baseline` should
   print `=== experiment OK ===` and append a TSV row ending in `OK` to
   `results.tsv`.

## Experiment loop

Repeat indefinitely:

1. **Pick an idea**. Either from the seed list below or your own. Feel free to pursue ideas you gave up on earlier if you reach a bottleneck.
2. **Edit** files under `src/point_add/` to implement it.
3. **Build**: `cargo build --release 2>&1 | tail -20`.
   - If it fails to compile, either fix immediately (if the fix is obvious
     and small) or `git checkout -- src/point_add/` and pick a different
     idea. Do not leave the tree broken.
4. **Run**: `cargo run --release -- --note "short description of the idea"`
   — `main.rs` automatically appends a TSV row to `results.tsv` with
   timestamp, commit, toffoli, clifford, qubits, ops, correct, and your note.
   Both `OK` and `FAIL` runs log a row.
5. **Decide**: read the last row of `results.tsv` (or the printed metrics).
   - If `correct == OK` AND `toffoli < best_toffoli` AND qubits constraint met:
     - `git add -A && git commit -m "<short desc>: toffoli <old> → <new>"`
     - Update your in-memory `best_toffoli`.
   - Else:
     - `git checkout -- src/point_add/` to revert. The TSV row stays;
       it's part of the research log.
6. Go to 1.

Never `git reset --hard` across multiple commits — only revert the current
in-progress edit. Keep every accepted commit.

## results.tsv format

Columns (tab-separated), written automatically by `main.rs`:
```
timestamp    commit    toffoli    clifford    qubits    ops    correct    notes
```
`main.rs` appends one row per `cargo run --release` invocation. The `notes`
column is whatever you pass via `--note "..."`. Tabs and newlines in the
note are stripped. Always pass a note — future-you needs it to interpret
the row.

## Idea seeds

- (Roetteler Naehrig Svore Lauter 2017 — Quantum resource estimates for ECDLP, https://arxiv.org/abs/1706.06752)                                                                                    
- (Litinski 2023 — How to compute a 256-bit elliptic curve private key with only 50 million Toffoli gates, https://arxiv.org/abs/2306.08585)
- (Häner Roetteler Soeken 2020 — Improved quantum circuits for elliptic curve discrete logarithms (eprint), https://eprint.iacr.org/2020/077.pdf)                                                   
- (Häner Roetteler Soeken 2020 — Improved quantum circuits for elliptic curve discrete logarithms (arXiv), https://arxiv.org/abs/2001.09580)                                                        
- (Gidney 2019 — Windowed quantum arithmetic, https://arxiv.org/abs/1905.07682)                                                                                                                     
- (Ragavan Gidney 2025 — Optimized circuits for windowed modular arithmetic, https://arxiv.org/abs/2502.17325)                                                                                      
- (Cuccaro Draper Kutin Moulton 2004 — A new quantum ripple-carry addition circuit, https://arxiv.org/abs/quant-ph/0410184)                                                                         
- (Remaud et al. 2024 — Optimizing T and CNOT gates in quantum ripple-carry adders and comparators, https://arxiv.org/abs/2401.17921)                                                               
- (Banegas et al. 2022 — Reducing the depth of quantum FLT-based inversion circuit, https://eprint.iacr.org/2022/463.pdf)
  
When you get stuck, try your own ideas and do more research.

## Go big

We're trying to **win** against SOTA (2.1M–2.7M Toffoli target).
You are authorized to make high-ambition changes: rewrite primitives
from scratch, swap in a new algorithm, do thousand-line refactors in a single
commit, and spend full sessions on moonshots that might not pan out.
The only immovable constraint is the harness contract.

Key: Micro optimizations are explicitly BANNED. They are a waste of time. Changes must be lead to large gains.

## Rules of thumb

- If a run takes longer than 5 minutes, something is wrong — kill and revert.
- Cliffords are free compared to Toffolis (~100× cheaper in fault tolerance).
  Do not optimize Cliffords at the cost of Toffolis.
- X/Z gates are not counted at all. Abuse them.
- Correctness is non-negotiable. A 0-Toffoli circuit that fails correctness
  is worth nothing. Run `cargo run --release` after every edit.

## Stop conditions

Keep iterating until one of:
- You hit the zenodo low-qubit target (2.7M Toffoli @ ≤1175 qubits).
- The user interrupts.

What is key is that you do not stop until you meet the above. The world is your oyster. LOOP FOREVER.

Never stop. I'm going to bed. Keep working. There is no budget, there are no constraints; we just need a high state of the art.                                                