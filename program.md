# quantum_ecc research loop

You are an autonomous research agent optimizing a reversible quantum circuit
for secp256k1 point addition. Your job is to iteratively reduce the
**average executed Toffoli count** while keeping the circuit correct and
within a qubit budget. Run continuously. Do not pause for human confirmation.

## Scope of edits

- You may ONLY modify `src/point_add.rs`.
- You may NOT modify `src/main.rs`, `src/builder.rs`, `src/circuit.rs`,
  `src/sim.rs`, `src/weierstrass_elliptic_curve.rs`, `Cargo.toml`,
  `Cargo.lock`, `rust-toolchain`, or anything else.
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
2. `qubits` (peak live) must be ≤ **3700** (≈ current baseline).
   Prefer to reduce qubits over time; never exceed the current best's
   qubit count by more than 5% unless the Toffoli win is >10%.
3. `cargo build --release` must succeed with no warnings introduced by your
   edits beyond those already present on the baseline.

### BANNED: dirty ancilla frees

`Builder::assert_zero_and_free(q)` is an ASSERTION, not a deallocation
primitive. It is the gate-level equivalent of `assert q == |0⟩`; calling
it on a qubit that is not genuinely back to |0⟩ is a correctness bug, not
an optimization.

**Rule**: every ancilla you allocate must be returned to |0⟩ by an
explicit inverse gate sequence BEFORE you call `assert_zero_and_free`.
No exceptions.

**Why this matters**: the simulator's `R` gate unconditionally zeros its
target, so a "dirty free" can silently pass classical tests on definite
inputs. But the ancilla is entangled with the computation state, and on
the superposition inputs Shor's algorithm actually uses, that entanglement
destroys the interference pattern the algorithm depends on. The circuit
would be worthless.

**How to apply**: if you compute `anc = f(x, y)` into a fresh ancilla,
you must emit the inverse of `f` (same gates in reverse order, since all
our gates are self-inverse) before freeing `anc`. The standard pattern is
compute / use / uncompute. The harness's `strict_apply` + forward∘reverse
identity check will catch every violation and fail the run. Do not try
to find loopholes; there are none.

**Do NOT** try to route around the assertion by using raw `builder.r(q)`,
never freeing the qubit, using `X` to flip a dirty ancilla back to 0, or
any other dodge. The forward∘reverse identity check catches all of these.

An optimization that appears to reduce Toffolis by skipping uncomputation
is not an optimization — it's a broken circuit that the harness will
reject. Spend your effort finding genuinely better algorithms (see seeds).

### Tie-breakers (when Toffoli counts are within ~0.5%)
- Lower peak qubits.
- Lower total Clifford.

## Baseline (starting point, commit `main`)

```
avg executed Toffoli  : 55399434
avg executed Clifford : 146099482
emitted ops           : 226847650
qubits                : 2574
```

Reference targets (zenodo `zkp_ecc` Pareto frontier, for calibration —
these are aspirational, not required):

| Variant | Toffoli | Qubits |
|---|---|---|
| low-qubit | 2,700,000 | 1,175 |
| low-gate  | 2,100,000 | 1,425 |

You are ~20–25× above these. There is substantial room.

## Setup

On first run only:
1. `git checkout -b autoresearch/<YYYY-MM-DD>` — work on a dated branch.
2. Read `src/point_add.rs`, `src/builder.rs`, and the module doc at the top
   of `point_add.rs` (steps 1–12 of the point-add algorithm).
3. Skim `src/circuit.rs` for the `Op` IR and `src/sim.rs` for how gates
   are counted (in particular `sim.rs:102` — `executed_shots` semantics).
4. Verify the baseline runs: `bash scripts/bench.sh` should print one TSV
   line ending in `OK`. Append it to `results.tsv` as the starting row.

## Experiment loop

Repeat indefinitely:

1. **Pick an idea**. Either from the seed list below or your own.
2. **Edit** `src/point_add.rs` to implement it.
3. **Build**: `cargo build --release 2>&1 | tail -20`.
   - If it fails to compile, either fix immediately (if the fix is obvious
     and small) or `git checkout -- src/point_add.rs` and pick a different
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
     - `git checkout -- src/point_add.rs` to revert. The TSV row stays;
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

Cheap / local:
- **Eliminate redundant mod-reductions**: our `mod_add_qq` always does
  sub-p + cond-add-p. Many chained adds can defer reduction until the end.
- **Share the `cmod_add_qq` scratch register** across iterations of a
  multiplication instead of alloc/free per iteration.
- **Cuccaro → Draper / QFT adder**: different adder, different Toffoli cost.
- **Replace `mod_neg_inplace` + `mod_add_qq`** in `mod_sub_qq` with a
  direct Cuccaro subtraction (no extra mod_neg round trip).

Medium:
- **Windowed multiplication**: process y in windows of w bits, precomputing
  multiples of x. Cuts the number of conditional adds by a factor of w at
  the cost of a lookup table. See Roetteler et al. 2017 §4.
- **Montgomery-domain arithmetic**: do the whole EC-add in Montgomery form
  so `mod_mul` is cheaper than the schoolbook shift-and-add we use now.
- **Better inverse**: Bernstein–Yang "safegcd" is O(n²) like Kaliski but
  with smaller constants and simpler control flow. Or: skip inverses
  entirely by working in projective (Jacobian) coordinates — but that
  changes the algorithm significantly.

Structural:
- **Fold the two Kaliski calls** in the forward + uncompute halves by
  keeping the inverse in an ancilla register instead of recomputing.
- **Fuse step 5 (`Py -= λ·dx`) with step 3 (`lam += dy/dx`)**: since
  Py ends up zero, there's structure to exploit.
- **Classical conditioning**: many CCX/CCZ in our current code run
  unconditionally. If a CCX's control can be decided classically for a
  given input, wrap it in `push_condition` / `pop_condition` and it stops
  costing Toffolis in `sim.stats`. Look for branches where one side is
  always classical.

Structure research in sweeps: pick one axis (e.g., "replace Kaliski with
safegcd"), implement it, measure, revert or keep. Don't try two ideas at
once — you can't attribute the result.

## Rules of thumb

- If a run takes longer than 5 minutes, something is wrong — kill and revert.
- Cliffords are free compared to Toffolis (~100× cheaper in fault tolerance).
  Do not optimize Cliffords at the cost of Toffolis.
- X/Z gates are not counted at all. Abuse them.
- The simulator's `R` gate unconditionally zeros a qubit — you can use it
  to "free" dirty ancillas without proper uncompute. `Builder::free_qubit`
  already does this. But remember that what you compute must still be
  correct as a mathematical function; you cannot use R to paper over bugs.
- Correctness is non-negotiable. A 0-Toffoli circuit that fails correctness
  is worth nothing. Run `bash scripts/bench.sh` after every edit.

## Stop conditions

Keep iterating until one of:
- You hit the zenodo low-qubit target (2.7M Toffoli @ ≤1175 qubits).
- You get stuck: 10 consecutive experiments with no improvement.
  In that case, try a structurally different idea (switch category in the
  seed list). Do not pause for human input.
- The user interrupts.
