# Transpile route: Qarton reference -> harness Op stream (WORK IN PROGRESS)

Goal: instead of hand-porting, lower the open-source Schrottenloher reference
circuit (Qarton, AGPL: gitlab.inria.fr/capsule/qarton-projects/ec-point-addition)
into the harness Op format, then polish to beat 2.48e9.

## Confirmed (this session)
- Reference builds: `build_circuit.py --gate_efficient` -> 1445 q x 1,865,521 Toffoli
  = 2.70e9 (leaderboard 2.48e9 = this + polish).
- `qarton.binary_operations.and_ccx.AndGate.replace_by_ccx = True` makes AND-compute
  a single CCX; AndGateUncompute stays as h;measure;cz(ctrl=measured);reset (=MBUC).
- With that flag, `circuit.iterate_basic_gates()` flattens the WHOLE circuit to
  {x,cx,ccx,h,measure,cz,reset,swap,z} -- ALL map to harness Ops:
    ccx->CCX, cx->CX, x->X, z->Z, swap->Swap,
    h;measure;reset  ->  Hmr(q, cbit)
    cz(targets, ctrl=measured q) -> cz_if(targets, cbit)
- MBUC pattern confirmed in the stream: `h(a); measure(a); cz([i,j], ctrl=a); reset(a)`.
- transpile_validate.py: transpiles 127-bit IPModMul -> 2,126,647 harness ops
  (CCX 465572, CX 945187, X 354066, Z 138276, HMR 116184, CZIF 47046, SWAP 60316).

## THE BUG TO FIX (blocks validation)
Qarton uses a UNIFIED position space for qubits AND classical bits
(nbr_bits=858 > nbr_qubits=606 for the 127-bit IPModMul; positions 606..857 are
classical bits). The basis-state Op-sim currently treats every position as a
qubit, so even `x` is not preserved. Must:
  1. Separate qubit positions from classical-bit positions (use input/output
     signatures + track which positions are classical).
  2. Model `measure` writing to its classical-bit position; `reset` clearing the
     qubit; correctly map the measured-bit control of cz to the harness cbit.
  3. Some measures have NO following cz/reset (h=measure=116184 but reset=47874,
     cz=47046) -- these are a different measurement use (likely clean ancilla
     release) and must map to R (measure-reset, discard) not Hmr+cz_if.

## Remaining after the bug (the multi-day tail)
- Validate transpiler on 127-bit IPModMul (xout==x, yout==x*y, phase==0, clean).
- Adapt the BENCHMARK INTERFACE: harness adds a CLASSICAL-register point (ox,oy)
  to a quantum point (tx,ty); Qarton's circuit adds a quantum/windowed point.
  Either modify point_add to take a classical operand, or load ox,oy into quantum
  ancillas (classically-conditioned X), run the q-q point add, unload.
- Map to the 4-register contract (tx,ty quantum; ox,oy classical bits).
- Scale to secp256k1 (256-bit), emit ~12M ops, embed in build() (Rust) for a
  valid submission, pass all 9024 shots.
- AGPL: output is AGPL-derived; submission must comply.

## Honest status
NOT beaten. 2.48e9 unbeaten; current real submission on main = 9.69e9.
The transpile route is de-risked to "transpiler runs on the real circuit;
debugging the measurement/classical-bit model"; the rest is sustained engineering.

## UPDATE: real obstacle in the iterate_basic_gates lowering
Two correction mechanisms, both controlled by measured CLASSICAL bits:
- cz: 47046, each (2 quantum targets, 1 measured control)  [AND-uncompute MBUC]
- z : 138276, each (1 quantum target, 1 measured control)  [classical-controlled Z phase]
The control positions are classical-bit positions (>= nbr_qubits, e.g. 645 with
nbr_qubits=606; nbr_bits=858). BUT iterate_basic_gates' `measure` op only carries
the measured QUBIT, not its output classical-bit position -- so the measure->cbit
dataflow (which correction reads which measurement) is LOST in this lowering.
=> A faithful transpile cannot use iterate_basic_gates alone; it needs a lowering
that preserves classical dataflow (e.g. via quantum_simulator.decompose_operation,
which threads measurement outputs, or a custom traversal of the op tree that keeps
classical bit ids). That is the next real step, and it is non-trivial.

Bug also found: controlled gates (all `z`, all `cz`) must NOT be emitted as bare
Z/CZ -- they carry measured controls. Earlier transpile_validate.py dropped z's
control (138k gates) -> garbage. Fix requires the cbit dataflow above.

## RESOLVED: Qarton measure -> cbit mapping
The classical bit IS the measured qubit's own position. Verified on 127-bit IPModMul:
- every `h(q)` is immediately followed by `measure(q)` (116184/116184) -> X-basis measure.
- ALL correction/feedforward controls are measure-target positions (subset check = True):
    cz corrections: 47046  (2 qubit targets, ctrl = measured pos)
    z  corrections: 138276 (1 qubit target,  ctrl = measured pos)
    cx feedforward: 99611  (control = measured pos)   <-- classical feedforward!
    ccx feedforward: 16146 (a control = measured pos)
- positions >= nbr_qubits (e.g. 645 with nbr_qubits=606) are just high-INDEX ancilla
  qubits; 606 is peak-simultaneous-live, not max index.
Mapping to harness: h(q);measure(q);...;reset(q) -> Hmr(q,c); every later gate whose
control is position q (until reset) is conditioned on classical bit c
(cz/z -> cz_if / conditioned-Z; cx/ccx -> classically-conditioned X/CX).
sim2.py implements this (classical[pos]=outcome; cval() reads it for controls).

## NEXT separate bug (not the mapping): basis-state sim still off
With the mapping applied, sim2.py still fails (xout!=x, ~73 dirty ancillas). The
mapping is correct; the remaining issue is matching the harness sim's measurement+
feedforward semantics EXACTLY for a single shot (collapsed-value vs reset-to-0 of a
measured position; ordering of phase kickback vs feedforward). sim.rs Hmr sets the
qubit to 0, classical bit = rng, phase ^= qubit&rng. Need to diff sim2 against sim.rs
gate-by-gate (or against Qarton's own simulator) to find the divergent op. This is the
next scoped step.

## KEY FINDING: iterate_basic_gates is architecturally insufficient for transpile
- Qarton's own simulator (qc.simulate, amplitude-accurate) CONFIRMS IPModMul(127-bit)
  computes (x, x*y mod p) exactly -> the reference circuit is correct; the oracle works.
  (Note: degenerate tiny inputs like x=5 trip the approximate circuit's assertions;
   must test with full-size random inputs.)
- WHY my basis-state transpile sim fails: only 47,874 resets for 116,184 measures, and
  424,564 "writes to a measured position" -> positions are measured, FREED, and REUSED
  as fresh qubits WITHOUT an explicit reset. iterate_basic_gates carries no qubit
  alloc/free/measure-destination lifecycle, so "cx controlled by a measured position"
  (99k) are mostly FALSE POSITIVES (position long-since reused as a fresh qubit).
- => A faithful transpile cannot use iterate_basic_gates. It must use a LIFECYCLE-AWARE
  lowering that threads measurement outcomes, resets, and qubit reallocations -- i.e.
  go through qarton.circuit.quantum_simulator's operation processing (decompose_operation
  / apply_*), which already tracks this, or a custom traversal that records each measure's
  destination and each qubit free/realloc.

## Next concrete step
Build the transpiler on top of quantum_simulator's op processing (it handles measure/
reset/realloc correctly) rather than iterate_basic_gates, emitting harness Ops as it
goes. Validate each component vs qc.simulate (the confirmed oracle).

## UPDATE: flattened IPModMul is salvageable with ordered positions + status tracking
The earlier "iterate_basic_gates is architecturally insufficient" conclusion was too
broad for `IPModMul`. The actual bugs in the flattened validator were:
- It converted `IntIntervalTuple` through `to_set()`, destroying the intentional
  interleaved register order used by Qarton remaps/dialog registers.
- It treated measured positions as classical until reset only; Qarton can promote a
  projected position back to quantum when a later gate structurally uses it as a
  quantum target.

`flatten_status_validate.py` now mirrors the Qarton quantum/classical status rules
and preserves tuple order. Validated:
- 61-bit `IPModMul(2^61-1)`: every top-level segment passes
  (`ToBitVector`, `ApplyBitVector`, swap, inverse, `Remap`), final phase=0/dirty=0.
- 127-bit `IPModMul(2^127-1)`: flattened 2,290,705 gates
  (`465,572` CCX) and passed 1/1 random full-size trial with phase=0/dirty=0.

This reopens a concrete route: transpile Qarton's reusable `IPModMul` and
`IPModMul.inverse()` into harness ops, then wire the benchmark's arbitrary
classical-offset affine formula around those blocks instead of trying to use
Qarton's fixed-window outer point-add circuit directly.
