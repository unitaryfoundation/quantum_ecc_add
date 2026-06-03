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
