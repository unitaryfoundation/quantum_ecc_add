#!/usr/bin/env python3
"""Validate Qarton flattened streams with explicit quantum/classical status.

The earlier flattened-stream validator treated a measured position as classical
until an explicit reset. Qarton can reuse projected positions without a reset:
the memory manager turns a classical bit back into a quantum bit when a later
gate structurally promotes it. This script mirrors those status transitions.
"""

from __future__ import annotations

import random
import os
from collections import Counter
from dataclasses import dataclass
from typing import Any

from qarton.binary_operations.and_ccx import AndGate

AndGate.replace_by_ccx = True

from point_add.gcd import IPModMul
from qarton.circuit.basic_gates import (
    CCX_GATE,
    CCZ_GATE,
    CX_GATE,
    CZ_GATE,
    H_GATE,
    MEASURE,
    RESET,
    SWAP_GATE,
    X_GATE,
    Z_GATE,
    BasicGate,
)
from qarton.circuit.circuit import (
    BasicOperation,
    ConditionalOperation,
    ControlledOperation,
    IteratedOperation,
)


P_BITS = int(os.environ.get("QARTON_VALIDATE_BITS", "127"))
P = 2**P_BITS - 1


def positions(reg: Any) -> list[int]:
    s = reg.positions
    return s.to_list() if hasattr(s, "to_list") else list(s)


def interval_positions(s: Any) -> list[int]:
    return s.to_list() if hasattr(s, "to_list") else list(s)


def gate_tuple(g: Any) -> tuple[Any, list[int], tuple[int, ...]]:
    t = g.targets
    targets = t.to_list() if hasattr(t, "to_list") else list(t)
    controls = tuple(g.controls) if isinstance(g, ControlledOperation) else ()
    return g.op, targets, controls


def flatten_one(operation: Any) -> list[Any]:
    """Flatten one Qarton operation, preserving classical controls."""
    out: list[Any] = []
    if isinstance(operation, ConditionalOperation):
        raise RuntimeError(f"conditional operation unsupported in segment mode: {operation}")
    if isinstance(operation, ControlledOperation):
        controls = operation.controls
        iterations = 1
        true_op: Any = operation
    elif isinstance(operation, IteratedOperation):
        controls = ()
        iterations = operation.iterations
        true_op = BasicOperation(operation.op, operation.targets)
    else:
        controls = ()
        iterations = 1
        true_op = operation

    for _ in range(iterations):
        op, targets = operation.op, operation.targets
        if isinstance(op, BasicGate):
            if controls:
                out.append(ControlledOperation(true_op.op, targets, controls))
            else:
                out.append(true_op)
        else:
            for tmp in op.iterate_basic_gates():
                op1 = tmp.op
                true_targets = targets.get_at_positions(tmp.targets)
                true_controls = (
                    controls
                    if isinstance(tmp, BasicOperation)
                    else tuple(targets[_t] for _t in tmp.controls) + controls
                )
                if true_controls:
                    out.append(ControlledOperation(op1, true_targets, true_controls))
                else:
                    out.append(BasicOperation(op1, true_targets))
    return out


@dataclass
class StatusSim:
    qval: dict[int, int]
    cval: dict[int, int]
    quantum: set[int]
    phase: int = 0
    pending_h: set[int] | None = None

    def __post_init__(self) -> None:
        if self.pending_h is None:
            self.pending_h = set()

    def val(self, p: int) -> int:
        if p in self.quantum:
            return self.qval.get(p, 0)
        return self.cval.get(p, 0)

    def set_value(self, p: int, v: int) -> None:
        if p in self.quantum:
            self.qval[p] = v & 1
        else:
            self.cval[p] = v & 1

    def promote(self, p: int) -> None:
        if p not in self.quantum:
            self.qval[p] = self.cval.pop(p, 0)
            self.quantum.add(p)

    def demote_measured(self, p: int, r: int) -> None:
        self.quantum.discard(p)
        self.qval[p] = 0
        self.cval[p] = r & 1

    def reset(self, p: int) -> None:
        self.quantum.discard(p)
        self.qval[p] = 0
        self.cval.pop(p, None)

    def controls_on(self, controls: tuple[int, ...]) -> bool:
        return all(self.val(c) for c in controls)

    def apply(self, op: Any, t: list[int], controls: tuple[int, ...], rng: random.Random) -> None:
        active = self.controls_on(controls)

        if op == X_GATE:
            if active:
                self.set_value(t[0], self.val(t[0]) ^ 1)
        elif op == Z_GATE:
            if active:
                self.phase ^= self.val(t[0])
        elif op == CX_GATE:
            c, q = t
            if c in self.quantum:
                self.promote(q)
            if active:
                self.set_value(q, self.val(q) ^ self.val(c))
        elif op == CZ_GATE:
            if t[0] in self.quantum and t[1] not in self.quantum:
                self.promote(t[1])
            if active:
                self.phase ^= self.val(t[0]) & self.val(t[1])
        elif op == CCX_GATE:
            c1, c2, q = t
            if c1 in self.quantum or c2 in self.quantum:
                self.promote(q)
            if active:
                self.set_value(q, self.val(q) ^ (self.val(c1) & self.val(c2)))
        elif op == CCZ_GATE:
            if active:
                self.phase ^= self.val(t[0]) & self.val(t[1]) & self.val(t[2])
        elif op == SWAP_GATE:
            a, b = t
            av, bv = self.val(a), self.val(b)
            aq, bq = a in self.quantum, b in self.quantum
            if aq:
                self.quantum.remove(a)
            if bq:
                self.quantum.remove(b)
            self.cval.pop(a, None)
            self.cval.pop(b, None)
            self.qval[a], self.qval[b] = bv, av
            if bq:
                self.quantum.add(a)
            else:
                self.cval[a] = bv
            if aq:
                self.quantum.add(b)
            else:
                self.cval[b] = av
        elif op == H_GATE:
            # In this codebase H is used as the first half of H+measure.
            # Promotion clears any stale measured-control mapping.
            self.promote(t[0])
            assert self.pending_h is not None
            self.pending_h.add(t[0])
        elif op == MEASURE:
            q = t[0]
            assert self.pending_h is not None
            if q not in self.pending_h:
                raise RuntimeError(f"measure without pending H at {q}")
            self.pending_h.remove(q)
            if active:
                r = rng.choices([0, 1], weights=[0.5, 0.5], k=1)[0]
                self.phase ^= self.qval.get(q, 0) & r
                self.demote_measured(q, r)
        elif op == RESET:
            for q in t:
                self.reset(q)
        else:
            raise RuntimeError(f"unmapped op {op} targets={t} controls={controls}")


def run_once(qc: Any, gates: list[Any], x: int, y: int, seed: int) -> tuple[bool, Counter[str], dict[str, Any]]:
    xpos = positions(qc.input_signature[0])
    ypos = positions(qc.input_signature[1])
    qval: dict[int, int] = {}
    for i, p in enumerate(xpos):
        if (x >> i) & 1:
            qval[p] = 1
    for i, p in enumerate(ypos):
        if (y >> i) & 1:
            qval[p] = 1
    sim = StatusSim(qval=qval, cval={}, quantum=set(xpos) | set(ypos))
    rng = random.Random(seed)
    counts: Counter[str] = Counter()
    for g in gates:
        op, targets, controls = gate_tuple(g)
        counts[getattr(op, "name", None) or str(op)] += 1
        sim.apply(op, targets, controls, rng)

    xout = sum(sim.val(p) << i for i, p in enumerate(xpos))
    yout = sum(sim.val(p) << i for i, p in enumerate(ypos))
    outputs = set(xpos) | set(ypos)
    dirty = [p for p in set(sim.qval) | set(sim.cval) if sim.val(p) and p not in outputs]
    expect = (x * y) % P
    detail = {
        "x_ok": xout == x,
        "y_ok": yout == expect,
        "phase": sim.phase,
        "dirty": len(dirty),
        "dirty_sample": dirty[:20],
        "x": x,
        "y": y,
        "xout": xout,
        "yout": yout,
        "expect": expect,
        "quantum_live": len(sim.quantum),
        "classical_live": len([p for p, v in sim.cval.items() if v]),
    }
    return detail["x_ok"] and detail["y_ok"] and sim.phase == 0 and not dirty, counts, detail


def read_bits(sim: StatusSim, pos: list[int]) -> int:
    return sum(sim.val(p) << i for i, p in enumerate(pos))


def bv_bits(v: Any) -> list[int]:
    return [int(v[i]) for i in range(len(v))]


def segment_debug(qc: Any, x: int, y: int) -> None:
    from point_add.gcd_functions import apply_bitvector, to_bitvector

    xpos = positions(qc.input_signature[0])
    ypos = positions(qc.input_signature[1])
    qval: dict[int, int] = {}
    for i, p in enumerate(xpos):
        if (x >> i) & 1:
            qval[p] = 1
    for i, p in enumerate(ypos):
        if (y >> i) & 1:
            qval[p] = 1
    sim = StatusSim(qval=qval, cval={}, quantum=set(xpos) | set(ypos))
    rng = random.Random(0)
    dialog = to_bitvector(P, x)
    expect = (x * y) % P

    for idx, top in enumerate(qc.get_operations()):
        gates = flatten_one(top)
        for g in gates:
            op, targets, controls = gate_tuple(g)
            sim.apply(op, targets, controls, rng)

        print(f"segment {idx}: gates={len(gates)} op={top.op}", flush=True)
        if idx == 0:
            d_abs = interval_positions(
                top.targets.get_at_positions(top.op.output_signature[0].positions)
            )
            got = [sim.val(p) for p in d_abs]
            want = bv_bits(dialog)
            print(f"  dialog_ok={got == want} ones={sum(got)}/{sum(want)}", flush=True)
        elif idx == 1:
            sig = top.op.output_signature
            y_abs = interval_positions(top.targets.get_at_positions(sig[1].positions))
            tmp_abs = interval_positions(top.targets.get_at_positions(sig[2].positions))
            print(
                f"  apply_y_zero={read_bits(sim, y_abs) == 0} "
                f"tmp_ok={read_bits(sim, tmp_abs) == expect}",
                flush=True,
            )
        elif idx == 2:
            print(f"  after_swap_y_ok={read_bits(sim, ypos) == expect}", flush=True)
        elif idx == 5:
            outputs = set(xpos) | set(ypos)
            dirty = [p for p in set(sim.qval) | set(sim.cval) if sim.val(p) and p not in outputs]
            print(
                f"  final x_ok={read_bits(sim, xpos) == x} "
                f"y_ok={read_bits(sim, ypos) == expect} phase={sim.phase} dirty={len(dirty)}",
                flush=True,
            )


def main() -> None:
    qc = IPModMul(P, gate_efficient=False, special_prime=False)
    print(
        f"flatten-status {P_BITS}-bit IPModMul bits={qc.nbr_bits()} "
        f"qubits={qc.nbr_qubits()}",
        flush=True,
    )
    gates = list(qc.iterate_basic_gates())
    print(f"flattened gates={len(gates)}", flush=True)
    if os.environ.get("QARTON_SEGMENT_DEBUG") == "1":
        x = int(os.environ.get("QARTON_X", "1231849834227650594"))
        y = int(os.environ.get("QARTON_Y", "1032789163651940088"))
        segment_debug(qc, x % P, y % P)
        return
    total = 1
    ok = 0
    aggregate: Counter[str] = Counter()
    for trial in range(total):
        x = random.randrange(max(1, 1 << (P_BITS - 2)), P)
        y = random.randrange(0, P)
        good, counts, detail = run_once(qc, gates, x, y, seed=trial)
        aggregate.update(counts)
        ok += int(good)
        if not good:
            print(f"trial {trial}: {detail}", flush=True)
    print("ops", dict(aggregate), flush=True)
    print(f"correct: {ok}/{total}", flush=True)


if __name__ == "__main__":
    main()
