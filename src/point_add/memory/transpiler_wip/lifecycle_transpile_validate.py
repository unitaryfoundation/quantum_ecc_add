#!/usr/bin/env python3
"""Lifecycle-aware Qarton lowering experiment.

This mirrors Qarton's QuantumSimulator decomposition instead of using
iterate_basic_gates, because iterate_basic_gates loses the lifetime of measured
positions that are later reused as qubits.
"""

from __future__ import annotations

import random
from collections import Counter
from collections.abc import Iterator
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
from qarton.circuit.bit_vector import BitVector
from qarton.circuit.circuit import (
    BasicOperation,
    Circuit,
    ConditionalOperation,
    ControlledOperation,
    IteratedOperation,
    Operation,
)
from qarton.circuit.util import IntIntervalTuple


P = 2**127 - 1


def positions(reg: Any) -> list[int]:
    s = reg.positions
    return s.to_list() if hasattr(s, "to_list") else list(s.to_set())


def gate_name(op: Any) -> str:
    return getattr(op, "name", None) or str(op)


@dataclass
class Lowered:
    ops: list[tuple]
    qubit_controls: set[int]
    measured: set[int]
    control_values: dict[int, int]
    q: dict[int, int]
    phase: int = 0

    def bit(self, p: int) -> int:
        if p in self.control_values:
            return self.control_values[p]
        return self.q.get(p, 0)

    def controls_on(self, controls: tuple[int, ...]) -> bool:
        return all(self.bit(c) for c in controls)

    def emit_x(self, t: int, controls: tuple[int, ...] = ()) -> None:
        if not self.controls_on(controls):
            return
        if len(controls) == 0:
            self.ops.append(("X", t))
        elif len(controls) == 1:
            self.ops.append(("CX", controls[0], t))
        elif len(controls) == 2:
            self.ops.append(("CCX", controls[0], controls[1], t))
        else:
            raise RuntimeError(f"unsupported x controls: {controls}")
        self.q[t] = self.q.get(t, 0) ^ 1

    def emit_phase(self, targets: tuple[int, ...], controls: tuple[int, ...] = ()) -> None:
        all_controls = targets + controls
        if not self.controls_on(all_controls):
            return
        if len(targets) == 1 and len(controls) == 0:
            self.ops.append(("Z", targets[0]))
        elif len(targets) == 2 and len(controls) == 0:
            self.ops.append(("CZ", targets[0], targets[1]))
        elif len(targets) == 3 and len(controls) == 0:
            self.ops.append(("CCZ", targets[0], targets[1], targets[2]))
        else:
            # For validation, keep symbolic op shape. Harness lowering can expand
            # this later with classically-conditioned phase ops.
            self.ops.append(("PHASE", targets, controls))
        self.phase ^= 1

    def emit_hmr(self, q: int, rng: random.Random) -> None:
        r = rng.randint(0, 1)
        self.ops.append(("HMR", q, r))
        self.phase ^= self.q.get(q, 0) & r
        self.q[q] = 0
        self.control_values[q] = r
        self.measured.add(q)


class Lowering:
    def __init__(self, qc: Circuit, rng: random.Random, initial: dict[int, int]):
        self.qc = qc
        self.rng = rng
        self.lowered = Lowered(
            ops=[],
            qubit_controls=set(),
            measured=set(),
            control_values={},
            q=dict(initial),
        )

    def decompose(self, operation: Operation) -> Iterator[BasicOperation]:
        if isinstance(operation, ConditionalOperation):
            data: list[Any] = []
            for t, d in zip(operation.controls, operation.datatypes):
                sub = BitVector(self.lowered.bit(i) for i in tuple(t))
                data.append(d.from_register(sub))
            sub_circuit = operation.true_op(*data)
            yield from self.decompose(BasicOperation(sub_circuit, operation.targets))
        elif isinstance(operation, IteratedOperation):
            for _ in range(operation.iterations):
                yield from self.decompose(BasicOperation(operation.op, operation.targets))
        elif isinstance(operation, ControlledOperation):
            if self.lowered.controls_on(operation.controls):
                yield from self.decompose(BasicOperation(operation.op, operation.targets))
        else:
            op, targets = operation
            if isinstance(op, BasicGate):
                yield operation
            else:
                for new_op in op.get_operations():
                    yield from self.decompose(new_op.map(targets))

    def run(self) -> Lowered:
        for op in self.qc.get_operations():
            for bop in self.decompose(op):
                self.apply_basic(bop)
        return self.lowered

    def apply_basic(self, bop: BasicOperation) -> None:
        op, targets_raw = bop
        targets = tuple(targets_raw)
        if op == X_GATE:
            self.lowered.emit_x(targets[0])
        elif op == CX_GATE:
            self.lowered.emit_x(targets[1], (targets[0],))
        elif op == CCX_GATE:
            self.lowered.emit_x(targets[2], (targets[0], targets[1]))
        elif op == SWAP_GATE:
            a, b = targets
            self.lowered.ops.append(("SWAP", a, b))
            self.lowered.q[a], self.lowered.q[b] = self.lowered.q.get(b, 0), self.lowered.q.get(a, 0)
        elif op == H_GATE:
            self.lowered.ops.append(("H", targets[0]))
        elif op == MEASURE:
            self.lowered.emit_hmr(targets[0], self.rng)
        elif op == RESET:
            for t in targets:
                self.lowered.ops.append(("RESET", t))
                self.lowered.control_values.pop(t, None)
                self.lowered.measured.discard(t)
                self.lowered.q[t] = 0
        elif op == Z_GATE:
            self.lowered.emit_phase((targets[0],))
        elif op == CZ_GATE:
            self.lowered.emit_phase((targets[0], targets[1]))
        elif op == CCZ_GATE:
            self.lowered.emit_phase((targets[0], targets[1], targets[2]))
        else:
            raise RuntimeError(f"unmapped basic op {op} targets={targets}")


def main() -> None:
    qc = IPModMul(P, gate_efficient=False, special_prime=False)
    xpos = positions(qc.input_signature[0])
    ypos = positions(qc.input_signature[1])
    print(f"lifecycle-lowering 127-bit IPModMul bits={qc.nbr_bits()} qubits={qc.nbr_qubits()}")

    total = 20
    ok = 0
    kinds = Counter()
    for trial in range(total):
        x = random.randrange(1 << 126, P)
        y = random.randrange(0, P)
        initial = {}
        for i, p in enumerate(xpos):
            if (x >> i) & 1:
                initial[p] = 1
        for i, p in enumerate(ypos):
            if (y >> i) & 1:
                initial[p] = 1
        lowered = Lowering(qc, random.Random(trial), initial).run()
        kinds.update(o[0] for o in lowered.ops)
        xout = sum(lowered.q.get(p, 0) << i for i, p in enumerate(xpos))
        yout = sum(lowered.q.get(p, 0) << i for i, p in enumerate(ypos))
        outset = set(xpos) | set(ypos)
        dirty = [p for p, v in lowered.q.items() if v and p not in outset]
        good = xout == x and yout == (x * y) % P and lowered.phase == 0 and not dirty
        if good:
            ok += 1
        elif trial < 3:
            print(
                f"trial {trial}: x={xout == x} y={yout == (x*y)%P} "
                f"phase={lowered.phase} dirty={len(dirty)}"
            )
    print("ops", dict(kinds))
    print(f"correct: {ok}/{total}")


if __name__ == "__main__":
    main()
