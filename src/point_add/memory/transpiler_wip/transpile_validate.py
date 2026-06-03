#!/usr/bin/env python3
"""Transpile a Qarton circuit's flattened gate stream into the harness Op model
(CCX/CX/X/Z/SWAP/CZ + Hmr/cz_if) and validate it end-to-end on a basis-state
Op-simulator: prove it computes x*y mod p with clean ancillas and zero phase
(random measurement outcomes test the MBUC phase corrections)."""
import random, sys

from qarton.binary_operations.and_ccx import AndGate
AndGate.replace_by_ccx = True
from point_add.gcd import IPModMul
from qarton.circuit.circuit import ControlledOperation

P = 2**127 - 1
N = 127

def positions(reg):
    # ordered list of qubit positions for a register (LSB first)
    s = reg.positions
    if hasattr(s, "to_list"): lst = s.to_list()
    else: lst = list(s.to_set())
    return lst

def gate_tuple(g):
    nm = getattr(g.op, "name", None) or str(g.op)
    t = g.targets
    tl = t.to_list() if hasattr(t, "to_list") else list(t.to_set())
    c = list(g.controls) if isinstance(g, ControlledOperation) else []
    return nm, tl, c

def transpile(qc):
    """Return (ops, nbits_positions). ops are tuples in the harness model."""
    ops = []
    # classical-bit allocator (separate id space)
    next_cbit = [0]
    # which positions are currently 'measured' -> classical bit id
    measured = {}  # pos -> cbit id
    it = qc.iterate_basic_gates()
    pending_h = {}  # pos -> True if an h was just applied (expect measure next)
    for g in it:
        nm, t, c = gate_tuple(g)
        if nm == "x":
            ops.append(("X", t[0]))
        elif nm == "z":
            ops.append(("Z", t[0]))
        elif nm == "cx":
            ops.append(("CX", t[0], t[1]))
        elif nm == "ccx":
            ops.append(("CCX", t[0], t[1], t[2]))
        elif nm == "swap":
            ops.append(("SWAP", t[0], t[1]))
        elif nm == "h":
            pending_h[t[0]] = True
        elif nm == "measure":
            q = t[0]
            assert pending_h.pop(q, False), "measure without preceding h at %d" % q
            cb = next_cbit[0]; next_cbit[0] += 1
            measured[q] = cb
            ops.append(("HMR", q, cb))   # H+measure->cb+reset
        elif nm == "reset":
            measured.pop(t[0], None)     # Hmr already reset; end classical window
        elif nm == "cz":
            # CZ on targets t, possibly controlled by a measured bit (the MBUC fix)
            if c and c[0] in measured:
                ops.append(("CZIF", t[0], t[1], measured[c[0]]))
            elif c:
                # cz controlled by a quantum qubit -> CCZ
                ops.append(("CCZ", t[0], t[1], c[0]))
            else:
                ops.append(("CZ", t[0], t[1]))
        else:
            raise RuntimeError("unmapped gate: " + nm)
    return ops, next_cbit[0]

def simulate(ops, init_qubits, nqpos, seed):
    rng = random.Random(seed)
    q = dict(init_qubits)  # pos -> 0/1 (default 0)
    def gq(p): return q.get(p, 0)
    cb = {}
    phase = 0
    for op in ops:
        k = op[0]
        if k == "X": q[op[1]] = gq(op[1]) ^ 1
        elif k == "Z": phase ^= gq(op[1])
        elif k == "CX": q[op[2]] = gq(op[2]) ^ gq(op[1])
        elif k == "CCX": q[op[3]] = gq(op[3]) ^ (gq(op[1]) & gq(op[2]))
        elif k == "SWAP":
            a,b=op[1],op[2]; q[a],q[b]=gq(b),gq(a)
        elif k == "CZ": phase ^= gq(op[1]) & gq(op[2])
        elif k == "CCZ": phase ^= gq(op[1]) & gq(op[2]) & gq(op[3])
        elif k == "HMR":
            qb, c = op[1], op[2]
            r = rng.randint(0,1)
            cb[c] = r
            phase ^= gq(qb) & r
            q[qb] = 0
        elif k == "CZIF":
            phase ^= gq(op[1]) & gq(op[2]) & cb.get(op[3],0)
    return q, phase

def main():
    qc = IPModMul(P, gate_efficient=False, special_prime=False)
    xpos = positions(qc.input_signature[0])
    ypos = positions(qc.input_signature[1])
    print("transpiling 127-bit IPModMul... qubits=%d" % qc.nbr_qubits())
    ops, ncb = transpile(qc)
    from collections import Counter
    kinds = Counter(o[0] for o in ops)
    print("transpiled ops:", len(ops), dict(kinds), "classical bits:", ncb)
    ok = 0; bad = 0; phasebad = 0; ancbad = 0
    for trial in range(40):
        x = random.randrange(1, P); y = random.randrange(0, P)
        init = {}
        for i,p in enumerate(xpos):
            if (x>>i)&1: init[p]=1
        for i,p in enumerate(ypos):
            if (y>>i)&1: init[p]=1
        q, phase = simulate(ops, init, qc.nbr_qubits(), seed=trial)
        xout = sum(((q.get(p,0))<<i) for i,p in enumerate(xpos))
        yout = sum(((q.get(p,0))<<i) for i,p in enumerate(ypos))
        # all non-input/output qubits must be 0
        outset = set(xpos)|set(ypos)
        anc_clean = all(v==0 for p,v in q.items() if p not in outset)
        expect = (x*y)%P
        if phase!=0: phasebad+=1
        if not anc_clean: ancbad+=1
        if yout==expect and xout==x and phase==0 and anc_clean: ok+=1
        else: bad+=1
    print("results over 40 random (x,y) with random measurements:")
    print("  correct (yout==x*y, xout==x, phase==0, ancillas clean):", ok)
    print("  bad:", bad, " phase!=0:", phasebad, " ancilla-not-clean:", ancbad)
    print("VERDICT:", "TRANSPILER VALIDATED on real 127-bit IPModMul" if ok==40 else "needs fixing")

main()
