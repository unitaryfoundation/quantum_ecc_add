# Phase-0 validation. Algorithm + compressor truth table derived from the AGPL
# reference: A. Schrottenloher, gitlab.inria.fr/capsule/qarton-projects/ec-point-addition
# (arXiv:2606.02235). This is an analysis/validation reimplementation.
#!/usr/bin/env python3
"""Phase-0 go/no-go: standalone classical model of the Schrottenloher dialog/Bezout
split-EEA modular inversion (no Qarton dependency). Validates on real secp256k1 inputs
and measures the failure rate at the paper's ITERATIONS_VAR=2.4."""
import random
from math import ceil, log2, sqrt

P = (1 << 256) - (1 << 32) - 977  # secp256k1 field prime (odd)

# compressor truth tables (from compressor.py), keys/vals as tuples of bits
_F = {
 (0,0,0,0,0,0):(0,0,1,0,1),(0,0,0,0,1,0):(0,0,1,0,0),(0,0,0,0,1,1):(0,0,1,1,1),
 (0,0,1,0,0,0):(0,0,0,0,1),(0,0,1,0,1,0):(0,0,0,0,0),(0,0,1,0,1,1):(0,0,0,1,1),
 (0,0,1,1,0,0):(1,1,1,1,1),(0,0,1,1,1,0):(0,0,1,1,0),(0,0,1,1,1,1):(1,1,1,0,1),
 (1,0,0,0,0,0):(1,0,0,0,1),(1,0,0,0,1,0):(1,0,0,0,0),(1,0,0,0,1,1):(1,0,0,1,1),
 (1,0,1,0,0,0):(1,0,1,0,1),(1,0,1,0,1,0):(1,0,1,0,0),(1,0,1,0,1,1):(1,0,1,1,1),
 (1,0,1,1,0,0):(1,1,0,1,1),(1,0,1,1,1,0):(1,0,0,1,0),(1,0,1,1,1,1):(1,1,0,0,1),
 (1,1,0,0,0,0):(0,1,1,0,0),(1,1,0,0,1,0):(0,1,1,0,1),(1,1,0,0,1,1):(0,1,1,1,0),
 (1,1,1,0,0,0):(0,1,0,0,0),(1,1,1,0,1,0):(0,1,0,0,1),(1,1,1,0,1,1):(0,1,0,1,0),
 (1,1,1,1,0,0):(1,1,1,1,0),(1,1,1,1,1,0):(0,1,1,1,1),(1,1,1,1,1,1):(1,1,1,0,0),
}
_R = {v:k for k,v in _F.items()}
def compress(t):   return _F[tuple(t)]
def uncompress(t): return list(_R[tuple(t)])

ITERATIONS_VAR = 2.4
U_PAD_VAR = 2.3

def to_bitvector(u, v):
    assert u % 2 == 1
    n = max(u.bit_length(), v.bit_length())
    it = ceil((1.413*n + ITERATIONS_VAR*sqrt(n)) / 3) * 3
    upad = ceil(U_PAD_VAR*sqrt(n))
    ng = (it//3)*5
    g = []
    for i in range(it//3):
        g += list(compress([0,0,0,0,0,0]))
    for i in range(it):
        b0 = v & 1
        b1 = 1 if u > v else 0
        small = g[(i//3)*5:(i//3+1)*5]
        d = uncompress(small)
        if d[2*(i%3)] != 0 or d[2*(i%3)+1] != 0: return None  # assert
        d[2*(i%3)] = b0
        d[2*(i%3)+1] = b0 & b1
        g[(i//3)*5:(i//3+1)*5] = list(compress(d))
        if u.bit_length() >= max(n - i*0.5*(3-log2(3)) + upad, 0): return None
        if v.bit_length() >= max(n - i*0.5*(3-log2(3)) + upad, 0): return None
        if b0 & b1: u, v = v, u
        if b0:     v -= u
        v >>= 1
    if v != 0 or u != 1: return None
    return g

def apply_bitvector(x, y, d, p):
    u, v = x, y
    nsteps = len(d)//5*3
    inv2 = pow(2, -1, p)
    for i in reversed(range(nsteps)):
        small = d[(i//3)*5:(i//3+1)*5]
        dec = uncompress(small)
        b0 = dec[2*(i%3)]; b0b1 = dec[2*(i%3)+1]
        v = (v*2) % p
        if b0: v = (v+u) % p
        if b0b1: u, v = v, u
    return u, v

def main():
    random.seed(12345)
    N = 10000
    fails = 0; wrong = 0; ok = 0
    for _ in range(N):
        dx = random.randrange(1, P)
        dy = random.randrange(0, P)
        d = to_bitvector(P, dx)
        if d is None:
            fails += 1
            continue
        _, lam = apply_bitvector(dy, 0, d, P)
        expected = (dy * pow(dx, -1, P)) % P
        if lam == expected:
            ok += 1
        else:
            wrong += 1
    n = 256
    it = ceil((1.413*n + ITERATIONS_VAR*sqrt(n)) / 3) * 3
    print(f"secp256k1 p, n={n}, iterations={it}, dialog_bits={(it//3)*5} (~{(it//3)*5/n:.2f}n)")
    print(f"inputs tested: {N}")
    print(f"  correct (lam == dy*dx^-1): {ok}")
    print(f"  GCD assertion failures   : {fails}")
    print(f"  WRONG results            : {wrong}")
    print(f"failure rate = {fails+wrong}/{N} = {(fails+wrong)/N:.5f}  (target ~2^-13.3 = {2**-13.3:.2e})")
    print("VERDICT:", "GO -- algorithm validated, math correct, failures only from iteration budget" if wrong==0 else "NO-GO -- wrong results, model misunderstood")

if __name__=="__main__": main()
