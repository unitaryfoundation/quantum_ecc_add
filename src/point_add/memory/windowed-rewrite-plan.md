---
name: windowed-rewrite-plan
description: Port plan + Phase-0 validation for the Schrottenloher split-EEA low-qubit point-add (the real SOTA architecture)
metadata:
  type: project
---

# Rewrite to the SOTA low-qubit architecture (Schrottenloher 2026)

**Why:** the ecdsa.fail leaderboard best is ~2.48e9 (1434 qubits x 1.73M Toffoli);
our Kaliski circuit (committed best 9.69e9, 2711 q) is structurally ~2x too wide
and CANNOT reach it by tweaks. The leaders implement the published construction:

- Paper: Andre Schrottenloher, "Optimized Point Addition Circuits for ECDL",
  arXiv:2606.02235 (Jun 2026). Design principles fully described.
- Open-source reference (Python/Qarton, AGPL): gitlab.inria.fr/capsule/qarton-projects/ec-point-addition
  (cloned to /tmp/qarton-ref during analysis). Qarton lib: gitlab.inria.fr/capsule/qarton
- Target numbers (paper Table 1): space-opt secp256k1 = 1192 q / 2.37M T;
  gate-opt = 1446 q / 1.85M T. Leaderboard = gate-opt + polish.

## The architecture (dialog / Bezout split EEA, "Khattar et al." method)

Modular inverse WITHOUT storing full forward+backward GCD state:
1. **Forward Euclidean (`to_bitvector`)**: binary GCD on (u,v)=(p, x), ~402 iters
   (= ceil((1.413n + 2.4*sqrt(n))/3)*3). Each iter: b0=v&1, b1=u>v; record (b0, b0&b1);
   if b0&b1 swap(u,v); if b0 v-=u; v>>=1. Ends cleanly at (u,v)=(1,0) -- no backward
   pass, no stored (r,s). Emits a **dialog** = compressed op bitstring (~2.12-2.62n bits).
2. **Compressor**: packs 3 iters' (b0,b0&b1) pairs (each in {00,10,11}) -> 5 bits.
   13-gate SAT-synthesized circuit + explicit truth table (in compressor.py).
3. **Bezout/apply (`apply_bitvector`)**: replay dialog reversed on (z,0) with mod-p
   double/cond-add/cond-swap -> (0, x*z). The REVERSE (`apply_bitvector_reverse`) gives
   (0, x^-1 * z). This folds the modular multiply INTO the inversion.
4. `IPModMul` = in-place (x,y)->(x, y*x): ToBitVector(x) -> ApplyBitVector(d,y,0)=(0,xy)
   -> swap -> uncompute dialog. `.inverse()` = (x,y)->(x, y*x^-1) = the division used
   for lambda = dy/dx.
5. Point-add formula (Gouzien et al. / [GRLGS23]): lambda=(y1-y2)/(x1-x2);
   x3=lambda^2-x1-x2; y3=lambda(x1-x3)-y1. Uses 2 in-place modmuls + 1 modular square
   (special-prime optimized). For THIS benchmark the "window" is a single classical
   offset (no table lookup needed) -- simpler than the full windowed ECDLP circuit.

Space win: between forward and Bezout you hold only the dialog (~2.12n) + (r,s)(2n),
never the full (u,v,r,s)+history. => ~4.12-4.36n qubits vs our ~5.5n.

## Phase 0 (DONE -- GO)
`memory/validate_split_eea.py` is a standalone classical model (no Qarton dep).
Results on secp256k1 p, 10000 random inputs, ITERATIONS_VAR=2.4:
- apply_bitvector(z,0,dialog_of_x) == x*z mod p: 10000/10000 CORRECT
- GCD assertion failures: 0/10000 (paper target ~2^-13.3 ~ 1.2e-4)
- from_bitvector reconstructs (p,x) exactly.
NOTE: at 2.4 the per-inversion fail rate ~1e-4; over 9024 shots E[fails]~0.9, so to
RELIABLY pass all 9024 the iteration budget must be raised (ITERATIONS_VAR ~3-4) --
a tunable qubit/Toffoli-vs-failure knob. This is the sanctioned approximate-correctness
regime (the SOTA itself fails ~2^-13.3); NOT the phase-island gaming we rejected.

## Port plan (Phases 1-5, multi-session Rust build into src/point_add/)
Keep current build() as default; add new architecture behind an env flag until it
validates, then switch the default.
- **P1 Compressor**: port the 13-gate compress/uncompress (CX/CCX/X) + Swapper/Absorber.
  Unit-test vs the truth table. Small, self-contained. START HERE.
- **P2 Forward dialog circuit (ToBitVector)**: per-iter b0/b1 (reuse our Kaliski
  step0/step2 logic: with_eq_zero, with_gt), cond-swap (cswap), cond-sub, shift;
  Absorb (b0,b0&b1) into the compressed dialog. Validate gcd reaches (1,0) + dialog
  matches classical to_bitvector on 9024.
- **P3 Bezout/apply circuit (ApplyBitVector)**: per-iter mod-p double + cond-add +
  cond-swap, replaying dialog. Reuse our mod_add_qq/mod_double/Solinas. Validate
  x*z / x^-1*z on 9024.
- **P4 IPModMul wiring** + point-add formula (lambda, x3, y3) using our modmul/square +
  pseudo-Mersenne reduction. Validate full point-add on 9024.
- **P5 Optimize**: Gidney venting adder (pays off HERE -- no per-iter GCD adds),
  CDKM/Gidney hybrid, MSB-only comparators, drive qubit budget to ~1446. Tune
  ITERATIONS_VAR for reliable 9024-pass.

## Reusable from current Rust code (~40%)
harness/Fiat-Shamir flow; pseudo-Mersenne Solinas reduction; Karatsuba/schoolbook
modmul; measurement-based uncompute (hmr/cz_if); mod_double; controlled mod-add;
the per-iter b0/b1/cswap/sub/shift logic from kaliski_iteration (close to ToBitVector).

## Honest status
Phase 0 validated; full port is a genuine multi-day engineering build with a hard
9024-shot validation gate, but it is a PORT of published+open-source code, not
research. This is the only path to a leaderboard-competitive score.
