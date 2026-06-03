//! Windowed/split-EEA architecture port (see memory/windowed-rewrite-plan.md).
//!
//! Phase 1: the **compressor** — packs 3 GCD-iteration bit-pairs (each in
//! {00,10,11}) into 5 bits, so the "dialog" costs ~2.12n instead of ~3n qubits.
//! 13-gate SAT-synthesized circuit from Schrottenloher (arXiv:2606.02235),
//! `compressor.py` (AGPL reference). This module is currently behind no build
//! path; it is the first validated brick of the rewrite.

#![allow(dead_code)]

use crate::circuit::QubitId;
use super::B;

/// A gate in the compressor, indexing into a 6-qubit register `q[0..6]`.
#[derive(Clone, Copy, Debug)]
pub enum CompGate {
    X(usize),
    Cx(usize, usize),
    Ccx(usize, usize, usize),
}

/// The 13-gate compressor sequence (maps a 6-bit input of the form
/// (00|10|11)^3 to a 5-bit output in q[0..5]; q[5] returns to 0 = clean ancilla).
pub fn compressor_gates() -> [CompGate; 13] {
    use CompGate::*;
    [
        Cx(1, 0), Cx(3, 2), Cx(5, 4),
        Cx(0, 2), Cx(5, 3), X(4),
        Ccx(1, 3, 5), Cx(1, 4), X(2),
        Ccx(3, 4, 5), Ccx(4, 5, 1), Ccx(2, 5, 0), Ccx(0, 1, 5),
    ]
}

/// Emit the compressor onto a builder, operating on the 6 given qubits.
pub fn emit_compressor(b: &mut B, q: &[QubitId; 6]) {
    for g in compressor_gates() {
        match g {
            CompGate::X(t) => b.x(q[t]),
            CompGate::Cx(c, t) => b.cx(q[c], q[t]),
            CompGate::Ccx(a, c, t) => b.ccx(q[a], q[c], q[t]),
        }
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Phase 2a: classical oracle for the split-EEA inversion (the function the
// forward dialog + Bezout circuits must reproduce). Ported from gcd_functions.py.
// Used as the test oracle for the circuit phases. Uses 512-bit ints to avoid
// overflow in the mod-p double/add of the Bezout pass.
// ───────────────────────────────────────────────────────────────────────────

type U512 = ruint::Uint<512, 8>;

const ITERATIONS_VAR: f64 = 2.4;
const U_PAD_VAR: f64 = 2.3;

/// 6-bit (form (00|10|11)^3) -> 5-bit compress, as the 27-entry table.
const TT27: [([u8; 6], [u8; 5]); 27] = [
    ([0,0,0,0,0,0],[0,0,1,0,1]), ([0,0,0,0,1,0],[0,0,1,0,0]), ([0,0,0,0,1,1],[0,0,1,1,1]),
    ([0,0,1,0,0,0],[0,0,0,0,1]), ([0,0,1,0,1,0],[0,0,0,0,0]), ([0,0,1,0,1,1],[0,0,0,1,1]),
    ([0,0,1,1,0,0],[1,1,1,1,1]), ([0,0,1,1,1,0],[0,0,1,1,0]), ([0,0,1,1,1,1],[1,1,1,0,1]),
    ([1,0,0,0,0,0],[1,0,0,0,1]), ([1,0,0,0,1,0],[1,0,0,0,0]), ([1,0,0,0,1,1],[1,0,0,1,1]),
    ([1,0,1,0,0,0],[1,0,1,0,1]), ([1,0,1,0,1,0],[1,0,1,0,0]), ([1,0,1,0,1,1],[1,0,1,1,1]),
    ([1,0,1,1,0,0],[1,1,0,1,1]), ([1,0,1,1,1,0],[1,0,0,1,0]), ([1,0,1,1,1,1],[1,1,0,0,1]),
    ([1,1,0,0,0,0],[0,1,1,0,0]), ([1,1,0,0,1,0],[0,1,1,0,1]), ([1,1,0,0,1,1],[0,1,1,1,0]),
    ([1,1,1,0,0,0],[0,1,0,0,0]), ([1,1,1,0,1,0],[0,1,0,0,1]), ([1,1,1,0,1,1],[0,1,0,1,0]),
    ([1,1,1,1,0,0],[1,1,1,1,0]), ([1,1,1,1,1,0],[0,1,1,1,1]), ([1,1,1,1,1,1],[1,1,1,0,0]),
];

fn compress6(d: [u8; 6]) -> [u8; 5] {
    TT27.iter().find(|(k, _)| *k == d).expect("invalid compress input").1
}
fn uncompress5(c: [u8; 5]) -> [u8; 6] {
    TT27.iter().find(|(_, v)| *v == c).expect("invalid uncompress input").0
}

fn iters_for(n: usize) -> usize {
    let nf = n as f64;
    (((1.413 * nf + ITERATIONS_VAR * nf.sqrt()) / 3.0).ceil() as usize) * 3
}

/// Forward Euclidean: returns the dialog garbage bits, or None on failure
/// (GCD didn't finish / overflowed the budget) — the approximate-correctness fail.
pub fn to_bitvector(u_in: U512, v_in: U512, n: usize) -> Option<Vec<u8>> {
    if u_in.bit(0) == false { return None; } // u must be odd
    let it = iters_for(n);
    let upad = (U_PAD_VAR * (n as f64).sqrt()).ceil() as f64;
    let ng = (it / 3) * 5;
    let mut g = vec![0u8; ng];
    for i in 0..(it / 3) {
        let c = compress6([0, 0, 0, 0, 0, 0]);
        g[5 * i..5 * i + 5].copy_from_slice(&c);
    }
    let (mut u, mut v) = (u_in, v_in);
    let two = U512::from(2u64);
    for i in 0..it {
        let b0 = if v.bit(0) { 1u8 } else { 0 };
        let b1 = if u > v { 1u8 } else { 0 };
        let mut small: [u8; 5] = g[(i / 3) * 5..(i / 3) * 5 + 5].try_into().unwrap();
        let mut dec = uncompress5(small);
        if dec[2 * (i % 3)] != 0 || dec[2 * (i % 3) + 1] != 0 { return None; }
        dec[2 * (i % 3)] = b0;
        dec[2 * (i % 3) + 1] = b0 & b1;
        small = compress6(dec);
        g[(i / 3) * 5..(i / 3) * 5 + 5].copy_from_slice(&small);
        let bound = (n as f64 - (i as f64) * 0.5 * (3.0 - 3.0_f64.log2()) + upad).max(0.0);
        if (u.bit_len() as f64) >= bound || (v.bit_len() as f64) >= bound { return None; }
        if (b0 & b1) == 1 { std::mem::swap(&mut u, &mut v); }
        if b0 == 1 { v -= u; }
        v >>= 1;
    }
    if v != U512::ZERO || u != U512::from(1u64) { return None; }
    let _ = two;
    Some(g)
}

fn addmod(a: U512, b: U512, p: U512) -> U512 { let s = a + b; if s >= p { s - p } else { s } }

/// Bezout "apply": apply_bitvector(z,0,dialog_of_x) = (·, x·z mod p) [multiply by x].
pub fn apply_bitvector(x: U512, y: U512, d: &[u8], p: U512) -> (U512, U512) {
    let (mut u, mut v) = (x % p, y % p);
    let nsteps = d.len() / 5 * 3;
    for i in (0..nsteps).rev() {
        let small: [u8; 5] = d[(i / 3) * 5..(i / 3) * 5 + 5].try_into().unwrap();
        let dec = uncompress5(small);
        let b0 = dec[2 * (i % 3)];
        let b0b1 = dec[2 * (i % 3) + 1];
        v = addmod(v, v, p);
        if b0 == 1 { v = addmod(v, u, p); }
        if b0b1 == 1 { std::mem::swap(&mut u, &mut v); }
    }
    (u, v)
}

#[cfg(test)]
mod oracle_tests {
    use super::*;

    fn secp_p() -> U512 {
        // 2^256 - 2^32 - 977
        (U512::from(1u64) << 256) - (U512::from(1u64) << 32) - U512::from(977u64)
    }

    // simple deterministic LCG for reproducible "random" inputs
    fn lcg(state: &mut u128) -> u128 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state
    }
    fn rand_mod(p: U512, st: &mut u128) -> U512 {
        let mut limbs = [0u64; 8];
        for k in 0..4 { limbs[k] = lcg(st) as u64; }
        U512::from_limbs(limbs) % p
    }

    #[test]
    fn split_eea_inversion_correct_and_low_failure() {
        let p = secp_p();
        let n = 256;
        let mut st: u128 = 0x1234_5678_9abc_def0;
        let (mut ok, mut fails) = (0u32, 0u32);
        let total = 2000;
        for _ in 0..total {
            let dx = (rand_mod(p, &mut st) % (p - U512::from(1u64))) + U512::from(1u64); // nonzero
            let dy = rand_mod(p, &mut st);
            // invert dx mod p: forward GCD on (p, dx)
            match to_bitvector(p, dx, n) {
                None => { fails += 1; }
                Some(dialog) => {
                    // apply_bitvector(dy,0) = dx*dy mod p (multiply); verify directly
                    let (_, prod) = apply_bitvector(dy, U512::ZERO, &dialog, p);
                    let expect = (dx * dy) % p;
                    assert_eq!(prod, expect, "split-EEA multiply wrong");
                    ok += 1;
                }
            }
        }
        // every completed case must be exactly correct; failures must be rare
        assert_eq!(ok + fails, total);
        assert!(fails * 1000 < total, "failure rate too high: {}/{}", fails, total);
        eprintln!("split-EEA oracle: {} correct, {} approximate-failures / {}", ok, fails, total);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simulate(state: [u8; 6]) -> [u8; 6] {
        let mut s = state;
        for g in compressor_gates() {
            match g {
                CompGate::X(t) => s[t] ^= 1,
                CompGate::Cx(c, t) => s[t] ^= s[c],
                CompGate::Ccx(a, c, t) => s[t] ^= s[a] & s[c],
            }
        }
        s
    }

    // The 27 valid (input -> output) pairs from compressor.py `_FUNCTION`.
    // input is 6 bits, output is the low 5 bits; bit 5 must return to 0.
    const TRUTH: [([u8; 6], [u8; 5]); 27] = [
        ([0,0,0,0,0,0],[0,0,1,0,1]), ([0,0,0,0,1,0],[0,0,1,0,0]), ([0,0,0,0,1,1],[0,0,1,1,1]),
        ([0,0,1,0,0,0],[0,0,0,0,1]), ([0,0,1,0,1,0],[0,0,0,0,0]), ([0,0,1,0,1,1],[0,0,0,1,1]),
        ([0,0,1,1,0,0],[1,1,1,1,1]), ([0,0,1,1,1,0],[0,0,1,1,0]), ([0,0,1,1,1,1],[1,1,1,0,1]),
        ([1,0,0,0,0,0],[1,0,0,0,1]), ([1,0,0,0,1,0],[1,0,0,0,0]), ([1,0,0,0,1,1],[1,0,0,1,1]),
        ([1,0,1,0,0,0],[1,0,1,0,1]), ([1,0,1,0,1,0],[1,0,1,0,0]), ([1,0,1,0,1,1],[1,0,1,1,1]),
        ([1,0,1,1,0,0],[1,1,0,1,1]), ([1,0,1,1,1,0],[1,0,0,1,0]), ([1,0,1,1,1,1],[1,1,0,0,1]),
        ([1,1,0,0,0,0],[0,1,1,0,0]), ([1,1,0,0,1,0],[0,1,1,0,1]), ([1,1,0,0,1,1],[0,1,1,1,0]),
        ([1,1,1,0,0,0],[0,1,0,0,0]), ([1,1,1,0,1,0],[0,1,0,0,1]), ([1,1,1,0,1,1],[0,1,0,1,0]),
        ([1,1,1,1,0,0],[1,1,1,1,0]), ([1,1,1,1,1,0],[0,1,1,1,1]), ([1,1,1,1,1,1],[1,1,1,0,0]),
    ];

    #[test]
    fn compressor_matches_truth_table() {
        for (inp, out) in TRUTH.iter() {
            let r = simulate(*inp);
            assert_eq!(&r[0..5], &out[..], "compress output wrong for input {:?}", inp);
            assert_eq!(r[5], 0, "ancilla (q[5]) not returned to 0 for input {:?}", inp);
        }
    }

    #[test]
    fn compressor_is_injective_on_valid_inputs() {
        // the 27 compressed outputs must be distinct (so uncompress is well-defined)
        let mut outs: Vec<[u8; 5]> = TRUTH.iter().map(|(_, o)| *o).collect();
        outs.sort();
        outs.dedup();
        assert_eq!(outs.len(), 27, "compressed outputs not injective");
    }
}
