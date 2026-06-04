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

/// The inverse compressor (reversed gate sequence): maps a 5-bit compressed
/// value (in q[0..5], q[5]=0) back to its 6-bit (00|10|11)^3 form.
fn compressor_inverse_gates() -> Vec<CompGate> {
    let mut g = compressor_gates().to_vec();
    g.reverse();
    g
}

/// Swapper(i): on a 6-qubit work register `w[0..6]` holding a compressed dialog
/// chunk (w[0..5] = compressed, w[5]=0 ancilla) and a 2-qubit pair `bb`, swap
/// `bb` with the (b0,b0&b1) pair at position `i ∈ {0,1,2}` of the chunk.
/// Sequence: uncompress (Compressor⁻¹) → swap bb↔(w[2i],w[2i+1]) → compress.
/// `Absorber` is the same circuit used when bb's target slot is known to be 0
/// (so bb is consumed into the dialog, no output) — identical gates.
fn swapper_gate_sim(i: usize, bb: [u8; 2], comp5: [u8; 5]) -> ([u8; 2], [u8; 5]) {
    // 6-qubit work register: [c0,c1,c2,c3,c4, anc=0]
    let mut w = [comp5[0], comp5[1], comp5[2], comp5[3], comp5[4], 0u8];
    let mut b = bb;
    let apply = |w: &mut [u8; 6], gates: &[CompGate]| {
        for g in gates {
            match *g {
                CompGate::X(t) => w[t] ^= 1,
                CompGate::Cx(c, t) => w[t] ^= w[c],
                CompGate::Ccx(a, c, t) => w[t] ^= w[a] & w[c],
            }
        }
    };
    apply(&mut w, &compressor_inverse_gates()); // now w[0..6] = uncompressed 6 bits
    std::mem::swap(&mut b[0], &mut w[2 * i]);
    std::mem::swap(&mut b[1], &mut w[2 * i + 1]);
    apply(&mut w, &compressor_gates()); // recompress; w[5] back to 0
    debug_assert_eq!(w[5], 0);
    (b, [w[0], w[1], w[2], w[3], w[4]])
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

// ───────────────────────────────────────────────────────────────────────────
// Phase 2b: the forward dialog (Euclidean) CIRCUIT, in our harness, built from
// exact reversible primitives (cuccaro_sub/add, cswap, compressor) so a clean
// basis-state op-sim is a correct judge. Matches the to_bitvector oracle.
// ───────────────────────────────────────────────────────────────────────────
use super::{cswap as kal_cswap, cuccaro_add, cuccaro_sub};

/// flag ^= (u > v), using a temp (n+1) register; u,v restored. Exact.
fn emit_gt(b: &mut B, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    let n = u.len();
    let tmp = b.alloc_qubits(n + 1); // holds v then v-u (extended)
    for j in 0..n {
        b.cx(v[j], tmp[j]);
    }
    let uhi = b.alloc_qubit(); // u extended high bit = 0
    let mut u_ext: Vec<QubitId> = u.to_vec();
    u_ext.push(uhi);
    let c_in = b.alloc_qubit();
    cuccaro_sub(b, &u_ext, &tmp, c_in); // tmp = v - u mod 2^(n+1); tmp[n]=borrow=(v<u)=(u>v)
    b.cx(tmp[n], flag);
    cuccaro_add(b, &u_ext, &tmp, c_in); // restore tmp = v
    b.free(c_in);
    b.free(uhi);
    for j in 0..n {
        b.cx(v[j], tmp[j]); // tmp -> 0
    }
    b.free_vec(&tmp);
}

/// Absorb (b0, b0b1) into a 5-bit compressed dialog chunk at position pos∈{0,1,2}.
/// chunk position pos must currently be (0,0) in uncompressed form. Uses 1 ancilla.
fn emit_absorb(b: &mut B, bb0: QubitId, bb1: QubitId, chunk: &[QubitId; 5], pos: usize) {
    let anc = b.alloc_qubit(); // 6th bit of the work reg, = 0
    let w = [chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], anc];
    // Compressor inverse (reversed gate list)
    for g in compressor_gates().iter().rev() {
        match *g {
            CompGate::X(t) => b.x(w[t]),
            CompGate::Cx(c, t) => b.cx(w[c], w[t]),
            CompGate::Ccx(a, c, t) => b.ccx(w[a], w[c], w[t]),
        }
    }
    // swap (bb0,bb1) into decompressed positions
    b.swap(bb0, w[2 * pos]);
    b.swap(bb1, w[2 * pos + 1]);
    // Compressor forward
    for g in compressor_gates() {
        match g {
            CompGate::X(t) => b.x(w[t]),
            CompGate::Cx(c, t) => b.cx(w[c], w[t]),
            CompGate::Ccx(a, c, t) => b.ccx(w[a], w[c], w[t]),
        }
    }
    b.free(anc);
}

/// Forward dialog: runs the binary-GCD on (u,v), recording the compressed dialog,
/// leaving (u,v) = (1,0). `dialog` must be 5*(iters/3) qubits, all 0 on entry.
/// Mirrors the to_bitvector oracle: per iter compute b0=v&1, b0b1=b0&(u>v);
/// apply (if b0b1 swap u,v); (if b0 v-=u); v>>=1; then absorb (b0,b0b1) into dialog.
fn forward_dialog(b: &mut B, u: &[QubitId], v: &[QubitId], dialog: &[QubitId], iters: usize) {
    let n = u.len();
    for i in 0..iters {
        let b0 = b.alloc_qubit();
        b.cx(v[0], b0); // b0 = v is odd
        let gt = b.alloc_qubit();
        emit_gt(b, u, v, gt); // gt = (u > v)
        let b0b1 = b.alloc_qubit();
        b.ccx(b0, gt, b0b1); // b0b1 = b0 & (u>v)
        emit_gt(b, u, v, gt); // uncompute gt -> 0
        b.free(gt);

        // apply ops BEFORE absorbing (absorb zeroes b0,b0b1)
        for j in 0..n {
            kal_cswap(b, b0b1, u[j], v[j]); // if b0b1: swap(u,v)
        }
        // if b0: v -= u  (load tmp = b0 & u, v -= tmp, unload)
        let tmp = b.alloc_qubits(n);
        for j in 0..n {
            b.ccx(b0, u[j], tmp[j]);
        }
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, &tmp, v, c_in); // v -= tmp
        b.free(c_in);
        for j in 0..n {
            b.ccx(b0, u[j], tmp[j]); // unload
        }
        b.free_vec(&tmp);
        // v >>= 1 (v[0] is now 0): rotate down
        for j in 0..n - 1 {
            b.swap(v[j], v[j + 1]);
        }

        // absorb (b0,b0b1) at chunk i/3, position i%3 (init chunk on i%3==0)
        let chunk: [QubitId; 5] = [
            dialog[5 * (i / 3)], dialog[5 * (i / 3) + 1], dialog[5 * (i / 3) + 2],
            dialog[5 * (i / 3) + 3], dialog[5 * (i / 3) + 4],
        ];
        if i % 3 == 0 {
            // initialize chunk to compress([0,0,0,0,0,0]) = [0,0,1,0,1]
            b.x(chunk[2]);
            b.x(chunk[4]);
        }
        emit_absorb(b, b0, b0b1, &chunk, i % 3);
        b.free(b0b1); // now 0
        b.free(b0); // now 0
    }
}

#[cfg(test)]
mod fwd_tests {
    use super::*;
    use crate::circuit::OperationType;
    use std::collections::HashMap;

    fn simulate(ops: &[crate::circuit::Op], init: &HashMap<u64, u8>) -> HashMap<u64, u8> {
        let mut q = init.clone();
        let g = |q: &HashMap<u64, u8>, id: u64| *q.get(&id).unwrap_or(&0);
        for op in ops {
            match op.kind {
                OperationType::CCX => {
                    let v = g(&q, op.q_control1.0) & g(&q, op.q_control2.0);
                    *q.entry(op.q_target.0).or_insert(0) ^= v;
                }
                OperationType::CX => {
                    let v = g(&q, op.q_control1.0);
                    *q.entry(op.q_target.0).or_insert(0) ^= v;
                }
                OperationType::X => {
                    *q.entry(op.q_target.0).or_insert(0) ^= 1;
                }
                OperationType::Swap => {
                    let a = g(&q, op.q_control1.0);
                    let b = g(&q, op.q_target.0);
                    q.insert(op.q_control1.0, b);
                    q.insert(op.q_target.0, a);
                }
                OperationType::R | OperationType::Hmr => {
                    // reset-on-free: a clean free leaves the qubit at 0
                    q.insert(op.q_target.0, 0);
                }
                OperationType::AppendToRegister | OperationType::Register => {}
                other => panic!("forward_dialog emitted non-exact op {:?}", other),
            }
        }
        q
    }

    #[test]
    fn forward_dialog_matches_oracle() {
        let n = 12usize;
        let p: u64 = 4093; // 12-bit prime
        let iters = iters_for(n);
        for v0 in [1u64, 2, 7, 100, 1234, 4000] {
            let b = &mut super::super::B::new();
            let u = b.alloc_qubits(n);
            let v = b.alloc_qubits(n);
            let dialog = b.alloc_qubits(5 * (iters / 3));
            forward_dialog(b, &u, &v, &dialog, iters);
            let mut init: HashMap<u64, u8> = HashMap::new();
            for j in 0..n {
                if (p >> j) & 1 == 1 { init.insert(u[j].0, 1); }
                if (v0 >> j) & 1 == 1 { init.insert(v[j].0, 1); }
            }
            let q = simulate(&b.ops, &init);
            let g = |id: u64| *q.get(&id).unwrap_or(&0);
            let uout: u64 = (0..n).map(|j| (g(u[j].0) as u64) << j).sum();
            let vout: u64 = (0..n).map(|j| (g(v[j].0) as u64) << j).sum();
            let dialog_bits: Vec<u8> = dialog.iter().map(|qb| g(qb.0)).collect();
            let expect = to_bitvector(U512::from(p), U512::from(v0), n)
                .expect("oracle converges");
            assert_eq!(uout, 1, "u != 1 for v0={}", v0);
            assert_eq!(vout, 0, "v != 0 for v0={}", v0);
            assert_eq!(dialog_bits, expect, "dialog mismatch for v0={}", v0);
        }
    }
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
    fn swapper_matches_reference_semantics() {
        // reference Swapper.dummy_classical_function (compressor.py)
        for i in 0..3 {
            for bb in [[0u8, 0], [1, 0], [1, 1]] {
                for (inp, comp) in TRUTH.iter() {
                    // `comp` is a valid compressed value; its uncompression is `inp`
                    let dec = *inp; // uncompress(comp) == inp by construction of TRUTH
                    let mut expect_dec = dec;
                    let new_bb = [dec[2 * i], dec[2 * i + 1]];
                    expect_dec[2 * i] = bb[0];
                    expect_dec[2 * i + 1] = bb[1];
                    let expect_comp = super::compress6(expect_dec);
                    let (got_bb, got_comp) = super::swapper_gate_sim(i, bb, *comp);
                    assert_eq!(got_bb, new_bb, "swapper bb out wrong i={} bb={:?} comp={:?}", i, bb, comp);
                    assert_eq!(got_comp, expect_comp, "swapper comp out wrong i={} bb={:?}", i, bb);
                }
            }
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
