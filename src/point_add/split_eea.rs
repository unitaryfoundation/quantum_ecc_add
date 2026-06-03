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
