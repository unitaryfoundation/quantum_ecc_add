//! More detailed HMR diagnostics for the specialized bulk-prefix step.
//!
//! We already know the specialized step now matches the generic step in HMR
//! count at iter 0, but the operand order still diverges. This file records the
//! first few divergences and checks whether the same problem persists at later
//! bulk iterations.

use crate::circuit::{OperationType, QubitId};

use super::{B, N, SECP256K1_P, kaliski_iteration, kaliski_iteration_bulk_prefix3};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HmrOp {
    pub q: QubitId,
}

fn extract_hmrs(ops: &[crate::circuit::Op]) -> Vec<HmrOp> {
    ops.iter()
        .filter(|op| matches!(op.kind, OperationType::Hmr))
        .map(|op| HmrOp { q: op.q_target })
        .collect()
}

fn build_generic(iter_idx: usize) -> Vec<HmrOp> {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m = b.alloc_qubit();
    let f = b.alloc_qubit();
    kaliski_iteration(&mut b, SECP256K1_P, &u, &v, &r, &s, m, f, iter_idx);
    extract_hmrs(&b.ops)
}

fn build_special(iter_idx: usize) -> Vec<HmrOp> {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m = b.alloc_qubit();
    kaliski_iteration_bulk_prefix3(&mut b, SECP256K1_P, &u, &v, &r, &s, m, iter_idx);
    extract_hmrs(&b.ops)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report(iter_idx: usize) {
        let g = build_generic(iter_idx);
        let s = build_special(iter_idx);
        eprintln!("--- iter {} ---", iter_idx);
        eprintln!("generic count     : {}", g.len());
        eprintln!("specialized count : {}", s.len());
        let common = g.len().min(s.len());
        let mut diffs = 0usize;
        for i in 0..common {
            if g[i] != s[i] {
                diffs += 1;
                if diffs <= 8 {
                    eprintln!("diff[{}]: generic={:?} special={:?}", i, g[i], s[i]);
                }
            }
        }
        eprintln!("total differing positions in common prefix: {}", diffs);
    }

    #[test]
    fn compare_hmr_sequences_multiple_iters() {
        eprintln!("=== detailed HMR diffs ===");
        for &iter_idx in &[0usize, 1, 2, 3, 7, 15, 31] {
            report(iter_idx);
        }
        eprintln!("==========================");
        assert!(true);
    }
}
