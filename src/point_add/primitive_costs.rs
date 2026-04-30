//! Measure the exact Toffoli cost of each modular arithmetic primitive in
//! isolation. Test-only; emits numbers via eprintln for the planner.
//!
//! We don't need these for live correctness — just for honest cost accounting
//! so we can sanity-check SOTA reachability.

#![cfg(test)]

use super::{
    mod_add_qb, mod_add_qc, mod_add_qq, mod_double_inplace_fast, mod_halve_inplace_fast,
    mod_mul_add_into_acc_schoolbook,
    mod_mul_sub_qq, mod_mul_write_into_zero_acc_schoolbook, mod_neg_inplace_fast, mod_sub_qb, N,
    SECP256K1_P,
};
use super::{B, QubitId};
use crate::circuit::OperationType;

enum ShiftUndoForCost {
    Doubles(usize),
    Chunk(usize, Vec<QubitId>, QubitId, QubitId),
}

fn shift_tmp_up_for_sparse_const_cost(
    b: &mut B,
    tmp: &[QubitId],
    p: alloy_primitives::U256,
    mut delta: usize,
    undo: &mut Vec<ShiftUndoForCost>,
) {
    while delta >= 22 {
        let (spill, flag_inv, ovf) = super::mod_shift_left_by_k(b, tmp, p, 22);
        undo.push(ShiftUndoForCost::Chunk(22, spill, flag_inv, ovf));
        delta -= 22;
    }
    if delta >= 12 {
        let (spill, flag_inv, ovf) = super::mod_shift_left_by_k(b, tmp, p, delta);
        undo.push(ShiftUndoForCost::Chunk(delta, spill, flag_inv, ovf));
    } else if delta > 0 {
        for _ in 0..delta {
            mod_double_inplace_fast(b, tmp, p);
        }
        undo.push(ShiftUndoForCost::Doubles(delta));
    }
}

fn undo_sparse_const_shifts_for_cost(
    b: &mut B,
    tmp: &[QubitId],
    p: alloy_primitives::U256,
    undo: Vec<ShiftUndoForCost>,
) {
    for item in undo.into_iter().rev() {
        match item {
            ShiftUndoForCost::Doubles(k) => {
                for _ in 0..k {
                    mod_halve_inplace_fast(b, tmp, p);
                }
            }
            ShiftUndoForCost::Chunk(k, spill, flag_inv, ovf) => {
                super::mod_shift_right_by_k(b, tmp, p, k, spill, flag_inv, ovf);
            }
        }
    }
}

fn mul_by_const_acc_chunked_shifts_for_cost(
    b: &mut B,
    x: &[QubitId],
    c: alloy_primitives::U256,
    acc: &[QubitId],
    p: alloy_primitives::U256,
) {
    let n = x.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n {
        b.cx(x[i], tmp[i]);
    }
    let mut positions = Vec::new();
    for i in 0..256 {
        if super::bit(c, i) {
            positions.push(i);
        }
    }
    let mut undo = Vec::new();
    let mut cur = 0usize;
    for pos in positions {
        shift_tmp_up_for_sparse_const_cost(b, &tmp, p, pos - cur, &mut undo);
        cur = pos;
        mod_add_qq(b, acc, &tmp, p);
    }
    undo_sparse_const_shifts_for_cost(b, &tmp, p, undo);
    for i in 0..n {
        b.cx(x[i], tmp[i]);
    }
    b.free_vec(&tmp);
}

fn count_ccx(ops: &[crate::circuit::Op]) -> usize {
    ops.iter()
        .filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
        .count()
}

fn new_builder_with_reg(n: usize) -> (B, Vec<QubitId>) {
    let mut b = B::new();
    let r = b.alloc_qubits(n);
    (b, r)
}

#[test]
fn cost_mul_write_schoolbook_n256() {
    let mut b = B::new();
    let p = SECP256K1_P;
    let acc = b.alloc_qubits(N);
    let x = b.alloc_qubits(N);
    let y = b.alloc_qubits(N);
    let start = b.ops.len();
    mod_mul_write_into_zero_acc_schoolbook(&mut b, &acc, &x, &y, p);
    let end = b.ops.len();
    let ccx = count_ccx(&b.ops[start..end]);
    eprintln!("mod_mul_write_into_zero_acc_schoolbook(n=256): {} CCX", ccx);
}

#[test]
fn cost_mul_add_schoolbook_n256() {
    let mut b = B::new();
    let p = SECP256K1_P;
    let acc = b.alloc_qubits(N);
    let x = b.alloc_qubits(N);
    let y = b.alloc_qubits(N);
    let start = b.ops.len();
    mod_mul_add_into_acc_schoolbook(&mut b, &acc, &x, &y, p);
    let end = b.ops.len();
    let ccx = count_ccx(&b.ops[start..end]);
    eprintln!("mod_mul_add_into_acc_schoolbook(n=256): {} CCX", ccx);
}

#[test]
fn cost_mul_sub_qq_n256() {
    let mut b = B::new();
    let p = SECP256K1_P;
    let acc = b.alloc_qubits(N);
    let x = b.alloc_qubits(N);
    let y = b.alloc_qubits(N);
    let start = b.ops.len();
    mod_mul_sub_qq(&mut b, &acc, &x, &y, p);
    let end = b.ops.len();
    let ccx = count_ccx(&b.ops[start..end]);
    eprintln!("mod_mul_sub_qq(n=256): {} CCX", ccx);
}

#[test]
fn cost_sub_qb_n256() {
    let mut b = B::new();
    let p = SECP256K1_P;
    let acc = b.alloc_qubits(N);
    let bits = b.alloc_bits(N);
    let start = b.ops.len();
    mod_sub_qb(&mut b, &acc, &bits, p);
    let end = b.ops.len();
    let ccx = count_ccx(&b.ops[start..end]);
    eprintln!("mod_sub_qb(n=256): {} CCX", ccx);
}

#[test]
fn cost_neg_inplace_fast_n256() {
    let (mut b, r) = new_builder_with_reg(N);
    let p = SECP256K1_P;
    let start = b.ops.len();
    mod_neg_inplace_fast(&mut b, &r, p);
    let end = b.ops.len();
    let ccx = count_ccx(&b.ops[start..end]);
    eprintln!("mod_neg_inplace_fast(n=256): {} CCX", ccx);
}
#[test]
fn cost_squaring_sub_n256() {
    use super::*;
    use crate::circuit::OperationType;
    fn count_ccx(ops: &[crate::circuit::Op]) -> usize {
        ops.iter().filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ)).count()
    }
    let mut b = B::new();
    let p = SECP256K1_P;
    let acc = b.alloc_qubits(N);
    let x = b.alloc_qubits(N);
    let start = b.ops.len();
    // mod_mul_sub_qq with same register is a squaring
    mod_mul_sub_qq(&mut b, &acc, &x, &x, p);
    let end = b.ops.len();
    let ccx = count_ccx(&b.ops[start..end]);
    eprintln!("squaring via mod_mul_sub_qq: {} CCX", ccx);
}

#[test]
fn fermat_fixed_chain_inversion_floor_misses_sota_by_order() {
    // Branchless inversion by Fermat/exponentiation is the obvious way to avoid
    // Euclidean branch histories.  But even an unrealistically optimal addition
    // chain for an exponent near 2^256 needs at least 255 modular
    // square/multiply layers (each layer can at most double the exponent).  With
    // the measured current n=256 modular-square floor, this is already tens of
    // millions of CCX per inverse before any Bennett cleanup, scratch pressure,
    // or the second point-add denominator.  So fixed-sequence exponentiation is
    // not the missing SOTA-shaped DIV/IMUL primitive.
    let mut b = B::new();
    let p = SECP256K1_P;
    let acc = b.alloc_qubits(N);
    let x = b.alloc_qubits(N);
    let start = b.ops.len();
    mod_mul_sub_qq(&mut b, &acc, &x, &x, p);
    let square_ccx = count_ccx(&b.ops[start..]);
    let chain_layer_lower_bound = 255usize;
    let inv_floor = square_ccx * chain_layer_lower_bound;
    println!("METRIC fermat_inv_square_floor_ccx={square_ccx}");
    println!("METRIC fermat_inv_chain_floor_ccx={inv_floor}");
    eprintln!(
        "Fermat inversion floor: square_ccx={square_ccx}, layers>={chain_layer_lower_bound}, inv_floor={inv_floor}"
    );
    assert!(inv_floor > 30_000_000);
}

#[test]
fn cost_halve_double_n256() {
    let mut b = B::new();
    let p = SECP256K1_P;
    let v = b.alloc_qubits(N);
    let start = b.ops.len();
    mod_halve_inplace_fast(&mut b, &v, p);
    let mid = b.ops.len();
    mod_double_inplace_fast(&mut b, &v, p);
    let end = b.ops.len();
    let halve_ccx = count_ccx(&b.ops[start..mid]);
    let double_ccx = count_ccx(&b.ops[mid..end]);
    eprintln!("mod_halve_inplace_fast(n=256): {} CCX", halve_ccx);
    eprintln!("mod_double_inplace_fast(n=256): {} CCX", double_ccx);
}

#[test]
fn chunked_shift_prescaler_reopens_small_scale_absorption_win() {
    // Scale absorption deletes a ~iters-long halve/double correction loop if we
    // initialize Kaliski with 2^iters*x.  The constants are sparse for secp256k1,
    // e.g. 2^404 = 2^148(2^32+977), so try a custom constant multiplier that
    // jumps between sparse set-bit positions with the Solinas k-bit shifter
    // instead of walking through every intermediate double.  This beats the old
    // mixed prescaler locally and is just below the correction-loop cost for the
    // current pair1/pair2 iteration counts, making scale absorption a small but
    // real env-gated integration candidate.
    use super::*;
    let p = SECP256K1_P;
    let x = B::new();
    drop(x);
    for &(iters, label) in &[(404usize, "pair1"), (401usize, "pair2")] {
        let scale = pow_mod_2_k(p, iters);
        let mut b = B::new();
        let src = b.alloc_qubits(N);
        let acc = b.alloc_qubits(N);
        let start = b.ops.len();
        mul_by_const_acc_exact_adds_fast_shifts(&mut b, &src, scale, &acc, p, false);
        let mixed_ccx = count_ccx(&b.ops[start..]);

        let mut b = B::new();
        let src = b.alloc_qubits(N);
        let acc = b.alloc_qubits(N);
        let start = b.ops.len();
        mul_by_const_acc_chunked_shifts_for_cost(&mut b, &src, scale, &acc, p);
        let chunked_ccx = count_ccx(&b.ops[start..]);

        let mut b = B::new();
        let v = b.alloc_qubits(N);
        let start = b.ops.len();
        for _ in 0..iters {
            if label == "pair1" {
                mod_halve_inplace_fast(&mut b, &v, p);
            } else {
                mod_double_inplace_fast(&mut b, &v, p);
            }
        }
        let correction_loop_ccx = count_ccx(&b.ops[start..]);
        let projected_delta = 2isize * chunked_ccx as isize - correction_loop_ccx as isize;
        eprintln!(
            "{label} scale prescaler: mixed_ccx={mixed_ccx}, chunked_ccx={chunked_ccx}, correction_loop_ccx={correction_loop_ccx}, projected_delta={projected_delta}"
        );
        println!("METRIC scale_absorb_{label}_mixed_prescale_ccx={mixed_ccx}");
        println!("METRIC scale_absorb_{label}_chunked_prescale_ccx={chunked_ccx}");
        println!("METRIC scale_absorb_{label}_correction_loop_ccx={correction_loop_ccx}");
        println!("METRIC scale_absorb_{label}_chunked_projected_delta={projected_delta}");
        assert!(chunked_ccx < mixed_ccx / 2, "chunked sparse shifts should strongly improve the local prescaler");
        assert!(projected_delta < 0, "chunked compute+uncompute should beat the deleted correction loop locally");
    }
}

#[test]
fn profile_point_add_by_phase() {
    use std::collections::HashMap;
    use crate::circuit::OperationType;
    let mut b = B::new();
    let p = SECP256K1_P;
    let n = 256;
    let tx = b.alloc_qubits(n);
    let ty = b.alloc_qubits(n);
    let ox = b.alloc_bits(n);
    let oy = b.alloc_bits(n);
    super::build_standard_point_add(&mut b, &tx, &ty, &ox, &oy, p);

    let mut phase_ccx: HashMap<&str, usize> = HashMap::new();
    let mut current_phase: &str = "(none)";
    let trans = &b.phase_transitions;
    let mut ti = 0;
    for (idx, op) in b.ops.iter().enumerate() {
        while ti < trans.len() && trans[ti].0 <= idx {
            current_phase = trans[ti].1;
            ti += 1;
        }
        if matches!(op.kind, OperationType::CCX | OperationType::CCZ) {
            *phase_ccx.entry(current_phase).or_insert(0) += 1;
        }
    }

    let mut entries: Vec<_> = phase_ccx.into_iter().collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1));
    let mut total = 0usize;
    eprintln!("\n=== Point Add Toffoli Profile by Phase ===");
    for (phase, ccx) in &entries {
        total += ccx;
        eprintln!("{:>10} {}", ccx, phase);
    }
    eprintln!("{:>10} TOTAL", total);
}
