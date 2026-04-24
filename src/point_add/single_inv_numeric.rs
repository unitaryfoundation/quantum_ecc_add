//! Classical numeric replays of single-inversion point-add strategies.
//!
//! Discipline: every strategy here is a *reversible-schedule simulation*.
//! Each `replay_strategy_X()` mirrors, step by step, what a quantum
//! scaffold would do — tracked as a U256 per "register". The end state
//! MUST match the reference `WeierstrassEllipticCurve::add`, across 200
//! random curve-point pairs. A strategy that doesn't pass 200/200 is
//! dead, full stop.
//!
//! No strategy graduates to reversible code in `mod.rs` unless it
//! passes and its op count beats the current 2-Kaliski scaffold.
//!
//! See `single_inv_plan.md` for the prose spec / register tables.

#![cfg(test)]

use alloy_primitives::U256;

use super::SECP256K1_P;

fn sub_mod(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        (a - b) % p
    } else {
        p - ((b - a) % p)
    }
}

fn neg_mod(a: U256, p: U256) -> U256 {
    sub_mod(U256::ZERO, a, p)
}

// ─────────────────────────────────────────────────────────────────────
// Reference: the established, classically-correct formula (sanity).
// ─────────────────────────────────────────────────────────────────────
pub fn single_inv_add(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let dx = sub_mod(px, qx, p);
    let dy = sub_mod(py, qy, p);
    let inv_dx = dx.inv_mod(p).expect("dx nonzero");
    let lam = dy.mul_mod(inv_dx, p);
    let lam2 = lam.mul_mod(lam, p);
    let rx = sub_mod(sub_mod(lam2, px, p), qx, p);
    let ry = sub_mod(lam.mul_mod(sub_mod(qx, rx, p), p), qy, p);
    (rx, ry)
}

// ─────────────────────────────────────────────────────────────────────
// Kaliski raw output convention (settled in 21b87fd):
//   inside a `with_kal_inv_raw(v)` body the inv_raw register holds
//     -v^{-1} * 2^{2n-1} mod p
//   (sign negative, scale 2^{2n-1} because iters = 2n-1 and Kaliski
//    doubles r every iter; the final positive-ation is skipped).
// ─────────────────────────────────────────────────────────────────────
fn kaliski_body_inv_raw(v: U256, p: U256) -> U256 {
    const TWO_2N_MINUS_1_EXP: u64 = 2 * 256 - 1;
    let two = U256::from(2);
    let scale = two.pow_mod(U256::from(TWO_2N_MINUS_1_EXP), p);
    let inv_v = v.inv_mod(p).expect("v nonzero");
    neg_mod(inv_v, p).mul_mod(scale, p)
}

fn pow2_mod(e: i64, p: U256) -> U256 {
    let two = U256::from(2);
    if e >= 0 {
        two.pow_mod(U256::from(e as u64), p)
    } else {
        two.pow_mod(U256::from((-e) as u64), p)
            .inv_mod(p)
            .expect("2 invertible")
    }
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY A
//
// Plan (from single_inv_plan.md §4):
//   1. tx = dx, ty = dy
//   2. a := tx * ty                            (fresh register)
//   3. run Kaliski on a, body entered with inv_raw_a = -a^{-1}*2^{2n-1}
//   4. inside body, lam := ty^2 * inv_raw_a * 2^{-(2n-1)} = -λ
//      (implemented as: dy_sq then mul by inv_raw then halve 2n-1 times)
//   5. Rx fold into tx: tx -= λ²; +3Qx; neg; +Qx  → tx = Rx - Qx.
//   6. ty += lam · tx; then ty -= Qy.          ← This is where the
//      ty=dy contamination has to clear. Expected: Ry mismatch.
// ─────────────────────────────────────────────────────────────────────
pub fn replay_strategy_a(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let mut tx = sub_mod(px, qx, p);
    let mut ty = sub_mod(py, qy, p);

    // step 2: a = tx·ty
    let a = tx.mul_mod(ty, p);

    // step 3: enter kaliski body
    let inv_raw_a = kaliski_body_inv_raw(a, p);

    // step 4: lam = ty^2 · inv_raw_a · 2^{-(2n-1)} = -λ.
    let scale_back = pow2_mod(-(2 * 256 - 1), p);
    let lam = ty
        .mul_mod(ty, p)
        .mul_mod(inv_raw_a, p)
        .mul_mod(scale_back, p);
    // sanity: lam == -dy/dx
    let lam_expected = neg_mod(ty.mul_mod(tx.inv_mod(p).unwrap(), p), p);
    assert_eq!(lam, lam_expected, "strategy A: lam must equal -λ");

    // step 5: Rx fold.  Goal: tx = Rx when done.
    //   tx := dx - λ² ; tx += 2Qx ; tx := -tx
    //      → tx = -(dx - λ² + 2Qx) = λ² - dx - 2Qx
    //          = λ² - (Px - Qx) - 2Qx = λ² - Px - Qx = Rx. ✓
    let lam2 = lam.mul_mod(lam, p);
    tx = sub_mod(tx, lam2, p);
    tx = tx.add_mod(qx.mul_mod(U256::from(2), p), p);
    tx = neg_mod(tx, p);

    // At this point tx = Rx. For Ry we need (Rx - Qx):
    let rx_minus_qx = sub_mod(tx, qx, p);

    // step 6: ty += lam·(Rx - Qx) - Qy.  lam = -λ, so
    //   lam · (Rx - Qx) = (-λ)(Rx - Qx) = λ(Qx - Rx) = Ry + Qy.
    //   ty = dy + (Ry + Qy) = (Py - Qy) + Ry + Qy = Py + Ry.
    //   ty -= Qy → ty = Py + Ry - Qy.
    // That is the predicted contamination.
    ty = ty.add_mod(lam.mul_mod(rx_minus_qx, p), p);
    ty = sub_mod(ty, qy, p);

    (tx, ty)
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY B2
//
// Plan (from §5, the version that only runs ONE Kaliski pass AND does
//       Rx/Ry computation inside the Kaliski body so lam can be
//       uncomputed via mul2-inverse before exit):
//
//   Body over tx = dx, with inv_raw = -dx^{-1} · 2^{2n-1}.
//   1. lam := ty · inv_raw                    ← lam = -dy·dx^{-1}·2^{2n-1}
//   2. halve lam (2n-1) times                 ← lam = -λ
//   3. ty += lam · tx                         ← ty := dy - λ·dx = 0 ✅
//      (at this point ty is zero and Py is gone from the state)
//   4. tx ← Rx - Qx via fold (uses lam² only)
//   5. ty += lam · tx                         ← ty := 0 + (-λ)(Rx-Qx)
//                                                  = λ(Qx - Rx) = Ry + Qy
//   6. ty -= Qy                               ← ty := Ry ✅
//   7. Now reverse-mul and reverse-halve and reverse-mul to put lam
//      back at 0 using only current state.
//
//   Step 7 details: we're still inside the Kaliski body (inv_raw, tx=dx,
//   ty=Ry, lam=-λ all live; kal_state also live). To free lam (=-λ):
//     - un-fold tx: reverse the Rx fold so tx returns to dx.
//     - subtract `ty_before_Ry_step` back from ty. We don't have it, but
//       we can recompute it: it was (Ry + Qy) prior to the -Qy. So
//       `ty += Qy` → ty = Ry + Qy. Then `ty -= lam · tx` → ty = 0. ✅
//     - But we need ty = dy at body exit for the Kaliski_backward pass
//       to close out correctly? Actually no — Kaliski_backward only
//       reverses the Kaliski state, it doesn't touch ty. ty can be
//       anything at body exit.
//     - So we can actually just LEAVE ty = Ry, skip the reverse steps
//       for ty, and only uncompute lam.
//     - To uncompute lam: reverse halves (2n-1 doublings) brings lam
//       back to `-dy·dx^{-1}·2^{2n-1} = -dy·inv_raw`. Then the inverse
//       of `lam += ty·inv_raw` is `lam -= ty·inv_raw`. But here ty has
//       changed from dy to Ry! So we need to do uncomputation BEFORE ty
//       is changed — i.e. before step 3.
//
// Revised Strategy B2:
//   1. lam := ty · inv_raw = -dy·dx^{-1}·2^{2n-1}
//   2. halve lam (2n-1) times  → lam = -λ
//   3. Take a Bennett snapshot of λ:
//        lam2 := 0
//        lam2 := lam · 1 (just cx-copy)
//   Actually cx-copy of lam into a fresh register is free in qubits
//   (just n qubits, which is what we'd pay anyway). Let's do that.
//        lam_out := lam   (fresh register, cx-copy, classical)
//   4. reverse: double lam (2n-1) times → lam = -dy·inv_raw
//   5. lam -= ty · inv_raw       → lam = 0, freeable ✅
//   6. But now we still have lam_out = -λ live. Proceed with it.
//   7. tx ← Rx - Qx fold using lam_out² (= λ²).
//   8. ty += lam_out · tx        → ty = dy + (-λ)(Rx-Qx)
//                                    = dy + λ(Qx-Rx)
//                                    = (Py-Qy) + (Ry+Qy) = Py + Ry. ✗
//
// SAME OBSTRUCTION. Ry replaces Qy but Py stays. The classical replay
// will catch this as an Ry mismatch.
//
// The ONLY way to avoid the Py contamination is to zero ty BEFORE
// computing Ry. That's what the live 2-Kaliski scaffold does via
// pair1_mul2 (ty := 0) then mul3_between_pair (ty := λ(Rx-Qx)) then
// pair2_cleanup (ty -= Qy). So inside a single Kaliski we MUST:
//   - zero ty inside the body (step 3 above), then
//   - compute Ry into the now-zero ty WITHOUT needing λ outside the body.
// That requires the Rx fold to happen inside the body too, and lam to
// STILL be live when we compute Ry.
//
// Final revised Strategy B2 (what we actually run):
//   body:
//     1. lam := ty · inv_raw                   (lam = -dy·dx^{-1}·2^{2n-1})
//     2. halve lam (2n-1)×                     (lam = -λ)
//     3. ty += lam · tx                        (ty = 0)
//     4. Rx fold in tx using lam²              (tx = Rx - Qx)
//     5. ty += lam · tx                        (ty = (-λ)(Rx-Qx) = Ry+Qy)
//     6. ty -= Qy                              (ty = Ry)
//     7. un-fold tx (reverse of step 4)         (tx = dx, lam² available)
//     8. double lam (2n-1)×                     (lam = -dy·dx^{-1}·2^{2n-1})
//     9. lam -= ty · inv_raw                    (...but ty=Ry not dy!)
//
//   Step 9 again: we need ty = dy to undo step 1, and ty = Ry now.
//   Ry - dy relation: Ry = λ(Qx-Rx) - Qy = λ(Qx - Rx) - Qy
//                     = -λ(Rx - Qx) - Qy.
//   No usable relation between Ry and dy without knowing λ and Rx.
//
//   The escape: uncompute lam via an ancilla BEFORE ty changes.
//
//   Fresh final Strategy B2:
//     body:
//      1. lam := ty · inv_raw                  (lam = -dy·inv_raw)
//      2. halve lam (2n-1)×                    (lam = -λ)
//      3. lam_copy(n) := 0 ; lam_copy ^= lam   (cx-copy; lam_copy = -λ)
//      4. double lam (2n-1)×                   (lam = -dy·inv_raw again)
//      5. lam -= ty · inv_raw                  (lam = 0, free it) ✅
//      6. ty += lam_copy · tx                  (ty = dy - λ·dx = 0) ✅
//      7. Rx fold in tx using lam_copy²        (tx = Rx - Qx)
//      8. ty += lam_copy · tx                  (ty = Ry + Qy)
//      9. ty -= Qy                             (ty = Ry)
//     10. un-fold tx                            (tx = dx briefly?)
//     11. Reverse lam_copy: …
//
//   At end we must free lam_copy. Its current value is -λ. To uncompute
//   without needing dy: compute `lam_copy -= lam_recomputed`, where
//   lam_recomputed = -dy·dx^{-1}·? ... no, dy is gone.
//
//   But we CAN recompute -λ using the current live registers:
//   tx = dx (after step 10), Qy, ty = Ry, Qx, Px (from tx+Qx).
//     Ry + Qy = λ(Qx - Rx) = -λ(Rx - Qx)
//     So λ = -(Ry + Qy) / (Rx - Qx).
//   That requires inverting (Rx - Qx). We already inverted dx. Can we
//   get (Rx - Qx)^{-1} from dx^{-1}?
//     Rx - Qx = λ² - dx - 2Qx, where λ = dy·dx^{-1}
//            = dy²/dx² - dx - 2Qx
//            = (dy² - dx³ - 2Qx·dx²) / dx²
//     So (Rx - Qx) = (dy² - dx·(dx² + 2Qx·dx)) / dx²
//              1 / (Rx - Qx) = dx² / (dy² - dx·(dx² + 2Qx·dx))
//   That's a second inversion (of dy² - dx·(dx² + 2Qx·dx)) unless we
//   have other structure. Dead for single-Kaliski.
//
// *** Therefore the classical replay below implements the REVISED
// Strategy B2 and expects Ry to fail, pinpointing the obstruction at
// step 11. ***
// ─────────────────────────────────────────────────────────────────────
pub fn replay_strategy_b2(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let mut tx = sub_mod(px, qx, p);
    let mut ty = sub_mod(py, qy, p);
    let dx_original = tx;

    // enter kaliski body on tx = dx.
    let inv_raw = kaliski_body_inv_raw(dx_original, p);

    // (1) lam := ty · inv_raw
    let mut lam = ty.mul_mod(inv_raw, p);
    // (2) halve lam 2n-1 times
    lam = lam.mul_mod(pow2_mod(-(2 * 256 - 1), p), p);
    // Now lam = -λ.

    // (3) cx-copy lam into lam_copy (reversibly, 0 Toffoli, just tracking)
    let lam_copy = lam;
    // (4) double lam 2n-1 times → lam = -dy·inv_raw again
    lam = lam.mul_mod(pow2_mod(2 * 256 - 1, p), p);
    // (5) lam -= ty · inv_raw  (ty is still dy here)
    lam = sub_mod(lam, ty.mul_mod(inv_raw, p), p);
    assert_eq!(lam, U256::ZERO, "lam should be zero after uncompute");

    // (6) ty += lam_copy · tx = dy + (-λ)(dx) = 0.
    ty = ty.add_mod(lam_copy.mul_mod(tx, p), p);
    assert_eq!(ty, U256::ZERO, "ty should be zero after step 6");

    // (7) Rx fold in tx using lam_copy².  Goal tx := Rx.
    let lam_sq = lam_copy.mul_mod(lam_copy, p);
    tx = sub_mod(tx, lam_sq, p);
    tx = tx.add_mod(qx.mul_mod(U256::from(2), p), p);
    tx = neg_mod(tx, p); // tx = Rx

    // For Ry we need Rx - Qx:
    let rx_minus_qx = sub_mod(tx, qx, p);

    // (8) ty += lam_copy · (Rx - Qx) = (-λ)(Rx - Qx) = λ(Qx - Rx) = Ry + Qy.
    //     But ty is 0 here (step 6), so ty := Ry + Qy.
    ty = ty.add_mod(lam_copy.mul_mod(rx_minus_qx, p), p);
    // (9) ty -= Qy → ty = Ry.
    ty = sub_mod(ty, qy, p);

    // ANCILLA LEAK: lam_copy still = -λ. Must be freed. In classical
    // replay we don't have to zero it — but as a reversible scaffold it
    // would need to. The open question is whether `lam_copy` can be
    // uncomputed from {tx=Rx, ty=Ry, ox=Qx, oy=Qy, dx_original} alone.
    // In the classical trace we just note the leak and return.
    let _ = dx_original;
    let _ = lam_copy;
    (tx, ty)
}

// ─────────────────────────────────────────────────────────────────────
// STRATEGY C — Montgomery batch: invert a single product w that yields
// both 1/dx and 1/(Rx-Qx).
//
// From §6: `dx² · (Rx - Qx) = dy² - dx² · (Px + Qx)`, so
//          `dx · (Rx - Qx) = (dy² - dx² · (Px + Qx)) / dx`.
// We can compute the numerator `u = dy² - dx² · (Px + Qx)` classically
// from live registers {dx, dy, Px_const? no Px quantum, Qx}. Px is
// quantum. But Px = dx + Qx. So:
//   u = dy² - dx²·((dx + Qx) + Qx) = dy² - dx³ - 2·Qx·dx²
// That's a polynomial in dx, dy, Qx. Fully computable from the live
// quantum registers.
//
// Letting `w = u / dx = (dy² - dx³ - 2·Qx·dx²) / dx = dy²/dx - dx² - 2·Qx·dx`
// we see w = dx·(Rx - Qx). We need w^{-1}. Then:
//   1/dx        = w^{-1} · (Rx - Qx) = w^{-1} · u/dx ...  hmm still dx
//   1/(Rx - Qx) = w^{-1} · dx
//
// We don't know Rx - Qx yet (need λ), so we can't compute w directly
// as `dx · (Rx - Qx)`. We compute it as `u / dx = (dy² - dx³ - 2Qx·dx²)/dx
//                                      = dy²·dx^{-1} - dx² - 2Qx·dx`.
// Dividing by dx needs 1/dx. Circular.
//
// Alternative: just take `w' = dx · u = dx·dy² - dx⁴ - 2Qx·dx³`. Then
// w'^{-1} = 1/(dx·u). Product of dx and u, inverted once. Then:
//   1/dx = u · w'^{-1}
//   1/u  = dx · w'^{-1}
// We get 1/dx from one inversion of a degree-4 polynomial in dx plus
// degree-2 in dy. That's two quantum muls (u, w' = dx·u) then inversion
// then two more muls to recover 1/dx. Five muls total before λ — more
// expensive than just inverting dx directly.
//
// But we ALSO get 1/u for free. And 1/u tells us λ² via:
//   u = dy² - dx²(Px+Qx)  =>  dy² = u + dx²(Px+Qx)
//   λ² = dy²/dx² = (u + dx²(Px+Qx))/dx² = u/dx² + Px + Qx
//   λ² - Px - Qx = u/dx².
//   So Rx = u/dx² and Rx - Qx = u/dx² - Qx. Hmm that's u·dx^{-2}.
//
// Simpler cleaner algebra:
//   Rx = λ² - Px - Qx
//   dx² · Rx = dy² - dx²·Px - dx²·Qx = dy² - dx²·(Px + Qx)
//   So Rx = (dy² - dx²·(Px+Qx)) / dx².
//   Let v = dy² - dx²·(Px+Qx). Then Rx = v/dx².
//   And for Ry:
//     Ry = λ(Qx - Rx) - Qy.
//     dx·Ry = dy·(Qx - Rx) - dx·Qy
//     dx³·Ry = dx²·dy·(Qx - Rx) - dx³·Qy
//            = dy·(dx²·Qx - dx²·Rx) - dx³·Qy
//            = dy·(dx²·Qx - v) - dx³·Qy
//   So Ry = (dy·(dx²·Qx - v) - dx³·Qy) / dx³.
//   Let w = dx³. Then:
//     Rx = v · w^{-1} · dx                 (since v/dx² = v·dx/dx³)
//     Ry = (dy·(dx²·Qx - v) - w·Qy) · w^{-1}
//
//   One inversion (of w = dx³), both outputs computable.
//
// This is the cleanest Strategy C candidate: invert dx³, recover
// (Rx, Ry) via muls on known-live quantum values. Counts later.
// ─────────────────────────────────────────────────────────────────────
pub fn replay_strategy_c(px: U256, py: U256, qx: U256, qy: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    let dx = sub_mod(px, qx, p);
    let dy = sub_mod(py, qy, p);

    // v = dy² - dx²·(Px + Qx)
    let dx2 = dx.mul_mod(dx, p);
    let dx3 = dx2.mul_mod(dx, p);
    let dy2 = dy.mul_mod(dy, p);
    let px_plus_qx = px.add_mod(qx, p);
    let v = sub_mod(dy2, dx2.mul_mod(px_plus_qx, p), p);

    // one inversion
    let w = dx3;
    let w_inv = w.inv_mod(p).expect("dx != 0");

    // Rx = v · w^{-1} · dx = v · (dx · w^{-1})
    let dx_winv = dx.mul_mod(w_inv, p);
    let rx = v.mul_mod(dx_winv, p);

    // Ry = (dy·(dx²·Qx - v) - w·Qy) · w^{-1}
    let dx2_qx = dx2.mul_mod(qx, p);
    let core = sub_mod(dx2_qx, v, p);
    let numer = sub_mod(dy.mul_mod(core, p), w.mul_mod(qy, p), p);
    let ry = numer.mul_mod(w_inv, p);

    (rx, ry)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

    fn curve() -> WeierstrassEllipticCurve {
        WeierstrassEllipticCurve {
            modulus: SECP256K1_P,
            a: U256::from(0),
            b: U256::from(7),
            gx: U256::from_str_radix(
                "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
                16,
            )
            .unwrap(),
            gy: U256::from_str_radix(
                "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
                16,
            )
            .unwrap(),
            order: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
                16,
            )
            .unwrap(),
        }
    }

    fn rand_u256(rng: &mut u64) -> U256 {
        let mut limbs = [0u64; 4];
        for l in &mut limbs {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *l = *rng;
        }
        U256::from_limbs(limbs) % SECP256K1_P
    }

    fn each_trial<F>(mut body: F)
    where
        F: FnMut(U256, U256, U256, U256, U256, U256),
    {
        let c = curve();
        let mut rng = 0xdead_beef_cafe_f00du64;
        let mut n = 0usize;
        let mut tried = 0usize;
        while n < 200 && tried < 2000 {
            tried += 1;
            let k1 = rand_u256(&mut rng);
            let k2 = rand_u256(&mut rng);
            let (px, py) = c.mul(c.gx, c.gy, k1);
            let (qx, qy) = c.mul(c.gx, c.gy, k2);
            if (px.is_zero() && py.is_zero())
                || (qx.is_zero() && qy.is_zero())
                || px == qx
            {
                continue;
            }
            let (rx_ref, ry_ref) = c.add(px, py, qx, qy);
            body(px, py, qx, qy, rx_ref, ry_ref);
            n += 1;
        }
        assert_eq!(n, 200, "needed 200 random valid trials");
    }

    #[test]
    fn reference_formula_sanity() {
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = single_inv_add(px, py, qx, qy);
            assert_eq!(rx, rx_ref);
            assert_eq!(ry, ry_ref);
        });
    }

    #[test]
    fn strategy_a_rx_ok_ry_contaminated_by_dy() {
        // Strategy A is predicted DEAD. Specifically:
        //   Rx: should match reference (200/200).
        //   Ry: should equal `ref_Ry + dy` on every trial (the exact
        //       contamination predicted in the plan doc:
        //       ty_final = Py + Ry - Qy = Ry + dy).
        let mut rx_ok = 0;
        let mut ry_off_by_dy = 0;
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = replay_strategy_a(px, py, qx, qy);
            if rx == rx_ref {
                rx_ok += 1;
            }
            let dy = sub_mod(py, qy, SECP256K1_P);
            let ry_expected_with_bug = ry_ref.add_mod(dy, SECP256K1_P);
            if ry == ry_expected_with_bug {
                ry_off_by_dy += 1;
            }
        });
        eprintln!("Strategy A: rx ok {rx_ok}/200 ; ry == ref+dy in {ry_off_by_dy}/200");
        assert_eq!(rx_ok, 200, "Strategy A Rx must match all 200");
        assert_eq!(
            ry_off_by_dy, 200,
            "Strategy A Ry should be off by exactly +dy on every trial"
        );
    }

    #[test]
    fn strategy_b2_passes_but_leaks_lam_copy() {
        // Strategy B2 with the lam_copy ancilla trick: Rx and Ry both
        // pass, BUT lam_copy (= -λ) is a leaked ancilla we don't yet
        // know how to reversibly zero without a second inversion.
        // Cost accounting below ignores that leak.
        let mut rx_ok = 0;
        let mut ry_ok = 0;
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = replay_strategy_b2(px, py, qx, qy);
            if rx == rx_ref {
                rx_ok += 1;
            }
            if ry == ry_ref {
                ry_ok += 1;
            }
        });
        eprintln!(
            "Strategy B2: rx matches {rx_ok}/200, ry matches {ry_ok}/200 (lam_copy leaked)"
        );
        assert_eq!(rx_ok, 200);
        assert_eq!(ry_ok, 200);
    }

    /// Falsification test for B2: at Kaliski body exit we have the
    /// registers {tx=Rx (or Rx-Qx), ty=Ry, dx_orig=dx, inv_raw=-dx⁻¹·2^{2n-1}}
    /// plus the classical constants ox=Qx, oy=Qy.
    /// Can -λ be expressed as a small-cost polynomial combination of these?
    ///
    /// Enumerate a catalogue of candidate expressions (each one corresponds
    /// to a tiny mul/add sequence, ~1–3 muls of cost ~150k each). If ANY
    /// candidate equals -λ on 200 random trials, B2 is alive. If NONE do,
    /// B2 is dead.
    #[test]
    fn strategy_b2_lam_copy_uncompute_falsification() {
        let p = SECP256K1_P;
        // Candidate functions: take (rx, ry, qx, qy, dx, dx_inv) → guess for -λ.
        type Cand = fn(U256, U256, U256, U256, U256, U256, U256) -> U256;
        let neg = |a: U256| sub_mod(U256::ZERO, a, p);
        let _ = neg;

        // Polynomial in {rx, ry, qx, qy, dx, dx_inv, dx_inv^2}. dx_inv^2 is
        // "free" relative to adding a new inversion: it's one extra mul.
        let candidates: &[(&str, Cand)] = &[
            ("(Ry+Qy)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                ry.add_mod(qy, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("-(Ry+Qy)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                sub_mod(U256::ZERO, ry.add_mod(qy, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P), SECP256K1_P)
            }),
            ("(Qy-Ry)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                sub_mod(qy, ry, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("(Ry-Qy)·dx_inv", |_rx, ry, _qx, qy, _dx, dx_inv, _| {
                sub_mod(ry, qy, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("(Rx-Qx)·dx_inv", |rx, _ry, qx, _qy, _dx, dx_inv, _| {
                sub_mod(rx, qx, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
            ("dx_inv_sq", |_rx, _ry, _qx, _qy, _dx, _dx_inv, dx_inv_sq| dx_inv_sq),
            ("(Ry+Qy)·dx_inv_sq", |_rx, ry, _qx, qy, _dx, _dx_inv, dx_inv_sq| {
                ry.add_mod(qy, SECP256K1_P).mul_mod(dx_inv_sq, SECP256K1_P)
            }),
            ("(Rx-Qx)·dx_inv_sq", |rx, _ry, qx, _qy, _dx, _dx_inv, dx_inv_sq| {
                sub_mod(rx, qx, SECP256K1_P).mul_mod(dx_inv_sq, SECP256K1_P)
            }),
            ("dx·dx_inv_sq", |_rx, _ry, _qx, _qy, dx, _dx_inv, dx_inv_sq| {
                dx.mul_mod(dx_inv_sq, SECP256K1_P)
            }),
            ("(Rx+Qx)·dx_inv", |rx, _ry, qx, _qy, _dx, dx_inv, _| {
                rx.add_mod(qx, SECP256K1_P).mul_mod(dx_inv, SECP256K1_P)
            }),
        ];

        let mut hits = vec![0usize; candidates.len()];
        let mut total = 0usize;
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            total += 1;
            let dx = sub_mod(px, qx, p);
            let dy = sub_mod(py, qy, p);
            let lam = dy.mul_mod(dx.inv_mod(p).unwrap(), p);
            let neg_lam = sub_mod(U256::ZERO, lam, p);
            let dx_inv = dx.inv_mod(p).unwrap();
            let dx_inv_sq = dx_inv.mul_mod(dx_inv, p);
            for (i, (_name, f)) in candidates.iter().enumerate() {
                let got = f(rx_ref, ry_ref, qx, qy, dx, dx_inv, dx_inv_sq);
                if got == neg_lam {
                    hits[i] += 1;
                }
            }
        });
        let mut any_matched_all = false;
        for (i, (name, _)) in candidates.iter().enumerate() {
            eprintln!("  candidate {name}: {}/{} matches -λ", hits[i], total);
            if hits[i] == total {
                any_matched_all = true;
            }
        }
        assert!(
            !any_matched_all,
            "B2 rescued: a low-cost polynomial expression in {{Rx, Ry, Qx, Qy, dx, dx⁻¹, dx⁻²}} \
             matches -λ on all trials; re-examine the obstruction"
        );
        eprintln!("B2 falsification: no low-cost polynomial in {{Rx, Ry, Qx, Qy, dx, dx⁻¹, dx⁻²}} equals -λ. B2 DEAD.");
    }

    #[test]
    fn strategy_c_passes_200() {
        // Strategy C: invert w = dx³, recover both Rx and Ry from it.
        // Everything is classical-reversible — only question is cost.
        each_trial(|px, py, qx, qy, rx_ref, ry_ref| {
            let (rx, ry) = replay_strategy_c(px, py, qx, qy);
            assert_eq!(rx, rx_ref);
            assert_eq!(ry, ry_ref);
        });
    }
}
