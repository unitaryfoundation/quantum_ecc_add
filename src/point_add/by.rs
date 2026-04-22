//! Bernstein-Yang divsteps2 : classical test harness and (later) quantum
//! implementation.
//!
//! Ref: Bernstein & Yang 2019, "Fast constant-time gcd computation and
//! modular inversion" (TCHES 2019(3)).  https://gcd.cr.yp.to/safegcd-20190413.pdf
//!
//! ## divstep (δ, f, g)
//!
//! ```text
//! if δ > 0 and g is odd:   (1 − δ, g, (g − f) / 2)
//! elif         g is odd:   (1 + δ, f, (g + f) / 2)
//! else:                    (1 + δ, f, g / 2)
//! ```
//!
//! ## Invariants
//!
//! - `f` always odd.
//! - `|f|, |g| ≤ max(|f₀|, |g₀|)` throughout.
//! - After `N ≥ safegcd(n)` iters with `gcd(f₀, g₀) = 1` and `n`-bit inputs,
//!   `f_N = ±1` and `g_N = 0`.
//!
//! ## Coefficient tracking (uniform `2^k` scaling)
//!
//! Track integers (`U`, `V`, `Q`, `R`) satisfying
//!
//! ```text
//! 2^k · f_k = U · f₀ + V · g₀
//! 2^k · g_k = Q · f₀ + R · g₀
//! ```
//!
//! Per-iteration updates (so both trackers gain a factor of 2 every step):
//!
//! - Case A (δ > 0 ∧ g odd): `(U, V, Q, R) ← (2Q, 2R, Q − U, R − V)`, `δ ← 1 − δ`.
//! - Case B (g odd, δ ≤ 0):  `(U, V, Q, R) ← (2U, 2V, Q + U, R + V)`, `δ ← 1 + δ`.
//! - Case C (g even):        `(U, V, Q, R) ← (2U, 2V, Q,     R    )`, `δ ← 1 + δ`.
//!
//! Recovery of `value^{-1}  mod p` from `(f₀, g₀) = (p, value)`:
//! at termination with `f_N ∈ {±1}`, `g_N = 0`, taking mod `p`:
//!
//! ```text
//! 2^N · f_N ≡ V_N · value  (mod p)
//! value^{-1} ≡ sign(f_N) · V_N · 2^{−N}  (mod p)
//! ```
//!
//! Safegcd iteration bound: `N_n = ⌈(49n + 80) / 17⌉`.
//! For `n = 256` bits, `N_256 = 743`.

pub fn safegcd_iters(n_bits: usize) -> usize {
    // ceil((49 * n + 80) / 17)
    (49 * n_bits + 80 + 16) / 17
}

/// Classical one-step-at-a-time (w = 1) divsteps2.
///
/// Tracks `(f, g)` as signed `i128` (caller: ensure inputs fit signed-127),
/// and `(U, V, Q, R)` as `u128` reduced mod `p`. Parity decisions read the
/// low bit of `g` as a signed integer.
///
/// Returns `(delta_final, f_final, g_final, U, V, Q, R)` with coefficients
/// in `[0, p)`.
pub fn classical_divsteps2_i128(
    n_iters: usize,
    delta_init: i64,
    f_init: i128,
    g_init: i128,
    p: u128,
) -> (i64, i128, i128, u128, u128, u128, u128) {
    assert!(f_init & 1 == 1, "f must be odd");
    assert!(p > 2 && p & 1 == 1, "p must be an odd modulus");
    assert!(p < (1u128 << 127), "p must fit so a+b doesn't wrap in u128");
    let addm = |a: u128, b: u128| -> u128 {
        let s = a + b;
        if s >= p { s - p } else { s }
    };
    let subm = |a: u128, b: u128| -> u128 {
        if a >= b { a - b } else { p - (b - a) }
    };

    let mut delta = delta_init;
    let mut f = f_init;
    let mut g = g_init;
    let mut uu: u128 = 1;
    let mut vv: u128 = 0;
    let mut qq: u128 = 0;
    let mut rr: u128 = 1;

    for _ in 0..n_iters {
        let g_odd = (g & 1) != 0;
        if delta > 0 && g_odd {
            // Case A: (f, g) ← (g, (g − f) / 2), δ ← 1 − δ.
            let nf = g;
            let ng = (g - f) >> 1; // g − f is even (odd − odd).
            let nu = addm(qq, qq);
            let nv = addm(rr, rr);
            let nq = subm(qq, uu);
            let nr = subm(rr, vv);
            delta = 1 - delta;
            f = nf; g = ng;
            uu = nu; vv = nv; qq = nq; rr = nr;
        } else if g_odd {
            // Case B: (f, g) ← (f, (g + f) / 2), δ ← 1 + δ.
            let ng = (g + f) >> 1;
            let nu = addm(uu, uu);
            let nv = addm(vv, vv);
            let nq = addm(qq, uu);
            let nr = addm(rr, vv);
            delta = 1 + delta;
            g = ng;
            uu = nu; vv = nv; qq = nq; rr = nr;
        } else {
            // Case C: (f, g) ← (f, g / 2), δ ← 1 + δ.
            let ng = g >> 1;
            let nu = addm(uu, uu);
            let nv = addm(vv, vv);
            // Q, R unchanged.
            delta = 1 + delta;
            g = ng;
            uu = nu; vv = nv;
        }
    }
    (delta, f, g, uu, vv, qq, rr)
}

pub fn pow_mod_u128(mut base: u128, mut exp: u128, p: u128) -> u128 {
    assert!(p < (1u128 << 63), "pow_mod_u128 requires p < 2^63 to avoid mul overflow");
    base %= p;
    let mut r: u128 = 1 % p;
    while exp > 0 {
        if exp & 1 == 1 { r = (r * base) % p; }
        exp >>= 1;
        if exp > 0 { base = (base * base) % p; }
    }
    r
}

pub fn gcd_u128(a: u128, b: u128) -> u128 {
    if b == 0 { a } else { gcd_u128(b, a % b) }
}

/// Modular inverse via classical B-Y.
///
/// For `gcd(value, p) == 1`, returns `Some(value^{-1} mod p)`.
/// Otherwise returns `None`.
pub fn classical_by_modinv_i128(value: u128, p: u128) -> Option<u128> {
    if value == 0 { return None; }
    let bits = 128 - p.leading_zeros() as usize;
    let n_iters = safegcd_iters(bits);
    let (_, f_final, g_final, _uu, vv, _qq, _rr) =
        classical_divsteps2_i128(n_iters, 1, p as i128, value as i128, p);
    if g_final != 0 { return None; }

    // value^{-1} ≡ sign(f_final) · V · 2^{-N}  (mod p)
    // Compute 2^{-1} mod p = 2^{p-2} mod p (Fermat), then raise to N.
    let two_inv = pow_mod_u128(2, p - 2, p);
    let two_inv_n = pow_mod_u128(two_inv, n_iters as u128, p);
    let v_scaled = (vv * two_inv_n) % p;

    match f_final {
        1  => Some(v_scaled),
        -1 => Some(if v_scaled == 0 { 0 } else { p - v_scaled }),
        _  => None,
    }
}

/// Env-gated: exhaustive verification against Fermat modinv over small primes.
/// Run with `BY_TEST=1`.
pub fn run_classical_test() {
    let primes: &[u128] = &[
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
        71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
        149, 151, 251, 257, 509, 1009, 65537, 1_000_003, 2_147_483_647,
        (1u128 << 61) - 1, // Mersenne M61
    ];
    let mut total: u64 = 0;
    let mut pass: u64 = 0;
    let mut first_fail: Option<(u128, u128, Option<u128>, u128)> = None;
    for &p in primes {
        let bound = p.min(400);
        for val in 1..bound {
            if gcd_u128(val, p) != 1 { continue; }
            total += 1;
            let expected = pow_mod_u128(val, p - 2, p);
            let got = classical_by_modinv_i128(val, p);
            if got == Some(expected) {
                pass += 1;
            } else if first_fail.is_none() {
                first_fail = Some((p, val, got, expected));
            }
        }
    }
    eprintln!("BY classical w=1: {}/{} pass", pass, total);
    if let Some((p, val, got, expected)) = first_fail {
        eprintln!("  first fail: p={} val={} got={:?} expected={}", p, val, got, expected);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_primes_exhaustive() {
        // Exhaustively verify BY modinv against Fermat for p up to 257.
        let primes: &[u128] = &[
            3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
            67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
            137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
            199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
        ];
        for &p in primes {
            for val in 1..p {
                if gcd_u128(val, p) != 1 { continue; }
                let expected = pow_mod_u128(val, p - 2, p);
                let got = classical_by_modinv_i128(val, p);
                assert_eq!(got, Some(expected),
                    "p={} val={}: BY got {:?}, Fermat got {}", p, val, got, expected);
            }
        }
    }

    #[test]
    fn larger_primes_spot_check() {
        // Spot-check 1000 values for each of a handful of larger primes.
        let primes: &[u128] = &[
            1009, 65537, 1_000_003, 2_147_483_647,
            (1u128 << 61) - 1, // Mersenne M61
        ];
        for &p in primes {
            let step = (p / 1000).max(1);
            let mut tested = 0u32;
            let mut val: u128 = 1;
            while tested < 1000 && val < p {
                if gcd_u128(val, p) == 1 {
                    let expected = pow_mod_u128(val, p - 2, p);
                    let got = classical_by_modinv_i128(val, p);
                    assert_eq!(got, Some(expected),
                        "p={} val={}: BY got {:?}, Fermat got {}", p, val, got, expected);
                    tested += 1;
                }
                val = val.wrapping_add(step);
            }
        }
    }

    #[test]
    fn convergence_check() {
        // After safegcd_iters bits for coprime (f₀, g₀), expect f = ±1, g = 0.
        let p: u128 = 2_147_483_647; // Mersenne M31
        let bits = 128 - p.leading_zeros() as usize;
        let n = safegcd_iters(bits);
        for val in [1u128, 2, 3, 7, 17, 1000, 123_456, p - 1, p / 2].iter() {
            let (_, f, g, _, _, _, _) =
                classical_divsteps2_i128(n, 1, p as i128, *val as i128, p);
            assert_eq!(g, 0, "p={} val={}: g did not converge to 0 (got {})", p, val, g);
            assert!(f == 1 || f == -1, "p={} val={}: f = {} (expected ±1)", p, val, f);
        }
    }
}
