pub mod final_exp;
pub mod miller;

use mathlib::U1024;

use crate::{
    algebra::fields::Fp6,
    models::{sextic_twist::STPoint as G2Projective, short_weierstrass::SWPoint as G1Affine},
};

/// Computes the Tate pairing of a G1 affine point and a G2 projective point.
///
/// The result is the pairing value in the Fp6 target field computed by first
/// running the Miller loop and then applying the final exponentiation.
///
/// # Examples
///
/// ```ignore
/// use mathlib::U1024;
/// use crate::models::{sextic_twist::STPoint as G2Projective, short_weierstrass::SWPoint as G1Affine};
/// use crate::protocol::pairing::tate_pairing;
///
/// let p: G1Affine = /* construct or obtain a G1 affine point */;
/// let q: G2Projective = /* construct or obtain a G2 projective point */;
/// let r_order: U1024 = /* curve subgroup order */;
/// let final_exp_val: U1024 = /* final exponentiation exponent */;
///
/// let result = tate_pairing(&p, &q, r_order, final_exp_val);
/// // `result` is an `Fp6` element representing the pairing value.
/// ```
pub fn tate_pairing<'a>(
    p: &G1Affine<'a>,
    q: &G2Projective<'a>,
    r_order: U1024,
    final_exp_val: U1024,
) -> Fp6<'a> {
    // 1. Miller Loop
    let f = miller::miller_loop(p, q, r_order);

    // 2. Final Exponentiation
    final_exp::final_exponentiation(&f, &final_exp_val)
}
