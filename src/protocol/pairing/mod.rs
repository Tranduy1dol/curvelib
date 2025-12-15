pub mod final_exp;
pub mod miller;

use mathlib::U1024;

use crate::{
    algebra::fields::Fp6,
    models::{sextic_twist::STPoint as G2Projective, short_weierstrass::SWPoint as G1Affine},
};

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
