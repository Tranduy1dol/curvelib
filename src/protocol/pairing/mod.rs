pub mod final_exp;
pub mod miller;

use mathlib::U1024;

use crate::algebra::fields::fp6::Fp6;
use crate::models::sextic_twist::G2Point;
use crate::models::short_weierstrass::SWPoint;

pub fn tate_pairing<'a>(
    p: &SWPoint<'a>,
    q: &G2Point<'a>,
    r_order: U1024,
    final_exp_val: U1024,
) -> Fp6<'a> {
    // 1. Miller Loop
    let f = miller::miller_loop(p, q, r_order);

    // 2. Final Exponentiation
    final_exp::final_exponentiation(&f, &final_exp_val)
}
