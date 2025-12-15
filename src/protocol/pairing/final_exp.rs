use mathlib::U1024;

use crate::algebra::fields::fp6::Fp6;
use crate::traits::field::Field;

pub fn final_exponentiation<'a>(f: &Fp6<'a>, exponent: &U1024) -> Fp6<'a> {
    // Perform f ^ exponent using Square-and-Multiply
    // Note: f is in Fp6.

    let params = f.c0.c0.params;
    // xi generation removed as it's not used in standard mul

    let mut res = Fp6::one(params);
    let mut base = *f;

    let num_bits = exponent.bits();

    for i in 0..num_bits {
        if exponent.bit(i) {
            res = res * base;
        }
        base = base * base;
    }

    res
}
